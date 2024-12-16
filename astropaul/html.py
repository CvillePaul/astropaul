import collections
from datetime import datetime
import os
import platform

import astropy.units as u
import dominate
import dominate.tags as tags
import dominate.util as util
import itables
import pandas as pd

import astropaul.priority as pr
import astropaul.targetlistcreator as tlc


def dataframe_to_datatable(
    table: pd.DataFrame, 
    table_name: str = "table", 
    caption: str = "", 
    table_options: dict = None, 
    buttons: list = None
):
    default_options = {
        "connected": True,
        "paging": False,
        "maxBytes": 0,
        "maxColumns": 0,
        "autoWidth": False,
        "layout": {"topStart": None, "topEnd": None, "bottomStart": None, "bottomEnd": None},
        "classes": "compact cell-border hover",
    }
    if buttons:
        default_options["buttons"] = buttons
    if table_options is None:
        table_options = {}
    if caption == "":
        caption = f"Created {datetime.now().astimezone().isoformat(sep=" ", timespec="seconds")} on {platform.node()}"
    table_name = table_name.replace(" ", "_")
    html = itables.to_html_datatable(df=table, caption=caption, table_id=table_name, **{**default_options, **table_options})
    html += f"""
        <style>
        td {{text-align: center}}
        #{table_name} th {{
            white-space: normal;
            word-wrap: break-word;
            text-align: center;
        }}
        caption {{
            text-align: left;
            font-style: italic;
        }}
        </style>
        """
    return html


def render_observing_pages(tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], dir: str = "html") -> str:
    # wipe out contents of dir
    # for file in os.listdir(dir):
    #     file_path = os.path.join(dir, file)
    #     os.unlink(file_path)

    # make overall summary page that links to all the other pages
    with dominate.document(title="Observing Files") as d:
        border = 1
        padding = 10
        d += tags.h1(tl.name)
        # output information about the observing session
        d += tags.h2("Observing session details:")
        t = tags.table(border=border, cellpadding=padding)
        with t:
            tags.tr(
                tags.td("Site"),
                tags.td(pl.session.site_info),
            )
            tags.tr(
                tags.td("Total Time"),
                tags.td(f"{pl.session.total_time.to(u.hour).value:.1f} hours"),
            )
            time_range = pl.session.time_range
            tags.tr(
                tags.td("Observing Segments"),
                tags.td(len(pl.segments)),
                # tags.td(util.raw(f"{len(pl.segments)} total<br>{time_range[0].iso[:19]} start<br>{time_range[1].iso[:19]} end")),
            )
            tags.tr(
                tags.td("UTC Range"),
                tags.td(f"{time_range[0].iso[:19]} through {time_range[1].iso[:19]}"),
            )
        d += t
        d += tags.br()
        # output information about the targets
        d += tags.h2(f"Target List ({len(tl.target_list)} targets)", tags.a("Link to table", href=f"{tl.name}.html"))
        target_types = collections.Counter(tl.target_list["Target Type"])
        t = tags.table(border=border, cellpadding=padding)
        with t:
            for type, count in target_types.items():
                tags.tr(
                    tags.td(count),
                    tags.td(type),
                )
        d += t
        d += tags.br()
        # output information about other files
        if other_files:
            d += tags.h2("Other Files")
            for file_type, contents in other_files.items():
                with dominate.document(title=file_type) as od:
                    od += tags.pre(contents)
                    with open(f"{dir}/{file_type}.html", "w") as f:
                        f.write(od.render())
                d += tags.span(f"{file_type}: ", tags.a("Contents", href=f"{file_type}.html"))
                d += tags.br()
            d += tags.br()
                
        # output information about target priorities
        d += tags.h2(f"Priorities ({pl.interval} interval)")
        t = tags.table(border=border, cellpadding=padding)
        with t:
            tags.tr(tags.th("Start UTC"), tags.th("Finish UTC"), tags.th(), tags.th())
            for subsegments in pl.segments:
                beg = subsegments[0][0]
                end = subsegments[-1][1]
                tags.tr(
                    tags.td(beg.iso[:19]),
                    tags.td(end.iso[:19]),
                    tags.td(tags.a("Numerical", href=f"Numerical Priorities {beg.iso[:10]}.html")),
                    tags.td(tags.a("Categorical", href=f"Categorical Priorities {beg.iso[:10]}.html")),
                )
        d += t
        with open(f"{dir}/index.html", "w") as f:
            f.write(d.render())

    # make page for target list
    tltl = tl.target_list.copy()
    tltl["Target Name"] = [f'<a href="targets/{target_name}.html">{target_name}</a>' for target_name in tltl["Target Name"]]
    title = tl.name
    with dominate.document(title=title) as d:
        d += tags.h1(title, style="text-align: center")
        buttons = [
            {
                "extend": "columnToggle",
                "columns": [f"{col}:title" for col in secondary_cols],
                "text": group,
                "redrawCalculations": True,
            }
            for group, (_, secondary_cols) in tl.column_groups.items()
            if secondary_cols
        ]
        d += util.raw(
            dataframe_to_datatable(
                tltl,
                title,
                table_options={
                    "showIndex": False,
                    "search": True,
                    "layout": {
                        "topStart": {
                            "searchBuilder": True,
                            "buttons": buttons,
                        }
                    },
                },
            )
        )
        with open(f"{dir}/{title}.html", "w") as f:
            f.write(d.render())

    # make a page for each target's priority scores
    for target, target_table_list in pl.target_tables.items():
        for target_table in target_table_list:
            tt = target_table.copy()
            start_utc = f"{tt.index[0]:%Y-%m-%d}"
            tt.index = [f"{time:%H:%M}" for time in tt.index]
            title = f"{start_utc} Target Scores for {target}"
            with dominate.document(title=title) as d:
                d += tags.h1("Target Scores for", tags.a(target, href=f"../targets/{target}.html"), style="text-align: center")
                d += tags.p(util.raw(dataframe_to_datatable(tt, "Target_Scores", table_options={"sort": False})))
                with open(f"{dir}/target scores/Target Scores {target} {start_utc}.html", "w") as f:
                    f.write(d.render())

    # make pages for numerical priorities
    if pl.numerical_priorities:
        for priority_table in pl.numerical_priorities:
            pt = priority_table.copy()
            start_utc = f"{pt.index[0]:%Y-%m-%d}"
            pt.index = [f"{time:%H:%M}" for time in pt.index]
            threshold = pl.category_bins[-2]
            for col in pt.columns:
                pt[col] = [
                    val if val < threshold else f'<span style="background-color: #EBF4FA">{val:.3f}</span>' for val in pt[col]
                ]
            pt.columns = [f'<a href="target scores/Target Scores {target} {start_utc}.html">{target}</a>' for target in pt.columns]
            title = f"{start_utc} Numerical Priorities"
            with dominate.document(title=title) as d:
                d += tags.h1(title, style="text-align: center")
                d += util.raw(dataframe_to_datatable(pt, "Numerical_Priority", table_options={"sort": False}))
                with open(f"{dir}/Numerical Priorities {start_utc}.html", "w") as f:
                    f.write(d.render())

    # make pages for categorical priorities
    if pl.categorical_priorities:
        for categories_table in pl.categorical_priorities:
            ct = categories_table.copy()  # make a copy we can alter for formatting purposes
            start_utc = f"{ct.index[0]:%Y-%m-%d}"
            ct.index = [f"{time:%H:%M}" for time in ct.index]
            highlight_value = "* * *"
            for col in ct.columns:
                ct[col] = [
                    val if val != highlight_value else f'<span style="background-color: #EBF4FA">{val}</span>'
                    for val in ct[col]
                ]
            list_targets = tl.target_list
            # not all targets on list have nonzero priority, so we need a list of targets for this particular segment
            segment_targets = list_targets[list_targets["Target Name"].isin(ct.columns)]
            # add some helpful informational rows to the top of the chart
            ct.loc["RA"] = [segment_targets[segment_targets["Target Name"] == col]["RA"].values[0][:8] for col in ct.columns]
            ct.loc["Dec"] = [segment_targets[segment_targets["Target Name"] == col]["Dec"].values[0][:9] for col in ct.columns]
            ct.loc["Teff"] = [
                f"{segment_targets[segment_targets["Target Name"] == col]["Teff"].values[0]:.0f}" for col in ct.columns
            ]
            ct.loc["RV Standard"] = [segment_targets[segment_targets["Target Name"] == col]["RV Standard"].values[0] for col in ct.columns]
            new_ordering = [*ct.index[-4:], *ct.index[:-4]]
            ct = ct.reindex(new_ordering)
            # elide the target names to first 4 digits
            ct.columns = [col[:8] + "..." for col in ct.columns]
            # generate the chart
            title = f"{start_utc} Categorical Priorities"
            with dominate.document(title=title) as d:
                d += tags.h1(title, style="text-align: center")
                d += util.raw(dataframe_to_datatable(ct, "Categorical_Priority", table_options={"sort": False}))
                with open(f"{dir}/Categorical Priorities {start_utc}.html", "w") as f:
                    f.write(d.render())

    # make pages for each individual target
    columns = ["RA", "Dec", "Vmag", "Teff", "ra", "dec"]
    for _, row in tl.target_list.iterrows():
        target_name = row["Target Name"]
        with dominate.document(title=target_name) as d:
            d += tags.h1(target_name, style="text-align: center")
            with tags.table(border=border, cellpadding=padding) as t:
                tags.tr([tags.th(item) for item in columns])
                tags.tr([tags.td(value) for value in row[columns]])
                d += t
            other_tables = [
                table_name
                for table_name, table in tl.other_lists.items()
                if "Target Name" in list(tl.other_lists[table_name].columns) + [tl.other_lists[table_name].index.name]
            ]
            table_styles = [
                {"selector": "td", "props": [("padding", "10px"), ("border", "1px solid black"), ("text-align", "center")]}
            ]
            for other_table in other_tables:
                d += tags.h2(f"{other_table} Entries")
                ot = tl.other_lists[other_table]
                if ot.index.name == "Target Name":
                    if target_name in ot.index:
                        entries = ot.loc[[target_name]].reset_index(drop=True)
                    else:
                        entries = pd.DataFrame()
                else:
                    entries = ot[ot["Target Name"] == target_name].drop("Target Name", axis=1)
                if not entries.empty:
                    d += util.raw(entries.style.hide(axis="index").set_table_styles(table_styles).to_html())
            with open(f"{dir}/targets/{target_name}.html", "w") as f:
                f.write(d.render())
