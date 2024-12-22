import collections
from datetime import datetime
import pathlib
import platform
import shutil
import textwrap

import astropy.units as u
import dominate
import dominate.tags as tags
import dominate.util as util
import itables
import pandas as pd

import astropaul.priority as pr
import astropaul.targetlistcreator as tlc


def dataframe_to_datatable(
    table: pd.DataFrame, table_name: str = "table", caption: str = "", table_options: dict = None, buttons: list = None
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
    html += textwrap.dedent(
        f"""
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
    )
    return html


def render_observing_pages(tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], dir: str = "html") -> str:
    # wipe out contents of dir
    try:
        shutil.rmtree(dir)
    except FileNotFoundError:
        pass  # don't care if directory didn't already exist
    pathlib.Path(f"{dir}/targets").mkdir(parents=True)
    pathlib.Path(f"{dir}/target scores").mkdir(parents=True)

    horizontal_space = tags.span(style="display: inline-block; width: 20px;")
    keybinding_script = textwrap.dedent(
        """
        document.addEventListener("keydown", function(event) {
            if (event.key === "ArrowRight") { document.getElementById("nextDay").click(); }
            if (event.key === "ArrowLeft")  { document.getElementById("prevDay").click(); }
            if (event.key === "s")          { document.getElementById("summary").click(); }
            if (event.key === "c")          { document.getElementById("categorical").click(); }
            if (event.key === "n")          { document.getElementById("numerical").click(); }
            if (event.key === "t")          { document.getElementById("targets").click(); }
            if (event.key === "j")          { document.getElementById("nextTarget").click(); }
            if (event.key === "k")          { document.getElementById("prevTarget").click(); }
        });
    """
    )
    # make overall summary page that links to all the other pages
    with dominate.document(title=tl.name) as d:
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
                tags.td("Available Time"),
                tags.td(f"{pl.session.total_time.to(u.hour).value:.1f} hours in {len(pl.segments)} observing segments"),
            )
            time_range = pl.session.time_range
            tags.tr(
                tags.td("UTC Range"),
                tags.td(f"{time_range[0].iso[:19]} through {time_range[1].iso[:19]}"),
            )
        d += t
        d += tags.br()
        # output information about the targets
        d += tags.h2(
            f"Target List ({len(tl.target_list)} targets)", tags.a("View target table", href=f"{tl.name}.html", id="targets")
        )
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
            first_iteration = True
            for subsegments in pl.segments:
                beg = subsegments[0][0]
                end = subsegments[-1][1]
                if first_iteration:
                    numerical_id = "numerical"
                    categorical_id = "categorical"
                    first_iteration = False
                else:
                    numerical_id = "numerical_other"
                    categorical_id = "categorical_other"
                tags.tr(
                    tags.td(beg.iso[:19]),
                    tags.td(end.iso[:19]),
                    tags.td(tags.a("Numerical", href=f"Numerical Priorities {beg.iso[:10]}.html", id=numerical_id)),
                    tags.td(tags.a("Categorical", href=f"Categorical Priorities {beg.iso[:10]}.html", id=categorical_id)),
                )
        d += t
        d += tags.script(util.raw(keybinding_script))
        with open(f"{dir}/index.html", "w") as f:
            f.write(d.render())

    # make page for target list
    tltl = tl.target_list.copy()
    tltl["Target Name"] = [f'<a href="targets/{target_name}.html">{target_name}</a>' for target_name in tltl["Target Name"]]
    title = tl.name
    with dominate.document(title=title) as d:
        d.head += tags.a("Summary Page", href="index.html", id="summary")
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
        d += tags.script(util.raw(keybinding_script))
        with open(f"{dir}/{title}.html", "w") as f:
            f.write(d.render())

    # make pages for numerical priorities
    if pl.numerical_priorities:
        start_times = [f"{pt.index[0]:%Y-%m-%d}" for pt in pl.numerical_priorities]
        for i, priority_table in enumerate(pl.numerical_priorities):
            pt = priority_table.copy()
            start_utc = start_times[i]
            pt.index = [f"{time:%H:%M}" for time in pt.index]
            # apply background highlighting to cells above threshold
            threshold = pl.category_bins[-2]
            for col in pt.columns:
                pt[col] = [
                    (
                        0
                        if val == 0
                        else f"{val:.3f}" if val < threshold else f'<span style="background-color: #EBF4FA">{val:.3f}</span>'
                    )
                    for val in pt[col]
                ]
            # turn column headings into links to target score pages
            all_targets = pt.columns.to_list() # save list of target names, in column order, before changing them into links
            pt.columns = [
                (
                    f'<a href="target scores/Target Scores {target} {start_utc}.html", id="targets">{target}</a>'
                    if target == pt.columns[0]
                    else f'<a href="target scores/Target Scores {target} {start_utc}.html">{target}</a>'
                )
                for target in pt.columns
            ]
            title = f"Numerical Priorities, {start_utc} UTC"
            with dominate.document(title=title) as d:
                d += tags.h1(title, style="text-align: center")
                d += util.raw(dataframe_to_datatable(pt, "Numerical_Priority", table_options={"sort": False}))
                if i > 0:
                    d.head += tags.a("<-Prev", href=f"Numerical Priorities {start_times[i - 1]}.html", id="prevDay")
                else:
                    d.head += tags.span("<-Prev")
                d.head += horizontal_space
                if i < len(start_times) - 1:
                    d.head += tags.a("Next->", href=f"Numerical Priorities {start_times[i + 1]}.html", id="nextDay")
                else:
                    d.head += tags.span("Next->")
                d.head += horizontal_space
                d.head += tags.a("Categorical Priorities", href=f"Categorical Priorities {start_utc}.html", id="categorical")
                d.head += horizontal_space
                d.head += tags.a("Summary Page", href="index.html", id="summary")
                d += tags.script(util.raw(keybinding_script))
                with open(f"{dir}/Numerical Priorities {start_utc}.html", "w") as f:
                    f.write(d.render())

                # now write target scores pages for each column of the numerical priorities table
                for target in all_targets:
                    tt = pl.target_tables[target][i].copy()
                    tt.index = [f"{time:%H:%M}" for time in tt.index]
                    title = f"{start_utc} Target Scores for {target}"
                    with dominate.document(title=title) as d:
                        d += tags.h1(
                            tags.a(target, href=f"../targets/{target}.html", id="targets"),
                            tags.span(style="display: inline-block; width: 10px;"),
                            f"Target Scores for {start_utc} UTC",
                            style="text-align: center",
                        )
                        d += tags.p(util.raw(dataframe_to_datatable(tt, "Target_Scores", table_options={"sort": False})))
                        if i > 0:
                            d.head += tags.a("<-Prev Day", href=f"Target Scores {target} {start_times[i - 1]}.html", id="prevDay")
                        else:
                            d.head += tags.span("<-Prev Day")
                        d.head += horizontal_space
                        if i < len(start_times) - 1:
                            d.head += tags.a("Next Day->", href=f"Target Scores {target} {start_times[i + 1]}.html", id="nextDay")
                        else:
                            d.head += tags.span("Next Day->")
                        d.head += horizontal_space
                        target_index = all_targets.index(target)
                        if target_index == 0:
                            d.head += tags.span("<-Prev Target")
                        else:
                            d.head += tags.a(
                                "<-Prev Target", href=f"Target Scores {all_targets[target_index - 1]} {start_utc}.html", id="prevTarget"
                            )
                        d.head += horizontal_space
                        if target_index == len(all_targets) - 1:
                            d.head += tags.span("Next Target->")
                        else:
                            d.head += tags.a(
                                "Next Target->", href=f"Target Scores {all_targets[target_index + 1]} {start_utc}.html", id="nextTarget"
                            )
                        d.head += horizontal_space
                        d.head += tags.a("Numerical Priorities", href=f"../Numerical Priorities {start_utc}.html", id="numerical")
                        d.head += horizontal_space
                        d.head += tags.a("Categorical Priorities", href=f"../Categorical Priorities {start_utc}.html", id="categorical")
                        d.head += horizontal_space
                        d.head += tags.a("Summary Page", href="../index.html", id="summary")
                        d += tags.script(util.raw(keybinding_script))
                        with open(f"{dir}/target scores/Target Scores {target} {start_utc}.html", "w") as f:
                            f.write(d.render())

    # make pages for categorical priorities
    if pl.categorical_priorities:
        start_times = [f"{pt.index[0]:%Y-%m-%d}" for pt in pl.categorical_priorities]
        for i, categories_table in enumerate(pl.categorical_priorities):
            ct = categories_table.copy()  # make a copy we can alter for formatting purposes
            if len(ct.columns) == 0:
                continue  # skip cases where there were no targets passing criteria
            start_utc = start_times[i]
            ct.index = [f"{time:%H:%M}" for time in ct.index]
            highlight_value = pl.category_labels[-1]
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
            num_added_rows = 2
            if "PEPSI exp_time" in segment_targets.columns:
                ct.loc["Teff"] = [
                    f"{segment_targets[segment_targets["Target Name"] == col]["Teff"].values[0]:.0f}" for col in ct.columns
                ]
                num_added_rows += 1
            if "RV Standard" in segment_targets.columns:
                ct.loc["RV Standard"] = [
                    segment_targets[segment_targets["Target Name"] == col]["RV Standard"].values[0] for col in ct.columns
                ]
                num_added_rows += 1
            new_ordering = [*ct.index[-num_added_rows:], *ct.index[:-num_added_rows]]
            ct = ct.reindex(new_ordering)
            # elide the target names to first 4 digits
            ct.columns = [col[:8] + "..." for col in ct.columns]
            # generate the chart
            title = f"Categorical Priorities, {start_utc} UTC"
            with dominate.document(title=title) as d:
                d += tags.h1(title, style="text-align: center")
                d += util.raw(dataframe_to_datatable(ct, "Categorical_Priority", table_options={"sort": False}))
                if i > 0:
                    d.head += tags.a("<-Prev", href=f"Categorical Priorities {start_times[i - 1]}.html", id="prevDay")
                else:
                    d.head += tags.span("<-Prev")
                d.head += horizontal_space
                if i < len(start_times) - 1:
                    d.head += tags.a("Next->", href=f"Categorical Priorities {start_times[i + 1]}.html", id="nextDay")
                else:
                    d.head += tags.span("Next->")
                d.head += horizontal_space
                d.head += tags.a("Numerical Priorities", href=f"Numerical Priorities {start_utc}.html", id="numerical")
                d.head += horizontal_space
                d.head += tags.a("Summary Page", href="index.html", id="summary")
                d += tags.script(util.raw(keybinding_script))
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
                ot = ot.sort_values(["Target Name"] + ot.columns[0:3].to_list())
                if ot.index.name == "Target Name":
                    if target_name in ot.index:
                        entries = ot.loc[[target_name]].reset_index(drop=True)
                    else:
                        entries = pd.DataFrame()
                else:
                    entries = ot[ot["Target Name"] == target_name].drop("Target Name", axis=1)
                if not entries.empty:
                    d += util.raw(entries.style.hide(axis="index").set_table_styles(table_styles).to_html())
            d.head += horizontal_space
            d.head += tags.a("Summary Page", href="../index.html", id="summary")
            d += tags.script(util.raw(keybinding_script))
            d.footer += tags.p(
                f"Created {datetime.now().astimezone().isoformat(sep=" ", timespec="seconds")} on {platform.node()}",
                style="text-align: left; font-style: italic;",
            )
            with open(f"{dir}/targets/{target_name}.html", "w") as f:
                f.write(d.render())
