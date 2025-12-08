import collections
from datetime import datetime
from pathlib import Path
import platform
import shutil
import textwrap

import astropy.units as u
from astropy.time import Time
import dominate
import dominate.tags as tags
import dominate.util as util
import itables
import pandas as pd

from astropaul.database import html_path, resources_path
import astropaul.priority as pr
import astropaul.targetlistcreator as tlc


def default_table_css(table_id: str) -> str:
    return textwrap.dedent(
        f"""
        <style>
        #{table_id} th {{
            white-space: normal;
            word-wrap: break-word;
            text-align: center;
        }}
        #{table_id} td {{
            text-align: center;
            word-wrap: normal;
        }}
        caption {{
            text-align: left;
            font-style: italic;
        }}
        button.exportButton,
        a.exportButton {{
            padding: 0px 4px !important;
            font-size: 12px !important;
        }}
        </style>
        """
    )


def dataframe_to_datatable(
    table: pd.DataFrame,
    table_id: str = "table",
    table_options: dict = None,
    buttons: list = None,
    column_defs: dict = None,
    table_css: str = None,
):
    default_options = {
        "connected": True,
        "paging": False,
        "maxBytes": 0,
        "maxColumns": 0,
        "autoWidth": True,
        "style": "width: auto; float: left; caption-side:bottom;",
        "layout": {"topStart": None, "topEnd": None, "bottomStart": None, "bottomEnd": None},
        "classes": "compact cell-border hover",
    }
    if buttons:
        default_options["buttons"] = buttons
    if not table_options:
        table_options = {}
    if not column_defs:
        column_defs = {}
    table_id = table_id.replace(" ", "_")
    html = itables.to_html_datatable(
        df=table,
        table_id=table_id,
        **{**default_options, **table_options},
        columnDefs=column_defs,
    )
    if table_css:
        html += table_css
    else:
        html += default_table_css(table_id)
    return html


def horizontal_space():
    return tags.span(style="display: inline-block; width: 20px;")


def keybinding_script():
    return textwrap.dedent(
        """
        document.addEventListener("keydown", function(event) {
            const el = event.target;
            if (
                el.tagName === 'INPUT' ||
                el.tagName === 'TEXTAREA' ||
                el.isContentEditable
            ) {
                return;  // Ignore keys typed into text fields or editable elements
            }
            if (event.key === "ArrowRight") { document.getElementById("nextDay").click(); }
            if (event.key === "ArrowLeft")  { document.getElementById("prevDay").click(); }
            if (event.key === "s")          { document.getElementById("summary").click(); }
            if (event.key === "c")          { document.getElementById("categorical").click(); }
            if (event.key === "n")          { document.getElementById("numerical").click(); }
            if (event.key === "t")          { document.getElementById("targets").click(); }
            if (event.key === "j")          { document.getElementById("nextTarget").click(); }
            if (event.key === "k")          { document.getElementById("prevTarget").click(); }
            if (event.key === "h")          { document.getElementById("home").click(); }
        });
    """
    )


def standard_footer():
    return tags.p(
        f"Created {datetime.now().astimezone().isoformat(sep=" ", timespec="seconds")} on {platform.node()}",
        style="text-align: left; font-style: italic;",
    )


def cell_border():
    return 1


def cell_padding():
    return 10


def make_summary_page(tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], dir: str = "html") -> None:
    with dominate.document(title=tl.name) as d:
        d.head += tags.a("Home", href="../index.html", id="home")
        d += tags.h1(tl.name)
        if pl:
            # output information about the observing session
            d += tags.h2("Observing session details:")
            t = tags.table(border=cell_border(), cellpadding=cell_padding())
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
                tags.tr(
                    tags.td("JD Range"),
                    tags.td(f"{time_range[0].jd:.4f} through {time_range[1].jd:.4f}"),
                )
            d += t
        d += tags.br()
        # output information about the targets
        d += tags.h2(
            f"Target List ({len(tl.target_list)} targets)",
            tags.a("View target table", href=f"{tl.name}.html", id="targets"),
        )
        target_types = collections.Counter(tl.target_list["Target Type"])
        t = tags.table(border=cell_border(), cellpadding=cell_padding())
        with t:
            for type, count in target_types.items():
                tags.tr(
                    tags.td(count),
                    tags.td(type),
                )
        d += t
        d += tags.br()
        if tl.list_criteria and len(tl.list_criteria) > 0:
            d += tags.h2("List Criteria")
            criteria = tags.ul()
            for line in str(tl.list_criteria).split("\n"):
                criteria += tags.li(tags.pre(line))
            d += criteria

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

        if pl:
            # output information about target priorities
            d += tags.h2(f"Priorities ({pl.interval} interval)")
            t = tags.table(border=cell_border(), cellpadding=cell_padding())
            t.set_attribute("style", "text-align: center;")
            with t:
                columns = [
                    tags.th("Start UTC/JD"),
                    tags.th("Finish UTC/JD"),
                    tags.th("Observing Time"),
                    tags.th("Moon Illumination"),
                    tags.th("Priorities"),
                    tags.th("Priorities"),
                ]
                for label in pl.category_labels[::-1]:
                    columns.append(tags.th(f"Category: {label if label else "(None)"}"))
                tags.tr(*columns)
                first_iteration = True
                for subsegments, illumination, members in zip(pl.segments, pl.moon_illumination, pl.segment_category_members):
                    beg = subsegments[0][0]
                    end = subsegments[-1][1]
                    if first_iteration:
                        numerical_id = "numerical"
                        categorical_id = "categorical"
                        first_iteration = False
                    else:
                        numerical_id = "numerical_other"
                        categorical_id = "categorical_other"
                    beg_cell = tags.td(beg.iso[:19])
                    beg_cell.add(tags.br())
                    beg_cell.add(f"{beg.jd:.4f}")
                    end_cell = tags.td(end.iso[:19])
                    end_cell.add(tags.br())
                    end_cell.add(f"{end.jd:.4f}")
                    columns = [
                        beg_cell,
                        end_cell,
                        tags.td(f"{(end - beg).to(u.hour):.1f}"),
                        tags.td(f"{illumination:.2f}"),
                        tags.td(tags.a("Numerical", href=f"Numerical Priorities {beg.iso[:10]}.html", id=numerical_id)),
                        tags.td(tags.a("Categorical", href=f"Categorical Priorities {beg.iso[:10]}.html", id=categorical_id)),
                    ]
                    for label in pl.category_labels[::-1]:
                        label_members = members[label]
                        content = tags.span(len(label_members))
                        content.set_attribute("title", ", ".join(label_members))
                        columns.append(tags.td(content))
                    tags.tr(*columns)
            d += t
            # write list of members of each category
            d += tags.h2("Category Members")
            d += tags.p("Highest category attained by a target across entire observing session")
            t = tags.table(border=cell_border(), cellpadding=cell_padding())
            with t:
                tags.tr(tags.th("Category"), tags.th("Count"), tags.th("Members"))
                for label in pl.category_labels[::-1]:
                    members = sorted(pl.overall_category_members[label])
                    tags.tr(tags.td(label if label else "(None)"), tags.td(len(members)), tags.td(", ".join(members)))
            d += t
        # write out summary page
        d += tags.script(util.raw(keybinding_script()))
        d.footer += standard_footer()

        with open(f"{dir}/index.html", "w") as f:
            f.write(d.render())


def make_target_list(tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], dir: str = "html") -> None:
    tltl = tl.target_list.copy()
    tltl["Target Name"] = [
        f'<a href="../Targets/targets/{target_name}.html">{target_name}</a>' for target_name in tltl["Target Name"]
    ]
    title = tl.name
    with dominate.document(title=title) as d:
        d.head += tags.a("Summary Page", href="index.html", id="summary")
        d += tags.h1(title, style="text-align: center")
        d.head += horizontal_space()
        d.head += tags.a("Home", href="../index.html", id="home")
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
        col_indexes_to_hide = [
            tltl.columns.get_loc(col)
            for _, secondary_cols in tl.column_groups.values()
            if secondary_cols
            for col in secondary_cols
        ]
        column_defs = [{"targets": col_indexes_to_hide, "visible": False}]
        d += util.raw(
            dataframe_to_datatable(
                tltl,
                table_options={
                    "showIndex": False,
                    "layout": {
                        "topStart": {
                            "search": True,
                            "buttons": buttons,
                        },
                        "topEnd": {
                            "info": True,
                            "searchBuilder": True,
                        },
                        "bottomStart": {
                            "buttons": [
                                {"extend": "csvHtml5", "className": "exportButton"},
                                {"extend": "copyHtml5", "className": "exportButton"},
                            ],
                        },
                    },
                    "autoWidth": False,
                },
                column_defs=column_defs,
            )
        )
        d += tags.script(util.raw(keybinding_script()))
        d.footer += standard_footer()
        with open(f"{dir}/{title}.html", "w") as f:
            f.write(d.render())


def make_numerical_scores_pages(
    tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], dir: str = "html"
) -> None:
    # make pages for numerical priorities
    if pl and pl.numerical_priorities:
        start_times = [f"{pt.index[0]:%Y-%m-%d}" for pt in pl.numerical_priorities]
        threshold = pl.category_bins[-2]
        above_threshold_style = "background-color: #EBF4FA;"
        not_observable_style = "color: #DDDDDD;"
        for i, priority_table in enumerate(pl.numerical_priorities):
            pt = priority_table.copy()
            start_utc = start_times[i]
            pt.index = [f"{time:%H:%M}" for time in pt.index]
            # apply background highlighting to cells above threshold
            for col in pt.columns:
                target_scores = pl.target_tables[col][i]
                if "Altitude Priority" in target_scores.columns:
                    altitudes = target_scores["Altitude Priority"].values
                else:
                    altitudes = [1] * len(target_scores)
                styled_scores = []
                for score, altitude in zip(pt[col], altitudes):
                    style = "text-align: center;"
                    precision = 0 if score == 0 else 3
                    if score > threshold:
                        style += above_threshold_style
                    if altitude == 0:
                        style += not_observable_style
                    styled_scores.append(f'<span style="{style}">{score:.{precision}f}</span>')
                pt[col] = styled_scores
            # turn column headings into links to target score pages
            all_targets = pt.columns.to_list()  # save list of target names, in column order, before changing them into links
            pt = pt.transpose()  # flip so targets are rows and columns are timeslots
            target_columns = [
                "Target Type",
                "Target Source",
                "RA",
                "Dec",
                "RA HMS",
                "Dec DMS",
                "Vmag",
                "Teff",
                "Num Ephemerides",
                "Num DSSI Observations",
                "Num Speckle Detections",
                "Num PEPSI Observations",
            ]
            pt = tl.target_list[["Target Name"] + target_columns].merge(
                pt, left_on="Target Name", right_index=True, how="right"
            )
            pt.reset_index(drop=True, inplace=True)
            first_target = pt.iloc[0]["Target Name"]
            pt["Target Name"] = [  # add target info columns to allow filtering of rows
                (
                    f'<div style=text-align:left><a href="target scores/Target Scores {target} {start_utc}.html", id="targets">{target}</a></div>'
                    if target == first_target
                    else f'<div style=text-align:left><a href="target scores/Target Scores {target} {start_utc}.html">{target}</a></div>'
                )
                for target in pt["Target Name"]
            ]
            col_indexes_to_hide = [pt.columns.get_loc(col) for col in target_columns]
            buttons = [
                {
                    "extend": "columnToggle",
                    "columns": [f"{col}:title" for col in target_columns],
                    "text": "Toggle Target Info",
                    "redrawCalculations": True,
                }
            ]
            pt.columns = [""] + list(pt.columns[1:])
            title = f"Numerical Priorities, {start_utc} UTC"
            with dominate.document(title=title) as d:
                d += tags.h1(title, style="text-align: center")
                table_options = {
                    "sort": False,
                    "layout": {
                        "topStart": {
                            "buttons": buttons,
                        },
                        "topEnd": {
                            "info": True,
                            "searchBuilder": True,
                        },
                        "bottomStart": {
                            "buttons": [
                                {"extend": "csvHtml5", "className": "exportButton"},
                                {"extend": "copyHtml5", "className": "exportButton"},
                            ],
                        },
                    },
                }
                column_defs = [{"targets": col_indexes_to_hide, "visible": False}]
                d += util.raw(
                    dataframe_to_datatable(pt, "Numerical_Priority", table_options=table_options, column_defs=column_defs)
                )
                d += util.raw("<style>#Numerical_Priority th {writing-mode:sideways-lr;}</style>")

                if i > 0:
                    d.head += tags.a("<-Prev", href=f"Numerical Priorities {start_times[i - 1]}.html", id="prevDay")
                else:
                    d.head += tags.span("<-Prev")
                d.head += horizontal_space()
                if i < len(start_times) - 1:
                    d.head += tags.a("Next->", href=f"Numerical Priorities {start_times[i + 1]}.html", id="nextDay")
                else:
                    d.head += tags.span("Next->")
                d.head += horizontal_space()
                d.head += tags.a("Categorical Priorities", href=f"Categorical Priorities {start_utc}.html", id="categorical")
                d.head += horizontal_space()
                d.head += tags.a("Summary Page", href="index.html", id="summary")
                d.head += horizontal_space()
                d.head += tags.a("Home", href="../index.html", id="home")
                d += tags.script(util.raw(keybinding_script()))
                d.footer += standard_footer()
                with open(f"{dir}/Numerical Priorities {start_utc}.html", "w") as f:
                    f.write(d.render())

            # now write target scores pages for each column of the numerical priorities table
            for target in all_targets:
                tt = pl.target_tables[target][i].copy()
                tt.index = [f"{time:%H:%M}" for time in tt.index]
                title = f"{start_utc} Target Scores for {target}"
                table_options = {
                    "sort": False,
                    "layout": {
                        "bottomStart": {
                            "buttons": [
                                {"extend": "csvHtml5", "className": "exportButton"},
                                {"extend": "copyHtml5", "className": "exportButton"},
                            ],
                        },
                    },
                }
                with dominate.document(title=title) as d:
                    d += tags.h1(
                        tags.a(target, href=f"../../Targets/targets/{target}.html", id="targets"),
                        tags.span(style="display: inline-block; width: 10px;"),
                        f"Target Scores for {start_utc} UTC",
                        style="text-align: center",
                    )
                    d += tags.p(util.raw(dataframe_to_datatable(tt, "Target_Scores", table_options=table_options)))
                    if i > 0:
                        d.head += tags.a("<-Prev Day", href=f"Target Scores {target} {start_times[i - 1]}.html", id="prevDay")
                    else:
                        d.head += tags.span("<-Prev Day")
                    d.head += horizontal_space()
                    if i < len(start_times) - 1:
                        d.head += tags.a("Next Day->", href=f"Target Scores {target} {start_times[i + 1]}.html", id="nextDay")
                    else:
                        d.head += tags.span("Next Day->")
                    d.head += horizontal_space()
                    target_index = all_targets.index(target)
                    if target_index == 0:
                        d.head += tags.span("<-Prev Target")
                    else:
                        d.head += tags.a(
                            "<-Prev Target",
                            href=f"Target Scores {all_targets[target_index - 1]} {start_utc}.html",
                            id="prevTarget",
                        )
                    d.head += horizontal_space()
                    if target_index == len(all_targets) - 1:
                        d.head += tags.span("Next Target->")
                    else:
                        d.head += tags.a(
                            "Next Target->",
                            href=f"Target Scores {all_targets[target_index + 1]} {start_utc}.html",
                            id="nextTarget",
                        )
                    d.head += horizontal_space()
                    d.head += tags.a("Numerical Priorities", href=f"../Numerical Priorities {start_utc}.html", id="numerical")
                    d.head += horizontal_space()
                    d.head += tags.a(
                        "Categorical Priorities", href=f"../Categorical Priorities {start_utc}.html", id="categorical"
                    )
                    d.head += horizontal_space()
                    d.head += tags.a("Summary Page", href="../index.html", id="summary")
                    d.head += horizontal_space()
                    d.head += tags.a("Home", href="../../index.html", id="home")
                    d += tags.script(util.raw(keybinding_script()))
                    with open(f"{dir}/target scores/Target Scores {target} {start_utc}.html", "w") as f:
                        f.write(d.render())


def make_categorical_scores(tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], dir: str = "html") -> None:
    # make pages for categorical priorities
    if pl and pl.categorical_priorities:
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
            ct.loc["RA"] = [
                segment_targets[segment_targets["Target Name"] == col]["RA HMS"].values[0][:8] for col in ct.columns
            ]
            ct.loc["Dec"] = [
                segment_targets[segment_targets["Target Name"] == col]["Dec DMS"].values[0][:9] for col in ct.columns
            ]
            ct.loc["Vmag"] = [
                f"{segment_targets[segment_targets["Target Name"] == col]["Vmag"].values[0]:.1f}" for col in ct.columns
            ]
            num_added_rows = 3
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
            if "Num SIDE" in segment_targets.columns:
                ct.loc["Num SIDE"] = [
                    segment_targets[segment_targets["Target Name"] == col]["Num SIDE"].values[0] for col in ct.columns
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
                d.head += horizontal_space()
                if i < len(start_times) - 1:
                    d.head += tags.a("Next->", href=f"Categorical Priorities {start_times[i + 1]}.html", id="nextDay")
                else:
                    d.head += tags.span("Next->")
                d.head += horizontal_space()
                d.head += tags.a("Numerical Priorities", href=f"Numerical Priorities {start_utc}.html", id="numerical")
                d.head += horizontal_space()
                d.head += tags.a("Summary Page", href="index.html", id="summary")
                d += tags.script(util.raw(keybinding_script()))
                d.head += horizontal_space()
                d.head += tags.a("Home", href="../index.html", id="home")
                with open(f"{dir}/Categorical Priorities {start_utc}.html", "w") as f:
                    f.write(d.render())


def make_target_pages(tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], dir: str = "html") -> None:
    columns_to_skip = {}
    # do special things if certain resources are present
    pepsi_table = "PEPSI Observations"
    if pepsi_table in tl.other_lists:
        pepsi_observations = tl.other_lists[pepsi_table]
        columns_to_skip[pepsi_table] = ["PEPSI SpectrumPlot"]
        spectra_dir = Path(dir) / "spectrum_plots"
        spectra_dir.mkdir(exist_ok=True)
        cd_stats = collections.defaultdict(set)
        for (target_name, cd), rows in pepsi_observations.sort_values("Mid JD").groupby(["Target Name", "Cross Disperser"]):
            name = f"{pepsi_table} for {target_name}"
            with dominate.document(name=name) as d:
                d.title = name
                d.head += tags.link(rel="stylesheet", href="https://unpkg.com/flickity@2/dist/flickity.min.css")
                # d.head += tags.link(rel="stylesheet", href="https://unpkg.com/flickity-fullscreen@2/fullscreen.css")
                d.head += tags.script(src="https://unpkg.com/flickity@2/dist/flickity.pkgd.min.js")
                # d.head += tags.script(src="https://unpkg.com/flickity-fullscreen@2/fullscreen.js")
                carousel = tags.div(cls="carousel")
                for _, row in rows.iterrows():
                    cell = tags.div(cls="flikity-cell") 
                    cell += tags.img(src=Path(row["PEPSI SpectrumPlot"]).name)
                    carousel += cell
                d += carousel
                d += util.raw("</div>")
                d += tags.style(".flickity-cell {width: 100%}")
                d += tags.script("""
                    var flkty = new Flickity('.carousel', {
                        wrapAround: false,
                        pageDots: false,
                    });

                    flkty.on( 'change', function( index ) {
                        flkty.select(index, false, true);
                    }); 
                    flkty.focus();
                    """)
                with open(f"{spectra_dir}/{target_name} - CD{cd[0]}.html", "w") as f:
                    f.write(d.render())
            cd_stats[target_name].add(cd[0])

    list_memberships_name = "List Memberships"
    other_tables = {
        table_name: table.copy()
        for table_name, table in tl.other_lists.items()
        if "Target Name" in list(table.columns) + [table.index.name]
        and table_name != list_memberships_name  # skip this because it's handled elsewhere
        # and len(table[table["Target Name"] == target_name]) > 0 # skip empty tables
    }
    # if spectral observations are present, add: links to plot in table, page w/ carousel of all plots
    # TODO: this code permanently alters the observations table, it should instead modify a copy of it, leaving original intact
    resources_dir = resources_path()
    spectroscopy_observations = [("PEPSI Observations", "PEPSI SpectrumPlot")] # (table, plot column name)
    for table_name, column_name in spectroscopy_observations:
        if table_name in other_tables:
            # add the link in the table
            table = other_tables[table_name]
            if not "Plot" in table.columns:
                for plot_resource in table[column_name]:
                    shutil.copy(resources_dir / plot_resource, spectra_dir)
                table["Plot"] = [
                    f'<a href="file:///C:{spectra_dir / Path(file).name}" target="_blank">Plot</a>'
                    for file in table[column_name]
                ]
                # observations.drop(column_name, axis=1, inplace=True)                    

    # make pages for each individual target
    Path(dir / "targets").mkdir(exist_ok=True)
    for _, row in tl.target_list.iterrows():
        target_name = row["Target Name"]
        with dominate.document(title=target_name) as d:
            d += tags.h1(f"{target_name} Target Details", style="text-align: center")
            with tags.table(cellpadding=cell_padding()) as t:
                tags.tr([tags.td("Target Type"), tags.td(row["Target Type"])])
                tags.tr([tags.td("Target Source"), tags.td(row["Target Source"])])
                # treat List Memberships specially - wrap contents as comma separated string instead of dataframe
                if list_memberships_name in tl.other_lists.keys():
                    list_memberships = tl.other_lists[list_memberships_name]
                    lists = list_memberships[list_memberships["Target Name"] == target_name]["List"].to_list()
                    if lists:
                        tags.tr([tags.td(list_memberships_name), tags.td(", ".join([list[5:] for list in lists]))])
                d += t
            d += tags.div(style="margin-bottom: 40px;")
            with tags.table(border=cell_border(), cellpadding=cell_padding()) as t:
                desired_columns = ["RA", "Dec", "RA HMS", "Dec DMS", "Vmag", "Teff", "Distance"]
                columns = [item for item in desired_columns if item in tl.target_list.columns]
                tags.tr([tags.th(item) for item in columns])
                tags.tr([tags.td(value) for value in row[columns]])
                d += t

            # add all the associated tables as grids
            for table_name, other_table in other_tables.items():
                for column_name in columns_to_skip.get(table_name, []):
                    other_table.drop(column_name, axis=1, inplace=True)
                if other_table.index.name == "Target Name":
                    if target_name in ot.index:
                        entries = other_table.loc[[target_name]].reset_index(drop=True)
                    else:
                        entries = pd.DataFrame()
                else:
                    entries = other_table[other_table["Target Name"] == target_name].drop("Target Name", axis=1)
                entries = entries.reset_index(drop=True)
                if not entries.empty:
                    d += tags.h2(f"{table_name} ({len(entries)})")
                    if table_name == pepsi_table and target_name in cd_stats:
                        for cd in sorted(cd_stats[target_name]):
                            d += tags.a(f"CD{cd} Plots",href=f"../{spectra_dir.name}/{target_name} - CD{cd}.html", target="_blank")
                            d += tags.br()
                    table_options = {
                        "autowidth": False,
                        "float": "left",
                        "layout": {
                            "topEnd": None,
                            "bottomStart": {
                                "buttons": [
                                    {"extend": "csvHtml5", "className": "exportButton"},
                                    {"extend": "copyHtml5", "className": "exportButton"},
                                ],
                            },
                        },
                    }
                    if "ID" in entries.columns:
                        entries.drop("ID", axis=1, inplace=True)
                    d += util.raw(dataframe_to_datatable(entries, table_name, table_options))
                # else:
                #     d += tags.span("(Empty Table)")

            d.head += horizontal_space()
            d.head += tags.a("Summary Page", href="../index.html", id="summary")
            d.head += horizontal_space()
            d.head += tags.a("Home", href="../../index.html", id="home")
            d += tags.script(util.raw(keybinding_script()))
            d.footer += standard_footer()
            with open(f"{dir}/targets/{target_name}.html", "w") as f:
                f.write(d.render())


def clear_directory(dir: str) -> None:
    # wipe out contents of dir
    dir_cleared = False
    fail_count = 0
    while not dir_cleared:
        try:
            shutil.rmtree(dir)
            dir_cleared = True
        except FileNotFoundError:
            break  # don't care if directory didn't already exist
        except PermissionError as e:
            # try again many times if permission error while removing directory (due to Dropbox file sync lock issue)
            fail_count += 1
            if fail_count >= 50:
                raise PermissionError(f"Tried {fail_count} times", e)


def render_observing_pages(tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], subdir: str = "html") -> str:
    dir = html_path() / subdir
    clear_directory(dir)
    Path(dir).mkdir()
    if pl:
        Path(dir / "target scores").mkdir(parents=True)

    make_summary_page(tl, pl, other_files, dir)
    make_target_list(tl, pl, other_files, dir)
    make_numerical_scores_pages(tl, pl, other_files, dir)
    make_categorical_scores(tl, pl, other_files, dir)
    if subdir == "targets":
        make_target_pages(tl, pl, other_files, dir)


from playwright.sync_api import sync_playwright

# import asyncio


def html_to_pdf(input_html_path, output_pdf_path):
    with open(input_html_path, "r") as f:
        html_content = f.read()
    beg = html_content.find("<head>")
    end = html_content.find("</head>")
    if beg > 0 and end > 0 and end > beg:
        html_content = html_content[0:beg] + html_content[end + 7 :]
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        # Letter size in pixels at 96 DPI: 8.5 x 11 inches → 816 x 1056
        page.set_viewport_size({"width": 816, "height": 1056})
        page.set_content(html_content)
        page.evaluate(
            """
            const table = document.querySelector('table.dataTable'); 
            if (table) {
                table.style.fontSize = 'small';
            }
        """
        )
        page.pdf(
            path=output_pdf_path,
            scale=0.8,
            format="Letter",
            margin={"left": "0.4 in"},
            landscape=True,
            print_background=True,
        )
        browser.close()
