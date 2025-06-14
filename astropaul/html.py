import collections
from datetime import datetime
import pathlib
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

import astropaul.priority as pr
import astropaul.targetlistcreator as tlc


def dataframe_to_datatable(
    table: pd.DataFrame, table_id: str = "table", table_options: dict = None, buttons: list = None, column_defs: dict = None
):
    default_options = {
        "connected": True,
        "paging": False,
        # "maxBytes": 0,
        # "maxColumns": 0,
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
        df=table, table_id=table_id, **{**default_options, **table_options}, columnDefs=column_defs
    )
    html += textwrap.dedent(
        f"""
        <style>
        #{table_id} th {{
            white-space: normal;
            word-wrap: break-word;
            text-align: center;
        }}
        #{{table_id}} td {{
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
    return html


def horizontal_space():
    return tags.span(style="display: inline-block; width: 20px;")


def keybinding_script():
    return textwrap.dedent(
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
            tags.a("View target table", href=f"{tl.name} Target List.html", id="targets"),
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
            for criterion in tl.list_criteria:
                criteria += tags.li(tags.code(criterion))
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
    tltl["Target Name"] = [f'<a href="targets/{target_name}.html">{target_name}</a>' for target_name in tltl["Target Name"]]
    title = f"{tl.name} Target List"
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
                table_options = {
                    "sort": False,
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
                d += util.raw(dataframe_to_datatable(pt, "Numerical_Priority", table_options=table_options))
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
                d += tags.script(util.raw(keybinding_script()))
                d.footer += standard_footer()
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
                        d += tags.p(util.raw(dataframe_to_datatable(tt, "Target_Scores", table_options=table_options)))
                        if i > 0:
                            d.head += tags.a(
                                "<-Prev Day", href=f"Target Scores {target} {start_times[i - 1]}.html", id="prevDay"
                            )
                        else:
                            d.head += tags.span("<-Prev Day")
                        d.head += horizontal_space()
                        if i < len(start_times) - 1:
                            d.head += tags.a(
                                "Next Day->", href=f"Target Scores {target} {start_times[i + 1]}.html", id="nextDay"
                            )
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
                        d.head += tags.a(
                            "Numerical Priorities", href=f"../Numerical Priorities {start_utc}.html", id="numerical"
                        )
                        d.head += horizontal_space()
                        d.head += tags.a(
                            "Categorical Priorities", href=f"../Categorical Priorities {start_utc}.html", id="categorical"
                        )
                        d.head += horizontal_space()
                        d.head += tags.a("Summary Page", href="../index.html", id="summary")
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
                with open(f"{dir}/Categorical Priorities {start_utc}.html", "w") as f:
                    f.write(d.render())


def make_target_pages(tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], dir: str = "html") -> None:
    # make pages for each individual target
    desired_columns = ["RA", "Dec", "RA HMS", "Dec DMS", "Vmag", "Teff", "Distance"]
    columns = [item for item in desired_columns if item in tl.target_list.columns]
    
    for _, row in tl.target_list.iterrows():
        target_name = row["Target Name"]
        with dominate.document(title=target_name) as d:
            d += tags.h1(f"{target_name} Target Details", style="text-align: center")
            with tags.table(border=cell_border(), cellpadding=cell_padding()) as t:
                tags.tr([tags.th(item) for item in columns])
                tags.tr([tags.td(value) for value in row[columns]])
                d += t
            other_tables = [
                table_name
                for table_name, table in tl.other_lists.items()
                if "Target Name" in list(tl.other_lists[table_name].columns) + [tl.other_lists[table_name].index.name]
            ]
            for other_table in other_tables:
                ot = tl.other_lists[other_table]
                # ot = ot.sort_values(["Target Name"] + ot.columns[0:3].to_list())
                if ot.index.name == "Target Name":
                    if target_name in ot.index:
                        entries = ot.loc[[target_name]].reset_index(drop=True)
                    else:
                        entries = pd.DataFrame()
                else:
                    entries = ot[ot["Target Name"] == target_name].drop("Target Name", axis=1)
                d += tags.h2(f"{other_table} ({len(entries)})")
                entries = entries.reset_index(drop=True)
                if not entries.empty:
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
                    d += util.raw(dataframe_to_datatable(entries, other_table, table_options))
                else:
                    d += tags.span("(Empty Table)")

            d.head += horizontal_space()
            d.head += tags.a("Summary Page", href="../index.html", id="summary")
            d += tags.script(util.raw(keybinding_script()))
            d.footer += standard_footer()
            with open(f"{dir}/targets/{target_name}.html", "w") as f:
                f.write(d.render())


def render_observing_pages(tl: tlc.TargetList, pl: pr.PriorityList, other_files: dict[str, str], dir: str = "html") -> str:
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
            if fail_count > 50:
                raise PermissionError(f"Tried {fail_count} times", e)
    pathlib.Path(f"{dir}/targets").mkdir(parents=True)
    if pl:
        pathlib.Path(f"{dir}/target scores").mkdir(parents=True)

    make_summary_page(tl, pl, other_files, dir)
    make_target_list(tl, pl, other_files, dir)
    make_numerical_scores_pages(tl, pl, other_files, dir)
    make_categorical_scores(tl, pl, other_files, dir)
    make_target_pages(tl, pl, other_files, dir)


from playwright.sync_api import sync_playwright


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
        page.set_content(html_content)
        page.evaluate(
            """
            const table = document.querySelector('table.dataTable'); 
            if (table) {
                table.style.width = 'auto';
                table.style.maxWidth = '100%'; 
                table.style.fontSize = 'smaller';
            }
        """
        )
        page.pdf(
            path=output_pdf_path,
            format="Letter",
            landscape=True,
            print_background=True,
        )
        browser.close()
