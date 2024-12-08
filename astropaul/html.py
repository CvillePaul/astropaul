from datetime import datetime
import platform

import itables
import pandas as pd


def dataframe_to_datatable(table: pd.DataFrame, table_name: str = "table", caption: str = "", table_options: dict = None):
    default_options = {
        "connected": True,
        "paging": False,
        "maxBytes": 0,
        "maxColumns": 0,
        "autoWidth": False,
        "layout": {"topStart": None, "topEnd": None, "bottomStart": None, "bottomEnd": None},
        "classes": "compact cell-border hover",
    }
    if table_options is None:
        table_options = {}
    if caption == "":
        caption = f"Created {datetime.now().astimezone().isoformat(sep=" ", timespec="seconds")} on {platform.node()}"
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
