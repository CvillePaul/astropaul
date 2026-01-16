from contextlib import contextmanager
from pathlib import Path

from sqlite3 import Connection


def base_path() -> Path:
    # TODO: change to use a relative path for portability
    return Path("/Users/User/Dropbox/Astro/")


def database_path() -> Path:
    return base_path() / "Data/astropaul.db"


def resources_path() -> Path:
    return base_path() / "Resources/"


def html_path() -> Path:
    return base_path() / "Observing Files/"


@contextmanager
def database_connection(path: str = None):
    if not path:
        path = database_path()
    conn = Connection(path)
    try:
        yield conn
    finally:
        conn.close()


def string_to_db_style(name: str) -> str:
    """Converts a directory or CSV column name to its database equivalent"""
    return name.lower().replace(" ", "_")


def db_style_to_string(name: str) -> str:
    """Converts a name in database style back to a more pretty, human form"""
    answer = name.replace("_", " ")
    answer = answer.title()
    for string in ["Jd", "Utc", "Id", "Dssi", "Ra", "Hms", "Dms", "Pepsi", "Rv", "Pm"]:
        answer = answer.replace(string, string.upper())
    special_cases = [("Zorroalopeke", "ZorroAlopeke"), ("Spectrumplot", "SpectrumPlot"), ("Hrcam", "HRCam")]
    for original, new in special_cases:
        answer = answer.replace(original, new)
    return answer
