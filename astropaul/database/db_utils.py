from contextlib import contextmanager
from pathlib import Path

from sqlite3 import Connection


def database_path() -> Path:
    return Path("/Users/User/Dropbox/Astro/Data/astropaul.db")


@contextmanager
def database_connection(path: str = None):
    if not path:
        path = database_path()
    conn = Connection(path)
    try:
        yield conn
    finally:
        conn.close()
