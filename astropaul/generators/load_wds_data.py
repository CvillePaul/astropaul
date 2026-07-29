from pathlib import Path

from astropy.coordinates import SkyCoord, Angle
from astropy.table import Column, QTable
import astropy.units as u
import pandas as pd

from astropaul.database import base_path, database_connection, database_path


def load_wds_data(out_dir: Path = None, wds_file: Path = None, verbose: bool = True) -> None:

    out_dir = out_dir or database_path().parent
    wds_file = wds_file or base_path() / "Files" / "Washington Double Star Catalog" / "wdsweb_summ2.ecsv"

    target_query = """
        select target_name, ra, dec 
        from targets
        where target_type in ("QuadEB", "SextEB")
        ;"""
    with database_connection() as conn:
        targets = pd.read_sql(target_query, conn)
    target_coords = SkyCoord(ra=targets["ra"], dec=targets["dec"], unit=(u.deg, u.deg))

    wds = QTable.read(wds_file)
    wds_coords = SkyCoord(ra=wds["RA"], dec=wds["Dec"], unit=(u.deg, u.deg))

    # find matches for each target in the WDS
    target_names = [""] * len(wds)
    idx, d2d, _ = target_coords.match_to_catalog_sky(wds_coords)
    for i, (wds_index, distance) in enumerate(zip(idx, d2d)):
        if distance < 1 * u.arcmin:
            target_names[wds_index] = targets.loc[i, "target_name"]
    wds.rename_column("ID", "WDS ID")
    wds.add_column(Column(target_names, name="Target Name"), index=0)
    wds = wds[wds["Target Name"] != ""]

    # write out WDS ID mapping
    out_file = out_dir / "Catalog Members" / "WDS.csv"
    catalog_members = pd.DataFrame()
    catalog_members["Target Name"] = wds["Target Name"]
    catalog_members["Catalog Name"] = "Washington Double Star"
    catalog_members["Catalog ID"] = wds["WDS ID"]
    catalog_members.to_csv(out_file, index=False)

    # write out matched WDS entries
    out_file = out_dir / "WDS" / "WDS.csv"
    wds.write(out_file, overwrite=True)


    if verbose:
        print(f"{len(targets)} targets had {len(wds)} WDS matches, written to {out_file}")


def process_wds_file(filename: Path) -> None:
    """Read in the text format version of the WDS and save as an Astropy QTable ecsv format in the same directory"""
    error_lines = []
    with open(filename, "r", newline="\n") as f:
        for _ in range(5):
            f.readline()
        rows = []
        row = {}
        for line in f.readlines():
            try:
                line = line.strip()  # remove EOL character
                row["ID"] = line[0:10]
                row["Discoverer"] = line[10:17]
                row["Components"] = line[17:22].strip()
                row["Num Obs"] = int(line[33:37])
                row["Separation"] = (float(line[46:51]) + float(line[51:57])) / 2 * u.arcsec
                row["Position Angle"] = int(line[42:45])
                maga_text = line[58:63]
                if not maga_text.strip() == ".":
                    row["Mag A"] = float(maga_text)
                    magb_text = line[63:69]
                    if not magb_text.strip() == ".":
                        row["Mag B"] = float(magb_text)
                        row["Mag Avg"] = (row["Mag A"] + row["Mag B"]) / 2
                    else:
                        row["Mag Avg"] = row["Mag A"]
                else:
                    row["Mag A"] = 0
                    row["Mag B"] = 0
                    row["Mag Avg"] = 0
                spectral_text = line[70:79]
                if not spectral_text.strip() == "":
                    row["Spectral Type"] = spectral_text
                notes_text = line[107:111]
                if not notes_text.strip() == "":
                    row["Note"] = notes_text
                ra_text = line[112:114].strip()
                if ra_text and ra_text != ".":
                    row["RA"] = Angle(f"{ra_text}h{line[114:116]}m{line[116:121]}s").to(u.deg)
                dec_text = line[121:124].strip()
                if dec_text and dec_text != ".":
                    row["Dec"] = Angle(f"{dec_text}d{line[124:126]}m{line[126:130]}s").to(u.deg)
                pmra_text = line[80:84]
                row["PM RA"] = u.Quantity(0 if pmra_text.strip() == "" else int(pmra_text), unit=u.mas / u.yr)
                pmdec_text = line[84:88]
                row["PM Dec"] = u.Quantity(0 if pmdec_text.strip() == "" else int(pmdec_text), unit=u.mas / u.yr)
                rows.append(row)
                row = {}
            except Exception as e:
                error_lines.append(line)
                continue
    print(f"{len(error_lines)} errors and {len(rows)} lines successfully processed")
    print("Error lines")
    for line in error_lines:
        print(line)

    wds = QTable(
        rows,
        names=[
            "ID",
            "Discoverer",
            "Components",
            "Num Obs",
            "Separation",
            "Position Angle",
            "Mag A",
            "Mag B",
            "Mag Avg",
            "Spectral Type",
            "Note",
            "RA",
            "Dec",
            "PM RA",
            "PM Dec",
        ],
    )
    wds["Separation"].info.format = ".2f"
    wds["Mag A"].info.format = ".2f"
    wds["Mag B"].info.format = ".2f"
    wds["Mag Avg"].info.format = ".2f"
    wds["RA"].info.format = "12.8f"
    wds["Dec"].info.format = "13.8f"
    wds.write(filename.parent / "wdsweb_summ2.csv", overwrite=True)
    wds.write(filename.parent / "wdsweb_summ2.ecsv", overwrite=True)
