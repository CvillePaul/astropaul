from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from astropaul.database import database_connection
import astropaul.lbt as lbt


def make_pepsi_resources(fits_dir: Path, resources_dir: Path):
    # make dictionary of all available fits files by name
    fits_files = {
        **{file.name: file for file in fits_dir.glob("20*/*.bwl")},
        **{file.name: file for file in fits_dir.glob("20*/*.nor")},
    }

    spectrum_resource_name = "PEPSI Spectrum"
    spectrum_sql_name = "pepsi_spectra"
    spectrum_dir = Path(spectrum_resource_name)
    if not (resources_dir / spectrum_dir).exists():
        (resources_dir / spectrum_dir).mkdir()

    spectrumplot_resource_name = "PEPSI SpectrumPlot"
    spectrumplot_sql_name = "pepsi_spectrumplots"
    spectrumplot_dir = Path(spectrumplot_resource_name)
    if not (resources_dir / spectrumplot_dir).exists():
        (resources_dir / spectrumplot_dir).mkdir()

    with database_connection() as conn:
        observations = pd.read_sql("select * from pepsi_observations;", conn)
        cursor = conn.cursor()

        spectrum_resources = pd.DataFrame()
        spectrum_resources["id"] = observations["id"]
        spectrum_resources[spectrum_sql_name] = [
            str(spectrum_dir / f"{spectrum_resource_name}.{id.replace(":", "_")}.spec") for id in observations["id"]
        ]
        spectra = observations["spectrum_file"].apply(lambda file: lbt.read_pepsi_file(fits_files[file]))
        [
            spectrum.write(resources_dir / file, format="tabular-fits", overwrite=True)
            for spectrum, file in zip(spectra, spectrum_resources[spectrum_sql_name])
        ]
        spectrum_table_name = f"resource_{spectrum_sql_name}"
        spectrum_resources.to_sql(spectrum_table_name, conn, if_exists="replace", index=False)
        cursor.execute(f"delete from table_metadata where table_name = '{spectrum_table_name}';")
        cursor.execute(
            f"""
            insert into table_metadata 
            (table_name, column_name, value_type, value) 
            values ('{spectrum_table_name}', NULL, 'Associated Table', 'PEPSI Observations')
            ;"""
        )
        print(f"Created {len(spectrum_resources)} {spectrum_resource_name} objects")

        spectrumplot_resources = pd.DataFrame()
        spectrumplot_resources["id"] = observations["id"]
        spectrumplot_resources[spectrumplot_sql_name] = [
            str(spectrumplot_dir / f"{spectrumplot_resource_name}.{id.replace(":", "_")}.png") for id in observations["id"]
        ]
        for spectrum, file in zip(spectra, spectrumplot_resources[spectrumplot_sql_name]):
            fig, ax = lbt.plot_pepsi_spectrum(spectrum)
            ax.set_ylim(-1, 0.5)
            fig.savefig(resources_dir / file)
            plt.close()
        spectrumplot_table_name = f"resource_{spectrumplot_sql_name}"
        spectrumplot_resources.to_sql(spectrumplot_table_name, conn, if_exists="replace", index=False)
        cursor.execute(f"delete from table_metadata where table_name = '{spectrumplot_table_name}';")
        cursor.execute(
            f"""
            insert into table_metadata 
            (table_name, column_name, value_type, value) 
            values ('{spectrumplot_table_name}', NULL, 'Associated Table', 'PEPSI Observations')
            ;"""
        )
        print(f"Created {len(spectrumplot_resources)} {spectrumplot_resource_name} objects")

        conn.commit()
