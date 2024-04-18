"""script will:
- generate plots with various zoom-levels for h5-files
- NOTE: let the generator- and converter-example run before

CLI-Version of this is:
shepherd-data plot file_or_dir [-s start_time] [-e end_time]
                               [-w plot_width] [-h plot_height]
                               [--multiplot]
"""

from pathlib import Path

import shepherd_data as shp

if __name__ == "__main__":
    with shp.Reader(Path("./hrv_sawtooth_10min.h5")) as db:
        db.plot_to_file()
        db.plot_to_file(0, 500)
        db.plot_to_file(0, 80)

    with shp.Reader(Path("./jogging_10m_ivcurve.h5")) as db:
        db.plot_to_file()
        db.plot_to_file(0, 100)
        db.plot_to_file(0, 10)
        db.plot_to_file(0, 1)
        db.plot_to_file(0, 0.2)
        db.plot_to_file(0.199, 0.201)

    # multiplot
    files = [
        Path("./jogging_10m_ivcurve.h5"),
        Path("./jogging_10m_isc_voc.h5"),
        Path("./jogging_10m_ivsample_voc.h5"),
        Path("./jogging_10m_ivsample_opt.h5"),
    ]
    data = []
    for file in files:
        with shp.Reader(file, verbose=False) as db:
            date = db.generate_plot_data(0, 1)
            date["name"] = file.stem[12:]  # cut away "jogging_10m_"
            data.append(date)
            print(db.energy())
    shp.Reader.multiplot_to_file(data, Path(files[0].stem[:11]))
