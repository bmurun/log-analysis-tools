import numpy as np
import os
import click
import sys
import multiprocessing
from common.data_parser.gsof_data_parser import GsofDataParser
from common.plot_utils.gsof_plot_generator import GsofPlotGenerator

def process_mcap_file(file_path, output_root_dir):
    file_name, ext = os.path.splitext(os.path.basename(file_path))
    click.echo(click.style(f"Processing mcap: {file_name + ext}", fg="green", bold=True))

    output_fqp = os.path.join(output_root_dir, file_name)
    output_fqp = os.path.join(output_fqp, "gsof_plots")

    if not os.path.exists(output_fqp):
        os.makedirs(output_fqp)
    
    data_parser = GsofDataParser(file_path)

    pvt_plot_output_file_path = os.path.join(output_fqp, file_name + "_gsof_pvt_plots.html")
    rtk_status_plot_output_file_path = os.path.join(output_fqp, file_name + "_gsof_rtk_status_plots.html")
    sky_plot_output_file_path = os.path.join(output_fqp, file_name + "_gsof_sky_plots.html")
    sat_plot_output_file_path = os.path.join(output_fqp, file_name + "_gsof_sat_plots.html")
    noise_plot_output_file_path = os.path.join(output_fqp, file_name + "_gsof_noise_plots.html")
    dop_plot_output_file_path = os.path.join(output_fqp, file_name + "_gsof_dop_plots.html")

    plot_generator = GsofPlotGenerator(data_parser)
    plot_generator.generate_pvt_plots(pvt_plot_output_file_path)
    plot_generator.generate_rtk_status_plots(rtk_status_plot_output_file_path)
    plot_generator.generate_sky_plots(sky_plot_output_file_path)
    plot_generator.generate_satellite_plots(sat_plot_output_file_path)
    plot_generator.generate_noise_plots(noise_plot_output_file_path)
    plot_generator.generate_dop_plots(dop_plot_output_file_path)

@click.command()
@click.option("-l", "--log-path", required=False, type=str, help="Path to .mcap")
@click.option(
    "-d",
    "--dir-path",
    required=False,
    type=str,
    help="Path to a root directory containing .mcap files",
)
@click.option(
    "-m", "--multiprocess", required=False, is_flag=True, help="Enable multiprocessing"
)
def generate_gsof_plots(log_path, dir_path, multiprocess=False):

    if log_path and dir_path:
        click.echo(click.style("Multiple options are not supported!", fg="red"))
        sys.exit()

    file_paths = []

    if log_path:
        file_paths.append(log_path)
    elif dir_path:

        for subdir, dirs, files in os.walk(dir_path):
            for file in files:
                fqp = subdir + os.sep + file

                if fqp.endswith(".mcap"):
                    file_paths.append(fqp)

    output_folder_name = "mcap_plots"
    dir_name = dir_path if dir_path else os.path.dirname(log_path)
    output_root_dir = os.path.join(dir_name, output_folder_name)

    if not os.path.exists(output_root_dir):
        os.makedirs(output_root_dir)
    
    output_root_dir = os.path.abspath(output_root_dir)

    if multiprocess:
        click.echo(click.style("Multiprocessing enabled!", fg="green", bold=True))
        with multiprocessing.Pool() as pool:
            pool.starmap(process_mcap_file, [(fp, output_root_dir) for fp in file_paths])
    else:
        for file_path in file_paths:
            process_mcap_file(file_path, output_root_dir)
    
if __name__ == "__main__":
    generate_gsof_plots()