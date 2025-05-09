import numpy as np
import os
import click
import sys
import multiprocessing
from common.data_parser.data_parser import DataParser
from common.plot_utils.plot_generator import PlotGenerator

def process_mcap_file(file_path, output_fqp):
    file_name, ext = os.path.splitext(os.path.basename(file_path))
    click.echo(click.style(f"Processing mcap: {file_name + ext}", fg="green", bold=True))

    data_parser = DataParser(file_path)
    se_output_file_path = os.path.join(output_fqp, file_name + "_state_plots.html")
    imu_output_file_path = os.path.join(output_fqp, file_name + "_imu_plots.html")
    nis_output_file_path = os.path.join(output_fqp, file_name + "_nis_metrics_plots.html")

    plot_generator = PlotGenerator(data_parser)
    plot_generator.generate_state_plots(se_output_file_path)
    # plot_generator.generate_imu_plots(imu_output_file_path)
    # plot_generator.generate_nis_metrics_plots(nis_output_file_path)

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
def generate_mcap_plots(log_path, dir_path, multiprocess=False):

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
    output_fqp = os.path.join(dir_name, output_folder_name)

    if not os.path.exists(output_fqp):
        os.makedirs(output_fqp)
    
    output_fqp = os.path.abspath(output_fqp)

    if multiprocess:
        click.echo(click.style("Multiprocessing enabled!", fg="green", bold=True))
        with multiprocessing.Pool() as pool:
            pool.starmap(process_mcap_file, [(fp, output_fqp) for fp in file_paths])
    else:
        for file_path in file_paths:
            process_mcap_file(file_path, output_fqp)
    
if __name__ == "__main__":
    generate_mcap_plots()