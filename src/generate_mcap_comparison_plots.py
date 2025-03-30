import numpy as np
import os
import click
import sys
import multiprocessing
from functools import partial
from common.data_parser.data_parser import DataParser
from common.plot_utils.comparison_plot_generator import ComparisonPlotGenerator

def parse_file(file_path):
    file_name, _ = os.path.splitext(os.path.basename(file_path))
    click.echo(click.style(f"Parsing mcap: {file_name}", fg="green", bold=True))
    data_parser = DataParser(file_path)
    return file_name, data_parser.mcap_data

@click.command()
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
def generate_mcap_comparison_plots(dir_path, multiprocess=False):

    file_paths = []

    for subdir, dirs, files in os.walk(dir_path):
        for file in files:
            fqp = subdir + os.sep + file

            if fqp.endswith(".mcap"):
                file_paths.append(fqp)

    output_folder_name = "mcap_comparison_plots"
    dir_name = dir_path
    output_fqp = os.path.join(dir_name, output_folder_name)

    if not os.path.exists(output_fqp):
        os.makedirs(output_fqp)
    
    output_fqp = os.path.abspath(output_fqp)
    
    data_container = {}

    if multiprocess:
        click.echo(click.style("Multiprocessing enabled!", fg="green", bold=True))
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.map(parse_file, file_paths)
        data_container = {file_name: mcap_data for file_name, mcap_data in results}
    else:
        for file_path in file_paths:
            file_name, file_extension = os.path.splitext(os.path.basename(file_path))
            click.echo(click.style("Input mcap: " + file_name, fg="green", bold=True))
            data_parser = DataParser(file_path)
            data_container[file_name] = data_parser.mcap_data
    
    plot_generator = ComparisonPlotGenerator(data_container)
    plot_generator.generate_plots(output_fqp)

if __name__ == "__main__":
    generate_mcap_comparison_plots()