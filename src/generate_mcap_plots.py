import numpy as np
import os
import click
import sys
from common.data_parser.data_parser import DataParser
from common.plot_utils.plot_generator import PlotGenerator


@click.command()
@click.option("-l", "--log-path", required=False, type=str, help="Path to .mcap")
@click.option(
    "-d",
    "--dir-path",
    required=False,
    type=str,
    help="Path to a root directory containing .mcap files",
)
def generate_mcap_plots(log_path, dir_path):

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
    
    for file_path in file_paths:
        file_name, file_extension = os.path.splitext(os.path.basename(file_path))
        click.echo(click.style("Input mcap: " + file_name, fg="green", bold=True))

        data_parser = DataParser(file_path)
        output_file_path = os.path.join(output_fqp, file_name + "_state_plots.html")
        state_plot_generator = PlotGenerator(data_parser).generate_state_plots(output_file_path)



if __name__ == "__main__":
    generate_mcap_plots()