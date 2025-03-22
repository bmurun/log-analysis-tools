import numpy as np
import os
import click
import sys
import multiprocessing
import itertools
from mcap.reader import make_reader
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message
from bokeh.plotting import figure, output_file, save
from bokeh.layouts import gridplot
from bokeh.plotting import figure
from bokeh.models import Label, CrosshairTool, Span, Circle, BoxAnnotation
from bokeh.palettes import Category20, Category10
from bokeh.io import output_file, show, save
from bokeh.models import ColumnDataSource, CustomJS, HoverTool
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn

def create_bokeh_figure(title, x_label, y_label):

    p = figure(max_width=500, height=250,
                    title=title,
                    tools="pan,wheel_zoom,box_zoom,reset",
                    active_scroll="wheel_zoom")

    width = Span(
        dimension="width",
        line_dash="dashed",
        line_color="grey",
        line_width=1,
        line_alpha=0.5,
    )
    height = Span(
        dimension="height",
        line_dash="dashed",
        line_color="grey",
        line_width=1,
        line_alpha=0.5,
    )
    p.add_tools(CrosshairTool(overlay=[width, height]))
    p.xaxis.axis_label = x_label
    p.yaxis.axis_label = y_label
    return p

def get_mcap_data(log_path):
    if not os.path.isfile(log_path):
        click.echo(click.style("Input file doesn't exist!", fg="red"))
        sys.exit()

    topics = [
        "/tracked_pose",
        "/odometry/gps",
        "/odometry/pose_in_map",
        "/debug/localization/kinematic_state"
    ]

    topics_misc = [
        "/debug/localization/is_ego_indoors",
        "/lvx_client/gsof/ins_solution_49"
    ]

    data_dict = {topic: {
        "timestamp_s": [],
        "pos_x": [],
        "pos_y": [],
        "pos_z": [],
        "pos_x_stdev": [],
        "pos_y_stdev": [],
        "pos_z_stdev": []
    } for topic in topics}

    data_dict["/debug/localization/is_ego_indoors"] = {
        "timestamp_s": [],
        "is_ego_indoors": []
    }

    data_dict["/lvx_client/gsof/ins_solution_49"] = {
        "timestamp_s": [],
        "gnss_status": []
    }

    with open(log_path, "rb") as f:
        reader = make_reader(f)

        for schema, channel, message in reader.iter_messages(topics=topics + topics_misc):

            topic_name = channel.topic
            msg_str = schema.name
            msg_type = get_message(msg_str)
            data = deserialize_message(message.data, msg_type)
            data_dict[topic_name]["timestamp_s"].append(data.header.stamp.sec + data.header.stamp.nanosec / 1e9)

            if (topic_name == "/tracked_pose"):

                pos = data.pose.position
                data_dict[topic_name]["pos_x"].append(pos.x)
                data_dict[topic_name]["pos_y"].append(pos.y)
                data_dict[topic_name]["pos_z"].append(pos.z)

            elif (topic_name == "/debug/localization/is_ego_indoors"):
                data_dict["/debug/localization/is_ego_indoors"]["timestamp_s"].append(data.header.stamp.sec + data.header.stamp.nanosec / 1e9)
                data_dict["/debug/localization/is_ego_indoors"]["is_ego_indoors"].append(data.data)

            elif (topic_name == "/lvx_client/gsof/ins_solution_49"):
                data_dict["/lvx_client/gsof/ins_solution_49"]["timestamp_s"].append(data.header.stamp.sec + data.header.stamp.nanosec / 1e9)
                data_dict["/lvx_client/gsof/ins_solution_49"]["gnss_status"].append(data.status.gnss)
            else:
                pos = data.pose.pose.position
                cov = data.pose.covariance
                data_dict[topic_name]["pos_x"].append(pos.x)
                data_dict[topic_name]["pos_y"].append(pos.y)
                data_dict[topic_name]["pos_z"].append(pos.z)

                data_dict[topic_name]["pos_x_stdev"].append(np.sqrt(cov[0]))
                data_dict[topic_name]["pos_y_stdev"].append(np.sqrt(cov[7]))
                data_dict[topic_name]["pos_z_stdev"].append(np.sqrt(cov[14]))

    min_timestamp = np.inf
    # Subtract the first timestamp from all
    for k, v in data_dict.items():
        try:
            min_timestamp_v = np.min(v["timestamp_s"])
            min_timestamp = np.minimum(min_timestamp, min_timestamp_v)
        except:
            continue

    for k, v in data_dict.items():
        try:
            v["timestamp_s"] -= min_timestamp
        except:
            continue

    return data_dict

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

    data_container = {}
    for file_path in file_paths:
        file_name, file_extension = os.path.splitext(os.path.basename(file_path))
        data_container[file_name] = {}

    for file_path in file_paths:
        # pool = multiprocessing.Pool() Do it later
        file_name, file_extension = os.path.splitext(os.path.basename(file_path))
        click.echo(click.style("Input ulog: " + file_name, fg="green", bold=True))
        data_container = get_mcap_data(file_path)

        file_path = os.path.join(output_fqp, file_name + "_position_plots.html")
        line_colors = line_colors = itertools.cycle(Category10[10])

        output_file(filename=file_path, title="Position Plots")

        x_plot = create_bokeh_figure("X Position", "Timestamp", "X [m]")
        y_plot = create_bokeh_figure("X Position", "Timestamp", "Y [m]")
        z_plot = create_bokeh_figure("X Position", "Timestamp", "Z [m]")

        x_stdev_plot = create_bokeh_figure("X Position Stdev", "Timestamp", "Stdev [m]")
        y_stdev_plot = create_bokeh_figure("X Position Stdev", "Timestamp", "Stdev [m]")
        z_stdev_plot = create_bokeh_figure("X Position Stdev", "Timestamp", "Stdev [m]")

        # /odometry/gps
        odometry_gps = data_container["/odometry/gps"]
        odometry_gps_timestamp = odometry_gps["timestamp_s"]
        odometry_gps_x = odometry_gps["pos_x"]
        odometry_gps_y = odometry_gps["pos_y"]
        odometry_gps_z = odometry_gps["pos_z"]
        odometry_gps_x_stdev = odometry_gps["pos_x_stdev"]
        odometry_gps_y_stdev = odometry_gps["pos_y_stdev"]
        odometry_gps_z_stdev = odometry_gps["pos_z_stdev"]

        # /localization/kinematic_state
        ekf = data_container["/debug/localization/kinematic_state"]
        ekf_timestamp = ekf["timestamp_s"][1:]
        ekf_x = ekf["pos_x"][1:]
        ekf_y = ekf["pos_y"][1:]
        ekf_z = ekf["pos_z"][1:]

        ekf_x_stdev = ekf["pos_x_stdev"][1:]
        ekf_y_stdev = ekf["pos_y_stdev"][1:]
        ekf_z_stdev = ekf["pos_z_stdev"][1:]

        # /odometry/pose_in_map
        pose_in_map = data_container["/odometry/pose_in_map"]
        pose_in_map_timestamp = pose_in_map["timestamp_s"]
        pose_in_map_x = pose_in_map["pos_x"]
        pose_in_map_y = pose_in_map["pos_y"]
        pose_in_map_z = pose_in_map["pos_z"]
        pose_in_map_x_stdev = pose_in_map["pos_x_stdev"]
        pose_in_map_y_stdev = pose_in_map["pos_y_stdev"]
        pose_in_map_z_stdev = pose_in_map["pos_z_stdev"]

        # /tracked_pose
        tracked_pose = data_container["/tracked_pose"]
        tracked_pose_timestamp = tracked_pose["timestamp_s"]
        tracked_pose_x = tracked_pose["pos_x"]
        tracked_pose_y = tracked_pose["pos_y"]
        tracked_pose_z = tracked_pose["pos_z"]

        x_plot.line(odometry_gps_timestamp, odometry_gps_x, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
        x_plot.line(tracked_pose_timestamp, tracked_pose_x, alpha=0.8, color="red", line_width=2, legend_label="/tracked_pose")
        x_plot.line(pose_in_map_timestamp, pose_in_map_x, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
        x_plot.line(ekf_timestamp, ekf_x, alpha=0.8, color="cyan", line_width=2, legend_label="/output_kinematic_state")
        y_plot.line(odometry_gps_timestamp, odometry_gps_y, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
        y_plot.line(tracked_pose_timestamp, tracked_pose_y, alpha=0.8, color="red", line_width=2, legend_label="/tracked_pose")
        y_plot.line(pose_in_map_timestamp, pose_in_map_y, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
        y_plot.line(ekf_timestamp, ekf_y, alpha=0.8, color="cyan", line_width=2, legend_label="/output_kinematic_state")
        z_plot.line(odometry_gps_timestamp, odometry_gps_z, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
        z_plot.line(tracked_pose_timestamp, tracked_pose_z, alpha=0.8, color="red", line_width=2, legend_label="/tracked_pose")
        z_plot.line(pose_in_map_timestamp, pose_in_map_z, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
        z_plot.line(ekf_timestamp, ekf_z, alpha=0.8, color="cyan", line_width=2, legend_label="/output_kinematic_state")

        x_stdev_plot.line(odometry_gps_timestamp, odometry_gps_x_stdev, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
        x_stdev_plot.line(pose_in_map_timestamp, pose_in_map_x_stdev, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
        x_stdev_plot.line(ekf_timestamp, ekf_x_stdev, alpha=0.8, color="cyan", line_width=2, legend_label="/output_kinematic_state")
        y_stdev_plot.line(odometry_gps_timestamp, odometry_gps_y_stdev, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
        y_stdev_plot.line(pose_in_map_timestamp, pose_in_map_y_stdev, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
        y_stdev_plot.line(ekf_timestamp, ekf_y_stdev, alpha=0.8, color="cyan", line_width=2, legend_label="/output_kinematic_state")
        z_stdev_plot.line(odometry_gps_timestamp, odometry_gps_z_stdev, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
        z_stdev_plot.line(pose_in_map_timestamp, pose_in_map_z_stdev, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
        z_stdev_plot.line(ekf_timestamp, ekf_z_stdev, alpha=0.8, color="cyan", line_width=2, legend_label="/output_kinematic_state")

        # Post mumbo jumbo
        x_plot.legend.location = "top_right"
        x_plot.legend.click_policy = "hide"
        y_plot.legend.location = "top_right"
        y_plot.legend.click_policy = "hide"
        z_plot.legend.location = "top_right"
        z_plot.legend.click_policy = "hide"
        x_stdev_plot.legend.location = "top_right"
        x_stdev_plot.legend.click_policy = "hide"
        y_stdev_plot.legend.location = "top_right"
        y_stdev_plot.legend.click_policy = "hide"
        z_stdev_plot.legend.location = "top_right"
        z_stdev_plot.legend.click_policy = "hide"

        save(gridplot([
            [x_plot, y_plot, z_plot],
            [x_stdev_plot, y_stdev_plot, z_stdev_plot]
        ]))

if __name__ == "__main__":
    generate_mcap_plots()