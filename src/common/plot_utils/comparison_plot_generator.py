import os
import itertools

from bokeh.plotting import figure, output_file, save
from bokeh.layouts import gridplot
from bokeh.plotting import figure
from bokeh.models import Label, CrosshairTool, Span, Circle, BoxAnnotation
from bokeh.palettes import Category20, Category10
from bokeh.io import output_file, show, save
from bokeh.models import ColumnDataSource, CustomJS, HoverTool
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn
from scipy.spatial.transform import Rotation as R

import numpy as np


class ComparisonPlotGenerator:

    def __init__(self, data_container) -> None:
        self.data_container = data_container

    def create_bokeh_figure(self, title, x_label, y_label):
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
    
    def generate_plots(self, output_file_name):
        
        filename = os.path.join(output_file_name, "comparison_plots.html")
        output_file(filename=filename, title="Comparison Plots")

        line_colors = itertools.cycle(Category10[10])

        odometry_gps_x_plot = self.create_bokeh_figure("/odometry/gps - X Position", "Timestamp", "X [m]")
        odometry_gps_y_plot = self.create_bokeh_figure("/odometry/gps - Y Position", "Timestamp", "Y [m]")
        odometry_gps_z_plot = self.create_bokeh_figure("/odometry/gps - Z Position", "Timestamp", "Z [m]")

        odometry_gps_x_stdev_plot = self.create_bokeh_figure("/odometry/gps - X Position Stdev", "Timestamp", "Stdev [m]")
        odometry_gps_y_stdev_plot = self.create_bokeh_figure("/odometry/gps - Y Position Stdev", "Timestamp", "Stdev [m]")
        odometry_gps_z_stdev_plot = self.create_bokeh_figure("/odometry/gps - Z Position Stdev", "Timestamp", "Stdev [m]")

        ekf_x_plot = self.create_bokeh_figure("/localization/kinematic_state - X Position", "Timestamp", "X [m]")
        ekf_y_plot = self.create_bokeh_figure("/localization/kinematic_state - Y Position", "Timestamp", "Y [m]")
        ekf_z_plot = self.create_bokeh_figure("/localization/kinematic_state - Z Position", "Timestamp", "Z [m]")

        ekf_x_stdev_plot = self.create_bokeh_figure("/localization/kinematic_state - X Position Stdev", "Timestamp", "Stdev [m]")
        ekf_y_stdev_plot = self.create_bokeh_figure("/localization/kinematic_state - Y Position Stdev", "Timestamp", "Stdev [m]")
        ekf_z_stdev_plot = self.create_bokeh_figure("/localization/kinematic_state - Z Position Stdev", "Timestamp", "Stdev [m]")

        ekf_vx_plot = self.create_bokeh_figure("/localization/kinematic_state - X velocity", "Timestamp", "Vx [m/s]")
        ekf_vy_plot = self.create_bokeh_figure("/localization/kinematic_state - Y velocity", "Timestamp", "Vy [m/s]")
        ekf_vz_plot = self.create_bokeh_figure("/localization/kinematic_state - Z velocity", "Timestamp", "Vz [m/s]")

        ekf_vx_stdev_plot = self.create_bokeh_figure("/localization/kinematic_state - X Velocity Stdev", "Timestamp", "Stdev [m/s]")
        ekf_vy_stdev_plot = self.create_bokeh_figure("/localization/kinematic_state - Y Velocity Stdev", "Timestamp", "Stdev [m/s]")
        ekf_vz_stdev_plot = self.create_bokeh_figure("/localization/kinematic_state - Z Velocity Stdev", "Timestamp", "Stdev [m/s]")

        ekf_roll_plot = self.create_bokeh_figure("/localization/kinematic_state - Roll", "Timestamp", "Roll [deg]")
        ekf_pitch_plot = self.create_bokeh_figure("/localization/kinematic_state - Pitch", "Timestamp", "Pitch [deg]")
        ekf_yaw_plot = self.create_bokeh_figure("/localization/kinematic_state - Yaw", "Timestamp", "Yaw [deg]")


        pose_in_map_x_plot = self.create_bokeh_figure("/odometry/pose_in_map - X Position", "Timestamp", "X [m]")
        pose_in_map_y_plot = self.create_bokeh_figure("/odometry/pose_in_map - Y Position", "Timestamp", "Y [m]")
        pose_in_map_z_plot = self.create_bokeh_figure("/odometry/pose_in_map - Z Position", "Timestamp", "Z [m]")

        pose_in_map_x_stdev_plot = self.create_bokeh_figure("/odometry/pose_in_maps - X Position Stdev", "Timestamp", "Stdev [m]")
        pose_in_map_y_stdev_plot = self.create_bokeh_figure("/odometry/pose_in_maps - Y Position Stdev", "Timestamp", "Stdev [m]")
        pose_in_map_z_stdev_plot = self.create_bokeh_figure("/odometry/pose_in_maps - Z Position Stdev", "Timestamp", "Stdev [m]")

        tracked_pose_x_plot = self.create_bokeh_figure("/tracked_pose - X Position", "Timestamp", "X [m]")
        tracked_pose_y_plot = self.create_bokeh_figure("/tracked_pose - Y Position", "Timestamp", "Y [m]")
        tracked_pose_z_plot = self.create_bokeh_figure("/tracked_pose - Z Position", "Timestamp", "Z [m]")

        for (k, mcap_data), color in zip(self.data_container.items(), line_colors):

            is_odometry_gps_available = False
            is_pose_in_map_available = False
            is_tracked_pose_available = False
            is_debug_ekf_available = False
            is_ekf_available = False
            is_debug_ekf_is_ego_indoors_available = False
            is_ekf_is_ego_indoors_available = False
            is_ins_solution_49_available = False

            if "/odometry/gps" in mcap_data.keys():
                odometry_gps = mcap_data["/odometry/gps"]
                odometry_gps_timestamp = odometry_gps["timestamp"]
                odometry_gps_x = odometry_gps["pose"]["pose"]["position"]["x"]
                odometry_gps_y = odometry_gps["pose"]["pose"]["position"]["y"]
                odometry_gps_z = odometry_gps["pose"]["pose"]["position"]["z"]
                odometry_gps_x_stdev = np.sqrt(odometry_gps["pose"]["covariance"][:, 0])
                odometry_gps_y_stdev = np.sqrt(odometry_gps["pose"]["covariance"][:, 7])
                odometry_gps_z_stdev = np.sqrt(odometry_gps["pose"]["covariance"][:, 14])
                is_odometry_gps_available = True

                odometry_gps_x_plot.line(odometry_gps_timestamp, odometry_gps_x, alpha=0.8, color=color, line_width=2, legend_label=k)
                odometry_gps_y_plot.line(odometry_gps_timestamp, odometry_gps_y, alpha=0.8, color=color, line_width=2, legend_label=k)
                odometry_gps_z_plot.line(odometry_gps_timestamp, odometry_gps_z, alpha=0.8, color=color, line_width=2, legend_label=k)

                odometry_gps_x_stdev_plot.line(odometry_gps_timestamp, odometry_gps_x_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                odometry_gps_y_stdev_plot.line(odometry_gps_timestamp, odometry_gps_y_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                odometry_gps_z_stdev_plot.line(odometry_gps_timestamp, odometry_gps_z_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)

            if "/odometry/pose_in_map" in mcap_data.keys():
                pose_in_map = mcap_data["/odometry/pose_in_map"]
                pose_in_map_timestamp = pose_in_map["timestamp"]
                pose_in_map_x = pose_in_map["pose"]["pose"]["position"]["x"]
                pose_in_map_y = pose_in_map["pose"]["pose"]["position"]["y"]
                pose_in_map_z = pose_in_map["pose"]["pose"]["position"]["z"]
                pose_in_map_x_stdev = np.sqrt(pose_in_map["pose"]["covariance"][:, 0])
                pose_in_map_y_stdev = np.sqrt(pose_in_map["pose"]["covariance"][:, 7])
                pose_in_map_z_stdev = np.sqrt(pose_in_map["pose"]["covariance"][:, 14])
                is_pose_in_map_available = True

                pose_in_map_x_plot.line(pose_in_map_timestamp, pose_in_map_x, alpha=0.8, color=color, line_width=2, legend_label=k)
                pose_in_map_y_plot.line(pose_in_map_timestamp, pose_in_map_y, alpha=0.8, color=color, line_width=2, legend_label=k)
                pose_in_map_z_plot.line(pose_in_map_timestamp, pose_in_map_z, alpha=0.8, color=color, line_width=2, legend_label=k)

                pose_in_map_x_stdev_plot.line(pose_in_map_timestamp, pose_in_map_x_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                pose_in_map_y_stdev_plot.line(pose_in_map_timestamp, pose_in_map_y_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                pose_in_map_z_stdev_plot.line(pose_in_map_timestamp, pose_in_map_z_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)

            if "/tracked_pose" in mcap_data.keys():
                tracked_pose = mcap_data["/tracked_pose"]
                tracked_pose_timestamp = tracked_pose["timestamp"]
                tracked_pose_x = tracked_pose["pose"]["position"]["x"]
                tracked_pose_y = tracked_pose["pose"]["position"]["y"]
                tracked_pose_z = tracked_pose["pose"]["position"]["z"]
                is_tracked_pose_available = True

                tracked_pose_x_plot.line(tracked_pose_timestamp, tracked_pose_x, alpha=0.8, color=color, line_width=2, legend_label=k)
                tracked_pose_y_plot.line(tracked_pose_timestamp, tracked_pose_y, alpha=0.8, color=color, line_width=2, legend_label=k)
                tracked_pose_z_plot.line(tracked_pose_timestamp, tracked_pose_z, alpha=0.8, color=color, line_width=2, legend_label=k)

            if "/debug/localization/kinematic_state" in mcap_data.keys():
                debug_ekf = mcap_data["/debug/localization/kinematic_state"]
                debug_ekf_timestamp = debug_ekf["timestamp"][1:]
                debug_ekf_x = debug_ekf["pose"]["pose"]["position"]["x"][1:]
                debug_ekf_y = debug_ekf["pose"]["pose"]["position"]["y"][1:]
                debug_ekf_z = debug_ekf["pose"]["pose"]["position"]["z"][1:]
                debug_ekf_x_stdev = np.sqrt(debug_ekf["pose"]["covariance"][1:, 0])
                debug_ekf_y_stdev = np.sqrt(debug_ekf["pose"]["covariance"][1:, 7])
                debug_ekf_z_stdev = np.sqrt(debug_ekf["pose"]["covariance"][1:, 14])

                debug_ekf_vx = debug_ekf["twist"]["twist"]["linear"]["x"][1:]
                debug_ekf_vy = debug_ekf["twist"]["twist"]["linear"]["y"][1:]
                debug_ekf_vz = debug_ekf["twist"]["twist"]["linear"]["z"][1:]
                debug_ekf_vx_stdev = np.sqrt(debug_ekf["twist"]["covariance"][1:, 0])
                debug_ekf_vy_stdev = np.sqrt(debug_ekf["twist"]["covariance"][1:, 7])
                debug_ekf_vz_stdev = np.sqrt(debug_ekf["twist"]["covariance"][1:, 14])

                debug_ekf_quat = debug_ekf["pose"]["pose"]["orientation"]
                debug_ekf_quat["w"] = debug_ekf_quat["w"][1:]
                debug_ekf_quat["x"] = debug_ekf_quat["x"][1:]
                debug_ekf_quat["y"] = debug_ekf_quat["y"][1:]
                debug_ekf_quat["z"] = debug_ekf_quat["z"][1:]

                debug_ekf_roll_deg = np.zeros(debug_ekf_timestamp.shape)
                debug_ekf_pitch_deg = np.zeros(debug_ekf_timestamp.shape)
                debug_ekf_yaw_deg = np.zeros(debug_ekf_timestamp.shape)

                for i in range(debug_ekf_timestamp.shape[0]):
                    w = debug_ekf_quat["w"][i]
                    x = debug_ekf_quat["x"][i]
                    y = debug_ekf_quat["y"][i]
                    z = debug_ekf_quat["z"][i]

                    Rot = R.from_quat(
                        [w, x, y, z], scalar_first=True
                    )

                    rpy = Rot.as_euler("xyz", degrees = True)
                    debug_ekf_roll_deg[i] = rpy[0]
                    debug_ekf_pitch_deg[i] = rpy[1]
                    debug_ekf_yaw_deg[i] = rpy[2]
                is_debug_ekf_available = True

                ekf_x_plot.line(debug_ekf_timestamp, debug_ekf_x, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_y_plot.line(debug_ekf_timestamp, debug_ekf_y, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_z_plot.line(debug_ekf_timestamp, debug_ekf_z, alpha=0.8, color=color, line_width=2, legend_label=k)

                ekf_x_stdev_plot.line(debug_ekf_timestamp, debug_ekf_x_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_y_stdev_plot.line(debug_ekf_timestamp, debug_ekf_y_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_z_stdev_plot.line(debug_ekf_timestamp, debug_ekf_z_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)

                ekf_vx_plot.line(debug_ekf_timestamp, debug_ekf_vx, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_vy_plot.line(debug_ekf_timestamp, debug_ekf_vy, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_vz_plot.line(debug_ekf_timestamp, debug_ekf_vz, alpha=0.8, color=color, line_width=2, legend_label=k)

                ekf_vx_stdev_plot.line(debug_ekf_timestamp, debug_ekf_vx_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_vy_stdev_plot.line(debug_ekf_timestamp, debug_ekf_vy_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_vz_stdev_plot.line(debug_ekf_timestamp, debug_ekf_vz_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)

                ekf_roll_plot.line(debug_ekf_timestamp, debug_ekf_roll_deg, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_pitch_plot.line(debug_ekf_timestamp, debug_ekf_pitch_deg, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_yaw_plot.line(debug_ekf_timestamp, debug_ekf_yaw_deg, alpha=0.8, color=color, line_width=2, legend_label=k)

            if "/localization/kinematic_state" in mcap_data.keys():
                ekf = mcap_data["/localization/kinematic_state"]
                ekf_timestamp = ekf["timestamp"][1:]
                ekf_x = ekf["pose"]["pose"]["position"]["x"][1:]
                ekf_y = ekf["pose"]["pose"]["position"]["y"][1:]
                ekf_z = ekf["pose"]["pose"]["position"]["z"][1:]
                ekf_x_stdev = np.sqrt(ekf["pose"]["covariance"][1:, 0])
                ekf_y_stdev = np.sqrt(ekf["pose"]["covariance"][1:, 7])
                ekf_z_stdev = np.sqrt(ekf["pose"]["covariance"][1:, 14])

                ekf_vx = ekf["twist"]["twist"]["linear"]["x"][1:]
                ekf_vy = ekf["twist"]["twist"]["linear"]["y"][1:]
                ekf_vz = ekf["twist"]["twist"]["linear"]["z"][1:]
                ekf_vx_stdev = np.sqrt(ekf["twist"]["covariance"][1:, 0])
                ekf_vy_stdev = np.sqrt(ekf["twist"]["covariance"][1:, 7])
                ekf_vz_stdev = np.sqrt(ekf["twist"]["covariance"][1:, 14])

                ekf_quat = ekf["pose"]["pose"]["orientation"]
                ekf_quat["w"] = ekf_quat["w"][1:]
                ekf_quat["x"] = ekf_quat["x"][1:]
                ekf_quat["y"] = ekf_quat["y"][1:]
                ekf_quat["z"] = ekf_quat["z"][1:]
                ekf_roll_deg = np.zeros(ekf_timestamp.shape)
                ekf_pitch_deg = np.zeros(ekf_timestamp.shape)
                ekf_yaw_deg = np.zeros(ekf_timestamp.shape)

                for i in range(ekf_timestamp.shape[0]):
                    w = ekf_quat["w"][i]
                    x = ekf_quat["x"][i]
                    y = ekf_quat["y"][i]
                    z = ekf_quat["z"][i]

                    Rot = R.from_quat(
                        [w, x, y, z], scalar_first=True
                    )

                    rpy = Rot.as_euler("xyz", degrees = True)
                    ekf_roll_deg[i] = rpy[0]
                    ekf_pitch_deg[i] = rpy[1]
                    ekf_yaw_deg[i] = rpy[2]
                
                is_ekf_available = True

                ekf_x_plot.line(ekf_timestamp, ekf_x, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_y_plot.line(ekf_timestamp, ekf_y, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_z_plot.line(ekf_timestamp, ekf_z, alpha=0.8, color=color, line_width=2, legend_label=k)

                ekf_x_stdev_plot.line(ekf_timestamp, ekf_x_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_y_stdev_plot.line(ekf_timestamp, ekf_y_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_z_stdev_plot.line(ekf_timestamp, ekf_z_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)

                ekf_vx_plot.line(ekf_timestamp, ekf_vx, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_vy_plot.line(ekf_timestamp, ekf_vy, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_vz_plot.line(ekf_timestamp, ekf_vz, alpha=0.8, color=color, line_width=2, legend_label=k)

                ekf_vx_stdev_plot.line(ekf_timestamp, ekf_vx_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_vy_stdev_plot.line(ekf_timestamp, ekf_vy_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_vz_stdev_plot.line(ekf_timestamp, ekf_vz_stdev, alpha=0.8, color=color, line_width=2, legend_label=k)

                ekf_roll_plot.line(ekf_timestamp, ekf_roll_deg, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_pitch_plot.line(ekf_timestamp, ekf_pitch_deg, alpha=0.8, color=color, line_width=2, legend_label=k)
                ekf_yaw_plot.line(ekf_timestamp, ekf_yaw_deg, alpha=0.8, color=color, line_width=2, legend_label=k)
            
            if "/debug/localization/is_ego_indoors" in mcap_data.keys():
                debug_ekf_is_ego_indoors = mcap_data["/debug/localization/is_ego_indoors"]
                debug_ekf_is_ego_indoors_flag = debug_ekf_is_ego_indoors["data"]

                if (is_debug_ekf_available):
                    debug_ekf_is_ego_indoors_timestamp = debug_ekf_timestamp
                    is_debug_ekf_is_ego_indoors_available = True
                elif (is_ekf_available):
                    debug_ekf_is_ego_indoors_timestamp= ekf_timestamp
                    is_debug_ekf_is_ego_indoors_available = True
                else:
                    is_debug_ekf_is_ego_indoors_available = False


            if "/localization/is_ego_indoors" in mcap_data.keys():
                ekf_is_ego_indoors = mcap_data["/localization/is_ego_indoors"]
                ekf_is_ego_indoors_flag = ekf_is_ego_indoors["data"]

                if (is_debug_ekf_available):
                    ekf_is_ego_indoors_timestamp = debug_ekf_timestamp
                    is_ekf_is_ego_indoors_available = True
                elif (is_ekf_available):
                    ekf_is_ego_indoors_timestamp= ekf_timestamp
                    is_ekf_is_ego_indoors_available = True
                else:
                    is_ekf_is_ego_indoors_available = False
            
            if "/lvx_client/gsof/ins_solution_49" in mcap_data.keys():
                ins_solution_49 = mcap_data["/lvx_client/gsof/ins_solution_49"]
                ins_solution_49_timestamp = ins_solution_49["timestamp"]
                ins_solution_49_gnss_status = ins_solution_49["status"]["gnss"]
                is_ins_solution_49_available = True

        # Post mumbo jumbo
        odometry_gps_x_plot.legend.location = "top_right"
        odometry_gps_x_plot.legend.click_policy = "hide"
        odometry_gps_y_plot.legend.location = "top_right"
        odometry_gps_y_plot.legend.click_policy = "hide"
        odometry_gps_z_plot.legend.location = "top_right"
        odometry_gps_z_plot.legend.click_policy = "hide"

        odometry_gps_x_stdev_plot.legend.location = "top_right"
        odometry_gps_x_stdev_plot.legend.click_policy = "hide"
        odometry_gps_y_stdev_plot.legend.location = "top_right"
        odometry_gps_y_stdev_plot.legend.click_policy = "hide"
        odometry_gps_z_stdev_plot.legend.location = "top_right"
        odometry_gps_z_stdev_plot.legend.click_policy = "hide"

        ekf_x_plot.legend.location = "top_right"
        ekf_x_plot.legend.click_policy = "hide"
        ekf_y_plot.legend.location = "top_right"
        ekf_y_plot.legend.click_policy = "hide"
        ekf_z_plot.legend.location = "top_right"
        ekf_z_plot.legend.click_policy = "hide"
        ekf_x_stdev_plot.legend.location = "top_right"
        ekf_x_stdev_plot.legend.click_policy = "hide"
        ekf_y_stdev_plot.legend.location = "top_right"
        ekf_y_stdev_plot.legend.click_policy = "hide"
        ekf_z_stdev_plot.legend.location = "top_right"
        ekf_z_stdev_plot.legend.click_policy = "hide"

        ekf_vx_plot.legend.location = "top_right"
        ekf_vx_plot.legend.click_policy = "hide"
        ekf_vy_plot.legend.location = "top_right"
        ekf_vy_plot.legend.click_policy = "hide"
        ekf_vz_plot.legend.location = "top_right"
        ekf_vz_plot.legend.click_policy = "hide"

        ekf_vx_stdev_plot.legend.location = "top_right"
        ekf_vy_stdev_plot.legend.click_policy = "hide"
        ekf_vz_stdev_plot.legend.location = "top_right"
        ekf_vx_stdev_plot.legend.click_policy = "hide"
        ekf_vy_stdev_plot.legend.location = "top_right"
        ekf_vz_stdev_plot.legend.click_policy = "hide"

        ekf_roll_plot.legend.location = "top_right"
        ekf_roll_plot.legend.click_policy = "hide"
        ekf_pitch_plot.legend.location = "top_right"
        ekf_pitch_plot.legend.click_policy = "hide"
        ekf_yaw_plot.legend.location = "top_right"
        ekf_yaw_plot.legend.click_policy = "hide"

        pose_in_map_x_plot.legend.location = "top_right"
        pose_in_map_x_plot.legend.click_policy = "hide"
        pose_in_map_y_plot.legend.location = "top_right"
        pose_in_map_y_plot.legend.click_policy = "hide"
        pose_in_map_z_plot.legend.location = "top_right"
        pose_in_map_z_plot.legend.click_policy = "hide"

        pose_in_map_x_stdev_plot.legend.location = "top_right"
        pose_in_map_x_stdev_plot.legend.click_policy = "hide"
        pose_in_map_y_stdev_plot.legend.location = "top_right"
        pose_in_map_y_stdev_plot.legend.click_policy = "hide"
        pose_in_map_z_stdev_plot.legend.location = "top_right"
        pose_in_map_z_stdev_plot.legend.click_policy = "hide"

        tracked_pose_x_plot.legend.location = "top_right"
        tracked_pose_x_plot.legend.click_policy = "hide"
        tracked_pose_y_plot.legend.location = "top_right"
        tracked_pose_y_plot.legend.click_policy = "hide"
        tracked_pose_z_plot.legend.location = "top_right"
        tracked_pose_z_plot.legend.click_policy = "hide"

        save(gridplot([
            [odometry_gps_x_plot, odometry_gps_y_plot, odometry_gps_z_plot],
            [tracked_pose_x_plot, tracked_pose_y_plot, tracked_pose_z_plot],
            [pose_in_map_x_plot, pose_in_map_y_plot, pose_in_map_z_plot],
            [ekf_x_plot, ekf_y_plot, ekf_z_plot],
            [odometry_gps_x_stdev_plot, odometry_gps_y_stdev_plot, odometry_gps_z_stdev_plot],
            [pose_in_map_x_stdev_plot, pose_in_map_y_stdev_plot, pose_in_map_z_stdev_plot],
            [ekf_x_stdev_plot, ekf_y_stdev_plot, ekf_z_stdev_plot],
            [ekf_vx_plot, ekf_vy_plot, ekf_vz_plot],
            [ekf_vx_stdev_plot, ekf_vy_stdev_plot, ekf_vz_stdev_plot],
            [ekf_roll_plot, ekf_pitch_plot, ekf_yaw_plot]
        ]))
