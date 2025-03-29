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


class PlotGenerator:

    def __init__(self, data_parser) -> None:
        self.data_parser = data_parser

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

    def generate_state_plots(self, output_file_name):

        output_file(filename=output_file_name, title="Plots")

        x_plot = self.create_bokeh_figure("X Position", "Timestamp", "X [m]")
        y_plot = self.create_bokeh_figure("Y Position", "Timestamp", "Y [m]")
        z_plot = self.create_bokeh_figure("Z Position", "Timestamp", "Z [m]")

        x_stdev_plot = self.create_bokeh_figure("X Position Stdev", "Timestamp", "Stdev [m]")
        y_stdev_plot = self.create_bokeh_figure("Y Position Stdev", "Timestamp", "Stdev [m]")
        z_stdev_plot = self.create_bokeh_figure("Z Position Stdev", "Timestamp", "Stdev [m]")

        vx_plot = self.create_bokeh_figure("X Velocity", "Timestamp", "Vx [m]")
        vy_plot = self.create_bokeh_figure("Y Velocity", "Timestamp", "Vy [m]")
        vz_plot = self.create_bokeh_figure("Z Velocity", "Timestamp", "Vz [m]")

        vx_stdev_plot = self.create_bokeh_figure("X Velocity Stdev", "Timestamp", "Stdev [m/s]")
        vy_stdev_plot = self.create_bokeh_figure("Y Velocity Stdev", "Timestamp", "Stdev [m/s]")
        vz_stdev_plot = self.create_bokeh_figure("Z Velocity Stdev", "Timestamp", "Stdev [m/s]")

        roll_plot = self.create_bokeh_figure("Roll", "Timestamp", "Roll [deg]")
        pitch_plot = self.create_bokeh_figure("Pitch", "Timestamp", "Pitch [deg]")
        yaw_plot = self.create_bokeh_figure("Yaw", "Timestamp", "Yaw [deg]")

        status_plot = self.create_bokeh_figure("Ego Indoors and GNSS Status", "Timestamp", "")

        is_odometry_gps_available = False
        is_pose_in_map_available = False
        is_tracked_pose_available = False
        is_debug_ekf_available = False
        is_ekf_available = False
        is_debug_ekf_is_ego_indoors_available = False
        is_ekf_is_ego_indoors_available = False
        is_ins_solution_49_available = False

        if "/odometry/gps" in self.data_parser.mcap_data.keys():
            odometry_gps = self.data_parser.mcap_data["/odometry/gps"]
            odometry_gps_timestamp = odometry_gps["timestamp"]
            odometry_gps_x = odometry_gps["pose"]["pose"]["position"]["x"]
            odometry_gps_y = odometry_gps["pose"]["pose"]["position"]["y"]
            odometry_gps_z = odometry_gps["pose"]["pose"]["position"]["z"]
            odometry_gps_x_stdev = np.sqrt(odometry_gps["pose"]["covariance"][:, 0])
            odometry_gps_y_stdev = np.sqrt(odometry_gps["pose"]["covariance"][:, 7])
            odometry_gps_z_stdev = np.sqrt(odometry_gps["pose"]["covariance"][:, 14])
            is_odometry_gps_available = True

        if "/odometry/pose_in_map" in self.data_parser.mcap_data.keys():
            pose_in_map = self.data_parser.mcap_data["/odometry/pose_in_map"]
            pose_in_map_timestamp = pose_in_map["timestamp"]
            pose_in_map_x = pose_in_map["pose"]["pose"]["position"]["x"]
            pose_in_map_y = pose_in_map["pose"]["pose"]["position"]["y"]
            pose_in_map_z = pose_in_map["pose"]["pose"]["position"]["z"]
            pose_in_map_x_stdev = np.sqrt(pose_in_map["pose"]["covariance"][:, 0])
            pose_in_map_y_stdev = np.sqrt(pose_in_map["pose"]["covariance"][:, 7])
            pose_in_map_z_stdev = np.sqrt(pose_in_map["pose"]["covariance"][:, 14])
            is_pose_in_map_available = True

        if "/tracked_pose" in self.data_parser.mcap_data.keys():
            tracked_pose = self.data_parser.mcap_data["/tracked_pose"]
            tracked_pose_timestamp = tracked_pose["timestamp"]
            tracked_pose_x = tracked_pose["pose"]["position"]["x"]
            tracked_pose_y = tracked_pose["pose"]["position"]["y"]
            tracked_pose_z = tracked_pose["pose"]["position"]["z"]
            is_tracked_pose_available = True

        if "/debug/localization/kinematic_state" in self.data_parser.mcap_data.keys():
            debug_ekf = self.data_parser.mcap_data["/debug/localization/kinematic_state"]
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

        if "/localization/kinematic_state" in self.data_parser.mcap_data.keys():
            ekf = self.data_parser.mcap_data["/localization/kinematic_state"]
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
        
        if "/debug/localization/is_ego_indoors" in self.data_parser.mcap_data.keys():
            debug_ekf_is_ego_indoors = self.data_parser.mcap_data["/debug/localization/is_ego_indoors"]
            debug_ekf_is_ego_indoors_flag = debug_ekf_is_ego_indoors["data"]

            if (is_debug_ekf_available):
                debug_ekf_is_ego_indoors_timestamp = debug_ekf_timestamp
                is_debug_ekf_is_ego_indoors_available = True
            elif (is_ekf_available):
                debug_ekf_is_ego_indoors_timestamp= ekf_timestamp
                is_debug_ekf_is_ego_indoors_available = True
            else:
                is_debug_ekf_is_ego_indoors_available = False


        if "/localization/is_ego_indoors" in self.data_parser.mcap_data.keys():
            ekf_is_ego_indoors = self.data_parser.mcap_data["/localization/is_ego_indoors"]
            ekf_is_ego_indoors_flag = ekf_is_ego_indoors["data"]

            if (is_debug_ekf_available):
                ekf_is_ego_indoors_timestamp = debug_ekf_timestamp
                is_ekf_is_ego_indoors_available = True
            elif (is_ekf_available):
                ekf_is_ego_indoors_timestamp= ekf_timestamp
                is_ekf_is_ego_indoors_available = True
            else:
                is_ekf_is_ego_indoors_available = False
        
        if "/lvx_client/gsof/ins_solution_49" in self.data_parser.mcap_data.keys():
            ins_solution_49 = self.data_parser.mcap_data["/lvx_client/gsof/ins_solution_49"]
            ins_solution_49_timestamp = ins_solution_49["timestamp"]
            ins_solution_49_gnss_status = ins_solution_49["status"]["gnss"]
            is_ins_solution_49_available = True


        # ============== Odometry GPS ==============
        if is_odometry_gps_available:
            x_plot.line(odometry_gps_timestamp, odometry_gps_x, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
            y_plot.line(odometry_gps_timestamp, odometry_gps_y, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
            z_plot.line(odometry_gps_timestamp, odometry_gps_z, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
            x_stdev_plot.line(odometry_gps_timestamp, odometry_gps_x_stdev, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
            y_stdev_plot.line(odometry_gps_timestamp, odometry_gps_y_stdev, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")
            z_stdev_plot.line(odometry_gps_timestamp, odometry_gps_z_stdev, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps")

        # ============== Pose in Map ==============
        if is_pose_in_map_available:
            x_plot.line(pose_in_map_timestamp, pose_in_map_x, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            y_plot.line(pose_in_map_timestamp, pose_in_map_y, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            z_plot.line(pose_in_map_timestamp, pose_in_map_z, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            x_stdev_plot.line(pose_in_map_timestamp, pose_in_map_x_stdev, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            y_stdev_plot.line(pose_in_map_timestamp, pose_in_map_y_stdev, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            z_stdev_plot.line(pose_in_map_timestamp, pose_in_map_z_stdev, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")

        # ============== Tracked Pose ==============
        if is_tracked_pose_available:
            x_plot.line(tracked_pose_timestamp, tracked_pose_x, alpha=0.8, color="red", line_width=2, legend_label="/tracked_pose")
            y_plot.line(tracked_pose_timestamp, tracked_pose_y, alpha=0.8, color="red", line_width=2, legend_label="/tracked_pose")
            z_plot.line(tracked_pose_timestamp, tracked_pose_z, alpha=0.8, color="red", line_width=2, legend_label="/tracked_pose")

        # ============== Replay EKF ==============
        if is_debug_ekf_available:
            x_plot.line(debug_ekf_timestamp, debug_ekf_x, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            y_plot.line(debug_ekf_timestamp, debug_ekf_y, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            z_plot.line(debug_ekf_timestamp, debug_ekf_z, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            x_stdev_plot.line(debug_ekf_timestamp, debug_ekf_x_stdev, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            y_stdev_plot.line(debug_ekf_timestamp, debug_ekf_y_stdev, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            z_stdev_plot.line(debug_ekf_timestamp, debug_ekf_z_stdev, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")

            vx_plot.line(debug_ekf_timestamp, debug_ekf_vx, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            vy_plot.line(debug_ekf_timestamp, debug_ekf_vy, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            vz_plot.line(debug_ekf_timestamp, debug_ekf_vz, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            vx_stdev_plot.line(debug_ekf_timestamp, debug_ekf_vx_stdev, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            vy_stdev_plot.line(debug_ekf_timestamp, debug_ekf_vy_stdev, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            vz_stdev_plot.line(debug_ekf_timestamp, debug_ekf_vz_stdev, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")

            roll_plot.line(debug_ekf_timestamp, debug_ekf_roll_deg, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            pitch_plot.line(debug_ekf_timestamp, debug_ekf_pitch_deg, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            yaw_plot.line(debug_ekf_timestamp, debug_ekf_yaw_deg, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")

        # ============== Online EKF ==============
        if is_ekf_available:
            x_plot.line(ekf_timestamp, ekf_x, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            y_plot.line(ekf_timestamp, ekf_y, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            z_plot.line(ekf_timestamp, ekf_z, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            x_stdev_plot.line(ekf_timestamp, ekf_x_stdev, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            y_stdev_plot.line(ekf_timestamp, ekf_y_stdev, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            z_stdev_plot.line(ekf_timestamp, ekf_z_stdev, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")

            vx_plot.line(ekf_timestamp, ekf_vx, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            vy_plot.line(ekf_timestamp, ekf_vy, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            vz_plot.line(ekf_timestamp, ekf_vz, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            vx_stdev_plot.line(ekf_timestamp, ekf_vx_stdev, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            vy_stdev_plot.line(ekf_timestamp, ekf_vy_stdev, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            vz_stdev_plot.line(ekf_timestamp, ekf_vz_stdev, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")

            roll_plot.line(ekf_timestamp, ekf_roll_deg, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            pitch_plot.line(ekf_timestamp, ekf_pitch_deg, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            yaw_plot.line(ekf_timestamp, ekf_yaw_deg, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")


        # ============== Misc ==============
        if is_debug_ekf_is_ego_indoors_available:
            status_plot.line(debug_ekf_is_ego_indoors_timestamp, debug_ekf_is_ego_indoors_flag, alpha=0.8, color="orange", line_width=2, legend_label="/debug/is_ekf_indoors")

        if is_ekf_is_ego_indoors_available:
            status_plot.line(ekf_is_ego_indoors_timestamp, ekf_is_ego_indoors_flag, alpha=0.8, color="cyan", line_width=2, legend_label="/is_ekf_indoors")

        if is_ins_solution_49_available:
            status_plot.line(ins_solution_49_timestamp, ins_solution_49_gnss_status, alpha=0.8, color="blue", line_width=2, legend_label="/gnss_status")


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

        vx_plot.legend.location = "top_right"
        vx_plot.legend.click_policy = "hide"
        vy_plot.legend.location = "top_right"
        vy_plot.legend.click_policy = "hide"
        vz_plot.legend.location = "top_right"
        vz_plot.legend.click_policy = "hide"
        vx_stdev_plot.legend.location = "top_right"
        vx_stdev_plot.legend.click_policy = "hide"
        vy_stdev_plot.legend.location = "top_right"
        vy_stdev_plot.legend.click_policy = "hide"
        vz_stdev_plot.legend.location = "top_right"
        vz_stdev_plot.legend.click_policy = "hide"

        roll_plot.legend.location = "top_right"
        roll_plot.legend.click_policy = "hide"
        pitch_plot.legend.location = "top_right"
        pitch_plot.legend.click_policy = "hide"
        yaw_plot.legend.location = "top_right"
        yaw_plot.legend.click_policy = "hide"

        status_plot.legend.location = "top_right"
        status_plot.legend.click_policy = "hide"

        save(gridplot([
            [x_plot, y_plot, z_plot],
            [x_stdev_plot, y_stdev_plot, z_stdev_plot],
            [vx_plot, vy_plot, vz_plot],
            [vx_stdev_plot, vy_stdev_plot, vz_stdev_plot],
            [roll_plot, pitch_plot, yaw_plot],
            [status_plot]
        ]))



        



                                                                      







