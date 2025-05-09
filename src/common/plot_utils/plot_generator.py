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

    def create_bokeh_figure(self, title, x_label, y_label, width=500, height=250):
        p = figure(max_width=width, height=height,
                            title=title,
                            tools="pan,wheel_zoom,box_zoom,reset",
                            active_scroll="wheel_zoom")
        w = Span(
            dimension="width",
            line_dash="dashed",
            line_color="grey",
            line_width=1,
            line_alpha=0.5,
        )
        h = Span(
            dimension="height",
            line_dash="dashed",
            line_color="grey",
            line_width=1,
            line_alpha=0.5,
        )
        p.add_tools(CrosshairTool(overlay=[w, h]))
        p.xaxis.axis_label = x_label
        p.yaxis.axis_label = y_label

        return p
    
    def generate_imu_plots(self, output_file_name):
        output_file(filename=output_file_name, title="IMU Plots")

        ax_plot = self.create_bokeh_figure("X Acceleration", "Timestamp", "Ax [m^2]")
        ay_plot = self.create_bokeh_figure("Y Acceleration", "Timestamp", "Ay [m^2]")
        az_plot = self.create_bokeh_figure("Z Specific Force", "Timestamp", "Az [m/s^2]")

        wx_plot = self.create_bokeh_figure("X Angular Velocity", "Timestamp", "wx [deg/s]")
        wy_plot = self.create_bokeh_figure("Y Angular Velocity", "Timestamp", "wy [deg/s]")
        wz_plot = self.create_bokeh_figure("Z Angular Velocity", "Timestamp", "wz [deg/s]")

        is_lvx_imu_available = False
        is_microstrain_imu_available = False

        if "/lvx_client/imu_raw" in self.data_parser.mcap_data.keys():
            lvx_imu = self.data_parser.mcap_data["/lvx_client/imu_raw"]
            lvx_imu_timestamp = lvx_imu["timestamp"]

            lvx_imu_ax = lvx_imu["linear_acceleration"]["x"]
            lvx_imu_ay = lvx_imu["linear_acceleration"]["y"]
            lvx_imu_az = lvx_imu["linear_acceleration"]["z"]

            lvx_imu_wx = np.rad2deg(lvx_imu["angular_velocity"]["x"])
            lvx_imu_wy = np.rad2deg(lvx_imu["angular_velocity"]["y"])
            lvx_imu_wz = np.rad2deg(lvx_imu["angular_velocity"]["z"])

            is_lvx_imu_available = True

        if "/microstrain/imu/data" in self.data_parser.mcap_data.keys():
            ms_imu = self.data_parser.mcap_data["/microstrain/imu/data"]
            ms_imu_timestamp = ms_imu["timestamp"]
            ms_imu_ax = ms_imu["linear_acceleration"]["x"]
            ms_imu_ay = ms_imu["linear_acceleration"]["y"]
            ms_imu_az = ms_imu["linear_acceleration"]["z"]

            ms_imu_wx = np.rad2deg(ms_imu["angular_velocity"]["x"])
            ms_imu_wy = np.rad2deg(ms_imu["angular_velocity"]["y"])
            ms_imu_wz = np.rad2deg(ms_imu["angular_velocity"]["z"])

            is_microstrain_imu_available = True
        
        if (is_lvx_imu_available):
            ax_plot.line(lvx_imu_timestamp, lvx_imu_ax, alpha=0.8, color="blue", line_width=2, legend_label="/lvx/imu_raw")
            ay_plot.line(lvx_imu_timestamp, lvx_imu_ay, alpha=0.8, color="blue", line_width=2, legend_label="/lvx/imu_raw")
            az_plot.line(lvx_imu_timestamp, lvx_imu_az, alpha=0.8, color="blue", line_width=2, legend_label="/lvx/imu_raw")

            wx_plot.line(lvx_imu_timestamp, lvx_imu_wx, alpha=0.8, color="blue", line_width=2, legend_label="/lvx/imu_raw")
            wy_plot.line(lvx_imu_timestamp, lvx_imu_wy, alpha=0.8, color="blue", line_width=2, legend_label="/lvx/imu_raw")
            wz_plot.line(lvx_imu_timestamp, lvx_imu_wz, alpha=0.8, color="blue", line_width=2, legend_label="/lvx/imu_raw")

        if (is_microstrain_imu_available):
            ax_plot.line(ms_imu_timestamp, ms_imu_ax, alpha=0.8, color="red", line_width=2, legend_label="/microstrain/imu")
            ay_plot.line(ms_imu_timestamp, ms_imu_ay, alpha=0.8, color="red", line_width=2, legend_label="/microstrain/imu")
            az_plot.line(ms_imu_timestamp, ms_imu_az, alpha=0.8, color="red", line_width=2, legend_label="/microstrain/imu")

            wx_plot.line(ms_imu_timestamp, ms_imu_wx, alpha=0.8, color="red", line_width=2, legend_label="/microstrain/imu")
            wy_plot.line(ms_imu_timestamp, ms_imu_wy, alpha=0.8, color="red", line_width=2, legend_label="/microstrain/imu")
            wz_plot.line(ms_imu_timestamp, ms_imu_wz, alpha=0.8, color="red", line_width=2, legend_label="/microstrain/imu")


        ax_plot.legend.location = "top_right"
        ax_plot.legend.click_policy = "hide"
        ay_plot.legend.location = "top_right"
        ay_plot.legend.click_policy = "hide"
        az_plot.legend.location = "top_right"
        az_plot.legend.click_policy = "hide"
        wx_plot.legend.location = "top_right"
        wx_plot.legend.click_policy = "hide"
        wy_plot.legend.location = "top_right"
        wy_plot.legend.click_policy = "hide"
        wz_plot.legend.location = "top_right"
        wz_plot.legend.click_policy = "hide"

        save(gridplot([
            [ax_plot, ay_plot, az_plot],
            [wx_plot, wy_plot, wz_plot],
        ]))

    def generate_state_plots(self, output_file_name):

        output_file(filename=output_file_name, title="State Estimate Plots")

        ne_plot = self.create_bokeh_figure("Northing-Easting", "Easting [m]", "Northing [m]", width=500, height=500)

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

        localization_mode_plot = self.create_bokeh_figure("Localization Mode", "Timestamp", "Mode")

        is_input_odometry_available = False
        is_input_imu_available = False

        is_odometry_gps_available = False
        is_odometry_gps_raw_cov_available = False
        is_input_gnss_available = False

        is_pose_in_map_available = False
        is_tracked_pose_available = False
        is_cartographer_available = False
        is_debug_ekf_available = False
        is_ekf_available = False
        is_debug_ekf_is_ego_indoors_available = False
        is_ekf_is_ego_indoors_available = False
        is_ins_solution_49_available = False
        is_cartographer_convergence_status_available = False
        is_localization_mode_available = False
        is_test_available = False

        ## This is to make the initial position start at (0, 0, 0) instead of relative to the map origin
        max_offset_x = -np.inf
        max_offset_y = -np.inf
        max_offset_z = -np.inf

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

            max_offset_x = np.maximum(max_offset_x, odometry_gps_x[0])
            max_offset_y = np.maximum(max_offset_y, odometry_gps_y[0])
            max_offset_z = np.maximum(max_offset_z, odometry_gps_z[0])

        if "/odometry/gps_raw_cov" in self.data_parser.mcap_data.keys():
            odometry_gps_raw_cov = self.data_parser.mcap_data["/odometry/gps_raw_cov"]
            odometry_gps_raw_cov_timestamp = odometry_gps_raw_cov["timestamp"]
            odometry_gps_raw_cov_x = odometry_gps_raw_cov["pose"]["pose"]["position"]["x"]
            odometry_gps_raw_cov_y = odometry_gps_raw_cov["pose"]["pose"]["position"]["y"]
            odometry_gps_raw_cov_z = odometry_gps_raw_cov["pose"]["pose"]["position"]["z"]
            odometry_gps_raw_cov_x_stdev = np.sqrt(odometry_gps_raw_cov["pose"]["covariance"][:, 0])
            odometry_gps_raw_cov_y_stdev = np.sqrt(odometry_gps_raw_cov["pose"]["covariance"][:, 7])
            odometry_gps_raw_cov_z_stdev = np.sqrt(odometry_gps_raw_cov["pose"]["covariance"][:, 14])
            is_odometry_gps_raw_cov_available = True

            max_offset_x = np.maximum(max_offset_x, odometry_gps_raw_cov_x[0])
            max_offset_y = np.maximum(max_offset_y, odometry_gps_raw_cov_y[0])
            max_offset_z = np.maximum(max_offset_z, odometry_gps_raw_cov_z[0])

        if "/localization/ekf/input/gnss" in self.data_parser.mcap_data.keys():
            input_gnss = self.data_parser.mcap_data["/localization/ekf/input/gnss"]
            input_gnss_timestamp = input_gnss["timestamp"]
            input_gnss_x = input_gnss["pose"]["pose"]["position"]["x"]
            input_gnss_y = input_gnss["pose"]["pose"]["position"]["y"]
            input_gnss_z = input_gnss["pose"]["pose"]["position"]["z"]
            input_gnss_x_stdev = np.sqrt(input_gnss["pose"]["covariance"][:, 0])
            input_gnss_y_stdev = np.sqrt(input_gnss["pose"]["covariance"][:, 7])
            input_gnss_z_stdev = np.sqrt(input_gnss["pose"]["covariance"][:, 14])
            is_input_gnss_available = True

            max_offset_x = np.maximum(max_offset_x, input_gnss_x[0])
            max_offset_y = np.maximum(max_offset_y, input_gnss_y[0])
            max_offset_z = np.maximum(max_offset_z, input_gnss_z[0])

        if "/localization/ekf/input/odometry" in self.data_parser.mcap_data.keys():
            input_odometry = self.data_parser.mcap_data["/localization/ekf/input/odometry"]
            input_odometry_timestamp = input_odometry["timestamp"]
            
            input_odometry_twist_x = input_odometry["twist"]["twist"]["linear"]["x"]
            input_odometry_twist_x_stdev = input_odometry["twist"]["covariance"][:, 0]
            is_input_odometry_available = True

        if "/localization/ekf/input/imu" in self.data_parser.mcap_data.keys():
            input_imu = self.data_parser.mcap_data["/localization/ekf/input/imu"]
            input_imu_timestamp = input_imu["timestamp"]
            input_imu_orientation = input_imu["orientation"]
            
            input_imu_roll_deg = np.zeros(input_imu_timestamp.shape)
            input_imu_pitch_deg = np.zeros(input_imu_timestamp.shape)
            input_imu_yaw_deg = np.zeros(input_imu_timestamp.shape)

            for i in range(input_imu_timestamp.shape[0]):
                w = input_imu_orientation["w"][i]
                x = input_imu_orientation["x"][i]
                y = input_imu_orientation["y"][i]
                z = input_imu_orientation["z"][i]

                Rot = R.from_quat(
                    [w, x, y, z], scalar_first=True
                )

                rpy = Rot.as_euler("xyz", degrees = True)
                input_imu_roll_deg[i] = rpy[0]
                input_imu_pitch_deg[i] = rpy[1]
                input_imu_yaw_deg[i] = rpy[2]

            is_input_imu_available = True

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

            max_offset_x = np.maximum(max_offset_x, pose_in_map_x[0])
            max_offset_y = np.maximum(max_offset_y, pose_in_map_y[0])
            max_offset_z = np.maximum(max_offset_z, pose_in_map_z[0])

        # Cartographer
        if "/tracked_pose" in self.data_parser.mcap_data.keys():
            tracked_pose = self.data_parser.mcap_data["/tracked_pose"]
            tracked_pose_timestamp = tracked_pose["timestamp"]
            tracked_pose_x = tracked_pose["pose"]["position"]["x"]
            tracked_pose_y = tracked_pose["pose"]["position"]["y"]
            tracked_pose_z = tracked_pose["pose"]["position"]["z"]
            is_tracked_pose_available = True

            max_offset_x = np.maximum(max_offset_x, tracked_pose_x[0])
            max_offset_y = np.maximum(max_offset_y, tracked_pose_y[0])
            max_offset_z = np.maximum(max_offset_z, tracked_pose_z[0])

        if "/localization/ekf/input/cartographer" in self.data_parser.mcap_data.keys():
            cartographer = self.data_parser.mcap_data["/localization/ekf/input/cartographer"]
            cartographer_timestamp = cartographer["timestamp"]
            cartographer_x = cartographer["pose"]["pose"]["position"]["x"]
            cartographer_y = cartographer["pose"]["pose"]["position"]["y"]
            cartographer_z = cartographer["pose"]["pose"]["position"]["z"]

            max_offset_x = np.maximum(max_offset_x, cartographer_x[0])
            max_offset_y = np.maximum(max_offset_y, cartographer_y[0])
            max_offset_z = np.maximum(max_offset_z, cartographer_z[0])

            cartographer_orientation = cartographer["pose"]["pose"]["orientation"]


            cartographer_roll_deg = np.zeros(cartographer_timestamp.shape)
            cartographer_pitch_deg = np.zeros(cartographer_timestamp.shape)
            cartographer_yaw_deg = np.zeros(cartographer_timestamp.shape)

            for i in range(cartographer_timestamp.shape[0]):
                w = cartographer_orientation["w"][i]
                x = cartographer_orientation["x"][i]
                y = cartographer_orientation["y"][i]
                z = cartographer_orientation["z"][i]

                Rot = R.from_quat(
                    [w, x, y, z], scalar_first=True
                )

                rpy = Rot.as_euler("xyz", degrees = True)
                cartographer_roll_deg[i] = rpy[0]
                cartographer_pitch_deg[i] = rpy[1]
                cartographer_yaw_deg[i] = rpy[2]

            is_cartographer_available = True
            

        if "/debug/localization/kinematic_state" in self.data_parser.mcap_data.keys():
            debug_ekf = self.data_parser.mcap_data["/debug/localization/kinematic_state"]
            debug_ekf_timestamp = debug_ekf["timestamp"][1:]
            debug_ekf_x = debug_ekf["pose"]["pose"]["position"]["x"][1:]
            debug_ekf_y = debug_ekf["pose"]["pose"]["position"]["y"][1:]
            debug_ekf_z = debug_ekf["pose"]["pose"]["position"]["z"][1:]
            debug_ekf_x_stdev = np.sqrt(debug_ekf["pose"]["covariance"][1:, 0])
            debug_ekf_y_stdev = np.sqrt(debug_ekf["pose"]["covariance"][1:, 7])
            debug_ekf_z_stdev = np.sqrt(debug_ekf["pose"]["covariance"][1:, 14])

            max_offset_x = np.maximum(max_offset_x, debug_ekf_x[0])
            max_offset_y = np.maximum(max_offset_y, debug_ekf_y[0])
            max_offset_z = np.maximum(max_offset_z, debug_ekf_z[0])

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

            max_offset_x = np.maximum(max_offset_x, ekf_x[0])
            max_offset_y = np.maximum(max_offset_y, ekf_y[0])
            max_offset_z = np.maximum(max_offset_z, ekf_z[0])

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

        try:
            if "/localization/debug/cartographer_convergence_status" in self.data_parser.mcap_data.keys():
                cartographer_convergence_status = self.data_parser.mcap_data["/localization/debug/cartographer_convergence_status"]
                cartographer_convergence_status_timestamp = cartographer_convergence_status["timestamp"]
                is_cartographer_reliable = cartographer_convergence_status["is_cartographer_reliable"]

                is_cartographer_convergence_status_available = True
        except:
            pass

        if "/output/test_cartographer_convergence_status_from_gnss" in self.data_parser.mcap_data.keys():
            test_cartographer_convergence_status_from_gnss = self.data_parser.mcap_data["/output/test_cartographer_convergence_status_from_gnss"]
            test_cartographer_convergence_status_from_gnss_timestamp = test_cartographer_convergence_status_from_gnss["timestamp"]
            test_cartographer_convergence_status_from_gnss_is_reliable = test_cartographer_convergence_status_from_gnss["is_cartographer_reliable"]

            is_test_available = True

        if "/localization/debug/localization_mode" in self.data_parser.mcap_data.keys():
            localization_mode = self.data_parser.mcap_data["/localization/debug/localization_mode"]
            localization_mode_timestamp = localization_mode["timestamp"]
            localization_mode_mode = localization_mode["mode"]
            is_localization_mode_available = True


        # ============== Odometry ================
        if is_input_odometry_available:
            vx_plot.line(input_odometry_timestamp, input_odometry_twist_x, alpha=0.8, color="grey", line_width=2, legend_label="/input/odometry")

        if is_input_imu_available:
            roll_plot.line(input_imu_timestamp, input_imu_roll_deg, alpha=0.8, color="grey", line_width=2, legend_label="/input/imu")
            pitch_plot.line(input_imu_timestamp, input_imu_pitch_deg, alpha=0.8, color="grey", line_width=2, legend_label="/input/imu")
            yaw_plot.line(input_imu_timestamp, input_imu_yaw_deg, alpha=0.8, color="grey", line_width=2, legend_label="/input/imu")

        # ============== GPS ==============
        if is_input_gnss_available:

            ne_plot.circle(input_gnss_y, input_gnss_x, alpha=0.8, radius=1, radius_units="screen", color="chartreuse", legend_label="/input/gnss")

            # x_plot.line(input_gnss_timestamp, input_gnss_x - max_offset_x, alpha=0.8, color="chartreuse", line_width=2, legend_label="/input/gnss")
            # y_plot.line(input_gnss_timestamp, input_gnss_y - max_offset_y, alpha=0.8, color="chartreuse", line_width=2, legend_label="/input/gnss")
            # z_plot.line(input_gnss_timestamp, input_gnss_z - max_offset_z, alpha=0.8, color="chartreuse", line_width=2, legend_label="/input/gnss")

            x_plot.circle(input_gnss_timestamp, input_gnss_x - max_offset_x, alpha=0.8, radius=1, radius_units="screen", color="chartreuse", legend_label="/input/gnss")
            y_plot.circle(input_gnss_timestamp, input_gnss_y - max_offset_y, alpha=0.8, radius=1, radius_units="screen", color="chartreuse", legend_label="/input/gnss")
            z_plot.circle(input_gnss_timestamp, input_gnss_z - max_offset_z, alpha=0.8, radius=1, radius_units="screen", color="chartreuse", legend_label="/input/gnss")

            x_stdev_plot.line(input_gnss_timestamp, input_gnss_x_stdev, alpha=0.8, color="chartreuse", line_width=2, legend_label="/input/gnss")
            y_stdev_plot.line(input_gnss_timestamp, input_gnss_y_stdev, alpha=0.8, color="chartreuse", line_width=2, legend_label="/input/gnss")
            z_stdev_plot.line(input_gnss_timestamp, input_gnss_z_stdev, alpha=0.8, color="chartreuse", line_width=2, 
            legend_label="/input/gnss")

        if is_odometry_gps_available:

            ne_plot.circle(odometry_gps_y, odometry_gps_x, alpha=0.8, radius=1, radius_units="screen", color="blueviolet", legend_label="/odometry/gps")

            x_plot.line(odometry_gps_timestamp, odometry_gps_x - max_offset_x, alpha=0.8, color="blueviolet", line_width=2, legend_label="/odometry/gps")
            y_plot.line(odometry_gps_timestamp, odometry_gps_y - max_offset_y, alpha=0.8, color="blueviolet", line_width=2, legend_label="/odometry/gps")
            z_plot.line(odometry_gps_timestamp, odometry_gps_z - max_offset_z, alpha=0.8, color="blueviolet", line_width=2, legend_label="/odometry/gps")
            x_stdev_plot.line(odometry_gps_timestamp, odometry_gps_x_stdev, alpha=0.8, color="blueviolet", line_width=2, legend_label="/odometry/gps")
            y_stdev_plot.line(odometry_gps_timestamp, odometry_gps_y_stdev, alpha=0.8, color="blueviolet", line_width=2, legend_label="/odometry/gps")
            z_stdev_plot.line(odometry_gps_timestamp, odometry_gps_z_stdev, alpha=0.8, color="blueviolet", line_width=2, legend_label="/odometry/gps")

        if is_odometry_gps_raw_cov_available:

            ne_plot.circle(odometry_gps_raw_cov_y, odometry_gps_raw_cov_x, alpha=0.8, radius=1, radius_units="screen", color="blue", legend_label="/odometry/gps_raw_cov")

            x_plot.line(odometry_gps_raw_cov_timestamp, odometry_gps_raw_cov_x - max_offset_x, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps_raw_cov")
            y_plot.line(odometry_gps_raw_cov_timestamp, odometry_gps_raw_cov_y - max_offset_y, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps_raw_cov")
            z_plot.line(odometry_gps_raw_cov_timestamp, odometry_gps_raw_cov_z - max_offset_z, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps_raw_cov")
            x_stdev_plot.line(odometry_gps_raw_cov_timestamp, odometry_gps_raw_cov_x_stdev, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps_raw_cov")
            y_stdev_plot.line(odometry_gps_raw_cov_timestamp, odometry_gps_raw_cov_y_stdev, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps_raw_cov")
            z_stdev_plot.line(odometry_gps_raw_cov_timestamp, odometry_gps_raw_cov_z_stdev, alpha=0.8, color="blue", line_width=2, legend_label="/odometry/gps_raw_cov")            

        # ============== Pose in Map ==============
        if is_pose_in_map_available:

            ne_plot.circle(pose_in_map_y, pose_in_map_x, alpha=0.8, radius=1, radius_units="screen", color="green", legend_label="/odometry/pose_in_map")

            x_plot.line(pose_in_map_timestamp, pose_in_map_x - max_offset_x, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            y_plot.line(pose_in_map_timestamp, pose_in_map_y - max_offset_y, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            z_plot.line(pose_in_map_timestamp, pose_in_map_z - max_offset_z, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            x_stdev_plot.line(pose_in_map_timestamp, pose_in_map_x_stdev, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            y_stdev_plot.line(pose_in_map_timestamp, pose_in_map_y_stdev, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")
            z_stdev_plot.line(pose_in_map_timestamp, pose_in_map_z_stdev, alpha=0.8, color="green", line_width=2, legend_label="/odometry/pose_in_map")

        # ============== Cartographer ==============
        if is_cartographer_available:

            ne_plot.circle(cartographer_y, cartographer_x, alpha=0.8, radius=1, radius_units="screen", color="red", legend_label="/input/cartographer")

            x_plot.line(cartographer_timestamp, cartographer_x - max_offset_x, alpha=0.8, color="red", line_width=2, legend_label="/input/cartographer")
            y_plot.line(cartographer_timestamp, cartographer_y - max_offset_y, alpha=0.8, color="red", line_width=2, legend_label="/input/cartographer")
            z_plot.line(cartographer_timestamp, cartographer_z - max_offset_z, alpha=0.8, color="red", line_width=2, legend_label="/input/cartographer")

            roll_plot.line(cartographer_timestamp, cartographer_roll_deg, alpha=0.8, color="red", line_width=2, legend_label="/input/cartographer")
            pitch_plot.line(cartographer_timestamp, cartographer_pitch_deg, alpha=0.8, color="red", line_width=2, legend_label="/input/cartographer")
            yaw_plot.line(cartographer_timestamp, cartographer_yaw_deg, alpha=0.8, color="red", line_width=2, legend_label="/input/cartographer")

        if is_tracked_pose_available:

            ne_plot.circle(tracked_pose_y, tracked_pose_x, alpha=0.8, radius=1, radius_units="screen", color="maroon", legend_label="/tracked_pose")

            x_plot.line(tracked_pose_timestamp, tracked_pose_x - max_offset_x, alpha=0.8, color="maroon", line_width=2, legend_label="/tracked_pose")
            y_plot.line(tracked_pose_timestamp, tracked_pose_y - max_offset_y, alpha=0.8, color="maroon", line_width=2, legend_label="/tracked_pose")
            z_plot.line(tracked_pose_timestamp, tracked_pose_z - max_offset_z, alpha=0.8, color="maroon", line_width=2, legend_label="/tracked_pose")

        # ============== Replay EKF ==============
        if is_debug_ekf_available:

            ne_plot.circle(debug_ekf_y, debug_ekf_x, alpha=0.8, radius=1, radius_units="screen", color="orange", legend_label="/debug/ekf")

            x_plot.line(debug_ekf_timestamp, debug_ekf_x - max_offset_x, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            y_plot.line(debug_ekf_timestamp, debug_ekf_y - max_offset_y, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
            z_plot.line(debug_ekf_timestamp, debug_ekf_z - max_offset_z, alpha=0.8, color="orange", line_width=2, legend_label="/debug/ekf")
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
            ne_plot.circle(ekf_y, ekf_x, alpha=0.8, radius=1, radius_units="screen", color="orange", legend_label="/ekf")

            x_plot.line(ekf_timestamp, ekf_x - max_offset_x, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            y_plot.line(ekf_timestamp, ekf_y - max_offset_y, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
            z_plot.line(ekf_timestamp, ekf_z - max_offset_z, alpha=0.8, color="cyan", line_width=2, legend_label="/ekf")
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

        if is_cartographer_convergence_status_available:
            status_plot.line(cartographer_convergence_status_timestamp, is_cartographer_reliable, alpha=0.8, color="red", line_width=2, legend_label="/cartographer_convergence_status")

        if is_test_available:
            status_plot.line(test_cartographer_convergence_status_from_gnss_timestamp, test_cartographer_convergence_status_from_gnss_is_reliable, alpha=0.8, color="purple", line_width=2, legend_label="/gnss_cartographer_convergence_status")

        if is_localization_mode_available:
            localization_mode_plot.line(localization_mode_timestamp, localization_mode_mode, alpha=0.8, color="blue", line_width=2, legend_label="/localization_mode")

        # Post mumbo jumbo 

        ne_plot.legend.location = "top_right"
        ne_plot.legend.click_policy = "hide"

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

        localization_mode_plot.legend.location = "top_right"
        localization_mode_plot.legend.click_policy = "hide"

        save(gridplot([
            [ne_plot],
            [x_plot, y_plot, z_plot],
            [x_stdev_plot, y_stdev_plot, z_stdev_plot],
            [vx_plot, vy_plot, vz_plot],
            [vx_stdev_plot, vy_stdev_plot, vz_stdev_plot],
            [roll_plot, pitch_plot, yaw_plot],
            [localization_mode_plot, status_plot]
        ]))

    def generate_nis_metrics_plots(self, output_file_name):

        output_file(filename=output_file_name, title="NIS Metrics Plots")

        num_pos_meas_plot = self.create_bokeh_figure("Position Number of Measurements Fused", "Timestamp", "Obs X [m]")

        x_obs_plot = self.create_bokeh_figure("X Position Observations", "Timestamp", "Obs X [m]")
        y_obs_plot = self.create_bokeh_figure("Y Position Observations", "Timestamp", "Obs Y [m]")
        z_obs_plot = self.create_bokeh_figure("Z Position Observations", "Timestamp", "Obs Z [m]")

        x_inn_plot = self.create_bokeh_figure("X Position Innovations", "Timestamp", "Inn X [m]")
        y_inn_plot = self.create_bokeh_figure("Y Position Innovations", "Timestamp", "Inn Y [m]")
        z_inn_plot = self.create_bokeh_figure("Z Position Innovations", "Timestamp", "Inn Z [m]")

        pos_mahalanobis_dist_plot = self.create_bokeh_figure("Position Mahalanobis Distances", "Timestamp", "d [m]")
        pos_test_statistics_plot = self.create_bokeh_figure("Position Test Statistics", "Timestamp", "Bool")

        if "/localization/ekf/debug/estimator_aid_src_cartographer_pose" in self.data_parser.mcap_data.keys():

            estimator_aid_src_cartographer_pose = self.data_parser.mcap_data["/localization/ekf/debug/estimator_aid_src_cartographer_pose"]
            
            timestamps = estimator_aid_src_cartographer_pose["timestamp"]
            observations = estimator_aid_src_cartographer_pose["observation"]
            innovations = estimator_aid_src_cartographer_pose["innovation"]
            mahalanobis_distances = estimator_aid_src_cartographer_pose["mahalanobis_distance"]
            mahalanobis_thresholds = estimator_aid_src_cartographer_pose["mahalanobis_threshold"]
            is_outliers = estimator_aid_src_cartographer_pose["is_outlier"]
            innovation_rejecteds = estimator_aid_src_cartographer_pose["innovation_rejected"]

            # Batching arrays
            num_pos_meas_times = []
            num_pos_meas_counts = []

            x_obs_times, x_obs_values = [], []
            y_obs_times, y_obs_values = [], []
            z_obs_times, z_obs_values = [], []

            x_inn_times, x_inn_values = [], []
            y_inn_times, y_inn_values = [], []
            z_inn_times, z_inn_values = [], []

            maha_times, maha_values, maha_thresh_values = [], [], []

            test_stat_times, test_stat_is_outlier, test_stat_inn_rejected = [], [], []

            # Loop over each filter update
            for timestamp, obs_list, inn_list, d_list, thresh_list, outlier_list, reject_list in zip(
                timestamps, observations, innovations, mahalanobis_distances, mahalanobis_thresholds, is_outliers, innovation_rejecteds
            ):
                # Number of measurements per timestamp
                num_pos_meas_times.append(timestamp)
                num_pos_meas_counts.append(len(obs_list))

                for obs, inn, d, thresh, outlier, reject in zip(obs_list, inn_list, d_list, thresh_list, outlier_list, reject_list):
                    # Observations
                    x_obs_times.append(timestamp)
                    x_obs_values.append(obs[0])

                    y_obs_times.append(timestamp)
                    y_obs_values.append(obs[1])

                    z_obs_times.append(timestamp)
                    z_obs_values.append(obs[2])

                    # Innovations
                    x_inn_times.append(timestamp)
                    x_inn_values.append(inn[0])

                    y_inn_times.append(timestamp)
                    y_inn_values.append(inn[1])

                    z_inn_times.append(timestamp)
                    z_inn_values.append(inn[2])

                    # Mahalanobis distance
                    maha_times.append(timestamp)
                    maha_values.append(d)
                    maha_thresh_values.append(thresh)

                    # Test statistics
                    test_stat_times.append(timestamp)
                    test_stat_is_outlier.append(int(outlier))
                    test_stat_inn_rejected.append(int(reject))

            # Scatter all batched points at once
            num_pos_meas_plot.scatter(num_pos_meas_times, num_pos_meas_counts, size=3, color="blue", legend_label="cartographer_pose")

            x_obs_plot.scatter(x_obs_times, x_obs_values, size=3, color="blue", legend_label="cartographer_pose")
            y_obs_plot.scatter(y_obs_times, y_obs_values, size=3, color="blue", legend_label="cartographer_pose")
            z_obs_plot.scatter(z_obs_times, z_obs_values, size=3, color="blue", legend_label="cartographer_pose")

            x_inn_plot.scatter(x_inn_times, x_inn_values, size=6, color="blue", legend_label="cartographer_pose")
            y_inn_plot.scatter(y_inn_times, y_inn_values, size=6, color="blue", legend_label="cartographer_pose")
            z_inn_plot.scatter(z_inn_times, z_inn_values, size=6, color="blue", legend_label="cartographer_pose")

            pos_mahalanobis_dist_plot.scatter(maha_times, maha_values, size=3, color="blue", legend_label="cartographer_pose")
            pos_mahalanobis_dist_plot.scatter(maha_times, maha_thresh_values, size=3, color="blue", marker="x", legend_label="cartographer_pose_threshold")

            pos_test_statistics_plot.scatter(test_stat_times, test_stat_is_outlier, size=6, color="blue", marker="circle", legend_label="cartographer_pose")
            pos_test_statistics_plot.scatter(test_stat_times, test_stat_inn_rejected, size=6, color="blue", marker="x", legend_label="cartographer_pose Rejected")


            x_obs_plot.legend.location = "top_right"
            x_obs_plot.legend.click_policy = "hide"
            y_obs_plot.legend.location = "top_right"
            y_obs_plot.legend.click_policy = "hide"
            z_obs_plot.legend.location = "top_right"
            z_obs_plot.legend.click_policy = "hide"

            x_inn_plot.legend.location = "top_right"
            x_inn_plot.legend.click_policy = "hide"
            y_inn_plot.legend.location = "top_right"
            y_inn_plot.legend.click_policy = "hide"
            z_inn_plot.legend.location = "top_right"
            z_inn_plot.legend.click_policy = "hide"

            pos_mahalanobis_dist_plot.legend.location = "top_right"
            pos_mahalanobis_dist_plot.legend.click_policy = "hide"
            pos_test_statistics_plot.legend.location = "top_right"
            pos_test_statistics_plot.legend.click_policy = "hide"
            num_pos_meas_plot.legend.location = "top_right"
            num_pos_meas_plot.legend.click_policy = "hide"



        save(gridplot([
            [num_pos_meas_plot],
            [x_obs_plot, x_inn_plot],
            [y_obs_plot, y_inn_plot],
            [z_obs_plot, z_inn_plot],
            [pos_mahalanobis_dist_plot, pos_test_statistics_plot]
        ]))

