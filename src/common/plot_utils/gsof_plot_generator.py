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
import copy
import pymap3d as pm
import math
import itertools

import os


class GsofPlotGenerator:

    def __init__(self, data_parser) -> None:
        self.data_parser = data_parser

        self.ref_lla = None
        self.ref_lla_initialized = False

    def _lla2ned(self, lat_deg, lon_deg, alt_m):

        n_m, e_m, d_m = 0, 0, 0

        if not self.ref_lla_initialized:
            self.ref_lla = [lat_deg, lon_deg, alt_m]
            self.ref_lla_initialized = True
            return n_m, e_m, d_m

        e_m, n_m, d_m = pm.geodetic2enu(
            lat_deg,
            lon_deg,
            alt_m,
            self.ref_lla[0],
            self.ref_lla[1],
            self.ref_lla[2],
            ell=pm.Ellipsoid.from_name("wgs84"),
        )

        d_m = -d_m

        return n_m, e_m, d_m


    def create_bokeh_figure(self, title, x_label, y_label, width=500, height=250):
        p = figure(max_width=width, height=height,
                            title=title,
                            tools="pan,wheel_zoom,box_zoom,reset,save",
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
    
    def generate_sky_plots(self, output_file_name):

        output_file(filename=output_file_name, title="IMU Plots")
            
        constellation_ids = ["0", "2", "3", "5"]
        constellation_names = ["GPS", "GLONASS", "Galileo", "BeiDou"]

        plts = []

        for i in range(len(constellation_ids)):

            p = figure(
                title=constellation_names[i] + " Sky Plot",
                width=600,
                height=600,
                x_range=(-90, 90),
                y_range=(-90, 90),
                x_axis_type=None,
                y_axis_type=None,
                match_aspect=True,
                # min_border=0,
                tools="pan,wheel_zoom,reset,save",
                active_scroll="wheel_zoom",
            )

            # Draw horizon
            p.circle(
                x=[0],
                y=[0],
                radius=90,
                fill_alpha=0.0,
                line_color="black",
                line_dash="dashed",
            )

            # Draw elevation circles
            for elev_line in [0, 15, 30, 45, 60, 75]:
                # Convert elevation to radial distance
                r = 90 - elev_line

                p.circle(
                    x=[0],
                    y=[0],
                    radius=r,
                    fill_color=None,
                    line_color="gray",
                    line_dash="dashed",
                )

                p.text(
                    x=[0],
                    y=[r],
                    text=[f"{elev_line}°"],
                    text_baseline="top",
                    text_align="center",
                    text_color="gray",
                    text_font_size="16pt",
                )
            

            # Draw radial lines every 30 degrees
            for az_line in range(0, 361, 30):
                # Convert az_line (in degrees) to radians
                theta = math.radians(az_line)

                # We'll draw from the center (0, 0) out to the edge at r=90 (the horizon)
                r = 90

                # Convert polar -> Cartesian
                x_end = r * math.sin(theta)
                y_end = r * math.cos(theta)

                # Draw the radial line (spoke)
                p.line(
                    x=[0, x_end],
                    y=[0, y_end],
                    line_color="black",
                    line_width=1,
                    line_dash="dotted",
                )

                # Label each spoke just beyond its end
                # 1.1 * r => place label 10% beyond the circle
                x_label = 1.1 * x_end
                y_label = 1.1 * y_end

                p.text(
                    x=[x_label],
                    y=[y_label],
                    text=[f"{az_line}°"],
                    text_align="center",
                    text_baseline="middle",
                    text_font_size="16pt",
                )

            for sat_name, sat_info in self.data_parser.gnss_data[constellation_ids[i]].items():

                az = np.array(sat_info["azimuth"])
                el = np.array(sat_info["elevation"])

                # ts = np.array(sat_info["timestamp"])
                ts = np.array(np.arange(len(az)))

                # Convert to radians
                theta = np.deg2rad(az)
                # Convert elevation -> r (0 at zenith, 90 at horizon)
                r = 90 - el

                x = r * np.sin(theta)
                y = r * np.cos(theta)

                # Create a size array that grows with time
                tmin, tmax = ts.min(), ts.max()
                if tmax > tmin:
                    size_vals = 8 + 12 * (ts - tmin) / (tmax - tmin)
                else:
                    # If all timestamps are the same, just use a constant size
                    size_vals = np.full_like(ts, 10)

                # Create a ColumnDataSource so we can map color to timestamp
                source = ColumnDataSource(
                    data=dict(
                        x=x,
                        y=y,
                        size=size_vals,
                    )
                )

                # Draw a line for the trail
                p.line(x="x", y="y", line_width=2, source=source, color="gray", alpha=0.5)

                # Scatter points for each satellite
                p.circle(x="x", y="y", size="size", alpha=0.8, source=source, color="blue")

                x_label = x[-1]
                y_label = y[-1]

                p.text(
                    x=[x_label],
                    y=[y_label],
                    text=[sat_name],
                    text_color="blue",
                    text_font_size="16pt",
                    text_baseline="middle",
                    text_align="left",
                    x_offset=10,  # small shift right so text doesn't overlap the point
                    y_offset=0,
                )

            plts.append(p)

        save(gridplot([
            plts
        ]))

    def generate_dop_plots(self, output_file_name):
        output_file(filename=output_file_name, title="DOP Plots")

        a = 1

        hdop = self.data_parser.pdop_info["hdop"]
        vdop = self.data_parser.pdop_info["vdop"]

        timestamp = np.arange(len(hdop))

        p_hdop = self.create_bokeh_figure("HDOP", "Timestamp", "")
        p_vdop = self.create_bokeh_figure("VDOP", "Timestamp", "")

        p_hdop.line(timestamp, hdop, alpha=0.8, color="blue", line_width=2, legend_label="pdop_info")
        p_vdop.line(timestamp, vdop, alpha=0.8, color="blue", line_width=2, legend_label="pdop_info")

        p_hdop.legend.location = "top_right"
        p_hdop.legend.click_policy = "hide"
        p_vdop.legend.location = "top_right"
        p_vdop.legend.click_policy = "hide"   

        save(gridplot([
            [p_hdop, p_vdop],
        ]))

    def generate_noise_plots(self, output_file_name):
        output_file(filename=output_file_name, title="Noise Plots")

        p_north_stdev = self.create_bokeh_figure("North Stdev", "Timestamp", "Stdev [m]")
        p_east_stdev = self.create_bokeh_figure("East Stdev", "Timestamp", "Stdev [m]")
        p_down_stdev = self.create_bokeh_figure("Down Stdev", "Timestamp", "Stdev [m]")

        first_ts = self.data_parser.ins_solution_rms["timestamp"][0]
        last_ts = self.data_parser.ins_solution_rms["timestamp"][-1]

        p_north_stdev.line(self.data_parser.ins_solution_rms["timestamp"], self.data_parser.ins_solution_rms["north_stdev"], alpha=0.8, color="green", line_width=2, legend_label="INSSolutionRms")
        p_east_stdev.line(self.data_parser.ins_solution_rms["timestamp"], self.data_parser.ins_solution_rms["east_stdev"], alpha=0.8, color="green", line_width=2, legend_label="INSSolutionRms")
        p_down_stdev.line(self.data_parser.ins_solution_rms["timestamp"], self.data_parser.ins_solution_rms["down_stdev"], alpha=0.8, color="green", line_width=2, legend_label="INSSolutionRms")

        len_pos_sigma = len(self.data_parser.position_sigma["north_stdev"])

        self.data_parser.position_sigma["timestamp"] = np.linspace(first_ts, last_ts, len_pos_sigma)

        p_north_stdev.line(self.data_parser.position_sigma["timestamp"], self.data_parser.position_sigma["north_stdev"], alpha=0.8, color="blue", line_width=2, legend_label="PositionSigma12")
        p_east_stdev.line(self.data_parser.position_sigma["timestamp"], self.data_parser.position_sigma["east_stdev"], alpha=0.8, color="blue", line_width=2, legend_label="PositionSigma12")
        p_down_stdev.line(self.data_parser.position_sigma["timestamp"], self.data_parser.position_sigma["down_stdev"], alpha=0.8, color="blue", line_width=2, legend_label="PositionSigma12")

        len_navsat = len(self.data_parser.navsat["north_stdev"])

        self.data_parser.navsat["timestamp"] = np.linspace(first_ts, last_ts, len_navsat)

        p_north_stdev.line(self.data_parser.navsat["timestamp"], self.data_parser.navsat["north_stdev"], alpha=0.8, color="red", line_width=2, legend_label="NavSat")
        p_east_stdev.line(self.data_parser.navsat["timestamp"], self.data_parser.navsat["east_stdev"], alpha=0.8, color="red", line_width=2, legend_label="NavSat")
        p_down_stdev.line(self.data_parser.navsat["timestamp"], self.data_parser.navsat["down_stdev"], alpha=0.8, color="red", line_width=2, legend_label="NavSat")

        p_north_stdev.legend.location = "top_right"
        p_north_stdev.legend.click_policy = "hide"
        p_east_stdev.legend.location = "top_right"
        p_east_stdev.legend.click_policy = "hide"
        p_down_stdev.legend.location = "top_right"
        p_down_stdev.legend.click_policy = "hide"

        p_gnss_status = self.create_bokeh_figure("GNSS Status", "Timestamp", "Status")

        p_gnss_status.line(self.data_parser.ins_solution_rms["timestamp"], self.data_parser.ins_solution_rms["gnss_status"], alpha=0.8,  color="red", line_width=2, legend_label="INSSolutionRms")

        p_gnss_status_ins_solution = self.create_bokeh_figure("GNSS Status", "Timestamp", "Status")

        p_gnss_status_ins_solution.line(self.data_parser.ins_solution["timestamp"], self.data_parser.ins_solution["gnss_status"], alpha=0.8,  color="red", line_width=2, legend_label="INSSolution")

        save(gridplot([
            [p_north_stdev, p_east_stdev, p_down_stdev],
            [p_gnss_status, p_gnss_status_ins_solution]
        ]))


    def generate_satellite_plots(self, output_file_name):
        output_file(filename=output_file_name, title="Satellite Plots")

        p_num_sat = self.create_bokeh_figure("Number of Satellites", "Timestamp", "Num", height=600)
        num_sat = self.data_parser.gnss_data["num_of_sat"]
        num_sat_timestamp = np.arange(len(num_sat))
        p_num_sat.line(num_sat_timestamp, num_sat, alpha=0.8, color="blue", line_width=2, legend_label="")

        p_num_sat.legend.location = "top_right"
        p_num_sat.legend.click_policy = "hide"

        constellation_ids = ["0", "2", "3", "5"]
        constellation_names = ["GPS", "GLONASS", "Galileo", "BeiDou"]

        azimuth_plts = []
        elevation_plts = []
        snr_l1_plts = []
        snr_l2_plts = []
        snr_l5_plts = []
        used_in_position_plts = []
        used_in_rtk_plts = []

        for i in range(len(constellation_ids)):

            p_az = self.create_bokeh_figure(constellation_names[i] + " Azimuth", "Timestamp", "Azimuth [deg]", height=600)
            p_el = self.create_bokeh_figure(constellation_names[i] + " Elevation", "Timestamp", "Elevation [deg]", height=600)

            p_snr_l1 = self.create_bokeh_figure(constellation_names[i] + " SNR L1", "Timestamp", "[dB]", height=600)
            p_snr_l2 = self.create_bokeh_figure(constellation_names[i] + " SNR L2", "Timestamp", "[dB]", height=600)
            p_snr_l5 = self.create_bokeh_figure(constellation_names[i] + " SNR L5", "Timestamp", "[dB]", height=600)
            
            p_pos = self.create_bokeh_figure(constellation_names[i] + " Pos", "Timestamp", "Used [bool]", height=600)
            p_rtk = self.create_bokeh_figure(constellation_names[i] + " RTK", "Timestamp", "Used [bool]", height=600)

            line_colors = itertools.cycle(Category20[20])

            for (sat_name, sat_data), color in zip(self.data_parser.gnss_data[constellation_ids[i]].items(), line_colors):

                azimuth = sat_data["azimuth"]
                elevation = sat_data["elevation"]
                snr_l1 = sat_data["snr_l1"]
                snr_l2 = sat_data["snr_l2"]
                snr_l5 = sat_data["snr_l5"]
                is_used_in_position = sat_data["is_used_in_position"]
                is_used_in_rtk = sat_data["is_used_in_rtk"]

                timestamp = np.arange(len(azimuth))

                p_az.line(timestamp, azimuth, alpha=0.8, color=color, line_width=2, legend_label=sat_name)
                p_el.line(timestamp, elevation, alpha=0.8, color=color, line_width=2, legend_label=sat_name)

                p_snr_l1.line(timestamp, snr_l1, alpha=0.8, color=color, line_width=2, legend_label=sat_name)
                p_snr_l2.line(timestamp, snr_l2, alpha=0.8, color=color, line_width=2, legend_label=sat_name)
                p_snr_l5.line(timestamp, snr_l5, alpha=0.8, color=color, line_width=2, legend_label=sat_name)

                p_pos.line(timestamp, is_used_in_position, alpha=0.8, color=color, line_width=2, legend_label=sat_name)
                p_rtk.line(timestamp, is_used_in_rtk, alpha=0.8, color=color, line_width=2, legend_label=sat_name)


            p_az.legend.location = "top_right"
            p_az.legend.click_policy = "hide"
            p_el.legend.location = "top_right"
            p_el.legend.click_policy = "hide"
            p_snr_l1.legend.location = "top_right"
            p_snr_l1.legend.click_policy = "hide"
            p_snr_l2.legend.location = "top_right"
            p_snr_l2.legend.click_policy = "hide"
            p_snr_l5.legend.location = "top_right"
            p_snr_l5.legend.click_policy = "hide"
            p_pos.legend.location = "top_right"
            p_pos.legend.click_policy = "hide"
            p_rtk.legend.location = "top_right"
            p_rtk.legend.click_policy = "hide"

            azimuth_plts.append(p_az)
            elevation_plts.append(p_el)
            snr_l1_plts.append(p_snr_l1)
            snr_l2_plts.append(p_snr_l2)
            snr_l5_plts.append(p_snr_l5)
            used_in_position_plts.append(p_pos)
            used_in_rtk_plts.append(p_rtk)

        save(gridplot([
            [p_num_sat],
            azimuth_plts,
            elevation_plts,
            snr_l1_plts,
            snr_l2_plts,
            snr_l5_plts,
            used_in_position_plts,
            used_in_rtk_plts
        ]))
    
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

    def generate_gsof_plots(self, output_file_name):

        output_file(filename=output_file_name, title="GSOF Plots")

        n_plot = self.create_bokeh_figure("North Position", "Timestamp", "North [m]")
        e_plot = self.create_bokeh_figure("East Position", "Timestamp", "East [m]")
        d_plot = self.create_bokeh_figure("Down Position", "Timestamp", "Down [m]")

        is_navsat_available = False

        if "/lvx_client/navsat" in self.data_parser.mcap_data.keys():
            lvx_client_navsat = self.data_parser.mcap_data["/lvx_client/navsat"]
            lvx_client_navsat_lat = lvx_client_navsat["latitude"]
            lvx_client_navsat_lon = lvx_client_navsat["longitude"]
            lvx_client_navsat_alt = lvx_client_navsat["altitude"]

            lvx_client_navsat_position_cov = lvx_client_navsat["position_covariance"]



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

