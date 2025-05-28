from bokeh.plotting import figure, output_file, save
from bokeh.layouts import gridplot
from bokeh.models import Label, CrosshairTool, Span, Circle, BoxAnnotation
from bokeh.palettes import Category20, Category10
from bokeh.io import show
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, Div
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn

from scipy.spatial.transform import Rotation as R

import numpy as np
import copy
import pymap3d as pm
import math
import itertools
import os


class GsofPlotGenerator:
    # Plot dimensions
    DEFAULT_WIDTH = 500
    DEFAULT_HEIGHT = 250
    SQUARE_PLOT_SIZE = 500
    TALL_PLOT_HEIGHT = 600
    EXTRA_TALL_PLOT_HEIGHT = 400
    
    # Constellation mapping
    CONSTELLATION_IDS = ["0", "2", "3", "5"]
    CONSTELLATION_NAMES = ["GPS", "GLONASS", "Galileo", "BeiDou"]
    
    # Sky plot parameters
    SKY_PLOT_SIZE = 600
    HORIZON_RADIUS = 90
    ELEVATION_CIRCLES = [0, 15, 30, 45, 60, 75]
    AZIMUTH_LINES = list(range(0, 361, 30))

    def __init__(self, data_parser) -> None:
        self.data_parser = data_parser

        self.ref_lla = None
        self.ref_lla_initialized = False

        # There are a lot of GSOF messages with no proper header timestamp
        # Interpolate from the INS solution message which is safe
        # to assume that it will always exist
        self.first_ts = self.data_parser.ins_solution["timestamp"][0]
        self.last_ts = self.data_parser.ins_solution["timestamp"][1]

    def lla2ned(self, lat_deg, lon_deg, alt_m):

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


    def create_bokeh_figure(self, title, x_label, y_label, width=None, height=None):
        width = width or self.DEFAULT_WIDTH
        height = height or self.DEFAULT_HEIGHT
        
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
    
    def generate_pvt_plots(self, output_file_name):

        output_file(filename=output_file_name, title="PVT Plots")

        ne_plot = self.create_bokeh_figure("Northing-Easting", "Easting [m]", "Northing [m]", width=self.SQUARE_PLOT_SIZE, height=self.SQUARE_PLOT_SIZE)

        x_plot = self.create_bokeh_figure("North Position", "Timestamp", "X [m]")
        y_plot = self.create_bokeh_figure("East Position", "Timestamp", "Y [m]")
        z_plot = self.create_bokeh_figure("Down Position", "Timestamp", "Z [m]")

        gnss_status_plot = self.create_bokeh_figure("GNSS Status\n0: No Fix, 1: GPS SPS, 2: DGNSS\n3: GPS PPS, 4: RTK Fixed, 5: RTK Float\n6: Dead Reckoning ", "Timestamp", "Status")

        x_stdev_plot = self.create_bokeh_figure("X Position Stdev", "Timestamp", "Stdev [m]")
        y_stdev_plot = self.create_bokeh_figure("Y Position Stdev", "Timestamp", "Stdev [m]")
        z_stdev_plot = self.create_bokeh_figure("Z Position Stdev", "Timestamp", "Stdev [m]")

        vx_plot = self.create_bokeh_figure("X Velocity", "Timestamp", "Vx [m]")
        vy_plot = self.create_bokeh_figure("Y Velocity", "Timestamp", "Vy [m]")
        vz_plot = self.create_bokeh_figure("Z Velocity", "Timestamp", "Vz [m]")

        vx_stdev_plot = self.create_bokeh_figure("X Velocity Stdev", "Timestamp", "Stdev [m/s]")
        vy_stdev_plot = self.create_bokeh_figure("Y Velocity Stdev", "Timestamp", "Stdev [m/s]")
        vz_stdev_plot = self.create_bokeh_figure("Z Velocity Stdev", "Timestamp", "Stdev [m/s]")

        is_code_pvt_available = False

        ins_solution_timestamp = self.data_parser.ins_solution["timestamp"]
        ins_solution_gnss_status = self.data_parser.ins_solution["gnss_status"]
        ins_solution_lat_deg = self.data_parser.ins_solution["latitude"]
        ins_solution_lon_deg = self.data_parser.ins_solution["longitude"]
        ins_solution_lat_height = self.data_parser.ins_solution["altitude"]
        ins_solution_north = []
        ins_solution_east = []
        ins_solution_down = []

        for lat, lon, alt in zip(ins_solution_lat_deg, ins_solution_lon_deg, ins_solution_lat_height):
            n_m, e_m, d_m = self.lla2ned(lat, lon, alt)
            ins_solution_north.append(n_m)
            ins_solution_east.append(e_m)
            ins_solution_down.append(d_m)


        if hasattr(self.data_parser, "code_lat_lon_ht"):
            code_lat_long_ht = self.data_parser.code_lat_lon_ht
            code_lat_long_ht_timestamp = code_lat_long_ht["timestamp"]
            code_lat_long_lat_deg = code_lat_long_ht["lat_deg"]
            code_lat_long_lon_deg = code_lat_long_ht["lon_deg"]
            code_lat_long_height = code_lat_long_ht["height"]
            code_lat_lon_position_type = code_lat_long_ht["position_type"]

            code_lat_long_north = []
            code_lat_long_east = []
            code_lat_long_down = []

            for lat, lon, alt in zip(code_lat_long_lat_deg, code_lat_long_lon_deg, code_lat_long_height):
                n_m, e_m, d_m = self.lla2ned(lat, lon, alt)
                code_lat_long_north.append(n_m)
                code_lat_long_east.append(e_m)
                code_lat_long_down.append(d_m)

            is_code_pvt_available = True

        
        if (is_code_pvt_available):
            init_timestamp = np.minimum(ins_solution_timestamp[0], code_lat_long_ht_timestamp[0])
            # Subtract the first timestamp from both
            ins_solution_timestamp -= init_timestamp
            code_lat_long_ht_timestamp -= init_timestamp
        else:
            ins_solution_timestamp -= ins_solution_timestamp[0]
        

        ne_plot.circle(ins_solution_east, ins_solution_north, alpha=0.8, radius=1, radius_units="screen", color="red", legend_label="/gsof/ins_solution")

        x_plot.circle(ins_solution_timestamp, ins_solution_north, alpha=0.8, radius=1, radius_units="screen", color="red", legend_label="/gsof/ins_solution")

        y_plot.circle(ins_solution_timestamp, ins_solution_east, alpha=0.8, radius=1, radius_units="screen", color="red", legend_label="/gsof/ins_solution")

        z_plot.circle(ins_solution_timestamp, ins_solution_down, alpha=0.8, radius=1, radius_units="screen", color="red", legend_label="/gsof/ins_solution")

        gnss_status_plot.circle(ins_solution_timestamp, ins_solution_gnss_status, alpha=0.8, radius=1, radius_units="screen", color="red", legend_label="/gsof/ins_solution")

        if (is_code_pvt_available):
            ne_plot.circle(code_lat_long_east, code_lat_long_north, alpha=0.8, radius=1, radius_units="screen", color="blue", legend_label="/gsof/code_lat_lon_ht")
        
            x_plot.circle(code_lat_long_ht_timestamp, code_lat_long_north, alpha=0.8, radius=1, radius_units="screen", color="blue", legend_label="/gsof/code_lat_lon_ht")

            y_plot.circle(code_lat_long_ht_timestamp, code_lat_long_east, alpha=0.8, radius=1, radius_units="screen", color="blue", legend_label="/gsof/code_lat_lon_ht")

            z_plot.circle(code_lat_long_ht_timestamp, code_lat_long_down, alpha=0.8, radius=1, radius_units="screen", color="blue", legend_label="/gsof/code_lat_lon_ht")

            gnss_status_plot.circle(code_lat_long_ht_timestamp, code_lat_lon_position_type, alpha=0.8, radius=1, radius_units="screen", color="blue", legend_label="/gsof/code_lat_lon_ht")


        ne_plot.legend.location = "top_right"
        ne_plot.legend.click_policy = "hide"

        x_plot.legend.location = "top_right"
        x_plot.legend.click_policy = "hide"   
        y_plot.legend.location = "top_right"
        y_plot.legend.click_policy = "hide"  
        z_plot.legend.location = "top_right"
        z_plot.legend.click_policy = "hide"  

        gnss_status_plot.legend.location = "top_right"
        gnss_status_plot.legend.click_policy = "hide"  

        save(gridplot([
            [ne_plot],
            [x_plot, y_plot, z_plot],
            [gnss_status_plot]
        ]))

        return

    
    def generate_rtk_status_plots(self, output_file_name):

        output_file(filename=output_file_name, title="RTK Status")

        num_of_svs_plot = self.create_bokeh_figure("Number of Satellites Used", "Timestamp", "Num")
        init_num_plot = self.create_bokeh_figure("Number of Inits", "Timestamp", "Num")
        position_flags_1_plot = self.create_bokeh_figure("Position Flags 1", "Timestamp", "Bool", height=self.EXTRA_TALL_PLOT_HEIGHT)
        position_flags_2_plot = self.create_bokeh_figure("Position Flags 2", "Timestamp", "Bool", height=self.EXTRA_TALL_PLOT_HEIGHT)

        network_solution_plot = self.create_bokeh_figure("Network Solution", "Timestamp", "Bool")
        correction_age_plot = self.create_bokeh_figure("Correction Age", "Timestamp", "Age [s]")
        rtk_fix_plot = self.create_bokeh_figure("RTK Fix", "Timestamp", "Bool")
        rtk_condition_plot = self.create_bokeh_figure("RTK Condition", "Timestamp", "Bool", height=self.EXTRA_TALL_PLOT_HEIGHT)
        position_fix_type_plot = self.create_bokeh_figure("Position Fix Type", "Timestamp", "Fix Type")

        base_valid_plot = self.create_bokeh_figure("Reference Station Valid", "Timestamp", "Bool")
        base_id_plot = self.create_bokeh_figure("Reference Station ID", "Timestamp", "ID")
        base_lat_plot = self.create_bokeh_figure("Reference Station Latitude", "Timestamp", "Latitude [deg]")
        base_lon_plot = self.create_bokeh_figure("Reference Station Latitude", "Timestamp", "Longitude [deg]")
        base_alt_plot = self.create_bokeh_figure("Reference Station Altitude WGS84", "Timestamp", "Altitude [m]")

        base_north_plot = self.create_bokeh_figure("Reference Station North", "Timestamp", "North [m]")
        base_east_plot = self.create_bokeh_figure("Reference Station East", "Timestamp", "East [m]")
        base_down_plot = self.create_bokeh_figure("Reference Station Down", "Timestamp", "Down [m]")

        is_position_time_info_available = False
        is_position_type_info_available = False
        is_received_base_info_available = False

        if hasattr(self.data_parser, "position_time_info"):

            position_time_info = self.data_parser.position_time_info
            position_time_info_timestamp = position_time_info["timestamp"]
            position_time_info_num_svs = position_time_info["number_space_vehicles_used"]
            position_time_info_init_num = position_time_info["init_num"]

            new_position = position_time_info["new_position"]
            clock_fix = position_time_info["clock_fix"]
            horizontal_coords = position_time_info["horizontal_coords"]
            height_calculated = position_time_info["height_calculated"]
            least_squares = position_time_info["least_squares"]
            filtered_L1 = position_time_info["filtered_L1"]

            is_differential = position_time_info["is_differential"]
            uses_phase = position_time_info["uses_phase"]
            is_RTK_fixed = position_time_info["is_RTK_fixed"]
            omniSTAR = position_time_info["omniSTAR"]
            static_constraint = position_time_info["static_constraint"]
            network_RTK = position_time_info["network_RTK"]
            dithered_RTK = position_time_info["dithered_RTK"]
            beacon_DGNSS = position_time_info["beacon_DGNSS"]

            is_position_time_info_available = True

        if hasattr(self.data_parser, "position_type_info"):
            position_type_info = self.data_parser.position_type_info
            position_type_info_timestamp = position_type_info["timestamp"]
            position_type_info_is_network_solution = position_type_info["is_network_solution"]
            position_type_info_is_rtk_fix = position_type_info["is_rtk_fix"]

            position_type_info_correction_age = position_type_info["correction_age"]
            position_type_info_position_fix_type = position_type_info["position_fix_type"]

            rtk_cond_new_position_computed = position_type_info["rtk_cond_new_position_computed"]
            rtk_cond_no_synced_pair = position_type_info["rtk_cond_no_synced_pair"]
            rtk_cond_insuff_dd_meas = position_type_info["rtk_cond_insuff_dd_meas"]
            rtk_cond_ref_pos_unavailable = position_type_info["rtk_cond_ref_pos_unavailable"]
            rtk_cond_failed_integer_ver_with_fix_sol = position_type_info["rtk_cond_failed_integer_ver_with_fix_sol"]
            rtk_cond_sol_res_rms_exceeds = position_type_info["rtk_cond_sol_res_rms_exceeds"]
            rtk_cond_pdop_exceeds = position_type_info["rtk_cond_pdop_exceeds"]

            is_position_type_info_available = True
        
        if hasattr(self.data_parser, "received_base_info"):
            received_base_info = self.data_parser.received_base_info
            received_base_info_timestamp = received_base_info["timestamp"]

            base_valid = received_base_info["base_valid"]
            base_id = received_base_info["base_id"]

            base_lat_deg = received_base_info["base_lat_deg"]
            base_lon_deg = received_base_info["base_lon_deg"]
            base_height_m = received_base_info["base_height_m"]

            base_north_m = []
            base_east_m = []
            base_down_m = []

            for lat, lon, alt in zip(base_lat_deg, base_lon_deg, base_height_m):
                n_m, e_m, d_m = self.lla2ned(lat, lon, alt)

                base_north_m.append(n_m)
                base_east_m.append(e_m)
                base_down_m.append(d_m)

            is_received_base_info_available = True

        if (is_position_time_info_available):

            num_of_svs_plot.line(position_time_info_timestamp, position_time_info_num_svs, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/position_time_info")
            init_num_plot.line(position_time_info_timestamp, position_time_info_init_num, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/position_time_info")

            position_flags_1_plot.line(position_time_info_timestamp, new_position, alpha=0.8, color="blue", line_width=2, legend_label="new_position")
            position_flags_1_plot.line(position_time_info_timestamp, clock_fix, alpha=0.8, color="red", line_width=2, legend_label="clock_fix")
            position_flags_1_plot.line(position_time_info_timestamp, horizontal_coords, alpha=0.8, color="green", line_width=2, legend_label="horizontal_coords")
            position_flags_1_plot.line(position_time_info_timestamp, height_calculated, alpha=0.8, color="orange", line_width=2, legend_label="height_calculated")
            position_flags_1_plot.line(position_time_info_timestamp, least_squares, alpha=0.8, color="purple", line_width=2, legend_label="least_squares")
            position_flags_1_plot.line(position_time_info_timestamp, filtered_L1, alpha=0.8, color="cyan", line_width=2, legend_label="filtered_L1")

            position_flags_2_plot.line(position_time_info_timestamp, is_differential, alpha=0.8, color="blue", line_width=2, legend_label="is_differential")
            position_flags_2_plot.line(position_time_info_timestamp, uses_phase, alpha=0.8, color="red", line_width=2, legend_label="uses_phase")
            position_flags_2_plot.line(position_time_info_timestamp, is_RTK_fixed, alpha=0.8, color="green", line_width=2, legend_label="is_RTK_fixed")
            position_flags_2_plot.line(position_time_info_timestamp, omniSTAR, alpha=0.8, color="orange", line_width=2, legend_label="omniSTAR")
            position_flags_2_plot.line(position_time_info_timestamp, static_constraint, alpha=0.8, color="purple", line_width=2, legend_label="static_constraint")
            position_flags_2_plot.line(position_time_info_timestamp, network_RTK, alpha=0.8, color="cyan", line_width=2, legend_label="network_RTK")
            position_flags_2_plot.line(position_time_info_timestamp, dithered_RTK, alpha=0.8, color="brown", line_width=2, legend_label="dithered_RTK")
            position_flags_2_plot.line(position_time_info_timestamp, beacon_DGNSS, alpha=0.8, color="magenta", line_width=2, legend_label="beacon_DGNSS")


        if (is_position_type_info_available):
            correction_age_plot.line(position_type_info_timestamp, position_type_info_correction_age, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/position_type_info")
            rtk_fix_plot.line(position_type_info_timestamp, position_type_info_is_rtk_fix, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/position_type_info")

            network_solution_plot.line(position_type_info_timestamp, position_type_info_is_network_solution, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/position_type_info")

            rtk_condition_plot.line(position_type_info_timestamp, rtk_cond_new_position_computed, alpha=0.8, color="blue", line_width=2, legend_label="rtk_cond_new_position_computed")
            rtk_condition_plot.line(position_type_info_timestamp, rtk_cond_no_synced_pair, alpha=0.8, color="red", line_width=2, legend_label="rtk_cond_no_synced_pair")
            rtk_condition_plot.line(position_type_info_timestamp, rtk_cond_insuff_dd_meas, alpha=0.8, color="green", line_width=2, legend_label="rtk_cond_insuff_dd_meas")
            rtk_condition_plot.line(position_type_info_timestamp, rtk_cond_ref_pos_unavailable, alpha=0.8, color="orange", line_width=2, legend_label="rtk_cond_ref_pos_unavailable")
            rtk_condition_plot.line(position_type_info_timestamp, rtk_cond_failed_integer_ver_with_fix_sol, alpha=0.8, color="purple", line_width=2, legend_label="rtk_cond_failed_integer_ver_with_fix_sol")
            rtk_condition_plot.line(position_type_info_timestamp, rtk_cond_sol_res_rms_exceeds, alpha=0.8, color="brown", line_width=2, legend_label="rtk_cond_sol_res_rms_exceeds")
            rtk_condition_plot.line(position_type_info_timestamp, rtk_cond_pdop_exceeds, alpha=0.8, color="cyan", line_width=2, legend_label="rtk_cond_pdop_exceeds")

            position_fix_type_plot.line(position_type_info_timestamp, position_type_info_position_fix_type, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/position_type_info")

        if (is_received_base_info_available):
            base_valid_plot.line(received_base_info_timestamp, base_valid, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/received_base_info")

            base_id_plot.line(received_base_info_timestamp, base_id, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/received_base_info")

            base_lat_plot.line(received_base_info_timestamp, base_lat_deg, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/received_base_info")

            base_lon_plot.line(received_base_info_timestamp, base_lon_deg, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/received_base_info")

            base_alt_plot.line(received_base_info_timestamp, base_height_m, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/received_base_info")

    
            base_north_plot.line(received_base_info_timestamp, base_north_m, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/received_base_info")

            base_east_plot.line(received_base_info_timestamp, base_east_m, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/received_base_info")

            base_down_plot.line(received_base_info_timestamp, base_down_m, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/received_base_info")


        num_of_svs_plot.legend.location = "top_right"
        num_of_svs_plot.legend.click_policy = "hide"  
        init_num_plot.legend.location = "top_right"
        init_num_plot.legend.click_policy = "hide"  
        position_flags_1_plot.legend.location = "top_right"
        position_flags_1_plot.legend.click_policy = "hide"  
        position_flags_2_plot.legend.location = "top_right"
        position_flags_2_plot.legend.click_policy = "hide"  

        correction_age_plot.legend.location = "top_right"
        correction_age_plot.legend.click_policy = "hide"  
        rtk_fix_plot.legend.location = "top_right"
        rtk_fix_plot.legend.click_policy = "hide"  
        network_solution_plot.legend.location = "top_right"
        network_solution_plot.legend.click_policy = "hide"  
        rtk_condition_plot.legend.location = "top_right"
        rtk_condition_plot.legend.click_policy = "hide"  
        position_fix_type_plot.legend.location = "top_right"
        position_fix_type_plot.legend.click_policy = "hide"  
        base_valid_plot.legend.location = "top_right"
        base_valid_plot.legend.click_policy = "hide"  
        base_id_plot.legend.location = "top_right"
        base_id_plot.legend.click_policy = "hide" 
        base_lat_plot.legend.location = "top_right"
        base_lat_plot.legend.click_policy = "hide" 
        base_lon_plot.legend.location = "top_right"
        base_lon_plot.legend.click_policy = "hide" 
        base_alt_plot.legend.location = "top_right"
        base_alt_plot.legend.click_policy = "hide" 
        base_north_plot.legend.location = "top_right"
        base_north_plot.legend.click_policy = "hide" 
        base_east_plot.legend.location = "top_right"
        base_east_plot.legend.click_policy = "hide" 
        base_down_plot.legend.location = "top_right"
        base_down_plot.legend.click_policy = "hide" 

        save(gridplot([
            [Div(text="<h1 style='margin-left: 2em;'>Position Time Info (1)</h1>")],
            [num_of_svs_plot, init_num_plot],
            [position_flags_1_plot, position_flags_2_plot],
            [Div(text="<h1 style='margin-left: 2em;'>Position Type Info (38)</h1>")],
            [network_solution_plot],
            [rtk_fix_plot, rtk_condition_plot],
            [correction_age_plot, position_fix_type_plot],
            [Div(text="<h1 style='margin-left: 2em;'>Reference Station (35)</h1>")],
            [base_valid_plot, base_id_plot, base_alt_plot],
            [base_lat_plot, base_lon_plot, base_alt_plot],
            [base_north_plot, base_east_plot, base_down_plot],
        ]))

    def generate_sky_plots(self, output_file_name):

        output_file(filename=output_file_name, title="IMU Plots")

        plts = []

        for i in range(len(self.CONSTELLATION_IDS)):

            p = figure(
                title=self.CONSTELLATION_NAMES[i] + " Sky Plot",
                width=self.SKY_PLOT_SIZE,
                height=self.SKY_PLOT_SIZE,
                x_range=(-self.HORIZON_RADIUS, self.HORIZON_RADIUS),
                y_range=(-self.HORIZON_RADIUS, self.HORIZON_RADIUS),
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
                radius=self.HORIZON_RADIUS,
                fill_alpha=0.0,
                line_color="black",
                line_dash="dashed",
            )

            # Draw elevation circles
            for elev_line in self.ELEVATION_CIRCLES:
                # Convert elevation to radial distance
                r = self.HORIZON_RADIUS - elev_line

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
            for az_line in self.AZIMUTH_LINES:
                # Convert az_line (in degrees) to radians
                theta = math.radians(az_line)

                # We'll draw from the center (0, 0) out to the edge at r=90 (the horizon)
                r = self.HORIZON_RADIUS

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

            for sat_name, sat_info in self.data_parser.gnss_data[self.CONSTELLATION_IDS[i]].items():

                az = np.array(sat_info["azimuth"])
                el = np.array(sat_info["elevation"])

                # ts = np.array(sat_info["timestamp"])
                ts = np.array(np.arange(len(az)))

                # Convert to radians
                theta = np.deg2rad(az)
                # Convert elevation -> r (0 at zenith, 90 at horizon)
                r = self.HORIZON_RADIUS - el

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

        p_num_sat = self.create_bokeh_figure("Number of Satellites", "Timestamp", "Num", height=self.TALL_PLOT_HEIGHT)
        num_sat = self.data_parser.gnss_data["num_of_sat"]
        num_sat_timestamp = np.arange(len(num_sat))
        p_num_sat.line(num_sat_timestamp, num_sat, alpha=0.8, color="blue", line_width=2, legend_label="")

        p_num_sat.legend.location = "top_right"
        p_num_sat.legend.click_policy = "hide"

        azimuth_plts = []
        elevation_plts = []
        snr_l1_plts = []
        snr_l2_plts = []
        snr_l5_plts = []
        used_in_position_plts = []
        used_in_rtk_plts = []

        for i in range(len(self.CONSTELLATION_IDS)):

            p_az = self.create_bokeh_figure(self.CONSTELLATION_NAMES[i] + " Azimuth", "Timestamp", "Azimuth [deg]", height=self.TALL_PLOT_HEIGHT)
            p_el = self.create_bokeh_figure(self.CONSTELLATION_NAMES[i] + " Elevation", "Timestamp", "Elevation [deg]", height=self.TALL_PLOT_HEIGHT)

            p_snr_l1 = self.create_bokeh_figure(self.CONSTELLATION_NAMES[i] + " SNR L1", "Timestamp", "[dB]", height=self.TALL_PLOT_HEIGHT)
            p_snr_l2 = self.create_bokeh_figure(self.CONSTELLATION_NAMES[i] + " SNR L2", "Timestamp", "[dB]", height=self.TALL_PLOT_HEIGHT)
            p_snr_l5 = self.create_bokeh_figure(self.CONSTELLATION_NAMES[i] + " SNR L5", "Timestamp", "[dB]", height=self.TALL_PLOT_HEIGHT)
            
            p_pos = self.create_bokeh_figure(self.CONSTELLATION_NAMES[i] + " Pos", "Timestamp", "Used [bool]", height=self.TALL_PLOT_HEIGHT)
            p_rtk = self.create_bokeh_figure(self.CONSTELLATION_NAMES[i] + " RTK", "Timestamp", "Used [bool]", height=self.TALL_PLOT_HEIGHT)

            line_colors = itertools.cycle(Category20[20])

            for (sat_name, sat_data), color in zip(self.data_parser.gnss_data[self.CONSTELLATION_IDS[i]].items(), line_colors):

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



        ne_plot = self.create_bokeh_figure("Northing-Easting", "Easting [m]", "Northing [m]", width=self.SQUARE_PLOT_SIZE, height=self.SQUARE_PLOT_SIZE)

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

