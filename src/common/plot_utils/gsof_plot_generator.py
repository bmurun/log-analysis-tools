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
        
        # Create hover tool for displaying (x, y) coordinates and source
        hover = HoverTool(tooltips=[
            ("Source", "@source"),
            ("X", "@x"),
            ("Y", "@y")
        ])
        
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
        p.add_tools(hover)
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
        

        # INS Solution plots with ColumnDataSource
        ins_ne_source = ColumnDataSource(data=dict(x=ins_solution_east, y=ins_solution_north, source=["/gsof/ins_solution"]*len(ins_solution_east)))
        ne_plot.scatter('x', 'y', source=ins_ne_source, alpha=0.8, size=2.5, color="red", legend_label="/gsof/ins_solution")

        ins_x_source = ColumnDataSource(data=dict(x=ins_solution_timestamp, y=ins_solution_north, source=["/gsof/ins_solution"]*len(ins_solution_timestamp)))
        x_plot.scatter('x', 'y', source=ins_x_source, alpha=0.8, size=2.5, color="red", legend_label="/gsof/ins_solution")

        ins_y_source = ColumnDataSource(data=dict(x=ins_solution_timestamp, y=ins_solution_east, source=["/gsof/ins_solution"]*len(ins_solution_timestamp)))
        y_plot.scatter('x', 'y', source=ins_y_source, alpha=0.8, size=2.5, color="red", legend_label="/gsof/ins_solution")

        ins_z_source = ColumnDataSource(data=dict(x=ins_solution_timestamp, y=ins_solution_down, source=["/gsof/ins_solution"]*len(ins_solution_timestamp)))
        z_plot.scatter('x', 'y', source=ins_z_source, alpha=0.8, size=2.5, color="red", legend_label="/gsof/ins_solution")

        ins_status_source = ColumnDataSource(data=dict(x=ins_solution_timestamp, y=ins_solution_gnss_status, source=["/gsof/ins_solution"]*len(ins_solution_timestamp)))
        gnss_status_plot.scatter('x', 'y', source=ins_status_source, alpha=0.8, size=2.5, color="red", legend_label="/gsof/ins_solution")

        if (is_code_pvt_available):
            # Code lat lon ht plots with ColumnDataSource
            code_ne_source = ColumnDataSource(data=dict(x=code_lat_long_east, y=code_lat_long_north, source=["/gsof/code_lat_lon_ht"]*len(code_lat_long_east)))
            ne_plot.scatter('x', 'y', source=code_ne_source, alpha=0.8, size=2.5, color="blue", legend_label="/gsof/code_lat_lon_ht")
        
            code_x_source = ColumnDataSource(data=dict(x=code_lat_long_ht_timestamp, y=code_lat_long_north, source=["/gsof/code_lat_lon_ht"]*len(code_lat_long_ht_timestamp)))
            x_plot.scatter('x', 'y', source=code_x_source, alpha=0.8, size=2.5, color="blue", legend_label="/gsof/code_lat_lon_ht")

            code_y_source = ColumnDataSource(data=dict(x=code_lat_long_ht_timestamp, y=code_lat_long_east, source=["/gsof/code_lat_lon_ht"]*len(code_lat_long_ht_timestamp)))
            y_plot.scatter('x', 'y', source=code_y_source, alpha=0.8, size=2.5, color="blue", legend_label="/gsof/code_lat_lon_ht")

            code_z_source = ColumnDataSource(data=dict(x=code_lat_long_ht_timestamp, y=code_lat_long_down, source=["/gsof/code_lat_lon_ht"]*len(code_lat_long_ht_timestamp)))
            z_plot.scatter('x', 'y', source=code_z_source, alpha=0.8, size=2.5, color="blue", legend_label="/gsof/code_lat_lon_ht")

            code_status_source = ColumnDataSource(data=dict(x=code_lat_long_ht_timestamp, y=code_lat_lon_position_type, source=["/gsof/code_lat_lon_ht"]*len(code_lat_long_ht_timestamp)))
            gnss_status_plot.scatter('x', 'y', source=code_status_source, alpha=0.8, size=2.5, color="blue", legend_label="/gsof/code_lat_lon_ht")


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

        # Create a second altitude plot to avoid repeated layout child error
        base_alt_plot2 = self.create_bokeh_figure("Reference Station Altitude WGS84", "Timestamp", "Altitude [m]")

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
            # Position time info plots with ColumnDataSource
            num_svs_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=position_time_info_num_svs, source=["/gsof/position_time_info"]*len(position_time_info_timestamp)))
            num_of_svs_plot.line('x', 'y', source=num_svs_source, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/position_time_info")
            
            init_num_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=position_time_info_init_num, source=["/gsof/position_time_info"]*len(position_time_info_timestamp)))
            init_num_plot.line('x', 'y', source=init_num_source, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/position_time_info")

            # Position flags 1 plots with ColumnDataSource
            new_pos_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=new_position, source=["new_position"]*len(position_time_info_timestamp)))
            position_flags_1_plot.line('x', 'y', source=new_pos_source, alpha=0.8, color="blue", line_width=2, legend_label="new_position")
            
            clock_fix_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=clock_fix, source=["clock_fix"]*len(position_time_info_timestamp)))
            position_flags_1_plot.line('x', 'y', source=clock_fix_source, alpha=0.8, color="red", line_width=2, legend_label="clock_fix")
            
            horiz_coords_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=horizontal_coords, source=["horizontal_coords"]*len(position_time_info_timestamp)))
            position_flags_1_plot.line('x', 'y', source=horiz_coords_source, alpha=0.8, color="green", line_width=2, legend_label="horizontal_coords")
            
            height_calc_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=height_calculated, source=["height_calculated"]*len(position_time_info_timestamp)))
            position_flags_1_plot.line('x', 'y', source=height_calc_source, alpha=0.8, color="orange", line_width=2, legend_label="height_calculated")
            
            least_sq_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=least_squares, source=["least_squares"]*len(position_time_info_timestamp)))
            position_flags_1_plot.line('x', 'y', source=least_sq_source, alpha=0.8, color="purple", line_width=2, legend_label="least_squares")
            
            filtered_l1_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=filtered_L1, source=["filtered_L1"]*len(position_time_info_timestamp)))
            position_flags_1_plot.line('x', 'y', source=filtered_l1_source, alpha=0.8, color="cyan", line_width=2, legend_label="filtered_L1")

            # Position flags 2 plots with ColumnDataSource
            is_diff_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=is_differential, source=["is_differential"]*len(position_time_info_timestamp)))
            position_flags_2_plot.line('x', 'y', source=is_diff_source, alpha=0.8, color="blue", line_width=2, legend_label="is_differential")
            
            uses_phase_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=uses_phase, source=["uses_phase"]*len(position_time_info_timestamp)))
            position_flags_2_plot.line('x', 'y', source=uses_phase_source, alpha=0.8, color="red", line_width=2, legend_label="uses_phase")
            
            is_rtk_fixed_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=is_RTK_fixed, source=["is_RTK_fixed"]*len(position_time_info_timestamp)))
            position_flags_2_plot.line('x', 'y', source=is_rtk_fixed_source, alpha=0.8, color="green", line_width=2, legend_label="is_RTK_fixed")
            
            omnistar_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=omniSTAR, source=["omniSTAR"]*len(position_time_info_timestamp)))
            position_flags_2_plot.line('x', 'y', source=omnistar_source, alpha=0.8, color="orange", line_width=2, legend_label="omniSTAR")
            
            static_const_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=static_constraint, source=["static_constraint"]*len(position_time_info_timestamp)))
            position_flags_2_plot.line('x', 'y', source=static_const_source, alpha=0.8, color="purple", line_width=2, legend_label="static_constraint")
            
            network_rtk_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=network_RTK, source=["network_RTK"]*len(position_time_info_timestamp)))
            position_flags_2_plot.line('x', 'y', source=network_rtk_source, alpha=0.8, color="cyan", line_width=2, legend_label="network_RTK")
            
            dithered_rtk_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=dithered_RTK, source=["dithered_RTK"]*len(position_time_info_timestamp)))
            position_flags_2_plot.line('x', 'y', source=dithered_rtk_source, alpha=0.8, color="brown", line_width=2, legend_label="dithered_RTK")
            
            beacon_dgnss_source = ColumnDataSource(data=dict(x=position_time_info_timestamp, y=beacon_DGNSS, source=["beacon_DGNSS"]*len(position_time_info_timestamp)))
            position_flags_2_plot.line('x', 'y', source=beacon_dgnss_source, alpha=0.8, color="magenta", line_width=2, legend_label="beacon_DGNSS")


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

            # Add the same data to the second altitude plot
            base_alt_plot2.line(received_base_info_timestamp, base_height_m, alpha=0.8, color="blue", line_width=2, legend_label="/gsof/received_base_info")

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

        # Set legend for the second altitude plot
        base_alt_plot2.legend.location = "top_right"
        base_alt_plot2.legend.click_policy = "hide"

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
            [base_lat_plot, base_lon_plot, base_alt_plot2],
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
                p.scatter(x="x", y="y", size="size", alpha=0.8, source=source, color="blue")

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

        # DOP plots with ColumnDataSource
        hdop_source = ColumnDataSource(data=dict(x=timestamp, y=hdop, source=["pdop_info"]*len(timestamp)))
        p_hdop.line('x', 'y', source=hdop_source, alpha=0.8, color="blue", line_width=2, legend_label="pdop_info")
        
        vdop_source = ColumnDataSource(data=dict(x=timestamp, y=vdop, source=["pdop_info"]*len(timestamp)))
        p_vdop.line('x', 'y', source=vdop_source, alpha=0.8, color="blue", line_width=2, legend_label="pdop_info")

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

        # Noise plots with ColumnDataSource
        ins_north_source = ColumnDataSource(data=dict(x=self.data_parser.ins_solution_rms["timestamp"], y=self.data_parser.ins_solution_rms["north_stdev"], source=["INSSolutionRms"]*len(self.data_parser.ins_solution_rms["timestamp"])))
        p_north_stdev.line('x', 'y', source=ins_north_source, alpha=0.8, color="green", line_width=2, legend_label="INSSolutionRms")
        
        ins_east_source = ColumnDataSource(data=dict(x=self.data_parser.ins_solution_rms["timestamp"], y=self.data_parser.ins_solution_rms["east_stdev"], source=["INSSolutionRms"]*len(self.data_parser.ins_solution_rms["timestamp"])))
        p_east_stdev.line('x', 'y', source=ins_east_source, alpha=0.8, color="green", line_width=2, legend_label="INSSolutionRms")
        
        ins_down_source = ColumnDataSource(data=dict(x=self.data_parser.ins_solution_rms["timestamp"], y=self.data_parser.ins_solution_rms["down_stdev"], source=["INSSolutionRms"]*len(self.data_parser.ins_solution_rms["timestamp"])))
        p_down_stdev.line('x', 'y', source=ins_down_source, alpha=0.8, color="green", line_width=2, legend_label="INSSolutionRms")

        len_pos_sigma = len(self.data_parser.position_sigma["north_stdev"])

        self.data_parser.position_sigma["timestamp"] = np.linspace(first_ts, last_ts, len_pos_sigma)

        # Position sigma plots with ColumnDataSource
        pos_sigma_north_source = ColumnDataSource(data=dict(x=self.data_parser.position_sigma["timestamp"], y=self.data_parser.position_sigma["north_stdev"], source=["PositionSigma12"]*len(self.data_parser.position_sigma["timestamp"])))
        p_north_stdev.line('x', 'y', source=pos_sigma_north_source, alpha=0.8, color="blue", line_width=2, legend_label="PositionSigma12")
        
        pos_sigma_east_source = ColumnDataSource(data=dict(x=self.data_parser.position_sigma["timestamp"], y=self.data_parser.position_sigma["east_stdev"], source=["PositionSigma12"]*len(self.data_parser.position_sigma["timestamp"])))
        p_east_stdev.line('x', 'y', source=pos_sigma_east_source, alpha=0.8, color="blue", line_width=2, legend_label="PositionSigma12")
        
        pos_sigma_down_source = ColumnDataSource(data=dict(x=self.data_parser.position_sigma["timestamp"], y=self.data_parser.position_sigma["down_stdev"], source=["PositionSigma12"]*len(self.data_parser.position_sigma["timestamp"])))
        p_down_stdev.line('x', 'y', source=pos_sigma_down_source, alpha=0.8, color="blue", line_width=2, legend_label="PositionSigma12")

        len_navsat = len(self.data_parser.navsat["north_stdev"])

        self.data_parser.navsat["timestamp"] = np.linspace(first_ts, last_ts, len_navsat)

        # NavSat plots with ColumnDataSource
        navsat_north_source = ColumnDataSource(data=dict(x=self.data_parser.navsat["timestamp"], y=self.data_parser.navsat["north_stdev"], source=["NavSat"]*len(self.data_parser.navsat["timestamp"])))
        p_north_stdev.line('x', 'y', source=navsat_north_source, alpha=0.8, color="red", line_width=2, legend_label="NavSat")
        
        navsat_east_source = ColumnDataSource(data=dict(x=self.data_parser.navsat["timestamp"], y=self.data_parser.navsat["east_stdev"], source=["NavSat"]*len(self.data_parser.navsat["timestamp"])))
        p_east_stdev.line('x', 'y', source=navsat_east_source, alpha=0.8, color="red", line_width=2, legend_label="NavSat")
        
        navsat_down_source = ColumnDataSource(data=dict(x=self.data_parser.navsat["timestamp"], y=self.data_parser.navsat["down_stdev"], source=["NavSat"]*len(self.data_parser.navsat["timestamp"])))
        p_down_stdev.line('x', 'y', source=navsat_down_source, alpha=0.8, color="red", line_width=2, legend_label="NavSat")

        p_north_stdev.legend.location = "top_right"
        p_north_stdev.legend.click_policy = "hide"
        p_east_stdev.legend.location = "top_right"
        p_east_stdev.legend.click_policy = "hide"
        p_down_stdev.legend.location = "top_right"
        p_down_stdev.legend.click_policy = "hide"

        p_gnss_status = self.create_bokeh_figure("GNSS Status", "Timestamp", "Status")

        gnss_status_rms_source = ColumnDataSource(data=dict(x=self.data_parser.ins_solution_rms["timestamp"], y=self.data_parser.ins_solution_rms["gnss_status"], source=["INSSolutionRms"]*len(self.data_parser.ins_solution_rms["timestamp"])))
        p_gnss_status.line('x', 'y', source=gnss_status_rms_source, alpha=0.8, color="red", line_width=2, legend_label="INSSolutionRms")

        p_gnss_status_ins_solution = self.create_bokeh_figure("GNSS Status", "Timestamp", "Status")

        gnss_status_ins_source = ColumnDataSource(data=dict(x=self.data_parser.ins_solution["timestamp"], y=self.data_parser.ins_solution["gnss_status"], source=["INSSolution"]*len(self.data_parser.ins_solution["timestamp"])))
        p_gnss_status_ins_solution.line('x', 'y', source=gnss_status_ins_source, alpha=0.8, color="red", line_width=2, legend_label="INSSolution")

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
    

