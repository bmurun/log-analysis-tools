import numpy as np
import os
import sys
import click
import bisect

from mcap.reader import make_reader
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message

from common.data_parser.topics import MCAP_ROS2_GSOF_TOPICS

class GsofDataParser:
    def __init__(self, log_path: str) -> None:
        if not os.path.isfile(log_path):
            click.echo(click.style("Input file doesn't exist!", fg="red"))
            return
        self.mcap_data = {}

        self.position_sigma = {
            "timestamp" : [],
            "north_stdev": [],
            "east_stdev": [],
            "down_stdev": []
        }

        self.gnss_data = {
            "timestamp": [],
            "num_of_sat": [],
            "0" : {}, # GPS
            "2": {},  # GLONASS
            "3": {},  # Galileo
            "5": {},  # BeiDou
        }

        self.ins_solution = {
            "timestamp": [],
            "gnss_status": [],
            "imu_alignment": [],
            "latitude": [],
            "longitude": [],
            "altitude": [],
            "roll": [],
            "pitch": [],
            "yaw": [],
            "north_vel": [],
            "east_vel": [],
            "down_vel": [],
            "gps_time_ms": [],
            "gps_week_number": [],
        }

        self.ins_solution_rms = {
            "timestamp": [],
            "north_stdev": [],
            "east_stdev": [],
            "down_stdev": [],
            "gnss_status": [],
        }

        self.navsat = {
            "timestamp": [],
            "latitude": [],
            "longitude": [],
            "altitude": [],
            "north_stdev": [],
            "east_stdev": [],
            "down_stdev": []
        }


        self.pdop_info = {
            "timestamp": [],
            "hdop": [],
            "pdop": [],
            "vdop": [],
            "tdop": []
        }

        self.received_base_info = {
            "timestamp": [],
            "base_valid": [],
            "base_id": [],
            "base_lat_deg": [],
            "base_lon_deg": [],
            "base_height_m": [],
        }

        self.position_time_info = {
            "timestamp": [],
            "gps_time_ms": [],
            "gps_week_number": [],
            "number_space_vehicles_used": [],
            "init_num": [], # Increments with each initialization (modulo 256)
            "new_position": [], # position flags 1
            "clock_fix": [],
            "horizontal_coords": [],
            "height_calculated": [],
            "least_squares": [],
            "filtered_L1": [],
            "is_differential": [],
            "uses_phase": [],
            "is_RTK_fixed": [],
            "omniSTAR": [],
            "static_constraint": [],
            "network_RTK": [],
            "dithered_RTK": [],
            "beacon_DGNSS": [],
            "network_flags_new_base_station_available": [],
            "network_flags_rtcm_not_available": [],
            "network_flags_rtcm_not_full_cycle": [],
            "network_flags_rtcm_insufficient": [],
            "network_flags_rtcm_good": [],
            "network_flags_outside_of_geofence": [],
            "network_flags_outside_of_rtk_range_base": [],
            "network_flags_xfill_operation": [],
            "network_flags_rtx_position_flag": [],
            "network_flags_rtx_xfill_is_down": [],
            "network_flags2_xfill_ready": [],
            "network_flags2_rtx_fast": [],
            "network_flags2_rtx_offset": [],
        }

        self.position_type_info = {
            "timestamp": [],
            "error_scale": [],
            "is_network_solution": [],
            "is_rtk_fix": [],
            "init_integrity_1": [],
            "init_integrity_2": [],
            "rtk_cond_new_position_computed": [],
            "rtk_cond_no_synced_pair": [],
            "rtk_cond_insuff_dd_meas": [],
            "rtk_cond_ref_pos_unavailable": [],
            "rtk_cond_failed_integer_ver_with_fix_sol": [],
            "rtk_cond_sol_res_rms_exceeds": [],
            "rtk_cond_pdop_exceeds": [],
            "correction_age": [],
            "position_fix_type": [],
            "network_flags_new_base_station_available": [],
            "network_flags_rtcm_not_available": [],
            "network_flags_rtcm_not_full_cycle": [],
            "network_flags_rtcm_insufficient": [],
            "network_flags_rtcm_good": [],
            "network_flags_outside_of_geofence": [],
            "network_flags_outside_of_rtk_range_base": [],
            "network_flags_xfill_operation": [],
            "network_flags_rtx_position_flag": [],
            "network_flags_rtx_xfill_is_down": [],
            "network_flags2_xfill_ready": [],
            "network_flags2_rtx_fast": [],
            "network_flags2_rtx_offset": [],
            "network_flags2_cmrxe_received": [],
            "network_flags2_rtx_wet_area": [],
        }

        self.code_lat_lon_ht = {
            "timestamp": [],
            "position_type": [],
            "lat_deg": [],
            "lon_deg": [],
            "height": [],
            "gps_week_number": [],
            "gps_time_ms": []
        }

        self._parse_mcap(log_path)
        self._process_data()

    def _parse_mcap(self, log_path) -> None:

        with open(log_path, "rb") as f:
            reader = make_reader(f)

            for schema, channel, message in reader.iter_messages(topics=MCAP_ROS2_GSOF_TOPICS):

                topic_name = channel.topic
                msg_str = schema.name
                msg_type = get_message(msg_str)
                data = deserialize_message(message.data, msg_type)

                if topic_name == "/lvx_client/gsof/all_sv_detailed_info_34":
                    self.gnss_data["num_of_sat"].append(len(data.sv_info))

                    sec= data.header.stamp.sec
                    nanosec = data.header.stamp.nanosec
                    timestamp = sec + nanosec / 1e9

                    self.gnss_data["timestamp"].append(timestamp)

                    for i in range(len(data.sv_info)):
                        sat_data = data.sv_info[i]

                        prn = sat_data.prn
                        sv_system = sat_data.sv_system

                        if (sv_system == 1): # SBAS
                            continue

                        azimuth = sat_data.azimuth
                        elevation = sat_data.elevation

                        snr_l1 = float(sat_data.snr_l1) / 4
                        snr_l2 = float(sat_data.snr_l2) / 4
                        snr_l5 = float(sat_data.snr_l5) / 4

                        sv_flags1 = sat_data.sv_flags1

                        is_used_in_position = (sv_flags1 & (1 << 6)) >> 6
                        is_used_in_rtk = (sv_flags1 & (1 << 7)) >> 7

                        if str(prn) not in self.gnss_data[str(sv_system)].keys():
                            self.gnss_data[str(sv_system)][str(prn)] = {
                                "azimuth": [],
                                "elevation": [],
                                "snr_l1": [],
                                "snr_l2": [],
                                "snr_l5": [],
                                "is_used_in_position": [],
                                "is_used_in_rtk": [],
                            }

                        self.gnss_data[str(sv_system)][str(prn)]["azimuth"].append(azimuth)
                        self.gnss_data[str(sv_system)][str(prn)]["elevation"].append(elevation)
                        self.gnss_data[str(sv_system)][str(prn)]["snr_l1"].append(snr_l1)
                        self.gnss_data[str(sv_system)][str(prn)]["snr_l2"].append(snr_l2)
                        self.gnss_data[str(sv_system)][str(prn)]["snr_l5"].append(snr_l5)
                        self.gnss_data[str(sv_system)][str(prn)]["is_used_in_position"].append(is_used_in_position)
                        self.gnss_data[str(sv_system)][str(prn)]["is_used_in_rtk"].append(is_used_in_rtk)
                elif topic_name == "/lvx_client/gsof/position_sigma_info_12":

                    north_stdev = data.sigma_north 
                    east_stdev = data.sigma_east 
                    down_stdev = data.sigma_up

                    sec= data.header.stamp.sec
                    nanosec = data.header.stamp.nanosec
                    timestamp = sec + nanosec / 1e9

                    self.position_sigma["timestamp"].append(timestamp)

                    self.position_sigma["north_stdev"].append(north_stdev)
                    self.position_sigma["east_stdev"].append(east_stdev)
                    self.position_sigma["down_stdev"].append(down_stdev)
                elif topic_name == "/lvx_client/gsof/lband_status_info_40":
                    pass
                elif topic_name == "/lvx_client/gsof/ins_solution_49":

                    sec = data.header.stamp.sec
                    nanosec = data.header.stamp.nanosec
                    timestamp = sec + nanosec / 1e9

                    self.ins_solution["timestamp"].append(timestamp)
                    self.ins_solution["gnss_status"].append(data.status.gnss)
                    self.ins_solution["imu_alignment"].append(data.status.imu_alignment)
                    self.ins_solution["latitude"].append(data.lla.latitude)
                    self.ins_solution["longitude"].append(data.lla.longitude)
                    self.ins_solution["altitude"].append(data.lla.altitude)
                    self.ins_solution["roll"].append(data.roll)
                    self.ins_solution["pitch"].append(data.pitch)
                    self.ins_solution["yaw"].append(data.heading)
                    self.ins_solution["north_vel"].append(data.velocity.north)
                    self.ins_solution["east_vel"].append(data.velocity.east)
                    self.ins_solution["down_vel"].append(data.velocity.down)

                    self.ins_solution["gps_time_ms"].append(data.gps_time.time)
                    self.ins_solution["gps_week_number"].append(data.gps_time.week)

                elif topic_name == "/lvx_client/navsat":

                    sec = data.header.stamp.sec
                    nanosec = data.header.stamp.nanosec
                    timestamp = sec + nanosec / 1e9

                    self.navsat["timestamp"].append(timestamp)
                    self.navsat["latitude"].append(data.latitude)
                    self.navsat["longitude"].append(data.longitude)
                    self.navsat["altitude"].append(data.altitude)
                    self.navsat["north_stdev"].append(np.sqrt(data.position_covariance[0]))
                    self.navsat["east_stdev"].append(np.sqrt(data.position_covariance[4]))
                    self.navsat["down_stdev"].append(np.sqrt(data.position_covariance[8]))
                elif topic_name == "/lvx_client/gsof/pdop_info_9":
                    sec = data.header.stamp.sec
                    nanosec = data.header.stamp.nanosec
                    timestamp = sec + nanosec / 1e9
                    self.pdop_info["timestamp"].append(timestamp)
                    self.pdop_info["hdop"].append(data.horiziontal_dop)
                    self.pdop_info["pdop"].append(data.position_dop)
                    self.pdop_info["vdop"].append(data.vertical_dop)
                    self.pdop_info["tdop"].append(data.time_dop)
                elif topic_name == "/lvx_client/gsof/ins_solution_rms_50":
                    sec = data.header.stamp.sec
                    nanosec = data.header.stamp.nanosec
                    timestamp = sec + nanosec / 1e9
                    self.ins_solution_rms["timestamp"].append(timestamp)
                    self.ins_solution_rms["north_stdev"].append(data.pos_rms_error.north)
                    self.ins_solution_rms["east_stdev"].append(data.pos_rms_error.east)
                    self.ins_solution_rms["down_stdev"].append(data.pos_rms_error.down)
                    self.ins_solution_rms["gnss_status"].append(data.status.gnss)

                elif topic_name == "/lvx_client/gsof/position_time_info_1":
                    # https://receiverhelp.trimble.com/alloy-gnss/en-us/gsof-messages-flags.html#Position%20flags%201
                    # sec = data.header.stamp.sec
                    # nanosec = data.header.stamp.nanosec
                    # timestamp = sec + nanosec / 1e9
                    timestamp = 0
                    self.position_time_info["timestamp"].append(timestamp)
                    self.position_time_info["gps_time_ms"].append(data.gps_time.time)
                    self.position_time_info["gps_week_number"].append(data.gps_time.week)
                    self.position_time_info["number_space_vehicles_used"].append(data.number_space_vehicles_used)
                    self.position_time_info["init_num"].append(data.init_num)

                    new_position = bool(data.position_flags_1 & (1 << 0))  # Bit 0
                    clock_fix = bool(data.position_flags_1 & (1 << 1))     # Bit 1
                    horizontal_coords = bool(data.position_flags_1 & (1 << 2))  # Bit 2
                    height_calculated = bool(data.position_flags_1 & (1 << 3))  # Bit 3
                    least_squares = bool(data.position_flags_1 & (1 << 5))      # Bit 5
                    filtered_L1 = bool(data.position_flags_1 & (1 << 7))        # Bit 7

                    self.position_time_info["new_position"].append(new_position)
                    self.position_time_info["clock_fix"].append(clock_fix)
                    self.position_time_info["horizontal_coords"].append(horizontal_coords)
                    self.position_time_info["height_calculated"].append(height_calculated)
                    self.position_time_info["least_squares"].append(least_squares)
                    self.position_time_info["filtered_L1"].append(filtered_L1)

                    is_differential = bool(data.position_flags_2 & (1 << 0))       # Bit 0
                    uses_phase = bool(data.position_flags_2 & (1 << 1))            # Bit 1
                    is_RTK_fixed = bool(data.position_flags_2 & (1 << 2))          # Bit 2
                    omniSTAR = bool(data.position_flags_2 & (1 << 3))              # Bit 3
                    static_constraint = bool(data.position_flags_2 & (1 << 4))     # Bit 4
                    network_RTK = bool(data.position_flags_2 & (1 << 5))           # Bit 5
                    dithered_RTK = bool(data.position_flags_2 & (1 << 6))          # Bit 6
                    beacon_DGNSS = bool(data.position_flags_2 & (1 << 7))          # Bit 7

                    self.position_time_info["is_differential"].append(is_differential)
                    self.position_time_info["uses_phase"].append(uses_phase)
                    self.position_time_info["is_RTK_fixed"].append(is_RTK_fixed)
                    self.position_time_info["omniSTAR"].append(omniSTAR)
                    self.position_time_info["static_constraint"].append(static_constraint)
                    self.position_time_info["network_RTK"].append(network_RTK)
                    self.position_time_info["dithered_RTK"].append(dithered_RTK)
                    self.position_time_info["beacon_DGNSS"].append(beacon_DGNSS)


                elif topic_name == "/lvx_client/gsof/received_base_info_35":
                    sec = data.header.stamp.sec
                    nanosec = data.header.stamp.nanosec
                    timestamp = sec + nanosec / 1e9

                    self.received_base_info["timestamp"].append(timestamp)
                    self.received_base_info["base_id"].append(data.id)
                    self.received_base_info["base_lat_deg"].append(np.rad2deg(data.llh.latitude))
                    self.received_base_info["base_lon_deg"].append(np.rad2deg(data.llh.longitude))
                    self.received_base_info["base_height_m"].append(data.llh.height)

                    base_valid = bool(data.flags & (1 << 3))
                    self.received_base_info["base_valid"].append(base_valid)
                elif topic_name == "/lvx_client/gsof/position_type_info_38":
                    # https://receiverhelp.trimble.com/alloy-gnss/en-us/gsof-messages-position-type.html?tocpath=Output%20Messages%7CGSOF%20messages%7C_____24
                    # sec = data.header.stamp.sec
                    # nanosec = data.header.stamp.nanosec
                    # timestamp = sec + nanosec / 1e9
                    # Header timestamp is not populated
                    timestamp = 0

                    self.position_type_info["timestamp"].append(timestamp)
                    self.position_type_info["error_scale"].append(data.error_scale)
                    self.position_type_info["is_network_solution"].append(bool(data.solution_flags & (1 << 0)))
                    self.position_type_info["is_rtk_fix"].append(bool(data.solution_flags & (1 << 1)))
                    self.position_type_info["init_integrity_1"].append(bool(data.solution_flags & (1 << 2)))
                    self.position_type_info["init_integrity_2"].append(bool(data.solution_flags & (1 << 3)))

                    # RTK condition is a value
                    self.position_type_info["rtk_cond_new_position_computed"].append(data.rtk_condition == 0)
                    self.position_type_info["rtk_cond_no_synced_pair"].append(data.rtk_condition == 1)
                    self.position_type_info["rtk_cond_insuff_dd_meas"].append(data.rtk_condition == 2)
                    self.position_type_info["rtk_cond_ref_pos_unavailable"].append(data.rtk_condition == 3)
                    self.position_type_info["rtk_cond_failed_integer_ver_with_fix_sol"].append(data.rtk_condition == 4)
                    self.position_type_info["rtk_cond_sol_res_rms_exceeds"].append(data.rtk_condition == 5)
                    self.position_type_info["rtk_cond_pdop_exceeds"].append(data.rtk_condition == 6)

                    self.position_type_info["correction_age"].append(data.correction_age)
                    self.position_type_info["position_fix_type"].append(data.position_fix_type)

                    # Network Flags (Field 12) - Parse as bitmap
                    network_flags = data.network_flags
                    
                    # Bit 0: New physical base station available
                    self.position_type_info["network_flags_new_base_station_available"].append(bool(network_flags & (1 << 0)))
                    
                    # Bits 2,1: RTCM v3 Network messages status (4 combinations)
                    rtcm_bits = (network_flags >> 1) & 0x03  # Extract bits 2,1
                    self.position_type_info["network_flags_rtcm_not_available"].append(rtcm_bits == 0)      # 0,0
                    self.position_type_info["network_flags_rtcm_not_full_cycle"].append(rtcm_bits == 1)    # 0,1  
                    self.position_type_info["network_flags_rtcm_insufficient"].append(rtcm_bits == 2)      # 1,0
                    self.position_type_info["network_flags_rtcm_good"].append(rtcm_bits == 3)              # 1,1
                    
                    # Bit 3: GeoFence option is enabled and unit is outside Geofence area
                    self.position_type_info["network_flags_outside_of_geofence"].append(bool(network_flags & (1 << 3)))
                    
                    # Bit 4: RTK Range limiting is enabled and unit is too far from the base
                    self.position_type_info["network_flags_outside_of_rtk_range_base"].append(bool(network_flags & (1 << 4)))
                    
                    # Bit 5: xFill operation
                    self.position_type_info["network_flags_xfill_operation"].append(bool(network_flags & (1 << 5)))
                    
                    # Bit 6: RTX position flag. 1 = RTX position, 0 = Not RTX position
                    self.position_type_info["network_flags_rtx_position_flag"].append(bool(network_flags & (1 << 6)))
                    
                    # Bit 7: RTX/xFill link is down. 1 = link is down, 0 = not applicable
                    self.position_type_info["network_flags_rtx_xfill_is_down"].append(bool(network_flags & (1 << 7)))

                    # Network Flags2 (Field 13) - Parse as bitmap
                    network_flags2 = data.network_flags2
                    
                    # Bit 0: xFill is ready to propagate RTK positions (or is already running)
                    self.position_type_info["network_flags2_xfill_ready"].append(bool(network_flags2 & (1 << 0)))
                    
                    # Bit 1: RTX solution is RTX Fast
                    self.position_type_info["network_flags2_rtx_fast"].append(bool(network_flags2 & (1 << 1)))
                    
                    # Bit 2: xFill-RTX offset (from RTK) is known to an acceptable accuracy to propagate RTK
                    self.position_type_info["network_flags2_rtx_offset"].append(bool(network_flags2 & (1 << 2)))
                    
                    # Bit 3: If set to 1, indicates that CMRxe is being received
                    self.position_type_info["network_flags2_cmrxe_received"].append(bool(network_flags2 & (1 << 3)))
                    
                    # Bit 4: If set, indicates RTX is in a "wet" area
                    self.position_type_info["network_flags2_rtx_wet_area"].append(bool(network_flags2 & (1 << 4)))

                elif topic_name == "/lvx_client/gsof/code_lat_long_ht_62":
                    sec = data.header.stamp.sec
                    nanosec = data.header.stamp.nanosec
                    timestamp = sec + nanosec / 1e9

                    self.code_lat_lon_ht["timestamp"].append(timestamp)
                    self.code_lat_lon_ht["position_type"].append(data.position_type)
                    self.code_lat_lon_ht["gps_week_number"].append(data.gps_time.week)
                    self.code_lat_lon_ht["gps_time_ms"].append(data.gps_time.time)
                    self.code_lat_lon_ht["lat_deg"].append(np.rad2deg(data.latitude))
                    self.code_lat_lon_ht["lon_deg"].append(np.rad2deg(data.longitude))
                    self.code_lat_lon_ht["height"].append(data.height)
                else:
                    pass
                

    def _process_data(self):
        """ Post processing: Interpolate or validate timestamps and normalize to start from zero
        """

        # First interpolate PositionTypeInfo timestamp given INS solution timestamp
        def compute_gps_total_ms(week, ms):
            return week * 7 * 24 * 3600 * 1000 + ms
        
        ins_gps_total_ms = [
            compute_gps_total_ms(week, ms)
            for week, ms in zip(self.ins_solution["gps_week_number"], self.ins_solution["gps_time_ms"])
        ]

        ins_timestamps = self.ins_solution["timestamp"]

        sorted_indices = sorted(range(len(ins_gps_total_ms)), key=lambda i: ins_gps_total_ms[i])
        ins_gps_total_ms_sorted = [ins_gps_total_ms[i] for i in sorted_indices]
        ins_timestamps_sorted = [ins_timestamps[i] for i in sorted_indices]

        def interpolate_timestamp(target_ms):
            idx = bisect.bisect_left(ins_gps_total_ms_sorted, target_ms)
            if idx == 0:
                return ins_timestamps_sorted[0]
            elif idx >= len(ins_gps_total_ms_sorted):
                return ins_timestamps_sorted[-1]
            else:
                x0, x1 = ins_gps_total_ms_sorted[idx - 1], ins_gps_total_ms_sorted[idx]
                y0, y1 = ins_timestamps_sorted[idx - 1], ins_timestamps_sorted[idx]
                return y0 + (y1 - y0) * (target_ms - x0) / (x1 - x0)

        # Step 4: Apply interpolation for all entries in self.position_time_info
        self.position_time_info["timestamp"] = []
        for week, ms in zip(self.position_time_info["gps_week_number"], self.position_time_info["gps_time_ms"]):
            target_ms = compute_gps_total_ms(week, ms)
            interp_ts = interpolate_timestamp(target_ms)
            self.position_time_info["timestamp"].append(interp_ts)

        # Now use PositionTimeInfo to extract PositionTypeInfo timestamp
        # use the uses_phase from PositionTimeInfo and position_fix_type from PositionTypeInfo

        # First generate base timestamp assuming uniform sampling over the same time period for PositionTypeInfo
        num_pts = len(self.position_type_info["position_fix_type"])
        start_time = self.position_time_info["timestamp"][0]
        end_time = self.position_time_info["timestamp"][-1]
        uniform_timestamps = np.linspace(start_time, end_time, num_pts)

        # Find the transition times (0→1) from PositionTimeInfo's uses_phase
        position_time_info_uses_phase = self.position_time_info["uses_phase"]
        position_time_info_timestamp = np.array(self.position_time_info["timestamp"])
        transition_times = position_time_info_timestamp[np.where(np.diff(np.array(position_time_info_uses_phase, dtype=int)) == 1)[0] + 1]

        # Find where position_fix_type transitions TO 29 (not just where it IS 29)
        position_fix_type_arr = np.array(self.position_type_info["position_fix_type"])
        transitions_to_29_indices = np.where((position_fix_type_arr[1:] == 29) & (position_fix_type_arr[:-1] != 29))[0] + 1

        # Create a mapping of transition indices to transition times
        transition_mapping = {}
        min_transitions = min(len(transitions_to_29_indices), len(transition_times))
        for i in range(min_transitions):
            transition_mapping[transitions_to_29_indices[i]] = transition_times[i]

        # First pass: create sparse array with known good timestamps at transition points
        position_type_info_timestamp = [None] * len(self.position_type_info["position_fix_type"])
        for i in transition_mapping:
            position_type_info_timestamp[i] = transition_mapping[i]
        
        # Second pass: interpolate between known good points
        known_indices = sorted(transition_mapping.keys())
        
        if len(known_indices) > 0:
            # Handle start to first known point
            start_time = self.position_time_info["timestamp"][0]
            for i in range(known_indices[0]):
                if known_indices[0] > 0:
                    position_type_info_timestamp[i] = start_time + (transition_mapping[known_indices[0]] - start_time) * i / known_indices[0]
                else:
                    position_type_info_timestamp[i] = transition_mapping[known_indices[0]]
            
            # Handle between known points
            for j in range(len(known_indices) - 1):
                start_idx = known_indices[j]
                end_idx = known_indices[j + 1]
                start_time = transition_mapping[start_idx]
                end_time = transition_mapping[end_idx]
                
                for i in range(start_idx + 1, end_idx):
                    progress = (i - start_idx) / (end_idx - start_idx)
                    position_type_info_timestamp[i] = start_time + (end_time - start_time) * progress
            
            # Handle after last known point
            last_idx = known_indices[-1]
            end_time = self.position_time_info["timestamp"][-1]
            for i in range(last_idx + 1, len(position_type_info_timestamp)):
                if i > last_idx:
                    # Simple extrapolation assuming constant rate
                    position_type_info_timestamp[i] = transition_mapping[last_idx] + (end_time - transition_mapping[last_idx]) * (i - last_idx) / (len(position_type_info_timestamp) - 1 - last_idx)
                else:
                    position_type_info_timestamp[i] = transition_mapping[last_idx]
        else:
            # Fallback: no transitions found, use uniform distribution
            position_type_info_timestamp = list(uniform_timestamps)

        # Update the dictionary
        self.position_type_info["timestamp"] = position_type_info_timestamp

        
        # Interpolate/extrapolate timestamp for messages that don't have
        # header timestamp properly logged
        self.first_ts = self.ins_solution["timestamp"][0]
        self.last_ts = self.ins_solution["timestamp"][-1]

        # Check each data source and update time range or interpolate invalid timestamps
        data_sources = [
            ("code_lat_lon_ht", self.code_lat_lon_ht),
            ("pdop_info", self.pdop_info),
            ("position_sigma", self.position_sigma),
            ("position_time_info", self.position_time_info),
            ("position_type_info", self.position_type_info),
            ("received_base_info", self.received_base_info),
            ("ins_solution_rms", self.ins_solution_rms),
            ("navsat", self.navsat),
            ("gnss_data", self.gnss_data)
        ]

        for name, data_dict in data_sources:
            if "timestamp" in data_dict and len(data_dict["timestamp"]) > 0:
                timestamp_array = data_dict["timestamp"]
                timestamp_valid = timestamp_array[-1] > 0

                if timestamp_valid:
                    self.first_ts = np.minimum(self.first_ts, timestamp_array[0])
                    self.last_ts = np.maximum(self.last_ts, timestamp_array[-1])
                else:
                    # Interpolate timestamps for invalid data
                    data_dict["timestamp"] = np.linspace(self.first_ts, self.last_ts, len(timestamp_array))

        # Normalize all timestamps to start from zero
        self._normalize_timestamps_to_zero()

    def _normalize_timestamps_to_zero(self):
        """Subtract first_ts from all timestamps to start from zero."""
        data_sources = [
            self.position_sigma,
            self.gnss_data,
            self.ins_solution,
            self.ins_solution_rms,
            self.navsat,
            self.pdop_info,
            self.received_base_info,
            self.position_time_info,
            self.position_type_info,
            self.code_lat_lon_ht
        ]

        for data_dict in data_sources:
            if "timestamp" in data_dict and len(data_dict["timestamp"]) > 0:
                data_dict["timestamp"] = data_dict["timestamp"] - self.first_ts

        # Update the time range
        self.last_ts = self.last_ts - self.first_ts
        self.first_ts = 0.0

    def _convert_lists_to_numpy(self, data_dict):
        """
        Recursively converts all list values in a dict to numpy arrays.
        """
        for key, value in data_dict.items():
            if isinstance(value, list):
                # Check if list contains dicts (nested fields)
                if len(value) > 0 and isinstance(value[0], dict):
                    for v in value:
                        self._convert_lists_to_numpy(v)
                else:
                    data_dict[key] = np.array(value)
            elif isinstance(value, dict):
                self._convert_lists_to_numpy(value)



