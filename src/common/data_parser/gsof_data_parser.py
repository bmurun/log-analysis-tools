import numpy as np
import os
import sys
import click

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
            "latitude": [],
            "longitude": [],
            "altitude": [],
            "roll": [],
            "pitch": [],
            "yaw": [],
            "north_vel": [],
            "east_vel": [],
            "down_vel": []
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

                    north_stdev = data.sigma_east 
                    east_stdev = data.sigma_north 
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
                    self.ins_solution["latitude"].append(data.lla.latitude)
                    self.ins_solution["longitude"].append(data.lla.longitude)
                    self.ins_solution["altitude"].append(data.lla.altitude)
                    self.ins_solution["roll"].append(data.roll)
                    self.ins_solution["pitch"].append(data.pitch)
                    self.ins_solution["yaw"].append(data.heading)
                    self.ins_solution["north_vel"].append(data.velocity.north)
                    self.ins_solution["east_vel"].append(data.velocity.north)
                    self.ins_solution["down_vel"].append(data.velocity.down)
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
                else:
                    pass

    def _process_data(self):
        """ Post processing
        """

        for topic, data in self.mcap_data.items():
            
            if topic in MCAP_ROS2_GSOF_TOPICS:
                if "header" in data:
                    sec_list = data.header.stamp.sec
                    nanosec_list = data.header.stamp.nanosec
                    timestamp_list = np.array(sec_list) + np.array(nanosec_list) / 1e9
                    data["timestamp"] = timestamp_list
                    del data["header"]
                if "stamp" in data:
                    sec_list = data["stamp"]["sec"]
                    nanosec_list = data["stamp"]["nanosec"]
                    timestamp_list = np.array(sec_list) + np.array(nanosec_list) / 1e9
                    data["timestamp"] = timestamp_list
                    del data["stamp"]

                self._convert_lists_to_numpy(data)

        min_timestamp = np.inf
        # Subtract the first timestamp from all
        for k, v in self.mcap_data.items():
            try:
                min_timestamp_v = np.min(v["timestamp"])
                min_timestamp = np.minimum(min_timestamp, min_timestamp_v)
            except:
                continue

        for k, v in self.mcap_data.items():
            try:
                v["timestamp"] -= min_timestamp
            except:
                continue

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



