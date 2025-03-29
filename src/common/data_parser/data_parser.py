import numpy as np
import os
import sys
import click

from mcap.reader import make_reader
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message

from common.data_parser.topics import MCAP_ROS2_TOPICS

class DataParser:
    def __init__(self, log_path: str) -> None:
        if not os.path.isfile(log_path):
            click.echo(click.style("Input file doesn't exist!", fg="red"))
            return
        self.mcap_data = {}
        self._parse_mcap(log_path)
        self._process_data()

    def _parse_mcap(self, log_path) -> None:

        with open(log_path, "rb") as f:
            reader = make_reader(f)

            for schema, channel, message in reader.iter_messages(topics=MCAP_ROS2_TOPICS):

                topic_name = channel.topic
                msg_str = schema.name
                msg_type = get_message(msg_str)
                data = deserialize_message(message.data, msg_type)
                
                if topic_name not in self.mcap_data:
                    self.mcap_data[topic_name] = {}


                self._populate_fields(topic_name, data, self.mcap_data[topic_name])

    def _populate_fields(self, topic_name, msg, storage):
        for field_name in msg.get_fields_and_field_types().keys():
            value = getattr(msg, field_name)

            if hasattr(value, 'get_fields_and_field_types'):
                if field_name not in storage:
                    storage[field_name] = {}
                self._populate_fields(topic_name, value, storage[field_name])

            elif isinstance(value, list) and len(value) > 0 and hasattr(value[0], 'get_fields_and_field_types'):
                if field_name not in storage:
                    storage[field_name] = []
                nested_list = []
                for v in value:
                    nested_storage = {}
                    self._populate_fields(topic_name, v, nested_storage)
                    nested_list.append(nested_storage)
                storage[field_name].append(nested_list)

            else:
                if field_name not in storage:
                    storage[field_name] = []
                storage[field_name].append(value)

    def _process_data(self):
        """ Post processing
        """

        for topic, data in self.mcap_data.items():
            if "header" in data:
                sec_list = data["header"]["stamp"]["sec"]
                nanosec_list = data["header"]["stamp"]["nanosec"]
                timestamp_list = np.array(sec_list) + np.array(nanosec_list) / 1e9
                data["timestamp"] = timestamp_list
                del data["header"]

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



