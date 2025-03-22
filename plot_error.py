import matplotlib.pyplot as plt
import numpy as np

from mcap.reader import make_reader
from rclpy.serialization import deserialize_message
from rosidl_runtime_py.utilities import get_message


mcap_path_0 = "/home/av_ws/output_eval/config_0/run_0/run_0_0.mcap"
mcap_path_1 = "/home/av_ws/output_eval/config_1/run_0/run_0_0.mcap"
mcap_path_2 = "/home/av_ws/output_eval/config_2/run_0/run_0_0.mcap"

x_pos = []
y_pos = []
z_pos = []

gps_type = []
is_indoor = []

# Need
# /tracked_pose (Cartographer)
# /odometry/gps (GPS)
# /odometry/pose_in_map (IO transition)
# /debug/localization/kinematic_state

def get_data(mcap_path):

    topics = [
        "/tracked_pose",
        "/odometry/gps",
        "/odometry/pose_in_map",
        "/debug/localization/kinematic_state"
    ]

    # Initialize a dict with topic names and empty lists
    data_dict = {topic: {
        "timestamp_s": [],
        "pos_x": [],
        "pos_y": [],
        "pos_z": [],
        "pos_x_stdev": [],
        "pos_y_stdev": [],
        "pos_z_stdev": []
    } for topic in topics}

    with open(mcap_path, "rb") as f:
        reader = make_reader(f)
        for schema, channel, message in reader.iter_messages(topics=topics):
            
            topic_name = channel.topic
            msg_str = schema.name
            msg_type = get_message(msg_str)
            data = deserialize_message(message.data, msg_type)

            data_dict[topic_name]["timestamp_s"].append(message.publish_time / 1e9)

            if (topic_name == "/tracked_pose"):

                pos = data.pose.position
                data_dict[topic_name]["pos_x"].append(pos.x)
                data_dict[topic_name]["pos_y"].append(pos.y)
                data_dict[topic_name]["pos_z"].append(pos.z)
            
            else:
                pos = data.pose.pose.position
                cov = data.pose.covariance
                data_dict[topic_name]["pos_x"].append(pos.x)
                data_dict[topic_name]["pos_y"].append(pos.y)
                data_dict[topic_name]["pos_z"].append(pos.z)

                data_dict[topic_name]["pos_x_stdev"].append(np.sqrt(cov[0]))
                data_dict[topic_name]["pos_y_stdev"].append(np.sqrt(cov[7]))
                data_dict[topic_name]["pos_z_stdev"].append(np.sqrt(cov[14]))

    return data_dict
        
log_data = {}
log_data["config_0"] = get_data(mcap_path_0)
log_data["config_1"] = get_data(mcap_path_1)
log_data["config_2"] = get_data(mcap_path_2)

log_data["config_1"]["diff_tracked_pos_x"] = []
log_data["config_1"]["diff_tracked_pos_y"] = []
log_data["config_1"]["diff_tracked_pos_z"] = []

for i in range(len(log_data["config_1"]["/tracked_pose"]["timestamp_s"])):
    time_1 = log_data["config_1"]["/tracked_pose"]["timestamp_s"][i]
    posx_1 = log_data["config_1"]["/tracked_pose"]["pos_x"][i]
    posy_1 = log_data["config_1"]["/tracked_pose"]["pos_y"][i]
    posz_1 = log_data["config_1"]["/tracked_pose"]["pos_z"][i]

    ind = np.argmin(np.abs(np.array(log_data["config_0"]["/tracked_pose"]["timestamp_s"]) - time_1))

    diff_x = np.abs(log_data["config_0"]["/tracked_pose"]["pos_x"][ind] - posx_1)
    diff_y = np.abs(log_data["config_0"]["/tracked_pose"]["pos_y"][ind] - posy_1)
    diff_z = np.abs(log_data["config_0"]["/tracked_pose"]["pos_z"][ind] - posz_1)

    log_data["config_1"]["diff_tracked_pos_x"].append(diff_x)
    log_data["config_1"]["diff_tracked_pos_y"].append(diff_y)
    log_data["config_1"]["diff_tracked_pos_z"].append(diff_z)

log_data["config_2"]["diff_tracked_pos_x"] = []
log_data["config_2"]["diff_tracked_pos_y"] = []
log_data["config_2"]["diff_tracked_pos_z"] = []

for i in range(len(log_data["config_2"]["/tracked_pose"]["timestamp_s"])):
    time_1 = log_data["config_2"]["/tracked_pose"]["timestamp_s"][i]
    posx_1 = log_data["config_2"]["/tracked_pose"]["pos_x"][i]
    posy_1 = log_data["config_2"]["/tracked_pose"]["pos_y"][i]
    posz_1 = log_data["config_2"]["/tracked_pose"]["pos_z"][i]

    ind = np.argmin(np.abs(np.array(log_data["config_0"]["/tracked_pose"]["timestamp_s"]) - time_1))

    diff_x = np.abs(log_data["config_0"]["/tracked_pose"]["pos_x"][ind] - posx_1)
    diff_y = np.abs(log_data["config_0"]["/tracked_pose"]["pos_y"][ind] - posy_1)
    diff_z = np.abs(log_data["config_0"]["/tracked_pose"]["pos_z"][ind] - posz_1)

    log_data["config_2"]["diff_tracked_pos_x"].append(diff_x)
    log_data["config_2"]["diff_tracked_pos_y"].append(diff_y)
    log_data["config_2"]["diff_tracked_pos_z"].append(diff_z)

    
##################################################
fig, ax = plt.subplots(3, 1, sharex=True)
fig.suptitle("/tracked_pose Diff vs. config_0", y=0.98)

ax[0].set_ylabel("X Abs Diff [m]")
ax[1].set_ylabel("Y Abs Diff [m]")
ax[2].set_ylabel("Z Abs Diff [m]")
ax[2].set_xlabel("Time [s]")
ax[0].grid(True)
ax[1].grid(True)
ax[2].grid(True)

colors = ["blue", "green", "red"]
alphas = [0.9, 1.0, 0.6]
linewidths = [1.8, 1.5, 2.5]

for i, (key, data) in enumerate(log_data.items()):

    if "diff_tracked_pos_x" in data.keys():
        ax[0].plot(data["/tracked_pose"]["timestamp_s"], data["diff_tracked_pos_x"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i], label=key)
        ax[1].plot(data["/tracked_pose"]["timestamp_s"], data["diff_tracked_pos_y"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
        ax[2].plot(data["/tracked_pose"]["timestamp_s"], data["diff_tracked_pos_z"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])


lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.92), ncol=4)
fig.tight_layout(rect=[0, 0, 1, 0.9])
fig.savefig("diff_tracked_pose.png", dpi=400, bbox_inches='tight')




###################### GPS ##########################
log_data["config_1"]["diff_gps_pos_x"] = []
log_data["config_1"]["diff_gps_pos_y"] = []
log_data["config_1"]["diff_gps_pos_z"] = []

for i in range(len(log_data["config_1"]["/odometry/gps"]["timestamp_s"])):
    time_1 = log_data["config_1"]["/odometry/gps"]["timestamp_s"][i]
    posx_1 = log_data["config_1"]["/odometry/gps"]["pos_x"][i]
    posy_1 = log_data["config_1"]["/odometry/gps"]["pos_y"][i]
    posz_1 = log_data["config_1"]["/odometry/gps"]["pos_z"][i]

    ind = np.argmin(np.abs(np.array(log_data["config_0"]["/odometry/gps"]["timestamp_s"]) - time_1))

    diff_x = np.abs(log_data["config_0"]["/odometry/gps"]["pos_x"][ind] - posx_1)
    diff_y = np.abs(log_data["config_0"]["/odometry/gps"]["pos_y"][ind] - posy_1)
    diff_z = np.abs(log_data["config_0"]["/odometry/gps"]["pos_z"][ind] - posz_1)

    log_data["config_1"]["diff_gps_pos_x"].append(diff_x)
    log_data["config_1"]["diff_gps_pos_y"].append(diff_y)
    log_data["config_1"]["diff_gps_pos_z"].append(diff_z)

log_data["config_2"]["diff_gps_pos_x"] = []
log_data["config_2"]["diff_gps_pos_y"] = []
log_data["config_2"]["diff_gps_pos_z"] = []

for i in range(len(log_data["config_2"]["/odometry/gps"]["timestamp_s"])):
    time_1 = log_data["config_2"]["/odometry/gps"]["timestamp_s"][i]
    posx_1 = log_data["config_2"]["/odometry/gps"]["pos_x"][i]
    posy_1 = log_data["config_2"]["/odometry/gps"]["pos_y"][i]
    posz_1 = log_data["config_2"]["/odometry/gps"]["pos_z"][i]

    ind = np.argmin(np.abs(np.array(log_data["config_0"]["/odometry/gps"]["timestamp_s"]) - time_1))

    diff_x = np.abs(log_data["config_0"]["/odometry/gps"]["pos_x"][ind] - posx_1)
    diff_y = np.abs(log_data["config_0"]["/odometry/gps"]["pos_y"][ind] - posy_1)
    diff_z = np.abs(log_data["config_0"]["/odometry/gps"]["pos_z"][ind] - posz_1)

    log_data["config_2"]["diff_gps_pos_x"].append(diff_x)
    log_data["config_2"]["diff_gps_pos_y"].append(diff_y)
    log_data["config_2"]["diff_gps_pos_z"].append(diff_z)

    
##################################################
fig, ax = plt.subplots(3, 1, sharex=True)
fig.suptitle("/odometry/gps Diff vs. config_0", y=0.98)

ax[0].set_ylabel("X Abs Diff [m]")
ax[1].set_ylabel("Y Abs Diff [m]")
ax[2].set_ylabel("Z Abs Diff [m]")
ax[2].set_xlabel("Time [s]")
ax[0].grid(True)
ax[1].grid(True)
ax[2].grid(True)

colors = ["blue", "green", "red"]
alphas = [0.9, 1.0, 0.6]
linewidths = [1.8, 1.5, 2.5]

for i, (key, data) in enumerate(log_data.items()):

    if "diff_tracked_pos_x" in data.keys():
        ax[0].plot(data["/odometry/gps"]["timestamp_s"], data["diff_gps_pos_x"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i], label=key)
        ax[1].plot(data["/odometry/gps"]["timestamp_s"], data["diff_gps_pos_y"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
        ax[2].plot(data["/odometry/gps"]["timestamp_s"], data["diff_gps_pos_z"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])


lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.92), ncol=4)
fig.tight_layout(rect=[0, 0, 1, 0.9])
fig.savefig("diff_gps_pos.png", dpi=400, bbox_inches='tight')


#################### EKF ##############
log_data["config_1"]["diff_ekf_pos_x"] = []
log_data["config_1"]["diff_ekf_pos_y"] = []
log_data["config_1"]["diff_ekf_pos_z"] = []

for i in range(len(log_data["config_1"]["/debug/localization/kinematic_state"]["timestamp_s"])):
    time_1 = log_data["config_1"]["/debug/localization/kinematic_state"]["timestamp_s"][i]
    posx_1 = log_data["config_1"]["/debug/localization/kinematic_state"]["pos_x"][i]
    posy_1 = log_data["config_1"]["/debug/localization/kinematic_state"]["pos_y"][i]
    posz_1 = log_data["config_1"]["/debug/localization/kinematic_state"]["pos_z"][i]

    ind = np.argmin(np.abs(np.array(log_data["config_0"]["/debug/localization/kinematic_state"]["timestamp_s"]) - time_1))

    diff_x = np.abs(log_data["config_0"]["/debug/localization/kinematic_state"]["pos_x"][ind] - posx_1)
    diff_y = np.abs(log_data["config_0"]["/debug/localization/kinematic_state"]["pos_y"][ind] - posy_1)
    diff_z = np.abs(log_data["config_0"]["/debug/localization/kinematic_state"]["pos_z"][ind] - posz_1)

    log_data["config_1"]["diff_ekf_pos_x"].append(diff_x)
    log_data["config_1"]["diff_ekf_pos_y"].append(diff_y)
    log_data["config_1"]["diff_ekf_pos_z"].append(diff_z)

log_data["config_2"]["diff_ekf_pos_x"] = []
log_data["config_2"]["diff_ekf_pos_y"] = []
log_data["config_2"]["diff_ekf_pos_z"] = []

for i in range(len(log_data["config_2"]["/debug/localization/kinematic_state"]["timestamp_s"])):
    time_1 = log_data["config_2"]["/debug/localization/kinematic_state"]["timestamp_s"][i]
    posx_1 = log_data["config_2"]["/debug/localization/kinematic_state"]["pos_x"][i]
    posy_1 = log_data["config_2"]["/debug/localization/kinematic_state"]["pos_y"][i]
    posz_1 = log_data["config_2"]["/debug/localization/kinematic_state"]["pos_z"][i]

    ind = np.argmin(np.abs(np.array(log_data["config_0"]["/debug/localization/kinematic_state"]["timestamp_s"]) - time_1))

    diff_x = np.abs(log_data["config_0"]["/debug/localization/kinematic_state"]["pos_x"][ind] - posx_1)
    diff_y = np.abs(log_data["config_0"]["/debug/localization/kinematic_state"]["pos_y"][ind] - posy_1)
    diff_z = np.abs(log_data["config_0"]["/debug/localization/kinematic_state"]["pos_z"][ind] - posz_1)

    log_data["config_2"]["diff_ekf_pos_x"].append(diff_x)
    log_data["config_2"]["diff_ekf_pos_y"].append(diff_y)
    log_data["config_2"]["diff_ekf_pos_z"].append(diff_z)

    
##################################################
fig, ax = plt.subplots(3, 1, sharex=True)
fig.suptitle("/localization/kinematic_state Diff vs. config_0", y=0.98)

ax[0].set_ylabel("X Abs Diff [m]")
ax[1].set_ylabel("Y Abs Diff [m]")
ax[2].set_ylabel("Z Abs Diff [m]")
ax[2].set_xlabel("Time [s]")
ax[0].grid(True)
ax[1].grid(True)
ax[2].grid(True)

colors = ["blue", "green", "red"]
alphas = [0.9, 1.0, 0.6]
linewidths = [1.8, 1.5, 2.5]

for i, (key, data) in enumerate(log_data.items()):

    if "diff_tracked_pos_x" in data.keys():
        ax[0].plot(data["/debug/localization/kinematic_state"]["timestamp_s"], data["diff_ekf_pos_x"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i], label=key)
        ax[1].plot(data["/debug/localization/kinematic_state"]["timestamp_s"], data["diff_ekf_pos_y"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
        ax[2].plot(data["/debug/localization/kinematic_state"]["timestamp_s"], data["diff_ekf_pos_z"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])


lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.92), ncol=4)
fig.tight_layout(rect=[0, 0, 1, 0.9])
fig.savefig("diff_ekf_pos.png", dpi=400, bbox_inches='tight')




# ##################################################
# fig, ax = plt.subplots(3, 1, sharex=True)
# fig.suptitle("/odometry/gps", y=0.98)

# ax[0].set_ylabel("X Pos Diff [m]")
# ax[1].set_ylabel("Y Pos Diff [m]")
# ax[2].set_ylabel("Z Pos Diff [m]")
# ax[2].set_xlabel("Time [s]")
# ax[0].grid(True)
# ax[1].grid(True)
# ax[2].grid(True)

# ax[0].plot(data["/odometry/gps"]["timestamp_s"], data["/odometry/gps"]["pos_x"], "-", label=key, color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
# ax[1].plot(data["/odometry/gps"]["timestamp_s"], data["/odometry/gps"]["pos_y"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
# ax[2].plot(data["/odometry/gps"]["timestamp_s"], data["/odometry/gps"]["pos_z"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])


# for i, (key, data) in enumerate(log_data.items()):

#     ax[0].plot(data["/odometry/gps"]["timestamp_s"], data["/odometry/gps"]["pos_x"], "-", label=key, color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
#     ax[1].plot(data["/odometry/gps"]["timestamp_s"], data["/odometry/gps"]["pos_y"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
#     ax[2].plot(data["/odometry/gps"]["timestamp_s"], data["/odometry/gps"]["pos_z"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])


# lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
# lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

# fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.92), ncol=4)
# fig.tight_layout(rect=[0, 0, 1, 0.9])
# fig.savefig("odometry_gps.png", dpi=400, bbox_inches='tight')

# ##################################################
# fig, ax = plt.subplots(3, 1, sharex=True)
# fig.suptitle("/odometry/pose_in_map", y=0.98)

# ax[0].set_ylabel("X Pos [m]")
# ax[1].set_ylabel("Y Pos [m]")
# ax[2].set_ylabel("Z Pos [m]")
# ax[2].set_xlabel("Time [s]")
# ax[0].grid(True)
# ax[1].grid(True)
# ax[2].grid(True)

# for i, (key, data) in enumerate(log_data.items()):

#     ax[0].plot(data["/odometry/pose_in_map"]["timestamp_s"], data["/odometry/pose_in_map"]["pos_x"], "-", label=key, color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
#     ax[1].plot(data["/odometry/pose_in_map"]["timestamp_s"], data["/odometry/pose_in_map"]["pos_y"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
#     ax[2].plot(data["/odometry/pose_in_map"]["timestamp_s"], data["/odometry/pose_in_map"]["pos_z"], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])


# lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
# lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

# fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.92), ncol=4)
# fig.tight_layout(rect=[0, 0, 1, 0.9])
# fig.savefig("odometry_pose_in_map.png", dpi=400, bbox_inches='tight')

# ##################################################
# fig, ax = plt.subplots(3, 1, sharex=True)
# fig.suptitle("/debug/localization/kinematic_state", y=0.98)

# ax[0].set_ylabel("X Pos [m]")
# ax[1].set_ylabel("Y Pos [m]")
# ax[2].set_ylabel("Z Pos [m]")
# ax[2].set_xlabel("Time [s]")
# ax[0].grid(True)
# ax[1].grid(True)
# ax[2].grid(True)

# for i, (key, data) in enumerate(log_data.items()):

#     ax[0].plot(data["/debug/localization/kinematic_state"]["timestamp_s"][1:], data["/debug/localization/kinematic_state"]["pos_x"][1:], "-", label=key, color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
#     ax[1].plot(data["/debug/localization/kinematic_state"]["timestamp_s"][1:], data["/debug/localization/kinematic_state"]["pos_y"][1:], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])
#     ax[2].plot(data["/debug/localization/kinematic_state"]["timestamp_s"][1:], data["/debug/localization/kinematic_state"]["pos_z"][1:], "-", color=colors[i], alpha=alphas[i], linewidth=linewidths[i])


# lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
# lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]

# fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.92), ncol=4)
# fig.tight_layout(rect=[0, 0, 1, 0.9])
# fig.savefig("localization_kinematic_state.png", dpi=400, bbox_inches='tight')
