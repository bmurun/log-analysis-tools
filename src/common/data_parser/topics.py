MCAP_ROS2_TOPICS = [
    "/tracked_pose",
    "/odometry/gps",
    "/odometry/pose_in_map",
    "/debug/localization/kinematic_state",
    "/localization/kinematic_state",
    "/microstrain/imu/data",
    "/lvx_client/imu_raw",
    "/lvx_client/gsof/ins_solution_rms_50",
    "/lvx_client/gsof/ins_solution_49",
    "/localization/is_ego_indoors",
    "/debug/localization/is_ego_indoors",
    "/localization/debug/cartographer_convergence_status",
    # New interface
    "/localization/ekf/input/odometry",
    "/localization/ekf/input/cartographer",
    "/localization/ekf/input/imu",
    "/localization/ekf/input/gnss",
    "/localization/debug/localization_mode",
    "/odometry/gps_raw_cov",
    # Test
    "/output/test_cartographer_convergence_status_from_gnss"
]

MCAP_ROS2_NIS_TOPICS = [
    "/localization/ekf/debug/estimator_aid_src_cartographer_pose",
    "/localization/ekf/debug/estimator_aid_src_gnss_pose",
    "/localization/ekf/debug/estimator_aid_src_imu_pose",
    "/localization/ekf/debug/estimator_aid_src_kiss_odometry_pose",
    "/localization/ekf/debug/estimator_aid_src_odometry_twist",
]