#!/usr/bin/env python3
from launch import LaunchDescription
from launch.actions import ExecuteProcess
from ament_index_python.packages import get_package_share_directory
import os

def generate_launch_description():

    ridgeback_urdf = os.path.join(
        get_package_share_directory("ridgeback_description"),
        "urdf", "ridgeback.urdf"
    )

    req_string = (
        f'sdf_filename: "{ridgeback_urdf}" '
        'name: "ridgeback" '
        'pose: { position: {x: 0, y: 0, z: 0.25}, orientation: {x: 0, y: 0, z: 0, w: 1} }'
    )

    return LaunchDescription([
        # Launch Gazebo Harmonic empty world
        ExecuteProcess(
            cmd=["gz", "sim", "-r", "empty.sdf"],
            output="screen"
        ),

        # Spawn the Ridgeback URDF
        ExecuteProcess(
            cmd=[
                "gz", "service",
                "-s", "/world/empty/create",
                "--reqtype", "gz.msgs.EntityFactory",
                "--reptype", "gz.msgs.Boolean",
                "--timeout", "1000",
                "--req", req_string
            ],
            output="screen"
        ),
    ])
