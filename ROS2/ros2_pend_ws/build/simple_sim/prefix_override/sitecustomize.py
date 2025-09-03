import sys
if sys.prefix == '/usr':
    sys.real_prefix = sys.prefix
    sys.prefix = sys.exec_prefix = '/home/bauer/REPOS/Robotics/ROS2/ros2_pend_ws/install/simple_sim'
