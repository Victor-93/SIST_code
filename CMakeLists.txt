cmake_minimum_required(VERSION 3.28)
project(code)

set(CMAKE_CXX_STANDARD 14)

include_directories(trans_compete)
include_directories(trans_three)

add_executable(code
        trans_compete/G_x.cpp
        trans_compete/G_x.h
        trans_compete/qsidd.cpp
        trans_compete/SIDD_1R.cpp
        trans_compete/SIDD_1R.h
        trans_compete/SIDD_2R.cpp
        trans_compete/SIDD_2R.h
        trans_compete/SIDD_3R.cpp
        trans_compete/SIDD_3R.h
        trans_compete/SIDD_4R.cpp
        trans_compete/SIDD_4R.h
        trans_compete/SIDD_Base.cpp
        trans_compete/SIDD_Base.h
        trans_compete/stat_1R.cpp
        trans_compete/stat_1R.h
        trans_three/G_x.cpp
        trans_three/G_x.h
        trans_three/qsidd.cpp
        trans_three/SIDD_1R.cpp
        trans_three/SIDD_1R.h
        trans_three/SIDD_2R.cpp
        trans_three/SIDD_2R.h
        trans_three/SIDD_3R.cpp
        trans_three/SIDD_3R.h
        trans_three/SIDD_4R.cpp
        trans_three/SIDD_4R.h
        trans_three/SIDD_Base.cpp
        trans_three/SIDD_Base.h
        trans_three/stat_1R.cpp
        trans_three/stat_1R.h)
