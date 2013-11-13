
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -W -Wall -Wextra -std=c99 -O3")

cmake_minimum_required (VERSION 2.6)
project (bilu-linial)

include_directories (.)

add_executable (bilu-linial bilu-linial.c nauty24r2/nauty.o nauty24r2/nautil.o nauty24r2/naugraph.o nauty24r2/gtools.o)

target_link_libraries(bilu-linial lapacke)

