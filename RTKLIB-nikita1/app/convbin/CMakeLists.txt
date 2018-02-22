cmake_minimum_required(VERSION 3.2)

project (convbin)

add_subdirectory(./../../src "%{CMAKE_BINARY_DIR}/rtklib")

add_executable(convbin convbin.c)

target_link_libraries(convbin rtklib )
