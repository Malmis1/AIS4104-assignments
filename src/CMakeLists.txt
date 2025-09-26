find_package(Eigen3 CONFIG REQUIRED)

add_executable(test main.cpp)

target_link_libraries(test PRIVATE Eigen3::Eigen)
