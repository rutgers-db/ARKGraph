set(CMAKE_CXX_STANDARD 11)



add_executable(groundtruth calculate_groundtruth.cc)
target_link_libraries(groundtruth UTIL kgrapha ${Boost_LIBRARIES} )




add_executable(exp exp.cc)
target_link_libraries(exp UTIL kgrapha ${Boost_LIBRARIES} )
