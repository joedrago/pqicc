cmake_minimum_required(VERSION 3.5)
project(pqicc)

add_subdirectory(liblcms2)

include_directories(liblcms2/include)
add_executable(pqicc
  main.c
  md5.c
)

target_link_libraries(pqicc lcms2)
