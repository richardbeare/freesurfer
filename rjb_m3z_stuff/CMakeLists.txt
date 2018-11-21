project(rjb_m3z_stuff)

include_directories(${FS_INCLUDE_DIRS})

add_executable(rjb_m3z_stuff rjb_m3z_stuff.cpp)
add_help(rjb_m3z_stuff rjb_m3z_stuff.help.xml)
target_link_libraries(rjb_m3z_stuff utils)

#add_test_script(NAME mri_warp_convert_test SCRIPT test.py DEPENDS mri_warp_convert)

install(TARGETS rjb_m3z_stuff DESTINATION bin)
