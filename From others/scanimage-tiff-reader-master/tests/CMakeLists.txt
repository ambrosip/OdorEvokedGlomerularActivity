set(tgt read_frame_interval)
add_executable(${tgt} ${tgt}.c)
target_add_scanimage_tiff_reader(${tgt})
add_test(test_${tgt} ${tgt})