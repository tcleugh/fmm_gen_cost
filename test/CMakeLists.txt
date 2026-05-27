# Portable, absolute data paths injected into the test binaries so they work
# regardless of the cwd they're invoked from (mirrors the polymatch_test
# convention; see CLAUDE.md "Path convention").
set(FMM_TEST_DATA_DIR     "${CMAKE_SOURCE_DIR}/test/data")
set(FMM_EXAMPLE_DATA_DIR  "${CMAKE_SOURCE_DIR}/example/data")

add_executable(algorithm_test algorithm_test.cpp
        $<TARGET_OBJECTS:CORE>
        $<TARGET_OBJECTS:ALGORITHM>)
target_link_libraries(algorithm_test ${GDAL_LIBRARIES} ${Boost_LIBRARIES})

add_executable(fmm_test fmm_test.cpp
        $<TARGET_OBJECTS:MM_OBJ>
        $<TARGET_OBJECTS:CORE>
        $<TARGET_OBJECTS:CONFIG>
        $<TARGET_OBJECTS:ALGORITHM>
        $<TARGET_OBJECTS:UTIL>
        $<TARGET_OBJECTS:IO>
        $<TARGET_OBJECTS:NETWORK>
        $<TARGET_OBJECTS:FMM_OBJ>)
target_link_libraries(fmm_test ${GDAL_LIBRARIES} ${Boost_LIBRARIES}
        ${OpenMP_CXX_LIBRARIES} ${OSMIUM_LIBRARIES})
target_compile_definitions(fmm_test PRIVATE
        FMM_TEST_DATA_DIR="${FMM_TEST_DATA_DIR}")

add_executable(network_graph_test network_graph_test.cpp
        $<TARGET_OBJECTS:CORE>
        $<TARGET_OBJECTS:CONFIG>
        $<TARGET_OBJECTS:ALGORITHM>
        $<TARGET_OBJECTS:UTIL>
        $<TARGET_OBJECTS:IO>
        $<TARGET_OBJECTS:NETWORK>)
target_link_libraries(network_graph_test ${GDAL_LIBRARIES} ${Boost_LIBRARIES}
        ${OpenMP_CXX_LIBRARIES} ${OSMIUM_LIBRARIES})
target_compile_definitions(network_graph_test PRIVATE
        FMM_TEST_DATA_DIR="${FMM_TEST_DATA_DIR}")

add_executable(network_test network_test.cpp
        $<TARGET_OBJECTS:CORE>
        $<TARGET_OBJECTS:CONFIG>
        $<TARGET_OBJECTS:ALGORITHM>
        $<TARGET_OBJECTS:UTIL>
        $<TARGET_OBJECTS:IO>
        $<TARGET_OBJECTS:NETWORK>)
target_link_libraries(network_test ${GDAL_LIBRARIES} ${OpenMP_CXX_LIBRARIES}
${OSMIUM_LIBRARIES})
target_compile_definitions(network_test PRIVATE
        FMM_TEST_DATA_DIR="${FMM_TEST_DATA_DIR}")

add_executable(weightmatch_test weightmatch_test.cpp
        $<TARGET_OBJECTS:MM_OBJ>
        $<TARGET_OBJECTS:CORE>
        $<TARGET_OBJECTS:CONFIG>
        $<TARGET_OBJECTS:ALGORITHM>
        $<TARGET_OBJECTS:UTIL>
        $<TARGET_OBJECTS:IO>
        $<TARGET_OBJECTS:NETWORK>
        $<TARGET_OBJECTS:WEIGHTMATCH_OBJ>)
target_link_libraries(weightmatch_test ${GDAL_LIBRARIES} ${Boost_LIBRARIES}
        ${OpenMP_CXX_LIBRARIES} ${OSMIUM_LIBRARIES})
target_compile_definitions(weightmatch_test PRIVATE
        FMM_TEST_DATA_DIR="${FMM_TEST_DATA_DIR}"
        FMM_EXAMPLE_DATA_DIR="${FMM_EXAMPLE_DATA_DIR}")

add_executable(polymatch_test polymatch_test.cpp
        $<TARGET_OBJECTS:MM_OBJ>
        $<TARGET_OBJECTS:CORE>
        $<TARGET_OBJECTS:CONFIG>
        $<TARGET_OBJECTS:ALGORITHM>
        $<TARGET_OBJECTS:UTIL>
        $<TARGET_OBJECTS:IO>
        $<TARGET_OBJECTS:NETWORK>
        $<TARGET_OBJECTS:WEIGHTMATCH_OBJ>
        $<TARGET_OBJECTS:POLYMATCH_OBJ>)
target_link_libraries(polymatch_test ${GDAL_LIBRARIES} ${Boost_LIBRARIES}
        ${OpenMP_CXX_LIBRARIES} ${OSMIUM_LIBRARIES})
target_compile_definitions(polymatch_test PRIVATE
        POLYMATCH_FIXTURE_DIR="${CMAKE_SOURCE_DIR}/test/data/polymatch"
        FMM_REAL_EXAMPLE_DIR="${CMAKE_SOURCE_DIR}/test/data/polymatch/real_example_area")

add_custom_target(tests
	DEPENDS algorithm_test network_test network_graph_test fmm_test weightmatch_test polymatch_test)
