cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)
project(ctbn-test-executables)

include_directories(${CTBN_INCLUDE_DIRS})

set(TEST_EXECUTABLES clique_unit_test decision exactfbquery factoredunifquery gibbsauxprequery gibbsprequery gibbsquery glpktest importancequery meanfieldquery sem structsearch structurescore unit_testing_graphic)

set(AUX_EXECUTABLES makequeryinput makedrugtst maketoroid)

set(DATA_FILES unit_test_traj.graph unit_test_traj.eps)

add_custom_command(OUTPUT queryinput.data
	COMMAND ./makequeryinput > queryinput.data
	DEPENDS makequeryinput
	)
add_custom_command(OUTPUT queryinput-toroid.data
	COMMAND ./maketoroid > queryinput-toroid.data
	DEPENDS maketoroid
	)
add_custom_command(OUTPUT drug.ctbn
	COMMAND ./makedrugtst > drug.ctbn
	DEPENDS makedrugtst
	)
add_custom_command(OUTPUT drug.ctbndyn
	COMMAND ./makedrugtst -Dsavedyn 1 > drug.ctbndyn
	DEPENDS makedrugtst
	)
add_custom_target(datafilestst ALL
	DEPENDS drug.ctbn drug.ctbndyn queryinput.data queryinput-toroid.data)

foreach(data_file ${DATA_FILES})
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${data_file}
	    ${CMAKE_CURRENT_BINARY_DIR}/${data_file} COPYONLY)
endforeach()

foreach(aux_executable ${AUX_EXECUTABLES})
    add_executable(${aux_executable} ${aux_executable}.cc)
    target_link_libraries(${aux_executable} ${CTBN_LIBRARIES})
endforeach()

foreach(test_executable ${TEST_EXECUTABLES})
    add_executable(${test_executable} ${test_executable}.cc)
    target_link_libraries(${test_executable} ${CTBN_LIBRARIES})
    add_test(${test_executable} ${test_executable})
endforeach()
