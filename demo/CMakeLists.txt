cmake_minimum_required(VERSION 2.6.4 FATAL_ERROR)
project(ctbn-demo-executables)

include_directories(${CTBN_INCLUDE_DIRS})
enable_testing()

set(DEMO_EXECUTABLES clique_example makechain gen_fulldata gen_partialdata learnparams sampleprocess using_graphic gibbssampquery exactquery isquery makeqinput makedrug)

add_custom_command(OUTPUT queryinput.data
		COMMAND makeqinput > queryinput.data
		DEPENDS makeqinput
		)

add_custom_command(OUTPUT drug.ctbn
		COMMAND makedrug > drug.ctbn
		DEPENDS makedrug
		)

add_custom_command(OUTPUT drug.ctbndyn
		COMMAND makedrug -Dsavedyn 1 > drug.ctbndyn
		DEPENDS makedrug
		)
add_custom_target(datafiles ALL DEPENDS drug.ctbn queryinput.data drug.ctbndyn)

foreach(demo_executable ${DEMO_EXECUTABLES})
    add_executable(${demo_executable} ${demo_executable}.cc)
    target_link_libraries(${demo_executable} ${CTBN_LIBRARIES})
endforeach()

