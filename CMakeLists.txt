#project name
PROJECT(tpmsfem)

#requierment of camke version
cmake_minimum_required(VERSION 3.0.0)

# Suppress MSVC CRT warnings.
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    #add_definitions(/Za)
    add_definitions(/MP) # multi-core builds
endif()

include(cmake/DownloadProject.cmake)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#add include
 
 set(TPMSFEM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

 file(GLOB INCLUDE_DIRS ${TPMSFEM_INCLUDE_DIR}/*.h ${TPMSFEM_INCLUDE_DIR}/*.cpp)

 add_library( tpmsfem ${INCLUDE_DIRS})
 target_include_directories(tpmsfem PUBLIC ${TPMSFEM_INCLUDE_DIR})

#tpmsgen

download_project(PROJ    tpmsGen
		 		 SOURCE_DIR   ${CMAKE_CURRENT_LIST_DIR}/external/tpmsGen
				 DOWNLOAD_DIR ${CMAKE_CURRENT_LIST_DIR}/external/.cache/tpmsGen
				 GIT_REPOSITORY      https://github.com/swannyPeng/tpmsGen.git
                 GIT_TAG             master			
		)
set( TPMSGEN_DIR "${PROJECT_SOURCE_DIR}/external/tpmsGen")
add_subdirectory("${TPMSGEN_DIR}" "tpmsGen")
target_include_directories(tpmsfem PUBLIC "${TPMSGEN_DIR}/include")
target_include_directories(tpmsfem PUBLIC ${TPMSGEN_DIR}/external/eigen)
target_link_libraries(tpmsfem INTERFACE tpmsgen)

#tetgen
set(TETGEN_DIR "${CMAKE_CURRENT_LIST_DIR}/external/tetgen")
if(NOT TARGET tetgen)
    download_project( PROJ tetgen
    				 SOURCE_DIR   ${CMAKE_CURRENT_LIST_DIR}/external/tetgen
					 DOWNLOAD_DIR ${CMAKE_CURRENT_LIST_DIR}/external/.cache/tetgen
					 GIT_REPOSITORY https://github.com/libigl/tetgen.git
					 GIT_TAG        d2dcc33cb8551f16e302c8130ce12fa52992cd09
    	)

    add_subdirectory("${TETGEN_DIR}" "tetgen")
    target_include_directories(tpmsfem PUBLIC ${TETGEN_DIR})
    target_link_libraries(tpmsfem INTERFACE tetgen)
  endif()




#Add your project files
file(GLOB HEADFILES *.h)
file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME}_bin ${HEADFILES} ${SRCFILES} )
#target_include_directories(${PROJECT_NAME}_bin PUBLIC ${HEXFILL_INCLUDE_DIR} )
target_link_libraries(${PROJECT_NAME}_bin  tpmsfem)