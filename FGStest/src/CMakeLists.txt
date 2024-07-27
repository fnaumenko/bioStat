cmake_minimum_required(VERSION 3.16.0)

set(PROJECT_NAME FGStest)
project(${PROJECT_NAME})

add_compile_definitions(
	MY_DEBUG
	_FEATURES
	_FEATR_SCORE
	_DUP_OUTPUT		# duplicated output
	)

# Source groups

if(SUB_PROJECT STREQUAL "true")		# sub-project
	message(">> ${PROJECT_NAME} is builded as a sub-project")
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../build)
else()								# standalone project
	message(">> ${PROJECT_NAME} is builded as a standalone project")
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY build)
	set(BASIC_DIR ../../../bioinfoBasic)
	# a standalone project is supported for development purposes 
	# and assumes a manual copying of the ZLIB.dll
	if(MSVC)
		set(ZLIB_DIR /ProgramData/zlib1.3)
		set(ZLIB_LIBRARY_DEBUG ${ZLIB_DIR}/Debug/zlibd.lib)
		set(ZLIB_LIBRARY_RELEASE ${ZLIB_DIR}/Release/zlib.lib)
		set(ZLIB_INCLUDE_DIR ${ZLIB_DIR})
	endif()
	find_package(ZLIB REQUIRED)	
endif()

# ZLIB
if(ZLIB_FOUND)		# find_package() is called in main project CMakeList
    add_compile_definitions(_ZLIB _BAM)
    include_directories(${ZLIB_INCLUDE_DIR})
endif()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/${BASIC_DIR} )

set(Header_Files
    ${BASIC_DIR}/common.h
    ${BASIC_DIR}/ChromData.h
    ${BASIC_DIR}/DataReader.h
    ${BASIC_DIR}/CrossRgns.h
    ${BASIC_DIR}/Options.h
    ${BASIC_DIR}/Features.h
    ${BASIC_DIR}/TxtFile.h
    Main.h
)
source_group("Header Files" FILES ${Header_Files})

set(Source_Files
    ${BASIC_DIR}/common.cpp
    ${BASIC_DIR}/ChromData.cpp
    ${BASIC_DIR}/DataReader.cpp
    ${BASIC_DIR}/Options.cpp
    ${BASIC_DIR}/Features.cpp
    ${BASIC_DIR}/TxtFile.cpp
    Main.cpp
)
if(ZLIB_FOUND)
	set(BAM_DIR ${BASIC_DIR}/bam)
	set(Source_Files ${Source_Files}
	    ${BAM_DIR}/BamReader.cpp
		${BAM_DIR}/BGZF.cpp
	)
endif()

source_group("Source Files" FILES ${Source_Files})

set(ALL_FILES
    ${Header_Files}
    ${Source_Files}
)

# Target

add_executable(${PROJECT_NAME} ${ALL_FILES})

if(ZLIB_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE ZLIB::ZLIB)
	# zlib.dll is copied only in this subproject
	add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different ${ZLIB_DIR}/Release/zlib.dll ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/
		)		
endif()

set(ROOT_NAMESPACE ${PROJECT_NAME})

set_target_properties(${PROJECT_NAME} PROPERTIES
    VS_GLOBAL_KEYWORD "Win32Proj"
)
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    set_target_properties(${PROJECT_NAME} PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION_DEBUG   "FALSE"
        INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
    )
elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
    set_target_properties(${PROJECT_NAME} PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION_DEBUG   "FALSE"
        INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
    )
endif()

# Compile definitions
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Debug>:"
            "_DEBUG"
        ">"
        "$<$<CONFIG:Release>:"
            "NDEBUG"
        ">"
        "_CONSOLE;"
        "UNICODE;"
        "_UNICODE"
        "_CRT_SECURE_NO_WARNINGS"
    )
elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Debug>:"
            "_DEBUG;"
            ""
        ">"
        "$<$<CONFIG:Release>:"
            "NDEBUG"
        ">"
        "WIN32;"
        "_CONSOLE;"
        "_CRT_SECURE_NO_WARNINGS"
    )
endif()

# Compile and link options
if(MSVC)
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>:
                /Oi;
                /Gy
            >
            /permissive-;
            /sdl;
            /W3;
            ${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
            ${DEFAULT_CXX_EXCEPTION_HANDLING};
            /Y-
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Debug>:
                /JMC;
                /Od;
                /Oi-
            >
            $<$<CONFIG:Release>:
                /Oi;
                /Gy
            >
            /permissive-;
            /std:c++17;
            /sdl;
            /W3;
            ${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
            ${DEFAULT_CXX_EXCEPTION_HANDLING};
             /Y-
       )
    endif()
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_link_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Debug>:
                /INCREMENTAL
            >
            $<$<CONFIG:Release>:
                /OPT:REF;
                /OPT:ICF;
                /INCREMENTAL:NO
            >
            /DEBUG;
            /SUBSYSTEM:CONSOLE
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
        target_link_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Debug>:
                /INCREMENTAL
            >
            $<$<CONFIG:Release>:
                /OPT:REF;
                /OPT:ICF;
                /INCREMENTAL:NO
            >
            /DEBUG;
            /SUBSYSTEM:CONSOLE
        )
    endif()
endif()
