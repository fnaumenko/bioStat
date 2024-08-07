cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)

set(CMAKE_SYSTEM_VERSION 10.0 CACHE STRING "" FORCE)
set(PROJECT_NAME biostat)
project(${PROJECT_NAME} CXX)

# Set target arch type if empty. Visual studio solution generator provides it.
if(NOT CMAKE_VS_PLATFORM_NAME)
    set(CMAKE_VS_PLATFORM_NAME "x64")
endif()
message("${CMAKE_VS_PLATFORM_NAME} architecture in use")

if(NOT ("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64"
     OR "${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86"))
    message(FATAL_ERROR "${CMAKE_VS_PLATFORM_NAME} arch is not supported!")
endif()

# Global configuration types
set(CMAKE_CONFIGURATION_TYPES
    "Debug"
    "Release"
    CACHE STRING "" FORCE
)

# Global compiler options
if(MSVC)
    # remove default flags provided with CMake for MSVC
    set(CMAKE_CXX_FLAGS "")
    set(CMAKE_CXX_FLAGS_DEBUG "")
    set(CMAKE_CXX_FLAGS_RELEASE "")
	add_compile_options(/EHsc)
endif()

# Global linker options
if(MSVC)
    # remove default flags provided with CMake for MSVC
    set(CMAKE_EXE_LINKER_FLAGS "")
    set(CMAKE_MODULE_LINKER_FLAGS "")
    set(CMAKE_SHARED_LINKER_FLAGS "")
    set(CMAKE_STATIC_LINKER_FLAGS "")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS_DEBUG "${CMAKE_MODULE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_STATIC_LINKER_FLAGS_DEBUG "${CMAKE_STATIC_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS_RELEASE "${CMAKE_MODULE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS}")
    set(CMAKE_STATIC_LINKER_FLAGS_RELEASE "${CMAKE_STATIC_LINKER_FLAGS}")
endif()

# Nuget packages function stub.
function(use_package TARGET PACKAGE VERSION)
    message(WARNING "No implementation of use_package. Create yours. "
                    "Package \"${PACKAGE}\" with version \"${VERSION}\" "
                    "for target \"${TARGET}\" is ignored!")
endfunction()

# Common utils
include(CMake/Utils.cmake)

# Additional Global Settings(add specific info there)
include(CMake/GlobalSettingsInclude.cmake OPTIONAL)

# Use solution folders feature
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Sub-projects
unset(SUB_PROJECT CACHE)
set(SUB_PROJECT "true" CACHE STRING "set when projects below are biuld as a subprojects")
unset(BASIC_DIR CACHE)
set(BASIC_DIR ../../bioinfoBasic CACHE STRING "set common basic folder")

# ZLIB for sub-projects
if(MSVC)
	set(ZLIB_DIR /ProgramData/zlib1.3)
	set(ZLIB_LIBRARY_DEBUG ${ZLIB_DIR}/Debug/zlibd.lib)
	set(ZLIB_LIBRARY_RELEASE ${ZLIB_DIR}/Release/zlib.lib)
	set(ZLIB_INCLUDE_DIR ${ZLIB_DIR})
endif()
find_package(ZLIB REQUIRED)
if(NOT ZLIB_FOUND)
	message(WARNING ">> ZLIB not found!")
endif()
# message(">> ZLIB_INCLUDE_DIRS: ${ZLIB_INCLUDE_DIRS}")

add_subdirectory(src ${CMAKE_BINARY_DIR}/build)
add_subdirectory(bioCC/src ${CMAKE_BINARY_DIR}/bioCC)
add_subdirectory(callDist/src ${CMAKE_BINARY_DIR}/callDist)
add_subdirectory(fqStatN/src ${CMAKE_BINARY_DIR}/fqStatN)
add_subdirectory(vAlign/src ${CMAKE_BINARY_DIR}/vAlign)
add_subdirectory(FGStest/src ${CMAKE_BINARY_DIR}/FGStest)

# add_subdirectory(../callDistTest/src ${CMAKE_BINARY_DIR}/callDistTest)


