cmake_minimum_required(VERSION 3.15)
project(PDE-06 LANGUAGES CXX)

# Set standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Source files
file(GLOB_RECURSE PROJECT_FILES "src/*.cpp" "src/*.h" "src/*.hpp")
add_executable(${PROJECT_NAME} ${PROJECT_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE src)

# MPI compiler.
# find_package(MPI REQUIRED)
# set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")

# Boost.
# find_package(Boost 1.72.0 REQUIRED
#   COMPONENTS filesystem iostreams serialization
#   HINTS ${BOOST_DIR} $ENV{BOOST_DIR} $ENV{mkBoostPrefix})
# message(STATUS "Using the Boost-${Boost_VERSION} configuration found at ${Boost_DIR}")
# message(STATUS)
# include_directories(${Boost_INCLUDE_DIRS})

# DealII
find_package(deal.II REQUIRED HINTS ${DEAL_II_DIR} $ENV{DEAL_II_DIR} $ENV{mkDealiiPrefix})
message(STATUS "Found deal.II version ${DEAL_II_VERSION}")

deal_ii_initialize_cached_variables()
deal_ii_setup_target(${PROJECT_NAME})

# Platform-specific settings
if (WIN32)
	message(STATUS "Configuring for Windows...")

	# Optional: Set Windows-specific compile definitions
	target_compile_definitions(${PROJECT_NAME} PRIVATE PLATFORM_WINDOWS)

	# Optional: Set MSVC flags
	if (MSVC)
		target_compile_options(${PROJECT_NAME} PRIVATE /W4)
	endif()

	if (CMAKE_GENERATOR MATCHES "Visual Studio")
		# Set startup project
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
		# Group by folder automatically in VS filters
		source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src PREFIX "Source Files" FILES ${PROJECT_FILES})
	endif()
elseif (UNIX)
	message(STATUS "Configuring for Linux...")

	target_compile_definitions(${PROJECT_NAME} PRIVATE PLATFORM_LINUX)

	# GCC/Clang warnings
	target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
endif()
