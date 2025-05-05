# CMakeList.txt: Top-level CMake project file. Do global configuration and include subprojects here.
#
cmake_minimum_required (VERSION 3.8)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable Hot Reload for MSVC compilers if supported. (https://docs.tizen.org/application/dotnet/guides/hot-reloading/overview/)
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("codiciTesi")

# Includere i sottoprogetti.
add_subdirectory ("codiciTesi")
