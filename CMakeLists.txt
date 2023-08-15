cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(FMM)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/dmolony3/SlicerFractionalMyocardialMass")
set(EXTENSION_CATEGORY "Cardiac")
set(EXTENSION_CONTRIBUTORS "David Molony (Georgia Heart Institute.)")
set(EXTENSION_DESCRIPTION "Extension for calculating the fractional myocardial mass and the myocardial mass at risk")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/dmolony3/SlicerFractionalMyocardialMass/main/FMM.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/dmolony3/SlicerFractionalMyocardialMass/main/Media/Mesh.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(FMM)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
