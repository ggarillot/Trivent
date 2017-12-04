########################################################
# cmake file for building Marlin example Package
# @author Gerald Grenier, IPNL (modifications)
# @author Jan Engels, Desy IT (original)
########################################################


### LIBRARY #################################################################


# include directories
include_directories( ./include )


FILE(GLOB_RECURSE headers "include/*.h" "include/*.hpp")
FILE(GLOB_RECURSE sources "src/*.cc")


add_library( ${PROJECT_NAME} SHARED ${sources} ${headers} )


install( TARGETS ${PROJECT_NAME} LIBRARY DESTINATION lib )
install( DIRECTORY ${PROJECT_SOURCE_DIR}/include DESTINATION . FILES_MATCHING PATTERN "*.h")



# include directories
#INCLUDE_DIRECTORIES( ./include )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.hh" )

# add library
#AUX_SOURCE_DIRECTORY( ./src library_sources )
#ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
#INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )


