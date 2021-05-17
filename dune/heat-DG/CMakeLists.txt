#install headers
install(FILES IncludesDUNE.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/heat)
add_subdirectory("operators")
add_subdirectory("heat")
