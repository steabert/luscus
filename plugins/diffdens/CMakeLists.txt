if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  project(diffdens)
  add_executable(diffdens diffdens.c open_luscus_file.c)
  add_definitions(-DWINDOWS)

  file(APPEND ${CONFIG_DIR}/calculation.rc "  libpath=${CONFIG_DIR} type=orbital extension=lus description=\"density difference\" plugin_name=run_densdif.bat ext_out=lus extrafile=yes extraargsym=-a extraargval=-1.0 extraargdes=weight  \n")
  file(APPEND ${CONFIG_DIR}/calculation.rc "  libpath=${CONFIG_DIR} type=orbital extension=lus description=\"density sum\" plugin_name=run_densdif.bat ext_out=lus extrafile=yes extraargsym=-a extraargval=1.0 extraargdes=weight  \n")

  install(TARGETS diffdens DESTINATION ${CONFIG_DIR})
  file(COPY run_densdif.bat DESTINATION ${CONFIG_DIR})
#  target_link_libraries()
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  project(diffdens.exe)
  add_executable(diffdens.exe diffdens.c open_luscus_file.c)
  add_definitions(-DLINUX)

  file(APPEND ${TMP_CONFIG_DIR}/calculation.rc "  libpath=${CMAKE_INSTALL_PREFIX}/bin type=orbital extension=lus description=\"density difference\" plugin_name=run_densdif.sh ext_out=lus extrafile=yes extraargsym=-a extraargval=-1.0 extraargdes=weight  \n")
  file(APPEND ${TMP_CONFIG_DIR}/calculation.rc "  libpath=${CMAKE_INSTALL_PREFIX}/bin type=orbital extension=lus description=\"density sum\" plugin_name=run_densdif.sh ext_out=lus extrafile=yes extraargsym=-a extraargval=1.0 extraargdes=weight  \n")

  install(TARGETS diffdens.exe DESTINATION bin)
  file(COPY run_densdif.sh DESTINATION bin)
#  target_link_libraries(gamess2lus.exe -lm)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")


