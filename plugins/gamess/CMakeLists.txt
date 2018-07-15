#set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_PREFIX_PATH} PARENT_SCOPE)

#set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} PARENT_SCOPE)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  project(gamess2lus)
  add_executable(gamess2lus gamess2lus.c mystring.c)
  target_link_libraries(gamess2lus -lm)
  file(APPEND ${TMP_CONFIG_DIR}\\plugin.rc "  libpath=${CONFIG_DIR} extension=out description=\"GAMESS output\" forward=gamess2lus.exe  \n")
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  project(gamess2lus.exe)
  add_executable(gamess2lus.exe gamess2lus.c mystring.c)
  target_link_libraries(gamess2lus.exe -lm)
  file(APPEND ${TMP_CONFIG_DIR}/plugin.rc "  libpath=${CMAKE_INSTALL_PREFIX}/bin extension=out description=\"GAMESS output\" forward=gamess2lus.exe  \n")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  install(TARGETS gamess2lus DESTINATION ${CONFIG_DIR})
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  install(TARGETS gamess2lus.exe DESTINATION bin)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

