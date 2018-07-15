INCLUDE(FindPkgConfig)
INCLUDE(CheckFunctionExists)
cmake_minimum_required(VERSION 2.8.8)
project(luscus)

set(luscus_VERSION_MAJOR 0)
set(luscus_VERSION_MINOR 8)
set(luscus_VERSION_REVISION 5)

check_function_exists(strcasestr HAVE_STRCASESTR)
check_function_exists(strcasecmp HAVE_STRCASECMP)

configure_file("${PROJECT_SOURCE_DIR}/luscus.h.in" "${PROJECT_BINARY_DIR}/luscus.h")

#if($HAVE_STRCASESTR)
#  message(STATUS "found strcasestr!")
#endif($HAVE_STRCASESTR)
#if(HAVE_STRCASECMP)
#  message(STATUS "found strcasecmp!")
#endif(HAVE_STRCASECMP)

set(USE_GTK3 ON)

# CMAKE_INSTALL_PREFIX

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_definitions(-DLINUX)
  pkg_check_modules(GTK2 QUIET gtk+-2.0 gl glu)
#  if (EXISTS ${CMAKE_PREFIX_PATH})
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT MATCHES 1)
    set(CONFIG_DIR "/etc/luscus")
#    message(status " CMAKE_PREFIX_PATH NOT DEFINED!")    # DEBUG
  else ()
    set(CONFIG_DIR "$ENV{HOME}/.luscus")
#    message(status " CMAKE_PREFIX_PATH DEFINED!")        # DEBUG
  endif ()
  set(TMP_CONFIG_DIR ${CMAKE_CURRENT_BINARY_DIR}/luscusrc)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_definitions(-DWINDOWS)
  pkg_check_modules(GTK2 QUIET gtk+-2.0)
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT MATCHES 1)
    set(CONFIG_DIR "$ENV{ALLUSERSPROFILE}\\luscus")
  else ()
#  if (EXISTS ${CMAKE_PREFIX_PATH})
    set(CONFIG_DIR "$ENV{LOCALAPPDATA}\\luscus")
  endif ()
  set(TMP_CONFIG_DIR ${CMAKE_CURRENT_BINARY_DIR}\\luscusrc)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions(-DLINUX)
  pkg_check_modules(GTK2 QUIET gtk+-2.0 gl glu)
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT MATCHES 1)
    set(CONFIG_DIR "/etc/luscus")
#  if (EXISTS ${CMAKE_PREFIX_PATH})
  else()
    set(CONFIG_DIR "$ENV{HOME}/.luscus")
  endif ()
  set(TMP_CONFIG_DIR ${CMAKE_CURRENT_BINARY_DIR}/luscusrc)
else()
  message(FATAL_ERROR "ERROR: Can not determine operating system. Your operating system is most probably unsupported.")
endif()

#--------------- DEBUG ----------------
#message(status " CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT  = ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT}")
#message(status " CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")
#message(status " CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
#message(status " config dir = ${CONFIG_DIR}")
#--------------------------------------

file(MAKE_DIRECTORY ${TMP_CONFIG_DIR})

if(GTK2_FOUND)
  set(USE_GTK3 OFF)
#  message(STATUS "Found GTK2")
endif(GTK2_FOUND) 

if(NOT USE_GTK3)
  pkg_check_modules(GTKGLEXT QUIET gtkglext-1.0)
  if(GTKGLEXT_FOUND) 
    set(USE_GTKGLEXT ON)
  endif(GTKGLEXT_FOUND) 
endif(NOT USE_GTK3)

if(USE_GTK3)
  pkg_check_modules(GL QUIET gl)
  if(NOT GL_FOUND)
    message(FATAL_ERROR "ERROR: GL development libraries are missing. Install development version (usually a -dev package) of gl.")
  endif()
  pkg_check_modules(GLU QUIET glu)
  if(NOT GLU_FOUND)
    message(FATAL_ERROR "ERROR: GLU development libraries are missing. Install development version (usually a -dev package) of glu.")
  endif()
  pkg_check_modules(GTK3 QUIET gtk+-3.0)
  if(NOT GTK3_FOUND)
    message(FATAL_ERROR "ERROR: GTK development libraries are missing. Install development version (usually a -dev package) of gl, glu, and gtkglext. If you have a choice, use GTK2 and gtkglext, else try GTK3.")
  endif(NOT GTK3_FOUND)
endif(USE_GTK3)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_definitions(-DLINUX)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  add_definitions(-DWINDOWS)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions(-DLINUX)
else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  message(FATAL_ERROR "ERROR: Can not determine operating system. Your operating system is most probably unsupported.")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

#option(COMPILE_PLUGINS ON)
option(INSTALL_BABEL_PLUGINS ON)

if(HAVE_STRCASESTR)
add_executable(luscus main.c gtk_gui.c gv_menubar.c gv_system.c gv_notebook.c backup_file.c findsym.c gv_about.c gv_atoms.c gv_geom.c gv_fragments.c gv_gtk_atom_list.c gv_gtk_multiview.c gv_notebook.c luscus_gtk_util.c mystring.c screenshot.c gveps.c gvgrp.c findsym.c gv_gtk_help.c gv_gtk_xyz_editor.c pixeldata.c write_file.c read_file.c surface.c gv_gtk_make_graph.c luscus_gtk_handle.c mcube.c vstring.c)
else(HAVE_STRCASESTR)
add_executable(luscus main.c gtk_gui.c gv_menubar.c gv_system.c gv_notebook.c backup_file.c findsym.c gv_about.c gv_atoms.c gv_geom.c gv_fragments.c gv_gtk_atom_list.c gv_gtk_multiview.c gv_notebook.c luscus_gtk_util.c mystring.c screenshot.c gveps.c gvgrp.c findsym.c gv_gtk_help.c gv_gtk_xyz_editor.c pixeldata.c write_file.c read_file.c surface.c gv_gtk_make_graph.c luscus_gtk_handle.c mcube.c vstring.c stdext.c)
endif(HAVE_STRCASESTR)


if(USE_GTK3)
  add_definitions(-DGTK3 -D_GNU_SOURCE ${GTK3_CFLAGS})
  target_link_libraries(luscus ${GTK3_LDFLAGS} -lm -lX11)
elseif(USE_GTKGLEXT)
#  add_definitions(-g -pg -DEBUG -DGTK2 -DGTK_OLD -DGTK_GLEXT -D_GNU_SOURCE ${GTKGLEXT_CFLAGS})
  add_definitions(-DGTK2 -DGTK_OLD -DGTK_GLEXT -D_GNU_SOURCE ${GTKGLEXT_CFLAGS})
  target_link_libraries(luscus ${GTKGLEXT_LDFLAGS} -lm)
else(USE_GTK3)
  add_definitions(-DGTK2 -D_GNU_SOURCE ${GTK2_CFLAGS})
  target_link_libraries(luscus ${GTK2_LDFLAGS} -lm -lX11)
endif(USE_GTK3)

include_directories("${PROJECT_BINARY_DIR}")

#--------------------------------------------------
# THIS IS OLD PLUG-IN HANDLING

#if(NOT EXISTS ${CONFIG_DIR})
#  file(MAKE_DIRECTORY ${CONFIG_DIR})
#endif(NOT EXISTS ${CONFIG_DIR})
#
#if(EXISTS ${CONFIG_DIR}/plugin.rc)
#  file(REMOVE ${CONFIG_DIR}/plugin.rc)
##  file(WRITE ${CONFIG_DIR}/plugin.rc "  ")
#  message(WARNING " Found old config file: ${CONFIG_DIR}/plugin.rc; it will be removed and replaced with new one!")
#endif(EXISTS ${CONFIG_DIR}/plugin.rc)
#
#if(EXISTS ${CONFIG_DIR}/calculation.rc)
#  file(REMOVE ${CONFIG_DIR}/calculation.rc)
##  file(WRITE ${CONFIG_DIR}/calculation.rc "  ")
#  message(WARNING " Found old config file: ${CONFIG_DIR}/calculation.rc; it will be removed and replaced with new one!")
#endif(EXISTS ${CONFIG_DIR}/calculation.rc)
#
##--------------------------------------------------

file(GLOB config_files RELATIVE ${PROJECT_SOURCE_DIR}/config_files config_files/*)

foreach(each_file ${config_files})
  set(SOUR_CONF_FIL "${PROJECT_SOURCE_DIR}/config_files/${each_file}")
#  set(DEST_CONF_FIL "${TMP_CONFIG_DIR}/${each_file}")
  file(COPY ${SOUR_CONF_FIL} DESTINATION ${TMP_CONFIG_DIR})
#  add_custom_command(TARGET luscus POST_BUILD
#                     COMMAND ${CMAKE_COMMAND} -E copy_if_different ${SOUR_CONF_FIL} ${DEST_CONF_FIL})
endforeach(each_file ${config_files})

add_subdirectory(plugins/gamess)
add_subdirectory(plugins/gvorb)
add_subdirectory(plugins/molden)
add_subdirectory(plugins/xyz)
add_subdirectory(plugins/diffdens)

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  find_program(HAS_BABEL NAME babel PATHS ${COMPILER_LOCATION})
  if(HAS_BABEL)
    message(status " found babel; installing babel plugins!")
    add_subdirectory(plugins_babel/cif)
    add_subdirectory(plugins_babel/g03)
    add_subdirectory(plugins_babel/g98)
    add_subdirectory(plugins_babel/mopac)
    add_subdirectory(plugins_babel/pdb)
  endif(HAS_BABEL)
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#add_custom_command(TARGET luscus POST_BUILD
#                   COMMAND ${CMAKE_COMMAND} -E copy_directory ${TMP_CONFIG_DIR} ${CONFIG_DIR})

install(TARGETS luscus  DESTINATION bin)
file(GLOB config_files ${TMP_CONFIG_DIR}/*)
install(FILES ${config_files} DESTINATION ${CONFIG_DIR})

#install config files
install(FILES ${PROJECT} )

