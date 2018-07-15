##----INSTALL(PROGRAMS scr1 ${CMAKE_CURRENT_BINARY_DIR}/scr2  DESTINATION bin)

INSTALL(PROGRAMS babel_cif_in.sh  babel_cif_out.sh DESTINATION bin)
file(APPEND ${TMP_CONFIG_DIR}/plugin.rc "  libpath=${CMAKE_INSTALL_PREFIX}/bin extension=cif description=\"Crystallographic Information File\"  forward=babel_cif_out.sh  backward=babel_cif_in.sh \n")

#file(COPY babel_cif_in.sh  babel_cif_out.sh DESTINATION ${TMP_CONFIG_DIR})
#file(APPEND ${TMP_CONFIG_DIR}/plugin.rc "  libpath=${TMP_CONFIG_DIR} extension=cif description=\"Crystallographic Information File\"  forward=babel_cif_out.sh  backward=babel_cif_in.sh \n")
