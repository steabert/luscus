INSTALL(PROGRAMS babel_mopac_in.sh  babel_mopac_out.sh DESTINATION bin)
#file(COPY babel_mopac_in.sh  babel_mopac_out.sh DESTINATION bin)
file(APPEND ${TMP_CONFIG_DIR}/plugin.rc "  libpath=${CMAKE_INSTALL_PREFIX}/bin  extension=out  description=mopac  forward=babel_mopac_out.sh  backward=babel_mopac_in.sh  \n")
