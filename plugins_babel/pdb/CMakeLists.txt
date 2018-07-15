INSTALL(PROGRAMS  babel_pdb_in.sh  babel_pdb_out.sh  DESTINATION bin)
#file(COPY babel_pdb_in.sh  babel_pdb_out.sh  DESTINATION bin)
file(APPEND ${TMP_CONFIG_DIR}/plugin.rc "  libpath=${CMAKE_INSTALL_PREFIX}/bin  extension=pdb  description=\"Protein Data Bank\"  forward=babel_pdb_out.sh  backward=babel_pdb_in.sh  \n")
