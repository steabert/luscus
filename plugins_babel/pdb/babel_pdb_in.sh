#!/bin/bash

if [ $1 = "--version" ]; then
  echo babel_pdb_in version 1.0.0
  exit 0
fi

$HOME/.luscus/gv2xyz.exe $1
#this is for linux only
babel -i xyz ${1%.*}.xyz -o pdb ${1%.*}.pdb
rm -f ${1%.*}.xyz
