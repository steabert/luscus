#!/bin/bash

if [ $1 = "--version" ]; then
  echo babel_cif_in version 1.0.0
  exit 0
fi

$HOME/.luscus/gv2xyz.exe $1
#this is for linux only
babel -i xyz ${1%.*}.xyz -o cif ${1%.*}.cif
rm -f ${1%.*}.xyz
