#!/bin/bash

if [ $1 = "--version" ]; then
  echo babel_mopac_in version 1.0.0
  exit 0
fi

$HOME/.luscus/gv2xyz.exe $1
#this is for linux only
babel -i xyz ${1%.*}.xyz -o mopin ${1%.*}.dat
rm -f ${1%.*}.xyz
