#!/bin/bash

if [ $1 = "--version" ]; then
  echo babel_g03_out version 1.0.0
  exit 0
fi

babel -i g03 $1 -o xyz ${1%.*}.xyz
#works for linux only
$HOME/.luscus/xyz2gv.exe ${1%.*}.xyz
rm -f ${1%.*}.xyz
