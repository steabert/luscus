#!/bin/bash

if [ $1 = "--version" ]; then
  echo babel_g03_in version 1.0.0
  exit 0
fi

$HOME/.luscus/gv2xyz.exe $1
#this is for linux only
babel -i xyz ${1%.*}.xyz -o gau ${1%.*}.com
rm -f ${1%.*}.xyz
