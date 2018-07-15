#!/bin/bash

#echo running density diff!
#echo 'arg0 =' $1
#echo 'diff weight =' $2
#echo 'file 2 =' $3
#echo 'file 1 =' $4
if [ -n "$LUSCUS_DIR" ]; then
  RUNDIR=$LUSCUS_DIR
elif [ -d $HOME/.luscus ]; then
  RUNDIR="$HOME/.luscus"
else
  RUNDIR="."
fi
$RUNDIR/diffdens.exe $4 $3 $2
sleep 1
