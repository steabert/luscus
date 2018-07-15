@echo off

if defined LUSCUS_DIR (
  set RUNDIR=%LUSCUS_DIR%\luscus
) else (
  if defined LOCALAPPDATA (
    set RUNDIR=%LOCALAPPDATA%\luscus
  ) else (
    set RUNDIR=.
  )
)

echo rundir = %RUNDIR%

%RUNDIR%\diffdens.exe %4 %3 %2

