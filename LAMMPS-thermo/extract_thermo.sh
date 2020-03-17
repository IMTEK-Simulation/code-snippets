#! /bin/sh
# Extracts thermo output from lammps log.
# Requires a SINGLE thermo section within the log.
# Does NOT extract multiple thermo sections, i.e.
# subsequent minimization and equilibration.
# 1st and 2nd positional arguments are optional and
# allow to specify input log file and output text file.
LOGFILE="lammps.log"
OUTFILE="thermo.out"
if [ -n "${1}" ] ; then
  LOGFILE="${1}"
  if [ -n "${2}" ] ; then
    OUTFILE="${2}"
  fi
fi
# first sed cuts out thermo section from log
# second sed removes all lines added by colvars
# other outputs (i.e. fix shake stats)  can be filtered in a similar manner if necessary
# head removes last line (might be incomplete)
if [ -f "${LOGFILE}" ]; then
  cat "${LOGFILE}" | sed -n '/^Step/,/^Loop time/p' | sed '/^colvars:/d' | head -n-1 > "${OUTFILE}"
else
  echo "${LOGFILE} does not exist!" >&2
  exit 1
fi
