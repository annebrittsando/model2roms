#!/bin/bash
#
#  Give the job a name
#PBS -N "model2roms_KINO"
#
#  Specify the project the job belongs to
#PBS -A nn9253k
#
#PBS -q normal
#
#PBS -l mppwidth=1
#PBS -l mppmem=2000MB
#PBS -l mppnppn=16
#
# Request 1 hour of walltime
#PBS -l walltime=196:00:00
#
#  Send me an email on  a=abort, b=begin, e=end
#PBS -m abe
#
#  Use this email address (check that it is correct):
#PBS -M trond.kristiansen@imr.no
#
#  Write the standard output of the job to file 'mpijob.out' (optional)
#PBS -o  model2roms_KINO.out
#
#  Write the standard error of the job to file 'mpijob.err' (optional)
#PBS -e  model2roms_KINO.err
#

#  Make sure I am in the correct directory
cd /work/users/trondk/KINO/model2roms
export MPLCONFIGDIR=${pwd}

export TMP=`pwd`
module unload notur
aprun -B python main.py > model2roms_output_24032015.log
