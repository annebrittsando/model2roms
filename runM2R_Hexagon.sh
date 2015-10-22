#!/bin/bash                                                                                                                                                                                               
#                                                                                                                                                                                                         
#  Give the job a name                                                                                                                                                                                    
#PBS -N "model2roms_NS8KM"                                                                                                                                                                                 
#                                                                                                                                                                                                         
#  Specify the project the job belongs to                                                                                                                                                                 
#PBS -A imr                                                                                                                                                                                               
#PBS -q normal                                                                                                                                                                                            
#PBS -l mppwidth=1,walltime=1:00:00                                                                                                                                                                      
#PBS -l mppmem=1000MB                                                                                                                                                                                     
echo "#PBS -l mppnppn=16"
#                                                                                                                                                                                                         
#  Send me an email on  a=abort, b=begin, e=end                                                                                                                                                           
#PBS -m abe                                                                                                                                                                                               
#                                                                                                                                                                                                         
#  Use this email address (check that it is correct):                                                                                                                                                     
#PBS -M trond.kristiansen@imr.no                                                                                                                                                                          
#                                                                                                                                                                                                         
#  Write the standard output of the job to file 'mpijob.out' (optional)                                                                                                                                   
#PBS -o  model2roms_NS8KM.out                                                                                                                                                                              
#                                                                                                                                                                                                         
#  Write the standard error of the job to file 'mpijob.err' (optional)                                                                                                                                    
#PBS -e  model2roms_NS8KM.err                                                                                                                                                                              
#                                                                                                                                                                                                         

#  Make sure I am in the correct directory                                                                                                                                                                
cd /work/users/trondk/NS8km/model2roms
module load python

export MPLCONFIGDIR=${pwd}
export TMP=`pwd`
export PYTHON_EGG_CACHE=/work/users/trondk/NS8km/model2roms

aprun -B python main.py > model2roms_output_19102015.log