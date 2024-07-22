#!/bin/bash

# First line expands aliases (e.g. PATH, EDITOR, etc.). Second line loads in
# aliases/functions from the .bashrc file
shopt -s expand_aliases
source ~/.bash_profile

#-----------------------------------------------------------------------------
# Is this the correct executable to run?
#-----------------------------------------------------------------------------
code_stem=`basename "$0" .sh`
printf "%s\n\n%*s\n\n%s\n" \
       "Using code stem:" \
       $(( `echo $code_stem | wc -c` + 20 | bc )) "${RED}$code_stem${NORMAL}" \
       "Is this correct? (y/n)"
read response
if [ "$response" != "y" ]; then echo "Exiting..."; exit; fi

# Build the code
#---------------
printf "\n%s\n\n%s\n\n" "Building code in:" "    ${RED}`pwd`${NORMAL}"
make $code_stem
#-----------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Parameter list
#------------------------------------------------------------------------------
output_flags_and_exit=1;
# t1=(0.9, 0.8, 0.7, 0.6, 0.5);
# t1=(0.88);
# t2=0.63
# s=0.7
t1=(0.98);
t2=0.53
s=0.15
counter=2

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Run it over whatever set of parameters we may want
#------------------------------------------------------------------------------
# Loop over the type of solver
for t in ${t1[@]}
do

		dir="./mesh$counter"
		mkdir -p "$dir/data"

    # Command line flags
    #-------------------
    # NOTE: ALWAYS have a space before the end of the cmd_line_flag string
    cmd_line_flags="";
    cmd_line_flags+="--t1 $t "
		cmd_line_flags+="--t2 $t2 "
    cmd_line_flags+="--s $s "
    cmd_line_flags+="--counter $counter "
		cmd_line_flags+="--dir $dir "


    ((counter++))
    # Split by the delimiter "--"
    #----------------------------
    printf "\n%s\n" "Using command-line flags:"
    if [ $output_flags_and_exit == 1 ]
    then
	printf "\n%s\n\n" "       $RED$cmd_line_flags$NORMAL"
    else
	temp1=`echo $cmd_line_flags | sed 's/--/    --/'`
	echo -ne "$RED"
	echo "$temp1" | sed 's/--/\'$'\n    --/g'
	echo -e "$NORMAL"
    fi

    # Start the timers
    #-----------------
    start=`date +%s`;
    start_time=`date +%T`

    # Output some stuff
    #------------------
    echo " "
    echo "============================================="
    echo "Start time: $start_time"
    echo "============================================="
    echo " "

    # Run the code
    #-------------
    # Run the code as usual
    ./$code_stem $cmd_line_flags

    # End the timer stuff
    #--------------------
    end=`date +%s`
    end_time=`date +%T`

    # Get the runtime
    #----------------
    run_time=$( echo "$end - $start" | bc -l )

    # How long did it take?
    #----------------------
    echo "============================================="
    echo "Start time: $start_time"
    echo "End time: $end_time"
    echo "Total runtime [sec]: $run_time"
    echo "============================================="
done
#------------------------------------------------------------------------------
