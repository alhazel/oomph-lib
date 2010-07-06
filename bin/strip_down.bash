#! /bin/bash

#-------------------------------------------------------
# Script to strip down oomph-lib disribution to bare
# essentials to allow test-driving of installation
# process.
#-------------------------------------------------------
home_dir=`pwd`


#====================================================================
# A few helper functions
#====================================================================
# A little function 'borrowed' from the tecplot installation script...
OptionPrompt() 
{ 
 printf "%s " "$1" 
}

# Another little function 'borrowed' from the tecplot installation script...
OptionRead()
{
 read Opt
 if test "$Opt" = "" ; then
  Opt=$1
 fi
 echo $Opt
}
echo " " 
echo "======================================================== " 
echo " " 
echo "WARNING: This script strips out most of the distribution "
echo "         and retains only the bare minimum (the build"
echo "         machinery, src/generic and src/poisson) required"
echo "         to test-drive the installation procedure."
echo "         Files that are not required will be DELETED!"
echo " "
echo "======================================================== " 
echo " " 
echo " "
echo "We are currently in: " $home_dir
echo " " 
echo " "
echo " Do you want to strip down this distribution [y/n -- default: n]"
reply=`OptionRead`
if test "$reply" != "y" -a "$reply" != "Y" ; then 
    echo "Not doing it..."
    exit 0
else
    echo "Doing it..."
fi


# Backup directory with overwritten config files
if [ -e backed_up_stripped_out_files ]
then
    echo " "
    echo "ERROR: Directory "
    echo " " 
    echo "       backed_up_stripped_out_files"
    echo " "
    echo "already exists. Please delete it and try again."
    echo " "
    exit 1
else
    mkdir backed_up_stripped_out_files
fi


# Fake delete/copy commands
del="echo Deleting:  `ls -d`"
cop="echo I would overwrite these files: `ls -d`"

# Real delete/copy commands
del="rm -rf "
cop="cp "

# Generic top level stuff
#========================

# Update makefiles
$cop Makefile.am backed_up_stripped_out_files/Makefile.am.top_level 
$cop $home_dir/config/stripped_down_files/Makefile.am.top_level Makefile.am
$cop Makefile.am backed_up_stripped_out_files/Makfile.am.bin
$cop $home_dir/config/stripped_down_files/Makefile.am.bin bin/Makefile.am


# Get rid of user src and drivers
echo "Deleting user_src and user-drivers"
$del user_src
$del user_drivers

# Kill private directory
echo "Deleting private"
$del private

# Get rid of entire doc directory
echo "Deleting doc"
$del doc



# Self tests: Only keep analyse stage
#====================================
echo "Wiping most of self_test"
cd self_test
mv analyse_self_tests ../tmp_analyse_self_tests
$del *
mv ../tmp_analyse_self_tests analyse_self_test
$cop Makefile.am backed_up_stripped_out_files/Makfile.am.self_test
$cop $home_dir/config/stripped_down_files/Makefile.am.self_test Makefile.am



# Demo drivers
#=============
cd $home_dir
echo "Wiping most of demo drivers"
cd demo_drivers


# Kill everything apart from Poisson and mpi
#-------------------------------------------
mv poisson ../tmp_poisson
mv mpi ../tmp_mpi
$del *
mv ../tmp_poisson poisson
mv ../tmp_mpi mpi
$cop Makefile.am backed_up_stripped_out_files/Makfile.am.demo_drivers
$cop $home_dir/config/stripped_down_files/Makefile.am.demo_drivers Makefile.am


# Within Poisson get rid of everything apart from adaptive flux bc
#-----------------------------------------------------------------
cd poisson
mv two_d_poisson_flux_bc_adapt ../tmp_two_d_poisson_flux_bc_adapt
$del *
mv ../tmp_two_d_poisson_flux_bc_adapt two_d_poisson_flux_bc_adapt 
$cop Makefile.am backed_up_stripped_out_files/Makfile.am.demo_drivers.poisson
$cop $home_dir/config/stripped_down_files/Makefile.am.demo_drivers.poisson Makefile.am
cd ..


# Within mpi get rid of everything apart from adaptive flux bc
#-------------------------------------------------------------
cd mpi
mv distribution ../tmp_distribution
$del *
mv ../tmp_distribution distribution
$cop Makefile.am backed_up_stripped_out_files/Makfile.am.demo_drivers.mpi
$cop $home_dir/config/stripped_down_files/Makefile.am.demo_drivers.mpi Makefile.am

cd distribution
mv two_d_poisson_flux_bc_adapt ../tmp_two_d_poisson_flux_bc_adapt
$del * 
mv ../tmp_two_d_poisson_flux_bc_adapt two_d_poisson_flux_bc_adapt
$cop Makefile.am backed_up_stripped_out_files/Makfile.am.demo_drivers.mpi.distribution
$cop $home_dir/config/stripped_down_files/Makefile.am.demo_drivers.mpi.distribution Makefile.am




# Sources
#========
cd $home_dir/src
echo "Wiping most of src"

# Keep only Poisson, generic and meshes
#--------------------------------------
mv generic ../tmp_generic
mv poisson ../tmp_poisson
mv meshes ../tmp_meshes
$del *
mv ../tmp_generic generic
mv ../tmp_poisson poisson
mv ../tmp_meshes  meshes
$cop Makefile.am backed_up_stripped_out_files/Makfile.am.src
$cop $home_dir/config/stripped_down_files/Makefile.am.src Makefile.am


# Keep only rectangular quad meshes
#----------------------------------
cd meshes
mkdir ../tmp_meshes
cp rectangular_quadmesh.template.h simple_rectangular_quadmesh.template.h ../tmp_meshes
cp rectangular_quadmesh.template.cc simple_rectangular_quadmesh.template.cc ../tmp_meshes
$del *.h *.cc
mv ../tmp_meshes/* .
rmdir ../tmp_meshes
$cop mesh_names.list backed_up_stripped_out_files/mesh_names.list
$cop $home_dir/config/stripped_down_files/mesh_names.list .
$cop Makefile.am backed_up_stripped_out_files/Makefile.am.meshes
$cop $home_dir/config/stripped_down_files/Makefile.am.meshes Makefile.am





# Update configure scripts
#=========================
echo "Updating configure scripts"

cd $home_dir
cd config/configure.ac_scripts


$cop private.dir_list backed_up_stripped_out_files/private.dir_list 
$cop $home_dir/config/stripped_down_files/private.dir_list .

$cop private_user_drivers.dir_list backed_up_stripped_out_files/private_user_drivers.dir_list 
$cop $home_dir/config/stripped_down_files/private_user_drivers.dir_list .

$cop private_user_src.dir_list backed_up_stripped_out_files/private_user_src.dir_list 
$cop $home_dir/config/stripped_down_files/private_user_src.dir_list .

$cop user_drivers.dir_list backed_up_stripped_out_files/user_drivers.dir_list 
$cop $home_dir/config/stripped_down_files/user_drivers.dir_list .

$cop user_src.dir_list backed_up_stripped_out_files/user_src.dir_list 
$cop $home_dir/config/stripped_down_files/user_src.dir_list .

$cop core.dir_list backed_up_stripped_out_files/core.dir_list 
$cop $home_dir/config/stripped_down_files/core.dir_list .

$cop doc.dir_list backed_up_stripped_out_files/doc.dir_list 
$cop $home_dir/config/stripped_down_files/doc.dir_list .



