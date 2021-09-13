#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=4

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation

# Validation for dense complex matrix
#-----------------------------------------
echo "Running dense complex matrix alidation "
../dense_complex_matrix_test > OUTPUT
echo "done"
echo " " >> validation.log
echo "Dense complex matrix validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/dense_complex_matrix_test.dat.gz  \
         OUTPUT >> validation.log
fi

mv OUTPUT OUTPUT_dense_complex_matrix_test
#-----------------------------------------

# Validation for compressed row complex matrix
#-----------------------------------------
echo "Running compressed row complex matrix validation "
../cr_complex_matrix_test > OUTPUT
echo "done"
echo " " >> validation.log
echo "Compressed row complex matrix validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/cr_complex_matrix_test.dat.gz  \
         OUTPUT >> validation.log
fi

mv OUTPUT OUTPUT_cr_complex_matrix_test
#-----------------------------------------

# Validation for compressed column complex matrix
#-----------------------------------------
echo "Running compressed column complex matrix validation "
../cc_complex_matrix_test > OUTPUT
echo "done"
echo " " >> validation.log
echo "Compressed column complex matrix validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/cc_complex_matrix_test.dat.gz  \
         OUTPUT >> validation.log
fi

mv OUTPUT OUTPUT_cc_complex_matrix_test
#-----------------------------------------

# Validation for complex eigensolver
#-----------------------------------------
echo "Running complex eigensolver validation "
../eigensolver_test > OUTPUT
echo "done"
echo " " >> validation.log
echo "Complex eigensolver validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/eigensolver_test.dat.gz  \
         OUTPUT >> validation.log
fi

mv OUTPUT OUTPUT_eigensolver_test
#-----------------------------------------

# Append log to main validation log
cat validation.log >> ../../../../validation.log

cd ..

#######################################################################

#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10