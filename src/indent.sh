#!/bin/bash

#
#        FILE: indent.sh
#      AUTHOR: Paul Bartholomew
# DESCRIPTION: Use emacs to indent files
#

for f in ./*.f90
do
    echo "Indenting ${f}"
	emacs -batch ${f} --eval '(indent-region (point-min) (point-max) nil)' -f save-buffer 2> /dev/null
done

rm *.f90~
exit
