#!/bin/bash

## We must be under sudo
#if [ `id -u` -ne 0 ]; 
#then
#        echo "You must execute this script as root!"
#        exit 1
#fi

METASEQR_HOME=/media/HD4/Fleming/dev/metaseqR
METASEQR_BUILD=/media/HD4/Fleming/dev/metaseqR-build
REPO_HOME=/var/www/Rrepos

#rm $METASEQR_HOME/inst/doc/metaseqr-pdf.tex
rm $METASEQR_HOME/inst/doc/metaseqr-html.html
#rm $METASEQR_HOME/inst/doc/metaseqr-html.md
#rm $METASEQR_HOME/inst/doc/metaseqr-pdf.log
#rm $METASEQR_HOME/inst/doc/metaseqr-pdf.out
rm $METASEQR_HOME/inst/doc/metaseqr-pdf.pdf
#rm $METASEQR_HOME/inst/doc/metaseqr-pdf.aux
rm $METASEQR_HOME/inst/doc/.build.timestamp
#rm -r $METASEQR_HOME/man
rm -r $REPO_HOME/*

Rscript -e "remove.packages('metaseqR')"
Rscript -e "source('roxybuild.R')"

if [ -d $METASEQR_BUILD ]
then
	break
else
	mkdir -p $METASEQR_BUILD
fi

rsync -r --exclude=.svn $METASEQR_HOME $METASEQR_BUILD

#R CMD build --resave-data $METASEQR_BUILD/metaseqR
#R CMD check --no-build-vignettes $METASEQR_BUILD/metaseqR
