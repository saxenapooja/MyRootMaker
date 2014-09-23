#!/bin/bash

for x in work.*/jobs/job_*.txt; do
    if grep SUCCESS $x > /dev/null; then
	NUM=`echo $x | sed 's/.*job_\([0-9]\+\)\.txt/\1/'`
	BLA=`dirname $x`
	rm -f `dirname $BLA`/output/job_$NUM/cmssw.log.gz
    fi
done
