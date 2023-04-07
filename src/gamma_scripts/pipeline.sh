#!/bin/bash
# Runs with $ sh quickrun.sh # # # ...
#
# Compile code.
#
# make
#
# Execute tests.
#
# echo -e "version\tn\tmax\tl2norm\ttime" 

GENS=20000
BATCH=25

#mkdir results_"$BATCH"
#./make_db.py $GENS $BATCH
#./view_db.py $BATCH
#mv best"$BATCH".h5 results_"$BATCH"/
./draw_full.py results_"$BATCH"/best"$BATCH".h5 $BATCH
