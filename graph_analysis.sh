#!/bin/bash
# set -x

dir=$PWD/../Results

# rm -rf $dir/*.log $dir/g*.json
rm -rf $dir
# mkdir -p $dir
vec_size="10000"

for players in 2
do
    sleep 1
    for rounds in 1
    do
        for party in $(seq 1 $players)
        do
            logdir=$dir/$players\_PC/$vec_size\_Nodes/TestRun/
            mkdir -p $logdir
            log=$logdir/round\_$rounds\_party\_$party.log
            tplog=$logdir/round\_$rounds\_party\_0.log
            if [ $party -eq 1 ]; then
                ./benchmarks/mpa_emgraph -p $party --localhost -l 100.0 -v $vec_size -i 10 -n $players 2>&1 | cat > $log &
            else
                ./benchmarks/mpa_emgraph -p $party --localhost -l 100.0 -v $vec_size -i 10 -n $players 2>&1 | cat > $log &
            fi
            codes[$party]=$!
        done

        ./benchmarks/mpa_emgraph -p 0 --localhost -l 100.0 -v $vec_size -i 10 -n $players 2>&1 | cat > $tplog &
        codes[0]=$!
        for party in $(seq 0 $players)
        do
            wait ${codes[$party]} || return 1
        done
    done
    python3 /code/pythonScripts/getAggStat.py $dir/$players\_PC/$vec_size\_Nodes/TestRun/
done
