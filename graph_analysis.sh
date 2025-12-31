#!/bin/bash
set -e

# graph_analysis.sh
# Runs a set of benchmarks (multiple graph sizes) using run.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SH="$SCRIPT_DIR/run.sh"

if [ ! -x "$RUN_SH" ]; then
    echo "Error: run.sh not found or not executable at $RUN_SH"
    exit 1
fi

SIZES=("1000:10000" "10000:100000" "10000:1000000" "100000:10000000")

run_benchmark() {
    local bench="$1"
    shift
    for s in "${SIZES[@]}"; do
        verts=${s%%:*}
        edges=${s##*:}
        echo "Running $bench with --num-vert $verts --num-edge $edges"
        if [ "$bench" = "add_vertices" ] || [ "$bench" = "add_edges" ]; then
            $RUN_SH "$bench" --num-vert "$verts" --num-edge "$edges" --num-inserts 1000 --latency 0
        else
            $RUN_SH "$bench" --num-vert "$verts" --num-edge "$edges" --latency 0
        fi
        echo "Finished $bench $verts/$edges"
        echo "----------------------------------------"
    done
}

echo "Starting graph analysis benchmarks"

run_benchmark graphiti_init
run_benchmark grapharo_init

# Run grapharo_init with varying client counts
echo "Running grapharo_init with varying client counts..."
for clients in 2 5 10 15; do
    echo "Running grapharo_init with --num-clients $clients --num-vert 100000 --num-edge 1000000"
    $RUN_SH grapharo_init --num-vert 100000 --num-edge 1000000 --num-clients "$clients" --latency 0
    echo "Finished grapharo_init with $clients clients"
    echo "----------------------------------------"
done

run_benchmark add_vertices
run_benchmark add_edges
run_benchmark delete
run_benchmark data_change

echo "All benchmarks queued/completed. Check Results/ for logs."
