#include <io/netmp.h>
#include <graphdb/offline_evaluator.h>
#include <graphdb/online_evaluator.h>
#include <utils/circuit.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <omp.h>

#include "utils.h"
#include "graphutils.h"

using namespace graphdb;
using json = nlohmann::json;
namespace bpo = boost::program_options;

common::utils::Circuit<Ring> generateCircuit(int nP, int pid, DistributedDaglist dist_daglist, DistributedDaglist new_vertices) {

    int nC = dist_daglist.num_clients;
    int nV = dist_daglist.nV;
    int nE = dist_daglist.nE;
    auto VSizes = dist_daglist.VSizes;
    auto ESizes = dist_daglist.ESizes;
    size_t vec_size = nV + nE;

    std::cout << "Generating circuit" << std::endl;
    
    common::utils::Circuit<Ring> circ;

    // Initialize all daglist field values
    std::vector<std::vector<wire_t>> vertex_src_values(nC);
    std::vector<std::vector<wire_t>> vertex_dst_values(nC);
    std::vector<std::vector<wire_t>> vertex_isV_values(nC);
    std::vector<std::vector<wire_t>> vertex_data_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigs_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigv_values(nC);
    std::vector<std::vector<wire_t>> vertex_sigd_values(nC);

    std::vector<std::vector<wire_t>> edge_src_values(nC);
    std::vector<std::vector<wire_t>> edge_dst_values(nC);
    std::vector<std::vector<wire_t>> edge_isV_values(nC);
    std::vector<std::vector<wire_t>> edge_data_values(nC);
    std::vector<std::vector<wire_t>> edge_sigs_values(nC);
    std::vector<std::vector<wire_t>> edge_sigv_values(nC);
    std::vector<std::vector<wire_t>> edge_sigd_values(nC);

    for (int i = 0; i < nC; ++i) {
        std::vector<wire_t> subg_vertex_src_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_dst_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_isV_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_data_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigs_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigv_values(VSizes[i]);
        std::vector<wire_t> subg_vertex_sigd_values(VSizes[i]);
        
        std::vector<wire_t> subg_edge_src_values(ESizes[i]);
        std::vector<wire_t> subg_edge_dst_values(ESizes[i]);
        std::vector<wire_t> subg_edge_isV_values(ESizes[i]);
        std::vector<wire_t> subg_edge_data_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigs_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigv_values(ESizes[i]);
        std::vector<wire_t> subg_edge_sigd_values(ESizes[i]);

        for (int j = 0; j < VSizes[i]; ++j){
            subg_vertex_src_values[j] = circ.newInputWire();
            subg_vertex_dst_values[j] = circ.newInputWire();
            subg_vertex_isV_values[j] = circ.newInputWire();
            subg_vertex_data_values[j] = circ.newInputWire();
            subg_vertex_sigs_values[j] = circ.newInputWire();
            subg_vertex_sigv_values[j] = circ.newInputWire();
            subg_vertex_sigd_values[j] = circ.newInputWire();
        }

        for (int j = 0; j < ESizes[i]; ++j){
            subg_edge_src_values[j] = circ.newInputWire();
            subg_edge_dst_values[j] = circ.newInputWire();
            subg_edge_isV_values[j] = circ.newInputWire();
            subg_edge_data_values[j] = circ.newInputWire();
            subg_edge_sigs_values[j] = circ.newInputWire();
            subg_edge_sigv_values[j] = circ.newInputWire();
            subg_edge_sigd_values[j] = circ.newInputWire();
        }

        vertex_src_values[i] = subg_vertex_src_values;
        vertex_dst_values[i] = subg_vertex_dst_values;
        vertex_isV_values[i] = subg_vertex_isV_values;
        vertex_data_values[i] = subg_vertex_data_values;
        vertex_sigs_values[i] = subg_vertex_sigs_values;
        vertex_sigv_values[i] = subg_vertex_sigv_values;
        vertex_sigd_values[i] = subg_vertex_sigd_values;

        edge_src_values[i] = subg_edge_src_values;
        edge_dst_values[i] = subg_edge_dst_values;
        edge_isV_values[i] = subg_edge_isV_values;
        edge_data_values[i] = subg_edge_data_values;
        edge_sigs_values[i] = subg_edge_sigs_values;
        edge_sigv_values[i] = subg_edge_sigv_values;
        edge_sigd_values[i] = subg_edge_sigd_values;
    
    }

    int add_nV = new_vertices.nV;
    auto add_VSizes = new_vertices.VSizes;
  
    // Initialize vertices to be added.
    std::vector<std::vector<wire_t>> new_vertex_src_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_dst_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_isV_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_data_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_sigs_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_sigv_values(nC);
    std::vector<std::vector<wire_t>> new_vertex_sigd_values(nC);

    for (int i = 0; i < nC; ++i) {
        std::vector<wire_t> subg_new_vertex_src_values(add_VSizes[i]);
        std::vector<wire_t> subg_new_vertex_dst_values(add_VSizes[i]);
        std::vector<wire_t> subg_new_vertex_isV_values(add_VSizes[i]);
        std::vector<wire_t> subg_new_vertex_data_values(add_VSizes[i]);
        
        for (int j = 0; j < add_VSizes[i]; ++j){
            subg_new_vertex_src_values[j] = circ.newInputWire();
            subg_new_vertex_dst_values[j] = circ.newInputWire();
            subg_new_vertex_isV_values[j] = circ.newInputWire();
            subg_new_vertex_data_values[j] = circ.newInputWire();
        }

        new_vertex_src_values[i] = subg_new_vertex_src_values;
        new_vertex_dst_values[i] = subg_new_vertex_dst_values;
        new_vertex_isV_values[i] = subg_new_vertex_isV_values;
        new_vertex_data_values[i] = subg_new_vertex_data_values;

    }

    // zero wire for assigning constant value to wire
    auto zero_wire = circ.newInputWire();

    // Compute per-client delta values (prefix sums of new vertices)
    std::vector<size_t> delta(nC);
    delta[0] = 0;
    for (int i = 1; i < nC; ++i) {
        delta[i] = delta[i-1] + add_VSizes[i-1];
    }

    // Assign positions for existing vertices
    
    std::vector<std::vector<wire_t>> updated_vertex_sigs_values(nC);
    std::vector<std::vector<wire_t>> updated_vertex_sigv_values(nC);
    std::vector<std::vector<wire_t>> updated_vertex_sigd_values(nC);

    for (int i = 0; i < nC; ++i) {
        updated_vertex_sigs_values[i].resize(VSizes[i]);
        updated_vertex_sigv_values[i].resize(VSizes[i]);
        updated_vertex_sigd_values[i].resize(VSizes[i]);

        for (int k = 0; k < VSizes[i]; ++k) {
            updated_vertex_sigs_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    vertex_sigs_values[i][k], delta[i]);
            updated_vertex_sigv_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    vertex_sigv_values[i][k], delta[i]);
            updated_vertex_sigd_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    vertex_sigd_values[i][k], delta[i]);
        }
    }
    
    // Assign positions to new vertices
    for (int i = 0; i < nC; ++i) {
        new_vertex_sigs_values[i].resize(add_VSizes[i]);
        new_vertex_sigv_values[i].resize(add_VSizes[i]);
        new_vertex_sigd_values[i].resize(add_VSizes[i]);

        for (int k = 0; k < add_VSizes[i]; ++k) {
            new_vertex_sigv_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    updated_vertex_sigv_values[i][VSizes[i]-1], Ring(k));
            new_vertex_sigd_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    updated_vertex_sigd_values[i][VSizes[i]-1], Ring(k));
        }   
    }

    for (int i = 0; i < nC - 1; ++i) {
        for (int k = 0; k < add_VSizes[i]; ++k) {
            new_vertex_sigs_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    updated_vertex_sigs_values[i+1][0], Ring(k*-1));
        } 
    }
    for (int k = 0; k < add_VSizes[nC - 1]; ++k) {
        new_vertex_sigs_values[nC - 1][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    zero_wire, Ring(vec_size + delta[nC - 1] + k));
    }

    // Assign positions of edges
    // Assign delta values to wires
    std::vector<wire_t> delta_wires(vec_size);
    int index = 0;
    for (size_t i = 0 ; i < nC; ++i) {
        for (size_t j = 0; j < VSizes[i]; ++j) {
            delta_wires[index] = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, delta[i]);
            index++;
        }
    }
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < ESizes[i]; ++j) {
            // Initialize wires to zero (don't set input)
            delta_wires[index] = circ.newInputWire();
            index++;
        }
    }

    // Flatten position maps
    std::vector<wire_t> sigv(vec_size);
    std::vector<wire_t> sigd(vec_size);
    index = 0;
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < VSizes[i]; ++j) {
            sigv[index] = vertex_sigv_values[i][j];
            sigd[index] = vertex_sigd_values[i][j];
            index++;
        }
    }
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < ESizes[i]; ++j) {
            sigv[index] = edge_sigv_values[i][j];
            sigd[index] = edge_sigd_values[i][j];
            index++;
        }
    }

    // Generate permutation for shuffle
    // Here we just pass identity permutations
    std::vector<int> base_perm(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        base_perm[i] = static_cast<int>(i);
    }
    std::vector<std::vector<int>> permutation;
    permutation.push_back(base_perm);
    if (pid == 0) {
        for (int i = 1; i < nP; ++i) {
            permutation.push_back(base_perm);
        }
    }

    // Propagate
    auto prop_delta = addSubCircPropagate(circ, sigd, delta_wires, nV, permutation, true);

    // Reorder back to vertex order
    auto delta_v = addSubCircPermList(circ, sigv, {prop_delta}, permutation)[0];

    // Update sigd; vector includes dummy values for first nV entries (corresponding to vertices)
    std::vector<wire_t> updated_edge_sigd_values(vec_size);
    for (size_t i = nV; i < vec_size; ++i) {
        updated_edge_sigd_values[i] = circ.addGate(common::utils::GateType::kAdd, sigd[i], delta_v[i]);
    }

    // Update sigv and sigs
    std::vector<std::vector<wire_t>> updated_edge_sigs_values(nC);
    std::vector<std::vector<wire_t>> updated_edge_sigv_values(nC);

    for (int i = 0; i < nC; ++i) {
        updated_edge_sigs_values[i].resize(ESizes[i]);
        updated_edge_sigv_values[i].resize(ESizes[i]);

        for (int k = 0; k < ESizes[i]; ++k) {
            updated_edge_sigs_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    edge_sigs_values[i][k], delta[i]);
            updated_edge_sigv_values[i][k] = 
                circ.addConstOpGate(common::utils::GateType::kConstAdd, 
                                    edge_sigv_values[i][k], delta[i]);
        }
    }

    // Set outputs
    index = 0;
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < VSizes[i]; ++j) {
            circ.setAsOutput(vertex_src_values[i][j]);
            circ.setAsOutput(vertex_dst_values[i][j]);
            circ.setAsOutput(vertex_isV_values[i][j]);
            circ.setAsOutput(vertex_data_values[i][j]);
            circ.setAsOutput(updated_vertex_sigs_values[i][j]);
            circ.setAsOutput(updated_vertex_sigv_values[i][j]);
            circ.setAsOutput(updated_vertex_sigd_values[i][j]);
            index++;
        }
        for (size_t j = 0; j < add_VSizes[i]; ++j) {
            circ.setAsOutput(new_vertex_src_values[i][j]);
            circ.setAsOutput(new_vertex_dst_values[i][j]);
            circ.setAsOutput(new_vertex_isV_values[i][j]);
            circ.setAsOutput(new_vertex_data_values[i][j]);
            circ.setAsOutput(new_vertex_sigs_values[i][j]);
            circ.setAsOutput(new_vertex_sigv_values[i][j]);
            circ.setAsOutput(new_vertex_sigd_values[i][j]);
        }
    }
    for (size_t i = 0; i < nC; ++i) {
        for (size_t j = 0; j < ESizes[i]; ++j) {
            circ.setAsOutput(edge_src_values[i][j]);
            circ.setAsOutput(edge_dst_values[i][j]);
            circ.setAsOutput(edge_isV_values[i][j]);
            circ.setAsOutput(edge_data_values[i][j]);
            circ.setAsOutput(updated_edge_sigs_values[i][j]);
            circ.setAsOutput(updated_edge_sigv_values[i][j]);
            circ.setAsOutput(updated_edge_sigd_values[index]);
            index++;
        }
    }

    return circ;
}