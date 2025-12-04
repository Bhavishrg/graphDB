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

common::utils::Circuit<Ring> generateCircuit(int nP, int pid, DistributedDaglist dist_daglist, DistributedDaglist new_edges) {

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

    int add_nE = new_edges.nE;
    auto add_ESizes = new_edges.ESizes;

    // Initialize edges to be added.
    std::vector<std::vector<wire_t>> new_edge_src_values(nC);
    std::vector<std::vector<wire_t>> new_edge_dst_values(nC);
    std::vector<std::vector<wire_t>> new_edge_isV_values(nC);
    std::vector<std::vector<wire_t>> new_edge_data_values(nC);
    std::vector<std::vector<wire_t>> new_edge_sigs_values(nC);
    std::vector<std::vector<wire_t>> new_edge_sigv_values(nC);
    std::vector<std::vector<wire_t>> new_edge_sigd_values(nC);
    std::vector<std::vector<wire_t>> new_edge_overall_sigs_values(nC);
    std::vector<std::vector<wire_t>> new_edge_overall_sigd_values(nC);

    for (int i = 0; i < nC; ++i) {
        std::vector<wire_t> subg_new_edge_src_values(add_ESizes[i]);
        std::vector<wire_t> subg_new_edge_dst_values(add_ESizes[i]);
        std::vector<wire_t> subg_new_edge_isV_values(add_ESizes[i]);
        std::vector<wire_t> subg_new_edge_data_values(add_ESizes[i]);
        // Here sigs, sigd are for within each client's list of added edges
        std::vector<wire_t> subg_new_edge_sigs_values(add_ESizes[i]);
        std::vector<wire_t> subg_new_edge_sigv_values(add_ESizes[i]);
        std::vector<wire_t> subg_new_edge_sigd_values(add_ESizes[i]);
        
        for (int j = 0; j < add_ESizes[i]; ++j){
            subg_new_edge_src_values[j] = circ.newInputWire();
            subg_new_edge_dst_values[j] = circ.newInputWire();
            subg_new_edge_isV_values[j] = circ.newInputWire();
            subg_new_edge_data_values[j] = circ.newInputWire();
            subg_new_edge_sigs_values[j] = circ.newInputWire();
            subg_new_edge_sigv_values[j] = circ.newInputWire();
            subg_new_edge_sigd_values[j] = circ.newInputWire();
        }

        new_edge_src_values[i] = subg_new_edge_src_values;
        new_edge_dst_values[i] = subg_new_edge_dst_values;
        new_edge_isV_values[i] = subg_new_edge_isV_values;
        new_edge_data_values[i] = subg_new_edge_data_values;
        new_edge_sigs_values[i] = subg_new_edge_sigs_values;
        new_edge_sigv_values[i] = subg_new_edge_sigv_values;
        new_edge_sigd_values[i] = subg_new_edge_sigd_values;

    }

    // Create V^Out_i and V^In_i inputs per party (both length |V|)
    std::vector<std::vector<wire_t>> Vout(nC);
    std::vector<std::vector<wire_t>> Vin(nC);
    for (int i = 0; i < nC; ++i) {
        std::vector<wire_t> Vout_party(nV);
        std::vector<wire_t> Vin_party(nV);
        for (size_t j = 0; j < nV; ++j) {
            Vout_party[j] = circ.newInputWire();
            Vin_party[j] = circ.newInputWire();
        }
        Vout[i] = Vout_party;
        Vin[i] = Vin_party;
    }

    // zero wire - input 0
    auto zero_wire = circ.newInputWire();

    // Compute aggregated V_in, V_out
    std::vector<wire_t> Vout_agg(nV);
    std::vector<wire_t> Vin_agg(nV);
    for (size_t j = 0; j < nV; ++j) {
        // start with party 0's value then add others
        wire_t acc_out = Vout[0][j];
        for (int p = 1; p < nC; ++p) {
            acc_out = circ.addGate(common::utils::GateType::kAdd, acc_out, Vout[p][j]);
        }
        Vout_agg[j] = acc_out;

        wire_t acc_in = Vin[0][j];
        for (int p = 1; p < nC; ++p) {
            acc_in = circ.addGate(common::utils::GateType::kAdd, acc_in, Vin[p][j]);
        }
        Vin_agg[j] = acc_in;
    }

    // Compute cumulative offsets OffIn[j] and OffOut[j] = sum_{k=0}^{j-1} Vout_agg[k]
    std::vector<wire_t> OffIn(nV);
    std::vector<wire_t> OffOut(nV);
    OffIn[0] = zero_wire;
    OffOut[0] = zero_wire;
    for (size_t j = 1; j < nV; ++j) {
        OffIn[j] = circ.addGate(common::utils::GateType::kAdd, OffIn[j - 1], Vin_agg[j - 1]);
        OffOut[j] = circ.addGate(common::utils::GateType::kAdd, OffOut[j - 1], Vout_agg[j - 1]);
    }
    
    // For each party i, compute indicator arrays
    std::vector<std::vector<wire_t>> IIn(nC);
    std::vector<std::vector<wire_t>> IOut(nC);
    for (int i = 0; i < nC; ++i) {
        IIn[i].resize(nV);
        IOut[i].resize(nV);
        for (size_t j = 0; j < nV; ++j) {
            // IIn: 1 - Eqz(Vin[i][j])  --> neg_eq0 = -Eqz(Vin); IIn = neg_eq0 + 1
            auto eq0 = circ.addGate(common::utils::GateType::kEqz, Vin[i][j]);
            auto neg_eq0 = circ.addConstOpGate(common::utils::GateType::kConstMul, eq0, Ring(-1));
            IIn[i][j] = circ.addConstOpGate(common::utils::GateType::kConstAdd, neg_eq0, Ring(1));

            // IOut: 1 - Eqz(Vout[i][j] - 1)
            auto v_minus_1 = circ.addConstOpGate(common::utils::GateType::kConstAdd, Vout[i][j], Ring(-1));
            auto eq1 = circ.addGate(common::utils::GateType::kEqz, v_minus_1);
            auto neg_eq1 = circ.addConstOpGate(common::utils::GateType::kConstMul, eq1, Ring(-1));
            IOut[i][j] = circ.addConstOpGate(common::utils::GateType::kConstAdd, neg_eq1, Ring(1));
        }
    }
    
    // Precompute absolute vertex indices for each party's existing vertices
    std::vector<size_t> party_vertex_offset(nC);
    party_vertex_offset[0] = 0;
    for (int i = 1; i < nC; ++i) {
        party_vertex_offset[i] = party_vertex_offset[i-1] + VSizes[i-1];
    }

    // For each vertex, compute data_e = sigs + OffOut[j] 
    std::vector<wire_t> data_e(nV);
    for (int i = 0; i < nC; ++i) {
        for (size_t j = 0; j < VSizes[i]; ++j) {
            size_t abs_vertex_idx = party_vertex_offset[i] + j;
            data_e[abs_vertex_idx] = circ.addGate(common::utils::GateType::kAdd, 
                vertex_sigs_values[i][j], OffOut[abs_vertex_idx]);   
        }
    }

    // Compute new edge source group keys
    std::vector<std::vector<wire_t>> new_edge_src_group(nC);
    for (int i = 0; i < nC; ++i) {
        if (add_ESizes[i] == 0) {
            new_edge_src_group[i] = {};
            continue;
        }
        std::vector<wire_t> new_edge_src_group_party(add_ESizes[i]);
        new_edge_src_group_party[0] = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, 1);
        for (size_t j = 1; j < add_ESizes[i]; ++j) {
            auto temp = 
                circ.addGate(common::utils::GateType::kSub, new_edge_src_values[i][j-1], new_edge_src_values[i][j]);
            temp = circ.addGate(common::utils::GateType::kEqz, temp);
            temp = circ.addConstOpGate(common::utils::GateType::kConstMul, temp, Ring(-1));
            new_edge_src_group_party[j] = circ.addConstOpGate(common::utils::GateType::kConstAdd, temp, Ring(1));
        }
        new_edge_src_group[i] = std::move(new_edge_src_group_party);
    }

    // Prepare permutations (two separate permutations: one for T1, one for T2)
    // Permutations set as identity for benchmarking
    std::vector<std::vector<int>> permutation;
    
    // T1 permutation
    std::vector<int> t1_perm(nV);
    for (size_t i = 0; i < nV; ++i) {
        t1_perm[i] = i;
    }
    
    // T2 permutation
    std::vector<int> t2_perm(add_nE);
    for (size_t i = 0; i < add_nE; ++i) {
        t2_perm[i] = i;
    }
    
    // For party 0, we need all parties' permutations
    // For other parties, just their own permutation
    if (pid == 0) {
        // Party 0 stores all parties' permutations
        // For simplicity, using same permutation for all parties
        permutation.push_back(t1_perm);  // T1 permutation
        permutation.push_back(t2_perm);  // T2 permutation
    } else {
        // Other parties store their own permutations
        permutation.push_back(t1_perm);  // T1 permutation
        permutation.push_back(t2_perm);  // T2 permutation
    }

    // Group-wise propagate and index gates
    // Compute sigs in overall graph for new edges
    for (int i = 0; i < nC; ++i) {
        if (add_ESizes[i] == 0) {
            new_edge_overall_sigs_values[i] = {};
            continue;
        }
        std::vector<wire_t> subg_new_edge_sigs_values(add_ESizes[i]);
        auto [prop_out_key, prop_out_v] = 
            circ.addGroupwisePropagateGate(IOut[i], data_e, new_edge_src_group[i], permutation);
        auto [out_ind, ind_output_key, ind_output_v] = 
            circ.addGroupwiseIndexGate(new_edge_src_group[i], prop_out_v, permutation);
        for (int j = 0; j < add_ESizes[i]; ++j) {
            subg_new_edge_sigs_values[j] = 
                circ.addGate(common::utils::GateType::kAdd, prop_out_v[j], out_ind[j]);
        }
        new_edge_overall_sigs_values[i] = subg_new_edge_sigs_values;
    }

    // Reorder new edges to destination order
    // Combine components of edge lists into payloads; only need dst and sigs
    std::vector<std::vector<std::vector<wire_t>>> payloads(nC);
    for (int i = 0; i < nC; ++i) {
        payloads[i].reserve(2);
        payloads[i].push_back(new_edge_dst_values[i]);
        payloads[i].push_back(new_edge_sigs_values[i]);
    }

    // PermList gates
    std::vector<std::vector<std::vector<wire_t>>> payloads_out(nC);
    for (int i = 0; i < nC; ++i) {
        if (add_ESizes[i] == 0) {
            payloads_out[i] = {};
        } else {
            payloads_out[i] = addSubCircPermList(circ, new_edge_sigd_values[i], payloads[i], permutation);
        }
    }
    
    // Compute new edge dest group keys
    std::vector<std::vector<wire_t>> new_edge_dest_group(nC);
    for (int i = 0; i < nC; ++i) {
        if (add_ESizes[i] == 0) {
            new_edge_dest_group[i] = {};
            continue;
        }
        std::vector<wire_t> new_edge_dest_group_party(add_ESizes[i]);
        new_edge_dest_group_party[0] = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, 1);
        for (size_t j = 1; j < add_ESizes[i]; ++j) {
            auto temp = 
                circ.addGate(common::utils::GateType::kSub, payloads_out[i][0][j-1], payloads_out[i][0][j]);
            temp = circ.addGate(common::utils::GateType::kEqz, temp);
            temp = circ.addConstOpGate(common::utils::GateType::kConstMul, temp, Ring(-1));
            new_edge_dest_group_party[j] = circ.addConstOpGate(common::utils::GateType::kConstAdd, temp, Ring(1));
        }
        new_edge_dest_group[i] = std::move(new_edge_dest_group_party);
    }

    // Compute overall sigd for new edges
    for (int i = 0; i < nC; ++i) {
        if (add_ESizes[i] == 0) {
            new_edge_overall_sigd_values[i] = {};
            continue;
        }
        std::vector<wire_t> subg_new_edge_sigd_values(add_ESizes[i]);
        auto [prop_out_key, prop_out_v] = 
            circ.addGroupwisePropagateGate(IIn[i], data_e, new_edge_dest_group[i], permutation);
        auto [out_ind, ind_output_key, ind_output_v] = 
            circ.addGroupwiseIndexGate(new_edge_dest_group[i], prop_out_v, permutation);
        for (int j = 0; j < add_ESizes[i]; ++j) {
            subg_new_edge_sigd_values[j] = 
                circ.addGate(common::utils::GateType::kAdd, prop_out_v[j], out_ind[j]);
        }
        new_edge_overall_sigd_values[i] = subg_new_edge_sigd_values;
    }

    // Reorder overall_sigd_values back to source order
    for (int i = 0; i < nC; ++i) {
        if (add_ESizes[i] == 0) {
            continue;
        }
        payloads[i][0] = new_edge_overall_sigd_values[i];
        payloads_out[i] = addSubCircPermList(circ, payloads_out[i][1], payloads[i], permutation);
    }

    // Flatten (original) G
    std::vector<wire_t> src(nV + nE);
    std::vector<wire_t> dst(nV + nE);
    std::vector<wire_t> isV(nV + nE);
    std::vector<wire_t> data(nV + nE);
    std::vector<wire_t> sigs(nV + nE);
    std::vector<wire_t> sigv(nV + nE);
    std::vector<wire_t> sigd(nV + nE);

    int index = 0; 

    for (int i = 0; i < nC; ++i) {
        for (int j = 0; j < VSizes[i]; ++j) {
            src[index] = vertex_src_values[i][j];
            dst[index] = vertex_dst_values[i][j];
            isV[index] = vertex_isV_values[i][j];
            data[index] = vertex_data_values[i][j];
            sigs[index] = vertex_sigs_values[i][j];
            sigv[index] = vertex_sigv_values[i][j];
            sigd[index] = vertex_sigd_values[i][j];
            index++;
        }
    }

    for (int i = 0; i < nC; ++i) {
        for (int j = 0; j < ESizes[i]; ++j) {
            src[index] = edge_src_values[i][j];
            dst[index] = edge_dst_values[i][j];
            isV[index] = edge_isV_values[i][j];
            data[index] = edge_data_values[i][j];
            sigs[index] = edge_sigs_values[i][j];
            sigv[index] = edge_sigv_values[i][j];
            sigd[index] = edge_sigd_values[i][j];
            index++;
        }
    }

    // Add in dummy entries for edges to Vout_agg, Vin_agg
    std::vector<wire_t> Gout(vec_size);
    std::vector<wire_t> Gin(vec_size);

    for (int i = 0; i < nV; ++i) {
        Gout[i] = Vout_agg[i];
        Gin[i] = Vin_agg[i];
    }
    for (int i = 0; i < nE; ++i) {
        Gout[nV + i] = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, 0);
        Gin[nV + i] = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, 0);
    }
    
    // Reorder G, Gout, Gin to source order 
    std::vector<std::vector<wire_t>> G_payload(5);
    G_payload[0] = sigs;
    G_payload[1] = sigv;
    G_payload[2] = sigd;
    G_payload[3] = Gout;
    G_payload[4] = Gin;

    std::vector<std::vector<wire_t>> G_reordered = addSubCircPermList(circ, sigs, G_payload, permutation);

    // Update sigs of existing entries using a running wire_t accumulator
    std::vector<wire_t> updated_sigs(vec_size);
    wire_t Gout_acc = zero_wire;
    for (int i = 0; i < vec_size; ++i) {
        updated_sigs[i] = circ.addGate(common::utils::GateType::kAdd, G_reordered[0][i], Gout_acc);
        Gout_acc = circ.addGate(common::utils::GateType::kAdd, Gout_acc, G_reordered[3][i]);
    }

    // Reorder G, Gin to destination order
    G_payload[0] = updated_sigs;
    G_payload[1] = G_reordered[1];
    G_payload[2] = G_reordered[2];
    G_payload[3] = G_reordered[4]; // Don't need Gout anymore

    // drop the last element from G_payload
    G_payload.pop_back();
    G_reordered.pop_back();

    G_reordered = addSubCircPermList(circ, G_reordered[2], G_payload, permutation);

    // Update sigd of existing entries
    std::vector<wire_t> updated_sigd(vec_size);
    wire_t Gin_acc = zero_wire;
    for (int i = 0; i < vec_size; ++i) {
        updated_sigd[i] = circ.addGate(common::utils::GateType::kAdd, G_reordered[2][i], Gin_acc);
        Gin_acc = circ.addGate(common::utils::GateType::kAdd, Gin_acc, G_reordered[3][i]);
    }

    // Reorder G to vertex order
    G_payload[0] = G_reordered[0];
    G_payload[1] = G_reordered[1];
    G_payload[2] = updated_sigd;

    G_payload.pop_back();
    G_reordered.pop_back();

    G_reordered = addSubCircPermList(circ, G_reordered[1], G_payload, permutation);

    // Update sigv and set outputs
    index = 0;
    for (int i = 0; i < nV; ++i) { // Vertices
        circ.setAsOutput(src[i]);
        circ.setAsOutput(dst[i]);
        circ.setAsOutput(isV[i]);
        circ.setAsOutput(data[i]);
        circ.setAsOutput(G_reordered[0][i]); // sigs
        circ.setAsOutput(G_reordered[1][i]); // sigv
        circ.setAsOutput(G_reordered[2][i]); // sigd
        index++;
    }
    wire_t index_wire;
    for (int i = 0; i < nC; ++i) { // New edges
        for (int j = 0; j < add_ESizes[i]; ++j){
            circ.setAsOutput(new_edge_src_values[i][j]);
            circ.setAsOutput(new_edge_dst_values[i][j]);
            circ.setAsOutput(new_edge_isV_values[i][j]);
            circ.setAsOutput(new_edge_data_values[i][j]);
            circ.setAsOutput(new_edge_overall_sigs_values[i][j]); // sigs
            index_wire = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, Ring(index));
            circ.setAsOutput(index_wire); // sigv
            circ.setAsOutput(payloads_out[i][0][j]); // sigd
            index++;
        }
    }

    for (int i = 0; i < nE; ++i) { // Existing edges
        circ.setAsOutput(src[nV + i]);
        circ.setAsOutput(dst[nV + i]);
        circ.setAsOutput(isV[nV + i]);
        circ.setAsOutput(data[nV + i]);
        circ.setAsOutput(G_reordered[0][nV + i]); // sigs
        index_wire = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, Ring(index));
        circ.setAsOutput(index_wire); // sigv
        circ.setAsOutput(G_reordered[2][nV + i]); // sigd
        index++;
    }
    
    return circ;
}