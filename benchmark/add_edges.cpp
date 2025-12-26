#include <io/netmp.h>
#include <graphdb/offline_evaluator.h>
#include <graphdb/online_evaluator.h>
#include <utils/circuit.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <omp.h>

#include "utils.h"
#include "graphutils.h"

using namespace graphdb;
using json = nlohmann::json;
namespace bpo = boost::program_options;

void printDaglistInfo(const DistributedDaglist& dist_daglist, const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << title << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int nC = dist_daglist.num_clients;
    int total_vertices = 0, total_edges = 0;
    int total_new_edges = 0;
    
    for (int i = 0; i < nC; ++i) {
        total_vertices += dist_daglist.VSizes[i];
        total_edges += dist_daglist.ESizes[i];
        total_new_edges += dist_daglist.InsertE[i].size();
    }
    
    std::cout << "Total Vertices: " << total_vertices << std::endl;
    std::cout << "Total Edges: " << total_edges << std::endl;
    std::cout << "New Edges to Insert: " << total_new_edges << std::endl;
    std::cout << "Total Entries: " << (total_vertices + total_edges) << std::endl;
    std::cout << std::endl;
    
    // Print distribution per client
    std::cout << "Distribution per Client:" << std::endl;
    for (int i = 0; i < nC; ++i) {
        std::cout << "  Client " << i << ": " << dist_daglist.VSizes[i] << " vertices, "
                  << dist_daglist.ESizes[i] << " edges, "
                  << dist_daglist.InsertE[i].size() << " new edges" << std::endl;
    }
    std::cout << std::endl;
    
    // Print sample vertices
    std::cout << "Sample Vertices (first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | Data | isV | sigv | sigs | sigd" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    int count = 0;
    for (int c = 0; c < nC && count < 10; ++c) {
        for (size_t i = 0; i < dist_daglist.VSizes[c] && count < 10; ++i) {
            const auto& v = dist_daglist.VertexLists[c][i];
            std::cout << "  " << std::setw(2) << count << " | "
                      << std::setw(3) << v.src << " | "
                      << std::setw(3) << v.dst << " | "
                      << std::setw(4) << v.data << " | "
                      << std::setw(3) << v.isV << " | "
                      << std::setw(4) << v.sigv << " | "
                      << std::setw(4) << v.sigs << " | "
                      << std::setw(4) << v.sigd << std::endl;
            count++;
        }
    }
    std::cout << std::endl;
    
    // Print sample edges
    std::cout << "Sample Existing Edges (first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | Data | isV | sigv | sigs | sigd" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    count = 0;
    for (int c = 0; c < nC && count < 10; ++c) {
        for (size_t i = 0; i < dist_daglist.ESizes[c] && count < 10; ++i) {
            const auto& e = dist_daglist.EdgeLists[c][i];
            std::cout << "  " << std::setw(2) << count << " | "
                      << std::setw(3) << e.src << " | "
                      << std::setw(3) << e.dst << " | "
                      << std::setw(4) << e.data << " | "
                      << std::setw(3) << e.isV << " | "
                      << std::setw(4) << e.sigv << " | "
                      << std::setw(4) << e.sigs << " | "
                      << std::setw(4) << e.sigd << std::endl;
            count++;
        }
    }
    std::cout << std::endl;
    
    // Print new edges to insert
    std::cout << "New Edges to Insert (first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | Data | isV | sigs | sigv | sigd" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    count = 0;
    for (int c = 0; c < nC && count < 10; ++c) {
        for (size_t i = 0; i < dist_daglist.InsertE[c].size() && count < 10; ++i) {
            const auto& e = dist_daglist.InsertE[c][i];
            std::cout << "  " << std::setw(2) << count << " | "
                      << std::setw(3) << e.src << " | "
                      << std::setw(3) << e.dst << " | "
                      << std::setw(4) << e.data << " | "
                      << std::setw(3) << e.isV << " | "
                      << std::setw(4) << e.sigs << " | "
                      << std::setw(4) << e.sigv << " | "
                      << std::setw(4) << e.sigd << std::endl;
            count++;
        }
    }
    std::cout << std::string(60, '=') << std::endl;
}

void printOutputs(const std::vector<Ring>& outputs, const DistributedDaglist& dist_daglist) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "OUTPUT: Graph After Edge Insertion" << std::endl;
    std::cout << std::string(60, '=') << std::endl;
    
    int nC = dist_daglist.num_clients;
    int nV = dist_daglist.nV;
    size_t expected_outputs = 0;
    for (int c = 0; c < nC; ++c) {
        expected_outputs += dist_daglist.VSizes[c] * 7; // vertices
        expected_outputs += dist_daglist.ESizes[c] * 7; // existing edges
        expected_outputs += dist_daglist.InsertE[c].size() * 7; // new edges
    }
    
    if (outputs.size() < expected_outputs) {
        std::cout << "Warning: Received " << outputs.size() << " outputs, expected " << expected_outputs << std::endl;
        return;
    }
    
    std::cout << "Updated Vertices (up to first 10):" << std::endl;
    std::cout << "  ID | Src | Dst | isV | Data | sigs | sigv | sigd" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    
    // The output structure is: ALL vertices (flattened), then for each client: existing edges, then new edges
    size_t output_idx = 0;
    int display_count = 0;
    
    // Print vertices (all flattened, not per-client)
    size_t total_vertices = 0;
    for (int c = 0; c < nC; ++c) {
        total_vertices += dist_daglist.VSizes[c];
    }
    
    size_t limit = std::min(total_vertices, size_t(10));
    for (size_t i = 0; i < limit && output_idx + 6 < outputs.size(); ++i) {
        std::cout << "  " << std::setw(2) << i << " | "
                  << std::setw(3) << outputs[output_idx + 0] << " | "
                  << std::setw(3) << outputs[output_idx + 1] << " | "
                  << std::setw(3) << outputs[output_idx + 2] << " | "
                  << std::setw(4) << outputs[output_idx + 3] << " | "
                  << std::setw(4) << outputs[output_idx + 4] << " | "
                  << std::setw(4) << outputs[output_idx + 5] << " | "
                  << std::setw(4) << outputs[output_idx + 6] << std::endl;
        output_idx += 7;
    }
    
    // Skip remaining vertices
    output_idx = total_vertices * 7;
    
    std::cout << "\nUpdated Edges (up to first 10, including new edges):" << std::endl;
    std::cout << "  ID | Src | Dst | isV | Data | sigs | sigv | sigd" << std::endl;
    std::cout << "  " << std::string(70, '-') << std::endl;
    
    // Now iterate through clients for edges
    display_count = 0;
    for (int c = 0; c < nC && display_count < 10; ++c) {
        // Existing edges
        size_t edge_limit = std::min(static_cast<size_t>(dist_daglist.ESizes[c]), size_t(10 - display_count));
        for (size_t i = 0; i < edge_limit && output_idx + 6 < outputs.size(); ++i) {
            std::cout << "  " << std::setw(2) << display_count << " | "
                      << std::setw(3) << outputs[output_idx + 0] << " | "
                      << std::setw(3) << outputs[output_idx + 1] << " | "
                      << std::setw(3) << outputs[output_idx + 2] << " | "
                      << std::setw(4) << outputs[output_idx + 3] << " | "
                      << std::setw(4) << outputs[output_idx + 4] << " | "
                      << std::setw(4) << outputs[output_idx + 5] << " | "
                      << std::setw(4) << outputs[output_idx + 6] << std::endl;
            output_idx += 7;
            display_count++;
        }
        output_idx += (dist_daglist.ESizes[c] - edge_limit) * 7;
        
        // New edges
        if (display_count < 10) {
            edge_limit = std::min(dist_daglist.InsertE[c].size(), size_t(10 - display_count));
            for (size_t i = 0; i < edge_limit && output_idx + 6 < outputs.size(); ++i) {
                std::cout << "  " << std::setw(2) << display_count << " | "
                          << std::setw(3) << outputs[output_idx + 0] << " | "
                          << std::setw(3) << outputs[output_idx + 1] << " | "
                          << std::setw(3) << outputs[output_idx + 2] << " | "
                          << std::setw(4) << outputs[output_idx + 3] << " | "
                          << std::setw(4) << outputs[output_idx + 4] << " | "
                          << std::setw(4) << outputs[output_idx + 5] << " | "
                          << std::setw(4) << outputs[output_idx + 6] << " *NEW*" << std::endl;
                output_idx += 7;
                display_count++;
            }
            output_idx += (dist_daglist.InsertE[c].size() - edge_limit) * 7;
        } else {
            output_idx += dist_daglist.InsertE[c].size() * 7;
        }
    }
    
    std::cout << std::string(60, '=') << std::endl;
}

common::utils::Circuit<Ring> generateCircuit(int nP, int pid, DistributedDaglist dist_daglist) {

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

    int add_nE = 0;
    std::vector<int> add_ESizes(nC, 0);
    for (int i = 0; i < nC; ++i) {
        add_ESizes[i] = dist_daglist.InsertE[i].size();
        add_nE += add_ESizes[i];
    }

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

    auto zero_wire = circ.addGate(kSub, vertex_src_values[0][0], vertex_src_values[0][0]);

    // Compute aggregated V_in, V_out
    std::vector<wire_t> Vout_agg(nV);
    std::vector<wire_t> Vin_agg(nV);
    for (size_t j = 0; j < nV; ++j) {
        // start with client 0's value then add others
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
    OffIn[0] = Vin_agg[0];
    OffOut[0] = zero_wire;
    for (size_t j = 1; j < nV; ++j) {
        OffIn[j] = circ.addGate(common::utils::GateType::kAdd, OffIn[j - 1], Vin_agg[j]);
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
    // Initialize data_e with zero wires to avoid uninitialized memory
    std::vector<wire_t> data_e(nV, zero_wire);
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

    // Prepare permutations for groupwise operations
    // Need three separate permutation structures: T1 (for IOut/IIn), T2 (for new_edge groups), and reverse
    // Permutations set as identity for benchmarking
    std::vector<std::vector<std::vector<int>>> permutation_t1(nC);  // For IOut/IIn (size nV)
    std::vector<std::vector<std::vector<int>>> permutation_t2(nC);  // For new edge groups (size add_ESizes[i])
    std::vector<std::vector<std::vector<int>>> permutation1(nC);    // For index operations

    for (size_t i = 0; i < nC; ++i) {
        // T1 permutation (size nV for vertices)
        std::vector<int> t1_perm(nV);
        for (size_t j = 0; j < nV; ++j) {
            t1_perm[j] = j;
        }

        // T2 permutation (size add_ESizes[i] for new edges)
        std::vector<int> t2_perm(add_ESizes[i]);
        for (size_t j = 0; j < add_ESizes[i]; ++j) {
            t2_perm[j] = j;
        }
        
        // For party 0, we need all parties' permutations
        // For other parties, just their own permutation
        if (pid == 0) {
            // Party 0 stores all parties' permutations
            permutation_t1[i].push_back(t1_perm);
            permutation_t2[i].push_back(t2_perm);
            for (int p = 1; p < nP; ++p) {
                permutation_t1[i].push_back(t1_perm);
                permutation_t2[i].push_back(t2_perm);
            }
        } else {
            // Other parties store their own permutations
            permutation_t1[i].push_back(t1_perm);
            permutation_t2[i].push_back(t2_perm);
        }

        permutation1[i].push_back(t2_perm);
        if (pid == 0) {
            for (int p = 1; p < nP; ++p) {
                permutation1[i].push_back(t2_perm);
            }
        }
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
            circ.addGroupwisePropagateSubcircuit(IOut[i], data_e, new_edge_src_group[i], permutation_t1[i], permutation_t2[i], permutation_t2[i], pid);
        auto [out_ind, ind_output_key, ind_output_v] = 
            circ.addGroupwiseIndexSubcircuit(new_edge_src_group[i], prop_out_v, permutation1[i], pid);
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
            payloads_out[i] = circ.addSubCircPermList(new_edge_sigd_values[i], payloads[i], permutation1[i]);
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
            circ.addGroupwisePropagateSubcircuit(IIn[i], data_e, new_edge_dest_group[i], permutation_t1[i], permutation_t2[i], permutation_t2[i], pid);
        auto [out_ind, ind_output_key, ind_output_v] = 
            circ.addGroupwiseIndexSubcircuit(new_edge_dest_group[i], prop_out_v, permutation1[i], pid);
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
        payloads_out[i] = circ.addSubCircPermList(payloads_out[i][1], payloads[i], permutation1[i]);
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

    // Update position maps for vertices
    
    std::vector<wire_t> updated_sigs(vec_size);
    std::vector<wire_t> updated_sigd(vec_size);
    for (int i = 0; i < nV; ++i) {
        updated_sigs[i] = circ.addGate(common::utils::GateType::kAdd, sigs[i], OffOut[i]);
        updated_sigd[i] = circ.addGate(common::utils::GateType::kAdd, sigd[i], OffIn[i]);
    }

    std::vector<int> base_perm2(vec_size);
    for (size_t i = 0; i < vec_size; ++i) {
        base_perm2[i] = static_cast<int>(i);
    }
    std::vector<std::vector<int>> permutation2;
    permutation2.push_back(base_perm2);
    if (pid == 0) {
        for (int i = 1; i < nP; ++i) {
            permutation2.push_back(base_perm2);
        }
    }

    // Propagate OffOut
    std::vector<wire_t> Gout(vec_size);
    for (int i = 0; i < nV; ++i) {
        Gout[i] = OffOut[i];
    }
    for (int i = 0; i < nE; ++i) {
        Gout[nV + i] = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, 0);
    }
    auto prop_out = circ.addSubCircPropagate(sigs, Gout, nV, permutation2);
    
    // Reorder sigv to source order
    auto sigv_s = circ.addSubCircPermList(sigs, {sigv}, permutation2)[0];

    // Reorder to vertex order
    auto out_v = circ.addSubCircPermList(sigv_s, {prop_out}, permutation2)[0];

    // Update sigs for existing edges
    for (int i = nV; i < vec_size; ++i) {
        updated_sigs[i] = circ.addGate(common::utils::GateType::kAdd, sigs[i], out_v[i]);
    }

    // Propagate OffIn
    std::vector<wire_t> Gin(vec_size);
    for (int i = 0; i < nV; ++i) {
        Gin[i] = OffIn[i];
    }
    for (int i = 0; i < nE; ++i) {
        Gin[nV + i] = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, 0);
    }
    auto prop_in = circ.addSubCircPropagate(sigs, Gin, nV, permutation2, true);

    // Reorder to vertex order
    auto in_v = circ.addSubCircPermList(sigv, {prop_in}, permutation2)[0];

    // Update sigs for existing edges
    for (int i = nV; i < vec_size; ++i) {
        updated_sigd[i] = circ.addGate(common::utils::GateType::kAdd, sigd[i], in_v[i]);
    }
    
    index = 0;
    for (int i = 0; i < nV; ++i) { // Vertices
        circ.setAsOutput(src[i]);
        circ.setAsOutput(dst[i]);
        circ.setAsOutput(isV[i]);
        circ.setAsOutput(data[i]);
        circ.setAsOutput(updated_sigs[i]); 
        circ.setAsOutput(sigv[i]);
        circ.setAsOutput(updated_sigd[i]); 
        index++;
    }

    wire_t index_wire;
    int edge_index = 0;
    for (int i = 0; i < nC; ++i) { // Edges
        for (int j = 0; j < ESizes[i]; ++j) { // Existing edges
            circ.setAsOutput(edge_src_values[i][j]);
            circ.setAsOutput(edge_dst_values[i][j]);
            circ.setAsOutput(edge_isV_values[i][j]);
            circ.setAsOutput(edge_data_values[i][j]);
            circ.setAsOutput(updated_sigs[nV + edge_index]); // sigs
            index_wire = circ.addConstOpGate(common::utils::GateType::kConstAdd, zero_wire, Ring(index));
            circ.setAsOutput(index_wire); // sigv
            circ.setAsOutput(updated_sigd[nV + edge_index]); // sigd
            index++;
            edge_index++;
        }
        for (int j = 0; j < add_ESizes[i]; ++j){ // New edges
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
    
    return circ;
}

void benchmark(const bpo::variables_map& opts) {

    bool save_output = false;
    std::string save_file;
    if (opts.count("output") != 0) {
        save_output = true;
        save_file = opts["output"].as<std::string>();
    }

    auto num_vert = opts["num-vert"].as<size_t>();
    auto num_edge = opts["num-edge"].as<size_t>();
    auto vec_size = num_vert + num_edge;
    auto nP = opts["num-parties"].as<int>();
    auto nC = opts["num-clients"].as<int>();
    auto latency = opts["latency"].as<double>();
    auto pid = opts["pid"].as<size_t>();
    auto threads = opts["threads"].as<size_t>();
    auto seed = opts["seed"].as<size_t>();
    auto repeat = opts["repeat"].as<size_t>();
    auto port = opts["port"].as<int>();
    auto use_pking = opts["use-pking"].as<bool>();

    omp_set_nested(1);
    if (nP < 10) { omp_set_num_threads(nP); }
    else { omp_set_num_threads(10); }

    std::cout << "Starting benchmarks" << std::endl;

    std::string net_config = opts.count("net-config") ? opts["net-config"].as<std::string>() : "";
    std::shared_ptr<io::NetIOMP> network = createNetwork(pid, nP, latency, port,
                                                          opts["localhost"].as<bool>(),
                                                          net_config);

    json output_data;
    output_data["details"] = {{"num_parties", nP},
                              {"num_clients", nC},
                              {"num_vert", num_vert},
                              {"num_edge", num_edge},
                              {"operation", "insert_edges"},
                              {"latency (ms)", latency},
                              {"pid", pid},
                              {"threads", threads},
                              {"seed", seed},
                              {"repeat", repeat}};
    output_data["benchmarks"] = json::array();

    std::cout << "--- Details ---" << std::endl;
    for (const auto& [key, value] : output_data["details"].items()) {
        std::cout << key << ": " << value << std::endl;
    }
    std::cout << std::endl;

    // Generate test graph
    Ring nV = static_cast<Ring>(num_vert);
    Ring nE = static_cast<Ring>(num_edge);
    
    std::cout << "============================\n" << std::endl;
    std::cout << "Generating random inputs " << std::endl;
    std::cout << "Generating scale-free graph with nV=" << nV << ", nE=" << nE << " (seed=" << seed << ")" << std::endl;
    auto edges = generate_scale_free(nV, nE, seed);
    std::cout << "Generated " << edges.size() << " edges" << std::endl;
    
    std::cout << "Building daglist..." << std::endl;
    auto daglist = build_daglist(nV, edges);
    std::cout << "Built daglist with " << daglist.size() << " entries" << std::endl;
    
    // Distribute daglist across clients
    std::cout << "Distributing daglist across " << nC << " clients..." << std::endl;
    auto dist_daglist = distribute_daglist(daglist, nC);
    
    // Generate random edges to insert (insert 20% of current edges)
    Ring num_inserts = static_cast<Ring>(nE * 0.2);
    std::cout << "Generating random edges to insert: " << num_inserts << " edges..." << std::endl;
    dist_daglist = generate_random_edges_to_insert(dist_daglist, num_inserts, seed);


    StatsPoint start(*network);
    network->sync();

    auto circ = generateCircuit(nP, pid, dist_daglist).orderGatesByLevel();
    network->sync();

    std::cout << "--- Circuit ---" << std::endl;
    std::cout << circ << std::endl;
    std::cout << "Party " << pid << ": circuit levels = " << circ.gates_by_level.size()
              << ", total gates = " << circ.num_gates
              << ", total outputs = " << circ.outputs.size() << std::endl;
    
    std::unordered_map<common::utils::wire_t, int> input_pid_map;
    for (const auto& g : circ.gates_by_level[0]) {
        if (g->type == common::utils::GateType::kInp) {
            input_pid_map[g->out] = 1;
        }
    }

    std::cout << "Starting preprocessing" << std::endl;
    StatsPoint preproc_start(*network);
    int latency_us = static_cast<int>(latency * 1000);
    OfflineEvaluator off_eval(nP, pid, network, circ, threads, seed, latency_us, use_pking);
    auto preproc = off_eval.run(input_pid_map);
    std::cout << "Preprocessing complete" << std::endl;
    network->sync();
    StatsPoint preproc_end(*network);

    std::cout << "Setting inputs" << std::endl;
    OnlineEvaluator eval(nP, pid, network, std::move(preproc), circ, threads, seed, latency_us, use_pking);
    
    std::unordered_map<common::utils::wire_t, Ring> inputs;
    
    // Collect all input wires owned by this party
    std::vector<common::utils::wire_t> input_wires;
    for (const auto& [wire, owner] : input_pid_map) {
        if (owner == static_cast<int>(pid)) {
            input_wires.push_back(wire);
        }
    }
    
    // Sort to ensure consistent ordering
    std::sort(input_wires.begin(), input_wires.end());
    
    std::cout << "Setting inputs for party " << pid << std::endl;
    
    // Only party 1 sets inputs
    std::vector<Ring> graph_input_values;

    if (pid == 1) {

        // Print distribution info
        std::cout << "\n=== Daglist Distribution ===" << std::endl;
        for (int i = 0; i < nC; ++i) {
            std::cout << "Client " << i << ": " << dist_daglist.VSizes[i] << " vertices, "
                    << dist_daglist.ESizes[i] << " edges, " << dist_daglist.InsertE[i].size() << " to insert" << std::endl;
        }
        std::cout << "============================\n" << std::endl;
    
        std::vector<Ring> all_input_values;
        
    // Collect all vertex and edge fields for all clients
    for (int c = 0; c < nC; ++c) {
        for (size_t i = 0; i < dist_daglist.VSizes[c]; ++i) {
            all_input_values.push_back(dist_daglist.VertexLists[c][i].src);
            all_input_values.push_back(dist_daglist.VertexLists[c][i].dst);
            all_input_values.push_back(dist_daglist.VertexLists[c][i].isV);
            all_input_values.push_back(dist_daglist.VertexLists[c][i].data);
            all_input_values.push_back(dist_daglist.VertexLists[c][i].sigs);
            all_input_values.push_back(dist_daglist.VertexLists[c][i].sigv);
            all_input_values.push_back(dist_daglist.VertexLists[c][i].sigd);
        }

        for (size_t i = 0; i < dist_daglist.ESizes[c]; ++i) {
            all_input_values.push_back(dist_daglist.EdgeLists[c][i].src);
            all_input_values.push_back(dist_daglist.EdgeLists[c][i].dst);
            all_input_values.push_back(dist_daglist.EdgeLists[c][i].isV);
            all_input_values.push_back(dist_daglist.EdgeLists[c][i].data);
            all_input_values.push_back(dist_daglist.EdgeLists[c][i].sigs);
            all_input_values.push_back(dist_daglist.EdgeLists[c][i].sigv);
            all_input_values.push_back(dist_daglist.EdgeLists[c][i].sigd);
        }
    }

    // Collect new edges to be inserted
    for (int c = 0; c < nC; ++c) {
        for (size_t i = 0; i < dist_daglist.InsertE[c].size(); ++i) {
            all_input_values.push_back(dist_daglist.InsertE[c][i].src);
            all_input_values.push_back(dist_daglist.InsertE[c][i].dst);
            all_input_values.push_back(dist_daglist.InsertE[c][i].isV);
            all_input_values.push_back(dist_daglist.InsertE[c][i].data);
            all_input_values.push_back(dist_daglist.InsertE[c][i].sigs);
            all_input_values.push_back(dist_daglist.InsertE[c][i].sigv);
            all_input_values.push_back(dist_daglist.InsertE[c][i].sigd);
        }
    }

    // Collect VIn, VOut
    for (int c = 0; c < nC; ++c) {
        for (size_t j = 0; j < nV; ++j) {
            all_input_values.push_back(dist_daglist.VOut[c][j]);
            all_input_values.push_back(dist_daglist.VIn[c][j]);
        }
    }

                // Map collected values into circuit input wires (in order)
        size_t wire_idx = 0;
        for (size_t i = 0; i < all_input_values.size() && wire_idx < input_wires.size(); ++i) {
            inputs[input_wires[wire_idx++]] = all_input_values[i];
        }
        
        // Store for verification
        graph_input_values = all_input_values;
        
    }

    std::cout << "Total inputs set by party " << pid << ": " << inputs.size() << std::endl;
    
    if (pid == 1) {
        std::cout << "Party 1 setting " << inputs.size() << " actual input values" << std::endl;
    } else {
        std::cout << "Party " << pid << " setting " << inputs.size() << " empty inputs (participant in MPC)" << std::endl;
    }
    
    eval.setInputs(inputs);
    network->sync();
    
    std::cout << "Starting online evaluation" << std::endl;
    StatsPoint online_start(*network);
    for (size_t i = 0; i < circ.gates_by_level.size(); ++i) {
        eval.evaluateGatesAtDepth(i);
        network->sync();
        network->flush();
    }
    network->sync();
    StatsPoint online_end(*network);
    std::cout << "Online evaluation complete" << std::endl;

    std::cout << "Getting outputs..." << std::endl;
    network->flush();
    auto outputs = eval.getOutputs();
    network->sync();
    std::cout << "Number of outputs: " << outputs.size() << std::endl;
    
    // Print formatted inputs
    if (pid == 1) {
        printDaglistInfo(dist_daglist, "INPUT: Graph Before Edge Insertion");
    }
    
    // Print formatted outputs
    if (pid == 1 && outputs.size() > 0) {
        printOutputs(outputs, dist_daglist);
    }

    StatsPoint end(*network);

     auto preproc_rbench = preproc_end - preproc_start;
    auto online_rbench = online_end - online_start;
    auto total_rbench = end - start;
    output_data["benchmarks"].push_back(preproc_rbench);
    output_data["benchmarks"].push_back(online_rbench);
    output_data["benchmarks"].push_back(total_rbench);

    size_t pre_bytes_sent = 0;
    for (const auto& val : preproc_rbench["communication"]) {
        pre_bytes_sent += val.get<int64_t>();
    }
    size_t online_bytes_sent = 0;
    for (const auto& val : online_rbench["communication"]) {
        online_bytes_sent += val.get<int64_t>();
    }
    size_t total_bytes_sent = 0;
    for (const auto& val : total_rbench["communication"]) {
        total_bytes_sent += val.get<int64_t>();
    }

    std::cout << "preproc time: " << preproc_rbench["time"] << " ms" << std::endl;
    std::cout << "preproc sent: " << pre_bytes_sent << " bytes" << std::endl;
    std::cout << "online time: " << online_rbench["time"] << " ms" << std::endl;
    std::cout << "online sent: " << online_bytes_sent << " bytes" << std::endl;
    std::cout << "total time: " << total_rbench["time"] << " ms" << std::endl;
    std::cout << "total sent: " << total_bytes_sent << " bytes" << std::endl;
    std::cout << std::endl;

    output_data["stats"] = {{"peak_virtual_memory", peakVirtualMemory()},
                            {"peak_resident_set_size", peakResidentSetSize()}};

    std::cout << "--- Statistics ---" << std::endl;
    for (const auto& [key, value] : output_data["stats"].items()) {
        std::cout << key << ": " << value << std::endl;
    }
    std::cout << std::endl;

    if (save_output) {
        saveJson(output_data, save_file);
    }
}

bpo::options_description programOptions() {
    bpo::options_description desc("Following options are supported by config file too.");
    desc.add_options()
        ("num-parties,n", bpo::value<int>()->required(), "Number of parties.")
        ("num-clients", bpo::value<int>()->default_value(2), "Number of parties.")
        ("num-vert", bpo::value<size_t>()->default_value(1000), "Number of vertices in the graph.")
        ("num-edge", bpo::value<size_t>()->default_value(4000), "Number of edges in the graph.")
        ("num-payloads", bpo::value<size_t>()->default_value(1), "Number of payload vectors.")
        ("latency,l", bpo::value<double>()->default_value(0.5), "Network latency in ms.")
        ("pid,p", bpo::value<size_t>()->required(), "Party ID.")
        ("threads,t", bpo::value<size_t>()->default_value(6), "Number of threads (recommended 6).")
        ("seed", bpo::value<size_t>()->default_value(200), "Value of the random seed.")
        ("net-config", bpo::value<std::string>(), "Path to JSON file containing network details of all parties.")
        ("localhost", bpo::bool_switch(), "All parties are on same machine.")
        ("port", bpo::value<int>()->default_value(10000), "Base port for networking.")
        ("output,o", bpo::value<std::string>(), "File to save benchmarks.")
        ("repeat,r", bpo::value<size_t>()->default_value(1), "Number of times to run benchmarks.")
        ("use-pking", bpo::value<bool>()->default_value(true), "Use king party for reconstruction (true) or direct reconstruction (false).");
  return desc;
}

int main(int argc, char* argv[]) {
    auto prog_opts(programOptions());
    bpo::options_description cmdline("Benchmark secure edge insertion circuit.");
    cmdline.add(prog_opts);
    cmdline.add_options()(
      "config,c", bpo::value<std::string>(),
      "configuration file for easy specification of cmd line arguments")(
      "help,h", "produce help message");
    bpo::variables_map opts;
    bpo::store(bpo::command_line_parser(argc, argv).options(cmdline).run(), opts);
    if (opts.count("help") != 0) {
        std::cout << cmdline << std::endl;
        return 0;
    }
    if (opts.count("config") > 0) {
        std::string cpath(opts["config"].as<std::string>());
        std::ifstream fin(cpath.c_str());
        if (fin.fail()) {
            std::cerr << "Could not open configuration file at " << cpath << std::endl;
            return 1;
        }
        bpo::store(bpo::parse_config_file(fin, prog_opts), opts);
    }
    try {
        bpo::notify(opts);
        if (!opts["localhost"].as<bool>() && (opts.count("net-config") == 0)) {
            throw std::runtime_error("Expected one of 'localhost' or 'net-config'");
        }
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    try {
        benchmark(opts);
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << "\nFatal error" << std::endl;
        return 1;
    }
    return 0;
}
// usage: ./../run.sh add_edges --num-parties 2 --num-clients 2 --num-vert 1000 --num-edge 4000