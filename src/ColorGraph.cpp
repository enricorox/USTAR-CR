//
// Created by enrico on 30/09/24.
//

#include "ColorGraph.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm> // fix compilation on the server
#include <cassert>

using namespace std;

inline char complement(char c){
    switch(c){
        // upper
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        // lower
        case 'a': return 't';
        case 't': return 'a';
        case 'c': return 'g';
        case 'g': return 'c';
        // unknown
        default: return 'N';
    }
}

std::string reverse_complement(std::string sequence){
    for(auto &c: sequence) c = complement(c);
    std::reverse(sequence.begin(), sequence.end());
    return sequence;
}

void encode_RLE(std::vector<color_id_t> colors, std::vector<color_id_t> &values, std::vector<size_t> &counts){
    color_id_t prev = colors[0];
    size_t count = 1;

    // join multiple RLEs
    if(!counts.empty() && colors[0] == values.back()){
        count += counts.back();
        //values.erase(values.end() - 1);
        //counts.erase(counts.end() - 1);
        values.pop_back();
        counts.pop_back();
    }

    for(size_t i = 1; i < colors.size(); i++){
        if(prev == colors[i])
            count++;
        else{
            values.push_back(prev);
            counts.push_back(count);
            prev = colors[i];
            count = 1;
        }
    }
    values.push_back(prev);
    counts.push_back(count);
}

ColorGraph::ColorGraph(std::string sequences_file_name, std::string colors_file_name, int kmer_length) {
    this->sequences_file_name = std::move(sequences_file_name);
    this->colors_file_name = std::move(colors_file_name);
    this->kmer_length = kmer_length;

    if(kmer_length < 3){
        cerr << "Error ColorGraph(): kmer_length must be at least 3" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "* Building the graph..." << endl;
    build_graph();
    print_stats();
    cout << "* Computing a path cover..." << endl;
    compute_path_cover();
    cout << "* Path cover ready!" << endl;
}

void ColorGraph::compute_path_cover() {
    for(auto& node: nodes){
        auto &seed = node.second;
        auto seed_id = node.first;

        // skip visited nodes
        if(seed.is_visited())
            continue;

        // visit the node and build a path
        seed.visit();
        Path path(seed_id);

        // extend direct
        while(has_next(path))
            path.extend(next(path));

        // extend reverse
        path.reverse();
        while(has_next(path))
            path.extend(next(path));

        // collect paths
        paths.push_back(path);
    }
}

void ColorGraph::build_graph() {
    // reading and decoding colors
    vector<color_id_t> colors = decode_RLE_colors();

    ifstream sequences_file(sequences_file_name);
    if(!sequences_file.is_open()){
        cerr << "Error build_graph(): cannot open file " << sequences_file_name << endl;
        exit(EXIT_FAILURE);
    }

    cout << "** Reading sequences " << sequences_file_name << endl;

    node_id_t node_id = 0;
    string line;
    while(getline(sequences_file, line)){
        if(line[0] == '>') // got a def-line
            continue;

        // extract info
        long length = static_cast<long>(line.size());

        if(length < kmer_length){
            cerr << "Error build_graph(): sequence too short!" << endl;
            exit(EXIT_FAILURE);
        }

        long n_kmer = length - kmer_length + 1;
        vector<color_id_t> sequence_colors = vector(colors.begin() + tot_kmers, colors.begin() + tot_kmers + n_kmer);

        // add node to the graph
        nodes[node_id] = Node(line, sequence_colors);
        nodes_head[sequence_colors.front()].push_back(node_id);
        nodes_tail[sequence_colors.back()].push_back(node_id);

        // update counters
        tot_kmers += n_kmer;
        node_id++;
    }
    if(tot_kmers != static_cast<long>(colors.size())){
        cerr    << "Error build_graph(): wrong number of colors\n"
                << "    tot_kmers = " << tot_kmers << "\n"
                << "    colors.size() = " << colors.size() << endl;
        exit(EXIT_FAILURE);
    }
}

void ColorGraph::print_stats() const{
    cout << "\n";
    cout << "Graph stats:\n";
    cout << "   number of kmers:        " << tot_kmers << "\n";
    cout << "   number of sequences:    " << nodes.size() << "\n";
    cout << endl;
}

std::vector<color_id_t> ColorGraph::decode_RLE_colors() {
    ifstream colors_file(colors_file_name);
    if(!colors_file.is_open()){
        cerr << "Error build_graph(): cannot open file " << sequences_file_name << endl;
        exit(EXIT_FAILURE);
    }

    std::vector<color_id_t> colors;
    string line;
    ulong value = 0;
    ulong count = 0;
    ulong pos = 0;

    cout << "** Reading colors " << colors_file_name << "..." << endl;

    while(colors_file >> line){
        pos = line.find(':');
        if(pos == string::npos) // value
            count = 1;
        else // value:count
            count = stoul(line.substr(pos + 1, string::npos));

        value = stoul(line.substr(0, pos));

        for(ulong i = 0; i < count; i++){
            colors.push_back(value);
        }
    }
    return colors;
}

bool ColorGraph::has_next(Path path) {
    color_id_t color;
    auto tail = path.get_tail_node_id();

    // check tail orientation and get the correct color
    if(path.get_tail_orientation() == orientation_t::direct)
        color = nodes[tail].colors.back();
    else
        color = nodes[tail].colors.front();

    // start searching
    size_t max_attempts = nodes.size();
    for(size_t i = 0; i < max_attempts; i++) {
        // search in nodes heads (direct)
        if (!nodes_head[color].empty()) {
            next_node = nodes_head[color].front(); // [C----->

            if (nodes[next_node].is_visited()) { // break loops (and auto-loops)!
                nodes_head[color].pop_front();
                continue;
                // return has_next(path);
            }

            next_orientation = orientation_t::direct;
            return true;
        }

        // search in nodes tails (reverse)
        if (!nodes_tail[color].empty()) {
            next_node = nodes_tail[color].front();

            if (nodes[next_node].is_visited()) { // break loops (and auto-loops)!
                nodes_tail[color].pop_front();
                continue;
                // return has_next(path);
            }

            next_orientation = orientation_t::reverse;
            return true;
        }

        // no colors
        return false;
    }

    cerr << "Warning has_next(): max attempts reached!" << endl;

    return false; // should not reach here
}

// undefined behavior if has_next() is not called before
oriented_node_t ColorGraph::next(Path path) {
    nodes[next_node].visit();
    return oriented_node_t{next_node, next_orientation};
}

void ColorGraph::write_cover(std::string sequences_filename, std::string colors_filename) {
    cout << "** Writing sequences to " << sequences_filename << endl;
    ofstream sequences_file(sequences_filename);

    vector<color_id_t> values;
    vector<size_t> counts;
    for(auto &path: paths){
        for(size_t i = 0; i < path.length(); i++){
            node_id_t node_id = path.get_n_node_id(i);
            Node &node = nodes.at(node_id); //nodes[node_id];
            if(path.get_n_orientation(i) == orientation_t::reverse)
                node.reverse();
            encode_RLE(node.colors, values, counts);
            sequences_file << ">\n" << node.sequence << "\n";
        }
    }

   assert(values.size() == counts.size());

    cout << "** Writing colors to " << colors_filename << "\n";
    cout << "       number of runs: " << counts.size() << endl;
    ofstream colors_file(colors_filename);
    for(size_t i = 0; i < values.size(); i++) {
        colors_file << values[i]; // value (only for 1-runs)
        if(counts[i] != 1)
            colors_file << ":" << counts[i]; // value:count
        colors_file << "\n";
    }
}

ColorGraph::ColorGraph(const std::vector<std::string>& sequences, const std::vector<std::vector<color_id_t>>& colors,
                       int kmer_length) {
    this->sequences_file_name = "tmp";
    this->colors_file_name = "tmp";
    this->kmer_length = kmer_length;

    if(kmer_length < 3){
        cerr << "Error ColorGraph(): kmer_length must be at least 3" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "* Building the graph..." << endl;
    build_graph(sequences, colors);
    print_stats();
    cout << "* Computing a path cover..." << endl;
    compute_path_cover();
    cout << "* Path cover ready!" << endl;
}

void ColorGraph::build_graph(const std::vector<std::string> &sequences, const std::vector<std::vector<color_id_t>> &colors) {
    assert(sequences.size() == colors.size());

    node_id_t node_id = 0;
    for(size_t i = 0; i < sequences.size(); i++){
        const string &sequence = sequences[i];
        const vector<color_id_t> &sequence_colors = colors[i];
        long length = static_cast<long>(sequence.size());

        if(length < kmer_length){
            cerr << "Error build_graph(): sequence too short!" << endl;
            exit(EXIT_FAILURE);
        }

        long n_kmer = length - kmer_length + 1;

        // add node to the graph
        nodes[node_id] = Node(sequence, sequence_colors);
        nodes_head[sequence_colors.front()].push_back(node_id);
        nodes_tail[sequence_colors.back()].push_back(node_id);

        // update counters
        tot_kmers += n_kmer;
        node_id++;
    }
}

Node::Node(std::string sequence, std::vector<color_id_t> colors) {
    this->sequence = std::move(sequence);
    this->colors = std::move(colors);
    this->visited = false;
}

bool Node::is_visited() const {
    return visited;
}

void Node::visit() {
    visited = true;
}

void Node::reverse() {
    sequence = reverse_complement(sequence);
    std::reverse(colors.begin(), colors.end());
}

std::string Node::colors_to_string() {
    string colors_str;

    for(color_id_t color: colors)
        colors_str += to_string(color) + '\n';

    return colors_str;
}

Node::Node() = default;

Path::Path(node_id_t seed) {
    nodes.push_back(seed);
    orientations.push_back(orientation_t::direct);
}

void Path::extend(oriented_node_t node) {
    nodes.push_back(node.node);
    orientations.push_back(node.orientation);
}

void Path::reverse() {
    std::reverse(nodes.begin(), nodes.end());
    std::reverse(orientations.begin(), orientations.end());
    for(auto &orientation: orientations)
        if(orientation == orientation_t::reverse)
            orientation = orientation_t::direct;
        else
            orientation = orientation_t::reverse;
}

orientation_t Path::get_tail_orientation() {
    return orientations.back();
}

color_id_t Path::get_tail_node_id() {
    return nodes.back();
}

size_t Path::length() {
    return nodes.size();
}

orientation_t Path::get_n_orientation(size_t n) {
    return orientations.at(n);
}

color_id_t Path::get_n_node_id(size_t n) {
    return nodes.at(n);
}