//
// Created by enrico on 30/09/24.
//

#ifndef USTAR_COLORGRAPH_H
#define USTAR_COLORGRAPH_H

#include <vector>
#include <map>
#include <list>

typedef unsigned long node_id_t;
typedef unsigned int color_id_t;

inline char complement(char c);

std::string reverse_complement(std::string sequence);

void encode_RLE(std::vector<color_id_t> colors, std::vector<color_id_t> &values, std::vector<size_t> &counts);

// nodes are default-direct
// head [===> tail
class Node{
    bool visited = false;
public:
    Node(std::string sequence, std::vector<color_id_t> colors);

    Node();

    void visit();

    void reverse();

    std::string colors_to_string();

    bool is_visited() const;
    std::vector<color_id_t> colors;
    std::string sequence;
};

enum orientation_t{
    direct, reverse
};

struct oriented_node_t{
    node_id_t node;
    orientation_t orientation;
};

class Path{
    std::vector<node_id_t> nodes;
    std::vector<orientation_t> orientations;

public:
    Path(node_id_t seed);

    void extend(oriented_node_t node);

    void reverse();

    // std::string sequences();

    // std::string colors();

    color_id_t get_tail_node_id();

    orientation_t get_tail_orientation();

    size_t length();

    orientation_t get_n_orientation(size_t n);

    color_id_t get_n_node_id(size_t n);
};

class ColorGraph{
private:
    std::vector<color_id_t> values;
    std::vector<size_t> counts;
    std::string sequences_file_name;
    std::string colors_file_name;
    int kmer_length;
    long tot_kmers = 0;

    node_id_t next_node{};
    orientation_t next_orientation;

    std::map<node_id_t, Node> nodes;

    // nodes_head[3] contains all the nodes with 3 in the head: [3--->
    std::map<color_id_t, std::list<node_id_t>> nodes_head;
    // nodes_tail[3] contains all the nodes with 3 in the tail: [---3>
    std::map<color_id_t, std::list<node_id_t>> nodes_tail;

    std::vector<Path> paths;

    void build_graph();

    void build_graph(const std::vector<std::string> &sequences, const std::vector<std::vector<color_id_t>> &colors);

    bool has_next(Path path);

    oriented_node_t next(Path path);

    void print_stats() const;

    void compute_path_cover();

    std::vector<color_id_t> decode_RLE_colors();

public:
    ColorGraph(std::string sequences_file_name, std::string colors_file_name, int kmer_length);

    ColorGraph(const std::vector<std::string>& sequences, const std::vector<std::vector<color_id_t>>& colors, int kmer_length = 0);

    void write_colors(std::string colors_filename);

    void write_sequences(std::string sequences_filename);
};

#endif //USTAR_COLORGRAPH_H
