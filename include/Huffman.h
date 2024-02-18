#pragma once

struct quadruples {
    int v,w,x,y;
};

struct code {
    int x,y;
};

#include <vector>
#include <cstddef>

struct huffman_node {
    size_t child[2];

    huffman_node()
    {
        child[0] = child[1] = 0;
    }
    bool is_leave();
    bool is_branch();
};

#include <InputStream.h>

struct huffman {
    size_t root;
    std::vector<huffman_node>tree;

    huffman()
    {
        tree = std::vector<huffman_node>(2);
        root = 1;
    }
    bool check_validation();
    size_t insert(const char *s, size_t length);
    size_t decode(MainDataStream *stream, size_t *limit);
};

struct huffman_quadruples : huffman {
    std::vector<const quadruples*>data;

    huffman_quadruples(): huffman() {}
    bool check_validation();
    size_t insert(const char *s, size_t length, const quadruples *__data);
    const quadruples *decode(MainDataStream *stream, size_t *limit);
};

struct huffman_code : huffman {
    std::vector<const code*>data;
    int linbits;

    huffman_code(int __linbits): huffman() {
        linbits = __linbits;
    }
    bool check_validation();
    size_t insert(const char *s, size_t length, const code *__data);
    const code *decode(MainDataStream *stream, size_t *limit);
};

huffman_quadruples* get_quadruple_huffman(int index);
huffman_code* get_code_huffman(int index);