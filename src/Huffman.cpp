#include <Huffman.h>
#include <functional>
#include <cstring>
#include <cassert>

bool huffman_node::is_leave()
{
    return child[0] == 0 && child[1] == 0;
}
bool huffman_node::is_branch()
{
    return child[0] != 0 && child[1] != 0;
}

bool huffman::check_validation()
{
    std::function<bool(size_t)>dfs = [&](size_t i)->bool {
        auto &node = tree[i];
        if (node.is_leave()) {
            return true;
        } else if (node.is_branch()) {
            return dfs(node.child[0]) && dfs(node.child[1]);
        } else {
            return false;
        }
    };
    return dfs(root);
}
size_t huffman::insert(const char *s, size_t length)
{
    assert(strlen(s) == length);
    size_t node = root;
    for (size_t i = 0 ; i < length ; i++) {
        assert('0' <= s[i] && s[i] <= '1');
        if (tree[node].child[s[i] - '0'] == 0) {
            tree.push_back(huffman_node());
            tree[node].child[s[i] - '0'] = tree.size() - 1;
        }
        node = tree[node].child[s[i] - '0'];
    }
    assert(tree[node].is_leave());
    return node;
}
size_t huffman::decode(MainDataStream *stream, size_t *limit)
{
    size_t node = root;
    while (tree[node].is_branch())
        node = tree[node].child[stream->read_bits(1, limit)];
    return node;
}

bool huffman_quadruples::check_validation()
{
    if (!huffman::check_validation())
        return false;
    for (size_t i = 1 ; i < tree.size() ; i++)
        if (data[i] != NULL && !tree[i].is_leave())
            return false;
    return true;
}
size_t huffman_quadruples::insert(const char *s, size_t length, const quadruples *__data)
{
    size_t node = huffman::insert(s, length);
    while (data.size() <= node)
        data.push_back(NULL);
    data[node] = __data;
    return node;
}
const quadruples * huffman_quadruples::decode(MainDataStream *stream, size_t *limit)
{
    return data[huffman::decode(stream, limit)];
}

bool huffman_code::check_validation()
{
    if (!huffman::check_validation())
        return false;
    for (size_t i = 1 ; i < tree.size() ; i++)
        if (data[i] != NULL && !tree[i].is_leave())
            return false;
    return true;
}
size_t huffman_code::insert(const char *s, size_t length, const code *__data)
{
    size_t node = huffman::insert(s, length);
    while (data.size() <= node)
        data.push_back(NULL);
    data[node] = __data;
    return node;
}
const code * huffman_code::decode(MainDataStream *stream, size_t *limit)
{
    return data[huffman::decode(stream, limit)];
}