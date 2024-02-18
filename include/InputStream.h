#pragma once

#include <cstdio>
#include <cstdint>
#include <queue>
#include <exception>

#define INPUT_STEAM_BUFFER_SIZE 256

class InputException : std::exception {
public:
    virtual const char* what() const noexcept override {
        return "No more Input";
    }
};



class FileStream {
private:
    FILE *fp;
    uint8_t buffer[INPUT_STEAM_BUFFER_SIZE];
    size_t buffer_bits;
    size_t buffer_bptr;
    size_t frame_offset;

    void check_buffer(void);
    void clear_buffer(void);
    void fill_buffer(void);
public:
    FileStream(const char *filepath) {
        fp = fopen(filepath, "r");
        buffer_bits = buffer_bptr = 0;
        frame_offset = 0;
    }
    FileStream(FILE *__fp) {
        fp = __fp;
        buffer_bits = buffer_bptr = 0;
        frame_offset = 0;
    }

    uint32_t read_bits(size_t bits);
    void read_bytes(void *dst, size_t count);
    void skip_bytes(size_t count);
    bool meet_sync_word(void);
    size_t get_frame_offset(void);
    void set_offset(size_t offset);
    size_t get_offset(void);
};

struct MainDataFrame {
    size_t frame_size;
    uint8_t *main_data;
};

class MainDataStream {
private:
    std::deque<MainDataFrame> frame_list;
    size_t frames_size;
    size_t frame_offset;

    uint8_t buffer[INPUT_STEAM_BUFFER_SIZE];
    size_t buffer_bits;
    size_t buffer_bptr;

    void check_buffer(void);
    void clear_buffer(void);
    void fill_buffer(void);
public:
    MainDataStream() {
        frames_size = frame_offset = 0;
    }
    uint32_t read_bits(size_t bits, size_t *limit);
    void insert_frame(MainDataFrame frame);
    void sync_frame_data(size_t main_data_begin);
};