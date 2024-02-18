#include <InputStream.h>
#include <cassert>
#include <cstring>
#include <algorithm>

void FileStream::check_buffer()
{
    if (buffer_bptr == buffer_bits) {
        clear_buffer();
        fill_buffer();
    }
}

void FileStream::clear_buffer()
{
    buffer_bptr = buffer_bits = 0;
}

void FileStream::fill_buffer()
{
    buffer_bits = 8 * fread(buffer, 1, INPUT_STEAM_BUFFER_SIZE, fp);
    if (buffer_bits == 0)
        throw InputException();
}

uint32_t FileStream::read_bits(size_t bits) 
{
    assert(0 <= bits && bits <= 32);
    uint32_t res = 0;
    if (buffer_bptr / 8 == (buffer_bptr + bits) / 8) {
        check_buffer();
        // in a byte
        int index = buffer_bptr / 8;
        int offset = 8 - (bits + buffer_bptr) % 8;
        res = (buffer[index] >> offset) & ((1 << bits) - 1);
        buffer_bptr += bits;
    } else {
        // accross bytes
        // align up
        if (buffer_bptr % 8 != 0) {
            int index = buffer_bptr / 8;
            int readed = 8 - buffer_bptr % 8;
            res = (res << readed) | (buffer[index] & ((1 << readed) - 1));
            buffer_bptr += readed;
            bits -= readed;
        }
        // read full bytes
        while (bits >= 8) {
            check_buffer();
            int index = buffer_bptr / 8;
            res = (res << 8) | buffer[index];
            buffer_bptr += 8;
            bits -= 8; 
        }
        // read left bits
        if (bits > 0) {
            check_buffer();
            int index = buffer_bptr / 8;
            int offset = 8 - bits;
            res = (res << bits) | (buffer[index] >> offset);
            buffer_bptr += bits;
            bits -= bits;
        }
    }
    assert(buffer_bptr <= buffer_bits);
    return res;
}

void FileStream::read_bytes(void *dst, size_t count)
{
    // this must be done
    assert(buffer_bptr % 8 == 0);
    if (count <= (buffer_bits - buffer_bptr) / 8) {
        memcpy(dst, buffer + buffer_bptr / 8, count);
        buffer_bptr += count * 8;
        count = 0;
    } else {
        size_t __count = (buffer_bits - buffer_bptr) / 8;
        memcpy(dst, buffer + buffer_bptr / 8, __count);
        clear_buffer();
        dst = (char *)dst + __count;
        count -= __count;
    }
    if (count > 0) {
        if (fread(dst, count, 1, fp) != 1)
            throw InputException();
    }
}

void FileStream::skip_bytes(size_t count)
{
    assert(buffer_bptr % 8 == 0);
    if (count <= (buffer_bits - buffer_bptr) / 8) {
        buffer_bptr += count * 8;
        count = 0;
    } else {
        size_t __count = (buffer_bits - buffer_bptr) / 8;
        buffer_bptr = buffer_bits = 0;
        count -= __count;
    }
    if (count > 0) {
        fseek(fp, count, SEEK_CUR);
    }
}

bool FileStream::meet_sync_word(void)
{
    assert(buffer_bptr % 8 == 0);
    uint16_t syncword = 0;
    while ((syncword & 0xfff) != 0xfff) {
        try {
            syncword = syncword << 1 | read_bits(1);
        } catch (InputException& e) {
            return false;
        }
    }
    assert(buffer_bptr % 8 == 4);
    frame_offset = ftell(fp) - (buffer_bits - (buffer_bptr - 12)) / 8;
    return true;
}

size_t FileStream::get_frame_offset()
{
    return frame_offset;
}

void FileStream::set_offset(size_t offset)
{
    clear_buffer();
    fseek(fp, offset, SEEK_SET);
}

size_t FileStream::get_offset()
{
    assert(buffer_bptr % 8 == 0);
    return ftell(fp) - (buffer_bits - buffer_bptr) / 8;
}

void MainDataStream::check_buffer(void)
{
    if (buffer_bptr == buffer_bits) {
        clear_buffer();
        fill_buffer();
    }
}

void MainDataStream::clear_buffer(void)
{
    buffer_bits = buffer_bptr = 0;
}

void MainDataStream::fill_buffer(void)
{
    size_t count = 0;
    while (count < INPUT_STEAM_BUFFER_SIZE && !frame_list.empty()) {
        auto &frame = frame_list.front();
        size_t buffer_left = INPUT_STEAM_BUFFER_SIZE - count;
        size_t frame_left = frame.frame_size - frame_offset;
        size_t readed = std::min(buffer_left, frame_left);
        memcpy(
            buffer + count,
            frame.main_data + frame_offset,
            readed
        );
        count += readed;
        frame_offset += readed;
        if (frame_offset == frame.frame_size) {
            delete frame.main_data;
            frames_size -= frame.frame_size;
            frame_list.pop_front();
            frame_offset = 0;
        }
    }
    buffer_bits = 8 * count;
    if (buffer_bits == 0)
        throw InputException();
}

uint32_t MainDataStream::read_bits(size_t bits, size_t *limit)
{
    assert(0 <= bits && bits <= 32);
    if (limit != NULL && *limit < bits)
        throw InputException();
    else
        *limit -= bits;
    uint32_t res = 0;
    if (buffer_bptr / 8 == (buffer_bptr + bits) / 8) {
        check_buffer();
        // in a byte
        int index = buffer_bptr / 8;
        int offset = 8 - (bits + buffer_bptr) % 8;
        res = (buffer[index] >> offset) & ((1 << bits) - 1);
        buffer_bptr += bits;
    } else {
        // accross bytes
        // align up
        if (buffer_bptr % 8 != 0) {
            int index = buffer_bptr / 8;
            int readed = 8 - buffer_bptr % 8;
            res = (res << readed) | (buffer[index] & ((1 << readed) - 1));
            buffer_bptr += readed;
            bits -= readed;
        }
        // read full bytes
        while (bits >= 8) {
            check_buffer();
            int index = buffer_bptr / 8;
            res = (res << 8) | buffer[index];
            buffer_bptr += 8;
            bits -= 8; 
        }
        // read left bits
        if (bits > 0) {
            check_buffer();
            int index = buffer_bptr / 8;
            int offset = 8 - bits;
            res = (res << bits) | (buffer[index] >> offset);
            buffer_bptr += bits;
            bits -= bits;
        }
    }
    assert(buffer_bptr <= buffer_bits);
    return res;
}

void MainDataStream::insert_frame(MainDataFrame frame)
{
    frame_list.push_back(frame);
    frames_size += frame.frame_size;
}

void MainDataStream::sync_frame_data(size_t main_data_begin)
{
    size_t total_size = frames_size - frame_offset + buffer_bits / 8;
    assert(total_size >= main_data_begin);
    size_t target_offset = total_size - main_data_begin;
    if (target_offset < buffer_bits / 8) {
        buffer_bptr = target_offset * 8;
    } else {
        size_t bytes_left = target_offset - buffer_bits / 8;
        clear_buffer();
        while (bytes_left) {
            assert(!frame_list.empty());
            auto &frame = frame_list.front();
            size_t skipped = std::min(bytes_left, frame.frame_size - frame_offset);
            bytes_left -= skipped;
            frame_offset += skipped;
            if (frame_offset == frame.frame_size) {
                delete frame.main_data;
                frames_size -= frame.frame_size;
                frame_list.pop_front();
                frame_offset = 0;
            }
        }
    }
}