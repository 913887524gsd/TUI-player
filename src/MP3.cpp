#include <MP3Header.h>
#include <Huffman.h>
#include <InputStream.h>
#include <PlayStream.h>
#include <cmath>
#include <cassert>
#include <cstring>

void print_frame_header(FrameHeader *header)
{
    static int frame_counter;
    printf("Frame: %d, layer: %u crc: %u bitrate: %uk frequency: %u\n",
        ++frame_counter, 4 - header->layer, header->protection_bit, bitrate_table[header->layer][header->bitrate_index],
        sampling_frequency_table[header->sampling_frequency]);
    printf("padding: %u mode: %u mode extension: %u frame length: %u\n",
        header->padding_bit, header->mode, header->mode_extension, get_frame_length(header));
    printf("private bit: %u copyright: %u copy: %u emphasis: %u\n", 
        header->private_bit, header->copyright, header->copy, header->emphasis);
}

void print_XING_info(XingInfo *info)
{
    printf("frames: %u bytes: %u quality: %u\n", info->frames, info->bytes, info->quality);
}

FileStream *fstream;
MainDataStream *mstream;

void skip_ID3V2()
{
    ID3V2Header header;
    fstream->read_bytes(&header, sizeof(header));
    assert(memcmp(&header, ID3V2_IDENTIFIER, 3) == 0);
    fstream->skip_bytes(get_ID3V2_length(&header));
}

inline void skip_crc(FrameHeader *header)
{
    if (header->protection_bit == 0)
        fstream->skip_bytes(4);
}

inline void skip_frame_left_bytes(FrameHeader *header)
{
    fstream->skip_bytes(get_frame_length(header) - (fstream->get_offset() - fstream->get_frame_offset()));
}

bool read_frame_header(FrameHeader *header, bool read_syncword)
{
    if (read_syncword)
        header->syncword = fstream->read_bits(12);
    else
        header->syncword = 0xfff;
    header->ID = fstream->read_bits(1);
    assert(header->ID == 1); // MPEG-1
    header->layer = fstream->read_bits(2);
    assert(header->layer == layer3);
    header->protection_bit = fstream->read_bits(1);
    header->bitrate_index = fstream->read_bits(4);
    assert(0 < header->bitrate_index && header->bitrate_index < 15);
    header->sampling_frequency = fstream->read_bits(2);
    assert(header->sampling_frequency != 3);
    header->padding_bit = fstream->read_bits(1);
    header->private_bit = fstream->read_bits(1);
    header->mode = fstream->read_bits(2);
    header->mode_extension = fstream->read_bits(2);
    header->copyright = fstream->read_bits(1);
    header->copy = fstream->read_bits(1);
    header->emphasis = fstream->read_bits(2);
    return true;
}

bool read_Xing_header(XingInfo *info)
{
    char ID[4];
    uint32_t flag;
    
    fstream->read_bytes(ID, 4);
    // now just recognize CBR
    if (memcmp(ID, XING_ID_CBR, 4) != 0)
        return false;
    flag = fstream->read_bits(32);
    if (flag & XING_FRAME_FLAG) {
        info->frames = fstream->read_bits(32);
    }
    if (flag & XING_BYTES_FLAG) {
        info->bytes = fstream->read_bits(32);
    }
    if (flag & XING_TOC_FLAG) {
        // useless, just skip
        fstream->skip_bytes(100);
    }
    if (flag & XING_VBR_SCALE_FLAG) {
        info->quality = fstream->read_bits(32);
    }
    return true;
}

bool read_tag_frame(XingInfo *info)
{
    // first frame
    fstream->meet_sync_word();
    FrameHeader header;
    read_frame_header(&header, false);
    fwrite(&header, sizeof(header), 1, frame_sender);
    // skip side info
    fstream->skip_bytes(get_side_info_length(&header));
    skip_crc(&header);
    if (!read_Xing_header(info))
        return false;
    print_XING_info(info);
    // skip unused data
    skip_frame_left_bytes(&header);
    return true;
}

bool read_layer3_side_info(FrameHeader *header, Layer3SideInfo *sideinfo)
{
    int nch = header->mode == single_channel ? 1 : 2;

    sideinfo->main_data_begin = fstream->read_bits(9);
    // skip private bits;
    if (header->mode == single_channel)
        fstream->read_bits(5);
    else
        fstream->read_bits(3);
    for (int ch = 0 ; ch < nch ; ch++)
        for (int scfsi_band = 0 ; scfsi_band < 4 ; scfsi_band++)
            sideinfo->scfsi[ch][scfsi_band] = fstream->read_bits(1);
    for (int gr = 0 ; gr < 2 ; gr++)
        for (int ch = 0 ; ch < nch ; ch++) {
            Layer3ChannelInfo *cinfo = &sideinfo->cinfo[gr][ch];
            cinfo->part2_3_length = fstream->read_bits(12);
            cinfo->big_values = fstream->read_bits(9);
            cinfo->global_gain = fstream->read_bits(8);
            cinfo->scalefac_compress = fstream->read_bits(4);
            cinfo->window_switching_flag = fstream->read_bits(1);
            if (cinfo->window_switching_flag) {
                cinfo->block_type = fstream->read_bits(2);
                assert(cinfo->block_type != 0);
                cinfo->mixed_block_flag = fstream->read_bits(1);
                for (int region = 0 ; region < 2 ; region++)
                    cinfo->table_select[region] = fstream->read_bits(5);
                for (int window = 0 ; window < 3 ; window++)
                    cinfo->subblock_gain[window] = fstream->read_bits(3);
                if (cinfo->block_type == 2 && cinfo->mixed_block_flag == 0)
                    cinfo->region0_count = 8;
                else
                    cinfo->region0_count = 7;
                cinfo->region1_count = 36;
            } else {
                cinfo->block_type = 0;
                for (int region = 0 ; region < 3 ; region++)
                    cinfo->table_select[region] = fstream->read_bits(5);
                cinfo->region0_count = fstream->read_bits(4);
                cinfo->region1_count = fstream->read_bits(3);
            }
            cinfo->preflag = fstream->read_bits(1);
            cinfo->scalefac_scale = fstream->read_bits(1);
            cinfo->count1table_select = fstream->read_bits(1);
        }
    return true;
}

bool read_layer3_main_data_frame(FrameHeader *header)
{
    MainDataFrame frame;
    frame.frame_size = get_frame_length(header) - (fstream->get_offset() - fstream->get_frame_offset());
    frame.main_data = new uint8_t[frame.frame_size];
    fstream->read_bytes(frame.main_data, frame.frame_size);
    mstream->insert_frame(frame);
    return true;
}

code read_layer3_big_value_huffman(int index, size_t *limit)
{
    static int16_t buffer[2];
    huffman_code *huffman = get_code_huffman(index);
    const code* code = huffman->decode(mstream, limit);
    buffer[0] = code->x;
    buffer[1] = code->y;
    for (int i = 0 ; i < 2 ; i++) {
        if (buffer[i] == 15 && huffman->linbits > 0)
            buffer[i] += mstream->read_bits(huffman->linbits, limit);
        if (buffer[i] != 0 && mstream->read_bits(1, limit))
            buffer[i] = -buffer[i];
    }
    return {buffer[0], buffer[1]};
}

quadruples read_layer3_count1_huffman(int index, size_t *limit)
{
    static int16_t buffer[4];
    huffman_quadruples *huffman = get_quadruple_huffman(index);
    const quadruples* quadruples = huffman->decode(mstream, limit);
    buffer[0] = quadruples->v;
    buffer[1] = quadruples->w;
    buffer[2] = quadruples->x;
    buffer[3] = quadruples->y;
    for (int i = 0 ; i < 4 ; i++)
        if (buffer[i] != 0 && mstream->read_bits(1, limit))
            buffer[i] = -buffer[i];
    return {buffer[0], buffer[1], buffer[2], buffer[3]};
}

bool read_layer3_huffman_code(FrameHeader *header, Layer3ChannelInfo *cinfo, Layer3ChannelMainData *cdata, size_t *limit)
{
    const int *limit_table = scale_band_long_acm_table[header->sampling_frequency];
    int big_values = cinfo->big_values;
    int region_limit[2];
    int l = 0;
    
    if (cinfo->block_type == 2) {
        region_limit[0] = 36; // this can be calculated
        region_limit[1] = 576;
    } else {
        region_limit[0] = std::min(cinfo->region0_count + 1, 22);
        region_limit[1] = std::min(region_limit[0] + cinfo->region1_count + 1, 22);
        region_limit[0] = limit_table[region_limit[0]];
        region_limit[1] = limit_table[region_limit[1]];
    }

    for (; l < big_values * 2 ; l += 2) {
        int index;
        if (l < region_limit[0]) {
            index = cinfo->table_select[0];
        } else if (l < region_limit[1]) {
            index = cinfo->table_select[1];
        } else {
            index = cinfo->table_select[2];
        }
        auto tmp = read_layer3_big_value_huffman(index, limit);
        cdata->samples[l + 0] = tmp.x;
        cdata->samples[l + 1] = tmp.y;
    }
    
    while (*limit && l < 576) {
        int index = cinfo->count1table_select;
        auto tmp = read_layer3_count1_huffman(index, limit);
        cdata->samples[l + 0] = tmp.v;
        cdata->samples[l + 1] = tmp.w;
        cdata->samples[l + 2] = tmp.x;
        cdata->samples[l + 3] = tmp.y;
        l += 4;
    }

    while (l < 576)
        cdata->samples[l++] = 0;

    return true;
}

bool read_layer3_main_data(FrameHeader *header, Layer3SideInfo *sideinfo, Layer3MainData *maindata)
{
    int nch = header->mode == single_channel ? 1 : 2;
    for (int gr = 0 ; gr < 2 ; gr++)
        for (int ch = 0 ; ch < nch ; ch++) {
            Layer3ChannelInfo *cinfo = &sideinfo->cinfo[gr][ch];
            Layer3ChannelMainData *cdata = &maindata->cdata[gr][ch];
            size_t limit = cinfo->part2_3_length;
            size_t slen1 = slen1_table[cinfo->scalefac_compress];
            size_t slen2 = slen2_table[cinfo->scalefac_compress];
            if (cinfo->block_type == 2) {
                if (cinfo->mixed_block_flag) {
                    for (int sfb = 0 ; sfb < 8 ; sfb++)
                        cdata->scalefac_l[sfb] = mstream->read_bits(slen1, &limit);
                    for (int sfb = 3 ; sfb < 6 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++)
                            cdata->scalefac_s[sfb][window] = mstream->read_bits(slen1, &limit);
                    for (int sfb = 6 ; sfb < 12 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++)
                            cdata->scalefac_s[sfb][window] = mstream->read_bits(slen2, &limit);
                } else {
                    for (int sfb = 0 ; sfb < 6 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++)
                            cdata->scalefac_s[sfb][window] = mstream->read_bits(slen1, &limit);
                    for (int sfb = 6 ; sfb < 12 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++)
                            cdata->scalefac_s[sfb][window] = mstream->read_bits(slen2, &limit);
                }
            } else {
                for (int sfb = 0 ; sfb < 6 ; sfb++)
                    if (sideinfo->scfsi[ch][0] == 0 || gr == 0)
                        cdata->scalefac_l[sfb] = mstream->read_bits(slen1, &limit);
                    else
                        cdata->scalefac_l[sfb] = maindata->cdata[0][ch].scalefac_l[sfb];
                for (int sfb = 6 ; sfb < 11 ; sfb++)
                    if (sideinfo->scfsi[ch][1] == 0 || gr == 0)
                        cdata->scalefac_l[sfb] = mstream->read_bits(slen1, &limit);
                    else
                        cdata->scalefac_l[sfb] = maindata->cdata[0][ch].scalefac_l[sfb];
                for (int sfb = 11 ; sfb < 16 ; sfb++)
                    if (sideinfo->scfsi[ch][2] == 0 || gr == 0)
                        cdata->scalefac_l[sfb] = mstream->read_bits(slen2, &limit);
                    else
                        cdata->scalefac_l[sfb] = maindata->cdata[0][ch].scalefac_l[sfb];
                for (int sfb = 16 ; sfb < 21 ; sfb++)
                    if (sideinfo->scfsi[ch][3] == 0 || gr == 0)
                        cdata->scalefac_l[sfb] = mstream->read_bits(slen2, &limit);
                    else
                        cdata->scalefac_l[sfb] = maindata->cdata[0][ch].scalefac_l[sfb];
            }
            if (cinfo->block_type == 2) {
                if (cinfo->mixed_block_flag)
                    assert(cinfo->part2_3_length - limit == 17 * slen1 + 18 * slen2);
                else
                    assert(cinfo->part2_3_length - limit == 18 * slen1 + 18 * slen2);
            } else {
                size_t checker = 11 * slen1 + 10 * slen2;
                for (int scfsi_band = 0 ; scfsi_band < 4 ; scfsi_band++)
                    if (sideinfo->scfsi[ch][scfsi_band] && gr == 1) {
                        switch (scfsi_band) {
                            case 0: checker -= slen1 * 6; break;
                            case 1: checker -= slen1 * 5; break;
                            case 2: checker -= slen2 * 5; break;
                            case 3: checker -= slen2 * 5; break;
                            default:break;
                        }
                    }
                assert(cinfo->part2_3_length - limit == checker);
            }
            read_layer3_huffman_code(header, cinfo, cdata, &limit);
        }
    return true;
}

#define PI 3.141592653589793
#define SQRT2 1.414213562373095

void requantization(FrameHeader *header, Layer3SideInfo *sideinfo, Layer3MainData *maindata)
{
    static const float pow_4_3[17] = {
        0.0, 1.0, 2.5198420997897464, 4.3267487109222245,
        6.3496042078727974, 8.549879733383484, 10.902723556992836, 13.390518279406722,
        16.0, 18.720754407467133, 21.544346900318832, 24.463780996262468,
        27.47314182127996, 30.567350940369842, 33.74199169845321, 36.993181114957046,
        40.317473596635935
    };
    int nch = header->mode == single_channel ? 1 : 2;
    const int *long_table = scale_band_long_table[header->sampling_frequency];
    const int *short_table = scale_band_short_table[header->sampling_frequency];
    
    for (int gr = 0 ; gr < 2 ; gr++)
        for (int ch = 0 ; ch < nch ; ch++) {
            Layer3ChannelInfo *cinfo = &sideinfo->cinfo[gr][ch];
            Layer3ChannelMainData *cdata = &maindata->cdata[gr][ch];
            const float scalefac_multiplier = cinfo->scalefac_scale ? 1.0 : 0.5;
            
            int fl = 0;
            if (cinfo->block_type == 2) {
                if (cinfo->mixed_block_flag) {
                    for (int sfb = 0 ; sfb < 8 ; sfb++)
                        for (int i = 0 ; i < long_table[sfb] ; i++, fl++) {
                            if(cdata->samples[fl] == 0.0)
                                continue;
                            float v = cdata->samples[fl];
                            float absv = abs(v);
                            cdata->samples[fl] = (v > 0 ? 1.0 : -1.0) * // sign
                                ((int)absv > 16 ? pow(absv, 4.0 / 3.0) : pow_4_3[(int)absv]) * // is
                                pow(2, (cinfo->global_gain - 210) / 4.0) *
                                pow(2, -(scalefac_multiplier * (cdata->scalefac_l[sfb] + cinfo->preflag * pretab[sfb])));
                        }
                    for (int sfb = 3 ; sfb < 12 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++)
                            for (int i = 0 ; i < short_table[sfb] ; i++, fl++) {
                                if (cdata->samples[fl] == 0.0)
                                    continue;
                                float v = cdata->samples[fl];
                                float absv = abs(v);
                                cdata->samples[fl] = (v > 0 ? 1.0 : -1.0) * // sign
                                    ((int)absv > 16 ? pow(absv, 4.0 / 3.0) : pow_4_3[(int)absv]) * // is
                                    pow(2, (cinfo->global_gain - 210 - 8 * cinfo->subblock_gain[window]) / 4.0) *
                                    pow(2, -(scalefac_multiplier * cdata->scalefac_s[sfb][window]));
                            }
                    while (fl < 576)
                        cdata->samples[fl++] = 0;
                } else {
                    for (int sfb = 0 ; sfb < 12 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++)
                            for (int i = 0 ; i < short_table[sfb] ; i++, fl++) {
                                if (cdata->samples[fl] == 0.0)
                                    continue;
                                float v = cdata->samples[fl];
                                float absv = abs(v);
                                cdata->samples[fl] = (v > 0 ? 1.0 : -1.0) * // sign
                                    ((int)absv > 16 ? pow(absv, 4.0 / 3.0) : pow_4_3[(int)absv]) * // is
                                    pow(2, (cinfo->global_gain - 210 - 8 * cinfo->subblock_gain[window]) / 4.0) *
                                    pow(2, -(scalefac_multiplier * cdata->scalefac_s[sfb][window]));
                            }
                    while (fl < 576)
                        cdata->samples[fl++] = 0;
                }
            } else {
                for (int sfb = 0 ; sfb < 21 ; sfb++)
                    for (int i = 0 ; i < long_table[sfb] ; i++, fl++) {
                        if(cdata->samples[fl] == 0.0)
                            continue;
                        float v = cdata->samples[fl];
                        float absv = abs(v);
                        cdata->samples[fl] = (v > 0 ? 1.0 : -1.0) * // sign
                            ((int)absv > 16 ? pow(absv, 4.0 / 3.0) : pow_4_3[(int)absv]) * // is
                            pow(2, (cinfo->global_gain - 210) / 4.0) *
                            pow(2, -(scalefac_multiplier * (cdata->scalefac_l[sfb] + cinfo->preflag * pretab[sfb])));
                    }
                while (fl < 576)
                    cdata->samples[fl++] = 0;
            }
        }
}

void stereo_processing(FrameHeader *header, Layer3SideInfo *sideinfo, Layer3MainData *maindata)
{
    if (header->mode != joint_stereo || header->mode_extension == 0)
        return;
    const int *long_table = scale_band_long_table[header->sampling_frequency];
    const int *short_table = scale_band_short_table[header->sampling_frequency];
    for (int gr = 0 ; gr < 2 ; gr++) {
        Layer3ChannelInfo *cinfos = sideinfo->cinfo[gr];
        Layer3ChannelMainData *cdatas = maindata->cdata[gr];
        assert(cinfos[0].block_type == cinfos[1].block_type);
        if (cinfos[0].block_type == 2)
            assert(cinfos[0].mixed_block_flag == cinfos[1].mixed_block_flag);
        
        auto MS_Stereo = [&](int limit) {
            for (int i = 0 ; i < limit ; i++) {
                float M = cdatas[0].samples[i];
                float S = cdatas[1].samples[i];
                cdatas[0].samples[i] = (M + S) / SQRT2;
                cdatas[1].samples[i] = (M - S) / SQRT2;
            }
        };

        auto Intensity_Stereo = [&](int fl, int count, uint8_t is_pos) {
            if (is_pos == 7)
                return;
            float is_ratio = tan(is_pos * PI / 12.0);
            for (; count > 0 ; fl++, count--) {
                float base = cdatas[1].samples[fl];
                cdatas[0].samples[fl] = base * (is_ratio / (1 + is_ratio));
                cdatas[1].samples[fl] = base * (1 / (1 + is_ratio));
            }
        };

        if (header->mode_extension == 2) {
            MS_Stereo(576);
        } else {
            int all_zero = 575;
            while (cdatas[1].samples[all_zero] == 0.0)
                all_zero--;
            all_zero++;
            int fl = 0, align_up = 0;
            if (cinfos[0].block_type == 2) {
                if (cinfos[0].mixed_block_flag) {
                    for (int sfb = 0 ; sfb < 8 ; sfb++) {
                        if (fl >= all_zero) {
                            Intensity_Stereo(fl, long_table[sfb], cdatas[1].scalefac_l[sfb]);
                            if (!align_up)
                                align_up = 1, all_zero = fl;
                        }
                        fl += long_table[sfb];
                    }
                    for (int sfb = 3 ; sfb < 12 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++) {
                            if (fl >= all_zero) {
                                Intensity_Stereo(fl, short_table[sfb], cdatas[1].scalefac_s[sfb][window]);
                                if (!align_up)
                                    align_up = 1, all_zero = fl;
                            }
                            fl += short_table[sfb];
                        }
                } else {
                    for (int sfb = 0 ; sfb < 12 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++) {
                            if (fl >= all_zero) {
                                Intensity_Stereo(fl, short_table[sfb], cdatas[1].scalefac_s[sfb][window]);
                                if (!align_up)
                                    align_up = 1, all_zero = fl;
                            }
                            fl += short_table[sfb];
                        }
                }
            } else {
                for (int sfb = 0 ; sfb < 21 ; sfb++) {
                    if (fl >= all_zero) {
                        Intensity_Stereo(fl, long_table[sfb], cdatas[1].scalefac_l[sfb]);
                        if (!align_up)
                            align_up = 1, all_zero = fl;
                    }
                    fl += long_table[sfb];
                }
            }
            if (header->mode_extension == 3) {
                if (align_up) 
                    MS_Stereo(all_zero);
                else
                    MS_Stereo(576);// I don't know how to handle
            }
        }
    }
}

void reorder(FrameHeader *header, Layer3SideInfo *sideinfo, Layer3MainData *maindata, Layer3XR xr[2][2])
{
    int nch = header->mode == single_channel ? 1 : 2;
    const int *long_table = scale_band_long_table[header->sampling_frequency];
    const int *short_table = scale_band_short_table[header->sampling_frequency];
    
    for (int gr = 0 ; gr < 2 ; gr++)
        for (int ch = 0 ; ch < nch ; ch++) {
            Layer3ChannelInfo *cinfo = &sideinfo->cinfo[gr][ch];
            Layer3ChannelMainData *cdata = &maindata->cdata[gr][ch];
            Layer3XR *cxr = &xr[gr][ch];
            
            int fl = 0;
            if (cinfo->block_type == 2) {
                if (cinfo->mixed_block_flag) {
                    int l = 0;
                    int s[3] = {};
                    for (int sfb = 0 ; sfb < 8 ; sfb++)
                        for (int i = 0 ; i < long_table[sfb] ; i++, l++)
                            cxr->l[l / 18][l % 18] = cdata->samples[fl++];
                    for (int sfb = 3; sfb < 12 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++)
                            for (int i = 0 ; i < short_table[sfb] ; i++, s[window]++)
                                cxr->s[s[window] / 6][window][s[window] % 6] = cdata->samples[fl++];
                    for (int window = 0 ; window < 3 ; window++)
                        for (; s[window] < 180 ; s[window]++)
                            cxr->s[s[window] / 6][window][s[window] % 6] = 0;
                } else {
                    int s[3] = {};
                    for (int sfb = 0; sfb < 12 ; sfb++)
                        for (int window = 0 ; window < 3 ; window++)
                            for (int i = 0 ; i < short_table[sfb] ; i++, s[window]++)
                                cxr->s[s[window] / 6][window][s[window] % 6] = cdata->samples[fl++];
                    for (int window = 0 ; window < 3 ; window++)
                        for (; s[window] < 192 ; s[window]++)
                            cxr->s[s[window] / 6][window][s[window] % 6] = 0;
                }
            } else {
                int l = 0;
                for (int sfb = 0 ; sfb < 21 ; sfb++)
                    for (int i = 0 ; i < long_table[sfb] ; i++, l++)
                        cxr->l[l / 18][l % 18] = cdata->samples[fl++];
                for (; l < 576 ; l++)
                    cxr->l[l / 18][l % 18] = 0;
            }
        }
}

void alias_reduction(FrameHeader *header, Layer3SideInfo *sideinfo, Layer3XR xr[2][2])
{
    const static float cs[8] = {
        0.8574929257125443, 0.8817419973177052, 0.9496286491027328, 0.9833145924917902,
        0.9955178160675858, 0.9991605581781475, 0.9998991952444471, 0.9999931550702803,
    };
    const static float ca[8] = {
        -0.5144957554275266, -0.47173196856497235, -0.31337745420390184, -0.18191319961098118,
        -0.09457419252642066, -0.04096558288530405, -0.01419856857247115, -0.0036999746737600373,
    };
    int nch = header->mode == single_channel ? 1 : 2;
    for (int gr = 0 ; gr < 2 ; gr++)
        for (int ch = 0 ; ch < nch ; ch++) {
            Layer3ChannelInfo *cinfo = &sideinfo->cinfo[gr][ch];
            Layer3XR *cxr = &xr[gr][ch];
            int sb_max;
            if (cinfo->block_type != 2)
                sb_max = 32;
            else if (cinfo->mixed_block_flag)
                sb_max = 2;
            else
                sb_max = 0;
            for (int sb = 1 ; sb < sb_max ; sb++)
                for (int i = 0 ; i < 8 ; i++) {
                    int offset_x = sb * 18 - 1 - i;
                    int offset_y = sb * 18 + i;
                    float x = cxr->pcm[offset_x];
                    float y = cxr->pcm[offset_y];
                    cxr->pcm[offset_x] = x * cs[i] - y * ca[i];
                    cxr->pcm[offset_y] = y * cs[i] + x * ca[i];
                }
        }
}

// 0 -> n = 36
// 1 -> n = 12
static float IMDCT_cos[2][36][18];
static float WIN_0_sin[36];
static float WIN_1_sin[36];
static float WIN_2_sin[12];
static float WIN_3_sin[36];

__attribute__((constructor))
static void init_IMDCT_const()
{
    for (int i = 0 ; i < 36 ; i++)
        for (int k = 0 ; k < 18 ; k++)
            IMDCT_cos[0][i][k] = cos(PI / 72 * (2 * i + 1 + 18) * (2 * k + 1));
    for (int i = 0 ; i < 12 ; i++)
        for (int k = 0 ; k < 6 ; k++)
            IMDCT_cos[1][i][k] = cos(PI / 24 * (2 * i + 1 + 6) * (2 * k + 1));
    for (int i = 0 ; i < 36 ; i++)
        WIN_0_sin[i] = sin(PI / 36 * (i + 0.5));
    for (int i = 0 ; i < 18 ; i++)
        WIN_1_sin[i] = sin(PI / 36 * (i + 0.5));
    for (int i = 24 ; i < 30 ; i++)
        WIN_1_sin[i] = sin(PI / 12 * (i - 18 + 0.5));
    for (int i = 0 ; i < 12 ; i++)
        WIN_2_sin[i] = sin(PI / 12 * (i + 0.5));
    for (int i = 6 ; i < 12 ; i++)
        WIN_3_sin[i] = sin(PI / 12 * (i - 6 + 0.5));
    for (int i = 18 ; i < 36 ; i++)
        WIN_3_sin[i] = sin(PI / 36 * (i + 0.5));
}

void IMDCT(FrameHeader *header, Layer3SideInfo *sideinfo, Layer3XR xr[2][2])
{
    int nch = header->mode == single_channel ? 1 : 2;
    static union {
        float l[36];
        float s[3][12];
    } buffer, buffer2;
    static float prev_frame[2][32][18];
    for (int gr = 0 ; gr < 2 ; gr++)
        for (int ch = 0 ; ch < nch ; ch++) {
            Layer3ChannelInfo *cinfo = &sideinfo->cinfo[gr][ch];
            Layer3XR *cxr = &xr[gr][ch];
            
            for (int sb = 0 ; sb < 32 ; sb++) {
                int short_block = cinfo->block_type == 2 && (sb > 2 || !cinfo->mixed_block_flag);
                int n = short_block ? 12 : 36;
                int block_type = cinfo->block_type;
                if (cinfo->window_switching_flag && cinfo->mixed_block_flag && sb < 2)
                    block_type = 0;
                // IMDCT
                memset(&buffer, 0, sizeof(buffer));
                if (short_block) {
                    for (int i = 0 ; i < n ; i++)
                        for (int j = 0 ; j < 3 ; j++)
                            for (int k = 0 ; k < n / 2 ; k++)
                                buffer.s[j][i] += cxr->s[sb][j][k] * IMDCT_cos[1][i][k];
                } else {
                    for (int i = 0 ; i < n ; i++)
                        for (int k = 0 ; k < n / 2 ; k++)
                            buffer.l[i] += cxr->l[sb][k] * IMDCT_cos[0][i][k];
                }
                // window
                if (block_type == 0) {
                    assert(!short_block);
                    for (int i = 0 ; i < 36 ; i++)
                        buffer.l[i] *= WIN_0_sin[i];
                } else if (block_type == 1) {
                    assert(!short_block);
                    for (int i = 0 ; i < 18 ; i++)
                        buffer.l[i] *= WIN_1_sin[i];
                    for (int i = 24 ; i < 30 ; i++)
                        buffer.l[i] *= WIN_1_sin[i];
                    for (int i = 30 ; i < 36 ; i++)
                        buffer.l[i] = 0;
                } else if (block_type == 2) {
                    assert(short_block);
                    for (int j = 0 ; j < 3; j++)
                        for (int i = 0 ; i < 12 ; i++)
                            buffer2.s[j][i] = buffer.s[j][i] * WIN_2_sin[i];
                    for (int i = 0 ; i < 6 ; i++)
                        buffer.l[i] = 0;
                    for (int i = 6 ; i < 12 ; i++)
                        buffer.l[i] = buffer2.s[0][i - 6];
                    for (int i = 12 ; i < 18 ; i++)
                        buffer.l[i] = buffer2.s[0][i - 6] + buffer2.s[1][i - 12];
                    for (int i = 18 ; i < 24 ; i++)
                        buffer.l[i] = buffer2.s[1][i - 12] + buffer2.s[2][i - 18];
                    for (int i = 24 ; i < 30 ; i++)
                        buffer.l[i] = buffer2.s[2][i - 18];
                    for (int i = 30 ; i < 36 ; i++)
                        buffer.l[i] = 0; 
                } else {
                    assert(!short_block);
                    for (int i = 0 ; i < 6 ; i++)
                        buffer.l[i] = 0;
                    for (int i = 6 ; i < 12 ; i++)
                        buffer.l[i] *= WIN_3_sin[i];
                    for (int i = 18 ; i < 36 ; i++)
                        buffer.l[i] *= WIN_3_sin[i];
                }
                // overlap
                for (int i = 0 ; i < 18 ; i++) {
                    cxr->l[sb][i] = buffer.l[i] + prev_frame[ch][sb][i];
                    prev_frame[ch][sb][i] = buffer.l[i + 18];
                }
            }
        }
}

void frequency_inversion(FrameHeader *header, Layer3XR xr[2][2])
{
    int nch = header->mode == single_channel ? 1 : 2;
    static float buffer[32][18];
    for (int gr = 0 ; gr < 2 ; gr++)
        for (int ch = 0 ; ch < nch ; ch++) {
            for (int sb = 0 ; sb < 32 ; sb++)
                for (int i = 0 ; i < 18 ; i++)
                    buffer[sb][i] = xr[gr][ch].l[sb][i];
            for (int sb = 1 ; sb < 32 ; sb += 2)
                for (int i = 1 ; i < 18 ; i += 2)
                    buffer[sb][i] *= -1;
            for (int i = 0 ; i < 18 ; i++)
                for (int sb = 0 ; sb < 32 ; sb++)
                    xr[gr][ch].r[i][sb] = buffer[sb][i];
        }
}

static float synth_n[64][32];

__attribute__((constructor))
static void init_synth_filterbank_const()
{
    for (int k = 0 ; k < 64 ; k++)
        for (int i = 0 ; i < 32 ; i++)
            synth_n[k][i] = cos((16.0 + k) * (2 * i + 1) * PI / 64);
}

void synth_filterbank(FrameHeader *header, Layer3XR xr[2][2])
{
    int nch = header->mode == single_channel ? 1 : 2;

    static float v[2][1024];
    float s[32], u[512], w[512];
    
    for (int gr = 0 ; gr < 2 ; gr++)
        for (int ch = 0 ; ch < nch ; ch++) {
            Layer3XR *cxr = &xr[gr][ch];
            for (int t = 0 ; t < 18 ; t++) {
                for (int i = 1023 ; i >= 64 ; i--)
                    v[ch][i] = v[ch][i - 64];
                for (int i = 0 ; i < 64 ; i++) {
                    v[ch][i] = 0;
                    for (int k = 0 ; k < 32 ; k++)
                        v[ch][i] += synth_n[i][k] * cxr->r[t][k];
                }
                for (int i = 0 ; i < 8 ; i++)
                    for (int j = 0 ; j < 32 ; j++) {
                        u[i * 64 + j] = v[ch][i * 128 + j];
                        u[i * 64 + 32 + j] = v[ch][i * 128 + 96 + j];
                    }
                for (int i = 0 ; i < 512 ; i++)
                    w[i] = u[i] * synth_D[i];
                for (int j = 0 ; j < 32 ; j++) {
                    s[j] = 0;
                    for (int i = 0 ; i < 16 ; i++)
                        s[j] += w[j + 32 * i];
                    cxr->r[t][j] = s[j];
                }
            }
        }
}

void orchestrate(FrameHeader *header, Layer3XR xr[2][2], float *pcm)
{
    int nch = header->mode == single_channel ? 1 : 2;
    int pcm_p = 0;
    
    for (int gr = 0 ; gr < 2 ; gr++)
        for (int i = 0 ; i < 576 ; i++)
            for (int ch = 0 ; ch < nch ; ch++)
                pcm[pcm_p++] = xr[gr][ch].pcm[i];
}

void *MP3_fetcher(void *filepath)
{
    fstream = new FileStream((const char *)filepath);
    mstream = new MainDataStream();
    FrameHeader *header = new FrameHeader;
    Layer3SideInfo *sideinfo = new Layer3SideInfo;
    Layer3MainData *maindata = new Layer3MainData;
    auto xr = new Layer3XR[2][2];
    float *pcm = new float[576 * 4];
    skip_ID3V2();
    XingInfo XING_info;
    if (read_tag_frame(&XING_info) == false) {
        fprintf(stderr, "read_tag_frame error");
        exit(EXIT_FAILURE);
    }
    for (size_t frame = 1 ; frame <= XING_info.frames && fstream->meet_sync_word() ; frame++) {
        read_frame_header(header, false);
        print_frame_header(header);
        skip_crc(header);
        read_layer3_side_info(header, sideinfo);
        mstream->sync_frame_data(sideinfo->main_data_begin);
        read_layer3_main_data_frame(header);
        read_layer3_main_data(header, sideinfo, maindata);
        requantization(header, sideinfo, maindata);
        stereo_processing(header, sideinfo, maindata);
        reorder(header, sideinfo, maindata, xr);
        alias_reduction(header, sideinfo, xr);
        IMDCT(header, sideinfo, xr);
        frequency_inversion(header, xr);
        synth_filterbank(header, xr);
        orchestrate(header, xr, pcm);
        int nch = header->mode == single_channel ? 1 : 2;
        fwrite(pcm, 2 * nch * 576 * sizeof(float) , 1, frame_sender);
    }
    fclose(frame_sender);
    delete fstream; delete mstream;
    delete header; delete sideinfo; delete maindata; delete xr; delete pcm;
    return NULL;
}