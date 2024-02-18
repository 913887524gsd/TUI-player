#pragma once

#include <cstdio>

extern FILE *frame_sender, *frame_receivier;

void *MP3_fetcher(void *filepath);
void *play_stream(void *arg);