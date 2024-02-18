#include <PlayStream.h>
#include <MP3Header.h>
#include <cstdlib>
#include <alsa/asoundlib.h>

void *play_stream(void * __attribute__((unused))arg)
{
    FrameHeader header;

    snd_pcm_t *handle;
    snd_pcm_hw_params_t *hw = NULL;
    
    if (snd_pcm_open(&handle, "default", SND_PCM_STREAM_PLAYBACK, 0) < 0) {
        fprintf(stderr, "snd_pcm_open failed");
        exit(EXIT_FAILURE);
    }
    snd_pcm_hw_params_alloca(&hw);
    snd_pcm_hw_params_any(handle, hw);
    if (snd_pcm_hw_params_set_access(handle, hw, SND_PCM_ACCESS_RW_INTERLEAVED) < 0) {
        fprintf(stderr, "snd_pcm_hw_params_set_access failed");
        exit(EXIT_FAILURE);
    }
    if (snd_pcm_hw_params_set_format(handle, hw, SND_PCM_FORMAT_FLOAT_LE) < 0) {
        fprintf(stderr, "snd_pcm_hw_params_set_format failed");
        exit(EXIT_FAILURE);
    }
    if (fread(&header, sizeof(header), 1, frame_receivier) == 0) {
        fprintf(stderr, "fread failed");
        exit(EXIT_FAILURE);
    }
    unsigned nch = header.mode == single_channel ? 1 : 2;
    if (snd_pcm_hw_params_set_channels(handle, hw, nch) < 0) {
        fprintf(stderr, "snd_pcm_hw_params_set_channels failed");
        exit(EXIT_FAILURE);
    }
    unsigned sample_rate = sampling_frequency_table[header.sampling_frequency];
    if (snd_pcm_hw_params_set_rate_near(handle, hw, &sample_rate, NULL) < 0) {
        fprintf(stderr, "snd_pcm_hw_params_set_rate_near failed");
        exit(EXIT_FAILURE);
    }
    snd_pcm_uframes_t frames = 4;
    if (snd_pcm_hw_params_set_period_size_near(handle, hw, &frames, NULL) < 0) {
        fprintf(stderr, "snd_pcm_hw_params_set_period_size_near failed");
        exit(EXIT_FAILURE);
    }
    if (snd_pcm_hw_params(handle, hw) < 0) {
        fprintf(stderr, "snd_pcm_hw_params failed");
        exit(EXIT_FAILURE);
    }

    float pcm[576 * 4];

    while (fread(pcm, 2 * nch * 576 * sizeof(float), 1, frame_receivier) == 1) {
        int e = snd_pcm_writei(handle, pcm, 1152);
		if (e == -EPIPE)
			snd_pcm_recover(handle, e, 0);
    }
    fclose(frame_receivier);

    snd_pcm_drain(handle);
    snd_pcm_close(handle);
    return NULL;
}