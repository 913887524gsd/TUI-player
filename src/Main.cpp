#include <PlayStream.h>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <pthread.h>

FILE *frame_sender, *frame_receivier;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "usage: %s filepath\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    int fd[2];
    pipe(fd);
    
    frame_receivier = fdopen(fd[0], "r");
    frame_sender = fdopen(fd[1], "w");

    pthread_t tid[2];
    pthread_create(tid + 0, NULL, MP3_fetcher, argv[1]);
    pthread_create(tid + 1, NULL, play_stream, NULL);
    pthread_join(tid[0], NULL);
    pthread_join(tid[1], NULL);
}