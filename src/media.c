#include "../headers/media.h"



void save_snapshot(int *LATTICE, uint8_t* tmp,  int n, int i) {
    for(int i = 0; i < n*n; i++) { 
        if(LATTICE[i] == -1){ 
            tmp[i] = 0;
        }
        else { 
            tmp[i] = 255;
        }
    }
    char buffer[32]; // The filename buffer.
    snprintf(buffer, sizeof(char) * 32, "../output/snapshots/file%i.pgm", i);
    FILE *fp = fopen(buffer, "wb"); /* b - binary mode */
    fprintf(fp, "P5\n%d %d\n255\n", n, n); 
    fwrite(tmp, 1, n*n, fp); // n*n pixel values, each 1 byte long.
    fclose(fp);
}


void save_video(int* simspace, int n, int nsteps) { 
    int n2 = n*n;
    int k = 0;
    uchar *frame = malloc(sizeof(uchar) * n2 * 3);
    FILE *pipeout = popen("ffmpeg -y -f rawvideo -vcodec rawvideo -pix_fmt rgb24 -s 100x100 -r 25 -i - -f mp4 -q:v 5 -an -vcodec mpeg4 ../output/output.mp4", "w");     
    for(int i = 0; i < nsteps; i++) {
        int *simsspace_at_i = simspace + i*n2;
        for(int j = 0; j < n*n; j++) { 
            k = j * 3;
            if (simsspace_at_i[j] == -1) {
                frame[k + 0] = frame[k + 1] = frame[k + 2] = 0;
                // frame[j + 1] = 0;
                // frame[j + 2] = 0;
            } 
            else { 
                frame[k + 0] = frame[k + 1] = frame[k + 2] = 250;
                // frame[j + 1] = 255;
                // frame[j + 2] = 255;
            }
        }
        fwrite(frame, 1, n*n*3, pipeout);
    }
    fflush(pipeout);
    pclose(pipeout);
    free(frame); 
}
