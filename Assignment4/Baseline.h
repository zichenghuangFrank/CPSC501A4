//
//  Baseline.h
//  
//
//  Created by 黄子宬 on 2019/4/2.
//




#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct header_file
{
    char chunk_id[4];
    int chunk_size;
    char format[4];
    char subchunk1_id[4];
    int subchunk1_size;
    short int audio_format;
    short int num_channels;
    int sample_rate;            // sample_rate denotes the sampling rate.
    int byte_rate;
    short int block_align;
    short int bits_per_sample;
    char subchunk2_id[4];
    int subchunk2_size;            // subchunk2_size denotes the number of samples.
} header;

typedef struct header_file* header_p;

void wavRead()
{
    FILE * infile = fopen("GuitarDry.wav","rb");            // Open wave file in read mode
    header_p meta = (header_p)malloc(sizeof(header));    // header_p points to a header struct that contains the wave file metadata fields
    
    fread(meta, sizeof(header), 1, infile);
    printf(" Size of Header file is %d  bytes\n", sizeof(*meta));
    printf(" Sampling rate of the input wave file is %d Hz \n", meta->sample_rate);
    printf(" Bits per sample in wave file is %d \n", meta->bits_per_sample );
    int i=0;
    short int value=0;
    int input_signal[meta->subchunk2_size];
    if(strcmp(meta->subchunk2_id, "data")){
        while( fread(&value,sizeof(value),1,infile)>0)
        {
            input_signal[i] = value;
            printf("%d\n", input_signal[i]);
            i++;
        }
        fclose(infile);
        
    }
    return 0;
}
