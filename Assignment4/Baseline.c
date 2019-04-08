//
//  Baseline.c
//  
//
//  Created by 黄子宬 on 2019/4/2.
//

/*  HEADER FILES  ************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/*  CONSTANTS  ***************************************************************/
#define PI                3.14159265358979

/*  Test tone frequency in Hz  */
#define FREQUENCY         440.0

/*  Test tone duration in seconds  */
#define DURATION          2.0

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE   16

/*  Standard sample size in bytes  */
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC        1
#define STEREOPHONIC      2


/*  FUNCTION PROTOTYPES  *****************************************************/
void createTestTone(double frequency, double duration,
                    int numberOfChannels, char *filename);
void writeWaveFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile);

size_t fwriteIntLSB(int data, FILE *stream);

size_t fwriteShortLSB(short int data, FILE *stream);

void convolve(float x[], int N, float h[], int M, float y[], int P);

void writeWAV(int channels, int numberSamples, double outputRate, char *filename, int16_t *samples);

float intToFloat(int16_t x);
int16_t floatToInt(float x);

struct header_file
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
    long int subchunk2_size;            // subchunk2_size denotes the number of samples.
};
struct header_file meta;

int main(int argc, char **argv){
//
//    if(argc != 4){
//        printf("This program needs 3 command arguments.\n");
//        exit(-1);
//    }
    
    FILE *file1;
    short int *buffer1;
    unsigned long fileLen1;
    struct header_file meta1;
    //    Open file
    file1 = fopen("GuitarDry.wav", "rb");
    if (!file1)
    {
        fprintf(stderr, "Unable to open file %s", "GuitarDry.wav");
        return 0;
    }
    //    Read header
    fread(&meta1, sizeof(meta1), 1, file1);
    printf(" Size of Header file : %lu  bytes\n", sizeof(struct header_file));
    printf(" Sampling rate of the input wave file : %d Hz \n", meta1.sample_rate);
    printf(" Bits per sample in wave file : %d \n", meta1.bits_per_sample );
    printf(" Channels: %d\n", meta1.num_channels);
    
    //    Get file length
    fseek(file1, 0, SEEK_END);
    fileLen1=ftell(file1);
    printf("Total bytes for samples: %lu\n", fileLen1);
    int count1 = fileLen1/2;
    printf("Number of samples (16bits or 2bytes per sample): %d\n", count1);
    
    fseek(file1, 0, SEEK_SET);
    
    //    Allocate memory
    buffer1=(short int *)malloc(fileLen1+1);
    if (!buffer1)
    {
        fprintf(stderr, "Memory error!");
        fclose(file1);
        return 0;
    }
    
    //    Read file contents into buffer
    fread(buffer1, fileLen1, 1, file1);
    fclose(file1);
    short int sample1[count1];
    for(int i=0 ; i <= count1; i++)
    {
        sample1[i]=((short int*)buffer1)[i];
        //printf("buffer[%d] = %d\n", i, input_signal[i]);
    }
    
    free(buffer1);
    
    float* array1 = malloc(sizeof(float) * count1);
    for(int i=0;i<count1;i++){
        array1[i]=intToFloat(sample1[i]);
    }
    
    FILE *file2;
    short int *buffer2;
    unsigned long fileLen2;
    struct header_file meta2;
    //    Open file
    file2 = fopen("BIG HALL E001 M2S.wav", "rb");
    if (!file2)
    {
        fprintf(stderr, "Unable to open file %s", "BIG HALL E001 M2S.wav");
        return 0;
    }
    //    Read header
    fread(&meta2, sizeof(meta2), 1, file2);
    printf(" Size of Header file : %lu  bytes\n", sizeof(struct header_file));
    printf(" Sampling rate of the input wave file : %d Hz \n", meta2.sample_rate);
    printf(" Bits per sample in wave file : %d \n", meta2.bits_per_sample );
    printf(" Channels: %d\n", meta2.num_channels);
    printf(" byte_rate: %d\n", meta2.byte_rate);
    
    //    Get file length
    fseek(file2, 0, SEEK_END);
    fileLen2=ftell(file2);
    printf("Total bytes for samples: %lu\n", fileLen2);
    int count2 = fileLen2/2;
    printf("Number of samples (16bits or 2bytes per sample): %d\n", count2);
    
    fseek(file2, 0, SEEK_SET);
    
    //    Allocate memory
    buffer2=(short int *)malloc(fileLen2+1);
    if (!buffer2)
    {
        fprintf(stderr, "Memory error!");
        fclose(file2);
        return 0;
    }
    
    //    Read file contents into buffer
    fread(buffer2, fileLen2, 1, file1);
    fclose(file2);
    short int sample2[count2];
    for(int i=0 ; i <= count2; i++)
    {
        sample2[i]=((short int*)buffer2)[i];
        //printf("buffer[%d] = %d\n", i, input_signal[i]);
    }
    
    free(buffer2);
    
    float* array2 = malloc(sizeof(float) * count2);
    for(int i=0;i<count2;i++){
        array2[i]=intToFloat(sample2[i]);
    }
    
    clock_t start, end;
    double elapsed;
    start = clock();
    
    int outputFileSize = count1 + count2 - 1;
    float* output = malloc(sizeof(float) * outputFileSize);

    
    convolve(array1, count1, array2, count2, output, outputFileSize);

    int16_t* resultArray = malloc(sizeof(int16_t)*outputFileSize);

    for(int i=0;i<outputFileSize;i++){
        resultArray[i]=floatToInt(output[i] * 0.2);
    }

//    int channels = meta.num_channels;
//    int byteRate = meta.byte_rate;
//
//    printf("Channels are %d \n", channels);
//    printf("ByteRate are %d \n", byteRate);

    writeWAV(1, outputFileSize, 44100, "output.wav", resultArray);

    end = clock();
    elapsed = (double)(end - start);
    printf("Block used %lf seconds.\n", elapsed / CLOCKS_PER_SEC);
    
    
    return 0;
}

/******************************************************************************
 *
 *       function:       createTestTone
 *
 *       purpose:        Calculates and writes out a sine test tone to file
 *
 *       arguments:      frequency:  frequency of the test tone in Hz
 *                       duration:  length of the test tone in seconds
 *                       numberOfChannels:  number of audio channels
 *                       filename:  name of the file to create
 *
 *       internal
 *       functions:      writeWaveFileHeader, fwriteShortLSB
 *
 *       library
 *       functions:      ceil, pow, fopen, fprintf, sin, rint, fclose
 *
 ******************************************************************************/

void createTestTone(double frequency, double duration,
                    int numberOfChannels, char *filename)
{
    int i;
    
    /*  Calculate the number of sound samples to create,
     rounding upwards if necessary  */
    int numberOfSamples = (int)ceil(duration * SAMPLE_RATE);
    
    /*  Calculate the maximum value of a sample  */
    int maximumValue = (int)pow(2.0, (double)BITS_PER_SAMPLE - 1) - 1;
    
    /*  Open a binary output file stream for writing */
    FILE *outputFileStream = fopen(filename, "wb");
    if (outputFileStream == NULL) {
        fprintf(stderr, "File %s cannot be opened for writing\n", filename);
        return;
    }
    
    /*  Write the WAVE file header  */
    writeWaveFileHeader(numberOfChannels, numberOfSamples,
                        SAMPLE_RATE, outputFileStream);
    
    /*  Create the sine tone and write it to file  */
    /*  Since the frequency is fixed, the angular frequency
     and increment can be precalculated  */
    double angularFrequency = 2.0 * PI * frequency;
    double increment = angularFrequency / SAMPLE_RATE;
    for (i = 0; i < numberOfSamples; i++) {
        /*  Calculate the sine wave in the range -1.0 to + 1.0  */
        double value = sin(i * increment);
        
        /*  Convert the value to a 16-bit integer, with the
         range -maximumValue to + maximumValue.  The calculated
         value is rounded to the nearest integer  */
        short int sampleValue = rint(value * maximumValue);
        
        /*  Write out the sample as a 16-bit (short) integer
         in little-endian format  */
        fwriteShortLSB(sampleValue, outputFileStream);
        
        /*  If stereo output, duplicate the sample in the right channel  */
        if (numberOfChannels == STEREOPHONIC)
            fwriteShortLSB(sampleValue, outputFileStream);
    }
    
    /*  Close the output file stream  */
    fclose(outputFileStream);
}



/******************************************************************************
 *
 *       function:       writeWaveFileHeader
 *
 *       purpose:        Writes the header in WAVE format to the output file.
 *
 *       arguments:      channels:  the number of sound output channels
 *                       numberSamples:  the number of sound samples
 *                       outputRate:  the sample rate
 *                       outputFile:  the output file stream to write to
 *
 *       internal
 *       functions:      fwriteIntLSB, fwriteShortLSB
 *
 *       library
 *       functions:      ceil, fputs
 *
 ******************************************************************************/

void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;
    
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
    
    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;
    
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);
    
    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
    
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
    
    /*  Form container type  */
    fputs("WAVE", outputFile);
    
    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
    
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);
    
    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);
    
    /*  Number of channels  */
    fwriteShortLSB((short)channels, outputFile);
    
    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);
    
    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);
    
    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);
    
    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);
    
    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);
    
    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}



/******************************************************************************
 *
 *       function:       fwriteIntLSB
 *
 *       purpose:        Writes a 4-byte integer to the file stream, starting
 *                       with the least significant byte (i.e. writes the int
 *                       in little-endian form).  This routine will work on both
 *                       big-endian and little-endian architectures.
 *
 *       internal
 *       functions:      none
 *
 *       library
 *       functions:      fwrite
 *
 ******************************************************************************/
size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];
    
    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}



/******************************************************************************
 *
 *       function:       fwriteShortLSB
 *
 *       purpose:        Writes a 2-byte integer to the file stream, starting
 *                       with the least significant byte (i.e. writes the int
 *                       in little-endian form).  This routine will work on both
 *                       big-endian and little-endian architectures.
 *
 *       internal
 *       functions:      none
 *
 *       library
 *       functions:      fwrite
 *
 ******************************************************************************/

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];
    
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}

/******************************************************************************
 *
 *       function:       writeWAV
 *
 *       purpose:        write the data into a .wav file
 *
 ******************************************************************************/

void writeWAV(int channels, int numberSamples, double outputRate, char *filename, int16_t *samples){

    /* Open a binary output file stream for writing */
    FILE *outputFile = fopen(filename, "wb");
    if (outputFile == NULL){
        fprintf(stderr, "The file %s cannot be opened and written!", filename);
    }

    /* Write the WAVE file header */
    writeWaveFileHeader(channels, numberSamples, outputRate, outputFile);
    for(int i = 0; i < numberSamples; i++){
        fwriteShortLSB(samples[i], outputFile);
        if(meta.num_channels == STEREOPHONIC){
            fwriteShortLSB(samples[i], outputFile);
        }
    }
    fclose(outputFile);
}

/*****************************************************************************
 *
 *    Function:     convolve
 *
 *    Description:  Convolves two signals, producing an output signal.
 *                  The convolution is done in the time domain using the
 *                  "Input Side Algorithm" (see Smith, p. 112-115).
 *
 *    Parameters:   x[] is the signal to be convolved
 *                  N is the number of samples in the vector x[]
 *                  h[] is the impulse response, which is convolved with x[]
 *                  M is the number of samples in the vector h[]
 *                  y[] is the output signal, the result of the convolution
 *                  P is the number of samples in the vector y[].  P must
 *                       equal N + M - 1
 *
 *****************************************************************************/

void convolve(float x[], int N, float h[], int M, float y[], int P)
{
    int n, m;
    
    /*  Make sure the output buffer is the right size: P = N + M - 1  */
    if (P != (N + M - 1)) {
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;
    }
    
    /*  Clear the output buffer y[] to all zero values  */
    for (n = 0; n < P; n++)
        y[n] = 0.0;
    
    /*  Do the convolution  */
    /*  Outer loop:  process each input value x[n] in turn  */
    for (n = 0; n < N; n++) {
        /*  Inner loop:  process x[n] with each sample of h[]  */
        for (m = 0; m < M; m++)
            y[n+m] += x[n] * h[m];
    }
}

/******************************************************************************
 *
 *       function:       intToFloat
 *
 *       purpose:        convert the signed 16 bit integer to a float between -1.0~1.0
 *
 *
 ******************************************************************************/

float intToFloat(int16_t x){
    float y;
    y = (1.0*x/32768);
    return y;
}

/******************************************************************************
 *
 *       function:       floatToInt
 *
 *       purpose:        convert a float between -1.0~1.0 to signed 16 bit integer
 *
 *
 ******************************************************************************/

int16_t floatToInt(float x){
    int16_t y;
    y = x*32768;
    return y;
}
