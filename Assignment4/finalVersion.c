//
//  codeT5.c
//  
//
//  Created by Zicheng Huang on 2019/4/12.
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

#define SIZE        8
#define TWO_PI      (2.0 * PI)
#define SWAP(a,b)   tempr=(a);(a)=(b);(b)=tempr


/*  FUNCTION PROTOTYPES  *****************************************************/

int16_t* readWav(char* fileName, int* dataLength);

void writeWaveFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile);

size_t fwriteIntLSB(int data, FILE *stream);

size_t fwriteShortLSB(short int data, FILE *stream);

void writeWAV(int channels, int numberSamples, double outputRate, char *filename, int16_t *samples);

double intToDouble(int16_t x);
int16_t doubleToInt(double x);

double findMax(double arr[], int n);
double findMin(double arr[], int n);

/* FFT function */

void four1(double data[], int nn, int isign);

/*******************************************************************************/

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
    
    if(argc != 4){
        printf("This program needs 3 command arguments.\n");
        exit(-1);
    }
    
    char* inputfile = argv[1];
    char* IRfile = argv[2];
    char* outputfile = argv[3];
    int count1;
    int count2;

    clock_t start, end;
    double elapsed;
    start = clock();
    
    int16_t* sample1 = readWav(inputfile, &count1);
    int16_t* sample2 = readWav(IRfile, &count2);
    
    double* array1 = malloc(sizeof(double) * count1);
    for(int i=0;i<count1;i++){
        array1[i]=intToDouble(sample1[i]);
    }
    
    double* array2 = malloc(sizeof(double) * count2);
    for(int i=0;i<count2;i++){
        array2[i]=intToDouble(sample2[i]);
    }
    
    printf("-------------------------No.1---------------------------\n");
    
    // Calculate the length for using in the fourier function
    int outputFileSize = count1 + count2 - 1;
    
    printf("Old length: %d\n", outputFileSize);
    double k = log2(outputFileSize);
    int y = ceil(k);
    int length = pow(2, y);
    
    printf("New length: %d\n", length);
    
    int dlength = length * 2;
    
    double* complexDataArr1 = malloc(sizeof(double)*dlength);
    
    double* complexDataArr2 = malloc(sizeof(double)*dlength);
    
    printf("--------------------------No.2---------------------------\n");
    
    // Do the initialization for each complex array (signal is real and the imaginary parts are set to 0.0)
    register int m, n, p, q;
    
    for(m=0,n=0,p=0,q=0;n<count1;m+=2,n++,p+=2,q++){
        complexDataArr1[m] = array1[n];
        complexDataArr1[m+1] = 0.0;
        if(q < count2){
            complexDataArr2[p] = array2[q];
            complexDataArr2[p+1] = 0.0;
        }
    }
    
    for(m=2*count1, n=2*count2; n<dlength; m++, n++){
        if(m < dlength){
            complexDataArr1[m] = 0.0;
        }
        complexDataArr2[n] = 0.0;
    }
    
    printf("---------------------------No.3---------------------------\n");
   
    // Do FFT
    four1(complexDataArr1-1, length, 1);
    four1(complexDataArr2-1, length, 1);
    
    // Do complex calculation (Multiply the frequency spectrum by the frequency response)
    double* complexOutput = malloc(sizeof(double)*dlength);
    for(int i=0; i<dlength; i+=2){
        double tempArR1 = complexDataArr1[i];
        double tempArI1 = complexDataArr1[i+1];
        double tempArR2 = complexDataArr2[i];
        double tempArI2 = complexDataArr2[i+1];
 
        complexOutput[i] = tempArR1*tempArR2-tempArI1*tempArI2;
        complexOutput[i+1] = tempArR1*tempArI2+tempArI1*tempArR2;
    }
    
    // Inverse FFT (IFFT)
    four1(complexOutput-1, length, -1);
  
    /************************************************************************************/
  
    register int j, h;
    
    for(j=0, h=0; j<length; j++, h+=2) {
        complexOutput[h] = complexOutput[h] / ((double)(length));
        complexOutput[h+1] = complexOutput[h+1] / ((double)(length));
    }
 
    
    double* output = malloc(sizeof(double) * outputFileSize);
 
    for(j=0, h=0; j<outputFileSize; j++, h+=2){
        output[j] = complexOutput[h];
    }
    
    double* temp = malloc(sizeof(double) * outputFileSize);
    
    for(j=0; j < outputFileSize-1; j+=2){
        temp[j] = output[j];
        temp[j+1] = output[j+1];
    }
    
    if(j == outputFileSize-1){
        temp[outputFileSize-1] = output[outputFileSize-1];
    }

    // Do normalization
    double max = findMax(temp, outputFileSize);
    double min = findMin(temp, outputFileSize);
    
    int16_t* resultArray = malloc(sizeof(int16_t)*outputFileSize);
    
    double denominator = max - min;
    for(int i=0; i<outputFileSize; i++){
        output[i] = ((2*(output[i]-min)) / denominator) - 1;
    }
    
    // Transform the double into int16_t and then write into file
    for(int i=0; i<outputFileSize; i++){
        resultArray[i]= doubleToInt(output[i]);
    }
    
    writeWAV(MONOPHONIC, outputFileSize, SAMPLE_RATE, outputfile, resultArray);
    
    end = clock();
    elapsed = (double)(end - start);
    printf("CodeTuning Version5 Optimization program used %lf seconds.\n", elapsed / CLOCKS_PER_SEC);
    
    free(sample1);
    free(sample2);
    free(array1);
    free(array2);
    free(output);
    free(temp);
    free(resultArray);
    return 0;
}
/******************************************************************************
 *
 *       function:       readWav
 *
 *       purpose:        read the data into an int16_t array
 *
 ******************************************************************************/

int16_t* readWav(char* fileName, int* dataLength){
    FILE *file;
    short int *buffer;
    unsigned long fileLen;
    struct header_file meta;
    //    Open file
    file = fopen(fileName, "rb");
    if (!file)
    {
        fprintf(stderr, "Unable to open file %s", fileName);
        return 0;
    }
    //    Read header
    fread(&meta, sizeof(meta), 1, file);
    printf(" Size of Header file : %lu  bytes\n", sizeof(struct header_file));
    printf(" Sampling rate of the input wave file : %d Hz \n", meta.sample_rate);
    printf(" Bits per sample in wave file : %d \n", meta.bits_per_sample );
    printf(" Channels: %d\n", meta.num_channels);
    
    //    Get file length
    fseek(file, 0, SEEK_END);
    fileLen=ftell(file);
    printf("Total bytes for samples: %lu\n", fileLen);
    int count = fileLen/2;
    printf("Number of samples (16bits or 2bytes per sample): %d\n", count);
    
    fseek(file, 0, SEEK_SET);
    
    //    Allocate memory
    buffer=(short int *)malloc(fileLen+1);
    if (!buffer)
    {
        fprintf(stderr, "Memory error!");
        fclose(file);
        return 0;
    }
    
    //    Read file contents into buffer
    fread(buffer, fileLen, 1, file);
    fclose(file);

    int16_t* sample = malloc(sizeof(int16_t)*count);
    * dataLength = count;
    for(int i=0 ; i <= count; i++)
    {
        sample[i]=((short int*)buffer)[i];
    }
    
    free(buffer);

    return sample;
}
/******************************************************************************
 *
 *       function:       four1
 *
 *       purpose:        do fft
 *
 *
 ******************************************************************************/

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;
    
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j+1], data[i+1]);
        }
        m = nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    
    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j+1];
                tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
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

void writeWAV(int channels, int numberSamples, double outputRate, char *filename, int16_t *resultSamples){
    
    /* Open a binary output file stream for writing */
    FILE *outputFile = fopen(filename, "wb");
    if (outputFile == NULL){
        fprintf(stderr, "The file %s cannot be opened and written!", filename);
    }
    
    /* Write the WAVE file header */
    writeWaveFileHeader(channels, numberSamples, outputRate, outputFile);
    for(int i = 0; i < numberSamples; i++){
        /* Use the fwriteFunction to write the data into outputFile*/
        fwriteShortLSB(resultSamples[i], outputFile);
        if(meta.num_channels == STEREOPHONIC){
            fwriteShortLSB(resultSamples[i], outputFile);
        }
    }
    fclose(outputFile);
}

/******************************************************************************
 *
 *       function:       intToDouble
 *
 *       purpose:        convert the signed 16 bit integer to a double between -1.0~1.0
 *
 *
 ******************************************************************************/

double intToDouble(int16_t x){
    double y;
    y = (1.0*x/32767.0);
    return y;
}

/******************************************************************************
 *
 *       function:       doubleToInt
 *
 *       purpose:        convert a double between -1.0~1.0 to signed 16 bit integer
 *
 *
 ******************************************************************************/

int16_t doubleToInt(double x){
    int16_t y;
    y = x*32767.0;
    return y;
}

/******************************************************************************
 *
 *       function:       findMax
 *
 *       purpose:        find the maximum value of a given array
 *
 *
 ******************************************************************************/

double findMax(double arr[], int n){
    
    int i;

    // Initialize maximum element
    double max = arr[0];

    // Traverse array elements
    // from second and compare
    // every element with current max
    for (i = 1; i < n; i++){
        if (arr[i] > max){
            max = arr[i];
        }
    }
    return max;
}

/******************************************************************************
 *
 *       function:       findMin
 *
 *       purpose:        find the minimum value of a given array
 *
 *
 ******************************************************************************/

double findMin(double arr[], int n){
    
    int i;
    
    // Initialize maximum element
    double min = arr[0];
    
    // Traverse array elements
    // from second and compare
    // every element with current min
    for (i = 1; i < n; i++){
        if (arr[i] < min){
            min = arr[i];
        }
    }
    return min;
}
