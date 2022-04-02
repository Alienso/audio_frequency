#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <complex>
#include <valarray>
#include <cmath>
#include "portaudio.h"
#include "fft.hpp"

#define PA_SAMPLE_TYPE      paFloat32
#define FRAMES_PER_BUFFER   (1024)
#define SAMPLE_RATE 48000

static int gNumNoInputs = 0;
using namespace std;

float data[2*FRAMES_PER_BUFFER];
short int* song;

short int* readWAV(char* path);

void applyHann(float data[2*FRAMES_PER_BUFFER]){
	for (int i=0;i<2*FRAMES_PER_BUFFER;i++)
		data[i] = data[i] * pow(sin((M_PI*data[i])/(2*FRAMES_PER_BUFFER)),2);
	return;
}

int calculate(){
	int bufferSize = 2*FRAMES_PER_BUFFER;
	double magnitude[bufferSize/2];
	valarray<Complex> fftArray(bufferSize);

	/*cout<<"DATA:\n";
	for(int i=0;i<bufferSize;i++){
		cout<<data[i]<<" ";
	}
	cout<<"\n";*/

	for (int i=0;i<bufferSize;i++)
		fftArray[i] = Complex(data[i],0);
	
	fft(fftArray);

	/*cout<<"FFTDATA:\n";
	for(int i=0;i<bufferSize;i++){
		cout<<fftArray[i]<<" ";
	}
	cout<<"\n";*/

	for(int i=0; i<(bufferSize/2) - 1;++i){
		double real = fftArray[i].real();
		double im = fftArray[i].imag();
		magnitude[i] = sqrt(real*real + im*im);
	}

	/*cout<<"MAGNITUDE:\n";
	for(int i=0;i<bufferSize/2;i++){
		cout<<magnitude[i]<<" ";
	}
	cout<<"\n";*/

	double max = magnitude[0];
	int index = 0;
	for (int i=0; i<bufferSize/2;++i){
		if (10*log10(magnitude[i]*magnitude[i])>-200){
			if(magnitude[i]>max){
				max = magnitude[i];
				index = i;
			}
		}
	}
	//cout<<"INDEX: " << index <<"\n";

	return SAMPLE_RATE * index / bufferSize;
}

static int fuzzCallback( const void *inputBuffer, void *outputBuffer,
                        unsigned long framesPerBuffer,
                        const PaStreamCallbackTimeInfo* timeInfo,
                        PaStreamCallbackFlags statusFlags,
                        void *userData )
	{
		const float *in = (const float*)inputBuffer;
		float *out = (float*)outputBuffer;
		
		unsigned int i=0;
		if( inputBuffer == NULL ){
			for( i=0; i<framesPerBuffer; i++ ){
				*out++ = 0;
				*out++ = 0; 
			}
			gNumNoInputs += 1;
		}
		else{
			for( i=0; i<2*framesPerBuffer; i+=2 ){
				data[i] = *in++;
				data[i+1] = *in++;
				/**out++ = *in++;
				*out++ = *in++;*/
			}

			applyHann(data);

			int rez = calculate();
			//if(rez>50 && rez<2000)
				printf("%d\n",rez);
		}

		return paContinue;
	}

long song_offset = 0;
static int songCallback( const void *inputBuffer, void *outputBuffer,
                        unsigned long framesPerBuffer,
                        const PaStreamCallbackTimeInfo* timeInfo,
                        PaStreamCallbackFlags statusFlags,
                        void *userData )
	{
		const float *in = (const float*)inputBuffer;
		float *out = (float*)outputBuffer;
		
		unsigned int i=0;
		for( i=0; i<2*framesPerBuffer; i+=2 ){
			*out++ = song[song_offset*2*framesPerBuffer + i];
			*out++ = song[song_offset*2*framesPerBuffer + i+1];
		}
		song_offset++;
		return paContinue;
	}

void test();
int main(){
	
	//song = readWAV("song.wav");

	PaStreamParameters inputParameters, outputParameters;
	PaStream *stream;
	PaError err;
	
	err = Pa_Initialize();
	if (err != paNoError){
		printf("Pa_init");
		return -1;
	}
	
	inputParameters.device = Pa_GetDefaultInputDevice(); 
	if (inputParameters.device == paNoDevice) {
         fprintf(stderr,"Error: No default input device.\n");
         return -1;
     }
    inputParameters.channelCount = 1;      
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;
	
	outputParameters.device = Pa_GetDefaultOutputDevice(); 
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        return -1;
     }
     outputParameters.channelCount = 1;       
     outputParameters.sampleFormat = PA_SAMPLE_TYPE;
     outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;
	

	err = Pa_OpenStream(
            &stream,
            &inputParameters,
            &outputParameters,
            SAMPLE_RATE,
            FRAMES_PER_BUFFER,
            0,
            fuzzCallback,
            NULL );
     if( err != paNoError ) return -1;

    err = Pa_StartStream( stream );
    if( err != paNoError ) return -1;
 
    printf("Hit ENTER to stop program.\n");
    getchar();
    err = Pa_CloseStream( stream );
    if( err != paNoError ) return -1;
 
    printf("Finished. gNumNoInputs = %d\n", gNumNoInputs );
    Pa_Terminate();
 
	return 0;
}

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// WAVE PCM soundfile format (you can find more in https://ccrma.stanford.edu/courses/422/projects/WaveFormat/ )
typedef struct header_file
{
    char chunk_id[4];
    int chunk_size;
    char format[4];
    char subchunk1_id[4];
    int subchunk1_size;
    short int audio_format;
    short int num_channels;
    int sample_rate;			// sample_rate denotes the sampling rate.
    int byte_rate;
    short int block_align;
    short int bits_per_sample;
    char subchunk2_id[4];
    int subchunk2_size;			// subchunk2_size denotes the number of samples.
} header;

typedef struct header_file* header_p;




short int* readWAV(char* path){
	FILE * infile = fopen(path,"rb");		// Open wave file in read mode
	FILE * outfile = fopen("Output.wav","wb");		// Create output ( wave format) file in write mode

	int BUFSIZE = 512;					// BUFSIZE can be changed according to the frame size required (eg:512)
	int count = 0;						// For counting number of frames in wave file.
	short int buff16[BUFSIZE];				// short int used for 16 bit as input data format is 16 bit PCM audio
	header_p meta = (header_p)malloc(sizeof(header));	// header_p points to a header struct that contains the wave file metadata fields
	int nb;							// variable storing number of byes returned
	short int* song_data = (short int*)malloc(SAMPLE_RATE*3000*sizeof(short int));
	if (infile)
	{
		fread(meta, 1, sizeof(header), infile);
		fwrite(meta,1, sizeof(*meta), outfile);
		cout << " Size of Header file is "<<sizeof(*meta)<<" bytes" << endl;
		cout << " Sampling rate of the input wave file is "<< meta->sample_rate <<" Hz" << endl;
		cout << " Number of samples in wave file are " << meta->subchunk2_size << " samples" << endl;



		while (!feof(infile))
		{
			nb = fread(buff16,1,BUFSIZE,infile);		// Reading data in chunks of BUFSIZE
			//cout << nb <<endl;
			count++;					// Incrementing Number of frames

			for(int j = 0;j<nb;j++){
				song_data[count*BUFSIZE+j] = buff16[j];
			}
			//memcpy(song_data+count*BUFSIZ,buff16,nb);

			fwrite(buff16,1,nb,outfile);			// Writing read data into output file
		}
	cout << " Number of frames in the input wave file are " <<count << endl;
	//fwrite(song_data,1,(count-1)*BUFSIZE,outfile);
	//fwrite(song_data,1,nb,outfile);
	}
	return 0;
}

/*class Recording extends Thread {

    @Override
    public void run() {

        while () {

            if (true) {                   
                int bufferReadResult = audioInput.read(buffer, 0, bufferSizeInBytes); // record data from mic into buffer
                if (bufferReadResult > 0) {
                    calculate();
                }
            } 
        }
    }
}*/

/*public void calculate() {

    double[] magnitude = new double[bufferSizeInBytes / 2];

    //Create Complex array for use in FFT
    Complex[] fftTempArray = new Complex[bufferSizeInBytes];
    for (int i = 0; i < bufferSizeInBytes; i++) {
        fftTempArray[i] = new Complex(buffer[i], 0);
    }

    //Obtain array of FFT data
    final Complex[] fftArray = FFT.fft(fftTempArray);
    // calculate power spectrum (magnitude) values from fft[]
    for (int i = 0; i < (bufferSizeInBytes / 2) - 1; ++i) {

        double real = fftArray[i].re();
        double imaginary = fftArray[i].im();
        magnitude[i] = Math.sqrt(real * real + imaginary * imaginary);

    }

    // find largest peak in power spectrum
    double max_magnitude = magnitude[0];
    int max_index = 0;
    for (int i = 0; i < magnitude.length; ++i) {
        if (magnitude[i] > max_magnitude) {
            max_magnitude = (int) magnitude[i];
            max_index = i;
        }
    }
    double freq = 44100 * max_index / bufferSizeInBytes;//here will get frequency in hz like(17000,18000..etc)        

}*/

//-----------------------------------------------------

void plot_x_y(vector<float> &x,vector<int> &y){
	return;
}

void plot_x(vector<float> &data){
	const int acc = 16;
	sort(data.begin(),data.end());
	vector<int> k(acc);
	float max = data[acc-1];
	float min = data[0];
	float d = (max-min)/(acc*1.0);
	cout<<max<<" "<<min<<"\n";
	for(int i=0;i<k.size();i++){
		k[i]=0;
		for (int j=0;j<data.size();j++){
			if(data[j]<d*(i+1) && data[j]>d*i )
				k[i]++;
		}
		cout<<k[i]<<" ";
	}
	cout<<"\n";
	return;
}

void print_vec(vector<float> &x){
	for(int i=0;i<x.size();i++){
		cout<<x[i]<<' ';
	}
	cout<<'\n';
	return;
}

vector<float> abs_vec(CArray &data){
	vector<float> ret(data.size());
	for(int i=0;i<data.size();i++)
		ret[i] = abs(data[i]);
	return ret; 
}

void print_max(vector<float> &x){
	float max = 0;
	for(int i=0;i<x.size();i++){
		if (x[i]>max)
			max = x[i];
	}
	cout<<max<<'\n';
	return;
}


void test(){

	CArray x = {Complex(-0.03480425839330703,0),Complex(0.07910192950176387,0),Complex(0.7233322451735928,0),Complex(0.1659819820667019,0)};
	fft(x);
	for (int i=0;i<4;i++){
		cout<<x[i]<<" ";
	}
	cout<<"/n";
    /**  -------------------
    *  0.9336118983487516
    *  -0.7581365035668999 + 0.08688005256493803i
    *  0.44344407521182005
    *  -0.7581365035668999 - 0.08688005256493803i
    **/
}