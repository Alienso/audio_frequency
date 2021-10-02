#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "portaudio.h"
#include "fft.hpp"

#define PA_SAMPLE_TYPE      paFloat32
#define FRAMES_PER_BUFFER   (64)

static int gNumNoInputs = 0;
using namespace std;

/*static int fuzzCallback( const void *inputBuffer, void *outputBuffer,
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
			for( i=0; i<framesPerBuffer; i++ ){
				*out++ = *in++;
				*out++ = *in++;
			}
		}
		//printf("%f\n",*out);

		return paContinue;
	}*/
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
			CArray data(framesPerBuffer*2);
			for( i=0; i<framesPerBuffer; i++ ){
				/*data[i] = {*in++,0};
				data[i+1] = {*in++,0};*/
				*out++ = *in++;
				*out++ = *in++;
			}
			/*data.resize(data.size()/2);
			fft(data);
			vector<float> rez = abs_vec(data);
			print_vec(rez);
			plot_x(rez);
			char c;
			cin>>c;*/
			//plt::plot(rez);
			//plt::save("./basic.png");
			//print_vec(rez);
			//print_max(rez);
		}

		return paContinue;
	}

int main(){
	
	PaStreamParameters inputParameters, outputParameters;
	PaStream *stream;
	PaError err;
	
	err = Pa_Initialize();
	if (err != paNoError){
		printf("Pa_init");
		return -1;
	}
	
	inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
	if (inputParameters.device == paNoDevice) {
         fprintf(stderr,"Error: No default input device.\n");
         return -1;
     }
    inputParameters.channelCount = 1;       /* stereo input */
    inputParameters.sampleFormat = PA_SAMPLE_TYPE;
    inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
    inputParameters.hostApiSpecificStreamInfo = NULL;
	
	outputParameters.device = Pa_GetDefaultOutputDevice(); /* default output device */
    if (outputParameters.device == paNoDevice) {
        fprintf(stderr,"Error: No default output device.\n");
        return -1;
     }
     outputParameters.channelCount = 1;       /* stereo output */
     outputParameters.sampleFormat = PA_SAMPLE_TYPE;
     outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency;
    outputParameters.hostApiSpecificStreamInfo = NULL;
	
	err = Pa_OpenStream(
            &stream,
            &inputParameters,
            &outputParameters,
            44100,
            FRAMES_PER_BUFFER,
            0, /* paClipOff, */  /* we won't output out of range samples so don't bother clipping them */
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