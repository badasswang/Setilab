#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

//new stuff from pthread
#include <sched.h>
#include <pthread.h>
#include <unistd.h>
//new stuff
pthread_t *tid;
int num_threads;
int num_procs;
typedef struct args {
int idnum;
    signal *sig;
int filter_order;
int num_bands;
double *lb;
double *ub;
int num_threads;
int num_procs;;
}ARG;
double *band_power;
double Fc,bandwidth;
void usage() 
{
    printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads num_processors \n");
}

double avg_power(double *data, int num)
{
    int i;
    double ss;
    
    ss=0;
    for (i=0;i<num;i++) { 
	ss += data[i]*data[i];
    }
    
    return ss/num;
}

double max_of(double *data, int num)
{
    double m=data[0];
    int i;
    
    for (i=1;i<num;i++) { 
	if (data[i]>m) { m=data[i]; } 
    }
    return m;
}

double avg_of(double *data, int num)
{
    double s=0;
    int i;
    
    for (i=0;i<num;i++) { 
	s+=data[i];
    }
    return s/num;
}

void remove_dc(double *data, int num)
{
  int i;
  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (i=0;i<num;i++) {
    data[i] -= dc;
  }
}

void preprocess(signal *sig, int num_bands)
{
    double signal_power;
    Fc=(sig->Fs)/2;
    bandwidth = Fc/num_bands;
    remove_dc(sig->data,sig->num_samples);

    signal_power = avg_power(sig->data,sig->num_samples);

    printf("signal average power:     %lf\n", signal_power);
}

void *worker(void *input){
    ARG* args = (ARG*) input;
    signal *sig = args->sig;
    int filter_order = args->filter_order;
    int num_bands = args->num_bands;
    double *lb = args->lb;
    double *ub = args->ub;
    int num_threads = args->num_threads;
    int num_procs = args->num_procs;

   

    double filter_coeffs[filter_order+1];

    long myid = args->idnum;
    int blocksize=num_bands / num_threads;
    int mystart = myid * blocksize;
    int myend;
    int band;
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(myid%num_procs,&set);
    if (sched_setaffinity(0,sizeof(set),&set)<0) {//do it
        perror("Can't set affinity"); //hopefully doesn't fail
        exit(-1);
        } //copied from pthread

    if (myid == (num_threads-1)) {myend=num_bands;}
    else {myend = (myid+1) * blocksize;}
    //end new


    for (band=mystart;band<myend;band++){
        generate_band_pass(sig->Fs, 
		    	   band*bandwidth+0.0001, // keep within limits
			       (band+1)*bandwidth-0.0001,
			       filter_order, 
			       filter_coeffs);
	    hamming_window(filter_order,filter_coeffs);

	// Convolve
	    convolve_and_compute_power(sig->num_samples,
				                   sig->data,
				                   filter_order,
				                   filter_coeffs,
				                   &(band_power[band]));
    }

    pthread_exit(NULL);

}

int analyze_signal(double *band_power, int num_bands, double *lb, double *ub) {
   int band;
    // Pretty print results
    double max_band_power = max_of(band_power,num_bands);
    double avg_band_power = avg_of(band_power,num_bands);
    int i;
    int wow=0;

#define MAXWIDTH 40

#define THRESHOLD 2.0

#define ALIENS_LOW   50000.0
#define ALIENS_HIGH  150000.0

    *lb=*ub=-1;

    for (band=0;band<num_bands;band++) {//edited 
      double band_low = band*bandwidth+0.0001;
      double band_high = (band+1)*bandwidth-0.0001;
      
      printf("%5d %20lf to %20lf Hz: %20lf ", 
	     band, band_low, band_high, band_power[band]);
      
      for (i=0;i<MAXWIDTH*(band_power[band]/max_band_power);i++) {
	printf("*");
      }
      
      if ( (band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
	   (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) { 

	// band of interest

	if (band_power[band] > THRESHOLD * avg_band_power) { 
	  printf("(WOW)");
	  wow=1;
	  if (*lb<0) { *lb=band*bandwidth+0.0001; }
	  *ub = (band+1)*bandwidth-0.0001;
	} else {
	  printf("(meh)");
	}
      } else {
	printf("(meh)");
      }
      
      printf("\n");
    }
    
    return wow;

}

int main(int argc, char *argv[])
{
    signal *sig;
    double Fs;
    char sig_type;
    char *sig_file;
    int filter_order;
    int num_bands;
    double start, end;
    int i;

    if (argc!=8) { 
	usage();
	return -1;
    }
  

    sig_type = toupper(argv[1][0]);
    sig_file = argv[2];
    Fs = atof(argv[3]);
    filter_order = atoi(argv[4]);
    num_bands = atoi(argv[5]);
    num_threads = atoi(argv[6]);
    num_procs = atoi(argv[6]);
    if (num_threads > num_bands)
        num_threads = num_bands; 
    band_power = malloc(num_bands * sizeof(double));

    assert(Fs>0.0);
    assert(filter_order>0 && !(filter_order & 0x1));
    assert(num_bands>0);

    printf("type:     %s\n"
	   "file:     %s\n"
	   "Fs:       %lf Hz\n"
	   "order:    %d\n"
	   "bands:    %d\n",
	   sig_type=='T' ? "Text" : sig_type=='B' ? "Binary" : sig_type=='M' ? "Mapped Binary" : "UNKNOWN TYPE",
	   sig_file,
	   Fs,
	   filter_order,
	   num_bands);
    
    printf("Load or map file\n");
    
    switch (sig_type) {
	case 'T':
	    sig = load_text_format_signal(sig_file);
	    break;

	case 'B':
	    sig = load_binary_format_signal(sig_file);
	    break;

	case 'M':
	    sig = map_binary_format_signal(sig_file);
	    break;
	    
	default:
	    printf("Unknown signal type\n");
	    return -1;
    }
    
    if (!sig) { 
	printf("Unable to load or map file\n");
	return -1;
    }

    sig->Fs=Fs;
    preprocess(sig, num_bands); 
    //new stuff
   
    long rc;
 tid = (pthread_t *) malloc(sizeof(pthread_t)*num_threads);


 for (i=0;i<num_threads;i++) {
     ARG* myargs = (ARG*) malloc(sizeof(ARG)); //need to initialize, ask Shae, also do I want to make a new struct for each thread?
     myargs->idnum = i;   
     myargs->sig = sig;
     myargs->num_bands = num_bands;
     myargs->lb = &start;
        myargs->ub = &end;
        myargs->filter_order = filter_order;
        myargs->num_threads = num_threads;
        myargs->num_procs = num_procs;
        rc = pthread_create( &(tid[i]),
                NULL,
                worker,
                myargs);
        if (rc!=0) {
            perror("failed to start thread");
            exit(-1);
        }}

    end = get_seconds();
  
    

    for (i=0;i<num_threads;i++) {
         rc=pthread_join(tid[i],NULL);   // 
         if (rc!=0) { 
        perror("join failed");
        exit(-1);
        }
    }

//don't know how to implement this
   if (analyze_signal(band_power,num_bands,&start,&end)) { 
	printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n",start,end,(end+start)/2.0);
    } else {
	printf("no aliens\n");
    } 


    free_signal(sig);

    return 0;
}


    
