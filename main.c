#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <pthread.h>


#define NUMTHREADS 2
#define SIZE 64 


#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
	  unsigned long long int x;
	  __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
	  return x;
}
#elif defined(__x86_64__)


static __inline__ unsigned long long rdtsc(void)
{
	  unsigned hi, lo;
	  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#elif defined(__powerpc__)
static __inline__ unsigned long long rdtsc(void)
{
	  unsigned long long int result=0;
	  unsigned long int upper, lower,tmp;
	  __asm__ volatile(
	                      "0:                  \n"
	                      "\tmftbu   %0           \n"
	                      "\tmftb    %1           \n"
	                      "\tmftbu   %2           \n"
	                      "\tcmpw    %2,%0        \n"
	                      "\tbne     0b         \n"
	                      : "=r"(upper),"=r"(lower),"=r"(tmp)
	      		   );
	 result = upper;
	 result = result<<32;
	 result = result|lower;

	 return(result);
}
#endif

/***********************************************************************/
/* START: MT 19937******************************************************/
/***********************************************************************/

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
	    mt[0]= s & 0xffffffffUL;
	        for (mti=1; mti<N; mti++) {
			        mt[mti] = 
						    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
				        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
				        /* In the previous versions, MSBs of the seed affect   */
				        /* only MSBs of the array mt[].                        */
				        /* 2002/01/09 modified by Makoto Matsumoto             */
				        mt[mti] &= 0xffffffffUL;
					        /* for >32 bit machines */
					    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
	    int i, j, k;
	        init_genrand(19650218UL);
		    i=1; j=0;
		        k = (N>key_length ? N : key_length);
			    for (; k; k--) {
				    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
					      + init_key[j] + j; /* non linear */
				    mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
				    i++; j++;
				    if (i>=N) { mt[0] = mt[N-1]; i=1; }
					    if (j>=key_length) j=0;
						}
			        for (k=N-1; k; k--) {
					        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
							          - i; /* non linear */
						        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
							        i++;
								        if (i>=N) { mt[0] = mt[N-1]; i=1; }
									    }

				    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
	    unsigned long y;
	static unsigned long mag01[2]={0x0UL, MATRIX_A};
	    /* mag01[x] = x * MATRIX_A  for x=0,1 */

	    if (mti >= N) { /* generate N words at one time */
	    int kk;

	    if (mti == N+1)   /* if init_genrand() has not been called, */
		init_genrand(5489UL); /* a default initial seed is used */

	    	for (kk=0;kk<N-M;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		for (;kk<N-1;kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
		}
		y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
		mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

		mti = 0;
				}

		y = mt[mti++];

		    /* Tempering */
		    y ^= (y >> 11);
			y ^= (y << 7) & 0x9d2c5680UL;
			    y ^= (y << 15) & 0xefc60000UL;
				y ^= (y >> 18);

				    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
	    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
	    return genrand_int32()*(1.0/4294967295.0); 
	        /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
	    return genrand_int32()*(1.0/4294967296.0); 
	        /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
	    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
	        /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
	    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
	        return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/***********************************************************************/
/* END: MT 19937 *******************************************************/
/***********************************************************************/

double **A=NULL;
double *B =NULL;
double *Bholder = NULL;
double **C=NULL;
double clock_rate=2666700000.0; // change to 700000000.0 for Blue Gene/L
unsigned int matrix_size=8192;
int worldSize, myRank, aSize, sendTo;
int offsetting;
int partitionSize;
unsigned long rng_init_seeds[6]={0x0, 0x123, 0x234, 0x345, 0x456, 0x789};
unsigned long rng_init_length=6;
long long sendTotal, recvTotal, compTotal;


void * multiplyRows(void * dealWithRows){
	int compStart, compEnd, i, j, k;
	//The row offset number that each thread will deal with	needs to be multiplied by rows/#threads to get actual offset
	int offset = *((int *) dealWithRows);

	//Calculate the matrix part
	compStart = rdtsc();



	//Each thread is only going to do part of the multiplication
	for(i = 0; i < partitionSize/NUMTHREADS; i++){
		for(k = 0; k < partitionSize; k++){
			for(j = 0; j < SIZE; j++){
				//Do the multiplication of A and part of B
				C[i + offset*partitionSize/NUMTHREADS][k + offsetting] += A[i + offset * partitionSize / NUMTHREADS][j] * B[k*SIZE + j];
			   }
		} 

	}
	compEnd = rdtsc();
	compTotal += compEnd - compStart;
	
	if(offset != 0){
		pthread_exit(0);
	}
	return NULL;
}



void multiply(){

	MPI_Request sentRequests[2];
	MPI_Request recvRequests[2];
	int i, outCount; 
	int statusArr[2];
	int recieved[2], sent[2];
	int offsetHolder; 
	int done = 0;
	int threads;
	long long recvStart, sendStart;
	MPI_Status statInfo[2];
	pthread_t threadHandler[NUMTHREADS-1];
	int doRowOffset[NUMTHREADS];



	while(done < worldSize){
		//Make sure to listen before you talk
		//Post the recieves but don't use them yet, use them after you are done computing information
		//The tag represents the number it should be recieving, 0 for first one 1 for second which will be the first it actually recieves, 2 for the third computation
		
		Bholder = (double *) calloc(partitionSize * SIZE, sizeof(double));

		recvStart = rdtsc();
		MPI_Irecv(&offsetHolder, 1, MPI_INT, MPI_ANY_SOURCE, done+1, MPI_COMM_WORLD, &(recvRequests[0]));
		MPI_Irecv(Bholder, partitionSize*SIZE, MPI_DOUBLE, MPI_ANY_SOURCE, done+1, MPI_COMM_WORLD, &(recvRequests[1]));
		recvTotal += rdtsc() - recvStart;
		

		
		doRowOffset[0] = 0;
		for(threads = 1; threads < NUMTHREADS; threads++){
			doRowOffset[threads] = threads;
			pthread_create(&threadHandler[threads-1], NULL, multiplyRows, (void *) &doRowOffset[threads]);
		}
		//Call the calculation function for the "main" thread for this task
		multiplyRows( &doRowOffset);

		for(threads = 1; threads < NUMTHREADS; threads++){
			pthread_join(threadHandler[threads-1], NULL);
		}

		//Check pointing done starts at 1 and goes up to n
		/*
		if(myRank == 0){
			printf("Done with %d calculation\n", done+1);
		}
		*/

		//Check to see if we have done n -1 multiplications	
		if(done == worldSize - 1){
			//Clean up outstanding recieves
			MPI_Cancel(&(recvRequests[0]));
			MPI_Cancel(&(recvRequests[1]));
			break;
		}


		done++;

		sendStart = rdtsc();
		//Send the necessary data. Offset and data
		MPI_Isend(&offsetting, 1, MPI_INT, sendTo, done, MPI_COMM_WORLD, &(sentRequests[0]));
		MPI_Isend(B, SIZE * partitionSize, MPI_DOUBLE, sendTo, done, MPI_COMM_WORLD, &(sentRequests[1]));

		sent[0] = 0;
		sent[1] = 0;
		while(sent[0] == 0 || sent[1] == 0){
			MPI_Testsome(2, sentRequests, &outCount, statusArr, statInfo);
			
			for(i = 0; i < outCount; i++){
				if(statusArr[i] == 0){
					sent[0] = 1;
				}else{
					sent[1] = 1;
				}
			}
		}
		sendTotal += rdtsc() - sendStart;

		recvStart = rdtsc();
		recieved[0] = 0;
		recieved[1] = 0;
		//Checking to see if we recieved both messages
		while(recieved[0] == 0 || recieved[1] == 0){
			MPI_Testsome(2, recvRequests, &outCount, statusArr, MPI_STATUSES_IGNORE);
			
			for(i = 0; i < outCount; i++){
				if(statusArr[i] == 0){
					recieved[0] = 1;
					offsetting = offsetHolder;
				}else{
					recieved[1] = 1;

					//Done with previous B so free the memory
					free(B);

					B = Bholder;
					//B now points to Bholder so we don't want to free Bholder or that will get rid of whats in B
					Bholder = NULL;
				}
			}
		}
		recvTotal += rdtsc() - recvStart;



	}

	//Do the final free on B
	free(B);
}

int main(int argc, char * argv[]){
	//Initialize everything
	int i, j;

	long long start = rdtsc(), finish;
	long long execTime, recvMax, recvMin, sendMax, sendMin, sendAv, recvAv, compAv, compMax, compMin;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	
	rng_init_seeds[0] = myRank;
	init_by_array(rng_init_seeds, rng_init_length);	
	
	partitionSize = SIZE / worldSize;


	//Correctly get the offset for the block that this process is making
	

	if(myRank == worldSize - 1){
		sendTo = 0;
	}else{
		sendTo = myRank+1;
	}

	offsetting = partitionSize * myRank;

	aSize = partitionSize;

	A = (double **) calloc(partitionSize,  sizeof(double *));
	//B is in row-column-major form for easier sending
	B = (double *) calloc(SIZE * partitionSize, sizeof(double));
	C = (double **) calloc(partitionSize, sizeof(double *));


	//Start initilization
	for(i = 0; i < partitionSize; i++){
		A[i] = (double *) calloc(SIZE, sizeof(double));
		C[i] = (double *) calloc(SIZE, sizeof(double));
		for(j = 0; j < SIZE; j++){
			A[i][j] = genrand_res53(); //Set the random values for A
			//Set the random values for B but since its in row-column-major do some fancy manipulations
			B[i*SIZE + j] = genrand_res53(); 
			C[i][j] = 0;
		}
	}


	/*
	if(myRank == 0){
		printf("A:\n");
		for(i = 0; i < partitionSize; i++){
			for(j = 0; j < SIZE; j++){
				printf("%f", A[i][j]);
			}
			printf("\n");
		}
	}	
	printf("\n");
		for(j = 0; j < SIZE; j++){
			for(i = 0; i < partitionSize; i++){
				printf("B%d: %f", myRank, B[i*SIZE + j]);
			}
			printf("\n");
		}

*/
		//After sending the information we just calculated we can now go into our processing of data recieved


	multiply();

	//Check output of C after the multiplication can be uncommented to make sure its correct
/*	if(myRank == 0){
		printf("C: \n");
		for(i = 0; i < partitionSize; i++){
			for(j = 0; j < SIZE; j++){
				printf("%f", C[i][j]);
			}
			printf("\n");
		}
		
		printf("\n\n");
	}*/
	finish = rdtsc();
	execTime = finish - start;



	//Compute the timing information leave as cycles for now
	

	MPI_Allreduce(&sendTotal, &sendMax, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&sendTotal, &sendAv, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	sendAv = sendAv / worldSize;
	MPI_Allreduce(&sendTotal, &sendMin, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);

	MPI_Allreduce(&recvTotal, &recvMax, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&recvTotal, &recvAv, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	recvAv = recvAv / worldSize;
	MPI_Allreduce(&recvTotal, &recvMin, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
	
	MPI_Allreduce(&compTotal, &compMax, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&compTotal, &compAv, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	compAv = compAv / worldSize;
	MPI_Allreduce(&compTotal, &compMin, 1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);


	if(myRank == 0)
		printf("Execution time: %lf\nSend time max, min, av: %lf, %lf, %lf\nRecv time max, min, av: %lf, %lf, %lf\nComputation time max, min, av: %lf, %lf, %lf\n", execTime/clock_rate, sendMax/clock_rate, sendMin/clock_rate, sendAv/clock_rate, recvMax/clock_rate, recvMin/clock_rate, recvAv/clock_rate, compMax/clock_rate, compMin/clock_rate, compAv/clock_rate);


	MPI_Finalize();
	return (0);
}
