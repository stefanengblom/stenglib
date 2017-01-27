//
// g++ bw-test.cpp -o bw-test -fopenmp
//
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>


double timenow(void) {

  double the_time_now = 0.0;

  struct timeval now;


  // change this for windows if needed
  gettimeofday(&now, NULL);
  the_time_now = now.tv_sec*1000000.0+(now.tv_usec);

  return the_time_now;

}


int main(void) {

  int N = 100000000;

  int *a, *b;

  a = (int *) malloc(N*sizeof(int));
  b = (int *) malloc(N*sizeof(int));

  double tick, tack;

  const int nThreads = omp_get_max_threads();

    // warm up
    for (int j=0; j<N; ++j)
      a[j] = b[j];

    tick = timenow();

    for (int j=0; j<N; ++j)
      a[j] = b[j];

    tack = timenow();

    double serial_time = tack-tick;
    std::cout << "Serial copy test:" << (tack-tick)/1000000 
              << " sec" << std::endl;

  
  for (int i=1; i<=nThreads; ++i) {

    omp_set_num_threads(i);

    // warm up
    for (int j=0; j<N; ++j)
      a[j] = b[j];

    tick = timenow();

#pragma omp parallel for
    for (int j=0; j<N; ++j)
      a[j] = b[j];

    tack = timenow();

    std::cout << "Copy test:" << (tack-tick)/1000000 
              << " sec; Threads:" << omp_get_max_threads() 
              << " ; speed-up:" << serial_time/(tack-tick)
              << std::endl;

  }

  return 0;
}
