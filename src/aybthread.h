#ifndef AYBTHREAD_H
#define AYBTHREAD_H

#ifdef _OPENMP
	#include <omp.h>
#else
	static inline int omp_get_thread_num(void){
		return 0;
	}

	static inline int omp_get_max_threads(void){
		return 1;
	}

	static inline void omp_set_num_threads(int t){}
#endif

#endif /* AYBTHREAD_H */
