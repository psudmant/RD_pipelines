import numpy as np
cimport numpy as np
#cimport cython 


#@cython.boundscheck(False)
#@cython.wraparound(False)

def calc_gc_depth(np.ndarray[np.uint64_t, ndim=1] depths, np.ndarray[np.uint16_t, ndim=1] gc_content, np.ndarray[np.uint8_t, ndim=1] cp2_unmasked_bool, int max_n_gc):
	
	cdef int sum
	sum = 0
	cdef Py_ssize_t i
	cdef np.ndarray[np.uint64_t, ndim=1] sum_depth_at_gc = np.zeros(max_n_gc, dtype=np.uint64)
	cdef np.ndarray[np.int32_t, ndim=1] count_at_gc = np.zeros(max_n_gc, dtype=np.int32)
	
	for i in range(gc_content.shape[0]):

		if(cp2_unmasked_bool[i]==0):
			sum+=depths[i]
			sum_depth_at_gc[gc_content[i]]+=depths[i]
			count_at_gc[gc_content[i]]+=1
	
		#print gc_content[i]

	return sum,sum_depth_at_gc,count_at_gc
