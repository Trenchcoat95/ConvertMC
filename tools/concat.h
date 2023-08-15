// Template function to concatenate two vectors.



#ifndef concatVector_h
#define concatVector_h



#include <vector>

template <typename T> std::vector<T> concat(std::vector<T>& vectorA, std::vector<T>& vectorB) {
    std::vector<T> retval;
    retval.reserve( vectorA.size() +vectorB.size() );
    retval.insert(retval.end(), vectorA.begin(),vectorA.end() );
    retval.insert(retval.end(), vectorB.begin(),vectorB.end() );
	return retval;
}



#endif
