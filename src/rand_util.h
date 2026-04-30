#ifndef __RAND_UTIL_H__
#define __RAND_UTIL_H__


#include "point.h"


using namespace std;

class RandUtil{
public:
    static void get_random_direction(size_t dim, Point& rp);
    static double randUnif(double lo, double hi);
};

#endif /* __RAND_UTIL_H__ */
