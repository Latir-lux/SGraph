#ifndef TIME_HPP
#define TIME_HPP

#include <sys/time.h>

inline double get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);//��ȡ��ǰʱ��
    return tv.tv_sec + (tv.tv_usec / 1e6);//tv_sec���룬tv_usec��΢��
}

#endif
