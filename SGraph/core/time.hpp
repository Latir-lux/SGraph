#ifndef TIME_HPP
#define TIME_HPP

#include <sys/time.h>

inline double get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);//获取当前时间
    return tv.tv_sec + (tv.tv_usec / 1e6);//tv_sec是秒，tv_usec是微秒
}

#endif
