#ifndef ATOMIC_HPP
#define ATOMIC_HPP

#include <atomic>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define ____force_inline __attribute__((always_inline))
#define compiler_fence() asm volatile("" ::: "memory")

template <typename T>
inline typename std::enable_if<sizeof(T) == 1, bool>::type cas(T *ptr, T oldv, T newv)//std::enable_if技术来条件性地启用此函数模板。它只在sizeof(T) == 1（即类型T的大小为1字节）时起作用
{
    static_assert(sizeof(char) == 1);
    return __sync_bool_compare_and_swap((char *)ptr, *((char *)&oldv), *((char *)&newv));
}

template <typename T>
inline typename std::enable_if<sizeof(T) == 2, bool>::type cas(T *ptr, T oldv, T newv)
{
    static_assert(sizeof(short) == 2);
    return __sync_bool_compare_and_swap((short *)ptr, *((short *)&oldv), *((short *)&newv));
}

template <typename T>
inline typename std::enable_if<sizeof(T) == 4, bool>::type cas(T *ptr, T oldv, T newv)
{
    static_assert(sizeof(int) == 4);
    return __sync_bool_compare_and_swap((int *)ptr, *((int *)&oldv), *((int *)&newv));
}

template <typename T>
inline typename std::enable_if<sizeof(T) == 8, bool>::type cas(T *ptr, T oldv, T newv)
{
    static_assert(sizeof(long) == 8);
    return __sync_bool_compare_and_swap((long *)ptr, *((long *)&oldv), *((long *)&newv));
}

template <typename T>
inline typename std::enable_if<sizeof(T) == 16, bool>::type cas(T *ptr, T oldv, T newv)
{
    static_assert(sizeof(__int128_t) == 16);
    return __sync_bool_compare_and_swap((__int128_t *)ptr, *((__int128_t *)&oldv), *((__int128_t *)&newv));//__sync_bool_compare_and_swap函数是gcc提供的一个内建函数，用于实现原子操作，它会对比第一个参数的值和第二个参数的值，如果相等，则将第三个参数的值赋给第一个参数，这个过程是原子的。
}

template <class T>
inline void write_add(T *ptr, T val)//一直尝试原子写入，直至写入成功
{
    volatile T new_val, old_val;
    do
    {
        old_val = *ptr;
        new_val = old_val + val;
    } while (!cas(ptr, old_val, new_val));
}

template <class T>
inline void write_sub(T *ptr, T val)
{
    volatile T new_val, old_val;
    do
    {
        old_val = *ptr;
        new_val = old_val - val;
    } while (!cas(ptr, old_val, new_val));
}

template <typename T>
____force_inline inline uint64_t atomic_append(std::vector<T> &array, std::atomic_uint64_t &length, const T &data, std::mutex &mutex)
{
    uint64_t idx = length.fetch_add(1, std::memory_order_acquire);
    if (idx >= array.size())
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (idx >= array.size())
        {
            if (array.size() < 1024)
                array.resize(1024);
            else
                array.resize(array.size() * 2);
        }
        std::atomic_thread_fence(std::memory_order_release);
    }
    array[idx] = data;
    return idx;
}

____force_inline inline uint64_t atomic_append(std::vector<char> &array, std::atomic_uint64_t &length, const void *data, uint64_t datalen, std::mutex &mutex)
{
    uint64_t idx = length.fetch_add(datalen, std::memory_order_acquire);
    if (idx + datalen >= array.size())
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (idx + datalen >= array.size())
        {
            if (array.size() < 1024)
                array.resize(1024);
            while (idx + datalen >= array.size())
                array.resize(array.size() * 2);
        }
        std::atomic_thread_fence(std::memory_order_release);
    }
    memcpy(array.data() + idx, data, datalen);
    return idx;
}

#endif
