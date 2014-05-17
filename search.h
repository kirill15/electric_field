#ifndef SEARCH_H
#define SEARCH_H

#include "cstddef"

template <class T>
int binarySearch(T x, size_t size, T *array, size_t step = 1)
{
    size_t first = 0;
    size_t last = size * step; // Номер элемента в массиве, СЛЕДУЮЩЕГО ЗА последним

    while (first < last)
    {
        size_t mid = first + (last - first) / 2;

        if (x <= array[mid])
            last = mid;
        else
            first = mid + step;
    }
    if (last < size && array[last] > x)
        return last;
    else
        return -1;
}



#endif // SEARCH_H
