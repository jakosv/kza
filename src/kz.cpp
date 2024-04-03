#include "kza.h"
#include "KZGeneric.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <functional>

#ifdef TIMER
#include "timer.h"
#endif

#include <thread>
#include <barrier>


template<class ValueType, class SizeType>
class KZ: public KZGeneric<ValueType, SizeType> {
public:
    KZ(SizeType window_size, const ValueType *data, SizeType data_size) : 
            KZGeneric<ValueType, SizeType>(window_size, data, data_size)
    {}

private:
    using KZGeneric<ValueType, SizeType>::window_size;
    using KZGeneric<ValueType, SizeType>::data;

    // window length is 2*window_size+1
    inline SizeType window_left_bound(SizeType win_center)
    {
        return (win_center >= window_size) ? (win_center - window_size) : 0;
    }

    inline SizeType window_right_bound(SizeType win_center)
    {
        return ((win_center + window_size) >= data.size()) ?
                                            data.size() - 1 :
                                            win_center + window_size;
    }

    inline void perform_single_iteration(SizeType start_idx, SizeType end_idx)
    {
        for (SizeType time = start_idx; time <= end_idx; ++time)
            this->ans[time] = this->average(this->window_left_bound(time),
                                            this->window_right_bound(time));
    }
};

double *kz1d(const double *data, int data_size, int win_size, int iterations)
{
    KZ<double, int> kz_data(win_size, data, data_size);

    kz_data.perform_iterations(iterations);

    return kz_data.get_ans();
}

double *kz(const double *data, int dimention, const int *data_size,
           const int *window, int iterations)
{
    double *ans = NULL;

#ifdef TIMER
    Timer timer;
#endif

    switch (dimention) {
    case 1:
#ifdef TIMER
        timer.start();
#endif
        ans = kz1d(data, data_size[0], 0.5 * window[0], iterations); 
#ifdef TIMER
        timer.stop();
        timer.print_elapsed("kz(): ");
#endif
        break;
    default:
        fprintf(stderr, "kza: Too many dimensions\n");
        break;
    }

    return ans;
}

void kz_free(double *data)
{
    delete[] data;
}
