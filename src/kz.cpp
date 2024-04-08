#include "kza.h"
#include "KZGeneric.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <functional>
#include <concepts>

#ifdef TIMER
#include "timer.h"
#endif

#include <thread>
#include <barrier>


template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZ: public KZGeneric<ValueT, SizeT, WinSizeT> {
public:
    KZ(WinSizeT window_size, const ValueT *data, SizeT data_size) : 
        KZGeneric<ValueT, SizeT, WinSizeT>(window_size, data, data_size)
    {}

private:
    using KZGeneric<ValueT, SizeT, WinSizeT>::window_size;
    using KZGeneric<ValueT, SizeT, WinSizeT>::data;
    using KZGeneric<ValueT, SizeT, WinSizeT>::ans;

    // window length is 2*window_size+1
    inline SizeT window_left_bound(SizeT win_center)
    {
        return (win_center >= SizeT(window_size))
                ? (win_center - window_size)
                : 0;
    }

    inline SizeT window_right_bound(SizeT win_center)
    {
        return (size_t(win_center + window_size) >= data.size())
                ? data.size() - 1
                : win_center + window_size;
    }

    inline void perform_single_iteration(SizeT start_idx,
                                         SizeT end_idx)
    {
        for (SizeT time = start_idx; time <= end_idx; ++time)
            ans[time] = this->average(this->window_left_bound(time),
                                      this->window_right_bound(time));
    }
};

double *kz1d(const double *data, int data_size, int win_size, int iterations)
{
    KZ<double, size_t, unsigned short> kz_data(win_size, data, data_size);

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
        ans = kz1d(data, data_size[0], window[0]>>1, iterations); 
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
