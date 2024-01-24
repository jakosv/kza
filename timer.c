#include "timer.h"

#include <stdio.h>
#include <sys/time.h>

void timestamp(struct timeval *time)
{
    gettimeofday(time, NULL);
}

void print_timer(struct timeval *t_start, const char *msg)
{
    enum { usecs_in_sec = 1000000, usecs_in_msec = 1000 };
    struct timeval t_end;
    long elapsed_time, secs, msecs;

    gettimeofday(&t_end, NULL);
    elapsed_time = (t_end.tv_sec - t_start->tv_sec) * 1e6 + 
        t_end.tv_usec - t_start->tv_usec;
    secs = elapsed_time / usecs_in_sec;
    msecs = elapsed_time % usecs_in_sec;
    msecs /= usecs_in_msec;
    printf("%s: elapsed time %lds, %ldms (%ld usec)\n", 
           msg, secs, msecs, elapsed_time);
}
