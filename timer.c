#include "timer.h"

#include <sys/time.h>
#include <stdio.h>

void timestamp(struct timeval *time)
{
    gettimeofday(time, NULL);
}

void print_timer(struct timeval *t_start, const char *msg)
{
    enum { usecs_in_msec = 1000 };
    struct timeval t_end;
    long elapsed_time, msecs;

    timestamp(&t_end);
    elapsed_time = (t_end.tv_sec - t_start->tv_sec) * 1e6 + 
        t_end.tv_usec - t_start->tv_usec;
    msecs = elapsed_time / usecs_in_msec;
    printf("%s elapsed time %ldms\n", msg, msecs);
}
