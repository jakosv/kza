#ifndef TIMER_H_SENTRY
#define TIMER_H_SENTRY

#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif

void timestamp(struct timeval *time);
void print_timer(struct timeval *t_start, const char *msg);

#ifdef __cplusplus
}
#endif

#endif
