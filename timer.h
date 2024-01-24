#ifdef TIMER_H_SENTRY
#define TIMER_H_SENTRY

#include <sys/time.h>

void timestamp(struct timeval *time);
void print_timer(struct timeval *t_start, const char *msg);

#endif
