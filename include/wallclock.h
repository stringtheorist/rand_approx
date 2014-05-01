#ifdef __MIC__
#define WALLCLOCK(time) do {                                 \
      unsigned long val;                                       \
      volatile unsigned int a, d;                              \
      __asm__ __volatile__("rdtsc" : "=a" (a), "=d" (d) : );   \
      val = ((unsigned long) a)|(((unsigned long)d)<<32);      \
      (time) = val / 1090000000.;                              \
    } while(0)
#else
#define WALLCLOCK(time) do {                                 \
      unsigned long val;                                       \
      volatile unsigned int a, d;                              \
      __asm__ __volatile__("rdtsc" : "=a" (a), "=d" (d) : );   \
      val = ((unsigned long) a)|(((unsigned long)d)<<32);      \
      (time) = val / 3330000000.;                              \
    } while(0)
#endif
