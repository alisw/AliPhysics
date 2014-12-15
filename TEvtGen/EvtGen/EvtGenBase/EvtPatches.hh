#if defined(SUNCC) && defined(_XOPEN_SOURCE) && ( _XOPEN_SOURCE - 0 == 500 )
#ifndef _CLOCK_T
#define _CLOCK_T
typedef long            clock_t; /* relative time in a specified resolution */
#endif  /* ifndef _CLOCK_T */
#endif

#define _unused(x) ((void)x)
