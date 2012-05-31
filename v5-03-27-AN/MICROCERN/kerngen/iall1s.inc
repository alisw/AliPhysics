C    set a word of all ones
#if defined(CERNLIB_QMDOS) && !defined(CERNLIB_QFMSOFT)
      PARAMETER   (IALL11 = NOT(0))
#else
      PARAMETER   (IALL11 = -1)
#endif
