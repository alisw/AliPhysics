*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/04/01 15:02:54  mclareni
* Mathlib gen
*
*
*
* v107ri.inc
*
#if defined(CERNLIB_ALLIANT)||defined(CERNLIB_DECS)||defined(CERNLIB_HPUX)||defined(CERNLIB_IBMRT)||defined(CERNLIB_SGI)||defined(CERNLIB_SUN)
      MR0  = IEOR (ISHFT(MCGX,-15), MCGX)
      MR1  = ISHFT (MR0, 17)
      MCGX = IEOR (MR0, MR1)
      MCGN = 69069 * MCGN
      IRNDM2 = ISHFT (IEOR(MCGN,MCGX), -1)

#elif defined(CERNLIB_APOLLO)
      MR0  = XOR (RSHFT(MCGX,+15), MCGX)
      MR1  = LSHFT (MR0, 17)
      MCGX = XOR (MR0, MR1)
      MCGN = 69069 * MCGN
      IRNDM2 = RSHFT (XOR(MCGN,MCGX), +1)

#elif 1
      MR0  = IOR (ISHFT(MCGX,-15), MCGX)
      MR1  = ISHFT (MR0, 17)
      MCGX = IOR (MR0, MR1)
      MCGN = 69069 * MCGN
      IRNDM2 = ISHFT (IOR(MCGN,MCGX), -1)

#endif
