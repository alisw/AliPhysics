*
* $Id$
*
* $Log$
* Revision 1.1.1.1  1996/03/06 10:46:54  mclareni
* Zebra
*
*
#ifndef CERNLIB_ZEBRA_QTRACEQ_INC
#define CERNLIB_ZEBRA_QTRACEQ_INC
*
*       without check on capacity
*
* qtraceq.inc
*
#if !defined(CERNLIB_QTRHOLL)||defined(CERNLIB_A6M)
      MQTRAC(NQTRAC+1) = NAMESR
      NQTRAC = NQTRAC + 1
#else
      MQTRAC(NQTRAC+1) = NAMESR(1)
      MQTRAC(NQTRAC+2) = NAMESR(2)
      NQTRAC = NQTRAC + 2
#endif

#endif
