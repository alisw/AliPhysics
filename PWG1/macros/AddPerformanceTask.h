#ifndef ADDPERFORMANCETASK_H
#define ADDPERFORMANCETASK_H

#define TPCBIT(n)      (1 << (n))
#define SETTPCBIT(n,i) ((n) |= TPCBIT(i))
#define TSTTPCBIT(n,i) ((Bool_t)(((n) & TPCBIT(i)) != 0))
#define CLRTPCBIT(n,i) ((n) &= ~TPCBIT(i))

#define NTPCPERFORMANCE 5
#define NTPCCALIBRATION 0
const Int_t NTPCTASKS = NTPCPERFORMANCE+NTPCCALIBRATION;

Char_t *fgkTPCtaskClassName[NTPCTASKS] = {
"AliPerformanceEff"
,"AliPerformanceRes"
,"AliPerformanceTPC"
,"AliPerformanceDEdx"
,"AliPerformanceDCA"
};

const Char_t *fgkTPCtaskOpt[NTPCTASKS+1] = {
"EFF"
,"RES"
,"TPC"
,"DEDX"
,"DCA"
,"ALL"
};

Bool_t fHpt = kFALSE;  // activated with option "HPT"

enum TPCAnalysisMode {
  kTPCMode = 0,
  kTPCITSMode = 1,
  kTPCConstrMode = 2,
  kTPCInnerMode = 3
};


#endif
