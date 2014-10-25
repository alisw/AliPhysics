#ifndef LHAPDFfw_H
#define LHAPDFfw_H

// Forward declarations of signatures of Fortran
// intermediate wrapper functions.

#include "LHAPDF/FortranWrappers.h"

extern "C" {

  #define fgetprefixpath FC_FUNC(getprefixpath, GETPREFIXPATH)
  void fgetprefixpath(char*, int len);
  #define fgetindexpath FC_FUNC(getindexpath, GETINDEXPATH)
  void fgetindexpath(char*, int len);
  #define fgetdirpath FC_FUNC(getdirpath, GETDIRPATH)
  void fgetdirpath(char*, int len);

  #define finitlhapdf FC_FUNC(initlhapdf, INITLHAPDF)
  void finitlhapdf();

  #define fgetlhapdfversion FC_FUNC(getlhapdfversion, GETLHAPDFVERSION)
  void fgetlhapdfversion(char*, int len);

  #define fgetmaxnumsets FC_FUNC(getmaxnumsets, GETMAXNUMSETS)
  void fgetmaxnumsets(int* len);

  #define finitpdfset FC_FUNC(finitpdfset, FINITPDFSET)
  void finitpdfset(char*, int len);
  #define finitpdfsetbyname FC_FUNC(finitpdfsetbyname, FINITPDFSETBYNAME)
  void finitpdfsetbyname(char*, int len);
  #define finitpdf FC_FUNC(finitpdf, FINITPDF)
  void finitpdf(int*);
  #define fevolvepdf FC_FUNC(fevolvepdf, FEVOLVEPDF)
  void fevolvepdf(double*, double *, double*);
  #define fevolvepdfp FC_FUNC(fevolvepdfp, FEVOLVEPDFP)
  void fevolvepdfp(double*, double *, double*, int*, double*);
  #define fevolvepdfa FC_FUNC(fevolvepdfa, FEVOLVEPDFA)
  void fevolvepdfa(double*, double *, double *, double*);
  #define fevolvepdfphoton FC_FUNC(fevolvepdfphoton, FEVOLVEPDFPHOTON)
  void fevolvepdfphoton(double*, double *, double*, double*);
  #define fhasphoton FC_FUNC(fhasphoton, FHASPHOTON)
  void fhasphoton(int*);
  #define fnumberpdf FC_FUNC(fnumberpdf, FNUMBERPDF)
  void fnumberpdf(int*);
  #define falphaspdf FC_FUNC(falphaspdf, FALPHASPDF)
  void falphaspdf(double*, double *);
  #define fgetorderpdf FC_FUNC(fgetorderpdf, FGETORDERPDF)
  void fgetorderpdf(int*);
  #define fgetorderas FC_FUNC(fgetorderas, FGETORDERAS)
  void fgetorderas(int*);
  #define fgetdesc FC_FUNC(fgetdesc, FGETDESC)
  void fgetdesc();
  #define fgetqmass FC_FUNC(fgetqmass, FGETQMASS)
  void fgetqmass(int*, double*);
  #define fgetthreshold FC_FUNC(fgetthreshold, FGETTHRESHOLD)
  void fgetthreshold(int*, double*);
  #define fgetnf FC_FUNC(fgetnf, FGETNF)
  void fgetnf(int*);
  #define fgetlam4 FC_FUNC(fgetlam4, FGETLAM4)
  void fgetlam4(int*, double*);
  #define fgetlam5 FC_FUNC(fgetlam5, FGETLAM5)
  void fgetlam5(int*, double*);
  #define fgetxmin FC_FUNC(fgetxmin, FGETXMIN)
  void fgetxmin(int*, double*);
  #define fgetxmax FC_FUNC(fgetxmax, FGETXMAX)
  void fgetxmax(int*, double*);
  #define fgetq2min FC_FUNC(fgetq2min, FGETQ2MIN)
  void fgetq2min(int*, double*);
  #define fgetq2max FC_FUNC(fgetq2max, FGETQ2MAX)
  void fgetq2max(int*, double*);
  #define fgetminmax FC_FUNC(fgetminmax, FGETMINMAX)
  void fgetminmax(int*, double*, double*, double*, double*);
  #define fextrapolate FC_FUNC(fextrapolate, FEXTRAPOLATE)
  void fextrapolate();

  // v5 subroutines for multiple set initialization
  #define finitpdfsetm FC_FUNC(finitpdfsetm, FINITPDFSETM)
  void finitpdfsetm(int*, char*, int len);
  #define finitpdfsetbynamem FC_FUNC(finitpdfsetbynamem, FINITPDFSETBYNAMEM)
  void finitpdfsetbynamem(int*, char*, int len);
  #define finitpdfm FC_FUNC(finitpdfm, FINITPDFM)
  void finitpdfm(int*, int*);
  #define fevolvepdfm FC_FUNC(fevolvepdfm, FEVOLVEPDFM)
  void fevolvepdfm(int*, double*, double *, double*);
  #define fevolvepdfpm FC_FUNC(fevolvepdfpm, FEVOLVEPDFPM)
  void fevolvepdfpm(int*, double*, double *, double*, int*, double*);
  #define fevolvepdfam FC_FUNC(fevolvepdfam, FEVOLVEPDFAM)
  void fevolvepdfam(int*, double*, double *, double *, double*);
  #define fevolvepdfphotonm FC_FUNC(fevolvepdfphotonm, FEVOLVEPDFPHOTONM)
  void fevolvepdfphotonm(int*, double*, double *, double*, double*);
  #define fnumberpdfm FC_FUNC(fnumberpdfm, FNUMBERPDFM)
  void fnumberpdfm(int*, int*);
  #define falphaspdfm FC_FUNC(falphaspdfm, FALPHASPDFM)
  void falphaspdfm(int*, double*, double *);
  #define fgetorderpdfm FC_FUNC(fgetorderpdfm, FGETORDERPDFM)
  void fgetorderpdfm(int*, int*);
  #define fgetorderasm FC_FUNC(fgetorderasm, FGETORDERASM)
  void fgetorderasm(int*, int*);
  #define fgetdescm FC_FUNC(fgetdescm, FGETDESCM)
  void fgetdescm(int*);
  #define fgetqmassm FC_FUNC(fgetqmassm, FGETQMASSM)
  void fgetqmassm(int*, int*, double*);
  #define fgetthresholdm FC_FUNC(fgetthresholdm, FGETTHRESHOLDM)
  void fgetthresholdm(int*, int*, double*);
  #define fgetnfm FC_FUNC(fgetnfm, FGETNFM)
  void fgetnfm(int*, int*);
  #define fgetlam4m FC_FUNC(fgetlam4m, FGETLAM4M)
  void fgetlam4m(int*, int*, double*);
  #define fgetlam5m FC_FUNC(fgetlam5m, FGETLAM5M)
  void fgetlam5m(int*, int*, double*);
  #define fgetxminm FC_FUNC(fgetxminm, FGETXMINM)
  void fgetxminm(int*, int*, double*);
  #define fgetxmaxm FC_FUNC(fgetxmaxm, FGETXMAXM)
  void fgetxmaxm(int*, int*, double*);
  #define fgetq2minm FC_FUNC(fgetq2minm, FGETQ2MINM)
  void fgetq2minm(int*, int*, double*);
  #define fgetq2maxm FC_FUNC(fgetq2maxm, FGETQ2MAXM)
  void fgetq2maxm(int*, int*, double*);
  #define fgetminmaxm FC_FUNC(fgetminmaxm, FGETMINMAXM)
  void fgetminmaxm(int*, int*, double*, double*, double*, double*);
  #define fextrapolateon FC_FUNC(fextrapolateon, FEXTRAPOLATEON)
  void fextrapolateon();
  #define fextrapolateoff FC_FUNC(fextrapolateoff, FEXTRAPOLATEOFF)
  void fextrapolateoff();
  #define fsilent FC_FUNC(fsilent, FSILENT)
  void fsilent();
  #define flowkey FC_FUNC(flowkey, FLOWKEY)
  void flowkey();
  #define fdefaultverb FC_FUNC(fdefaultverb, FDEFAULTVERB)
  void fdefaultverb();
  #define fsetpdfpath FC_FUNC(fsetpdfpath, FSETPDFPATH)
  void fsetpdfpath(char*, int len);

  #define fsetlhaparm FC_FUNC(setlhaparm, SETLHAPARM)
  void fsetlhaparm(char*, int len);


}

#endif
