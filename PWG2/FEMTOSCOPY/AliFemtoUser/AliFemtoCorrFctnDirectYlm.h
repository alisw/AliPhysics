////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDirectYlm - Correlation function that is binned in Ylms    //
// directly. Provides a way to store the numerator and denominator            //
// in Ylms directly and correctly calculate the correlation                   //
// function from them.                                                        //
//                                                                            //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOCORRFCTNDIRECTYLM_H
#define ALIFEMTOCORRFCTNDIRECTYLM_H

#include <math.h>
#include <complex>
#include <TH1D.h>
#include <TH3D.h>
#include <TFile.h>
#include "AliFemtoCorrFctn.h"
#include "AliFemtoYlm.h"

using namespace std;

class AliFemtoCorrFctnDirectYlm: public AliFemtoCorrFctn {
 public:
  AliFemtoCorrFctnDirectYlm();
  AliFemtoCorrFctnDirectYlm(const char *name, int maxl, int ibin, double vmin, double vmax);
  AliFemtoCorrFctnDirectYlm(const AliFemtoCorrFctnDirectYlm& aCorrFctn);
  ~AliFemtoCorrFctnDirectYlm();

  AliFemtoCorrFctnDirectYlm& operator=(const AliFemtoCorrFctnDirectYlm& aCorrFctn);

  void AddRealPair(double *qvec, double weight=1.0);
  void AddMixedPair(double *qvec, double weight=1.0);

  void AddRealPair(double qout, double qside, double qlong, double weight=1.0);
  void AddMixedPair(double qout, double qside, double qlong, double weight=1.0);

  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();
  virtual TList* GetOutputList();

  void Write();

  void ReadFromFile(TFile *infile, const char *name, int maxl);

  TH1D *GetNumRealHist(int el, int em);
  TH1D *GetNumImagHist(int el, int em);

  TH1D *GetDenRealHist(int el, int em);
  TH1D *GetDenImagHist(int el, int em);

 private:
  double ClebschGordan(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm);
  double DeltaJ(double aJot1, double aJot2, double aJot);
  double WignerSymbol(double aJot1, double aEm1, double aJot2, double aEm2, double aJot, double aEm);
  
  void GetMtilde(complex<double>* aMat, double *aMTilde); 
  
  int  GetMaxJM() const;
  void GetElEmForIndex(int aIndex, double *aEl, double *aEm) const;
  void GetElEmForIndex(int aIndex, int *aEl, int *aEm) const;
  int  GetBin(int qbin, int ilmzero, int zeroimag, int ilmprim, int primimag);

  int  PackYlmVector(const double *invec, double *outvec);
  int  PackYlmMatrix(const double *inmat, double *outmat);

  int GetIndexForLM(int el, int em) const;

  void PackCovariances();
  void UnpackCovariances();

  TH1D **fnumsreal;            // Real parts of Ylm components of the numerator
  TH1D **fnumsimag;            // Imaginary parts of Ylm components of the numerator
  TH1D **fdensreal;            // Real parts of Ylm components of the denominator	    
  TH1D **fdensimag;            // Imaginary parts of Ylm components of the denominator

  TH1D *fbinctn;               // Bin occupation for the numerator
  TH1D *fbinctd;               // Bin occupation for the denominator

  TH3D *fcovnum;               // Numerator covariance matrix packed into TH3D
  TH3D *fcovden;               // Denominator covariance matrix packed into TH3D

  double *fcovmnum;            // Covariance matrix for the numerator
  double *fcovmden;            // Covariance matrix for the denominator

  int fMaxL;                  // l cut-off of the decomposition

  int    fMaxJM;               // number of l-m combinations
  double *fels;                // table of l's
  double *fems;                // table of m's
  int    *felsi;               // table of integer l's
  int    *femsi;               // table of integer m's

  complex<double> *fYlmBuffer; // buffer for ylm calculation
  double *factorials;         // Helper table of factorials

  double fSout;                // Save last calculated qout
  double fSside;               // Save last calculated qside
  double fSlong;               // Save last calculated qlong

};

#endif

