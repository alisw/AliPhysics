#ifndef AliMultGlauberNBDFitter_H
#define AliMultGlauberNBDFitter_H

#include <iostream>
#include "TNamed.h"
#include "AliVEvent.h"
//For Run Ranges functionality
#include <map>

using namespace std;
class AliMultGlauberNBDFitter : public TNamed {
  
public:
  //basic functionality
  AliMultGlauberNBDFitter();
  AliMultGlauberNBDFitter(const char * name, const char * title = "Glauber+NBD fitter");
  ~AliMultGlauberNBDFitter();
  
  //Master fitter function
  Double_t ProbDistrib(Double_t *x, Double_t *par);

  //Do Fit: where everything happens 
  Bool_t DoFit();
  
  //Set input characteristics: the 2D plot with Npart, Nanc
  Bool_t SetNpartNcollCorrelation(TH2 *hNpNc); 
  
  //Set main input to be fitted (the V0M distribution)
  Bool_t SetInputV0M(TH1 *hV0M);
  
  //Interface to get funtions if asked to
  TF1 *GetNBD();
  TF1 *GetGlauberNBD();
  
  //Helper
  Bool_t InitializeNpNc();
  
  //Interface to set vals
  void SetMu ( Double_t lVal ) {fMu = lVal;}
  void Setk ( Double_t lVal ) {fk = lVal;}
  void Setf ( Double_t lVal ) {ff = lVal;}
  void SetNorm ( Double_t lVal ) {fnorm = lVal;}
  
  //Interface to get vals
  Double_t GetMu ()   {return fMu;}
  Double_t Getk ()    {return fk;}
  Double_t Getf ()    {return ff;}
  Double_t GetNorm () {return fnorm;}
  
  void SetFitRange  (Double_t lMin, Double_t lMax);
  void SetFitOptions(TString lOpt);
  
  //void    Print(Option_t *option="") const;
  
private:
  //This function serves as the (analytical) NBD
  TF1 *fNBD;
  
  //This function is the key fitting function
  TF1 *fGlauberNBD;
  
  //Reference histo
  TH1D *fhNanc; //basic ancestor distribution
  TH2 *fhNpNc; //correlation between Npart and Ncoll
  TH1 *fhV0M; //basic ancestor distribution
  
  //Fitting utilities
  Bool_t ffChanged;
  Double_t fCurrentf;
  
  //Buffer for (Npart, Ncoll) pairs in memory
  Double_t *fNpart;
  Double_t *fNcoll;
  Long_t *fContent;
  Long_t fNNpNcPairs; //number of pairs to use
  Long_t fMaxNpNcPairs;
  
  //The actual output: mu, k, f, norm
  Double_t fMu;
  Double_t fk;
  Double_t ff;
  Double_t fnorm;
  
  TString fFitOptions; 
  
  ClassDef(AliMultGlauberNBDFitter, 1);
};
#endif
