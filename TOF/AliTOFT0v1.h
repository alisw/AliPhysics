#ifndef AliTOFT0v1_H
#define AliTOFT0v1_H

#include "TString.h"
#include "AliESDEvent.h"
#include "AliStack.h"

class AliESDtrack;
class AliTOFcalibHisto;

class AliTOFT0v1: public TObject {
public:
  
  AliTOFT0v1(AliESDEvent*) ;
  AliTOFT0v1(const AliTOFT0v1 & tzero);
  virtual ~AliTOFT0v1() ; // dtor
 
  void SetCalib(AliTOFcalibHisto *calib){fCalib = calib;};

  const char*   GetTZeroFile() const {return fT0File.Data();}   
  Double_t* DefineT0(Option_t *option); 
  Double_t* DefineT0RawCorrection(Option_t *option); 
  
  void      SetTimeResolution(Double_t timeresolution);// timeresolution in [s] e.g. for 120 ps -> 1.2e-10
  
  Double_t       GetTimeResolution(){return fTimeResolution;}
  
  void          SetTZeroFile(char* file) ;
  void          SetMomBounds(Float_t pLow, Float_t pUp) { fLowerMomBound=pLow; fUpperMomBound=pUp;} // momenta are expressed in [GeV/c]
  void          SetTimeCorr(Float_t timecorr) {fTimeCorr=timecorr;} //in ns!!!
  void          SetT0Offset(Float_t t0offset){fT0Offset=t0offset;} //in ns!!!
  void          SetLOffset(Float_t loffset){fT0Offset=loffset;}  //in m!!!
  Float_t       GetMomError(Int_t index, Float_t mom, Float_t texp);
  void  Print(Option_t* option) const ;
  Bool_t   operator == (const AliTOFT0v1 & tzero) const ;

 private:

  Bool_t AcceptTrack(AliESDtrack *track); /* accept track */
  Float_t GetSigmaToVertex(AliESDtrack *track); /* get sigma to vertex */

  AliTOFcalibHisto *fCalib;

  Double_t fTimeResolution;  // global time resolution used to calculate T0
  Float_t fTimeCorr;  // global time resolution used to calculate T0
  Float_t fLowerMomBound;   // momentum lower bound for selected primary tracks   
  Float_t fT0Offset;
  Float_t fLOffset;
  Float_t fUpperMomBound;   // momentum upper bound for selected primary tracks 
  Float_t fDeltaTfromMisallinement;

  Double_t fT0SigmaT0def[4];
  AliESDEvent* fEvent;      //evento per il quale si vuole calcolare il T0
  //AliStack* fStack;         //stack associata all'evento fEvent
  TString fT0File ;         // output file; it contains for time being only 3 histos 
  
  ClassDef(AliTOFT0v1,1);  // Calculate the time zero using TOF detector */
  
};

#endif 
