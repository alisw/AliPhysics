

#ifndef ALIEVENTPLANECORRELATIONS_H
#define ALIEVENTPLANECORRELATIONS_H

#include <TMath.h>
#include <THn.h>
#include <TProfile.h>
#include <AliQnCorrectionsQnVector.h>
#include "AliChargeOnePwrtRP.h"

#include <iostream>
using std::cout;
using std::endl;


//_________________________________________________________________________
class AliEventPlaneCorrelations : public TObject {

 public:
  AliEventPlaneCorrelations(TAxis* ax,TString epA, TString epB, TString epC="");
  ~AliEventPlaneCorrelations();

  void SetEventPlanes(AliQnCorrectionsQnVector* qvecA,AliQnCorrectionsQnVector* qvecB,AliQnCorrectionsQnVector* qvecC) {fEventPlanes[0]=qvecA;fEventPlanes[1]=qvecB;fEventPlanes[2]=qvecC;}
  TProfile* CorrelationProfile(Int_t icor, Int_t ih, Int_t icomp) const {return fEPcorrelation[icor][ih-1][icomp];}
  TString CorrelationName() const {return fCorrelationName;}

  void FillCorrelations(Float_t xvalue);

 private:
  TProfile* fEPcorrelation[3][AliChargeOnePwrtRP::Nharmonics][4];
  AliQnCorrectionsQnVector* fEventPlanes[3];
  TString fCorrelationName;

  AliEventPlaneCorrelations(const AliEventPlaneCorrelations &c);
  AliEventPlaneCorrelations& operator= (const AliEventPlaneCorrelations &c);

  ClassDef(AliEventPlaneCorrelations, 1);
};




#endif
