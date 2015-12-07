// Class for reduced fmd information information
// Author: Jaap Onderwaater
// 

#ifndef ALIREDUCEDFMDINFO_H
#define ALIREDUCEDFMDINFO_H

#include <TMath.h>
#include <TObject.h>


//_____________________________________________________________________
class AliReducedFMDInfo : public TObject {

  friend class AliAnalysisTaskReducedTreeMaker;                                             // friend analysis task which fills the object

 public:
  AliReducedFMDInfo();
  ~AliReducedFMDInfo();

  // getters
  Float_t Multiplicity()	    const {return fMultiplicity;}
  Short_t Id()	    	        const {return fId;}                                           // id is mapped to bin number for FMD histograms (fId = iEta*nPhi+iPhi for FMD1I, FMD2I and FMD3I;  fId=-(iEta*nPhi+iPhi) for FMD2O and FMD3O)

  // The following functions are based on the ring by ring 2d-histograms for corrected FMD information from the AliForwardMultiplicityTask in PWGLF
  Float_t Phi()     const  {return ( fId>0 ? TMath::Pi()/20. + TMath::Pi()/10.*PhiBin() :  TMath::Pi()/40. + TMath::Pi()/20.*PhiBin() );} // get phi for bin number of FMD histogram
  Float_t Eta()     const  {return -3.975+0.05*(EtaBin()-1);}                                             // get eta for bin number of FMD histogram
  Int_t PhiBin()    const  {return ( fId>0 ? (fId-1)%20 : (-fId-1)%40);}                                  // get bin number along phi-axis for FMD histogram
  Int_t EtaBin()    const  {return ( fId>0 ? (fId-PhiBin())/20 : (-fId-PhiBin())/40) ;}                   // get bin number along eta-axis for FMD histogram

  Double_t RawPhi(Int_t det, Int_t nEtaSlices) const;
  Double_t RawPhi(Int_t det) const;
  
 private:

  Float_t fMultiplicity;        
  Short_t fId;           
    
  AliReducedFMDInfo(const AliReducedFMDInfo &c);
  AliReducedFMDInfo& operator= (const AliReducedFMDInfo &c);

  ClassDef(AliReducedFMDInfo, 1);
};

//_______________________________________________________________________________
inline Double_t AliReducedFMDInfo::RawPhi(Int_t det) const
{

Double_t phi=-999.;
Int_t ich=Id();//-10240*det;

switch (det){
        case 0:
          phi = TMath::Pi()/10.*(Int_t (ich/512.));  // stepsize * strip number 
          phi+=11*TMath::Pi()/20.;        // offset
        break;
        case 1:
          phi = TMath::Pi()/10.*(Int_t (ich/512.));  // stepsize * strip number
          phi+=TMath::Pi()/20.;          // offset
        break;
        case 2:
          phi = TMath::Pi()/20.*(Int_t (ich/256.));  // stepsize * strip number
          phi+=TMath::Pi()/40.;          // offset
        break;
        case 3:
          phi = TMath::Pi()/10.*(20-(Int_t (ich/512.)));  // stepsize * strip number (negative slope)
          phi+=19*TMath::Pi()/20.;         // offset
        break;
        case 4:
          phi = TMath::Pi()/20.*(40-(Int_t (ich/256.)));  // stepsize * strip number (negative slope)
          phi+=39*TMath::Pi()/40.;         // offset
        break;
}

if(phi>TMath::Pi()*2) phi = phi-TMath::Pi()*2.; // domain of (0,2 pi]
if(phi<-TMath::Pi()*2) phi = phi+TMath::Pi()*2.;

return phi;

}


//_______________________________________________________________________________
inline Double_t AliReducedFMDInfo::RawPhi(Int_t det, Int_t nEtaSlices) const
{
  Double_t phi=-999.;
  Double_t delta=0.0;
  if(nEtaSlices==0) phi=RawPhi(det);
  else{
  Int_t id=Id();
  Int_t nClus=0;
    if(det==2||det==4){   // if FMD2O or FMD3O  n sectors is 40, otherwise 20  (FMD reduced coordinates are aligned)
      nClus=40;
      delta=2.*TMath::Pi()/40.;
    }
    else {
      nClus=20;
      delta=2.*TMath::Pi()/20.;
    }
    
    phi = delta*(id%nClus+0.5);
  }
  
  return phi;
}

#endif
