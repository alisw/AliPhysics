/*
***********************************************************
    Q-vector class for event plane correction framework and analysis
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
***********************************************************
*/

#ifndef ALIEVENTPLANEQVECTOR_H
#define ALIEVENTPLANEQVECTOR_H

#include <TClonesArray.h>
#include <TBits.h>
#include <TMath.h>
#include <TList.h>
#include <iostream>
#include <Rtypes.h>
#include <TArrayS.h>
#include <THn.h>
#include <TVector.h>
#include <TArrayI.h>
#include <TProfile.h>

const Int_t fgkEPMinHarmonics = 1;
const Int_t fgkEPMaxHarmonics = 6;
const Int_t fgkNHarmonics = 6;

//_________________________________________________________________________
class AliEventPlaneQvector : public TObject {

  
 public:
        
  AliEventPlaneQvector();
  ~AliEventPlaneQvector();

  enum EventPlaneStatus {
    kRaw=0,
    kEqualized,
    kRecentered,
    kAligned,
    kDiagonalized,
    kRescaled,
    kUndefined,
    kNMaxFlowFlags
  };

  // setters
  void SetQx(Int_t harmonic, Double_t qx) { if(harmonic>=fgkEPMinHarmonics && harmonic<=fgkEPMaxHarmonics) fQvector[harmonic-fgkEPMinHarmonics][0]=(Float_t) qx;}
  void SetQy(Int_t harmonic, Double_t qy) { if(harmonic>=fgkEPMinHarmonics && harmonic<=fgkEPMaxHarmonics) fQvector[harmonic-fgkEPMinHarmonics][1]=(Float_t) qy;}
  void SetBin(Int_t bin) { fBin=bin;}
  void SetEventPlaneStatus(Int_t harmonic, EventPlaneStatus status) { 
    if(harmonic>0 && harmonic<=fgkEPMaxHarmonics) 
      fEventPlaneStatus[harmonic-fgkEPMinHarmonics] |= (1<<status);
  }
  void UnsetEventPlaneStatus(Int_t harmonic, EventPlaneStatus status) { 
    if(harmonic>0 && harmonic<=fgkEPMaxHarmonics) 
      fEventPlaneStatus[harmonic-fgkEPMinHarmonics] |= (0<<status);
  }
  void SetMultiplicity(Float_t mult)   {fMultiplicity = mult;}
  void Set(AliEventPlaneQvector* qvec);
  void Add(AliEventPlaneQvector* qvec);
  void Add(Double_t phi, Double_t w=1.);
  void Normalize();
  void SetQoverM();
  void Reset();


  // getters
  Float_t Qx(Int_t harmonic) const { return (harmonic>=fgkEPMinHarmonics && harmonic<=fgkEPMaxHarmonics ? fQvector[harmonic-fgkEPMinHarmonics][0] : 0.0 );}
  Float_t Qy(Int_t harmonic) const { return (harmonic>=fgkEPMinHarmonics && harmonic<=fgkEPMaxHarmonics ? fQvector[harmonic-fgkEPMinHarmonics][1] : 0.0 );}
  Float_t Length(Int_t harmonic) const { return (harmonic>=fgkEPMinHarmonics && harmonic<=fgkEPMaxHarmonics ? TMath::Sqrt(Qx(harmonic)*Qx(harmonic)+Qy(harmonic)*Qy(harmonic)) : 0.0 );}
  Float_t QxNorm(Int_t harmonic) const { return  (harmonic>=fgkEPMinHarmonics && harmonic<=fgkEPMaxHarmonics ? Qx(harmonic)/Length(harmonic) : 0.0 );}
  Float_t QyNorm(Int_t harmonic) const { return  (harmonic>=fgkEPMinHarmonics && harmonic<=fgkEPMaxHarmonics ? Qy(harmonic)/Length(harmonic) : 0.0 );}
  Int_t Bin() const { return fBin;}
  UChar_t GetEventPlaneStatus(Int_t h) const {return (h>=fgkEPMinHarmonics && h<=fgkEPMaxHarmonics ? fEventPlaneStatus[h-fgkEPMinHarmonics] : 999);} 
  Bool_t  CheckEventPlaneStatus(Int_t h, EventPlaneStatus flag) const;
  Float_t  Multiplicity()   const {return fMultiplicity;}
  Double_t EventPlane(Int_t h) const;



 private:
  Float_t fQvector[fgkNHarmonics][2];     // Q vector components for n harmonics
  Short_t fBin;
  UChar_t fEventPlaneStatus[fgkNHarmonics];  // Bit maps for the event plane status (1 char per detector and per harmonic)
  Float_t fMultiplicity;

     
  AliEventPlaneQvector(const AliEventPlaneQvector &c);
  AliEventPlaneQvector& operator= (const AliEventPlaneQvector &c);


  ClassDef(AliEventPlaneQvector, 1);
};




//_______________________________________________________________________________
inline Bool_t AliEventPlaneQvector::CheckEventPlaneStatus(Int_t h, EventPlaneStatus flag) const {
  //
  // Check the status of the event plane for a given detector and harmonic
  //
  if(h<fgkEPMinHarmonics || h>fgkEPMaxHarmonics) return kFALSE;
  return (flag<kNMaxFlowFlags ? (fEventPlaneStatus[h-fgkEPMinHarmonics]&(1<<flag)) : kFALSE);
}



//_______________________________________________________________________________
inline Double_t AliEventPlaneQvector::EventPlane(Int_t harmonic) const
{
  //
  // Event plane for harmonic "harmonic"
  //
  if(harmonic<1 || harmonic>fgkEPMaxHarmonics) return -999.;
  if(TMath::Abs(Qx(harmonic))<1.E-6&&TMath::Abs(Qy(harmonic))<1.E-6) return 0.;
  return TMath::ATan2(Qy(harmonic), Qx(harmonic))/Double_t(harmonic);
}

//_______________________________________________________________________________
inline void AliEventPlaneQvector::Set(AliEventPlaneQvector* qvec) 
{
  //
  // Add a Q-vector
  //

  fMultiplicity=qvec->Multiplicity();
  fBin=qvec->Bin();
  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
    fQvector[ih-1][0]=qvec->Qx(ih);
    fQvector[ih-1][1]=qvec->Qy(ih);
    fEventPlaneStatus[ih-1]=qvec->GetEventPlaneStatus(ih);
  }


}


//_______________________________________________________________________________
inline void AliEventPlaneQvector::Reset() 
{
  //
  // Add a Q-vector
  //

  fMultiplicity=0.0;
  fBin=-1;
  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
    fQvector[ih-1][0]=0.0;
    fQvector[ih-1][1]=0.0;
    fEventPlaneStatus[ih-1]=0;
  }
}


//_______________________________________________________________________________
inline void AliEventPlaneQvector::Add(AliEventPlaneQvector* qvec) 
{
  //
  // Add a Q-vector
  //

  fMultiplicity+=qvec->Multiplicity();
  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
    fQvector[ih-1][0]+=qvec->Qx(ih);
    fQvector[ih-1][1]+=qvec->Qy(ih);
  }


}


//_______________________________________________________________________________
inline void AliEventPlaneQvector::Add(Double_t phi, Double_t w) 
{
  //
  // Add a Q-vector
  //

  fMultiplicity+=w;
  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
    fQvector[ih-1][0]+=TMath::Cos(ih*phi);
    fQvector[ih-1][1]+=TMath::Sin(ih*phi);
  }


}


//_______________________________________________________________________________
inline void AliEventPlaneQvector::SetQoverM() 
{
  //
  // Set Q/M
  //

  Double_t qx,qy;
  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
    qx = Qx(ih)/Multiplicity();
    qy = Qy(ih)/Multiplicity();
    SetQx(ih,qx);
    SetQy(ih,qy);
  }
}


//_______________________________________________________________________________
inline void AliEventPlaneQvector::Normalize() 
{
  //
  // Normalize Q-vector
  //

  Double_t qx,qy;
  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
    qx = QxNorm(ih);
    qy = QyNorm(ih);
    SetQx(ih,qx);
    SetQy(ih,qy);
  }
}

   
#endif
