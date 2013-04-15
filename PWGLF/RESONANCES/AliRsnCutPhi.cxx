#include "AliRsnCutPhi.h"

ClassImp(AliRsnCutPhi)

AliRsnCutPhi::AliRsnCutPhi() :
  AliRsnCut("cut", AliRsnTarget::kDaughter),
  fOption("")
{
   //Default constructor
  SetPhiRange(0.0, 360.0);
}

//_________________________________________________________________________________________________
AliRsnCutPhi::AliRsnCutPhi(const char *name, TString opt) :
  AliRsnCut(name, AliRsnTarget::kDaughter),
  fOption(opt.Data())
{
   //main constructor
  SetPhiRange(0.0, 360.0);
}

//_________________________________________________________________________________________________
AliRsnCutPhi::AliRsnCutPhi(const AliRsnCutPhi &copy) :
  AliRsnCut(copy),
  fOption(copy.fOption)
{
  //copy constructor
  SetPhiRange(copy.fPhiRange[0], copy.fPhiRange[1]);
}
//_________________________________________________________________________________________________
AliRsnCutPhi &AliRsnCutPhi::operator=(const AliRsnCutPhi &copy)
{
  //
  // operator =
  //
  AliRsnCut::operator=(copy);
  if (this == &copy)
    return *this;
  
  fOption=copy.fOption;
  SetPhiRange(copy.fPhiRange[0], copy.fPhiRange[1]);
  return (*this); 
}
//_________________________________________________________________________________________________
Bool_t AliRsnCutPhi::IsSelected(TObject *object)
{
  //
  // Checks if the track passes the phi cut
  //
  Bool_t accept = kFALSE;
  if (!TargetOK(object)) return accept;
  
  AliVTrack *vtrack = fDaughter->Ref2Vtrack();
  if (!vtrack) {
    AliError("Referenced daughter is not a track");
    return accept;
  }
  
  if (fOption.Contains("InTRD")) return IsInsideTRD(vtrack);
  if (fOption.Contains("OutTRD")) return IsOutsideTRD(vtrack);
  
  Double_t value = 90.0; 
  if (fOption.Contains("OuterTPC")) value = GetTrackPhi(vtrack, 278.0);
  if (fOption.Contains("InnerTOF")) value = GetTrackPhi(vtrack, 378.0);
  
  if ( (value>=fPhiRange[0]) && (value<=fPhiRange[1]) ) 
    accept = kTRUE;
  else 
    accept = kFALSE;
  
  return accept;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPhi::IsInsideTRD(AliVTrack *vtrack)
{
  //
  // Checks if track falls inside the TRD sectors
  // implemented for 2010 configuration only 
  // edge effects removed by tightening the phi cut by 5 deg 
  //
  Bool_t accept = kFALSE;
  if (!vtrack) {
    AliError("Referenced daughter is not a track");
    return accept;
  }
  Double_t value = GetTrackPhi(vtrack, 278.0);
  if ( ((value>=0.0) && (value<=35.0)) ||
       ((value>=135.0) && (value<=215.0)) ||
       ((value>=345.0) && (value<=360.0)) ) 
    accept = kTRUE;
  else 
    accept = kFALSE;
  return accept;
}
//_________________________________________________________________________________________________
Bool_t AliRsnCutPhi::IsOutsideTRD(AliVTrack *vtrack)
{
  //
  // Checks if track falls inside the TRD sectors 
  // implemented for 2010 configuration only 
  // edge effects removed by tightening the phi cut by 5 deg 
  //
  Bool_t accept = kFALSE;
  if (!vtrack) {
    AliError("Referenced daughter is not a track");
    return accept;
  }
  Double_t value = GetTrackPhi(vtrack, 278.0);
  if ( ((value>=45.0) && (value<=125.0)) ||
       ((value>=225.0) && (value<=335.0)) ) 
    accept = kTRUE;
  else 
    accept = kFALSE;
  return accept;
}

//----------------------------------------------------------------------------
Double_t AliRsnCutPhi::GetTrackPhi(AliVTrack * vtrack, Double_t radius = 0.0)
{
  //
  // Extract phi from vtrack object at radius r 
  // If r==0 (default), provides phi at vertex 
  //
  Double_t pos[3]={0.,0.,0.};
  Double_t phiOut = -999.0;
  
  if (!vtrack) {
    AliError("Invalid VTrack object");
    return phiOut;
  }
  if (radius==0.0){
    phiOut=vtrack->Phi()*TMath::RadToDeg();
  } else {
    AliExternalTrackParam etp; 
    etp.CopyFromVTrack(vtrack);
    if(etp.GetXYZAt(radius, 5., pos)){
      phiOut=TMath::ATan2(pos[1],pos[0])*TMath::RadToDeg();
      if (phiOut<0) phiOut+= (2*TMath::Pi()*TMath::RadToDeg());
    }
  }
  return phiOut;	
}
