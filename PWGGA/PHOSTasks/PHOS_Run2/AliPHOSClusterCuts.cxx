#include "TObject.h"
#include "AliLog.h"
#include "AliCaloPhoton.h"
#include "AliPHOSClusterCuts.h"

// Author: Daiki Sekihata (Hiroshima University)

ClassImp(AliPHOSClusterCuts)

//________________________________________________________________________
AliPHOSClusterCuts::AliPHOSClusterCuts(const char *name):
  fUseCPV(kFALSE),
  fUseDisp(kFALSE),
  fNsigmaCPV(-1),
  fNsigmaDisp(-1),
  fIsCore(kFALSE)
{
  //Constructor
  SetName(name);
  AliInfo(Form("Cut Parameter | CPV:%2.1f sigma , Dispersion:%2.1f sigma.",fNsigmaCPV,fNsigmaDisp));

}
//________________________________________________________________________
AliPHOSClusterCuts::~AliPHOSClusterCuts()
{
  //destructor


}
//________________________________________________________________________
Bool_t AliPHOSClusterCuts::AcceptPhoton(AliCaloPhoton *ph)
{
  if(!IsNeutral(ph))  return kFALSE;
  if(!AcceptDisp(ph)) return kFALSE;

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliPHOSClusterCuts::IsNeutral(AliCaloPhoton *ph)
{
  //be careful! criteria of CPV and Disp is opposite sign of inequality.
  Double_t sigma = 999;

  if(fUseCPV){
    sigma = ph->GetNsigmaCPV();
    if(sigma < fNsigmaCPV) return kFALSE;
  }

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliPHOSClusterCuts::AcceptDisp(AliCaloPhoton *ph)
{
  //be careful! criteria of CPV and Disp is opposite sign of inequality.
  Double_t sigma = 999;

  if(fUseDisp){
    if(fIsCore){
      sigma = ph->GetNsigmaCoreDisp();
      if(sigma > fNsigmaDisp) return kFALSE;
    }
    else{
      sigma = ph->GetNsigmaFullDisp();
      if(sigma > fNsigmaDisp) return kFALSE;
    }
  }

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliPHOSClusterCuts::AcceptChargedParticle(AliCaloPhoton *ph)
{
  return !IsNeutral(ph);
}
//________________________________________________________________________
