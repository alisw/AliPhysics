#include "TObject.h"
#include "TMath.h"
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
  fIsCore(kFALSE),
  fUseCoreEnergy(kFALSE),
  fMinDistBC(0.)
{
  //Constructor
  SetName(name);

}
//________________________________________________________________________
AliPHOSClusterCuts::~AliPHOSClusterCuts()
{
  //destructor


}
//________________________________________________________________________
Bool_t AliPHOSClusterCuts::AcceptPhoton(AliCaloPhoton *ph)
{
  if(!IsNeutral(ph))   return kFALSE;
  if(!AcceptDisp(ph))  return kFALSE;
  if(!IsFarFromBC(ph)) return kFALSE;

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
  Double_t e     = ph->Energy();
  if(fUseCoreEnergy) e = (ph->GetMomV2())->Energy();

  Int_t Ncell    = ph->GetNCells();
  Double_t M02   = ph->GetLambda2();

  if(fUseDisp){
    Double_t Nmax = -1.20 + 10.2 * TMath::Sqrt(e);
    Double_t M02max = TMath::Min(3.7 , TMath::Min(1.75 + 2.73 * e , 1.94 + 1.86/TMath::Sqrt(e)));

    if(Nmax < Ncell) return kFALSE;
    //if(M02max < M02) return kFALSE;

    Double_t Nmin   = -999;
    Double_t M02min = -999;
    if(e > 1.){//enegy > 1 GeV
      M02min = TMath::Min(0.483 + 0.037 * e , 1.2);
      Nmin   = -3.04 + 3.14 * TMath::Sqrt(e);
    }

    if(fIsCore){
      if(e < 1.5){//Ecluster < 1.5 GeV
        if(M02 < M02min || M02max < M02) return kFALSE;
        if(Ncell < Nmin || Nmax < Ncell) return kFALSE;
      }
      else{
        sigma = ph->GetNsigmaCoreDisp();
        if(sigma > fNsigmaDisp) return kFALSE;
      }
    }
    else{
      if(e < 1.5){//Ecluster < 1.5 GeV
        if(M02 < M02min || M02max < M02) return kFALSE;
        if(Ncell < Nmin || Nmax < Ncell) return kFALSE;
      }
      else{
        sigma = ph->GetNsigmaFullDisp();
        if(sigma > fNsigmaDisp) return kFALSE;
      }
    }
  }

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliPHOSClusterCuts::IsFarFromBC(AliCaloPhoton *ph)
{
  Double_t distance = ph->DistToBadfp();
  if(distance < fMinDistBC) return kFALSE;
  else return kTRUE;
}
//________________________________________________________________________
Bool_t AliPHOSClusterCuts::AcceptElectron(AliCaloPhoton *ph)
{
  if(IsNeutral(ph))    return kFALSE;
  if(!AcceptDisp(ph))  return kFALSE;
  if(!IsFarFromBC(ph)) return kFALSE;

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliPHOSClusterCuts::AcceptChargedParticle(AliCaloPhoton *ph)
{
  if(!IsFarFromBC(ph)) return kFALSE;
  else return !IsNeutral(ph);
}
//________________________________________________________________________

