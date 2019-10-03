/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//Utils for identified fragmentation function (IDFF) analysis
//Author: Xianguo Lu (xianguo.lu@cern.ch)

#include "AliIDFFUtils.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"

#include "TClonesArray.h"

ClassImp(AliIDFFUtils);

AliPIDResponse * AliIDFFUtils::fPid=0x0;

Int_t AliIDFFUtils::PDG2Type(const Int_t pdg)
{
  //
  //conversion from pdg code to local definition of particle type
  //

  Int_t itype = kNOTSELECTED;
  switch(pdg){
  case 11:
    itype = kELECTRON;
    break;
  case 211:
    itype = kPION;
    break;
  case 2212:
    itype = kPROTON;
    break;
  case 321:
    itype = kKAON;
    break;
  default:
    itype = kNOTSELECTED;
    break;
  }
  return itype;
}

THnSparseD *AliIDFFUtils::GetTHn(const TString name)
{
  //
  //get THnSparseD
  //

  const Int_t nvar = 11;
  //                                   0       1              2              3               4     5      6       7       8         9     10   
  const TString atitle[nvar]={"TrackEta","JetPt", "TrackTPCsig", "Log10TrackP", "Log10TrackPt",  "z",  "xi",  "pdg",  "comb",   "tof",  "tpc"};
  const Int_t nbins[nvar]   ={         4,     15,          1200,          Nx(),             50,   30,    60,      7,      7,        7,      7};
  const Double_t xmins[nvar]={         0,      5,             0,        Xmin(),         Xmin(),    0,     0,   -3.5,   -3.5,     -3.5,   -3.5};
  const Double_t xmaxs[nvar]={       0.9,     20,           200,        Xmax(),         Xmax(),  1.2,     6,    3.5,    3.5,      3.5,    3.5};

  THnSparseD * hh = new THnSparseD(name,"", nvar, nbins, xmins, xmaxs);
  for(Int_t ia=0; ia<nvar; ia++){
    hh->GetAxis(ia)->SetTitle(atitle[ia]);
  }

  TAxis * ax = 0x0;
  Int_t nb = 0;

  //0: eta bin
  ax = hh->GetAxis(0);
  nb = ax->GetNbins();
  const Double_t etas[]={0, 0.2, 0.4, 0.6, 0.9};
  ax->Set(nb, etas);

  return hh;
}

/*
Bool_t AliIDFFUtils::HMPIDAcceptance(const AliAODTrack *track)
{
  //
  //check HMPID acceptance
  //From S. Pochybova
  //

  Double_t tEta = TMath::Abs(track->Eta());
  Double_t tPhi = track->Phi();
  if(tPhi < 0.) tPhi += TMath::TwoPi();
  if(tPhi > TMath::TwoPi()) tPhi -= TMath::TwoPi();

  if(tEta > 0.46){
    return kFALSE;
  }

  if(tPhi < 0.08 || tPhi > 1.12){
    return kFALSE;
  }

  return kTRUE;
}

Bool_t AliIDFFUtils::HMPIDQA(const AliAODTrack *track)
{
  //
  //check HMPID PID quality
  //From S. Pochybova
  //

  if(track->GetHMPIDsignal() <= 0.){
    return kFALSE;
  }
  
  //check track-quality cuts
  //dist_(mip-trk)
  //track variables
  Float_t tX, tY, tTh, tPh;
  //mip variables
  Float_t mpX, mpY;
  Int_t mpQ, mpNph;
  track->GetHMPIDtrk(tX, tY, tTh, tPh);
  track->GetHMPIDmip(mpX, mpY, mpQ, mpNph);
  const Double_t dist = TMath::Sqrt((tX-mpX)*(tX-mpX)+(tY-mpY)*(tY-mpY));

  //taking the pass 2 case for the moment
  if(dist > 1){
    //Printf("Track did not pass the distance cut");
    return kFALSE;
  }

  //cut on charge
  //have to check if this also varies with pass
  //taking the pass 2 case for the moment
  if(mpQ < 120){
    //Printf("Track did not pass the MipQ cut");
    return kFALSE;
  }

  return kTRUE;
}

Int_t AliIDFFUtils::HMPIDType(const AliAODTrack * track)
{
  //
  //return the (locally defined) particle type judged by HMPID
  //From S. Pochybova
  //

  if(!HMPIDAcceptance(track))
    return kNOTACCEPTED;

  if(!HMPIDQA(track))
    return kNOINFO;

  Double_t nsigma[]={-999,-999,-999};
  nsigma[kPION]     = fPid->NumberOfSigmasHMPID( track, AliPID::kPion);
  nsigma[kKAON]     = fPid->NumberOfSigmasHMPID( track, AliPID::kKaon);
  nsigma[kPROTON]   = fPid->NumberOfSigmasHMPID( track, AliPID::kProton);

  const Double_t inclusion=2;
  const Double_t exclusion=3;

  if(TMath::Abs(nsigma[kPION])<inclusion && TMath::Abs(nsigma[kKAON])>exclusion && TMath::Abs(nsigma[kPROTON])>exclusion) return kPION;
  if(TMath::Abs(nsigma[kPION])>exclusion && TMath::Abs(nsigma[kKAON])<inclusion && TMath::Abs(nsigma[kPROTON])>exclusion) return kKAON;
  if(TMath::Abs(nsigma[kPION])>exclusion && TMath::Abs(nsigma[kKAON])>exclusion && TMath::Abs(nsigma[kPROTON])<inclusion) return kPROTON;

  return kNOTSELECTED;
}
*/


Int_t AliIDFFUtils::TPCType(const AliAODTrack * trackptr)
{
  //
  //return the (locally defined) particle type judged by TPC
  //use tofmode for TPC mode
  //

  const AliPIDResponse::EDetPidStatus tpcstatus =  fPid->CheckPIDStatus(AliPIDResponse::kTPC, trackptr);
  if(tpcstatus != AliPIDResponse::kDetPidOk)
    return kNOINFO;

  Double_t nsigma[]={-999,-999,-999, -999};
  nsigma[kPION]     = fPid->NumberOfSigmasTPC( trackptr, AliPID::kPion);
  nsigma[kKAON]     = fPid->NumberOfSigmasTPC( trackptr, AliPID::kKaon);
  nsigma[kPROTON]   = fPid->NumberOfSigmasTPC( trackptr, AliPID::kProton);
  nsigma[kELECTRON] = fPid->NumberOfSigmasTPC( trackptr, AliPID::kElectron);

  //so that the effective region is really low momentum
  const Double_t inclusion=5;
  const Double_t exclusion=5;

  //don't destroy TPC signal shape below 120
  const Double_t maxsig = 150;
  if(trackptr->GetTPCsignal()> maxsig){
    if(TMath::Abs(nsigma[kPION])<inclusion && TMath::Abs(nsigma[kKAON])>exclusion && TMath::Abs(nsigma[kPROTON])>exclusion && TMath::Abs(nsigma[kELECTRON])>exclusion) return kPION;
    if(TMath::Abs(nsigma[kPION])>exclusion && TMath::Abs(nsigma[kKAON])<inclusion && TMath::Abs(nsigma[kPROTON])>exclusion && TMath::Abs(nsigma[kELECTRON])>exclusion) return kKAON;
    if(TMath::Abs(nsigma[kPION])>exclusion && TMath::Abs(nsigma[kKAON])>exclusion && TMath::Abs(nsigma[kPROTON])<inclusion && TMath::Abs(nsigma[kELECTRON])>exclusion) return kPROTON;
    if(TMath::Abs(nsigma[kPION])>exclusion && TMath::Abs(nsigma[kKAON])>exclusion && TMath::Abs(nsigma[kPROTON])>exclusion && TMath::Abs(nsigma[kELECTRON])<inclusion) return kELECTRON;
  }

  return kNOTSELECTED;
}

Int_t AliIDFFUtils::TOFType(const AliAODTrack * trackptr, const Int_t tofmode)
{
  //
  //return the (locally defined) particle type judged by TOF
  //

  //check kTOFout, kTIME, mismatch 
  const AliPIDResponse::EDetPidStatus tofstatus =  fPid->CheckPIDStatus(AliPIDResponse::kTOF, trackptr);
  if(tofstatus != AliPIDResponse::kDetPidOk)
    return kNOINFO;

  Double_t nsigma[]={-999,-999,-999, -999};
  nsigma[kPION]     = fPid->NumberOfSigmasTOF( trackptr, AliPID::kPion);
  nsigma[kKAON]     = fPid->NumberOfSigmasTOF( trackptr, AliPID::kKaon);
  nsigma[kPROTON]   = fPid->NumberOfSigmasTOF( trackptr, AliPID::kProton);
  nsigma[kELECTRON] = fPid->NumberOfSigmasTOF( trackptr, AliPID::kElectron);

  Double_t inclusion=-999;
  Double_t exclusion=-999;
  if(tofmode == 1){
    inclusion = 2;
    exclusion = 2;
  }
  else if(tofmode == 2){
    inclusion = 2;
    exclusion = 3;
  }
  else if(tofmode == 3){
    inclusion = 3;
    exclusion = 3;
  }
  else if(tofmode == 4){
    inclusion = 3;
    exclusion = 4;
  }
  else{
    printf("AliIDFFUtils::TOFType bad tofmode ! %d\n", tofmode); exit(1);
  }

  const Bool_t kpassEle = kTRUE;
  /*
  const Double_t cutEle = 1.5; 
  //tofmode = 1x then require electron exclusion cut
  if( tofmode == 4 ){
    if(TMath::Abs(nsigma[kELECTRON])> cutEle ){
      kpassEle = kTRUE;
    }
    else{
      kpassEle = kFALSE;
    }
  }
  */

  //cut on electron for pion because the precision of pion is good and the contamination of electron can not be ignored
  //+1 exclusion sigma in electron/pion to enforce better purity, otherwise not only pion, but also kaon is bias for jet pt 5-10 GeV/c
  if(TMath::Abs(nsigma[kPION])<inclusion     && TMath::Abs(nsigma[kKAON])>exclusion && TMath::Abs(nsigma[kPROTON])>exclusion && kpassEle) return kPION;
  if(TMath::Abs(nsigma[kPION])>exclusion     && TMath::Abs(nsigma[kKAON])<inclusion && TMath::Abs(nsigma[kPROTON])>exclusion && kpassEle) return kKAON;
  if(TMath::Abs(nsigma[kPION])>exclusion     && TMath::Abs(nsigma[kKAON])>exclusion && TMath::Abs(nsigma[kPROTON])<inclusion && kpassEle) return kPROTON;
  if(TMath::Abs(nsigma[kPION])>exclusion+0.5 && TMath::Abs(nsigma[kKAON])>exclusion && TMath::Abs(nsigma[kPROTON])>exclusion && TMath::Abs(nsigma[kELECTRON])<inclusion) return kELECTRON;

  return kNOTSELECTED;
}

Int_t AliIDFFUtils::CombineTPCTOF(const Int_t ktpc, const Int_t ktof)
{
  //tpc and tof, if <0 only noinfo or notselected
  if(ktpc == ktof)
    return ktpc;
  else if(ktpc < 0 && ktof >= 0 )
    return ktof;
  else if(ktof < 0 && ktpc >= 0)
    return ktpc;
  else
    return kNOTACCEPTED;
}

void AliIDFFUtils::FillTHn(THnSparseD * hh, Double_t jetpt, const AliAODTrack * track,  AliAODEvent *aodevent, const Int_t tofmode) //AliMCEvent * mcevent)
{
  //
  //fill variables
  //

  Int_t mcpdg = kNOINFO;

  TClonesArray *tca = dynamic_cast<TClonesArray*>(aodevent->FindListObject(AliAODMCParticle::StdBranchName()));
  if(tca){
    const Int_t mclabel = TMath::Abs(track->GetLabel());

    AliAODMCParticle* gentrack = dynamic_cast<AliAODMCParticle*> (tca->At(mclabel));
    if(gentrack){
      //printf("null gentrack in AliIDFFUtils::FillTHn mclabel %d\n", mclabel);
      mcpdg = PDG2Type(TMath::Abs(gentrack->GetPdgCode()));
    }
  }
  
  //===========================================================================================================

  //use tofmode for tpcmode
  const Int_t ktpc = TPCType(track);
  const Int_t ktof = TOFType(track, tofmode);

  //fake kcomb by pretpc+pretof
  const Int_t kcomb = CombineTPCTOF(ktpc, ktof);

  //const Int_t khmpid = HMPIDType(track);

  //===========================================================================================================

  const Double_t eps = 1e-6;
  const Double_t tracketa = TMath::Abs(track->Eta());
  Double_t tracktpc = track->GetTPCsignal();
  if(tracktpc>=200)
    tracktpc = 200 - eps;
  const Double_t tracklogmom = TMath::Log10(track->P());
  const Double_t tracklogpt  = TMath::Log10(track->Pt());

  Double_t zz = -999;
  Double_t xi = -999;
  if(jetpt<-990){
    jetpt = zz = xi = eps;
  }
  //from Oliver
  else if(track->Pt()>(1-eps)*jetpt && track->Pt()<(1+eps)*jetpt){ // case z=1 : move entry to last histo bin <1
    zz  = 1-eps;
    xi = eps;
  }
  else if(jetpt>0){
    zz = track->Pt()/jetpt;
    xi = TMath::Log(1/zz);
  }

  //===========================================================================================================

  const Double_t vars[]={tracketa, jetpt, tracktpc, tracklogmom, tracklogpt, zz, xi, (Double_t)mcpdg, (Double_t)kcomb, (Double_t)ktof, (Double_t)ktpc};
  
  hh->Fill(vars);
}

//_________________________________________________________________________________
Bool_t AliIDFFUtils::TPCCutPIDN(const AliAODTrack * track)
{
  //
  //TPC Cut kPIDN
  //

  //assuming this cut is particle type independent, then it is fine
  //need further investigation
  if(track->GetTPCsignalN()<60){
    return kFALSE;
  }

  return kTRUE;
}

//_________________________________________________________________________________
Bool_t AliIDFFUtils::TPCCutMIGeo(const AliAODTrack * track, const AliVEvent* evt, TTreeStream * streamer)
{
  //
  //TPC Cut MIGeo
  //

  Short_t sign = track->Charge();
  Double_t xyz[50];
  Double_t pxpypz[50];
  Double_t cv[100];

  track->GetXYZ(xyz);
  track->GetPxPyPz(pxpypz);

  AliExternalTrackParam * par = new AliExternalTrackParam(xyz, pxpypz, cv, sign);
  AliESDtrack dummy;

  Double_t varGeom = dummy.GetLengthInActiveZone(par,3,236, evt->GetMagneticField(), 0,0);
  Double_t varNcr  = track->GetTPCClusterInfo(3,1);
  Double_t varNcls = track->GetTPCsignalN();

  Bool_t cutGeom     = varGeom > 1.*(130-5*TMath::Abs(1./track->Pt()));
  Bool_t cutNcr      = varNcr  > 0.85*(130-5*TMath::Abs(1./track->Pt()));
  Bool_t cutNcls     = varNcls > 0.7*(130-5*TMath::Abs(1./track->Pt()));

  Bool_t kout = cutGeom && cutNcr && cutNcls;

  if(streamer){
    Double_t dedx = track->GetTPCsignal();

    (*streamer)<<"tree"<<
      "param.="<< par<<

      "varGeom="<<varGeom<<
      "varNcr="<<varNcr<<
      "varNcls="<<varNcls<<

      "cutGeom="<<cutGeom<<
      "cutNcr="<<cutNcr<<
      "cutNcls="<<cutNcls<<

      "kout="<<kout<<
      "dedx="<<dedx<<

      "\n";
  }

  delete par;

  return kout;
}

