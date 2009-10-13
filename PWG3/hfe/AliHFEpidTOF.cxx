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

#include <TH2F.h>
#include <TList.h>
#include <TMath.h>

#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliPID.h"

#include "AliHFEpidTOF.h"
#include "AliHFEpidBase.h"


ClassImp(AliHFEpidTOF)
  
//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF(const Char_t *name):
  AliHFEpidBase(name)
  , fPID(0x0)
  , fQAList(0x0)
{
  //
  // Constructor
  //
}
//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF(const AliHFEpidTOF &c):
  AliHFEpidBase("")
  , fPID(0x0)
  , fQAList(0x0)
{  
  // 
  // Copy operator
  //

  c.Copy(*this);
}
//___________________________________________________________________
AliHFEpidTOF &AliHFEpidTOF::operator=(const AliHFEpidTOF &ref){
  //
  // Assignment operator
  //

  if(this != &ref){
    ref.Copy(*this);
  }

  return *this;
}
//___________________________________________________________________
AliHFEpidTOF::~AliHFEpidTOF(){
  //
  // Destructor
  //
  if(fPID) delete fPID;
  if(fQAList){
    fQAList->Delete();
    delete fQAList;
  }
}
//___________________________________________________________________
void AliHFEpidTOF::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEpidTOF &target = dynamic_cast<AliHFEpidTOF &>(ref);

  target.fPID = fPID;          
  target.fQAList = fQAList;

  AliHFEpidBase::Copy(ref);
}
//___________________________________________________________________
Bool_t AliHFEpidTOF::InitializePID(){
  //
  // InitializePID: TOF experts have to implement code here
  //
  return kTRUE;
}

//___________________________________________________________________
Int_t AliHFEpidTOF::IsSelected(AliHFEpidObject *vtrack)
{

  //
  // as of 22/05/2006 :
  // returns AliPID based on the ESD TOF PID decision
  // the ESD PID will be checked and if necessary improved 
  // in the mean time. Best of luck
  //
  // returns 10 (== kUnknown)if PID can not be assigned
  //
  if(vtrack->fAnalysisType == AliHFEpidObject::kESDanalysis){
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(vtrack->fRecTrack);
    if(!esdTrack) return 0;
    AliMCParticle *mcTrack = dynamic_cast<AliMCParticle *>(vtrack->fMCtrack);
    return MakePIDesd(esdTrack, mcTrack);
  } else {
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(vtrack->fRecTrack);
    if(!aodTrack) return 0;
    AliAODMCParticle *aodmc = dynamic_cast<AliAODMCParticle *>(vtrack->fMCtrack);
    return MakePIDaod(aodTrack, aodmc);
  }
}

//___________________________________________________________________
Int_t AliHFEpidTOF::MakePIDesd(AliESDtrack *track, AliMCParticle * /*mcTrack*/){
  Long_t status = 0;
  status = track->GetStatus(); 

  if(!(status & AliESDtrack::kTOFout)) return AliPID::kUnknown;
  
  (dynamic_cast<TH1F *>(fQAList->At(kHistTOFpidFlags)))->Fill(0.);

  Double_t tItrackL = track->GetIntegratedLength();
  Double_t tTOFsignal = track->GetTOFsignal();
  
  if(tItrackL > 0)
    (dynamic_cast<TH1F *>(fQAList->At(kHistTOFpidFlags)))->Fill(1.);

  if(tTOFsignal > 0)
    (dynamic_cast<TH1F *>(fQAList->At(kHistTOFpidFlags)))->Fill(2.);
  

  if(tItrackL <=0 || tTOFsignal <=0) return AliPID::kUnknown;

  (dynamic_cast<TH1F *>(fQAList->At(kHistTOFpidFlags)))->Fill(3.);
  (dynamic_cast<TH1F *>(fQAList->At(kHistTOFsignal)))->Fill(tTOFsignal/1000.);
  (dynamic_cast<TH1F *>(fQAList->At(kHistTOFlength)))->Fill(tItrackL);
  // get the TOF pid probabilities
  Double_t tESDpid[5] = {0., 0., 0., 0., 0.};
  Float_t tTOFpidSum = 0.;
  // find the largest PID probability
  track->GetTOFpid(tESDpid);
  Double_t tMAXpid = 0.;
  Int_t tMAXindex = -1;
  for(Int_t i=0; i<5; ++i){
    tTOFpidSum += tESDpid[i];
    if(tESDpid[i] > tMAXpid){
      tMAXpid = tESDpid[i];
      tMAXindex = i;
    }
  }
  
  Double_t p = track->GetOuterParam()->P();
  Double_t beta = (tItrackL/100.)/(TMath::C()*(tTOFsignal/1e12));
  
  if(TMath::Abs(tTOFpidSum - 1) > 0.01) return AliPID::kUnknown;
  else{
    // should be the same as AliPID flags
    
    (dynamic_cast<TH2F *>(fQAList->At(kHistTOFpid0+tMAXindex)))->Fill(beta, p);
    (dynamic_cast<TH2F *>(fQAList->At(kHistTOFpidBetavP)))->Fill(beta, p);
    return tMAXindex;
  }
}

//___________________________________________________________________
Int_t AliHFEpidTOF::MakePIDaod(AliAODTrack * /*aodTrack*/, AliAODMCParticle * /*mcTrack*/){
  AliError("AOD PID not yet implemented");
  return 0;
}

//___________________________________________________________________
void AliHFEpidTOF::AddQAhistograms(TList *qaList){
  //
  // Create QA histograms for TOF PID
  //

  fQAList = new TList;
  fQAList->SetName("fTOFqaHistos");
  fQAList->AddAt(new TH1F("hTOF_flags", "TOF flags;flags (see code for info);counts", 10, -0.25, 4.75), kHistTOFpidFlags);
  fQAList->AddAt(new TH2F("fTOFbeta_v_P_no","beta -v- P; beta;momentum [GeV/c]", 120, 0, 1.2, 200, 0, 20), kHistTOFpidBetavP);
  fQAList->AddAt(new TH1F("hTOF_signal", "TOF signal; TOF signal [ns];counts", 1000, 12, 50), kHistTOFsignal);
  fQAList->AddAt(new TH1F("hTOF_length", "TOF track length; length [cm];counts", 400, 300, 700), kHistTOFlength);
  fQAList->AddAt(new TH2F("hTOFpid_electron", "TOF reco electron; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5), kHistTOFpid0);
  fQAList->AddAt(new TH2F("hTOFpid_muon", "TOF reco muon; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5), kHistTOFpid1);
  fQAList->AddAt(new TH2F("hTOFpid_pion", "TOF reco pion; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5), kHistTOFpid2);
  fQAList->AddAt(new TH2F("hTOFpid_kaon", "TOF reco kaon; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5), kHistTOFpid3);
  fQAList->AddAt(new TH2F("hTOFpid_proton", "TOF reco proton; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5), kHistTOFpid4);

  qaList->AddLast(fQAList);
}
