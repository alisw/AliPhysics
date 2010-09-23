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
//
// Class for TOF PID
// Implements the abstract base class AliHFEpidBase
// IsInitialized() does the PID decision
// 
// Authors:
//   Markus Fasel  <M.Fasel@gsi.de>
//   Matus Kalisky <matus.kalisky@cern.ch>  (contact)
//

#include <TList.h>
#include <TMath.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>

#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliESDpid.h"

#include "AliHFEcollection.h"
#include "AliHFEpidTOF.h"
#include "AliHFEpidBase.h"


ClassImp(AliHFEpidTOF)
  
//___________________________________________________________________
AliHFEpidTOF::AliHFEpidTOF(const Char_t *name):
  AliHFEpidBase(name)
  , fPID(0x0)
  , fQAList(0x0)
  , fNsigmaTOF(3)
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
  , fNsigmaTOF(3)
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
  target.fNsigmaTOF = fNsigmaTOF;

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
  // returns 10 (== kUnknown)if PID can not be assigned for whatever reason
  //

  if(vtrack->fAnalysisType == AliHFEpidObject::kESDanalysis){
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack *>(vtrack->fRecTrack);
    if(!esdTrack) return 0;
    AliMCParticle *mcTrack = dynamic_cast<AliMCParticle *>(vtrack->fMCtrack);
    return MakePIDesdV3(esdTrack, mcTrack);
  } else {
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(vtrack->fRecTrack);
    if(!aodTrack) return 0;
    AliAODMCParticle *aodmc = dynamic_cast<AliAODMCParticle *>(vtrack->fMCtrack);
    return MakePIDaod(aodTrack, aodmc);
  }
}

//___________________________________________________________________
Int_t AliHFEpidTOF::MakePIDesd(AliESDtrack *track, AliMCParticle * /*mcTrack*/){
  //
  // Does particle identification as discribed in IsSelected
  //
  if(!fESDpid){
    AliError("No ESD PID object available");
    return kFALSE;
  }
  Long_t status = 0;
  status = track->GetStatus(); 

  if(!(status & AliESDtrack::kTOFout)) return 0;
  
  if(IsQAon()) fQAList->Fill("hTOF_flags", 0.);

  Double_t tItrackL = track->GetIntegratedLength();
  Double_t tTOFsignal = track->GetTOFsignal();
  
  if(IsQAon()){
    if(tItrackL > 0)
      fQAList->Fill("hTOF_flags", 1.);

    if(tTOFsignal > 0)
      fQAList->Fill("hTOF_flags", 2.);
  }
  

  if(tItrackL <=0 || tTOFsignal <=0) return 0;

  if(IsQAon()){
    fQAList->Fill("hTOF_flags", 3.);
    fQAList->Fill("hTOF_signal", tTOFsignal/1000.);
    fQAList->Fill("hTOF_length", tItrackL);
  }
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

  Int_t pdg = 0;

  TString specname;
  switch(tMAXindex){
    case 0:    pdg = 11; specname = "electron";  break;
    case 1:    pdg = 13; specname = "muon"; break;
    case 2:    pdg = 211; specname = "pion"; break;
    case 3:    pdg = 321; specname = "kaon"; break;
    case 4:    pdg = 2212; specname = "proton"; break;
    default:   pdg = 0;
  };

  
  Double_t p = track->GetOuterParam()->P();
  Double_t beta = (tItrackL/100.)/(TMath::C()*(tTOFsignal/1e12));
  
  // sanity check, should not be necessary
  if(TMath::Abs(tTOFpidSum - 1) > 0.01) return 0;

  Double_t nSigmas = fESDpid->NumberOfSigmasTOF(track, (AliPID::EParticleType)tMAXindex, 0.);
  if(TMath::Abs(nSigmas) > fNsigmaTOF) return 0;

  
  // should be the same as AliPID flags
  
  if(IsQAon()){
    TString histname = "hTOFpid_" + specname;
    fQAList->Fill(histname.Data(), beta, p);
    fQAList->Fill("fTOFbeta_v_P_no", beta, p);
  }
  //return tMAXindex;
  return pdg;
  
}
//__________________________________________________________________
Int_t AliHFEpidTOF::MakePIDesdV2(AliESDtrack *track, AliMCParticle * /*mcTrack*/){
  //
  // Computes the PID response based on TOF & T0 signal
  //
  
  if(!fESDpid){
    AliError("No ESD PID object available");
    return kFALSE;
  }
  Long_t status = 0;
  status = track->GetStatus(); 
  if(!(status & AliESDtrack::kTOFpid)) return 0;

  Double_t p = track->GetOuterParam()->P();  
  // track integrated times for 5 different hypothesis and T0 time
  Double_t times[5];
  track->GetIntegratedTimes(times);
  Double_t tItrackL = track->GetIntegratedLength();
  Double_t tTOFsignal = track->GetTOFsignal();
  Double_t t0 = fESDpid->GetTOFResponse().GetTimeZero();
  //printf("-D: tof: %f, T0: %f\n", tTOFsignal, t0);
  // suppress missing or wrong T0 information
  if(t0 > 999990.0) return 0;
  Double_t tof = tTOFsignal - t0;
  Double_t beta = (tItrackL/100.)/((TMath::C()*tof)/1e12);
  //if(IsQAon())fQAList->Fill("hTOFbetaV2all", p, beta);  

  const Int_t pdg[5] = {11, 13, 211, 321, 2212};
  const Double_t invMass[5] = {TDatabasePDG::Instance()->GetParticle(11)->Mass(),
			       TDatabasePDG::Instance()->GetParticle(13)->Mass(),
			       TDatabasePDG::Instance()->GetParticle(211)->Mass(),
			       TDatabasePDG::Instance()->GetParticle(321)->Mass(),
			       TDatabasePDG::Instance()->GetParticle(2212)->Mass()};


  // accepted beta bands as function of momentum - parameters
  // line: par[0]/p + par[1] + expected_tof
  const Double_t bMin[5][2] = {{0., -0.03}, {-0.005, -0.02}, {-0.005, -0.02}, {-0.02, -0.006}, {-0.03, -0.005}}; 
  const Double_t bMax[5][2] = {{0., 0.03}, {0.005, 0.02}, {0.005, 0.02}, {0.02, 0.006}, {0.03, 0.005}};	     
	     
  Int_t index = -1;
  Double_t rdiff = 1.;
  for(Int_t i=0; i<5; ++i){
    Double_t d = (TMath::Abs(times[i] - tof))/times[i];
    if(d < rdiff){
      rdiff = d;
      index = i;
    }
  }

  // stupid and unnecessary complicated - to be improved soon
  Double_t a = p/(invMass[index]);
  a *= a;
  Double_t betaMatch = TMath::Sqrt(a/(1+a));

  // check wheter the most probable match is within allowed region of beta for given momentum and species
  Double_t min = bMin[index][0]/p + bMin[index][1] + betaMatch;
  Double_t max = bMax[index][0]/p + bMax[index][1] + betaMatch;  

  // debug
  //printf("-D: p: %f, beta: %f, pdg: %i, min: %f, max: %f\n", p, beta, pdg[index], min, max);

  //
  // PID decision - can be simpler than the QA histograms above could indicate !
  // 
  
  // suppress nonsense
  if(beta < 0.2) return 0;  

  // 1) Simple version - protect electrons
  if(beta > (1+bMin[0][1]) && beta < (1+bMax[0][1])){
    //if(IsQAon())fQAList->Fill("hTOFbetaV2electron", p, beta);
    return 11;
  }
  else return 0;
 
  
  // NOT ACTIVE when (1) activated
  // 2) more complex version - return true PID of the particle based on the best TOF estimate
  // under development - still keep protecting electrons
  if(beta > (1+bMin[0][1]) && beta < (1+bMax[0][1])){
    if(IsQAon())fQAList->Fill("hTOFbetaV2_electron", p, beta);
    return 11;
  }
  // above 3 GeV/c the supression gets weak
  if(p > 3.0) return 0;
  if(beta > min && beta < max) {
    if(IsQAon())fQAList->Fill("hTOFbetaV2selected", p, beta);
    return pdg[index];
  }
  

  return 0;
}

//___________________________________________________________________
Int_t AliHFEpidTOF::MakePIDesdV3(AliESDtrack *track, AliMCParticle * /*mctrack*/){
  //
  // TOF PID based on n-Sigma cut
  // Selects Protons and Kaons via n-sigma cut up to 3 GeV/c
  // In addition histos for n-sigma before (all species) and after (only closest species) are filled
  //
  if(!fESDpid){
    AliError("No ESD pid Object available. Return");
    return 0;
  }
  if(!(track->GetStatus() & AliESDtrack::kTOFpid)) return 0;

  Double_t t0 = fESDpid->GetTOFResponse().GetTimeZero();
  Double_t p = track->GetOuterParam() ? track->GetOuterParam()->P() : track->P();

  // Fill before selection
  Double_t sigEle = fESDpid->NumberOfSigmasTOF(track, AliPID::kElectron, t0);
  //printf("-D: p: %f, t0: %f, nSigma: %f\n", p, t0, sigEle);
  Int_t pdg = 0;
  if(TMath::Abs(sigEle) < fNsigmaTOF)
    pdg = 11;
  if(IsQAon()){
    // Draw TPC cleanup
    Double_t nsigTPC = fESDpid->NumberOfSigmasTPC(track, AliPID::kElectron);
    if(nsigTPC > -3 && nsigTPC < 5) fQAList->Fill("hTOFsigmaTPCcleanup", p, sigEle);
    // Draw TOF signal
    Double_t hcontent[3] = {p, sigEle, 0};
    THnSparseF * hptr = dynamic_cast<THnSparseF *>(fQAList->Get("hTOFsigmaElectron"));
    hptr->Fill(hcontent);
    if(pdg == 11){
      hcontent[2] = 1;
      hptr->Fill(hcontent);
    }
  }
 
  return pdg;
}
//___________________________________________________________________
Double_t AliHFEpidTOF::Likelihood(const AliESDtrack *track, Int_t species, Float_t rsig){
  
  //gives probability for track to come from a certain particle species;
  //no priors considered -> probabilities for equal abundances of all species!
  //optional restriction to r-sigma bands of different particle species; default: rsig = 2. (see below)
  
  //IMPORTANT: Tracks which are judged to be outliers get negative likelihoods -> unusable for combination with further detectors!
  
  if(!track) return -1.;
  Bool_t outlier = kTRUE;
  // Check whether distance from the respective particle line is smaller than r sigma
  for(Int_t hypo = 0; hypo < AliPID::kSPECIES; hypo++){
    if(TMath::Abs(fESDpid->NumberOfSigmasTOF(track, (AliPID::EParticleType)hypo, 0.)) > rsig)
      outlier = kTRUE;
    else {
      outlier = kFALSE;
      break;
    }
  }
  if(outlier)
    return -2.;
  
  Double_t tofProb[5];
  
  track->GetTOFpid(tofProb);
  
  return tofProb[species];
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

  fQAList = new AliHFEcollection("TOFqaHistos", "Collection for TOF PID histograms");
  //fQAList->SetName("fTOFqaHistos");
  fQAList->CreateTH1F("hTOF_flags", "TOF flags;flags (see code for info);counts", 10, -0.25, 4.75);
  fQAList->CreateTH2F("fTOFbeta_v_P_no","beta -v- P; beta;momentum [GeV/c]", 120, 0, 1.2, 200, 0, 20);
  fQAList->CreateTH1F("hTOF_signal", "TOF signal; TOF signal [ns];counts", 1000, 12, 50);
  fQAList->CreateTH1F("hTOF_length", "TOF track length; length [cm];counts", 400, 300, 700);
  fQAList->CreateTH2F("hTOFpid_electron", "TOF reco electron; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5);
  fQAList->CreateTH2F("hTOFpid_muon", "TOF reco muon; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5);
  fQAList->CreateTH2F("hTOFpid_pion", "TOF reco pion; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5);
  fQAList->CreateTH2F("hTOFpid_kaon", "TOF reco kaon; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5);
  fQAList->CreateTH2F("hTOFpid_proton", "TOF reco proton; beta ; momentum [GeV/c]", 120, 0, 1.2, 200, 0, 5);
  //fQAList->CreateTH2F("hTOFbetaV2all", "TOF #beta vs p for all tracks; momentum [GeV/c]; beta ", 400, 0.1, 10., 1200, 0, 1.2, 0);
  //fQAList->CreateTH2F("hTOFbetaV2electron", "TOF #beta vs p for selected electron tracks; momentum [GeV/c]; beta", 400, 0.1, 10., 1200, 0, 1.2, 0);

  // histograms for sigma cut
  const Int_t kNdim= 3;
  Int_t nBins[kNdim] = {1000, 1400, 2};
  Double_t binMin[kNdim] = {0.1, -12., 0};
  Double_t binMax[kNdim] = {20, 12., 2};
  fQAList->CreateTHnSparse("hTOFsigmaElectron", "TOF N#sigma around the Electron Line for all tracks; p (GeV/c); N sigma; Selection Step", kNdim, nBins, binMin, binMax);
  fQAList->CreateTH2F("hTOFsigmaTPCcleanup", "TOF N#sigma around the Electron Line for all tracks after TPC cleanup; p (GeV/c); N sigma;", 100, 0.1, 20, 140, -12, 12);

  qaList->AddLast(fQAList->GetList());
}


