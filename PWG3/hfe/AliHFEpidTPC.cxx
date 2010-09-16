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
// Class for TPC PID
// Implements the abstract base class AliHFEpidBase
// 
// Class contains TPC specific cuts and QA histograms
// Two selection strategies are offered: Selection of certain value
// regions in the TPC dE/dx (by IsSelected), and likelihoods
//
// Authors: 
//
//   Markus Fasel <M.Fasel@gsi.de> 
//   Markus Heide <mheide@uni-muenster.de> 
//  
#include <TF1.h>
#include <TList.h>
#include <TMath.h>
#include <THnSparse.h>

#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliESDpid.h"

#include "AliHFEcollection.h"
#include "AliHFEpidTPC.h"

ClassImp(AliHFEpidTPC)

//___________________________________________________________________
AliHFEpidTPC::AliHFEpidTPC(const char* name) :
  // add a list here
  AliHFEpidBase(name)
  , fLineCrossingType(0)
  , fLineCrossingsEnabled(0)
  , fUpperSigmaCut(NULL)
  , fLowerSigmaCut(NULL)
  , fNsigmaTPC(3)
  , fRejectionEnabled(0)
  , fPID(NULL)
  , fQAList(NULL)
{
  //
  // default  constructor
  // 
  memset(fLineCrossingSigma, 0, sizeof(Double_t) * AliPID::kSPECIES);
  memset(fPAsigCut, 0, sizeof(Float_t) * 2);
  memset(fNAsigmaTPC, 0, sizeof(Float_t) * 2);
  fPID = new AliPID;
}

//___________________________________________________________________
AliHFEpidTPC::AliHFEpidTPC(const AliHFEpidTPC &ref) :
  AliHFEpidBase("")
  , fLineCrossingType(0)
  , fLineCrossingsEnabled(0)
  , fUpperSigmaCut(NULL)
  , fLowerSigmaCut(NULL)
  , fNsigmaTPC(2)
  , fRejectionEnabled(0)
  , fPID(NULL)
  , fQAList(NULL)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);
}

//___________________________________________________________________
AliHFEpidTPC &AliHFEpidTPC::operator=(const AliHFEpidTPC &ref){
  //
  // Assignment operator
  //
  if(this != &ref){
    ref.Copy(*this);
  } 
  return *this;
}
//___________________________________________________________________
void AliHFEpidTPC::Copy(TObject &o) const{
  //
  // Copy function 
  // called in copy constructor and assigment operator
  //
  AliHFEpidTPC &target = dynamic_cast<AliHFEpidTPC &>(o);

  target.fLineCrossingsEnabled = fLineCrossingsEnabled;
  target.fUpperSigmaCut = fUpperSigmaCut;
  target.fLowerSigmaCut = fLowerSigmaCut;
  target.fNsigmaTPC = fNsigmaTPC;
  target.fRejectionEnabled = fRejectionEnabled;
  target.fPID = new AliPID(*fPID);
  target.fQAList = new AliHFEcollection(*fQAList);
  memcpy(target.fLineCrossingSigma, fLineCrossingSigma, sizeof(Double_t) * AliPID::kSPECIES);
  memcpy(target.fPAsigCut, fPAsigCut, sizeof(Float_t) * 2);
  memcpy(target.fNAsigmaTPC, fNAsigmaTPC, sizeof(Float_t) * 2);
 
  AliHFEpidBase::Copy(target);
}

//___________________________________________________________________
AliHFEpidTPC::~AliHFEpidTPC(){
  //
  // Destructor
  //
  if(fPID) delete fPID;
  if(fQAList){
    delete fQAList;
  }
}

//___________________________________________________________________
Bool_t AliHFEpidTPC::InitializePID(){
  //
  // Add TPC dE/dx Line crossings
  //
  //AddTPCdEdxLineCrossing(AliPID::kKaon, 0.3, 0.018);
  //AddTPCdEdxLineCrossing(AliPID::kProton, 0.9, 0.054);
  return kTRUE;
}

//___________________________________________________________________
Int_t AliHFEpidTPC::IsSelected(AliHFEpidObject *track)
{
  //
  // For the TPC pid we use the 2-sigma band around the bethe bloch curve
  // for electrons
  // exclusion of the crossing points
  //
  if(track->fAnalysisType == AliHFEpidObject::kESDanalysis){
    AliESDtrack *esdTrack = NULL;
    if(!(esdTrack = dynamic_cast<AliESDtrack *>(track->fRecTrack))) return 0;
    AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(track->fMCtrack);
    return MakePIDesd(esdTrack, mctrack);
  }else{
    AliAODTrack *aodtrack = NULL;
    if(!(aodtrack = dynamic_cast<AliAODTrack *>(track->fRecTrack))) return 0;
    AliAODMCParticle *aodmctrack = dynamic_cast<AliAODMCParticle *>(track->fMCtrack);
    return MakePIDaod(aodtrack, aodmctrack);
  }
}
//___________________________________________________________________
Int_t AliHFEpidTPC::MakePIDesd(AliESDtrack *esdTrack, AliMCParticle *mctrack){
  //
  //  Doing TPC PID as explained in IsSelected for ESD tracks
  //
  if(!fESDpid){
    AliError("No ESD PID object available");
    return kFALSE;
  }
  AliDebug(1, "Doing TPC PID based on n-Sigma cut approach");
  Float_t nsigma = fESDpid->NumberOfSigmasTPC(esdTrack, AliPID::kElectron);
  if(IsQAon()) FillTPChistograms(esdTrack, mctrack, kFALSE);
  // exclude crossing points:
  // Determine the bethe values for each particle species
  Bool_t isLineCrossing = kFALSE;
  fLineCrossingType = 0;  // default value
  for(Int_t ispecies = 0; ispecies < AliPID::kSPECIES; ispecies++){
    if(ispecies == AliPID::kElectron) continue;
    if(!(fLineCrossingsEnabled & 1 << ispecies)) continue;
    if(TMath::Abs(fESDpid->NumberOfSigmasTPC(esdTrack, (AliPID::EParticleType)ispecies)) < fLineCrossingSigma[ispecies] && TMath::Abs(nsigma) < fNsigmaTPC){
      // Point in a line crossing region, no PID possible, but !PID still possible ;-)
      isLineCrossing = kTRUE;      
      fLineCrossingType = ispecies;
      break;
    }
  }
  if(isLineCrossing) return 0;

  // Check particle rejection
  if(HasParticleRejection()){
    Int_t reject = Reject(esdTrack);
    if(reject != 0) return reject;
  }

  // Check if we have an asymmetric sigma model set
  Int_t pdg = 0;
  if(fUpperSigmaCut || fLowerSigmaCut){
    pdg = CutSigmaModel(esdTrack) ? 11 : 0;
  } else { 
    // Perform Asymmetric n-sigma cut if required, else perform symmetric TPC sigma cut
    Float_t p = 0.;
    if(HasAsymmetricSigmaCut() && (p = esdTrack->P()) >= fPAsigCut[0] && p <= fPAsigCut[1]){ 
      if(nsigma >= fNAsigmaTPC[0] && nsigma <= fNAsigmaTPC[1]) pdg = 11; 
    } else {
      if(TMath::Abs(nsigma) < fNsigmaTPC ) pdg = 11;
    }
  }
  if(IsQAon() && pdg != 0) FillTPChistograms(esdTrack, mctrack, kTRUE);
  return pdg;
}

//___________________________________________________________________
Int_t AliHFEpidTPC::MakePIDaod(AliAODTrack * /*aodTrack*/, AliAODMCParticle * /*mctrack*/){
  AliError("AOD PID not yet implemented");
  return 0;
}

//___________________________________________________________________
Bool_t AliHFEpidTPC::CutSigmaModel(AliESDtrack *track){
  //
  // N SigmaCut using parametrization of the cuts
  //
  Bool_t isSelected = kTRUE;
  Float_t nsigma = fESDpid->NumberOfSigmasTPC(track, AliPID::kElectron);
  Double_t p =  track->GetInnerParam() ? track->GetInnerParam()->P() : track->P();
  if(fUpperSigmaCut && nsigma > fUpperSigmaCut->Eval(p)) isSelected = kFALSE;
  if(fLowerSigmaCut && nsigma < fLowerSigmaCut->Eval(p)) isSelected = kFALSE;
  return isSelected;
}

//___________________________________________________________________
Int_t AliHFEpidTPC::Reject(AliESDtrack *track){
  //
  // reject particles based on asymmetric sigma cut
  //
  Int_t pdc[AliPID::kSPECIES] = {11,13,211,321,2212};
  Double_t p = track->GetOuterParam() ? track->GetOuterParam()->P() : track->P();
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    if(!TESTBIT(fRejectionEnabled, ispec)) continue;
    // Particle rejection enabled
    if(p < fRejection[4*ispec] || p > fRejection[4*ispec+2]) continue;
    Double_t sigma = fESDpid->NumberOfSigmasTPC(track, static_cast<AliPID::EParticleType>(ispec));
    if(sigma >= fRejection[4*ispec+1] && sigma <= fRejection[4*ispec+3]) return pdc[ispec] * track->Charge();
  }
  return 0;
}

//___________________________________________________________________
void AliHFEpidTPC::AddTPCdEdxLineCrossing(Int_t species, Double_t sigma){
  //
  // Add exclusion point for the TPC PID where a dEdx line crosses the electron line
  // Stores line center and line sigma
  //
  if(species >= AliPID::kSPECIES){
    AliError("Species doesn't exist");
    return;
  }
  fLineCrossingsEnabled |= 1 << species;
  fLineCrossingSigma[species] = sigma;
}

//___________________________________________________________________
Double_t AliHFEpidTPC::Likelihood(const AliESDtrack *track, Int_t species, Float_t rsig)
{
  //gives probability for track to come from a certain particle species;
  //no priors considered -> probabilities for equal abundances of all species!
  //optional restriction to r-sigma bands of different particle species; default: rsig = 2. (see below)

  //IMPORTANT: Tracks which are judged to be outliers get negative likelihoods -> unusable for combination with further detectors!
  
  if(!track) return -1.;
  Bool_t outlier = kTRUE;
  // Check whether distance from the respective particle line is smaller than r sigma
  for(Int_t hypo = 0; hypo < AliPID::kSPECIES; hypo++){
    if(TMath::Abs(fESDpid->NumberOfSigmasTPC(track, (AliPID::EParticleType)hypo)) > rsig)
      outlier = kTRUE;
    else {
	    outlier = kFALSE;
	    break;
	  }
  }
  if(outlier)
    return -2.;

  Double_t tpcProb[5];

  track->GetTPCpid(tpcProb);

  return tpcProb[species];
}

//___________________________________________________________________
Double_t  AliHFEpidTPC::Suppression(const AliESDtrack *track, Int_t species)
{
  //ratio of likelihoods to be whatever species/to be an electron;
  //as a cross-check for possible particle type suppression compared to electrons
  const Double_t kVerySmall = 1e-10;
  if(!track) return -20;
  if((TMath::Abs(Likelihood(track,species) + 2) < kVerySmall)||(TMath::Abs(Likelihood(track,0) + 2 ) < kVerySmall))
    return -30;
  if(TMath::Abs(Likelihood(track,species)) < kVerySmall)
    return -10;
  if (TMath::Abs(Likelihood(track,0)) < kVerySmall)
    return 10.;
  else
    return TMath::Log10(Likelihood(track,species)/(Likelihood(track,0)));


}

//___________________________________________________________________
void AliHFEpidTPC::FillTPChistograms(const AliESDtrack *track, const AliMCParticle *mctrack, Bool_t stepSelected){
  // 
  // Fill the QA histogtrams
  //
  if(!track)
    return;
 
  Double_t tpcSignal = track->GetTPCsignal();
  Double_t p = track->GetInnerParam() ? track->GetInnerParam()->P() : track->P();
  Int_t species = -1;
  THnSparse *hptr = NULL;
  if(HasMCData() && mctrack){
    switch(TMath::Abs(mctrack->Particle()->GetPdgCode())){
      case 11:   
        species = AliPID::kElectron;
        if(!stepSelected){
          Double_t contentElHist[4];
          for(Int_t ispec = AliPID::kMuon; ispec < AliPID::kSPECIES; ispec++){
            contentElHist[0] = ispec;
            contentElHist[1] = p;
            contentElHist[2] = -Suppression(track, ispec);
            contentElHist[3] = Likelihood(track, ispec);
            hptr = dynamic_cast<THnSparseF *>(fQAList->Get("fHistTPCel"));
            hptr->Fill(contentElHist);
          }
        }
	      break;
      case 13:    species = AliPID::kMuon; break;
      case 211:   species = AliPID::kPion; break;
      case 321:   species = AliPID::kKaon; break;
      case 2212:  species = AliPID::kProton; break; 
      default:    species = -1; break;
    }
    if(!stepSelected){
      // Fill Probability Histogram
      Double_t contentProb[3] = {species , p, Likelihood(track, 0)};
      hptr = dynamic_cast<THnSparseF *>(fQAList->Get("fHistTPCprob"));
      hptr->Fill(contentProb);
      // Fill suppression Histogram
      if(species > 0 && species < AliPID::kSPECIES){
        Double_t contentSup[3] = {species, p, Suppression(track, species)};
        hptr = dynamic_cast<THnSparseF *>(fQAList->Get("fHistTPCsuppression"));
        hptr->Fill(contentSup);
      }
    }
  }
  
  // Fill signal histogram
  Double_t contentSignal[5] = {species, p, tpcSignal, fESDpid->NumberOfSigmasTPC(track, AliPID::kElectron), stepSelected ? 1 : 0};
  hptr = dynamic_cast<THnSparseF *>(fQAList->Get("fHistTPCsignal"));
  hptr->Fill(contentSignal);
}

//___________________________________________________________________
void AliHFEpidTPC::AddQAhistograms(TList *qaList){
  //
  // Create QA histograms for TPC PID
  //
  fQAList = new AliHFEcollection("fQAhistosTPC", "TPC QA histos");
  
  // First THnSparse we fill with the signal
  const Int_t kNdimSignal  = 5;
  Int_t nBins[kNdimSignal];
  Double_t binMin[kNdimSignal], binMax[kNdimSignal];
  nBins[0] = AliPID::kSPECIES + 1; binMin[0] = -1.; binMax[0] = AliPID::kSPECIES; // MC Species;
  nBins[1] = 1000; binMin[1] = 0.; binMax[1] = 20.; 
  nBins[2] = 6000; binMin[2] = 0.; binMax[2] = 600.;
  nBins[3] = 1400; binMin[3] = -12.; binMax[3] = 12.;
  nBins[4] = 2; binMin[4] = 0.; binMax[4] = nBins[4];  // Selected or not
  fQAList->CreateTHnSparse("fHistTPCsignal", "TPC signal; Species; p [GeV/c]; TPC Signal [a.u.]; Normalized TPC distance to the electron Line [#sigma]; Selection Status", kNdimSignal, nBins, binMin, binMax);

  const Int_t kNdimProbEl = 3;
  nBins[2] = 200; binMin[2] = 0.; binMax[2] = 1.;
  fQAList->CreateTHnSparse("fHistTPCprob", "TPC Likelihood to be an electron; Species; p [GeV/c]; TPC Likelihood [a.u.]", kNdimProbEl, nBins, binMin, binMax);

  const Int_t kNdimSuppression = 3;
  nBins[2] = 200; binMin[2] = -1.; binMax[2] = 5.8;  // log10 of TPC Likelihood(species)/Likelihood(elec) for species i neq electron
  fQAList->CreateTHnSparse("fHistTPCsuppression", "TPC non-electron Suppression; Species; p [GeV/c]; Suppression [a.u.]", kNdimSuppression, nBins, binMin, binMax);

  const Int_t kNdimEle = 4;
  nBins[0] = AliPID::kSPECIES - 1; binMin[0] = 1.; binMax[0] = AliPID::kSPECIES;
  nBins[2] = 100; binMin[2] = -1.; binMax[2] = 5.8;
  nBins[3] = 200; binMin[3] = 0.; binMax[3] = 1.; 
  fQAList->CreateTHnSparse("fHistTPCel", "TPC electron Histogram; Species; p [GeV/c]; Electron Enhancement:Electron Likelihood", kNdimEle, nBins, binMin, binMax);

  qaList->AddLast(fQAList->GetList());
}

