
//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

// $Id$

/// @file   AliDxHFECorrelationMC.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-11-02
/// @brief  Worker class for DxHFE correlation
///

#include "AliDxHFECorrelationMC.h"
#include "AliVParticle.h"
#include "AliLog.h"
//#include "AliAnalysisCuts.h"         // required dependency libANALYSISalice.so
//#include "AliFlowTrackSimple.h"      // required dependency libPWGflowBase.so
//#include "AliFlowCandidateTrack.h"   // required dependency libPWGflowTasks.so
//#include "AliCFContainer.h"          // required dependency libCORRFW.so
//#include "AliAODRecoDecayHF2Prong.h"   // libPWGHFvertexingHF
#include "TObjArray.h"
#include "AliDxHFECorrelation.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "AliReducedParticle.h"
#include "AliDxHFEParticleSelection.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

ClassImp(AliDxHFECorrelationMC)

AliDxHFECorrelationMC::AliDxHFECorrelationMC(const char* name)
  : AliDxHFECorrelation(name?name:"AliDxHFECorrelationMC")
  , fMCEventType(0)
{
  // default constructor
  // 
  //

}


AliDxHFECorrelationMC::~AliDxHFECorrelationMC()
{
  // destructor
  //
  //

}

THnSparse* AliDxHFECorrelationMC::DefineTHnSparse()
{
  //
  // Defines the THnSparse. 

  AliDebug(1, "Creating Corr THnSparse");
  // here is the only place to change the dimension
  static const int sizeEventdphi = 10;  
  InitTHnSparseArray(sizeEventdphi);
  const double pi=TMath::Pi();
  Double_t minPhi=GetMinPhi();
  Double_t maxPhi=GetMaxPhi();

  // TODO: Everything here needed for eventmixing? 
  // 			                        0        1     2      3      4     5       6      7    8      9
  // 			                      D0invmass  PtD0 PhiD0 PtbinD0 Pte  dphi    dEta OrigD0 origEl process
  int         binsEventdphi[sizeEventdphi] = {   200,   1000,  100,  21,   1000, 100,     100,   10,     14,  100 };
  double      minEventdphi [sizeEventdphi] = { 1.5648,     0,    0,   0,     0 , minPhi, -0.9, -1.5,   -1.5, -0.5 };
  double      maxEventdphi [sizeEventdphi] = { 2.1648,   100, 2*pi,  20,    100, maxPhi,  0.9,  8.5,   12.5, 99.5 };
  const char* nameEventdphi[sizeEventdphi] = {
    "D0InvMass",
    "PtD0",
    "PhiD0",
    "PtBinD0",
    "PtEl",
    "#Delta#Phi",
    "#Delta#eta", 
    "Origin D0", 
    "Origin Electron",
    "Original Process"
};

  TString name;
  name.Form("%s info", GetName());


  return CreateControlTHnSparse(name,sizeEventdphi,binsEventdphi,minEventdphi,maxEventdphi,nameEventdphi);


}


int AliDxHFECorrelationMC::FillParticleProperties(AliVParticle* tr, AliVParticle *as, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  AliReducedParticle *ptrigger=(AliReducedParticle*)tr;
  AliReducedParticle *assoc=(AliReducedParticle*)as;
  if (!ptrigger || !assoc) return -ENODATA;
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  if(AliDxHFECorrelation::GetTriggerParticleType()==kD){
    data[i++]=ptrigger->GetInvMass();
    data[i++]=ptrigger->Pt();
    data[i++]=ptrigger->Phi();
    data[i++]=ptrigger->GetPtBin(); 
    data[i++]=assoc->Pt();
  } 
  else{
    data[i++]=assoc->GetInvMass();
    data[i++]=assoc->Pt();
    data[i++]=assoc->Phi();
    data[i++]=assoc->GetPtBin(); 
    data[i++]=ptrigger->Pt();
  }
  data[i++]=AliDxHFECorrelation::GetDeltaPhi();
  data[i++]=AliDxHFECorrelation::GetDeltaEta();
  if(AliDxHFECorrelation::GetTriggerParticleType()==kD){
    data[i++]=ptrigger->GetOriginMother();
    data[i++]=assoc->GetOriginMother();
  }
  else {
    data[i++]=assoc->GetOriginMother();
    data[i++]=ptrigger->GetOriginMother();
  }
  data[i++]=fMCEventType;
  
  return i;
}

int AliDxHFECorrelationMC::Fill(const TObjArray* triggerCandidates, TObjArray* associatedTracks, const AliVEvent* pEvent)
{
  // TODO: Implement more on MC?? (Needed?)
  return AliDxHFECorrelation::Fill(triggerCandidates,associatedTracks,pEvent);
}

