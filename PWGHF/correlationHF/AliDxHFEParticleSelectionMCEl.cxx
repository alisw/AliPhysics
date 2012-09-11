// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/// @file   AliDxHFEParticleSelectionMCEl.cxx
/// @author Hege Erdal, Matthias Richter
/// @date   2012-07-19
/// @brief  MC El selection for D0-HFE correlation
///

#include "AliDxHFEParticleSelectionMCEl.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "THnSparse.h"
#include "AliAODMCParticle.h"
#include "TH1F.h"
#include "TAxis.h"
#include "AliAODTrack.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;
using std::vector;


/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionMCEl)

AliDxHFEParticleSelectionMCEl::AliDxHFEParticleSelectionMCEl(const char* opt)
  : AliDxHFEParticleSelectionEl(opt)
  , fMCTools()
  , fElectronProperties(NULL)
  , fWhichCut(NULL)
  , fOriginMother(0)
  , fResultMC(0)
{
  // constructor
  // 
  // 
  // 
  // 

  fMCTools.~AliDxHFEToolsMC();

  // TODO: argument scan, build tool options accordingly
  // e.g. set mc mode first/last, skip control histograms
  TString toolopt("pdg=11 mc-last");
  new (&fMCTools) AliDxHFEToolsMC(toolopt);
}

//at the moment not used, keep for now (to be used when plotting)
const char* AliDxHFEParticleSelectionMCEl::fgkPDGMotherBinLabels[]={
  "d",
  "u",
  "s",
  "c",
  "b",
  "gluon",
  "gamma",
  "#pi^{0}",
  "#eta",
  "proton",
  "others"
};

const char* AliDxHFEParticleSelectionMCEl::fgkTrackControlBinNames[]={
  "Pt",
  "Phi",
  "Eta", 
  "MC statistics electron"
  "Mother", 
  "PDG of selected electron", 
};

AliDxHFEParticleSelectionMCEl::~AliDxHFEParticleSelectionMCEl()
{
  // destructor
}

THnSparse* AliDxHFEParticleSelectionMCEl::DefineTHnSparse() const
{
  //
  // Defines the THnSparse. 

  const int thnSize = 6;
  const double Pi=TMath::Pi();
  TString name;
  name.Form("%s info", GetName());

  //     		       0    1      2     3       4         5
  // 	 	               Pt   Phi   Eta  'stat e'  mother "pdg e"
  int    thnBins[thnSize] = { 1000,  200, 500,    2,       14,      10 };
  double thnMin [thnSize] = {    0,    0, -1., -0.5,     -1.5,    -0.5 };
  double thnMax [thnSize] = {  100, 2*Pi,  1.,  1.5,     12.5,     9.5 };
 
  return CreateControlTHnSparse(name,thnSize,thnBins,thnMin,thnMax,fgkTrackControlBinNames);

}

int AliDxHFEParticleSelectionMCEl::DefineParticleProperties(AliVParticle* p, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  AliAODTrack *track=(AliAODTrack*)p;
  if (!track) return -ENODATA;
  int i=0;
  // TODO: this corresponds to the THnSparse dimensions which is available in the same class
  // use this consistently
  const int requiredDimension=6;
  if (dimension!=requiredDimension) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  data[i++]=track->Pt();
  data[i++]=track->Phi();
  data[i++]=track->Eta();
  data[i++]=fResultMC;     // stat electron (MC electron or not)
  data[i++]=fOriginMother; // at the moment not included background. Should expand
  data[i++]=1;             // PDG e - not sure if needed, maybe only needed as separate histo? 
  
  return i;
}

int AliDxHFEParticleSelectionMCEl::IsSelected(AliVParticle* p, const AliVEvent* pEvent)
{
  /// overloaded from AliDxHFEParticleSelection: check particle
  /// H: Have changed function. Now doing particle selection first, then run MC over 
  /// selected tracks. Could configure it to be configurable, but not sure if it
  /// is needed.  
  /// Result from normal track selection is returned, result from MC is stored in
  /// THnSparse. 

  int iResult=0;
  if (!p || !pEvent){
    return -EINVAL;
  }
  // step 1:
  // optional MC selection before the particle selection
  if (fMCTools.MCFirst() && (iResult=CheckMC(p, pEvent))==0) {
    // histograming?
    return iResult;
  }

  // step 2 or 1, depending on sequence:
  // normal particle selection
  iResult=AliDxHFEParticleSelectionEl::IsSelected(p, pEvent);
  if (fMCTools.MCFirst() || iResult==0) return iResult;

  // step 2, only executed if MC check is last
  // optional MC selection after the particle selection
  iResult=CheckMC(p, pEvent);
  // TODO: why do we need to store the result in a member?
  fResultMC=iResult;
 
  return iResult;
}

int AliDxHFEParticleSelectionMCEl::CheckMC(AliVParticle* p, const AliVEvent* pEvent)
{
  /// check if MC criteria are fulfilled
  if (!p || !pEvent){
    return -EINVAL;
  }
  fOriginMother=-1;
  int iResult=0;

  if (!fMCTools.IsInitialized() && (iResult=fMCTools.InitMCParticles(pEvent))<0) {
    // TODO: message? but has to be filtered in order to avoid message flood
    return 0; // no meaningful filtering on mc possible
  }
 
  if (fMCTools.RejectByPDG(p,false)) {
    // rejected by pdg
    // would also like to use this info? or not meaningful?
    // NEED TO CHANGE THIS!
    return 0;
  }

  int pdgMother=0;
  pdgMother=fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetFirstMother);

  // Particles considered HFE background. can be expanded
  // Should be created only once 
  // TODO: that needs to be configured once to avoid performance penalty
  vector<int> motherPDGs;
  motherPDGs.push_back(AliDxHFEToolsMC::kPDGpi0); 
  motherPDGs.push_back(AliDxHFEToolsMC::kPDGeta);
  motherPDGs.push_back(AliDxHFEToolsMC::kPDGgamma);
  motherPDGs.push_back(AliDxHFEToolsMC::kPDGJpsi);

  if(fMCTools.RejectByPDG(pdgMother,motherPDGs)){
    pdgMother=fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetOriginMother);
    fOriginMother=fMCTools.GetOriginMother();
  }
  else{
    //Could this be done in a more elegant way?
    switch(pdgMother){
    case(AliDxHFEToolsMC::kPDGpi0): fOriginMother=AliDxHFEToolsMC::kNrOrginMother; break;
    case(AliDxHFEToolsMC::kPDGeta): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+1; break;
    case(AliDxHFEToolsMC::kPDGgamma): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+2;break;
    case(AliDxHFEToolsMC::kPDGJpsi): fOriginMother=AliDxHFEToolsMC::kNrOrginMother+3;break;
    }
  }

  /*if (fMCTools.RejectByMotherPDG(p)) {
    // rejected by pdg of original mother
    // H: want pdg of origin process to be stored in THnSparse
    // Not sure if this is needed... Need to check, using AliDxHFEToolsMC, who
    // first mother are, and also what origin is. Use this info here. 
    return 0;
    }*/

  return 1;
}

void AliDxHFEParticleSelectionMCEl::Clear(const char* option)
{
  /// clear internal memory
  fMCTools.Clear(option);
}
