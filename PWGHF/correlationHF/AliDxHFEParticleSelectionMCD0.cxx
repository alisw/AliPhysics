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

/// @file   AliDxHFEParticleSelectionMCD0.cxx
/// @author Hege Erdal, Matthias Richter
/// @date   2012-07-19
/// @brief  MC D0 selection for D0-HFE correlation
///

#include "AliDxHFEParticleSelectionMCD0.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliVParticle.h"
#include "AliReducedParticle.h"
#include "TH1F.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionMCD0)

AliDxHFEParticleSelectionMCD0::AliDxHFEParticleSelectionMCD0(const char* opt)
  : AliDxHFEParticleSelectionD0(opt)
  , fMCTools()
  , fPDGnotMCD0(NULL)
  , fResultMC(0)
  , fOriginMother(0)
  , fUseKine(kFALSE)
{
  // constructor
  // 
  // 
  // 
  // 

  // TODO: argument scan, pass only relevant arguments to tools
  fMCTools.~AliDxHFEToolsMC();
  TString toolopt("pdg=421 mc-last");
  new (&fMCTools) AliDxHFEToolsMC(toolopt);
}

AliDxHFEParticleSelectionMCD0::~AliDxHFEParticleSelectionMCD0()
{
  // destructor
}

THnSparse* AliDxHFEParticleSelectionMCD0::DefineTHnSparse()
{
  //
  // Defines the THnSparse.

  // here is the only place to change the dimension
  const int thnSize2 = 6;
  InitTHnSparseArray(thnSize2);
  const double Pi=TMath::Pi();
  TString name;
  name.Form("%s info", GetName());

  // 			             0     1      2       3        4     5
  // 	 	                     Pt   Phi   Ptbin  D0InvMass  Eta   mother 
  int         thnBins [thnSize2] = {1000, 200,   15,     200,     500,    10  };
  double      thnMin  [thnSize2] = {   0,  0,     0,    1.5648,   -1.,  -1.5  };
  double      thnMax  [thnSize2] = { 100, 2*Pi,  14,    2.1648,    1.,   8.5  };
  const char* thnNames[thnSize2] = {
    "Pt",
    "Phi",
    "Ptbin", 
    "D0InvMass", 
    "Eta",
    "Mother of D0"  // Bin -1 = not MC truth D0, rest OK
  };

  // Add Histo displaying pdg of D0 candidates not passing MatchToMC()
  fPDGnotMCD0= new TH1F("fPDGnotMCD0","PDG of track not MC truth D0",1002,-2.5,999.5);
  AddControlObject(fPDGnotMCD0);

  return CreateControlTHnSparse(name,thnSize2,thnBins,thnMin,thnMax,thnNames);
}

int AliDxHFEParticleSelectionMCD0::FillParticleProperties(AliVParticle* p, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  AliAODTrack *track=(AliAODTrack*)p;
  if (!track) return -ENODATA;
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  data[i++]=track->Pt();
  data[i++]=track->Phi();
  data[i++]=AliDxHFEParticleSelectionMCD0::GetPtBin(); 
  data[i++]=AliDxHFEParticleSelectionMCD0::GetInvMass();
  data[i++]=track->Eta();
  data[i++]=fOriginMother; // at the moment not included background. Should expand

  return i;
}

int AliDxHFEParticleSelectionMCD0::IsSelected(AliVParticle* p, const AliVEvent* pEvent)
{
  /// overloaded from AliDxHFEParticleSelection: check particle
  /// H: Have changed function. Now doing particle selection first, then run MC over 
  /// selected tracks. Could configure it to be configurable, but not sure if it
  /// is needed.  
  /// result from normal track selection is returned, result from MC is stored in
  /// THnSparse. 

  int iResult=0;
  fOriginMother=-1;

  // step 1:
  // MC selection
  if (fMCTools.MCFirst() && (iResult=CheckMC(p, pEvent))==0) {
    // histograming?
    return iResult;
  }

  // step 2 or 1, depending on sequence:
  // normal particle selection
  iResult=AliDxHFEParticleSelectionD0::IsSelected(p, pEvent);
  if (fMCTools.MCFirst() || iResult==0) return iResult;

  // step 2, only executed if MC check is last
  // MC selection  - > Should maybe also distinguish between D0 and D0bar
  // result stored to be filled into THnSparse
  // TODO: strictly speaken the particles should be rejected
  // if not mc selected, however skip this for the moment, because of
  // the logic outside
  fResultMC=CheckMC(p, pEvent);

  return iResult;
}

int AliDxHFEParticleSelectionMCD0::CheckMC(AliVParticle* p, const AliVEvent* pEvent)
{
  /// check if MC criteria are fulfilled
  // Check both D0 and D0bar (for now only D0)

  if (!p || !pEvent){
    return -EINVAL;
  }
  int iResult=0;

  if (!fMCTools.IsInitialized() && (iResult=fMCTools.InitMCParticles(pEvent))<0) {
    // TODO: message? but has to be filtered in order to avoid message flood
    return 0; // no meaningful filtering on mc possible
  }

  AliAODRecoDecayHF2Prong *particle = dynamic_cast<AliAODRecoDecayHF2Prong*>(p);

  if(!particle) return 0;

  Int_t pdgDgD0toKpi[2]={AliDxHFEToolsMC::kPDGkaon,AliDxHFEToolsMC::kPDGpion};

  TClonesArray* fMCArray = dynamic_cast<TClonesArray*>(fMCTools.GetMCArray());
  if(!fMCArray) {cout << "no array" << endl; return -1;}

  // find associated MC particle for D0->Kpi
  Int_t MClabel=-9999;

  //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx). Checks both D0s and daughters
  MClabel=particle->MatchToMC(AliDxHFEToolsMC::kPDGD0,fMCArray,2,pdgDgD0toKpi); 
  
  //TODO: Need a different strategy!!!
  // ALSO: look at AliAnalysisTaskSED0Mass for tips
  if(MClabel<0){
    // Checking PDG of particle if not MC truth D0
    // TODO: done the right way??
    Int_t MCl = p->GetLabel();
    if(MCl<0) {
      fPDGnotMCD0->Fill(-2);
      return 0;
    }
    int pdgPart=-1;
    AliAODMCParticle* aodmcp=0;
    aodmcp=dynamic_cast<AliAODMCParticle*>(fMCArray->At(MCl));
    if (aodmcp)
      pdgPart=TMath::Abs(aodmcp->GetPdgCode());
    if (pdgPart<0){
      fPDGnotMCD0->Fill(-1);
      return 0;
    }
    else{
      fPDGnotMCD0->Fill(pdgPart);
    }
    fOriginMother=-1;
    return 0;
  }

  fMCTools.SetMClabel(MClabel);
  fMCTools.FindMotherPDG(p,AliDxHFEToolsMC::kGetOriginMother);
  fOriginMother=fMCTools.GetOriginMother();

  return 1;
}

void AliDxHFEParticleSelectionMCD0::Clear(const char* option)
{
  /// clear internal memory
  fMCTools.Clear(option);
}

AliVParticle *AliDxHFEParticleSelectionMCD0::CreateParticle(AliVParticle* track)
{
  //
  //Created object which contain variables needed for correlation. 
  //

  AliReducedParticle *part = new AliReducedParticle(track->Eta(), track->Phi(), track->Pt(),AliDxHFEParticleSelectionMCD0::GetInvMass(),AliDxHFEParticleSelectionMCD0::GetPtBin(), fOriginMother);

  return part;

}
