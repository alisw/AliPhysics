// $Id$

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

/// @file   AliDxHFEParticleSelectionD0.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  D0 selection for D0-HFE correlation
///

#include "AliDxHFEParticleSelectionD0.h"
#include "AliVParticle.h"
//#include "AliAnalysisCuts.h"       // required dependency libANALYSISalice.so
//#include "AliFlowTrackSimple.h"    // required dependency libPWGflowBase.so
//#include "AliFlowCandidateTrack.h" // required dependency libPWGflowTasks.so
//#include "AliCFContainer.h"        // required dependency libCORRFW.so
#include "AliAODRecoDecayHF2Prong.h" // libPWGHFvertexingHF
#include "AliRDHFCutsD0toKpi.h"
#include "TObjArray.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TString.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelectionD0)

AliDxHFEParticleSelectionD0::AliDxHFEParticleSelectionD0(const char* opt)
  : AliDxHFEParticleSelection("D0", opt)
  , fD0Properties(NULL)
  , fD0Daughter0(NULL)
  , fD0Daughter1(NULL)
  , fCuts(NULL)
  , fFillOnlyD0D0bar(0)
{
  // constructor
  // 
  // 
  // 
  // 
  TString strOption(opt);
  // TODO: one might need a proper argument parsing including
  // chopping whole string into individual arguments
  if (strOption.Contains("FillD0D0bar")) fFillOnlyD0D0bar=0;
  else if (strOption.Contains("FillOnlyD0")) fFillOnlyD0D0bar=1;
  else if (strOption.Contains("FillOnlyD0bar")) fFillOnlyD0D0bar=2;
}

AliDxHFEParticleSelectionD0::~AliDxHFEParticleSelectionD0()
{
  // destructor
  if (fD0Properties) {
    delete fD0Properties;
    fD0Properties=NULL;
  }
  if (fD0Daughter0) {
    delete fD0Daughter0;
    fD0Daughter0=NULL;
  }
  if (fD0Daughter1) {
    delete fD0Daughter1;
    fD0Daughter1=NULL;
  }

  // Note: external object deleted elsewhere  
  fCuts=NULL;
}

int AliDxHFEParticleSelectionD0::InitControlObjects()
{
  /// init the control objects, can be overloaded by childs which should
  /// call AliDxHFEParticleSelection::InitControlObjects() explicitly
  AliInfo("Setting up THnSparse");
  TString name;
  const int thnSize = 4;
  const double pi=TMath::Pi();

  // TODO: theta?
  // 			         0       1     2      3
  // 			      mass      Pt   Phi    Ptbin
  int    thnBins[thnSize] = {   200,   1000,  100,   100  };
  double thnMin [thnSize] = {  1.5648,   0,    0,     0   };
  double thnMax [thnSize] = {  2.1648, 100,  (2*pi), 100  };

  name.Form("%s info", GetName());
  std::auto_ptr<THnSparseF> D0Properties(new THnSparseF(name, name, thnSize, thnBins, thnMin, thnMax));

  if (D0Properties.get()==NULL) {
    return -ENOMEM;
  }
  int axis=0;
  D0Properties->GetAxis(axis++)->SetTitle("D0 inv mass");
  D0Properties->GetAxis(axis++)->SetTitle("Pt");
  D0Properties->GetAxis(axis++)->SetTitle("Phi"); 
  D0Properties->GetAxis(axis++)->SetTitle("Ptbin"); 

  fD0Properties=D0Properties.release();
  AddControlObject(fD0Properties);

  //Adding control objects for the daughters
  InitControlObjectsDaughters("pi information",0);
  InitControlObjectsDaughters("K information",1);

  return AliDxHFEParticleSelection::InitControlObjects();
}

int AliDxHFEParticleSelectionD0::InitControlObjectsDaughters(TString name, int daughter)
{
  //Setting up Control objects for the daughters.
  AliInfo("Setting up daughter THnSparse");

  const int thnSize2 = 5;
  const double Pi=TMath::Pi();
  // 			       0    1      2      3          4
  // 	 	               Pt   Phi   Ptbin  D0InvMass  Eta
  int    thnBins[thnSize2] = { 1000,  200, 21,     200,     500};
  double thnMin [thnSize2] = {    0,    0,  0,    1.5648,   -1.};
  double thnMax [thnSize2] = {  100, 2*Pi, 20,    2.1648,    1.};

  std::auto_ptr<THnSparseF> DaughterProperties(new THnSparseF(name, name, thnSize2, thnBins, thnMin, thnMax));

  if (DaughterProperties.get()==NULL) {
    return -ENOMEM;
  }
  int axis=0;
  DaughterProperties->GetAxis(axis++)->SetTitle("Pt");
  DaughterProperties->GetAxis(axis++)->SetTitle("Phi");
  DaughterProperties->GetAxis(axis++)->SetTitle("Ptbin"); 
  DaughterProperties->GetAxis(axis++)->SetTitle("D0InvMass"); 
  DaughterProperties->GetAxis(axis++)->SetTitle("Eta"); 

  if(daughter==0){ 
    fD0Daughter0=DaughterProperties.release();
    AddControlObject(fD0Daughter0);
  }
  
  if(daughter==1){
    fD0Daughter1=DaughterProperties.release();
    AddControlObject(fD0Daughter1);
  }
  return 0;
}

int AliDxHFEParticleSelectionD0::HistogramParticleProperties(AliVParticle* p, int selectionCode)
{

  /// histogram particle properties
  if (!p) return -EINVAL;

  // fill the common histograms
  AliDxHFEParticleSelection::HistogramParticleProperties(p, selectionCode);

  // no daughters to fill if 0 (= no candidate)
  if (selectionCode==0) return 0;

  AliAODRecoDecayHF2Prong* part=dynamic_cast<AliAODRecoDecayHF2Prong*>(p);

  if(!part) return 0;
  // Convention: 1. daughter is postive track, 2. = negative
  AliAODTrack *prongpos=(AliAODTrack*)part->GetDaughter(0);
  AliAODTrack *prongneg=(AliAODTrack*)part->GetDaughter(1);

  if(!prongpos || !prongneg) {
    return 0;
  }
 
  // Only D0s are filled 
  // TODO: Also include D0bar
  if ((selectionCode==1 || selectionCode==3) && fFillOnlyD0D0bar<2) {
    Double_t invmassD0 = part->InvMassD0();
    AliDebug(3,Form("pt %f,  nr pt bins %d",part->Pt(),fCuts->GetNPtBins()));
    Int_t ptbin=fCuts->PtBin(part->Pt());
    Double_t D0Properties[] = {invmassD0,part->Pt(),part->Phi(),ptbin};
    Double_t KProperties[]={prongneg->Pt(),prongneg->Phi(),ptbin, invmassD0,prongneg->Eta()};
    Double_t piProperties[]={prongpos->Pt(),prongpos->Phi(),ptbin,invmassD0,prongpos->Eta()};

    if(fD0Properties) fD0Properties->Fill(D0Properties);
    if(fD0Daughter0) fD0Daughter0->Fill(piProperties);
    if(fD0Daughter1) fD0Daughter1->Fill(KProperties);
  }

  return 0;
}

TObjArray* AliDxHFEParticleSelectionD0::Select(TObjArray* pTracks, const AliVEvent *pEvent)
{
  /// create selection, array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  if (!pTracks) return NULL;
  TObjArray* selectedTracks=new TObjArray;
  if (!selectedTracks) return NULL;
  TIter itrack(pTracks);
  TObject* pObj=NULL;
  while ((pObj=itrack())!=NULL) {
    AliVParticle* track=dynamic_cast<AliVParticle*>(pObj);
    if (!track) continue;
    int selectionCode=IsSelected(track,pEvent);
    HistogramParticleProperties(track, selectionCode);
    //TODO: Also add selection for D0bar

    // Add track if it is either defined as D0(selectionCode==1) or both 
    // D0bar and a D0 (selectionCode==3)
    if (! ((selectionCode==1 || selectionCode==3) && fFillOnlyD0D0bar<2)) continue;
    selectedTracks->Add(track);
  }
  return selectedTracks;
}

int AliDxHFEParticleSelectionD0::IsSelected(AliVParticle* p, const AliVEvent* pEvent)
{
  /// TODO: implement specific selection of D0 candidates
  /// Could also return values based on where where selection "failed"
  int selectionCode=0;

  AliAODRecoDecayHF2Prong *d0 = dynamic_cast<AliAODRecoDecayHF2Prong*>(p);
  if(d0->GetSelectionMap()) if(!d0->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
      AliDebug(1,"Skip D0 from Dstar");
      return 0; //skip the D0 from Dstar
    }

  // TODO: the cuts instance should be const but the function definition of
  // AliRDHFCuts::IsSelected does not allow this
  AliRDHFCuts* cuts=const_cast<AliRDHFCuts*>(fCuts);
  if (!cuts) {
    selectionCode=1;
  } else if(cuts->IsInFiducialAcceptance(d0->Pt(),d0->Y(421)) ) {
    //    if(cuts->IsSelected(d0,AliRDHFCuts::kTracks,pEvent))fNentries->Fill(6);       
    Int_t ptbin=cuts->PtBin(d0->Pt());
    if(ptbin==-1) {
      AliDebug(1,"Pt out of bounds");
      return 0;
    } //out of bounds

    // TODO: the aod pointer should also be const but the function definition of
    // AliRDHFCuts::IsSelected does not allow this
    AliAODEvent* aod=NULL;
    if (pEvent) aod=dynamic_cast<AliAODEvent*>(const_cast<AliVEvent*>(pEvent));
  
    // Selected. Return 0 (none), 1(D0), 2(D0bar) or 3 (both)
    selectionCode=cuts->IsSelected(d0,AliRDHFCuts::kAll,aod); 

    AliDebug(1,Form("Candidate is %d \n", selectionCode));
    TObjArray daughters;
    daughters.AddAt((AliAODTrack*)d0->GetDaughter(0),0);
    daughters.AddAt((AliAODTrack*)d0->GetDaughter(1),1);

    //check daughters
    if(!daughters.UncheckedAt(0) || !daughters.UncheckedAt(1)) {
      AliDebug(1,"at least one daughter not found!");
      daughters.Clear();
      return 0;
    }
  }

  return selectionCode;
}

void AliDxHFEParticleSelectionD0::SetCuts(TObject* cuts, int /*level*/)
{
  /// set cuts objects
  fCuts=dynamic_cast<AliRDHFCuts*>(cuts);
  if (!fCuts && cuts) {
    AliError(Form("cuts object is not of required type AliRDHFCuts but %s", cuts->ClassName()));
  }
}
