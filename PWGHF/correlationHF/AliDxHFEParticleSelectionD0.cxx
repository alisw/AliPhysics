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
#include "AliReducedParticle.h"
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
  , fD0InvMass(0.0)
  , fPtBin(-1)
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

const char* AliDxHFEParticleSelectionD0::fgkDgTrackControlBinNames[]={
  "Pt",
  "Phi",
  "Ptbin", 
  "D0InvMass", 
  "Eta"
};

int AliDxHFEParticleSelectionD0::InitControlObjects()
{
  /// init the control objects, can be overloaded by childs which should
  /// call AliDxHFEParticleSelection::InitControlObjects() explicitly

  fD0Properties=DefineTHnSparse();
  AddControlObject(fD0Properties);

  //Adding control objects for the daughters
  InitControlObjectsDaughters("pi information",0);
  InitControlObjectsDaughters("K information",1);

  return AliDxHFEParticleSelection::InitControlObjects();
}

THnSparse* AliDxHFEParticleSelectionD0::DefineTHnSparse()
{
  //
  // Defines the THnSparse. 

  // here is the only place to change the dimension
  const int thnSize2 = 5;
  InitTHnSparseArray(thnSize2);
  
  const double Pi=TMath::Pi();
  TString name;
  name.Form("%s info", GetName());

  // 			             0     1     2       3         4
  // 	 	                     Pt   Phi   Ptbin  D0InvMass  Eta  
  int         thnBins [thnSize2] = {1000, 200,  15,     200,     500 };
  double      thnMin  [thnSize2] = {  0,    0,   0,    1.5648,   -1. };
  double      thnMax  [thnSize2] = { 100, 2*Pi, 14,    2.1648,    1. };
  const char* thnNames[thnSize2] = {"Pt", "Phi","Ptbin","D0InvMass","Eta"};

  return CreateControlTHnSparse(name,thnSize2,thnBins,thnMin,thnMax,thnNames);
}

int AliDxHFEParticleSelectionD0::FillParticleProperties(AliVParticle* p, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  AliAODRecoDecayHF2Prong* track=dynamic_cast<AliAODRecoDecayHF2Prong*>(p);
  if (!track) return -ENODATA;
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  data[i++]=track->Pt();
  data[i++]=track->Phi();
  data[i++]=fPtBin;
  data[i++]=fD0InvMass;
  data[i++]=track->Eta();

  return i;
}

int AliDxHFEParticleSelectionD0::InitControlObjectsDaughters(TString name, int daughter)
{
  // Setting up Control objects for the daughters.
  // Move to ParticleSelection?? 
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

  for(int iLabel=0; iLabel< 5;iLabel++)
    DaughterProperties->GetAxis(iLabel)->SetTitle(fgkDgTrackControlBinNames[iLabel]);  

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
 
  fD0InvMass= part->InvMassD0();
  fPtBin=fCuts->PtBin(part->Pt());
  
  // TODO: avoid repeated allocation of the arrays
  Double_t KProperties[]={prongneg->Pt(),prongneg->Phi(),fPtBin, fD0InvMass,prongneg->Eta()};
  Double_t piProperties[]={prongpos->Pt(),prongpos->Phi(),fPtBin,fD0InvMass,prongpos->Eta()};


  // Fills only for D0 or both.. 
  if ((selectionCode==1 || selectionCode==3) && fFillOnlyD0D0bar<2) {

    if(fD0Properties && ParticleProperties()) {
      memset(ParticleProperties(), 0, GetDimTHnSparse()*sizeof(ParticleProperties()[0]));
      FillParticleProperties(p, ParticleProperties(), GetDimTHnSparse());
      fD0Properties->Fill(ParticleProperties());
    }
    if(fD0Daughter0) fD0Daughter0->Fill(piProperties);
    if(fD0Daughter1) fD0Daughter1->Fill(KProperties);
  }
  else{
    // If not D0 (or both), check for D0bar (actually now also checks for both, not sure if this is needed)
    if ((selectionCode==2 || selectionCode==3) && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)) {

      if(fD0Properties && ParticleProperties()) {
	memset(ParticleProperties(), 0, GetDimTHnSparse()*sizeof(ParticleProperties()[0]));
	FillParticleProperties(p, ParticleProperties(), GetDimTHnSparse());
	fD0Properties->Fill(ParticleProperties());
      }
      if(fD0Daughter0) fD0Daughter0->Fill(piProperties);
      if(fD0Daughter1) fD0Daughter1->Fill(KProperties);
    }
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
  selectedTracks->SetOwner();
  TIter itrack(pTracks);
  TObject* pObj=NULL;
  while ((pObj=itrack())!=NULL) {
    AliVParticle* track=dynamic_cast<AliVParticle*>(pObj);
    if (!track) continue;
    int selectionCode=IsSelected(track,pEvent);
    HistogramParticleProperties(track, selectionCode);

    // This should make sure the array only gets filled with D0, D0bar or both:

    // Add track if it is either defined as D0(selectionCode==1) or both 
    // D0bar and a D0 (selectionCode==3)
    if ((selectionCode==1 || selectionCode==3) && fFillOnlyD0D0bar<2) 
      selectedTracks->Add(CreateParticle(track));
    else{
      // Add track if it is either defined as D0bar(selectionCode==2) or both 
      // D0bar and a D0 (selectionCode==3)
      if ((selectionCode==2 || selectionCode==3) && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)) 
	selectedTracks->Add(CreateParticle(track));
    }
  }
  return selectedTracks;
}

int AliDxHFEParticleSelectionD0::IsSelected(AliVParticle* p, const AliVEvent* pEvent)
{
  /// TODO: implement specific selection of D0 candidates
  /// Could also return values based on where where selection "failed
  /// Selected. Return 0 (none), 1(D0), 2(D0bar) or 3 (both)

  int selectionCode=0;

  AliAODRecoDecayHF2Prong *d0 = dynamic_cast<AliAODRecoDecayHF2Prong*>(p);
  if (!d0) return 0;
  if(d0->GetSelectionMap()) if(!d0->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
      AliDebug(1,"Skip D0 from Dstar");
      return 0; //skip the D0 from Dstar
    }

  // TODO: the cuts instance should be const but the function definition of
  // AliRDHFCuts::IsSelected does not allow this
  AliRDHFCuts* cuts=const_cast<AliRDHFCuts*>(fCuts);
  if (!cuts) {
    selectionCode=0;
  } 
  else if(cuts->IsInFiducialAcceptance(d0->Pt(),d0->Y(421)) ) {

    Int_t ptbin=cuts->PtBin(d0->Pt());
    if(ptbin==-1) {
      AliDebug(1,"Pt out of bounds");
      return 0;
    } //out of bounds

    // TODO: the aod pointer should also be const but the function definition of
    // AliRDHFCuts::IsSelected does not allow this
    AliAODEvent* aod=NULL;
    if (pEvent) aod=dynamic_cast<AliAODEvent*>(const_cast<AliVEvent*>(pEvent));
  
    // Selected. Return 0 (none), 1 (D0), 2 (D0bar) or 3 (both)
    // check daughters before calling as there is unchecked code in
    // AliAODRecoDecayHF::HasBadDaughters called
    TObject* o=NULL;
    if ((o=d0->GetDaughter(0))!=NULL && dynamic_cast<AliAODTrack*>(o)!=NULL &&
	(o=d0->GetDaughter(1))!=NULL && dynamic_cast<AliAODTrack*>(o)!=NULL) {
    selectionCode=cuts->IsSelected(d0,AliRDHFCuts::kAll,aod); 

    AliDebug(1,Form("Candidate is %d \n", selectionCode));
    } else {
      AliDebug(1,"at least one daughter not found!");
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

AliVParticle *AliDxHFEParticleSelectionD0::CreateParticle(AliVParticle* track)
{
  //
  // Creates object containing only the variables needed for correlation
  // 

  AliReducedParticle *part = new AliReducedParticle(track->Eta(), track->Phi(), track->Pt(),fD0InvMass,fPtBin);

  return part;

}
