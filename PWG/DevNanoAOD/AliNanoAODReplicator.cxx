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

// $Id: AliNanoAODReplicator.cxx 56492 2012-05-15 18:42:47Z pcrochet $

// Implementation of a branch replicator 
// to produce nanoAODs.
//
//
// This replicator is in charge of replicating the tracks,vertices,headers
// branches of a standard AOD or ESD file into a nanoAODs 
// (AliAOD.Special.root)
// 
// The class was inspired by AliAODMuonReplicator
// 
// Author: Michele Floris, michele.floris@cern.ch


class AliAODVertex;
class AliAODRecoDecay;

#include "AliAODDimuon.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODTZERO.h"
#include "AliAODTrack.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAnalysisCuts.h"
#include "TF1.h"
#include "AliExternalTrackParam.h"
#include "AliAODv0.h"
#include "AliPIDResponse.h"
#include <iostream>
#include <cassert>
#include "TObjArray.h"
#include "AliAnalysisFilter.h"
#include "AliNanoAODTrack.h"

#include <TFile.h>
#include <TDatabasePDG.h>
#include <TString.h>
#include <TList.h>
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
#include "AliVertexerTracks.h"
#include "AliKFVertex.h"
#include "AliAODEvent.h"
#include "AliAnalysisFilter.h"

#include "AliNanoAODReplicator.h"
#include "TH1.h"
#include "TCanvas.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODCustomSetter.h"
#include "AliV0ReaderV1.h"
#include "AliAnalysisNanoAODCuts.h"

using std::cout;
using std::endl;

ClassImp(AliNanoAODReplicator)

//_____________________________________________________________________________
AliNanoAODReplicator::AliNanoAODReplicator() :
AliAODBranchReplicator(), 
  fTrackCuts(0), 
  fV0Cuts(0), 
  fCascadeCuts(0), 
  fConversionPhotonCuts(0),
  fMCParticleCuts(nullptr),
  fTracks(0x0), 
  fHeader(0x0), 
  fVertices(0x0), 
  fList(0x0),
  fMCParticles(0x0),
  fMCHeader(0x0),
  fMCMode(0),
  fLabelMap(),
  fParticleSelected(),
  fVarList(""),
  fVarListHeader(""),
  fVarListHeader_fTC(""),
  fCustomSetters(),
  fVzero(0x0),
  fAodZDC(0x0),
  fV0s(0x0),
  fCascades(0x0),
  fConversionPhotons(0x0),
  fSaveZDC(0),
  fSaveVzero(0),
  fSaveV0s(0),
  fSaveCascades(kFALSE),
  fSaveConversionPhotons(kFALSE),
  fPhotonFromDeltas(kFALSE),
  fDeltaAODBranchName(""),
  fInputArrayName(""),
  fOutputArrayName("tracks"),
  fKeepDaughters(),
  fClonedVertices()
  {
  // Default ctor. we need it to avoid instantiating a wrong mapping when reading from file
  }

AliNanoAODReplicator::AliNanoAODReplicator(const char* name, const char* title) :
  AliAODBranchReplicator(name,title), 
  fTrackCuts(0), 
  fV0Cuts(0), 
  fCascadeCuts(0), 
  fConversionPhotonCuts(0),
  fMCParticleCuts(nullptr),
  fTracks(0x0), 
  fHeader(0x0), 
  fVertices(0x0), 
  fList(0x0),
  fMCParticles(0x0),
  fMCHeader(0x0),
  fMCMode(0),
  fLabelMap(),
  fParticleSelected(),
  fVarList(""),
  fVarListHeader(""),
  fVarListHeader_fTC(""),
  fCustomSetters(),
  fVzero(0x0),
  fAodZDC(0x0),
  fV0s(0x0),
  fCascades(0x0),
  fConversionPhotons(0x0),
  fSaveZDC(0),
  fSaveVzero(0),
  fSaveV0s(0),
  fSaveCascades(kFALSE),
  fSaveConversionPhotons(kFALSE),
  fPhotonFromDeltas(kFALSE),
  fDeltaAODBranchName(""),
  fInputArrayName(""),
  fOutputArrayName("tracks"),
  fKeepDaughters(),
  fClonedVertices()
{
  // default ctor
}

//_____________________________________________________________________________
AliNanoAODReplicator::~AliNanoAODReplicator()
{
  // dtor
  delete fTrackCuts;
  delete fList;
}

//_____________________________________________________________________________
void AliNanoAODReplicator::SelectParticle(Int_t i)
{
  // taking the absolute values here, need to take care 
  // of negative daughter and mother
  // IDs when setting!
  
  if (!IsParticleSelected(TMath::Abs(i)))
  {
    fParticleSelected.Add(TMath::Abs(i),1);    
  }
}

//_____________________________________________________________________________
Bool_t AliNanoAODReplicator::IsParticleSelected(Int_t i)  
{
  // taking the absolute values here, need to take 
  // care with negative daughter and mother
  // IDs when setting!
  return (fParticleSelected.GetValue(TMath::Abs(i))==1);
}


//_____________________________________________________________________________
void AliNanoAODReplicator::CreateLabelMap(const AliAODEvent& source)
{  
  //
  // this should be called once all selections are done 
  // This method associates to the original index of the mc particle
  // (i) the new one (j). J runs only over particles which are
  // actually kept.
  //
  
  fLabelMap.Delete();
  
  TClonesArray* mcParticles = static_cast<TClonesArray*>(source.FindListObject(AliAODMCParticle::StdBranchName()));
  
  Int_t j(0);
  Int_t i(0); // We need i, we cannot rely on part->GetLabel, because some of the original mc particles are not kept in the stack, apparently
  
  TIter next(mcParticles);

  while ( next() )
  {
    if (IsParticleSelected(i))
    {
      fLabelMap.Add(i,j++);
      //      std::cout << i <<  "->" << j-1 << std::endl;
    }
    ++i;
  }  


}

//_____________________________________________________________________________
Int_t AliNanoAODReplicator::GetNewLabel(Int_t i) 
{
  // Gets the label from the new created Map
  // Call CreatLabelMap before
  // otherwise only 0 returned
  return fLabelMap.GetValue(TMath::Abs(i));
}


//_____________________________________________________________________________
void AliNanoAODReplicator::RelabelAODPhotonCandidates(AliAODConversionPhoton *PhotonCandidate) {
  // taken from PWGGA/GammaConvBase/AliV0ReaderV1.cxx

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel
  Bool_t AODLabelPos = kFALSE;
  Bool_t AODLabelNeg = kFALSE;

  TIter nextTRACK(fTracks);
  AliNanoAODTrack* track;

  while ((track = static_cast<AliNanoAODTrack*>(nextTRACK()))) {
    if (!AODLabelPos) {
      if (track->GetID() == PhotonCandidate->GetTrackLabelPositive()) {
        PhotonCandidate->SetMCLabelPositive(TMath::Abs(track->GetLabel()));
        // PhotonCandidate->SetLabelPositive(i);
        AODLabelPos = kTRUE;
      }
    }
    if (!AODLabelNeg) {
      if (track->GetID() == PhotonCandidate->GetTrackLabelNegative()) {
        PhotonCandidate->SetMCLabelNegative(TMath::Abs(track->GetLabel()));
        // PhotonCandidate->SetLabelNegative(i);
        AODLabelNeg = kTRUE;
      }
    }
    if (AODLabelNeg && AODLabelPos) {
      return;
    }
  }
  if (!AODLabelPos || !AODLabelNeg) {
    if (!AODLabelNeg) {
      PhotonCandidate->SetMCLabelNegative(-999999);
      PhotonCandidate->SetLabelNegative(-999999);
    }
    if (!AODLabelPos) {
      PhotonCandidate->SetMCLabelPositive(-999999);
      PhotonCandidate->SetLabelPositive(-999999);
    }
  }
}

//_____________________________________________________________________________
void AliNanoAODReplicator::FilterMC(const AliAODEvent& source)
{
  // Filter MC information

  AliAODMCHeader* mcHeader(0x0);
  TClonesArray* mcParticles(0x0);
  
  fParticleSelected.Delete();

  //  std::cout << "MC Mode: " << fMCMode << ", Tracks " << fTracks->GetEntries() << std::endl;
  
  mcHeader = static_cast<AliAODMCHeader*>(source.FindListObject(AliAODMCHeader::StdBranchName()));
  if (mcHeader) 
    *fMCHeader = *mcHeader;
  
  mcParticles = static_cast<TClonesArray*>(source.FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcParticles)
    return;
  
  if (fMCMode == 1)
  {
    // simple copy of input MC particles to ouput MC particles
    TIter nextMC(mcParticles);
    AliAODMCParticle* p;
    Int_t nmcout(0);
  
    while ( ( p = static_cast<AliAODMCParticle*>(nextMC()) ) )
      new ((*fMCParticles)[nmcout++]) AliAODMCParticle(*p);
  }
  else if (fMCMode == 2)
  {
    // keep all primaries
    TIter nextPart(mcParticles);
    AliAODMCParticle * prim = 0;
    Int_t iprim = 0;  // We need iprim, we cannot rely on part->GetLabel, because some of the original mc particles are not kept in the stack, apparently
    // also select all charged primaries 
    while ((prim = (AliAODMCParticle*) nextPart())) {
      if (fMCParticleCuts && fMCParticleCuts->IsSelected(prim)) {
        SelectParticle(iprim);
      }
      iprim++;
    }

    // loop on (kept) tracks to find their ancestors
    TIter nextTRACK(fTracks);
    AliNanoAODTrack* track;
    while ((track = static_cast<AliNanoAODTrack*>(nextTRACK()))) {
      Int_t label = TMath::Abs(track->GetLabel());
      while (label >= 0) {
        SelectParticle(label);
        AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(mcParticles
            ->UncheckedAt(label));
        if (!mother) {
          AliError(Form("Got a null mother ! Check that ! (label %d", label));  // FIXME: I think this error is not needed
          label = -1;
        } else {
          label = mother->GetMother();  // do not only keep particles which created a track, but all their mothers
        }
      }
    }

    // loop on (kept) v0 to find their ancestors
    TIter nextV0(fV0s);
    AliAODv0* v0;
    // Get the PDG codes we want to match to from the cut object
    std::vector<int> pdgCodesV0;
    if (fMCParticleCuts) {
      pdgCodesV0 = static_cast<AliAnalysisNanoAODMCParticleCuts*>(fMCParticleCuts)->GetKeepV0s();
    }
    while ((v0 = static_cast<AliAODv0*>(nextV0()))) {
      // Get the daughter labels
      for (Int_t i = 0; i < v0->GetNDaughters(); ++i) {
        AliNanoAODTrack *trk = (AliNanoAODTrack*) v0->GetDaughter(i);
        SelectParticle(trk->GetLabel());
      }
      // loop over all PDG codes we want to match the V0 to
      for (auto it : pdgCodesV0) {
        int label = v0->MatchToMC(TMath::Abs(it), mcParticles);
        while (label >= 0) {
          SelectParticle(label);
          AliAODMCParticle* mother = static_cast<AliAODMCParticle*>(mcParticles
              ->UncheckedAt(label));
          if (!mother) {
            AliError(Form("Got a null mother ! Check that ! (label %d", label));  // FIXME: I think this error is not needed
            label = -1;
          } else {
            label = mother->GetMother();  // do not only keep particles which created a track, but all their mothers
          }
        }
      }
    }

    // loop on (kept) cascades to find their ancestors
    std::vector<int> pdgCodesV0Cascade;
    if (fMCParticleCuts) {
      pdgCodesV0Cascade =
          static_cast<AliAnalysisNanoAODMCParticleCuts*>(fMCParticleCuts)
              ->GetKeepV0sCascades();
    }
    TIter nextCascade(fCascades);
    AliAODcascade* casc;

    while ((casc = static_cast<AliAODcascade*>(nextCascade()))) {
      AliNanoAODTrack *nTrackXi = (AliNanoAODTrack*) (casc->GetDaughter(1));
      AliNanoAODTrack *pTrackXi = (AliNanoAODTrack*) (casc->GetDaughter(0));
      AliNanoAODTrack *bachTrackXi =
          (AliNanoAODTrack*) (casc->GetDecayVertexXi()->GetDaughter(0));
      SelectParticle(nTrackXi->GetLabel());
      SelectParticle(pTrackXi->GetLabel());
      SelectParticle(bachTrackXi->GetLabel());
      for (auto it : pdgCodesV0Cascade) {
        // You might think this is wrong, its not ....
        int MamaLauda = casc->MatchToMC(TMath::Abs(it), mcParticles);
        auto mamaBachPart = ((AliAODMCParticle*) mcParticles->At(
            bachTrackXi->GetLabel()));
        int MamaBach = (mamaBachPart) ? mamaBachPart->GetMother() : -1;
        // Is it a brother from another mother
        if (MamaLauda > 0 && MamaLauda == MamaBach) {
          int label = MamaLauda;
          while (label >= 0) {
            SelectParticle (label);
            AliAODMCParticle* mother =
                static_cast<AliAODMCParticle*>(mcParticles->UncheckedAt(label));
            if (!mother) {
              AliError(
                  Form("Got a null mother ! Check that ! (label %d", label));  // FIXME: I think this error is not needed
              label = -1;
            } else {
              label = mother->GetMother();  // do not only keep particles which created a track, but all their mothers
            }
          }
        }
      }
    }

    // loop on (kept) photons to find their ancestors
    TIter nextPhoton(fConversionPhotons);
    AliAODConversionPhoton* phot;

    while ((phot = static_cast<AliAODConversionPhoton*>(nextPhoton()))) {
      RelabelAODPhotonCandidates(phot);
      const int labelPos = phot->GetMCLabelPositive();
      const int labelNeg = phot->GetMCLabelNegative();
      AliAODMCParticle * mcPartPos = (AliAODMCParticle*) mcParticles
          ->UncheckedAt(labelPos);
      AliAODMCParticle * mcPartNeg = (AliAODMCParticle*) mcParticles
          ->UncheckedAt(labelNeg);
      SelectParticle(labelPos);
      SelectParticle(labelNeg);

      if (mcPartPos && mcPartNeg) {
        if (mcPartPos->GetMother() > -1
            && (mcPartNeg->GetMother() == mcPartPos->GetMother())) {
          int label = mcPartPos->GetMother();
          while (label >= 0) {
            SelectParticle(label);
            AliAODMCParticle* mother =
                static_cast<AliAODMCParticle*>(mcParticles->UncheckedAt(label));
            if (!mother) {
              AliError(
                  Form("Got a null mother ! Check that ! (label %d", label));  // FIXME: I think this error is not needed
              label = -1;
            } else {
              label = mother->GetMother();  // do not only keep particles which created a track, but all their mothers
            }
          }
        }
      }
    }

    // Actual filtering and label remapping (shamelessly taken for the implementation of AliAODHandler::StoreMCParticles)
    CreateLabelMap(source);
  
    TIter nextMC(mcParticles);
    AliAODMCParticle* p;
    Int_t nmc(0);  // We need nmc, we cannot rely on part->GetLabel, because some of the original mc particles are not kept in the stack, apparently
    Int_t nmcout(0);
  
    while ( ( p = static_cast<AliAODMCParticle*>(nextMC()) ) )
    {
      AliAODMCParticle c(*p);

      if ( IsParticleSelected(nmc) )
      {
        // 
        Int_t d0 =  p->GetDaughterLabel(0);
        Int_t d1 =  p->GetDaughterLabel(1);
        Int_t m =   p->GetMother();
      
        // other than for the track labels, negative values mean
        // no daughter/mother so preserve it
      
        if(d0<0 && d1<0)
        {
          // no first daughter -> no second daughter
          // nothing to be done
          // second condition not needed just for sanity check at the end
          c.SetDaughter(0,d0);
          c.SetDaughter(1,d1);
        } 
        else if(d1 < 0 && d0 >= 0) 
        {
          // Only one daughter
          // second condition not needed just for sanity check at the end
          if(IsParticleSelected(d0))
            {
              c.SetDaughter(0,GetNewLabel(d0));
            } else 
            {
              c.SetDaughter(0,-1);
            }
          c.SetDaughter(1,d1);
        }
        else if (d0 > 0 && d1 > 0 )
        {
          // we have two or more daughters loop on the stack to see if they are
          // selected
          Int_t d0tmp = -1;
          Int_t d1tmp = -1;
          for (int id = d0; id<=d1;++id)
            {
              if (IsParticleSelected(id))
            {
              if(d0tmp==-1)
                {
                  // first time
                  d0tmp = GetNewLabel(id);
                  d1tmp = d0tmp; // this is to have the same schema as on the stack i.e. with one daugther d0 and d1 are the same 
                }
              else d1tmp = GetNewLabel(id);
            }
            }
          c.SetDaughter(0,d0tmp);
          c.SetDaughter(1,d1tmp);
        } 
        else 
        {
          AliFatal(Form("Unxpected indices %d %d",d0,d1));
        }
        
        if ( m < 0 )
        {
          c.SetMother(m);
        } 
        else 
        {
          if (IsParticleSelected(m)) 
            {
              c.SetMother(GetNewLabel(m));              
            }
          // else // FIXME: re-enable this checj. Sometimes it gets here. Still to be understood why
          //   {
          //     //              AliError(Form("PROBLEM Mother not selected %d", m));              
          //   }
        }
        
        new ((*fMCParticles)[nmcout++]) AliAODMCParticle(c);
      }
      
      ++nmc;        
    } //closes loop over MC particles
    
    // now remap the tracks...
  
    TIter nextTrack(fTracks);
    AliNanoAODTrack* t;
    //      std::cout << "Remapping tracks" << std::endl;
  
    while ( ( t = dynamic_cast<AliNanoAODTrack*>(nextTrack()) ) )
    {
      t->SetLabel(GetNewLabel(t->GetLabel()));
    }
  }
  
  AliDebug(1,Form("input mc %d output mc %d",
                  mcParticles ? mcParticles->GetEntries() : 0,
                  fMCParticles ? fMCParticles->GetEntries() : 0));
}

// //_____________________________________________________________________________
TList* AliNanoAODReplicator::GetList() const
{
  // return (and build if not already done) our internal list of managed objects
  
  if (!fList)
    {
      // sanity checks
      if (fSaveConversionPhotons) {
        // check if id field is in fVarList
        AliNanoAODTrackMapping::GetInstance(fVarList);
        if (AliNanoAODTrackMapping::GetInstance()->GetVarIndex("ID") == -1)
          AliFatal("Conversion Photons requested but field 'id' missing in track variables");
      }
      
      fList = new TList;
      fList->SetOwner(kTRUE);

      fTracks = new TClonesArray("AliNanoAODTrack");
      fTracks->SetName(fOutputArrayName.Data());
      fList->Add(fTracks);

      Int_t numberOfHeaderParam = 0;
      Int_t numberOfHeaderParamInt = 0;
      for (Int_t i=0; i < fVarListHeader.Length(); i++){
          if (fVarListHeader.Data()[i] == ',') numberOfHeaderParam++;
      }
      
      if (fVarListHeader.Contains("T0Spread"))
        numberOfHeaderParam+=3;
      
      if (fVarListHeader.CompareTo("")==0) 
        numberOfHeaderParam = 2;
      else 
        numberOfHeaderParam = numberOfHeaderParam+1 ;

      numberOfHeaderParamInt = AliNanoAODHeader::GetIntParameters(fVarListHeader);
      numberOfHeaderParam -= numberOfHeaderParamInt;
      fHeader = new AliNanoAODHeader(numberOfHeaderParam, numberOfHeaderParamInt);

      fHeader->SetName("header");
      fList->Add(fHeader);    
        
      if(fSaveVzero){
          fVzero = new AliAODVZERO();
          fList->Add(fVzero);
      }
        
      if(fSaveZDC){
          fAodZDC = new AliAODZDC();
          fList->Add(fAodZDC);
      }

      if(fSaveV0s){
          fV0s = new TClonesArray("AliAODv0",2);
          fV0s->SetName("v0s");
          fList->Add(fV0s);
      }
      
      if(fSaveCascades){
          fCascades = new TClonesArray("AliAODcascade",2);
          fCascades->SetName("cascades");
          fList->Add(fCascades);
      }
      
      if(fSaveConversionPhotons){
          fConversionPhotons = new TClonesArray("AliAODConversionPhoton",2);
          fConversionPhotons->SetName("conversionphotons");
          fList->Add(fConversionPhotons);
      }

      fVertices = new TClonesArray("AliAODVertex",2);
      fVertices->SetName("vertices");    
      fList->Add(fVertices);
    
      if ( fMCMode > 0 )
      {
        fMCHeader = new AliAODMCHeader;    
        fMCParticles = new TClonesArray("AliAODMCParticle",1000);
        fMCParticles->SetName(AliAODMCParticle::StdBranchName());
        fList->Add(fMCHeader);
        fList->Add(fMCParticles);
      }
    }
  return fList;
}

AliAODVertex* AliNanoAODReplicator::CloneAndStoreVertex(AliAODVertex* toClone)
{
  // Clone vertex if not yet cloned. Update list of to store daughter objects.
  
  if (fClonedVertices.find(toClone) != fClonedVertices.end())
    return fClonedVertices[toClone];
  
  AliAODVertex* copiedVertex = new((*fVertices)[fVertices->GetEntriesFast()]) AliAODVertex(*toClone);
  copiedVertex->SetUniqueID(0); // avoid reusing old unique ID which confuses TRef
  for (int nD = 0; nD<copiedVertex->GetNDaughters(); nD++)
    fKeepDaughters[copiedVertex].push_back(copiedVertex->GetDaughter(nD));
  copiedVertex->RemoveDaughters();
  fClonedVertices[toClone] = copiedVertex;
  return copiedVertex;
}

//_____________________________________________________________________________
void AliNanoAODReplicator::ReplicateAndFilter(const AliAODEvent& source)
{
  // Replicate (and filter if filters are there) the relevant parts we're interested in AODEvent
  
  fTracks->Clear("C");
  
  assert(fVertices!=0x0);
  fVertices->Clear("C");
  
  if (fV0s)
    fV0s->Clear("C");
    
  if (fCascades)
    fCascades->Clear("C");
  
  if (fConversionPhotons)
    fConversionPhotons->Clear("C");
  
  if (fMCMode > 0){
    if(!fMCHeader) {
      AliFatal(Form("fMCMode = %d, but MC header not found", fMCMode));
    }
    fMCHeader->Reset();
    if(!fMCParticles){
      AliFatal(Form("fMCMode = %d, but MC particles not found", fMCMode));
    }
    fMCParticles->Clear("C");
  }
  
  fKeepDaughters.clear();
  fClonedVertices.clear();

  fHeader->SetMapFiredTriggerClasses(fVarListHeader_fTC);

  // Set custom variables in the header if the callback is set
  for (std::list<AliNanoAODCustomSetter*>::iterator it = fCustomSetters.begin(); it != fCustomSetters.end(); ++it)
    (*it)->SetNanoAODHeader(&source, fHeader, fVarListHeader);

  // keep here only *primary* vertices
  TIter nextV(source.GetVertices());
  AliAODVertex* v;
  Int_t nvertices(0);
  while ( ( v = static_cast<AliAODVertex*>(nextV()) ) )
  {
    if (v->GetType()!=AliAODVertex::kPrimary && v->GetType()!=AliAODVertex::kMainSPD && 
      v->GetType()!=AliAODVertex::kPileupSPD && v->GetType()!=AliAODVertex::kMainTPC && v->GetType()!=AliAODVertex::kPrimaryTPC)
      continue;

    AliAODVertex* tmp = v->CloneWithoutRefs();
    AliAODVertex* copiedVertex = new((*fVertices)[nvertices++]) AliAODVertex(*tmp);
    copiedVertex->SetNContributors(v->GetNContributors()); 
    copiedVertex->SetUniqueID(0); // avoid reusing old unique ID which confuses TRef
    delete tmp;
    fClonedVertices[v] = copiedVertex;
  }
  
  if(fSaveVzero==1){
      fVzero->Clear("C");
      AliAODVZERO *vzeroAOD(0x0);
      vzeroAOD = static_cast<AliAODVZERO*>(source.GetVZEROData());
      *fVzero = *vzeroAOD;
  }
  
  if(fSaveZDC){
      fAodZDC->Clear("C");
      AliAODZDC *aodZDC(0x0);
      aodZDC = static_cast<AliAODZDC*>(source.GetZDCData());
      *fAodZDC = *aodZDC;
  }
  
  if (fSaveCascades) {
    TIter nextC(const_cast<AliAODEvent&>(source).GetCascades());
    AliAODcascade* cascade;
    Int_t n = 0;
    while ( ( cascade = static_cast<AliAODcascade*>(nextC()) ) )
    {
      if (fCascadeCuts && !fCascadeCuts->IsSelected(cascade))
        continue;
  
      // bachelor track and xi vertex
      AliAODVertex* copiedXi = CloneAndStoreVertex(cascade->GetDecayVertexXi());
      // store additional daughter AODVertex if needed
      for (int nD = 0; nD<cascade->GetDecayVertexXi()->GetNDaughters(); nD++) {
        auto vertex = dynamic_cast<AliAODVertex*> (cascade->GetDecayVertexXi()->GetDaughter(nD));
        if (vertex != nullptr)
          CloneAndStoreVertex(vertex);
      }
      
      // v0 vertex and tracks
      AliAODVertex* copiedV0Vertex = CloneAndStoreVertex(cascade->GetSecondaryVtx());
      
      // NOTE we don't have AliAODcascade::SetDecayVertexXi so have to use copy constructor here
      //AliAODcascade* nanoCascade = new((*fCascades)[n++]) AliAODcascade(*cascade); 
      //nanoCascade->SetDecayVertexXi(copiedXi);
      const Double_t momBach[] = { cascade->MomBachX(), cascade->MomBachY(), cascade->MomBachZ() };
      AliAODcascade* nanoCascade = new((*fCascades)[n++]) AliAODcascade(copiedXi, cascade->ChargeXi(), cascade->DcaXiDaughters(), cascade->DcaXiToPrimVertex(), 
                                                                        cascade->DcaBachToPrimVertex(), (const Double_t*) momBach, *cascade);
      nanoCascade->SetSecondaryVtx(copiedV0Vertex);
    }
  }  
  
  if(fSaveV0s){
    TIter nextV(source.GetV0s());
    AliAODv0* v;
    Int_t nV0s = 0;
    while ( ( v = static_cast<AliAODv0*>(nextV()) ) )
    {
      if (fV0Cuts && !fV0Cuts->IsSelected(v))
        continue;

      AliAODv0* nanoV0 = new((*fV0s)[nV0s++]) AliAODv0(*v);
      AliAODVertex* copiedVertex = CloneAndStoreVertex(v->GetSecondaryVtx());
      nanoV0->SetSecondaryVtx(copiedVertex);
    }
  }
  
  // Fix parent association. Not clear if this is needed downstream, though.
  for (auto it = fClonedVertices.begin(); it != fClonedVertices.end(); it++) {
    if (it->first->GetParent() == nullptr)
      continue;
    auto parent = fClonedVertices.find((AliAODVertex*) (it->first->GetParent()));
    if (parent != fClonedVertices.end())
      it->second->SetParent(parent->second);
    else
      it->second->SetParent(0x0);
  }

  // Tracks
  Int_t entries = -1;
  TClonesArray* particleArray = 0x0;

  if(!fInputArrayName.IsNull()){
    particleArray = static_cast<TClonesArray*> (source.FindListObject(fInputArrayName.Data()));
    entries = particleArray->GetEntries();
  }else{
    entries = source.GetNumberOfTracks();
  }

  // Photons
  std::list<Int_t> trackIDs; // tracks which should be kept as they are referred to
  if (fSaveConversionPhotons) {
    TClonesArray* gammaArray = nullptr;
    Int_t nConvPhotons = 0;
    static AliV0ReaderV1* photonReader = (AliV0ReaderV1*) AliAnalysisManager::GetAnalysisManager()->GetTask("ConvGammaAODProduction");
    if (photonReader) {
      if (photonReader->AreAODsRelabeled() != kFALSE)
        AliFatal("This requires not relabled tracks");
      gammaArray = photonReader->GetReconstructedGammas();
    } else if (fPhotonFromDeltas) {
      gammaArray = dynamic_cast<TClonesArray*>(source.FindListObject(fDeltaAODBranchName.Data()));
    } else {
      AliFatal("No Photon reader and no branch name specified");
    }
    for (int iGamma = 0; iGamma < gammaArray->GetEntriesFast(); ++iGamma) {
      auto *photonCandidate = dynamic_cast<AliAODConversionPhoton *>(gammaArray->At(iGamma));
      if (fConversionPhotonCuts && !fConversionPhotonCuts->IsSelected(photonCandidate))
        continue;
      trackIDs.push_back(photonCandidate->GetTrackLabelPositive());
      trackIDs.push_back(photonCandidate->GetTrackLabelNegative());
      auto copiedPhoton = new((*fConversionPhotons)[nConvPhotons++]) AliAODConversionPhoton(*photonCandidate);
      copiedPhoton->SetV0Index(-1); // related V0 is not stored
    }
  }
  
  std::map<TObject*, AliNanoAODTrack*> trackAssociation;
  
  // Tracks
  Int_t ntracks(0);
  for(Int_t j=0; j<entries; j++) {
    AliVTrack *track = 0x0;
    if (particleArray) track = (AliVTrack*)particleArray->At(j);
    else track = (AliVTrack*)source.GetTrack(j);

    AliAODTrack *aodtrack = (AliAODTrack*) track;

    Bool_t selected = kFALSE;
    if (!fTrackCuts || fTrackCuts->IsSelected(aodtrack)) 
      selected = kTRUE;
    
    // store tracks needed for V0s
    for (std::map<AliAODVertex*, std::vector<TObject*> >::iterator it = fKeepDaughters.begin(); it != fKeepDaughters.end(); it++) {
      if (std::find(it->second.begin(), it->second.end(), aodtrack) != it->second.end())
        selected = kTRUE;
    }
    
    // store tracks needed for conversions
    if (std::find(trackIDs.begin(), trackIDs.end(), aodtrack->GetID()) != trackIDs.end())
      selected = kTRUE;
    
    if (!selected)
      continue;

    AliNanoAODTrack* nanoTrack = new((*fTracks)[ntracks++]) AliNanoAODTrack (aodtrack, fVarList);

    for (std::list<AliNanoAODCustomSetter*>::iterator it = fCustomSetters.begin(); it != fCustomSetters.end(); ++it)
      (*it)->SetNanoAODTrack(aodtrack, nanoTrack);
    
    trackAssociation[aodtrack] = nanoTrack;
  }
  
  // Replace references to stored tracks. 
  // NOTE this has to respect the order in which they were stored (e.g. for a V0 the first daugther needs to be the positive one).
  for (std::map<AliAODVertex*, std::vector<TObject*> >::iterator it = fKeepDaughters.begin(); it != fKeepDaughters.end(); it++) {
    for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      //Printf("Vertex %p Track %p", it->first, *it2);
      auto track = dynamic_cast<AliAODTrack*> (*it2);
      auto vertex = dynamic_cast<AliAODVertex*> (*it2);
      if (track != nullptr && trackAssociation.find(*it2) != trackAssociation.end())
        it->first->AddDaughter(trackAssociation[*it2]);
      else if (vertex != nullptr && fClonedVertices.find(vertex) != fClonedVertices.end())
        it->first->AddDaughter(fClonedVertices[vertex]);
      else {
        Printf("Dumping useful information before abort.");
        it->first->Dump();
        (*it2)->Dump();
        AliFatal("You requested to store an AliAODVertex for which the daughter object is not available. Aborting.");
      }
    }
  }
  
  AliDebug(1,Form("tracks=%d vertices=%d", fTracks->GetEntries(),fVertices->GetEntries())); 
  
  // Finally, deal with MC information, if needed
  if ( fMCMode > 0 ) {
    FilterMC(source);      
  }
}

void AliNanoAODReplicator::Terminate()
{
}
