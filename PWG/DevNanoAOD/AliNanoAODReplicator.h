#ifndef ALINANOAODREPLICATO_H
#define ALINANOAODREPLICATO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id: AliAODMuonReplicator.h 56492 2012-05-15 18:42:47Z pcrochet $


#ifndef ALIDAODBRANCHREPLICATOR_H
#  include "AliAODBranchReplicator.h"
#endif
#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif

#include <iostream>

/* #ifndef AliAOD3LH_H */
/* #include "AliAOD3LH.h" */
/* #endif */

//#include "AliAOD3LH.h" 

//
// Implementation of a branch replicator 
// to produce nano AOD.
//
// Author: Michele Floris, michele.floris@cern.ch

class AliAnalysisCuts;
class TClonesArray;
class AliAODMCHeader;
class AliAODVZERO;
class AliAODTZERO;
class AliPIDResponse;
class AliESDv0;
class TArrayI;
class AliAODv0;  
class TRefArray;
class AliAODRecoDecay;
class AliAODRecoDecayHF;
class AliAODRecoDecayHF2Prong;
class AliVertexerTracks;

class AliESDVertex;
class AliESDtrack;
class AliVEvent;
class AliAODVertex;
class AliVertexerTracks;
class AliESDv0; 
class AliAODv0; 
class AliAODHeader;
class AliNanoAODHeader;
class AliAnalysisTaskSE;
class AliNanoAODTrack;
class AliAODTrack;
class AliNanoAODCustomSetter;

class TH1F;

class AliNanoAODReplicator : public AliAODBranchReplicator
{
 public:
  
  AliNanoAODReplicator(const char* name="AliNanoAODReplicator", 
			  const char* title="Branch Replicator for special tracks related branches",
			  const char * varlist="pt,eta",
			  AliAnalysisCuts* trackCut=0x0,
			  Int_t mcMode=0
			);
  
  virtual ~AliNanoAODReplicator();
  
  virtual TList* GetList() const ; // FIXME: This is declared const in the interface
  
  virtual void ReplicateAndFilter(const AliAODEvent& source);	

  virtual void Terminate();

  const char * GetVarList() { return fVarList; }
  void  SetVarList (const char * var) { fVarList = var;}
  const char * GetVarListHeader() { return fVarListHeader; }
  void  SetVarListHeader (const char * var) { fVarListHeader = var;}

  // Call backs for custom variables
  AliNanoAODCustomSetter * GetCustomSetter() { return fCustomSetter; }
  void  SetCustomSetter (AliNanoAODCustomSetter * var) { fCustomSetter = var;  }


 private:

  void SelectParticle(Int_t i);
  Bool_t IsParticleSelected(Int_t i);
  void CreateLabelMap(const AliAODEvent& source);
  Int_t GetNewLabel(Int_t i);
  void FilterMC(const AliAODEvent& source);
 

 private:
  
  AliAnalysisCuts* fTrackCut; // decides which tracks to keep
  mutable TClonesArray* fTracks; //! internal array of arrays of NanoAOD tracks
  mutable AliNanoAODHeader* fHeader; //! internal array of headers
  Int_t fNTracksVariables; //! Number of variables in the array
 
  mutable TClonesArray* fVertices; //! internal array of vertices
 
  mutable TList* fList; //! internal list of managed objects (fVertices and fTracks)
  
  mutable TClonesArray* fMCParticles; //! internal array of MC particles
  mutable AliAODMCHeader* fMCHeader; //! internal array of MC header
  Int_t fMCMode; // MC filtering switch (0=none=no mc information,1=normal=simple copy,>=2=aggressive=filter out : keep only particles leading to tracks and trheir relatives + all charged primaries)

  TExMap fLabelMap; //! for MC label remapping (in case of aggressive filtering)
  TExMap fParticleSelected; //! List of selected MC particles
			
  TString fVarList; // list of variables to be filterered
  TString fVarListHeader; // list of variables to be filtered (header)

  AliNanoAODCustomSetter * fCustomSetter;  // Setter class for custom variables

 private:

  
  AliNanoAODReplicator(const AliNanoAODReplicator&);
  AliNanoAODReplicator& operator=(const AliNanoAODReplicator&);
  
  ClassDef(AliNanoAODReplicator,1) // Branch replicator for ESD to muon AOD.
};

#endif
