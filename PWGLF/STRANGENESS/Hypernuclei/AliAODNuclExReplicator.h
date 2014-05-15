#ifndef ALIAODNUCLEXREPLICATOR_H
#define ALIAODNUCLEXREPLICATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

#ifndef ALIDAODBRANCHREPLICATOR_H
#  include "AliAODBranchReplicator.h"
#endif
#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif


//
// Implementation of a branch replicator 
// to produce aods with only few branches.
//
// Authors: S. Bufalino (stefania.bufalino@cern.ch)
//          R. Lea      (ramona.lea@cern.ch)
//         
// Based on AliAODMuonReplicator.h (L. Aphecetche (Subatech))

//class AliAnalysisCuts;
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
class AliAODRecoDecayLF;
class AliAODRecoDecayLF2Prong;
class AliVertexerTracks;
class AliAODHeader;

class AliESDVertex;
class AliESDtrack;
class AliVEvent;
class AliAODVertex;
class AliVertexerTracks;
class AliESDv0; 
class AliAODv0; 

class TH1F;

// TODO add new constructor for the 3 prong case (it will write an AliAODRecoDecayLF3Prong

class AliAODNuclExReplicator : public AliAODBranchReplicator
{
 public:
  
  AliAODNuclExReplicator(const char* name="AliAODNuclExReplicator", 
			 const char* title="Branch Replicator for muon related branches",
			 /* AliAnalysisCuts* trackCut=0x0, */
			 /* AliAnalysisCuts* vertexCut=0x0, */
			 Int_t mcMode=0,
			 Int_t nsigmaTrk1=3, Int_t partType1 = 2,
			 Int_t nsigmaTrk2=3, Int_t partType2 = 7
			 ); 
  //  Int_t partType; // 0 = e; 1 = mu; 2 = pi; 3 = K; 4= p; 5 = d; 6 =t ; 7 = 3He; 8=4He;

  virtual ~AliAODNuclExReplicator();
  
  virtual TList* GetList() const;
  
  virtual void ReplicateAndFilter(const AliAODEvent& source);	
  
  virtual void Terminate();

 private:

  // TO DO : Implemet MC
  // void FilterMC(const AliAODEvent& source);
  void SelectParticle(Int_t i);
  Bool_t IsParticleSelected(Int_t i);
  void CreateLabelMap(const AliAODEvent& source);
  Int_t GetNewLabel(Int_t i);
 

  Double_t fBzkG;     // z componenent of field in kG

  Double_t fCosMin       ;
  Double_t fDCAtracksMin ;
  Double_t fRmax         ;
  Double_t fRmin         ;
  Double_t fDNmin        ;
  Double_t fDPmin        ;

 private:
  
  mutable AliAODHeader* fHeader; //! internal header object
  
  //  AliAnalysisCuts* fVertexCut; // decides which vertices to keep
  mutable TClonesArray* fVertices; //! internal array of vertices
  
  mutable TClonesArray* fNuclei; //! internal array of nuclei tracks
  
  //  AliAnalysisCuts* fTrackCut; // decides which tracks to keep
  mutable TClonesArray* fSecondaryVerices; //! internal array of secondary vertices canditates
   
  mutable TClonesArray* fDaughterTracks; //! internal SV daughter tracks



  mutable TList* fList; //! internal list of managed objects (fVertices and fTracks)
  
  mutable TClonesArray* fMCParticles; //! internal array of MC particles
  mutable AliAODMCHeader* fMCHeader; //! internal array of MC header
  Int_t fMCMode; // MC filtering switch (0=none=no mc information,1=normal=simple copy,>=2=aggressive=filter out)

  TExMap fLabelMap; //! for MC label remapping (in case of aggressive filtering)
  TExMap fParticleSelected; //! List of selected MC particles

  Bool_t fReplicateHeader; //! whether or not the replicate the AOD Header

  Int_t fnSigmaTrk1;
  Int_t fnSigmaTrk2;
  Int_t fpartType1;
  Int_t fpartType2;

			

  Bool_t fSecVtxWithKF; // if kTRUE use KF vertexer, else AliVertexerTracks
  AliVertexerTracks* fVertexerTracks; // vertexer, to compute secondary vertices
  AliESDVertex *fV1;   // primary vertex
  
  Int_t fAODMapSize;  // size of fAODMap 
  Int_t *fAODMap;     //[fAODMapSize] map between index and ID for AOD tracks
  
  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray,Double_t &dispersion,Bool_t useTRefArray=kTRUE) const;
  
  AliAODRecoDecayLF2Prong* Make2Prong(TObjArray *twoTrackArray,const AliAODEvent &evento,
  				      AliAODVertex *secVert,Double_t dca);

  /* AliAODRecoDecayLF2Prong* Make2Prong(TObjArray *twoTrackArray,AliAODEvent *evento, */
  /* 				      AliAODVertex *secVert,Double_t dca); */
  

  void AddDaughterRefs(AliAODVertex *v, const AliAODEvent &event, const TObjArray *trkArray) const;
  /* void AddDaughterRefs(AliAODVertex *v, AliAODEvent *event, const TObjArray *trkArray) const; */
  void AddRefs(AliAODVertex *v,AliAODRecoDecayLF *rd, const AliAODEvent &event, const TObjArray *trkArray) const;
  /* void AddRefs(AliAODVertex *v,AliAODRecoDecayLF *rd, AliAODEvent *event, const TObjArray *trkArray) const; */

 private:

  //  AliPIDResponse  *fPIDResponse;                  //! PID response object
  
  AliAODNuclExReplicator(const AliAODNuclExReplicator&);
  AliAODNuclExReplicator& operator=(const AliAODNuclExReplicator&);
  
  ClassDef(AliAODNuclExReplicator,5) // Branch replicator for ESD to muon AOD.
    };

#endif
