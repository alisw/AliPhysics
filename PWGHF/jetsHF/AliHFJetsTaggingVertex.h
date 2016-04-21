#ifndef ALIHFJETSTAGGINGVERTEX_H
#define ALIHFJETSTAGGINGVERTEX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* mailto: ycorrale@cern.ch */

/* Manager class for the HF b-jets analysis */

//--Root--
class TClonesArray;
class TObjArray;

//--AliRoot--
#include "AliPID.h"

class AliAODEvent;
class AliAODVertex;
class AliEmcalJet;
class AliESDtrack;
class AliESDVertex;

//--AliHFJetsClass--
#include "AliHFJetsUtils.h"
#include "AliHFJetsTagging.h"
#include "AliRDHFJetsCutsVertex.h"

//-------------------------------------------------------------------------------------

class AliHFJetsTaggingVertex : public AliHFJetsTagging {
  
public:
  
  // Constructors
  AliHFJetsTaggingVertex();
  
  AliHFJetsTaggingVertex(const char *name);
  
  // Destructor
  virtual ~AliHFJetsTaggingVertex();
  
  // Assignment operator
  AliHFJetsTaggingVertex &operator=(const AliHFJetsTaggingVertex &corr);
  
  Int_t FindVertices(const AliEmcalJet *jet,
                     TClonesArray      *fTrackArrayIn,
                     AliAODEvent       *aodEvent,
                     AliESDVertex      *primaryESDVertex,
                     Double_t           magZkG,
                     TClonesArray      *arrayVtxHF,
                     map_int_bool      *mapV0gTrks,
                     vctr_pair_dbl_int &vecVtxDisp);
  
  AliAODVertex *ReconstructSecondaryVertex(TObjArray     *trkArray,
                                           AliESDVertex  *v1,
                                           Double_t       magzkG,
                                           Double_t      &sigmaVtx) const;
  
  void     SetCuts(AliRDHFJetsCutsVertex *cuts);
  
  void     GetVtxPxy(AliAODVertex *vtx, Double_t *pxyzSum);
  
  Double_t GetVertexInvariantMass(AliAODVertex *vtx,
                                  Double_t massParticle = AliPID::ParticleMass(AliPID::kPion));
  
  Int_t    GetNTracksFromCommonVertex(AliAODVertex *vtx,
                                      const TClonesArray *mcPart,
                                      Int_t    &mcVert,
                                      Double_t &xVtxMC,
                                      Double_t &yVtxMC,
                                      Int_t    &nFromBandD,
                                      Int_t    &nFromD,
                                      Int_t    &nFromPromptD);

 private:
  
  AliRDHFJetsCutsVertex *fCutsHFjets;  // jet cut object
  
  TObjArray             *fTrackArray;  //! track array
  
  ClassDef(AliHFJetsTaggingVertex, 2);
};

//-------------------------------------------------------------------------------------
inline void AliHFJetsTaggingVertex::SetCuts(AliRDHFJetsCutsVertex *cuts) {

  if (fCutsHFjets)
    delete fCutsHFjets;
  
  fCutsHFjets = (AliRDHFJetsCutsVertex *)cuts->Clone("fCutsHFjets");
}

#endif
