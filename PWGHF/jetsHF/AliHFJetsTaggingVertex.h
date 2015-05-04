#ifndef ALIHFJETSTAGGINGVERTEX_H
#define ALIHFJETSTAGGINGVERTEX_H

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

// ******************************************
// Manager class for HF jet analysis   
// Author: andrea.rossi@cern.ch, elena.bruna@cern.ch,min.jung.kweon@cern.ch,linus.feldkamp@cern.ch, s.lapointe@cern.ch
// *******************************************

#include "TClonesArray.h"
#include "AliAODJet.h"
#include "AliEmcalJet.h"
#include "AliAODEvent.h"
#include "AliRDHFJetsCutsVertex.h"
#include "AliHFJetsTagging.h"

class AliAODVertex;
class AliESDVertex;
class AliVertexerTracks;
class AliESDtrack;
class AliVEvent;
class AliHFJetstagging;

class AliHFJetsTaggingVertex : public AliHFJetsTagging {
 public: 

  AliHFJetsTaggingVertex();
  AliHFJetsTaggingVertex(const char* name);
  ~AliHFJetsTaggingVertex();
  // Assignment operator
  AliHFJetsTaggingVertex& operator=(const AliHFJetsTaggingVertex& corr);

  Int_t FindVertices(const AliEmcalJet *jet, AliAODTrack **fAODTrackInfoP, AliAODTrack **fAODTrackInfoN,TClonesArray *fTrackArrayIn, AliAODEvent* aod, AliESDVertex* v1, Double_t magzkG ,TClonesArray *arrVertices, Double_t *arrDispersion);
  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray, AliESDVertex* v1, Double_t magzkG ,Double_t &dispersion) const;

  void SetCuts(AliRDHFJetsCutsVertex *cuts){delete fCutsHFjets;fCutsHFjets=(AliRDHFJetsCutsVertex*)cuts->Clone("fCutsHFjets");}
  
  Double_t GetVertexInvariantMass(AliAODVertex *vtx,Double_t massParticle=0.138);
  void GetVtxPxy(AliAODVertex *vtx,Double_t *pxyzSum);
  Int_t GetNTracksFromCommonVertex(AliAODVertex *vtx,const TClonesArray *mcPart,Int_t &mcVert,Double_t &xVtxMC,Double_t &yVtxMC,Int_t &nfromBandD,Int_t &nfromD,Int_t &nfromPromptD);

 private:
  AliRDHFJetsCutsVertex *fCutsHFjets;  // jet cut object 
  TObjArray *fTrackArray; //! track array 
   
  
 private:
   
  ClassDef(AliHFJetsTaggingVertex,1);

};
#endif
