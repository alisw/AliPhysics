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

//////////////////////////////////////////////////////////////////////
// AliCFParticleGenCut implementation
// This class is designed to handle 
// particle selection at generated level.
//
// added support for MC in AOD tree (2008-11-04)
// added a bool flag for the alternative (standard MC) vs (AOD MC).
//
// author : R. Vernet (renaud.vernet@cern.ch)
//////////////////////////////////////////////////////////////////////


#ifndef ALICFPARTICLEGENCUTS_H
#define ALICFPARTICLEGENCUTS_H

#include "AliCFCutBase.h"

class AliMCEvent;
class TObject;
class AliMCParticle;
class AliStack;
class TList;
class TH1F;
class TH2F;
class TBits;
class TArrayF;
class TDecayChannel;
class AliVParticle;
class AliVEvent;
class AliAODMCParticle;


class AliCFParticleGenCuts : public AliCFCutBase
{
 public :
  AliCFParticleGenCuts() ;
  AliCFParticleGenCuts           (const Char_t* name, const Char_t* title) ;
  AliCFParticleGenCuts           (const AliCFParticleGenCuts& c) ;
  AliCFParticleGenCuts& operator=(const AliCFParticleGenCuts& c) ;
  virtual ~AliCFParticleGenCuts() { };
  virtual Bool_t IsSelected(TObject* obj) ;
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  virtual void   SetEvtInfo(TObject* mcEvent) ;
  void    SetAODMC(Bool_t flag) {fIsAODMC=flag;}

  Bool_t IsPrimaryCharged(AliVParticle *mcPart);
  Bool_t IsPrimary(AliMCParticle    *mcPart) ;
  Bool_t IsPrimary(AliAODMCParticle *mcPart) ;
  //static checkers
  static Bool_t IsCharged(AliVParticle *mcPart);
  static Bool_t IsA(AliMCParticle    *mcPart, Int_t pdg, Bool_t abs=kFALSE);
  static Bool_t IsA(AliAODMCParticle *mcPart, Int_t pdg, Bool_t abs=kFALSE);

  void SetRequireIsCharged   () {fRequireIsCharged  =kTRUE; fRequireIsNeutral  =kFALSE;}
  void SetRequireIsNeutral   () {fRequireIsNeutral  =kTRUE; fRequireIsCharged  =kFALSE;}
  void SetRequireIsPrimary   () {fRequireIsPrimary  =kTRUE; fRequireIsSecondary=kFALSE;}
  void SetRequireIsSecondary () {fRequireIsSecondary=kTRUE; fRequireIsPrimary  =kFALSE;}
  void SetRequirePdgCode     (Int_t pdg)            {fRequirePdgCode=kTRUE; fPdgCode=pdg;}
  void SetProdVtxRangeX    (Double32_t xmin, Double32_t xmax) {fProdVtxXMin   =xmin; fProdVtxXMax   =xmax;}
  void SetProdVtxRangeY    (Double32_t ymin, Double32_t ymax) {fProdVtxYMin   =ymin; fProdVtxYMax   =ymax;}
  void SetProdVtxRangeZ    (Double32_t zmin, Double32_t zmax) {fProdVtxZMin   =zmin; fProdVtxZMax   =zmax;}
  void SetDecayVtxRangeX   (Double32_t xmin, Double32_t xmax) {fDecayVtxXMin  =xmin; fDecayVtxXMax  =xmax;}
  void SetDecayVtxRangeY   (Double32_t ymin, Double32_t ymax) {fDecayVtxYMin  =ymin; fDecayVtxYMax  =ymax;}
  void SetDecayVtxRangeZ   (Double32_t zmin, Double32_t zmax) {fDecayVtxZMin  =zmin; fDecayVtxZMax  =zmax;}
  void SetDecayLengthRange (Double32_t rmin, Double32_t rmax) {fDecayLengthMin=rmin; fDecayLengthMax=rmax;}
  void SetDecayRxyRange    (Double32_t rmin, Double32_t rmax) {fDecayRxyMin   =rmin; fDecayRxyMax   =rmax;}
  void SetDecayChannel     (TDecayChannel* dc) {fDecayChannel = dc ;}

  enum { 
    kCutCharge,       // ischarged cut
    kCutPrimSec,      // isprimary cut
    kCutPDGCode,      // PDG code  cut
    kCutProdVtxXMin,  // production vertex cut
    kCutProdVtxXMax,  // production vertex cut
    kCutProdVtxYMin,  // production vertex cut
    kCutProdVtxYMax,  // production vertex cut
    kCutProdVtxZMin,  // production vertex cut
    kCutProdVtxZMax,  // production vertex cut
    kCutDecVtxXMin,   // decay vertex cut
    kCutDecVtxXMax,   // decay vertex cut
    kCutDecVtxYMin,   // decay vertex cut
    kCutDecVtxYMax,   // decay vertex cut
    kCutDecVtxZMin,   // decay vertex cut
    kCutDecVtxZMax,   // decay vertex cut
    kCutDecLgthMin,   // decay length cut
    kCutDecLgthMax,   // decay length cut
    kCutDecRxyMin,    // transverse decay length cut
    kCutDecRxyMax,    // transverse decay length cut
    kCutDecayChannel, // decay channel reuired
    kNCuts,           // number of single selections
    kNStepQA=2        // number of QA steps (before/after the cuts)
  };

 private:
  Bool_t fIsAODMC ;       // flag for standard MC or MC from AOD tree
  AliVEvent* fMCInfo ;    // pointer to the MC event information
  Bool_t     fRequireIsCharged;   // require charged particle
  Bool_t     fRequireIsNeutral;   // require neutral particle
  Bool_t     fRequireIsPrimary;   // require primary particle
  Bool_t     fRequireIsSecondary; // require secondary particle
  Bool_t     fRequirePdgCode;     // require check of the PDG code
  Int_t      fPdgCode ;           // particle PDG code
  Double32_t fProdVtxXMin;        // min X of particle production vertex
  Double32_t fProdVtxYMin;        // min Y of particle production vertex
  Double32_t fProdVtxZMin;        // min Z of particle production vertex
  Double32_t fProdVtxXMax;        // max X of particle production vertex
  Double32_t fProdVtxYMax;        // max Y of particle production vertex
  Double32_t fProdVtxZMax;        // max Z of particle production vertex
  Double32_t fDecayVtxXMin;       // min X of particle decay vertex
  Double32_t fDecayVtxYMin;       // min Y of particle decay vertex
  Double32_t fDecayVtxZMin;       // min Z of particle decay vertex
  Double32_t fDecayVtxXMax;       // max X of particle decay vertex
  Double32_t fDecayVtxYMax;       // max Y of particle decay vertex
  Double32_t fDecayVtxZMax;       // max Z of particle decay vertex
  Double32_t fDecayLengthMin;     // min decay length (absolute)
  Double32_t fDecayLengthMax;     // max decay length (absolute)
  Double32_t fDecayRxyMin;        // min decay length in transverse plane wrt (0,0,0)
  Double32_t fDecayRxyMax;        // max decay length in transverse plane wrt (0,0,0)
  TDecayChannel* fDecayChannel;   // decay channel 

  //QA histos
  TH1F*    fhCutStatistics;        // Histogram: statistics of what cuts the tracks did not survive
  TH2F*    fhCutCorrelation;	   // Histogram: 2d statistics plot
  TH1F*    fhQA[kNCuts][kNStepQA]; // QA Histograms
  TArrayF* fCutValues;             // array of cut values
  TBits* fBitmap ;                 // stores single selection decisions

  void SelectionBitMap(AliMCParticle*    obj); // for MC got from Kinematics
  void SelectionBitMap(AliAODMCParticle* obj); // for MC got from AOD
  void FillHistograms(TObject* obj, Bool_t afterCuts);
  void AddQAHistograms(TList *qaList) ;
  void DefineHistograms();

  ClassDef(AliCFParticleGenCuts,2);
};

#endif
