#ifndef ALIJETAODFILLUNITARRAYTRACKS_H
#define ALIJETAODFILLUNITARRAYTRACKS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet Fill Unit Array 
// Called by ESD Reader for jet analysis
// Author: Magali Estienne (magali.estienne@subatech.in2p3.fr)
//---------------------------------------------------------------------

#include "AliJetFillUnitArray.h"

class AliJetReader;
class AliJetAODReader;
class AliEMCALGeometry;

class AliJetAODFillUnitArrayTracks : public AliJetFillUnitArray
{
 public: 
  AliJetAODFillUnitArrayTracks();
  AliJetAODFillUnitArrayTracks(AliAODEvent *fAOD);
  virtual ~AliJetAODFillUnitArrayTracks();
  
  // Setter
  void SetHadCorrector(AliJetHadronCorrection* const corr) {fHadCorr = corr;}
  void SetApplyMIPCorrection(Bool_t const val)             {fApplyMIPCorrection = val;}
  void SetAOD(AliAODEvent* const aod)                      {fAOD = aod;}
  void SetGrid0(AliJetGrid* const grid0)                   {fGrid0 = grid0;}
  void SetGrid1(AliJetGrid* const grid1)                   {fGrid1 = grid1;}
  void SetGrid2(AliJetGrid* const grid2)                   {fGrid2 = grid2;}
  void SetGrid3(AliJetGrid* const grid3)                   {fGrid3 = grid3;}
  void SetGrid4(AliJetGrid* const grid4)                   {fGrid4 = grid4;}

  // Getter
  Int_t GetHadCorrection()  const {return fApplyMIPCorrection;}
  Int_t GetMult()           const {return fNTracks;}
  Int_t GetMultCut()        const {return fNTracksCut;}

  // Other
  void Exec(Option_t* const option);

 protected:
  Int_t                     fNumUnits;           // Number of units in the unit object array (same as num towers in EMCAL)
  AliJetHadronCorrection   *fHadCorr;            // Pointer to Hadron Correction Object
  Bool_t                    fApplyMIPCorrection; // Apply MIP or not ? Exclusive with fApplyFractionHadronicCorrection

  AliAODEvent              *fAOD;                // AOD output Event
  AliJetGrid               *fGrid0;              // Grid used for dead zones definition
  AliJetGrid               *fGrid1;              // Grid used for dead zones definition
  AliJetGrid               *fGrid2;              // Grid used for dead zones definition
  AliJetGrid               *fGrid3;              // Grid used for dead zones definition
  AliJetGrid               *fGrid4;              // Grid used for dead zones definition

 private:
  AliJetAODFillUnitArrayTracks(const AliJetAODFillUnitArrayTracks &det);
  AliJetAODFillUnitArrayTracks &operator=(const AliJetAODFillUnitArrayTracks &det);
  void InitParameters();

  ClassDef(AliJetAODFillUnitArrayTracks,1) // Fill Unit Array with tpc and/or emcal information
};

#endif
