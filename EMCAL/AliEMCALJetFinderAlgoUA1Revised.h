
#ifndef ALIEMCALJETFINDERALGOUA1REVISED_H
#define ALIEMCALJETFINDERALGOUA1REVISED_H

//THIS IS SARAH'S REVISED UA1 CODE WITH CHANGES FOR ETA/PHI ITERATION INCLUDED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  *  * See cxx source for full Copyright notice     */

/* $Id$ */

//*--Author: Sarah Blyth (LBL)
//*--Based on UA1 jet algorithm from LUND JETSET called from EMC-erj


#include "TTask.h"
#include "AliEMCALJetFinderInput.h"
#include "AliEMCALJetFinderOutput.h"
#include "AliEMCALJetFinderAlgo.h"
#include "AliEMCALJetFinderAlgoUA1Unit.h"
#include "AliEMCALJet.h"
#include "AliEMCALHadronCorrectionv1.h"


class AliEMCALJetFinderAlgoUA1Revised : public AliEMCALJetFinderAlgo
{

public:
  AliEMCALJetFinderAlgoUA1Revised();
  ~AliEMCALJetFinderAlgoUA1Revised();
  void InitUnitArray();                      
  void SetNumUnits(Int_t numUnits)           {fNumUnits = numUnits;}
  Int_t GetNumUnits() const                  {return fNumUnits;}
  void SetJetESeed(Float_t eSeed)            {fESeed = eSeed;}
  Float_t GetJetESeed() const                {return fESeed;}
  void SetConeRad(Float_t coneRad)           {fConeRad = coneRad;}
  Float_t GetConeRad() const                 {return fConeRad;}
  void SetJetEMin(Float_t jetEMin)           {fJetEMin = jetEMin;}
  Float_t GetJetEMin() const                 {return fJetEMin;}
  void SetEtMin(Float_t etMin)               {fEtMin = etMin;}
  Float_t GetEtMin() const                   {return fEtMin;}
  void SetMinMove(Float_t minMove)           {fMinMove = minMove;}
  void SetMaxMove(Float_t maxMove)           {fMaxMove = maxMove;}
  void SetBGMaxMove(Float_t bgMaxMove)       {fBGMaxMove = bgMaxMove;}
  void SetPtCut(Float_t ptCut)               {fPtCut = ptCut;}
  Float_t GetPtCut() const                   {return fPtCut;}
  void SetHadronCorrection(AliEMCALHadronCorrectionv1 *hadCorr) {fHadCorr = hadCorr;}
  AliEMCALHadronCorrectionv1* GetHadronCorrection() const {return fHadCorr;}
  void SetJetFindingParameters(Int_t numUnits, Float_t eSeed, Float_t coneRad, Float_t jetEMin, Float_t etMin, 
                               Float_t minMove, Float_t maxMove, Float_t bgMaxMove); 
  void SetJetFindingParameters(Int_t numUnits, Float_t eSeed, Float_t coneRad, Float_t jetEMin, Float_t etMin);
  void FillUnitArray(AliEMCALJetFinderAlgoUA1FillUnitFlagType_t flag);
  void Sort(AliEMCALJetFinderAlgoUA1Unit* unit,Int_t integer);
  void FindBG();
  void FindJetEtaPhi(Int_t counter);
  void FindJetEnergy();
  void StoreJetInfo();
  void FindJets();
  void QS(AliEMCALJetFinderAlgoUA1Unit *unit, Int_t left, Int_t right);

protected:
  AliEMCALJetFinderAlgoUA1Unit   *fUnit; //Array of JetFinder Unit objects (treated as the cells)
  AliEMCALHadronCorrectionv1 *fHadCorr; //Pointer to Hadron Correction Object
  AliEMCALJet       fJet;              //Jet object
  Int_t             fNumIter;          //Number of iterations for entire algorithm
  Int_t             fNumUnits;         //Number of units in the unit object array (same as num towers in EMCAL)
  Float_t           fESeed;            //Minimum energy a cell must have to be considered a jet seed
  Float_t           fConeRad;          //Size of cone radius 
  Float_t           fJetEMin;          //Minimum energy a cluster must have to be considered a jet
  Float_t           fEtMin;            //Minimum cell Et cut
  Float_t           fMinMove;          //Minimum move of jet centre from its previous position
  Float_t           fMaxMove;          //Maximum move allowed of jet centre from its initiator cell position
  Float_t           fBGMaxMove;        //Maximum allowed change in background energy between iterations 
  Float_t           fPtCut;            //Pt cut for tracks to minimise background contribution

  Float_t           fEBGTotal;         //Total background energy
  Float_t           fEBGTotalOld;      //Old total background energy
  Float_t           fEBGAve;           //Average background energy
  Float_t           fEnergy;           //Energy 
  Float_t           fJetEta;           //Jet eta value
  Float_t           fJetPhi;           //Jet phi value
  Float_t           fEtaInit;          //Jet initiate cell eta
  Float_t           fPhiInit;          //Jet initiate cell phi
  Float_t           fEtaB;             //Value of jet eta Before
  Float_t           fPhiB;             //Value of jet phi Before
  Float_t           fJetESum;          //Sum of weighted jet energy
  Float_t           fJetEtaSum;        //Sum of weighted jet eta
  Float_t           fJetPhiSum;        //Sum of weighted jet phi
  Float_t           fDEta;             //Offset of unit from jet eta value
  Float_t           fDPhi;             //Offset of unit from jet phi value
  Float_t           fDistP;            //Distance of new jet axis position from Previous position
  Float_t           fDistI;            //Distance of new jet axis position from Initiator cell position
  Float_t           fTempE;            //Temp E for comparing with JetEMin
  Float_t           fRad;              //Distance of cell from jet cone centre
  Int_t             fNumInCone;        //Number of units in the jet cone
  Int_t             fNumJets;          //Number of jets in an event
  Bool_t            fArrayInitialised; //To check that array of units is initialised

  ClassDef(AliEMCALJetFinderAlgoUA1Revised,2)

};
#endif

