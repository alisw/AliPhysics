#ifndef ALIEMCALJETFINDER_H
#define ALIEMCALJETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */
//*-- Author:
//*-- Andreas Morsch (CERN)

#include <TTask.h>
#include "AliEMCALJet.h"

class TClonesArray;
class AliEMCALHadronCorrection;

class TH2F;

class AliEMCALJetFinder : public TTask {
 public:
    AliEMCALJetFinder();
    AliEMCALJetFinder(const char* name, const char *title);
    virtual ~AliEMCALJetFinder();
    virtual void Find(Int_t ncell, Int_t ncell_tot, Float_t etc[30000], 
		      Float_t etac[30000], Float_t phic[30000],
		      Float_t min_move, Float_t max_move, Int_t mode,
		      Float_t prec_bg, Int_t ierror);
    virtual void Find();
    virtual void FindTracksInJetCone();
    virtual void Test();
    virtual void  BuildTrackFlagTable();
    virtual Int_t SetTrackFlag(Float_t radius, Int_t pdgCode, Double_t charge);
    // Geometry
    virtual void SetCellSize(Float_t eta, Float_t phi);
    // Parameters
    virtual void SetDebug(Int_t flag = 0) {fDebug = flag;}
    virtual void SetConeRadius(Float_t par);
    virtual void SetEtSeed(Float_t par);
    virtual void SetMinJetEt(Float_t par);
    virtual void SetMinCellEt(Float_t par);
    virtual void SetPtCut(Float_t par = 1.);
    virtual void SetMomentumSmearing(Bool_t flag = kFALSE) {fSmear = flag;}
    virtual void SetEfficiencySim(Bool_t flag = kFALSE)    {fEffic = flag;}
    virtual void SetSamplingFraction(Float_t par = 12.9) {fSamplingF = par;}
    virtual void SetIncludeK0andN(Bool_t flag = kFALSE) {fK0N = flag;}
    // Correction of hadronic energy
    virtual void SetHadronCorrector(AliEMCALHadronCorrection* corr)
	{fHadronCorrector = corr;}
    virtual void SetHadronCorrection(Int_t flag = 1) {fHCorrection = flag;}
    // Results
    virtual Int_t   Njets();
    virtual Float_t JetEnergy(Int_t);
    virtual Float_t JetPhiL(Int_t);
    virtual Float_t JetPhiW(Int_t);
    virtual Float_t JetEtaL(Int_t);  
    virtual Float_t JetEtaW(Int_t);
    virtual TH2F* GetLego() {return fLego;}
    // I/O
    virtual void FillFromHits(Int_t flag = 0);
    virtual void FillFromHitFlaggedTracks(Int_t flag = 0);
    virtual void FillFromDigits(Int_t flag = 0);
    virtual void FillFromTracks(Int_t flag = 0, Int_t ich = 0);
    virtual void AddJet(const AliEMCALJet& jet);
    virtual void WriteJets();
    virtual void ResetJets();
    virtual TClonesArray* Jets() {return fJets;}
 private:
    virtual void BookLego();
    virtual void DumpLego();
    virtual void ResetMap();
    virtual Float_t PropagatePhi(Float_t pt, Float_t charge, Bool_t& curls);
 protected:
    TClonesArray*                  fJets;            //! List of Jets
    TH2F*                          fLego;            //! Lego Histo
    AliEMCALJet*                   fJetT[10];        //! Jet temporary storage
    AliEMCALHadronCorrection*      fHadronCorrector; //! Pointer to hadronic correction
    Int_t                          fHCorrection;     //  Hadron correction flag
    Int_t                          fDebug;           //! Debug flag
    Float_t                        fConeRadius;      //  Cone radius
    Float_t                        fPtCut;           //  Pt cut on charged tracks
    Float_t                        fEtSeed;          //  Min. Et for seed
    Float_t                        fMinJetEt;        //  Min Et of jet
    Float_t                        fMinCellEt;       //  Min Et in one cell
    Float_t                        fSamplingF;
    Bool_t                         fSmear;           //  Flag for momentum smearing
    Bool_t                         fEffic;           //  Flag for efficiency simulation
    Bool_t                         fK0N;             //  Flag for efficiency simulation
    Int_t                          fNjets;           //! Number of Jets
    Float_t                        fDeta;            //! eta cell size 
    Float_t                        fDphi;            //! phi cell size
    Int_t                          fNcell;           //! number of cells
    Int_t                          fNtot;            //! total number of cells
    Int_t                          fNbinEta;         //! number of cells in eta
    Int_t                          fNbinPhi;         //! number of cells in phi
    Float_t                        fEtCell[30000];   //! Cell Energy
    Float_t                        fEtaCell[30000];  //! Cell eta
    Float_t                        fPhiCell[30000];  //! Cell phi
    Int_t*                         fTrackList;       //! List of selected tracks
    Int_t                          fNt;              //! number of tracks
    Float_t*                       fPtT;             //! Pt   of tracks in jet cone
    Float_t*                       fEtaT;            //! Eta  of tracks in jet cone
    Float_t*                       fPhiT;            //! Phi  of tracks in jet cone
    ClassDef(AliEMCALJetFinder,2)        // JetFinder for EMCAL
}
;
#endif // ALIEMCALJetFinder_H


