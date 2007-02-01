#ifndef ALIEMCALJETFINDER_H
#define ALIEMCALJETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id$ */
//*-- Author:
//*-- Andreas Morsch (CERN)

#include <TFile.h>
#include <TString.h>
#include <TTask.h>

class TClonesArray;
class TH2F;
class TH1F;
class TCanvas;
class TList;

#include "AliEMCALJet.h"

class AliEMCALHadronCorrection;

class AliEMCALJetFinder : public TTask {
  friend class AliEMCALJetMicroDst; //PH Temporary solution
 public:
    AliEMCALJetFinder();
    AliEMCALJetFinder(const char* name, const char *title);
    virtual ~AliEMCALJetFinder();

    AliEMCALJetFinder (const AliEMCALJetFinder&);
    AliEMCALJetFinder & operator = (const AliEMCALJetFinder & ) {
      Fatal("operator =", "not implemented") ;
      return *this ;
    }

    virtual void  Init();
    virtual void  Find(Int_t ncell, Int_t ncelltot, Float_t etc[30000], 
		      Float_t etac[30000], Float_t phic[30000],
		      Float_t minmove, Float_t maxmove, Int_t mode,
		      Float_t precbg, Int_t ierror);
    virtual void  Find();
    virtual void  FindChargedJet();
    virtual void  FindTracksInJetCone();
    virtual void  Test();
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
    virtual void SetSamplingFraction(Float_t par) {fSamplingF = par;}
    virtual void SetEnergyWeightingFlag(Bool_t flag) {fWeightingMethod = flag;}
    virtual void SetEMCALWeight(Float_t par) {fEMCALWeight = par;}
    virtual void SetTrackWeight(Float_t par) {fTrackWeight = par;}
    virtual void SetIncludeK0andN(Bool_t flag = kFALSE) {fK0N = flag;}
    // Correction of hadronic energy
    virtual void SetHadronCorrector(AliEMCALHadronCorrection* corr)
	{fHadronCorrector = corr;}
    virtual void SetHadronCorrection(Int_t flag = 1) {fHCorrection = flag;}
    // PAI
    void SetWriteKey(Bool_t flag = kFALSE) {fWrite = flag;}
    void SetMode(Int_t mode = 0) {fMode = mode;}
    void SetMinMove(Float_t minMove = 0.05) {fMinMove = minMove;}
    void SetMaxMove(Float_t maxMove = 0.15) {fMaxMove = maxMove;}
    void SetPrecBg (Float_t precBg = 0.035) {fPrecBg = precBg;}
    void SetParametersForBgSubtraction
    (Int_t mode=0, Float_t minMove=0.05, Float_t maxMove=0.15, Float_t precBg=0.035);
    //    virtual void Print(Option_t* option="") const;    // *MENU*
    void  SetRandomBg(Bool_t flag) {fRandomBg = flag;}
    Bool_t GetWriteKey() const {return fWrite;}
  //AliEMCALJet* GetJetT() {return fJetT[0];}
    AliEMCALJet* GetJetT(Int_t n = 0) {return fJetT[n];}
    virtual void DrawHistsForTuning(Int_t mode=0);           // *MENU*
    virtual void PrintParameters(Int_t mode=0);              // *MENU*
    virtual const Char_t* GetFileNameForParameters(const char* dir="RES/");

    // Access to Results
    virtual Int_t Njets() const;
    virtual Float_t JetEnergy (Int_t count) const;
    virtual Float_t JetPhiL  (Int_t count) const;
    virtual Float_t JetPhiW  (Int_t count) const ;
    virtual Float_t JetEtaL (Int_t count) const ;  
    virtual Float_t JetEtaW (Int_t count) const;
    TH2F*   GetLego() const {return fLego;}
    TH2F*   GetLegoB() const {return fLegoB;}
    TH2F*   GetLegoEMCAL() const {return fhLegoEMCAL;}
    TH2F*   GetLegoTracks() const {return fhLegoTracks;}
    TH2F*   GethEff() const {return fhEff;}
    TH1F*   GetCellEt() const {return fhCellEt;}
    TH1F*   GetCellEMCALEt() const {return fhCellEMCALEt;}
    TH1F*   GetTrackPt() const {return fhTrackPt;}
    TH1F*   GetTrackPtBcut() const {return fhTrackPtBcut;}
    TList*  GetHistsList() const {return fHistsList;}
    Int_t   GetNChTpc() const {return fNChTpc;}
    Bool_t  GetEnergyWeightingFlag() const {return fWeightingMethod ;}
    Float_t GetEMCALWeight() const {return fEMCALWeight;}
    Float_t GetTrackWeight() const {return fTrackWeight;}
    void    DrawLego(const char *opt="lego");         // *MENU*
    void    DrawLegoEMCAL(const char *opt="lego");    // *MENU*
    void    DrawLegos();                          // *MENU*
    void    DrawLegoBackground(const char *opt="lego"); // *MENU*
    Bool_t  IsThisPartonsOrDiQuark(Int_t pdg);
    // I/O
    virtual void SetOutputFileName(const char* name) {fOutFileName = name;}
    virtual void FillFromHits(Int_t flag = 0);
    virtual void FillFromHitFlaggedTracks(Int_t flag = 0);
    virtual void FillFromDigits(Int_t flag = 0);
    virtual void FillFromTracks(Int_t flag = 0, Int_t ich = 0);
    virtual void FillFromParticles();
    virtual void FillFromPartons();

    virtual void SaveBackgroundEvent(const char *name="");
    virtual void InitFromBackground();
    virtual void AddJet(const AliEMCALJet& jet);
    virtual void WriteJets();
    virtual void ResetJets();
    virtual TClonesArray* Jets() const {return fJets;}
    const char* GetNameOfVariant();

    virtual Bool_t  IsFolder() const;
    virtual void Browse(TBrowser* b);



    
 protected:
    TString fBGFileName;				// file name for background
    Float_t			   fEMCALWeight;	// EMCal energy weighting
    Float_t			   fTrackWeight;	// Track energy weighting
    Bool_t                         fRandomBg;        //  Flag for Random Background 
    Bool_t                         fWrite;           // Key for writing
    Bool_t                         fWeightingMethod; // Key for writing
    TClonesArray*                  fJets;            //! List of Jets
    TH2F*                          fLego;            //! Lego Histo
    TH2F*                          fLegoB;           //! Lego Histo Backg
    TH2F*                          fhLegoTracks;     //! Lego for Tracks
    TH2F*                          fhLegoEMCAL;      //! Lego for EMCAL itself
    TH2F*                          fhLegoHadrCorr;   //! Lego for hadron correction
    TH2F*                          fhEff;            //! Hist. for controling eff.
    TH1F*                          fhCellEt;         //! Et distr. for cells from fLego
    TH1F*                          fhCellEMCALEt;    //! Et distr. for cells from fLegoEMCAL
    TH1F*                          fhTrackPt;        //! Pt distr. for charge particles
    TH1F*                          fhTrackPtBcut;    //! Pt distr. for charge particles + cut due to magnetic field
    TH1F*                          fhChPartMultInTpc;//! Ch. part. multiplicity in TPC acceptance
    TH1F*                          fhSinTheta;       //! sin(theta)
    TCanvas*                       fC1;              //! first canvas for drawing
    TList*                         fHistsList;       //! List of hists - 4-mar-2002
    AliEMCALJet*                   fJetT[10];        //! Jet temporary storage
    AliEMCALHadronCorrection*      fHadronCorrector; //! Pointer to hadronic correction
    Int_t                          fHCorrection;     //  Hadron correction flag
    Int_t                          fDebug;           //! Debug flag
    Int_t                          fBackground;      //! Background flag
    Float_t                        fConeRadius;      //  Cone radius
    Float_t                        fPtCut;           //  Pt cut on charged tracks
    Float_t                        fEtSeed;          //  Min. Et for seed
    Float_t                        fMinJetEt;        //  Min Et of jet
    Float_t                        fMinCellEt;       //  Min Et in one cell
    Float_t                        fSamplingF;       //  Sampling Fraction
    Bool_t                         fSmear;           //  Flag for momentum smearing
    Bool_t                         fEffic;           //  Flag for efficiency simulation
    Bool_t                         fK0N;             //  Flag for efficiency simulation
    Int_t                          fNjets;           //! Number of Jetsp
    Float_t                        fDeta;            //! eta cell size 
    Float_t                        fDphi;            //! phi cell size
    Int_t                          fNcell;           //! number of cells
    Int_t                          fNtot;            //! total number of cells
    Int_t                          fNbinEta;         //! number of cells in eta
    Int_t                          fNbinPhi;         //! number of cells in phi
    Float_t                        fEtaMin;          //! minimum eta  
    Float_t                        fEtaMax;          //! maximum eta
    Float_t                        fPhiMin;          //! minimun phi
    Float_t                        fPhiMax;          //! maximum phi
    Float_t                        fEtCell[30000];   //! Cell Energy
    Float_t                        fEtaCell[30000];  //! Cell eta
    Float_t                        fPhiCell[30000];  //! Cell phi
    Int_t                          fNt;              //! number of tracks
    Int_t                          fNChTpc;          //! number of ch.part in TPC

    Int_t                          fNtS;             //! number of tracks selected
    Int_t*                         fTrackList;       //! List of selected tracks
    Float_t*                       fPtT;             //! Pt   of tracks 
    Float_t*                       fEtaT;            //! Eta  of tracks
    Float_t*                       fPhiT;            //! Phi  of tracks
    Int_t*                         fPdgT;            //! PDG code of tracks
 
    Int_t                          fNtB;             //! number of tracks in Bg
    Int_t*                         fTrackListB;      //! List of selected tracks in Bg
    Float_t*                       fPtB;             //! Pt   of tracks in Bg
    Float_t*                       fEtaB;            //! Eta  of tracks in Bg
    Float_t*                       fPhiB;            //! Phi  of tracks in Bg
    Int_t*                         fPdgB;            //! PDG  of tracks in Bg

    // parameter for jet_finder_ua1
    Int_t                          fMode;            // key for BG subtraction
    Float_t                        fMinMove;         // min cone move 
    Float_t                        fMaxMove;         // max cone move
    Float_t                        fPrecBg;          // max value of change for BG (in %)
    Int_t                          fError;           // error variables 

    const char*                    fOutFileName;     //! Output file name
    TFile*                         fOutFile;         //! Output file
    TFile*                         fInFile;          //! Output file
    Int_t                          fEvent;           //! Processed event
 private:
    virtual void BookLego();
    Float_t WeightedJetEnergy(Float_t eta, Float_t phi);
    Float_t EMCALConeEnergy(Float_t eta, Float_t phi);
    Float_t TrackConeEnergy(Float_t eta, Float_t phi);
    virtual void DumpLego();
    virtual void ResetMap();
    virtual void RearrangeParticlesMemory(Int_t npart);
 public:
    virtual Float_t PropagatePhi(Float_t pt, Float_t charge, Bool_t& curls);

    ClassDef(AliEMCALJetFinder,5)                    // JetFinder for EMCAL
}
;
#endif // ALIEMCALJetFinder_H
