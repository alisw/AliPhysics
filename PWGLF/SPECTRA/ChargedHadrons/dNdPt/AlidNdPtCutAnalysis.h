#ifndef ALIDNDPTCUTANALYSIS_H
#define ALIDNDPTCUTANALYSIS_H

//------------------------------------------------------------------------------
// AlidNdPtCutAnalysis class to determine
// cuts to be used for dNdPt analysis.
//
// Author: J.Otwinowski 04/11/2008
// Modified: J.Otwinowski 21/10/2015
//------------------------------------------------------------------------------

class iostream;
class TFile;
class TCint;
class TProfile;
class TFolder;
class TObjArray;
class TString;
class THnSparse;

class AliESDtrackCuts;
class AliVertexerTracks;
class AliESD;
class AliESDfriend;
class AliESDfriendTrack;

#include "AlidNdPt.h"

class AlidNdPtCutAnalysis : public AlidNdPt {
public :
  AlidNdPtCutAnalysis();
  AlidNdPtCutAnalysis(Char_t* name, Char_t* title);
  ~AlidNdPtCutAnalysis();

  // Init data members
  virtual void Init();

  // Process events
  virtual void Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0);

  // Merge output objects (needed by PROOF)
  virtual Long64_t Merge(TCollection* const list);

  // Set and get sparse histogram filling flag
  void SetFillSparseHisto(Bool_t bFillSparseHisto) {fFillSparseHisto = bFillSparseHisto;};
  Bool_t GetFillSparseHisto() {return fFillSparseHisto;};

  // Analyse output histograms
  virtual void Analyse();

  // Export objects to folder
  virtual TFolder *ExportToFolder(TObjArray * const array=0);

  // Get analysis folder
  TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderdNdPtAnalysis",TString title = "Analysed dNdPt histograms");

  // Fill histograms
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack) const;

  // Getters
  THnSparseF *GetEventCount()   const {return fEventCount;}
  THnSparseF *GetRecEventHist() const   {return fRecEventHist;}
  THnSparseF *GetMCEventHist()  const   {return fMCEventHist;}
  THnSparseF *GetRecMCEventHist() const {return fRecMCEventHist;}

  //
  THnSparseF *GetRecMCTrackHist() const {return fRecMCTrackHist;}

private:

  // histogram fill flag
  Bool_t fFillSparseHisto;

  // analysis folder
  TFolder *fAnalysisFolder; // folder for analysed histograms

  //
  // THnSparse event histograms
  //
  THnSparseF *fEventCount; //-> trig, trig + vertex

  THnSparseF *fRecEventHist;   //-> Xv:Yv:Zv:ResZv:Mult
  THnSparseF *fMCEventHist;    //-> mcXv:mcYv:mcZv
  THnSparseF *fRecMCEventHist; //-> Xv-mcXv:Yv-mcYv:Zv-mcZv:Mult

  //
  // THnSparse track histograms
  //
  //THnSparseF *fRecMCTrackHist; //-> nClust:chi2PerClust:nClust/nFindableClust:DCAy:DCAz:eta:phi:pt:kinkIdx:isPrim:polarity
  THnSparseF *fRecMCTrackHist; //-> nCrossRows:chi2PerClust:nCrossRows/nFindableClust:fracSharedClust:DCAy:DCAz:eta:phi:pt:hasStrangeMother:isFromMaterial:isPrim:charge

  // THx event histograms

  TH2D *fTrigVertCount;         //-> trig, trig + vertex

  TH3D *fRecEventXYZ;           //-> Xv:Yv:Zv
  TH3D *fRecEventXYMult;        //-> Xv:Yv:Mult
  TH3D *fRecEventXZMult;        //-> Xv:Zv:Mult
  TH3D *fRecEventYZMult;        //-> Yv:Zv:Mult
  TH3D *fRecEventZResZMult;     //-> Zv:ResZv:Mult

  TH3D *fMCEventXYZ;            //-> mcXv:mcYv:mcZv

  TH3D *fRecMCEventDXDYDZ;      //-> Xv-mcXv:Yv-mcYv:Zv-mcZv
  TH3D *fRecMCEventDXDYMult;    //-> Xv-mcXv:Yv-mcYv:Mult
  TH3D *fRecMCEventDXDZMult;    //-> Xv-mcXv:Zv-mcZv:Mult
  TH3D *fRecMCEventDYDZMult;    //-> Yv-mcYv:Zv-mcZv:Mult

  // THx track histograms

  TH3D *fCrossRowsEtaPt;            // nCrossRows:eta:pt
  TH3D *fChi2PerClustEtaPt;         // chi2PerClust:eta:pt
  TH3D *fCrossRowsOverFindEtaPt;    // nCrossRows/nFindableClust:eta:pt
  TH3D *fFracSharedClustEtaPt;      // fracSharedClust:eta:pt
  TH3D *fDCAyEtaPt;                 // DCAy:eta:pt
  TH3D *fDCAzEtaPt;                 // DCAz:eta:pt

  TH3D *fCrossRowsPhiPt;            // nCrossRows:phi:pt
  TH3D *fChi2PerClustPhiPt;         // chi2PerClust:phi:pt
  TH3D *fCrossRowsOverFindPhiPt;    // nCrossRows/nFindableClust:phi:pt
  TH3D *fFracSharedClustPhiPt;      // fracSharedClust:phi:pt
  TH3D *fDCAyPhiPt;                 // DCAy:phi:pt
  TH3D *fDCAzPhiPt;                 // DCAz:phi:pt

  TH3D *fCrossRowsEtaPhi;           // nCrossRows:eta:phi
  TH3D *fChi2PerClustEtaPhi;        // chi2PerClust:eta:phi
  TH3D *fCrossRowsOverFindEtaPhi;   // nCrossRows/nFindableClust:eta:phi
  TH3D *fFracSharedClustEtaPhi;     // fracSharedClust:eta:phi
  TH3D *fDCAyEtaPhi;                // DCAy:eta:phi
  TH3D *fDCAzEtaPhi;                // DCAz:eta:phi

  TH3D *fDCAyEtaPtMCPrim;             // DCAy:eta:pt for primary particles
  TH3D *fDCAzEtaPtMCPrim;             // DCAz:eta:pt for primary particles
  TH3D *fDCAyPhiPtMCPrim;             // DCAy:phi:pt for primary particles
  TH3D *fDCAzPhiPtMCPrim;             // DCAz:phi:pt for primary particles
  TH3D *fDCAyEtaPhiMCPrim;            // DCAy:eta:phi for primary particles
  TH3D *fDCAzEtaPhiMCPrim;            // DCAz:eta:phi for primary particles

  TH3D *fDCAyEtaPtMCSec;             // DCAy:eta:pt for secondary particles
  TH3D *fDCAzEtaPtMCSec;             // DCAz:eta:pt for secondary particles
  TH3D *fDCAyPhiPtMCSec;             // DCAy:phi:pt for secondary particles
  TH3D *fDCAzPhiPtMCSec;             // DCAz:phi:pt for secondary particles
  TH3D *fDCAyEtaPhiMCSec;            // DCAy:eta:phi for secondary particles
  TH3D *fDCAzEtaPhiMCSec;            // DCAz:eta:phi for secondary particles

  TH3D *fDCAyEtaPtSecMCDecays;     // DCAy:eta:pt for secondary particles from decays in MC
  TH3D *fDCAzEtaPtSecMCDecays;     // DCAz:eta:pt for secondary particles from decays in MC
  TH3D *fDCAyPhiPtSecMCDecays;     // DCAy:phi:pt for secondary particles from decays in MC
  TH3D *fDCAzPhiPtSecMCDecays;     // DCAz:phi:pt for secondary particles from decays in MC
  TH3D *fDCAyEtaPhiSecMCDecays;    // DCAy:eta:phi for secondary particles from decays in MC
  TH3D *fDCAzEtaPhiSecMCDecays;    // DCAz:eta:phi for secondary particles from decays in MC

  TH3D *fDCAyEtaPtSecMCDecaysK0s;     // DCAy:eta:pt for secondary particles from decays K0s in MC
  TH3D *fDCAzEtaPtSecMCDecaysK0s;     // DCAz:eta:pt for secondary particles from decays K0s in MC
  TH3D *fDCAyPhiPtSecMCDecaysK0s;     // DCAy:phi:pt for secondary particles from decays K0s in MC
  TH3D *fDCAzPhiPtSecMCDecaysK0s;     // DCAz:phi:pt for secondary particles from decays K0s in MC
  TH3D *fDCAyEtaPhiSecMCDecaysK0s;    // DCAy:eta:phi for secondary particles from decays K0s in MC
  TH3D *fDCAzEtaPhiSecMCDecaysK0s;    // DCAz:eta:phi for secondary particles from decays K0s in MC

  TH3D *fDCAyEtaPtSecMCDecaysLambda;     // DCAy:eta:pt for secondary particles from decays Lambda in MC
  TH3D *fDCAzEtaPtSecMCDecaysLambda;     // DCAz:eta:pt for secondary particles from decays Lambda in MC
  TH3D *fDCAyPhiPtSecMCDecaysLambda;     // DCAy:phi:pt for secondary particles from decays Lambda in MC
  TH3D *fDCAzPhiPtSecMCDecaysLambda;     // DCAz:phi:pt for secondary particles from decays Lambda in MC
  TH3D *fDCAyEtaPhiSecMCDecaysLambda;    // DCAy:eta:phi for secondary particles from decays Lambda in MC
  TH3D *fDCAzEtaPhiSecMCDecaysLambda;    // DCAz:eta:phi for secondary particles from decays Lambda in MC

  TH3D *fDCAyEtaPtSecMCMaterial;   // DCAy:eta:pt for secondary particles from material in MC
  TH3D *fDCAzEtaPtSecMCMaterial;   // DCAz:eta:pt for secondary particles from material in MC
  TH3D *fDCAyPhiPtSecMCMaterial;   // DCAy:phi:pt for secondary particles from material in MC
  TH3D *fDCAzPhiPtSecMCMaterial;   // DCAz:phi:pt for secondary particles from material in MC
  TH3D *fDCAyEtaPhiSecMCMaterial;  // DCAy:eta:phi for secondary particles from material in MC
  TH3D *fDCAzEtaPhiSecMCMaterial;  // DCAz:eta:phi for secondary particles from material in MC

  AlidNdPtCutAnalysis(const AlidNdPtCutAnalysis&); // not implemented
  AlidNdPtCutAnalysis& operator=(const AlidNdPtCutAnalysis&); // not implemented

  ClassDef(AlidNdPtCutAnalysis,3);
};

#endif
