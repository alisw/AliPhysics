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
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, const Float_t centralityF) const;

  // Getters
  THnSparseF *GetEventCount()   const {return fEventCount;}
  THnSparseF *GetRecEventHist() const   {return fRecEventHist;}
  THnSparseF *GetMCEventHist()  const   {return fMCEventHist;}
  THnSparseF *GetRecMCEventHist() const {return fRecMCEventHist;}

  //
  THnSparseF *GetRecMCTrackHist() const {return fRecMCTrackHist;}


  // THx event histograms getters

  TH2D *GetTrigVertCount() const {return fTrigVertCount;}
  TH3D *GetRecEventXYZ() const {return fRecEventXYZ;}
  TH3D *GetRecEventXYMult() const {return fRecEventXYMult;}
  TH3D *GetRecEventXZMult() const {return fRecEventXZMult;}
  TH3D *GetRecEventYZMult()  const{return fRecEventYZMult;}
  TH3D *GetRecEventZResZMult() const {return fRecEventZResZMult;}

  TH3D *GetMCEventXYZ() const {return fMCEventXYZ;}

  TH3D *GetRecMCEventDXDYDZ() const {return fRecMCEventDXDYDZ;}
  TH3D *GetRecMCEventDXDYMult() const {return fRecMCEventDXDYMult;}
  TH3D *GetRecMCEventDXDZMult() const {return fRecMCEventDXDZMult;}
  TH3D *GetRecMCEventDYDZMult() const {return fRecMCEventDYDZMult;}

  // THx track histograms

  TH3D *GetCrossRowsEtaPt() const {return fCrossRowsEtaPt;}
  TH3D *GetChi2PerClustEtaPt() const {return fChi2PerClustEtaPt;}
  TH3D *GetCrossRowsOverFindEtaPt() const {return fCrossRowsOverFindEtaPt;}
  TH3D *GetFracSharedClustEtaPt() const {return fFracSharedClustEtaPt;}
  TH3D *GetDCAyEtaPt() const {return fDCAyEtaPt;}
  TH3D *GetDCAzEtaPt() const {return fDCAzEtaPt;}
  TH3D *GetDCAyCenPt() const {return fDCAyCenPt;}

  TH3D *GetCrossRowsPhiPt() const {return fCrossRowsPhiPt;} 
  TH3D *GetChi2PerClustPhiPt() const {return fChi2PerClustPhiPt;}
  TH3D *GetCrossRowsOverFindPhiPt() const {return fCrossRowsOverFindPhiPt;}
  TH3D *GetFracSharedClustPhiPt() const {return fFracSharedClustPhiPt;}
  TH3D *GetDCAyPhiPt() const {return fDCAyPhiPt;}
  TH3D *GetDCAzPhiPt() const {return fDCAzPhiPt;}

  TH3D *GetCrossRowsEtaPhi() const {return fCrossRowsEtaPhi;}
  TH3D *GetChi2PerClustEtaPhi() const {return fChi2PerClustEtaPhi;}
  TH3D *GetCrossRowsOverFindEtaPhi() const {return fCrossRowsOverFindEtaPhi;}
  TH3D *GetFracSharedClustEtaPhi() const {return fFracSharedClustEtaPhi;}
  TH3D *GetDCAyEtaPhi() const {return fDCAyEtaPhi;}
  TH3D *GetDCAzEtaPhi() const {return fDCAzEtaPhi;}

  TH3D *GetDCAyEtaPtMCPrim() const {return fDCAyEtaPtMCPrim;}
  TH3D *GetDCAzEtaPtMCPrim() const {return fDCAzEtaPtMCPrim;}
  TH3D *GetDCAyPhiPtMCPrim() const {return fDCAyPhiPtMCPrim;}
  TH3D *GetDCAzPhiPtMCPrim() const {return fDCAzPhiPtMCPrim;}
  TH3D *GetDCAyEtaPhiMCPrim() const {return fDCAyEtaPhiMCPrim;}
  TH3D *GetDCAzEtaPhiMCPrim() const {return fDCAzEtaPhiMCPrim;}
  TH3D *GetDCAyCenPtMCPrim() const {return fDCAyCenPtMCPrim;}

  TH3D *GetDCAyEtaPtMCSec() const {return fDCAyEtaPtMCSec;}
  TH3D *GetDCAzEtaPtMCSec() const {return fDCAzEtaPtMCSec;}
  TH3D *GetDCAyPhiPtMCSec() const {return fDCAyPhiPtMCSec;}
  TH3D *GetDCAzPhiPtMCSec() const {return fDCAzPhiPtMCSec;}
  TH3D *GetDCAyEtaPhiMCSec() const {return fDCAyEtaPhiMCSec;}
  TH3D *GetDCAzEtaPhiMCSec() const {return fDCAzEtaPhiMCSec;}
  TH3D *GetDCAyCenPtMCSec() const {return fDCAyCenPtMCSec;}

  TH3D *GetDCAyEtaPtMCSecDecays() const {return fDCAyEtaPtMCSecDecays;}
  TH3D *GetDCAzEtaPtMCSecDecays() const {return fDCAzEtaPtMCSecDecays;}
  TH3D *GetDCAyPhiPtMCSecDecays() const {return fDCAyPhiPtMCSecDecays;}
  TH3D *GetDCAzPhiPtMCSecDecays() const {return fDCAzPhiPtMCSecDecays;}
  TH3D *GetDCAyEtaPhiMCSecDecays() const {return fDCAyEtaPhiMCSecDecays;}
  TH3D *GetDCAzEtaPhiMCSecDecays() const {return fDCAzEtaPhiMCSecDecays;}
  TH3D *GetDCAyCenPtMCSecDecays() const {return fDCAyCenPtMCSecDecays;}

  TH3D *GetDCAyEtaPtMCSecDecaysK0s() const {return fDCAyEtaPtMCSecDecaysK0s;}
  TH3D *GetDCAzEtaPtMCSecDecaysK0s() const {return fDCAzEtaPtMCSecDecaysK0s;}
  TH3D *GetDCAyPhiPtMCSecDecaysK0s() const {return fDCAyPhiPtMCSecDecaysK0s;}
  TH3D *GetDCAzPhiPtMCSecDecaysK0s() const {return fDCAzPhiPtMCSecDecaysK0s;} 
  TH3D *GetDCAyEtaPhiMCSecDecaysK0s() const {return fDCAyEtaPhiMCSecDecaysK0s;}
  TH3D *GetDCAzEtaPhiMCSecDecaysK0s() const {return fDCAzEtaPhiMCSecDecaysK0s;}
  TH3D *GetDCAyCenPtMCSecDecaysK0s() const {return fDCAyCenPtMCSecDecaysK0s;}

  TH3D *GetDCAyEtaPtMCSecDecaysLambda() const {return fDCAyEtaPtMCSecDecaysLambda;}
  TH3D *GetDCAzEtaPtMCSecDecaysLambda() const {return fDCAzEtaPtMCSecDecaysLambda;}
  TH3D *GetDCAyPhiPtMCSecDecaysLambda() const {return fDCAyPhiPtMCSecDecaysLambda;}
  TH3D *GetDCAzPhiPtMCSecDecaysLambda() const {return fDCAzPhiPtMCSecDecaysLambda;}
  TH3D *GetDCAyEtaPhiMCSecDecaysLambda() const {return fDCAyEtaPhiMCSecDecaysLambda;}
  TH3D *GetDCAzEtaPhiMCSecDecaysLambda() const {return fDCAzEtaPhiMCSecDecaysLambda;}
  TH3D *GetDCAyCenPtMCSecDecaysLambda() const {return fDCAyCenPtMCSecDecaysLambda;}

  TH3D *GetDCAyEtaPtMCSecMaterial() const {return fDCAyEtaPtMCSecMaterial;} 
  TH3D *GetDCAzEtaPtMCSecMaterial() const {return fDCAzEtaPtMCSecMaterial;}
  TH3D *GetDCAyPhiPtMCSecMaterial() const {return fDCAyPhiPtMCSecMaterial;}
  TH3D *GetDCAzPhiPtMCSecMaterial() const {return fDCAzPhiPtMCSecMaterial;}
  TH3D *GetDCAyEtaPhiMCSecMaterial() const {return fDCAyEtaPhiMCSecMaterial;}
  TH3D *GetDCAzEtaPhiMCSecMaterial() const {return fDCAzEtaPhiMCSecMaterial;}
  TH3D *GetDCAyCenPtMCSecMaterial() const {return fDCAyCenPtMCSecMaterial;} 

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
  TH3D *fDCAyCenPt;                 // DCAy:centrality:pt

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
  TH3D *fDCAyCenPtMCPrim;             // DCAy:centrality:pt for primary particles

  TH3D *fDCAyEtaPtMCSec;             // DCAy:eta:pt for secondary particles
  TH3D *fDCAzEtaPtMCSec;             // DCAz:eta:pt for secondary particles
  TH3D *fDCAyPhiPtMCSec;             // DCAy:phi:pt for secondary particles
  TH3D *fDCAzPhiPtMCSec;             // DCAz:phi:pt for secondary particles
  TH3D *fDCAyEtaPhiMCSec;            // DCAy:eta:phi for secondary particles
  TH3D *fDCAzEtaPhiMCSec;            // DCAz:eta:phi for secondary particles
  TH3D *fDCAyCenPtMCSec;             // DCAy:centrality:pt for secondary particles

  TH3D *fDCAyEtaPtMCSecDecays;     // DCAy:eta:pt for secondary particles from decays in MC
  TH3D *fDCAzEtaPtMCSecDecays;     // DCAz:eta:pt for secondary particles from decays in MC
  TH3D *fDCAyPhiPtMCSecDecays;     // DCAy:phi:pt for secondary particles from decays in MC
  TH3D *fDCAzPhiPtMCSecDecays;     // DCAz:phi:pt for secondary particles from decays in MC
  TH3D *fDCAyEtaPhiMCSecDecays;    // DCAy:eta:phi for secondary particles from decays in MC
  TH3D *fDCAzEtaPhiMCSecDecays;    // DCAz:eta:phi for secondary particles from decays in MC
  TH3D *fDCAyCenPtMCSecDecays;     // DCAy:centrality:pt for secondary particles from decays in MC

  TH3D *fDCAyEtaPtMCSecDecaysK0s;     // DCAy:eta:pt for secondary particles from decays K0s in MC
  TH3D *fDCAzEtaPtMCSecDecaysK0s;     // DCAz:eta:pt for secondary particles from decays K0s in MC
  TH3D *fDCAyPhiPtMCSecDecaysK0s;     // DCAy:phi:pt for secondary particles from decays K0s in MC
  TH3D *fDCAzPhiPtMCSecDecaysK0s;     // DCAz:phi:pt for secondary particles from decays K0s in MC
  TH3D *fDCAyEtaPhiMCSecDecaysK0s;    // DCAy:eta:phi for secondary particles from decays K0s in MC
  TH3D *fDCAzEtaPhiMCSecDecaysK0s;    // DCAz:eta:phi for secondary particles from decays K0s in MC
  TH3D *fDCAyCenPtMCSecDecaysK0s;     // DCAy:centrality:pt for secondary particles from decays K0s in MC

  TH3D *fDCAyEtaPtMCSecDecaysLambda;     // DCAy:eta:pt for secondary particles from decays Lambda in MC
  TH3D *fDCAzEtaPtMCSecDecaysLambda;     // DCAz:eta:pt for secondary particles from decays Lambda in MC
  TH3D *fDCAyPhiPtMCSecDecaysLambda;     // DCAy:phi:pt for secondary particles from decays Lambda in MC
  TH3D *fDCAzPhiPtMCSecDecaysLambda;     // DCAz:phi:pt for secondary particles from decays Lambda in MC
  TH3D *fDCAyEtaPhiMCSecDecaysLambda;    // DCAy:eta:phi for secondary particles from decays Lambda in MC
  TH3D *fDCAzEtaPhiMCSecDecaysLambda;    // DCAz:eta:phi for secondary particles from decays Lambda in MC
  TH3D *fDCAyCenPtMCSecDecaysLambda;     // DCAy:centrality:pt for secondary particles from decays Lambda in MC

  TH3D *fDCAyEtaPtMCSecMaterial;   // DCAy:eta:pt for secondary particles from material in MC
  TH3D *fDCAzEtaPtMCSecMaterial;   // DCAz:eta:pt for secondary particles from material in MC
  TH3D *fDCAyPhiPtMCSecMaterial;   // DCAy:phi:pt for secondary particles from material in MC
  TH3D *fDCAzPhiPtMCSecMaterial;   // DCAz:phi:pt for secondary particles from material in MC
  TH3D *fDCAyEtaPhiMCSecMaterial;  // DCAy:eta:phi for secondary particles from material in MC
  TH3D *fDCAzEtaPhiMCSecMaterial;  // DCAz:eta:phi for secondary particles from material in MC
  TH3D *fDCAyCenPtMCSecMaterial;   // DCAy:centrality:pt for secondary particles from material in MC

  AlidNdPtCutAnalysis(const AlidNdPtCutAnalysis&); // not implemented
  AlidNdPtCutAnalysis& operator=(const AlidNdPtCutAnalysis&); // not implemented

  ClassDef(AlidNdPtCutAnalysis,4);
};

#endif
