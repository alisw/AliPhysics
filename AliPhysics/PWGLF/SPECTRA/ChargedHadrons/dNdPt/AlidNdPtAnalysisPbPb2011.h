
#ifndef AlidNdPtAnalysisPbPb2011_H
#define AlidNdPtAnalysisPbPb2011_H


//------------------------------------------------------------------------------
// AlidNdPtAnalysisPbPb2011 class used for dNdPt analysis.in PbPb collision 
// 
// Author: J.Otwinowski 04/11/2008 
// last change: 2012-09-04 by P. Luettig
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
class AlidNdPtHelper;
class AliTriggerAnalysis;

#include "AlidNdPt.h"
#include "TObjString.h"

class AlidNdPtAnalysisPbPb2011 : public AlidNdPt {
public :
  AlidNdPtAnalysisPbPb2011(); 
  AlidNdPtAnalysisPbPb2011(Char_t* name, Char_t* title);
  ~AlidNdPtAnalysisPbPb2011();

  // Init data members
  virtual void Init();

  // Process events
  virtual void Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0);

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* const list);

  // Analyse output histograms 
  virtual void Analyse();

  // Export objects to folder
  virtual TFolder *ExportToFolder(TObjArray * const array=0);

  // Get analysis folder
  TFolder* GetAnalysisFolder() const {return fAnalysisFolder;}

  // Fill control histograms
  void SetHistogramsOn(Bool_t histOn=kTRUE) {fHistogramsOn = histOn;}
  Bool_t IsHistogramsOn() const {return fHistogramsOn;}
    
  // Create folder for analysed histograms
  TFolder *CreateFolder(TString folder = "folderdNdPtAnalysis",TString title = "Analysed dNdPt histograms");
  
  // Set binning for Histograms (if not set default binning is used)
  void SetBinsMult(Int_t nbins, Double_t* edges) { Printf("[I] Setting Mult Bins"); fMultNbins = nbins; fBinsMult = CloneArray(nbins,edges); }
  void SetBinsPt(Int_t nbins, Double_t* edges) { Printf("[I] Setting pT Bins"); fPtNbins = nbins; fBinsPt = CloneArray(nbins,edges); }
  void SetBinsPtCorr(Int_t nbins, Double_t* edges) { Printf("[I] Setting pTcorr Bins"); fPtCorrNbins = nbins; fBinsPtCorr = CloneArray(nbins,edges); }
  void SetBinsEta(Int_t nbins, Double_t* edges) { Printf("[I] Setting Eta Bins"); fEtaNbins = nbins; fBinsEta = CloneArray(nbins,edges); }
  void SetBinsZv(Int_t nbins, Double_t* edges) { Printf("[I] Setting Zv Bins"); fZvNbins = nbins; fBinsZv= CloneArray(nbins,edges); }
  void SetBinsCentrality(Int_t nbins, Double_t* edges) { Printf("[I] Setting Cent Bins"); fCentralityNbins = nbins; fBinsCentrality = CloneArray(nbins,edges); }

  // Fill histograms
  void FillHistograms(AliESDtrack *const esdTrack, AliStack *const stack, AlidNdPtHelper::TrackObject trackObj, Float_t centralityF);
  void FillHistograms(AliStack *const stack, Int_t label, AlidNdPtHelper::TrackObject trackObj, Float_t centralityF);
  void FillHistograms(TObjArray *const allChargedTracks,Int_t *const labelsAll,Int_t multAll,Int_t *const labelsAcc,Int_t multAcc,Int_t *const labelsRec,Int_t multRec,  Float_t centralityF);

  // Getters
  THnSparseF *GetTrackPtCorrelationMatrix()   const {return fTrackPtCorrelationMatrix;}
  //
  THnSparseF *GetGenEventMatrix() const {return fGenEventMatrix;}
  THnSparseF *GetTriggerEventMatrix() const {return fTriggerEventMatrix;}
  THnSparseF *GetRecEventMatrix() const {return fRecEventMatrix;}
  // 
  THnSparseF *GetGenTrackEventMatrix() const {return fGenTrackEventMatrix;}
  THnSparseF *GetTriggerTrackEventMatrix() const {return fTriggerTrackEventMatrix;}
  THnSparseF *GetRecTrackEventMatrix() const {return fRecTrackEventMatrix;}
  //
  THnSparseF *GetGenTrackMatrix() const {return fGenTrackMatrix;}
  THnSparseF *GetGenPrimTrackMatrix() const {return fGenPrimTrackMatrix;}
  THnSparseF *GetRecPrimTrackMatrix() const {return fRecPrimTrackMatrix;}

  THnSparseF *GetRecTrackMatrix() const {return fRecTrackMatrix;}
  THnSparseF *GetRecSecTrackMatrix() const {return fRecSecTrackMatrix;}
  THnSparseF *GetRecMultTrackMatrix() const {return fRecMultTrackMatrix;}
  //
  // control histograms
  //
  THnSparseF *GetMCEventHist1() const {return fMCEventHist1;}
  THnSparseF *GetRecEventHist1() const {return fRecEventHist1;}
  THnSparseF *GetRecEventHist2() const {return fRecEventHist2;}


  THnSparseF *GetRecMCEventHist1() const {return fRecMCEventHist1;}
  THnSparseF *GetRecMCTrackHist1() const {return fRecMCTrackHist1;}

  THnSparseF *GetRecMCEventHist2() const {return fRecMCEventHist2;}

  THnSparseF *GetMCTrackHist1(Int_t i) const {return fMCTrackHist1[i];}
  THnSparseF *GetMCPrimTrackHist1(Int_t i) const {return fMCPrimTrackHist1[i];}
  THnSparseF *GetMCPrimTrackHist2(Int_t i) const {return fMCPrimTrackHist2[i];}
  THnSparseF *GetMCSecTrackHist1(Int_t i) const {return fMCSecTrackHist1[i];}

  THnSparseF *GetRecTrackHist1(Int_t i) const {return fRecTrackHist1[i];}
  THnSparseF *GetRecTrackHist2(Int_t i) const {return fRecTrackHist2[i];}
  THnSparseF *GetRecTrackMultHist1(Int_t i) const {return fRecTrackMultHist1[i];}


  THnSparseF *GetMCMultRecTrackHist1() const {return fMCMultRecTrackHist1;}

  THnSparseF *GetRecTrackHist3() const {return fRecTrackHist3;}
  
  TString GetCentralityEstimator() const {return fCentralityEstimator; }
  void SetCentralityEstimator(TString centEst="V0M") { fCentralityEstimator = centEst; }

private:

  // analysis folder 
  TFolder *fAnalysisFolder; // folder for analysed histograms
  Bool_t fHistogramsOn; // switch on/off filling of control histograms 

  // 
  // correlation matrices (histograms)
  //
  // rec. track pt vs true track pt correlation matrix for given eta
  THnSparseF *fTrackPtCorrelationMatrix; //-> Pt:mcPt:mcEta:centrality

  //
  // event level correction 
  //
  // all generated
  THnSparseF *fGenEventMatrix; //-> mcZv:multMB:centrality 

  // trigger bias corrections (fTriggerEventMatrix / fGenEventMatrix)
  THnSparseF *fTriggerEventMatrix; //-> mcZv:multMB:centrality

  // event vertex rec. eff correction (fRecEventMatrix / fTriggerEventMatrix)
  THnSparseF *fRecEventMatrix; //-> mcZv:multMB:centrality 

  // track-event level correction 
  THnSparseF *fGenTrackEventMatrix; //-> mcZv:mcPt:mcEta:centrality

  // trigger bias corrections (fTriggerTrackEventMatrix / fGenTrackEventMatrix)
  THnSparseF *fTriggerTrackEventMatrix; //-> mcZv:mcPt:mcEta:centrality

  // event vertex rec. corrections (fRecTrackEventMatrix / fTriggerTrackEventMatrix)
  THnSparseF *fRecTrackEventMatrix; //-> mcZv:Pt:mcEta:centrality

  //
  // track level correction 
  //
  // track rec. efficiency correction (fRecPrimTrackMatrix / fGenPrimTrackMatrix)
  THnSparseF *fGenTrackMatrix; //-> mcZv:mcPt:mcEta:centrality
  THnSparseF *fGenPrimTrackMatrix; //-> mcZv:mcPt:mcEta:centrality
  THnSparseF *fRecPrimTrackMatrix; //-> mcZv:mcPt:mcEta:centrality
  // secondary track contamination correction (fRecSecTrackMatrix / fRecTrackMatrix)
  THnSparseF *fRecTrackMatrix;    //-> mcZv:mcPt:mcEta:centrality
  THnSparseF *fRecSecTrackMatrix; //-> mcZv:mcPt:mcEta:centrality
  // multiple rec. track corrections (fRecMultTrackMatrix / fRecTrackMatrix)
  THnSparseF *fRecMultTrackMatrix; //-> mcZv:Pt:mcEta:centrality

  //
  // ESD and MC control analysis histograms
  //
  // THnSparse event histograms
  THnSparseF *fMCEventHist1;  //-> mcXv:mcYv:mcZv:centrality
  THnSparseF *fRecEventHist1; //-> Xv:Yv:Zv:centrality
  THnSparseF *fRecEventHist2; //-> Zv:multMB:mult:centrality
  THnSparseF *fRecMCEventHist1; //-> Xv-mcXv:Yv-mcYv:Zv-mcZv:centrality
  THnSparseF *fRecMCEventHist2; //-> Xv-mcXv:Zv-mcZv:mult:centrality

  // [0] - after charged track selection, [1] - after acceptance cuts, [2] - after esd track cuts
  THnSparseF *fMCTrackHist1[AlidNdPtHelper::kCutSteps];     //-> mcPt:mcEta:mcPhi:centrality
  THnSparseF *fMCPrimTrackHist1[AlidNdPtHelper::kCutSteps]; //-> mcPt:mcEta:pid:mech:mother:centrality
  THnSparseF *fMCPrimTrackHist2[AlidNdPtHelper::kCutSteps]; //-> pdg:mech:mother:centrality
  THnSparseF *fMCSecTrackHist1[AlidNdPtHelper::kCutSteps];  //-> mcPt:mcEta:pid:mech:mother:centrality

  THnSparseF *fRecTrackHist1[AlidNdPtHelper::kCutSteps];     //-> Pt:Eta:Phi:centrality
  THnSparseF *fRecTrackHist2[AlidNdPtHelper::kCutSteps];     //-> Zv:Pt:Eta:centrality
  THnSparseF *fRecTrackMultHist1[AlidNdPtHelper::kCutSteps]; //-> Pt:mult:centrality
  THnSparseF *fRecMCTrackHist1; //-> mcPt:mcEta:(Pt-mcPt)/mcPt:(Eta-mcEta):centrality

  //multple reconstructed tracks
  THnSparseF *fMCMultRecTrackHist1; //-> mcPt:mcEta:pid:centrality
  // track control histograms
  THnSparseF *fRecTrackHist3;  //-> nclust:chi2:Pt:Eta:Phi:centrality

  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis object;
  TString fCentralityEstimator;     // use centrality can be "VOM" (default), "FMD", "TRK", "TKL", "CL0", "CL1", "V0MvsFMD", "TKLvsV0M", "ZEMvsZDC"
   
  //binning for THNsparse
  Int_t fMultNbins;
  Int_t fPtNbins;
  Int_t fPtCorrNbins;
  Int_t fEtaNbins;
  Int_t fZvNbins;
  Int_t fCentralityNbins;
  Double_t* fBinsMult; //[fMultNbins]
  Double_t* fBinsPt; //[fPtNbins]
  Double_t* fBinsPtCorr; //[fPtCorrNbins]
  Double_t* fBinsEta; //[fEtaNbins]
  Double_t* fBinsZv; //[fZvNbins]
  Double_t* fBinsCentrality; //[fCentralityNbins]
  
  Bool_t fIsInit;
  

  AlidNdPtAnalysisPbPb2011(const AlidNdPtAnalysisPbPb2011&); // not implemented
  AlidNdPtAnalysisPbPb2011& operator=(const AlidNdPtAnalysisPbPb2011&); // not implemented  

  ClassDef(AlidNdPtAnalysisPbPb2011,6);
};

#endif
