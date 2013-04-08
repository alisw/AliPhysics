// *************************************************************************
// * Task for Fragmentation Function Analysis in PWG4 Jet Task Force Train *
// *************************************************************************

#ifndef ALIANALYSISTASKFRAGMENTATIONFUNCTIONN_H
#define ALIANALYSISTASKFRAGMENTATIONFUNCTIONN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class AliESDEvent;
class AliAODEvent;
class AliAODJet;
class AliAODExtension;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class THnSparse; 
class TRandom3;
class TArrayS;
class AliAODTrack;

#include "AliAnalysisTaskSE.h"
#include "TAxis.h"
#include "THnSparse.h"
  
class AliAnalysisTaskFragmentationFunction : public AliAnalysisTaskSE {

 public:
  
  //----------------------------------------
  class AliFragFuncHistos : public TObject
  {
    
    public:
    
    AliFragFuncHistos(const char* name = "FFhistos", 
		      Int_t nJetPt = 0, Float_t jetPtMin = 0, Float_t jetPtMax = 0,
		      Int_t nPt = 0, Float_t ptMin = 0, Float_t ptMax = 0,
		      Int_t nXi = 0, Float_t xiMin = 0, Float_t xiMax = 0,
		      Int_t nZ  = 0, Float_t zMin  = 0, Float_t zMax  = 0);
    
    AliFragFuncHistos(const AliFragFuncHistos& copy);
    AliFragFuncHistos& operator=(const AliFragFuncHistos &o);
    virtual ~AliFragFuncHistos();
    
    virtual void DefineHistos();
    virtual void FillFF(Float_t trackPt, Float_t jetPt,Bool_t incrementJetPt, Float_t norm = 0, Bool_t scaleStrangeness = kFALSE, Float_t scaleFacStrangeness = 1.);

    virtual void AddToOutput(TList* list) const;

  private:

    Int_t   fNBinsJetPt; // FF histos bins
    Float_t fJetPtMin;   // FF histos limits
    Float_t fJetPtMax;   // FF histos limits
    Int_t   fNBinsPt;    // FF histos bins
    Float_t fPtMin;      // FF histos limits
    Float_t fPtMax;      // FF histos limits
    Int_t   fNBinsXi;    // FF histos bins
    Float_t fXiMin;      // FF histos limits
    Float_t fXiMax;      // FF histos limits
    Int_t   fNBinsZ;     // FF histos bins
    Float_t fZMin;       // FF histos limits
    Float_t fZMax;       // FF histos limits
  
    TH2F*   fh2TrackPt;   //! FF: track transverse momentum 
    TH2F*   fh2Xi;        //! FF: xi 
    TH2F*   fh2Z;         //! FF: z  
    TH1F*   fh1JetPt;     //! jet pt 

    TString fNameFF;      // histo names prefix
    
    ClassDef(AliFragFuncHistos, 1);
  };
  
  //----------------------------------------
  class AliFragFuncQAJetHistos : public TObject
  {

  public:
 
    AliFragFuncQAJetHistos(const char* name = "QAJethistos",
		Int_t nPt  = 0, Float_t ptMin  = 0, Float_t ptMax  = 0,
		Int_t nEta = 0, Float_t etaMin = 0, Float_t etaMax = 0,
		Int_t nPhi = 0, Float_t phiMin = 0, Float_t phiMax = 0);
      
    AliFragFuncQAJetHistos(const AliFragFuncQAJetHistos& copy);
    AliFragFuncQAJetHistos& operator=(const AliFragFuncQAJetHistos &o);
    virtual ~AliFragFuncQAJetHistos();
    virtual void DefineHistos();
    virtual void FillJetQA(Float_t eta, Float_t phi, Float_t pt);
    virtual void AddToOutput(TList* list) const;

  private:
    
    Int_t   fNBinsPt;    // jet QA histos bins
    Float_t fPtMin;      // jet QA histos limits
    Float_t fPtMax;      // jet QA histos limits
    Int_t   fNBinsEta;   // jet QA histos bins
    Float_t fEtaMin;     // jet QA histos limits
    Float_t fEtaMax;     // jet QA histos limits
    Int_t   fNBinsPhi;   // jet QA histos bins
    Float_t fPhiMin;     // jet QA histos limits
    Float_t fPhiMax;     // jet QA histos limits
    
    TH2F*   fh2EtaPhi;   //! jet phi vs eta 
    TH1F*   fh1Pt;       //! jet transverse momentum 
    TString fNameQAJ;    // histo names prefix
    
    ClassDef(AliFragFuncQAJetHistos, 1);
  };
  
  //----------------------------------------
  class AliFragFuncQATrackHistos : public TObject
  {

  public:

    AliFragFuncQATrackHistos(const char* name = "QATrackhistos", 
		  Int_t nPt  = 0, Float_t ptMin  = 0, Float_t ptMax  = 0,
		  Int_t nEta = 0, Float_t etaMin = 0, Float_t etaMax = 0,
		  Int_t nPhi = 0, Float_t phiMin = 0, Float_t phiMax = 0, 
		  Float_t ptThresh = 0);
    
    AliFragFuncQATrackHistos(const AliFragFuncQATrackHistos& copy);
    AliFragFuncQATrackHistos& operator=(const AliFragFuncQATrackHistos &o);
    virtual ~AliFragFuncQATrackHistos();
    virtual void DefineHistos();
    virtual void FillTrackQA(Float_t eta, Float_t phi, Float_t pt, Bool_t weightPt = kFALSE, Float_t norm = 0., Bool_t scaleStrangeness = kFALSE, Float_t scaleFacStrangeness = 1.);
    virtual void AddToOutput(TList* list) const;

  private:
    
    Int_t   fNBinsPt;    // track QA histos bins in pt
    Float_t fPtMin;      // track QA histos limits in pt
    Float_t fPtMax;      // track QA histos limits in pt
    Int_t   fNBinsEta;   // track QA histos bins in eta
    Float_t fEtaMin;     // track QA histos limits in eta
    Float_t fEtaMax;     // track QA histos limits in eta
    Int_t   fNBinsPhi;   // track QA histos bins in phi
    Float_t fPhiMin;     // track QA histos limits in phi
    Float_t fPhiMax;     // track QA histos limits in phi

    Float_t fHighPtThreshold; //  high pt track phi vs eta distribution

    TH2F*   fh2EtaPhi;        //! track phi vs eta 
    TH1F*   fh1Pt;            //! track transverse momentum 
    TH2F*   fh2HighPtEtaPhi;  //! phi vs eta for high pt (>fgHighPtThreshold) tracks
    TH2F*   fh2PhiPt;         //! track phi vs pt

    TString fNameQAT;         // histo names prefix
    
    ClassDef(AliFragFuncQATrackHistos, 2);
  };
  

  AliAnalysisTaskFragmentationFunction(); 
  AliAnalysisTaskFragmentationFunction(const char *name);
  AliAnalysisTaskFragmentationFunction(const  AliAnalysisTaskFragmentationFunction &copy);
  AliAnalysisTaskFragmentationFunction& operator=(const  AliAnalysisTaskFragmentationFunction &o);
  virtual ~AliAnalysisTaskFragmentationFunction();
  
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   LocalInit() {Init();}

  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t* );
  virtual Bool_t Notify();

  virtual void   SetNonStdFile(char* c){fNonStdFile = c;} 

  virtual void   SetTrackTypeGen(Int_t i){fTrackTypeGen = i;}
  virtual void   SetJetTypeGen(Int_t i){fJetTypeGen = i;}
  virtual void   SetJetTypeRecEff(Int_t i){fJetTypeRecEff = i;}

  virtual void   SetBranchRecBackClusters(const char* c){fBranchRecBckgClusters = c;}
  virtual void   SetBranchGenJets(const char* c){fBranchGenJets = c;}
  virtual void   SetBranchRecJets(const char* c){fBranchRecJets = c;}
  virtual void   SetBranchEmbeddedJets(const char* c){fBranchEmbeddedJets = c;}

  virtual void   SetTrackCuts(Float_t trackPt = 0.15, Float_t trackEtaMin = -0.9, Float_t trackEtaMax = 0.9, 
			      Float_t trackPhiMin = 0., Float_t trackPhiMax = 2*TMath::Pi())
  {fTrackPtCut = trackPt; fTrackEtaMin = trackEtaMin; fTrackEtaMax = trackEtaMax; 
    fTrackPhiMin = trackPhiMin; fTrackPhiMax = trackPhiMax;}

  virtual void   UseExtraTracks()        { fUseExtraTracks =  1;}
  virtual void   UseExtraonlyTracks()    { fUseExtraTracks = -1;}

  virtual void   UseExtraTracksBgr()     { fUseExtraTracksBgr =  1;}
  virtual void   UseExtraonlyTracksBgr() { fUseExtraTracksBgr = -1;}

  virtual void   SetCutFractionPtEmbedded(Float_t cut = 0) { fCutFractionPtEmbedded = cut; }
  virtual void   SetUseEmbeddedJetAxis(Bool_t b = kTRUE)   { fUseEmbeddedJetAxis = b; }
  virtual void   SetUseEmbeddedJetPt(Bool_t  b = kTRUE)    { fUseEmbeddedJetPt   = b; }

  virtual void   UseAODInputJets(Bool_t b) {fUseAODInputJets = b;}  
  virtual void   SetFilterMask(UInt_t i) {fFilterMask = i;}
  virtual void   UsePhysicsSelection(Bool_t b) {fUsePhysicsSelection = b;}
  virtual void   SetEventSelectionMask(UInt_t i){fEvtSelectionMask = i;}
  virtual void   SetEventClass(Int_t i){fEventClass = i;}
  virtual void   SetMaxVertexZ(Float_t z){fMaxVertexZ = z;}
  virtual void   SetJetCuts(Float_t jetPt = 5., Float_t jetEtaMin = -0.5, Float_t jetEtaMax = 0.5, 
			    Float_t jetPhiMin = 0., Float_t jetPhiMax = 2*TMath::Pi())
  {fJetPtCut = jetPt; fJetEtaMin = jetEtaMin; fJetEtaMax = jetEtaMax; 
    fJetPhiMin = jetPhiMin; fJetPhiMax = jetPhiMax;}

  virtual void   SetFFRadius(Float_t r = 0.4) { fFFRadius = r; }
  virtual void   SetFFMinLTrackPt(Float_t pt = -1) { fFFMinLTrackPt = pt; }
  virtual void   SetFFMaxTrackPt(Float_t pt = -1) { fFFMaxTrackPt = pt; }
  virtual void   SetFFMinNTracks(Int_t nTracks = 0) { fFFMinnTracks = nTracks; }
  virtual void   SetFFBckgRadius(Float_t r = 0.7) { fFFBckgRadius = r; }
  virtual void   SetBckgMode(Bool_t bg = 1) { fBckgMode = bg; }
  virtual void   SetBckgType(Int_t bg0 = 0, Int_t bg1 = 0,Int_t bg2 = 0, Int_t bg3 = 0, Int_t bg4 = 0) 
  { fBckgType[0] = bg0; fBckgType[1] = bg1; fBckgType[2] = bg2; fBckgType[3] = bg3; fBckgType[4] = bg4; }
  virtual void   SetQAMode(Int_t qa = 3)      {fQAMode = qa;}
  virtual void   SetFFMode(Int_t ff = 1)      {fFFMode = ff;}
  virtual void   SetEffMode(Int_t eff = 1)    {fEffMode = eff;}
  virtual void   SetJSMode(Int_t js = 1)      {fJSMode = js;}

  static  void   SetProperties(TH1* h,const char* x, const char* y);
  static  void   SetProperties(TH1* h,const char* x, const char* y,const char* z);
  static  void   SetProperties(THnSparse* h,const Int_t dim, const char** labels);

  void   SetHighPtThreshold(Float_t pt = 5.) { fQATrackHighPtThreshold = pt; }

  void   SetFFHistoBins(Int_t nJetPt = 245, Float_t jetPtMin = 5, Float_t jetPtMax = 250, 
			Int_t nPt = 200, Float_t ptMin = 0., Float_t ptMax = 200., 
			Int_t nXi = 70, Float_t xiMin = 0., Float_t xiMax = 7.,
			Int_t nZ = 22,  Float_t zMin = 0.,  Float_t zMax = 1.1)
  { fFFNBinsJetPt = nJetPt; fFFJetPtMin = jetPtMin; fFFJetPtMax = jetPtMax; 
    fFFNBinsPt = nPt; fFFPtMin = ptMin; fFFPtMax = ptMax;
    fFFNBinsXi = nXi; fFFXiMin = xiMin; fFFXiMax = xiMax;
    fFFNBinsZ  = nZ;  fFFZMin  = zMin;  fFFZMax  = zMax; }

  void  SetQAJetHistoBins(Int_t nPt = 300, Float_t ptMin = 0., Float_t ptMax = 300.,
			  Int_t nEta = 20, Float_t etaMin = -1.0, Float_t etaMax = 1.0,
			  Int_t nPhi = 60, Float_t phiMin = 0., Float_t phiMax = 2*TMath::Pi())
    { fQAJetNBinsPt = nPt; fQAJetPtMin = ptMin; fQAJetPtMax = ptMax;
      fQAJetNBinsEta = nEta; fQAJetEtaMin = etaMin; fQAJetEtaMax = etaMax;
      fQAJetNBinsPhi = nPhi; fQAJetPhiMin = phiMin; fQAJetPhiMax = phiMax; }
  
  void  SetQATrackHistoBins(Int_t nPt = 200, Float_t ptMin = 0., Float_t ptMax = 200.,
			    Int_t nEta = 20, Float_t etaMin = -1.0, Float_t etaMax = 1.0,
			    Int_t nPhi = 60, Float_t phiMin = 0., Float_t phiMax = 2*TMath::Pi())
  { fQATrackNBinsPt = nPt; fQATrackPtMin = ptMin; fQATrackPtMax = ptMax;
    fQATrackNBinsEta = nEta; fQATrackEtaMin = etaMin; fQATrackEtaMax = etaMax;
    fQATrackNBinsPhi = nPhi; fQATrackPhiMin = phiMin; fQATrackPhiMax = phiMax; }
  


  Float_t  GetFFRadius() const { return fFFRadius; }
  Float_t  GetFFMinLTrackPt() const { return fFFMinLTrackPt; }
  Float_t  GetFFMaxTrackPt() const { return fFFMaxTrackPt; }
  Float_t  GetFFMinNTracks() const { return fFFMinnTracks; }
  Float_t  GetFFBckgRadius() const { return fFFBckgRadius; }
  void	   GetJetTracksTrackrefs(TList* l, const AliAODJet* j, const Double_t minPtL, const Double_t maxPt, Bool_t& isBadPt);
  void	   GetJetTracksPointing(TList* in, TList* out, const AliAODJet* j, const Double_t r, Double_t& sumPt, const Double_t minPtL, const Double_t maxPt, Bool_t& isBadPt);  
  void     GetTracksOutOfNJets(Int_t nCases, TList* in, TList* out, TList* jets, Double_t& pt);
  void     GetTracksOutOfNJetsStat(Int_t nCases, TList* in, TList* out, TList* jets, Double_t& pt, Double_t &normFactor);
  void     GetTracksTiltedwrpJetAxis(Float_t alpha, TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius, Double_t& sumPt);
  void     GetTracksTiltedwrpJetAxisWindow(Float_t alpha, TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius, Double_t& sumPt, Double_t &normFactor);

  void     AssociateGenRec(TList* tracksAODMCCharged,TList* tracksRec, TArrayI& indexAODTr,TArrayI& indexMCTr,TArrayS& isRefGen,TH2F* fh2PtRecVsGen);

  void     FillSingleTrackHistosRecGen(AliFragFuncQATrackHistos* trackQAGen, AliFragFuncQATrackHistos* trackQARec, TList* tracksGen, 
				       const TArrayI& indexAODTr, const TArrayS& isRefGen, const Bool_t scaleStrangeness = kFALSE);


  void     FillJetTrackHistosRec(AliFragFuncHistos* histRec,  AliAODJet* jet, 
				 TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, const TArrayI& indexAODTr,
				 const TArrayS& isRefGen, TList* jetTrackListTR = 0, const Bool_t scaleStrangeness = kFALSE,
				 Bool_t fillJS = kFALSE, TProfile* hProNtracksLeadingJet = 0, TProfile** hProDelRPtSum = 0, TProfile* hProDelR80pcPt = 0);


  Float_t  CalcJetArea(const Float_t etaJet, const Float_t rc) const;
  void     GetClusterTracksOutOf1Jet(AliAODJet* jet, TList* outputlist, Double_t &normFactor);
  void     GetClusterTracksMedian(TList* outputlist, Double_t &normFactor);

  void     FillBckgHistos(Int_t type, TList* inputtracklist, TList* inputjetlist, AliAODJet* jet, 
			  AliFragFuncHistos* ffbckghistocuts,AliFragFuncQATrackHistos* qabckghistos,TH1F* fh1Mult = 0); 
 
  Double_t GetMCStrangenessFactor(const Double_t pt);
  void FillJetShape(AliAODJet* jet, TList* list,  TProfile* hProNtracksLeadingJet, TProfile** hProDelRPtSum, TProfile* hProDelR80pcPt=0, Double_t dPhiUE=0, Double_t normUE = 0, Bool_t scaleStrangeness = kFALSE);


  // Consts
  enum {kTrackUndef=0, kTrackAOD, kTrackAODQualityCuts, kTrackAODCuts, 
	kTrackAODExtra, kTrackAODExtraonly, kTrackAODExtraCuts, kTrackAODExtraonlyCuts, 
	kTrackKineAll, kTrackKineCharged, kTrackKineChargedAcceptance, 
	kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance, kTrackAODMCChargedSecS, kTrackAODMCChargedSecNS, kTrackAOCMCChargedPrimAcceptance};
  enum {kJetsUndef=0, kJetsRec, kJetsRecAcceptance, kJetsGen, kJetsGenAcceptance, kJetsKine, kJetsKineAcceptance,kJetsEmbedded};
  enum {kBckgNone=0, kBckgPerp, kBckgOutLJ, kBckgOut2J, kBckgClusters, kBckgClustersOutLeading, kBckgOut3J, kBckgOutAJ, kBckgOutLJStat, 
	kBckgOut2JStat, kBckgOut3JStat, kBckgOutAJStat,  kBckgASide, kBckgASideWindow, kBckgPerpWindow, kBckgPerp2, kBckgPerp2Area};

 
 protected:
  
  Int_t   GetListOfTracks(TList* list, Int_t type);
  Int_t	  GetListOfJets(TList* list, Int_t type);
  Int_t   GetListOfBckgJets(TList *list, Int_t type);

  AliESDEvent* fESD;      // ESD event
  AliAODEvent* fAOD;      // AOD event
  AliAODEvent* fAODJets;  // AOD event with jet branch (case we have AOD both in input and output)
  AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD
  TString       fNonStdFile; // name of delta aod file to catch the extension
 
 
  TString fBranchRecJets;         // branch name for reconstructed jets
  TString fBranchRecBckgClusters; // branch name for reconstructed background clusters 
  TString fBranchGenJets;         // branch name for generated jets
  TString fBranchEmbeddedJets;    // branch name for embedded jets

  Int_t   fTrackTypeGen;        // type of generated tracks
  Int_t   fJetTypeGen;          // type of generated jets

  Int_t   fJetTypeRecEff;       // type of jets used for filling reconstruction efficiency histos

  Bool_t  fUseAODInputJets;     // take jets from in/output - only relevant if AOD event both in input AND output and we want to use output
  UInt_t  fFilterMask;	        // filter bit for selected tracks
  Bool_t  fUsePhysicsSelection; // switch for event selection
  UInt_t  fEvtSelectionMask;    // trigger class selection
  Int_t   fEventClass;          // centrality class selection
  Float_t fMaxVertexZ;          // maximum abs(z) position of primiary vertex [cm]

  // track cuts
  Float_t fTrackPtCut;    // track transverse momentum cut
  Float_t fTrackEtaMin;   // track eta cut
  Float_t fTrackEtaMax;   // track eta cut
  Float_t fTrackPhiMin;   // track phi cut
  Float_t fTrackPhiMax;   // track phi cut
  
  Int_t   fUseExtraTracks;         // +/- 1: embedded extra/extra only tracks, default: 0 (ignore extra tracks)
  Int_t   fUseExtraTracksBgr;      // +/- 1: background: use embedded extra/extra only tracks, default: 0 (ignore extra tracks)
  Float_t fCutFractionPtEmbedded;  // cut on ratio of embedded pt found in jet
  Bool_t  fUseEmbeddedJetAxis;     // use axis of embedded jet for FF
  Bool_t  fUseEmbeddedJetPt;       // use axis of embedded jet for FF

  // jet cuts
  Float_t fJetPtCut;      // jet transverse momentum cut
  Float_t fJetEtaMin;     // jet eta cut
  Float_t fJetEtaMax;     // jet eta cut
  Float_t fJetPhiMin;     // jet phi cut
  Float_t fJetPhiMax;     // jet phi cut

  Float_t fFFRadius;        // if radius > 0 construct FF from tracks within cone around jet axis, otherwise use trackRefs  
  Float_t fFFMinLTrackPt;   // reject jets with leading track with pt smaller than this value
  Float_t fFFMaxTrackPt;    // reject jets containing any track with pt larger than this value
  Int_t   fFFMinnTracks;    // reject jets with less tracks than this value
  Float_t fFFBckgRadius;    // compute background outside cone of this radius around jet axes
  Bool_t  fBckgMode;        // Set background subtraction mode
  Int_t   fBckgType[5];     // Set background subtraction mode
  Int_t   fQAMode;          // QA mode: 0x00=0 none, 0x01=1 track qa, 0x10=2 track qa, 0x11=3 both
  Int_t   fFFMode;          // fragmentation function mode
  Int_t   fEffMode;         // efficiency mode
  Int_t   fJSMode;          // jet shape mode

  Float_t fAvgTrials;       // average number of trials per event
  
  TList* fTracksRecCuts;           //! reconstructed tracks after cuts
  TList* fTracksGen;               //! generated tracks 
  TList* fTracksAODMCCharged;      //! AOD MC tracks 
  TList* fTracksAODMCChargedSecNS; //! AOD MC tracks - secondaries (non-strangeness) 
  TList* fTracksAODMCChargedSecS;  //! AOD MC tracks - secondaries (from strangeness)
  TList* fTracksRecQualityCuts;    //! reconstructed tracks after quality cuts, no acceptance/pt cut

  
  TList* fJetsRec;        //! jets from reconstructed tracks
  TList* fJetsRecCuts;    //! jets from reonstructed tracks after jet cuts 
  TList* fJetsGen;        //! jets from generated tracks
  TList* fJetsRecEff;     //! jets used for reconstruction efficiency histos 
  TList* fJetsEmbedded;   //! jets used for embedding

  TList* fBckgJetsRec;      //! jets from reconstructed tracks
  TList* fBckgJetsRecCuts;  //! jets from reonstructed tracks after jet cuts
  TList* fBckgJetsGen;      //! jets from generated tracks
 
  
  AliFragFuncQATrackHistos* fQATrackHistosRecCuts;  //! track QA: reconstructed tracks after cuts
  AliFragFuncQATrackHistos* fQATrackHistosGen;      //! track QA: generated tracks
  
  AliFragFuncQAJetHistos*  fQAJetHistosRec;             //! jet QA: jets from reconstructed tracks
  AliFragFuncQAJetHistos*  fQAJetHistosRecCuts;         //! jet QA: jets from reconstructed tracks after jet cuts 
  AliFragFuncQAJetHistos*  fQAJetHistosRecCutsLeading;  //! jet QA: leading jet from reconstructed tracks after jet cuts 
  AliFragFuncQAJetHistos*  fQAJetHistosGen;             //! jet QA: jets from generated tracks  
  AliFragFuncQAJetHistos*  fQAJetHistosGenLeading;      //! jet QA: leading jet from generated tracks  
  AliFragFuncQAJetHistos*  fQAJetHistosRecEffLeading;   //! jet QA: leading jet used for reconstruction efficiency histos  
  

  AliFragFuncHistos*  fFFHistosRecCuts;         //! FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFHistosGen;             //! FF generated tracks after cuts 

  Float_t  fQATrackHighPtThreshold;       // track QA high transverse momentum threshold
  
  // histogram bins  

  Int_t   fFFNBinsJetPt;    // FF histos bins
  Float_t fFFJetPtMin;      // FF histos limits
  Float_t fFFJetPtMax;      // FF histos limits

  Int_t   fFFNBinsPt;       // FF histos bins
  Float_t fFFPtMin;         // FF histos limits
  Float_t fFFPtMax;         // FF histos limits

  Int_t   fFFNBinsXi;       // FF histos bins
  Float_t fFFXiMin;         // FF histos limits
  Float_t fFFXiMax;         // FF histos limits

  Int_t   fFFNBinsZ;        // FF histos bins
  Float_t fFFZMin;          // FF histos limits
  Float_t fFFZMax;          // FF histos limits

  Int_t   fQAJetNBinsPt;    // jet QA histos bins
  Float_t fQAJetPtMin;      // jet QA histos limits
  Float_t fQAJetPtMax;      // jet QA histos limits
  
  Int_t   fQAJetNBinsEta;   // jet QA histos bins
  Float_t fQAJetEtaMin;     // jet QA histos limits
  Float_t fQAJetEtaMax;     // jet QA histos limits
  
  Int_t   fQAJetNBinsPhi;   // jet QA histos bins
  Float_t fQAJetPhiMin;     // jet QA histos limits
  Float_t fQAJetPhiMax;     // jet QA histos limits
  
  Int_t   fQATrackNBinsPt;  // track QA histos bins
  Float_t fQATrackPtMin;    // track QA histos limits
  Float_t fQATrackPtMax;    // track QA histos limits
  
  Int_t   fQATrackNBinsEta; // track QA histos bins
  Float_t fQATrackEtaMin;   // track QA histos limits
  Float_t fQATrackEtaMax;   // track QA histos limits
  
  Int_t   fQATrackNBinsPhi; // track QA histos bins
  Float_t fQATrackPhiMin;   // track QA histos limits
  Float_t fQATrackPhiMax;   // track QA histos limits
  
  // Histograms
  TList	*fCommonHistList;         // List of common histos
  
  TH1F  *fh1EvtSelection;         //! event cuts 
  TH1F	*fh1VertexNContributors;  //! NContributors to prim vertex
  TH1F	*fh1VertexZ;              //! prim vertex z distribution
  TH1F	*fh1EvtMult;              //! number of reconstructed tracks after cuts 
  TH1F	*fh1EvtCent;              //! centrality percentile 

  TProfile* fh1Xsec;              //! pythia cross section and trials
  TH1F*     fh1Trials;            //! sum of trials
  TH1F*     fh1PtHard;            //! pt hard of the event
  TH1F*     fh1PtHardTrials;      //! pt hard of the event

  TH1F  *fh1nRecJetsCuts;         //! number of jets from reconstructed tracks per event 
  TH1F  *fh1nGenJets;             //! number of jets from generated tracks per event
  TH1F  *fh1nRecEffJets;          //! number of jets for reconstruction eff per event
  TH1F  *fh1nEmbeddedJets;        //! number of embedded jets per event

  TH1F  *fh1nRecBckgJetsCuts;     //! number of jets from reconstructed tracks per event
  TH1F  *fh1nGenBckgJets;         //! number of jets from generated tracks per event
  TH2F  *fh2PtRecVsGenPrim;       //! association rec/gen MC: rec vs gen pt, primaries 
  TH2F  *fh2PtRecVsGenSec;        //! association rec/gen MC: rec vs gen pt, secondaries 

  // tracking efficiency / secondaries
  
  AliFragFuncQATrackHistos* fQATrackHistosRecEffGen;      //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffRec;      //! tracking efficiency: reconstructed primaries
  AliFragFuncQATrackHistos* fQATrackHistosSecRecNS;       //! reconstructed secondaries (non-strangeness)
  AliFragFuncQATrackHistos* fQATrackHistosSecRecS;        //! reconstructed secondaries (strange mothers)
  AliFragFuncQATrackHistos* fQATrackHistosSecRecSsc;      //! reconstructed secondaries (strange mothers) - scale factor

  AliFragFuncHistos*  fFFHistosRecEffRec;                 //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosSecRecNS;                  //! secondary contamination: FF reconstructed secondaries (non-strangeness)
  AliFragFuncHistos*  fFFHistosSecRecS;                   //! secondary contamination: FF reconstructed secondaries (strange mothers)
  AliFragFuncHistos*  fFFHistosSecRecSsc;                 //! secondary contamination: FF reconstructed secondaries (strange mothers) - scale factor


  // Background
  TH1F  *fh1BckgMult0; //! background multiplicity
  TH1F  *fh1BckgMult1; //! background multiplicity
  TH1F  *fh1BckgMult2; //! background multiplicity
  TH1F  *fh1BckgMult3; //! background multiplicity
  TH1F  *fh1BckgMult4; //! background multiplicity

  // embedding
  TH1F* fh1FractionPtEmbedded;         //! ratio embedded pt in rec jet to embedded jet pt 
  TH1F* fh1IndexEmbedded;              //! index embedded jet matching to leading rec jet 
  TH2F*	fh2DeltaPtVsJetPtEmbedded;     //! delta pt rec - embedded jet
  TH2F*	fh2DeltaPtVsRecJetPtEmbedded;  //! delta pt rec - embedded jet
  TH1F* fh1DeltaREmbedded;             //! delta R  rec - embedded jet


  AliFragFuncQATrackHistos* fQABckgHisto0RecCuts;  //! track QA: reconstructed tracks after cuts
  AliFragFuncQATrackHistos* fQABckgHisto0Gen;      //! track QA: generated tracks
  AliFragFuncQATrackHistos* fQABckgHisto1RecCuts;  //! track QA: reconstructed tracks after cuts
  AliFragFuncQATrackHistos* fQABckgHisto1Gen;      //! track QA: generated tracks
  AliFragFuncQATrackHistos* fQABckgHisto2RecCuts;  //! track QA: reconstructed tracks after cuts
  AliFragFuncQATrackHistos* fQABckgHisto2Gen;      //! track QA: generated tracks
  AliFragFuncQATrackHistos* fQABckgHisto3RecCuts;  //! track QA: reconstructed tracks after cuts
  AliFragFuncQATrackHistos* fQABckgHisto3Gen;      //! track QA: generated tracks
  AliFragFuncQATrackHistos* fQABckgHisto4RecCuts;  //! track QA: reconstructed tracks after cuts
  AliFragFuncQATrackHistos* fQABckgHisto4Gen;      //! track QA: generated tracks
  
  AliFragFuncHistos*  fFFBckgHisto0RecCuts;       //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto0Gen;           //! Bckg (outside leading jet or 2 jets or more) FF generated tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto1RecCuts;       //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto1Gen;           //! Bckg (outside leading jet or 2 jets or more) FF generated tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto2RecCuts;       //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto2Gen;           //! Bckg (outside leading jet or 2 jets or more) FF generated tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto3RecCuts;       //! Bckg (outside leading jet or 3 jets or more) FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto3Gen;           //! Bckg (outside leading jet or 3 jets or more) FF generated tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto4RecCuts;       //! Bckg (outside leading jet or 4 jets or more) FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto4Gen;           //! Bckg (outside leading jet or 4 jets or more) FF generated tracks after cuts 

  AliFragFuncHistos*  fFFBckgHisto0RecEffRec;     //! Bckg (outside leading jet or 2 jets or more) FF reconstructed primaries after cuts 
  AliFragFuncHistos*  fFFBckgHisto0SecRecNS;      //! secondary contamination: FF reconstructed secondaries (non-strangeness)
  AliFragFuncHistos*  fFFBckgHisto0SecRecS;       //! secondary contamination: FF reconstructed secondaries (strange mothers)
  AliFragFuncHistos*  fFFBckgHisto0SecRecSsc;     //! secondary contamination: FF reconstructed secondaries (strange mothers) - scale factor

  TProfile* fProNtracksLeadingJet;          //! jet shape 
  TProfile* fProDelR80pcPt;                 //! jet shape 
  TProfile* fProDelRPtSum[5];               //! jet shape 

  TProfile* fProNtracksLeadingJetGen;       //! jet shape 
  TProfile* fProDelR80pcPtGen;              //! jet shape 
  TProfile* fProDelRPtSumGen[5];            //! jet shape 

  TProfile* fProNtracksLeadingJetBgrPerp2;  //! jet shape 
  TProfile* fProDelRPtSumBgrPerp2[5];       //! jet shape 

  TProfile* fProNtracksLeadingJetRecPrim;   //! jet shape 
  TProfile* fProDelR80pcPtRecPrim;          //! jet shape 
  TProfile* fProDelRPtSumRecPrim[5];        //! jet shape 

  TProfile* fProNtracksLeadingJetRecSecNS;  //! jet shape 
  TProfile* fProDelRPtSumRecSecNS[5];       //! jet shape 

  TProfile* fProNtracksLeadingJetRecSecS;   //! jet shape 
  TProfile* fProDelRPtSumRecSecS[5];        //! jet shape 

  TProfile* fProNtracksLeadingJetRecSecSsc; //! jet shape 
  TProfile* fProDelRPtSumRecSecSsc[5];      //! jet shape 
  

  TRandom3*                   fRandom;          // TRandom3 for background estimation 

  ClassDef(AliAnalysisTaskFragmentationFunction, 12);
};

#endif
