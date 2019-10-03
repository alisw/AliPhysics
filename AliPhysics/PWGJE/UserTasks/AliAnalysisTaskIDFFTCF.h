// *************************************************************************
// * Task for Fragmentation Function Analysis in PWG4 Jet Task Force Train *
// *************************************************************************

#ifndef ALIANALYSISTASKIDFFTCFN_H
#define ALIANALYSISTASKIDFFTCFN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class AliESDEvent;
class AliAODEvent;
class AliAODJets;
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
class AliAODMCParticle;

#include "AliAnalysisTaskSE.h"
#include "TAxis.h"
#include "THnSparse.h"
#include <TTreeStream.h>
  
class AliAnalysisTaskIDFFTCF : public AliAnalysisTaskSE {

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
    virtual void FillFF(Float_t trackPt, Float_t trackEta, Float_t jetPt, 
			Bool_t incrementJetPt, Float_t norm = 0, Bool_t scaleStrangeness = kFALSE, Float_t scaleFacStrangeness = 1.);

    virtual void SetLogPt(Bool_t doLog = kTRUE) { fLogPt = doLog; }

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
    Bool_t  fLogPt;      // FF track log pt bins   

    TH2F*   fh2TrackPt;   //! FF: track transverse momentum 
    TH2F*   fh2Xi;        //! FF: xi 
    TH2F*   fh2Z;         //! FF: z  
    TH1F*   fh1JetPt;     //! jet pt 

    TH3F*   fh3TrackPtVsEta;  //! FF: track transverse momentum vs track eta 
    TH3F*   fh3TrackPVsEta;   //! FF: track momentum vs track eta 

    TString fNameFF;      // histo names prefix
    
    ClassDef(AliFragFuncHistos, 1);
  };

  //----------------------------------------
  class AliFragFuncHistosMult : public TObject
  {
    
    public:
    
    AliFragFuncHistosMult(const char* name = "FFhistos", 
			  Int_t nJetPt = 0, Float_t jetPtMin = 0, Float_t jetPtMax = 0,
			  Int_t nPt = 0, Float_t ptMin = 0, Float_t ptMax = 0,
			  Int_t nXi = 0, Float_t xiMin = 0, Float_t xiMax = 0,
			  Int_t nZ  = 0, Float_t zMin  = 0, Float_t zMax  = 0,
			  Int_t nMult = 0, Float_t multMin  = 0, Float_t multMax  = 0);
    
    AliFragFuncHistosMult(const AliFragFuncHistosMult& copy);
    AliFragFuncHistosMult& operator=(const AliFragFuncHistosMult &o);
    virtual ~AliFragFuncHistosMult();
    
    virtual void DefineHistos();
    virtual void FillFF(Float_t trackPt, Float_t trackEta, Int_t mult, Float_t jetPt,  
			Bool_t incrementJetPt, Float_t norm = 0, Bool_t scaleStrangeness = kFALSE, Float_t scaleFacStrangeness = 1.);

    virtual void SetLogPt(Bool_t doLog = kTRUE) { fLogPt = doLog; }

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
    Bool_t  fLogPt;      // FF track log pt bins   
    Int_t   fNBinsMult;  // FF histos bins
    Float_t fMultMin;    // FF histos limits
    Float_t fMultMax;    // FF histos limits

  
    TH3F*   fh3TrackPt;   //! FF: track transverse momentum 
    TH3F*   fh3Xi;        //! FF: xi 
    TH3F*   fh3Z;         //! FF: z  
    TH2F*   fh2JetPt;     //! jet pt 

    TH3F*   fh3TrackPtVsEta;  //! FF: track transverse momentum vs track eta 
    TH3F*   fh3TrackPVsEta;   //! FF: track momentum vs track eta 

    TString fNameFF;      // histo names prefix
    
    ClassDef(AliFragFuncHistosMult, 1);
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
  
  enum TPCCUTMODE{
    kPIDNone = 0, 
    kPIDN,
    kMIGeo
  };
  static Bool_t fkDump;  //=1: enable debug streamer; =0 : not.

  AliAnalysisTaskIDFFTCF(); 
  AliAnalysisTaskIDFFTCF(const char *name);
  AliAnalysisTaskIDFFTCF(const  AliAnalysisTaskIDFFTCF &copy);
  AliAnalysisTaskIDFFTCF& operator=(const  AliAnalysisTaskIDFFTCF &o);
  virtual ~AliAnalysisTaskIDFFTCF();
  
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

  virtual void   SetBranchGenJets(const char* c){fBranchGenJets = c;}
  virtual void   SetBranchRecJets(const char* c){fBranchRecJets = c;}

  virtual void   SetTrackCuts(Float_t trackPt = 0.15, Float_t trackEtaMin = -0.9, Float_t trackEtaMax = 0.9, 
			      Float_t trackPhiMin = 0., Float_t trackPhiMax = 2*TMath::Pi())
  {fTrackPtCut = trackPt; fTrackEtaMin = trackEtaMin; fTrackEtaMax = trackEtaMax; 
    fTrackPhiMin = trackPhiMin; fTrackPhiMax = trackPhiMax;}

  virtual void   UseAODInputJets(Bool_t b) {fUseAODInputJets = b;}  
  virtual void   SetFilterMask(UInt_t i) {fFilterMask = i;}
  virtual void   UsePhysicsSelection(Bool_t b) {fUsePhysicsSelection = b;}
  virtual void   SetEventSelectionMask(UInt_t i){fEvtSelectionMask = i;}
  virtual void   SetEventClass(Int_t i){fEventClass = i;}
  virtual void   SetMaxVertexZ(Float_t z){fMaxVertexZ = z;}
  virtual void   RejectPileupEvents(Bool_t b){fRejectPileup = b;}
  virtual void   UseLeadingJet(Bool_t b){fLeadingJets = b;}

  virtual void   SetJetCuts(Float_t jetPt = 5., Float_t jetEtaMin = -0.5, Float_t jetEtaMax = 0.5, 
			    Float_t jetPhiMin = 0., Float_t jetPhiMax = 2*TMath::Pi())
  {fJetPtCut = jetPt; fJetEtaMin = jetEtaMin; fJetEtaMax = jetEtaMax; 
    fJetPhiMin = jetPhiMin; fJetPhiMax = jetPhiMax;}

  virtual void   SetFFRadius(Float_t r = 0.4) { fFFRadius = r; }
  virtual void   SetFFMinLTrackPt(Float_t pt = -1) { fFFMinLTrackPt = pt; }
  virtual void   SetFFMaxTrackPt(Float_t pt = -1) { fFFMaxTrackPt = pt; }
  virtual void   SetFFMinNTracks(Int_t nTracks = 0) { fFFMinnTracks = nTracks; }
  virtual void   SetQAMode(Int_t qa = 3)      {fQAMode = qa;}
  virtual void   SetFFMode(Int_t ff = 1)      {fFFMode = ff;}
  virtual void   SetEffMode(Int_t eff = 1)    {fEffMode = eff;}
  virtual void   SetRespMode(Int_t eff = 1)   {fRespMode = eff;}

  static  void   SetProperties(TH1* h,const char* x, const char* y);
  static  void   SetProperties(TH1* h,const char* x, const char* y,const char* z);
  static  void   SetProperties(THnSparse* h, Int_t dim, const char** labels);

  void SetTPCCutMode(Int_t mode){ fTPCCutMode = mode; }
  Int_t GetTPCCutMode(){return fTPCCutMode; }

  void SetTOFCutMode(Int_t mode){ fTOFCutMode = mode; }
  Int_t GetTOFCutMode(){return fTOFCutMode; }

  void   SetHighPtThreshold(Float_t pt = 5.) { fQATrackHighPtThreshold = pt; }

  void   SetFFHistoBins(Int_t nJetPt = 245, Float_t jetPtMin = 5, Float_t jetPtMax = 250, 
			Int_t nPt = 200, Float_t ptMin = 0., Float_t ptMax = 200., 
			Int_t nXi = 70, Float_t xiMin = 0., Float_t xiMax = 7.,
			Int_t nZ = 22,  Float_t zMin = 0.,  Float_t zMax = 1.1,
			Int_t nMult = 20, Float_t multMin = 0., Float_t multMax = 20)
  { fFFNBinsJetPt = nJetPt; fFFJetPtMin = jetPtMin; fFFJetPtMax = jetPtMax; 
    fFFNBinsPt = nPt; fFFPtMin = ptMin; fFFPtMax = ptMax;
    fFFNBinsXi = nXi; fFFXiMin = xiMin; fFFXiMax = xiMax;
    fFFNBinsZ = nZ; fFFZMin = zMin; fFFZMax = zMax;
    fFFNBinsMult  = nMult;  fFFMultMin  = multMin;  fFFMultMax  = multMax; }

  void  SetFFLogPt(Bool_t doLog = kTRUE) { fFFLogPt = doLog; }

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

  void	   GetJetTracksTrackrefs(TList* l, const AliAODJet* j, Double_t minPtL, Double_t maxPt, Bool_t& isBadPt);
  void	   GetJetTracksPointing(TList* in, TList* out, const AliAODJet* j, Double_t r, Double_t& sumPt, Double_t minPtL, Double_t maxPt, Bool_t& isBadPt);  

  void     AssociateGenRec(TList* tracksAODMCCharged,TList* tracksRec, TArrayI& indexAODTr,TArrayI& indexMCTr,TArrayS& isRefGen,TH2F* fh2PtRecVsGen);

  void     FillSingleTrackHistosRecGen(AliFragFuncQATrackHistos* trackQAGen, AliFragFuncQATrackHistos* trackQARec, TList* tracksGen, 
				       const TArrayI& indexAODTr, const TArrayS& isRefGen, Int_t pdg = 0, 
				       Bool_t scaleGFL = kFALSE, Bool_t scaleStrangeness = kFALSE);


  void     FillJetTrackHistosRec(AliFragFuncHistos* histRec,  AliAODJet* jet, 
				 TList* jetTrackList, const TList* tracksGen, const TList* tracksRec, const TArrayI& indexAODTr,
				 const TArrayS& isRefGen, TList* jetTrackListTR = 0, Int_t pdg = 0, 
				 Bool_t scaleGFL = kFALSE, Bool_t scaleStrangeness = kFALSE);


  Float_t  CalcJetArea(Float_t etaJet, Float_t rc) const;
 
  void     BookQAHistos(TList* list = 0, AliFragFuncQATrackHistos** rec = 0, TString strTitRec = "", AliFragFuncQATrackHistos** gen = 0, TString strTitGen = "",
			AliFragFuncQATrackHistos** sec = 0, TString strTitSec = "");

  void     BookFFHistos(TList* list, AliFragFuncHistos** rec = 0, TString strTitRec = "", AliFragFuncHistos** gen = 0, TString strTitGen = "",
			AliFragFuncHistos** sec = 0, TString strTitSec = "");

  void     BookFFHistosM(TList* list, AliFragFuncHistosMult** rec = 0, TString strTitRec = "", 
			 AliFragFuncHistosMult** gen = 0, TString strTitGen = "",
			 AliFragFuncHistosMult** sec = 0, TString strTitSec = "");

  Double_t  TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc);
  Double_t  TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc);
  Double_t  GetMCStrangenessFactorCMS(AliAODMCParticle* daughter);
  void      FillResponse(AliAODJet* genJet, AliAODJet* recJet);
  void      GetTracksTiltedwrpJetAxis(Float_t alpha, TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius, Double_t& sumPt);

  // Consts
  enum {kTrackUndef=0, kTrackAOD, kTrackAODQualityCuts, kTrackAODCuts,  
	kTrackKineAll, kTrackKineCharged, kTrackKineChargedAcceptance, 
	kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance, kTrackAODMCChargedSec, kTrackAOCMCChargedPrimAcceptance};
  enum {kJetsUndef=0, kJetsRec, kJetsRecAcceptance, kJetsGen, kJetsGenAcceptance, kJetsKine, kJetsKineAcceptance};

 
 protected:
  
  Int_t   GetListOfTracks(TList* list, Int_t type);
  Int_t	  GetListOfJets(TList* list, Int_t type);

  AliESDEvent* fESD;      // ESD event
  AliAODEvent* fAOD;      // AOD event
  AliAODEvent* fAODJets;  // AOD event with jet branch (case we have AOD both in input and output)
  AliAODExtension  *fAODExtension; //! where we take the jets from can be input or output AOD
  TString       fNonStdFile; // name of delta aod file to catch the extension
 
 
  TString fBranchRecJets;         // branch name for reconstructed jets
  TString fBranchGenJets;         // branch name for generated jets

  Int_t   fTrackTypeGen;        // type of generated tracks
  Int_t   fJetTypeGen;          // type of generated jets

  Int_t   fJetTypeRecEff;       // type of jets used for filling reconstruction efficiency histos

  Bool_t  fUseAODInputJets;     // take jets from in/output - only relevant if AOD event both in input AND output and we want to use output
  UInt_t  fFilterMask;	        // filter bit for selected tracks
  Bool_t  fUsePhysicsSelection; // switch for event selection
  UInt_t  fEvtSelectionMask;    // trigger class selection
  Int_t   fEventClass;          // centrality class selection
  Float_t fMaxVertexZ;          // maximum abs(z) position of primiary vertex [cm]
  Bool_t  fRejectPileup;        // SPD pileup rejection
  Bool_t  fLeadingJets;         // leading/all jets


  Int_t fTPCCutMode;      //mode for cutting TPC for good dE/dx
  Int_t fTOFCutMode;      //mode for cutting TOF
  TTreeStream * fStream; //debug streamer
  TTree * fTree;         //tree of streamer

  // track cuts
  Float_t fTrackPtCut;    // track transverse momentum cut
  Float_t fTrackEtaMin;   // track eta cut
  Float_t fTrackEtaMax;   // track eta cut
  Float_t fTrackPhiMin;   // track phi cut
  Float_t fTrackPhiMax;   // track phi cut
  

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
  Int_t   fQAMode;          // QA mode: 0x00=0 none, 0x01=1 track qa, 0x10=2 track qa, 0x11=3 both
  Int_t   fFFMode;          // fragmentation function mode
  Int_t   fEffMode;         // efficiency mode
  Int_t   fRespMode;        // response mode

  Float_t fAvgTrials;       // average number of trials per event
  
  TList* fTracksRecCuts;           //! reconstructed tracks after cuts
  TList* fTracksGen;               //! generated tracks 
  TList* fTracksAODMCCharged;      //! AOD MC tracks 
  TList* fTracksAODMCChargedSec;   //! AOD MC tracks - secondaries 
  TList* fTracksRecQualityCuts;    //! reconstructed tracks after quality cuts, no acceptance/pt cut

  TList* fJetsRec;        //! jets from reconstructed tracks
  TList* fJetsRecCuts;    //! jets from reonstructed tracks after jet cuts 
  TList* fJetsGen;        //! jets from generated tracks
  TList* fJetsRecEff;     //! jets used for reconstruction efficiency histos 

   
  AliFragFuncQATrackHistos* fQATrackHistosRecCuts;  //! track QA: reconstructed tracks after cuts
  AliFragFuncQATrackHistos* fQATrackHistosGen;      //! track QA: generated tracks
  
  AliFragFuncQAJetHistos*  fQAJetHistosRec;             //! jet QA: jets from reconstructed tracks
  AliFragFuncQAJetHistos*  fQAJetHistosRecCuts;         //! jet QA: jets from reconstructed tracks after jet cuts 
  AliFragFuncQAJetHistos*  fQAJetHistosRecCutsLeading;  //! jet QA: leading jet from reconstructed tracks after jet cuts 
  AliFragFuncQAJetHistos*  fQAJetHistosGen;             //! jet QA: jets from generated tracks  
  AliFragFuncQAJetHistos*  fQAJetHistosGenLeading;      //! jet QA: leading jet from generated tracks  
  AliFragFuncQAJetHistos*  fQAJetHistosRecEffLeading;   //! jet QA: leading jet used for reconstruction efficiency histos  
  

  AliFragFuncHistos*  fFFHistosRecCutsInc;       //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosRecCutsIncPi;     //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosRecCutsIncPro;    //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosRecCutsIncK;      //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosRecCutsIncEl;     //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosRecCutsIncMu;     //! inclusive FF (all jets) 

  AliFragFuncHistosMult*  fFFHistosMRecCutsInc;       //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMRecCutsIncPi;     //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMRecCutsIncPro;    //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMRecCutsIncK;      //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMRecCutsIncEl;     //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMRecCutsIncMu;     //! inclusive FF (all jets) 

  AliFragFuncHistos*  fFFHistosRecLeadingTrack; //! FF reconstructed tracks after cuts: leading track pt / jet pt (all jets)

  AliFragFuncHistos*  fFFHistosGenInc;       //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosGenIncPi;     //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosGenIncPro;    //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosGenIncK;      //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosGenIncEl;     //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosGenIncMu;     //! inclusive FF (all jets) 
  AliFragFuncHistos*  fFFHistosGenLeadingTrack; //! FF reconstructed tracks after cuts: leading track pt / jet pt (all jets)

  AliFragFuncHistosMult*  fFFHistosMRecIncMatch;       //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMGenIncMatch;       //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMGenIncMatchPi;     //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMGenIncMatchPro;    //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMGenIncMatchK;      //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMGenIncMatchEl;       //! inclusive FF (all jets) 
  AliFragFuncHistosMult*  fFFHistosMGenIncMatchMu;     //! inclusive FF (all jets) 


  Float_t  fQATrackHighPtThreshold;       // track QA high transverse momentum threshold
  
  THnSparseD * fTHnIDFF;                //! tracks in jets
  THnSparseD * fTHnIncl;                //! inclusive tracks

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

  Int_t   fFFNBinsMult;     // FF histos bins
  Float_t fFFMultMin;       // FF histos limits
  Float_t fFFMultMax;       // FF histos limits

  Bool_t fFFLogPt;          // fill FF histos with track log pt

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

  TH2F  *fh2PtRecVsGenPrim;       //! association rec/gen MC: rec vs gen pt, primaries 
  TH2F  *fh2PtRecVsGenSec;        //! association rec/gen MC: rec vs gen pt, secondaries 

  TH1F* hJetSpecIncRec;           //! jet spec
  TH1F* hJetSpecIncRecUEsub;      //! jet spec
  TH1F* hJetSpecIncGen;           //! jet spec
  TH1F* hJetSpecIncGenUEsub;      //! jet spec

  TH2F* h2UERec;                  //! UE
  TH2F* h2UEGen;                  //! UE

  // tracking efficiency / secondaries
  
  AliFragFuncQATrackHistos* fQATrackHistosRecEffGen;      //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffRec;      //! tracking efficiency: reconstructed primaries
  AliFragFuncQATrackHistos* fQATrackHistosSecRec;         //! reconstructed secondaries
  AliFragFuncQATrackHistos* fQATrackHistosSecRecSSc;      //! reconstructed secondaries

  AliFragFuncQATrackHistos* fQATrackHistosRecEffGenPi;     //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffGenPro;    //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffGenK;      //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffGenEl;     //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffGenMu;     //! tracking efficiency: generated primaries 

  AliFragFuncQATrackHistos* fQATrackHistosRecEffRecPi;       //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffRecPro;      //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffRecK;        //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffRecEl;       //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffRecMu;       //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffRecProGFL;   //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffRecKGFL;     //! tracking efficiency: generated primaries 

  AliFragFuncQATrackHistos* fQATrackHistosSecRecPi;       //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecPro;      //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecK;        //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecEl;       //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecMu;       //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecProGFL;   //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecKGFL;     //! tracking efficiency: generated primaries 

  AliFragFuncQATrackHistos* fQATrackHistosSecRecPiSSc;       //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecProSSc;      //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecKSSc;        //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecElSSc;       //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecMuSSc;       //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecProGFLSSc;   //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosSecRecKGFLSSc;     //! tracking efficiency: generated primaries 

  AliFragFuncHistos*  fFFHistosRecEffRec;                 //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosSecRec;                    //! secondary contamination: FF reconstructed secondaries 
  AliFragFuncHistos*  fFFHistosSecRecSSc;                 //! secondary contamination: FF reconstructed secondaries 

  AliFragFuncHistos*  fFFHistosRecEffRecPi;               //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosRecEffRecPro;              //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosRecEffRecK;                //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosRecEffRecEl;               //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosRecEffRecMu;               //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosRecEffRecProGFL;           //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosRecEffRecKGFL;             //! tracking efficiency: FF reconstructed primaries

  AliFragFuncHistos*  fFFHistosSecRecPi;                  //! secondary contamination: FF reconstructed secondaries 
  AliFragFuncHistos*  fFFHistosSecRecPro;                 //! secondary contamination: FF reconstructed secondaries 
  AliFragFuncHistos*  fFFHistosSecRecK;                   //! secondary contamination: FF reconstructed secondaries 
  AliFragFuncHistos*  fFFHistosSecRecEl;                  //! secondary contamination: FF reconstructed secondaries 
  AliFragFuncHistos*  fFFHistosSecRecMu;                  //! secondary contamination: FF reconstructed secondaries 
  AliFragFuncHistos*  fFFHistosSecRecProGFL;              //! secondary contamination: FF reconstructed secondaries 
  AliFragFuncHistos*  fFFHistosSecRecKGFL;                //! secondary contamination: FF reconstructed secondaries 

  AliFragFuncHistos*  fFFHistosSecRecPiSSc;            //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosSecRecProSSc;           //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosSecRecKSSc;             //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosSecRecElSSc;            //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosSecRecMuSSc;            //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosSecRecProGFLSSc;        //! tracking efficiency: FF reconstructed primaries
  AliFragFuncHistos*  fFFHistosSecRecKGFLSSc;          //! tracking efficiency: FF reconstructed primaries

  THnSparse* fhnResponseJetPt;
  THnSparse* fhnRespJetPtHistG;
  THnSparse* fhnRespJetPtHistR;
  THnSparse* fhnRespJetPtHistM;
  
  THnSparse* fhnResponseZ;
  THnSparse* fhnRespZHistG; 
  THnSparse* fhnRespZHistR;  
  THnSparse* fhnRespZHistM;  
  THnSparse* fhnRespZHistMPrim;  

  THnSparse* fhnResponsePt;
  THnSparse* fhnRespPtHistG; 
  THnSparse* fhnRespPtHistR;  
  THnSparse* fhnRespPtHistM;  
  THnSparse* fhnRespPtHistMPrim;  

  
  TRandom3*                   fRandom;          // TRandom3 for background estimation 

  TH2F  *fh2EtaPhiUnm;        //!  
  TH1F  *fh1AreaUnm;          //!  
  TH1F  *fh1AreaM;            //!  

  ClassDef(AliAnalysisTaskIDFFTCF, 1);
};

#endif
