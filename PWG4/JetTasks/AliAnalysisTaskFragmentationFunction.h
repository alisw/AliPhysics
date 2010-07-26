/*************************************************************************
 * Task for Fragmentation Function Analysis in PWG4 Jet Task Force Train *
 *************************************************************************/

#ifndef ALIANALYSISTASKFRAGMENTATIONFUNCTION_H
#define ALIANALYSISTASKFRAGMENTATIONFUNCTION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class AliESDEvent;
class AliAODEvent;
class TList;
class TH1F;
class TH2F;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFragmentationFunction : public AliAnalysisTaskSE {

 public:

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
    virtual void FillFF(Float_t trackPt, Float_t jetPt,Bool_t incrementJetPt);
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
  
    TH2F*  fh2TrackPt;   //! FF: track transverse momentum 
    TH2F*  fh2Xi;        //! FF: xi 
    TH2F*  fh2Z;         //! FF: z  
    TH1F*  fh1JetPt;     //! jet pt 

    TString fName;       // histo names prefix
    
    ClassDef(AliFragFuncHistos, 1);
  };
  

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
    TString fName;       // histo names prefix
    
    ClassDef(AliFragFuncQAJetHistos, 1);
  };
  

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
    virtual void FillTrackQA(Float_t eta, Float_t phi, Float_t pt);
    virtual void AddToOutput(TList* list) const;

  private:
    
    Int_t   fNBinsPt;    // track QA histos bins
    Float_t fPtMin;      // track QA histos limits
    Float_t fPtMax;      // track QA histos limits
    Int_t   fNBinsEta;   // track QA histos bins
    Float_t fEtaMin;     // track QA histos limits
    Float_t fEtaMax;     // track QA histos limits
    Int_t   fNBinsPhi;   // track QA histos bins
    Float_t fPhiMin;     // track QA histos limits
    Float_t fPhiMax;     // track QA histos limits

    Float_t fHighPtThreshold; //  high pt track phi vs eta distribution

    TH2F*	fh2EtaPhi;       //! track phi vs eta 
    TH1F*	fh1Pt;           //! track transverse momentum 
    TH2F*	fh2HighPtEtaPhi; //! phi vs eta for high pt (>fgHighPtThreshold) tracks

    TString fName;               // histo names prefix
    
    ClassDef(AliFragFuncQATrackHistos, 1);
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
  
  virtual void   SetTrackTypeGen(Int_t i){fTrackTypeGen = i;}
  virtual void   SetJetTypeGen(Int_t i){fJetTypeGen = i;}

  virtual void   SetBranchGenJets(const char* c){fBranchGenJets = c;}
  virtual void   SetBranchRecJets(const char* c){fBranchRecJets = c;}

  virtual void   SetTrackCuts(Float_t trackPt = 0.15, Float_t trackEtaMin = -0.9, Float_t trackEtaMax = 0.9, Float_t trackPhiMin = 0., Float_t trackPhiMax = 2*TMath::Pi())
  {fTrackPtCut = trackPt; fTrackEtaMin = trackEtaMin; fTrackEtaMax = trackEtaMax; fTrackPhiMin = trackPhiMin; fTrackPhiMax = trackPhiMax;}
  virtual void   SetFilterMask(UInt_t i) {fFilterMask = i;}
  virtual void   SetJetCuts(Float_t jetPt = 5., Float_t jetEtaMin = -0.5, Float_t jetEtaMax = 0.5, Float_t jetPhiMin = 0., Float_t jetPhiMax = 2*TMath::Pi())
  {fJetPtCut = jetPt; fJetEtaMin = jetEtaMin; fJetEtaMax = jetEtaMax; fJetPhiMin = jetPhiMin; fJetPhiMax = jetPhiMax;}
  virtual void   SetDijetCuts(Float_t deltaPhi = 0., Float_t invMassMin = -1., Float_t invMassMax = -1., Float_t cdfCut = -1., Float_t eMeanMin = -1., Float_t eMeanMax = -1., Float_t eFraction = -1.)
  {fDijetDeltaPhiCut = deltaPhi; fDijetInvMassMin = invMassMin; fDijetInvMassMax = invMassMax; fDijetCDFcut = cdfCut; fDijetEMeanMin = eMeanMin; fDijetEMeanMax = eMeanMax; fDijetEFraction = eFraction;}
  virtual void   SetFFRadius(Float_t r = 0.4) { fFFRadius = r; }

  static  void   SetProperties(TH1* h,const char* x, const char* y);
  static  void   SetProperties(TH2* h,const char* x, const char* y,const char* z);

  void   SetHighPtThreshold(Float_t pt = 5.) { fQATrackHighPtThreshold = pt; }


  void   SetFFHistoBins(Int_t nJetPt = 55, Float_t jetPtMin = 5, Float_t jetPtMax = 60, 
			Int_t nPt = 70, Float_t ptMin = 0., Float_t ptMax = 70., 
			Int_t nXi = 70, Float_t xiMin = 0., Float_t xiMax = 7.,
			Int_t nZ = 22,  Float_t zMin = 0.,  Float_t zMax = 1.1)
  { fFFNBinsJetPt = nJetPt; fFFJetPtMin = jetPtMin; fFFJetPtMax = jetPtMax; 
    fFFNBinsPt = nPt; fFFPtMin = ptMin; fFFPtMax = ptMax;
    fFFNBinsXi = nXi; fFFXiMin = xiMin; fFFXiMax = xiMax;
    fFFNBinsZ  = nZ;  fFFZMin  = zMin;  fFFZMax  = zMax; }
  
  void  SetQAJetHistoBins(Int_t nPt = 200, Float_t ptMin = 0., Float_t ptMax = 200.,
			  Int_t nEta = 20, Float_t etaMin = -1.0, Float_t etaMax = 1.0,
			  Int_t nPhi = 60, Float_t phiMin = 0., Float_t phiMax = 2*TMath::Pi())
    { fQAJetNBinsPt = nPt; fQAJetPtMin = ptMin; fQAJetPtMax = ptMax;
      fQAJetNBinsEta = nEta; fQAJetEtaMin = etaMin; fQAJetEtaMax = etaMax;
      fQAJetNBinsPhi = nPhi; fQAJetPhiMin = phiMin; fQAJetPhiMax = phiMax; }
  
  void  SetQATrackHistoBins(Int_t nPt = 70, Float_t ptMin = 0., Float_t ptMax = 70.,
			    Int_t nEta = 20, Float_t etaMin = -1.0, Float_t etaMax = 1.0,
			    Int_t nPhi = 60, Float_t phiMin = 0., Float_t phiMax = 2*TMath::Pi())
  { fQATrackNBinsPt = nPt; fQATrackPtMin = ptMin; fQATrackPtMax = ptMax;
    fQATrackNBinsEta = nEta; fQATrackEtaMin = etaMin; fQATrackEtaMax = etaMax;
    fQATrackNBinsPhi = nPhi; fQATrackPhiMin = phiMin; fQATrackPhiMax = phiMax; }
  

  Float_t		GetFFRadius() const { return fFFRadius; }
  void			GetJetTracksTrackrefs(TList* l, AliAODJet* j);
  void			GetJetTracksPointing(TList* in, TList* out, AliAODJet* j, const Double_t r, Double_t& pt);  
  
  
 private:
    
  // Consts
  
  enum {kTrackUndef=0, kTrackAOD, kTrackAODCuts, kTrackKineAll, kTrackKineCharged, kTrackKineChargedAcceptance, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance};
  enum {kJetsUndef=0, kJetsRec, kJetsRecAcceptance, kJetsGen, kJetsGenAcceptance, kJetsKine, kJetsKineAcceptance};
  
  
  Int_t   GetListOfTracks(TList* list, Int_t type);
  Int_t	  GetListOfJets(TList* list, Int_t type);
  
  AliESDEvent* fESD;      // ESD event
  AliAODEvent* fAOD;      // AOD event
  AliMCEvent*  fMCEvent;  // MC event
  
  TString fBranchRecJets;         // branch name for reconstructed jets
  TString fBranchGenJets;         // branch name for generated jets
  
  Int_t fTrackTypeGen;      // type of generated tracks
  Int_t fJetTypeGen;        // type of generated jets

  UInt_t fFilterMask;	    // filter bit for selected tracks
	
  // track cuts
  Float_t fTrackPtCut;      // track transverse momentum cut
  Float_t fTrackEtaMin;     // track eta cut
  Float_t fTrackEtaMax;     // track eta cut
  Float_t fTrackPhiMin;     // track phi cut
  Float_t fTrackPhiMax;     // track phi cut
  
  // jet cuts
  Float_t fJetPtCut;        // jet transverse momentum cut
  Float_t fJetEtaMin;       // jet eta cut
  Float_t fJetEtaMax;       // jet eta cut
  Float_t fJetPhiMin;       // jet phi cut
  Float_t fJetPhiMax;       // jet phi cut

  
  Float_t fFFRadius;        // if radius > 0 construct FF from tracks within cone around jet axis, otherwise use trackRefs  
  
  // dijet cuts
  Float_t fDijetDeltaPhiCut;  // should be comment here 
  Float_t fDijetInvMassMin;   // should be comment here 
  Float_t fDijetInvMassMax;   // should be comment here 
  Float_t fDijetCDFcut;       // should be comment here  
  Float_t fDijetEMeanMin;     // should be comment here  
  Float_t fDijetEMeanMax;     // should be comment here  
  Float_t fDijetEFractionCut; // should be comment here  
  Float_t fDijetEFraction;    // should be comment here  
  
  TList* fTracksRec;      //! reconstructed tracks
  TList* fTracksRecCuts;  //! reconstructed tracks after cuts
  TList* fTracksGen;      //! generated tracks 
  
  TList* fJetsRec;        //! jets from reconstructed tracks
  TList* fJetsRecCuts;    //! jets from reonstructed tracks after jet cuts 
  TList* fJetsGen;        //! jets from generated tracks
  
  
  AliFragFuncQATrackHistos* fQATrackHistosRec;      //! track QA: reconstructed tracks
  AliFragFuncQATrackHistos* fQATrackHistosRecCuts;  //! track QA: reconstructed tracks after cuts
  AliFragFuncQATrackHistos* fQATrackHistosGen;      //! track QA: generated tracks
  
  AliFragFuncQAJetHistos*  fQAJetHistosRec;             //! jet QA: jets from reconstructed tracks
  AliFragFuncQAJetHistos*  fQAJetHistosRecCuts;         //! jet QA: jets from reconstructed tracks after jet cuts 
  AliFragFuncQAJetHistos*  fQAJetHistosRecCutsLeading;  //! jet QA: leading jet from reconstructed tracks after jet cuts 
  AliFragFuncQAJetHistos*  fQAJetHistosGen;             //! jet QA: jets from generated tracks  
  AliFragFuncQAJetHistos*  fQAJetHistosGenLeading;      //! jet QA: leading jet from generated tracks  
  
  AliFragFuncHistos*  fFFHistosRecCuts;         //! FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFHistosRecLeading;      //! FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFHistosRecLeadingTrack; //! FF reconstructed tracks after cuts: leading track pt / jet pt
  AliFragFuncHistos*  fFFHistosGen;             //! FF generated tracks after cuts 
  AliFragFuncHistos*  fFFHistosGenLeading;      //! FF generated tracks after cuts: all generated tracks pt / leading track pt  
  AliFragFuncHistos*  fFFHistosGenLeadingTrack; //! FF generated tracks after cuts: leading track pt / jet pt

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
  TH1F  *fh1nRecJetsCuts;         //! number of jets from reconstructed tracks per event 
  TH1F  *fh1nGenJets;             //! number of jets from generated tracks per event

  ClassDef(AliAnalysisTaskFragmentationFunction, 2);
};

#endif
