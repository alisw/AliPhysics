// *************************************************************************
// * Task for Fragmentation Function Analysis in PWG4 Jet Task Force Train *
// *************************************************************************

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
class TH3F;
class TProfile;
class THnSparse; 
class TRandom3;
class TArrayS;

#include "AliAnalysisTaskSE.h"

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
    virtual void FillFF(Float_t trackPt, Float_t jetPt,Bool_t incrementJetPt, Float_t norm = 0);
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
    virtual void FillTrackQA(Float_t eta, Float_t phi, Float_t pt, Bool_t weightPt = kFALSE, Float_t norm = 0.);
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
  
  //----------------------------------------
  class AliFragFuncIntraJetHistos : public TObject
  {

  public:
    
    AliFragFuncIntraJetHistos(const char* name = "IntraJethistos", 
	     Int_t nJetPt    = 0, Float_t jetPtMin    = 0, Float_t jetPtMax    = 0,
	     Int_t nPt       = 0, Float_t ptMin       = 0, Float_t ptMax       = 0,
	     Int_t nZ        = 0, Float_t zMin        = 0, Float_t zMax        = 0,
	     Int_t nCosTheta = 0, Float_t costhetaMin = 0, Float_t costhetaMax = 0,
	     Int_t nTheta    = 0, Float_t thetaMin    = 0, Float_t thetaMax    = 0,
	     Int_t nJt       = 0, Float_t jtMin       = 0, Float_t jtMax       = 0);
    AliFragFuncIntraJetHistos(const AliFragFuncIntraJetHistos& copy);
    AliFragFuncIntraJetHistos& operator=(const AliFragFuncIntraJetHistos &o);
    virtual ~AliFragFuncIntraJetHistos();
    
    virtual void DefineHistos();
    virtual void FillIntraJet(const TLorentzVector* trackV, const TLorentzVector* jetV, Float_t norm = 0);
    virtual void AddToOutput(TList* list) const;

  private:

    Int_t   fNBinsJetPt;    // IntraJet histos bins in jet pt
    Float_t fJetPtMin;      // IntraJet histos limits in jet pt
    Float_t fJetPtMax;      // IntraJet histos limits in jet pt
    Int_t   fNBinsPt;       // IntraJet histos bins in pt
    Float_t fPtMin;         // IntraJet histos limits in pt
    Float_t fPtMax;         // IntraJet histos limits in pt
    Int_t   fNBinsZ;        // IntraJet histos bins in z
    Float_t fZMin;          // IntraJet histos limits in z
    Float_t fZMax;          // IntraJet histos limits in z
    Int_t   fNBinsJt;       // IntraJet histos bins in jt
    Float_t fJtMin;         // IntraJet histos limits in jt
    Float_t fJtMax;         // IntraJet histos limits in jt
    Int_t   fNBinsTheta;    // IntraJet histos bins in theta
    Float_t fThetaMin;      // IntraJet histos limits in theta
    Float_t fThetaMax;      // IntraJet histos limits in theta
    Int_t   fNBinsCosTheta; // IntraJet histos bins in cos(theta)
    Float_t fCosThetaMin;   // IntraJet histos limits in cos(theta)
    Float_t fCosThetaMax;   // IntraJet histos limits in cos(theta)
  
    TH2F*   fh2CosTheta;    //! IntraJet: cos(theta) distribution
    TH2F*   fh2PtZ;       //! IntraJet: pt vs z distribution

    TH3F*   fh3ThetaZ;      //! IntraJet: theta, z, jet pt
    TH3F*   fh3JtTheta;     //! IntraJet: jt, theta, jet pt
    TH3F*   fh3JtZ;         //! IntraJet: jt, z, jet pt

    TString fNameIJ;         // histo names prefix
    
    ClassDef(AliFragFuncIntraJetHistos, 1);
  };

  //----------------------------------------
  class AliFragFuncDiJetHistos : public TObject
  {

    public:
    
    AliFragFuncDiJetHistos(const char* name = "DiJetHistos", Int_t kindSlices = 0,
	     Int_t nJetinvMass = 0, Float_t jetInvMassMin = 0, Float_t jetInvMassMax = 0,
	     Int_t nJetPt = 0, Float_t jetPtMin = 0, Float_t jetPtMax = 0,
	     Int_t nPt = 0, Float_t ptMin = 0, Float_t ptMax = 0,
	     Int_t nXi = 0, Float_t xiMin = 0, Float_t xiMax = 0,
	     Int_t nZ  = 0, Float_t zMin  = 0, Float_t zMax  = 0);
    AliFragFuncDiJetHistos(const AliFragFuncDiJetHistos& copy);
    AliFragFuncDiJetHistos& operator=(const AliFragFuncDiJetHistos &o);
    virtual ~AliFragFuncDiJetHistos();
    
    virtual void DefineDiJetHistos();
    virtual void FillDiJetFF(Int_t jetType, Float_t trackPt, Float_t jetPt, Double_t jetBin, Bool_t incrementJetPt);
    virtual void AddToOutput(TList* list) const;
    
    private:

    Int_t   fKindSlices;      // DJ kind of slices
    Int_t   fNBinsJetInvMass; // FF histos bins
    Float_t fJetInvMassMin;   // FF histos limits
    Float_t fJetInvMassMax;   // FF histos limits
    Int_t   fNBinsJetPt;      // FF histos bins
    Float_t fJetPtMin;        // FF histos limits
    Float_t fJetPtMax;        // FF histos limits
    Int_t   fNBinsPt;         // FF histos bins
    Float_t fPtMin;           // FF histos limits
    Float_t fPtMax;           // FF histos limits
    Int_t   fNBinsXi;         // FF histos bins
    Float_t fXiMin;           // FF histos limits
    Float_t fXiMax;           // FF histos limits
    Int_t   fNBinsZ;          // FF histos bins
    Float_t fZMin;            // FF histos limits
    Float_t fZMax;            // FF histos limits

    TH2F*   fh2TrackPtJet1; //! FF dijet : track transverse momentum of jet 1 
    TH2F*   fh2TrackPtJet2; //! FF dijet : track transverse momentum of jet 2 
    TH2F*   fh2TrackPtJet;  //! FF dijet : track transverse momentum of jets 1 and 2   
    TH1F*   fh1Jet1Pt;      //! jet 1 pt 
    TH1F*   fh1Jet2Pt;      //! jet 2 pt 
    TH1F*   fh1JetPt;       //! jet 1 and 2 pt 
    TH2F*   fh2Xi1;         //! FF dijet : xi of jet 1 
    TH2F*   fh2Xi2;         //! FF dijet : xi of jet 2
    TH2F*   fh2Xi;          //! FF dijet : xi of jet 1 and 2 
    TH2F*   fh2Z1;          //! FF dijet : z of jet 1   
    TH2F*   fh2Z2;          //! FF dijet : z of jet 2 
    TH2F*   fh2Z;           //! FF dijet : z of jet 1 and 2 
    TH2F*   fh2Pt1;         //! FF dijet : z of jet 1   
    TH2F*   fh2Pt2;         //! FF dijet : z of jet 2 
    TH2F*   fh2Pt;          //! FF dijet : z of jet 1 and 2 

    TString fNameDJ;        // histo names prefix
    
    ClassDef(AliFragFuncDiJetHistos, 1);
  };

  //----------------------------------------
  class AliFragFuncQADiJetHistos : public TObject
  {
    public:
    AliFragFuncQADiJetHistos(const char* name = "QADiJetHistos", Int_t kindSlices = 0,
	     Int_t nInvMass = 0, Float_t invMassMin = 0, Float_t invMassMax = 0,
	     Int_t nJetPt = 0,   Float_t jetPtMin = 0, Float_t jetPtMax = 0,
	     Int_t nDeltaPhi = 0, Float_t deltaPhiMin = 0, Float_t deltaPhiMax = 0,
	     Int_t nDeltaEta = 0, Float_t deltaEtaMin = 0, Float_t deltaEtaMax = 0,
	     Int_t nDeltaPt  = 0, Float_t deltaPtMin  = 0, Float_t deltaPtMax  = 0,
	     Int_t nInBal    = 0, Float_t inBalMin  = 0, Float_t inBalMax  = 0);
    AliFragFuncQADiJetHistos(const AliFragFuncQADiJetHistos& copy);
    AliFragFuncQADiJetHistos& operator=(const AliFragFuncQADiJetHistos &o);
    virtual ~AliFragFuncQADiJetHistos();
    
    virtual void DefineQADiJetHistos();
    virtual void FillDiJetQA(Double_t invMass, Double_t deltaPhi, Double_t deltaEta, Double_t deltaPt, Double_t inBal, Double_t jetBin);
    virtual void AddToOutput(TList* list) const;
    
    private:
    
    Int_t   fKindSlices;       // DJ kind of slices
    Int_t   fNBinsJetInvMass;  // FF histos bins in jet invariant mass
    Float_t fJetInvMassMin;    // FF histos limits in jet invariant mass
    Float_t fJetInvMassMax;    // FF histos limits in jet invariant mass
    Int_t   fNBinsJetPt;       // FF histos bins in jet pt
    Float_t fJetPtMin;         // FF histos limits in jet pt
    Float_t fJetPtMax;         // FF histos limits in jet pt
    Int_t   fNBinsDeltaPhi;    // FF histos bins in jet delta phi
    Float_t fDeltaPhiMin;      // FF histos limits in jet delta phi
    Float_t fDeltaPhiMax;      // FF histos limits in jet delta phi
    Int_t   fNBinsDeltaEta;    // FF histos bins in jet delta eta
    Float_t fDeltaEtaMin;      // FF histos limits in jet delta eta
    Float_t fDeltaEtaMax;      // FF histos limits in jet delta eta
    Int_t   fNBinsDeltaPt;     // FF histos bins in jet delta pt
    Float_t fDeltaPtMin;       // FF histos limits in jet delta pt
    Float_t fDeltaPtMax;       // FF histos limits in jet delta pt
    Int_t   fNBinsInBal;       // FF histos bins in jet delta pt
    Float_t fInBalMin;         // FF histos limits in pt inbalance
    Float_t fInBalMax;         // FF histos limits in pt inbalance

    TH2F*   fh2InvMass;        // FF dijet invariant mass histos
    TH2F*   fh2DeltaPhi;       // FF dijet delta phi histos
    TH2F*   fh2DeltaEta;       // FF dijet delta eta histos
    TH2F*   fh2DeltaPt;        // FF dijet delta pt histos
    TH2F*   fh2InBal;        // FF dijet delta pt histos

    TString fNameQADJ;         // histo names prefix
    
    ClassDef(AliFragFuncQADiJetHistos, 1);
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
  
  virtual void   SetTrackTypeGen(Int_t i){fTrackTypeGen = i;}
  virtual void   SetJetTypeGen(Int_t i){fJetTypeGen = i;}
  virtual void   SetJetTypeRecEff(Int_t i){fJetTypeRecEff = i;}

  virtual void   SetBranchRecBackJets(const char* c){fBranchRecBackJets = c;}
  virtual void   SetBranchRecBackClusters(const char* c){fBranchRecBckgClusters = c;}
  virtual void   SetBranchGenJets(const char* c){fBranchGenJets = c;}
  virtual void   SetBranchRecJets(const char* c){fBranchRecJets = c;}

  virtual void   SetTrackCuts(Float_t trackPt = 0.15, Float_t trackEtaMin = -0.9, Float_t trackEtaMax = 0.9, 
			      Float_t trackPhiMin = 0., Float_t trackPhiMax = 2*TMath::Pi())
  {fTrackPtCut = trackPt; fTrackEtaMin = trackEtaMin; fTrackEtaMax = trackEtaMax; 
    fTrackPhiMin = trackPhiMin; fTrackPhiMax = trackPhiMax;}
  virtual void   SetFilterMask(UInt_t i) {fFilterMask = i;}
  virtual void   UsePhysicsSelection(Bool_t b) {fUsePhysicsSelection = b;}
  virtual void   SetEventClass(Int_t i){fEventClass = i;}
  virtual void   SetMaxVertexZ(Float_t z){fMaxVertexZ = z;}
  virtual void   SetJetCuts(Float_t jetPt = 5., Float_t jetEtaMin = -0.5, Float_t jetEtaMax = 0.5, 
			    Float_t jetPhiMin = 0., Float_t jetPhiMax = 2*TMath::Pi())
  {fJetPtCut = jetPt; fJetEtaMin = jetEtaMin; fJetEtaMax = jetEtaMax; 
    fJetPhiMin = jetPhiMin; fJetPhiMax = jetPhiMax;}
  virtual void   SetDiJetCuts(Int_t cutType = 1, Float_t deltaPhiCut = 0.,  
			      Float_t cdfCut = 0.5, Float_t ptFractionCut = 0.6)
  {fDiJetCut = cutType; fDiJetDeltaPhiCut = deltaPhiCut;  
    fDiJetCDFCut = cdfCut; fDiJetPtFractionCut = ptFractionCut;}
  virtual void   SetKindSlices(Int_t slice = 1) {fDiJetKindBins = slice;}

  virtual void   SetFFRadius(Float_t r = 0.4) { fFFRadius = r; }
  virtual void   SetFFBckgRadius(Float_t r = 0.7) { fFFBckgRadius = r; }
  virtual void   SetBckgMode(Bool_t bg = 1) { fBckgMode = bg; }
  virtual void   SetBckgType(Int_t bg0 = 0, Int_t bg1 = 1,Int_t bg2 = 2, Int_t bg3 = 3, Int_t bg4 = 4) 
  { fBckgType[0] = bg0; fBckgType[1] = bg1; fBckgType[2] = bg2; fBckgType[3] = bg3; fBckgType[4] = bg4; }
  virtual void   SetIJMode(Int_t ij = 1)      {fIJMode = ij;}
  virtual void   SetQAMode(Int_t qa = 3)      {fQAMode = qa;}
  virtual void   SetFFMode(Int_t ff = 1)      {fFFMode = ff;}
  virtual void   SetDJMode(Int_t dj = 1)      {fDJMode = dj;}
  virtual void   SetEffMode(Int_t eff = 1)    {fEffMode = eff;}
  virtual void   SetPhiCorrMode(Int_t pc = 1) {fPhiCorrMode = pc;}
  virtual void   SetBckgSubMethod(Int_t bg = 0) {fBckgSubMethod = bg;}

  virtual void   UseRecEffRecJetPtBins(Bool_t useRec = kFALSE)  { fUseRecEffRecJetPtBins = useRec; }
  virtual void   UseResponseRecJetPtBins(Bool_t useRec = kTRUE) { fUseResponseRecJetPtBins = useRec; }

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
  
  void  SetPhiCorrHistoBins(Int_t nPt = 200, Float_t ptMin = 0., Float_t ptMax = 200.,
			         Int_t nEta = 1, Float_t etaMin = -0.9, Float_t etaMax = 0.9,
			         Int_t nPhi = 64, Float_t phiMin = -3.2, Float_t phiMax = 3.2)
  { fPhiCorrNBinsPt = nPt; fPhiCorrPtMin = ptMin; fPhiCorrPtMax = ptMax;
    fPhiCorrNBinsEta = nEta; fPhiCorrEtaMin = etaMin; fPhiCorrEtaMax = etaMax;
    fPhiCorrNBinsPhi = nPhi; fPhiCorrPhiMin = phiMin; fPhiCorrPhiMax = phiMax; }

  void   SetIJHistoBins(Int_t nJetPt = 23, Float_t jetPtMin = 5, Float_t jetPtMax = 120, Int_t nPt = 120, Float_t ptMin = 0., Float_t ptMax = 120., 
			Int_t nZ = 22,  Float_t zMin = 0.,  Float_t zMax = 1.1,	Int_t nCosTheta = 100,  Float_t costhetaMin = 0.,  Float_t costhetaMax = 1.,
			Int_t nTheta = 200,  Float_t thetaMin = -0.5,  Float_t thetaMax = 1.5, Int_t nJt = 25,  Float_t jtMin = 0.,  Float_t jtMax = 5.)
  { fIJNBinsJetPt = nJetPt; fIJJetPtMin = jetPtMin; fIJJetPtMax = jetPtMax; fIJNBinsPt = nPt; fIJPtMin = ptMin; fIJPtMax = ptMax;
    fIJNBinsZ = nZ; fIJZMin = zMin; fIJZMax = zMax;fIJNBinsCosTheta  = nCosTheta;  fIJCosThetaMin  = costhetaMin;  fIJCosThetaMax  = costhetaMax;
    fIJNBinsTheta  = nTheta;  fIJThetaMin  = thetaMin;  fIJThetaMax  = thetaMax; fIJNBinsJt  = nJt;  fIJJtMin  = jtMin;  fIJJtMax  = jtMax; }

  void SetDiJetHistoBins(Int_t nJetInvMass = 245, Float_t jetInvMassMin = 5, Float_t jetInvMassMax = 250, Int_t nJetPt = 245, Float_t jetPtMin = 5, 
			 Float_t jetPtMax = 250, Int_t nPt = 200, Float_t ptMin = 0., Float_t ptMax = 200., Int_t nXi = 70, Float_t xiMin = 0., 
			 Float_t xiMax = 7., Int_t nZ = 22,  Float_t zMin = 0.,  Float_t zMax = 1.1)
  {
    fDiJetNBinsJetInvMass = nJetInvMass; fDiJetJetInvMassMin = jetInvMassMin; fDiJetJetInvMassMax = jetInvMassMax; fDiJetNBinsJetPt = nJetPt; fDiJetJetPtMin = jetPtMin; 
    fDiJetJetPtMax = jetPtMax; fDiJetNBinsPt = nPt; fDiJetPtMin = ptMin; fDiJetPtMax = ptMax; fDiJetNBinsXi = nXi; fDiJetXiMin = xiMin; 
    fDiJetXiMax = xiMax; fDiJetNBinsZ = nZ; fDiJetZMin = zMin; fDiJetZMax = zMax;
  }

  void SetQADiJetHistoBins(Int_t nInvMass = 245, Float_t invMassMin = 5., Float_t invMassMax = 250., Int_t nJetPt = 245, Float_t jetPtMin = 5, Float_t jetPtMax = 250, 
			   Int_t nDeltaPhi = 100, Float_t deltaPhiMin = 0., Float_t deltaPhiMax = TMath::Pi(), Int_t nDeltaEta = 22, Float_t deltaEtaMin = 0., 
			   Float_t deltaEtaMax = 1.1, Int_t nDeltaPt = 100, Float_t deltaPtMin = 0., Float_t deltaPtMax = 100.,
			   Int_t nInBal = 22, Float_t inBalMin = -1.1, Float_t inBalMax = 1.1)
  {
    fQADiJetNBinsInvMass = nInvMass; fQADiJetInvMassMin = invMassMin; fQADiJetInvMassMax = invMassMax; fQADiJetNBinsJetPt = nJetPt; fQADiJetJetPtMin = jetPtMin; 
    fQADiJetJetPtMax = jetPtMax; fQADiJetNBinsDeltaPhi = nDeltaPhi; fQADiJetDeltaPhiMin = deltaPhiMin; fQADiJetDeltaPhiMax = deltaPhiMax; fQADiJetNBinsDeltaEta = nDeltaEta; 
    fQADiJetDeltaEtaMin = deltaEtaMin; fQADiJetDeltaEtaMax = deltaEtaMax; fQADiJetNBinsDeltaPt = nDeltaPt; fQADiJetDeltaPtMin = deltaPtMin; fQADiJetDeltaPtMax = deltaPtMax;
    fQADiJetNBinsInBal = nInBal; fQADiJetInBalMin = inBalMin; fQADiJetInBalMax = inBalMax;
  }

  Float_t  GetFFRadius() const { return fFFRadius; }
  Float_t  GetFFBckgRadius() const { return fFFBckgRadius; }
  void	   GetJetTracksTrackrefs(TList* l, const AliAODJet* j);
  void	   GetJetTracksPointing(TList* in, TList* out, const AliAODJet* j, const Double_t r, Double_t& pt);  
  void     GetTracksOutOfNJets(Int_t nCases, TList* in, TList* out, TList* jets, Double_t& pt);
  void     GetTracksOutOfNJetsStat(Int_t nCases, TList* in, TList* out, TList* jets, Double_t& pt, Double_t &normFactor);
  void     GetTracksTiltedwrpJetAxis(Float_t alpha, TList* inputlist, TList* outputlist, AliAODJet* jet, Double_t radius, Double_t& sumPt);
  void     GetTracksTiltedwrpJetAxisWindow(Float_t alpha, TList* inputlist, TList* outputlist, AliAODJet* jet, Double_t radius, Double_t& sumPt, Double_t &normFactor);

  Double_t GetDiJetBin(Double_t invMass, Double_t leadingJetPt, Double_t eMean, Int_t kindSlices); // function to find which bin fill
  Double_t InvMass(const AliAODJet* jet1, const AliAODJet* jet2);
  void     AssociateGenRec(TList* tracksAODMCCharged,TList* tracksRec, TArrayI& indexAODTr,TArrayI& indexMCTr,TArrayS& isGenPrim);
  void     FillSingleTrackRecEffHisto(AliFragFuncQATrackHistos* trackQAGen, AliFragFuncQATrackHistos* trackQARec, TList* tracksGen, const TArrayI& indexAODTr, const TArrayS& isGenPrim);
  void     FillJetTrackRecEffHisto(TObject* histGen,TObject* histRec,Double_t jetPtGen,Double_t jetPtRec, TList* jetTrackList, const TList* tracksGen,
				   const TArrayI& indexAODTr,const TArrayS& isGenPrim, const Bool_t useRecJetPt);

  void     FillSingleTrackResponse(THnSparse* hnResponse, TList* tracksGen, TList* tracksRec, const TArrayI& indexAODTr, const TArrayS& isGenPrim);

  void     FillJetTrackResponse(THnSparse* hnResponsePt, THnSparse* hnResponseZ, THnSparse* hnResponseXi, 
				Double_t jetPtGen, Double_t jetPtRec, TList* jetTrackList,
				const TList* tracksGen, TList* tracksRec, const TArrayI& indexAODTr, const TArrayS& isGenPrim,const Bool_t useRecJetPt);
  
  // 


  Float_t  CalcJetArea(const Float_t etaJet, const Float_t rc) const;
  void     GetClusterTracksOutOf1Jet(AliAODJet* jet, TList* outputlist, Double_t &normFactor);
  void     GetClusterTracksMedian(TList* outputlist, Double_t &normFactor);

  void     FillBckgHistos(Int_t type, TList* inputtracklist, TList* inputjetlist, AliAODJet* jet, 
			  Float_t leadTrackPt, TLorentzVector* leadTrackV, AliFragFuncHistos* ffbckghistocuts,
			  AliFragFuncHistos* ffbckghistoleading,AliFragFuncIntraJetHistos* ijbckghistocuts, 
			  AliFragFuncIntraJetHistos* ijbckghistoleading,AliFragFuncQATrackHistos* qabckghistos);    
  AliAODJet* GetAODBckgSubJet(AliAODJet* jet, Int_t method);

  // Consts
  
  enum {kTrackUndef=0, kTrackAOD, kTrackAODQualityCuts, kTrackAODCuts, kTrackKineAll, kTrackKineCharged, kTrackKineChargedAcceptance, 
	kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance};
  enum {kJetsUndef=0, kJetsRec, kJetsRecAcceptance, kJetsGen, kJetsGenAcceptance, kJetsKine, kJetsKineAcceptance};
  enum {kBckgPerp=0, kBckgOutLJ, kBckgOut2J, kBckgClusters, kBckgClustersOutLeading, kBckgOut3J, kBckgOutAJ, kBckgOutLJStat, 
	kBckgOut2JStat, kBckgOut3JStat, kBckgOutAJStat,  kBckgASide, kBckgASideWindow, kBckgPerpWindow};

 
 private:
  
  Int_t   GetListOfTracks(TList* list, Int_t type);
  Int_t	  GetListOfJets(TList* list, Int_t type);
  Int_t   GetListOfBckgJets(TList *list, Int_t type);

  AliESDEvent* fESD;      // ESD event
  AliAODEvent* fAOD;      // AOD event
  //AliMCEvent*  fMCEvent;  // MC event
  
  TString fBranchRecJets;         // branch name for reconstructed jets
  TString fBranchRecBackJets;     // branch name for reconstructed background jets
  TString fBranchRecBckgClusters; // branch name for reconstructed background clusters 
  TString fBranchGenJets;         // branch name for generated jets
  
  Int_t   fTrackTypeGen;        // type of generated tracks
  Int_t   fJetTypeGen;          // type of generated jets

  Int_t   fJetTypeRecEff;       // type of jets used for filling reconstruction efficiency histos

  UInt_t  fFilterMask;	        // filter bit for selected tracks
  Bool_t  fUsePhysicsSelection; // switch for event selection
  Int_t   fEventClass;          // event class to be looked at for this instace of the task
  Float_t fMaxVertexZ;          // maximum abs(z) position of primiary vertex [cm]

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

  // di-jet cuts
  Int_t   fDiJetCut;            // dijet cut selection
  Float_t fDiJetDeltaPhiCut;    // delta phi cut value
  Float_t fDiJetPtFractionCut;  // fraction of pt cut value
  Float_t fDiJetCDFCut;         // cdf cut value

  Int_t   fDiJetKindBins;       // type of bins: invmass, etleading, emean

  Float_t fFFRadius;        // if radius > 0 construct FF from tracks within cone around jet axis, otherwise use trackRefs  
  Float_t fFFBckgRadius;    // compute background outside cone of this radius around jet axes
  Bool_t  fBckgMode;        // Set background subtraction mode
  Int_t   fBckgType[5];     // Set background subtraction mode
  Int_t   fIJMode;          // Set intrajet mode
  Int_t   fQAMode;          // QA mode: 0x00=0 none, 0x01=1 track qa, 0x10=2 track qa, 0x11=3 both
  Int_t   fFFMode;          // fragmentation function mode
  Int_t   fDJMode;          // dijet mode: 0x00=0 none, 0x01=1 di-jet, 0x10=2 di-jet qa, 0x11=3 both
  Int_t   fEffMode;         // efficiency mode
  Int_t   fPhiCorrMode;     // track phi correlation mode

  Bool_t  fUseRecEffRecJetPtBins;   // bin track reconstruction efficiency in reconstructed/generated jet pt bins 
  Bool_t  fUseResponseRecJetPtBins; // bin track response matrix in reconstructed/generated jet pt bins 

  Float_t fAvgTrials;       // average number of trials per event
  
  TList* fTracksRec;            //! reconstructed tracks
  TList* fTracksRecCuts;        //! reconstructed tracks after cuts
  TList* fTracksGen;            //! generated tracks 
  TList* fTracksAODMCCharged;   //! AOD MC tracks 
  TList* fTracksRecQualityCuts; //! reconstructed tracks after quality cuts, no acceptance/pt cut

  
  TList* fJetsRec;        //! jets from reconstructed tracks
  TList* fJetsRecCuts;    //! jets from reonstructed tracks after jet cuts 
  TList* fJetsGen;        //! jets from generated tracks
  TList* fJetsRecEff;     //! jets used for reconstruction efficiency histos 
  TList* fBckgJetsRec;      //! jets from reconstructed tracks
  TList* fBckgJetsRecCuts;  //! jets from reonstructed tracks after jet cuts
  TList* fBckgJetsGen;      //! jets from generated tracks
 
  
  AliFragFuncQATrackHistos* fQATrackHistosRec;      //! track QA: reconstructed tracks
  AliFragFuncQATrackHistos* fQATrackHistosRecCuts;  //! track QA: reconstructed tracks after cuts
  AliFragFuncQATrackHistos* fQATrackHistosGen;      //! track QA: generated tracks
  
  AliFragFuncQAJetHistos*  fQAJetHistosRec;             //! jet QA: jets from reconstructed tracks
  AliFragFuncQAJetHistos*  fQAJetHistosRecCuts;         //! jet QA: jets from reconstructed tracks after jet cuts 
  AliFragFuncQAJetHistos*  fQAJetHistosRecCutsLeading;  //! jet QA: leading jet from reconstructed tracks after jet cuts 
  AliFragFuncQAJetHistos*  fQAJetHistosGen;             //! jet QA: jets from generated tracks  
  AliFragFuncQAJetHistos*  fQAJetHistosGenLeading;      //! jet QA: leading jet from generated tracks  
  AliFragFuncQAJetHistos*  fQAJetHistosRecEffLeading;   //! jet QA: leading jet used for reconstruction efficiency histos  
  

  AliFragFuncHistos*  fFFHistosRecCuts;         //! FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFHistosRecLeading;      //! FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFHistosRecLeadingTrack; //! FF reconstructed tracks after cuts: leading track pt / jet pt
  AliFragFuncHistos*  fFFHistosGen;             //! FF generated tracks after cuts 
  AliFragFuncHistos*  fFFHistosGenLeading;      //! FF generated tracks after cuts: all generated tracks pt / leading track pt  
  AliFragFuncHistos*  fFFHistosGenLeadingTrack; //! FF generated tracks after cuts: leading track pt / jet pt


  AliFragFuncIntraJetHistos*  fIJHistosRecCuts;         //! IJ reconstructed tracks after cuts 
  AliFragFuncIntraJetHistos*  fIJHistosRecLeading;      //! IJ reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncIntraJetHistos*  fIJHistosRecLeadingTrack; //! IJ reconstructed tracks after cuts: leading track pt / jet pt
  AliFragFuncIntraJetHistos*  fIJHistosGen;             //! IJ generated tracks after cuts 
  AliFragFuncIntraJetHistos*  fIJHistosGenLeading;      //! IJ generated tracks after cuts: all generated tracks pt / leading track pt  
  AliFragFuncIntraJetHistos*  fIJHistosGenLeadingTrack; //! IJ generated tracks after cuts: leading track pt / jet pt

  AliFragFuncDiJetHistos* fFFDiJetHistosRecCuts;         //! DiJet FF reconstructed tracks after cuts
  AliFragFuncDiJetHistos* fFFDiJetHistosRecLeading;      //! DiJet FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt
  AliFragFuncDiJetHistos* fFFDiJetHistosRecLeadingTrack; //! DiJet FF reconstructed tracks after cuts: leading track pt / jet pt

  AliFragFuncDiJetHistos* fFFDiJetHistosGen;             //! DiJet FF generated tracks after cuts 
  AliFragFuncDiJetHistos* fFFDiJetHistosGenLeading;      //! DiJet FF generated tracks after cuts: all generated tracks pt / leading track pt 
  AliFragFuncDiJetHistos* fFFDiJetHistosGenLeadingTrack; //! DiJet FF generated tracks after cuts: leading track pt / jet pt 

  AliFragFuncQADiJetHistos* fQADiJetHistosRecCuts;       //! Dijet QA : reconstructed tracks after cuts
  AliFragFuncQADiJetHistos* fQADiJetHistosGen;           //! DiJet QA: jets from generated tracks  

  AliFragFuncQATrackHistos* fPhiCorrHistosJetArea;        //! tracks in area of leading jet (phi = phi_jet - phi_track, eta = eta_track)
  AliFragFuncQATrackHistos* fPhiCorrHistosTransverseArea; //! tracks in area transverse region (shift of phi by 90�)
  AliFragFuncQATrackHistos* fPhiCorrHistosAwayArea;       //! tracks in area in away region (shift of phi by 180�)

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
  
  Int_t   fIJNBinsJetPt;    // IJ histos bins
  Float_t fIJJetPtMin;      // IJ histos limits
  Float_t fIJJetPtMax;      // IJ histos limits

  Int_t   fIJNBinsPt;       // IJ histos bins
  Float_t fIJPtMin;         // IJ histos limits
  Float_t fIJPtMax;         // IJ histos limits

  Int_t   fIJNBinsZ;        // IJ histos bins
  Float_t fIJZMin;          // IJ histos limits
  Float_t fIJZMax;          // IJ histos limits

  Int_t   fIJNBinsCosTheta; // IJ histos bins
  Float_t fIJCosThetaMin;   // IJ histos limits
  Float_t fIJCosThetaMax;   // IJ histos limits

  Int_t   fIJNBinsTheta;    // IJ histos bins
  Float_t fIJThetaMin;      // IJ histos limits
  Float_t fIJThetaMax;      // IJ histos limits
  
  Int_t   fIJNBinsJt;       // IJ histos bins
  Float_t fIJJtMin;         // IJ histos limits
  Float_t fIJJtMax;         // IJ histos limits

  Int_t   fDiJetNBinsJetInvMass; // FF dijet histos bins
  Float_t fDiJetJetInvMassMin;   // FF dijet histos limits
  Float_t fDiJetJetInvMassMax;   // FF dijet histos limits
  Int_t   fDiJetNBinsJetPt;      // FF dijet histos bins
  Float_t fDiJetJetPtMin;        // FF dijet histos limits
  Float_t fDiJetJetPtMax;        // FF dijet histos limits
  Int_t   fDiJetNBinsPt;         // FF dijet histos bins
  Float_t fDiJetPtMin;           // FF dijet histos limits
  Float_t fDiJetPtMax;           // FF dijet histos limits
  Int_t   fDiJetNBinsXi;         // FF dijet histos bins
  Float_t fDiJetXiMin;           // FF dijet histos limits
  Float_t fDiJetXiMax;           // FF dijet histos limits
  Int_t   fDiJetNBinsZ;          // FF dijet histos bins
  Float_t fDiJetZMin;            // FF dijet histos limits
  Float_t fDiJetZMax;            // FF dijet histos limits

  Int_t   fQADiJetNBinsInvMass;  // dijet QA histos bins
  Float_t fQADiJetInvMassMin;    // dijet QA histos limits
  Float_t fQADiJetInvMassMax;    // dijet QA histos limits

  Int_t   fQADiJetNBinsJetPt;    // dijet QA histos bins
  Float_t fQADiJetJetPtMin;      // dijet QA histos limits
  Float_t fQADiJetJetPtMax;      // dijet QA histos limits

  Int_t   fQADiJetNBinsDeltaPhi; // dijet QA histos bins
  Float_t fQADiJetDeltaPhiMin;   // dijet QA histos limits
  Float_t fQADiJetDeltaPhiMax;   // dijet QA histos limits

  Int_t   fQADiJetNBinsDeltaEta; // dijet QA histos bins
  Float_t fQADiJetDeltaEtaMin;   // dijet QA histos limits
  Float_t fQADiJetDeltaEtaMax;   // dijet QA histos limits

  Int_t   fQADiJetNBinsDeltaPt;  // dijet QA histos bins
  Float_t fQADiJetDeltaPtMin;    // dijet QA histos limits
  Float_t fQADiJetDeltaPtMax;    // dijet QA histos limits

  Int_t   fQADiJetNBinsInBal;  // dijet QA histos bins
  Float_t fQADiJetInBalMin;    // dijet QA histos limits
  Float_t fQADiJetInBalMax;    // dijet QA histos limits

  // phi correlation
  Int_t   fPhiCorrNBinsPt;  // track related to jet histos bins 
  Float_t fPhiCorrPtMin;    // track related to jet histos limits
  Float_t fPhiCorrPtMax;    // track related to jet histos limits
  
  Int_t   fPhiCorrNBinsEta; // track related to jet histos bins
  Float_t fPhiCorrEtaMin;   // track related to jet histos limits
  Float_t fPhiCorrEtaMax;   // track related to jet histos limits
  
  Int_t   fPhiCorrNBinsPhi; // track related to jet histos bins
  Float_t fPhiCorrPhiMin;   // track related to jet histos limits
  Float_t fPhiCorrPhiMax;   // track related to jet histos limits

  // Histograms
  TList	*fCommonHistList;         // List of common histos
  
  TH1F  *fh1EvtSelection;         //! event cuts 
  TH1F	*fh1VertexNContributors;  //! NContributors to prim vertex
  TH1F	*fh1VertexZ;              //! prim vertex z distribution
  TH1F	*fh1EvtMult;              //! number of reconstructed tracks after cuts 
  TH1F	*fh1EvtCent;              //! centrality percentile 

  TH2F  *fh2TrackPtVsDCAXY;       //! track pt vs DCA 
  TH2F  *fh2TrackPtVsDCAZ;        //! track pt vs DCA


  TProfile* fh1Xsec;              //! pythia cross section and trials
  TH1F*     fh1Trials;            //! sum of trials
  TH1F*     fh1PtHard;            //! pt hard of the event
  TH1F*     fh1PtHardTrials;      //! pt hard of the event

  TH1F  *fh1nRecJetsCuts;         //! number of jets from reconstructed tracks per event 
  TH1F  *fh1nGenJets;             //! number of jets from generated tracks per event
  TH1F  *fh1nRecEffJets;          //! number of jets for reconstruction eff per event
  TH1F  *fh1nRecBckgJetsCuts;     //! number of jets from reconstructed tracks per event
  TH1F  *fh1nGenBckgJets;         //! number of jets from generated tracks per event
  TH2F  *fh2PtRecVsGenPrim;       //! association rec/gen MC: rec vs gen pt 

  // tracking efficiency 
  
  AliFragFuncQATrackHistos* fQATrackHistosRecEffGen;      //! tracking efficiency: generated primaries 
  AliFragFuncQATrackHistos* fQATrackHistosRecEffRec;      //! tracking efficiency: reconstructed primaries

  AliFragFuncHistos*  fFFHistosRecEffGen;                 //! tracking efficiency: FF generated primaries  
  AliFragFuncHistos*  fFFHistosRecEffRec;                 //! tracking efficiency: FF reconstructed primaries

  // momentum resolution 
  THnSparse* fhnResponseSinglePt;    //! single track response pt
  THnSparse* fhnResponseJetTrackPt;  //! jet track response pt 
  THnSparse* fhnResponseJetZ;        //! jet track response z 
  THnSparse* fhnResponseJetXi;       //! jet track response xi


  // Background
  TH1F  *fh1OutLeadingMult;       //! background multiplicity outside leading jet
  TH1F  *fh1OutLeadingStatMult;   //! background multiplicity outside leading jet (stat)
  TH1F  *fh1PerpMult;             //! background multiplicity perpendicular to the leading jet
  TH1F  *fh1ASideMult;            //! background multiplicity perpendicular to the leading jet
  TH1F  *fh1ASideWindowMult;      //! background multiplicity perpendicular to the leading jet
  TH1F  *fh1PerpWindowMult;       //! background multiplicity perpendicular to the leading jet
  TH1F  *fh1Out2JetsMult;         //! background multiplicity outside 2 jets
  TH1F  *fh1Out3JetsMult;         //! background multiplicity outside 3 jets
  TH1F  *fh1MedianClustersMult;   //! background multiplicity median cluster
  TH1F  *fh1OutClustersMult;      //! background multiplicity clusters outside leading jet

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
  AliFragFuncHistos*  fFFBckgHisto0RecLeading;    //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFBckgHisto0Gen;           //! Bckg (outside leading jet or 2 jets or more) FF generated tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto0GenLeading;    //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFBckgHisto1RecCuts;       //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto1RecLeading;    //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFBckgHisto1Gen;           //! Bckg (outside leading jet or 2 jets or more) FF generated tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto1GenLeading;    //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFBckgHisto2RecCuts;       //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto2RecLeading;    //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFBckgHisto2Gen;           //! Bckg (outside leading jet or 2 jets or more) FF generated tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto2GenLeading;    //! Bckg (outside leading jet or 2 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFBckgHisto3RecCuts;       //! Bckg (outside leading jet or 3 jets or more) FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto3RecLeading;    //! Bckg (outside leading jet or 3 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFBckgHisto3Gen;           //! Bckg (outside leading jet or 3 jets or more) FF generated tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto3GenLeading;    //! Bckg (outside leading jet or 3 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt
  AliFragFuncHistos*  fFFBckgHisto4RecCuts;       //! Bckg (outside leading jet or 4 jets or more) FF reconstructed tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto4RecLeading;    //! Bckg (outside leading jet or 4 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt  
  AliFragFuncHistos*  fFFBckgHisto4Gen;           //! Bckg (outside leading jet or 4 jets or more) FF generated tracks after cuts 
  AliFragFuncHistos*  fFFBckgHisto4GenLeading;    //! Bckg (outside leading jet or 4 jets or more) FF reconstructed tracks after cuts: all reconstructed tracks pt / leading track pt


  AliFragFuncIntraJetHistos*  fIJBckgHisto0RecCuts;    //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto0RecLeading; //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto0Gen;        //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto0GenLeading; //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto1RecCuts;    //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto1RecLeading; //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto1Gen;        //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto1GenLeading; //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto2RecCuts;    //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto2RecLeading; //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto2Gen;        //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto2GenLeading; //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto3RecCuts;    //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto3RecLeading; //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto3Gen;        //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto3GenLeading; //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto4RecCuts;    //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto4RecLeading; //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto4Gen;        //!
  AliFragFuncIntraJetHistos*  fIJBckgHisto4GenLeading; //!

  TRandom3*                   fRandom;          // TRandom3 for background estimation 
  Int_t                       fBckgSubMethod;   // Bckg method: 1 = leading jet excluded, 2 = 2 most energetic jets excluded        

  ClassDef(AliAnalysisTaskFragmentationFunction, 10);
};

#endif
