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
class TProfile;
class THnSparse; 

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
    virtual void FillTrackQA(Float_t eta, Float_t phi, Float_t pt);
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

    TString fNameQAT;         // histo names prefix
    
    ClassDef(AliFragFuncQATrackHistos, 1);
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
    virtual void FillIntraJet(const TLorentzVector* trackV, const TLorentzVector* jetV);
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
  
    TH2F*  fh2Theta;         //! IntraJet: theta distribution
    TH2F*  fh2CosTheta;      //! IntraJet: cos(theta) distribution
    TH2F*  fh2Jt;            //! IntraJet: jt distribution
    TH2F*  fh2PtvsZ;         //! IntraJet: pt vs z distribution

    THnSparseF* fhnIntraJet; //! IntraJet
    Int_t fnDim;             // HnSparseF dimensions 

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
	     Int_t nDeltaPt  = 0, Float_t deltaPtMin  = 0, Float_t deltaPtMax  = 0);
    AliFragFuncQADiJetHistos(const AliFragFuncQADiJetHistos& copy);
    AliFragFuncQADiJetHistos& operator=(const AliFragFuncQADiJetHistos &o);
    virtual ~AliFragFuncQADiJetHistos();
    
    virtual void DefineQADiJetHistos();
    virtual void FillDiJetQA(Double_t invMass, Double_t deltaPhi, Double_t deltaEta, Double_t deltaPt, Double_t jetBin);
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

    TH2F*   fh2InvMass;        // FF dijet invariant mass histos
    TH2F*   fh2DeltaPhi;       // FF dijet delta phi histos
    TH2F*   fh2DeltaEta;       // FF dijet delta eta histos
    TH2F*   fh2DeltaPt;        // FF dijet delta pt histos

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

  virtual void   SetBranchGenJets(const char* c){fBranchGenJets = c;}
  virtual void   SetBranchRecJets(const char* c){fBranchRecJets = c;}

  virtual void   SetTrackCuts(Float_t trackPt = 0.15, Float_t trackEtaMin = -0.9, Float_t trackEtaMax = 0.9, 
			      Float_t trackPhiMin = 0., Float_t trackPhiMax = 2*TMath::Pi())
  {fTrackPtCut = trackPt; fTrackEtaMin = trackEtaMin; fTrackEtaMax = trackEtaMax; 
    fTrackPhiMin = trackPhiMin; fTrackPhiMax = trackPhiMax;}
  virtual void   SetFilterMask(UInt_t i) {fFilterMask = i;}
  virtual void   UsePhysicsSelection(Bool_t b) {fUsePhysicsSelection = b;}
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

  static  void   SetProperties(TH1* h,const char* x, const char* y);
  static  void   SetProperties(TH2* h,const char* x, const char* y,const char* z);
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
  
  void   SetIJHistoBins(Int_t nJetPt = 245, Float_t jetPtMin = 5, Float_t jetPtMax = 250, Int_t nPt = 200, Float_t ptMin = 0., Float_t ptMax = 200., 
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
			   Float_t deltaEtaMax = 1.1, Int_t nDeltaPt = 100, Float_t deltaPtMin = 0., Float_t deltaPtMax = 100.)
  {
    fQADiJetNBinsInvMass = nInvMass; fQADiJetInvMassMin = invMassMin; fQADiJetInvMassMax = invMassMax; fQADiJetNBinsJetPt = nJetPt; fQADiJetJetPtMin = jetPtMin; 
    fQADiJetJetPtMax = jetPtMax; fQADiJetNBinsDeltaPhi = nDeltaPhi; fQADiJetDeltaPhiMin = deltaPhiMin; fQADiJetDeltaPhiMax = deltaPhiMax; fQADiJetNBinsDeltaEta = nDeltaEta; 
    fQADiJetDeltaEtaMin = deltaEtaMin; fQADiJetDeltaEtaMax = deltaEtaMax; fQADiJetNBinsDeltaPt = nDeltaPt; fQADiJetDeltaPtMin = deltaPtMin; fQADiJetDeltaPtMax = deltaPtMax;
  }

  Float_t  GetFFRadius() const { return fFFRadius; }
  void	   GetJetTracksTrackrefs(TList* l, const AliAODJet* j);
  void	   GetJetTracksPointing(TList* in, TList* out, const AliAODJet* j, const Double_t r, Double_t& pt);  
  Double_t GetDiJetBin(Double_t invMass, Double_t leadingJetPt, Double_t eMean, Int_t kindSlices); // function to find which bin fill
  Double_t InvMass(const AliAODJet* jet1, const AliAODJet* jet2);
  void     AssociateGenRec(TList* tracksAODMCCharged,TList* tracksRec, TArrayI& indexAODTr,TArrayI& indexMCTr,TArrayS& isGenPrim);
  void     FillSingleTrackRecEffHisto(THnSparse* histo, TList* tracksGen, const TList* tracksRec, const TArrayI& indexAODTr, const TArrayS& isGenPrim);
  void     FillJetTrackRecEffHisto(THnSparse* histo,Double_t jetPhi,Double_t jetEta,Double_t jetPtGen,Double_t jetPtRec, TList* jetTrackList, TList* tracksGen,
				   const TArrayI& indexAODTr,const TArrayS& isGenPrim);

    
  // Consts
  
  enum {kTrackUndef=0, kTrackAOD, kTrackAODQualityCuts, kTrackAODCuts, kTrackKineAll, kTrackKineCharged, kTrackKineChargedAcceptance, 
	kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance};
  enum {kJetsUndef=0, kJetsRec, kJetsRecAcceptance, kJetsGen, kJetsGenAcceptance, kJetsKine, kJetsKineAcceptance};
 
 
 private:
  
  Int_t   GetListOfTracks(TList* list, Int_t type);
  Int_t	  GetListOfJets(TList* list, Int_t type);
  
  AliESDEvent* fESD;      // ESD event
  AliAODEvent* fAOD;      // AOD event
  //AliMCEvent*  fMCEvent;  // MC event
  
  TString fBranchRecJets;         // branch name for reconstructed jets
  TString fBranchGenJets;         // branch name for generated jets
  
  Int_t   fTrackTypeGen;  // type of generated tracks
  Int_t   fJetTypeGen;    // type of generated jets

  Int_t   fJetTypeRecEff; // type of jets used for filling reconstruction efficiency histos

  UInt_t  fFilterMask;	  // filter bit for selected tracks
  Bool_t  fUsePhysicsSelection; // switch for event selection
	
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

  // Histograms
  TList	*fCommonHistList;         // List of common histos
  
  TH1F  *fh1EvtSelection;         //! event cuts 
  TH1F	*fh1VertexNContributors;  //! NContributors to prim vertex
  TH1F	*fh1VertexZ;              //! prim vertex z distribution
  TH1F	*fh1EvtMult;              //! number of reconstructed tracks after cuts 

  TProfile* fh1Xsec;              //! pythia cross section and trials
  TH1F*     fh1Trials;            //! sum of trials
  TH1F*     fh1PtHard;            //! pt hard of the event
  TH1F*     fh1PtHardTrials;      //! pt hard of the event

  TH1F  *fh1nRecJetsCuts;         //! number of jets from reconstructed tracks per event 
  TH1F  *fh1nGenJets;             //! number of jets from generated tracks per event
  TH1F  *fh1nRecEffJets;          //! number of jets for reconstruction eff per event

  // tracking efficiency 

  THnSparseF *fhnSingleTrackRecEffHisto; //! track reconstruction efficiency 
  THnSparseF *fhnJetTrackRecEffHisto;    //! reconstruction efficiency jet tracks 

  ClassDef(AliAnalysisTaskFragmentationFunction, 6);
};

#endif
