//-*- Mode: C++ -*-
#ifndef ALIANALYSEUE_H
#define ALIANALYSEUE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class  for transverse regions analysis
//---------------------------------------------
////////////////////////////////////////////////

// --- ROOT system ---
#include <TObject.h> 

#if __GNUC__ >= 3
using namespace std;
#endif

class AliAnalysisTaskUE;
class AliAODEvent;
class AliAODTrack;
class AliGenPythiaEventHeader;
class AliHistogramsUE;
class AliMCEvent;
class TVector3;

class AliAnalyseUE : public TObject {

 public: 

  AliAnalyseUE();                                         //constructor
  AliAnalyseUE(const AliAnalyseUE & g);                   //copy constructor
  AliAnalyseUE & operator = (const AliAnalyseUE & g);     //assignment operator
  virtual ~AliAnalyseUE();                                //virtual destructor

  void 		AnalyseMC(TVector3 *jetVect, AliMCEvent *mcEvent, AliGenPythiaEventHeader  *pythiaGenHeader, Int_t conePosition, Bool_t useAliStack, Bool_t constrainDistance, Double_t minDistance);
  Bool_t   	AnaTypeSelection(TVector3 *jetVect);

  TList*	CreateHistograms(Int_t bins, Double_t min, Double_t max);

  void          FillLeadingJet( Double_t  w);
 
  void          FillRegions(Bool_t isNorm2Area, TVector3 *jetVect);

  void          FillTrials(const char *namex, Double_t  w);
  
  void          FillVertex(Double_t  w);
  
  void          FillXsec(const char *namex, Double_t  w);

  void 		FindMaxMinRegions(TVector3 *jetVect, Int_t conePosition);
  
  TList*        GetHistograms();
  
  TVector3	GetOrderedClusters(TString aodBranch, Bool_t chargedJets, Double_t chJetPtMin);

  void          Initialize(AliAnalysisTaskUE& tmp);
  void          Initialize(Int_t anaType, AliAODEvent* aod,Double_t coneRadius, Int_t debug, Int_t filterBit, Double_t 
  jet1EtaCut, Double_t jet2DeltaPhiCut, Double_t jet2RatioPtCut, Double_t jet3PtCut, Int_t ordering, Int_t regionType,Bool_t simulateChJetPt, Double_t trackEtaCut, Double_t trackPtCut, Bool_t useChargeHadrons, Bool_t useChPartJet, Bool_t usePositiveCharge, Bool_t useSingleCharge);

  Bool_t        VertexSelection(AliAODEvent *value, Int_t tracks, Double_t zed);
  void		WriteSettings();

  // Various setters when you do not want to initialize members from AliAnalysisTaskUE
  void          SetAnaTopology(Int_t value)		{ fAnaType = value; }
  void          SetAOD(AliAODEvent *value)		{ fkAOD = value; }
  void          SetConeRadius(Double_t value)		{ fConeRadius = value; }
  void          SetDebug(Int_t value)			{ fDebug = value; }
  void          SetFilterBit(Int_t value)		{ fFilterBit = value; }
  void          SetJet1EtaCut(Double_t value)		{ fJet1EtaCut = value; }
  void          SetJet2DeltaPhiCut(Double_t value)	{ fJet2DeltaPhiCut = value; }
  void          SetJet2RatioPtCut(Double_t value)	{ fJet2RatioPtCut = value; }
  void          SetJet3PtCut(Double_t value)		{ fJet3PtCut = value; }
  void          SetOrdering(Int_t value)		{ fOrdering = value; }
  void          SetRegionType(Int_t value)		{ fRegionType = value; }
  void          SetSimulateChJetPt(Bool_t value)	{ fSimulateChJetPt = value; }
  void          SetTrackEtaCut(Double_t value)		{ fTrackEtaCut = value; }
  void          SetTrackPtCut(Double_t value)		{ fTrackPtCut = value; }
  void          SetUseChargeHadrons(Bool_t value)	{ fUseChargeHadrons = value; }
  void          SetUseChPartJet(Bool_t value)	        { fUseChPartJet = value; }
  void          SetUsePositiveCharge(Bool_t value)	{ fUsePositiveCharge = value; }
  void          SetUseSingleCharge(Bool_t value)	{ fUseSingleCharge = value; }

 private:

  void          FillAvePartPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin );
  void          FillMultRegion( Double_t leadingE, Double_t nTrackPtmax, Double_t nTrackPtmin, Double_t ptMin );
  void          FillSumPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin );
  TObjArray*  	FindChargedParticleJets( Double_t chJetPtMin);
  Int_t 	IsTrackInsideRegion(TVector3 *jetVect, TVector3 *partVect, Int_t conePosition);
  void  	QSortTracks(TObjArray &a, Int_t first, Int_t last);
  void  	SetRegionArea(TVector3 *jetVect);
  TObjArray*    SortChargedParticles();     
  virtual const Bool_t	TrackSelected(AliAODTrack* part);
  virtual const Bool_t	TrackMCSelected(Double_t charge, Double_t pT, Double_t eta, Int_t pdgCode);

  
    //AliAnalysisTaskUE    fTaskUE;        //  current instance of the analysis-task
    const AliAODEvent*   fkAOD;             //! AOD Event 
    Int_t          fDebug;           //  Debug flag

    
    // For MC
    Bool_t         fSimulateChJetPt;      // Naive simulation of charged jet Pt from original Jet in MC Header
    
    // Cuts UE analysis
    Int_t          fAnaType;              // Analysis type on jet topology: 
    Double_t       fAreaReg;              // Area of the region To be used as normalization factor when filling histograms
    Double_t       fConeRadius;           // if selected Cone-like region type, set Radius (=0.7 default)
    UInt_t         fFilterBit;            // Select tracks from an specific track cut (default 0xFF all track selected)
    Int_t	   fRegionType;           // 1 = transverse regions (default)
                                          // 2 = cone regions   
    Bool_t         fUseChargeHadrons;     // Only use charge hadrons
    Bool_t         fUseChPartJet;         // Use "Charged Particle Jet" instead of jets from AOD see FindChargedParticleJets()

    Bool_t         fUsePositiveCharge;    //If Single type of charge used then set which one (=kTRUE default positive)
    Bool_t         fUseSingleCharge;      //Make analysis for a single type of charge (=kFALSE default)
    
    Int_t          fOrdering;             //  Pt and multiplicity summation ordering:
    
    // Jet cuts 
    Double_t      fJet1EtaCut;       // |jet1 eta| < fJet1EtaCut   (fAnaType = 1,2,3)
    Double_t      fJet2DeltaPhiCut;  // |Jet1.Phi - Jet2.Phi| < fJet2DeltaPhiCut (fAnaType = 2,3)
    Double_t      fJet2RatioPtCut;   // Jet2.Pt/Jet1Pt > fJet2RatioPtCut  (fAnaType = 2,3)
    Double_t      fJet3PtCut;        // Jet3.Pt < fJet3PtCut  (fAnaType = 3)

    // track cuts
    Double_t      fTrackEtaCut;      // Eta cut on tracks in the regions (fRegionType=1)
    Double_t      fTrackPtCut;       // Pt cut of tracks in the regions
  
    AliHistogramsUE* fHistos;	     // Pointer to histogram class	

    //to fill the different regions
    Double_t	  fSumPtRegionPosit;	// Sum pT in positive region
    Double_t      fSumPtRegionNegat;	// Sum pT in negative region
    Double_t      fSumPtRegionForward;	// Sum pT in forward region
    Double_t      fSumPtRegionBackward;	// Sum pT in backward region
    Double_t      fMaxPartPtRegion;	// Max part pt in region
    Int_t         fNTrackRegionPosit;	// Tracks in positive region
    Int_t         fNTrackRegionNegat;	// Tracks in negative region 
    Int_t         fNTrackRegionForward;	// Track in forward region
    Int_t         fNTrackRegionBackward;// Tracks in backward region

    //Store analysis settings
    TTree* 	  fSettingsTree;	// To store analysis settings
    
    ClassDef(AliAnalyseUE,1)
};
#endif
