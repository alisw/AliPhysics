#ifndef ALIANALYSISTASKUE_H
#define ALIANALYSISTASKUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTask.h"

class AliESDEvent;
class AliAODEvent;
class TH1F;
class TH2F;
class TH1I;
class TProfile;
class TVector3;
class TTree;

class  AliAnalysisTaskUE : public AliAnalysisTask
  {
  public:
    AliAnalysisTaskUE(const char* name="AliAnalysisTaskUE");
    virtual           ~AliAnalysisTaskUE() {if ( fListOfHistos ) delete fListOfHistos; }
    
    // Implementation of interace methods
    virtual     Bool_t Notify();
    virtual     void   ConnectInputData(Option_t *);
    virtual     void   CreateOutputObjects();
    virtual     void   Exec(Option_t *option);
    virtual     void   Terminate(Option_t *);

    //Select the trigger
    void SelectTrigger(Int_t trig) { fTrigger = trig; }

    //  Setters
    virtual     void   SetDebugLevel( Int_t level )  { fDebug = level; }
    void   SetPtRangeInHist( Int_t bin, Double_t min, Double_t max ) {
      fBinsPtInHist = bin; 
      fMinJetPtInHist = min; 
      fMaxJetPtInHist = max; 
    }

    // Read deltaAODs
    void   ReadDeltaAOD()                   { fDeltaAOD = kTRUE; }
    void   SelectDeltaAODBranch(const char* val)     { fDeltaAODBranch = val;   }
    void   SelectAODBranch(const char* val)     { fAODBranch = val;   }

    // Setters for MC
    void  SetUseMCBranch(){fUseMCParticleBranch = kTRUE;}
    void  SetConstrainDistance(Bool_t val1, Double_t val2){ fMinDistance = val2; fConstrainDistance = val1;}
    void  SetSimulateChJetPt(){fSimulateChJetPt = kTRUE;}
    void  SetUseAODMCParticle(){fUseAliStack = kFALSE;}
    
    //Setters for Events QA
    void  SetZVertex( Double_t val ) { fZVertex = val; }
    void  SetTracksInVertex( Int_t val ){ fnTracksVertex = val; }
    
    // Stters for UE Analysis
    void   SetAnaTopology( Int_t val )    { fAnaType = val;    }
    void   SetRegionType( Int_t val )     { fRegionType = val; }
    void   SetUseChPartJet( Int_t val )   { fUseChPartJet = val; }
    void   SetUseChargeHadrons( Bool_t val ){ fUseChargeHadrons = val; }
    void   SetPtSumOrdering( Int_t val ) { fOrdering = val;   }
    void   SetFilterBit( UInt_t val )     { fFilterBit = val;  }
    void   SetJetsOnFly( Bool_t val )     { fJetsOnFly = val;  }
    void   SetConeRadius( Double_t val )  { fConeRadius = val; }
    void   SetConePosition(Int_t val)     { fConePosition= val; }
    void   SetUseSingleCharge()  { fUseSingleCharge = kTRUE; } 
    void   SetUseNegativeChargeType()        { fUsePositiveCharge = kFALSE; }
    void   SetDoNotNormalizeQuantities()  { fIsNorm2Area = kFALSE; }
    // Jet cuts
    void   SetPtMinChPartJet( Double_t val )  { fChJetPtMin = val; }
    void   SetJet1EtaCut( Double_t val )      { fJet1EtaCut = val; }
    void   SetJet2DeltaPhiCut( Double_t val ) { fJet2DeltaPhiCut = val; }
    void   SetJet2RatioPtCut( Double_t val )  { fJet2RatioPtCut = val; }
    void   SetJet3PtCut( Double_t val )       { fJet3PtCut = val; }
    // track cuts
    void   SetTrackPtCut( Double_t val )  { fTrackPtCut = val; }
    void   SetTrackEtaCut( Double_t val ) { fTrackEtaCut = val; }
    
  private:
    AliAnalysisTaskUE(const  AliAnalysisTaskUE &det);
    AliAnalysisTaskUE&   operator=(const  AliAnalysisTaskUE &det);
    void   AnalyseUE();
    Int_t   IsTrackInsideRegion(TVector3 *jetVect, TVector3 *partVect);
    void   CreateHistos();
    void   SetRegionArea(TVector3 *jetVect);
    void   FillSumPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin );
    void   FillAvePartPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin );
    void   FillMultRegion( Double_t leadingE, Double_t nTrackPtmax, Double_t nTrackPtmin, Double_t ptMin );
    TObjArray*  FindChargedParticleJets();
    TObjArray*  SortChargedParticles();
    void   QSortTracks(TObjArray &a, Int_t first, Int_t last);
    void   WriteSettings();
    
    Int_t      fTrigger;         //Trigger flag as defined in AliAnalysisHelperJetTasks.h
    Int_t      fDebug;           //  Debug flag
    Bool_t      fDeltaAOD;        //  Read jets from delta AOD 
    TString     fDeltaAODBranch;  //  Jet branch name from delta AOD
    TString     fAODBranch;       //  Jet branch name from standard AOD
    TClonesArray*  fArrayJets;       //!  Array of Jets from delta AOD

    AliAODEvent*  fAOD;             //! AOD Event 
    AliAODEvent*  fAODjets;         //! AOD Event for reconstructed on the fly (see ConnectInputData()
    TList*  fListOfHistos;    //  Output list of histograms
    
    // Config
    Int_t   fBinsPtInHist;     //  # bins for Pt histos range
    Double_t   fMinJetPtInHist;   //  min Jet Pt value for histo range
    Double_t   fMaxJetPtInHist;   //  max Jet Pt value for histo range
    Bool_t     fIsNorm2Area;      // Apply Area Normalization to collected observables
    
    // For MC
    Bool_t fUseMCParticleBranch;  // Run Over mcparticles branch in AOD
    Bool_t fConstrainDistance;    // Constrain Distance between rec jet and pyth
    Double_t fMinDistance;  // Minimum distance between rec jet and pyth
    Bool_t fSimulateChJetPt; // Naive simulation of charged jet Pt from original Jet in MC Header
    Bool_t fUseAliStack;     // Use AliSatck for particle info otherwise "mcparticles" branch in AOD
    
    // Cuts Events type
    Int_t fnTracksVertex;  // QA tracks pointing to principal vertex (= 3 default) 
    Double_t fZVertex;     // Position of Vertex in Z direction
    
    // Cuts 
    Int_t   fAnaType;          // Analysis type on jet topology: 
    //     1=inclusive  (default) 
    //     2=back to back inclusive
    //     3=back to back exclusive
    //     4=Pt max (max Pt track in region)
    //     5=gama jet (back to back) ???
    //  Minimum bias
    //     31 = Semi jet (charged leading particle jets)
    //     32 = Random jetcone  ?
    //     33 = Swiss chees   ?
    
    // UE analysis is conducted in different type of regions
    // Transverse are those like defined in: R. Field Acta Physica Polonica B. Vol 36 No. 2 pg 167 (2005) 
    // Cone regions like defined in: Phys. Rev. D 70, 072002 (2004)
    Int_t   fRegionType;       // 1 = transverse regions (default)
                               // 2 = cone regions   
    Double_t   fConeRadius;       // if selected Cone-like region type, set Radius (=0.7 default)
    Int_t   fConePosition;     // This parameter set how will cone center in transversal zone will be set
                               //    1 : To be used in any jet topology (default value)
                               //        eta_cone = eta_leadingjet
                               //        phi_cone = phi_leadingjet + - 90
                               //    2 : To be used in multiple jet topology (code will cry otherwise)
                               //        eta_cone = (eta_leadingjet + eta_subleadingjet)/2
                               //        phi_cone = phi_leadingjet + - 90
    Double_t   fAreaReg;       // Area of the region To be used as normalization factor when filling histograms
                               // if fRegionType = 2 not always it is included within eta range
    Bool_t   fUseChPartJet;     // Use "Charged Particle Jet" instead of jets from AOD see FindChargedParticleJets()
    Bool_t   fUseChPartMaxPt;   // Use "Charged Particle with max Pt" instead of any jets to define forward region

    Bool_t   fUseChargeHadrons;   // Only use charge hadrons
    
    // Theoreticians ask for tools charge-aware
    // especially those which are related to multiplicity and MC-tunings
    // see arXiv:hep-ph/0507008v3
    Bool_t   fUseSingleCharge;     //Make analysis for a single type of charge (=kFALSE default)
    Bool_t   fUsePositiveCharge;   //If Single type of charge used then set which one (=kTRUE default positive)
    Int_t   fOrdering;         //  Pt and multiplicity summation ordering:
    //     1=CDF-like -independent sorting according quantity to be scored: Double sorting- (default)
    //       if Pt summation will be scored take Pt minimum between both zones and 
    //          fill Pt Max. and Min. histog. accordingly
    //       if Multiplicity summation will be scored take Mult. minimum between both zones and 
    //          fill Mult Max and Min histog. accordingly
    //       Bib:
    //     2=Marchesini-like (Only Pt sorting: Single sorting)
    //          sort only according Pt summation scored, find minimum between both zones and
    //          fill Pt and Multiplicity Max and Min summation histog. following only this criterium
    //       Bib: Phys. Rev. D 38, 3419 (1988)
    //     3=Nameless pt per track single sorting
    //          sort according to pt per track scored in each transverse zone 
    //          lowest values indicates minimum zone.   
    //     4=User Selection sorting (NOTE: USER must implement it within cxx)
    
    UInt_t   fFilterBit;        // Select tracks from an specific track cut (default 0xFF all track selected)
    Bool_t   fJetsOnFly;        // if jets are reconstructed on the fly from AOD tracks (see ConnectInputData() )
    
    // Jet cuts 
    Double_t   fChJetPtMin;       // Min Pt for charged Particle Jet
    Double_t   fJet1EtaCut;       // |jet1 eta| < fJet1EtaCut   (fAnaType = 1,2,3)
    Double_t   fJet2DeltaPhiCut;  // |Jet1.Phi - Jet2.Phi| < fJet2DeltaPhiCut (fAnaType = 2,3)
    Double_t   fJet2RatioPtCut;   // Jet2.Pt/Jet1Pt > fJet2RatioPtCut  (fAnaType = 2,3)
    Double_t   fJet3PtCut;        // Jet3.Pt < fJet3PtCut  (fAnaType = 3)
    // track cuts
    Double_t   fTrackPtCut;       // Pt cut of tracks in the regions
    Double_t   fTrackEtaCut;      // Eta cut on tracks in the regions (fRegionType=1)
    Double_t   fAvgTrials;        // average trials used to fill the fh1Triasl histogram in case we do not have trials on a event by event basis
    	  	 
    // Histograms    ( are owned by fListOfHistos TList )
    TH1F*  fhNJets;                  //!
    TH1F*  fhEleadingPt;             //!
    
    TH1F*  fhMinRegPtDist;           //!
    TH1F*  fhRegionMultMin;          //!
    TH1F*  fhMinRegAvePt;            //!
    TH1F*  fhMinRegSumPt;            //!
    TH1F*  fhMinRegMaxPtPart;        //!
    TH1F*  fhMinRegSumPtvsMult;      //!
    
    TH2F*  fhdNdEtaPhiDist;          //!
    TH2F*  fhFullRegPartPtDistVsEt;  //!
    TH2F*  fhTransRegPartPtDistVsEt; //!
    
    TH1F*  fhRegionSumPtMaxVsEt;     //!
    TH1I*  fhRegionMultMax;          //!
    TH1F*  fhRegionMultMaxVsEt;      //!
    TH1F*  fhRegionSumPtMinVsEt;     //!
    TH1F*  fhRegionMultMinVsEt;      //!
    TH1F*  fhRegionAveSumPtVsEt;     //!
    TH1F*  fhRegionDiffSumPtVsEt;    //!
    
    TH1F*  fhRegionAvePartPtMaxVsEt; //!
    TH1F*  fhRegionAvePartPtMinVsEt; //!
    TH1F*  fhRegionMaxPartPtMaxVsEt; //!
    
    //TH1F*  fhRegForwardSumPtVsEt;    //!
    //TH1F*  fhRegForwardMultVsEt;     //!
    //TH1F*  fhRegBackwardSumPtVsEt;   //!
    //TH1F*  fhRegBackwardMultVsEt;    //!
    TH2F*  fhRegForwardMult;         //!
    TH2F*  fhRegForwardSumPtvsMult;  //!
    TH2F*  fhRegBackwardMult;        //!
    TH2F*  fhRegBackwardSumPtvsMult; //!
    TH2F*  fhRegForwardPartPtDistVsEt; //!
    TH2F*  fhRegBackwardPartPtDistVsEt; //!
    TH2F*  fhRegTransMult;         //!
    TH2F*  fhRegTransSumPtVsMult;    //!
    TH2F*  fhMinRegSumPtJetPtBin;    //!
    TH2F*  fhMaxRegSumPtJetPtBin;    //!
    TH1F*  fhVertexMult;

    //        TH2F*  fhValidRegion; //! test to be canceled
    
    TProfile* fh1Xsec;               //!
    TH1F*  fh1Trials;                //!
    
    TTree* fSettingsTree;            //! Fast Settings saving

    ClassDef( AliAnalysisTaskUE, 4); // Analysis task for Underlying Event analysis
  };

#endif

    
