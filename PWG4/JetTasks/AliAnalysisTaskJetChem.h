// Analysis task for jet chemistry analysis
// Oliver Busch, o.busch@gsi.de
// UE task & CDF jet finder based on UE task by Arian Abrahantes Quintana and Enesto Lopez


#ifndef ALIANALYSISTASKJETCHEM_H
#define ALIANALYSISTASKJETCHEM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h" 
#include "AliPID.h"

//class AliAnalysisTaskSE; 
class AliESDEvent;
class AliAODEvent;
class TH1F;
class TH2F;
class TH1I;
class TProfile;
class TVector3;
class TTree;
class AliAODJet;
class AliAODTrack;
class AliPID;
class TDatabasePDG;

class  AliAnalysisTaskJetChem : public AliAnalysisTaskSE
  {
  public:
    AliAnalysisTaskJetChem(const char* name="AliAnalysisTaskJetChem");
    virtual           ~AliAnalysisTaskJetChem() {;}
    
    // Implementation of interface methods
    virtual     Bool_t UserNotify();
    virtual     void   UserCreateOutputObjects();
    virtual     void   UserExec(Option_t *option);
    virtual     void   Terminate(Option_t *);
    
    //  Setters
    virtual     void   SetDebugLevel( Int_t level )  { fDebug = level; }
 
    // Read deltaAODs
    void   ReadDeltaAOD()                       { fDeltaAOD = kTRUE; }
    void   SelectDeltaAODBranch(const char* val){ fDeltaAODBranch = val;   }
    void   SelectAODBranch(const char* val)     { fAODBranch = val;   }
    void   SelectDeltaAODBranchMC(const char* val){ fDeltaAODBranchMC = val;   }
    void   SelectAODBranchMC(const char* val)     { fAODBranchMC = val;   }

    void   SetJetsOnFly( Bool_t val )           { fJetsOnFly = val;  }

    // use internal jet finder 

    void   SetUseLOConeJets( )              { fUseLOConeJets = kTRUE; }
    void   SetUseLOConeMCJets()             { fUseLOConeMCJets = kTRUE;  fUsePythiaJets = kFALSE;}
    void   SetUsePythiaJets()               { fUseLOConeMCJets = kFALSE; fUsePythiaJets = kTRUE; }
    void   SetConeRadius( Double_t val )    { fConeRadius = val; }
    void   SetTrackPtCutJF( Double_t val )  { fTrackPtCutJF = val; }

    void   SetFilterBitJF( UInt_t val )     { fFilterBitJF = val; }
    void   SetRequireITSRefitJF()           { fRequireITSRefitJF = kTRUE; }
    void   SetRejectK0TracksJF()            { fRejectK0TracksJF  = kTRUE; }

    // Jet cuts
    void   SetJetPtCut( Double_t val )    { fJetPtCut = val;  }
    void   SetJetEtaCut( Double_t val )   { fJetEtaCut = val; }

    // track cuts
    void   SetFilterBit( UInt_t val )     { fFilterBit = val;  }
    void   SetTrackPtCut( Double_t val )  { fTrackPtCut = val; }
    void   SetTrackEtaCut( Double_t val ) { fTrackEtaCut = val; }

    // K0 cuts 

    void SetUseOnFlyV0s()                 { fUseOnFlyV0s = kTRUE; }
    void SetCutnSigdEdx( Double_t val)    { fCutnSigdEdx = val; }


    // Setters for MC
    void  SetUseAODMCTracksForUE(){ fUseAODMCTracksForUE = kTRUE;}

    
  private:
    AliAnalysisTaskJetChem(const  AliAnalysisTaskJetChem &det);
    AliAnalysisTaskJetChem&   operator=(const  AliAnalysisTaskJetChem &det);
    
    void   AnalyseEvent();
    Int_t  IsTrackInsideRegion(const AliAODJet* aodjetVect, const TVector3 *partVect);

    void   FillPIDhisto(TH1F* hist,Int_t pdg,Float_t weight=1);
    TH1F*  CreatePIDhisto(const char* name);
    TH1F*  CreatePythiaIDhisto(const char* name);
    void   FillPythiaIDhisto(TH1F* h, const Int_t PID);
    void   CreateHistos();
    void   FillSumPtRegion( Double_t leadingE, Double_t ptMax, Double_t ptMin );
    void   FillMultRegion(Double_t leadingE, Double_t nTrackPtmax, Double_t nTrackPtmin);
    TClonesArray* FindChargedParticleJets();
    TClonesArray* FindChargedParticleJetsMC();
    TClonesArray* GetPythiaJets();
    void   QSortTracks(TObjArray &a, Int_t first, Int_t last);
    void   WriteSettings();
    
    Bool_t IsK0InvMass(const Double_t mass) const;
    Bool_t IsLambdaInvMass(const Double_t mass) const;
    Bool_t IsAcceptedDCAK0(/*const Double_t dca*/) const;
    Bool_t IsAcceptedDCALambda(/*const Double_t dca*/) const;
    Bool_t IsAccepteddEdx(const Double_t mom, const Double_t dEdx, AliPID::EParticleType n, const Double_t cutnSig) const;
    void   CheckV0s(AliAODJet* jetVect, Int_t maxPtRegionIndex, Bool_t& foundK0);
    void   CheckMCParticles(AliAODJet* jetVect, Int_t maxPtRegionIndex,Bool_t& isK0event);
    Double_t AssociateV0MC(const TVector3* V0Mom,const Int_t pdg);
    void CompLeadingJets(AliAODJet* jetLeadingAOD,AliAODJet* jetLeadingMC,const Int_t pythiaPID,
			 const Bool_t foundK0AID,const Bool_t foundK0MC);
    Double_t GetJetRadius(const AliAODJet* jet, const Double_t energyFrac);
    Bool_t IsTrackFromK0(const Int_t indexTrackCheck);

    Bool_t IsQuarkHardScatteringEvent(const Int_t PID);
    Bool_t IsGluonHardScatteringEvent(const Int_t PID);
    Bool_t IsDiffractiveEvent(const Int_t PID);
    Int_t  GetPythiaProcessID();

    void GetJetTracksResum(TList* list, AliAODJet* jet, const Double_t radius);
    void GetJetTracksTrackrefs(TList* list, AliAODJet* jet);
    void FillReferenceFF(AliAODJet* jet);
    void FillReferencePlotsTracks(); 

    Int_t   fDebug;                 //  Debug flag

    Bool_t      fDeltaAOD;          //  Read jets from delta AOD 
    TString     fDeltaAODBranch;    //  Jet branch name from delta AOD
    TString     fAODBranch;         //  Jet branch name from standard AOD
    TString     fDeltaAODBranchMC;  //  MC Jet branch name from delta AOD
    TString     fAODBranchMC;       //  MC Jet branch name from standard AOD

    TClonesArray*  fArrayJetsAOD;   //  Array of AOD Jets from delta AOD
    TClonesArray*  fArrayJetsMC;    //  Array of MC  Jets from delta AOD
    

    AliAODEvent*  fAOD;             //! AOD Event 
    AliAODEvent*  fAODjets;         //! AOD Event for reconstructed on the fly (see ConnectInputData()
    TList*  fListOfHistos;          //  Output list of histograms
    Bool_t   fJetsOnFly;            // if jets are reconstructed on the fly from AOD tracks 

    Bool_t   fUseLOConeJets;        // Use LO cone finder instead of jets from AOD
    Bool_t   fUseLOConeMCJets;      // Use LO cone finder on aodmcparticles
    Bool_t   fUsePythiaJets;        // use pythia internal jet finder output

    Double_t fConeRadius;          // if selected Cone-like region type, set Radius (=0.7 default)
    Double_t fTrackPtCutJF;        // track lower pt for JF     
    UInt_t   fFilterBitJF;         // track filter for JetFinder 
    Bool_t   fRequireITSRefitJF;   // additional ITS refit requirement in JF 
    Bool_t   fRejectK0TracksJF;    // exclude tracks from K0 decay in JF 


    // Jet cuts 
    Double_t   fJetPtCut;          // Min Pt for charged Particle Jet
    Double_t   fJetEtaCut;         // |jet1 eta| < fJet1EtaCut   (fAnaType = 1,2,3)

    // track cuts
    UInt_t   fFilterBit;           // Select tracks from an specific track cut (default 0xFF all track selected)
    Double_t fTrackPtCut;          // Pt cut of tracks in the regions
    Double_t fTrackEtaCut;         // Eta cut on tracks in the regions (fRegionType=1)

    // K0 cuts 
    Bool_t fUseOnFlyV0s;           // on-the-fly V0s versus vertex track 'offline' V0s
    Double_t fCutnSigdEdx;         // TPC dEdx cut 
 

    // For MC
    Bool_t fUseAODMCTracksForUE;   // use aodmcparticles branch to determine max/min transverse region
    Double_t fAreaReg;             // Area of the region To be used as normalization factor when filling histograms
    Double_t fAvgTrials;           // average trials used to fill the fh1Triasl histogram in case we do not have trials on a event by event basis
    

    // Histograms    ( are owned by fListOfHistos TList )
    TH1F*  fhPrimVertexNCont;        //!
    TH1F*  fhPrimVertexRho;          //!
    TH1F*  fhPrimVertexZ;            //!
    TH1F*  fhNJets;                  //!
    TH1F*  fhNJetsMC;                //!
    TH1F*  fhLeadingEta;             //!
    TH2F*  fhLeadingNTracksVsEta;    //!
    TH2F*  fhLeadingPtVsEta;         //!
    TH1F*  fhLeadingPhi;             //!  
    TH1F*  fhLeadingPt;              //!
    TH1F*  fhLeadingPtDiffr;         //!
    TH1F*  fhLeadingEtaMC;           //!
    TH1F*  fhLeadingPhiMC;           //!  
    TH1F*  fhLeadingPtMC;            //!
    TH1F*  fhLeadingPtMCDiffr;       //!
    TH2F*  fhPhiEtaTracksNoCut;      //!
    TH1F*  fhPtTracksNoCut;          //!
    TH2F*  fhPhiEtaTracks;           //!
    TH1F*  fhPtTracks;               //!
    TH1F*  fhTrackMult;              //!

    TH1F*  fhEtaMCTracks;            //!
    TH1F*  fhPhiMCTracks;            //!
    TH1F*  fhPtMCTracks;             //!
    TH2F*  fhnTracksVsPtLeading;     //!
  
    TH1F*  fhdNdEtaPhiDist;          //!
   
    TH1F*  fhRegionSumPtMaxVsEt;     //!
    TH1F*  fhRegionMultMaxVsEt;      //!
    TH1F*  fhRegionSumPtMinVsEt;     //!
    TH1F*  fhRegionMultMinVsEt;      //!


    TH1F* fhNV0s;                     //!
    TH1F* fhV0onFly;                  //!
    TH1F* fhV0DCADaughters;           //!
    TH1F* fhV0Radius;                 //!
    TH1F* fhV0DCAToVertex;            //!
    TH1F* fhV0DCAToVertexK0;          //!
    
    TH1F* fhV0InvMassK0;               //!
    TH1F* fhV0InvMassK0JetEvt;         //!       
    TH1F* fhV0InvMassLambda;           //!
    TH1F* fhV0InvMassAntiLambda;       //!
    TH1F* fhV0InvMassLambdaJetEvt;     //!
    TH1F* fhV0InvMassAntiLambdaJetEvt; //!

    TH2F* fhdROpanK0VsPt;             //!
    TH1F* fhdPhiJetV0;                //!
    TH1F* fhdPhiJetK0;                //!
    TH1F* fhdRJetK0;                  //!

    TH1F* fhdNdptV0;                  //!
    TH1F* fhdNdptK0;                  //!
    TH2F* fhPtVsEtaK0;                //!

    TH1F* fhV0InvMassK0DCA;              //!
    TH1F* fhV0InvMassK0DCAdEdx;          //!
    TH1F* fhV0InvMassK0DCAPID;           //!
    TH1F* fhV0InvMassLambdaDCAdEdx;      //!
    TH1F* fhV0InvMassAntiLambdaDCAdEdx;  //!
    TH1F* fhdNdptK0DCA;                  //!
    TH1F* fhdNdptK0DCAdEdx;              //!

    TH1F* fhV0InvMassK0Min;          //!
    TH1F* fhV0InvMassLambdaMin;      //!
    TH1F* fhV0InvMassAntiLambdaMin;  //!

    TH1F* fhV0InvMassK0Max;          //!
    TH1F* fhV0InvMassLambdaMax;      //!
    TH1F* fhV0InvMassAntiLambdaMax;  //!

    TH1F* fhV0InvMassK0Jet;          //!
    TH1F* fhV0InvMassLambdaJet;      //!
    TH1F* fhV0InvMassAntiLambdaJet;  //!

    TH1F* fhV0InvMassK0Lambda;       //!

    TH1F* fhdNdptK0JetEvt;           //!

    TH1F* fhdNdzK0;                   //!
    TH1F* fhdNdzK05to10;              //!
    TH1F* fhdNdzK010to20;             //!
    TH1F* fhdNdzK020to30;             //!
    TH1F* fhdNdzK030to40;             //!
    TH1F* fhdNdzK040to60;             //!

    TH1F* fhdNdxiK0;                 //!

    TH1F* fhdNdzLambda;              //!
    TH1F* fhdNdzAntiLambda;          //!

    TH1F* fhdNdzK0Max;               //!
    TH1F* fhdNdxiK0Max;              //!
    TH1F* fhdNdzLambdaMax;           //!
    TH1F* fhdNdxiLambdaMax;          //!
    
    TH1F* fhdNdptK0Max;              //!
    TH1F* fhdNdptLambdaMax;          //!
    
    TH1F* fhdNdzK0Min;               //!
    TH1F* fhdNdxiK0Min;              //!
    TH1F* fhdNdzLambdaMin;           //!
    TH1F* fhdNdxiLambdaMin;          //!
    
    TH1F* fhdNdptK0Min;              //!
    TH1F* fhdNdptLambdaMin;          //!

    TH1F* fhdNdzK0Jet;               //!
    TH1F* fhdNdxiK0Jet;              //!
    TH1F* fhdNdzLambdaJet;           //!
    TH1F* fhdNdxiLambdaJet;          //!

    TH1F* fhdNdptK0Jet;              //!
    TH1F* fhdNdptLambdaJet;          //!
    
    TH2F* fhdEdxVsMomV0;             //!
    TH2F* fhdEdxVsMomV0pidEdx;       //! 
    TH2F* fhdEdxVsMomV0piPID;        //! 

    TH1F* fhdPhiJetK0MC;             //!
    TH1F* fhdRJetK0MC;               //!
    TH1F* fhdRV0MC;                  //!

    TH1F* fhdNdptchPiMCMax;               //!
    TH1F* fhdNdptK0MCMax;                 //!
    TH1F* fhdNdptchKMCMax;                //!
    TH1F* fhdNdptpMCMax;                  //!
    TH1F* fhdNdptpBarMCMax;               //!
    TH1F* fhdNdptLambdaMCMax;             //!
    TH1F* fhdNdptLambdaBarMCMax;          //!

    TH1F* fhdNdptchPiMCMin;               //!
    TH1F* fhdNdptK0MCMin;                 //!
    TH1F* fhdNdptchKMCMin;                //!
    TH1F* fhdNdptpMCMin;                  //!
    TH1F* fhdNdptpBarMCMin;               //!
    TH1F* fhdNdptLambdaMCMin;             //!
    TH1F* fhdNdptLambdaBarMCMin;          //!
    TH1F* fhdNdptOmegaMCMin;              //!
    TH1F* fhdNdptOmegaBarMCMin;           //!

    TH1F* fhdNdptchPiMCJet;               //!
    TH1F* fhdNdptK0MCJet;                 //!
    TH1F* fhdNdptchKMCJet;                //!
    TH1F* fhdNdptpMCJet;                  //!
    TH1F* fhdNdptpBarMCJet;               //!
    TH1F* fhdNdptLambdaMCJet;             //!
    TH1F* fhdNdptLambdaBarMCJet;          //!


    // kine tree 
    TH1F* fhPIDMC;                    //!
    TH1F* fhPIDMC_quarkEv;            //!
    TH1F* fhPIDMC_gluonEv;            //!
    TH1F* fhPIDMCAll;                 //!
    TH1F* fhPIDMCMin;                 //!
    TH1F* fhPIDMCJet;                 //!
 
    TH1F* fhPIDMCMotherK0;            //!
    TH1F* fhPIDMCGrandMotherK0;       //!
    TH1F* fhPIDMCMotherChK;           //!
    TH1F* fhPIDMCMotherK0Trans;       //!
    TH1F* fhPIDMCGrandMotherK0Trans;  //!
    TH1F* fhPIDMCMotherChKTrans;      //!

    TH1F* fhdNdptgammaMC;             //!
    TH1F* fhdNdptchPiMC;              //!
    TH1F* fhdNdptpi0MC;               //!
    TH1F* fhdNdptK0MC;                //!
    TH1F* fhdNdptchKMC;               //!
    TH1F* fhdNdptpMC;                 //!
    TH1F* fhdNdptpBarMC;              //!
    TH1F* fhdNdptLambdaMC;            //!
    TH1F* fhdNdptLambdaBarMC;         //!
    TH1F* fhdNdptOmegaMC;             //!
    TH1F* fhdNdptOmegaBarMC;          //!

    TH1F* fhdNdxiMC;                 //!
    TH1F* fhdNdxiK0MC;               //!
    TH1F* fhdNdxiK0MCJet;            //!

    TH1F* fhdNdzK0MC;                //!
    TH1F* fhdNdzK0MCJet;             //!
    TH1F* fhdNdptK0MCJetEvt;         //!
      
    TH2F* fhnJetsAODvsMC;             //!
    TH2F* fhLeadingPtAODvsMC;         //!
    TH2F* fhLeadingEtaAODvsMC;        //!
    TH2F* fhLeadingPhiAODvsMC;        //!
    TH2F* fhnTracksLeadingAODvsMC;    //!
    
    TH1F* fhLeadingdRAODMC;            //!
    TH2F* fhLeadingPtAODvsMCdRcut;     //!
    TH2F* fhdnTracksVsdPtLeadingAODMC; //!

    TH2F* fhnTracksJetVsPtAOD;         //!
    TH2F* fhnTracksJetVsPtAODquarkEv;  //!
    TH2F* fhRadiusJetVsPtAOD;          //!
    TH2F* fhnTracksJetVsPtMC;          //!
    TH2F* fhnTracksJetVsPtMCquarkEv;   //!
    TH2F* fhRadiusJetVsPtMC;           //!

    TH2F* fhnTracksJetVsPtMCK0;         //!
    TH2F* fhnTracksJetVsPtMCK0quarkEv;  //!
    TH2F* fhRadiusJetVsPtMCK0;          //!

    TH2F* fhnTracksJetVsPtAODK0;         //!
    TH2F* fhnTracksJetVsPtAODK0quarkEv;  //!
    TH2F* fhRadiusJetVsPtAODK0;          //!
    TH2F* fhnTracksJetVsPtAODpKch;       //!
    TH2F* fhRadiusJetVsPtAODpKch;        //!

    TH1F* fhPythiaProcess;           //!
    TH1F* fhPythiaProcessK0;         //!
    TH1F* fhPythiaProcessKch;        //!
    TH1F* fhPythiaProcessp;          //!
    TH1F* fhPythiaProcesspbar;       //!

    TH1F* fhdNdzJets5to10;  //!
    TH1F* fhdNdzJets10to20; //!
    TH1F* fhdNdzJets20to30; //!
    TH1F* fhdNdzJets30to40; //!
    TH1F* fhdNdzJets40to60; //!

    TH1F* fhdNdxiJets5to10;  //!
    TH1F* fhdNdxiJets10to20; //!
    TH1F* fhdNdxiJets20to30; //!
    TH1F* fhdNdxiJets30to40; //!
    TH1F* fhdNdxiJets40to60; //!

    TH1F* fhdNdptTracksJetPt5to10;  //!
    TH1F* fhdNdptTracksJetPt10to20; //!
    TH1F* fhdNdptTracksJetPt20to30; //!
    TH1F* fhdNdptTracksJetPt30to40; //!
    TH1F* fhdNdptTracksJetPt40to60; //!


    TProfile* fh1Xsec;               //!
    TH1F*  fh1Trials;                //!
    
    //TTree* fSettingsTree;            //! Fast Settings saving
    
    TDatabasePDG* fpdgdb;            //!

    enum PythiaPIDHistoBin{kPythiaPIDP11Bin=1, kPythiaPIDP12Bin=3, kPythiaPIDP13Bin=5, kPythiaPIDP28Bin=7, 
			   kPythiaPIDP53Bin=9, kPythiaPIDP68Bin=11, kPythiaPIDP92Bin=13, kPythiaPIDP93Bin=15, 
			   kPythiaPIDP94Bin=17,kPythiaPIDP95Bin=19, kPythiaPIDPOtherBin=21}; 


    enum PIDHistoBin{kPDGpm311Bin=48,kPDG333Bin=49,kPDGpm313Bin=50,kPDGp323Bin=51,kPDGm323Bin=52,
		     kPDGNeutrinoBin=53,kPDGCharmedBaryonBin=54,kPDGQuarkBin=55,kPDGDiQuarkBin=56};
   
    ClassDef( AliAnalysisTaskJetChem, 1); // Analysis task for jet chemistry analysis 
  };

#endif

    
