/// \class AliAnalysisTaskMKBase
/// \brief Base class providing functionality for derived AnalaysisTasks
///
/// This class does some basic event analysis and provides a framework and 
/// functionality to derived analysis tasks. The idea is not duplicate code
/// which is needed by many tasks and keep the derived tasks very simple and 
/// clean.
/// 
/// In addition, this class provides (static) functionality to easily create and fill THnSparseD.
///
/// \author Michael Linus Knichel <michael.linus.knichel@cern.ch>, CERN
/// \date Mar 8, 2019

#ifndef AliAnalysisTaskMKBase_H
#define AliAnalysisTaskMKBase_H

#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"
#include "AlidNdPtTools.h"
#include "AliEventCuts.h"

class AliESDtrackCuts;
class AliInputEventHandler;
class AliAnalysisManager;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class AliESDtrack;
class AliEventCuts;
class AliMCParticle;

class AliAnalysisTaskMKBase : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskMKBase();
                                AliAnalysisTaskMKBase(const char *name);
        virtual                 ~AliAnalysisTaskMKBase();

        virtual void           UserCreateOutputObjects();
        virtual void           UserExec(Option_t* option);
        virtual void           Terminate(Option_t* option);
        virtual void           FinishTaskOutput();        
        
        virtual void           SetESDtrackCutsM(AliESDtrackCuts* cut) { fESDtrackCutsM = cut; }
        virtual Bool_t         SetESDtrackCuts(Int_t i, AliESDtrackCuts* cut) { if (i<10) fESDtrackCuts[i] = cut; return (i<10); }        
        virtual void           SetTriggerMaskRequired(UInt_t trigger) { fTriggerMaskRequired = trigger; }
        virtual void           SetTriggerMaskRejected(UInt_t trigger) { fTriggerMaskRejected = trigger; }        
        virtual void           SetUseEventCuts(Bool_t use = kTRUE) { fUseEventCuts = use; }
        
        AliESDtrackCuts*       GetESDtrackCutsM() { return fESDtrackCutsM; }
        AliESDtrackCuts*       GetESDtrackCuts(Int_t i) { return (i < 10) ? fESDtrackCuts[i] : 0; }
        
        static Long64_t        FillHist(THnSparseD* s, Double_t x1, Double_t x2=0, Double_t x3=0, Double_t x4=0, Double_t x5=0, Double_t x6=0, Double_t x7 =0, Double_t x8 =0, Double_t x9 =0, Double_t x10 =0, Double_t x11 =0, Double_t x12 =0) { return AlidNdPtTools::FillHist(s, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12); }
        static Int_t           AddAxis(const char* label, Int_t nbins, Double_t xmin, Double_t xmax, const char* option = 0) { return AlidNdPtTools::AddAxis(label, nbins, xmin, xmax, option); }
        static Int_t           AddAxis(const char* label, const char* title, Int_t nbins, Double_t xmin, Double_t xmax, const char* option = 0) { return AlidNdPtTools::AddAxis(label, title, nbins, xmin, xmax, option); }
        static Int_t           AddAxis(const char* label, Int_t nbins, Double_t* xbins, const char* option = 0) { return AlidNdPtTools::AddAxis(label, nbins, xbins, option); }
        static Int_t           AddAxis(const char* label, const char* title, Int_t nbins, Double_t* xbins, const char* option = 0) { return AlidNdPtTools::AddAxis(label, title, nbins, xbins, option = 0); }
        static Int_t           AddAxis(const char* label, const char* title, const char* option) { return AlidNdPtTools::AddAxis(label, title, option); }
        static Int_t           AddAxis(const char* label, const char* option) { return AlidNdPtTools::AddAxis(label, option); }
        static Int_t           AddAxis(const char* option) { return AlidNdPtTools::AddAxis(option); }                 
        static THnSparseD*     CreateHist(const char* name) { return AlidNdPtTools::CreateHist(name); }
        static void            ResetHist() { AlidNdPtTools::ResetHist(); }
        static TH1D*           CreateLogHist(const char* name, const char* title) { return AlidNdPtTools::CreateLogHist(name, title); }
        static TH1D*           CreateLogHist(const char* name) { return AlidNdPtTools::CreateLogHist(name); }
        
        static AliAnalysisTaskMKBase* AddTaskMKBase(const char* name = "TaskMKBase", const char* outfile = 0);

    protected:
        
        virtual void          Log(const char* name) { Log(fLogHist,name); }        
        virtual void          Err(const char* name) { Log(fLogErr,name); }        
        virtual void          LogEvent(const char* name) { Log(fLogEvent,name); }
        
        virtual void          Log(const char* name, Int_t n)    { Log(fLogHist,name,n); }
        virtual void          Log(const char* name, Double_t n) { Log(fLogHist,name,n); }
        
        virtual void          Log(TH1D* h, const char* name) { if (h) h->Fill(name,1); }
        virtual void          Log(TH1D* h, const char* name, Int_t n)    { TString s(name); s+=n; Log(h,s.Data()); }
        virtual void          Log(TH1D* h, const char* name, Double_t n) { TString s(name); s+=n; Log(h,s.Data()); }
        
        virtual Bool_t        ReadEvent(); // read the event info
        virtual Bool_t        ReadMCEvent(); // read the mc event info
 
        virtual void          AddOutput() {}; // add histograms to fOutputList

        virtual void          AnaEvent() {}; //called for every event, both data and mc (to be implemented in derived class)
        virtual void          AnaMCEvent() {}; //called for every mc event (to be implemented in derived class)
        virtual void          AnaDATAEvent() {}; //called only every mc event (to be implemented in derived class)
        virtual void          AnaTrack() {}; //called for every track, both data and mc (to be implemented in derived class)
        virtual void          AnaMCTrack() {}; //called for every track with MC info (to be implemented in derived class)
        virtual void          AnaMCParticle() {}; //called for every MC Particle (to be implemented in derived class)
        
        //         
        virtual Bool_t          InitEvent();   // loads event-related properties, this need to be called in AnaEvent
        virtual Bool_t          InitMCEvent(); //load mc event-related properties, this needs to be called in anamcevent
        virtual Bool_t          InitTrack();   //initializes track related quantities, this needs ot be called in AnaTrack
        virtual Bool_t          InitTrackCuts(); //check all track cuts and set corresponding variables
        virtual Bool_t          InitMCTrack();
        virtual Bool_t          InitMCParticle();
        virtual Bool_t          InitVZERO();        
        
        virtual Bool_t          InitEventMult();   //initialize multiplicity specific variables, requires corresponding task
        virtual Bool_t          InitEventCent();   //initialize (old) centrality specific variables, requires corresponding task
        
        virtual Bool_t          InitTrackIP();  //initialize inner params 
        virtual Bool_t          InitTrackTPC();  //initialize inner params tpc
        
        virtual void            LoopOverAllTracks();    // loops over all tracks in the event, calls AnaTrack() and AnaMCTrack() for each track
        virtual void            LoopOverAllParticles(); // loops over all MC particles in the event, calls AnaMCParticle() for each particle
        
        virtual void          FillTriggerLog();   // fill the trigger histogram
        virtual void          CheckEvent();      // check if the event is ok
        
        // event related properties
        AliAnalysisManager*     fAnalysisManager;       //!<! analysis manager
        AliInputEventHandler*   fInputEventHandler;     //!<! input event handler
        UInt_t                  fEventSelected;         //!<! AliInputEventHandler::IsEventSelected()
        AliVEvent*              fEvent;                 //!<! input event       
        AliESDEvent*            fESD;                   //!<! input ESD event
        AliAODEvent*            fAOD;                   //!<! input AOD event
        AliMCEvent*             fMC;                    //!<! MC event
        AliESDVZERO*            fVZERO;                 //!<! VZERO information
        AliStack*               fMCStack;               //!<! MC stack
        AliHeader*              fMCHeader;              //!<! MC header
        AliGenEventHeader*      fMCGenHeader;           //!<! MC gen event header
        Int_t                   fNTracksESD;            //!<! number of esd tracks in the event
        Int_t                   fNTracksAcc;            //!<! number of accepted trackswith trackcuts
        Bool_t                  fIsMC;                  //!<! do we have an MC event?
        Int_t                   fRunNumber;             //!<! run n
        TString                 fRunNumberString;       //!<! run number as string
        TString                 fFiredTriggerClasses;   //!<! all trigger classes as string
        UInt_t                  fEventSpecie;           //!<! event specie
        Double_t                fOldCentPercentileV0M;  //!<! centrality percentile from old framework
        Double_t                fMultPercentileV0M;     //!<! centrality/multiplicity percentile from new framework
        Bool_t                  fEventCutsPassed;       //!<! accepted by AliEventCuts?
        Double_t                fZv;                    //!<! z vertex position
        Double_t                fMCzv;                  //!<! mc truth z vertex position                
        Int_t                   fMultMB;                //!<! MinBias Multiplicity (no of contributers to vertex)
        Double_t                fMultV0M;               //!<! v0a + v0c
        Double_t                fMCb;                   //!<! impact parameter in MC
        Int_t                   fMCnPrim;               //!<! number of primary particles according to mc
        Int_t                   fMCnPrimV0M;            //!<! number of primary particles in the v0 acceptance
        Int_t                   fMCnTracks;             //!<! number of "tracks" i.e. particles in MCevent
        Bool_t                  fIsTrigger;             //!<! is event triggered?
        Bool_t                  fHasVertex;             //!<! has the event a vertex?
        Bool_t                  fIsIncompleteDAQ;           //!<! incomplete daq event
        Bool_t                  fIsSPDClusterVsTrackletBG;  //!<! spd cluster vs tracklet bg
        Bool_t                  fIsFirstEventInChunk;       //!<! first event in chunk
        Bool_t                  fIsPileUpMV;                //!<! is pileup mv
        Bool_t                  fIsOutOfBunchPileUp;        //!<! is out of bunch pileup
        Bool_t                  fIsPileUpEvent;             //!<! is pileup event
        Bool_t                  fIsPileUpSPD;               //!<! is pileup from spd
        Bool_t                  fNOTIsVertexSelected2013pA; //!<! vertex 2013 pA rejected
        Bool_t                  fIsPileupFromSPD508;        //!<! is pileup from spd(5,0,8)        
        
        // track related properties
        AliESDtrack*            fESDTrack;              //!<! current esd track
        Double_t                fPt;                    //!<! track pT
        Double_t                fEta;                   //!<! track Eta
        Double_t                fPhi;                   //!<! track Phi
        Float_t                 fDCA[2];                //!<! impact parameter (DCA)
        Float_t                 fDCACov[3];             //!<! impat parameter (DCA) covariance
        Double_t                fDCAr;                  //!<! impact parameter (DCA) in xy-direction
        Double_t                fDCAz;                  //!<! impact parameter (DCA) in z-direction
        Double_t                fSigma1Pt2;             //!<! sigma(1/pT)**2
        Double_t                fSigma1Pt;              //!<! sigma(1/pT)
        Double_t                fSigned1Pt;             //!<! signed 1/pT
        
        AliMCParticle*          fMCParticle;            //!<! mc particle
        Int_t                   fMCLabel;               //!<! mc label
        Double_t                fMCPt;                  //!<! mc pt
        Double_t                fMCEta;                 //!<! mc eta
        Double_t                fMCPhi;                 //!<! mc phi
        Bool_t                  fMCisPrim;              //!<! is physical primary?
        Bool_t                  fMCisSec;               //!<! is secondary?
        Bool_t                  fMCisSecDecay;          //!<! is secondary from decay?
        Bool_t                  fMCisSecMat;            //!<! is secondary from material?
        Int_t                   fMCPrimSec;             //!<! status of mc track: 0=prim, 1=decay 2=material
        AlidNdPtTools::ParticleType   fMCParticleType;  //!<! which particle is it
        AlidNdPtTools::ProductionType fMCProdcutionType;//!<! production mechanism (prim,material,decay)
        Int_t                   fMCPDGCode;             //!<! PDG code
        Short_t                 fMCCharge;              //!<! charge in units of 1/3e
        Double_t                fMCQ;                   //!<! charge in units of e
        Bool_t                  fMCIsCharged;           //!<! charged particle
        Short_t                 fMCChargeSign;          //!<! Sign of the charge
        
        const AliExternalTrackParam*  fInnerP;          //!<! innerparams
        const AliExternalTrackParam*  fTPCinnerP;       //!<! tpc inner params
        Double_t                fPtInner;               //!<! inner param pt
        Double_t                fEtaInner;              //!<! inner param eta
        Double_t                fPhiInner;              //!<! inner param phi
        Double_t                fPtInnerTPC;            //!<! tpc inner pt
        Double_t                fEtaInnerTPC;           //!<! tpc inner eta
        Double_t                fPhiInnerTPC;           //!<! tpc inner phi
        Float_t                 fDCATPC[2];             //!<! TPC impact parameter (DCA)
        Float_t                 fDCACovTPC[3];          //!<! TPC impat parameter (DCA) covariance
        Double_t                fDCArTPC;               //!<! TPC impact parameter (DCA) in xy-direction 
        Double_t                fDCAzTPC;               //!<! TPC impact parameter (DCA) in z-direction        

        AliEventCuts            fEventCuts;             /// event cuts        
        Bool_t                  fUseEventCuts;          /// use event cuts?
        AliESDtrackCuts*        fESDtrackCutsM;         //-> trackcuts used for mult estimate
        Bool_t                  fAcceptTrackM;          //-> is track accepted by fESDtrackCutsM
        AliESDtrackCuts*        fESDtrackCuts[10];      //-> several track cuts that can be used in the analysis
        Bool_t                  fAcceptTrack[10];       //-> is track accepted by fESDtrackCuts[10]
        TString                 fMultEstimator;         // mult/cent estimator
        TString                 fCentEstimator;         // old cent estimator
        UInt_t                  fTriggerMaskRequired;   // only events with this trigger mask are accepted
        UInt_t                  fTriggerMaskRejected;   // reject events with this trigger mask

        TList*                  fOutputList;            //->  output list
        TH1D*                   fLogHist;               //->  generic log histogram
        TH1D*                   fLogErr;                //->  histogram with errors that should not occur
        TH1D*                   fLogEvent;              //->  event related logs
        TH1D*                   fRunHist;               //->  distribution of events in runs before event selection
        TH1D*                   fRunHistSelected;       //->  distribution of events in runs after event selection
        TH1D*                   fTrigInfo;              //->  all the trigger strings
        TH1D*                   fTrigHist;              //->  AliVEvent trigger classes
        
    private:
        AliAnalysisTaskMKBase(const AliAnalysisTaskMKBase&); // not implemented
        AliAnalysisTaskMKBase& operator=(const AliAnalysisTaskMKBase&); // not implemented
        
    /// \cond CLASSIMP      
    ClassDef(AliAnalysisTaskMKBase, 3);
    /// \endcond
    
};

#endif