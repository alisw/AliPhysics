#ifndef ALIANALYSISTASKSE_H
#define ALIANALYSISTASKSE_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliAnalysisTask.h"
#include "AliVEvent.h"

class AliVfriendEvent;
class AliVEventHandler;
class AliAODEvent;
class AliAODHeader;
class AliTOFHeader;
class AliAODVZERO;
class AliAODTracklets;
class AliAODCaloCells;
class AliAODCaloTrigger;
class AliMCEvent;
class AliMCEventHandler;
class AliInputEventHandler;
class AliMultiInputEventHandler;
class AliAnalysisCuts;
class AliESDfriend;
class AliEventTag;
class AliTrackSelectionFactory;
class AliVTrackSelection;

class TTree;
class TList;


class AliAnalysisTaskSE : public AliAnalysisTask
{
 public:
    AliAnalysisTaskSE();
    AliAnalysisTaskSE(const char* name);
    AliAnalysisTaskSE(const AliAnalysisTaskSE& obj);
    AliAnalysisTaskSE& operator=(const AliAnalysisTaskSE& other);
    virtual ~AliAnalysisTaskSE() {;}
    // Implementation of interface methods
    virtual void   ConnectInputData(Option_t *option = "");
    virtual void   CreateOutputObjects();
    virtual void   Exec(Option_t* option);
    virtual void   SetDebugLevel(Int_t level) {fDebug = level;}
    virtual void   Init() {;}
    virtual Bool_t Notify();
    // To be implemented by user
    virtual void   UserCreateOutputObjects()  {;}
    virtual void   UserExec(Option_t* /*option*/) {;}
    virtual void   UserExecMix(Option_t */*option*/) {;}
    virtual Bool_t UserNotify() {return kTRUE;}
    virtual void   NotifyRun()  {;}
    
    // Helpers for adding branches to the AOD
    virtual void   AddAODBranch(const char* cname, void* addobj, const char *fname="");
    // Event Selection
    virtual void   SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kMB) {fOfflineTriggerMask = offlineTriggerMask;}
    // Loading the declared input branches
    void           LoadBranches() const;
 // Getters
    virtual Int_t         DebugLevel() const  {return fDebug;     }
    virtual AliVEvent*    InputEvent() const  {return fInputEvent;}
    virtual AliVfriendEvent* ESDfriend()  const  {return fESDfriend; }
    virtual AliAODEvent*  AODEvent()   const  {return fOutputAOD; }
    virtual TTree*        OutputTree() const  {return fTreeA;     }
    virtual AliMCEvent*   MCEvent()    const  {return fMCEvent;   }
    virtual Long64_t      Entry()      const  {return fEntry;     }
    virtual const AliEventTag *EventTag() const;
    virtual const char*   CurrentFileName();
    virtual Bool_t        IsStandardAOD() const;
    virtual TList*        GetQAHistos()   const {return fHistosQA;}
    virtual Bool_t        IsEventInBinZero() { return kFALSE;}
    virtual UInt_t        GetCollisionCandidates() const { return fOfflineTriggerMask;}

    void SetTrackSelectionFactory(AliTrackSelectionFactory *factory)  { fTrackSelectionFactory = factory; }
    void SetTrackSelection(AliVTrackSelection *sel)                   { fTrackSelection = sel; }

 protected:
    void ConnectMultiHandler();
    void DisconnectMultiHandler();

    TObjArray *GetAcceptedTracks();

  protected:
    Int_t                 fDebug;           //  Debug flag
    // IO
    Int_t                 fEntry;           //  Current entry in the chain
    AliVEvent*            fInputEvent;      //! VEvent Input
    AliVfriendEvent*      fESDfriend;       //! ESD friend
    AliVEventHandler* fInputHandler;    //! Input Handler
    AliAODEvent*          fOutputAOD;       //! AOD out 
    AliMCEvent*           fMCEvent;         //! MC
    TTree*                fTreeA;           //  AOD output Tree
    Int_t                 fCurrentRunNumber;//! Current run number
    // Output histos for QA
    TList*                fHistosQA;        //! Output histos for QA
    // Provisions for replication
    static AliVHeader*      fgAODHeader;        //! Header for replication
    static AliTOFHeader*    fgTOFHeader;        //! TOFHeader for replication
    static AliAODVZERO*     fgAODVZERO;         //! VZERO for replication
    static TClonesArray*    fgAODTracks;        //! Tracks for replication
    static TClonesArray*    fgAODVertices;      //! Vertices for replication
    static TClonesArray*    fgAODV0s;           //! V0s for replication
    static TClonesArray*    fgAODPMDClusters;   //! PMDClusters for replication
    static TClonesArray*    fgAODJets;          //! Jets for replication
    static TClonesArray*    fgAODFMDClusters;   //! FMDClusters for replication
    static TClonesArray*    fgAODCaloClusters;  //! CaloClusters for replication
    static AliAODCaloTrigger* fgAODEMCALTrigger; //! Emcal Trigger for replication
    static AliAODCaloTrigger* fgAODPHOSTrigger;  //! Phos Trigger for replication
    static TClonesArray*    fgAODMCParticles;   //! MC Particles for replicatio
    static AliAODTracklets* fgAODTracklets;     //! Tracklets for replication
    static AliAODCaloCells* fgAODEmcalCells;    //! Emcal Cell replication
    static AliAODCaloCells* fgAODPhosCells;     //! Phos  Cell replication
    static TClonesArray*    fgAODDimuons;       //! Dimuons replication
    static TClonesArray*    fgAODHmpidRings;    //! HMPID replication
    // Event Selection
    UInt_t fOfflineTriggerMask;   //  Task processes collision candidates only
    // Event Mixing
    AliMultiInputEventHandler *fMultiInputHandler;  //! pointer to multihandler
    AliInputEventHandler      *fMCEventHandler;     //! pointer to MCEventHandler
    AliTrackSelectionFactory  *fTrackSelectionFactory; /// track selection factory
    AliVTrackSelection        *fTrackSelection;        /// track selection
    ClassDef(AliAnalysisTaskSE, 5); // Analysis task for standard jet analysis
};
 
#endif
