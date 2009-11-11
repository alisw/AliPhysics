#ifndef ALIANALYSISTASKSE_H
#define ALIANALYSISTASKSE_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliAnalysisTask.h"
class AliVEvent;
class AliAODEvent;
class AliAODHeader;
class AliAODTracklets;
class AliAODCaloCells;
class AliMCEvent;
class AliInputEventHandler;

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
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Exec(Option_t* option);
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    virtual void Init() {;}
    virtual Bool_t Notify();
    // To be implemented by user
    virtual void UserCreateOutputObjects()  {;}
    virtual void UserExec(Option_t* /*option*/) {;}
    virtual Bool_t UserNotify() {return kTRUE;}
    virtual void   NotifyRun()  {;}

    // Helpers for adding branches to the AOD
   virtual void AddAODBranch(const char* cname, void* addobj, const char *fname="");
// Getters
    virtual Int_t        DebugLevel()  {return fDebug;     }
    virtual AliVEvent*   InputEvent()  {return fInputEvent;}
    virtual AliAODEvent* AODEvent()    {return fOutputAOD; }
    virtual TTree*       OutputTree()  {return fTreeA;     }
    virtual AliMCEvent*  MCEvent()     {return fMCEvent;   }
    virtual Long64_t     Entry()       {return fEntry;     }
    virtual const char*  CurrentFileName();
    virtual Bool_t       IsStandardAOD() const;
    virtual TList*       GetQAHistos()   const {return fHistosQA;}
  protected:
    Int_t                 fDebug;           //  Debug flag
    // IO
    Int_t                 fEntry;           //  Current entry in the chain
    AliVEvent*            fInputEvent;      //! VEvent Input
    AliInputEventHandler* fInputHandler;    //! Input Handler
    AliAODEvent*          fOutputAOD;       //! AOD out 
    AliMCEvent*           fMCEvent;         //! MC
    TTree*                fTreeA;           //  AOD output Tree
    Int_t                 fCurrentRunNumber;//! Current run number
    // Output histos for QA
    TList*                fHistosQA;        //! Output histos for QA
    // Provisions for replication
    static AliAODHeader*    fgAODHeader;        //! Header for replication
    static TClonesArray*    fgAODTracks;        //! Tracks for replication
    static TClonesArray*    fgAODVertices;      //! Vertices for replication
    static TClonesArray*    fgAODV0s;           //! V0s for replication
    static TClonesArray*    fgAODPMDClusters;   //! PMDClusters for replication
    static TClonesArray*    fgAODJets;          //! Jets for replication
    static TClonesArray*    fgAODFMDClusters;   //! FMDClusters for replication
    static TClonesArray*    fgAODCaloClusters;  //! CaloClusters for replication
    static TClonesArray*    fgAODMCParticles;   //! MC Particles for replicatio
    static AliAODTracklets* fgAODTracklets;     //! Tracklets for replication
    static AliAODCaloCells* fgAODEmcalCells;    //! Emcal Cell replication
    static AliAODCaloCells* fgAODPhosCells;     //! Phos  Cell replication
	
    ClassDef(AliAnalysisTaskSE, 1); // Analysis task for standard jet analysis
};
 
#endif
