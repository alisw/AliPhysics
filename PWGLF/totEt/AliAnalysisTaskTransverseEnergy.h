//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for Et tasks
//  - reconstruction and MonteCarlo output
//
//*-- Authors: Oystein Djuvsland (Bergen)
//_________________________________________________________________________//
#ifndef ALIANALYSISTASKTRANSVERSEENERGY_H
#define ALIANALYSISTASKTRANSVERSEENERGY_H

#include "AliAnalysisTaskSE.h"

class AliCentrality;
class AliAnalysisEtSelectionHandler;
class AliESDtrackCuts;
class AliESDEvent;
class TH2F;

class AliAnalysisTaskTransverseEnergy : public AliAnalysisTaskSE
{

public:

    /** Constructor */
    AliAnalysisTaskTransverseEnergy(const char* name, Bool_t isMc);

    /** Destructro */
    virtual ~AliAnalysisTaskTransverseEnergy();

    AliESDtrackCuts* GetTPCITSTrackCuts(){return (AliESDtrackCuts*) fOutputList->FindObject("fEsdTrackCuts");}
    AliESDtrackCuts* GetTPCOnlyTrackCuts(){return (AliESDtrackCuts*) fOutputList->FindObject("fEsdTrackCutsTPCOnly");}
    AliESDtrackCuts* GetITSTrackCuts(){return (AliESDtrackCuts*) fOutputList->FindObject("fEsdTrackCutsITS");}

    Int_t LoadPhysicsSelection(TString name);
    
    void SetMcData(Bool_t isMc = true) { fIsMc = isMc; }

protected:

    /** Check if the physics selection is still valid, if not load new */
    Int_t CheckPhysicsSelection(Int_t runNumber); // check if the current physics selection is valid, if not load new
    
    /** Check if the event is a physics event */
    Bool_t IsPhysicsSelected() const;
    
    /** Get the centrality object */
    AliCentrality* GetCentralityObject();
    
    /** The ESD event */
    AliESDEvent *fESDEvent; //The ESD event
  
    TString       fMCConfigFile;        // the name of the ConfigFile
    TString       fRecoConfigFile;        // the name of the ConfigFile

    TH2F *fHistEtRecvsEtMC; // Rec vs MC histo 
    TH2F *fHistEtRecOverEtMC; // Rec over MC histo 
    TH2F *fHistDiffEtRecEtMCOverEtMC; // Rec - MC over MC histo 

    AliESDtrackCuts* fEsdtrackCutsITSTPC; // track cuts ITS&TPC
    AliESDtrackCuts* fEsdtrackCutsTPC; // track cuts TPC
    AliESDtrackCuts* fEsdtrackCutsITS; // track cuts ITS

    TList *fOutputList; //output list
    
    TString fPhysSelTaskName; // If we need to access the physics selection task
    TString fCentSelTaskName; // If we need to access the centrality selection task
    
    Bool_t fIsMc; // Are we analysing MC data

    Bool_t fUsingDefaultSelection; // Are we using the default physics selection

private:

    Int_t fCurrentRunNum; // The current run number
    
    AliAnalysisEtSelectionHandler* fSelectionHandler; //! Taking care of loading the correct selections
    AliAnalysisTaskTransverseEnergy();
  //Declare it private to avoid compilation warning
    AliAnalysisTaskTransverseEnergy & operator = (const AliAnalysisTaskTransverseEnergy &);//assignment
    AliAnalysisTaskTransverseEnergy(const AliAnalysisTaskTransverseEnergy &) ; //copy constructor

    ClassDef(AliAnalysisTaskTransverseEnergy, 1)
};

#endif // ALIANALYSISTASKTRANSVERSEENERGY_H
