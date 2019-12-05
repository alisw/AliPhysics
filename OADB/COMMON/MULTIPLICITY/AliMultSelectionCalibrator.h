#ifndef ALIMULTSELECTIONCALIBRATOR_H
#define ALIMULTSELECTIONCALIBRATOR_H

#include <iostream>
#include "TNamed.h"
#include "AliVEvent.h"
//For Run Ranges functionality
#include <map>

using namespace std;
class AliESDEvent;
class AliMultSelectionCalibrator : public TNamed {
    
public:
    
    //Constructors/Destructor
    AliMultSelectionCalibrator();
    AliMultSelectionCalibrator(const char * name, const char * title = "Multiplicity Calibration Class");
    ~AliMultSelectionCalibrator();
    
    //void    Print(Option_t *option="") const;
    
    //_________________________________________________________________________
    //Interface: steering functions to be used in calibration macro
  
    //Set Filenames
    void SetInputFile ( TString lFile ) { fInputFileName = lFile.Data(); } 
    void SetBufferFile ( TString lFile ) { fBufferFileName = lFile.Data(); } 
    void SetOutputFile ( TString lFile ) { fOutputFileName = lFile.Data(); }
    //Set Boundaries to find
    void SetBoundaries ( Long_t lNB, Double_t *lB ){
        if ( lNB<2 || lNB > 1e+6 ){
            cout<<"Please make sure you are using a reasonable number of boundaries!"<<endl;
            lNB = -1;
        }
        lDesiredBoundaries = lB;
        lNDesiredBoundaries = lNB;
    }
    
    //Task Configuration: trigger selection
    //This is in addition to the "IsTriggered" functionality. 
    void SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes trigType) { fTrigType = trigType; fCheckTriggerType=kTRUE; }
    
    void SetFiredTriggerString(TString lData) { fFiredTrigString = lData.Data(); } 
    
    //Run Ranges Interface
    Long_t GetNRunRanges() const {return fNRunRanges; }
    void AddRunRange ( Int_t lFirst, Int_t lLast, AliMultSelection *lMultSelProvided );
    
    //Getter for event selection criteria
    AliMultSelectionCuts * GetEventCuts() const { return fMultSelectionCuts;     }
    
    //Getter for MultSelection object
    AliMultSelection * GetMultSelection() const { return fSelection;     }
    
    //Setter (warning: not a copy...) 
    void SetMultSelection(AliMultSelection *lMultSelProvided ){ fSelection = lMultSelProvided; }

    //Getter for MultInput object
    AliMultInput * GetMultInput() const { return fInput;     }
    
    //Setter for golden run
    void SetRunToUseAsDefault ( Int_t lRunNumber ) { fRunToUseAsDefault = lRunNumber; }
    
    //Getter for golden run
    Int_t GetRunToUseAsDefault() const { return fRunToUseAsDefault; } 
    
    //Getter for golden run
    void SetMaxEventsPerRun(Long_t lVal) { fMaxEventsPerRun = lVal; }
    
    //Configure standard input
    void SetupStandardInput();
    
    //Filter only flag
    void SetFilterOnly(Bool_t lOpt = kTRUE){ fPrefilterOnly = lOpt; }
    
    //Master Function in this Class: To be called once filenames are set
    Bool_t Calibrate();
    
    //Helper
    Float_t MinVal( Float_t A, Float_t B );
    
private:
    AliMultInput     *fInput;     //Object for all input
    AliMultSelection *fSelection; //(current) transient pointer object

    //Calibration Boundaries to locate
    Double_t *lDesiredBoundaries;
    Long_t   lNDesiredBoundaries;
    
    Int_t fRunToUseAsDefault; //Give preference for this run to be the default
    
    Long_t fMaxEventsPerRun; //Implemented to get a grip on huge runs
    
    Bool_t fCheckTriggerType; 
    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type to calibrate
    Bool_t fPrefilterOnly; //stop before calibrating stuff
    TString fFiredTrigString; //select on fired trigger string if desired
    
    //Run Ranges map - master storage
    Long_t fNRunRanges;
    std::map<int, int> fRunRangesMap;
    Int_t fFirstRun[1000];
    Int_t fLastRun[1000]; 
    
    TList *fMultSelectionList; // List of AliMultSelection objects to be used per run period 
    
    TString fInputFileName;  // Filename for TTree object for calibration purposes
    TString fBufferFileName; // Filename for TTree object (buffer file)
    TString fOutputFileName; // Filename for calibration OADB output
    
    // Object for storing event selection configuration
    AliMultSelectionCuts *fMultSelectionCuts;
    
    // TList object for storing histograms
    TList *fCalibHists; 

    ClassDef(AliMultSelectionCalibrator, 2);
    //(this classdef is only for bookkeeping, class will not usually
    // be streamed according to current workflow except in very specific
    // tests!) 
    //2 - Adjustments of extra event selections
};
#endif
