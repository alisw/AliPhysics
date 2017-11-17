#ifndef ALIMULTSELECTIONCALIBRATOR_H
#define ALIMULTSELECTIONCALIBRATOR_H

#include <iostream>
#include "TNamed.h"
#include "TArrayD.h"
#include "TArrayI.h"
//For Run Ranges functionality
#include <map>

using namespace std;
class AliESDEvent;
class AliMultSelectionCalibrator : public TNamed {
public:

    AliMultSelectionCalibrator(const char * name="AliMultSelectionCalibrator", const char * title = "Multiplicity Calibration Class");
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
        lDesiredBoundaries.Set(lNB, lB);
    }

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

    //Configure standard input
    void SetupStandardInput();

    //Master Function in this Class: To be called once filenames are set
    Bool_t Calibrate();

    //Helper
    Float_t MinVal( Float_t A, Float_t B );

private:
    AliMultInput     *fInput;     //Object for all input
    AliMultSelection *fSelection; //(current) transient pointer object

    //Calibration Boundaries to locate
    TArrayD lDesiredBoundaries;

    Int_t fRunToUseAsDefault; //Give preference for this run to be the default

    //Run Ranges map - master storage
    Long_t fNRunRanges;
    std::map<int, int> fRunRangesMap;
    TArrayI fFirstRun;
    TArrayI fLastRun;

    TList *fMultSelectionList; // List of AliMultSelection objects to be used per run period

    TString fInputFileName;  // Filename for TTree object for calibration purposes
    TString fBufferFileName; // Filename for TTree object (buffer file)
    TString fOutputFileName; // Filename for calibration OADB output

    // Object for storing event selection configuration
    AliMultSelectionCuts *fMultSelectionCuts;

    // TList object for storing histograms
    TList *fCalibHists;

    ClassDef(AliMultSelectionCalibrator, 3);
    //(this classdef is only for bookkeeping, class will not usually
    // be streamed according to current workflow except in very specific
    // tests!)
    //2 - Adjustments of extra event selections
};
#endif
