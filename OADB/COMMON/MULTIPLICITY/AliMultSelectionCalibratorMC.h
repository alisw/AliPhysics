#ifndef ALIMULTSELECTIONCALIBRATORMC_H
#define ALIMULTSELECTIONCALIBRATORMC_H

#include <iostream>
#include "TNamed.h"

using namespace std;
class AliESDEvent;
class AliMultSelectionCalibratorMC : public TNamed {
    
public:
    
    //Constructors/Destructor
    AliMultSelectionCalibratorMC();
    AliMultSelectionCalibratorMC(const char * name, const char * title = "Multiplicity Calibration Class");
    ~AliMultSelectionCalibratorMC();
    
    //void    Print(Option_t *option="") const;
    
    //_________________________________________________________________________
    //Interface: steering functions to be used in calibration macro
  
    //Set Filenames
    void SetInputFileData ( TString lFile ) { fInputFileNameData = lFile.Data(); }
    void SetInputFileOADB ( TString lFile ) { fInputFileNameOADB = lFile.Data(); }
    void SetInputFileMC   ( TString lFile ) { fInputFileNameMC = lFile.Data(); } 
    
    void SetBufferFileData    ( TString lFile ) { fBufferFileNameData = lFile.Data(); }
    void SetBufferFileMC      ( TString lFile ) { fBufferFileNameMC   = lFile.Data(); } 
    void SetOutputFile    ( TString lFile ) { fOutputFileName = lFile.Data(); }
    //Set Boundaries to find
    void SetBoundaries ( Long_t lNB, Double_t *lB ){
        if ( lNB<2 || lNB > 1e+6 ){
            cout<<"Please make sure you are using a reasonable number of boundaries!"<<endl;
            lNB = -1;
        }
        lDesiredBoundaries = lB;
        lNDesiredBoundaries = lNB;
    }
    
    //Getter for event selection criteria
    AliMultSelectionCuts * GetEventCuts() const { return fMultSelectionCuts;     }
    
    //Getter for MultSelection object
    //WARNING - should not be used in this case! 
    AliMultSelection * GetMultSelection() const { return fSelection;     }

    //Getter for MultInput object
    AliMultInput * GetMultInput() const { return fInput;     }
    
    //Setter for golden run
    void SetRunToUseAsDefault ( Int_t lRunNumber ) { fRunToUseAsDefault = lRunNumber; }
    
    //Getter for golden run
    Int_t GetRunToUseAsDefault() const { return fRunToUseAsDefault; }
    
    //Configure standard input
    void SetupStandardInput();
    
    //Switch to configure <Ntracklet> fit type
    void SetUseQuadraticMapping(Bool_t lOpt){ fkUseQuadraticMapping=lOpt; } 
    
    //Master Function in this Class: To be called once filenames are set
    Bool_t Calibrate();
    
    //Helper
    Float_t MinVal( Float_t A, Float_t B );
    
private:
    AliMultInput     *fInput;     //Object for all input
    AliMultSelection *fSelection; //Object for all estimators

    //Calibration Boundaries to locate
    Double_t *lDesiredBoundaries;
    Long_t   lNDesiredBoundaries;
    
    //Run to use as default for scaling factor in this period
    Int_t fRunToUseAsDefault;
    
    TString fInputFileNameData;  // Filename for TTree object for calibration purposes
    TString fInputFileNameOADB;  // Filename for TTree object for calibration purposes
    TString fInputFileNameMC;    // Filename for TTree object for calibration purposes
    TString fBufferFileNameData; // Filename for TTree object (buffer file)
    TString fBufferFileNameMC;   // Filename for TTree object (buffer file)
    TString fOutputFileName;     // Filename for calibration OADB output (MC)
    
    //Configuration
    Bool_t fkUseQuadraticMapping; //switch to toggle quadratic <Ntracklets> fits
    
    // Object for storing event selection configuration
    AliMultSelectionCuts *fMultSelectionCuts;
    
    // TList object for storing histograms
    TList *fCalibHists; 

    ClassDef(AliMultSelectionCalibratorMC, 1);
};
#endif
