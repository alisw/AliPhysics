#ifndef ALIMULTSELECTIONCALIBRATOR_H
#define ALIMULTSELECTIONCALIBRATOR_H

#include <iostream>
#include "TNamed.h"

using namespace std;
class AliESDEvent;
class AliMultSelectionCalibrator : public TNamed {
    
public:
    
    //Constructors/Destructor
    AliMultSelectionCalibrator();
    AliMultSelectionCalibrator(const char * name, const char * title = "Multiplicity Calibration Class");
    ~AliMultSelectionCalibrator();
    
    //void    Print(Option_t *option="") const;
    
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
    
    TString fInputFileName;  // Filename for TTree object for calibration purposes
    TString fBufferFileName; // Filename for TTree object (buffer file)
    TString fOutputFileName; // Filename for calibration OADB output
    
    // Object for storing event selection configuration
    AliMultSelectionCuts *fMultSelectionCuts;
    
    // TList object for storing histograms
    TList *fCalibHists; 

    ClassDef(AliMultSelectionCalibrator, 1);
};
#endif
