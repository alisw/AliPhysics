#ifndef AliAnalysisTaskEmeanCalib_H
#define AliAnalysisTaskEmeanCalib_H
#include "AliAnalysisTaskSE.h"

class AliESDEvent;
class AliPHOSCalibData;
class AliPHOSGeometry;

class AliAnalysisTaskEmeanCalib : public AliAnalysisTaskSE {
    
public:
    
    AliAnalysisTaskEmeanCalib();
    AliAnalysisTaskEmeanCalib(const char *name);
    virtual ~AliAnalysisTaskEmeanCalib() {}
    
    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *);
    
    void SetCalibrations(AliPHOSCalibData* c);
    
protected:
    
    void SetGeometry();
    void SetMisalignment();
    
private:
    
    AliAnalysisTaskEmeanCalib(const AliAnalysisTaskEmeanCalib&);
    AliAnalysisTaskEmeanCalib& operator=(const AliAnalysisTaskEmeanCalib&);
    
    AliESDEvent* fESDEvent;             // ESD event
    TList *fOutput;                     //< List of histograms for data
    
    AliPHOSGeometry*  fPHOSGeo;         //!
    AliPHOSCalibData* fCalibData;       //!
    Int_t fRunNumber;                   //!
    
    TH2F* hmpt_diff;                    //! two clusters in different modules
    TH2F* hmpt[5];                      //! [0]-[3] for modules 1-4,
    TH2F* hmpt_orig[5];                 //! [4] - all modules.
    
    TH1F* hCellAmplitudes[5][64][56];   //!
    TH2F* hCellMeanAmps[5];             //!
    
    Float_t fCC[5][64][56];             //! (Re)calibrations ~1
    
    ClassDef(AliAnalysisTaskEmeanCalib, 1);
    
};

#endif
