#ifndef AliPmdCustomCalibTask_h
#define AliPmdCustomCalibTask_h


/*------------------------------------------------------------------------
  .  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
  .		     Satyajit Jena, IIT Bombay
  .		     sjena@cern.ch
  .		     3/8/2011
  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
  ------------------------------------------------------------------------*/

class TH1F;
class TH2F;
class TList;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliPmdCustomCalibTask : public AliAnalysisTaskSE {
 public:
    AliPmdCustomCalibTask();
    AliPmdCustomCalibTask(const char *name);
    virtual ~AliPmdCustomCalibTask();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    
    void SetSmn(Int_t i ) {fSmn = i; }


 private:
    
    TList *fOutput;              //! Output list
    Int_t fSmn;                  //! Module under Consideration 
    TH1F *fHistAdcPre;           //! ADC spectrum for PRE for entire module 
    TH1F *fHistAdcCpv;           //! ADC spectrum for CPV for entire module
    TH2F *fHistClusterXYPre;     //! cluster XY for PRE for entire module
    TH2F *fHistClusterXYCpv;     //! cluster XY for CPV for entir module
    
    TH1F *fHistAdcPreRC[48][96]; //! cell-wise ADC PRE
    TH1F *fHistAdcCpvRC[48][96]; //! cell-wise ADC CPV
       
    AliPmdCustomCalibTask(const AliPmdCustomCalibTask&); // not implemented
    AliPmdCustomCalibTask& operator=(const AliPmdCustomCalibTask&); // not implemented
    
    ClassDef(AliPmdCustomCalibTask, 1); // example of analysis
};

#endif

