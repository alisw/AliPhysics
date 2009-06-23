#ifndef AliT0MultiplicityTask_cxx
#define AliT0MultiplicityTask_cxx

class TList;
class TH1F;
class TH2F;
class TClonesArray;


#include "AliAnalysisTaskSE.h"

class AliT0MultiplicityTask: public AliAnalysisTaskSE
{
  public:
   AliT0MultiplicityTask ();
   AliT0MultiplicityTask (const char* name);
    virtual ~AliT0MultiplicityTask() {}
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
  
  private:
    TList* fListOfHistos;
    TH1F* fOrA; // Or A
    TH1F* fOrC; // Or C
    TH1F* fMean; // multiplicity histogram
    TH1F* fVertex; // multiplicity histogram
    TH1F* fTime;
    TH1F* fAmp;
    TH1F* fTotalMult;
    TH2F* fMultRecTotal;
    TH2F* fMultRecRealA;
    TH2F* fMultRecRealC;
    TH1F*  fPrim;
    TH1F* fRef;
    TH1F* fRec;
    
  AliT0MultiplicityTask (const AliT0MultiplicityTask&); // not implemented
  AliT0MultiplicityTask& operator=(const AliT0MultiplicityTask&); // not implemented

    ClassDef(AliT0MultiplicityTask, 1); // example of analysis 
};

#endif
