#ifndef AliAnalysisTaskPythiaNuclei_h
#define AliAnalysisTaskPythiaNuclei_h
 
#include "AliAnalysisTaskSE.h"

class AliInputEventHandler;
class AliMCEvent;
class TH1I;
class TH2F;
class TList;

class AliAnalysisTaskPythiaNuclei: public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskPythiaNuclei();
    AliAnalysisTaskPythiaNuclei(const Char_t* name);
    virtual ~AliAnalysisTaskPythiaNuclei();

    virtual void UserCreateOutputObjects();
    virtual void Init(); 
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

  protected:
	
    AliMCEvent*           fMcEvent;    //!<! MC event                    
    AliInputEventHandler* fMcHandler;  //!<! MCEventHandler             
    TList*                fOutputList; //!<! output list for histograms

    TH1I*                 fCounter;    //!<!
    TH2F*                 fYPtA;       //!<!
    TH2F*                 fYPtM;       //!<!
    TH2F*                 fYPA;        //!<!
    TH2F*                 fYPM;        //!<!

    ClassDef(AliAnalysisTaskPythiaNuclei, 1) 
  };
 
#endif
