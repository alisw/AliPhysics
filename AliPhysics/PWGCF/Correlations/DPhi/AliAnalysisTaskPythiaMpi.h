#ifndef ALIANALYSISTASKPYTHIAMPI_H
#define ALIANALYSISTASKPYTHIAMPI_H
 

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"


class TH1I;
class TParticle;
class AliStack;
class AliVVertex;
class AliVParticle;

class AliAnalysisTaskPythiaMpi: public AliAnalysisTaskSE
{
 public:

  
  AliAnalysisTaskPythiaMpi();
  AliAnalysisTaskPythiaMpi(const Char_t* name);
  AliAnalysisTaskPythiaMpi(const AliAnalysisTaskPythiaMpi& c);
  AliAnalysisTaskPythiaMpi& operator= (const AliAnalysisTaskPythiaMpi& c);
  virtual ~AliAnalysisTaskPythiaMpi();
  

    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init(); 
    virtual void LocalInit();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);


 protected:
	
    AliMCEvent*              fMcEvent;    //! MC event                    
    AliInputEventHandler*    fMcHandler;  //! MCEventHandler             
 
    TList*        fOutputList; //! output list for histograms

    ClassDef(AliAnalysisTaskPythiaMpi, 2) 
};
 
#endif
