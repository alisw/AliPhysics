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
 
    TList*        fOutputList;
    TH1I*         fHistEvents;
    TH1F*         fHistPt; //pT distribution
    TH1F*         fHistEta; //eta distribution
    TH1F*         fHistMpi; //MPIs distribution
    TH2F*         fHistMultMpi; //Multiplicity distribution vs MPIs
    TH2F*         fHistdNdetaMpi; //dNdEta vs MPIs

    ClassDef(AliAnalysisTaskPythiaMpi, 1) 
};
 
#endif
