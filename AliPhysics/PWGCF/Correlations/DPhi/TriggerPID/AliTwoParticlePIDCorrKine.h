#ifndef AliTwoParticlePIDCorrKine_H
#define AliTwoParticlePIDCorrKine_H

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"




class TH1F;
class TList;
class TParticle;
class AliStack;
class AliVEvent;
class AliVVertex;
class AliVParticle;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif


class AliTwoParticlePIDCorrKine : public AliAnalysisTaskSE {

 public:


  AliTwoParticlePIDCorrKine();
  AliTwoParticlePIDCorrKine(const Char_t* name);
  
  virtual ~AliTwoParticlePIDCorrKine();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  virtual void UserCreateOutputObjects();
 
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  protected:
    AliVEvent*              fEvent;    //! MC event                       
    AliInputEventHandler*    fMcHandler;  //! MCEventHandler                 
    

  TH1F  *fHistEventsProcessed;   //! histo for monitoring the number of events processed slot 1
  TList       *fOutputList; //! Output list
  TH1F       *fHistZvtx;
  TH1F        *fHistPt; //!Pt spectrum
  TH1F        *fHistImpact;
Double_t fZvtxLim;
TString fCentralityFrom;
TString fCentralityEstimator;
	
    AliTwoParticlePIDCorrKine(const AliTwoParticlePIDCorrKine&); // not implemented
    AliTwoParticlePIDCorrKine& operator=(const AliTwoParticlePIDCorrKine&); // not implemented

  ClassDef(AliTwoParticlePIDCorrKine,1)
    
};

#endif
