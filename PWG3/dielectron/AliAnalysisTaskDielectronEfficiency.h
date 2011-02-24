#ifndef ALIANALYSISTASKDIELECTRONEFFICIENCY_H
#define ALIANALYSISTASKDIELECTRONEFFICIENCY_H

/* $Id$ */ 

//#####################################################
//#                                                   # 
//#  Analysis Task for Event Mixing for dielectron    #
//#                                                   #
//#     J.Wiechula (Jens.Wiechula@cern.ch)            #
//#                                                   #
//#####################################################

#include <AliAnalysisTask.h>

class TDatabasePDG;

class AliESDtrackCuts;
class AliKineTrackCuts;
class AliDielectronHistos;
class AliVEvent;
class AliStack;

class AliAnalysisTaskDielectronEfficiency : public AliAnalysisTask {

public:
  AliAnalysisTaskDielectronEfficiency();
  AliAnalysisTaskDielectronEfficiency(const char *name);
  virtual ~AliAnalysisTaskDielectronEfficiency();
  
  void SetupDefaultCuts(Int_t type=0);
  
	virtual void ConnectInputData(Option_t *);
	virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *);

//   virtual Int_t Merge(TCollection *list);
  //getters
  AliESDtrackCuts  *GetESDTrackCuts()   const {return fESDtrackCuts;}
  AliKineTrackCuts *GetKineCutsLeg()    const {return fKineCutsLegs;}
  AliKineTrackCuts *GetKineCutsMother() const {return fKineCutsMother;}
  //
  Int_t GetIdMother() const {return fIdMCMother;}

  //setters
  void SetIdMother(Int_t id) {fIdMCMother=id;}
  void SetIdDaughters(Int_t idPositive, Int_t idNegative) {fIdMCDaughterP=idPositive; fIdMCDaughterN=idNegative;}
  
private:

  AliVEvent           *fInputEvent;     //! Input event
  AliDielectronHistos *fHist;           //! Histogram container
  //cut objects
  AliESDtrackCuts     *fESDtrackCuts;   //  ESD track cuts
  AliKineTrackCuts    *fKineCutsLegs;   //  MC cuts on Legs
  AliKineTrackCuts    *fKineCutsMother; //  MC cuts on Mother (Id see below)
  
  Int_t                fIdMCMother;   //  MC Id of mother particle of interest (eg. Jpsi=443)
  Int_t                fIdMCDaughterP; //  MC Id of legs without sign
  Int_t                fIdMCDaughterN; //  MC Id of legs without sign
  //
  TDatabasePDG        *fPDG;           //! PDG database

  void FillPlots(AliVEvent *event);
  void FillMCInfo(AliStack * const pStack);

  AliAnalysisTaskDielectronEfficiency(const AliAnalysisTaskDielectronEfficiency &c);
  AliAnalysisTaskDielectronEfficiency& operator= (const AliAnalysisTaskDielectronEfficiency &c);
  
	ClassDef(AliAnalysisTaskDielectronEfficiency, 1);
};
#endif
