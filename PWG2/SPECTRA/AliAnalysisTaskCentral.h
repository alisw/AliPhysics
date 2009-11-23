/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * $Id$
 */

//  *************************************
//  * analysis task for azimuthal isotropic   *
//  * expansion in highly central collisions     *
//  * author: Cristian Andrei                          *
//  *         acristian@niham.nipne.ro              *
//  * ************************************


#ifndef ALI_ANALYSIS_TASK_CENTRAL_H
#define ALI_ANALYSIS_TASK_CENTRAL_H

#include "AliAnalysisTask.h"

class TH1F;
class TObject;
class TObjArray;
class TList;

class AliESDEvent;
class AliMCEvent;
class AliCFContainer;


class AliAnalysisTaskCentral : public AliAnalysisTask {
 public:
  AliAnalysisTaskCentral(const char *name="AliAnalysisTaskCentral");
  virtual ~AliAnalysisTaskCentral();
  
  void SetCuts(Int_t const  no, TObjArray* const array) {fCutsList[no] = array;} //used to set the cuts to the Task

  void SendEvent(TObject *obj) const;  //used to send the MCEvent to the cuts that need it (i.e MC IsPrimary)

  Bool_t CheckCuts(Int_t no, TObject *obj) const; //used to check if a track/particle is selected

  void SetSimulation(Bool_t type) {fSim = type;} // set to kTRUE if running on simulated data

  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliAnalysisTaskCentral(const AliAnalysisTaskCentral& ref);
  AliAnalysisTaskCentral& operator=(const AliAnalysisTaskCentral& ref);

  void InitCuts();  //initialize cuts

  AliESDEvent *fESD;   //ESD object
  AliMCEvent  *fMC;    //MC Object

  TH1F        *fNoEvt;  //Number of events processed

  AliCFContainer *fCFContainerPi; // CF Container used to calc/apply eff - Pions
  AliCFContainer *fCFContainerK; // CF Container used to calc/apply eff - Kaons
  AliCFContainer *fCFContainerP; // CF Container used to calc/apply eff - Protons

  Bool_t fSim; // kTRUE = running on simulated data (look at MC Truth too)

  TObjArray *fCutsList[10];  //list containing the cuts

  TList *fOutList;  //list containing the output objects


  ClassDef(AliAnalysisTaskCentral, 1);
};

#endif
