#ifndef ALIANALYSISTASK_H
#define ALIANALYSISTASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 31/05/2006

//==============================================================================
//   AliAnalysysTask - Class representing a basic analysis task. Any
// user-defined task should derive from it and implement the Exec() virtual
// method.
//==============================================================================

#ifndef ROOT_TTask
#include "TTask.h"
#endif

#ifndef ROOT_TObjArray
#include "TObjArray.h"
#endif

class TClass;
class AliAnalysisDataSlot;
class AliAnalysisDataContainer;

class AliAnalysisTask : public TTask {
 public:
  enum EAnalysisTaskFlags {
    kTaskUsed    = BIT(14),
    kTaskZombie  = BIT(15),
    kTaskChecked = BIT(16)
  };   
  
  AliAnalysisTask();
  AliAnalysisTask(const char *name, const char *title);
  AliAnalysisTask(const AliAnalysisTask &task); 
  virtual ~AliAnalysisTask();
  
  // Assignment
  AliAnalysisTask& operator=(const AliAnalysisTask &task);
  // Conect inputs/outputs to data containers (by AliAnalysisModule)
  Bool_t                    ConnectInput(Int_t islot, AliAnalysisDataContainer *cont);
  Bool_t                    ConnectOutput(Int_t islot, AliAnalysisDataContainer *cont);
  // Check connectivity
  Bool_t                    AreSlotsConnected();
  // Check if data for all inputs is ready
  void                      CheckNotify(Bool_t init=kFALSE);
  // Check if there are illegal circular dependencies
  Bool_t                    CheckCircularDeps();
  // Getters
  Int_t                     GetNinputs() const  {return fNinputs;}
  Int_t                     GetNoutputs() const {return fNoutputs;}
  TObject                  *GetPublishedData() const {return fPublishedData;}
  AliAnalysisDataSlot      *GetInputSlot(Int_t islot) const  {return (AliAnalysisDataSlot*)fInputs->At(islot);}
  AliAnalysisDataSlot      *GetOutputSlot(Int_t islot) const {return (AliAnalysisDataSlot*)fOutputs->At(islot);}
  TClass                   *GetInputType(Int_t islot) const;
  TClass                   *GetOutputType(Int_t islot) const;
  // === USE THIS TO RETREIVE INPUT DATA AND STATICALLY CAST IT TO THE DECLARED TYPE
  TObject                  *GetInputData(Int_t islot) const;
  Bool_t                    IsOutputReady(Int_t islot) const {return fOutputReady[islot];}
  Bool_t                    IsChecked() const  {return TObject::TestBit(kTaskChecked);}
  Bool_t                    IsInitialized() const  {return fInitialized;}
  Bool_t                    IsReady() const  {return fReady;}
  Bool_t                    IsUsed() const   {return TObject::TestBit(kTaskUsed);}
  Bool_t                    IsZombie() const {return TObject::TestBit(kTaskZombie);}
  void                      PrintTask(Option_t *option="all", Int_t indent=0) const;
  void                      PrintContainers(Option_t *option="all", Int_t indent=0) const;
  void                      SetChecked(Bool_t flag=kTRUE) {TObject::SetBit(kTaskChecked,flag);}
  void                      SetUsed(Bool_t flag=kTRUE);
  void                      SetZombie(Bool_t flag=kTRUE) {TObject::SetBit(kTaskZombie,flag);}
  // Main task execution 
  //=== IMPLEMENT THIS !!! ==============================================
  virtual void              Exec(Option_t *option) = 0;
  //=====================================================================
  Bool_t                    HasExecuted() const {return fHasExecuted;}
  //=====================================================================
  // === OVERLOAD THIS IF YOU WANT TO DO SOMETHING WITH THE OUTPUT
  virtual void              Terminate(Option_t *option="");
  //=====================================================================
  
 protected:
  // Define the input/output slots (called by user in the ctor of the derived class)
  //=== CALL IN THE CONSTRUCTOR OF DERIVED CLASS TO DEFINE INPUTS/OUTPUTS ===
  void                      DefineInput(Int_t islot, TClass *type);
  void                      DefineOutput(Int_t islot, TClass *type);
  //=====================================================================
  
  //=====================================================================
  // === OVERLOAD THIS TO CONNECT TREE BRANCHES AT INPUT SLOTS. YOU
  // SHOULD DEFINE HERE ALSO THE OBJECTS TO BE CONNECTED TO YOUR OUTPUTS
  virtual void              Init(Option_t *option="");
  //=====================================================================
  
  // Post output data (called by Exec() when data is ready)
  //=== CALL IN EXEC() FOR EACH OUTPUT WHEN READY ===
  Bool_t                    PostData(Int_t iout, TObject *data, Option_t *option="");
  //=====================================================================
  
  // === USE THIS FIRST IN YOUR Init() TO CHECH IF A BRANCH IS ALREADY CONNECTED
  // TO SOME ADDRESS.
  char                     *GetBranchAddress(Int_t islot, const char *branch) const;
  // === CALL THIS AFTERWARDS IN Init() IF THE BRANCH ADDRESS IS NOT YET SET
  Bool_t                    SetBranchAddress(Int_t islot, const char *branch, void *address) const;
  //=====================================================================
  // === CALL THIS IN INIT IF THE OUTPUT IS TO BE WRITTEN AT OUTPUT IOUT
  void                      OpenFile(Int_t iout, const char *name, Option_t *option) const;
  
  Bool_t                    fReady;      // Flag if the task is ready
  Bool_t                    fInitialized; // True if Init() was called
  Int_t                     fNinputs;    // Number of inputs
  Int_t                     fNoutputs;   // Number of outputs
  Bool_t                   *fOutputReady; //[fNoutputs] Flags for output readyness
  TObject                  *fPublishedData; // !published data
  TObjArray                *fInputs;     // Array of input slots
  TObjArray                *fOutputs;    // Array of output slots
  
  ClassDef(AliAnalysisTask,2)  // Class describing an analysis task
};
#endif
