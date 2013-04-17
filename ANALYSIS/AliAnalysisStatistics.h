#ifndef ALIANALYSISSTATISTICS_H
#define ALIANALYSISSTATISTICS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 20/12/2010

//==============================================================================
//   AliAnalysisStatistics - Class holding statistics information regarding the
//      number of events processed, failed and accepted.
//==============================================================================

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class TObjArray;
class TStopwatch;

class AliAnalysisStatistics : public TNamed {

protected:
  Long64_t                    fNinput;            // Total number of input events
  Long64_t                    fNprocessed;        // Number of events processed
  Long64_t                    fNfailed;           // Number of events for which reading failed
  Long64_t                    fNaccepted;         // Number of events that passed filtering criteria
  UInt_t                      fOfflineMask;       // Offline mask used for accepted events
  Int_t                       fMaxTasks;          // Allocated size for the task timing arrays
  Int_t                       fNtasks;            // Number of tasks
  Int_t                       fCurrentTask;       // Current task being timed
  Double_t                   *fTaskTimeReal;      //[fNtasks] Cumulated CPU time per task
  Double_t                   *fTaskTimeCPU;       //[fNtasks] Cumulated CPU time per task
  TObjArray                  *fTaskNames;         // Task names
  TStopwatch                 *fTaskTimer;         //! Stopwatch for task timing
  
public:
  AliAnalysisStatistics() : TNamed(),fNinput(0),fNprocessed(0),fNfailed(0),
    fNaccepted(0),fOfflineMask(0), fMaxTasks(0),fNtasks(0), fCurrentTask(-1),
    fTaskTimeReal(0), fTaskTimeCPU(0), fTaskNames(0), fTaskTimer(0) {}
  AliAnalysisStatistics(const char *name) 
                          : TNamed(name,""),fNinput(0),fNprocessed(0),fNfailed(0),
    fNaccepted(0),fOfflineMask(0), fMaxTasks(0),fNtasks(0), fCurrentTask(-1),
    fTaskTimeReal(0), fTaskTimeCPU(0), fTaskNames(0), fTaskTimer(0) {}
  AliAnalysisStatistics(const AliAnalysisStatistics &other);
  virtual ~AliAnalysisStatistics() {}
  
  AliAnalysisStatistics& operator=(const AliAnalysisStatistics &other);
  // Update methods
  void                        AddInput(Int_t nevents=1)     {fNinput += nevents;}
  void                        AddProcessed(Int_t nevents=1) {fNprocessed += nevents;}
  void                        AddFailed(Int_t nevents=1)    {fNfailed += nevents;}
  void                        AddAccepted(Int_t nevents=1)  {fNaccepted += nevents;}
  // Getters
  Long64_t                    GetNinput()     const         {return fNinput;}
  Long64_t                    GetNprocessed() const         {return fNprocessed;}
  Long64_t                    GetNfailed()    const         {return fNfailed;}
  Long64_t                    GetNaccepted()  const         {return fNaccepted;}
  UInt_t                      GetOfflineMask() const        {return fOfflineMask;}
  static const char          *GetMaskAsString(UInt_t mask);
  Int_t                       GetNtasks() const             {return fNtasks;}
  const char                 *GetTaskName(Int_t itask) const;
  Double_t                    GetRealTime(Int_t itask) const {return fTaskTimeReal[itask];}
  Double_t                    GetCPUTime(Int_t itask) const {return fTaskTimeCPU[itask];}
  
  void                        SetOfflineMask(UInt_t mask)   {fOfflineMask = mask;}
  virtual Long64_t            Merge(TCollection* list);
  virtual void                Print(const Option_t *option="") const;
  // Task timing
  void                        StartTimer(Int_t itask, const char *name, const char *classname = "");
  void                        StopTimer();

  ClassDef(AliAnalysisStatistics,2)  // Class holding the processed events statistics
};
#endif
