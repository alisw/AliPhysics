/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */
// Author: Andrei Gheata, 20/12/2010

//==============================================================================
// AliAnalysisStatistics - basic class for storing statistics for the processed
//   events. The object is mergeable and can be used for general purpose. In case
//   a AliAnalysisTaskStat is used, this will set the global statistics object
//   to the analysis manager and will update it for the accepted events.
//==============================================================================

#include "AliAnalysisStatistics.h"

#include "Riostream.h"
#include "TObjArray.h"
#include "TStopwatch.h"

#include "AliVEvent.h"

using std::cout;
using std::endl;
ClassImp(AliAnalysisStatistics)

//______________________________________________________________________________
AliAnalysisStatistics::AliAnalysisStatistics(const AliAnalysisStatistics &other)
      :TNamed(other),
       fNinput(other.fNinput),
       fNprocessed(other.fNprocessed),
       fNfailed(other.fNfailed),
       fNaccepted(other.fNaccepted),
       fOfflineMask(other.fOfflineMask),
       fMaxTasks(other.fMaxTasks),
       fNtasks(other.fNtasks),
       fCurrentTask(other.fCurrentTask),
       fTaskTimeReal(0),
       fTaskTimeCPU(0),
       fTaskNames(0),
       fTaskTimer(0)
{
// Copy constructor.
  if (fNtasks) {
    fTaskTimer = new TStopwatch();
    fTaskTimeReal = new Double_t[fMaxTasks];
    memset(fTaskTimeReal, 0, fMaxTasks*sizeof(Double_t));
    memcpy(fTaskTimeReal, other.fTaskTimeReal, fNtasks*sizeof(Double_t));
    fTaskTimeCPU = new Double_t[fMaxTasks];
    memset(fTaskTimeCPU, 0, fMaxTasks*sizeof(Double_t));
    memcpy(fTaskTimeCPU, other.fTaskTimeCPU, fNtasks*sizeof(Double_t));
    fTaskNames = new TObjArray(fMaxTasks);
    for (Int_t i=0; i<fNtasks; i++) fTaskNames->AddAt(new TObjString(other.GetTaskName(i)), i);
  }
}

//______________________________________________________________________________
AliAnalysisStatistics &AliAnalysisStatistics::operator=(const AliAnalysisStatistics &other)
{
// Assignment.
  if (&other == this) return *this;
  fNinput       = other.fNinput;
  fNprocessed   = other.fNprocessed;
  fNfailed      = other.fNfailed;
  fNaccepted    = other.fNaccepted;
  fOfflineMask  = other.fOfflineMask;
  fMaxTasks     = other.fMaxTasks;
  fNtasks       = other.fNtasks;
  fCurrentTask  = other.fCurrentTask;
  fTaskTimeReal = 0;
  fTaskTimeCPU  = 0;
  fTaskNames    = 0;
  fTaskTimer   = 0;
  if (fNtasks) {
    fTaskTimer = new TStopwatch();
    fTaskTimeReal = new Double_t[fMaxTasks];
    memset(fTaskTimeReal, 0, fMaxTasks*sizeof(Double_t));
    memcpy(fTaskTimeReal, other.fTaskTimeReal, fNtasks*sizeof(Double_t));
    fTaskTimeCPU = new Double_t[fMaxTasks];
    memset(fTaskTimeCPU, 0, fMaxTasks*sizeof(Double_t));
    memcpy(fTaskTimeCPU, other.fTaskTimeCPU, fNtasks*sizeof(Double_t));
    fTaskNames = new TObjArray(fMaxTasks);
    for (Int_t i=0; i<fNtasks; i++) fTaskNames->AddAt(new TObjString(other.GetTaskName(i)), i);
  }  
  return *this;
}

//______________________________________________________________________________
void AliAnalysisStatistics::StartTimer(Int_t itask, const char *name, const char *classname)
{
// Measure the CPU time done by this task in the interval
  if (!fTaskTimer) {
    // Create the arrays with timings with the initial size
    if (!fMaxTasks) fMaxTasks = 100;
    fTaskTimer = new TStopwatch();
    fTaskTimeReal = new Double_t[fMaxTasks];
    memset(fTaskTimeReal, 0, fMaxTasks*sizeof(Double_t));
    fTaskTimeCPU = new Double_t[fMaxTasks];
    memset(fTaskTimeCPU, 0, fMaxTasks*sizeof(Double_t));
    fTaskNames = new TObjArray(fMaxTasks);
  } else {
  // Stop the timer if it was timing some task
    StopTimer();
  }  
  
  if (fNtasks<itask+1) {
  // Double the array size
    if (itask>=fMaxTasks) {
      Int_t newsize = TMath::Max(2*fMaxTasks, itask+1);
      Double_t *taskTimeReal = new Double_t[newsize];
      memset(taskTimeReal, 0, newsize*sizeof(Double_t));
      memcpy(taskTimeReal, fTaskTimeReal, fMaxTasks*sizeof(Double_t));
      delete [] fTaskTimeReal;
      fTaskTimeReal = taskTimeReal;
      Double_t *taskTimeCPU = new Double_t[newsize];
      memset(taskTimeCPU, 0, newsize*sizeof(Double_t));
      memcpy(taskTimeCPU, fTaskTimeCPU, fMaxTasks*sizeof(Double_t));
      delete [] fTaskTimeCPU;
      fTaskTimeCPU = taskTimeCPU;
      fMaxTasks = newsize;
    }  
    fNtasks = itask+1;
  }
  // Start the timer for the new task
  fCurrentTask = itask;
  if (!fTaskNames->At(fCurrentTask)) {
    TString sname = name;
    if (strlen(classname)) sname += Form("(%s)", classname);
    fTaskNames->AddAt(new TObjString(sname), fCurrentTask);
  }  
  fTaskTimer->Start(kTRUE);  
}   

//______________________________________________________________________________
void AliAnalysisStatistics::StopTimer()
{
// Stop the current task timing.
  if (fCurrentTask>=0) {
    fTaskTimer->Stop();
    fTaskTimeReal[fCurrentTask] += fTaskTimer->RealTime();
    fTaskTimeCPU[fCurrentTask]  += fTaskTimer->CpuTime();
    fCurrentTask = -1;
  }
}  

//______________________________________________________________________________
Long64_t AliAnalysisStatistics::Merge(TCollection* list)
{
// Merge statistics objets from list on top of this.
  TIter next(list);
  AliAnalysisStatistics *current;
  Long64_t count = 1;
  while ((current = (AliAnalysisStatistics*)next())) {
    fNinput     += current->GetNinput();
    fNprocessed += current->GetNprocessed();
    fNfailed    += current->GetNfailed();
    fNaccepted  += current->GetNaccepted();
    for (Int_t i=0; i<fNtasks; i++) {
      fTaskTimeReal[i] += current->GetRealTime(i);
      fTaskTimeCPU[i] += current->GetCPUTime(i);
    }   
  }
  return count;
}

//______________________________________________________________________________
void AliAnalysisStatistics::Print(const Option_t *) const
{
// Print info about the processed statistics.
  cout << "### Input events                 : " << fNinput << endl;
  cout << "### Processed events w/o errors  : " << fNprocessed << endl;
  cout << "### Failed events                : " << fNfailed << endl;
  cout << "### Accepted events for mask: " << GetMaskAsString(fOfflineMask) << ": " << fNaccepted << endl;
  if (fNtasks) {
    cout << "Timing per task:" <<endl;
    TString s;
    for (Int_t i=0; i<fNtasks; i++) {
      s = Form("%03d:  real: %9.2f   cpu: %9.2f   => %s", i,fTaskTimeReal[i], fTaskTimeCPU[i], GetTaskName(i));
      cout << s << endl;
    }
  }  
}

//______________________________________________________________________________
const char *AliAnalysisStatistics::GetMaskAsString(UInt_t mask)
{
// Returns a string corresponding to the offline mask.
   static TString smask;
   smask = "ALL EVT.";
   if (!mask) return smask.Data();
   smask.Clear();
   if (mask & AliVEvent::kMB)   smask = "MB";
   if (mask & AliVEvent::kMUON) {
      if (!smask.IsNull()) smask += " | ";
      smask += "MUON";
   }
   if (mask & AliVEvent::kHighMult) {
      if (!smask.IsNull()) smask += " | ";
      smask += "HighMult";
   }
   if (mask & AliVEvent::kUserDefined) {
      if (!smask.IsNull()) smask += " | ";
      smask += "UserDefined";
   }
   if (mask ==  AliVEvent::kAny) smask = "ANY";
   return smask.Data();
}

//______________________________________________________________________________
const char *AliAnalysisStatistics::GetTaskName(Int_t itask) const
{
// Returns task name
  if (!fTaskNames || !fTaskNames->At(itask)) return 0;
  return fTaskNames->At(itask)->GetName();
}
  
