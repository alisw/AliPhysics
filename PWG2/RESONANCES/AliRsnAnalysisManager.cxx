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

////////////////////////////////////////////////////////////////////////////////
//
//  This is the uppermost level in the analysis task and is built as a 
//  separate object since up to this level the execution of the analysis 
//  is not directly related to the real analysis task structure.
//  An analysis task for resonances will just configure and initialize an
//  instance of this class, which will do all the work, and will share with it
//  the TList containing all the outputs.
//
//  This object collects all AliRsnPair and AliRsnMonitor objects which do
//  the computations and fill histograms, each properly initialized with 
//  outputs (histos, ntuples) and cuts on whatever is needed.
//  
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TH1.h>
#include <TROOT.h>

#include "AliLog.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliRsnEvent.h"
#include "AliRsnPairFunctions.h"
#include "AliRsnPairNtuple.h"
#include "AliRsnMonitorFunctions.h"
#include "AliRsnMonitorNtuple.h"

#include "AliRsnAnalysisManager.h"


ClassImp(AliRsnAnalysisManager)

//_____________________________________________________________________________
AliRsnAnalysisManager::AliRsnAnalysisManager(const char *name) :
   TNamed(name, ""),
   fAddUsageHist(kFALSE),
   fList(0x0),
   fPairs(0),
   fMonitors(0)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnAnalysisManager::AliRsnAnalysisManager(const AliRsnAnalysisManager& copy) :
   TNamed(copy),
   fAddUsageHist(copy.fAddUsageHist),
   fList(copy.fList),
   fPairs(copy.fPairs),
   fMonitors(copy.fMonitors)
{
//
// Copy constructor
//
}

//_____________________________________________________________________________
AliRsnAnalysisManager& AliRsnAnalysisManager::operator=(const AliRsnAnalysisManager& copy)
{
//
// Assignment operator
//

   TNamed::operator=(copy);

   fAddUsageHist = copy.fAddUsageHist;
   fList = copy.fList;
   fPairs = copy.fPairs;
   fMonitors = copy.fMonitors;

   return (*this);
}

//_____________________________________________________________________________
void AliRsnAnalysisManager::Add(AliRsnPair *pair)
{
//
// Adds a new pair to the list.
//

   AliDebug(AliLog::kDebug + 2, "<-");

   if (!pair) {
      AliWarning(Form("AliRsnPairManager is %p. Skipping ...", pair));
      return;
   }

   AliDebug(AliLog::kDebug + 1, Form("Adding %s [%d]...", pair->GetName(), fPairs.GetEntries()));
   fPairs.Add(pair);

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnAnalysisManager::Add(AliRsnMonitor *monitor)
{
//
// Adds a new monitor to the list.
//

   AliDebug(AliLog::kDebug + 2, "<-");

   if (!monitor) {
      AliWarning(Form("AliRsnPairManager is %p. Skipping ...", monitor));
      return;
   }

   AliDebug(AliLog::kDebug + 1, Form("Adding %s [%d]...", monitor->GetName(), fMonitors.GetEntries()));
   fMonitors.Add(monitor);

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnAnalysisManager::Print(Option_t* /*dummy*/) const
{
//
// Overload of the TObject::Print() method
//

   AliInfo(Form("\t======== Analysis Manager %s ========", GetName()));
   PrintArray();
}

//_____________________________________________________________________________
void AliRsnAnalysisManager::PrintArray() const
{
//
// Calls the "Print" method of all included pair managers
//

   AliDebug(AliLog::kDebug + 2, "<-");

   AliRsnPair *pair = 0;
   TObjArrayIter nextPair(&fPairs);
   while ((pair = (AliRsnPair*)nextPair())) pair->Print();
   
   AliRsnMonitor *monitor = 0;
   TObjArrayIter nextMonitor(&fMonitors);
   while ((monitor = (AliRsnMonitor*)nextMonitor())) monitor->Print();

   AliDebug(AliLog::kDebug + 2, "->");
}

//_____________________________________________________________________________
void AliRsnAnalysisManager::InitAllPairs(TList *list)
{
//
// Initialize all pair managers, and put all the TList of histograms
// generated by each one into a unique final output TList
//

   Int_t i = 0;

   // pairs
   AliRsnPair   *pair = 0;
   TObjArrayIter nextPair(&fPairs);
   while ((pair = (AliRsnPair*)nextPair())) {
      AliDebug(AliLog::kDebug + 1, Form("InitAllPairs of the PairManager(%s) [%d] ...", pair->GetName(), i++));
      pair->Init(GetName(), list);

      // add a counter for used/unused events for each pair
      if (fAddUsageHist) {
         TH1I *hPairUsed = new TH1I(Form("%s_%s_USED", GetName(), pair->GetName()), "Used events for pair", 2, 0, 2);
         list->Add(hPairUsed);
      }
   }
   
   // monitors
   i = 0;
   AliRsnMonitor *monitor = 0;
   TObjArrayIter  nextMonitor(&fMonitors);
   while ((monitor = (AliRsnMonitor*)nextMonitor())) {
      AliDebug(AliLog::kDebug + 1, Form("InitAllPairs of the PairManager(%s) [%d] ...", monitor->GetName(), i++));
      monitor->Init(GetName(), list);
   }

   // set list pointer
   fList = list;
   fList->Print();
}

//_____________________________________________________________________________
void AliRsnAnalysisManager::ProcessAll(AliRsnEvent *ev0, AliRsnEvent *ev1, Bool_t pureMC)
{
//
// Loop on all pairs/monitors stored here 
// and process all candidate daughters.
//

   AliDebug(AliLog::kDebug + 2, "<-");
   
   // don't do anything if the output list isn't initialized
   if (!fList) return;
   
   // if second event is NULL, we assume that analysis is on single event
   Bool_t sameEvent = kFALSE;
   if (!ev1) {
      ev1 = ev0;
      sameEvent = kTRUE;
   }
   
   // if MC reference is absent, cannot process pure MC
   if (pureMC && (!ev0->GetRefMC() || !ev1->GetRefMC())) {
      AliError("Cannot process a pure MC analysis without MC references!");
      return;
   }

   // count total number of candidates per event
   // (sum of tracks, V0s and cascades)
   Int_t nTot[2];
   if (pureMC) {
      nTot[0] = ev0->GetRefMC()->GetNumberOfTracks();
      nTot[1] = ev1->GetRefMC()->GetNumberOfTracks();
   } else {
      nTot[0] = ev0->GetAbsoluteSum();
      nTot[1] = ev1->GetAbsoluteSum();
   }
   
   // if there are only monitors and no pairs
   // disable the inner loop by settin counter to zero
   if (fPairs.IsEmpty()) nTot[1] = 0;

   // variables
   Int_t           i0, i1, i;
   AliRsnDaughter  daughter0, daughter1;
   AliRsnPair     *pair = 0x0;
   AliRsnMonitor  *monitor = 0x0;
   TObjArrayIter   nextPair(&fPairs);
   TObjArrayIter   nextMonitor(&fMonitors);

   // reset all counters which tell us
   // how many entries were added now
   while ((pair = (AliRsnPair*)nextPair())) {
      pair->ResetCount();
      pair->GetMother()->SetRefEvent(daughter0.GetOwnerEvent());
   }

   // external loop
   for (i0 = 0; i0 < nTot[0]; i0++) {
      
      // try to assign first track
      // and check global cuts
      // in case of ESD pure MC, skip not physical primaries
      if (pureMC) {
         if (ev0->IsESD()) if (!ev0->GetRefMCESD()->Stack()->IsPhysicalPrimary(i0)) continue;
         if (!ev0->SetDaughterMC(daughter0, i0)) continue;
      } else {
         if (!ev0->SetDaughterAbs(daughter0, i0)) continue;
      }
      // when the assigment is unsuccessful, this i known in internal status flag
      if (!daughter0.IsOK() || !daughter0.GetRef()) continue;
      
      // process all single-track monitors
      nextMonitor.Reset();
      while ((monitor = (AliRsnMonitor*)nextMonitor())) 
         if (monitor->Fill(&daughter0)) 
            monitor->Compute();

      // internal loop
      // starts from next track or first, depending if not mixing or yes
      for (i1 = (sameEvent ? i0 + 1 : 0); i1 < nTot[1]; i1++) {

         // try to assign first track
         // and check global cuts
         // in case of ESD pure MC, skip not physical primaries
         if (pureMC) {
            if (ev1->IsESD()) if (!ev1->GetRefMCESD()->Stack()->IsPhysicalPrimary(i1)) continue;
            if (!ev1->SetDaughterMC(daughter1, i1)) continue;
         } else {
            if (!ev1->SetDaughterAbs(daughter1, i1)) continue;
         }
         // when the assigment is unsuccessful, this i known in internal status flag
         if (!daughter1.IsOK() || !daughter1.GetRef()) continue;

         // loop over all pairs and make computations
         i = 0;
         nextPair.Reset();
         while ((pair = (AliRsnPair*)nextPair())) {
            
            // debug message
            AliDebug(AliLog::kDebug + 1, Form("ProcessAllPairs of the AnalysisManager(%s) [%d] ...", pair->GetName(), i++));

            // process the two tracks
            // the Fill() method returns kTRUE if the two daughters
            // do match the pair definition of each AliRsnPait object,
            // and the pair is processed only in this case
            if (pair->Fill(&daughter0, &daughter1, kTRUE)) {
               pair->Compute();
            } else if (pair->Fill(&daughter1, &daughter0, kFALSE)) {
               pair->Compute();
            }
         }
      }
   }

   // update all count histograms counters
   if (fAddUsageHist) {
      nextPair.Reset();
      while ((pair = (AliRsnPair*)nextPair())) {
         TH1I *hist = (TH1I*)fList->FindObject(Form("%s_%s_USED", GetName(), pair->GetName()));
         if (!hist) continue;
         if (pair->GetCount() > 0) hist->Fill(1); else hist->Fill(0);
      }
   }

   AliDebug(AliLog::kDebug + 2, "->");
}
