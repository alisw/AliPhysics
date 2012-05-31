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
#include "TCollection.h"

#include "AliVEvent.h"

ClassImp(AliAnalysisStatistics)

//______________________________________________________________________________
AliAnalysisStatistics::AliAnalysisStatistics(const AliAnalysisStatistics &other)
      :TNamed(other),
       fNinput(other.fNinput),
       fNprocessed(other.fNprocessed),
       fNfailed(other.fNfailed),
       fNaccepted(other.fNaccepted),
       fOfflineMask(other.fOfflineMask)
{
// Copy constructor.
}

//______________________________________________________________________________
AliAnalysisStatistics &AliAnalysisStatistics::operator=(const AliAnalysisStatistics &other)
{
// Assignment.
  if (&other == this) return *this;
  fNinput     = other.fNinput;
  fNprocessed = other.fNprocessed;
  fNfailed    = other.fNfailed;
  fNaccepted  = other.fNaccepted;
  fOfflineMask = other.fOfflineMask;
  return *this;
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
    current++;
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
   
