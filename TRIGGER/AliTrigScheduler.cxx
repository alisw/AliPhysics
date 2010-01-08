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
// Author: Andrei Gheata, 05/01/2010

#include "AliTrigScheduler.h"

#include <TMath.h>
#include <TList.h>
#include <TObjArray.h>
#include "AliTrigScheduledEntry.h"

ClassImp(AliTrigScheduledGroup)

//==============================================================================
//
//   AliTrigScheduledGroup - A group of scheduled entries that will simply be
//                           fired-up sequentially. The group delay in global time
//                           units is the latest start time of the contained
//                           entries. A group has a priority assigned by the
//                           owner scheduler object. Groups are fired-up according
//                           a programable sequence.
//
//==============================================================================

//______________________________________________________________________________
AliTrigScheduledGroup::AliTrigScheduledGroup()
                      :TNamed(),
                       fPriority(0),
                       fDelay(0),
                       fEntries(NULL)
{
// I/O constructor. 
}

//______________________________________________________________________________
AliTrigScheduledGroup::AliTrigScheduledGroup(const char *name, Int_t priority)
                      :TNamed(name,""),
                       fPriority(priority),
                       fDelay(0),
                       fEntries(new TObjArray())
{
// Constructor.
}

//______________________________________________________________________________
AliTrigScheduledGroup::~AliTrigScheduledGroup()
{
// Destructor.
  if (fEntries) delete fEntries;
}  

//______________________________________________________________________________
void AliTrigScheduledGroup::AddEntry(AliTrigScheduledEntry *entry)
{
// Add a scheduled entry to the group. There is no check if an entry was added twice !
  if (!fEntries) fEntries = new TObjArray();
  fEntries->Add(entry);
}

//______________________________________________________________________________
void AliTrigScheduledGroup::FireUp(Int_t time)
{
// Fire-up all entries in the group.
  Int_t nentries = GetNentries();
  AliTrigScheduledEntry *entry;
  for (Int_t i=0; i<nentries; i++) {
    entry = (AliTrigScheduledEntry*)fEntries->At(i);
    entry->FireUp(time);
  }
}    

//______________________________________________________________________________
Int_t AliTrigScheduledGroup::GetNentries() const
{
// Get number of scheduled entries in the group.
  return (fEntries)?fEntries->GetEntriesFast():0;
}

//______________________________________________________________________________
void AliTrigScheduledGroup::Print(Option_t *option) const
{
// Print the group content.
  Int_t nentries = GetNentries();
  printf("Group: %s containing %d entries.\n", GetName(), nentries);
  TString opt(option);
  opt.ToUpper();
  // Check if details are requested.
  if (!opt.Contains("D")) return;
  for (Int_t i=0; i<nentries; i++) {
    printf("   %d: ", i);
    fEntries->At(i)->Print(option);
  }   
}   

//______________________________________________________________________________
void AliTrigScheduledGroup::RemoveEntry(AliTrigScheduledEntry *entry)
{
// Remove an entry.
   if (!fEntries) return;
   fEntries->RecursiveRemove(entry);
   fEntries->Compress();
}
   
ClassImp(AliTrigScheduledSequence)

//==============================================================================
//
//   AliTrigScheduledSequence - A programable group sequence. Scheduled groups
//                      are owned and controlled by a trigger scheduler. They
//                      are fired-up in such a sequence. A sequence supports some
//                      default modes but can also be programed manually.
//
//==============================================================================

//______________________________________________________________________________
AliTrigScheduledSequence::AliTrigScheduledSequence()
                :TNamed(),
                 fScheduler(NULL),
                 fNgroups(0),
                 fType(AliTrigScheduledSequence::kDefault),
                 fArray(NULL)
{
// I/O constructor
}

//______________________________________________________________________________
AliTrigScheduledSequence::AliTrigScheduledSequence(const char *name, AliTrigScheduler *scheduler)
                :TNamed(name,""),
                 fScheduler(scheduler),
                 fNgroups(scheduler->GetNgroups()),
                 fType(AliTrigScheduledSequence::kDefault),
                 fArray(NULL)
{
// Constructor
  if (fNgroups) {
    fArray = new Int_t[fNgroups];
    for (Int_t i=0; i<fNgroups; i++) fArray[i] = i;
  }
}    
   
//______________________________________________________________________________
AliTrigScheduledSequence::~AliTrigScheduledSequence()
{
// Destructor.
  if (fArray) delete [] fArray;
}  

//______________________________________________________________________________
void AliTrigScheduledSequence::Print(Option_t *) const
{
// Print the sequence.
  printf("Sequence: %s scheduled by: %s\n", GetName(), fScheduler->GetName());
  printf("   type: %d  sequence: ", (Int_t)fType);
  for (Int_t i=0; i<fNgroups; i++) printf("%d ", fArray[i]);
  printf("\n");
}   

//______________________________________________________________________________
void AliTrigScheduledSequence::SortArray(Int_t *array, Bool_t increasing)
{
// Sort the internal sequence array.
  Int_t *ind = new Int_t[fNgroups];
  TMath::Sort(fNgroups, array, ind, !increasing);
  memcpy(fArray, ind, fNgroups*sizeof(Int_t));
  delete [] ind;
}
  
//______________________________________________________________________________
void AliTrigScheduledSequence::Sort(ESortingType type, Int_t *sequence)
{
  // Sort the group sequence according a predefined type. The sequence
  // custom input array is considered only for kCustom type.
  Int_t i;
  Int_t *array = 0;
  switch (type) {
    case AliTrigScheduledSequence::kDefault:
      // Just ID permutation
      for (i=0; i<fNgroups; i++) fArray[i] = i;
      break;
    case AliTrigScheduledSequence::kTimeInc:
      // Sort by increasing start time
      array = new Int_t[fNgroups];
      for (i=0; i<fNgroups; i++) {
        array[i] = fScheduler->GetScheduledGroup(i)->GetDelay();
      }  
      SortArray(array, kTRUE);
      break;
    case AliTrigScheduledSequence::kTimeDec:
      // Sort by decreasing start time
      array = new Int_t[fNgroups];
      for (i=0; i<fNgroups; i++) {
        array[i] = fScheduler->GetScheduledGroup(i)->GetDelay();
      }  
      SortArray(array, kFALSE);
      break;
    case AliTrigScheduledSequence::kPriorityInc:
      // Sort by increasing priority
      array = new Int_t[fNgroups];
      for (i=0; i<fNgroups; i++) {
        array[i] = fScheduler->GetScheduledGroup(i)->GetPriority();
      }  
      SortArray(array, kTRUE);
      break;
    case AliTrigScheduledSequence::kPriorityDec:
      // Sort by decreasing priority
      array = new Int_t[fNgroups];
      for (i=0; i<fNgroups; i++) {
        array[i] = fScheduler->GetScheduledGroup(i)->GetPriority();
      }  
      SortArray(array, kFALSE);
      break;
    case AliTrigScheduledSequence::kCustom:
      if (!sequence) {
        Error("Sort", "Sequence array must be provided for custom type");
        return;
      }
      memcpy(fArray, sequence, fNgroups*sizeof(Int_t));
  }
  if (array) delete [] array;  
}            

ClassImp(AliTrigScheduler)

//==============================================================================
//   AliTrigScheduler - Device response function scheduler. Every device has a
//                      scheduler, but the same scheduler can replay responses of
//                      several devices. A scheduler holds groups of scheduled 
//                      entries. The groups can be replayed in programable
//                      sequences. A default group and sequence are always created.
//==============================================================================

//______________________________________________________________________________
AliTrigScheduler::AliTrigScheduler()
                 :TNamed(),
                  fNgroups(0),
                  fGroups(NULL),
                  fSequences(NULL),
                  fCurrentSequence(NULL)
{
// I/O constructor.
}

//______________________________________________________________________________
AliTrigScheduler::AliTrigScheduler(const char *name)
                 :TNamed(name,""),
                  fNgroups(0),
                  fGroups(new TObjArray()),
                  fSequences(new TObjArray()),
                  fCurrentSequence(NULL)
{
// Constructor.
  AddGroup("default");
  AddSequence("default");
}

//______________________________________________________________________________
AliTrigScheduler::~AliTrigScheduler()
{
// Destructor.
  if (fGroups) {fGroups->Delete(); delete fGroups;}
  if (fSequences) {fSequences->Delete(); delete fSequences;}
}

//______________________________________________________________________________
void AliTrigScheduler::AddScheduledEntry(AliTrigScheduledEntry *entry, const char *togroup)
{
// Add a scheduled entry to a given group.
  AliTrigScheduledGroup *group = GetScheduledGroup(togroup);
  if (!group) {
    Error("AddScheduledEntry", "Group %s does not exist in scheduler %s", togroup, GetName());
    return;
  }
  group->AddEntry(entry);
}

//______________________________________________________________________________
AliTrigScheduledGroup *AliTrigScheduler::AddGroup(const char *groupname)
{
// Add a group to the list of groups.
  if (!fGroups) fGroups = new TObjArray();
  if (fGroups->FindObject(groupname)) {
    Error("AddGroup", "Scheduler %s contains already a group named: %s", GetName(), groupname);
    return NULL;
  }
  AliTrigScheduledGroup *group = new AliTrigScheduledGroup(groupname);
  fGroups->Add(group);
  return group;
}

//______________________________________________________________________________
AliTrigScheduledGroup *AliTrigScheduler::AddGroup(AliTrigScheduledGroup *group)
{
// Add a group to the list of groups.
  if (!fGroups) fGroups = new TObjArray();
  if (fGroups->FindObject(group->GetName())) {
    Error("AddGroup", "Scheduler %s contains already a group named: %s", GetName(), group->GetName());
    return NULL;
  }
  fGroups->Add(group);
  return group;
}

//______________________________________________________________________________
AliTrigScheduledSequence *AliTrigScheduler::AddSequence(const char *seqname, AliTrigScheduledSequence::ESortingType type, Int_t *sequence)
{
// Add a sequence to the scheduler. Becomes the current sequence.
  if (!fSequences) fSequences = new TObjArray();
  if (fSequences->FindObject(seqname)) {
    Error("AddSequence", "Scheduler %s contains already a sequence named: %s", GetName(), seqname);
    return NULL;
  }
  AliTrigScheduledSequence *seq = new AliTrigScheduledSequence(seqname, (AliTrigScheduler*)this); 
  seq->Sort(type, sequence);
  fCurrentSequence = seq;
  return seq;
}

//______________________________________________________________________________
void AliTrigScheduler::FireUp(Int_t time)
{
// Fire-up groups in the order given by the current sequence.
  if (!fCurrentSequence) Fatal("FireUp", "No scheduled sequence booked for scheduler: %s", GetName());
  Int_t *sequence = fCurrentSequence->GetArray();
  AliTrigScheduledGroup *group;
  for (Int_t i=0; i<fNgroups; i++) {
    group = GetScheduledGroup(sequence[i]);
    group->FireUp(time);
  }
}

//______________________________________________________________________________
AliTrigScheduledGroup *AliTrigScheduler::GetScheduledGroup(Int_t i) const
{
// Get i-th registered group (scheduling order does not matter, only group addition order).
  return (AliTrigScheduledGroup*)fGroups->At(i);
}

//______________________________________________________________________________
AliTrigScheduledGroup *AliTrigScheduler::GetScheduledGroup(const char *name) const
{
// Get a scheduled group by name.
  return (AliTrigScheduledGroup*)fGroups->FindObject(name);
}

//______________________________________________________________________________
void AliTrigScheduler::SetGroupPriority(const char *groupname, Int_t priority)
{
// Set the priority of a group.
  AliTrigScheduledGroup *group = GetScheduledGroup(groupname);
  if (!group) {
    Error("SetGroupPriority", "Scheduler %s has no group named: %s", GetName(), groupname);
    return;
  }
  group->SetPriority(priority);
}
  
