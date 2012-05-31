#ifndef ALITRIGSCHEDULER_H
#define ALITRIGSCHEDULER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 27/07/2009

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

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

class TObjArray;
class AliTrigScheduledEntry;
class AliTrigScheduler;

//______________________________________________________________________________
class AliTrigScheduledGroup : public TNamed {

private:
  Int_t                     fPriority;     // Group priority
  Int_t                     fDelay;        // Group max. delay.
  TObjArray                *fEntries;      // List of scheduled entries

private:
  AliTrigScheduledGroup(const AliTrigScheduledGroup &other);
  AliTrigScheduledGroup &operator=(const AliTrigScheduledGroup &other);

public:  
  AliTrigScheduledGroup();
  AliTrigScheduledGroup(const char *name, Int_t priority=0);
  virtual ~AliTrigScheduledGroup();
  
  void                      AddEntry(AliTrigScheduledEntry *entry);
  void                      FireUp(Int_t time);
  TObjArray                *GetScheduledEntries() const {return fEntries;}
  Int_t                     GetNentries() const;
  Int_t                     GetPriority() const {return fPriority;}
  Int_t                     GetDelay()    const {return fDelay;}
  virtual void              Print(Option_t *option) const;
  void                      RemoveEntry(AliTrigScheduledEntry *entry);  
  void                      SetPriority(Int_t priority) {fPriority = priority;}
    
  ClassDef(AliTrigScheduledGroup, 1) // Groups of scheduled response functions
};   

//==============================================================================
//
//   AliTrigScheduledSequence - A programable group sequence. Scheduled groups
//                      are owned and controlled by a trigger scheduler. They
//                      are fired-up in such a sequence. A sequence supports some
//                      default modes but can also be programed manually.
//
//==============================================================================

//______________________________________________________________________________
class AliTrigScheduledSequence : public TNamed {

public:
enum ESortingType {
  kDefault       = 0,
  kTimeInc       = 1,
  kTimeDec       = 2,
  kPriorityInc   = 3,
  kPriorityDec   = 4,
  kCustom        = 5
};  

private:
  AliTrigScheduler         *fScheduler;    // Scheduler to which the sequence applies
  Int_t                     fNgroups;      // Number of groups
  ESortingType              fType;         // Sorting type
  Int_t                    *fArray;        //[fNgroups] Array specifying the sequence

private:
  AliTrigScheduledSequence(const AliTrigScheduledSequence &other);
  AliTrigScheduledSequence &operator=(const AliTrigScheduledSequence &other);
  void                      SortArray(Int_t *array, Bool_t increasing);

public:
  AliTrigScheduledSequence();
  AliTrigScheduledSequence(const char *name, AliTrigScheduler *scheduler);
  virtual ~AliTrigScheduledSequence();
  
  Int_t                    *GetArray() const {return fArray;}
  Int_t                     GetNgroups() const {return fNgroups;}
  AliTrigScheduler         *GetScheduler() const {return fScheduler;}
  ESortingType              GetSortingType() const {return fType;}
  virtual void              Print(Option_t *option) const;
  void                      Sort(ESortingType type, Int_t *sequence=0);
  
  ClassDef(AliTrigScheduledSequence, 1)  // Class for a scheduled group sequence
};

//==============================================================================
//
//   AliTrigScheduler - Device response function scheduler. Every device has a
//                      scheduler, but the same scheduler can replay responses of
//                      several devices. A scheduler holds groups of scheduled 
//                      entries. The groups can be replayed in programable
//                      sequences. A default group and sequence are always created.
//
//==============================================================================

class AliTrigScheduledGroup;
class AliTrigScheduledSequence;

//______________________________________________________________________________
class AliTrigScheduler : public TNamed {

private:
  Int_t                     fNgroups;      // Number of scheduled groups (at least one)
  TObjArray                *fGroups;       // List of groups of response functions
  TObjArray                *fSequences;    // List of group replay sequences
  AliTrigScheduledSequence *fCurrentSequence; // Current group replay sequence

private:
  AliTrigScheduler(const AliTrigScheduler &other);
  AliTrigScheduler &operator=(const AliTrigScheduler &other);

public:
  AliTrigScheduler();
  AliTrigScheduler(const char *name);
  virtual ~AliTrigScheduler();
  
  void                      AddScheduledEntry(AliTrigScheduledEntry *entry, const char *togroup="default");
  AliTrigScheduledGroup    *AddGroup(const char *groupname);
  AliTrigScheduledGroup    *AddGroup(AliTrigScheduledGroup *group);
  AliTrigScheduledSequence *AddSequence(const char *seqname, AliTrigScheduledSequence::ESortingType type=AliTrigScheduledSequence::kDefault, 
                                        Int_t *sequence = 0);
  void                      FireUp(Int_t time);
  AliTrigScheduledSequence *GetCurrentSequence() const {return fCurrentSequence;}
  Int_t                     GetNgroups() const {return fNgroups;}
  TObjArray                *GetScheduledGroups() const {return fGroups;}
  AliTrigScheduledGroup    *GetScheduledGroup(Int_t i) const;
  AliTrigScheduledGroup    *GetScheduledGroup(const char *name) const;
  void                      SetCurrentSequence(AliTrigScheduledSequence *seq) {fCurrentSequence = seq;}
  void                      SetGroupPriority(const char *groupname, Int_t priority);
     
  ClassDef(AliTrigScheduler,1)  // Trigger scheduler class
};
#endif
