/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTTriggerDomain.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Nov 2008
/// @brief  Implementation of the AliHLTTriggerDomain class.
///
/// The trigger domain class is the set of HLT raw data block types that should
/// be readout and sent to HLTOUT.

#include "AliHLTTriggerDomain.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTReadoutList.h"
#include "Riostream.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "AliDAQ.h"

ClassImp(AliHLTTriggerDomain)


AliHLTTriggerDomain::AliHLTTriggerDomain() :
  TObject(), fEntries(AliHLTDomainEntry::Class(), 10)
{
  // Default constructor.
}


AliHLTTriggerDomain::AliHLTTriggerDomain(const char* list) :
  TObject(), fEntries(AliHLTDomainEntry::Class(), 10)
{
  // Constructs the domain from a list of entries.
  
  TString lst = list;
  TObjArray* entries = lst.Tokenize(",");
  for (Int_t i = 0; i < entries->GetEntriesFast(); i++)
  {
    TString entry = static_cast<TObjString*>(entries->UncheckedAt(i))->GetString();
    TObjArray* domainStrings = entry.Tokenize(":");
    if (domainStrings->GetEntriesFast() <= 0 or domainStrings->GetEntriesFast() > 3)
    {
      Error("AliHLTTriggerDomain",
            "The domain string must contain 1, 2 or 3 fields separated by a ':'."
           );
      delete domainStrings;
      continue;
    }
    
    bool inclusiveEntry = true;
    TString typeString = "*******";
    if (domainStrings->GetEntriesFast() >= 1)
    {
      typeString = static_cast<TObjString*>(domainStrings->UncheckedAt(0))->GetString();
      if (typeString.Length() > 0)
      {
        if (typeString[0] == '+')
        {
          inclusiveEntry = true;
          typeString.Remove(0, 1);
        }
        if (typeString[0] == '-')
        {
          inclusiveEntry = false;
          typeString.Remove(0, 1);
        }
      }
    }
    TString originString = "***";
    if (domainStrings->GetEntriesFast() >= 2)
    {
      originString = static_cast<TObjString*>(domainStrings->UncheckedAt(1))->GetString();
    }
    bool usespec = false;
    UInt_t spec = 0;
    if (domainStrings->GetEntriesFast() == 3)
    {
      TString specString = static_cast<TObjString*>(domainStrings->UncheckedAt(2))->GetString();
      char* error = NULL;
      spec = UInt_t( strtoul(specString.Data(), &error, 0) );
      if (error == NULL or *error != '\0')
      {
        Error("AliHLTTriggerDomain",
              "The last field of the domain string must be a number, but we received '%s'.",
              specString.Data()
             );
      }
      else
      {
        usespec = true;
      }
    }
    
    if (usespec)
    {
      if (inclusiveEntry)
        Add(typeString.Data(), originString.Data(), spec);
      else
        Remove(typeString.Data(), originString.Data(), spec);
    }
    else
    {
      if (inclusiveEntry)
        Add(typeString.Data(), originString.Data());
      else
        Remove(typeString.Data(), originString.Data());
    }
    
    delete domainStrings;
  }
  delete entries;
}


AliHLTTriggerDomain::AliHLTTriggerDomain(const AliHLTReadoutList& list) :
  TObject(), fEntries(AliHLTDomainEntry::Class(), 10)
{
  // Constructor creates a trigger domain from a readout list.
  // See header file for more details.
  
  Add(list);
}


AliHLTTriggerDomain::AliHLTTriggerDomain(const AliHLTTriggerDomain& domain) :
  TObject(domain),
  fEntries(AliHLTDomainEntry::Class(), domain.fEntries.GetEntriesFast())
{
  // Copy constructor performs a deep copy.
  // See header file for more details.
  
  for (Int_t i = 0; i < domain.fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* entry = static_cast<const AliHLTDomainEntry*>( domain.fEntries[i] );
    new (fEntries[fEntries.GetEntriesFast()]) AliHLTDomainEntry(*entry);
  }
}


AliHLTTriggerDomain::~AliHLTTriggerDomain()
{
  // Default destructor.
}


void AliHLTTriggerDomain::Add(const AliHLTReadoutList& list)
{
  // Adds the readout list to the trigger domain.
  // See header file for more details.
  
  Int_t detId[AliDAQ::kNDetectors] = {
      AliHLTReadoutList::kITSSPD, AliHLTReadoutList::kITSSDD, AliHLTReadoutList::kITSSSD,
      AliHLTReadoutList::kTPC, AliHLTReadoutList::kTRD, AliHLTReadoutList::kTOF,
      AliHLTReadoutList::kHMPID, AliHLTReadoutList::kPHOS, AliHLTReadoutList::kCPV,
      AliHLTReadoutList::kPMD, AliHLTReadoutList::kMUONTRK, AliHLTReadoutList::kMUONTRG,
      AliHLTReadoutList::kFMD, AliHLTReadoutList::kT0, AliHLTReadoutList::kV0,
      AliHLTReadoutList::kZDC, AliHLTReadoutList::kACORDE, AliHLTReadoutList::kTRG,
      AliHLTReadoutList::kEMCAL, AliHLTReadoutList::kDAQTEST, AliHLTReadoutList::kHLT
    };
  
  for (Int_t deti = 0; deti < AliDAQ::kNDetectors; deti++)
  {
    if (list.DetectorEnabled(detId[deti]))
    {
      Add("DAQRDOUT", AliDAQ::OnlineName(deti));
    }
    else
    {
      for (Int_t i = 0; i < AliDAQ::NumberOfDdls(deti); i++)
      {
        Int_t ddlId = AliDAQ::DdlID(deti, i);
        if (list.IsDDLEnabled(ddlId)) Add("DAQRDOUT", AliDAQ::OnlineName(deti), ddlId);
      }
    }
  }
}

void AliHLTTriggerDomain::Add(const AliHLTDomainEntry& entry)
{
  // Adds a new domain entry to the trigger domain.
  // See header file for more details.
  
  AliHLTDomainEntry intersect;
  bool alreadyInSet = false;
  
  // Get the initial size of the fEntries array since we might add things to the
  // end during the calculation.
  Int_t count = fEntries.GetEntriesFast();
  
  // Go through each entry that is already in fEntries and see if we can remove
  // it because it will become redundant, or if we need to patch exclusion entries
  // by adding inclusive intersects, or if we do not even need to add the new entry
  // because it is already part of the trigger domain.
  for (Int_t i = 0; i < count; i++)
  {
    AliHLTDomainEntry* ientry = static_cast<AliHLTDomainEntry*>(fEntries[i]);
    if (ientry->Inclusive())
    {
      if (entry.SubsetOf(*ientry))
      {
        alreadyInSet = true;
      }
      else if (ientry->SubsetOf(entry))
      {
        ientry->SetBit(14, true);  // mark for removal.
      }
    }
    else
    {
      if (ientry->SubsetOf(entry))
      {
        ientry->SetBit(14, true);  // mark for removal.
      }
      else if (entry.SubsetOf(*ientry))
      {
        alreadyInSet = false;
      }
      else if (ientry->IntersectWith(entry, intersect))
      {
        MarkForDeletionSubsetsOf(intersect, count);
        new (fEntries[fEntries.GetEntriesFast()]) AliHLTDomainEntry(kFALSE, intersect);
      }
    }
  }
  
  // Check if we need to add the new entry.
  if (not alreadyInSet)
  {
    MarkForDeletionSubsetsOf(entry, count);
    new (fEntries[fEntries.GetEntriesFast()]) AliHLTDomainEntry(kFALSE, entry);
  }
  RemoveMarkedEntries();
}


void AliHLTTriggerDomain::Add(const AliHLTComponentDataType& datatype)
{
  // Adds a new domain entry with the given data type to the trigger domain.
  // But the data block specification is set to the any matching wild card.
  // See header file for more details.
  
  Add(AliHLTDomainEntry(datatype));
}


void AliHLTTriggerDomain::Add(const char* blocktype, const char* origin)
{
  // Adds a new domain entry with the given data type and origin to the trigger domain.
  // But the data block specification is set to the any matching wild card.
  // See header file for more details.
  
  Add(AliHLTDomainEntry(blocktype, origin));
}


void AliHLTTriggerDomain::Add(const AliHLTComponentDataType& datatype, UInt_t spec)
{
  // Adds a new domain entry to the trigger domain with the data type and data block
  // specification bits.
  // See header file for more details.
  
  Add(AliHLTDomainEntry(datatype, spec));
}


void AliHLTTriggerDomain::Add(const char* blocktype, const char* origin, UInt_t spec)
{
  // Adds a new domain entry to the trigger domain with the given data type, origin
  // and data block specification bits.
  // See header file for more details.
  
  Add(AliHLTDomainEntry(blocktype, origin, spec));
}


void AliHLTTriggerDomain::Remove(const AliHLTReadoutList& list)
{
  // Removes the entries in the readout list from the trigger domain that are enabled.
  // See header file for more details.
  
  Int_t detId[AliDAQ::kNDetectors] = {
      AliHLTReadoutList::kITSSPD, AliHLTReadoutList::kITSSDD, AliHLTReadoutList::kITSSSD,
      AliHLTReadoutList::kTPC, AliHLTReadoutList::kTRD, AliHLTReadoutList::kTOF,
      AliHLTReadoutList::kHMPID, AliHLTReadoutList::kPHOS, AliHLTReadoutList::kCPV,
      AliHLTReadoutList::kPMD, AliHLTReadoutList::kMUONTRK, AliHLTReadoutList::kMUONTRG,
      AliHLTReadoutList::kFMD, AliHLTReadoutList::kT0, AliHLTReadoutList::kV0,
      AliHLTReadoutList::kZDC, AliHLTReadoutList::kACORDE, AliHLTReadoutList::kTRG,
      AliHLTReadoutList::kEMCAL, AliHLTReadoutList::kDAQTEST, AliHLTReadoutList::kHLT
    };
  
  for (Int_t deti = 0; deti < AliDAQ::kNDetectors; deti++)
  {
    if (list.DetectorEnabled(detId[deti]))
    {
      Remove("DAQRDOUT", AliDAQ::OnlineName(deti));
    }
    else
    {
      for (Int_t i = 0; i < AliDAQ::NumberOfDdls(deti); i++)
      {
        Int_t ddlId = AliDAQ::DdlID(deti, i);
        if (list.IsDDLEnabled(ddlId)) Remove("DAQRDOUT", AliDAQ::OnlineName(deti), ddlId);
      }
    }
  }
}


void AliHLTTriggerDomain::Remove(const AliHLTDomainEntry& entry)
{
  // Removes the given domain entry from the trigger domain.
  // See header file for more details.

  AliHLTDomainEntry intersect;
  bool addToExcludeSet = false;
  
  // Get the initial size of the fEntries array since we might add things to the
  // end during the calculation.
  Int_t count = fEntries.GetEntriesFast();
  
  // We need to go through all existing entries and see if they need to be removed
  // because they would become redundant when we add the new 'entry' to the end of
  // the fEntries list. We also need to check if the new entry needs to be added
  // at all because the trigger domain might already not contain those entries.
  // Lastly, some intersection entries might need to be added to patch up existing
  // inclusive trigger domain entries (rules / patterns).
  for (Int_t i = 0; i < count; i++)
  {
    AliHLTDomainEntry* ientry = static_cast<AliHLTDomainEntry*>(fEntries[i]);
    if (ientry->Inclusive())
    {
      if (ientry->SubsetOf(entry))
      {
        ientry->SetBit(14, true);  // mark for removal.
      }
      else if (entry.SubsetOf(*ientry))
      {
        addToExcludeSet = true;
      }
      else if (ientry->IntersectWith(entry, intersect))
      {
        new (fEntries[fEntries.GetEntriesFast()]) AliHLTDomainEntry(kTRUE, intersect);
      }
    }
    else
    {
      if (entry.SubsetOf(*ientry))
      {
        addToExcludeSet = false;
      }
      else if (ientry->SubsetOf(entry))
      {
        ientry->SetBit(14, true);  // mark for removal.
      }
    }
  }
  
  // Check if we need to add the new entry.
  if (addToExcludeSet)
  {
    MarkForDeletionSubsetsOf(entry, count);
    new (fEntries[fEntries.GetEntriesFast()]) AliHLTDomainEntry(kTRUE, entry);
  }
  RemoveMarkedEntries();
}


void AliHLTTriggerDomain::Remove(const AliHLTComponentDataType& datatype)
{
  // Removes the domain entries that have the given data type from the trigger domain.
  // See header file for more details.
  
  Remove(AliHLTDomainEntry(datatype));
}


void AliHLTTriggerDomain::Remove(const char* blocktype, const char* origin)
{
  // Removes the domain entries that have the given data type and origin from the
  // trigger domain.
  // See header file for more details.
  
  Remove(AliHLTDomainEntry(blocktype, origin));
}


void AliHLTTriggerDomain::Remove(const AliHLTComponentDataType& datatype, UInt_t spec)
{
  // Removes the domain entries that have the given data type and data block
  // specification bits from the trigger domain.
  // See header file for more details.
  
  Remove(AliHLTDomainEntry(datatype, spec));
}


void AliHLTTriggerDomain::Remove(const char* blocktype, const char* origin, UInt_t spec)
{
  // Removes the domain entries that have the given data type, origin and data
  // block specification bits from the trigger domain.
  // See header file for more details.
  
  Remove(AliHLTDomainEntry(blocktype, origin, spec));
}


bool AliHLTTriggerDomain::Contains(const AliHLTDomainEntry& entry) const
{
  // Checks to see if the given domain entry is part of the trigger domain set.
  // See header file for more details.

  // Simply go through the whole list of fEntries and for each entry see if the
  // given domain entry 'entry' being checked matches. If there is a match then
  // update the result depending on the entry type. i.e. set to false if the entry
  // in fEntries is an exclusion and set to true if it is an inclusion.
  bool result = false;
  for (Int_t i = 0; i < fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* ientry = static_cast<const AliHLTDomainEntry*>(fEntries[i]);
    if (ientry->Inclusive())
    {
      if (*ientry == entry) result = true;
    }
    else
    {
      if (entry.SubsetOf(*ientry)) result = false;
    }
  }
  return result;
}


bool AliHLTTriggerDomain::IncludeInReadout(const AliHLTComponentBlockData* block) const
{
  // Checks to see if the given data block is part of the trigger domain set and
  // should be readout.
  // See header file for more details.

  // Same algorithm as for Contains() but applied directly to the data block
  // descriptor structure.
  AliHLTDomainEntry blockEntry(block->fDataType, block->fSpecification);
  bool result = false;
  for (Int_t i = 0; i < fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* entry = static_cast<const AliHLTDomainEntry*>(fEntries[i]);
    if (entry->Inclusive())
    {
      if (*entry == block) result = true;
    }
    else
    {
      if (blockEntry.SubsetOf(*entry)) result = false;
    }
  }
  return result;
}


void AliHLTTriggerDomain::Clear(Option_t* option)
{
  // Clears the trigger domain (Removes all entries).
  
  fEntries.Clear(option);
}


void AliHLTTriggerDomain::Print(Option_t* /*option*/) const
{
  // Prints the trigger domain entries in the order that they are applied.
  // See header file for more details.

  cout << "Trigger domain rules (applied in order of first to last):" << endl;
  for (Int_t i = 0; i < fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* entry = static_cast<const AliHLTDomainEntry*>( fEntries[i] );
    if (entry->Inclusive())
    {
      cout << "Include ";
    }
    else
    {
      cout << "Exclude ";
    }
    entry->Print();
  }
  if (fEntries.GetEntriesFast() == 0)
  {
    cout << "(empty)" << endl;
  }
}


AliHLTTriggerDomain& AliHLTTriggerDomain::operator = (const AliHLTTriggerDomain& domain)
{
  // Assignment operator performs a deep copy.
  // See header file for more details.
  
  if (this == &domain) return *this;
  TObject::operator = (domain);
  fEntries.Clear();
  for (Int_t i = 0; i < domain.fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* entry = static_cast<const AliHLTDomainEntry*>( domain.fEntries[i] );
    new (fEntries[fEntries.GetEntriesFast()]) AliHLTDomainEntry(*entry);
  }
  return *this;
}


AliHLTTriggerDomain& AliHLTTriggerDomain::operator |= (const AliHLTTriggerDomain& domain)
{
  // This operator performs the set union.
  // See header file for more details.
  
  // Note that we partition the fEntries array into 3 regions for this calculation.
  //   - 0..entriesCount-1 : contains the initial entries of this trigger domain.
  //   - entriesCount..startOfIntersects-1 : is space reserved for the new entries
  //       from 'domain'.
  //   - startOfIntersects..fEntries.GetEntriesFast()-1 : This will grow as we add
  //       all the new domain intersections created during the calculation.
  //
  // Get the number of entries now before we start adding more entries from 'domain'.
  Int_t count = fEntries.GetEntriesFast();
  // Mark the start location for new intersection entries.
  Int_t startOfIntersects = count + domain.fEntries.GetEntriesFast();
  Int_t newIndex = startOfIntersects;

  // Allocate and initialise a single block of memory so that we do not call new twice.
  bool* buffer = new bool[startOfIntersects];
  for (Int_t i = 0; i < startOfIntersects; i++) buffer[i] = false;
  bool* removeThisEntry = buffer;
  bool* removeDomainEntry = buffer + count;

  AliHLTDomainEntry intersect;
  
  // The idea behind this algorithm is that we need to add all inclusion domain
  // entries from 'domain' to this object that will not be redundant, but for
  // the exclusion entries we patch the fEntries rule set by adding the appropriate
  // intersections to the end of fEntries.
  for (Int_t i = 0; i < domain.fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* newEntry = static_cast<const AliHLTDomainEntry*>( domain.fEntries[i] );
    for (Int_t j = 0; j < count; j++)
    {
      const AliHLTDomainEntry* currentEntry = static_cast<const AliHLTDomainEntry*>( fEntries[j] );
      if (currentEntry->Inclusive() and newEntry->Inclusive())
      {
        // If either entry is a subset of the other then we do not need to add
        // both, so make sure to remove the one that is redundant.
        if (newEntry->SubsetOf(*currentEntry))
        {
          removeDomainEntry[i] = true;
        }
        else if (currentEntry->SubsetOf(*newEntry))
        {
          removeThisEntry[j] = true;
        }
      }
      else
      {
        if (newEntry->IntersectWith(*currentEntry, intersect))
        {
          // We can remove all intersections that were already added that will
          // become redundant when this intersection is added to fEntries.
          MarkForDeletionSubsetsOf(intersect, startOfIntersects);
          
          // Make the new intersection entry an exclusion if the newEntry and
          // currentEntry flags are the same.
          bool exclude = newEntry->Exclusive() == currentEntry->Exclusive();
          new (fEntries[newIndex++]) AliHLTDomainEntry(exclude, intersect);
          
          // We can also remove entries that are subsets of another entry in the
          // opposite list, since they will be redundant when everything is merged
          // together. For example, remove entry x from fEntries if it is a subset
          // of entry y in domain.fEntries.
          if (currentEntry->IdenticalTo(intersect)) removeThisEntry[j] = true;
          if (newEntry->IdenticalTo(intersect)) removeDomainEntry[i] = true;
        }
      }
    }
  }

  MergeEntries(removeThisEntry, count, removeDomainEntry, startOfIntersects, domain);
  delete [] buffer;
  Optimise();
  return *this;
}


AliHLTTriggerDomain& AliHLTTriggerDomain::operator ^= (const AliHLTTriggerDomain& domain)
{
  // This operator performs the set union, less the set intersect (something like and xor).
  // See header file for more details.
  
  // Note that we partition the fEntries array into 3 regions for this calculation.
  //   - 0..entriesCount-1 : contains the initial entries of this trigger domain.
  //   - entriesCount..startOfIntersects-1 : is space reserved for the new entries
  //       from 'domain'.
  //   - startOfIntersects..fEntries.GetEntriesFast()-1 : This will grow as we add
  //       all the new domain intersections created during the calculation.
  //
  // Get the number of entries now before we start adding more entries from 'domain'.
  Int_t count = fEntries.GetEntriesFast();
  // Mark the start location for new intersection entries.
  Int_t startOfIntersects = count + domain.fEntries.GetEntriesFast();
  Int_t newIndex = startOfIntersects;

  // Allocate and initialise a single block of memory so that we do not call new twice.
  bool* buffer = new bool[startOfIntersects];
  for (Int_t i = 0; i < startOfIntersects; i++) buffer[i] = false;
  bool* removeThisEntry = buffer;
  bool* removeDomainEntry = buffer + count;

  AliHLTDomainEntry intersect;
  
  // This algorithm is similar to the case for the set union (operator |=), except
  // that we make sure to remove from the trigger domain all parts where the entries
  // from fEntries and domain.fEntries intersect.
  // This is done by adding the intersections to the end of fEntries such that they
  // effectively remove those overlapping trigger domain entries when calculating
  // IncludeInReadout() or Contains().
  for (Int_t i = 0; i < domain.fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* newEntry = static_cast<const AliHLTDomainEntry*>( domain.fEntries[i] );
    for (Int_t j = 0; j < count; j++)
    {
      const AliHLTDomainEntry* currentEntry = static_cast<const AliHLTDomainEntry*>( fEntries[j] );
      if (newEntry->IntersectWith(*currentEntry, intersect))
      {
        // We can remove all intersections that were already added that will
        // become redundant when this intersection is added to fEntries.
        MarkForDeletionSubsetsOf(intersect, startOfIntersects);
        
        // Make the new intersection entry an exclusion if the newEntry and
        // currentEntry flags are the same.
        bool exclude = newEntry->Exclusive() == currentEntry->Exclusive();
        new (fEntries[newIndex++]) AliHLTDomainEntry(exclude, intersect);
        
        // We can also remove entries that are subsets of another entry in the
        // opposite list, since they will be redundant when everything is merged
        // together. For example, remove entry x from fEntries if it is a subset
        // of entry y in domain.fEntries.
        if (currentEntry->IdenticalTo(intersect)) removeThisEntry[j] = true;
        if (newEntry->IdenticalTo(intersect)) removeDomainEntry[i] = true;
      }
    }
  }

  MergeEntries(removeThisEntry, count, removeDomainEntry, startOfIntersects, domain);
  delete [] buffer;
  Optimise();
  return *this;
}


AliHLTTriggerDomain& AliHLTTriggerDomain::operator -= (const AliHLTTriggerDomain& domain)
{
  // This operator performs the set difference.
  // See header file for more details.
  
  // Mark the number of entries in fEntries now before we start adding more
  // entries from 'domain' or intersections.
  Int_t startOfIntersects = fEntries.GetEntriesFast();
  Int_t newIndex = startOfIntersects;
  
  AliHLTDomainEntry intersect;
  
  // To compute the set difference we need to remove all all parts that overlap
  // with 'domain'. i.e. we need to find all the intersects between the domain
  // entries in fEntries and those in domain.fEntries, and add the intersects
  // to the fEntries list, such that they will cancel or remove the overlapping
  // parts of the two trigger domains.
  for (Int_t i = 0; i < domain.fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* checkEntry = static_cast<const AliHLTDomainEntry*>( domain.fEntries[i] );
    if (checkEntry->Inclusive())
    {
      // For inclusive entries we need to find the overlaps with the inclusive
      // entries in fEntries and add exclusive entries that will remove that
      // part of the trigger domain set.
      for (Int_t j = 0; j < startOfIntersects; j++)
      {
        AliHLTDomainEntry* currentEntry = static_cast<AliHLTDomainEntry*>( fEntries[j] );
        
        // We only need to consider the case where both entries are inclusive,
        // since an exclusion in fEntries already eliminates those data blocks
        // from the trigger domain set.
        if (currentEntry->Exclusive()) continue;
        
        if (checkEntry->IntersectWith(*currentEntry, intersect))
        {
          // We can remove all intersections that were already added that will
          // become redundant when this intersection is added to fEntries.
          MarkForDeletionSubsetsOf(intersect, startOfIntersects);
          
          new (fEntries[newIndex++]) AliHLTDomainEntry(kTRUE, intersect);
          if (currentEntry->IdenticalTo(intersect))
          {
            currentEntry->SetBit(14, true);
          }
        }
      }
    }
    else
    {
      // For an exclusive entry in 'domain' we need to find the intersections with
      // all of fEntries and re-apply these with the same exclude flags.
      for (Int_t j = 0; j < startOfIntersects; j++)
      {
        AliHLTDomainEntry* currentEntry = static_cast<AliHLTDomainEntry*>( fEntries[j] );
        if (checkEntry->IntersectWith(*currentEntry, intersect))
        {
          // We can remove all intersections that were already added that will
          // become redundant when this intersection is added to fEntries.
          MarkForDeletionSubsetsOf(intersect, startOfIntersects);
          
          new (fEntries[newIndex++]) AliHLTDomainEntry(currentEntry->Exclusive(), intersect);
        }
      }
    }
  }

  RemoveMarkedEntries();
  Optimise();
  return *this;
}


AliHLTTriggerDomain AliHLTTriggerDomain::operator ~ () const
{
  // Performs a set complement of the trigger domain.
  
  // The set complement is calculated by creating a new trigger domain which
  // accepts all possible data blocks, and then apply all the trigger domain
  // entries (rules / patterns) from top to bottom, but apply them with the
  // opposite meaning. For example, this->fEntries contains an inclusive domain
  // entry then remove it from the new trigger domain 'result', but if it is
  // an exclusion then add it.
  AliHLTTriggerDomain result;
  result.Add(kAliHLTAnyDataType);
  for (Int_t i = 0; i < fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* entry = static_cast<const AliHLTDomainEntry*>( fEntries[i] );
    if (entry->Inclusive())
    {
      result.Remove(*entry);
    }
    else
    {
      result.Add(*entry);
    }
  }
  return result;
}


AliHLTTriggerDomain AliHLTTriggerDomain::operator & (const AliHLTTriggerDomain& domain) const
{
  // This operator finds the set intersect.
  // See header file for more details.
  
  AliHLTTriggerDomain result;
  Int_t newIndex = 0;
  AliHLTDomainEntry intersect;
  
  // To find the set intersect we need to compare each entry in 'domain' to those
  // of fEntries. For each inclusive entry in 'domain' we need to add to the result
  // the intersect between it and each entry of fEntries, with the same exclude flag
  // value as the domain entry from fEntries.
  // However, in principle, for the exclusion entries in 'domain' we just add them
  // to the result, since those entries do not form part of the 'domain' trigger
  // domain set, so they should not form part of the result (remember any data block
  // must be contained in both trigger domains for a set intersect).
  // In actual fact we just add the intersect of the exclusion entries in 'domain'
  // with those of fEntries to the result. This has the same overall effect, but
  // makes sure that all exclusion entries are always subsets of inclusion entries.
  for (Int_t i = 0; i < domain.fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* checkEntry = static_cast<const AliHLTDomainEntry*>( domain.fEntries[i] );
    if (checkEntry->Inclusive())
    {
      for (Int_t j = 0; j < fEntries.GetEntriesFast(); j++)
      {
        AliHLTDomainEntry* currentEntry = static_cast<AliHLTDomainEntry*>( fEntries[j] );
        if (checkEntry->IntersectWith(*currentEntry, intersect))
        {
          // We can remove all entries that were already added to the result that
          // will become redundent because they are subsets of the new entry.
          result.MarkForDeletionSubsetsOf(intersect, 0);
          
          new (result.fEntries[newIndex++]) AliHLTDomainEntry(currentEntry->Exclusive(), intersect);
        }
      }
    }
    else
    {
      for (Int_t j = 0; j < fEntries.GetEntriesFast(); j++)
      {
        AliHLTDomainEntry* currentEntry = static_cast<AliHLTDomainEntry*>( fEntries[j] );
        if (checkEntry->IntersectWith(*currentEntry, intersect))
        {
          // We can remove all entries that were already added to the result that
          // will become redundant because they are subsets of the new entry.
          result.MarkForDeletionSubsetsOf(intersect, 0);
          
          new (result.fEntries[newIndex++]) AliHLTDomainEntry(kTRUE, intersect);
        }
      }
    }
  }

  result.RemoveMarkedEntries();
  result.Optimise();
  return result;
}


AliHLTTriggerDomain::operator AliHLTReadoutList () const
{
  // Typecast operator which constructs a readout list from the trigger domain.
  
  AliHLTReadoutList result;
  for (Int_t deti = 0; deti < AliDAQ::kNDetectors; deti++)
  {
    for (Int_t i = 0; i < AliDAQ::NumberOfDdls(deti); i++)
    {
      Int_t ddlId = AliDAQ::DdlID(deti, i);
      if (Contains(AliHLTDomainEntry("DAQRDOUT", AliDAQ::OnlineName(deti), ddlId)))
      {
        result.EnableDDLBit(ddlId);
      }
    }
  }
  return result;
}


void AliHLTTriggerDomain::MergeEntries(
    const bool* removeThisEntry, Int_t entriesCount,
    const bool* removeDomainEntry, Int_t startOfIntersects,
    const AliHLTTriggerDomain& domain
  )
{
  // Merges the entries in this trigger domain with the ones in 'domain', while
  // removing all entries that were marked for removal.
  // See header file for more information.
  
  bool anythingRemoved = false;
  
  // Remember this method is used at the end of the calculation of the binary operators
  // and that fEntries is expected to be partitioned into 3 regions.
  //   - 0..entriesCount-1 : contains the original (initial) entries of this trigger domain.
  //   - entriesCount..startOfIntersects-1 : is space reserved for the new entries
  //       from the given trigger domain 'domain' being processed.
  //   - startOfIntersects..fEntries.GetEntriesFast()-1 : contains all new domain entry
  //       intersection created and added to fEntries.
  //
  // First we need to remove all entries marked for removal from the original entries.
  for (Int_t i = 0; i < entriesCount; i++)
  {
    if (removeThisEntry[i])
    {
      fEntries.RemoveAt(i);
      anythingRemoved = true;
    }
  }
  
  // Now we copy over all the new entries from 'domain' which were not marked for removal
  // and indicate anythingRemoved = true since there will now be gaps in the clones array
  // that need to be compressed away later.
  for (Int_t i = 0; i < domain.fEntries.GetEntriesFast(); i++)
  {
    if (removeDomainEntry[i])
    {
      anythingRemoved = true;
    }
    else
    {
      const AliHLTDomainEntry* newEntry = static_cast<const AliHLTDomainEntry*>( domain.fEntries[i] );
      new (fEntries[entriesCount+i]) AliHLTDomainEntry(*newEntry);
    }
  }
  
  // Finally remove all new intersection entries that were marked for removal by
  // the MarkForDeletionSubsetsOf method.
  for (Int_t i = startOfIntersects; i < fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* ientry = static_cast<const AliHLTDomainEntry*>( fEntries[i] );
    if (ientry->TestBit(14))
    {
      fEntries.RemoveAt(i);
      anythingRemoved = true;
    }
  }
  if (anythingRemoved) fEntries.Compress();
}


void AliHLTTriggerDomain::MarkForDeletionSubsetsOf(const AliHLTDomainEntry& entry, Int_t min)
{
  // Marks for deletion all the entries in this trigger domain that are subsets
  // of the given entry.
  // See header file for more information.

  AliHLTDomainEntry intersect;
  for (Int_t i = min; i < fEntries.GetEntriesFast(); i++)
  {
    AliHLTDomainEntry* ientry = static_cast<AliHLTDomainEntry*>( fEntries[i] );
    if (ientry->TestBit(14)) continue;
    if (ientry->SubsetOf(entry))
    {
      ientry->SetBit(14, true);
    }
  }
}


void AliHLTTriggerDomain::RemoveMarkedEntries()
{
  // Removes all entries in this trigger domain which were marked for removal.
  // See header file for more information.
  
  bool anythingRemoved = false;
  for (Int_t i = 0; i < fEntries.GetEntriesFast(); i++)
  {
    const AliHLTDomainEntry* ientry = static_cast<const AliHLTDomainEntry*>( fEntries[i] );
    if (ientry->TestBit(14))
    {
      fEntries.RemoveAt(i);
      anythingRemoved = true;
    }
  }
  if (anythingRemoved) fEntries.Compress();
}


void AliHLTTriggerDomain::Optimise()
{
  // Removes redundant trigger domain entries from the trigger domain.
  // See header file for more information.

  AliHLTDomainEntry intersect;
  
  // Check that the first entry is not and exclusion which would be redundent.
  if (fEntries.GetEntriesFast() == 0) return;
  AliHLTDomainEntry* firstEntry = static_cast<AliHLTDomainEntry*>( fEntries[0] );
  if (firstEntry->Exclusive()) firstEntry->SetBit(14, true);
  
  for (Int_t i = 1; i < fEntries.GetEntriesFast(); i++)
  {
    AliHLTDomainEntry* ientry = static_cast<AliHLTDomainEntry*>( fEntries[i] );
    
    // For the i'th entry in fEntries, compare it in reverse order with all other
    // entries that are before it and look for redundant ones, i.e. that are subsets
    // of the i'th entry.
    for (Int_t j = i-1; j >= 0; j--)
    {
      AliHLTDomainEntry* jentry = static_cast<AliHLTDomainEntry*>( fEntries[j] );
      if (jentry->TestBit(14)) continue;
      // Find entries that intersect
      if (jentry->SubsetOf(*ientry))
      {
        // jentry is a subset of ientry so it is redundant because for all values
        // ientry will override jentry when calling IncludeInReadout.
        jentry->SetBit(14, true);
      }
      else if (*ientry == *jentry)
      {
        // If intersecting entries have opposite exclude flags then search no further,
        // we know that we will need this entry for correct behaviour of IncludeInReadout.
        if (ientry->Inclusive() == jentry->Exclusive()) goto processNextEntry;
        
        if (ientry->SubsetOf(*jentry))
        {
          ientry->SetBit(14, true);
          goto processNextEntry;
        }
      }
    }
    
    // If we got to this point then we hit the top of the trigger domain rules
    // (pattern matching) list without hitting any and overlapping entries.
    // So now we need to check if ientry is an exclusion. If it is, then it is
    // redundant and we can mark it for removal.
    if (ientry->Exclusive()) ientry->SetBit(14, true);
    
    processNextEntry: ;
  }
  
  RemoveMarkedEntries();
}

