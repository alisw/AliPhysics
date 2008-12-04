#ifndef ALIHLTTRIGGERDOMAIN_H
#define ALIHLTTRIGGERDOMAIN_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTDomainEntry.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Nov 2008
/// @brief  Declaration of the AliHLTTriggerDomain class used to store the set of data block types to readout.

#include "TObject.h"
#include "TClonesArray.h"
#include "AliHLTDataTypes.h"

class AliHLTDomainEntry;

/**
 * \class AliHLTTriggerDomain
 * The trigger domain class is the set of HLT raw data block types that should
 * be readout and sent to HLTOUT.
 * It is implemented as a list of domain entries, where each domain entry is
 * like a rule, or pattern to match against. When trying to decide if a given
 * data block falls within the trigger domain, i.e. is part of the readout, the
 * domain entries are applied one after the other from top to bottom. Two kinds
 * of domain entries, inclusive and exclusive are possible, which indicate if a
 * data block type or range, is part of the trigger domain or not. As we process
 * the domain entries we update our decision of whether the data block is part
 * of the trigger domain or not. If the domain entry is an inclusion then we update
 * the decision to true, if it is an exclusion then we update the decision to false.
 * The value of the result after applying the last domain entry then indicates
 * if the data block is part of the trigger domain or not.
 * In this way we can specify trigger domains as sets to arbitrary complexity
 * and manipulate them as mathematical sets accordingly.
 *
 * The other feature of the AliHLTTriggerDomain class is that it overloads the
 * following operators to provide set like behaviour:
 *   |  this provides the set union operation (The + operator does the same).
 *   &  this provides the set intersect operation.
 *   -  this provides the set difference operation.
 *   ^  this provides an exclusive or (xor) operation. i.e. given two sets A and B
 *        the result C is given by:
 *        C = {x : x elf A and not x elf B, or x elf B and not x elf A}
 *        where 'elf' means "is an element of".
 *   ~  this returns the set complement.
 * These operators then allow expressions to be formed from trigger domain objects
 * which behave like sets.
 */
class AliHLTTriggerDomain : public TObject
{
 public:
  
  /**
   * Default constructor.
   */
  AliHLTTriggerDomain();
  
  /**
   * Copy constructor performs a deep copy.
   * \param domain  The domain entry to copy from.
   */
  AliHLTTriggerDomain(const AliHLTTriggerDomain& domain);
  
  /**
   * Default destructor.
   */
  virtual ~AliHLTTriggerDomain();
  
  /**
   * Adds the given entry to this trigger domain as an inclusive entry.
   * Existing entries are modified as required to optimise the trigger domain
   * rule / pattern matching list.
   * \param entry  The domain entry object to add.
   * \note The entry.Exclusive() flag is ignored and is treated as if it was kFALSE.
   */
  void Add(const AliHLTDomainEntry& entry);
  
  /**
   * Adds the given data type to the trigger domain such that all data blocks
   * that match this type will form part of the trigger domain.
   * \param datatype  The data block type and origin to match.
   */
  void Add(const AliHLTComponentDataType& datatype);
  
  /**
   * Adds the given data type and origin to the trigger domain such that all data
   * blocks that match will form part of this trigger domain.
   * \param blocktype  The data block type string of the data block that must match.
   *    The value of kAliHLTAnyDataTypeID can be used to specify the 'any' type
   *    wild card value.
   * \param origin  The origin of the data block, such as the detector name, that
   *    must match. The value of kAliHLTDataOriginAny can be used to specify the
   *    'any' origin wild card value.
   */
  void Add(const char* blocktype, const char* origin);
  
  /**
   * Adds the given data type with particular specification bits to the trigger
   * domain, such that all data blocks that match these will form part of this
   * trigger domain.
   * \param datatype  The data block type and origin that must match.
   * \param spec  The specification bits that must match.
   */
  void Add(const AliHLTComponentDataType& datatype, UInt_t spec);
  
  /**
   * Adds the given data type, origin and specification bits of data blocks that
   * should form part of this trigger domain.
   * \param blocktype  The data block type string of the data block that must match.
   *    The value of kAliHLTAnyDataTypeID can be used to specify the 'any' type
   *    wild card value.
   * \param origin  The origin of the data block, such as the detector name, that
   *    must match. The value of kAliHLTDataOriginAny can be used to specify the
   *    'any' origin wild card value.
   * \param spec  The specification bits that must match.
   */
  void Add(const char* blocktype, const char* origin, UInt_t spec);
  
  /**
   * Removes or modifies all entries from the trigger domain, such that data blocks
   * that match the given domain entry will not form part of this trigger domain.
   * Existing entries are modified as required to optimise the trigger domain
   * rule / pattern matching list.
   * \param entry  The domain entry object to indicating values that should be removed.
   * \note The entry.Exclusive() flag is ignored and is treated as if it was kTRUE.
   */
  void Remove(const AliHLTDomainEntry& entry);
  
  /**
   * Removes the given data type from the trigger domain, such that all data blocks
   * that match this type will not form part of the trigger domain.
   * \param datatype  The data block type and origin that must match the blocks not
   *    forming part of this trigger domain.
   */
  void Remove(const AliHLTComponentDataType& datatype);
  
  /**
   * Removes the given data type and origin from the trigger domain, such that all
   * data blocks that match these will not form part of the trigger domain.
   * \param blocktype  The data block type string that must match the data blocks
   *    not forming part of this trigger domain.
   *    The value of kAliHLTAnyDataTypeID can be used to specify the 'any' type
   *    wild card value.
   * \param origin  The origin string, such as the detector name, that must match
   *    the data blocks not forming part of this trigger domain.
   *    The value of kAliHLTDataOriginAny can be used to specify the 'any' origin
   *    wild card value.
   */
  void Remove(const char* blocktype, const char* origin);
  
  /**
   * Removes the given data type with given specification bit from the trigger
   * domain, such that all data blocks that match these will not form part of the
   * trigger domain.
   * \param datatype  The data block type and origin that must match the blocks
   *    not forming part of this trigger domain.
   * \param spec  The specification bits that must match for the blocks that do
   *    not form part of this trigger domain.
   */
  void Remove(const AliHLTComponentDataType& datatype, UInt_t spec);
  
  /**
   * Removes the given data type, origin and specification from the trigger domain,
   * such that all data blocks that match these will not form part of the trigger
   * domain.
   * \param blocktype  The data block type string that must match the data blocks
   *    not forming part of this trigger domain.
   *    The value of kAliHLTAnyDataTypeID can be used to specify the 'any' type
   *    wild card value.
   * \param origin  The origin string, such as the detector name, that must match
   *    the data blocks not forming part of this trigger domain.
   *    The value of kAliHLTDataOriginAny can be used to specify the 'any' origin
   *    wild card value.
   * \param spec  The specification bits that must match for the blocks that do
   *    not form part of this trigger domain.
   */
  void Remove(const char* blocktype, const char* origin, UInt_t spec);
  
  /**
   * This checks to see if the given entry (or class of entries, if the entry uses
   * wild card values) is part of the trigger domain.
   * \param entry  This is the entry to check for.
   * \return  true if data blocks that match the entry are part of this
   *    trigger domain and false otherwise.
   * \note If the block contains the 'any' wild card values for the data type
   *    origin or specification, then the inclusive domains are treated
   *    optimistically and the exclusive domains pessimistically. This means that
   *    the wild card values are assumed to fall within the trigger domain for the
   *    optimistic case, but fall outside the domain for the pessimistic case.
   */
  bool Contains(const AliHLTDomainEntry& entry) const;
  
  /**
   * This checks to see if the given data block should be included in the HLT readout.
   * \param block  The data block descriptor to check.
   * \return  true if data block forms part of this trigger domain and should
   *    be part of the readout and false otherwise.
   * \note If the block contains the 'any' wild card values for the data type
   *    or origin, then the inclusive domains are treated optimistically and the
   *    exclusive domains pessimistically. This means that the wild card values
   *    are assumed to fall within the trigger domain for the optimistic case,
   *    but fall outside the domain for the pessimistic case.
   */
  bool IncludeInReadout(const AliHLTComponentBlockData* block) const;
  
  /**
   * This checks to see if the given data block should not be included in the
   * HLT readout.
   * \param block  The data block descriptor to check.
   * \return  true if data block does not form part of this trigger domain and
   *    false otherwise.
   */
  bool ExcludeFromReadout(const AliHLTComponentBlockData* block) const
  {
    return ! IncludeInReadout(block);
  }
  
  /**
   * This method removes all entries in the trigger domain list, giving us and
   * empty trigger domain set.
   * \param  option  This is passed onto the internal fEntries TClonesArray.
   * The method is inherited from TObject.
   */
  virtual void Clear(Option_t* option = "");
  
  /**
   * Prints all the domain entries in this trigger domain in the order in which
   * they are applied and if they are inclusive or exclusive.
   * \param  option  This is not used by this method.
   * The method is inherited from TObject.
   */
  virtual void Print(Option_t* option = "") const;

  /**
   * Assignment operator performs a deep copy.
   * \param domain  The domain entry to copy from.
   * \return  A reference to this object.
   */
  AliHLTTriggerDomain& operator = (const AliHLTTriggerDomain& domain);
  
  /**
   * This operator adds all domain entries in <i>domain</i> to this trigger domain
   * in such a way, so as to effectively perform a set union.
   * \param domain  The domain object on the right hand side of the operator.
   * \return  A reference to this object.
   */
  AliHLTTriggerDomain& operator |= (const AliHLTTriggerDomain& domain);
  
  /**
   * This operator adds all domain entries in <i>domain</i> that do not exist in
   * this trigger domain, but removes all entries that do exist, effectively
   * performing an exclusive or (xor) operation.
   * i.e. given two sets A and B the result C is given by:
   *    C = {x : x elf A and not x elf B, or x elf B and not x elf A}
   * where 'elf' means "is an element of".
   * \param domain  The domain object on the right hand side of the operator.
   * \return  A reference to this object.
   */
  AliHLTTriggerDomain& operator ^= (const AliHLTTriggerDomain& domain);
  
  /**
   * This operator removes all domain entries from this trigger domain that do
   * not also exist in <i>domain</i>, effectively performing a set intersect.
   * \param domain  The domain object on the right hand side of the operator.
   * \return  A reference to this object.
   */
  AliHLTTriggerDomain& operator &= (const AliHLTTriggerDomain& domain)
  {
    return this->operator = (*this & domain);
  }
  
  /**
   * This operator performs the same operation as the '|=' operator.
   * \param domain  The domain object on the right hand side of the operator.
   * \return  A reference to this object.
   */
  AliHLTTriggerDomain& operator += (const AliHLTTriggerDomain& domain)
  {
    return operator |= (domain);
  }
  
  /**
   * This operator removes all domain entries from this trigger domain that exisit
   * in <i>domain</i>, effectively implementing a set difference.
   * \param domain  The domain object on the right hand side of the operator.
   * \return  A reference to this object.
   */
  AliHLTTriggerDomain& operator -= (const AliHLTTriggerDomain& domain);
  
  /**
   * This operator returns the set complement of the trigger domain.
   * \return  The complement of this trigger domain, such that any data block that
   *    returns true for AliHLTTriggerDomain::IncludeInReadout() for this trigger
   *    domain, will return false for the same method call in the returned object.
   */
  AliHLTTriggerDomain operator ~ () const;
  
  /**
   * This operator performs a set union between this trigger domain and <i>domain</i>.
   * \param domain  The domain object on the right hand side of the operator.
   * \return  a trigger domain object R, such that for each data block D, we will have
   *    R.IncludeInReadout(D) == this->IncludeInReadout(D) or domain.IncludeInReadout(D)
   */
  AliHLTTriggerDomain operator | (const AliHLTTriggerDomain& domain) const
  {
    AliHLTTriggerDomain result = *this;
    return result.operator |= (domain);
  }
  
  /**
   * This operator performs an exclusive or (xor) like operation between this trigger
   * domain and <i>domain</i>.
   * \param domain  The domain object on the right hand side of the operator.
   * \return  a trigger domain object R, such that for each data block D, we will have
   *    R.IncludeInReadout(D) == this->IncludeInReadout(D) xor domain.IncludeInReadout(D)
   */
  AliHLTTriggerDomain operator ^ (const AliHLTTriggerDomain& domain) const
  {
    AliHLTTriggerDomain result = *this;
    return result.operator ^= (domain);
  }
  
  /**
   * This operator performs a set intersect operation between this trigger domain
   * and <i>domain</i>.
   * \param domain  The domain object on the right hand side of the operator.
   * \return  a trigger domain object R, such that for each data block D, we will have
   *    R.IncludeInReadout(D) == this->IncludeInReadout(D) and domain.IncludeInReadout(D)
   */
  AliHLTTriggerDomain operator & (const AliHLTTriggerDomain& domain) const;
  
  /**
   * This operator performs the same operation as the '|' operator.
   * \param domain  The domain object on the right hand side of the operator.
   * \return  a trigger domain object R, such that for each data block D, we will have
   *    R.IncludeInReadout(D) == this->IncludeInReadout(D) or domain.IncludeInReadout(D)
   */
  AliHLTTriggerDomain operator + (const AliHLTTriggerDomain& domain) const
  {
    AliHLTTriggerDomain result = *this;
    return result.operator += (domain);
  }
  
  /**
   * This operator implements the set difference between this trigger domain and
   * <i>domain</i>.
   * \param domain  The domain object on the right hand side of the operator.
   * \return  a trigger domain object R, such that for each data block D, we will have
   *    R.IncludeInReadout(D) == this->IncludeInReadout(D) and not domain.IncludeInReadout(D)
   */
  AliHLTTriggerDomain operator - (const AliHLTTriggerDomain& domain) const
  {
    AliHLTTriggerDomain result = *this;
    return result.operator -= (domain);
  }
  
 private:
  
  /**
   * This method merges the domain entries from <i>domain</i> by copying them into
   * fEntries, but only the ones not marked for removal in <i>removeDomainEntry</i>.
   * Any entries that were marked for removal in fEntries are also removed.
   * \param removeThisEntry  Flags which indicate if the corresponding fEntries[i]
   *    should be removed.
   * \param entriesCount The number of entries in <i>removeThisEntry</i>.
   * \param removeDomainEntry  Flags which indicate if the corresponding domain.fEntries[i]
   *    was marked for removal or not. If marked for removal then it will not be copied
   *    into this trigger domain. The size of the array is given by domain.GetEntriesFast().
   * \param startOfIntersects  This is the start location of the new intersection domain
   *    entries that were added to fEntries. i.e. fEntries[startOfIntersects] is the
   *    first new intersect entry.
   */
  void MergeEntries(
      const bool* removeThisEntry, Int_t entriesCount,
      const bool* removeDomainEntry, Int_t startOfIntersects,
      const AliHLTTriggerDomain& domain
    );
  
  /**
   * Goes throught the list of domain entries in fEntries from the first entry
   * indicated by 'min' to the end of the list and marks for deletion all entries
   * in fEntries that are subsets of 'entry'.
   * The entries are marked by setting the 14'th bit in fBits with a call to
   * AliHLTDomainEntry::SetBit(14, true).
   * \param entry  The entry that should be the super set of the entries we mark
   *    for removal.
   * \param min  This is the first entry we consider, all the way up to
   *    fEntries.GetEntriesFast() - 1.
   */
  void MarkForDeletionSubsetsOf(const AliHLTDomainEntry& entry, Int_t min);
  
  /**
   * Removes all entries in this trigger domain which were marked for removal.
   * These are all domain entries that have the 14'th bit set in their fBits field
   * with a call to AliHLTDomainEntry::SetBit(14, true).
   */
  void RemoveMarkedEntries();
  
  /**
   * Removes any redundant trigger domain entries from the fEntries list.
   * Entries that are subsets of each other are removed. Also exclusive entries
   * that are not subsets of any inclusive entry are also removed, because we
   * implicitly assume a data block does not form part of the trigger domain,
   * unless explicitly included with an inclusive domain entry. So these kinds
   * of entries are redundant.
   */
  void Optimise();
  
  TClonesArray fEntries;  /// The list of domain entries used to decide if a data block forms part of trigger domain set.
  
  ClassDef(AliHLTTriggerDomain, 1) // This is a list of internal HLT data block types which should be forwarded for readout.

};

#endif // ALIHLTTRIGGERDOMAIN_H

