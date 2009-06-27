#ifndef ALIHLTDOMAINENTRY_H
#define ALIHLTDOMAINENTRY_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTDomainEntry.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   20 Nov 2008
/// @brief  Declaration of the AliHLTDomainEntry class used to store identifying information about HLT data blocks.

#include "TObject.h"
#include "AliHLTDataTypes.h"

class TString;

/**
 * \class AliHLTDomainEntry
 * The AliHLTDomainEntry class is used to store information identifying a particular
 * HLT internal data block, or set of data blocks using wild card values. This
 * class is used by AliHLTTriggerDomain to store a list of data block classes
 * that should be readout by the HLT. The information identifying a data block is
 * the following:
 *  - the data block type
 *  - the data block's origin (detector name)
 *  - the data block's specification (detector specific bits)
 * Several useful operators and methods are defined to help manipulate this
 * information in the AliHLTTriggerDomain class.
 */
class AliHLTDomainEntry : public TObject
{
 public:
  
  /**
   * Default constructor.
   */
  AliHLTDomainEntry();
  
  /**
   * Copy constructor performs a deep copy.
   * \param domain  The domain entry to copy from.
   */
  AliHLTDomainEntry(const AliHLTDomainEntry& domain);
  
  /**
   * This constructs a domain entry with a particular data type and an 'any' wild
   * card value for specification indicating any specification will match.
   * \param type  The data block type and origin to use.
   */
  AliHLTDomainEntry(const AliHLTComponentDataType& type);
  
  /**
   * This constructs a domain entry with a particular data type and origin. The
   * specification is marked as an 'any' wild card value, indicating any data
   * block specification will match.
   * \param blocktype  The data block type string of the data block. The value
   *    kAliHLTAnyDataTypeID can be used to specify the any type wild card value.
   * \param origin  The origin of the data block, such as the detector name.
   *    The value kAliHLTDataOriginAny can be used to specify the any origin
   *    wild card value.
   */
  AliHLTDomainEntry(const char* blocktype, const char* origin);
  
  /**
   * This constructs a domain entry with a particular data type, origin and
   * specification.
   * \param type  The data block type and origin to use.
   * \param spec  The data block specification to use.
   */
  AliHLTDomainEntry(const AliHLTComponentDataType& type, UInt_t spec);
  
  /**
   * This constructs a domain entry with a particular data type, origin and
   * specification.
   * \param blocktype  The data block type string of the data block. The value
   *    kAliHLTAnyDataTypeID can be used to specify the any type wild card value.
   * \param origin  The origin of the data block, such as the detector name.
   *    The value kAliHLTDataOriginAny can be used to specify the any origin
   *    wild card value.
   * \param spec  The data block specification to use.
   */
  AliHLTDomainEntry(const char* blocktype, const char* origin, UInt_t spec);
  
  /**
   * The constructor deep copies the domain entry but overrides the exclude flag.
   * \param exclude  The new exclude flag value to use. If 'true' then the entry
   *    forms part of a trigger domain exclusion rule.
   * \param domain  The domain entry to copy from.
   */
  AliHLTDomainEntry(Bool_t exclude, const AliHLTDomainEntry& domain);
  
  /**
   * This constructs a domain entry with a particular data type and exclude flag
   * value, but an 'any' wild card value is used for the data block specification.
   * \param exclude  The exclude flag value to use. If 'true' then the entry forms
   *    part of a trigger domain exclusion rule.
   * \param type  The data block type and origin to use.
   */
  AliHLTDomainEntry(Bool_t exclude, const AliHLTComponentDataType& type);
  
  /**
   * This constructs a domain entry with a particular data type, origin and exclusion
   * value. The specification is marked as an 'any' wild card value, indicating any
   * data block specification will match.
   * \param exclude  The exclude flag value to use. If 'true' then the entry forms
   *    part of a trigger domain exclusion rule.
   * \param blocktype  The data block type string of the data block. The value
   *    kAliHLTAnyDataTypeID can be used to specify the any type wild card value.
   * \param origin  The origin of the data block, such as the detector name.
   *    The value kAliHLTDataOriginAny can be used to specify the any origin
   *    wild card value.
   */
  AliHLTDomainEntry(Bool_t exclude, const char* blocktype, const char* origin);
  
  /**
   * This constructs a domain entry with a particular exclude flag value, data type,
   * origin and specification.
   * \param exclude  The exclude flag value to use. If 'true' then the entry forms
   *    part of a trigger domain exclusion rule.
   * \param type  The data block type and origin to use.
   * \param spec  The data block specification to use.
   */
  AliHLTDomainEntry(Bool_t exclude, const AliHLTComponentDataType& type, UInt_t spec);
  
  /**
   * This constructs a domain entry with a particular exclude flag value, data type,
   * origin and specification.
   * \param exclude  The exclude flag value to use. If 'true' then the entry forms
   *    part of a trigger domain exclusion rule.
   * \param blocktype  The data block type string of the data block. The value
   *    kAliHLTAnyDataTypeID can be used to specify the any type wild card value.
   * \param origin  The origin of the data block, such as the detector name.
   *    The value kAliHLTDataOriginAny can be used to specify the any origin
   *    wild card value.
   * \param spec  The data block specification to use.
   */
  AliHLTDomainEntry(Bool_t exclude, const char* blocktype, const char* origin, UInt_t spec);
  
  /**
   * Default destructor.
   */
  virtual ~AliHLTDomainEntry();
  
  /**
   * Returns the value of the exclude flag.
   * \return  true if the domain entry is an exclusion and the matching data blocks
   *    should not form part of the trigger domain for readout.
   */
  Bool_t Exclusive() const { return fExclude; }
  
  /**
   * Sets the value of the exclude flag.
   * \param value  The value to set the flag to. If 'true' then the domain entry
   *    is an exclusion and the matching data blocks should not form part of the
   *    trigger domain for readout. If 'false' then the matching data blocks should
   *    form part of the readout.
   */
  void Exclusive(Bool_t value) { fExclude = value; }
  
  /**
   * Indicates if the domain entry is an inclusive rule.
   * \return  true if the domain entry is an inclusion and the matching data blocks
   *    should form part of the trigger domain for readout.
   */
  Bool_t Inclusive() const { return ! fExclude; }
  
  /**
   * Used to set if the domain entry should be an inclusion or exclusion.
   * \param value  The value to set. If 'true' then the domain entry is an inclusion
   *    and the matching data blocks should form part of the trigger domain for readout.
   *    If 'false' then the matching data blocks should not form part of the readout.
   */
  void Inclusive(Bool_t value) { fExclude = ! value; }
  
  /**
   * Returns the data type of the domain entry.
   * \return  The data type that data blocks are compared to.
   */
  const AliHLTComponentDataType& DataType() const { return fType; }
  
  /**
   * Indicates if the specification is used.
   * \return  true if the specification is used when matching data blocks, otherwise
   *     false, indicating that the specification is treated as a wild card value.
   */
  Bool_t IsValidSpecification() const { return fUseSpec; }
  
  /**
   * Returns the data block specification of the domain entry.
   * \return  The data block specification that data blocks are compared to.
   */
  UInt_t Specification() const { return fSpecification; }
  
  /**
   * The copy operator performs a deep copy.
   * \param domain  The domain entry to copy from.
   */
  AliHLTDomainEntry& operator = (const AliHLTDomainEntry& domain);
  
  /**
   * The comparison operator checks to see if two domain entries match.
   * \param rhs  The right hand side domain entry to compare to.
   * \return  true if the domain entries are identical or if they overlap (match)
   *    due to wild card values. False is returned if there is absolutely no
   *    overlap between this and the right hand side domain entries.
   */
  bool operator == (const AliHLTDomainEntry& rhs) const
  {
    return (fType == rhs.fType) && (fUseSpec && rhs.fUseSpec ? fSpecification == rhs.fSpecification : true);
  }
  
  /**
   * The comparison operator checks to see if two domain entries do not match.
   * \param rhs  The right hand side domain entry to compare to.
   * \return  true if the domain entries do not overlap (match) in any way, also
   *    after considering any wild card values. False is returned if the entries
   *    are identical or if they overlap due to wild card values.
   */
  bool operator != (const AliHLTDomainEntry& rhs) const
  {
    return ! this->operator == (rhs);
  }
  
  /**
   * The comparison operator checks to see if the data block matches the domain entry.
   * \note The data block's specification is treated as exact and never as a wild card
   *    'any' value. To be able to treat the specification as 'any', create a new
   *    AliHLTDomainEntry object with the
   *      \code AliHLTDomainEntry(const AliHLTComponentDataType& type) \endcode
   *    constructor, using the data blocks type for the <i>type</i> parameter.
   *    With the new AliHLTDomainEntry object one can make the required wild card comparison.
   * \param block  The data block to compare to.
   * \return  true if the data block matches the domain entry and false otherwise.
   */
  bool operator == (const AliHLTComponentBlockData* block) const
  {
    return (fType == block->fDataType) && (fUseSpec ? fSpecification == block->fSpecification : true);
  }
  
  /**
   * The comparison operator checks to see if the data block does not match the domain entry.
   * \note The data block's specification is treated as exact and never as a wild card
   *    'any' value. To be able to make the required comparison treating the specification
   *    as 'any' try the following code:
   *    \code
   *      AliHLTComponentBlockData* block;  // assumed initialised.
   *      AliHLTDomainEntry entryToCompareTo;  // assumed initialised.
   *      AliHLTDomainEntry newEntryForBlock(block->fDataType);
   *      bool comparisonResult = (entryToCompareTo == newEntryForBlock);
   *    \endcode
   * \param block  The data block to compare to.
   * \return  true if the data block matches the domain entry and false otherwise.
   */
  bool operator != (const AliHLTComponentBlockData* block) const
  {
    return ! this->operator == (block);
  }
  
  /**
   * This typecast operator returns the data type of the domain entry.
   * \return  Copy of the data block type structure.
   */
  operator AliHLTComponentDataType () const { return fType; }
  
  /**
   * Compares this domain entry to another to see if they are identical.
   * \param rhs  The domain entry to compare to.
   * \return  True if the two domain entries have the same data types, origins and
   *   specifications, character for character, ignoring wild card symantics.
   *   False is returned otherwise.
   * \note No comparison is done for the exclude flag.
   */
  bool IdenticalTo(const AliHLTDomainEntry& rhs) const;
  
  /**
   * Compares this domain entry is a subset of the given entry.
   * If we consider the possibility of wild card characters, then the domain entry
   * can be thought of as a set of possible data block entries. This operator
   * therefore effectively implements set logic.
   * \param rhs  The domain entry to compare to.
   * \return  True if the this domain entry is either identical to <i>rhs</i>, i.e.
   *   IdenticalTo(rhs) returns true, or if <i>rhs</i> can match to all data blocks
   *   that this domain entry can match to, but also <i>rhs</i> can match to other
   *   data blocks that this entry cannot match to.
   * \note No comparison is done for the exclude flag.
   */
  bool SubsetOf(const AliHLTDomainEntry& rhs) const;
  
  /**
   * Finds the set intersection between this domain entry and 'rhs', and puts the
   * intersection value into 'result'.
   * If we consider the possibility of wild card characters, then the domain entry
   * can be thought of as a set of possible data block entries. This operator
   * therefore effectively implements the set intersection.
   * \param rhs <i>[in]</i> The domain entry to compare to.
   * \param result <i>[out]</i>  The resulting intersect is written into this
   *    variable if this method returns true. The contents is not modified if
   *    there is no intersect and this method returns false.
   * \return true is returned if there is a intersect between the domain entries
   *    and false otherwise.
   */
  bool IntersectWith(const AliHLTDomainEntry& rhs, AliHLTDomainEntry& result) const;
  
  /**
   * Inherited from TObject. Prints the domain entry in the following format:<br>
   *  \<type\>:\<origin\>:\<specification\><br>
   * where<br>
   *  \<type\> is the 8 character data block type.<br>
   *  \<origin\> is the 4 character data block origin.<br>
   *  \<specification\> is the data block specification printed in hexadecimal format.<br>
   * The "\0" string is printed for NULL characters in the type and origin strings.
   * While "********" is printed for the 'any' data type wild card value, "****"
   * is printed for the 'any' origin wild card value and "**********" is printed
   * for the 'any' specification wild card.
   * \param option  If set to "noendl" then no end of line is printed.
   */
  virtual void Print(Option_t* option = "") const;
  
  /**
   * Converts the domain entry type, origin and specification into a string
   * representation.
   * \returns  A string in the format \<type\>:\<origin\>:\<specification\>
   */
  TString AsString() const;
  
 private:
  
  Bool_t fExclude;  /// Indicates if the domain entry is exclusive, indicating data blocks that should not be readout.
  Bool_t fUseSpec;  /// Indicates if the fSpecification field should be used. If not set then the specification is treated as an 'any' wild card value.
  AliHLTComponentDataType fType;  /// The HLT data block type.
  UInt_t fSpecification;  /// The data block specification to match.
  
  ClassDef(AliHLTDomainEntry, 1) // A data block type and possible specification entry, which forms part of a trigger domain.
};

#endif // ALIHLTDOMAINENTRY_H

