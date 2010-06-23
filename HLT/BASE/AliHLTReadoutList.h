//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTREADOUTLIST_H
#define ALIHLTREADOUTLIST_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTReadoutList.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Nov 2008
/// @brief  Declaration of the AliHLTReadoutList class used to handle AliHLTEventDDL structures.

#include "TNamed.h"
#include "AliHLTDataTypes.h"

/**
 * \class AliHLTReadoutList
 * This class is used as an interface or wrapper to the AliHLTEventDDL structure.
 * It makes it easy to manipulate the bits in this structure, which define what DDLs
 * should be readout by DAQ.
 * Several operators are also overloaded which are meant to be used in the trigger
 * menu specification for the AliHLTGlobalTrigger. It allows one to construct
 * expressions for the readout lists, which is necessary to be able to evaluate
 * or compose the final readout list, given multiple input readout lists received
 * from individual components that derive from AliHLTTrigger.
 * The operators implemented are:
 *  |  applies a bitwise or on the DDL bits.
 *  &  applies a bitwise and on the DDL bits.
 *  ^  applies a bitwise xor on the DDL bits.
 *  ~  applies a bitwise not on the DDL bits.
 *  -  unsets the bits in readout list A that are set in readout list B.
 *      This effectively applies A & (A ^ B).
 */
class AliHLTReadoutList : public TNamed
{
 public:
  
  /**
   * Identifiers for different detectors used by methods in AliHLTReadoutList.
   */
  enum EDetectorId
  {
    kITSSPD = 0x1 << 0,    /// ID for SPD detector
    kITSSDD = 0x1 << 1,    /// ID for SDD detector
    kITSSSD = 0x1 << 2,    /// ID for SSD detector
    kTPC = 0x1 << 3,       /// ID for TPC detector
    kTRD = 0x1 << 4,       /// ID for TRD detector
    kTOF = 0x1 << 5,       /// ID for TOF detector
    kHMPID = 0x1 << 6,     /// ID for HMPID detector
    kPHOS = 0x1 << 7,      /// ID for PHOS detector
    kCPV = 0x1 << 8,       /// ID for CPV detector
    kPMD = 0x1 << 9,       /// ID for PMD detector
    kMUONTRK = 0x1 << 10,  /// ID for MUON tracking chambers
    kMUONTRG = 0x1 << 11,  /// ID for MUON trigger detector
    kFMD = 0x1 << 12,      /// ID for FMD detector
    kT0 = 0x1 << 13,       /// ID for T0 detector
    kV0 = 0x1 << 14,       /// ID for V0 detector
    kZDC = 0x1 << 15,      /// ID for ZDC detector
    kACORDE = 0x1 << 16,   /// ID for ACORDE detector
    kTRG = 0x1 << 17,      /// ID for TRG detector
    kEMCAL = 0x1 << 18,    /// ID for EMCAL detector
    kDAQTEST = 0x1 << 19,  /// ID for DAQ_TEST detector
    kHLT = 0x1 << 30,      /// ID for HLT detector
    // kALLDET sets readout for all detectors except DAQ_TEST
    kALLDET = (kITSSPD | kITSSDD | kITSSSD | kTPC | kTRD | kTOF | kHMPID | kPHOS
               | kCPV | kPMD | kMUONTRK | kMUONTRG | kFMD | kT0 | kV0 | kZDC
               | kACORDE | kTRG | kEMCAL | kHLT)
  };
  
  /**
   * Default constructor.
   */
  AliHLTReadoutList();
  
  /**
   *  Constructor to select which detectors to enable for readout.
   * \param enabledDetectors  Detector bit field. Can be any values for
   *     EDetectorId or'ed together.
   */
  AliHLTReadoutList(Int_t enabledDetectors);
  
  /**
   * Constructor to select which detectors and DDLs to enable for readout.
   * \param enabledList The string format is a space separated list where
   *     each item is either a detector acronym name or DDL number.
   * Invalid sub-strings are simply ignored. The special ALL string is
   * equivalent to kALLDET for AliHLTReadoutList(Int_t enabledDetectors).
   */
  AliHLTReadoutList(const char* enabledList);
  
  /**
   * Constructor to create readout list from AliHLTEventDDL structure.
   * \param list  The AliHLTEventDDL structure from which to create this object.
   */
  AliHLTReadoutList(const AliHLTEventDDL& list);
  
  /**
   * The copy constructor performs a deep copy.
   * \param list  The readout list to copy from.
   */
  AliHLTReadoutList(const AliHLTReadoutList& list);
  
  /**
   * Default destructor.
   */
  virtual ~AliHLTReadoutList();
  
  /**
   * Checks if the readout list is empty, i.e. all DDLs are disabled.
   * \returns true if the readout list is empty and false otherwise.
   */
  bool Empty() const;
  
  /**
   * Disables all bits in the readout list.
   * \param  option  This parameter is ignored.
   * The method is inherited from TObject.
   */
  virtual void Clear(Option_t* option = "");
  
  /**
   * Enables a specific DDL bit in the readout list.
   * \param ddlId  The ID number of the DDL to enable.
   */
  void EnableDDLBit(Int_t ddlId)
  {
    SetDDLBit(ddlId, kTRUE);
  }
  
  /**
   * Disables a specific DDL bit in the readout list.
   * \param ddlId  The ID number of the DDL to disable.
   */
  void DisableDDLBit(Int_t ddlId)
  {
    SetDDLBit(ddlId, kFALSE);
  }
  
  /**
   * Fetches the bit value for a particular DDL in the readout list.
   * \param ddlId  The ID number of the DDL to fetch.
   * \return the bit value for the specified DDL.
   */
  Bool_t GetDDLBit(Int_t ddlId) const;
  
  /**
   * Sets the bit value for a particular DDL in the readout list.
   * \param ddlId  The ID number of the DDL to set.
   * \param state  The value to set the bit to.
   */
  void SetDDLBit(Int_t ddlId, Bool_t state);
  
  /**
   * Checks if a particular DDL is enabled for readout.
   * \param ddlId  The ID number of the DDL to check.
   * \return the if the DDL is enabled for readout.
   */
  bool IsDDLEnabled(Int_t ddlId) const
  {
    return GetDDLBit(ddlId) == kTRUE;
  }
  
  /**
   * Checks if a particular DDL is disabled for readout.
   * \param ddlId  The ID number of the DDL to check.
   * \return the if the DDL is disabled for readout.
   */
  bool IsDDLDisabled(Int_t ddlId) const
  {
    return GetDDLBit(ddlId) == kFALSE;
  }
  
  /**
   * Enables all DDLs for a particular detector or detectors.
   * \param detector  A bitmap of detectors to enable. Should be any values from
   *    EDetectorId that can be or'ed together for multiple detector selection.
   */
  void Enable(Int_t detector);
  
  /**
   * Disables all DDLs for a particular detector or detectors.
   * \param detector  A bitmap of detectors to disable. Should be any values from
   *    EDetectorId that can be or'ed together for multiple detector selection.
   */
  void Disable(Int_t detector);
  
  /**
   * Checks if a particular detector's DDLs are enabled for readout.
   * \param detector  A bitmap of detectors to check. Should be any values from
   *    EDetectorId that can be or'ed together for multiple detector selection.
   * \return true if all DDLs for the specified detectors are enabled for readout.
   */
  bool DetectorEnabled(Int_t ddlId) const;
  
  /**
   * Inherited from TObject. Prints the DDLs that will be readout according to
   * this readout list.
   * \param option  This is not used by this method.
   */
  virtual void Print(Option_t* option = "") const;
  
  /**
   * This typecast operator converts the readout list to the AliHLTEventDDL
   * structure format.
   * \return  Copy of the AliHLTEventDDL raw structure.
   */
  operator AliHLTEventDDL () const { return fReadoutList; }
  
  /**
   * This typecast operator converts the readout list to the AliHLTEventDDL
   * structure format.
   * \return  Reference to the AliHLTEventDDL raw structure.
   */
  operator AliHLTEventDDL& () { return fReadoutList; }

  /**
   * Access method to the binary buffer.
   * \return pointer to the binary buffer.
   */
  AliHLTEventDDL* Buffer() { return &fReadoutList; }

  /**
   * Access method to the binary buffer.
   * \return const pointer to the binary buffer.
   */
  const AliHLTEventDDL* Buffer() const { return &fReadoutList; }

  /**
   * Access to the size of the binary buffer.
   * \return size of the binary buffer
   */
  unsigned BufferSize() const { return sizeof(fReadoutList); }
  
  /**
   * Assignment operator performs a deep copy.
   * \param list  The readout list to copy from.
   * \return  A reference to this object.
   */
  AliHLTReadoutList& operator = (const AliHLTReadoutList& list);
  
  /**
   * This operator performs a bitwise inclusive or operation on all DDL bits
   * between this readout and <i>list</i>.
   * \param list  The right hand side readout list to operate on.
   * \return  A reference to this object.
   */
  AliHLTReadoutList& operator |= (const AliHLTReadoutList& list);

  /// same as operator |=
  AliHLTReadoutList& OrEq(const AliHLTReadoutList& list);
  
  /**
   * This operator performs a bitwise exclusive or (xor) operation on all DDL
   * bits between this readout and <i>list</i>.
   * \param list  The right hand side readout list to operate on.
   * \return  A reference to this object.
   */
  AliHLTReadoutList& operator ^= (const AliHLTReadoutList& list);

  /// same as operator ^=
  AliHLTReadoutList& XorEq(const AliHLTReadoutList& list);
  
  /**
   * This operator performs a bitwise and operation on all DDL bits between
   * this readout and <i>list</i>.
   * \param list  The right hand side readout list to operate on.
   * \return  A reference to this object.
   */
  AliHLTReadoutList& operator &= (const AliHLTReadoutList& list);

  /// same as operator &=
  AliHLTReadoutList& AndEq(const AliHLTReadoutList& list);
  
  /**
   * This operator performs the effective operation of "this and (this xor list)".
   * It removes all the DDLs specified in list from this readout list.
   * \param list  The right hand side readout list to operate on.
   * \return  A reference to this object.
   */
  AliHLTReadoutList& operator -= (const AliHLTReadoutList& list);
  
  /**
   * This operator performs a bitwise ones compliment on all DDL bits of this
   * readout list.
   * \return  The result of the unary operator.
   */
  AliHLTReadoutList operator ~ () const;
  
  /**
   * This operator performs a bitwise inclusive or operation on all DDL bits
   * between this readout and <i>list</i>.
   * \param list  The right hand side readout list to operate on.
   * \return  The result of the binary operator.
   */
  AliHLTReadoutList operator | (const AliHLTReadoutList& list) const
  {
    AliHLTReadoutList result = *this;
    return result.operator |= (list);
  }
  
  /**
   * This operator performs a bitwise exclusive or (xor) operation on all DDL
   * bits between this readout and <i>list</i>.
   * \param list  The right hand side readout list to operate on.
   * \return  The result of the binary operator.
   */
  AliHLTReadoutList operator ^ (const AliHLTReadoutList& list) const
  {
    AliHLTReadoutList result = *this;
    return result.operator ^= (list);
  }
  
  /**
   * This operator performs a bitwise and operation on all DDL bits between
   * this readout and <i>list</i>.
   * \param list  The right hand side readout list to operate on.
   * \return  The result of the binary operator.
   */
  AliHLTReadoutList operator & (const AliHLTReadoutList& list) const
  {
    AliHLTReadoutList result = *this;
    return result.operator &= (list);
  }
  
  /**
   * This operator performs the effective operation of "this and (this xor list)".
   * i.e. the set difference.
   * It removes all the DDLs specified in list from this readout list.
   * \param list  The right hand side readout list to operate on.
   * \return  The result of the binary operator.
   */
  AliHLTReadoutList operator - (const AliHLTReadoutList& list) const
  {
    AliHLTReadoutList result = *this;
    return result.operator -= (list);
  }
  
 private:
  
  /**
   * Decodes the word index and bit index within that word for the readout list structure.
   * \param ddlId <i>[in]</i>  The ID number of the DDL to decode.
   * \param wordIndex <i>[out]</i>  the word index of the word to modify or check
   *    within fReadoutList.fList
   * \param bitIndex <i>[out]</i>   the bit index of the bit to modify or check
   *    within the word pointed to by <i>wordIndex</i>.
   * \return  true if the ddlId was decoded and false if it was invalid.
   * \note We do not check extensively if the ddlId is invalid. Just simple checks
   *    are performed to see that we do not overflow the buffer fReadoutList.fList.
   */
  static bool DecodeDDLID(Int_t ddlId, Int_t& wordIndex, Int_t& bitIndex);
  
  AliHLTEventDDL fReadoutList; /// The DDL readout list structure.
  
  ClassDef(AliHLTReadoutList, 3) // Readout list object used for manipulating and storing an AliHLTEventDDL structure.

};

#endif // ALIHLTREADOUTLIST_H

