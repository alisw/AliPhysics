#ifndef ALIATTRIB_H
#define ALIATTRIB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"

class AliAttrib
{
 public:
  AliAttrib();                                 // Default constructor
  virtual ~AliAttrib();                        // Destructor
  AliAttrib(const AliAttrib& a);               // Copy constructor
  Int_t GetNgains() const;                     // Provide the number of specified gains
  Int_t GetNoffsets() const;                   // Provide the number of specified offsets
  Int_t GetNcalflags() const;                  // Provide the number of specified calib. flags
  Int_t GetNnames() const;                     // Provide the maximum number of specified names
  void SetGain(Double_t gain,Int_t j=1);       // Set gain of the j-th attribute slot
  Float_t GetGain(Int_t j=1) const;            // Provide gain of the j-th attribute slot
  void SetOffset(Double_t off,Int_t j=1);      // Set offset value of the j-th attribute slot
  Float_t GetOffset(Int_t j=1) const;          // Provide offset value of the j-th attribute slot
  Int_t GetGainFlag(Int_t j=1) const;          // Provide gain flag of the j-th attribute slot
  Int_t GetOffsetFlag(Int_t j=1) const;        // Provide offset flag of the j-th attribute slot
  void ResetGain(Int_t j=1);                   // Reset j-th gain value and flag
  void ResetOffset(Int_t j=1);                 // Reset j-th offset value and flag
  void DeleteCalibrations(Int_t mode=0);       // User selected delete of gains and/or offsets
  void SetDead(Int_t j=1);                     // Indicate j-th attribute slot as 'dead'
  void SetAlive(Int_t j=1);                    // Indicate j-th attribute slot as 'active'
  Int_t GetDeadValue(Int_t j=1) const;         // Return the 'dead flag' of the j-th attribute slot
  void SetEdgeOn(Int_t j=1);                   // Indicate j-th slot as 'detector edge module'
  void SetEdgeOff(Int_t j=1);                  // Indicate j-th slot as 'detector non-edge module'
  void IncreaseEdgeValue(Int_t j=1);           // Increase the edge value of the j-th slot by 1
  void DecreaseEdgeValue(Int_t j=1);           // Decrease the edge value of the j-th slot by 1
  void SetEdgeValue(Int_t val,Int_t j=1);      // Set a specific edge value for the j-th slot
  Int_t GetEdgeValue(Int_t j=1) const;         // Provide the edge value of the j-th slot
  virtual void List(Int_t j=0) const;          // Printout of attribute data
  virtual void Load(AliAttrib& a,Int_t j=0);   // Load j-th slot or all attributes of the input AliAttrib
  void SetSlotName(TString s,Int_t j=1);       // Set user defined name for the j-th slot
  TString GetSlotName(Int_t j=1) const;        // Provide user defined name for the j-th slot
  Int_t GetSlotIndex(TString name) const;      // Provide the slot index of the matching name

 protected:
  void SetCalFlags(Int_t gf,Int_t of,Int_t j); // Set flags for gain and/or offset settings
  TArrayF* fGains;                             // Gain values
  TArrayF* fOffsets;                           // Offset values
  TArrayI* fCalflags;                          // Flags to mark dead, edge, and gain/offset calibrated signals
  TObjArray* fNames;                           // User defined names for the various slots

 ClassDef(AliAttrib,2) // Generic handling of detector signal (calibration) attributes.
};
#endif
