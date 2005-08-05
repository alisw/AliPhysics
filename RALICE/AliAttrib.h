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
#include "TF1.h"

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
  void SetGain(Double_t gain,TString name);    // Set gain of the name-specified attribute slot
  Float_t GetGain(Int_t j=1) const;            // Provide gain of the j-th attribute slot
  Float_t GetGain(TString name) const;         // Provide gain of the name-specified attribute slot
  void SetOffset(Double_t off,Int_t j=1);      // Set offset value of the j-th attribute slot
  void SetOffset(Double_t off,TString name);   // Set offset value of the name specified attribute slot
  Float_t GetOffset(Int_t j=1) const;          // Provide offset value of the j-th attribute slot
  Float_t GetOffset(TString name) const;       // Provide offset value of the name-specified attribute slot
  Int_t GetGainFlag(Int_t j=1) const;          // Provide gain flag of the j-th attribute slot
  Int_t GetGainFlag(TString name) const;       // Provide gain flag of the name-specified attribute slot
  Int_t GetOffsetFlag(Int_t j=1) const;        // Provide offset flag of the j-th attribute slot
  Int_t GetOffsetFlag(TString name) const;     // Provide offset flag of the name-specified attribute slot
  void ResetGain(Int_t j=1);                   // Reset j-th gain value and flag
  void ResetGain(TString name);                // Reset name-specified gain value and flag
  void ResetOffset(Int_t j=1);                 // Reset j-th offset value and flag
  void ResetOffset(TString name);              // Reset name-specified offset value and flag
  void DeleteCalibrations(Int_t mode=0);       // User selected delete of gains and/or offsets
  void SetDead(Int_t j=1);                     // Indicate j-th attribute slot as 'dead'
  void SetDead(TString name);                  // Indicate name-specified attribute slot as 'dead'
  void SetAlive(Int_t j=1);                    // Indicate j-th attribute slot as 'active'
  void SetAlive(TString name);                 // Indicate name-specified attribute slot as 'active'
  Int_t GetDeadValue(Int_t j=1) const;         // Return the 'dead flag' of the j-th attribute slot
  Int_t GetDeadValue(TString name) const;      // Return the 'dead flag' of the name-specified attribute slot
  void SetEdgeOn(Int_t j=1);                   // Indicate j-th slot as 'detector edge module'
  void SetEdgeOn(TString name);                // Indicate name-spcified slot as 'detector edge module'
  void SetEdgeOff(Int_t j=1);                  // Indicate j-th slot as 'detector non-edge module'
  void SetEdgeOff(TString name);               // Indicate name-specified slot as 'detector non-edge module'
  void IncreaseEdgeValue(Int_t j=1);           // Increase the edge value of the j-th slot by 1
  void IncreaseEdgeValue(TString name);        // Increase the edge value of the name-specified slot by 1
  void DecreaseEdgeValue(Int_t j=1);           // Decrease the edge value of the j-th slot by 1
  void DecreaseEdgeValue(TString name);        // Decrease the edge value of the name-specified slot by 1
  void SetEdgeValue(Int_t val,Int_t j=1);      // Set a specific edge value for the j-th slot
  void SetEdgeValue(Int_t val,TString name);   // Set a specific edge value for the name-specified slot
  Int_t GetEdgeValue(Int_t j=1) const;         // Provide the edge value of the j-th slot
  Int_t GetEdgeValue(TString name) const;      // Provide the edge value of the name-specified slot
  virtual void List(Int_t j=0) const;          // Printout of attribute data
  virtual void List(TString name) const;       // Printout of name-specified attribute data
  virtual void Load(AliAttrib& a,Int_t j=0);   // Load j-th slot or all attributes of the input AliAttrib
  virtual void Load(AliAttrib& a,TString name);// Load name-specified slot attributes of the input AliAttrib
  void SetSlotName(TString s,Int_t j=1);       // Set user defined name for the j-th slot
  TString GetSlotName(Int_t j=1) const;        // Provide user defined name for the j-th slot
  Int_t GetSlotIndex(TString name) const;      // Provide the slot index of the matching name
  void SetCalFunction(TF1* f,Int_t j=1);       // Set calibration function of the j-th attribute slot
  void SetCalFunction(TF1* f,TString name);    // Set calibration function of the name-specified attribute slot
  void SetDecalFunction(TF1* f,Int_t j=1);     // Set de-calibration function of the j-th attribute slot
  void SetDecalFunction(TF1* f,TString name);  // Set de-calibration function of the name-specified attribute slot
  TF1* GetCalFunction(Int_t j=1) const;        // Get calibration function of the j-th attribute slot
  TF1* GetCalFunction(TString name) const;     // Get calibration function of the name-specified attribute slot
  TF1* GetDecalFunction(Int_t j=1) const;      // Get de-calibration function of the j-th attribute slot
  TF1* GetDecalFunction(TString name) const;   // Get de-calibration function of the name-specified attribute slot
  Int_t GetNcalfuncs() const;                  // Provide the number of calibration functions
  Int_t GetNdecalfuncs() const;                // Provide the number of de-calibration functions

 protected:
  void SetCalFlags(Int_t gf,Int_t of,Int_t j); // Set flags for gain and/or offset settings
  TArrayF* fGains;                             // Gain values
  TArrayF* fOffsets;                           // Offset values
  TArrayI* fCalflags;                          // Flags to mark dead, edge, and gain/offset calibrated signals
  TObjArray* fNames;                           // User defined names for the various slots
  TObjArray* fCalfuncs;                        // Explicit signal calibration functions
  TObjArray* fDecalfuncs;                      // Explicit signal de-calibration functions

 ClassDef(AliAttrib,4) // Generic handling of detector signal (calibration) attributes.
};
#endif
