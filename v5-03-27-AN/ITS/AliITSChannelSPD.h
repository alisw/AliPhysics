#ifndef ALIITSCHANNELSPD_H
#define ALIITSCHANNELSPD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
// AliITSChannelSPD declaration by P. Nilsson 2005
// AUTHOR/CONTACT: Paul.Nilsson@cern.ch
//
// Objects of this class are stored in a TObjArray and should be
// interpreted as "bad" channels, i.e. either noisy or dead channels
// depending on where they are stored
///////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliITSChannelSPD: public TObject {

 public:

  AliITSChannelSPD(void);    			            // Default constructor
  AliITSChannelSPD(Int_t column, Int_t row);    // Constructor for already existing "bad" channel
  AliITSChannelSPD(const AliITSChannelSPD &ch);             // Copy constructor
  virtual ~AliITSChannelSPD(void) { };                      // Default destructor
  AliITSChannelSPD& operator=(const AliITSChannelSPD &ch);  // Assignment operator
  Bool_t operator==(const AliITSChannelSPD &channel) const; // Equivalence operator

  // Getters and setters
  Int_t GetColumn(void) const { return fColumn; };          // Get column
  Int_t GetRow(void) const { return fRow; };                // Get row
  void SetColumn(Int_t c) { fColumn = c; };                 // Set column
  void SetRow(Int_t r) { fRow = r; };                       // Set row

 protected:

  Int_t fColumn;                                            // SPD column (real range [0,31], but not checked)
  Int_t fRow;                                               // SPD row (real range [0,255] but not checked)

  ClassDef(AliITSChannelSPD,1)
};

#endif
