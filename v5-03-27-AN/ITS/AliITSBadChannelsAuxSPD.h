#ifndef ALIITSBADCHANNELSAUXSPD_H
#define ALIITSBADCHANNELSAUXSPD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
// AliITSBadChannelsAuxSPD declaration by P. Nilsson 2005
// AUTHOR/CONTACT: Paul.Nilsson@cern.ch                  
//
// Auxiliary algorithms for the SPD
/////////////////////////////////////////////////////////


#include <TObjArray.h>
#include <TString.h>
#include "AliITSChannelSPD.h"
#include "AliITSdigitSPD.h"

class AliITSBadChannelsAuxSPD {

 public:

  AliITSBadChannelsAuxSPD(void);    			      // Default constructor
  virtual ~AliITSBadChannelsAuxSPD(void) { };                 // Default destructor

  // General algorithms
  Bool_t Diff(TObjArray *&in1, TObjArray *&in2, TObjArray *&out1, TObjArray *&out2) const; // Diff algorithm
  Bool_t Find(AliITSChannelSPD *&channel, TObjArray *&array) const; // Find a channel in the array
  Bool_t Find(AliITSdigitSPD *&digit, TObjArray *&array) const;     // Find a digit in the array
  Int_t GetNumberOfBadChannels(Int_t* &array, Int_t* &indexArray, Int_t size) const; // Get the number of bad channels

  // Converters
  AliITSdigitSPD* CreateDigitFromChannel(const AliITSChannelSPD *&channel) const; // Create a digit from a channel
  AliITSChannelSPD* CreateChannelFromDigit(const AliITSdigitSPD *&digit) const;   // Create a channel from a digit

  // Miscellanious
  Bool_t CreateHTMLReport(char *name, Int_t* &array, Int_t* &indexArray,    // Create an HTML report
			  Int_t indexSize, TString *buffer, Bool_t tags);

 protected:

  ClassDef(AliITSBadChannelsAuxSPD,1)
};

#endif
