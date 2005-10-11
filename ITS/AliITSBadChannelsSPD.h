#ifndef ALIITSBADCHANNELSSPD_H
#define ALIITSBADCHANNELSSPD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
// AliITSBadChannelsSPD declaration by P. Nilsson 2005
// AUTHOR/CONTACT: Paul.Nilsson@cern.ch
//
// The class is used by the AliITSPreprocessorSPD class to store the
// final noisy and dead channel objects in the calibration database for
// the SPD.
//
// (See the source file for more information)
///////////////////////////////////////////////////////////////////////////


#include "TObjArray.h"
#include "AliITSChannelSPD.h"

class AliITSBadChannelsSPD: public TObject {

 public:

  AliITSBadChannelsSPD(void);		                // Default constructor
  AliITSBadChannelsSPD(const AliITSBadChannelsSPD &bc); // Default copy constructor
  virtual ~AliITSBadChannelsSPD(void);                  // Default destructor
  AliITSBadChannelsSPD& operator=(const AliITSBadChannelsSPD& bc); // Assignment operator

  void Put(Int_t* &array, const Int_t &arraySize,       // Add new arrays to the collection
	   Int_t* &index, const Int_t &indexSize);
  Bool_t Get(Int_t* &array, Int_t* &index) const;       // Retrieve the stored arrays (true if non empty arrays)

  Int_t GetIndexArraySize(void) const                   // Return the size of the index array
    { return fIndexArraySize; };
  Int_t GetBadChannelsArraySize(void) const             // Return the size of the bad channels array
    { return fBadChannelsArraySize; };

  Int_t* CreateModuleArray(const Int_t module) const;   // Create an array with sequential data for a given module
  Int_t GetModuleArraySize(const Int_t module) const    // Return array size for a given module
    { return (2*fBadChannelsArray[fIndexArray[module]] + 1); };

  TObjArray* CreateModuleObjArray(const Int_t module) const; // Create a TObjArray with data for a given module
  Int_t GetModuleObjArraySize(const Int_t module) const      // Return TObjArray size for a given module
    { return (fBadChannelsArray[fIndexArray[module]]); };

 protected:

  Int_t fIndexArraySize;                                // Size of the index array
  Int_t fBadChannelsArraySize;                          // Size of the bad channels array
  Int_t *fIndexArray;                                   //[fIndexArraySize]
  Int_t *fBadChannelsArray;                             //[fBadChannelsArraySize]

  ClassDef(AliITSBadChannelsSPD,1)
};

#endif
