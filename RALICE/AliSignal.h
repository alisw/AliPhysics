#ifndef ALISIGNAL_H
#define ALISIGNAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include "TObject.h"
#include "TArrayF.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TArrayI.h"

#include "AliPosition.h"
#include "AliAttrib.h"
#include "AliObjMatrix.h"

class AliSignal : public TNamed,public AliPosition,public AliAttrib
{
 public:
  AliSignal();                                          // Default constructor
  virtual ~AliSignal();                                 // Destructor
  AliSignal(AliSignal& s);                              // Copy constructor
  virtual TObject* Clone(const char* name="");          // Make a deep copy and provide its pointer
  virtual void SetSignal(Double_t sig,Int_t j=1);       // Store j-th signal value
  virtual void AddSignal(Double_t sig,Int_t j=1);       // Add value to j-th signal value
  virtual Float_t GetSignal(Int_t j=1,Int_t mode=0);    // Provide j-th (corrected) signal value
  virtual void SetSignalError(Double_t dsig,Int_t j=1); // Store error on j-th signal value
  virtual Float_t GetSignalError(Int_t j=1);            // Provide error j-th signal value
  virtual void ResetSignals(Int_t mode=0);              // User selected reset of signal values and/or errors
  virtual void DeleteSignals(Int_t mode=0);             // User selected delete of signal values and/or errors
  virtual void Reset(Int_t mode=0);                     // Reset signal and position values and errors
  virtual void Data(TString f="car");                   // Print all signal info for coord. frame f
  virtual void List(Int_t j=0);                         // Print signal info for the j-th (or all) slot(s)
  Int_t GetNvalues();                                   // Provide the number of signal values
  Int_t GetNerrors();                                   // Provide the number of specified errors
  Int_t GetNwaveforms();                                // Provide the number of specified waveforms
  void SetWaveform(TH1F* waveform,Int_t j=1);           // Set the j-th waveform histogram
  TH1F* GetWaveform(Int_t j=1);                         // Provide pointer of the j-th waveform histogram 
  void ResetWaveform(Int_t j=1);                        // Reset the j-th waveform histogram 
  void DeleteWaveform(Int_t j=1);                       // Delete the j-th waveform histogram 
  Int_t GetNlinks(TObject* obj=0,Int_t j=0);            // Provide the number of links for the specified object
  void SetLink(TObject* obj,Int_t j=1,Int_t k=1);       // Link object to the j-th slot at position k
  void AddLink(TObject* obj,Int_t j=1);                 // Link object to the j-th slot at first free position
  TObject* GetLink(Int_t j=1,Int_t k=1);                // Provide pointer of the object linked to the j-th slot
  Int_t GetIndices(TObject* obj,TArrayI& js,TArrayI& ks);// Provide slot and position indices for linked objects
  Int_t GetIndices(TObject* obj,Int_t j,TArrayI& ks);   // Provide pos. indices for linked objects of j-th slot 
  Int_t GetIndices(TObject* obj,TArrayI& js,Int_t k);   // Provide slot indices for linked objects at pos. k 
  void ResetLink(Int_t j=1,Int_t k=1);                  // Reset the link(s) of the j-th slot 
  void ResetLinks(TObject* obj,Int_t j=0,Int_t k=0);    // Reset link(s) to the specified object for j-th slot
  void SetSwapMode(Int_t swap=1);                       // Set swapmode flag for the link storage
  Int_t GetSwapMode();                                  // Provide swapmode flag for the link storage

 protected:
  TArrayF* fSignals;                           // Signal values
  TArrayF* fDsignals;                          // Errors on signal values
  TObjArray* fWaveforms;                       // The 1D histograms containing the signal waveforms
  AliObjMatrix* fLinks;                        // Pointers of objects related to the various slots

 ClassDef(AliSignal,11) // Generic handling of (extrapolated) detector signals.
};
#endif
