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

class AliDevice;
class AliTrack;

class AliSignal : public TNamed,public AliPosition,public AliAttrib
{
 public:
  AliSignal();                                                  // Default constructor
  virtual ~AliSignal();                                         // Destructor
  AliSignal(const AliSignal& s);                                // Copy constructor
  virtual TObject* Clone(const char* name="") const;            // Make a deep copy and provide its pointer
  virtual void SetSignal(Double_t sig,Int_t j=1);               // Store signal value for the j-th slot
  virtual void SetSignal(Double_t sig,TString name);            // Store signal value for the name-specified slot
  virtual void AddSignal(Double_t sig,Int_t j=1);               // Add value to the signal of the j-th slot
  virtual void AddSignal(Double_t sig,TString name);            // Add value to the signal of the name-specified slot
  virtual Float_t GetSignal(Int_t j=1,Int_t mode=0) const;      // Provide j-th (corrected) signal value
  virtual Float_t GetSignal(TString name,Int_t mode=0) const;   // Provide name-specified (corrected) signal value
  virtual void SetSignalError(Double_t dsig,Int_t j=1);         // Store error on the signal of the j-th slot
  virtual void SetSignalError(Double_t dsig,TString name);      // Store error on the signal of the name-specified slot
  virtual Float_t GetSignalError(Int_t j=1) const;              // Provide error on the signal of the j-th slot
  virtual Float_t GetSignalError(TString name) const;           // Provide err. on the sig. of the name-specified slot
  virtual void ResetSignals(Int_t mode=0);                      // Selective reset of signal values and/or errors
  virtual void DeleteSignals(Int_t mode=0);                     // Selectie delete of signal values and/or errors
  virtual void Reset(Int_t mode=0);                             // Reset signal and position values and errors
  virtual void Data(TString f="car",TString u="rad") const;     // Print all signal info for frame f and ang units u
  virtual void List(Int_t j=0) const;                           // Print signal info for the j-th (all) slot(s)
  virtual void List(TString name) const;                        // Print signal info for the name-specified slot
  void ListWaveform(Int_t j=0) const;                           // Print info for the j-th (all) waveform(s)
  void ListTrack(Int_t j=0) const;                              // Print info for the j-th (all) assoc. track(s)
  Int_t GetNvalues() const;                                     // Provide the number of signal values
  Int_t GetNerrors() const;                                     // Provide the number of specified errors
  virtual Int_t GetNslots() const;                              // Provide the number of exising slots
  Int_t GetNwaveforms() const;                                  // Provide the number of specified waveforms
  void SetWaveform(TH1F* waveform,Int_t j=1);                   // Set the histogram for the j-th waveform
  TH1F* GetWaveform(Int_t j=1) const;                           // Pointer to the histo of the j-th waveform
  TH1F* GetWaveform(TString name) const;                        // Pointer to the waveform with the specified name
  Int_t GetWaveformIndex(TString name) const;                   // Index of the waveform with the specified name
  void ResetWaveform(Int_t j=1);                                // Reset the histo of the j-th waveform
  void ResetWaveform(TString name);                             // Reset the waveform histo with the specified name
  void DeleteWaveform(Int_t j=1);                               // Delete histo of the j-th waveform
  void DeleteWaveform(TString name);                            // Delete waveform histo with the specified name
  Int_t GetNlinks(TObject* obj=0,Int_t j=0) const;              // Number of links for the specified object
  Int_t GetNlinks(TObject* obj,TString name) const;             // Number of links for the specified object
  void SetLink(TObject* obj,Int_t j=1,Int_t k=1);               // Link object to the j-th slot at position k
  void SetLink(TObject* obj,TString name,Int_t k=1);            // Link object to the name-specified slot at pos. k
  void AddLink(TObject* obj,Int_t j=1);                         // Link obj to the j-th slot at 1st free position
  void AddLink(TObject* obj,TString name);                      // Link obj to the name-specified slot at 1st free pos.
  TObject* GetLink(Int_t j=1,Int_t k=1) const;                  // Pointer of the object linked to the j-th slot
  TObject* GetLink(TString name,Int_t k=1) const;               // Pointer of object linked to the name-specified slot
  Int_t GetIndices(TObject* obj,TArrayI& js,TArrayI& ks) const; // Slot and pos. indices for linked objects
  Int_t GetIndices(TObject* obj,Int_t j,TArrayI& ks) const;     // Pos. indices for linked objects of j-th slot 
  Int_t GetIndices(TObject* obj,TString name,TArrayI& ks) const;// Pos. indices for linked objects of name-spec. slot 
  Int_t GetIndices(TObject* obj,TArrayI& js,Int_t k) const;     // Slot indices for linked objects at pos. k 
  void ResetLink(Int_t j=1,Int_t k=1);                          // Reset the link(s) of the j-th slot 
  void ResetLink(TString name,Int_t k=1);                       // Reset the link(s) of the name-specified slot 
  void ResetLinks(TObject* obj,Int_t j=0,Int_t k=0);            // Reset link(s) to object obj for j-th slot
  void ResetLinks(TObject* obj,TString name,Int_t k=0);         // Reset link(s) to object obj for name-specified slot
  void SetSwapMode(Int_t swap=1);                               // Set swapmode flag for the link storage
  Int_t GetSwapMode() const;                                    // Provide swapmode flag for the link storage
  void SetDevice(TObject* dev);                                 // Store pointer to the device that owns this signal
  AliDevice* GetDevice() const;                                 // Provide pointer to the owning device 
  void AddTrack(AliTrack& t,Int_t mode=1);                      // Relate an AliTrack to this signal
  void RemoveTrack(AliTrack& t,Int_t mode=1);                   // Remove related AliTrack from this signal
  void RemoveTracks(Int_t mode=1);                              // Remove all related AliTracks from this signal
  Int_t GetNtracks(AliTrack* t=0) const;                        // Provide number of related AliTracks
  AliTrack* GetTrack(Int_t j) const;                            // Access to the related AliTrack number j
  AliTrack* GetIdTrack(Int_t id) const;                         // Access to the related AliTrack with ID=id
  Int_t GetSignalFlag(Int_t j=1) const;                         // Provide signal value flag of the j-th slot
  Int_t GetSignalFlag(TString name) const;                      // Provide signal value flag of the name specified slot
  Int_t GetErrorFlag(Int_t j=1) const;                          // Provide error value flag of the j-th slot
  Int_t GetErrorFlag(TString name) const;                       // Provide error value flag of the name specified slot

 protected:
  void SetSigFlags(Int_t is,Int_t ie,Int_t j); // Set flags for signal and/or error value settings
  TArrayF* fSignals;                           // Signal values
  TArrayF* fDsignals;                          // Errors on signal values
  TArrayI* fSigflags;                          // Flags to mark setting of signal and/or error values 
  TObjArray* fWaveforms;                       // The 1D histograms containing the signal waveforms
  AliObjMatrix* fLinks;                        // Pointers of objects related to the various slots
  TObject* fDevice;                            // Pointer to the device that owns this signal
  TObjArray* fTracks;                          // Pointers to associated tracks

 ClassDef(AliSignal,16) // Generic handling of (extrapolated) detector signals.
};
#endif
