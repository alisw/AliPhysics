#ifndef ALIDEVICE_H
#define ALIDEVICE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include "TPolyMarker3D.h"

#include "AliSignal.h"

class AliDevice : public AliSignal
{
 public:
  AliDevice();                                       // Default constructor
  virtual ~AliDevice();                              // Default destructor
  AliDevice(const AliDevice& dev);                   // Copy constructor
  virtual TObject* Clone(const char* name="") const; // Make a deep copy and provide its pointer
  void SetHitCopy(Int_t j);                          // (De)activate creation of private copies of hits
  Int_t GetHitCopy() const;                          // Provide HitCopy flag value      
  void AddHit(AliSignal& s);                         // Register an AliSignal object as a hit to this module
  void AddHit(AliSignal* s) { if (s) AddHit(*s); }
  void RemoveHit(AliSignal& s);                      // Remove AliSignal object as hit from this module
  void RemoveHit(AliSignal* s) { if (s) RemoveHit(*s); }
  void RemoveHits();                                 // Remove all AliSignals as hits from this module
  Int_t GetNhits() const;                            // Provide number of registered hits
  AliSignal* GetHit(Int_t j) const;                  // Access to the AliSignal registered as hit number j
  AliSignal* GetIdHit(Int_t id) const;               // Provide the hit with unique identifier "id"
  TObjArray* GetHits();                              // Provide the references to all the registered hits
  virtual void Reset(Int_t mode=0);                  // Reset registered hits and AliSignal attributes
  void ShowHit(Int_t j=0) const;                     // Show data of the j-th hit (j=0 means all hits)
  virtual void Data(TString f="car") const;          // Print device and all signal info for coord. frame f
  void GetExtremes(Float_t& vmin,Float_t& vmax,Int_t idx=1,TObjArray* hits=0) const; // Get min. and max. signal value
  void GetExtremes(Float_t& vmin,Float_t& vmax,TString name,TObjArray* hits=0) const;// Get min. and max. signal value
  TObjArray* SortHits(TString name,Int_t mode=-1,TObjArray* hits=0); // Sort hits by named signal value
  TObjArray* SortHits(Int_t idx=1,Int_t mode=-1,TObjArray* hits=0);  // Sort hits by indexed signal value
  void DisplayHits(TString name,Float_t scale=-1,TObjArray* hits=0,Int_t dp=0,Int_t mstyle=8,Int_t mcol=4);// Hit disp.
  void DisplayHits(Int_t idx=1,Float_t scale=-1,TObjArray* hits=0,Int_t dp=0,Int_t mstyle=8,Int_t mcol=4); // Hit disp.

 protected:
  Int_t fHitCopy;      // Flag to denote making private copies of added hits
  TObjArray* fHits;    // Array to hold the registered hits
  TObjArray* fOrdered; //! Temp. array to hold the ordered hits
  TObjArray* fMarkers; //! Temp. array to hold the 3D markers for the hit display

 ClassDef(AliDevice,5) // Signal (Hit) handling of a generic device.
};
#endif
