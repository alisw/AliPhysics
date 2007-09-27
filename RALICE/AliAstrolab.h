#ifndef ALIASTROLAB_H
#define ALIASTROLAB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

#include <math.h>

#include "TTask.h"
#include "TString.h"
#include "TRotMatrix.h"
#include "TObjArray.h"
#include "TArrayI.h"

#include "AliTimestamp.h"
#include "AliPosition.h"
#include "AliSignal.h"
 
class AliAstrolab : public TTask,public AliTimestamp
{
 public:
  AliAstrolab(const char* name="AliAstrolab",const char* title="Generic lab"); // Constructor
  virtual ~AliAstrolab();                                      // Destructor
  AliAstrolab(const AliAstrolab& t);                           // Copy constructor
  virtual TObject* Clone(const char* name="") const;           // Make a deep copy and provide its pointer
  void Data(Int_t mode=1,TString u="deg");                     // Lab info in angular units u
  void SetLabPosition(Ali3Vector& r);                          // Set lab position in terrestrial frame
  void SetLabPosition(Double_t l,Double_t b,TString u="deg");  // Set lab terrestrial position
  AliPosition GetLabPosition() const;                          // Provide the lab terrestrial position 
  void GetLabPosition(Double_t& l,Double_t& b,TString u="deg") const;// Provide the lab terrestrial position
  using AliTimestamp::GetLT;
  Double_t GetLT();  // Provide Local Time (LT) in fractional hours
  using AliTimestamp::GetLMST;
  Double_t GetLMST(); // Provide Local Mean Sidereal Time (LMST) in fractional hours
  using AliTimestamp::GetLAST;
  Double_t GetLAST(); // Provide Local Apparent Sidereal Time (LAST) in fractional hours
  using AliTimestamp::SetLT;
  void SetLT(Int_t y,Int_t m,Int_t d,Int_t hh,Int_t mm,Int_t ss,Int_t ns=0,Int_t ps=0); // Set specified LT
  void SetLT(Int_t y,Int_t d,Int_t s,Int_t ns=0,Int_t ps=0); // Set LT based on elapsed days, secs etc...
  Double_t ConvertAngle(Double_t a,TString in,TString out) const;       // Angular format conversions
  void PrintAngle(Double_t a,TString in,TString out,Int_t ndig=1) const;// Print angle in various formats
  void SetSignal(Ali3Vector* r,TString frame,TString mode,AliTimestamp* ts,Int_t jref=0,TString name=""); // Store a generic signal
  void SetSignal(Double_t a,Double_t d,TString s,Double_t e,TString mode,Int_t jref=0,TString name="");   // Store RA, decl. and time
  void SetSignal(Double_t a,Double_t d,TString mode,AliTimestamp* ts,Int_t jref=0,TString name="");       // Store RA, decl. and time
  AliSignal* GetSignal(Ali3Vector& r,TString frame,TString mode,AliTimestamp* ts,Int_t jref=0);// Provide stored signal data
  AliSignal* GetSignal(Ali3Vector& r,TString frame,TString mode,AliTimestamp* ts,TString name);// Provide stored signal data
  AliSignal* GetSignal(Double_t& a,Double_t& d,TString mode,AliTimestamp* ts,Int_t jref=0);    // Provide corrected RA and decl.
  AliSignal* GetSignal(Double_t& a,Double_t& d,TString mode,AliTimestamp* ts,TString name);    // Provide corrected RA and decl.
  AliSignal* GetSignal(Double_t& a,Double_t& d,TString s,Double_t e,TString mode,Int_t jref=0);// Provide corrected RA and decl.
  AliSignal* GetSignal(Double_t& a,Double_t& d,TString s,Double_t e,TString mode,TString name);// Provide corrected RA and decl.
  AliSignal* GetSignal(Int_t jref=0);                // Provide pointer to a stored signal object
  AliSignal* GetSignal(TString name);                // Provide pointer to a stored signal object
  void RemoveRefSignal(Int_t j,Int_t compress);      // Remove a stored reference signal object
  void RemoveRefSignal(TString name,Int_t compress); // Remove a stored reference signal object
  void PrintSignal(TString frame,TString mode,AliTimestamp* ts,Int_t ndig,Int_t jref=0); // Print stored signal data
  void PrintSignal(TString frame,TString mode,AliTimestamp* ts,Int_t ndig,TString name); // Print stored signal data
  void ListSignals(TString frame,TString mode,Int_t ndig=1); // List all stored signals
  Int_t GetSignalIndex(TString name); // Provide storage index of the signal with the specified name
  Double_t GetHourAngle(TString mode,AliTimestamp* ts,Int_t jref=0);// Provide the Local Hour Angle in degrees
  void SetLocalFrame(Double_t t1,Double_t p1,Double_t t2,Double_t p2,Double_t t3,Double_t p3); // Define local coordinate frame
  using AliTimestamp::GetDifference;
  Double_t GetDifference(Int_t jref,TString au,Double_t& dt,TString tu,Int_t mode=1,Int_t* ia=0,Int_t* it=0); // Provide space and time difference
  Double_t GetDifference(TString name,TString au,Double_t& dt,TString tu,Int_t mode=1);// Provide space and time difference
  TArrayI* MatchRefSignal(Double_t da,TString au,Double_t dt,TString tu,Int_t mode=1); // Provide space and time matching reference signals
 
 protected:
  AliPosition fLabPos;   // Position of the lab in the terrestrial longitude-latitude frame
  Double_t fToffset;     // Lab time offset in fractional hours w.r.t. UT
  AliSignal* fXsig;      // Signal entry for object or event studies
  TObjArray* fRefs;      // Array holding the reference signals
  TRotMatrix fB;         //! The frame bias matrix for conversion of ICRS to J2000 coordinates
  Int_t fBias;           //! Initialisation flag for fB values (0=uninitialised  1=initialised)
  TRotMatrix fP;         //! Matrix for precession correction  
  TRotMatrix fN;         //! Matrix for nutation correction  
  TRotMatrix fG;         //! Matrix for conversion of equatorial to galactic coordinates
  Int_t fGal;            //! Type indicator for fG values (0=uninitialised  1=B1950  2=J2000)
  TRotMatrix fE;         //! Matrix for conversion of equatorial to ecliptic coordinates
  TRotMatrix fH;         //! Matrix for conversion of equatorial to horizontal coordinates
  TRotMatrix fL;         //! Matrix for conversion of horizontal to local-frame coordinates
  TArrayI* fIndices;     //! Storage indices of the matching reference signals
  void SetBmatrix();                 // Set the frame bias matrix
  void SetPmatrix(AliTimestamp* ts); // Set precession matrix for Julian date jd w.r.t. J2000.
  void SetNmatrix(AliTimestamp* ts); // Set nutation matrix for Julian date jd w.r.t. J2000.
  void SetGmatrix(TString mode);     // Set the equatorial to galactic conversion matrix
  void SetEmatrix(AliTimestamp* ts); // Set the equatorial to ecliptic conversion matrix
  void SetHmatrix(AliTimestamp* ts); // Set the equatorial to horizontal conversion matrix
  void Precess(Ali3Vector& r,AliTimestamp* ts1,AliTimestamp* ts2); // Correct RA and decl. for earth's precession
  void Nutate(Ali3Vector& r,AliTimestamp* ts); // Correct RA and decl. for earth's nutation
 
 ClassDef(AliAstrolab,2) // Virtual lab to relate measurements with astrophysical phenomena
};
#endif
