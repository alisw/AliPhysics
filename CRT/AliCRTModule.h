#ifndef ALICRTMODULE_H
#define ALICRTMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: */
/////////////////////////////////
// CRT module geometry manager //
/////////////////////////////////

#include <TNamed.h>

class AliCRTModule : public TNamed {
public:
  AliCRTModule();
  AliCRTModule(const char* name, const char* title);
  AliCRTModule(const AliCRTModule& mod);
  virtual ~AliCRTModule();

  AliCRTModule& operator=(const AliCRTModule& mod);

  void SetScintillatorThickness(const Float_t thickness);
  void SetScintillatorWidth(const Float_t width);
  void SetScintillatorLenght(const Float_t length);

  void SetFrameThickness(const Float_t thickness);
  void SetFrameWidth(const Float_t width);
  void SetFrameLength(const Float_t length);

  void SetNumberOfColumns(const Int_t ncols);
  void SetNumberOfRows(const Int_t nrows);

  void SetZGap(const Float_t zgap);
  void SetXGap(const Float_t xgap);

  const Float_t ScintillatorThickness() const;
  const Float_t ScintillatorWidth() const;
  const Float_t ScintillatorLenght() const;

  const Float_t FrameThickness() const;
  const Float_t FrameWidth() const;
  const Float_t FrameLength() const;

  const Int_t NumberOfModules() const;
  const Int_t NumberOfColumns() const;
  const Int_t NumberOfRows() const;

  const Float_t ZGap() const;
  const Float_t XGap() const;

private:
  Float_t fScintillatorThickness; // Scintillator thickness
  Float_t fScintillatorWidth; // Scintillator width
  Float_t fScintillatorLength; // Scintillator length
  Float_t fFrameThickness; // Aluminium frame thickness
  Float_t fFrameWidth; // Aluminium frame width
  Float_t fFrameLength; // Aliuminium frame length
  Int_t fNColumns;//Number of modules per column per magnet face (z coordinate)
  Int_t fNRows; // Number of module rows per magnet face (x coordinate)
  Float_t fZGap; // Gap in Z betwen modules
  Float_t fXGap; // Gap in X betwen modules
  ClassDef(AliCRTModule, 1)// CRT module geometry manager
};

inline void AliCRTModule::SetScintillatorThickness(const Float_t thick)
{ fScintillatorThickness = thick; }

inline void AliCRTModule::SetScintillatorWidth(const Float_t width)
{ fScintillatorWidth = width; }

inline void AliCRTModule::SetScintillatorLenght(const Float_t length)
{ fScintillatorLength = length; }

inline void AliCRTModule::SetFrameThickness(const Float_t thick)
{ fFrameThickness = thick; }

inline void AliCRTModule::SetFrameWidth(const Float_t width)
{ fFrameWidth = width; }

inline void AliCRTModule::SetFrameLength(const Float_t length)
{ fFrameLength = length; }

inline void AliCRTModule::SetNumberOfColumns(const Int_t ncols)
{ fNColumns = ncols; }

inline void AliCRTModule::SetNumberOfRows(const Int_t nrows)
{ fNRows = nrows; }

inline void AliCRTModule::SetZGap(const Float_t zgap)
{ fZGap = zgap; }

inline void AliCRTModule::SetXGap(const Float_t xgap)
{ fXGap = xgap; }

inline const Float_t AliCRTModule::ScintillatorThickness() const
{ return fScintillatorThickness; }

inline const Float_t AliCRTModule::ScintillatorWidth() const
{ return fScintillatorWidth; }

inline const Float_t AliCRTModule::ScintillatorLenght() const
{ return fScintillatorLength; }

inline const Float_t AliCRTModule::FrameThickness() const
{ return fFrameThickness; }

inline const Float_t AliCRTModule::FrameWidth() const
{ return fFrameWidth; }

inline const Float_t AliCRTModule::FrameLength() const
{ return fFrameLength; }

inline const Int_t AliCRTModule::NumberOfModules() const
{ return fNColumns*fNRows; }

inline const Int_t AliCRTModule::NumberOfColumns() const
{ return fNColumns; }

inline const Int_t AliCRTModule::NumberOfRows() const
{ return fNRows; }

inline const Float_t AliCRTModule::ZGap() const
{ return fZGap; }

inline const Float_t AliCRTModule::XGap() const
{ return fXGap; }
#endif // ALICRTMODULE_H
