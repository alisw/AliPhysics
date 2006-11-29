#ifndef ALIACORDEMODULE_H
#define ALIACORDEMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: */
/////////////////////////////////
// ACORDE module geometry manager //
/////////////////////////////////

#include <TNamed.h>

class AliACORDEModule : public TNamed {
public:
  AliACORDEModule();
  AliACORDEModule(const char* name, const char* title);
  AliACORDEModule(const AliACORDEModule& mod);
  virtual ~AliACORDEModule();

  AliACORDEModule& operator=(const AliACORDEModule& mod);

  void SetScintillatorThickness(Float_t thickness);
  void SetScintillatorWidth(Float_t width);
  void SetScintillatorLenght(Float_t length);

  void SetFrameThickness(Float_t thickness);
  void SetFrameWidth(Float_t width);
  void SetFrameLength(Float_t length);

  void SetNumberOfColumns(Int_t ncols);
  void SetNumberOfRows(Int_t nrows);

  void SetZGap(Float_t zgap);
  void SetXGap(Float_t xgap);

  Float_t ScintillatorThickness() const;
  Float_t ScintillatorWidth() const;
  Float_t ScintillatorLenght() const;

  Float_t FrameThickness() const;
  Float_t FrameWidth() const;
  Float_t FrameLength() const;

  Int_t NumberOfModules() const;
  Int_t NumberOfColumns() const;
  Int_t NumberOfRows() const;

  Float_t ZGap() const;
  Float_t XGap() const;

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
  ClassDef(AliACORDEModule, 1)// ACORDE module geometry manager
};

typedef AliACORDEModule AliCRTModule; // for backward compatibility

inline void AliACORDEModule::SetScintillatorThickness(Float_t thick)
{ fScintillatorThickness = thick; }

inline void AliACORDEModule::SetScintillatorWidth(Float_t width)
{ fScintillatorWidth = width; }

inline void AliACORDEModule::SetScintillatorLenght(Float_t length)
{ fScintillatorLength = length; }

inline void AliACORDEModule::SetFrameThickness(Float_t thick)
{ fFrameThickness = thick; }

inline void AliACORDEModule::SetFrameWidth(Float_t width)
{ fFrameWidth = width; }

inline void AliACORDEModule::SetFrameLength(Float_t length)
{ fFrameLength = length; }

inline void AliACORDEModule::SetNumberOfColumns(Int_t ncols)
{ fNColumns = ncols; }

inline void AliACORDEModule::SetNumberOfRows(Int_t nrows)
{ fNRows = nrows; }

inline void AliACORDEModule::SetZGap(Float_t zgap)
{ fZGap = zgap; }

inline void AliACORDEModule::SetXGap(Float_t xgap)
{ fXGap = xgap; }

inline Float_t AliACORDEModule::ScintillatorThickness() const
{ return fScintillatorThickness; }

inline Float_t AliACORDEModule::ScintillatorWidth() const
{ return fScintillatorWidth; }

inline Float_t AliACORDEModule::ScintillatorLenght() const
{ return fScintillatorLength; }

inline Float_t AliACORDEModule::FrameThickness() const
{ return fFrameThickness; }

inline Float_t AliACORDEModule::FrameWidth() const
{ return fFrameWidth; }

inline Float_t AliACORDEModule::FrameLength() const
{ return fFrameLength; }

inline Int_t AliACORDEModule::NumberOfModules() const
{ return fNColumns*fNRows; }

inline Int_t AliACORDEModule::NumberOfColumns() const
{ return fNColumns; }

inline Int_t AliACORDEModule::NumberOfRows() const
{ return fNRows; }

inline Float_t AliACORDEModule::ZGap() const
{ return fZGap; }

inline Float_t AliACORDEModule::XGap() const
{ return fXGap; }
#endif // ALIACORDEMODULE_H
