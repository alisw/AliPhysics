#ifndef ALIMUONPAINTERMATRIX_H
#define ALIMUONPAINTERMATRIX_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterMatrix
/// \brief A matrix of AliMUONVPainter
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif
#ifndef ALIMUONATTPAINTER_H
#  include "AliMUONAttPainter.h"
#endif

class AliMUONVPainter;
class AliMUONVTrackerData;
class TObjArray;

class AliMUONPainterMatrix : public TObject
{
public:
  AliMUONPainterMatrix(const char* basename="", Int_t nx=1, Int_t ny=1);
  virtual ~AliMUONPainterMatrix();
  
  /// Adopt a painter in this matrix
  void Adopt(AliMUONVPainter* painter);
  
  using TObject::Clone;
  
  AliMUONPainterMatrix* Clone(const AliMUONAttPainter& attributes) const;

  void Connect(const char* sourceMethod, const char* destClassName, 
               void* destObject, const char* destMethod);
    
//  void ChangeAttributes(const AliMUONAttPainter& attributes);
  
  /// Get our attributes
  const AliMUONAttPainter& Attributes() const { return fAttributes; }
  
  /// Compute the data range for this matrix
  void ComputeDataRange();
  
  /// Get the data range for this matrix
  void GetDataRange(Double_t& dataMin, Double_t& dataMax) const;
  
  /// Get our name
  virtual const char* GetName() const { return Name().Data(); }

  static TString NameIt(const TString& basename, const AliMUONAttPainter& att);

  /// Matrix name
  virtual TString Name() const { return fName; }
  
  /// Base name (short name)
  virtual TString Basename() const { return fBasename; }
  
  void GetTypes(TObjArray& types) const;

  /// Number of painters to arrange in x-direction
  Int_t Nx() const { return fNx; }
  
  /// Number of painters to arrange in y-direction
  Int_t Ny() const { return fNy; }
  
  /// Get a painter
  AliMUONVPainter* Painter(Int_t index) const;

  /// Printout
  void Print(Option_t* opt="") const;
  
  AliMUONVTrackerData* Data() const;
  
  TString DataPattern() const;
  
  Int_t DataIndex() const;
  
  void SetData(const char* pattern, AliMUONVTrackerData* d, Int_t indexInData);
    
  /// Force a given data range for all painter groups belonging to this matrix
  void SetDataRange(Double_t min, Double_t max);
  
  void SetOutlined(const char* pattern, Bool_t value);

  void SetResponder(const char* pattern);
  
  /// Number of painters (should be <= Nx*Ny)
  Int_t Size() const;

  /// Normalize attributes
  AliMUONAttPainter Validate(const AliMUONAttPainter& att) const;
  
private:
  /// Not implemented
  AliMUONPainterMatrix(const AliMUONPainterMatrix& rhs);
  /// Not implemented
  AliMUONPainterMatrix& operator=(const AliMUONPainterMatrix& rhs);

  void UpdateAttributes();
  
private:
  TString fBasename; ///< base name of that matrix
  TString fName; ///< complete name
  Int_t fNx; ///< number of rows
  Int_t fNy; ///< number of columns
  TObjArray* fPainters; ///< painters in that matrix
  AliMUONAttPainter fAttributes; ///< attributes of our painter(s)
  
  ClassDef(AliMUONPainterMatrix,1) // Matrix of AliMUONVPainter
};

#endif
