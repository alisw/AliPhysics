#ifndef ALIMUONALIGNMENTRECORD_H
#define ALIMUONALIGNMENTRECORD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup rec
/// \class AliMUONAlignmentClusterRecord
/// \brief Class to store alignment local and global derivatives of muon spectrometer
//
// Author: Hugo Pereira Da Costa

#include <TClonesArray.h>
#include <TObject.h>
#include <cassert>

class AliMUONAlignmentClusterRecord:public TObject
{

public:

  /// constructor
  AliMUONAlignmentClusterRecord( void );

  /// destructor
  virtual ~AliMUONAlignmentClusterRecord( void )
  {}

  /// Print
  virtual void	Print(Option_t* = "") const;

  ///@name modifiers
  ///@{

  /// detetion element Id
  void SetDetElemId( Int_t value )
  {
    fDetElemId = value;
    fDetElemNumber = fDetElemId%100;
  }

  /// detection element number
  void SetDetElemNumber( Int_t value )
  { fDetElemNumber = value; }

  /// measurement
  void SetMeas( Double_t value )
  { fMeas = value; }

  /// error on measurement
  void SetSigma( Double_t value )
  { fSigma = value; }

  /// local derivative
  void SetLocalDerivative( Int_t index, Double_t value )
  {
    // todo: bound check ?
    fLocalDerivatives[index] = value;
  }

  /// global derivative
  void SetGlobalDerivative( Int_t index, Double_t value )
  {
    // todo: bound check ?
    fGlobalDerivatives[index] = value;
  }

  ///@}

  ///@name accessors
  ///@{

  /// detection element id
  Int_t GetDetElemId( void ) const
  { return fDetElemId; }

  /// detection element number
  Int_t GetDetElemNumber( void ) const
  { return fDetElemNumber; }

  /// measurement
  Double_t GetMeas( void ) const
  { return fMeas; }

  /// error on measurement
  Double_t GetSigma( void ) const
  { return fSigma; }

  /// local derivative
  Double_t GetLocalDerivative( Int_t index ) const
  {
    // todo: bound check ?
    return fLocalDerivatives[index];
  }

  /// global derivative
  Double_t GetGlobalDerivative( Int_t index ) const
  {
    // todo: bound check ?
    return fGlobalDerivatives[index];
  }

  /// @}

private:

  /// Detection element Id
  Int_t fDetElemId;

  /// Detection element number
  Int_t fDetElemNumber;

  /// measurement
  Double_t fMeas;

  /// error on measurement
  Double_t fSigma;

  /// local derivatives
  Double_t fLocalDerivatives[4];

  /// global derivatives
  Double_t fGlobalDerivatives[4];

ClassDef(AliMUONAlignmentClusterRecord, 1)

};

/// \ingroup rec
/// \class AliMUONAlignmentTrackRecord
/// \brief Class to store alignment local and global derivatives of muon spectrometer
//
// Author: Hugo Pereira Da Costa
class AliMUONAlignmentTrackRecord:public TObject
{

public:

  /// constructor
  AliMUONAlignmentTrackRecord( void );

  /// destructor
  virtual ~AliMUONAlignmentTrackRecord( void );

  /// copy constructor
  AliMUONAlignmentTrackRecord (const AliMUONAlignmentTrackRecord& );

  /// assignment operator
  AliMUONAlignmentTrackRecord& operator=(const AliMUONAlignmentTrackRecord& );

  /// Print
  virtual void	Print(Option_t* option = "") const;

  ///@name accessors
  //@{

  /// cluster records
  TClonesArray* GetClusterRecords( void ) const
  { return fClusterRecords; }

  /// number of records
  Int_t GetNRecords( void ) const
  { return fClusterCount; }

  /// return record for given index
  AliMUONAlignmentClusterRecord* GetRecord( Int_t index ) const
  { return fClusterRecords ? static_cast<AliMUONAlignmentClusterRecord*>(fClusterRecords->UncheckedAt( index )) : 0x0; }

  //@}

  ///@name modifiers
  //@{

  /// add cluster record
  void AddClusterRecord( const AliMUONAlignmentClusterRecord& );

  /// remove cluster record
  void RemoveClusterRecord( AliMUONAlignmentClusterRecord* );

  /// clear memory
  virtual void Clear( Option_t* ="" );

  //@}

private:

  /// default clonesarray size
  enum { fSize = 20 };

  /// alignment parameters at cluster
  TClonesArray* fClusterRecords;

  /// number of cluster records in array
  Int_t fClusterCount;

ClassDef(AliMUONAlignmentTrackRecord, 1)

};

#endif
