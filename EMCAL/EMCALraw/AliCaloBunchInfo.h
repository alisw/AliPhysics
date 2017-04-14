#ifndef ALICALOBUNCHINFO_H
#define ALICALOBUNCHINFO_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include "Rtypes.h"

//_________________________________________________________________________
/// \class AliCaloBunchInfo
/// \ingroup EMCALraw
/// \brief   Container of ALTRO information 
///
/// Container class to hold information  about ALTRO
/// bunches from the altro stream.
/// Used by the AliCaloRawAnalyzer classes
///
/// Each bunch has a start marker, ( fStartTimebin ) 
/// the number of ADC samples in the bunch fLength, and a pointer
/// to the last (fStartTimebin + fLength ) time bin of the bunch.
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
//_________________________________________________________________________      

class  AliCaloBunchInfo
{
 
public:

  AliCaloBunchInfo( UInt_t starttimebin, Int_t length,  const UShort_t * data );
  virtual ~AliCaloBunchInfo();

  AliCaloBunchInfo                ( const AliCaloBunchInfo  & rhs);
  AliCaloBunchInfo   & operator = ( const  AliCaloBunchInfo & rhs);

  UInt_t          GetStartBin() const { return fStartTimebin ; }
  Int_t           GetLength  () const { return fLength       ; }
  const UShort_t *GetData    () const { return fkData        ; }
  
 private:
  
  AliCaloBunchInfo();
  
  UInt_t fStartTimebin;   ///< Starttimebin as given by the ALTRO stream
  
  Int_t fLength;          ///< Length of the bunch
  
  const UShort_t *fkData; ///< Pointer to the last data enetry of the bunch (data is reversed with respect to fStartTimebin)
  
};

#endif //ALICALOBUNCHINFO_H
