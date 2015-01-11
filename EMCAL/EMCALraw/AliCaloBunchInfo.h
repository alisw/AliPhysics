#ifndef ALICALOBUNCHINFO_H
#define ALICALOBUNCHINFO_H

/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <p.t.hille@fys.uio.no>                *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to p.t.hille@fys.uio.no                             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "Rtypes.h"

// Container class to hold 
// information  about ALTRO
// Bunces from the altro stream.
// Each bunch has a start marker, ( fStartTimebin ) 
// the number of ADC samples in the bunch fLength, and a pointer
// to the last (fStartTimebin + fLength ) time bin of the bunch.
// 
class  AliCaloBunchInfo
{
 public:
  AliCaloBunchInfo( UInt_t starttimebin, Int_t length,  const UShort_t * data );
  virtual ~AliCaloBunchInfo();

  AliCaloBunchInfo( const AliCaloBunchInfo  & rhs);
  AliCaloBunchInfo   & operator = ( const  AliCaloBunchInfo & rhs);

  
  UInt_t GetStartBin( ) const { return fStartTimebin;};
  Int_t GetLength() const { return fLength; };
  const UShort_t *GetData() const { return fkData; };
  
 private:
  AliCaloBunchInfo();
  UInt_t fStartTimebin;   //Starttimebin as given by the ALTRO stream
  Int_t fLength;          //Length of the bunch
  const UShort_t *fkData; //Pointer to the last data enetry of the bunch (data is reversed with respect to fStartTimebin)
};



#endif
