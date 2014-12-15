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


#include "AliCaloBunchInfo.h"

// Container class to hold information 
// about  of bunches.
// Used by the
// AliCaloRawAnalyzer
// classes
AliCaloBunchInfo::AliCaloBunchInfo( UInt_t start, Int_t length, const UShort_t * data ) :  fStartTimebin(start),
											     fLength(length),
											     fkData(data)
{


}



AliCaloBunchInfo::~AliCaloBunchInfo()
{

  
}



AliCaloBunchInfo::AliCaloBunchInfo( const AliCaloBunchInfo  & rhs) :fStartTimebin( rhs.fStartTimebin ),
								       fLength(  rhs.fLength ),
								       fkData( rhs.fkData )
{
  
  
}




AliCaloBunchInfo&  AliCaloBunchInfo::operator = ( const  AliCaloBunchInfo & rhs)
{
  //This is just to get of compliation warning. Its not really needed
  if(this != & rhs) 
    {
      fStartTimebin  = rhs.fStartTimebin;
      fLength = rhs.fLength;
      fkData = rhs.fkData;	
    }
  return *this;
}

