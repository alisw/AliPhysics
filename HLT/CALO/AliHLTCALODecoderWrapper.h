#ifndef ALIHLTCALODECODERWRAPPER_H
#define ALIHLTCALODECODERWRAPPER_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
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

/*
class AliAltroRawStreamV3;
class AliCaloRawStreamV3;
class AliRawReaderMemory;
*/
 
#include "AliAltroRawStreamV3.h"
#include "AliCaloRawStreamV3.h"
#include "AliRawReaderMemory.h"

class AliHLTComponentBlockData;


class  AliHLTCALODecoderWrapper
{
 public:
  AliHLTCALODecoderWrapper();
  virtual ~AliHLTCALODecoderWrapper();
  void SetMemory( AliHLTComponentBlockData *dtaptr );
  inline bool  NextChannel          ( )       { return  fAltroRawStream->NextChannel();  };
  inline bool NextBunch             ( )       { return  fAltroRawStream->NextBunch();    };
  inline const UShort_t *GetSignals ( )       { return  fAltroRawStream->GetSignals();   };
  inline Int_t  GetHWAddress        ( ) const { return  fAltroRawStream->GetHWAddress();}; 
  inline Int_t  GetBunchLength      ( ) const { return  fAltroRawStream->GetBunchLength();  };
  inline UInt_t GetStartTimeBin     ( ) const { return  fAltroRawStream->GetEndTimeBin(); };
  inline UInt_t GetEndTimeBin       ( ) const { return  fAltroRawStream->GetStartTimeBin(); };

 private:
  AliAltroRawStreamV3 *fAltroRawStream;
  //  AliCaloRawStreamV3  *fCaloRawStream;
  AliRawReaderMemory  *fReaderMemory;
};

#endif
