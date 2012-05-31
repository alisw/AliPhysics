#ifndef ALIVMFT_H
#define ALIVMFT_H

//-------------------------------------------------------------------------
//     Base class for ESD and AOD MFT data
//     Author: Cvetan Cheshkov
//     cvetan.cheshkov@cern.ch 2/02/2011
//-------------------------------------------------------------------------

#include "TObject.h"

class AliVMFT : public TObject 
{
public:
  AliVMFT() { }
  AliVMFT(const AliVMFT& source);
  AliVMFT &operator=(const AliVMFT& source);

  virtual ~AliVMFT() { }


protected:  

    
  ClassDef(AliVMFT,1)
};

#endif
