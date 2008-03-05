#ifndef ALITRIGGERBCMASK_H
#define ALITRIGGERBCMASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This class represents the CTP bunch-crossing mask                         //
//                                                                           //
// The Mask contains name and 3565 bits for each bunch-crossing in an orbit  //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class AliTriggerBCMask : public TNamed {

public:
                          AliTriggerBCMask();
                          AliTriggerBCMask( TString & name);
                          AliTriggerBCMask( TString & name, TString & mask );
                          AliTriggerBCMask( const AliTriggerBCMask& mask );
               virtual   ~AliTriggerBCMask();
  AliTriggerBCMask&   operator=(const AliTriggerBCMask& mask);

           const UChar_t* GetFullMask () const {return fBCMask; }
		  Bool_t  GetMask(UShort_t index) const;
		    void  Print( const Option_t* ) const;
  
  enum {kNBytesPerBCMask = 446}; // Number of bytes to store the 3565 bits of BC mask

private:
                void   CreateMask(TString &/*mask*/) {} 

                UChar_t   fBCMask[kNBytesPerBCMask];         // Bunch cross mask (3565 bit)

   ClassDef( AliTriggerBCMask, 1 )  // Define a trigger bunch-crossing mask
};

#endif
