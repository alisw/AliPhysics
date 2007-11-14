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

class TNamed;


class AliTriggerBCMask : public TNamed {

public:
                          AliTriggerBCMask();
                          AliTriggerBCMask( TString & name, UChar_t *mask = NULL );
                          AliTriggerBCMask( const AliTriggerBCMask& mask );
               virtual   ~AliTriggerBCMask();
  AliTriggerBCMask&   operator=(const AliTriggerBCMask& mask);

                 UChar_t* GetFullMask () const {return &fBCMask; }
		  Bool_t  GetMask(UShort_t index) const;
  
  enum {kNBytesPerBCMask = 446}; // Number of bytes to store the 3565 bits of BC mask

private:
                UChar_t    fBCMask[kNBytesPerBCMask];         // Bunch cross mask (3565 bit)

   ClassDef( AliTriggerBCMask, 1 )  // Define a trigger bunch-crossing mask
};

#endif
