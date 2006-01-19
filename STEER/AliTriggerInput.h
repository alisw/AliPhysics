#ifndef ALITRIGGERINPUT_H
#define ALITRIGGERINPUT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//
//  Class to define a Trigger Input from an specific detector                                                                                           //
//
//
//                        name         description     id mask
//    Ej:
//      AliTriggerInput( "V0_MB_L0", "VO minimum bias", 0x01 );
//      AliTriggerInput( "V0_SC_L0", "VO semi central", 0x02 );
//      AliTriggerInput( "V0_C_L0",  "VO central",      0x04 );

//    The name must be globaly unique. Spaces are not allowed.
//    As convention should start with detector name then an id
//    and the trigger level (L0, L1, L2)
//
//    A maximun of 60 inputs trigger are allow.
//    So, the id mask should set only bit from the position 1 to 60.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class AliTriggerInput : public TNamed {

public:
                          AliTriggerInput() {}
                          AliTriggerInput( TString name, TString description, Long_t mask )
                                           : TNamed( name.Data(), description.Data() ),
                                             fMask( mask ),
                                             fValue( 0 ) {}
               virtual   ~AliTriggerInput() {}

  //  Setters
          virtual void    Set()   { fValue = fMask; }
          virtual void    Reset() { fValue = 0; }

  //  Getters
                Bool_t    Status() const   { return (Bool_t)fValue; }
                Long_t    GetValue() const { return fValue; }
                 Int_t    GetMask() const  { return fMask; }
  //             ULong_t    Hash() const { return TMath::Hash( GetName().Data() ); };

           virtual void    Print( const Option_t* opt ="" ) const;

protected:
                Long_t    fMask;        //  Trigger ID mask (1 bit)
                Long_t    fValue;       //  Trigger Signal (0 = false, > 1 = true = fMask )
     //          Int_t      fLevel;       //  Trigger Level (L0, L1, L2)

   ClassDef( AliTriggerInput, 1 )  // Define a Trigger Input
};


#endif
