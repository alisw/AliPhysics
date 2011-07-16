#ifndef ALITRIGGERDETECTOR_H
#define ALITRIGGERDETECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Base Class for Detector specific Trigger                                 //                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include <TNamed.h>
class TString;
class AliTriggerInput;


class AliTriggerDetector : public TNamed {

public:
                          AliTriggerDetector();
               virtual   ~AliTriggerDetector();
	       AliTriggerDetector(const AliTriggerDetector & de );

          virtual void    AssignInputs(const TObjArray& inputs);
          virtual void    CreateInputs();
          virtual void    Trigger();
  //  Setters
                  void    AddInput( TObject * input ) { fInputs.AddLast( input ); }
                  void    SetInput( TString & name );
                  void    SetInput( const char * name );
  //  Getters
             TObjArray*   GetInputs() { return &fInputs; }
                Long_t    GetMask() const { return fMask; }

       AliTriggerInput*   GetInput( TString & name ) {
                             return ((AliTriggerInput*)fInputs.FindObject( name.Data() ));
                          }
       AliTriggerInput*   GetInput( const char *  name ) {
                             return ((AliTriggerInput*)fInputs.FindObject( name ));
                          }
          virtual void    Print( const Option_t* opt ="" ) const;

protected:
                Long_t    fMask;      // Trigger Mask ( bitwise OR of trigger inputs )
             TObjArray    fInputs;    // Array of Triggers Inputs (AliTriggerInput class)

private:
	     AliTriggerDetector&   operator=(const AliTriggerDetector& de);

   ClassDef( AliTriggerDetector, 1 )  // Base Class for Detector specific Trigger
};

#endif
