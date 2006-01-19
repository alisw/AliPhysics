#ifndef ALITRIGGERCONDITION_H
#define ALITRIGGERCONDITION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//
//  Class to define a Trigger Condition  
//                                                                                                              //
//  A Trigger condition is defined from logical combination of trigger
//  inputs names (boolean expression)
//
///////////////////////////////////////////////////////////////////////////////

class TNamed;
class TObjArray;
class TString;

class AliTriggerCondition : public TNamed {

public:
                          AliTriggerCondition();
                          AliTriggerCondition( const AliTriggerCondition &cond );
                          AliTriggerCondition( TString & condition, TString & name,
                                               TString & description, Long_t mask );
               virtual   ~AliTriggerCondition() {}
   AliTriggerCondition&   operator=(const AliTriggerCondition& rhs);

                  void    Trigger( TObjArray & inputs );
                Bool_t    CheckInputs( TObjArray & inputs );
  //  Setters
                  void    Reset() { fStatus = kFALSE; }
  //  Getters
                Long_t    GetValue() const { return (fStatus) ? fClassMask : 0; }
                Bool_t    GetStatus() const { return fStatus; }
          virtual void    Print( const Option_t* opt ="" ) const;
protected:
                Long_t    fClassMask;   // UID "class mask" should set only 1 bit from the position 0 to 50
               TString    fCondition;   // Definition of the condition
                Bool_t    fStatus;      // true = Condition has been satisfied after Trigger

   ClassDef( AliTriggerCondition, 1 )  // Define a Trigger Condition
};

#endif
