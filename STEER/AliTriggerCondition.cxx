/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *

 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//
//  Class to define a Trigger Condition                                                                                                                 //
//
//  Ej                       Condition                        name   Description   class mask
//         inputs names  ___ _________ _________
//                          |         |         |
//   AliTriggerCondition("(T0_L0 & VZERO_MB_L0 & TRD_PRE_L0)", "MB", "Minimum Bias", 0x0100 );
//
//  A Trigger condition is defined from logical combination of trigger
//  inputs names (boolean expression), trigger inputs names must match
//  with the inputs defined in AliTriggerDetector classes
//
//      Allow operators:
//                &    =>  and
//                |    =>  or
//                !    =>  not
//
//    The name must be globally unique. Spaces are not allowed.
//
//    A maximun of 50 diffentes trigger signatures ("trigger classes" or conditions)
//    are allow to run simultaneously.
//    So, the "class mask" should set only 1 bit from the position 1 to 50.
//
//
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>

#include "AliLog.h"
#include "AliExpression.h"
#include "AliTriggerInput.h"
#include "AliTriggerCondition.h"

ClassImp( AliTriggerCondition )

//_____________________________________________________________________________
AliTriggerCondition::AliTriggerCondition() :
   TNamed(),
   fClassMask( 0 ),
   fCondition( "" ),
   fStatus( kFALSE )
{
   // Default ctor
}

//______________________________________________________________________________
AliTriggerCondition::AliTriggerCondition(const AliTriggerCondition &cond) :
   TNamed( cond ),
   fClassMask( cond.fClassMask ),
   fCondition( cond.fCondition ),
   fStatus( cond.fStatus )
{
   // AliTriggerCondition copy ctor.
}

//______________________________________________________________________________
AliTriggerCondition& AliTriggerCondition::operator=(const AliTriggerCondition& rhs)
{
   // AliTriggerCondition assignment operator.

   if (this != &rhs) {
      TNamed::operator=(rhs);
      fClassMask  = rhs.fClassMask;
      fCondition = rhs.fCondition;
   }
   return *this;
}

//_____________________________________________________________________________
AliTriggerCondition::AliTriggerCondition( TString & condition, TString & name,
                                          TString & description, ULong64_t mask ) :
   TNamed( name, description ),
   fClassMask( mask ),
   fCondition( condition ),
   fStatus( kFALSE )
{
   // Default Constructor

   // check the expression 
   AliExpression* exp = new AliExpression( fCondition );
   delete exp;
}


//_____________________________________________________________________________
Bool_t AliTriggerCondition::CheckInputs( TObjArray& inputs )
{
   // The "inputs" array should be the list of all possible inputs
   // so each input in the condition is checked to be present in the array
   // return false if one input is missing

   TString condition( fCondition );
   TObjArray* tokens = condition.Tokenize(" !&|()");

   Int_t ntokens = tokens->GetEntriesFast();
   for( Int_t i=0; i<ntokens; i++ ) {
      TObjString* iname = (TObjString*)tokens->At( i );
      Int_t nInputs = inputs.GetEntriesFast();
      Int_t j;
      for( j=0; j<nInputs; j++ ) {
         AliTriggerInput* in = (AliTriggerInput*)inputs.At( j );
         if( (iname->String()).CompareTo( in->GetName() ) == 0 ) break;
      }

      if( j >= nInputs ) {
         AliWarning( Form( "The trigger input (%s) is not available for Condition (%s)",
                      iname->String().Data(), GetName() ) );
         delete tokens;
         return kFALSE;
      }
   }

   delete tokens;
   return kTRUE;
}

//_____________________________________________________________________________
void AliTriggerCondition::Trigger( TObjArray& inputs )
{
   // Check if the inputs satify the expression condition 
   AliExpression* exp = new AliExpression( fCondition );
   fStatus = exp->Value( inputs );
   delete exp;
}

//_____________________________________________________________________________
void AliTriggerCondition::Print( const Option_t* ) const
{
   // Print
   cout << "Trigger Condition:" << endl;
   cout << "  Name:        " << GetName() << endl;
   cout << "  Description: " << GetTitle() << endl;
   cout << "  Condition:   " << fCondition << endl;
   cout << "  Class Mask:  " << "0x" << hex << fClassMask << endl;
   cout << "  Value:       " << "0x" << hex << GetValue() << dec << endl;
}
