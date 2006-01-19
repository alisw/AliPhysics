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
//                                                                           //
//  Base Class for Detector specific Trigger                                 //                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TString.h>
#include <TObjArray.h>
#include <TRandom.h>


#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"

#include "AliTriggerInput.h"
#include "AliTriggerDetector.h"



ClassImp( AliTriggerDetector )

//_____________________________________________________________________________
AliTriggerDetector::AliTriggerDetector() : TNamed()
{
   // Default constructor
   fMask    = 0;
}

//_____________________________________________________________________________
void AliTriggerDetector::CreateInputs()
{
   // Define the inputs to the Central Trigger Processor
   // This is a dummy version 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   TString name = GetName();
   fInputs.AddLast( new AliTriggerInput( name+"_TEST1_L0", "Dummy input 1", 0x01 ) );
   fInputs.AddLast( new AliTriggerInput( name+"_TEST2_L0", "Dummy input 2", 0x02 ) );
   fInputs.AddLast( new AliTriggerInput( name+"_TEST3_L0", "Dummy input 3", 0x04 ) );
}

//_____________________________________________________________________________
void AliTriggerDetector::Trigger()
{
   // This is a dummy version set all inputs in a random way

   AliWarning( Form( "Triggering dummy detector %s", GetName() ) );

   //  ********** Get Digits for the current event **********
   AliRunLoader* runLoader = gAlice->GetRunLoader();
   AliInfo( Form( "Event %d", runLoader->GetEventNumber() ) );
   
   TString loadername = GetName(); 
   loadername.Append( "Loader" );
   AliLoader * loader = runLoader->GetLoader( loadername.Data() );
   if( loader ) {
      loader->LoadDigits( "READ" );
      TTree* digits = loader->TreeD();
      // Do something with the digits !!!
      if( digits ) { 
//         digits->Print();
      }   
      loader->UnloadDigits();
   }
   //  ******************************************************

   // set all inputs in a random way
   Int_t nInputs = fInputs.GetEntriesFast();
   for( Int_t j=0; j<nInputs; j++ ) {
      AliTriggerInput* in = (AliTriggerInput*)fInputs.At( j );
      if( gRandom->Rndm() > 0.5 ) {
         in->Set();
         fMask |= in->GetValue();
      }
   }
}

//_____________________________________________________________________________
void AliTriggerDetector::SetInput( TString& name )
{
   // Set Input by name
   AliTriggerInput* in = ((AliTriggerInput*)fInputs.FindObject( name.Data() ));
   in->Set();
   fMask |= in->GetValue();
}

//_____________________________________________________________________________
void AliTriggerDetector::SetInput( const char * name )
{
   // Set Input by name
   AliTriggerInput* in = ((AliTriggerInput*)fInputs.FindObject( name ));
   in->Set();
   fMask |= in->GetValue();
}

//_____________________________________________________________________________
void AliTriggerDetector::SetInput( Int_t mask )
{
   // Set Input by mask
   Int_t nInputs = fInputs.GetEntriesFast();
   for( Int_t j=0; j<nInputs; j++ ) {
      AliTriggerInput* in = (AliTriggerInput*)fInputs.At( j );
      if( in->GetMask() == mask ) { 
         in->Set();
         fMask |= in->GetValue();
         break;
      }
   }
}

