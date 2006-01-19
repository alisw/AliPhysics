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
// This class for running the Central Trigger Processor                      //
//                                                                           //
//                                                                           //
//    Load Descriptors                                                       //
//    Make a list the trigger detectors involve from the descriptors         //
//    For the each event                                                     //
//           Run the Trigger for the each detector                           //
//           Get the inputs                                                  //
//           Check the condition classes                                     //
//           Create the class mask                                           //
//           Save result                                                     //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TStopwatch.h>

#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliModule.h"

#include "AliTriggerInput.h"
#include "AliTriggerDetector.h"
#include "AliTriggerCondition.h"
#include "AliTriggerDescriptor.h"
#include "AliCentralTrigger.h"

ClassImp( AliCentralTrigger )

//_____________________________________________________________________________
AliCentralTrigger::AliCentralTrigger() :
   TObject(),
   fClassMask(0)
{
   // Default constructor
//   LoadDescriptor("Pb-Pb");
}

//_____________________________________________________________________________
AliCentralTrigger::AliCentralTrigger( TString & descriptor ) :
   TObject(),
   fClassMask(0)
{
   // Default constructor
   LoadDescriptor( descriptor );
}

//_____________________________________________________________________________
AliCentralTrigger::~AliCentralTrigger()
{
   // Destructor
   fDescriptors.SetOwner();
   fDescriptors.Delete();
}

//_____________________________________________________________________________
Bool_t AliCentralTrigger::LoadDescriptor( TString & descriptor )
{
   // Load one or more pre-created Descriptors from database/file that match
   // with the input string 'descriptor'
   // Ej: "Pb-Pb" or "p-p-DIMUON CALIBRATION-CENTRAL-BARREL"

   // Delete any descriptor
   fDescriptors.Delete();

   // Load the selected descriptors
   TObjArray* desArray = descriptor.Tokenize( " " );
   Int_t ndes = desArray->GetEntriesFast();
   for( Int_t i=0; i<ndes; i++ ) {
      TObjString* val = (TObjString*)desArray->At( i );
      AliTriggerDescriptor* des = AliTriggerDescriptor::LoadDescriptor( val->String() );
      if( des ) 
         fDescriptors.AddLast( des );
      else
         AliWarning( Form( "Descriptor (%s) not found", val->String().Data() ) );
   }
   Bool_t desfound = kTRUE;
   if( fDescriptors.GetEntriesFast() == 0 ) desfound = kFALSE;

   delete desArray;

   return desfound;
}

//_____________________________________________________________________________
TString AliCentralTrigger::GetDetectors()
{
   // return TString with the detectors to be trigger
   // merging detectors from all descriptors

   TString result;

   Int_t ndes = fDescriptors.GetEntriesFast();
   for( Int_t i=0; i<ndes; i++ ) {
      TString detStr = ((AliTriggerDescriptor*)fDescriptors.At( i ))->GetDetectorCluster();
      TObjArray* det = detStr.Tokenize(" ");
      Int_t ndet = det->GetEntriesFast();
      for( Int_t j=0; j<ndet; j++ ) {
         if( result.Contains( ((TObjString*)det->At(j))->String() ) )continue;
         result.Append( " " );
         result.Append( ((TObjString*)det->At(j))->String() );
      }
   }

   return result;
}

//_____________________________________________________________________________
Bool_t AliCentralTrigger::RunTrigger( AliRunLoader* runLoader )
{
   // run the trigger
   
   if( fDescriptors.GetEntriesFast() == 0 ) {
      AliError("not trigger descriptor loaded, skipping trigger");
      return kFALSE;
   }

   TStopwatch stopwatch;
   stopwatch.Start();
   
   // Process each event
   for( Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++ ) {
      AliInfo( Form("\n ***** Processing event %d *****\n", iEvent) );
      runLoader->GetEvent( iEvent );       
      // Get detectors involve
      TString detStr = GetDetectors();
      TObjArray* detArray = runLoader->GetAliRun()->Detectors();
      // Reset Mask
      fClassMask = 0;
      TObjArray trgdetArray;
      for( Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++ ) {
         AliModule* det = (AliModule*) detArray->At( iDet );
         if( !det || !det->IsActive() ) continue;
         if( IsSelected(det->GetName(), detStr) ) {
            AliInfo( Form("triggering from digits for %s", det->GetName() ) );
            AliTriggerDetector* trgdet = det->CreateTriggerDetector();
            trgdet->CreateInputs();
            trgdetArray.AddLast( trgdet );
            TStopwatch stopwatchDet;
            stopwatchDet.Start();
            trgdet->Trigger();
            AliInfo( Form("Execution time for %s: R:%.2fs C:%.2fs",
                     det->GetName(), stopwatchDet.RealTime(), stopwatchDet.CpuTime() ) );
            // Get the inputs
            TObjArray* detInp = trgdet->GetInputs();
            for( Int_t i=0; i<detInp->GetEntriesFast(); i++ ) {
               fInputs.AddLast( detInp->At(i) );
            }
         }
      }

      // Check trigger conditions and create the trigger class mask
      CheckConditions();

      fInputs.Clear();
      
      // Clear trigger detectors
      trgdetArray.SetOwner();
      trgdetArray.Delete();

      if( (detStr.CompareTo( "ALL" ) != 0) && !detStr.IsNull() ) {
         AliError( Form("the following detectors were not found: %s",
                   detStr.Data()));
         return kFALSE;
      }

      // Write trigger ???? -> Event Header
      // Write();
      Print();

   } // end event loop
   return kTRUE;
}



//_____________________________________________________________________________
Long_t AliCentralTrigger::CheckConditions()
{
   // Check trigger conditions and create the trigger class mask

   Int_t ndes = fDescriptors.GetEntriesFast();
   for( Int_t i=0; i<ndes; i++ ) {
      TObjArray* condArray = ((AliTriggerDescriptor*)fDescriptors.At( i ))->GetTriggerConditions();
      Int_t ncond = condArray->GetEntriesFast();
      for( Int_t j=0; j<ncond; j++ ) {
         AliTriggerCondition* cond = (AliTriggerCondition*)condArray->At( j );
         if( !cond->CheckInputs( fInputs ) ) continue;
         cond->Trigger( fInputs );
    //     cond->Print();
         fClassMask |= cond->GetValue();
      }
   }
   return fClassMask;
}
//_____________________________________________________________________________
TObjArray* AliCentralTrigger::GetResultConditions()
{
   // return only the true conditions

   TObjArray* result = new TObjArray();

   Int_t ndes = fDescriptors.GetEntriesFast();
   for( Int_t i=0; i<ndes; i++ ) {
      TObjArray* condArray = ((AliTriggerDescriptor*)fDescriptors.At( i ))->GetTriggerConditions();
      Int_t ncond = condArray->GetEntriesFast();
      for( Int_t j=0; j<ncond; j++ ) {
         AliTriggerCondition* cond = (AliTriggerCondition*)condArray->At( j );
         if( cond->GetStatus() ) result->AddLast( cond );
      }
   }

   return result;
}


//_____________________________________________________________________________
void AliCentralTrigger::Print( const Option_t*  ) const
{
   // Print
   cout << "Central Trigger: " << endl;
   cout << "  Trigger Class Mask: 0x" << hex << fClassMask << dec << endl;

   Int_t ndes = fDescriptors.GetEntriesFast();
   for( Int_t i=0; i<ndes; i++ ) {
      AliTriggerDescriptor* des = (AliTriggerDescriptor*)fDescriptors.At( i );
      if( des ) des->Print();
   }
   cout << endl;
}


//////////////////////////////////////////////////////////////////////////////
// Helper method

//_____________________________________________________________________________
Bool_t AliCentralTrigger::IsSelected( TString detName, TString& detectors ) const
{
   // check whether detName is contained in detectors
   // if yes, it is removed from detectors

   // check if all detectors are selected
   if( (detectors.CompareTo("ALL") == 0 ) ||
        detectors.BeginsWith("ALL ") ||
        detectors.EndsWith(" ALL") ||
        detectors.Contains(" ALL ") ) {
      detectors = "ALL";
      return kTRUE;
   }

   // search for the given detector
   Bool_t result = kFALSE;
   if( (detectors.CompareTo( detName ) == 0) ||
        detectors.BeginsWith( detName+" " ) ||
        detectors.EndsWith( " "+detName ) ||
        detectors.Contains( " "+detName+" " ) ) {
      detectors.ReplaceAll( detName, "" );
      result = kTRUE;
   }

   // clean up the detectors string
   while( detectors.Contains("  ") )  detectors.ReplaceAll( "  ", " " );
   while( detectors.BeginsWith(" ") ) detectors.Remove( 0, 1 );
   while( detectors.EndsWith(" ") )   detectors.Remove( detectors.Length()-1, 1 );

   return result;
}
