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
// This class for running and define a Trigger Descriptor 
//
// A Trigger Descriptor define a trigger setup for specific runnign
// condition (Pb-Pb, p-p, p-A, Calibration, etc).
// It keep:
//    - cluster detector (List of detectors involved)
//    - List of conditions
//
// Descriptors could be create in advance and store in a file.
//
//   Example how to create a Trigger Descriptor:
//
//   AliTriggerDescriptor descrip( "TEST", "Test Descriptor" );
//
//   // Define a Cluster Detector
//   descrip.AddDetectorCluster( "VZERO ZDC MUON" );
//
//   // Define the trigger conditions (see AliTriggerCondition.cxx)
//   descrip.AddCondition( "VZERO_TEST1_L0 & MUON_SPlus_LPt_L0 & ZDC_TEST2_L0", // condition
//                         "VO1_M1_ZDC2",      // short name
//                         "Dummy",            // short description
//                          0x0100 );          // class mask (set one bit)
//
//   descrip.AddCondition( "VZERO_TEST2_L0 & MUON_SMinus_HPt_L0 & ZDC_TEST1_L0",
//                         "VO2_M3_ZDC1",
//                         "Dummy",
//                          0x0200 );
//
//   descrip.AddCondition( "VZERO_TEST3_L0 | MUON_Unlike_LPt_L0 | ZDC_TEST3_L0",
//                         "VO3_M1_ZDC3",
//                         "Dummy",
//                          0x0400 );
//   descrip.CheckInputsConditions("Config.C");
//   descrip.Print();
//
//   // save the descriptor to file 
//   // (default file name $ALICE_ROOT/data/triggerDescriptor.root)
//   descrip.WriteDescriptor(); or descrip.WriteDescriptor( filename );
//
///////////////////////////////////////////////////////////////////////////////

#include <TString.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TKey.h>
#include <TFile.h>

#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliModule.h"

#include "AliTriggerInput.h"
#include "AliTriggerDetector.h"
#include "AliTriggerCondition.h"
#include "AliTriggerDescriptor.h"


ClassImp(AliTriggerDescriptor)

//_____________________________________________________________________________
const char* AliTriggerDescriptor::fgkDetectorName[AliTriggerDescriptor::fgkNDetectors] =
             { "ITS", "TRD", "PHOS", "EMCAL", "MUON", "ZDC", "T0", "VZERO", "ACORDE", "TOF" };

const TString AliTriggerDescriptor::fgkDescriptorFileName("/data/triggerDescriptors.root");

//_____________________________________________________________________________
AliTriggerDescriptor::AliTriggerDescriptor():
  TNamed(),
  fDetectorCluster(""),
  fConditions()
{
}

//_____________________________________________________________________________
AliTriggerDescriptor::AliTriggerDescriptor( TString & name, TString & description ):
  TNamed( name, description ),
  fDetectorCluster(""),
  fConditions()
{
}

//_____________________________________________________________________________
AliTriggerDescriptor::AliTriggerDescriptor( const AliTriggerDescriptor& des ):
  TNamed( des ),
  fDetectorCluster( des.fDetectorCluster ),
  fConditions()
{
   // Copy constructor
   Int_t ncond = des.fConditions.GetEntriesFast();
   for( Int_t j=0; j<ncond; j++ ) {
      AddCondition( new AliTriggerCondition( *(AliTriggerCondition*)(des.fConditions.At( j )) ) );
   }
}

//______________________________________________________________________________
AliTriggerDescriptor& AliTriggerDescriptor::operator=(const AliTriggerDescriptor& des)
{
   // AliTriggerDescriptor assignment operator.

   if (this != &des) {
      TNamed::operator=(des);
      fDetectorCluster = des.fDetectorCluster;
      fConditions.Delete();
      Int_t ncond = des.fConditions.GetEntriesFast();
      for( Int_t j=0; j<ncond; j++ ) {
         AddCondition( new AliTriggerCondition( *(AliTriggerCondition*)(des.fConditions.At( j )) ) );
      }
   }
   return *this;
}

//_____________________________________________________________________________
Bool_t AliTriggerDescriptor::AddDetectorCluster( TString & cluster )
{
   // Add a List of Detectors to be read together (Detector Cluster)
   // Ej "TO VO ZDC MUON" or "ALL"

   TString olddet = fDetectorCluster;
   TString newdet = cluster;
   for( Int_t iDet = 0; iDet < fgkNDetectors; iDet++ ) {
      if( IsSelected( fgkDetectorName[iDet], newdet ) && !IsSelected( fgkDetectorName[iDet], olddet ) ) {
         // Add the detector
         fDetectorCluster.Append( " " );
         fDetectorCluster.Append( fgkDetectorName[iDet] );
      }
   }

   // check if there are no trigger detectors (E.g. TPC)
   if ((newdet.CompareTo("ALL") != 0) && !newdet.IsNull()) {
      AliError( Form("the following detectors are not trigger detectors: %s",
                newdet.Data() ) );
      return kFALSE;
   }

   return kTRUE;
}

//_____________________________________________________________________________
void AliTriggerDescriptor::AddCondition( TString & cond, TString & name, TString & description, ULong64_t mask   )
{
   // Add a new condition
   AliTriggerCondition* acond = new AliTriggerCondition( cond, name, description, mask );
   fConditions.AddLast( acond );
}


//_____________________________________________________________________________
AliTriggerDescriptor* AliTriggerDescriptor::LoadDescriptor( TString & descriptor, const char* filename )
{
   // Load one pre-created Descriptors from database/file that match
   // with the input string 'descriptor'
   // Ej: "Pb-Pb" or "p-p-DIMUON CALIBRATION-CENTRAL-BARREL"

   // Load the selected descriptor
   TString path;
   if( !filename[0] ) {
      path += gSystem->Getenv("ALICE_ROOT");
      path += fgkDescriptorFileName;
   }
   else
      path += filename;

   if( gSystem->AccessPathName( path.Data() ) ) {
      AliErrorGeneral( "AliTriggerDescriptor", Form( "file (%s) not found", path.Data() ) );
      return NULL;
   }

   TFile file( path.Data(), "READ" );
   AliTriggerDescriptor* des = (AliTriggerDescriptor*)(file.Get( descriptor.Data() ));

   file.Close();

   return des;
}

//_____________________________________________________________________________
TObjArray* AliTriggerDescriptor::GetAvailableDescriptors( const char* filename )
{
   // Return an array of descriptor in the file

   TString path;
   if( !filename[0] ) {
      path += gSystem->Getenv( "ALICE_ROOT" );
      path += fgkDescriptorFileName;
   }
   else
      path += filename;

   if( gSystem->AccessPathName( path.Data() ) ) {
      AliErrorGeneral( "AliTriggerDescriptor", Form( "file (%s) not found", path.Data() ) );
      return NULL;
   }

   TObjArray* desArray = new TObjArray();

   TFile file( path.Data(), "READ" );
   if( file.IsZombie() ) {
      AliErrorGeneral( "AliTriggerDescriptor", Form( "Error opening file (%s)", path.Data() ) );
      return NULL;
   }

   file.ReadAll();

   TKey* key;
   TIter next( file.GetListOfKeys() );
   while( (key = (TKey*)next()) ) {
      TObject* obj = key->ReadObj();
      if( obj->InheritsFrom( "AliTriggerDescriptor" ) ) {
         desArray->AddLast( obj );
      }
   }
   file.Close();

   return desArray;
}

//_____________________________________________________________________________
void AliTriggerDescriptor::WriteDescriptor( const char* filename )
{
   // Load one pre-created Descriptors from database/file that match
   // with the input string 'descriptor'
   // Ej: "Pb-Pb" or "p-p-DIMUON CALIBRATION-CENTRAL-BARREL"

   // Load the selected descriptor
   TString path;
   if( !filename[0] ) {
      path += gSystem->Getenv("ALICE_ROOT");
      path += fgkDescriptorFileName;
   }
   else
      path += filename;

   TFile file( path.Data(), "UPDATE" );
   if( file.IsZombie() ) {
      AliErrorGeneral( "AliTriggerDescriptor", 
                        Form( "Can't open file (%s)", path.Data() ) );
      return;
   }

   Bool_t result = (Write( GetName(), TObject::kOverwrite ) != 0);
   if( !result )
      AliErrorGeneral( "AliTriggerDescriptor",
                        Form( "Can't write entry to file <%s>!", path.Data() ) );
   file.Close();
}

//_____________________________________________________________________________
Bool_t AliTriggerDescriptor::CheckInputsConditions( TString& configfile )
{
   // To be used on the pre-creation of Descriptors to check if the
   // conditions have valid inputs names.
   //
   // Initiate detectors modules from a Config file
   // Ask to each active module present in the fDetectorCluster
   // to create a Trigger detector and retrive the inputs from it
   // to create a list of inputs.
   // Each condition in the descriptor is then checked agains 
   // the list of inputs


   if (!gAlice) {
      AliError( "no gAlice object. Restart aliroot and try again." );
      return kFALSE;
   }
   if (gAlice->Modules()->GetEntries() > 0) {
      AliError( "gAlice was already run. Restart aliroot and try again." );
      return kFALSE;
   }

   AliInfo( Form( "initializing gAlice with config file %s",
            configfile.Data() ) );
   StdoutToAliInfo( StderrToAliError(
      gAlice->Init( configfile.Data() );
   ););

   AliRunLoader* runLoader = gAlice->GetRunLoader();
   if( !runLoader ) {
      AliError( Form( "gAlice has no run loader object. "
                      "Check your config file: %s", configfile.Data() ) );
      return kFALSE;
   }

   // get the possible inputs to check the condition
   TObjArray inputs;
   TObjArray* detArray = runLoader->GetAliRun()->Detectors();

   TString detStr = fDetectorCluster;
   for( Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++ ) {
      AliModule* det = (AliModule*) detArray->At(iDet);
      if( !det || !det->IsActive() ) continue;
      if( IsSelected( det->GetName(), detStr ) ) {
         AliInfo( Form( "Creating inputs for %s", det->GetName() ) );
         AliTriggerDetector* dtrg = det->CreateTriggerDetector();
         dtrg->CreateInputs();
         TObjArray* detInp = dtrg->GetInputs();
         for( Int_t i=0; i<detInp->GetEntriesFast(); i++ ) {
            AliInfo( Form( "Adding input %s", ((AliTriggerInput*)detInp->At(i))->GetName() ) );
            inputs.AddLast( detInp->At(i) );
         }
      }
   }

   // check if the condition is compatible with the triggers inputs
   Int_t ncond = fConditions.GetEntriesFast();
   Bool_t check = kTRUE;
   ULong64_t mask = 0L;
   for( Int_t j=0; j<ncond; j++ ) {
      AliTriggerCondition* cond = (AliTriggerCondition*)(fConditions.At( j ));
      if( !(cond->CheckInputs( inputs )) ) check = kFALSE;
      else AliInfo( Form( "Condition (%s) inputs names OK, class mask (0x%Lx)",
                    cond->GetName(), cond->GetMask( ) ) );
      // check if condition mask is duplicated
      if( mask & cond->GetMask() ) {
         AliError( Form("Condition (%s). The class mask (0x%Lx) is ambiguous. It was previous defined",
                   cond->GetName(), cond->GetMask()  ) );
         check = kFALSE;
      }
      mask |= cond->GetMask();
   }

   return check;
}


//_____________________________________________________________________________
void AliTriggerDescriptor::Print( const Option_t*  ) const
{
   // Print
   cout << "Trigger Descriptor:"  << endl;
   cout << "  Name:             " << GetName() << endl; 
   cout << "  Description:      " << GetTitle() << endl;
   cout << "  Detector Cluster: " << fDetectorCluster << endl;

//   Int_t ninputs = fInputs->GetEntriesFast();
//   for( Int_t i=0; i<ninputs; i++ ) {
//      AliTriggerInput* in = (AliTriggerInput*)fInputs->At(i)
//      in->Print();
//   }

   Int_t ncond = fConditions.GetEntriesFast();
   for( Int_t i=0; i<ncond; i++ ) {
      AliTriggerCondition* in = (AliTriggerCondition*)fConditions.At(i);
      in->Print();
   }
   cout << endl;
}


//////////////////////////////////////////////////////////////////////////////
// Helper method

//_____________________________________________________________________________
Bool_t AliTriggerDescriptor::IsSelected( TString detName, TString& detectors ) const
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
