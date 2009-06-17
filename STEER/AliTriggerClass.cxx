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

///////////////////////////////////////////////////////////////////////////////
//
// This class which defines the trigger classes objects
//
//
///////////////////////////////////////////////////////////////////////////////
#include <Riostream.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliTriggerClass.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerDescriptor.h"
#include "AliTriggerCluster.h"
#include "AliTriggerPFProtection.h"
#include "AliTriggerBCMask.h"

ClassImp(AliTriggerClass)

//_____________________________________________________________________________
AliTriggerClass::AliTriggerClass():
  TNamed(),
  fClassMask(0),
  fIndex(0),
  fDescriptor(NULL),
  fCluster(NULL),
  fPFProtection(NULL),
  fMask(NULL),
  fPrescaler(0),
  fAllRare(kFALSE),
  fStatus(kFALSE)
{
  // Default constructor
}

//_____________________________________________________________________________
AliTriggerClass::AliTriggerClass( TString & name, UChar_t index,
				  AliTriggerDescriptor *desc, AliTriggerCluster *clus,
				  AliTriggerPFProtection *pfp, AliTriggerBCMask *mask,
				  UInt_t prescaler, Bool_t allrare) :
  TNamed( name, name ),
  fClassMask( 1ull << ULong64_t(index-1)),
  fIndex(index),
  fDescriptor( desc ),
  fCluster( clus ),
  fPFProtection( pfp ),
  fMask( mask ),
  fPrescaler( prescaler ),
  fAllRare( allrare ),
  fStatus(kFALSE)
{
  // Constructor
}

//_____________________________________________________________________________
AliTriggerClass::AliTriggerClass( AliTriggerConfiguration *config,
				  TString & name, UChar_t index,
				  TString &desc, TString &clus,
				  TString &pfp, TString &mask,
				  UInt_t prescaler, Bool_t allrare) :
  TNamed( name, name ),
  fClassMask( 1ull << ULong64_t(index-1)),
  fIndex(index),
  fDescriptor( NULL ),
  fCluster( NULL ),
  fPFProtection( NULL ),
  fMask( NULL ),
  fPrescaler( prescaler ),
  fAllRare( allrare ),
  fStatus(kFALSE)
{
  fDescriptor = (AliTriggerDescriptor*)config->GetDescriptors().FindObject(desc);
  fCluster = (AliTriggerCluster*)config->GetClusters().FindObject(clus);
  pfp.ReplaceAll("{","");
  pfp.ReplaceAll("}","");
  fPFProtection = (AliTriggerPFProtection*)config->GetPFProtections().FindObject(pfp);
  mask.ReplaceAll("{","");
  mask.ReplaceAll("}","");
  fMask = (AliTriggerBCMask*)config->GetMasks().FindObject(mask);
}

//_____________________________________________________________________________
AliTriggerClass::~AliTriggerClass() 
{ 
  // Destructor
}
//_____________________________________________________________________________
AliTriggerClass::AliTriggerClass( const AliTriggerClass& trclass ):
  TNamed( trclass ),
  fClassMask(trclass.fClassMask),
  fIndex(trclass.fIndex),
  fDescriptor(trclass.fDescriptor),
  fCluster(trclass.fCluster),
  fPFProtection(trclass.fPFProtection),
  fMask(trclass.fMask),
  fPrescaler(trclass.fPrescaler),
  fAllRare(trclass.fAllRare),
  fStatus(trclass.fStatus)
{
   // Copy constructor
}

//______________________________________________________________________________
AliTriggerClass& AliTriggerClass::operator=(const AliTriggerClass& trclass)
{
   // AliTriggerClass assignment operator.

   if (this != &trclass) {
      TNamed::operator=(trclass);
      fClassMask = trclass.fClassMask;
      fIndex=trclass.fIndex;
      fDescriptor = trclass.fDescriptor;
      fCluster = trclass.fCluster;
      fPFProtection = trclass.fPFProtection;
      fMask = trclass.fMask;
      fPrescaler = trclass.fPrescaler;
      fAllRare = trclass.fAllRare;
      fStatus = trclass.fStatus;
   }
   return *this;
}

//_____________________________________________________________________________
Bool_t AliTriggerClass::CheckClass(AliTriggerConfiguration* config) const
{
  // Check the existance of trigger inputs and functions
  // and the logic used.
  // Return false in case of wrong class
  // definition.

  if (!fClassMask) {
    AliError(Form("The class (%s) has invalid mask pattern !",GetName()));
    return kFALSE;
  }

  // check comsistency of index and mask

  if (!config->GetDescriptors().FindObject(fDescriptor)) {
    AliError(Form("The class (%s) contains invalid descriptor !",GetName()));
    return kFALSE;
  }
  else {
    if (!(fDescriptor->CheckInputsAndFunctions(config->GetInputs(),config->GetFunctions()))) {
      AliError(Form("The class (%s) contains bad descriptor !",GetName()));
      return kFALSE;
    }
  }

  if (!config->GetClusters().FindObject(fCluster)) {
    AliError(Form("The class (%s) contains invalid cluster !",GetName()));
    return kFALSE;
  }

  if (!config->GetPFProtections().FindObject(fPFProtection)) {
    AliError(Form("The class (%s) contains invalid past-future protection !",GetName()));
    return kFALSE;
  }

  if (!config->GetMasks().FindObject(fMask)) {
    AliError(Form("The class (%s) contains invalid BC mask !",GetName()));
    return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
void AliTriggerClass::Trigger( const TObjArray& inputs , const TObjArray& functions)
{
   // Check if the inputs satify the trigger class conditions
  fStatus = fDescriptor->Trigger(inputs,functions);
}

//_____________________________________________________________________________
Bool_t AliTriggerClass::IsActive( const TObjArray& inputs, const TObjArray& functions) const
{
   // Check if the inputs satify the trigger class conditions
  if (fDescriptor)
    return fDescriptor->IsActive(inputs,functions);

  return kFALSE;
}

//_____________________________________________________________________________
void AliTriggerClass::Print( const Option_t* ) const
{
   // Print
  cout << "Trigger Class:" << endl;
  cout << "  Name:         " << GetName() << endl;
  cout << "  ClassBit:     0x" << hex << fClassMask << dec << endl;
  cout << "  Index:        " <<  (UInt_t)fIndex <<  endl;
  cout << "  Descriptor:   " << fDescriptor->GetName() << endl;
  cout << "  Cluster:      " << fCluster->GetName() << endl;
  cout << "  PF Protection:" << fPFProtection->GetName() << endl;
  cout << "  BC Mask:      " << fMask->GetName() << endl;
  cout << "  Prescaler:    " << fPrescaler << endl;
  cout << "  AllRare:      " << fAllRare << endl;
  if (fStatus)
     cout << "   Class is fired      " << endl;
   else
     cout << "   Class is not fired  " << endl;
}
