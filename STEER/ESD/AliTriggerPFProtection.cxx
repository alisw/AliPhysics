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
//  Class to define a Trigger Past Future Protection
//
//                              name      description         INT1  INT2
//    Ej:
//  AliTriggerPFProtection sc( "BARREL", "BARREL DET Pb-Pb",  "SC","PE" );
//  sc.SetINTa("PE");  // Peripheral
//  sc.SetINTb("SC");  // Semicentral
//  sc.SetINT("PE");
//  sc.SetNa1( 5 );
//  sc.SetNa2( 5 );
//  sc.SetTa( 88 );
//  sc.SetNb1( 1 );
//  sc.SetNb2( 0 );
//  sc.SetTb( 88 );
//
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include <TObject.h>
#include <TString.h>
#include <TObjString.h>

#include "AliLog.h"
#include "AliTriggerPFProtection.h"

ClassImp( AliTriggerPFProtection )
//_____________________________________________________________________________
AliTriggerPFProtection::AliTriggerPFProtection() :
TNamed(), 
fINTa(), fINTb(),fINT(),
fNa1(0),fNa2(0),fTa(0),
fNb1(0),fNb2(0),fTb(0) 
{
for(Int_t i=0;i<12;i++)pfdef[i]=0;
}
AliTriggerPFProtection::AliTriggerPFProtection( TString & name) :
TNamed(name, name), 
fINTa(), fINTb(),fINT(),
fNa1(0),fNa2(0),fTa(0),
fNb1(0),fNb2(0),fTb(0)
{
for(Int_t i=0;i<12;i++)pfdef[i]=0;
}
AliTriggerPFProtection::AliTriggerPFProtection( TString & name,TString & inta, TString & intb, TString & interaction ) :
TNamed(name, name), 
fINTa(inta), fINTb(intb),fINT(interaction),
fNa1(0),fNa2(0),fTa(0),
fNb1(0),fNb2(0),fTb(0)
{
for(Int_t i=0;i<12;i++)pfdef[i]=0;
}
AliTriggerPFProtection::AliTriggerPFProtection(TString& name,UInt_t* pfdefin)
:TNamed(name,name),
fINTa(), fINTb(),fINT(),
fNa1(0),fNa2(0),fTa(0),
fNb1(0),fNb2(0),fTb(0)
{
for(Int_t i=0;i<12;i++)pfdef[i]=pfdefin[i];
}
//_____________________________________________________________________________
void AliTriggerPFProtection::Print( const Option_t* ) const
{
   // Print
  cout << "Trigger Past-Future Protection: " << endl;
  cout << "  Name:                          " << GetName() << endl;
  cout << "  Interaction_a:                 " << fINTa.Data() << endl;
  cout << "  Interaction_b:                 " << fINTb.Data() << endl;
  cout << "  Interaction:                   " << fINT.Data() << endl;
  cout << "  Na1: " << fNa1 << " Na2: " << fNa2 << " Ta: " << fTa << endl;
  cout << "  Nb1: " << fNb1 << " Nb2: " << fNb2 << " Tb: " << fTb << endl;
  if(1){
    cout << "pfdef: " << hex;
    for(Int_t i=0;i<12;i++)cout << " 0x" << pfdef[i];
    cout << dec << endl;
  }
}

//_____________________________________________________________________________
Bool_t AliTriggerPFProtection::CheckInteractions(TObjArray &interactions) const
{
  // Check if the interactions are valid
  {
    TString logic( GetINTa() );
    TObjArray* tokens = logic.Tokenize(" !&|()\t");

    Int_t ntokens = tokens->GetEntriesFast();
    for( Int_t i=0; i<ntokens; i++ ) {
      TObjString* iname = (TObjString*)tokens->At( i );

      if (!interactions.FindObject(iname->String().Data())) {
	AliError( Form( "The trigger interaction (%s) is not available for past-future protection (%s)",
			iname->String().Data(), GetName() ) );
	delete tokens;
	return kFALSE;
      }
    }
    delete tokens;
  }
  {
    TString logic( GetINTb() );
    TObjArray* tokens = logic.Tokenize(" !&|()\t");

    Int_t ntokens = tokens->GetEntriesFast();
    for( Int_t i=0; i<ntokens; i++ ) {
      TObjString* iname = (TObjString*)tokens->At( i );

      if (!interactions.FindObject(iname->String().Data())) {
	AliError( Form( "The trigger interaction (%s) is not available for past-future protection (%s)",
			iname->String().Data(), GetName() ) );
	delete tokens;
	return kFALSE;
      }
    }
    delete tokens;
  }
  {
    TString logic( GetINT() );
    TObjArray* tokens = logic.Tokenize(" !&|()\t");

    Int_t ntokens = tokens->GetEntriesFast();
    for( Int_t i=0; i<ntokens; i++ ) {
      TObjString* iname = (TObjString*)tokens->At( i );

      if (!interactions.FindObject(iname->String().Data())) {
	AliError( Form( "The trigger interaction (%s) is not available for past-future protection (%s)",
			iname->String().Data(), GetName() ) );
	delete tokens;
	return kFALSE;
      }
    }
    delete tokens;
  }
  return kTRUE;
}
