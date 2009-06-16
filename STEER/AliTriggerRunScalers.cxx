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

/* $Id: AliTriggerRunScalers.cxx 22322 2007-11-22 11:43:14Z cvetan $ */

///////////////////////////////////////////////////////////////////////////////
//
// Class to define the ALICE Trigger Scalers Record per Run.
//
//
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <Riostream.h>
#include <TObject.h>
#include <TClass.h>
#include <TArrayC.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TFile.h>

#include "AliLog.h"  
#include "AliTimeStamp.h"
#include "AliTriggerScalers.h"
#include "AliTriggerScalersRecord.h"
#include "AliTriggerRunScalers.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"

ClassImp( AliTriggerRunScalers )

//_____________________________________________________________________________
AliTriggerRunScalers::AliTriggerRunScalers():
  TObject(),
  fVersion(0),
  fRunNumber(0),
  fnClasses(0),
  fClassIndex(),                    
  fScalersRecord()
{
  // Default constructor
}

//_____________________________________________________________________________
void AliTriggerRunScalers::AddTriggerScalers( AliTriggerScalersRecord* scaler ) 
{ 
  //
  fScalersRecord.AddLast( scaler );
  if (!AliTriggerRunScalers::ConsistencyCheck()) AliErrorClass("Trigger counters not in the right order or decreasing!");
//  fScalersRecord.Sort(); 
}

//_____________________________________________________________________________
AliTriggerRunScalers* AliTriggerRunScalers::ReadScalers( TString & filename )
{
  
  if( gSystem->AccessPathName( filename.Data() ) ) {
     AliErrorClass( Form( "file (%s) not found", filename.Data() ) );
     return NULL;
  }

  ifstream *file = new ifstream ( filename.Data() );
  if (!*file) {
    AliErrorClass(Form("Error opening file (%s) !\n",filename.Data()));
    file->close();
    delete file;
    return NULL;
  }
  
  AliTriggerRunScalers* rScaler = new AliTriggerRunScalers();
  
  TString strLine;
  Bool_t verflag = kFALSE;
  Bool_t classflag = kFALSE;
  UChar_t nclass = 0;
  while (strLine.ReadLine(*file)) {
    if (strLine.BeginsWith("#")) continue;
    
    TObjArray *tokens = strLine.Tokenize(" \t");
    Int_t ntokens = tokens->GetEntriesFast();
    // 1st line, version, it is one number, 
    if (!verflag) {
      if (ntokens != 1) { 
        AliErrorClass( Form( "Error reading version number from (%s), line :%s\n", 
                              filename.Data() , strLine.Data() ) );  
        return NULL;
      }
  //    cout << "Version "<< ((TObjString*)tokens->At(0))->String().Atoi() << endl;
      rScaler->SetVersion( ((TObjString*)tokens->At(0))->String().Atoi() );
      verflag = kTRUE;
      delete tokens;
      continue;
    }
   
    // 2nd line, run number , number of classes, list of classes used in this partition

    if (!classflag) {
      if ( !((TObjString*)tokens->At(1))->String().IsDigit() ) {
        AliErrorClass( Form( "Error reading Run number from (%s)\n", filename.Data() )); 
      }
  //    cout << "Run Number " << ((TObjString*)tokens->At(0))->String().Atoi() << endl;
      rScaler->SetRunNumber( ((TObjString*)tokens->At(0))->String().Atoi() );
      nclass = (UChar_t)((TObjString*)tokens->At(1))->String().Atoi();
  //    cout << "Number of classes " << nclass << endl;
      rScaler->SetNumClasses( nclass );
      if ( nclass != ntokens - 2 ) {
        AliErrorClass( Form( "Error reading number of classes from (%s)\n", filename.Data() )); 
      }
      for (UChar_t i=0; i<nclass; ++i) {
        rScaler->SetClass( i, (Char_t)((TObjString*)tokens->At(2+i))->String().Atoi() );
      }
      classflag = kTRUE;
      delete tokens;
      continue;
    }
    
    // Records
    // Each record consists of 1+(number of classes in partition) lines
    //  1st line of record = time stamp (4 words)
    //                        1st word = ORBIT(24 bits)
    //                        2nd word = Period Counter (28 bit)
    //                        3rd word = seconds (32 bits) from epoch
    //                        4th word = microsecs (32 bits)
    //  other lines = 6 words of counters (L0 before,L0 after, ....)
    if (ntokens != 4) { 
      AliErrorClass( Form( "Error reading timestamp from (%s): line (%s)", 
                            filename.Data(), strLine.Data() )); 
      return NULL;
    }

    UInt_t orbit     = strtoul(((TObjString*)tokens->At(0))->String(), NULL, 10);
    UInt_t period    = strtoul(((TObjString*)tokens->At(1))->String(), NULL, 10);
    UInt_t seconds   = strtoul(((TObjString*)tokens->At(2))->String(), NULL, 10);
    UInt_t microSecs = strtoul(((TObjString*)tokens->At(3))->String(), NULL, 10);

    AliTriggerScalersRecord * rec = new AliTriggerScalersRecord();
    rec->SetTimeStamp( orbit, period, seconds, microSecs );
    TString strLine1;
    for (Int_t i=0; i<nclass; ++i) {
      strLine1.ReadLine(*file);
      if (strLine1.BeginsWith("#")) continue;
      TObjArray *tokens1 = strLine1.Tokenize(" \t");
      Int_t ntokens1 = tokens1->GetEntriesFast();
      if( ntokens1 != 6 ) {
        AliErrorClass( Form( "Error reading scalers from (%s): line (%s)", 
			     filename.Data(), strLine1.Data() ));
	delete rec;
	return rScaler;
      }
      UInt_t LOCB = strtoul(((TObjString*)tokens1->At(0))->String(), NULL, 10);
      UInt_t LOCA = strtoul(((TObjString*)tokens1->At(1))->String(), NULL, 10);
      UInt_t L1CB = strtoul(((TObjString*)tokens1->At(2))->String(), NULL, 10);
      UInt_t L1CA = strtoul(((TObjString*)tokens1->At(3))->String(), NULL, 10);
      UInt_t L2CB = strtoul(((TObjString*)tokens1->At(4))->String(), NULL, 10);
      UInt_t L2CA = strtoul(((TObjString*)tokens1->At(5))->String(), NULL, 10);

      rScaler->GetClass(i);
      rec->AddTriggerScalers( rScaler->GetClass(i),
                              LOCB, LOCA, L1CB,
                              L1CA, L2CB, L2CA );
      delete tokens1;
    } 
    rScaler->AddTriggerScalers( rec );
    
    delete tokens;     
  }
  file->close();
  delete file;

  return  rScaler; 
}
  
//_____________________________________________________________________________
Int_t  AliTriggerRunScalers::FindNearestScalersRecord( AliTimeStamp * stamp )
{
   // Find Trigger scaler record with the closest timestamp <= "stamp"
   // using a binary search. 
   // return the index in the array of records, if the timestamp 
   // is out of range return -1

   Int_t   base, position=-1, last, result = 0;
   AliTimeStamp *op2 = NULL;
   
   fScalersRecord.Sort();

   base = 0;
   last = fScalersRecord.GetEntriesFast();

   while (last >= base) {
      position = (base+last) / 2;
      cout << "pos " <<   position<< " base " <<   base << "last " <<   last << endl;
      AliTriggerScalersRecord* rec = (AliTriggerScalersRecord*)fScalersRecord.At(position);
      if( rec ) op2 = rec->GetTimeStamp();
      if( op2 && (result = stamp->Compare(op2)) == 0  )
         return position;  // exact match
      cout << "result " <<   result << " op2 " << op2 << " rec "<< rec << endl;
      if (!op2 || result < 0)
         last = position-1;
      else
         base = position+1;
      op2 = NULL;  
   }
   if( (position == 0 && result < 0) || position >= fScalersRecord.GetEntriesFast() ) 
    return -1;  // out of range
   else 
    return (result < 0 ) ? position-1 : position; // nearst < stamp   
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::ConsistencyCheck()
{

UInt_t LOCB_2, LOCA_2, L1CB_2, L1CA_2, L2CB_2, L2CA_2, LOCB_1, LOCA_1, L1CB_1, L1CA_1, L2CB_1, L2CA_1;
char *LOCBstat = "NULL", *LOCAstat = "NULL", *L1CBstat = "NULL", *L1CAstat = "NULL", *L2CBstat = "NULL", *L2CAstat = "NULL";
Int_t position = fScalersRecord.GetEntriesFast()-1;

if (position == 0) return 1;

AliTriggerScalersRecord* Scalers2 = (AliTriggerScalersRecord*)fScalersRecord.At(position);
AliTriggerScalersRecord* Scalers1 = (AliTriggerScalersRecord*)fScalersRecord.At(position-1);
if (Scalers2->Compare((AliTriggerScalersRecord*)fScalersRecord.At(position-1)) == -1) return 0;
else for( Int_t i=0; i<fnClasses; ++i ){

TObjArray* ScalersArray2 = (TObjArray*)Scalers2->GetTriggerScalers();
AliTriggerScalers* counters2 = (AliTriggerScalers*)ScalersArray2->At(i);
	LOCB_2 = counters2->GetLOCB();
	LOCA_2 = counters2->GetLOCA();
	L1CB_2 = counters2->GetL1CB();
	L1CA_2 = counters2->GetL1CA();
	L2CB_2 = counters2->GetL2CB();
	L2CA_2 = counters2->GetL2CA();

TObjArray* ScalersArray1 = (TObjArray*)Scalers1->GetTriggerScalers();
AliTriggerScalers* counters1 = (AliTriggerScalers*)ScalersArray1->At(i);
	LOCB_1 = counters1->GetLOCB();
	LOCA_1 = counters1->GetLOCA();
	L1CB_1 = counters1->GetL1CB();
	L1CA_1 = counters1->GetL1CA();
	L2CB_1 = counters1->GetL2CB();
	L2CA_1 = counters1->GetL2CA();

if ( LOCB_2 > LOCB_1 ) LOCBstat = "increase";
else if ( LOCB_2 < LOCB_1 && (LOCB_1 - LOCB_2) > 1000000000) LOCBstat = "overflow";
else return 0;

if ( LOCA_2 > LOCA_1 ) LOCAstat = "increase";
else if ( LOCA_2 < LOCA_1 && (LOCA_1 - LOCA_2) > 1000000000) LOCAstat = "overflow";
else return 0;

if ( L1CB_2 > L1CB_1 ) L1CBstat = "increase";
else if ( L1CB_2 < L1CB_1 && (L1CB_1 - L1CB_2) > 1000000000) L1CBstat = "overflow";
else return 0;

if ( L1CA_2 > L1CA_1 ) L1CAstat = "increase";
else if ( L1CA_2 < L1CA_1 && (L1CA_1 - L1CA_2) > 1000000000) L1CAstat = "overflow";
else return 0;

if ( L2CB_2 > L2CB_1 ) L2CBstat = "increase";
else if ( L2CB_2 < L2CB_1 && (L2CB_1 - L2CB_2) > 1000000000) L2CBstat = "overflow";
else return 0;

if ( L2CA_2 > L2CA_1 ) L2CAstat = "increase";
else if ( L2CA_2 < L2CA_1 && (L2CA_1 - L2CA_2) > 1000000000) L2CAstat = "overflow";
else return 0;


if ( (LOCB_2 - LOCB_1) < (LOCA_2 - LOCA_1) && LOCBstat == "increase" && LOCAstat == "increase" ) return 0;
else if ( (4294967295ul - LOCB_1 + LOCB_2 ) < (LOCA_2 - LOCA_1) && LOCBstat == "overflow" && LOCAstat == "increase" ) return 0;
else if ( (LOCB_2 - LOCB_1) < (4294967295ul - LOCA_1 + LOCA_2) && LOCBstat == "increase" && LOCAstat == "overflow" ) return 0;
else if ( (4294967295ul - LOCB_1 + LOCB_2 ) < (4294967295ul - LOCA_1 + LOCA_2) && LOCBstat == "overflow" && LOCAstat == "overflow" ) return 0;

if ( (LOCA_2 - LOCA_1) < (L1CB_2 - L1CB_1) && LOCAstat == "increase" && L1CBstat == "increase" ) return 0;
else if ( (4294967295ul - LOCA_1 + LOCA_2 ) < (L1CB_2 - L1CB_1) && LOCAstat == "overflow" && L1CBstat == "increase" ) return 0;
else if ( (LOCA_2 - LOCA_1) < (4294967295ul - L1CB_1 + L1CB_2) && LOCAstat == "increase" && L1CBstat == "overflow" ) return 0;
else if ( (4294967295ul - LOCA_1 + LOCA_2 ) < (4294967295ul - L1CB_1 + L1CB_2) && LOCAstat == "overflow" && L1CBstat == "overflow" ) return 0;

if ( (L1CB_2 - L1CB_1) < (L1CA_2 - L1CA_1) && L1CBstat == "increase" && L1CAstat == "increase" ) return 0;
else if ( (4294967295ul - L1CB_1 + L1CB_2 ) < (L1CA_2 - L1CA_1) && L1CBstat == "overflow" && L1CAstat == "increase" ) return 0;
else if ( (L1CB_2 - L1CB_1) < (4294967295ul - L1CA_1 + L1CA_2) && L1CBstat == "increase" && L1CAstat == "overflow" ) return 0;
else if ( (4294967295ul - L1CB_1 + L1CB_2 ) < (4294967295ul - L1CA_1 + L1CA_2) && L1CBstat == "overflow" && L1CAstat == "overflow" ) return 0;

if ( (L1CA_2 - L1CA_1) < (L2CB_2 - L2CB_1) && L1CAstat == "increase" && L2CBstat == "increase" ) return 0;
else if ( (4294967295ul- L1CA_1 + L1CA_2 ) < (L2CB_2 - L2CB_1) && L1CAstat == "overflow" && L2CBstat == "increase" ) return 0;
else if ( (L1CA_2 - L1CA_1) < (4294967295ul- L2CB_1 + L2CB_2) && L1CAstat == "increase" && L2CBstat == "overflow" ) return 0;
else if ( (4294967295ul- L1CA_1 + L1CA_2 ) < (4294967295ul- L2CB_1 + L2CB_2) && L1CAstat == "overflow" && L2CBstat == "overflow" ) return 0;

if ( (L2CB_2 - L2CB_1) < (L2CA_2 - L2CA_1) && L2CBstat == "increase" && L2CAstat == "increase" ) return 0;
else if ( (4294967295ul- L2CB_1 + L2CB_2 ) < (L2CA_2 - L2CA_1) && L2CBstat == "overflow" && L2CAstat == "increase" ) return 0;
else if ( (L2CB_2 - L2CB_1) < (4294967295ul- L2CA_1 + L2CA_2) && L2CBstat == "increase" && L2CAstat == "overflow" ) return 0;
else if ( (4294967295ul- L2CB_1 + L2CB_2 ) < (4294967295ul- L2CA_1 + L2CA_2) && L2CBstat == "overflow" && L2CAstat == "overflow" ) return 0;




cout<<"LOCB_1 =" << LOCB_1 <<"  LOCB_2="<<LOCB_2<<endl;
if ( L1CB_2 < L1CB_1 || L1CA_2 < L1CA_1 || L2CB_2 < L2CB_1 || L2CA_2 < L2CA_1) return 0;
else if ( LOCB_2 < LOCB_1 && ((LOCB_1 - LOCB_2) < 1000000000) ) return 0;//verify this number with Roman!
 
}

return 1;
}
//_____________________________________________________________________________
void AliTriggerRunScalers::Print( const Option_t* ) const
{
   // Print
  cout << "Trigger Scalers Record per Run: " << endl;
  cout << "  File version :" <<  fVersion << endl;            
  cout << "  Run Number :" <<  fRunNumber << endl;          
  cout << "  Number of Classes :" <<  (Int_t)fnClasses << endl;          
  cout << "    Classes ID:";
  for( Int_t i=0; i<fnClasses; ++i ) 
    cout << "  " << (Int_t)fClassIndex[i];       
  cout << endl; 
    
  for( Int_t i=0; i<fScalersRecord.GetEntriesFast(); ++i ) 
     ((AliTriggerScalersRecord*)fScalersRecord.At(i))->Print();
}

