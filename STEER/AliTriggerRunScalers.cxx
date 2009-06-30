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
// ReadScalers(): read the txt file (rXXXX.cnt) provided by CTP and creates array of records.
// ConsistencyCheck(): checks if the last scaler added to record is consistent, 
// i.e. time increases and valued are not decreasing in time.
// 
//
//////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <Riostream.h>

//#include <TObject.h>
//#include <TArrayC.h>
//#include <TFile.h>
#include <TClass.h>
#include <TSystem.h>
#include <TObjString.h>
#include <TObjArray.h>

#include "AliLog.h"  
#include "AliTimeStamp.h"
#include "AliTriggerScalers.h"
#include "AliTriggerScalersRecord.h"
#include "AliTriggerRunScalers.h"

ClassImp( AliTriggerRunScalers )

//_____________________________________________________________________________
AliTriggerRunScalers::AliTriggerRunScalers():
  TObject(),
  fVersion(0),
  fRunNumber(0),
  fnClasses(0),
  fClassIndex(),                    
  fScalersRecord(),
  fScalersRecordESD()
{
  // Default constructor
}

//_____________________________________________________________________________
void AliTriggerRunScalers::AddTriggerScalers( AliTriggerScalersRecord* scaler ) 
{ 
  // Add scaler and check consistency
  fScalersRecord.AddLast( scaler );
  if (!AliTriggerRunScalers::ConsistencyCheck()) AliErrorClass("Trigger counters not in the right order or decreasing!");
//  fScalersRecord.Sort(); 
}
//_____________________________________________________________________________

AliTriggerRunScalers::AliTriggerRunScalers(const AliTriggerRunScalers &run) :
 TObject(),
 fVersion(run.fVersion),
 fRunNumber(run.fRunNumber),
 fnClasses(run.fnClasses),
 fClassIndex(),                    
 fScalersRecord(),
 fScalersRecordESD()
{
// copy constructor
for (Int_t i = 0; i < run.fClassIndex.GetSize(); i++) {
    if (run.fClassIndex[i]) fClassIndex.AddAt(run.fClassIndex[i], i);
  }
for (Int_t i = 0; i < run.fScalersRecord.GetEntriesFast(); i++) {
    if (run.fScalersRecord[i]) fScalersRecord.Add(run.fScalersRecord[i]->Clone());
  }
for (Int_t i = 0; i < run.fScalersRecordESD.GetEntriesFast(); i++) {
    if (run.fScalersRecordESD[i]) fScalersRecordESD.Add(run.fScalersRecordESD[i]->Clone());
  }

}
//_____________________________________________________________________________
AliTriggerRunScalers &AliTriggerRunScalers::operator=(const AliTriggerRunScalers& run)
{
// assignment operator
if(&run == this) return *this;
((TObject *)this)->operator=(run);

fVersion = run.fVersion;
fRunNumber = run.fRunNumber;
fnClasses = run.fnClasses;

for (Int_t i = 0; i < run.fClassIndex.GetSize(); i++) {
    if (run.fClassIndex[i]) fClassIndex.AddAt(run.fClassIndex[i], i);
  }
for (Int_t i = 0; i < run.fScalersRecord.GetEntriesFast(); i++) {
    if (run.fScalersRecord[i]) fScalersRecord.Add(run.fScalersRecord[i]->Clone());
  }
for (Int_t i = 0; i < run.fScalersRecordESD.GetEntriesFast(); i++) {
    if (run.fScalersRecordESD[i]) fScalersRecordESD.Add(run.fScalersRecordESD[i]->Clone());
  }
return *this;
} 
//_____________________________________________________________________________
AliTriggerRunScalers* AliTriggerRunScalers::ReadScalers( TString & filename )
{
  // Read scalers from text file(.cnt) provided by CTP 
  // for given run and convert it to root format 
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

      UInt_t lOCB = strtoul(((TObjString*)tokens1->At(0))->String(), NULL, 10);
      UInt_t lOCA = strtoul(((TObjString*)tokens1->At(1))->String(), NULL, 10);
      UInt_t l1CB = strtoul(((TObjString*)tokens1->At(2))->String(), NULL, 10);
      UInt_t l1CA = strtoul(((TObjString*)tokens1->At(3))->String(), NULL, 10);
      UInt_t l2CB = strtoul(((TObjString*)tokens1->At(4))->String(), NULL, 10);
      UInt_t l2CA = strtoul(((TObjString*)tokens1->At(5))->String(), NULL, 10);

      rScaler->GetClass(i);
      rec->AddTriggerScalers( rScaler->GetClass(i),
                              lOCB, lOCA, l1CB,
                              l1CA, l2CB, l2CA );

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
Int_t  AliTriggerRunScalers::FindNearestScalersRecord( const AliTimeStamp *stamp ) const
{
   // Find Trigger scaler record with the closest timestamp <= "stamp"
   // using a binary search. 
   // return the index in the array of records, if the timestamp 
   // is out of range return -1

   Int_t   base, position=-1, last, result = 0;
   Int_t op2 = 0;
   
   //fScalersRecord.Sort();

   base = 0;
   last = fScalersRecord.GetEntriesFast();

   while (last >= base) {
      position = (base+last) / 2;
      cout << "pos " <<   position<< " base " <<   base << "last " <<   last << endl;
      AliTriggerScalersRecord* rec = (AliTriggerScalersRecord*)fScalersRecord.At(position);
      if( rec && rec->GetTimeStamp()) op2 = 1;
      if( op2 && (result = stamp->Compare(rec->GetTimeStamp())) == 0  )
         return position;  // exact match
      cout << "result " <<   result << " op2 " << op2 << " rec "<< rec << endl;
      if (!op2 || result < 0)
         last = position-1;
      else
         base = position+1;
      op2 = 0;  
   }
   if( (position == 0 && result < 0) || position >= fScalersRecord.GetEntriesFast() ) 
    return -1;  // out of range
   else 
    return (result < 0 ) ? position-1 : position; // nearst < stamp   
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::ConsistencyCheck() const
{
   //Check if counters are consistent(increase). Example: lOCB(n) < lOCB(n+1) and lOCB > lOCA
   UInt_t lOCBtwo, lOCAtwo, l1CBtwo, l1CAtwo, l2CBtwo, l2CAtwo, lOCBone, lOCAone, l1CBone, l1CAone, l2CBone, l2CAone;
   Bool_t increase0B=0, increase0A=0, increase1B=0,  increase1A=0, increase2B=0, increase2A=0;
   Bool_t overflow0B=0, overflow0A=0, overflow1B=0,  overflow1A=0, overflow2B=0, overflow2A=0;
   Int_t position = fScalersRecord.GetEntriesFast()-1;
   if (position == 0) return 1;

   AliTriggerScalersRecord* scalers2 = (AliTriggerScalersRecord*)fScalersRecord.At(position);
   AliTriggerScalersRecord* scalers1 = (AliTriggerScalersRecord*)fScalersRecord.At(position-1);
   if (scalers2->Compare((AliTriggerScalersRecord*)fScalersRecord.At(position-1)) == -1) return 0;
   else for( Int_t i=0; i<fnClasses; ++i ){

   TObjArray* scalersArray2 = (TObjArray*)scalers2->GetTriggerScalers();
   AliTriggerScalers* counters2 = (AliTriggerScalers*)scalersArray2->At(i);
	lOCBtwo = counters2->GetLOCB();
	lOCAtwo = counters2->GetLOCA();
	l1CBtwo = counters2->GetL1CB();
	l1CAtwo = counters2->GetL1CA();
	l2CBtwo = counters2->GetL2CB();
	l2CAtwo = counters2->GetL2CA();

   TObjArray* scalersArray1 = (TObjArray*)scalers1->GetTriggerScalers();
   AliTriggerScalers* counters1 = (AliTriggerScalers*)scalersArray1->At(i);
	lOCBone = counters1->GetLOCB();
	lOCAone = counters1->GetLOCA();
	l1CBone = counters1->GetL1CB();
	l1CAone = counters1->GetL1CA();
	l2CBone = counters1->GetL2CB();
	l2CAone = counters1->GetL2CA();

   UInt_t const max1 = 4294967295ul;  //32bit counters overflow after 4294967295
   UInt_t const max2 = 1000000000ul;  //when counters overflow they seem to be decreasing. Assume decrease cannot be smaller than max2.

   if ( lOCBtwo > lOCBone ) increase0B=1;
   else if ( lOCBtwo < lOCBone && (lOCBone - lOCBtwo) > max2) overflow0B=1;
   else return 0;

   if ( lOCAtwo > lOCAone ) increase0A=1;
   else if ( lOCAtwo < lOCAone && (lOCAone - lOCAtwo) > max2) overflow0A=1;
   else return 0;

   if ( l1CBtwo > l1CBone ) increase1B=1;
   else if ( l1CBtwo < l1CBone && (l1CBone - l1CBtwo) > max2) overflow1B=1;
   else return 0;

   if ( l1CAtwo > l1CAone ) increase1A=1;
   else if ( l1CAtwo < l1CAone && (l1CAone - l1CAtwo) > max2) overflow1A=1;
   else return 0;

   if ( l2CBtwo > l2CBone ) increase2B=1;
   else if ( l2CBtwo < l2CBone && (l2CBone - l2CBtwo) > max2) overflow2B=1;
   else return 0;

   if ( l2CAtwo > l2CAone ) increase2A=1;
   else if ( l2CAtwo < l2CAone && (l2CAone - l2CAtwo) > max2) overflow2A=1;
   else return 0;


   if ( (lOCBtwo - lOCBone) < (lOCAtwo - lOCAone) && increase0B && increase0A ) return 0;
   else if ( (max1 - lOCBone + lOCBtwo ) < (lOCAtwo - lOCAone) && overflow0B && increase0A ) return 0;
   else if ( (lOCBtwo - lOCBone) < (max1 - lOCAone + lOCAtwo) && increase0B && overflow0A ) return 0;
   else if ( (max1 - lOCBone + lOCBtwo ) < (max1 - lOCAone + lOCAtwo) && overflow0B && overflow0A ) return 0;
 
   if ( (lOCAtwo - lOCAone) < (l1CBtwo - l1CBone) && increase0A && increase1B ) return 0;
   else if ( (max1 - lOCAone + lOCAtwo ) < (l1CBtwo - l1CBone) && overflow0A && increase1B ) return 0;
   else if ( (lOCAtwo - lOCAone) < (max1 - l1CBone + l1CBtwo) && increase0A && overflow1B ) return 0;
   else if ( (max1 - lOCAone + lOCAtwo ) < (max1 - l1CBone + l1CBtwo) && overflow0A && overflow1B ) return 0;

   if ( (l1CBtwo - l1CBone) < (l1CAtwo - l1CAone) && increase1B && increase1A ) return 0;
   else if ( (max1 - l1CBone + l1CBtwo ) < (l1CAtwo - l1CAone) && overflow1B && increase1A ) return 0;
   else if ( (l1CBtwo - l1CBone) < (max1 - l1CAone + l1CAtwo) && increase1B && overflow1A ) return 0;
   else if ( (max1 - l1CBone + l1CBtwo ) < (max1 - l1CAone + l1CAtwo) && overflow1B && overflow1A ) return 0;

   if ( (l1CAtwo - l1CAone) < (l2CBtwo - l2CBone) && increase1A && increase2B ) return 0;
   else if ( (max1 - l1CAone + l1CAtwo ) < (l2CBtwo - l2CBone) && overflow1A && increase2B ) return 0;
   else if ( (l1CAtwo - l1CAone) < (max1 - l2CBone + l2CBtwo) && increase1A && overflow2B ) return 0;
   else if ( (max1 - l1CAone + l1CAtwo ) < (max1 - l2CBone + l2CBtwo) && overflow1A && overflow2B ) return 0;

   if ( (l2CBtwo - l2CBone) < (l2CAtwo - l2CAone) && increase2B && increase2A ) return 0;
   else if ( (max1 - l2CBone + l2CBtwo ) < (l2CAtwo - l2CAone) && overflow2B && increase2A ) return 0;
   else if ( (l2CBtwo - l2CBone) < (max1 - l2CAone + l2CAtwo) && increase2B && overflow2A ) return 0;
   else if ( (max1 - l2CBone + l2CBtwo ) < (max1 - l2CAone + l2CAtwo) && overflow2B && overflow2A ) return 0;

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

