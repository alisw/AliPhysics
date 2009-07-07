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
#include "AliTriggerScalersESD.h"
#include "AliTriggerScalersRecord.h"
#include "AliTriggerScalersRecordESD.h"
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
//______________________________________________________________________________
AliTriggerRunScalers::~AliTriggerRunScalers() 
{
 // Destructor
 fScalersRecord.SetOwner(); 
 fScalersRecord.Delete(); 
 fScalersRecordESD.SetOwner(); 
 fScalersRecordESD.Delete(); 
}
//_____________________________________________________________________________
void AliTriggerRunScalers::AddTriggerScalers( AliTriggerScalersRecord* scaler ) 
{ 
  // Add scaler and check consistency
  fScalersRecord.AddLast( scaler );
  if (AliTriggerRunScalers::ConsistencyCheck(fScalersRecord.GetEntriesFast()-1,kFALSE)){
   AliErrorClass("Trigger counters not in the right order or decreasing!");
  //  scaler->Print();
  //  fScalersRecord.Sort(); 
 }
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
Int_t AliTriggerRunScalers::ConsistencyCheck(Int_t position,Bool_t correctOverflow)
{
   //Check if counters are consistent(increase). Example: lOCB(n) < lOCB(n+1) and lOCB > lOCA
   // scalers coding 0,1,2,3,4,5=0b,0a,1b,1a,2b,2a
   // returns: 
   //         1 = decresing time 
   //         2 = too big jump in scalers, looks like some readings are missing
   //         3 = (level+1) > (level)
   if (position == 0){
      if(correctOverflow){
        AliErrorClass("position=0\n");
        return 1;
      }else return 0; // to work correctlu in AddScalers
   };
   UInt_t c2[6], c1[6];
   ULong64_t c64[6]; 
   Bool_t increase[6], overflow[6];  
   for(Int_t i=0;i<6;i++){increase[i]=0;overflow[i]=0;}
   ULong64_t const max1 = 4294967295ul;  //32bit counters overflow after 4294967295
   ULong64_t const max2 = 1000000000ul;  //when counters overflow they seem to be decreasing. Assume decrease cannot be smaller than max2.

   AliTriggerScalersRecord* scalers2 = (AliTriggerScalersRecord*)fScalersRecord.At(position);
   AliTriggerScalersRecord* scalers1 = (AliTriggerScalersRecord*)fScalersRecord.At(position-1);
   if (scalers2->Compare((AliTriggerScalersRecord*)fScalersRecord.At(position-1)) == -1) return 1;
   
   AliTriggerScalersRecordESD* recESD = 0;
   if(correctOverflow){
     recESD = new AliTriggerScalersRecordESD();
     recESD->SetTimeStamp(scalers2->GetTimeStamp());
   }
   for( Int_t ic=0; ic<fnClasses; ++ic ){
      TObjArray* scalersArray2 = (TObjArray*)scalers2->GetTriggerScalers();
      AliTriggerScalers* counters2 = (AliTriggerScalers*)scalersArray2->At(ic);
      counters2->GetAllScalers(c2);
      TObjArray* scalersArray1 = (TObjArray*)scalers1->GetTriggerScalers();
      AliTriggerScalers* counters1 = (AliTriggerScalers*)scalersArray1->At(ic);
      counters1->GetAllScalers(c1);
      for(Int_t i=0;i<5;i++){
         if ( c2[i] >= c1[i] ) increase[i]=1;
         else if ( c2[i] < c1[i] && (c1[i] - c2[i]) > max2) overflow[i]=1;
         else return 2;
      }
      for(Int_t i=0;i<5;i++){
         if ((c2[i] - c1[i]) < (c2[i+1] - c1[i+1]) && increase[i] && increase[i+1] ) {
                 if ( ((c2[i+1] - c1[i+1]) - (c2[i] - c1[i])) < 16 ) {AliWarningClass("Trigger scaler Level[i+1] > Level[i]. Diff < 16!");}
                 else return 3; }
         else if ( (max1 - c1[i]+c2[i]) < (c2[i+1] - c1[i+1]) && overflow[i] && increase[i+1] ) {
                 if ( ((c2[i+1] - c1[i+1]) - (max1 - c1[i]+c2[i])) < 16 ) {AliWarningClass("Trigger scaler Level[i+1] > Level[i]. Diff < 16!");}
                 else return 3; }
         else if ( (c2[i] - c1[i]) < (max1 - c1[i+1] + c2[i+1]) && increase[i] && overflow[i+1] ) {
                 if ( ((max1 - c1[i+1] + c2[i+1]) - (c2[i] - c1[i])) < 16 ) {AliWarningClass("Trigger scaler Level[i+1] > Level[i]. Diff < 16!");}
                 else return 3; }
         else if ( (max1 - c1[i] + c2[i] ) < (max1 - c1[i+1] + c2[i+1]) && overflow[i] && overflow[i+1] ) {
                 if ( ((max1 - c1[i+1] + c2[i+1]) - (max1 - c1[i] + c2[i] )) < 16 ) {AliWarningClass("Trigger scaler Level[i+1] > Level[i]. Diff < 16!");}
                 else return 3; }
      }
      if(correctOverflow){ 
        for(Int_t i=0;i<6;i++){ c64[i]=c2[i]+max1*overflow[i]; }
        AliTriggerScalersESD* s= new AliTriggerScalersESD(fClassIndex[ic],c64);
        recESD->AddTriggerScalers(s);
         }

 }
 if(correctOverflow)fScalersRecordESD.AddLast(recESD);
 return 0;
}
//____________________________________________________________________________
Int_t AliTriggerRunScalers::CorrectScalersOverflow()
{
 // Run over fScalersRecord, check overflow using CheckConsistency methos
 // and save corrected result in fScalersRecordESD.
 UInt_t c1[6];
 ULong64_t c64[6];
 AliTriggerScalersRecordESD* recESD = new AliTriggerScalersRecordESD();
 // add 0
 AliTriggerScalersRecord* scalers = (AliTriggerScalersRecord*)fScalersRecord.At(0);
 for( Int_t ic=0; ic<fnClasses; ++ic ){
    TObjArray* scalersArray = (TObjArray*)scalers->GetTriggerScalers();
    AliTriggerScalers* counters = (AliTriggerScalers*)scalersArray->At(ic);
    counters->GetAllScalers(c1);
    for(Int_t i=0; i<6; i++)c64[i]=c1[i];
    AliTriggerScalersESD* s= new AliTriggerScalersESD(fClassIndex[ic],c64);
    recESD->AddTriggerScalers(s);
 }
 fScalersRecordESD.AddLast(recESD);
 for(Int_t i=1;i<fScalersRecord.GetEntriesFast(); i++){
  if(ConsistencyCheck(i,kTRUE)){
    fScalersRecord.At(i)->Print();
    fScalersRecord.At(i-1)->Print();
    fScalersRecordESD.SetOwner(); 
    fScalersRecordESD.Delete(); 
    AliErrorClass("Inconsistent scalers, they will not be provided.\n");
    return 1;
  }
 }
 if(fScalersRecordESD.GetEntriesFast() != fScalersRecord.GetEntriesFast()){
    AliErrorClass("Internal error: #scalers ESD != #scalers \n");
    return 1;
 }
 return 0;
}
//_____________________________________________________________________________
AliTriggerScalersESD* AliTriggerRunScalers::GetScalersForEventClass(const AliTimeStamp* stamp,const Int_t classIndex) const
{
 // Find scalers for event for class in fScalersRecordESD
 // Assumes that fScalerRecord = fScalerRecordESD
 Int_t position = FindNearestScalersRecord(stamp);
 if ( position == -1 ) { 
  AliErrorClass("Event AliTimeStamp out of range!");
  return 0; 
 }
 // check also position=max
 AliTriggerScalersRecordESD* scalrec1 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(position);
 AliTriggerScalersRecordESD* scalrec2 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(position+1);
 TObjArray* scalers1 = (TObjArray*)scalrec1->GetTriggerScalers();
 TObjArray* scalers2 = (TObjArray*)scalrec2->GetTriggerScalers();
 cout << " Position = " << position << endl;
 if(scalers1->GetEntriesFast() != fnClasses){
  AliErrorClass("Internal error: #classes in RecordESD != fnClasses\n");
  return 0; 
 }

 AliTriggerScalersESD *s1,*s2;
 for ( Int_t ic=0; ic < (Int_t)fnClasses; ++ic ){

      s1 = (AliTriggerScalersESD*)scalers1->At(ic);
      s2 = (AliTriggerScalersESD*)scalers2->At(ic);

      Bool_t classfound = (s1->GetClassIndex() == classIndex) && (s2->GetClassIndex() == classIndex);
      if(classfound){
        ULong64_t max = 4294967295ul;
        AliTriggerScalersRecordESD* scalrec0 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(0);
        TObjArray* scalers0 = (TObjArray*)scalrec0->GetTriggerScalers();
        AliTriggerScalersESD *s0 = (AliTriggerScalersESD*)scalers0->At(ic);
	ULong64_t base[6],c1[6],c2[6],cint[6];
        ULong64_t orbit = max*(stamp->GetPeriod()) + stamp->GetOrbit();
	s0->GetAllScalers(base);
	s1->GetAllScalers(c1);
	s2->GetAllScalers(c2);
	ULong64_t orbit1 = max*(scalrec1->GetTimeStamp()->GetPeriod())+scalrec1->GetTimeStamp()->GetOrbit();
	ULong64_t orbit2 = max*(scalrec2->GetTimeStamp()->GetPeriod())+scalrec2->GetTimeStamp()->GetOrbit();
        for(Int_t i=0;i<6;i++){
	   Double_t slope=Double_t(c2[i]-c1[i])/Double_t(orbit2-orbit1);
	   cint[i]=ULong64_t(slope*(orbit-orbit1)) +c1[i] -base[i];
	}
	AliTriggerScalersESD* result = new AliTriggerScalersESD(classIndex,cint);
        return result;
      }
 }
 AliErrorClass(Form("Classindex %i not found.\n",classIndex));
 return 0;
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

