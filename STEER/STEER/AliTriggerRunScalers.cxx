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
#include <TMath.h>
#include <TGraphErrors.h>

#include "AliLog.h"  
#include "AliTimeStamp.h"
#include "AliTriggerScalers.h"
#include "AliTriggerScalersESD.h"
#include "AliTriggerScalersRecord.h"
#include "AliTriggerScalersRecordESD.h"
#include "AliTriggerRunScalers.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"
#include "AliTriggerBCMask.h"

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
  UInt_t* overflow[6];
  for(Int_t i=0; i<6; i++) {
     overflow[i] = new UInt_t[fnClasses];
     for(Int_t j=0; j<fnClasses; j++) overflow[i][j] = 0;
  }
  if (AliTriggerRunScalers::ConsistencyCheck(fScalersRecord.GetEntriesFast()-1,kFALSE,overflow)){
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
  fClassIndex.Set(run.fClassIndex.GetSize());
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
      AliDebug(1, Form("position= %d   base= %d    last= %d  ",position,base,last));
      AliTriggerScalersRecord* rec = (AliTriggerScalersRecord*)fScalersRecord.At(position);
      if( rec && rec->GetTimeStamp()) op2 = 1;
      if( op2 && (result = stamp->Compare(rec->GetTimeStamp())) == 0  )
         return position;  // exact match 
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
Int_t AliTriggerRunScalers::ConsistencyCheck(Int_t position,Bool_t correctOverflow, UInt_t** overflow)
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
      }else return 0; // to work correctly in AddScalers
   };
   UInt_t c2[6], c1[6];
   ULong64_t c64[6]; 
   Bool_t increase[6];  
   for(Int_t i=0;i<6;i++){increase[i]=0;}
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
      UChar_t iclass = counters2->GetClassIndex();
      counters2->GetAllScalers(c2);
      TObjArray* scalersArray1 = (TObjArray*)scalers1->GetTriggerScalers();
      AliTriggerScalers* counters1 = (AliTriggerScalers*)scalersArray1->At(ic);
      counters1->GetAllScalers(c1);
      for(Int_t i=0;i<5;i++){
         if ( c2[i] >= c1[i] ) increase[i]=1;
         else if ( c2[i] < c1[i] && (c1[i] - c2[i]) > max2) overflow[i][ic]++;
         else return 2;
      }
      for(Int_t i=0;i<5;i++){
         if ((c2[i] - c1[i]) < (c2[i+1] - c1[i+1]) && increase[i] && increase[i+1] ) {
                 if ( ((c2[i+1] - c1[i+1]) - (c2[i] - c1[i])) < 16 ) {AliWarningClass("Trigger scaler Level[i+1] > Level[i]. Diff < 16!");}
                 else return 3; }
         else if ( (max1 - c1[i]+c2[i]) < (c2[i+1] - c1[i+1]) && overflow[i][ic] && increase[i+1] ) {
                 if ( ((c2[i+1] - c1[i+1]) - (max1 - c1[i]+c2[i])) < 16 ) {AliWarningClass("Trigger scaler Level[i+1] > Level[i]. Diff < 16!");}
                 else return 3; }
         else if ( (c2[i] - c1[i]) < (max1 - c1[i+1] + c2[i+1]) && increase[i] && overflow[i+1][ic] ) {
                 if ( ((max1 - c1[i+1] + c2[i+1]) - (c2[i] - c1[i])) < 16 ) {AliWarningClass("Trigger scaler Level[i+1] > Level[i]. Diff < 16!");}
                 else return 3; }
         else if ( (max1 - c1[i] + c2[i] ) < (max1 - c1[i+1] + c2[i+1]) && overflow[i][ic] && overflow[i+1][ic] ) {
                 if ( ((max1 - c1[i+1] + c2[i+1]) - (max1 - c1[i] + c2[i] )) < 16 ) {AliWarningClass("Trigger scaler Level[i+1] > Level[i]. Diff < 16!");}
                 else return 3; }
      }
      if(correctOverflow){ 
        for(Int_t i=0;i<6;i++){ c64[i]=c2[i]+max1*overflow[i][ic]; }
        AliTriggerScalersESD* s= new AliTriggerScalersESD(iclass,c64);
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

 // Temporary fix for the OCDB entries written with v4-16-Release
 // where the wrong sorting was used
 fScalersRecord.Sort();
 UInt_t c1[6];
 ULong64_t c64[6];
 AliTriggerScalersRecordESD* recESD = new AliTriggerScalersRecordESD();
 // add 0
 AliTriggerScalersRecord* scalers = (AliTriggerScalersRecord*)fScalersRecord.At(0);
 recESD->SetTimeStamp(scalers->GetTimeStamp());
 for( Int_t ic=0; ic<fnClasses; ++ic ){
    TObjArray* scalersArray = (TObjArray*)scalers->GetTriggerScalers();
    AliTriggerScalers* counters = (AliTriggerScalers*)scalersArray->At(ic);
    counters->GetAllScalers(c1);
    UChar_t iclass = counters->GetClassIndex();
    for(Int_t i=0; i<6; i++)c64[i]=c1[i];
    AliTriggerScalersESD* s= new AliTriggerScalersESD(iclass,c64);
    recESD->AddTriggerScalers(s);
 }
 fScalersRecordESD.AddLast(recESD);

 UInt_t* overflow[6];
 for(Int_t i=0; i<6; i++) {
    overflow[i] = new UInt_t[fnClasses];
    for(Int_t j=0; j<fnClasses; j++) overflow[i][j] = 0;
 }

 for(Int_t i=1;i<fScalersRecord.GetEntriesFast(); i++){
  if(ConsistencyCheck(i,kTRUE,overflow)){
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
        ULong64_t max = 16777216ul;
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
const AliTriggerScalersRecordESD* AliTriggerRunScalers::GetScalersDeltaForEvent(const AliTimeStamp* stamp) const
{
 // Find scalers for event for class in fScalersRecordESD
 // Assumes that fScalerRecord = fScalerRecordESD
 Int_t position = FindNearestScalersRecord(stamp);
 if ( position == -1 ) { 
  AliErrorClass("Event AliTimeStamp out of range!");
  return 0; 
 }
 // check also position=max
 if (fScalersRecordESD.GetEntriesFast()<2) return 0;
 
 AliTriggerScalersRecordESD* scalrec1 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(position);
 AliTriggerScalersRecordESD* scalrec2 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(position+1);
 TObjArray* scalers1 = (TObjArray*)scalrec1->GetTriggerScalers();
 TObjArray* scalers2 = (TObjArray*)scalrec2->GetTriggerScalers();

 if(scalers1->GetEntriesFast() != fnClasses){
  AliErrorClass("Internal error: #classes in RecordESD != fnClasses\n");
  return 0; 
 }
 AliTriggerScalersRecordESD *scalrec = new AliTriggerScalersRecordESD();
 AliTriggerScalersESD *s1,*s2;
 for ( Int_t ic=0; ic < (Int_t)fnClasses; ++ic ){

  s1 = (AliTriggerScalersESD*)scalers1->At(ic);
  s2 = (AliTriggerScalersESD*)scalers2->At(ic);

  ULong64_t c1[6],c2[6],cint[6];
  s1->GetAllScalers(c1);
  s2->GetAllScalers(c2);
  for(Int_t i=0;i<6;i++){
     cint[i]=c2[i]-c1[i];
  }
  AliTriggerScalersESD* result = new AliTriggerScalersESD(s1->GetClassIndex(),cint);
  scalrec->AddTriggerScalers(result);
 }
 UInt_t max = 16777216ul;
 UInt_t orbit, period;
 
 if (scalrec2->GetTimeStamp()->GetOrbit() > scalrec1->GetTimeStamp()->GetOrbit()) {
  orbit = scalrec2->GetTimeStamp()->GetOrbit() - scalrec1->GetTimeStamp()->GetOrbit();
  period = scalrec2->GetTimeStamp()->GetPeriod() - scalrec1->GetTimeStamp()->GetPeriod();
 }
 else {
  orbit = max - (scalrec1->GetTimeStamp()->GetOrbit() - scalrec2->GetTimeStamp()->GetOrbit());
  period = scalrec2->GetTimeStamp()->GetPeriod() - scalrec1->GetTimeStamp()->GetPeriod() - 1;
 }
 
 AliTimeStamp *timestamp = new AliTimeStamp(orbit, period, 0);
 scalrec->SetTimeStamp(timestamp);
 delete timestamp;
 return scalrec;
}
//_____________________________________________________________________________
const AliTriggerScalersRecordESD* AliTriggerRunScalers::GetScalersDeltaForRun() const
{
 // Find scalers for event for class in fScalersRecordESD
 // Assumes that fScalerRecord = fScalerRecordESD
 if (fScalersRecordESD.GetEntriesFast()<2) return 0;

 AliTriggerScalersRecordESD* scalrec1 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(0);
 AliTriggerScalersRecordESD* scalrec2 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(fScalersRecord.GetEntriesFast()-1);
 TObjArray* scalers1 = (TObjArray*)scalrec1->GetTriggerScalers();
 TObjArray* scalers2 = (TObjArray*)scalrec2->GetTriggerScalers();

 if(scalers1->GetEntriesFast() != fnClasses){
  AliErrorClass("Internal error: #classes in RecordESD != fnClasses\n");
  return 0; 
 }
 AliTriggerScalersESD *s1,*s2;
 AliTriggerScalersRecordESD *scalrec = new AliTriggerScalersRecordESD();
 for ( Int_t ic=0; ic < (Int_t)fnClasses; ++ic ){

  s1 = (AliTriggerScalersESD*)scalers1->At(ic);
  s2 = (AliTriggerScalersESD*)scalers2->At(ic);

  ULong64_t c1[6],c2[6],cint[6];
  s1->GetAllScalers(c1);
  s2->GetAllScalers(c2);
  for(Int_t i=0;i<6;i++){
     cint[i]=c2[i]-c1[i];
  }
  AliTriggerScalersESD* result = new AliTriggerScalersESD(s1->GetClassIndex(),cint);
  scalrec->AddTriggerScalers(result);
 }
 UInt_t max = 16777216ul;
 UInt_t orbit, period;
 
 if (scalrec2->GetTimeStamp()->GetOrbit() > scalrec1->GetTimeStamp()->GetOrbit()) {
  orbit = scalrec2->GetTimeStamp()->GetOrbit() - scalrec1->GetTimeStamp()->GetOrbit();
  period = scalrec2->GetTimeStamp()->GetPeriod() - scalrec1->GetTimeStamp()->GetPeriod();
 }
 else {
  orbit = max - (scalrec1->GetTimeStamp()->GetOrbit() - scalrec2->GetTimeStamp()->GetOrbit());
  period = scalrec2->GetTimeStamp()->GetPeriod() - scalrec1->GetTimeStamp()->GetPeriod() - 1;
 }
 
 AliTimeStamp *timestamp = new AliTimeStamp(orbit, period, 0);
 scalrec->SetTimeStamp(timestamp);
 delete timestamp;
 return scalrec;
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::CalculateMu(Double_t &mu, Double_t &errmu, ULong64_t countsB, ULong64_t countsAC, UShort_t nB, UShort_t nAC, UInt_t orbits, Bool_t bkgCorr, Double_t triggerEff, Double_t errorEff) 
{
 
   if (nB!=0 && orbits!=0)  {
      Double_t pB = (Double_t)countsB/(nB*orbits); // probability for B trigger
      if (!bkgCorr || nAC==0 ) {
         mu = -log(1-pB)/triggerEff;
         errmu = TMath::Sqrt(pB/((1-pB)*nB*orbits) + mu*mu*errorEff*errorEff/(triggerEff*triggerEff)); //
         return kTRUE;
      }
      else 
      {
         Double_t pAC = (Double_t)countsAC/(nAC*orbits); // probability for AC trigger (background)
         mu = ( log(1.-pAC) - log(1.-pB) )/triggerEff;
         // error
         errmu =  TMath::Sqrt(pB/((1.-pB)*nB*orbits) + pAC/((1.-pAC)*nAC*orbits) /*- 2*TMath::Sqrt(pB*pAC/(nB*nAC*(1.-pB)*(1.-pAC)))/orbits*/ + mu*mu*errorEff*errorEff/(triggerEff*triggerEff)); // assume no correlation between B and AC rates, hence no cov term in error
         return kTRUE;
      }
   }
   return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::CalculateMu(Double_t &mu, Double_t &errmu, ULong64_t countsB, ULong64_t countsAC, ULong64_t beamB, UShort_t nB, UShort_t nAC, Bool_t bkgCorr, Double_t triggerEff, Double_t errorEff) 
{

   if (beamB!=0)  {
      Double_t pB = (Double_t)countsB/beamB; // probability for B trigger
      if (!bkgCorr || nAC==0 || nB==0) {
         mu = -log(1-pB)/triggerEff;
         errmu = TMath::Sqrt(pB/((1-pB)*beamB) + mu*mu*errorEff*errorEff/(triggerEff*triggerEff)); //
         return kTRUE;
      }
      else
      {
         Double_t pAC = (Double_t)countsAC/(nAC*beamB/nB); // probability for AC trigger (background)
         mu = ( log(1-pAC) - log(1-pB) )/triggerEff;
         // error
         errmu =  TMath::Sqrt(pB/((1-pB)*beamB) + pAC/((1-pAC)*nAC*beamB/nB) + mu*mu*errorEff*errorEff/(triggerEff*triggerEff));
         return kTRUE;
      }
   }
   return kFALSE;
}
//_____________________________________________________________________________
ULong64_t AliTriggerRunScalers::GetDeltaScaler(const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, Int_t classIndex, TString level)
{
  const AliTriggerScalersESD* scalers1 = scalRec1->GetTriggerScalersForClass(classIndex);
  const AliTriggerScalersESD* scalers2 = scalRec2->GetTriggerScalersForClass(classIndex);
  ULong64_t s1=0;
  ULong64_t s2=0;

  if (level == "l0b") {s1=scalers1->GetLOCB(); s2=scalers2->GetLOCB();}
  else if (level == "l0a") {s1=scalers1->GetLOCA(); s2=scalers2->GetLOCA();}
  else if (level == "l1b") {s1=scalers1->GetL1CB(); s2=scalers2->GetL1CB();}
  else if (level == "l1a") {s1=scalers1->GetL1CA(); s2=scalers2->GetL1CA();}
  else if (level == "l2b") {s1=scalers1->GetL2CB(); s2=scalers2->GetL2CB();}
  else if (level == "l2a") {s1=scalers1->GetL2CA(); s2=scalers2->GetL2CA();}
  else return 0;

  return s2-s1;
}
//_____________________________________________________________________________
Double_t AliTriggerRunScalers::GetDeltaTime(const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2)
{
  const AliTimeStamp* stamp1 = scalRec1->GetTimeStamp();
  const AliTimeStamp* stamp2 = scalRec2->GetTimeStamp();
  UInt_t orbit1 = stamp1->GetOrbit();
  UInt_t orbit2 = stamp2->GetOrbit();
  UInt_t period1 = stamp1->GetPeriod();
  UInt_t period2 = stamp2->GetPeriod();
  UInt_t max = 16777216;  // The period counter increases when 24 bit orbit counter overflow
  Double_t orbitSec = 89.1*1.e-6; // Length of 1 orbit in seconds
  return ((period2 - period1)*max + (orbit2-orbit1))*orbitSec;
}
//_____________________________________________________________________________
UInt_t AliTriggerRunScalers::GetDeltaOrbits(const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2)
{
  const AliTimeStamp* stamp1 = scalRec1->GetTimeStamp();
  const AliTimeStamp* stamp2 = scalRec2->GetTimeStamp();
  UInt_t orbit1 = stamp1->GetOrbit();
  UInt_t orbit2 = stamp2->GetOrbit();
  UInt_t period1 = stamp1->GetPeriod();
  UInt_t period2 = stamp2->GetPeriod();
  UInt_t max = 16777216;  // The period counter increases when 24 bit orbit counter overflow
  return (period2 - period1)*max + (orbit2-orbit1);
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::GetScalerRate(Double_t &rate, Double_t &error, const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, Int_t classIndex, TString level)
{
  if (level != "l0b" && level != "l0a" && level != "l1b" && level != "l1a" && level != "l2b" && level != "l2a") return kFALSE;
 
  ULong64_t scaler = GetDeltaScaler(scalRec1, scalRec2, classIndex, level);
  Double_t time = GetDeltaTime(scalRec1, scalRec2 );
  if (time==0.) return kFALSE;
  rate = (Double_t)scaler/time;
  error = (Double_t)sqrt(scaler)/time;
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::GetScalerRatePerBC(Double_t &rate, Double_t &error, const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, AliTriggerConfiguration* cfg, Int_t classIndex, TString level)
{
  if (level != "l0b" && level != "l0a" && level != "l1b" && level != "l1a" && level != "l2b" && level != "l2a") return kFALSE;
  const AliTriggerClass* trgclass = cfg->GetTriggerClass(classIndex);
  AliTriggerBCMask* bcMask = new AliTriggerBCMask();
  if (trgclass) bcMask = (AliTriggerBCMask*)trgclass->GetBCMask();
  Int_t nBC=0;
  if (TString(bcMask->GetName()).CompareTo("NONE")!=0){
     nBC = (UShort_t)bcMask->GetNUnmaskedBCs();
  }
  if (nBC<1) return kFALSE;
  ULong64_t scaler = GetDeltaScaler(scalRec1, scalRec2, classIndex, level);
  Double_t time = GetDeltaTime(scalRec1, scalRec2 );
  if (time==0.) return kFALSE;
  rate = (Double_t)scaler/time/nBC;
  error = (Double_t)sqrt(scaler)/time/nBC;
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::GetClassL2L0(Double_t &l2l0, Double_t &error, const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, Int_t classIndex)
{
  ULong64_t l0 = GetDeltaScaler(scalRec1, scalRec2, classIndex, "l0b");
  ULong64_t l2 = GetDeltaScaler(scalRec1, scalRec2, classIndex, "l2a");
  if (l0!=0) {
    l2l0 = (Double_t)l2/l0;
    error = (Double_t)sqrt(l2-l2*l2/l0)/l0;
    return kTRUE;
  }
  return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::GetMuFromClassScaler(Double_t &mu, Double_t &errmu, const char* className, const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, AliTriggerConfiguration* cfg, Bool_t colBCsFromFillScheme, Bool_t bkgCorr, Double_t triggerEff, Double_t errorEff)
{
  // className = the first part of the class name. For example CINT1 for the class CINT1-B-NOPF-ALL
  // colBCsFromFillScheme=kTRUE - use filling scheme and orbit counter
  // colBCsFromFillScheme=kFALSE - use cbeamb scaler to get number of colliding BCs
  // One can also switch background correction ON or OFF with bkgCorr
  TObjArray cint1bNames;
  TObjArray cint1acNames;
  TObjArray cbeambNames;
  //
  cint1bNames.Add(new TNamed(Form("%s-ABCE",className),NULL));
  cint1bNames.Add(new TNamed(Form("%s-B",className),NULL));
  cint1acNames.Add(new TNamed(Form("%s-AC",className),NULL));
  cbeambNames.Add(new TNamed("CBEAMB",NULL));
  cbeambNames.Add(new TNamed("CTRUE",NULL));
  //
  Int_t cint1bIndex=-1, cint1acIndex=-1, cbeambIndex=-1;
  TString nameString;
  //
  for (Int_t i=0; i<cfg->GetClasses().GetEntriesFast(); i++ ) {
     nameString = cfg->GetClassNameFromIndex(i);
     for (Int_t j=0; j<cint1bNames.GetEntriesFast(); j++ ) {
       if (nameString.BeginsWith(cint1bNames.At(j)->GetName())) cint1bIndex = i;
     }
     for (Int_t j=0; j<cint1acNames.GetEntriesFast(); j++ ) {
       if (nameString.BeginsWith(cint1acNames.At(j)->GetName())) cint1acIndex = i;
     }
     for (Int_t j=0; j<cbeambNames.GetEntriesFast(); j++ ) {
       if (nameString.BeginsWith(cbeambNames.At(j)->GetName())) cbeambIndex = i;
     }
     nameString.Clear();
  }
  //
  ULong64_t cint1b = 0;
  UShort_t nB=0;
  if (cint1bIndex!=-1) {
     cint1b=GetDeltaScaler(scalRec1, scalRec2, cint1bIndex, "l0b");
     const AliTriggerClass* cint1bClass = cfg->GetTriggerClass(cint1bIndex);
     AliTriggerBCMask* cint1bBCMask = new AliTriggerBCMask();
     if (cint1bClass) cint1bBCMask = (AliTriggerBCMask*)cint1bClass->GetBCMask();
     if (TString(cint1bBCMask->GetName()).CompareTo("NONE")!=0){
        nB = (UShort_t)cint1bBCMask->GetNUnmaskedBCs();
     }
  } else return kFALSE;
  ULong64_t cint1ac = 0;
  UShort_t nAC=0;
  if (cint1acIndex!=-1) {
     cint1ac=GetDeltaScaler(scalRec1, scalRec2, cint1acIndex, "l0b");
     AliTriggerClass* cint1acClass = cfg->GetTriggerClass(cint1acIndex);
     AliTriggerBCMask* cint1acBCMask = new AliTriggerBCMask();
     if (cint1acClass) cint1acBCMask = (AliTriggerBCMask*)cint1acClass->GetBCMask();
     if (TString(cint1acBCMask->GetName()).CompareTo("NONE")!=0){
        nAC = (UShort_t)cint1acBCMask->GetNUnmaskedBCs();
     }
  }
  ULong64_t cbeamb = 0;
  if (cbeambIndex!=-1) cbeamb=GetDeltaScaler(scalRec1, scalRec2, cbeambIndex, "l0b");
  UInt_t orbits = GetDeltaOrbits(scalRec1, scalRec2);
  //
  if (bkgCorr && (nB==0 || nAC==0 )) return kFALSE;
  //
  if (colBCsFromFillScheme) {
     if (CalculateMu(mu, errmu, cint1b, cint1ac, nB, nAC, orbits, bkgCorr, triggerEff, errorEff)) return kTRUE;
  }
  else {
     if (cint1b!=0 && cbeamb==0) return kFALSE;
     if (CalculateMu(mu, errmu, cint1b, cint1ac, cbeamb, nB, nAC, bkgCorr, triggerEff, errorEff)) return kTRUE;
  }
  //
  return kFALSE;
}
//_____________________________________________________________________________
ULong64_t AliTriggerRunScalers::GetDeltaScalerForRun(Int_t classIndex, TString level)
{
  if (fScalersRecordESD.GetEntriesFast()==0) {
     if (CorrectScalersOverflow()==1) return 0;
  }
  if (level != "l0b" && level != "l0a" && level != "l1b" && level != "l1a" && level != "l2b" && level != "l2a") return 0;

  if (fScalersRecordESD.GetEntriesFast()>1) { 
     const AliTriggerScalersRecordESD* scalRec1 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(0);
     const AliTriggerScalersRecordESD* scalRec2 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(fScalersRecordESD.GetEntriesFast()-1);
     return GetDeltaScaler(scalRec1, scalRec2, classIndex, level);
  }
  return 0;
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::GetScalerRateForRun(Double_t &rate, Double_t &error, Int_t classIndex, TString level)
{
  if (level != "l0b" && level != "l0a" && level != "l1b" && level != "l1a" && level != "l2b" && level != "l2a") return 0;
  if (fScalersRecordESD.GetEntriesFast()==0) {
     if (CorrectScalersOverflow()==1) return kFALSE;
  }
  if (fScalersRecordESD.GetEntriesFast()>=4) {
     const AliTriggerScalersRecordESD* scalRec1 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(0);
     const AliTriggerScalersRecordESD* scalRec2 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(fScalersRecordESD.GetEntriesFast()-1);
     if (GetScalerRate(rate, error, scalRec1, scalRec2, classIndex, level)) return kTRUE;
  }
  return kFALSE;
}
//_____________________________________________________________________________
Bool_t AliTriggerRunScalers::GetClassL2L0ForRun(Double_t &l2l0, Double_t &error, Int_t classIndex)
{
  if (fScalersRecordESD.GetEntriesFast()==0) {
     if (CorrectScalersOverflow()==1) return kFALSE;
  }
  if (fScalersRecordESD.GetEntriesFast()>=4) {
     const AliTriggerScalersRecordESD* scalRec1 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(0);
     const AliTriggerScalersRecordESD* scalRec2 = (AliTriggerScalersRecordESD*)fScalersRecordESD.At(fScalersRecordESD.GetEntriesFast()-1);  // Skip the first and last scaler readings for timing reasons

     if ( GetClassL2L0(l2l0, error, scalRec1, scalRec2, classIndex)) return kTRUE;
  }
  return kFALSE;
}
//_____________________________________________________________________________
TGraphErrors* AliTriggerRunScalers::GetGraphScalerRate(const char* className, TString level, AliTriggerConfiguration* cfg)
{
  Int_t classIndex = cfg->GetClassIndexFromName(className);
  if (classIndex == -1) return 0;
  if (fScalersRecordESD.GetEntriesFast()==0) {
     if (CorrectScalersOverflow()==1) return 0;
  }
  if (level != "l0b" && level != "l0a" && level != "l1b" && level != "l1a" && level != "l2b" && level != "l2a") return 0;
  Int_t nent = fScalersRecordESD.GetEntriesFast();
  Double_t* time = new Double_t[nent];
  Double_t* etime = new Double_t[nent];
  Double_t* rate = new Double_t[nent];
  Double_t* erate = new Double_t[nent];
  for (Int_t i=0;i<nent-1;i++) {
     if (i>0) time[i] = time[i-1]+GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(i-1), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i))/2.+GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(i), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i+1))/2.;
     else time[0] = GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(0), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(1))/2.;
     etime[i] = GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(i), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i+1))/2.;
     if (!GetScalerRate( rate[i], erate[i], (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i+1),  classIndex, level)) {rate[i]=-1; erate[i]=0.;}
  }
  TGraphErrors* graph = new TGraphErrors(nent-1, time, rate, etime, erate);  
  return graph;
}
//_____________________________________________________________________________
TGraphErrors* AliTriggerRunScalers::GetGraphScalerL2L0Ratio(const char* className, AliTriggerConfiguration* cfg)
{
  Int_t classIndex = cfg->GetClassIndexFromName(className);
  if (classIndex == -1) return 0;
  if (fScalersRecordESD.GetEntriesFast()==0) {
     if (CorrectScalersOverflow()==1) return 0;
  }
  Int_t nent = fScalersRecordESD.GetEntriesFast();
  Double_t* time = new Double_t[nent];
  Double_t* etime = new Double_t[nent];
  Double_t* ratio = new Double_t[nent];
  Double_t* eratio = new Double_t[nent];

  for (Int_t i=0;i<nent-1;i++) {
     if (i>0) time[i] = time[i-1]+GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(i-1), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i))/2.+GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(i), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i+1))/2.;
     else time[0] = GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(0), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(1))/2.;
     etime[i] = GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(i), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i+1))/2.;
     if (!GetClassL2L0(ratio[i], eratio[i], (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i+1),  classIndex)) {ratio[i]=-1; eratio[i]=0.;}
  }
  TGraphErrors* graph = new TGraphErrors(nent-1, time, ratio, etime, eratio);  
  return graph;
}
//_____________________________________________________________________________
TGraphErrors* AliTriggerRunScalers::GetGraphMu(AliTriggerConfiguration* cfg, const char* className, Bool_t colBCsFromFillScheme, Bool_t bkgCorr, Double_t triggerEff, Double_t errorEff)
{
  if (fScalersRecordESD.GetEntriesFast()==0) {
     if (CorrectScalersOverflow()==1) return 0;
  }
  Int_t nent = fScalersRecordESD.GetEntriesFast();
  Double_t* time = new Double_t[nent];
  Double_t* etime = new Double_t[nent];
  Double_t* mu = new Double_t[nent];
  Double_t* emu = new Double_t[nent];
  
  for (Int_t i=0;i<nent-1;i++) {
     Double_t m=0.;
     Double_t err=0.;
     if (i!=0) time[i] = time[i-1]+GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(i-1), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i))/2.+GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(i), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i+1))/2.;
     else time[0] = GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(0), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(1))/2.;
     etime[i] = GetDeltaTime((AliTriggerScalersRecordESD*)fScalersRecordESD.At(i), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i+1))/2.;
     if (GetMuFromClassScaler( m,err, className, (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i), (AliTriggerScalersRecordESD*)fScalersRecordESD.At(i+1), cfg, colBCsFromFillScheme, bkgCorr, triggerEff, errorEff)!=kFALSE) {mu[i]=m; emu[i]=err;}
     else {mu[i]=-1.; emu[i]=0.;}
  }
  TGraphErrors* graph = new TGraphErrors(nent-1, time, mu, etime, emu);  
  return graph;
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

