#ifndef ALITRIGGERRUNSCALERS_H
#define ALITRIGGERRUNSCALERS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTriggerRunScalers.h 22322 2007-11-22 11:43:14Z cvetan $ */

///////////////////////////////////////////////////////////////////////////////
//
//  Class to define a collection scalers per Run  
//
// 
//
//////////////////////////////////////////////////////////////////////////////
class TObject;
class TGraphErrors;
class AliTimeStamp;
class AliTriggerScalersESD;
class AliTriggerScalersRecord;
class AliTriggerScalersRecordESD;
class AliTriggerConfiguration;

#include "TArrayC.h"

class AliTriggerRunScalers : public TObject {

public:
                         AliTriggerRunScalers();
              virtual   ~AliTriggerRunScalers();
  //  Getters
                  Short_t    GetVersion()          const { return fVersion;       }            
                  ULong_t    GetRunNumber()        const { return fRunNumber;     }
                  UChar_t    GetNumClasses()       const { return fnClasses;      }
                   Char_t    GetClass( Int_t i )   const { return fClassIndex[i]; }
          const TObjArray*   GetScalersRecords()   const { return &fScalersRecord; } 
          const TObjArray*   GetScalersRecordsESD()   const { return &fScalersRecordESD; } 
  AliTriggerScalersRecord*   GetScalersRecord( Int_t index ) const { return (AliTriggerScalersRecord*)fScalersRecord.At(index); }
                    Int_t    FindNearestScalersRecord( const AliTimeStamp *stamp ) const;
     AliTriggerScalersESD*   GetScalersForEventClass(const AliTimeStamp* stamp,const Int_t classIndex) const;
     const AliTriggerScalersRecordESD*   GetScalersDeltaForEvent(const AliTimeStamp* stamp) const;
     const AliTriggerScalersRecordESD*   GetScalersDeltaForRun() const;

 // Analysis		    
                    Int_t    ConsistencyCheck(Int_t position,Bool_t correctOverflow, UInt_t** overflow);
		    Int_t    CorrectScalersOverflow();
		    Int_t    CheckRunScalers(){return (fScalersRecord.GetEntriesFast()==fScalersRecordESD.GetEntriesFast());}
  //  Setters
                     void    SetVersion( Short_t ver )       { fVersion = ver;   }            
                     void    SetRunNumber( ULong_t run )     { fRunNumber = run; }
                     void    SetNumClasses( UChar_t nclass ) { fnClasses = nclass; fClassIndex.Set(nclass); }
                     void    SetClass( UChar_t i, UChar_t index ) { fClassIndex[i]=index; }
                     void    AddTriggerScalers( AliTriggerScalersRecord* scal );
             virtual void    Print( const Option_t* opt ="" ) const;
     AliTriggerRunScalers( const AliTriggerRunScalers &run );
     AliTriggerRunScalers&    operator=(const AliTriggerRunScalers& run);
                                        
static AliTriggerRunScalers* ReadScalers( TString & filename );
	  static Bool_t    CalculateMu(Double_t &mu, Double_t &errmu, ULong64_t countsB, ULong64_t countsAC, UShort_t nB, UShort_t nAC, UInt_t orbits, Bool_t bkgCorr=kTRUE, Double_t triggerEff=1., Double_t errorEff=0.);
	  static Bool_t    CalculateMu(Double_t &mu, Double_t &errmu, ULong64_t countsB, ULong64_t countsAC, ULong64_t beamB, UShort_t nB, UShort_t nAC, Bool_t bkgCorr=kTRUE, Double_t triggerEff=1., Double_t errorEff=0.);
	  static ULong64_t    GetDeltaScaler(const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, Int_t classIndex, TString level);
	  static Double_t    GetDeltaTime(const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2);
	  static UInt_t    GetDeltaOrbits(const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2);
	  static Bool_t    GetScalerRate(Double_t &rate, Double_t &error, const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, Int_t classIndex, TString level);
	  static Bool_t    GetScalerRatePerBC(Double_t &rate, Double_t &error, const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, AliTriggerConfiguration* cfg, Int_t classIndex, TString level);
	  static Bool_t    GetClassL2L0(Double_t &l2l0, Double_t &error, const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, Int_t classIndex);
          static Bool_t    GetMuFromClassScaler(Double_t &mu, Double_t &errmu, const char* className, const AliTriggerScalersRecordESD* scalRec1, const AliTriggerScalersRecordESD* scalRec2, AliTriggerConfiguration* cfg, Bool_t colBCsFromFillScheme=kTRUE, Bool_t bkgCorr=kTRUE, Double_t triggerEff=1., Double_t errorEff=0.);
	        ULong64_t    GetDeltaScalerForRun(Int_t classIndex, TString level);
	         Bool_t    GetScalerRateForRun(Double_t &rate, Double_t &error, Int_t classIndex, TString level);
	         Bool_t    GetClassL2L0ForRun(Double_t &l2l0, Double_t &error, Int_t classIndex);
             TGraphErrors*   GetGraphScalerRate(const char* className, TString level, AliTriggerConfiguration* cfg);
             TGraphErrors*   GetGraphScalerL2L0Ratio(const char* className, AliTriggerConfiguration* cfg);
             TGraphErrors*   GetGraphMu(AliTriggerConfiguration* cfg, const char* className, Bool_t colBCsFromFillScheme=kTRUE, Bool_t bkgCorr=kTRUE, Double_t triggerEff=1., Double_t errorEff=0.);

private:
                  Short_t    fVersion;            // Version
                  ULong_t    fRunNumber;          // Run number
                  UChar_t    fnClasses;           // Number of trigger classes
                  TArrayC    fClassIndex;         // list of classes used in this partition
                TObjArray    fScalersRecord;      // Array of records (AliTriggerScalersRecord)
                TObjArray    fScalersRecordESD;   // Array of records with 64bit scalers (AliTriggerScalersRecordESD)

    

   ClassDef( AliTriggerRunScalers, 4 )  // Define a Run Trigger Scalers (Scalers)
};

#endif
