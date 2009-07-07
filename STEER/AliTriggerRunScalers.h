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
class AliTimeStamp;
class AliTriggerScalersESD;
class AliTriggerScalersRecord;
class AliTriggerScalersRecordESD;

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
 // Analysis		    
                    Int_t    ConsistencyCheck(Int_t position,Bool_t correctOverflow);
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
                                        
static AliTriggerRunScalers*  ReadScalers( TString & filename );
                                      
private:
                  Short_t    fVersion;            // Version
                  ULong_t    fRunNumber;          // Run number
                  UChar_t    fnClasses;           // Number of trigger classes
                  TArrayC    fClassIndex;         // list of classes used in this partition
                TObjArray    fScalersRecord;      // Array of records (AliTriggerScalersRecord)
                TObjArray    fScalersRecordESD;   // Array of records with 64bit scalers (AliTriggerScalersRecordESD)

    

   ClassDef( AliTriggerRunScalers, 2 )  // Define a Run Trigger Scalers (Scalers)
};

#endif
