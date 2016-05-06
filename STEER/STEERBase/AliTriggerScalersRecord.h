#ifndef ALITRIGGERSCALERSRECORD_H
#define ALITRIGGERSCALERSRECORD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTriggerScalersRecord.h 22322 2007-11-22 11:43:14Z cvetan $ */

///////////////////////////////////////////////////////////////////////////////
//
// Class to define the ALICE Trigger Scalers Record 
//
// Each record consists of 1 time stamp (4 words)  (AliTimeStamp)
// and an array with the scalers (AliTriggerScalers) for each trigger class 
// in partition  
//
//////////////////////////////////////////////////////////////////////////////
#include "AliTimeStamp.h"

class TObjArray;
class AliTriggerScalers;

class AliTriggerScalersRecord : public TObject {

public:
                            AliTriggerScalersRecord();
                 virtual   ~AliTriggerScalersRecord() { fScalers.SetOwner(); fScalers.Delete(); }
                 
                 
                    void    SetTimeStamp( UInt_t orbit, UInt_t period, UInt_t seconds, UInt_t microsecs );
		    void    SetTimeGroup(UInt_t tgr){fTimeGroup=tgr;};
                    void    AddTriggerScalers( AliTriggerScalers* scaler );
                    void    AddTriggerScalers( UChar_t classIndex, UInt_t LOCB, UInt_t LOCA,        
                                              UInt_t L1CB, UInt_t L1CA, UInt_t L2CB, UInt_t L2CA );
                    void    AddTriggerScalers( UChar_t classIndex, UInt_t LOCB, UInt_t LOCA,        
                                              UInt_t L1CB, UInt_t L1CA, UInt_t L2CB, UInt_t L2CA,
					      UInt_t LMCB, UInt_t LMCA);
                            
      const AliTimeStamp*   GetTimeStamp() const { return &fTimestamp; }
         const TObjArray*   GetTriggerScalers()  const { return  &fScalers; }
 const AliTriggerScalers*   GetTriggerScalersForClass( const Int_t classindex ) const;       
       AliTriggerScalers*   GetTriggerScalersForClassBinary( const Int_t classindex ) ;     
                   UInt_t   GetTimeGroup(){return fTimeGroup;}
          virtual Bool_t    IsSortable() const { return kTRUE; }
                                
           virtual Int_t    Compare( const TObject* obj ) const;
            virtual void    Print( const Option_t* opt ="" ) const;

       
             
     
private:  

            AliTimeStamp    fTimestamp;    // record timestamp
               TObjArray    fScalers;      // Array of scalers (AliTriggerScalers) 
	          UInt_t    fTimeGroup;    // Time group of record


                            AliTriggerScalersRecord( const AliTriggerScalersRecord &rec );
 AliTriggerScalersRecord&   operator=(const AliTriggerScalersRecord& rec);

   ClassDef( AliTriggerScalersRecord, 2 )  // Define a Record of Trigger Scalers 
};

#endif
