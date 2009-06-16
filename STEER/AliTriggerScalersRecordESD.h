#ifndef ALITRIGGERSCALERSRECORDESD_H
#define ALITRIGGERSCALERSRECORDESD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//
// Class to define the ALICE Trigger Scalers Record 
//
// Each record consists of 1 time stamp (4 words)  (AliTimeStamp)
// and an array with the scalers (AliTriggerScalers) for each trigger class 
// in partition  
//
//////////////////////////////////////////////////////////////////////////////
#include <TObjArray.h>
class TObjArray;
class AliTriggerScalersESD;

class AliTriggerScalersRecordESD : public TObject {

public:
                            AliTriggerScalersRecordESD();
                    
                            virtual   ~AliTriggerScalersRecordESD() { fScalers.SetOwner(); fScalers.Delete(); }
                 
                 
                    void    AddTriggerScalers( AliTriggerScalersESD* scaler );
                    void    AddTriggerScalers( UChar_t classIndex, ULong64_t LOCB, ULong64_t LOCA,        
                                              ULong64_t L1CB, ULong64_t L1CA, ULong64_t L2CB, ULong64_t L2CA );
                            
               TObjArray*   GetTriggerScalers()  { return  &fScalers; }
    AliTriggerScalersESD*   GetTriggerScalersForClass( Int_t classindex );       
                                
            virtual void    Print( const Option_t* opt ="" ) const;

               

private:  

                      TObjArray    fScalers;      // Array of scalers (AliTriggerScalersESD) 
                      AliTriggerScalersRecordESD( const AliTriggerScalersRecordESD &rec );
                      AliTriggerScalersRecordESD&   operator=(const AliTriggerScalersRecordESD& rec);      

   ClassDef( AliTriggerScalersRecordESD, 1 )  // Define a Record of Trigger Scalers 
};

#endif
