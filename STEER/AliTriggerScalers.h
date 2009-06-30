#ifndef ALITRIGGERSCALERS_H
#define ALITRIGGERSCALERS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTriggerScalers.h 22322 2007-11-22 11:43:14Z cvetan $ */

///////////////////////////////////////////////////////////////////////////////
//
//  Class to define the ALICE Trigger Scalers  
//
//  For each trigger class there are six scalers:
//
//    LOCB       L0 triggers before any vetos 
//    LOCA       L0 triggers after all vetos 
//    L1CB       L1 triggers before any vetos 
//    L1CA       L1 triggers after all vetos 
//    L2CB       L2 triggers before any vetos 
//    L2CA       L2 triggers after all vetos 
//
//////////////////////////////////////////////////////////////////////////////

class AliTriggerScalers : public TObject {

public:
                         AliTriggerScalers();
                         AliTriggerScalers(
                              UChar_t    classIndex, 
                               UInt_t    LOCB,        
                               UInt_t    LOCA,        
                               UInt_t    L1CB,        
                               UInt_t    L1CA,        
                               UInt_t    L2CB,        
                               UInt_t    L2CA       
                         );   
              virtual   ~AliTriggerScalers() {}
              
       virtual Bool_t    IsSortable() const { return kTRUE; }
        virtual Int_t    Compare( const TObject* obj ) const;
         virtual void    Print( const Option_t* opt ="" ) const;
               UInt_t    GetLOCB() const { return fLOCB; }
               UInt_t    GetLOCA() const { return fLOCA; }
               UInt_t    GetL1CB() const { return fL1CB; }
               UInt_t    GetL1CA() const { return fL1CA; }
               UInt_t    GetL2CB() const { return fL2CB; }
               UInt_t    GetL2CA() const { return fL2CA; }
	         void    GetAllScalers(UInt_t *scalers) const;
              UChar_t    GetClassIndex() const { return fClassIndex; }
private: 
   
              UChar_t    fClassIndex;      // class index 
               UInt_t    fLOCB;            //  L0 triggers before any vetos  (32 bits)
               UInt_t    fLOCA;            //  L0 triggers after all vetos   (32 bits)
               UInt_t    fL1CB;            //  L1 triggers before any vetos  (32 bits)
               UInt_t    fL1CA;            //  L1 triggers after all vetos   (32 bits)
               UInt_t    fL2CB;            //  L2 triggers before any vetos  (32 bits)
               UInt_t    fL2CA;            //  L2 triggers after all vetos   (32 bits)


                         AliTriggerScalers( const AliTriggerScalers &run );
    AliTriggerScalers&   operator=(const AliTriggerScalers& clus);

   ClassDef( AliTriggerScalers, 1 )  // Define a Run Trigger Scalers (Scalers)
};

#endif
