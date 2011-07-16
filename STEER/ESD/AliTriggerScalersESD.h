#ifndef ALITRIGGERSCALERSESD_H
#define ALITRIGGERSCALERSESD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

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

class AliTriggerScalersESD : public TObject {

public:
                         AliTriggerScalersESD();
                         AliTriggerScalersESD(
                                 UChar_t    classIndex,                                 
                               ULong64_t    LOCB,        
                               ULong64_t    LOCA,        
                               ULong64_t    L1CB,        
                               ULong64_t    L1CA,        
                               ULong64_t    L2CB,        
                               ULong64_t    L2CA     
                         );   
			 AliTriggerScalersESD(UChar_t xlassIndex,ULong64_t* s64);
              virtual   ~AliTriggerScalersESD() {}
         virtual void    Print( const Option_t* opt ="" ) const;

          AliTriggerScalersESD( const AliTriggerScalersESD &scal );
          AliTriggerScalersESD&   operator=(const AliTriggerScalersESD& scal);

               ULong64_t    GetLOCB() const { return fLOCB; }
               ULong64_t    GetLOCA() const { return fLOCA; }
               ULong64_t    GetL1CB() const { return fL1CB; }
               ULong64_t    GetL1CA() const { return fL1CA; }
               ULong64_t    GetL2CB() const { return fL2CB; }
               ULong64_t    GetL2CA() const { return fL2CA; }
 	            void    GetAllScalers(ULong64_t *scalers) const;

                 UChar_t    GetClassIndex() const { return fClassIndex; }
    
private:    
                 UChar_t    fClassIndex;      //  number of triggered classes        
               ULong64_t    fLOCB;            //  L0 triggers before any vetos  (64 bits)
               ULong64_t    fLOCA;            //  L0 triggers after all vetos   (64 bits)
               ULong64_t    fL1CB;            //  L1 triggers before any vetos  (64 bits)
               ULong64_t    fL1CA;            //  L1 triggers after all vetos   (64 bits)
               ULong64_t    fL2CB;            //  L2 triggers before any vetos  (64 bits)
               ULong64_t    fL2CA;            //  L2 triggers after all vetos   (64 bits)
                        
   ClassDef( AliTriggerScalersESD, 1 )  // Define a Run Trigger Scalers (Scalers)
};

#endif
