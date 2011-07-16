#ifndef ALITRIGGERCLASS_H
#define ALITRIGGERCLASS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This class represents the CTP class objects                               //
//                                                                           //
// The Class consists of Name, descriptor, mask, protection, index in the    //
// trigger mask                                                              //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>

class AliTriggerConfiguration;
class AliTriggerDescriptor;
class AliTriggerCluster;
class AliTriggerPFProtection;
class AliTriggerBCMask;

class AliTriggerClass : public TNamed {

public:
                          AliTriggerClass();
                          AliTriggerClass( TString & name, UChar_t index,
					   AliTriggerDescriptor *desc, AliTriggerCluster *clus,
					   AliTriggerPFProtection *pfp, AliTriggerBCMask *mask,
					   UInt_t prescaler, Bool_t allrare);
                          AliTriggerClass( AliTriggerConfiguration *config,
					   TString & name, UChar_t index,
					   TString &desc, TString &clus,
					   TString &pfp, TString &mask,
					   UInt_t prescaler, Bool_t allrare);

                          AliTriggerClass( const AliTriggerClass& trclass );
               virtual   ~AliTriggerClass();
  AliTriggerClass&   operator=(const AliTriggerClass& trclass);

                  void    Reset() { fStatus = kFALSE; }

             ULong64_t    GetValue() const { return (fStatus) ? fClassMask : 0; }
                Bool_t    GetStatus() const { return fStatus; }
               ULong64_t  GetMask() const { return fClassMask; }
    AliTriggerDescriptor* GetDescriptor() const { return fDescriptor; }
       AliTriggerCluster* GetCluster() const { return fCluster; }
        AliTriggerBCMask* GetBCMask() const { return fMask; }

                    void  Trigger( const TObjArray& inputs , const TObjArray& functions);
		    void  Print( const Option_t* ) const;

                  Bool_t  CheckClass(AliTriggerConfiguration *config) const;
		  Bool_t  IsActive( const TObjArray& inputs, const TObjArray& functions) const;
private:
	       ULong64_t  fClassMask;    // trigger mask (1<< (index-1))
	       	 UChar_t  fIndex;        // position of class in mask
    AliTriggerDescriptor* fDescriptor;   // pointer to the descriptor
       AliTriggerCluster* fCluster;      // pointer to the cluster
  AliTriggerPFProtection* fPFProtection; // pointer to the past-future protection
        AliTriggerBCMask* fMask;         // pointer to bunch-crossing mask
                  UInt_t  fPrescaler;    // Downscaling factor
                  Bool_t  fAllRare;      // All or Rare trigger
		  Bool_t  fStatus;       //! true = Condition has been satisfied after Trigger

  ClassDef( AliTriggerClass, 3 )  // Define a trigger class object
};

#endif
