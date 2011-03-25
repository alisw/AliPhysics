#ifndef ALIEMCALTRIGGERSTURAWSTREAM_H
#define ALIEMCALTRIGGERSTURAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*

 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#include <TObject.h>
#include <map>

class AliRawReader;

class AliEMCALTriggerSTURawStream: public TObject 
{
  public:
	         AliEMCALTriggerSTURawStream();
             AliEMCALTriggerSTURawStream(AliRawReader* rawReader);
    virtual ~AliEMCALTriggerSTURawStream();
  
    virtual void   Reset();
    virtual Bool_t ReadPayLoad();
	virtual void   DumpPayLoad(const Option_t *option = "ALL") const;

	virtual void                GetADC(Int_t iTRU, UInt_t ADC[]);
	virtual UInt_t   GetL1JetThreshold() const {return   fL1JetThreshold;}
	virtual UInt_t GetL1GammaThreshold() const {return fL1GammaThreshold;}
	
	virtual Int_t     GetNL0GammaPatch() const {return fNL0GammaPatch;}
	virtual Int_t     GetNL1GammaPatch() const {return fNL1GammaPatch;}
	virtual Int_t       GetNL1JetPatch() const {return fNL1JetPatch;}
	virtual Int_t          GetRawData() const {return fGetRawData;}
	
	virtual Bool_t     GetL0GammaPatch(const Int_t i, Int_t& x, Int_t& y) const;
	virtual Bool_t     GetL1GammaPatch(const Int_t i, Int_t& x, Int_t& y, Int_t& z) const;
	virtual Bool_t       GetL1JetPatch(const Int_t i, Int_t& x, Int_t& y) const;
		
private:
    
	AliEMCALTriggerSTURawStream(const AliEMCALTriggerSTURawStream& rhs);
    AliEMCALTriggerSTURawStream& operator = (const AliEMCALTriggerSTURawStream& rhs);

    AliRawReader* fRawReader;   // object for reading the raw data

	UInt_t              	             fL1JetThreshold;          //
	UInt_t              	           fL1GammaThreshold;          //
	UShort_t                          fL0GammaPatchIndex[3100];    //
	UShort_t                          fL1GammaPatchIndex[3100];    //
	UShort_t                            fL1JetPatchIndex[200];     //
	
	Int_t                                 fNL0GammaPatch;          //
	Int_t                                   fNL1JetPatch;          //
	Int_t                                 fNL1GammaPatch;          //
	
	Int_t                                     fGetRawData;         //
	
	UInt_t                                          fADC[32][96];  //
	
    ClassDef(AliEMCALTriggerSTURawStream,1)   // class for reading EMCAL STU DDL raw data
};

#endif
