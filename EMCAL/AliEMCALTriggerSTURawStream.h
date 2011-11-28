#ifndef ALIEMCALTRIGGERSTURAWSTREAM_H
#define ALIEMCALTRIGGERSTURAWSTREAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*

 
Author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
*/

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
//#include <map>

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
	virtual Int_t           GetRawData() const {return fGetRawData;}
	
	virtual Bool_t     GetL0GammaPatch(const Int_t i, Int_t& x, Int_t& y) const;
	virtual Bool_t     GetL1GammaPatch(const Int_t i, Int_t& x, Int_t& y, Int_t& z) const;
	virtual Bool_t       GetL1JetPatch(const Int_t i, Int_t& x, Int_t& y) const;
	
	virtual UInt_t              GetV0A()           const {return fV0A;}
	virtual UInt_t              GetV0C()           const {return fV0C;}
	virtual UInt_t              GetGA()            const {return fGA;}
	virtual UInt_t              GetGB()            const {return fGB;}
	virtual UInt_t              GetGC()            const {return fGC;}
	virtual UInt_t              GetJA()            const {return fJA;}
	virtual UInt_t              GetJB()            const {return fJB;}
	virtual UInt_t              GetJC()            const {return fJC;}
	virtual UInt_t              GetRegionEnable()  const {return fRegionEnable;}
	virtual UInt_t              GetFrameReceived() const {return fFrameReceived;}
	virtual UInt_t              GetFwVersion()     const {return fFwVersion;}
	
private:
    
	AliEMCALTriggerSTURawStream(const AliEMCALTriggerSTURawStream& rhs);
    AliEMCALTriggerSTURawStream& operator = (const AliEMCALTriggerSTURawStream& rhs);

    UShort_t GetThreshold(Short_t A, Short_t B, Short_t C, UShort_t V0A, UShort_t V0C);

    AliRawReader* fRawReader;   // object for reading the raw data

	UInt_t              	             fL1JetThreshold;          // L1 Jet Threshold
	UInt_t              	           fL1GammaThreshold;          // L1 Gamma Threshold
	UShort_t                          fL0GammaPatchIndex[3100];    // L0 Gamma Patch Index
	UShort_t                          fL1GammaPatchIndex[3100];    // L1 Gamma Patch Index
	UShort_t                            fL1JetPatchIndex[200];     // L1 Jet Patch Index
	                                                                  
	Int_t                                 fNL0GammaPatch;          // N L0 Gamma Patch
	Int_t                                   fNL1JetPatch;          // N L1 Jet Patch
	Int_t                                 fNL1GammaPatch;          // N L1 Gamma Patch
	                                                                  
	Int_t                                    fGetRawData;          // Get Raw Data
	                                                                  
	UInt_t                                          fADC[32][96];  // ADC
	                                                                  
	UInt_t                                          fV0A;          // V0A
    UInt_t                                          fV0C;          // V0C
    UInt_t                                           fGA;          // GA
    UInt_t                                           fGB;          // GB
    UInt_t                                           fGC;          // GC
    UInt_t                                           fJA;          // JA
    UInt_t                                           fJB;          // JB
    UInt_t                                           fJC;          // JC
    UInt_t                                 fRegionEnable;          // Region Enable
    UInt_t                                fFrameReceived;          // Frame Received
    UInt_t                                    fFwVersion;          // Fw Version

    ClassDef(AliEMCALTriggerSTURawStream,1)   // class for reading EMCAL STU DDL raw data
};

#endif
