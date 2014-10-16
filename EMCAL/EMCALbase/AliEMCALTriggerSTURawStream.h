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
	
	virtual void DecodeL1JetPatchIndexes(  const int i, UInt_t *word32, const int offset);
	virtual void DecodeL1GammaPatchIndexes(const int i, UInt_t *word32, const int offset);
	virtual void DecodeL0GammaPatchIndexes(             UInt_t *word32, const int offset);
    virtual void DecodeTRUADC(                          UInt_t *word32, const int offset);
	
	virtual void   DumpPayLoad(const Option_t *option = "ALL") const;

	virtual void                GetADC(Int_t iTRU, UInt_t ADC[]);
	virtual UInt_t   GetL1JetThreshold(const int i) const {return   fL1JetThreshold[i];}
	virtual UInt_t GetL1GammaThreshold(const int i) const {return fL1GammaThreshold[i];}
	
	virtual Int_t     GetNL0GammaPatch() const {return fNL0GammaPatch;}
	virtual Int_t     GetNL1GammaPatch(const int i) const {return fNL1GammaPatch[i];}
	virtual Int_t       GetNL1JetPatch(const int i) const {return fNL1JetPatch[i];}
	virtual Int_t           GetRawData() const {return fGetRawData;}
	
	virtual Bool_t     GetL0GammaPatch(const Int_t i, Int_t& x, Int_t& y) const;
	virtual Bool_t     GetL1GammaPatch(const Int_t i, const Int_t j, Int_t& x, Int_t& y, Int_t& z) const;
	virtual Bool_t       GetL1JetPatch(const Int_t i, const Int_t j, Int_t& x, Int_t& y) const;
	
	virtual UInt_t              GetV0A()           const {return fV0A;}
	virtual UInt_t              GetV0C()           const {return fV0C;}
	virtual UInt_t              GetG(int i, int j) const {return fG[i][j];}
	virtual UInt_t              GetJ(int i, int j) const {return fJ[i][j];}
	virtual UInt_t              GetRegionEnable()  const {return fRegionEnable;}
	virtual UInt_t              GetFrameReceived() const {return fFrameReceived;}
	virtual UInt_t              GetFwVersion()     const {return fFwVersion;}
	
private:
    
	AliEMCALTriggerSTURawStream(const AliEMCALTriggerSTURawStream& rhs);
    AliEMCALTriggerSTURawStream& operator = (const AliEMCALTriggerSTURawStream& rhs);

    UShort_t GetThreshold(Short_t A, Short_t B, Short_t C, UShort_t V0A, UShort_t V0C) const;

    AliRawReader* fRawReader;   // object for reading the raw data

	UInt_t              	             fL1JetThreshold[2];          // L1 Jet Threshold
	UInt_t              	           fL1GammaThreshold[2];          // L1 Gamma Threshold
	UShort_t                          fL0GammaPatchIndex[3100];       // L0 Gamma Patch Index
	UShort_t                          fL1GammaPatchIndex[3100][2];    // L1 Gamma Patch Index
	UShort_t                            fL1JetPatchIndex[200][2];     // L1 Jet Patch Index
	                                                                  
	Int_t                                 fNL0GammaPatch;             // N L0 Gamma Patch
	Int_t                                   fNL1JetPatch[2];          // N L1 Jet Patch
	Int_t                                 fNL1GammaPatch[2];          // N L1 Gamma Patch
	                                                                  
	Int_t                                    fGetRawData;          // Get Raw Data
	                                                                  
	UInt_t                                          fADC[32][96];  // ADC
	                                                                  
	UInt_t                                          fV0A;          // V0A
    UInt_t                                          fV0C;          // V0C
	UInt_t                                      fG[3][2];          // Gamma
	UInt_t                                      fJ[3][2];          // Jet
    UInt_t                                 fRegionEnable;          // Region Enable
    UInt_t                                fFrameReceived;          // Frame Received
    UInt_t                                    fFwVersion;          // Fw Version

    ClassDef(AliEMCALTriggerSTURawStream,2)   // class for reading EMCAL STU DDL raw data
};

#endif
