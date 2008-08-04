#ifndef ALIITSGAINSSDV2_H
#define ALIITSGAINSSDV2_H
 
//////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
// Modified: 08/07/2008
//                                          //
//////////////////////////////////////////////
#include "TObject.h"

class AliITSGainSSDv2 : public TObject {

 public:

    AliITSGainSSDv2();
    virtual ~AliITSGainSSDv2();
    AliITSGainSSDv2(const AliITSGainSSDv2 &source); // copy constructor
    AliITSGainSSDv2& operator=(const AliITSGainSSDv2 &source); // ass. op.

    void AddGainP(Int_t module, Int_t strip, Float_t value) { 
      fGain[module*2*fgkDefaultNStripsSSD+strip] = (UShort_t) (1000.*value);
    }       
    Float_t GetGainP(Int_t module, Int_t strip) {
      return ((Float_t) ( (Float_t) fGain[module*2*fgkDefaultNStripsSSD+strip] ) /1000.); 
    }

    void AddGainN(Int_t module, Int_t strip, Float_t value) { 
      fGain[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip] = 
	(UShort_t) (1000.*value);
    }       
    Float_t GetGainN(Int_t module, Int_t strip) {
      return ( (Float_t) ((Float_t) fGain[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip])/1000.); 
    }

 protected:

    static const Int_t fgkDefaultNModulesSSD = 1698;
    static const Int_t fgkDefaultNStripsSSD = 768;

    //static const Int_t fgkDefaultNModulesSSD; // Total numbers of SSD modules
    //static const Int_t fgkDefaultNStripsSSD; // Total numbers of SSD modules

    UShort_t fGain[2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD]; 
    //   UShort_t *fGain; 

 private:
    
    ClassDef(AliITSGainSSDv2,1) // Gain  class for SSD
      };
#endif
