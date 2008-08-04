#ifndef ALIITSNOISESSDV2_H
#define ALIITSNOISESSDV2_H
 
//////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
// Modified: 08/07/2008
//                                          //
//////////////////////////////////////////////
#include "TObject.h"

class AliITSNoiseSSDv2 : public TObject {

 public:

    AliITSNoiseSSDv2();
    virtual ~AliITSNoiseSSDv2();
    AliITSNoiseSSDv2(const AliITSNoiseSSDv2 &source); // copy constructor
    AliITSNoiseSSDv2& operator=(const AliITSNoiseSSDv2 &source); // ass. op.

    void AddNoiseP(Int_t module, Int_t strip, Float_t value) { 
      fNois[module*2*fgkDefaultNStripsSSD+strip] = value;
    }       
    Float_t GetNoiseP(Int_t module, Int_t strip) {
      return fNois[module*2*fgkDefaultNStripsSSD+strip]; 
    }

    void AddNoiseN(Int_t module, Int_t strip, Float_t value) { 
      fNois[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip] = value;
    }       
    Float_t GetNoiseN(Int_t module, Int_t strip) {
      return fNois[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip]; 
    }

 protected:

    static const Int_t fgkDefaultNModulesSSD = 1698;
    static const Int_t fgkDefaultNStripsSSD = 768;

    //   static const Int_t fgkDefaultNModulesSSD; // Total numbers of SSD modules
    //static const Int_t fgkDefaultNStripsSSD; // Total numbers of SSD modules

     Float_t fNois[2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD]; 
     //Float_t *fNois;

 private:
    
    ClassDef(AliITSNoiseSSDv2,1) // Noise  class for SSD
      };
#endif
