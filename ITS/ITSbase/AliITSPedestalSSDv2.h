#ifndef ALIITSPEDESTALSSDV2_H
#define ALIITSPEDESTALSSDV2_H
 
//////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
// Modified: 08/07/2008
//                                          //
//////////////////////////////////////////////
#include "TObject.h"

class AliITSPedestalSSDv2 : public TObject {

 public:

    AliITSPedestalSSDv2();
    virtual ~AliITSPedestalSSDv2();
    AliITSPedestalSSDv2(const AliITSPedestalSSDv2 &source); // copy constructor
    AliITSPedestalSSDv2& operator=(const AliITSPedestalSSDv2 &source); // ass. op.

    void AddPedestalP(Int_t module, Int_t strip, Float_t value) { 
      fPedestal[module*2*fgkDefaultNStripsSSD+strip] = value;
    }       
    Float_t GetPedestalP(Int_t module, Int_t strip) const {
      return fPedestal[module*2*fgkDefaultNStripsSSD+strip]; 
    }

    void AddPedestalN(Int_t module, Int_t strip, Float_t value) { 
      fPedestal[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip] = value;
    }       
    Float_t GetPedestalN(Int_t module, Int_t strip) const {
      return fPedestal[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip]; 
    }

 protected:

    static const Int_t fgkDefaultNModulesSSD = 1698;// # of default SSD modules
    static const Int_t fgkDefaultNStripsSSD = 768;// # of default SSD strips

    //   static const Int_t fgkDefaultNModulesSSD; // Total numbers of SSD modules
    //static const Int_t fgkDefaultNStripsSSD; // Total numbers of SSD modules

Float_t fPedestal[2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD];  // tot # of default SSD strips

 private:
    
    ClassDef(AliITSPedestalSSDv2,1) // Pedestal  class for SSD
      };
#endif
