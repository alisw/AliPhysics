#ifndef ALIITSPEDESTALSSD_H
#define ALIITSPEDESTALSSD_H
 
//////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
// Modified: 08/07/2008
//                                          //
//////////////////////////////////////////////
#include "TObject.h"

class AliITSPedestalSSD : public TObject {

 public:

    AliITSPedestalSSD();
    virtual ~AliITSPedestalSSD();
    AliITSPedestalSSD(const AliITSPedestalSSD &source); // copy constructor
    AliITSPedestalSSD& operator=(const AliITSPedestalSSD &source); // ass. op.

    void AddPedestalP(Int_t module, Int_t strip, Float_t value) { 
      fPedestal[module*2*fgkDefaultNStripsSSD+strip] = value;
    }       
    Float_t GetPedestalP(Int_t module, Int_t strip) {
      return fPedestal[module*2*fgkDefaultNStripsSSD+strip]; 
    }

    void AddPedestalN(Int_t module, Int_t strip, Float_t value) { 
      fPedestal[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip] = value;
    }       
    Float_t GetPedestalN(Int_t module, Int_t strip) {
      return fPedestal[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip]; 
    }

 protected:

    static const Int_t fgkDefaultNModulesSSD = 1698;
    static const Int_t fgkDefaultNStripsSSD = 768;

    //   static const Int_t fgkDefaultNModulesSSD; // Total numbers of SSD modules
    //static const Int_t fgkDefaultNStripsSSD; // Total numbers of SSD modules

Float_t fPedestal[2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD]; 

 private:
    
    ClassDef(AliITSPedestalSSD,2) // Pedestal  class for SSD
      };
#endif
