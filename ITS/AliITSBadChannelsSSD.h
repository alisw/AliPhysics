#ifndef ALIITSBADCHANNELSSSD_H
#define ALIITSBADCHANNELSSSD_H
 
//////////////////////////////////////////////
// Author: Enrico Fragiacomo
// Date: 23/08/2007
// Modified: 08/07/2008
//                                          //
//////////////////////////////////////////////
#include "TObject.h"

class AliITSBadChannelsSSD : public TObject {

 public:

    AliITSBadChannelsSSD();
    virtual ~AliITSBadChannelsSSD();
    AliITSBadChannelsSSD(const AliITSBadChannelsSSD &source); // copy constructor
    AliITSBadChannelsSSD& operator=(const AliITSBadChannelsSSD &source); // ass. op.

    void AddBadChannelP(Int_t module, Int_t strip, Char_t value) { 
      fBadChannels[module*2*fgkDefaultNStripsSSD+strip] = value;
    }       
    Char_t GetBadChannelP(Int_t module, Int_t strip) {
      return fBadChannels[module*2*fgkDefaultNStripsSSD+strip]; 
    }

    void AddBadChannelN(Int_t module, Int_t strip, Char_t value) { 
      fBadChannels[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip] = value;
    }       
    Char_t GetBadChannelN(Int_t module, Int_t strip) {
      return fBadChannels[module*2*fgkDefaultNStripsSSD+fgkDefaultNStripsSSD+strip]; 
    }

 protected:

    static const Int_t fgkDefaultNModulesSSD = 1698;
    static const Int_t fgkDefaultNStripsSSD = 768;
    //static const Int_t fgkDefaultNModulesSSD; // Total numbers of SSD modules
   //   static const Int_t fgkDefaultNStripsSSD; // Total numbers of SSD modules
   Char_t fBadChannels[2*fgkDefaultNModulesSSD*fgkDefaultNStripsSSD];
    //   Char_t *fBadChannels; 

 private:
   
   ClassDef(AliITSBadChannelsSSD,2) // BadChannels  class for SSD
     };

#endif
