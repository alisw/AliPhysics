#ifndef ALIITSBADCHANNELSSSD_H
#define ALIITSBADCHANNELSSSD_H
 
#include "TObjArray.h"
#include "TArrayI.h"

//////////////////////////////////////////////
// Response class for SSD                   //
//                                          //
//////////////////////////////////////////////
class AliITSBadChannelsSSD : public TObject {

 public:
    AliITSBadChannelsSSD();
    virtual ~AliITSBadChannelsSSD();

    void SetNBadPChannelsList(Int_t n) { fBadPChannelsList.Set(n); }
    void AddBadPChannel(Int_t c, Int_t n) { fBadPChannelsList.AddAt(n,c);}
    TArrayI GetBadPChannelsList() const {return fBadPChannelsList; }
    void SetNBadNChannelsList(Int_t n) { fBadNChannelsList.Set(n); }
    void AddBadNChannel(Int_t c, Int_t n) { fBadNChannelsList.AddAt(n,c);}
    TArrayI GetBadNChannelsList() const {return fBadNChannelsList; }
    //

    void SetMod(UShort_t mod) {fMod = mod;}
    UShort_t GetMod() { return fMod;}

 protected:

    UShort_t fMod;       // module number (from 0 to 1535). Needed by the preprocessor to 
    //   extract the information from the Detector Algorithm
    
    TArrayI  fBadNChannelsList;  // list of P side dead channels
    TArrayI  fBadPChannelsList;  // list of N side dead channels
    
 private:
    AliITSBadChannelsSSD(const AliITSBadChannelsSSD &source); // copy constructor
    AliITSBadChannelsSSD& operator=(const AliITSBadChannelsSSD &source); // ass. op.

    ClassDef(AliITSBadChannelsSSD,1) //Response class for SSD
};
#endif
