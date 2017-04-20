#ifndef AliVWeakResult_H
#define AliVWeakResult_H
#include <TNamed.h>
#include <TList.h>
#include <TH3F.h>
#include <TProfile.h>

//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// TObject to hold weak decay results
// This is a base class for AliV0Result and AliCascadeResult
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class AliVWeakResult : public TNamed {
    
public:
    //Simple constructor
    AliVWeakResult();
    
    //TNamed-inspired constructor
    AliVWeakResult(const char * name, const char * title = "Result");
    
    //Simple destructor
    ~AliVWeakResult();
    
    void Clear(Option_t* = "") {}; //dummy
    
    //Basic Functionality
    virtual Double_t GetMass() const { return 0; }
    virtual TString GetParticleName() const { return ""; }
    //Make sure to have getters available no matter what result 
    virtual Double_t GetCutMinRapidity() const { return -100; }
    virtual Double_t GetCutMaxRapidity() const { return -100; }
    
    //Main Result
    virtual TH3F* GetHistogram       ()       { return 0x0; }
    virtual TH3F* GetHistogramToCopy () const { return 0x0; }
    //Proton Profile
    virtual TProfile* GetProtonProfile       ()       { return 0x0; }
    virtual TProfile* GetProtonProfileToCopy () const { return 0x0; }
    //Feeddown matrix
    virtual TH3F* GetHistogramFeeddown       ()       { return 0x0; }
    virtual TH3F* GetHistogramFeeddownToCopy () const { return 0x0; }
    virtual Bool_t HasSameCuts( AliVWeakResult *lCompare, Bool_t lCheckdEdx = kTRUE ) { return kFALSE; }
    virtual void Print() {};
    
private:
    ClassDef(AliVWeakResult, 3)
    // 1 - original implementation
    // 2 - adjustments for general-purpose functionality
    // 3 - adjustments for proton profile
};
#endif
