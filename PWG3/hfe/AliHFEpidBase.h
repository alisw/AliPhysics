/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// Base Class for Detector PID Objects
// For more information see the implementation file
//
#ifndef ALIHFEPIDBASE_H
#define ALIHFEPIDBASE_H
 
 #ifndef ROOT_TNamed
 #include <TNamed.h>
 #endif

class TList;
class AliVParticle;
class AliMCParticle;

struct AliHFEpidObject{
    typedef enum{ 
      kESDanalysis,
      kAODanalysis
    }AnalysisType_t;
    AliVParticle *fRecTrack;    // Reconstructed track
    AliVParticle *fMCtrack;     // Monte Carlo track
    UChar_t fAnalysisType;      // Analysis Mode (ESD or AOD)
    AliHFEpidObject():fRecTrack(NULL), fMCtrack(NULL), fAnalysisType(kESDanalysis){}
};

class AliHFEpidBase : public TNamed{
  public:
    AliHFEpidBase(const Char_t *name);
    AliHFEpidBase(const AliHFEpidBase &c);
    AliHFEpidBase &operator=(const AliHFEpidBase &c);
    virtual ~AliHFEpidBase() {};
    // Framework functions that have to be implemented by the detector PID classes
    virtual Bool_t InitializePID() = 0;
    virtual Int_t IsSelected(AliHFEpidObject *track) = 0;
    virtual Bool_t HasQAhistos() const = 0;

    Int_t GetDebugLevel() const { return fDebugLevel; };
    Bool_t IsQAon() const { return TestBit(kQAon);};
    Bool_t HasMCData() const { return TestBit(kHasMCData); };

    void SetDebugLevel(Int_t debugLevel) { fDebugLevel = debugLevel; }; 
    inline void SetQAOn(TList *fQAlist);
    void SetHasMCData(Bool_t hasMCdata = kTRUE) { SetBit(kHasMCData,hasMCdata); };

  protected:
    void Copy(TObject &ref) const;
    virtual void AddQAhistograms(TList *){};
  private:
    enum{
      kQAon = BIT(14),
      kHasMCData = BIT(15)
    };

    Int_t fDebugLevel;              // Debug Level

    ClassDef(AliHFEpidBase, 1)      // Base class for detector Electron ID
};

//___________________________________________________________________
void AliHFEpidBase::SetQAOn(TList *qaList){
  //
  // Initialize QA for Detector PID class
  //
  if(HasQAhistos()){
    SetBit(kQAon, kTRUE);
    AddQAhistograms(qaList);
  }
}
#endif
