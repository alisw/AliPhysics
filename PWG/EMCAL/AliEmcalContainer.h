#ifndef AliEmcalContainer_H
#define AliEmcalContainer_H

//
// container with name, TClonesArray
//

class TLorentzVector;
class AliVEvent;
class AliNamedArrayI;

#include <TNamed.h>
#include <TClonesArray.h>

class AliEmcalContainer : public TNamed {
 public:
  enum RejectionReason {
    // General
    kNullObject = 1<<0,
    kPtCut = 1<<1,
    kAcceptanceCut = 1<<2,
    kBitMapCut = 1<<3,
    // leave bits 4-7 free for future implementations
    
    // AliParticleContainer
    kMCFlag = 1<<8,
    kMCGeneratorCut = 1<<9,
    kChargeCut = 1<<10,
    // leave bits 11-12 free for future implementations

    // AliClusterContainer
    kIsEMCalCut = 1<<13,
    kTimeCut = 1<<14,
    kEnergyCut = 1<<15,
    // leave bits 16-17 free for future implementations

    // AliJetContainer
    kAreaCut = 1<<18,
    kAreaEmcCut = 1<<19,
    kZLeadingChCut = 1<<20,
    kZLeadingEmcCut = 1<<21,
    kNEFCut = 1<<22,
    kMinLeadPtCut = 1<<23,
    kMaxTrackPtCut = 1<<24,
    kMaxClusterPtCut = 1<<25,
    kFlavourCut = 1<<26,
    kTagStatus = 1<<27
  };

  AliEmcalContainer();
  AliEmcalContainer(const char *name); 
  virtual ~AliEmcalContainer(){;}

  TClonesArray               *GetArray() const                      { return fClArray                   ; }
  const TString&              GetArrayName()                  const { return fClArrayName               ; }
  Int_t                       GetCurrentID()                  const { return fCurrentID-1               ; }
  Bool_t                      GetIsParticleLevel()            const { return fIsParticleLevel           ; }
  Int_t                       GetIndexFromLabel(Int_t lab)    const;
  Int_t                       GetNEntries()                   const { return fClArray->GetEntriesFast() ; }
  virtual void                GetMomentum(TLorentzVector &mom, Int_t i) const = 0;
  void                        ResetCurrentID(Int_t i=0)             { fCurrentID = i                    ; }
  virtual void                SetArray(AliVEvent *event);
  void                        SetArrayName(const char *n)           { fClArrayName = n                  ; }
  void                        SetIsParticleLevel(Bool_t b)          { fIsParticleLevel = b              ; }
  void                        SortArray()                           { fClArray->Sort()                  ; }
  UInt_t                      GetRejectionReason()            const { return fRejectionReason           ; }
  UInt_t                      TestRejectionReason(UInt_t rs)  const { return fRejectionReason & rs      ; }
  UShort_t                    GetRejectionReasonBitPosition() const;

 protected:
  TString                     fClArrayName;             // name of branch
  TString                     fClassName;               // name of the class in the TClonesArray
  Bool_t                      fIsParticleLevel;         // whether or not it is a particle level object collection
  TClonesArray               *fClArray;                 //!TClonesArray
  Int_t                       fCurrentID;               //!current ID for automatic loops
  AliNamedArrayI             *fLabelMap;                //!Label-Index map
  Double_t                    fVertex[3];               //!event vertex array
  UInt_t                      fRejectionReason;         //!reject reason bit map for the last call to an accept object function

 private:
  AliEmcalContainer(const AliEmcalContainer& obj); // copy constructor
  AliEmcalContainer& operator=(const AliEmcalContainer& other); // assignment

  ClassDef(AliEmcalContainer,3);
};
#endif
