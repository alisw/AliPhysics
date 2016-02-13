#ifndef AliEmcalContainer_H
#define AliEmcalContainer_H

//
// container with name, TClonesArray
//

class TLorentzVector;
class AliVEvent;
class AliNamedArrayI;
class AliVParticle;

#include <TNamed.h>
#include <TClonesArray.h>

class AliEmcalContainer : public TObject {
 public:
  enum RejectionReason {
    // General
    kNullObject = 1<<0,
    kPtCut = 1<<1,
    kAcceptanceCut = 1<<2,
    kMCLabelCut = 1<<3,
    kBitMapCut = 1<<4,
    // leave bits 4-6 free for future implementations
    
    // AliParticleContainer
    kNotHybridTrack = 1<<7,
    kMCFlag = 1<<8,
    kMCGeneratorCut = 1<<9,
    kChargeCut = 1<<10,
    kMinDistanceTPCSectorEdgeCut = 1<<11,

    // AliClusterContainer
    kIsEMCalCut = 1<<13,
    kTimeCut = 1<<14,
    kEnergyCut = 1<<15,
    kExoticCut = 1<<16,
    // leave bit 17 free for future implementations

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
    kTagStatus = 1<<27,
    kMinNConstituents = 1<<28
  };

  AliEmcalContainer();
  AliEmcalContainer(const char *name); 
  virtual ~AliEmcalContainer(){;}

  TClonesArray               *GetArray() const                      { return fClArray                   ; }
  const TString&              GetArrayName()                  const { return fClArrayName               ; }
  const TString&              GetClassName()                  const { return fClassName                 ; }
  Int_t                       GetCurrentID()                  const { return fCurrentID                 ; }
  Bool_t                      GetIsParticleLevel()            const { return fIsParticleLevel           ; }
  Int_t                       GetIndexFromLabel(Int_t lab)    const;
  Int_t                       GetNEntries()                   const { return fClArray->GetEntriesFast() ; }
  virtual Bool_t              GetMomentum(TLorentzVector &mom, Int_t i) = 0;
  virtual Bool_t              GetAcceptMomentum(TLorentzVector &mom, Int_t i) = 0;
  virtual Bool_t              GetNextMomentum(TLorentzVector &mom, Int_t i=-1) = 0;
  virtual Bool_t              GetNextAcceptMomentum(TLorentzVector &mom, Int_t i=-1) = 0;
  virtual Bool_t              AcceptObject(Int_t i) = 0;
  virtual Bool_t              AcceptObject(TObject* obj) = 0;
  void                        ResetCurrentID(Int_t i=-1)            { fCurrentID = i                    ; }
  virtual void                SetArray(AliVEvent *event);
  void                        SetArrayName(const char *n)           { fClArrayName = n                  ; }
  void                        SetBitMap(UInt_t m)                   { fBitMap = m                       ; }
  void                        SetIsParticleLevel(Bool_t b)          { fIsParticleLevel = b              ; }
  void                        SortArray()                           { fClArray->Sort()                  ; }
  UInt_t                      GetRejectionReason()            const { return fRejectionReason           ; }
  UInt_t                      TestRejectionReason(UInt_t rs)  const { return fRejectionReason & rs      ; }
  UShort_t                    GetRejectionReasonBitPosition() const;
  TClass*                     GetLoadedClass()                      { return fLoadedClass               ; }
  virtual void                NextEvent() {;}

  const char*                 GetName()                       const { return fName.Data()               ; }
  void                        SetName(const char* n)                { fName = n                         ; }

  static Bool_t               SamePart(const AliVParticle* part1, const AliVParticle* part2, Double_t dist = 1.e-4);

 protected:
  TString                     fName;                    // object name
  TString                     fClArrayName;             // name of branch
  TString                     fClassName;               // name of the class in the TClonesArray
  Bool_t                      fIsParticleLevel;         // whether or not it is a particle level object collection
  UInt_t                      fBitMap;                  // bitmap mask
  TClonesArray               *fClArray;                 //!TClonesArray
  Int_t                       fCurrentID;               //!current ID for automatic loops
  AliNamedArrayI             *fLabelMap;                //!Label-Index map
  Double_t                    fVertex[3];               //!event vertex array
  UInt_t                      fRejectionReason;         //!reject reason bit map for the last call to an accept object function
  TClass                     *fLoadedClass;             //!Class of teh objects contained in the TClonesArray

 private:
  AliEmcalContainer(const AliEmcalContainer& obj); // copy constructor
  AliEmcalContainer& operator=(const AliEmcalContainer& other); // assignment

  ClassDef(AliEmcalContainer,6);
};
#endif
