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
  AliEmcalContainer();
  AliEmcalContainer(const char *name); 
  virtual ~AliEmcalContainer(){;}

  TClonesArray               *GetArray()                          { return fClArray                   ; }
  const TString&              GetArrayName()                const { return fClArrayName               ; }
  Int_t                       GetCurrentID()                const { return fCurrentID-1               ; }
  Bool_t                      GetIsParticleLevel()          const { return fIsParticleLevel           ; }
  Int_t                       GetIndexFromLabel(Int_t lab)  const;
  Int_t                       GetNEntries()                 const { return fClArray->GetEntriesFast() ; }
  virtual void                GetMomentum(TLorentzVector &mom, Int_t i) const = 0;
  void                        ResetCurrentID(Int_t i=0)           { fCurrentID = i                    ; }
  virtual void                SetArray(AliVEvent *event);
  void                        SetArrayName(const char *n)         { fClArrayName = n                  ; }
  void                        SetIsParticleLevel(Bool_t b)        { fIsParticleLevel = b              ; }
  void                        SortArray()                         { fClArray->Sort()                  ; }

 protected:
  TString                     fClArrayName;             // name of branch
  TString                     fClassName;               // name of the class in the TClonesArray
  Bool_t                      fIsParticleLevel;         // whether or not it is a particle level object collection
  TClonesArray               *fClArray;                 //!TClonesArray
  Int_t                       fCurrentID;               //!current ID for automatic loops
  AliNamedArrayI             *fLabelMap;                //!Label-Index map
  Double_t                    fVertex[3];               //!event vertex array

 private:
  AliEmcalContainer(const AliEmcalContainer& obj); // copy constructor
  AliEmcalContainer& operator=(const AliEmcalContainer& other); // assignment

  ClassDef(AliEmcalContainer,3);
};
#endif
