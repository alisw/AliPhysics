#ifndef ALIANALYSISMUMUCUTREGISTRY_H
#define ALIANALYSISMUMUCUTREGISTRY_H

/**
 *
 * \class AliAnalysisMuMuCutRegistry
 *
 * \brief Container of AliAnalysisMuMuCutElement and AliAnalysisMuMuCutCombination objects.
 *
 * \author L. Aphecetche (Subatech)
 *
 */

#include "TObject.h"
#include "TString.h"
#include "TMethodCall.h"
#include "AliAnalysisMuMuCutElement.h"

class AliVEvent;
class AliAnalysisMuMuCutElementBar;
class AliAnalysisMuMuCutCombination;
class AliVParticle;
class AliVEventHandler;

class AliAnalysisMuMuCutRegistry : public TObject
{
public:
  AliAnalysisMuMuCutRegistry();
  virtual ~AliAnalysisMuMuCutRegistry();

  AliAnalysisMuMuCutElement* AddEventCut(TObject& cutClass,
                                           const char* cutMethodName,
                                           const char* cutMethodPrototype,
                                           const char* defaultParameters);

  AliAnalysisMuMuCutElement* AddTrackCut(TObject& cutClass,
                                           const char* cutMethodName,
                                           const char* cutMethodPrototype,
                                           const char* defaultParameters);

  AliAnalysisMuMuCutElement* AddTrackPairCut(TObject& cutClass,
                                             const char* cutMethodName,
                                             const char* cutMethodPrototype,
                                             const char* defaultParameters);

  AliAnalysisMuMuCutElement* AddTriggerClassCut(TObject& cutClass,
                                                const char* cutMethodName,
                                                const char* cutMethodPrototype,
                                                const char* defaultParameters);

  AliAnalysisMuMuCutElement* Not(const AliAnalysisMuMuCutElement& cutElement);

  AliAnalysisMuMuCutElement* AddCutElement(AliAnalysisMuMuCutElement* ce);

  Int_t AddCutCombination(const TObjArray& cutElements);

  Int_t AddCutCombination(AliAnalysisMuMuCutElement* ce1);
  Int_t AddCutCombination(AliAnalysisMuMuCutElement* ce1, AliAnalysisMuMuCutElement* ce2);
  Int_t AddCutCombination(AliAnalysisMuMuCutElement* ce1, AliAnalysisMuMuCutElement* ce2, AliAnalysisMuMuCutElement* ce3);
  Int_t AddCutCombination(AliAnalysisMuMuCutElement* ce1, AliAnalysisMuMuCutElement* ce2, AliAnalysisMuMuCutElement* ce3,
                          AliAnalysisMuMuCutElement* ce4);
  Int_t AddCutCombination(AliAnalysisMuMuCutElement* ce1, AliAnalysisMuMuCutElement* ce2, AliAnalysisMuMuCutElement* ce3,
                          AliAnalysisMuMuCutElement* ce4, AliAnalysisMuMuCutElement* ce5);
  Int_t AddCutCombination(AliAnalysisMuMuCutElement* ce1, AliAnalysisMuMuCutElement* ce2, AliAnalysisMuMuCutElement* ce3,
                          AliAnalysisMuMuCutElement* ce4, AliAnalysisMuMuCutElement* ce5, AliAnalysisMuMuCutElement* ce6);

  /// Get cut combinations of a given type
  const TObjArray* GetCutCombinations(AliAnalysisMuMuCutElement::ECutType type) const;
  TObjArray* GetCutCombinations(AliAnalysisMuMuCutElement::ECutType type);

  /// Get cut elements of a given type
  const TObjArray* GetCutElements(AliAnalysisMuMuCutElement::ECutType type) const;
  TObjArray* GetCutElements(AliAnalysisMuMuCutElement::ECutType type);

  virtual void Print(Option_t* opt="") const;

  Bool_t AlwaysTrue(const AliVEvent& /*event*/) const { return kTRUE; }
  void NameOfAlwaysTrue(TString& name) const { name="ALL"; }
  Bool_t AlwaysTrue(const AliVEventHandler& /*eventHandler*/) const { return kTRUE; }
  Bool_t AlwaysTrue(const AliVParticle& /*part*/) const { return kTRUE; }
  Bool_t AlwaysTrue(const AliVParticle& /*p1*/, const AliVParticle& /*p2*/) const { return kTRUE; }

private:

  /// not implemented on purpose
  AliAnalysisMuMuCutRegistry(const AliAnalysisMuMuCutRegistry& rhs);
  /// not implemented on purpose
  AliAnalysisMuMuCutRegistry& operator=(const AliAnalysisMuMuCutRegistry& rhs);

  AliAnalysisMuMuCutElement* CreateCutElement(AliAnalysisMuMuCutElement::ECutType expectedType,
                                              TObject& cutClass,
                                              const char* cutMethodName,
                                              const char* cutMethodPrototype,
                                              const char* defaultParameters);

private:

  mutable TObjArray* fCutElements; // cut elements
  mutable TObjArray* fCutCombinations; // cut combinations

  ClassDef(AliAnalysisMuMuCutRegistry,1) // storage for cut pointers
};

#endif
