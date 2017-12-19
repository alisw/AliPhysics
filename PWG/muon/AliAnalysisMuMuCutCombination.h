#ifndef ALIANALYSISMUMUCUTCOMBINATION_H
#define ALIANALYSISMUMUCUTCOMBINATION_H

/**
 *
 * \class AliAnalysisMuMuCutCombination
 *
 * \brief Defines a cut by combining several cut elements
 *
 * \author L. Aphecetche (Subatech)
 */

#include "TObject.h"
#include "TString.h"

class AliAnalysisMuMuCutElement;
class AliVEvent;
class AliVEventHandler;
class AliVParticle;

class AliAnalysisMuMuCutCombination : public TObject
{
public:
  AliAnalysisMuMuCutCombination();
  virtual ~AliAnalysisMuMuCutCombination();

  void Add(AliAnalysisMuMuCutElement* ce);

  Bool_t Pass(const AliVEventHandler& eventHandler) const;

  Bool_t Pass(const AliVParticle& particle) const;

  Bool_t Pass(const AliVParticle& p1, const AliVParticle& p2) const;

  Bool_t Pass(const TString& firedTriggerClasses, TString& acceptedTriggerClasses,
              UInt_t L0, UInt_t L1, UInt_t L2) const;

  const char* GetName() const { return fName.Data(); }

  Bool_t IsEventCutter() const { return fIsEventCutter; }
  Bool_t IsEventHandlerCutter() const { return fIsEventHandlerCutter; }
  Bool_t IsTrackCutter() const { return fIsTrackCutter; }
  Bool_t IsTrackPairCutter() const { return fIsTrackPairCutter; }
  Bool_t IsTriggerClassCutter() const { return fIsTriggerClassCutter; }

  void Print(Option_t* opt="") const;

  Bool_t IsEqual(const TObject* obj) const;

  Bool_t IsEqualForTrackCutter(const AliAnalysisMuMuCutCombination& other) const;

private:
  /// not implemented on purpose
  AliAnalysisMuMuCutCombination(const AliAnalysisMuMuCutCombination& rhs);
  /// not implemented on purpose
  AliAnalysisMuMuCutCombination& operator=(const AliAnalysisMuMuCutCombination& rhs);

private:
  TObjArray* fCuts; // array of cut elements that form this cut combination
  TString fName; // name of the combination
  Bool_t fIsEventCutter; // whether or not the combination cuts on event
  Bool_t fIsEventHandlerCutter; // whether or not the combination cuts on event handler
  Bool_t fIsTrackCutter; // whether or not the combination cuts on track
  Bool_t fIsTrackPairCutter; // whether or not the combination cuts on track pairs
  Bool_t fIsTriggerClassCutter; // whether or not the combination cuts on trigger class

  ClassDef(AliAnalysisMuMuCutCombination,1) // combination of 1 or more individual cuts
};

#endif
