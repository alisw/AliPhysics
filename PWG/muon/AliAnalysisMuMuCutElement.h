#ifndef ALIANALYSISMUMUCUTELEMENT_H
#define ALIANALYSISMUMUCUTELEMENT_H

/**
 *
 * \class AliAnalysisMuMuCutElement
 *
 * \brief Describes an elementary cut (either event cut, track cut, pair cut, or trigger class cut)
 *
 */

#include "TObject.h"
#include "TString.h"

#include <vector>

class TMethodCall;
class AliVEvent;
class AliVEventHandler;
class AliVParticle;

class AliAnalysisMuMuCutElement : public TObject
{
public:

  enum ECutType
  {
    kEvent=0, // a cut on event
    kTrack=1, // a cut on track
    kTrackPair=2, // a cut on track pair
    kTriggerClass=3, // a cut on trigger class
    kAny=4 // must be the last one
  };

  static const char* CutTypeName(ECutType type);

  AliAnalysisMuMuCutElement();

  AliAnalysisMuMuCutElement(ECutType expectedType,
                            TObject& cutObject,
                            const char* cutMethodName,
                            const char* cutMethodPrototype,
                            const char* defaultParameters);

  virtual ~AliAnalysisMuMuCutElement();

  virtual Bool_t IsValid() const { return (fCutMethod != 0x0); }

  const char* GetName() const { return fName.Data(); }

  virtual Bool_t Pass(const AliVEvent& event) const;

  virtual Bool_t Pass(const AliVEventHandler& eventHandler) const;

  virtual Bool_t Pass(const AliVParticle& particle) const;

  virtual Bool_t Pass(const AliVParticle& p1, const AliVParticle& p2) const;

  virtual Bool_t Pass(const TString& firedTriggerClasses, TString& acceptedTriggerClasses,
                      UInt_t L0, UInt_t L1, UInt_t L2) const;

  virtual void Print(Option_t* opt="") const;

  Bool_t IsEventCutter() const { return fIsEventCutter; }
  Bool_t IsEventHandlerCutter() const { return fIsEventHandlerCutter; }
  Bool_t IsTrackCutter() const { return fIsTrackCutter; }
  Bool_t IsTrackPairCutter() const { return fIsTrackPairCutter; }
  Bool_t IsTriggerClassCutter() const { return fIsTriggerClassCutter; }

  TObject* GetCutObject() const { return fCutObject; }

  const Long_t* GetCallParams() const { return &fCallParams[0]; }

  const char* GetCallMethodName() const;
  const char* GetCallMethodProto() const;

  Bool_t IsEqual(const TObject* obj) const;

private:

  void Init(ECutType type=kAny) const;

  Bool_t CallCutMethod(Long_t p) const;
  Bool_t CallCutMethod(Long_t p1, Long_t p2) const;

  Int_t CountOccurences(const TString& prototype, const char* search) const;

  /// not implemented on purpose
  AliAnalysisMuMuCutElement(const AliAnalysisMuMuCutElement& rhs);
  /// not implemented on purpose
  AliAnalysisMuMuCutElement& operator=(const AliAnalysisMuMuCutElement& rhs);

protected:
  TString fName; // name of the cut
  mutable Bool_t fIsEventCutter; // whether or not the cut is an event cutter
  mutable Bool_t fIsEventHandlerCutter; // whether or not the cut is an event handler cutter
  mutable Bool_t fIsTrackCutter; // whether or not the cut is a track cutter
  mutable Bool_t fIsTrackPairCutter; // whether or not the cut is a track pair cutter
  mutable Bool_t fIsTriggerClassCutter; // whether or not the cut is a trigger class cutter
private:

  TObject* fCutObject; // pointer (not owner) to the object doing the actual cut work
  TString fCutMethodName; // method (of fCutObject) to be called to do the cut
  TString fCutMethodPrototype; // prototype of the method to be called to do the cut
  TString fDefaultParameters; // default parameters of the cut method (might be empty)
  mutable Int_t fNofParams; // number of parameters
  mutable TMethodCall* fCutMethod; //! cut method

  mutable std::vector<Long_t> fCallParams; //! vector of parameters for the fCutMethod
  mutable std::vector<Double_t> fDoubleParams; //! temporary vector to hold the references

  ClassDef(AliAnalysisMuMuCutElement,1) // One piece of a cut combination
};

class AliAnalysisMuMuCutElementBar : public AliAnalysisMuMuCutElement
{
public:
  AliAnalysisMuMuCutElementBar();

  AliAnalysisMuMuCutElementBar(const AliAnalysisMuMuCutElement& ce);

  virtual ~AliAnalysisMuMuCutElementBar();

  Bool_t IsValid() const { return fCutElement && fCutElement->IsValid(); }

  Bool_t Pass(const AliVEvent& event) const { return !fCutElement->Pass(event); }

  Bool_t Pass(const AliVEventHandler& eventHandler) const { return !fCutElement->Pass(eventHandler); }

  Bool_t Pass(const AliVParticle& particle) const { return !fCutElement->Pass(particle); }

  Bool_t Pass(const AliVParticle& p1, const AliVParticle& p2) const { return !fCutElement->Pass(p1,p2); }

  Bool_t Pass(const TString& firedTriggerClasses, TString& acceptedTriggerClasses,
              UInt_t L0, UInt_t L1, UInt_t L2) const
  { return fCutElement->Pass(firedTriggerClasses,acceptedTriggerClasses,L0,L1,L2); }

  void Print(Option_t* opt="") const;

private:

  /// not implemented on purpose
  AliAnalysisMuMuCutElementBar(const AliAnalysisMuMuCutElementBar& rhs);
  /// not implemented on purpose
  AliAnalysisMuMuCutElementBar& operator=(const AliAnalysisMuMuCutElementBar& rhs);

  const AliAnalysisMuMuCutElement* fCutElement; // the cut element we're the negation of

  ClassDef(AliAnalysisMuMuCutElementBar,1) // opposite of cut element (i.e. !cutelement)
};

#endif
