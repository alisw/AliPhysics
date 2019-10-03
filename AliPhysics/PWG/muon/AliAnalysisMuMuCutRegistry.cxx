#include "AliAnalysisMuMuCutRegistry.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuCutRegistry
 * The cut elements and cut combinations are stored per type, i.e. there's a set for
 * event cuts/combinations, a set for track cuts/combinations, and a set for track pair cuts/combinations.
 *
 * To define a new cut use the AddEventCut, AddTrackCut, AndTrackPairCut and AddTriggerClassCut methods.
 *
 * To add an existing cut see AddCutElement.
 *
 * To define the negation of a cut, use the Not method.
 *
 * To add a new combination use one of the AddCutCombination methods, depending on the number
 * of cut element(s) the combination is made of.
 *
 * Note that what the sub-analysis are really concerned with are cut combinations (i.e. if
 * you fail to define any combination, nothing will be cut, whatever the number of cut elements
 * you've defined).
 *
 * This class also defines a few default control cut elements aptly named AlwaysTrue.
 *
 */

#include <utility>
#include "AliLog.h"
#include "TMethodCall.h"
#include "AliVEvent.h"
#include <set>
#include "AliAnalysisMuMuCutElement.h"
#include "AliAnalysisMuMuCutCombination.h"
#include "TObjArray.h"
#include "Riostream.h"
#include "TList.h"

ClassImp(AliAnalysisMuMuCutRegistry)

//_____________________________________________________________________________
AliAnalysisMuMuCutRegistry::AliAnalysisMuMuCutRegistry()
: TObject(),
fCutElements(0x0),
fCutCombinations(0x0)
{
  /// ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuCutRegistry::~AliAnalysisMuMuCutRegistry()
{
  /// dtor

  delete fCutElements;
  delete fCutCombinations;
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuCutRegistry::AddCutCombination(const TObjArray& cutElements)
{
  /// Add a cut combination composed of the cuts in the cutElements array.
  ///
  /// \return 1 in case of success, 0 if already there, -1 failure

  if ( cutElements.IsEmpty() ) return -1;

  AliAnalysisMuMuCutCombination* cutCombination = new AliAnalysisMuMuCutCombination;

  TIter next(&cutElements);
  AliAnalysisMuMuCutElement* ce;

  while ( ( ce = static_cast<AliAnalysisMuMuCutElement*>(next()) ) )
  {
    cutCombination->Add(ce);
  }

  if ( GetCutCombinations(AliAnalysisMuMuCutElement::kAny)->FindObject(cutCombination) )
  {
    delete cutCombination;
    return 0;
  }

  GetCutCombinations(AliAnalysisMuMuCutElement::kAny)->Add(cutCombination);

  if ( cutCombination->IsEventCutter() || cutCombination->IsEventHandlerCutter() )
  {
    GetCutCombinations(AliAnalysisMuMuCutElement::kEvent)->Add(cutCombination);
  }

  if ( cutCombination->IsTrackCutter() )
  {
    TObjArray* a = GetCutCombinations(AliAnalysisMuMuCutElement::kTrack);
    TIter nextCutComb(a);
    AliAnalysisMuMuCutCombination* other;
    Bool_t alreadyThere(kFALSE);

    while ( ( other = static_cast<AliAnalysisMuMuCutCombination*>(nextCutComb())) && !alreadyThere )
    {
        if ( cutCombination->IsEqualForTrackCutter(*other) )
        {
          alreadyThere = kTRUE;
        }
    }

    if (!alreadyThere)
    {
      a->Add(cutCombination);
    }
  }

  if ( cutCombination->IsTrackPairCutter() )
  {
    GetCutCombinations(AliAnalysisMuMuCutElement::kTrackPair)->Add(cutCombination);
  }

  if ( cutCombination->IsTriggerClassCutter() )
  {
    GetCutCombinations(AliAnalysisMuMuCutElement::kTriggerClass)->Add(cutCombination);
  }

  return 1;
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuCutRegistry::AddCutCombination(AliAnalysisMuMuCutElement* ce1)
{
  /// Convenience method to create a cut combination made of a single cut
  TObjArray cutElements;
  if ( ce1 ) cutElements.Add(ce1);
  return AddCutCombination(cutElements);
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuCutRegistry::AddCutCombination(AliAnalysisMuMuCutElement* ce1,
                                                    AliAnalysisMuMuCutElement* ce2)
{
  /// Convenience method to create a cut combination made of 2 cuts
  TObjArray cutElements;
  if ( ce1 ) cutElements.Add(ce1);
  if ( ce2 ) cutElements.Add(ce2);
  return AddCutCombination(cutElements);
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuCutRegistry::AddCutCombination(AliAnalysisMuMuCutElement* ce1,
                                                    AliAnalysisMuMuCutElement* ce2,
                                                    AliAnalysisMuMuCutElement* ce3)
{
  /// Convenience method to create a cut combination made of 3 cuts
  TObjArray cutElements;
  if ( ce1 ) cutElements.Add(ce1);
  if ( ce2 ) cutElements.Add(ce2);
  if ( ce3 ) cutElements.Add(ce3);
  return AddCutCombination(cutElements);
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuCutRegistry::AddCutCombination(AliAnalysisMuMuCutElement* ce1, AliAnalysisMuMuCutElement* ce2, AliAnalysisMuMuCutElement* ce3,
                        AliAnalysisMuMuCutElement* ce4)
{
  /// Convenience method to create a cut combination made of 4 cuts
  TObjArray cutElements;
  if ( ce1 ) cutElements.Add(ce1);
  if ( ce2 ) cutElements.Add(ce2);
  if ( ce3 ) cutElements.Add(ce3);
  if ( ce4 ) cutElements.Add(ce4);
  return AddCutCombination(cutElements);
}
//_____________________________________________________________________________
Int_t AliAnalysisMuMuCutRegistry::AddCutCombination(AliAnalysisMuMuCutElement* ce1, AliAnalysisMuMuCutElement* ce2, AliAnalysisMuMuCutElement* ce3,
                        AliAnalysisMuMuCutElement* ce4, AliAnalysisMuMuCutElement* ce5)
{
  /// Convenience method to create a cut combination made of 5 cuts

  TObjArray cutElements;
  if ( ce1 ) cutElements.Add(ce1);
  if ( ce2 ) cutElements.Add(ce2);
  if ( ce3 ) cutElements.Add(ce3);
  if ( ce4 ) cutElements.Add(ce4);
  if ( ce5 ) cutElements.Add(ce5);
  return AddCutCombination(cutElements);
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuCutRegistry::AddCutCombination(AliAnalysisMuMuCutElement* ce1, AliAnalysisMuMuCutElement* ce2, AliAnalysisMuMuCutElement* ce3,
                        AliAnalysisMuMuCutElement* ce4, AliAnalysisMuMuCutElement* ce5, AliAnalysisMuMuCutElement* ce6)
{
  /// Convenience method to create a cut combination made of 6 cuts

  TObjArray cutElements;
  if ( ce1 ) cutElements.Add(ce1);
  if ( ce2 ) cutElements.Add(ce2);
  if ( ce3 ) cutElements.Add(ce3);
  if ( ce4 ) cutElements.Add(ce4);
  if ( ce5 ) cutElements.Add(ce5);
  if ( ce6 ) cutElements.Add(ce6);
  return AddCutCombination(cutElements);
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElement*
AliAnalysisMuMuCutRegistry::CreateCutElement(AliAnalysisMuMuCutElement::ECutType type,
                                             TObject& cutClass,
                                             const char* cutMethodName,
                                             const char* cutMethodPrototype,
                                             const char* defaultParameters)
{
  /** Create a cut element. See the ctor of AliAnalysisMuMuCutElement for the meaning
   * of the parameters.
   */

  AliAnalysisMuMuCutElement* ce = new AliAnalysisMuMuCutElement(type,cutClass,cutMethodName,cutMethodPrototype,defaultParameters);

  AliAnalysisMuMuCutElement* added = AddCutElement(ce);

  if (!added)
  {
    delete ce;
  }

  return added;
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElement*
AliAnalysisMuMuCutRegistry::AddCutElement(AliAnalysisMuMuCutElement* ce)
{
  /// Add an existing cut element to the registry if it is valid

  if ( ce && ce->IsValid() )
  {
    if (!GetCutElements(AliAnalysisMuMuCutElement::kAny)->FindObject(ce))
    {
      GetCutElements(AliAnalysisMuMuCutElement::kAny)->Add(ce);
      if ( ce->IsEventCutter() || ce->IsEventHandlerCutter() )
      {
        GetCutElements(AliAnalysisMuMuCutElement::kEvent)->Add(ce);
      }
      else if ( ce->IsTrackCutter() )
      {
        GetCutElements(AliAnalysisMuMuCutElement::kTrack)->Add(ce);
      }
      else if ( ce->IsTrackPairCutter() )
      {
        GetCutElements(AliAnalysisMuMuCutElement::kTrackPair)->Add(ce);
      }
      else if ( ce->IsTriggerClassCutter() )
      {
        GetCutElements(AliAnalysisMuMuCutElement::kTriggerClass)->Add(ce);
      }
    }
    return ce;
  }
  return 0x0;
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElement* AliAnalysisMuMuCutRegistry::AddEventCut(TObject& cutClass,
                                       const char* cutMethodName,
                                       const char* cutMethodPrototype,
                                       const char* defaultParameters)
{
  /// Shortcut method to create a cut element of type kEvent
  return CreateCutElement(AliAnalysisMuMuCutElement::kEvent,cutClass,cutMethodName,
                          cutMethodPrototype,defaultParameters);
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElement* AliAnalysisMuMuCutRegistry::AddTrackCut(TObject& cutClass,
                                       const char* cutMethodName,
                                       const char* cutMethodPrototype,
                                       const char* defaultParameters)
{
  /// Shortcut method to create a cut element of type kTrack
  return CreateCutElement(AliAnalysisMuMuCutElement::kTrack,cutClass,cutMethodName,
                          cutMethodPrototype,defaultParameters);
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElement* AliAnalysisMuMuCutRegistry::AddTrackPairCut(TObject& cutClass,
                                        const char* cutMethodName,
                                        const char* cutMethodPrototype,
                                        const char* defaultParameters)
{
  /// Shortcut method to create a cut element of type kTrackPair
  return CreateCutElement(AliAnalysisMuMuCutElement::kTrackPair,cutClass,cutMethodName,
                          cutMethodPrototype,defaultParameters);
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElement* AliAnalysisMuMuCutRegistry::AddTriggerClassCut(TObject& cutClass,
                                              const char* cutMethodName,
                                              const char* cutMethodPrototype,
                                              const char* defaultParameters)
{
  /// Shortcut method to create a cut element of type kTriggerClass
  return CreateCutElement(AliAnalysisMuMuCutElement::kTriggerClass,cutClass,cutMethodName,
                          cutMethodPrototype,defaultParameters);
}

//_____________________________________________________________________________
const TObjArray* AliAnalysisMuMuCutRegistry::GetCutCombinations(AliAnalysisMuMuCutElement::ECutType type) const
{
  /// Get (and create if not already done) the array of cut combinations of the given type

  if (!fCutCombinations) return 0x0;

  return static_cast<TObjArray*>(fCutCombinations->At(type));
}

//_____________________________________________________________________________
TObjArray* AliAnalysisMuMuCutRegistry::GetCutCombinations(AliAnalysisMuMuCutElement::ECutType type)
{
  /// Get (and create if not already done) the array of cut combinations of the given type

  if (!fCutCombinations)
  {
    // the fCutCombinations array will be the owner of all the cut combinations

    Int_t N = AliAnalysisMuMuCutElement::kAny + 1;

    fCutCombinations = new TObjArray(N);
    fCutCombinations->SetOwner(kTRUE);

    for ( Int_t i = 0; i < N; ++i )
    {
      TObjArray* array = new TObjArray;
      array->SetOwner(kFALSE);
      if (i==AliAnalysisMuMuCutElement::kAny)
      {
        // only the first array, containing all the combinations
        // is the owner of the combinations
        // the other arrays are just pointing to those
        array->SetOwner(kTRUE);
      }
      fCutCombinations->AddAt(array,i);
    }
  }
  return static_cast<TObjArray*>(fCutCombinations->At(type));
}

//_____________________________________________________________________________
const TObjArray* AliAnalysisMuMuCutRegistry::GetCutElements(AliAnalysisMuMuCutElement::ECutType type) const
{
  /// Get the array of cut elements of the given type. Return 0x0 if the array does not exist yet

  if (!fCutElements) return 0x0;

  return static_cast<TObjArray*>(fCutElements->At(type));
}

//_____________________________________________________________________________
TObjArray* AliAnalysisMuMuCutRegistry::GetCutElements(AliAnalysisMuMuCutElement::ECutType type)
{
  /// Get (and create if not already done) the array of cut elements of the given type

  if (!fCutElements)
  {
    // owner of all the cut elements
    Int_t N = AliAnalysisMuMuCutElement::kAny + 1;
    fCutElements = new TObjArray(N);
    fCutElements->SetOwner(kTRUE);

    for ( Int_t i = 0; i < N; ++i )
    {
      TObjArray* array = new TObjArray;
      array->SetOwner(kFALSE);
      if (i == AliAnalysisMuMuCutElement::kAny )
      {
        // only the first array is the owner of the cuts
        // the other ones are just pointing to this one
        array->SetOwner(kTRUE);
      }
      fCutElements->AddAt(array,i);
    }
  }
  return static_cast<TObjArray*>(fCutElements->At(type));
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElement* AliAnalysisMuMuCutRegistry::Not(const AliAnalysisMuMuCutElement& cutElement)
{
  /// Create a cut which is the opposite of cutElement, and adds it.

  AliAnalysisMuMuCutElementBar* bar = new AliAnalysisMuMuCutElementBar(cutElement);

  AliAnalysisMuMuCutElement* added = AddCutElement(bar);

  if (!added)
  {
    delete bar;
  }

  return added;
}

//_____________________________________________________________________________
void AliAnalysisMuMuCutRegistry::Print(Option_t* opt) const
{
  /// Printout

  TString sopt(opt);
  sopt.ToUpper();

  std::cout << "++++ Cut combinations defined : " << std::endl;

  AliAnalysisMuMuCutElement::ECutType cutTypes[] = { AliAnalysisMuMuCutElement::kEvent, AliAnalysisMuMuCutElement::kTrack,
    AliAnalysisMuMuCutElement::kTrackPair, AliAnalysisMuMuCutElement::kTriggerClass };

  Int_t i(1);

  for ( Int_t iCutType = 0; iCutType < 4; ++iCutType )
  {
    if (GetCutElements(cutTypes[iCutType])->IsEmpty()) continue;
    std::cout << "  Cutting on " << AliAnalysisMuMuCutElement::CutTypeName(cutTypes[iCutType]) << std::endl;
    TIter next(GetCutCombinations(cutTypes[iCutType]));
    AliAnalysisMuMuCutCombination* cutCombination;

    while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(next())) )
    {
      std::cout << Form("    %4d ",i);
      cutCombination->Print(sopt.Data());
      ++i;
    }
  }

  if ( sopt.Contains("FULL") || sopt.Contains("ALL") )
  {
    std::cout << "++++ Individual cuts defined : " << std::endl;

    for ( Int_t iCutType = 0; iCutType < 4; ++iCutType )
    {
      if (GetCutElements(cutTypes[iCutType])->IsEmpty()) continue;
      std::cout << "  Cutting on " << AliAnalysisMuMuCutElement::CutTypeName(cutTypes[iCutType]) << std::endl;
      TIter nextCutRef(GetCutElements(cutTypes[iCutType]));
      AliAnalysisMuMuCutElement* ce;
      i = 1;
      while ( ( ce = static_cast<AliAnalysisMuMuCutElement*>(nextCutRef()) ) )
      {
        std::cout << Form("%4d ",i);
        ce->Print(sopt.Data());
        ++i;
      }
    }
  }
}
