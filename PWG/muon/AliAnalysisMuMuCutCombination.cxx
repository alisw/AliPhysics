#include "AliAnalysisMuMuCutCombination.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuCutCombination
 *
 * A cut combination is the real cut class that is used in the sub-analysis. It is composed
 * of one or several AliAnalysisMuMuCutElement.
 *
 * Unlike the cut elements, it can be of several types at the same time (i.e. it can
 * be an event cutter and a track cutter for instance). The work is done in the
 * Pass method(s).
 *
 */

#include "AliAnalysisMuMuCutElement.h"
#include "TList.h"
#include "Riostream.h"
#include "AliVEventHandler.h"
#include "AliLog.h"

ClassImp(AliAnalysisMuMuCutCombination)

//_____________________________________________________________________________
AliAnalysisMuMuCutCombination::AliAnalysisMuMuCutCombination()
: TObject(), fCuts(0x0), fName(""),
fIsEventCutter(kFALSE), fIsEventHandlerCutter(kFALSE),
fIsTrackCutter(kFALSE), fIsTrackPairCutter(kFALSE),
fIsTriggerClassCutter(kFALSE)
{
  /// Default ctor.
}

//_____________________________________________________________________________
AliAnalysisMuMuCutCombination::~AliAnalysisMuMuCutCombination()
{
  /// Dtor
  delete fCuts;
}

//_____________________________________________________________________________
void AliAnalysisMuMuCutCombination::Add(AliAnalysisMuMuCutElement* ce)
{
  /** Add a cut element to this combination, if the cut is not void and
   *  not already part of the combination
   */

  if (!ce || !ce->IsValid()) return;

  if (!fCuts)
  {
    fCuts = new TObjArray;
    fCuts->SetOwner(kFALSE);
    fIsEventCutter = ce->IsEventCutter();
    fIsEventHandlerCutter = ce->IsEventHandlerCutter();
    fIsTrackCutter = ce->IsTrackCutter();
    fIsTrackPairCutter = ce->IsTrackPairCutter();
    fIsTriggerClassCutter = ce->IsTriggerClassCutter();
  }

  if (!fCuts->FindObject(ce))
  {
    fCuts->Add(ce);
    fName += ce->GetName();

    fIsEventCutter = fIsEventCutter || ce->IsEventCutter();
    fIsEventHandlerCutter = fIsEventHandlerCutter || ce->IsEventHandlerCutter();
    fIsTrackCutter = fIsTrackCutter || ce->IsTrackCutter();
    fIsTrackPairCutter = fIsTrackPairCutter || ce->IsTrackPairCutter();
    fIsTriggerClassCutter = fIsTriggerClassCutter || ce->IsTriggerClassCutter();

  }

  // update the name

  if ( IsTrackPairCutter() )
  {
    if ( fName[0] == 's' )
    {
      fName[0] = 'p';
    }
    else if ( fName[0] == 'p')
    {
      // already ok
    }
    else
    {
      TString tmp = fName;
      fName = "p";
      fName += tmp;
    }
  }
  else if ( IsTrackCutter() )
  {
    if ( fName[0] == 's')
    {
      // already ok
    }
    else
    {
      TString tmp = fName;
      fName = "s";
      fName += tmp;
    }
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutCombination::IsEqualForTrackCutter(const AliAnalysisMuMuCutCombination& other) const
{
  /// whether or not we are the same combination as other, only considering the cuts
  /// or the given type

  for ( Int_t i = 0; i <= fCuts->GetLast(); ++i )
  {
    AliAnalysisMuMuCutElement* thisCuti = static_cast<AliAnalysisMuMuCutElement*>(fCuts->At(i));

    if ( thisCuti->IsTrackCutter() && !other.fCuts->FindObject(thisCuti)  )
    {
      return kFALSE;
    }
  }

  for ( Int_t i = 0; i <= other.fCuts->GetLast(); ++i )
  {
    AliAnalysisMuMuCutElement* otherCuti = static_cast<AliAnalysisMuMuCutElement*>(other.fCuts->At(i));

    if ( otherCuti->IsTrackCutter() && !fCuts->FindObject(otherCuti)  )
    {
      return kFALSE;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutCombination::IsEqual(const TObject* obj) const
{
  /// Whether or not we are the same cut combination as obj

  if ( obj->IsA() != AliAnalysisMuMuCutCombination::Class() )
  {
    return kFALSE;
  }

  const AliAnalysisMuMuCutCombination* other = static_cast<const AliAnalysisMuMuCutCombination*>(obj);

  if ( IsEventCutter() != other->IsEventCutter() ) return kFALSE;
  if ( IsTrackCutter() != other->IsTrackCutter() ) return kFALSE;
  if ( IsTrackPairCutter() != other->IsTrackPairCutter() ) return kFALSE;
  if ( IsTriggerClassCutter() != other->IsTriggerClassCutter() ) return kFALSE;

  if ( !fCuts && !other->fCuts ) return kTRUE;

  // no cuts, nothing to check further...

  if ( ( fCuts && !other->fCuts ) || ( !fCuts && other->fCuts ) ) return kFALSE;

  if ( fCuts->GetEntries() != other->fCuts->GetEntries() ) return kFALSE;

  // ok, looks similar so far, now have compute the set of cuts in common
  // (whatever the order) to see if they really are the same combination or not

  Int_t n1in2(0);
  Int_t n2in1(0);

  for ( Int_t i = 0; i <= fCuts->GetLast(); ++i )
  {
    AliAnalysisMuMuCutElement* thisCuti = static_cast<AliAnalysisMuMuCutElement*>(fCuts->At(i));

    if ( other->fCuts->FindObject(thisCuti) )
    {
      ++n1in2;
    }
  }

  for ( Int_t i = 0; i <= fCuts->GetLast(); ++i )
  {
    AliAnalysisMuMuCutElement* otherCuti = static_cast<AliAnalysisMuMuCutElement*>(other->fCuts->At(i));

    if ( fCuts->FindObject(otherCuti) )
    {
      ++n2in1;
    }
  }

  return (n1in2==n2in1 && n1in2==fCuts->GetLast()+1);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutCombination::Pass(const AliVEventHandler& eventHandler) const
{
  /// Whether or not the event handler is passing the cut

  if (!fCuts) return kFALSE;
  TIter next(fCuts);
  AliAnalysisMuMuCutElement* ce;

  const AliVEvent* event = eventHandler.GetEvent();

  Bool_t passEvent(kTRUE);
  Bool_t passEventHandler(kTRUE);

  while ( ( ce = static_cast<AliAnalysisMuMuCutElement*>(next()) ) )
  {
    if ( ce->IsEventCutter() && !ce->Pass(*event) )
    {
      passEvent = kFALSE;
    }

    if ( ce->IsEventHandlerCutter() && !ce->Pass(eventHandler) )
    {
      passEventHandler = kFALSE;
    }
  }

  if ( IsEventCutter() && IsEventHandlerCutter() )
  {
    return passEvent && passEventHandler;
  }

  if ( IsEventHandlerCutter() )
  {
    return passEventHandler;
  }

  if ( IsEventCutter() )
  {
    return passEvent;
  }

  return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutCombination::Pass(const AliVParticle& particle) const
{
  /// Whether or not the particle is passing the cut

  if (!fCuts) return kFALSE;
  TIter next(fCuts);
  AliAnalysisMuMuCutElement* ce;

  while ( ( ce = static_cast<AliAnalysisMuMuCutElement*>(next()) ) )
  {
    if (ce->IsTrackCutter() && !ce->Pass(particle))
    {
      return kFALSE;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutCombination::Pass(const AliVParticle& p1, const AliVParticle& p2) const
{
  /// Whether or not the particle pair is passing the cut

  if (!fCuts) return kFALSE;
  TIter next(fCuts);
  AliAnalysisMuMuCutElement* ce;

  while ( ( ce = static_cast<AliAnalysisMuMuCutElement*>(next()) ) )
  {
    if (ce->IsTrackPairCutter() && !ce->Pass(p1,p2))
    {
      return kFALSE;
    }
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutCombination::Pass(const TString& firedTriggerClasses,
                                           TString& acceptedTriggerClasses,
                                           UInt_t L0, UInt_t L1, UInt_t L2) const
{
  /** Whether or not the firedTriggerClasses pass the cut.
   * \param firedTriggerClasses (input) list of fired trigger classes (separated by space)
   * \param acceptedTriggerClasses (output) list of accepted classes
   * \param L0 (input, optional) level 0 trigger mask
   * \param L1 (input, optional) level 1 trigger mask
   * \param L2 (input, optional) level 2 trigger mask
   */

  if (!fCuts) return kFALSE;
  TIter next(fCuts);
  AliAnalysisMuMuCutElement* ce;
  Bool_t rv(kFALSE);

  // contrary to the other cut types, here we make a full loop on all
  // the cuts, as we need to give each cut a chance to update the acceptedTriggerClasses
  // string

  acceptedTriggerClasses = "";

  while ( ( ce = static_cast<AliAnalysisMuMuCutElement*>(next()) ) )
  {
    TString tmp;

    if (ce->IsTriggerClassCutter() && ce->Pass(firedTriggerClasses,tmp,L0,L1,L2))
    {
      acceptedTriggerClasses += tmp;
      acceptedTriggerClasses += " ";
      rv = kTRUE;
    }
  }

  return rv;
}

//_____________________________________________________________________________
void AliAnalysisMuMuCutCombination::Print(Option_t* opt) const
{
  /// Printout of the cut combination
  TString sopt(opt);
  sopt.ToUpper();

  if ( sopt.Contains("PTR"))
  {
    std::cout << Form("%s(%p) [",GetName(),this);
  }
  else
  {
    std::cout << Form("%s [",GetName());
  }
  if ( IsEventCutter() ) std::cout << " E";
  if ( IsEventHandlerCutter() ) std::cout << " EH";
  if ( IsTrackCutter() ) std::cout << " T";
  if ( IsTrackPairCutter() ) std::cout << " TP";
  if ( IsTriggerClassCutter() ) std::cout << " TC";
  std::cout << " ]" << std::endl;

  TIter next(fCuts);
  AliAnalysisMuMuCutElement* ce;

  while ( ( ce = static_cast<AliAnalysisMuMuCutElement*>(next()) ) )
  {
    std::cout << "            ";
    ce->Print(sopt.Data());
  }
}
