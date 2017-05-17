#include "AliAnalysisMuMuCutElement.h"

/**
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuCutElement
 *
 * An AliAnalysisMuMuCutElement is an interface/proxy to another method (see ctor)
 *
 * A cut has a type (event cut, track cut, etc..., see ECutType and the IsXXX methods),
 * a name, and offers a Pass method. A cut can be of one type only.
 * Generally a real cut is made of several cut elements,
 * see \ref AliAnalysisMuMuCutCombination
 *
 *  \author L. Aphecetche (Subatech)
 */

#include "TMethodCall.h"
#include "AliLog.h"
#include "Riostream.h"
#include "AliVParticle.h"

ClassImp(AliAnalysisMuMuCutElement)
ClassImp(AliAnalysisMuMuCutElementBar)

//_____________________________________________________________________________
AliAnalysisMuMuCutElement::AliAnalysisMuMuCutElement()
: TObject(), fName(""), fIsEventCutter(kFALSE), fIsEventHandlerCutter(kFALSE),
fIsTrackCutter(kFALSE), fIsTrackPairCutter(kFALSE), fIsTriggerClassCutter(kFALSE),
fCutObject(0x0), fCutMethodName(""), fCutMethodPrototype(""),
fDefaultParameters(""), fNofParams(0), fCutMethod(0x0), fCallParams(), fDoubleParams()
{
  /// Default ctor, leading to an invalid cut object
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElement::AliAnalysisMuMuCutElement(ECutType expectedType,
                                                     TObject& cutObject,
                                                     const char* cutMethodName,
                                                     const char* cutMethodPrototype,
                                                     const char* defaultParameters)
: TObject(), fName(""), fIsEventCutter(kFALSE), fIsEventHandlerCutter(kFALSE),
fIsTrackCutter(kFALSE), fIsTrackPairCutter(kFALSE), fIsTriggerClassCutter(kFALSE),
fCutObject(&cutObject), fCutMethodName(cutMethodName),
fCutMethodPrototype(cutMethodPrototype),fDefaultParameters(defaultParameters),
fNofParams(0), fCutMethod(0x0), fCallParams(), fDoubleParams()
{
  /**
   * Construct a cut, which is a proxy to another method of (most probably) another object
   *
   * \param expectedType the type of cut which is expected (this cut will be invalidated if
   * the actual method does not comply with this type)
   * \param cutObject the object which has the cut method
   * \param cutMethodName the name of the method of cutObject that should do the cut work
   * \param cutMethodPrototype the prototype (i.e. list of arguments) of cutMethod
   * \param defaultParameters values of the default parameters (if any) of the cutMethod
   *
   * See the \ref Init method for details.
   */

  Init(expectedType);
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElement::~AliAnalysisMuMuCutElement()
{
  /// Dtor
  delete fCutMethod;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutElement::CallCutMethod(Long_t p) const
{
  /// Call the cut method with one parameter
  if (!fCutMethod)
  {
    Init();
    if (!fCutMethod) return kFALSE;
  }

  fCallParams[0] = p;

  fCutMethod->SetParamPtrs(&fCallParams[0]);
  Long_t result;
  fCutMethod->Execute(fCutObject,result);

  return (result!=0);
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutElement::CallCutMethod(Long_t p1, Long_t p2) const
{
  /// Call the cut method with two parameters

  if (!fCutMethod)
  {
    Init();
    if (!fCutMethod) return kFALSE;
  }

  fCallParams[0] = p1;
  fCallParams[1] = p2;

  fCutMethod->SetParamPtrs(&fCallParams[0]);
  Long_t result;
  fCutMethod->Execute(fCutObject,result);

  return (result!=0);
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuCutElement::CountOccurences(const TString& prototype, const char* search) const
{
  /// Count the number of times "search" is found in prototype

  TObjArray* a = prototype.Tokenize(",");
  TObjString* str;
  TIter next(a);

  Int_t n(0);

  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    if  ( str->String().Contains(search) )
    {
      ++n;
    }
  }
  delete a;
  return n;
}

//_____________________________________________________________________________
const char* AliAnalysisMuMuCutElement::CutTypeName(ECutType type)
{
  /// Convert the enum into a string

  if ( type == kEvent )
  {
    return "Event";
  }
  if ( type == kTrack )
  {
    return "Track";
  }
  if ( type == kTrackPair )
  {
    return "TrackPair";
  }
  if ( type == kTriggerClass )
  {
    return "TriggerClass";
  }
  return "Any";
}

//_____________________________________________________________________________
const char* AliAnalysisMuMuCutElement::GetCallMethodName() const
{
  /// Return the cut method name
  return ( fCutMethod ? fCutMethod->GetMethodName() : "");
}

//_____________________________________________________________________________
const char* AliAnalysisMuMuCutElement::GetCallMethodProto() const
{
  /// Return the cut method prototype
  return ( fCutMethod ? fCutMethod->GetProto() : "");
}

//_____________________________________________________________________________
void AliAnalysisMuMuCutElement::Init(ECutType expectedType) const
{
  /** Build the non-persistent members (TMethodCalls)
    *
    * Each cut method XXXX must have
    * a companion method, called NameOfXXXX(TString& nameOfCut) which fills the nameOfCut
    * string with the name of cut (using the default parameter values if need be to distinguish
    * different cuts using the same method).
    * The cut method we're proxy-ing can not be of any kind. We basically only recognize the
    * following prototypes (more or less based on the Pass() method prototypes). The checks
    * performed below to ensure that are far from being bullet-proof though... You've been
    * warned !
    *
    * For event cutters (3 or them instead of just one because of the, imo,
    * improper event/eventhandler interfaces we currently have in AliRoot)
    *
    * Bool_t XXXX(const AliVEvent&, ...)
    * Bool_t XXXX(const AliVEventHandler&, ...)
    * Bool_t XXXX(const AliInputEventHandler&, ...)
    *
    * For track cutters :
    *
    * Bool_t XXXX(const AliVParticle&, ...)
    *
    * For track pairs cutters :
    *
    * Bool_t XXXX(const AliVParticle&, const AliVParticle&, ...)
    *
    * For trigger class cutters :
    *
    * Bool_t XXXX(const TString& firedTriggerClasses, TString& acceptedTriggerClasses)
    *
    * where ... stands for optional arguments (different from VEvent VParticle, etc...),
    * which can have default values
    *
    * Note that Root reflexion does not allow (yet?) to check for constness of the arguments,
    * so AliVEvent& and const AliVEvent& will be the same.
    *
    *
   */

  TString scutMethodPrototype(fCutMethodPrototype);

  // some basic checks first

  TObjArray* tmp = fCutMethodPrototype.Tokenize(",");
  fNofParams = tmp->GetEntries();
  delete tmp;

  Int_t nVEvent = CountOccurences(fCutMethodPrototype,"AliVEvent");
  Int_t nVEventHandler = CountOccurences(fCutMethodPrototype,"AliVEventHandler") + CountOccurences(fCutMethodPrototype,"AliInputEventHandler");
  Int_t nparticles = CountOccurences(fCutMethodPrototype,"AliVParticle");
  Int_t nstrings = CountOccurences(fCutMethodPrototype,"TString");

  if ( expectedType == kEvent && ( nVEvent == 0 && nVEventHandler == 0 ) )
  {
    AliError(Form("Cut not of the expected %s type : did not find required prototype arguments AliVEvent, AliVEventHandler or AliInputEventHandler",CutTypeName(kEvent)));
    return;
  }

  if ( expectedType == kTrack && ( nparticles != 1 ) )
  {
    AliError(Form("Cut not of the expected %s type : did not find the required prototype argument AliVParticle (one and only one required)",CutTypeName(kTrack)));
    return;
  }

  if ( expectedType == kTrackPair && ( nparticles != 2 ) )
  {
    AliError(Form("Cut not of the expected %s type : did not find the required prototype arguments AliVParticle (2 of them required)",CutTypeName(kTrackPair)));
    return;
  }

  if ( expectedType == kTriggerClass && ( nstrings != 2 ) )
  {
    AliError(Form("Cut not of the expected %stype : did not find the required prototype arguments TString& (2 of them required)",CutTypeName(kTriggerClass)));
    return;
  }

  // OK, at least the prototype seems to match what we require. Let's continue...

  scutMethodPrototype.ReplaceAll("  ","");

  fCutMethod = new TMethodCall;

  fCutMethod->InitWithPrototype(fCutObject->IsA(),fCutMethodName.Data(),scutMethodPrototype.Data());

  if (!fCutMethod->IsValid())
  {
    AliError(Form("Could not find method %s(%s) in class %s",fCutMethodName.Data(),
                  scutMethodPrototype.Data(),fCutObject->ClassName()));
    delete fCutMethod;
    fCutMethod=0x0;
    return;
  }

  TMethodCall nameOfMethod;

  TString prototype("TString&");

  Int_t nMainPar = 0;

  if ( scutMethodPrototype.Contains("AliVEvent") )
  {
    fIsEventCutter=kTRUE;
    ++nMainPar;
  }
  if ( scutMethodPrototype.Contains("AliInputEventHandler") || scutMethodPrototype.Contains("AliVEventHandler") )
  {
    fIsEventHandlerCutter=kTRUE;
    ++nMainPar;
  }

  if ( nMainPar > 1 )
  {
    AliError(Form("Got an invalid prototype %s (more than one main parameter)",scutMethodPrototype.Data()));
    delete fCutMethod;
    fCutMethod=0x0;
    return;
  }

  if ( nparticles == 1 )
  {
    fIsTrackCutter=kTRUE;
  }
  else if ( nparticles == 2 )
  {
    fIsTrackPairCutter=kTRUE;
  }
  else if ( nstrings == 2 )
  {
    fIsTriggerClassCutter = kTRUE;
  }

  nMainPar += nparticles;
  nMainPar += nstrings;

  if ( nMainPar > 2 )
  {
    AliError(Form("Got an invalid prototype %s (more than two main parameters)",scutMethodPrototype.Data()));
    delete fCutMethod;
    fCutMethod=0x0;
    return;
  }

  if ( nMainPar == 0 )
  {
    AliError(Form("Got an invalid prototype %s (no main parameter found)",scutMethodPrototype.Data()));
    delete fCutMethod;
    fCutMethod=0x0;
    return;
  }

  if ( !fIsTriggerClassCutter )
  {
    scutMethodPrototype.ReplaceAll("const AliVEvent&","");
    scutMethodPrototype.ReplaceAll("AliVEvent&","");
    scutMethodPrototype.ReplaceAll("const AliVParticle&","");
    scutMethodPrototype.ReplaceAll("const AliInputEventHandler&","");
    scutMethodPrototype.ReplaceAll("const AliVEventHandler&","");

    prototype += scutMethodPrototype;

    nameOfMethod.InitWithPrototype(fCutObject->IsA(),Form("NameOf%s",fCutMethodName.Data()),prototype);

    if (!nameOfMethod.IsValid())
    {
      AliError(Form("Could not find method NameOf%s(%s) in class %s",fCutMethodName.Data(),
                    prototype.Data(),fCutObject->ClassName()));
      delete fCutMethod;
      fCutMethod=0x0;
      return;
    }

    // Now check if we have some default parameters for the NameOf method
    // Note that the only supported types for those default parameters
    // are Int_t and Double_t (which must be then of the form const Double_t&,
    // note the reference).

    prototype.ReplaceAll("TString&","");

    TObjArray* paramTypes = prototype.Tokenize(",");
    TObjArray* paramValues = fDefaultParameters.Tokenize(",");

    fDoubleParams.resize(paramValues->GetEntries());

    Int_t nparams = paramValues->GetEntries();

    // first parameter is always the TString&, i.e. the "output" of the NameOf
    // method

    fCallParams.resize(nparams+nMainPar);

    if ( nMainPar == 2 )
    {
      fCallParams[0] = 0;
      fCallParams[1] = reinterpret_cast<Long_t>(&fName);
    }
    else
    {
      fCallParams[0] = reinterpret_cast<Long_t>(&fName);
    }

    for ( Int_t i = 0; i < nparams; ++i )
    {
      TString pValue = static_cast<TObjString*>(paramValues->At(i))->String();
      TString pType = static_cast<TObjString*>(paramTypes->At(i))->String();

      if ( pType.Contains("Double_t"))
      {
        fDoubleParams[i] = pValue.Atof();
        fCallParams[i+nMainPar] = reinterpret_cast<Long_t>(&fDoubleParams[i]);
      }
      else if ( pType.Contains("Int_t") )
      {
        fCallParams[i+nMainPar] = pValue.Atoi();
      }
      else
      {
        AliError(Form("Got a parameter of type %s which I don't exactly know how to deal with. Expect something bad to happen...",pType.Data()));
        fCallParams[i+nMainPar] = reinterpret_cast<Long_t>(&pValue);
      }
    }

    nameOfMethod.SetParamPtrs(&fCallParams[0+nMainPar-1]);

    nameOfMethod.Execute(fCutObject);

    delete paramTypes;
    delete paramValues;
  }
  else
  {
      // check whether we have the input bit masks as well
    Int_t nuint = CountOccurences(scutMethodPrototype,"UInt_t");

    Bool_t ok =
     ( nuint == 3 && nstrings == 2 && ( nuint + nstrings ) == fNofParams ) ||
    ( nuint == 0 && nstrings == 2 && nstrings == fNofParams );

    if (!ok)
    {
      AliError("Incorrect prototype for a trigger class cutter");
      delete fCutMethod;
      fCutMethod=0x0;
      return;
    }
  }

  // Final consistency ross-check

  if ( expectedType == kEvent && ! (fIsEventCutter || fIsEventHandlerCutter) )
  {
    AliError("No, it's not an event cutter. Invalidate");
    delete fCutMethod;
    fCutMethod=0x0;
  }

  if ( expectedType == kTrack && !fIsTrackCutter )
  {
    AliError("No, it's not a track cutter. Invalidate");
    delete fCutMethod;
    fCutMethod=0x0;
  }

  if ( expectedType == kTrackPair && !fIsTrackPairCutter )
  {
    AliError("No, it's not a track pair cutter. Invalidate");
    delete fCutMethod;
    fCutMethod=0x0;
  }

  if ( expectedType == kTriggerClass && !fIsTriggerClassCutter )
  {
    AliError("No, it's not a trigger class cutter. Invalidate");
    delete fCutMethod;
    fCutMethod=0x0;
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutElement::IsEqual(const TObject* obj) const
{
  /// Whether we're the same cut as obj

  if ( obj->IsA() != IsA() ) return kFALSE;

  const AliAnalysisMuMuCutElement* cut = static_cast<const AliAnalysisMuMuCutElement*>(obj);

  return ( fName == cut->fName &&
    fIsEventCutter == cut->fIsEventCutter &&
    fIsEventHandlerCutter == cut->fIsEventHandlerCutter &&
    fIsTrackCutter == cut->fIsTrackCutter &&
    fIsTrackPairCutter == cut->fIsTrackPairCutter &&
    fIsTriggerClassCutter == cut->fIsTriggerClassCutter &&
    fCutMethodName == cut->fCutMethodName &&
    fCutMethodPrototype == cut->fCutMethodPrototype &&
    fDefaultParameters == cut->fDefaultParameters &&
    fCutObject == cut->fCutObject /* CHECK: should we really impose object equality or class equality would be enough ? */
  );
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutElement::Pass(const AliVEvent& event) const
{
  /// Whether the event pass this cut
  return CallCutMethod(reinterpret_cast<Long_t>(&event));
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutElement::Pass(const AliInputEventHandler& eventHandler) const
{
  /// Whether the eventHandler pass this cut
  return CallCutMethod(reinterpret_cast<Long_t>(&eventHandler));
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutElement::Pass(const AliVParticle& part) const
{
  /// Whether the particle pass this cut
  return CallCutMethod(reinterpret_cast<Long_t>(&part));
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutElement::Pass(const AliVParticle& p1, const AliVParticle& p2) const
{
  /// Whether the particle pair pass this cut
  return CallCutMethod(reinterpret_cast<Long_t>(&p1),reinterpret_cast<Long_t>(&p2));
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuCutElement::Pass(const TString& firedTriggerClasses,
                                       TString& acceptedTriggerClasses,
                                       UInt_t L0, UInt_t L1, UInt_t L2) const
{
  /** Whether the firedTriggerClasses pass the cut.
   * \param firedTriggerClasses (input) list of fired trigger classes (separated by space)
   * \param acceptedTriggerClasses (output) list of accepted classes
   * \param L0 (input, optional) level 0 trigger mask
   * \param L1 (input, optional) level 1 trigger mask
   * \param L2 (input, optional) level 2 trigger mask
   */

  if (!fCutMethod)
  {
    Init();
    if (!fCutMethod) return kFALSE;
  }

  acceptedTriggerClasses = "";

  Long_t result;
  Long_t params[] = { reinterpret_cast<Long_t>(&firedTriggerClasses),
    reinterpret_cast<Long_t>(&acceptedTriggerClasses),
    L0,L1,L2 };

  fCutMethod->SetParamPtrs(params);
  fCutMethod->Execute(fCutObject,result);
  return (result!=0);
}

//_____________________________________________________________________________
void AliAnalysisMuMuCutElement::Print(Option_t* opt) const
{
  /// Printout of the cut information
  if ( !fCutMethod )
  {
    Init();
  }

  TString sopt(opt);
  sopt.ToUpper();

  if (sopt.Contains("PTR"))
  {
    std::cout << Form("Cut %s(%p) %s(%p)::%s(%s) [",
                    fName.Data(),this,
                    fCutObject->ClassName(),
                    fCutObject,
                    GetCallMethodName(),
                    GetCallMethodProto());
  }
  else
  {
    std::cout << Form("Cut %s %s::%s(%s) [",
                      fName.Data(),
                      fCutObject->ClassName(),
                      GetCallMethodName(),
                      GetCallMethodProto());
  }

  if ( IsEventCutter() ) std::cout << " E";
  if ( IsEventHandlerCutter() ) std::cout << " EH";
  if ( IsTrackCutter() ) std::cout << " T";
  if ( IsTrackPairCutter() ) std::cout << " TP";
  if ( IsTriggerClassCutter() ) std::cout << " TC";

  std::cout << " ]" << std::endl;
}

//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________

/**
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuCutElementBar
 *
 * \brief The negation of an AliAnalysisMuMuCutElement
 *
 *  \author L. Aphecetche (Subatech)
 */

//_____________________________________________________________________________
AliAnalysisMuMuCutElementBar::AliAnalysisMuMuCutElementBar() : AliAnalysisMuMuCutElement(),
fCutElement(0x0)
{
  /// default ctor
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElementBar::AliAnalysisMuMuCutElementBar(const AliAnalysisMuMuCutElement& ce)
: AliAnalysisMuMuCutElement(), fCutElement(&ce)
{
  /// ctor
  fIsEventCutter = ce.IsEventCutter();
  fIsEventHandlerCutter = ce.IsEventHandlerCutter();
  fIsTrackCutter = ce.IsTrackCutter();
  fIsTrackPairCutter = ce.IsTrackPairCutter();
  fName = TString::Format("NOT%s",ce.GetName());
}

//_____________________________________________________________________________
AliAnalysisMuMuCutElementBar::~AliAnalysisMuMuCutElementBar()
{
  /// dtor (nop as we're not the owner of fCutElement)
}

//_____________________________________________________________________________
void AliAnalysisMuMuCutElementBar::Print(Option_t* /*opt*/) const
{
  /// Printout of the cut information
  std::cout << Form("Cut %s(%p) : negation of %s(%p)",GetName(),this,fCutElement->GetName(),fCutElement)
  << std::endl;
}

