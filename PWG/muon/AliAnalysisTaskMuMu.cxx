#include "AliAnalysisTaskMuMu.h"

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliCounterCollection.h"
#include "AliHistogramCollection.h"
#include "AliInputEventHandler.h"
#include "AliLog.h" 
#include "AliMuonTrackCuts.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TFormula.h"
#include "TH1.h"
#include "TH2.h"
#include "THashList.h"
#include "TList.h"
#include "TMath.h"
#include "TObjString.h"
#include "TPaveText.h"
#include "TROOT.h"
#include <algorithm>
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliAODTrack.h"
#include "AliAODTZERO.h"
#include "AliESDTZERO.h"
#include "AliCodeTimer.h"

//
// AliAnalysisTaskMuMu : base class for mu pairs analysis 
//
// Mainly invariant mass (for J/psi and Upsilon) but also 
// some single control histograms.
//
// This base class contains common things for ESD-based
// and AOD-based analysis
//
// The output contains an AliHistogramCollection and
// an AliCounterCollection
//
// author: L. Aphecetche (Subatech)
//

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskMuMu)
ClassImp(AliAnalysisTaskMuMu::PairCut)

namespace
{
  Int_t GetNbins(Double_t xmin, Double_t xmax, Double_t xstep)
  {
    if ( TMath::AreEqualRel(xstep,0.0,1E-9) ) return 1;
    
    return TMath::Nint(TMath::Abs((xmax-xmin)/xstep));
  }
  
  TObjArray* GetMuonTriggerList()
  {
    TObjArray* a = new TObjArray;
    a->SetOwner(kTRUE);
    
    a->Add(new TObjString("CMUL"));
    a->Add(new TObjString("CMLL"));
    a->Add(new TObjString("C0MUL"));
    a->Add(new TObjString("CMSL"));
    a->Add(new TObjString("CMSH"));
    
    return a;
  }

  TObjArray* GetEmcalTriggerList()
  {
    TObjArray* a = new TObjArray;
    a->SetOwner(kTRUE);
    
    a->Add(new TObjString("CEMC"));
    
    return a;
  }

  TObjArray* GetMBTriggerList()
  {
    TObjArray* a = new TObjArray;
    a->SetOwner(kTRUE);
    
    a->Add(new TObjString("CINT7-S-"));
    a->Add(new TObjString("CINT8-S-"));
    
    return a;
  }

}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::PairCut::Print(Option_t* /*opt*/) const
{
  std::cout << Form("PAIR   CUT %20s SINGLE MASK %x PAIR MASK %x",GetName(),MaskForOneOrBothTrack(),MaskForTrackPair())
  << std::endl;
}

//_____________________________________________________________________________
AliAnalysisTaskMuMu::AliAnalysisTaskMuMu() : AliAnalysisTaskSE("AliAnalysisTaskMuMu"),
fOutput(0),
fHistogramCollection(0),
fEventCounters(0),
fMuonTrackCuts(0x0),
fPrecomputedTrackMasks(),
fIsFromESD(kFALSE),
fIsDynamicTriggerClasses(kFALSE),
fShouldSeparatePlusAndMinus(kFALSE),
fIsHistogrammingDisabled(kFALSE),
fBeamYear("pp"),
fCentralityLimits(),
fTriggerClasses(0),
fSingleTrackCutNames(0x0),
fPairTrackCutNames(0x0),
fCentralityNames(0x0),
fEventCutNames(0x0),
fUseBackgroundTriggers(kFALSE),
fTriggerInputBitMap()
{
  /// default ctor
}

//_____________________________________________________________________________
AliAnalysisTaskMuMu::AliAnalysisTaskMuMu(Bool_t fromESD, const char* beamYear, TArrayF* centralities)
: AliAnalysisTaskSE(Form("AliAnalysisTaskMuMu-from%s",fromESD ? "ESD":"AOD")),
fOutput(0),
fHistogramCollection(0),
fEventCounters(0),
fMuonTrackCuts(0x0),
fPrecomputedTrackMasks(),
fIsFromESD(fromESD),
fIsDynamicTriggerClasses(kTRUE),
fShouldSeparatePlusAndMinus(kFALSE),
fIsHistogrammingDisabled(kFALSE),
fBeamYear(beamYear),
fCentralityLimits(),
fTriggerClasses(new THashList),
fSingleTrackCutNames(0x0),
fPairTrackCutNames(0x0),
fCentralityNames(new TObjArray),
fEventCutNames(0x0),
fUseBackgroundTriggers(kFALSE),
fTriggerInputBitMap()
{
  /// Constructor
  /// The list of triggers to be considered will be updated on the fly
  /// (see method AddTriggerClasses)
  
  fTriggerClasses->SetOwner(kTRUE);
  
  DefineOutput(1, TList::Class());
  
  DefineCentralityClasses(centralities);  
}

//_____________________________________________________________________________
AliAnalysisTaskMuMu::AliAnalysisTaskMuMu(Bool_t fromESD, TList* triggerClasses, const char* beamYear, TArrayF* centralities)
: AliAnalysisTaskSE(Form("AliAnalysisTaskMuMu-from%s",fromESD ? "ESD":"AOD")),
fOutput(0),
fHistogramCollection(0),
fEventCounters(0),
fMuonTrackCuts(0x0),
fPrecomputedTrackMasks(),
fIsFromESD(fromESD),
fIsDynamicTriggerClasses(kFALSE),
fShouldSeparatePlusAndMinus(kFALSE),
fIsHistogrammingDisabled(kFALSE),
fBeamYear(beamYear),
fCentralityLimits(),
fTriggerClasses(new THashList),
fSingleTrackCutNames(0x0),
fPairTrackCutNames(0x0),
fCentralityNames(new TObjArray),
fEventCutNames(0x0),
fUseBackgroundTriggers(kFALSE),
fTriggerInputBitMap()
{
  /// Constructor with a predefined list of triggers to consider

  fTriggerClasses->SetOwner(kTRUE);
  
  DefineOutput(1, TList::Class());
  
  TObjString* tname;
  TIter next(triggerClasses);
  
  while ( ( tname = static_cast<TObjString*>(next()) ) )
  {
    fTriggerClasses->Add(new TObjString(*tname));
    if ( !beamYear ) 
    {
      if ( tname->String().BeginsWith("CMB") )
      {
        fBeamYear = "PbPb2010";
      }
      if ( tname->String().BeginsWith("CPBI") )
      {
        fBeamYear = "PbPb2011";
      }
    }
  }
  
  DefineCentralityClasses(centralities);  
}

//_____________________________________________________________________________
AliAnalysisTaskMuMu::~AliAnalysisTaskMuMu()
{
  /// dtor

  if (fOutput && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode()) 
  {
    delete fOutput;     
  }

  delete fMuonTrackCuts;

  delete fTriggerClasses;
  delete fSingleTrackCutNames;
  delete fPairTrackCutNames;
  delete fCentralityNames;
  delete fEventCutNames;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::AddTriggerClasses(const char* triggerlist, const char* sep)
{
  /// Given a list of trigger names (triggerlist) separated by sep
  /// add those which we don't know yet (selecting only the few ones
  /// we're interested in, e.g CINT* CMU* CMB*
  ///
  
  TString slist(triggerlist);
  
  TObjArray* a = slist.Tokenize(sep);
  TObjString* s;
  TIter next(a);
  
  while ( ( s = static_cast<TObjString*>(next()) ) )
  {
    TString trigger(s->String());
    Bool_t add(kFALSE);

    if ( trigger.Contains("WU") )
    {
      add = kFALSE; 
      continue;
    }

    // MB triggers
    if ( trigger.BeginsWith("CMB") && trigger.Contains("-B-") ) add = kTRUE; // MB PbPb 2010
    
    if ( trigger.BeginsWith("CPBI") && trigger.Contains("-B-") ) add = kTRUE; // MB PbPb 2011
    if ( trigger.BeginsWith("CVHN") && trigger.Contains("-B-") ) add = kTRUE; // Central PbPb 2011
    if ( trigger.BeginsWith("CVLN") && trigger.Contains("-B-") ) add = kTRUE; // Semi Central PbPb 2011

    if ( ( trigger.BeginsWith("CINT") /* || trigger.BeginsWith("CTRUE") */ ) && TriggerSBACECondition(trigger) )  add = kTRUE;
    
    // main muon triggers
    if ( ( trigger.BeginsWith("CMU") || trigger.BeginsWith("CMS") || trigger.BeginsWith("CML") ) &&
         TriggerSBACECondition(trigger) ) add = kTRUE;
    
    // muon triggers not in coincidence with anything else
    if ( trigger.BeginsWith("C0MUL") ) add = kTRUE;
    
    /*
    if ( trigger.BeginsWith("CMUP") ||
         trigger.BeginsWith("CCUP") )
    {
      // discard ultra-peripheral triggers
      add = kFALSE;
    }
    */
    
    // triggers in Monte-Carlo
    if ( trigger.BeginsWith("MB") ||
         trigger.BeginsWith("SC") ||
         ( trigger.BeginsWith("CE") && !trigger.BeginsWith("CEMC") ) ||
         trigger.BeginsWith("DMU") ||
         trigger.BeginsWith("DML") ||
         trigger.BeginsWith("MU") ||
        trigger.BeginsWith("CMSNGL") ||
        trigger.BeginsWith("CMLK") )
    {
      add = kTRUE; // for MC
    }
    
    if ( add && !fTriggerClasses->FindObject(trigger.Data()) )
    {
      AliDebug(1,Form("Adding %s to considered trigger classes",trigger.Data()));
      fTriggerClasses->Add(new TObjString(trigger));
    }
    
  }
  
  delete a;  
}
                                  
//_____________________________________________________________________________
void AliAnalysisTaskMuMu::AddPairCut(const char* cutName, UInt_t maskForOneOrBothTrack, UInt_t maskForTrackPair)
{
  /// Add a cut for a pair.
  /// maskForOneOrBothTrack is the mask of cuts that at least one track must satisfy
  /// maskForTrackPair is the mask of cuts that *both* tracks must satisfy.
  /// if maskForTrackPair is 0, then no (extra) condition is applied to the pair
  
  if ( !fPairTrackCutNames ) 
  {
    fPairTrackCutNames = new TObjArray;
    fPairTrackCutNames->SetOwner(kTRUE);
  }
  fPairTrackCutNames->Add(new AliAnalysisTaskMuMu::PairCut(cutName,maskForOneOrBothTrack,maskForTrackPair));
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::AddSingleCut(const char* name, UInt_t mask)
{
  /// Add a cut for single tracks
  if ( !fSingleTrackCutNames ) 
  {
    fSingleTrackCutNames = new TObjArray;
    fSingleTrackCutNames->SetOwner(kTRUE);
  }
  TObjString* oname = new TObjString(Form("s%s",name));
  oname->SetUniqueID(mask);
  fSingleTrackCutNames->Add(oname);
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::AddEventCut(const char* name, UInt_t mask)
{
  /// Add a cut at event level
  if ( !fEventCutNames )
  {
    fEventCutNames = new TObjArray;
    fEventCutNames->SetOwner(kTRUE);
  }
  TObjString* oname = new TObjString(Form("%s",name));
  oname->SetUniqueID(mask);
  fEventCutNames->Add(oname);
}


//_____________________________________________________________________________
void AliAnalysisTaskMuMu::AssertHistogramCollection(const char* physics, const char* triggerClassName)
{
  // insure that a given set of histogram is created
  if (!fHistogramCollection->Histo(physics,triggerClassName,DefaultCentralityName(),"Zvertex")) 
  {
    FillHistogramCollection(physics,triggerClassName);
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::BeautifyHistos()
{
  /// Put the titles, marker sizes, color, etc...
  
  TIter next(fHistogramCollection->CreateIterator());
  TH1* h;
  
  while ( ( h = static_cast<TH1*>(next()) ) )
  {
    TString name(h->GetName());
    
    if ( name.Contains("Plus") ) 
    {
      h->SetMarkerStyle(kFullCircle);
      h->SetMarkerColor(kBlue);
      h->SetLineColor(kBlue);        
    }
    else
    {
      h->SetMarkerStyle(kOpenCircle);
      h->SetMarkerColor(kRed);
      h->SetLineColor(kRed);        
    }      
  }
}

//_____________________________________________________________________________
const char* 
AliAnalysisTaskMuMu::CentralityName(Double_t centrality) const
{
  /// Get centrality name corresponding to the floating ^point value
  
  if ( centrality > 0 && centrality <= 100.0 ) 
  {
    return Form("CENT%02d",TMath::Nint(centrality));
  }
  else
  {
    return DefaultCentralityName();
  }
}

//_____________________________________________________________________________
Bool_t 
AliAnalysisTaskMuMu::CheckTriggerClass(const TString& toCheck,
                                       const TString& firedTriggerClasses,
                                       UInt_t l0Inputs) const
{
  // Check if the "toCheck" class (or logical combination of classes and L0 inputs)
  // are within the "firedTriggerClasses"
  
  Bool_t ok(kFALSE);
  
  TString tn(toCheck);
  TString comp("");
  
  if ( tn.Contains("(") || tn.Contains(")") || tn.Contains("&") || tn.Contains("|") )
  {
    comp=tn;
    tn.ReplaceAll("(","");
    tn.ReplaceAll(")","");
    tn.ReplaceAll("&","");
    tn.ReplaceAll("|","");
    TObjArray* a = tn.Tokenize(" ");
    TIter nextA(a);
    TObjString* an;
    while ( ( an = static_cast<TObjString*>(nextA()) ) )
    {
      //        AliDebug(1,Form("    an %s => %d",an->String().Data(),firedTriggerClasses.Contains(an->String().Data())));
      if ( an->String().BeginsWith("0") )
      {
        // that's an input
        UInt_t bit = GetTriggerInputBitMaskFromInputName(an->String().Data());
        comp.ReplaceAll(an->String().Data(),( (l0Inputs & (bit)) == bit) ? "1" : "0");
      }
      else
      {
        comp.ReplaceAll(an->String().Data(),Form("%d",firedTriggerClasses.Contains(an->String().Data())));
      }
    }
    delete a;
    
    TFormula formula("TriggerClassFormulaCheck", comp.Data());
    if (formula.Compile() > 0)
    {
      AliError(Form("Could not evaluate formula %s",comp.Data()));
    }
    else
    {
      ok = formula.Eval(0);
    }
  }
  else
  {
    ok = firedTriggerClasses.Contains(toCheck);
  }
  
  AliDebug(1,Form("tname %s => %d comp=%s",toCheck.Data(),ok,comp.Data()));
  
  return ok;
}

//_____________________________________________________________________________
void 
AliAnalysisTaskMuMu::CreateSingleHisto(const char* physics,
                                   const char* triggerClassName,
                                   const char* hname, const char* htitle, 
                                   Int_t nbinsx, Double_t xmin, Double_t xmax,
                                   Int_t nbinsy, Double_t ymin, Double_t ymax,
                                   Bool_t separatePlusAndMinus) const
{  
  /// Append histograms for single track to our histogram collection
  
  if ( separatePlusAndMinus ) 
  {
    const char* suffix[] = { "Plus", "Minus" };
    const char* symbol[] = { "+", "-" };

    for ( Int_t i = 0; i < 2; ++i )
    {
      TString shtitle(htitle);
      TString shname(hname);
      
      shtitle.ReplaceAll("#mu",Form("#mu^{%s}",symbol[i]));
      
      shname += suffix[i];
      
      CreateHisto(fSingleTrackCutNames,physics,triggerClassName,shname.Data(),shtitle.Data(),
                  nbinsx,xmin,xmax,nbinsy,ymin,ymax);
    }
  }
  else 
  {
    CreateHisto(fSingleTrackCutNames,physics,triggerClassName,hname,htitle,
                nbinsx,xmin,xmax,nbinsy,ymin,ymax);
  }
}

//_____________________________________________________________________________
void 
AliAnalysisTaskMuMu::CreatePairHisto(const char* physics,
                                 const char* triggerClassName,
                                 const char* hname, const char* htitle, 
                                 Int_t nbinsx, Double_t xmin, Double_t xmax,
                                 Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  /// Append histograms for track pairs to our histogram collection

  CreateHisto(fPairTrackCutNames,physics,triggerClassName,hname,htitle,
              nbinsx,xmin,xmax,nbinsy,ymin,ymax);
}

//_____________________________________________________________________________
void 
AliAnalysisTaskMuMu::CreateEventHisto(const char* physics,
                                  const char* triggerClassName,
                                  const char* hname, const char* htitle, 
                                  Int_t nbinsx, Double_t xmin, Double_t xmax,
                                  Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  /// Append histograms at the event level
  
  TIter next(fCentralityNames);
  TObjString* cent;
  
  while ( ( cent = static_cast<TObjString*>(next()) ) )
  {
    TH1* h(0x0);
    
    if ( nbinsy > 0 )
    {  
      h = new TH2F(hname,htitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
    }
    else
    {
      h = new TH1F(hname,htitle,nbinsx,xmin,xmax);
    }
    
    fHistogramCollection->Adopt(physics,triggerClassName,cent->String().Data(),h); 
  }
}

//_____________________________________________________________________________
void 
AliAnalysisTaskMuMu::CreateHisto(TObjArray* array,
                             const char* physics,
                             const char* triggerClassName,
                             const char* hname, const char* htitle, 
                             Int_t nbinsx, Double_t xmin, Double_t xmax,
                             Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  /// Create a bunch of histograms for all centralities
  /// FIXME: have a way to specify the histo precision (i.e. F vs I vs S ...)
  
  TIter next(array);
  TObjString* tcn;
  while ( ( tcn = static_cast<TObjString*>(next()) ) )
  {
    TIter nextCent(fCentralityNames);
    TObjString* cent;
    
    while ( ( cent = static_cast<TObjString*>(nextCent()) ) )
    {
      TH1* h(0x0);
    
      if ( nbinsy > 0 )
      {  
        h = new TH2F(hname,htitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
      }
      else
      {
        h = new TH1F(hname,htitle,nbinsx,xmin,xmax);
      }
    
      fHistogramCollection->Adopt(physics,triggerClassName,cent->String().Data(),tcn->String().Data(),h);
    }
  }      
}

//_____________________________________________________________________________
const char* 
AliAnalysisTaskMuMu::DefaultCentralityName() const
{
  /// Get default centrality name
  if ( !fBeamYear.Contains("pp") ) return "CENTX";
  else return "PP";
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::DefineCentralityClasses(TArrayF* centralities)
{
  /// Define the default centrality classes that will be used.
  
  if ( !fBeamYear.Contains("pp") ) 
  {
    if ( !centralities ) 
    {
      // default values
      fCentralityLimits.push_back(10.0);
      fCentralityLimits.push_back(30.0);
      fCentralityLimits.push_back(50.0);
      fCentralityLimits.push_back(80.0);
    }
    else
    {
      for ( Int_t i = 0; i < centralities->GetSize(); ++i ) 
      {
        fCentralityLimits.push_back(centralities->At(i));
      }
    }
  }
  
  for ( std::vector<double>::size_type i = 0; i < fCentralityLimits.size(); ++i )
  {
    Double_t limit = fCentralityLimits[i];
    fCentralityNames->Add(new TObjString(CentralityName(limit)));
  }
  
  fCentralityNames->Add(new TObjString(DefaultCentralityName()));
  fCentralityNames->SetOwner(kTRUE);
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::EAComputeTrackMasks(const AliVEvent& event)
{
  // compute the track masks for the event
  
  fPrecomputedTrackMasks.Reset();
  Int_t n = EAGetNumberOfMuonTracks(event);
  fPrecomputedTrackMasks.Set(n);
  
  if ( event.IsA() == AliESDEvent::Class() ) 
  {
    const AliESDEvent& esd = static_cast<const AliESDEvent&>(event);
    AliESDEvent& ncesd = const_cast<AliESDEvent&>(esd);
    
    for ( Int_t i = 0; i < n; ++i ) 
    {
      ComputeTrackMask(*(ncesd.GetMuonTrack(i)),i);
    }
  }
  else if ( event.IsA() == AliAODEvent::Class() ) 
  {
    const AliAODEvent& aod = static_cast<const AliAODEvent&>(event);
    for ( Int_t i = 0; i < n; ++i ) 
    {
      ComputeTrackMask(*(aod.GetTrack(i)),i);
    }
  }

}

//_____________________________________________________________________________
TString AliAnalysisTaskMuMu::EAGetFiredTriggerClasses(const AliVEvent& event) const
{
  // hack (this method should really be part of AliVEvent itself) FIXME
  
  if ( event.IsA() == AliESDEvent::Class() ) 
  {
    return static_cast<const AliESDEvent&>(event).GetFiredTriggerClasses();
  }
  else if ( event.IsA() == AliAODEvent::Class() ) 
  {
    return static_cast<const AliAODEvent&>(event).GetFiredTriggerClasses();    
  }
  else
  {
    AliError(Form("Unknown class for the event = %s",event.ClassName()));
    return "";
  }
}

//_____________________________________________________________________________
UInt_t AliAnalysisTaskMuMu::EAGetL0TriggerInputs(const AliVEvent& event) const
{
  // Get the list of level 0 trigger inputs
  
  if ( event.IsA() == AliESDEvent::Class() ) 
  {
    return static_cast<const AliESDEvent&>(event).GetHeader()->GetL0TriggerInputs();
  }
  else if ( event.IsA() == AliAODEvent::Class() ) 
  {
    return static_cast<const AliAODEvent&>(event).GetHeader()->GetL0TriggerInputs();
  }
  else
  {
    AliError(Form("Unknown class for the event = %s",event.ClassName()));
    return 0;
  }
  
  return 0;
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskMuMu::EAGetNumberOfMuonTracks(const AliVEvent& event) const
{
  // get the number of muon tracks from the event
  if ( event.IsA() == AliESDEvent::Class() ) 
  {
    Int_t n(0);
    const AliESDEvent& esd = static_cast<const AliESDEvent&>(event);
    for ( Int_t i = 0; i < esd.GetNumberOfMuonTracks(); ++i )
    {
      AliESDMuonTrack* muon = const_cast<AliESDEvent&>(esd).GetMuonTrack(i);
      if ( muon->ContainTrackerData() )
      {
        ++n;
      }
    }
    return n;
  }
  else if ( event.IsA() == AliAODEvent::Class() ) 
  {
    return static_cast<const AliAODEvent&>(event).GetNumberOfTracks();    
  }
  else
  {
    AliError(Form("Unknown class for the event = %s",event.ClassName()));
    return 0;
  }
}

//_____________________________________________________________________________
AliVParticle* AliAnalysisTaskMuMu::EAGetTrack(const AliVEvent& event, Int_t i) const
{
  // get the i-th track from the event
  if ( event.IsA() == AliESDEvent::Class() ) 
  {
    AliESDMuonTrack* track = static_cast<AliESDEvent&>(const_cast<AliVEvent&>(event)).GetMuonTrack(i);
    if ( track->ContainTrackerData() ) 
    {
      return track;
    }
    return 0x0;
  }
  else if ( event.IsA() == AliAODEvent::Class() ) 
  {
    AliAODTrack* part = static_cast<const AliAODEvent&>(event).GetTrack(i);
    if ( part->IsMuonTrack() )
    {
      return part;
    }
    return 0x0;
  }
  else
  {
    AliError(Form("Unknown class for the event = %s",event.ClassName()));
  }  
  return 0x0;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskMuMu::EAGetTrackDCA(const AliVParticle& track) const
{
  // Get track DCA
  
  Double_t xdca(1E9);
  Double_t ydca(xdca);
  
  if ( track.IsA() == AliAODTrack::Class() ) 
  {
    xdca = static_cast<const AliAODTrack&>(track).XAtDCA();
    ydca = static_cast<const AliAODTrack&>(track).YAtDCA();
  }
  else if ( track.IsA() == AliESDMuonTrack::Class() ) 
  {
    xdca = static_cast<const AliESDMuonTrack&>(track).GetNonBendingCoorAtDCA(); 
    ydca = static_cast<const AliESDMuonTrack&>(track).GetBendingCoorAtDCA();    
  }
  else
  {
    AliError(Form("Unknown track class: %s",track.ClassName()));
  }
  
  return TMath::Sqrt(xdca*xdca+ydca*ydca);
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskMuMu::EAGetTrackNormalizedChi2(const AliVParticle& track) const
{
  // get the chi2 per NDF of the track
  
  Double_t chi2(1E9);
  
  if ( track.IsA() == AliAODTrack::Class() ) 
  {
    chi2 = static_cast<const AliAODTrack&>(track).Chi2perNDF();
  }
  else if ( track.IsA() == AliESDMuonTrack::Class() ) 
  {
    chi2 = static_cast<const AliESDMuonTrack&>(track).GetNormalizedChi2();
  }
  else
  {
    AliError(Form("Unknown track class: %s",track.ClassName()));
  }
  
  return chi2;

}

//_____________________________________________________________________________
Double_t AliAnalysisTaskMuMu::EAGetTrackChi2MatchTrigger(const AliVParticle& track) const
{
  /// Get the Chi2 of the matching MCH - MTR
  
  if ( track.IsA() == AliAODTrack::Class() ) 
  {
    return static_cast<const AliAODTrack&>(track).GetChi2MatchTrigger();
  }
  else
  {
    return static_cast<const AliESDMuonTrack&>(track).GetChi2MatchTrigger();
  }
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskMuMu::GetTrackTheta(const AliVParticle& track) const
{
  // Get track theta (in radian)
  
  return TMath::ATan(EAGetTrackRabs(track)/AbsZEnd());
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskMuMu::EAGetTrackRabs(const AliVParticle& track) const
{
  // Get track DCA
  
  Double_t rabs(1E9);
  
  if ( track.IsA() == AliAODTrack::Class() ) 
  {
    rabs = static_cast<const AliAODTrack&>(track).GetRAtAbsorberEnd();
  }
  else if ( track.IsA() == AliESDMuonTrack::Class() ) 
  {
    rabs = static_cast<const AliESDMuonTrack&>(track).GetRAtAbsorberEnd(); 
  }
  else
  {
    AliError(Form("Unknown track class: %s",track.ClassName()));
  }
  
  return TMath::ATan(rabs)/AbsZEnd();
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::EAGetTZEROFlags(const AliVEvent& event, Bool_t& backgroundFlag, Bool_t& pileupFlag, Bool_t& satelliteFlag) const
{
  // get the TZERO decisions
  // return false if there's no tzero information in this event
  
  Bool_t rv(kFALSE);
  
  if ( event.IsA() == AliESDEvent::Class() ) 
  {
    const AliESDTZERO* tzero = static_cast<AliESDEvent&>(const_cast<AliVEvent&>(event)).GetESDTZERO();
    if ( tzero )
    {
      backgroundFlag = tzero->GetBackgroundFlag();
      pileupFlag = tzero->GetPileupFlag();
      satelliteFlag = tzero->GetSatellite();
      rv = kTRUE;
    }    
  }
  else if ( event.IsA() == AliAODEvent::Class() ) 
  {
    AliAODTZERO* tzero = static_cast<const AliAODEvent&>(event).GetTZEROData();
    if ( tzero )
    {
      backgroundFlag = tzero->GetBackgroundFlag();
      pileupFlag = tzero->GetPileupFlag();
      satelliteFlag = tzero->GetSatellite();
      rv = kTRUE;
    }
  }
  else
  {
    AliError(Form("Unknown class for the event = %s",event.ClassName()));
  }  
  
  return rv;

}

//_____________________________________________________________________________
Bool_t
AliAnalysisTaskMuMu::AtLeastOneMBTrigger(const TString& firedTriggerClasses) const
{
  // whether or not we have a least one MB trigger in the fired trigger classes
  static TObjArray* triggerList = GetMBTriggerList();
  TIter next(triggerList);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next())))
  {
    if ( firedTriggerClasses.Contains(str->String().Data())) return kTRUE;
  }
  
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t
AliAnalysisTaskMuMu::AtLeastOneMuonTrigger(const TString& firedTriggerClasses) const
{
  // whether or not we have a least one muon trigger in the fired trigger classes
  static TObjArray* triggerList = GetMuonTriggerList();
  TIter next(triggerList);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next())))
  {
    if ( firedTriggerClasses.Contains(str->String().Data())) return kTRUE;
  }
  
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t
AliAnalysisTaskMuMu::AtLeastOneEmcalTrigger(const TString& firedTriggerClasses) const
{
  // whether or not we have a least one emcal trigger in the fired trigger classes
  static TObjArray* triggerList = GetEmcalTriggerList();
  TIter next(triggerList);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next())))
  {
    if ( firedTriggerClasses.Contains(str->String().Data())) return kTRUE;
  }
  
  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::Fill(const char* eventtype, 
                           TObjString* tname, 
                           const char* centrality, 
                           float fcent, 
                           const AliVEvent& event)
{
  // Fill one set of histograms
  
  TString seventtype(eventtype);
  seventtype.ToLower();
  
  fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", seventtype.Data(), tname->GetName(), fCurrentRunNumber));

  if ( !IsHistogrammingDisabled() )
  {
    AssertHistogramCollection(eventtype,tname->String().Data());
    FillHistos(eventtype,tname->String().Data(),centrality,event);    
    fHistogramCollection->Histo(eventtype,tname->String().Data(),"Centrality")->Fill(fcent);
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillHistogramCollection(const char* physics, const char* triggerClassName)
{
  /// Actually create the histograms for phyics/triggerClassName
  
  AliDebug(1,Form("(%s,%s)",physics,triggerClassName));
  
  Double_t ptMin = 0;
  Double_t ptMax = 12*3;
  Int_t nbinsPt = GetNbins(ptMin,ptMax,0.5);
  Double_t pMin = 0;
  Double_t pMax = 100*3;
  Int_t nbinsP = GetNbins(pMin,pMax,2.0);
  Double_t etaMin = -5;
  Double_t etaMax = -2;
  Int_t nbinsEta = GetNbins(etaMin,etaMax,0.05);
  
  Double_t rapidityMin = -5;
  Double_t rapidityMax = -2;
  Int_t nbinsRapidity = GetNbins(rapidityMin,rapidityMax,0.05);
  
  CreateSingleHisto(physics,triggerClassName,"Chi2MatchTrigger","Chi2 Match Trigger",72,0,72);

  CreateSingleHisto(physics,triggerClassName,"EtaRapidityMu", "Eta distribution vs Rapidity for #mu", nbinsRapidity,rapidityMin,rapidityMax,nbinsEta,etaMin,etaMax, fShouldSeparatePlusAndMinus);

  CreateSingleHisto(physics,triggerClassName,"PtEtaMu", "P_{T} distribution vs Eta for #mu", nbinsEta,etaMin,etaMax, nbinsPt,ptMin,ptMax,fShouldSeparatePlusAndMinus);

  CreateSingleHisto(physics,triggerClassName,"PtRapidityMu", "P_{T} distribution vs Rapidity for #mu", nbinsRapidity,rapidityMin,rapidityMax, nbinsPt,ptMin,ptMax,fShouldSeparatePlusAndMinus);

  CreateSingleHisto(physics,triggerClassName,"PEtaMu", "P distribution for #mu",nbinsEta,etaMin,etaMax,nbinsP,pMin,pMax,fShouldSeparatePlusAndMinus);

  Double_t chi2min = 0;
  Double_t chi2max = 20;
  Int_t nbinchi2 = GetNbins(chi2min,chi2max,0.05);
  
  CreateSingleHisto(physics, triggerClassName, "Chi2Mu", "chisquare per NDF #mu", nbinchi2, chi2min, chi2max,fShouldSeparatePlusAndMinus);
    
  Double_t minvMin = 0;
  Double_t minvMax = 16;
  Int_t nMinvBins = GetNbins(minvMin,minvMax,0.025);
  
  CreatePairHisto(physics,triggerClassName,"Chi12","Chi2MatchTrigger of muon 1 vs muon 2",72,0,72,72,0,72);
  CreatePairHisto(physics,triggerClassName,"Rabs12","Rabs of muon 1 vs muon ",100,0,100,100,0,100);
  
  CreatePairHisto(physics,triggerClassName,"MinvUSPt", "#mu+#mu- inv. mass vs Pt",nbinsPt,ptMin,ptMax,nMinvBins,minvMin,minvMax);
  CreatePairHisto(physics,triggerClassName,"MinvUSRapidity", "#mu+#mu- inv. mass vs rapidity",nbinsRapidity,rapidityMin,rapidityMax,nMinvBins,minvMin,minvMax);
  CreatePairHisto(physics,triggerClassName,"USPtRapidity", "#mu+#mu- Pt vs Rapidity",nbinsRapidity,rapidityMin,rapidityMax,nbinsPt,ptMin,ptMax);
  
  CreatePairHisto(physics,triggerClassName,"MinvLSPt", "unlike sign inv. mass vs Pt",nbinsPt,ptMin,ptMax,nMinvBins,minvMin,minvMax);

  Double_t xmin = -35;
  Double_t xmax = +35;
  Int_t nbins = GetNbins(xmin,xmax,0.4);
  
  CreateEventHisto(physics,triggerClassName,"Zvertex","z vertex",nbins,xmin,xmax);  

  CreateEventHisto(physics,triggerClassName,"T0Zvertex","T0 zvertex",nbins,xmin,xmax);

  xmin = -5;
  xmax = 5;
  nbins = GetNbins(xmin,xmax,0.01);

  CreateEventHisto(physics,triggerClassName,"Xvertex","x vertex",nbins,xmin,xmax);  
  CreateEventHisto(physics,triggerClassName,"Yvertex","y vertex",nbins,xmin,xmax);  
  
//  CreateEventHisto(physics,triggerClassName,"YXvertex","y vs x vertex",nbins,xmin,xmax,nbins,xmin,xmax);  
  
//  if (!fIsFromESD)
//  {
//    
//    CreateEventHisto(physics,triggerClassName,"PileUpZvertex","pileup z vertex",nbins,xmin,xmax);  
//    
//    CreateEventHisto(physics,triggerClassName,"PileUpXvertex","pileup x vertex",nbins,xmin,xmax);  
//    CreateEventHisto(physics,triggerClassName,"PileUpYvertex","pileup y vertex",nbins,xmin,xmax);  
//    
//    CreateEventHisto(physics,triggerClassName,"PileUpYXvertex","pileup y vs x vertex",nbins,xmin,xmax,
//                     nbins,xmin,xmax);  
//    
//  }

  CreateEventHisto(physics,triggerClassName,"Nevents","number of events",2,-0.5,1.5);  

  xmin = 0;
  xmax = 3564;
  nbins = GetNbins(xmin,xmax,1.0);
  
  CreateEventHisto(physics,triggerClassName,"BCX","bunch-crossing ids",nbins,xmin-0.5,xmax-0.5);

  CreateEventHisto(physics,triggerClassName,"BCXD","bunch-crossing distances",nbins,xmin-0.5,xmax-0.5);

  xmin = -200;
  xmax = +200;
  nbins = GetNbins(xmin,xmax,1.0);
  
  xmin = 0;
  xmax = 150;
  nbins = GetNbins(xmin,xmax,2.0);
  
  CreateSingleHisto(physics,triggerClassName,"dcaP23Mu","#mu DCA vs P for 2-3 degrees;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax,fShouldSeparatePlusAndMinus);  

  CreateSingleHisto(physics,triggerClassName,"dcaP310Mu","#mu DCA vs P for 3-10 degrees;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax,fShouldSeparatePlusAndMinus);  

  CreateSingleHisto(physics,triggerClassName,"dcaPwPtCut23Mu","#mu DCA vs P for 2-3 degrees with Pt Cut;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax,fShouldSeparatePlusAndMinus);  

  CreateSingleHisto(physics,triggerClassName,"dcaPwPtCut310Mu","#mu DCA vs P for 3-10 degrees with Pt Cut;P (GeV);DCA (cm)",nbinsP,pMin,pMax,nbins,xmin,xmax,fShouldSeparatePlusAndMinus);  
  
  xmin = -30;
  xmax = +30;
  nbins = GetNbins(xmin,xmax,0.1);
  
  CreateEventHisto(physics,triggerClassName,"V02D","V0C+V0A versus V0A-V0C;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);

  CreateEventHisto(physics,triggerClassName,"V02DwT0BB","V0C+V0A versus V0A-V0C with T0 BB;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);

  CreateEventHisto(physics,triggerClassName,"V02DwT0BG","V0C+V0A versus V0A-V0C with T0 background flag on;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);

  CreateEventHisto(physics,triggerClassName,"V02DwT0PU","V0C+V0A versus V0A-V0C with T0 pile up flag on;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);

  CreateEventHisto(physics,triggerClassName,"V02DwT0SAT","V0C+V0A versus V0A-V0C with T0 satellite flag on;Time V0A - V0C (ns);Time V0A+V0C (ns)",nbins,xmin,xmax,nbins,xmin,xmax);

  /*
  CreateEventHisto(physics,triggerClassName,"T02D","(T0C+T0A)/2 versus (T0A-T0C)/2;Time (T0A-T0C)/2 (ns);Time (T0A+T0C)/2 (ns)",nbins,xmin,xmax,nbins,xmin,xmax);
   CreateEventHisto(physics,triggerClassName,"T0Flags","T0 flags",3,0,3);
   */

  
  
  TH1* h = new TH1F("Centrality","Centrality",12,-10,110);

  fHistogramCollection->Adopt(physics,triggerClassName,h);                                                              
  
  xmin = 0;
  xmax = 5000;
  nbins = GetNbins(xmin,xmax,10);
  
  CreateEventHisto(physics,triggerClassName,"Tracklets","Number of tracklets",nbins,xmin,xmax);
  
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillHistosForTrack(const char* physics,
                                             const char* triggerClassName, 
                                             const char* centrality,
                                             const AliVParticle& track,
                                             Int_t trackIndex)
{
  /// Fill histograms for one track
  
  TLorentzVector p(track.Px(),track.Py(),track.Pz(),
                   TMath::Sqrt(MuonMass2()+track.P()*track.P()));
  
  
  TString charge("");
  
  if ( ShouldSeparatePlusAndMinus() ) 
  { 
    if ( track.Charge() < 0 ) 
    {
      charge = "Minus";
    }
    else
    {
      charge = "Plus";
    }
  }
  
  UInt_t mask = GetTrackMask(trackIndex);
  
  Double_t dca = EAGetTrackDCA(track);
  
  Double_t theta = GetTrackTheta(track); 
  
  TIter next(fSingleTrackCutNames);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    Bool_t test = ( ( str->GetUniqueID() & mask ) == str->GetUniqueID() );
    
    if ( test ) 
    {
      Histo(physics,triggerClassName,centrality,str->String().Data(),"Chi2MatchTrigger")->Fill(EAGetTrackChi2MatchTrigger(track));

      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("EtaRapidityMu%s",charge.Data()))->Fill(p.Rapidity(),p.Eta());

      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PtEtaMu%s",charge.Data()))->Fill(p.Eta(),p.Pt());
      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PtRapidityMu%s",charge.Data()))->Fill(p.Rapidity(),p.Pt());
      
      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PEtaMu%s",charge.Data()))->Fill(p.Eta(),p.P());
      
      
      //      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("XYdcaMu%s",charge.Data()))->Fill(ydca,xdca);
      
      Histo(physics,triggerClassName,centrality,str->String().Data(),Form("Chi2Mu%s",charge.Data()))->Fill(EAGetTrackNormalizedChi2(track));
      
      if ( theta >= Deg2() && theta < Deg3() )         
      {
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaP23Mu%s",charge.Data()))->Fill(p.P(),dca);
        if ( p.Pt() > 2 )
        {
          Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaPwPtCut23Mu%s",charge.Data()))->Fill(p.P(),dca);
        }
      }
      else if ( theta >= Deg3() && theta < Deg10() )
      {
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaP310Mu%s",charge.Data()))->Fill(p.P(),dca);
        if ( p.Pt() > 2 )
        {
          Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaPwPtCut310Mu%s",charge.Data()))->Fill(p.P(),dca);
        }
      }
    }
  }
  
  
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillHistos(const char* physics, const char* triggerClassName, 
                                      const char* centrality, const AliVEvent& event)
{
  /// Fill histograms for /physics/triggerClassName/centrality
  
  Histo(physics,triggerClassName,centrality,"BCX")->Fill(1.0*event.GetBunchCrossNumber());
  
  Histo(physics,triggerClassName,centrality,"Nevents")->Fill(1.0);
  
  const AliVVertex* vertex = event.GetPrimaryVertex();
  
  if ( vertex ) 
  {
    Histo(physics,triggerClassName,centrality,"Xvertex")->Fill(vertex->GetX());
    Histo(physics,triggerClassName,centrality,"Yvertex")->Fill(vertex->GetY());
    Histo(physics,triggerClassName,centrality,"Zvertex")->Fill(vertex->GetZ());
  }

  if ( !fIsFromESD )
  {
    const AliAODTZERO* tzero = static_cast<const AliAODEvent&>(event).GetTZEROData();
  
    if (tzero)
    {
      Histo(physics,triggerClassName,centrality,"T0Zvertex")->Fill(tzero->GetT0VertexRaw());
    }
  }
  else
  {
    const AliESDTZERO* tzero = static_cast<const AliESDEvent&>(event).GetESDTZERO();
    
    if (tzero)
    {
      Histo(physics,triggerClassName,centrality,"T0Zvertex")->Fill(tzero->GetT0zVertex());
    }    
  }
  
  AliVVZERO* vzero = event.GetVZEROData();
    
  if (vzero)
  {
    Float_t v0a = vzero->GetV0ATime();
    Float_t v0c = vzero->GetV0CTime();
    
    Float_t x = v0a-v0c;
    Float_t y = v0a+v0c;
    
    Histo(physics,triggerClassName,centrality,"V02D")->Fill(x,y);
    
    Bool_t background,pileup,satellite;
    
    Bool_t tzero = EAGetTZEROFlags(event,background,pileup,satellite);
    
    if (tzero)
    { 
      if ( background ) 
      {
        Histo(physics,triggerClassName,centrality,"V02DwT0BG")->Fill(x,y);      
      }
      
      if ( pileup ) 
      {
        Histo(physics,triggerClassName,centrality,"V02DwT0PU")->Fill(x,y);
      }
      
      if ( satellite ) 
      {
        Histo(physics,triggerClassName,centrality,"V02DwT0SAT")->Fill(x,y);
      }
      
      if ( !background && !pileup && !satellite )
      {
        Histo(physics,triggerClassName,centrality,"V02DwT0BB")->Fill(x,y);        
      }
    }
  }
  
  /* FIXME : how to properly get multiplicity from AOD and ESD consistently ?
     is is doable at all ?
  Int_t ntracklets(0);  
  AliAODTracklets* tracklets = aod.GetTracklets();   
  if ( tracklets ) 
  {
    ntracklets = tracklets->GetNumberOfTracklets();
  }  
  Histo(physics,triggerClassName,centrality,"Tracklets")->Fill(ntracklets);
  */
  
  // Track loop
  
  Int_t nMuonTracks = EAGetNumberOfMuonTracks(event);
  
  for (Int_t i = 0; i < nMuonTracks; ++i) 
  {
    AliVParticle* tracki = EAGetTrack(event,i);
    
    if (!tracki) continue;
    
    FillHistosForTrack(physics,triggerClassName,centrality,*tracki,i);
    
    TLorentzVector pi(tracki->Px(),tracki->Py(),tracki->Pz(),
                      TMath::Sqrt(MuonMass2()+tracki->P()*tracki->P()));
    
    for (Int_t j = i+1; j < nMuonTracks; ++j) 
    {
      AliVParticle* trackj = EAGetTrack(event,j);
      
      if (!trackj) continue;
      
      TLorentzVector pj(trackj->Px(),trackj->Py(),trackj->Pz(),
                        TMath::Sqrt(MuonMass2()+trackj->P()*trackj->P()));
      
      pj += pi;
      
      TIter next(fPairTrackCutNames);
      AliAnalysisTaskMuMu::PairCut* str;
      
      UInt_t maski(0),maskj(0),maskij(0);
      
      GetPairMask(*tracki,*trackj,i,j,maski,maskj,maskij);
      
      while ( ( str = static_cast<AliAnalysisTaskMuMu::PairCut*>(next()) ) )
      {
        UInt_t singleTrackMask = str->MaskForOneOrBothTrack();
        UInt_t pairMask = str->MaskForTrackPair();
        
        Bool_t testi = ( ( maski & singleTrackMask ) == singleTrackMask ) ;
        Bool_t testj = ( ( maskj & singleTrackMask ) == singleTrackMask ) ;
        Bool_t testij(kTRUE);
        
        if (pairMask>0) testij = ( ( maskij & pairMask ) == pairMask ) ;
        
        if ( ( testi || testj ) && testij )
        {
          Histo(physics,triggerClassName,centrality,str->GetName(),"Chi12")->Fill(EAGetTrackChi2MatchTrigger(*tracki),
                                                                                  EAGetTrackChi2MatchTrigger(*trackj));
          
          Histo(physics,triggerClassName,centrality,str->GetName(),"Rabs12")->Fill(EAGetTrackRabs(*tracki),
                                                                                   EAGetTrackRabs(*trackj));
          
          if ( tracki->Charge() != trackj->Charge() )
          {
            Histo(physics,triggerClassName,centrality,str->GetName(),"MinvUSPt")->Fill(pj.Pt(),pj.M());            
            Histo(physics,triggerClassName,centrality,str->GetName(),"MinvUSRapidity")->Fill(pj.Rapidity(),pj.M());   
            Histo(physics,triggerClassName,centrality,str->GetName(),"USPtRapidity")->Fill(pj.Rapidity(),pj.Pt());            
          }
          else
          {
            Histo(physics,triggerClassName,centrality,str->GetName(),"MinvLSPt")->Fill(pj.Pt(),pj.M());                        
          }
        }
      }
    }
  } //track loop  
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FinishTaskOutput()
{
  /// prune empty histograms BEFORE mergin, in order to save some bytes...
  
  if ( fHistogramCollection )
  {
    fHistogramCollection->PruneEmptyHistograms();      
  }
}

//_____________________________________________________________________________
UInt_t AliAnalysisTaskMuMu::GetEventMask(const AliVEvent& event) const
{
  /// Compute the event mask
  
  /*
  
  kEventAll = all events
  kEventPS = physics selected events
  kEventTVX = events with 0TVX input present
  kEventV0AND = events with 0VBA and 0VBC present
  kEventV0UP = a check of events within a narrow square range of v0a+c vs v0a-c
  kEventZSPD = events with a vertex computed by SPD
  kEventZ7 = events with | zvertex | < 7cm
  kEventZ10 = events with | zvertex | < 10 cm
  kEventSD2 = events with 0SD2 input present (was for PbPb 2010)
  kEventMSL = events with 0MSL input present
   
   */
  
  UInt_t m(AliAnalysisTaskMuMu::kEventAll);
  
  if ( PassPhysicsSelection() ) m |= AliAnalysisTaskMuMu::kEventPS;
  
  UInt_t trigger = EAGetL0TriggerInputs(event); 
  
  UInt_t l0TVXBIT = GetTriggerInputBitMaskFromInputName("0TVX");
  
  if ( ( trigger & (l0TVXBIT) ) == l0TVXBIT )
  {
    m |= AliAnalysisTaskMuMu::kEventTVX;
  }
  
  UInt_t l0VBABIT = GetTriggerInputBitMaskFromInputName("0VBA");
  UInt_t l0VBCBIT = GetTriggerInputBitMaskFromInputName("0VBC");
  
  if ( ( ( trigger & (l0VBABIT ) ) == l0VBABIT ) && 
      ( ( trigger & (l0VBCBIT ) ) == l0VBCBIT ) )
  {
    m |= AliAnalysisTaskMuMu::kEventV0AND;
  }
  
  UInt_t l0MSLBIT = GetTriggerInputBitMaskFromInputName("0MSL");
  
  if ( ( trigger & (l0MSLBIT) ) == l0MSLBIT )
  {
    m |= AliAnalysisTaskMuMu::kEventMSL;
  }
  
  if ( fBeamYear == "PbPb2010" ) 
  {
    // consider only events with OSM2 fired
    UInt_t sd2 = GetTriggerInputBitMaskFromInputName("0SM2");
    if ( ( trigger & sd2 ) == sd2 )
    {
      m |= AliAnalysisTaskMuMu::kEventSD2;
    }
  }
  
  //  Bool_t hasPileUp = aod->IsPileupFromSPD(3,0.8);
  //  Bool_t hasPileUp2 = aod->IsPileupFromSPD(5,0.8);  
  //  Bool_t isIsolated = ( aod->GetBunchCrossNumber() > 1000 && aod->GetBunchCrossNumber() < 2900 );
  
  //  Bool_t hasPileUp2 = aod->IsPileupFromSPDInMultBins();
  
  //  TIter nextV(aod->GetVertices());
  //  AliAODVertex* v;
  //  while ( ( v = static_cast<AliAODVertex*>(nextV())) && !hasPileUp2 )
  //  {
  //    if ( v->GetType() == AliAODVertex::kPileupSPD ) hasPileUp2 = kTRUE;
  //  }
  
  //  Bool_t isIsolated = false;//( aod->GetClosestBunchCrossingDistance() > 10 );
  
  const AliVVertex* vertex = event.GetPrimaryVertex();
  
  if ( vertex->IsA() == AliAODVertex::Class() )
  {
    AliAODVertex* spdVertex = static_cast<const AliAODEvent&>(event).GetPrimaryVertexSPD();

    if ( spdVertex && spdVertex->GetNContributors() > 0 ) 
    {
      m |= AliAnalysisTaskMuMu::kEventZSPD;
    }
  }
  
  if ( TMath::Abs(vertex->GetZ()) < 10.0 )
  {
    m |= AliAnalysisTaskMuMu::kEventZ10;
  }
  
  if ( TMath::Abs(vertex->GetZ()) < 7.0 )
  {
    m |= AliAnalysisTaskMuMu::kEventZ7;
  }
  
  AliVVZERO* vzero = event.GetVZEROData();
  
  if (vzero)
  {
    Float_t v0a = vzero->GetV0ATime();
    Float_t v0c = vzero->GetV0CTime();
    
    Float_t x = v0a-v0c;
    Float_t y = v0a+v0c;
    
    if ( ( x > 6 && x < 10 ) && y > 20 )
    {
      m |= AliAnalysisTaskMuMu::kEventV0UP;
    }
  }
  
  return m;
}

//_____________________________________________________________________________
UInt_t AliAnalysisTaskMuMu::GetTriggerInputBitMaskFromInputName(const char* inputName) const
{
  // Get trigger input bit from its name
  // FIXME : this should really come directly from the trigger configuration
  // object, if only this one would be available in a more convenient
  // way than the OCDB (e.g. in RunBasedContainer ?)
  //

  if ( fTriggerInputBitMap.empty() )
  {
    // nothing given to us, use the bad bad hard-coded values !
    
    TString sInputName(inputName);

    if ( sInputName == "0SM2" ) return (1<<12);

    
    if ( sInputName == "0VBA" ) return (1<<0);
    if ( sInputName == "0VBC" ) return (1<<1);
    if ( sInputName == "0SMB" ) return (1<<2);
    if ( sInputName == "0TVX" ) return (1<<3);
    if ( sInputName == "0VGC" ) return (1<<4);
    if ( sInputName == "0VGA" ) return (1<<5);
    if ( sInputName == "0SH1" ) return (1<<6);
    if ( sInputName == "0SH2" ) return (1<<7);
    if ( sInputName == "0HPT" ) return (1<<8);
    if ( sInputName == "0AMU" ) return (1<<9);
    if ( sInputName == "0OB0" ) return (1<<10);
    if ( sInputName == "0ASL" ) return (1<<11);
    if ( sInputName == "0MSL" ) return (1<<12);
    if ( sInputName == "0MSH" ) return (1<<13);
    if ( sInputName == "0MUL" ) return (1<<14);
    if ( sInputName == "0MLL" ) return (1<<15);
    if ( sInputName == "0EMC" ) return (1<<16);
    if ( sInputName == "0PH0" ) return (1<<17);
    if ( sInputName == "0HWU" ) return (1<<18);
    if ( sInputName == "0LSR" ) return (1<<19);
    if ( sInputName == "0T0A" ) return (1<<20);
    if ( sInputName == "0BPA" ) return (1<<21);
    if ( sInputName == "0BPC" ) return (1<<22);
    if ( sInputName == "0T0C" ) return (1<<23);

    if ( sInputName == "1EJE" ) return (1<<0);
    if ( sInputName == "1EGA" ) return (1<<1);
    if ( sInputName == "1EJ2" ) return (1<<2);
    if ( sInputName == "1EG2" ) return (1<<3);
    if ( sInputName == "1PHL" ) return (1<<4);
    if ( sInputName == "1PHM" ) return (1<<5);
    if ( sInputName == "1PHH" ) return (1<<6);
    if ( sInputName == "1HCO" ) return (1<<8);
    if ( sInputName == "1HJT" ) return (1<<9);
    if ( sInputName == "1HSE" ) return (1<<10);
    if ( sInputName == "1DUM" ) return (1<<11);
    if ( sInputName == "1HQU" ) return (1<<12);
    if ( sInputName == "1H14" ) return (1<<13);
    if ( sInputName == "1ZMD" ) return (1<<14);
    if ( sInputName == "1ZMB" ) return (1<<16);
    if ( sInputName == "1ZED" ) return (1<<17);
    if ( sInputName == "1ZAC" ) return (1<<18);
    if ( sInputName == "1EJE" ) return (1<<19);
    
    AliError(Form("Don't know this input %s",inputName));
    
    return (1<<31);
  }
  else
  {
    std::map<std::string,int>::const_iterator it = fTriggerInputBitMap.find(inputName);
    if ( it != fTriggerInputBitMap.end() )
    {
      return ( 1 << it->second );
    }
    else
    {
      AliError(Form("Don't know this input %s",inputName));
      
      return (1<<31);
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::GetPairMask(const AliVParticle& t1, const AliVParticle& t2,
                                      Int_t trackIndex1, Int_t trackIndex2,
                                      UInt_t& mask1, UInt_t& mask2,
                                      UInt_t& mask12) const
{
  /// Get the mask of the track pair
  
  mask1 = GetTrackMask(trackIndex1);
  mask2 = GetTrackMask(trackIndex2);
  
  mask12 = mask1 & mask2;
  
  if ( PairRapidityCut(t1,t2) ) mask12 |= kPairRapidity;
}

//_____________________________________________________________________________
UInt_t AliAnalysisTaskMuMu::GetTrackMask(Int_t trackIndex) const
{
  /// Get the mask of all the cuts this track pass
  
  return static_cast<UInt_t>(fPrecomputedTrackMasks.At(trackIndex));
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::ComputeTrackMask(const AliVParticle& track, Int_t trackIndex)
{
  // Compute the track mask
  UInt_t m(kAll);
  
  UInt_t selectionMask = fMuonTrackCuts ? fMuonTrackCuts->GetSelectionMask(&track) : 0;
  
  if ( ( selectionMask & AliMuonTrackCuts::kMuThetaAbs ) == AliMuonTrackCuts::kMuThetaAbs ) m |= kRabs;
  
  Double_t angle = GetTrackTheta(track);
  
  if ( angle >= Deg2() && angle < Deg3() ) m |= kDeg23;
  
  if ( angle >= Deg3() && angle < Deg10() ) m |= kDeg310;
  
  if ( selectionMask & AliMuonTrackCuts::kMuEta ) m |= kEta;

  Double_t pt = track.Pt();
  
  if ( pt >  1.0 ) m |= kPt1;  
  if ( pt >  1.2 ) m |= kPt1dot2;  
  if ( pt >  1.5 ) m |= kPt1dot5;  
  if ( pt >  2.0 ) m |= kPt2;  

  if ( track.P() > 10.0 ) m |= kP10;  
  
  if ( pt < 4.0 ) m |= kBelowPt;
  
  if ( ( selectionMask & AliMuonTrackCuts::kMuMatchApt ) == AliMuonTrackCuts::kMuMatchApt ) m |= kMatched;
  
  if ( ( selectionMask & AliMuonTrackCuts::kMuMatchLpt ) == AliMuonTrackCuts::kMuMatchLpt ) m |= kMatchedLow;
  
  if ( ( selectionMask & AliMuonTrackCuts::kMuMatchHpt ) == AliMuonTrackCuts::kMuMatchHpt) m |= kMatchedHigh;
  
  if ( ( selectionMask & AliMuonTrackCuts::kMuTrackChiSquare ) ==  AliMuonTrackCuts::kMuTrackChiSquare ) m |= kChi2;
  
  if ( ( selectionMask & AliMuonTrackCuts::kMuPdca ) == AliMuonTrackCuts::kMuPdca ) m |= kDCA;
  
  if ( EAGetTrackChi2MatchTrigger(track) < 16.0 ) m |= kChi2MatchTrigger;
  
  fPrecomputedTrackMasks.SetAt(m,trackIndex);
}


//_____________________________________________________________________________
TH1* AliAnalysisTaskMuMu::Histo(const char* physics, const char* triggerClassName, const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(physics,triggerClassName,histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisTaskMuMu::Histo(const char* physics, const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(physics,histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisTaskMuMu::Histo(const char* physics,
                            const char* triggerClassName, 
                            const char* what,
                            const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(physics,triggerClassName,what,histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisTaskMuMu::Histo(const char* physics,
                            const char* triggerClassName, 
                            const char* cent,
                            const char* what,
                            const char* histoname)
{
  /// Get one histo back

  return fHistogramCollection ? fHistogramCollection->Histo(physics,triggerClassName,cent,what,histoname) : 0x0;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::IsPP() const
{
  // whether we're dealing with proton proton collisions
  return fBeamYear.Contains("pp");
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::MergeCentralities(AliHistogramCollection* histogramCollection)
{
  /// Merge CENT10 + CENT20 + ... into CENTMB
  
  TList* listA = histogramCollection->CreateListOfKeysA();
  TList* listB = histogramCollection->CreateListOfKeysB();
  TList* listC = histogramCollection->CreateListOfKeysC();
  TList* listD = histogramCollection->CreateListOfKeysD();
  
  if (!listA) 
  {
    AliErrorClass("listA=0x0");
    return;
  }
  
  if (!listB) 
  {
    AliErrorClass("listB=0x0");
    return;
  }
  
  if (!listC) 
  {
    AliErrorClass("listC=0x0");
    return;
  }
  
  if (!listD) 
  {
    AliErrorClass("listD=0x0");
    return;
  }
  
  for ( Int_t id = 0; id <= listD->GetLast(); ++id ) 
  {
    TString keyD = static_cast<TObjString*>(listD->At(id))->String();
    
    for ( Int_t ia = 0; ia <= listA->GetLast(); ++ia )
    {
      TString keyA = static_cast<TObjString*>(listA->At(ia))->String();
      
      for ( Int_t ib = 0; ib <= listB->GetLast(); ++ib ) 
      {
        TString keyB = static_cast<TObjString*>(listB->At(ib))->String();
        
        TList* list = new TList;
        list->SetOwner(kTRUE);
        
        AliHistogramCollection* hmerge(0x0);
        
        for ( Int_t ic = 0; ic <= listC->GetLast(); ++ic ) 
        {
          TString keyC = static_cast<TObjString*>(listC->At(ic))->String();
          
          if ( keyC != "CENTX" && keyC != "CENTMB" )
          {
            AliHistogramCollection* hc = histogramCollection->Project(keyA.Data(),keyB.Data(),keyC.Data(),keyD.Data());
            if (!hmerge) 
            {
              hmerge = hc;
            }
            else
            {
              list->Add(hc);
            }
          }
        }
        if (hmerge)
        {
          hmerge->Merge(list);
          TIter next(hmerge->CreateIterator());
          TH1* h;
          while ( ( h = static_cast<TH1*>(next()) ) )
          {
            histogramCollection->Adopt(keyA.Data(),keyB.Data(),"CENTMB",keyD.Data(),static_cast<TH1*>(h->Clone()));
          }
        }
        delete list;
        delete hmerge;
      }      
    }
  }
  
  delete listA;
  delete listB;
  delete listC;
  delete listD;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskMuMu::MuonMass2() const
{
  /// A usefull constant
  static Double_t m2 = 1.11636129640000012e-02; // using a constant here as the line below is a problem for CINT...
  //  static Double_t m2 = TDatabasePDG::Instance()->GetParticle("mu-")->Mass()*TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  return m2;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::NotifyRun()
{
  /// Called at each change of run 
  
  AliDebug(1,Form("Run %09d File %s",fCurrentRunNumber,CurrentFileName()));
 
  if (!fMuonTrackCuts)
  {
    fMuonTrackCuts = new AliMuonTrackCuts;
      
      fMuonTrackCuts->SetAllowDefaultParams(kTRUE );

    
    fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuEta |
                                  AliMuonTrackCuts::kMuThetaAbs |
                                  AliMuonTrackCuts::kMuPdca |
                                  AliMuonTrackCuts::kMuMatchApt |
                                  AliMuonTrackCuts::kMuMatchLpt |
                                  AliMuonTrackCuts::kMuMatchHpt |
                                  AliMuonTrackCuts::kMuTrackChiSquare);    
  }
  
  fMuonTrackCuts->SetRun(fInputHandler);
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::PairRapidityCut(const AliVParticle& t1, const AliVParticle& t2) const
{
  /// Whether the pair passes the rapidity cut
  
  TLorentzVector p1(t1.Px(),t1.Py(),t1.Pz(),TMath::Sqrt(MuonMass2()+t1.P()*t1.P()));
  TLorentzVector p2(t2.Px(),t2.Py(),t2.Pz(),TMath::Sqrt(MuonMass2()+t2.P()*t2.P()));
  
  TLorentzVector total(p1+p2);
  
  Double_t y = total.Rapidity();
  
  Bool_t ok = ( y < -2.5 && y > -4.0 );
  
  return ok;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::PassPhysicsSelection() const
{
  /// Whether the event pass the physics selection or not
  
  AliInputEventHandler* inputEventHandler = static_cast<AliInputEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

  Bool_t isPhysicsSelected = (inputEventHandler->IsEventSelected() & AliVEvent::kAny);
  
  /*
  ( inputEventHandler->IsEventSelected() & AliVEvent::kMUSH7 )
  || ( inputEventHandler->IsEventSelected() & AliVEvent::kMUL7 )
  || ( inputEventHandler->IsEventSelected() & AliVEvent::kMUS7 )
  || ( inputEventHandler->IsEventSelected() & AliVEvent::kMUU7 )
  || ( inputEventHandler->IsEventSelected() & AliVEvent::kMuonSingleLowPt8 )
  || ( inputEventHandler->IsEventSelected() & AliVEvent::kMuonSingleHighPt8 )
  || ( inputEventHandler->IsEventSelected() & AliVEvent::kMuonLikeLowPt8 )
  || ( inputEventHandler->IsEventSelected() & AliVEvent::kMuonUnlikeLowPt8 );

  if ( IsPP() )
  {
    isPhysicsSelected = isPhysicsSelected || ( inputEventHandler->IsEventSelected() & AliVEvent::kINT7 ) || ( inputEventHandler->IsEventSelected() & AliVEvent::kINT8 );
  }
  else
  {
    isPhysicsSelected = isPhysicsSelected ||  ( inputEventHandler->IsEventSelected() & AliVEvent::kAny );    
  }
  */
  
  return isPhysicsSelected;
}

//_____________________________________________________________________________
void 
AliAnalysisTaskMuMu::Print(Option_t* /*opt*/) const
{
  /// Print the definition of this analysis
  
  cout << ClassName() << " - " << GetName() << " - " << fBeamYear.Data() << endl;
  
  if  ( !fSingleTrackCutNames || !fSingleTrackCutNames )
  {
    cout << "No single track cut defined yet" << endl;    
  }
  else
  {
    TIter next(fSingleTrackCutNames);
    TObjString* str;
    
    while ( ( str = static_cast<TObjString*>(next()) ) )
    {
      cout << Form("SINGLE CUT %20s MASK %x",str->String().Data(),str->GetUniqueID()) << endl;
    }
  }
  
  if  ( !fPairTrackCutNames || !fPairTrackCutNames )
  {
    cout << "No track pair cut defined yet" << endl;    
  }
  else
  {
    TIter next2(fPairTrackCutNames);
    AliAnalysisTaskMuMu::PairCut* str;

    while ( ( str = static_cast<AliAnalysisTaskMuMu::PairCut*>(next2()) ) )
    {
      str->Print();
    }
  }
  
  if ( !fTriggerClasses || !fTriggerClasses->First() ) 
  {
    cout << "No trigger classes defined yet" << endl;
  }
  else
  {
    cout << "Trigger classes that will be considered:" << endl;
    TIter next(fTriggerClasses);
    TObjString* s;
    
    while ( ( s = static_cast<TObjString*>(next()) ) )
    {
      cout << s->String().Data() << endl;
    }
  }
}

//_____________________________________________________________________________
void
AliAnalysisTaskMuMu::Terminate(Option_t *)
{
  /// Called once at the end of the query
  /// Just a simple printout of the stat we analyse and how many histograms
  /// we got
  
  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  
  if (!fOutput) return;
  
  fHistogramCollection = dynamic_cast<AliHistogramCollection*>(fOutput->FindObject("mumu"));
  
  if (!fHistogramCollection)
  {
    AliError("Could not find back histogram collection in output...");
    return;
  }
  
  fHistogramCollection->PruneEmptyHistograms();

  UInt_t size2 = fHistogramCollection->EstimateSize();

  AliInfo(Form("size after prune histograms = %5.1f MB",size2/1024.0/1024.0));
  
  if ( !IsPP() && fCentralityLimits.size() > 1 ) 
  {
    MergeCentralities(fHistogramCollection);    
  }
    
  BeautifyHistos();

  fHistogramCollection->Print("*Minv*");
  
  fEventCounters = dynamic_cast<AliCounterCollection*>(fOutput->FindObject("eventCounters"));
  
  if (!fEventCounters)
  {
    AliError("Could not find back counters in output...");
    return;
  }
  
  fEventCounters->Print("trigger");
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::TriggerSBACECondition(const TString& triggerName) const
{
  // check the beam condition in the trigger name
  
  if ( triggerName.Contains("-S-") ) return kTRUE;

  if ( triggerName.Contains("-B-") ) return kTRUE;

  if ( fUseBackgroundTriggers )
  {
    if ( triggerName.Contains("-ACE-") ) return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::UserExec(Option_t* /*opt*/)
{
  /// Executed at each event
  
  AliVEvent* event = InputEvent();
  
  if (!event) return;
  
  Bool_t anyTrigger(kFALSE);
  
  if ( IsDynamicTriggerClasses() ) 
  {
    AddTriggerClasses(EAGetFiredTriggerClasses(*event)," ");
  }
  else
  {
    if (fTriggerClasses->Contains("ANY") )
    {
      anyTrigger = kTRUE;
    }
  }
  
  TString centrality(DefaultCentralityName());
  
  Float_t fcent(-1);
  
  AliCentrality* acent = event->GetCentrality();
  
  if (acent) fcent = acent->GetCentralityPercentile("V0M");
  
  Double_t cent(-100);
  
  if ( fcent > 0 )
  {
    for ( std::vector<double>::size_type i = 0 ; i < fCentralityLimits.size() && cent < 0 ; ++i )
    {
      if ( fcent < fCentralityLimits[i] ) 
      {
        cent = fCentralityLimits[i];
      }
    }
  }
  
  if ( cent > -1 ) 
  {
    centrality = CentralityName(cent);
  }
  
  int nmu = EAGetNumberOfMuonTracks(*event);
  
  EAComputeTrackMasks(*event);
  
  TIter next(fTriggerClasses);
  TObjString* tname;
  TString firedTriggerClasses(EAGetFiredTriggerClasses(*event));

  // first loop to count things not associated to a specific trigger
  TIter nextEventCut(fEventCutNames);
  TObjString* et;
  
  UInt_t mask = GetEventMask(*event);

  while ( ( et = static_cast<TObjString*>(nextEventCut()) ) )
  {
    Bool_t test = ( ( et->GetUniqueID() & mask ) == et->GetUniqueID() );
    
    if ( test )
    {
      fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "EVERYTHING", fCurrentRunNumber));

      if ( firedTriggerClasses == "" )
      {
        fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "EMPTY", fCurrentRunNumber));
      }
      
      if (nmu)
      {
        fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "ATLEASTONEMUONTRACK", fCurrentRunNumber));
      }
      
      if ( AtLeastOneMuonTrigger(firedTriggerClasses) )
      {
        fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "ATLEASTONEMUONTRIGGER", fCurrentRunNumber));
        if ( AtLeastOneMBTrigger(firedTriggerClasses) )
        {
          fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "ATLEASTONEMUONORMBTRIGGER", fCurrentRunNumber));
        }
      }

      if ( AtLeastOneMBTrigger(firedTriggerClasses) )
      {
        fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "ATLEASTONEMBTRIGGER", fCurrentRunNumber));
      }

//      if ( AtLeastOneEmcalTrigger(firedTriggerClasses) )
//      {
//        fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "ATLEASTONEEMCALTRIGGER", fCurrentRunNumber));
//
//        if ( AtLeastOneMBTrigger(firedTriggerClasses) )
//        {
//          fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "ATLEASTONEEMCALORMBTRIGGER", fCurrentRunNumber));
//        }
//      }

    }
  }

  UInt_t l0Inputs = EAGetL0TriggerInputs(*event);

  // second loop to count only the triggers we're interested in
  while ( ( tname = static_cast<TObjString*>(next()) ) )
  {
    if ( !CheckTriggerClass(tname->String(),firedTriggerClasses,l0Inputs) ) continue;
    
    nextEventCut.Reset();
    
    while ( ( et = static_cast<TObjString*>(nextEventCut()) ) )
    {
      Bool_t test = ( ( et->GetUniqueID() & mask ) == et->GetUniqueID() );
      
      if ( test ) 
      {
        Fill(et->String().Data(),tname,centrality,fcent,*event);
      }
    }
  }
  
  // Post output data.
  PostData(1, fOutput);      
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  OpenFile(1);
  
  fTriggerClasses->Print();  
  
  fOutput = new TList;
  fOutput->SetOwner(kTRUE);
  
  fHistogramCollection = new AliHistogramCollection("mumu");
  
  fOutput->Add(fHistogramCollection);  
  
  // initialize event counters
  fEventCounters = new AliCounterCollection("eventCounters");

  TIter nextEventCutName(fEventCutNames);
  TObjString* str;
  TString eventRubric;
  while ( ( str = static_cast<TObjString*>(nextEventCutName()) ) )
  {
    if ( eventRubric.Length() > 0 ) eventRubric += "/"; 
    eventRubric += str->String();
  }
  
  eventRubric += "/SCALER";
  
  fEventCounters->AddRubric("event", eventRubric.Data());
  
  fEventCounters->AddRubric("trigger", 100);
  
  fEventCounters->AddRubric("run", 1000000);
  
  fEventCounters->Init();
  
  fOutput->Add(fEventCounters);  
  
  // Post output data.
  PostData(1, fOutput);        
}
