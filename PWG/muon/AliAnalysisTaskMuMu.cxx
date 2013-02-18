#include "AliAnalysisTaskMuMu.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisMuMuBinning.h"
#include "AliAnalysisMuonUtility.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODTZERO.h"
#include "AliCentrality.h"
#include "AliCodeTimer.h"
#include "AliCounterCollection.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliESDTZERO.h"
#include "AliInputEventHandler.h"
#include "AliLog.h" 
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMergeableCollection.h"
#include "AliMuonEventCuts.h"
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
#include "TRegexp.h"
#include "TROOT.h"
#include <algorithm>
#include <cassert>

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
    a->Add(new TObjString("CINT7-B-"));
    a->Add(new TObjString("CINT8-S-"));
    
    return a;
  }
  
  TString GetMinvHistoName(const AliAnalysisMuMuBinning::Range& r)
  {
    return TString::Format("MinvUS%s",r.AsString().Data());
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
fHistogramCollection(0),
fEventCounters(0),
fMuonTrackCuts(0x0),
fPrecomputedTrackMasks(),
fIsFromESD(kFALSE),
fShouldSeparatePlusAndMinus(kFALSE),
fBeamYear("pp"),
fSingleTrackCutNames(0x0),
fPairTrackCutNames(0x0),
fCentralityNames(0x0),
fEventCutNames(0x0),
fUseBackgroundTriggers(kFALSE),
fTriggerInputBitMap(),
fBinning(0x0),
fHistogramToDisable(0x0),
fBinArray(0x0),
fHasMC(kFALSE),
fEventCuts(0x0)
{
  /// default ctor
}

//_____________________________________________________________________________
AliAnalysisTaskMuMu::AliAnalysisTaskMuMu(Bool_t fromESD, const char* beamYear, TArrayF* centralities)
: AliAnalysisTaskSE(Form("AliAnalysisTaskMuMu-from%s",fromESD ? "ESD":"AOD")),
fHistogramCollection(0),
fEventCounters(0),
fMuonTrackCuts(0x0),
fPrecomputedTrackMasks(),
fIsFromESD(fromESD),
fShouldSeparatePlusAndMinus(kFALSE),
fBeamYear(beamYear),
fSingleTrackCutNames(0x0),
fPairTrackCutNames(0x0),
fCentralityNames(new TObjArray),
fEventCutNames(0x0),
fUseBackgroundTriggers(kFALSE),
fTriggerInputBitMap(),
fBinning(0x0),
fHistogramToDisable(0x0),
fBinArray(0x0),
fHasMC(kFALSE),
fEventCuts(0x0)
{
  /// Constructor
  /// The list of triggers to be considered will be updated on the fly
  
  DefineOutput(1,AliMergeableCollection::Class());
  DefineOutput(2,AliCounterCollection::Class());
  DefineOutput(3,AliAnalysisMuMuBinning::Class());
  
  DefineDefaultBinning();

  DefineCentralityClasses(centralities);
}

//_____________________________________________________________________________
AliAnalysisTaskMuMu::AliAnalysisTaskMuMu(Bool_t fromESD, TList* triggerClasses, const char* beamYear, TArrayF* centralities)
: AliAnalysisTaskSE(Form("AliAnalysisTaskMuMu-from%s",fromESD ? "ESD":"AOD")),
fHistogramCollection(0),
fEventCounters(0),
fMuonTrackCuts(0x0),
fPrecomputedTrackMasks(),
fIsFromESD(fromESD),
fShouldSeparatePlusAndMinus(kFALSE),
fBeamYear(beamYear),
fSingleTrackCutNames(0x0),
fPairTrackCutNames(0x0),
fCentralityNames(new TObjArray),
fEventCutNames(0x0),
fUseBackgroundTriggers(kFALSE),
fTriggerInputBitMap(),
fBinning(0x0),
fHistogramToDisable(0x0),
fBinArray(0x0),
fHasMC(kFALSE),
fEventCuts(0x0)
{
  /// Constructor with a predefined list of triggers to consider

  DefineOutput(1,AliMergeableCollection::Class());
  DefineOutput(2,AliCounterCollection::Class());
  DefineOutput(3,AliAnalysisMuMuBinning::Class());
  
  TObjString* tname;
  TIter next(triggerClasses);
  TString tclasses;
  
  while ( ( tname = static_cast<TObjString*>(next()) ) )
  {
    if (tclasses.Length()>0)
    {
      tclasses += ",";
    }
    
    tclasses += tname->String();    
  }

  EventCuts()->SetTrigClassPatterns(tclasses);
  
  DefineDefaultBinning();
  
  DefineCentralityClasses(centralities);
}

//_____________________________________________________________________________
AliAnalysisTaskMuMu::~AliAnalysisTaskMuMu()
{
  /// dtor

  if (fHistogramCollection && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode())
  {
    delete fHistogramCollection;
  }

  if (fEventCounters && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode())
  {
    delete fEventCounters;
  }

  if (fBinning && ! AliAnalysisManager::GetAnalysisManager()->IsProofMode())
  {
    delete fBinning;
  }

  delete fMuonTrackCuts;

  delete fSingleTrackCutNames;
  delete fPairTrackCutNames;
  delete fCentralityNames;
  delete fEventCutNames;
  delete fHistogramToDisable;
  delete fBinArray;
  delete fEventCuts;
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::AddBin(const char* particle, const char* type,
                                 Double_t xmin, Double_t xmax,
                                 Double_t ymin, Double_t ymax)
{
  /// Add one bin
  fBinning->AddBin(particle,type,xmin,xmax,ymin,ymax);

  // invalidate cached bins, if any
  if (fBinArray)
  {
    delete fBinArray;
    fBinArray = 0x0;
  }
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
  if (!fHistogramCollection->Histo(Form("/%s/%s/%s/Zvertex",physics,triggerClassName,DefaultCentralityName())))
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
    else if ( name.Contains("Minus") )
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
void AliAnalysisTaskMuMu::ComputeTrackMask(const AliVParticle& track, Int_t trackIndex)
{
  // Compute the track mask
  UInt_t m(kAll);
  
  UInt_t selectionMask = fMuonTrackCuts ? fMuonTrackCuts->GetSelectionMask(&track) : 0;
  
  if ( ( selectionMask & AliMuonTrackCuts::kMuThetaAbs ) == AliMuonTrackCuts::kMuThetaAbs ) m |= kRabs;
  
  Double_t angle = AliAnalysisMuonUtility::GetThetaAbsDeg(&track);
  
  if ( angle >= 2.0 && angle < 3.0 ) m |= kDeg23;
  
  if ( angle >= 3.0 && angle < 10.0 ) m |= kDeg310;
  
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
  
  if ( AliAnalysisMuonUtility::GetChi2MatchTrigger(&track) < 16.0 ) m |= kChi2MatchTrigger;
  
  fPrecomputedTrackMasks.SetAt(m,trackIndex);
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::CreateMesh(const char* particle, const char* type1, const char* type2, Bool_t remove12)
{
  /// Create a 2d binning from 2 1d binning.
  /// WARNING : not fully tested yet
  
  if (!fBinning)
  {
    AliError("Cannot create a mesh as I have no bin at all !");
    return;
  }
  fBinning->CreateMesh(particle,type1,type2,remove12);
}

//_____________________________________________________________________________
void
AliAnalysisTaskMuMu::CreateMinvHistograms(const char* physics,
                                          const char* triggerClassName)
{
  /// Create invariant mass histograms
  
  Double_t minvMin = 0;
  Double_t minvMax = 16;
  Int_t nMinvBins = GetNbins(minvMin,minvMax,0.025);

  Int_t nMCMinvBins = GetNbins(minvMin,minvMax,0.1);

  TObjArray* bins = fBinning->CreateBinObjArray("psi","pt vs y,pt,y,phi");
  
  CreatePairHisto(physics,triggerClassName,"Pt","#mu+#mu- Pt distribution",
                  200,0,20);

//  CreatePairHisto(physics,triggerClassName,"BinFlowPt","#mu+#mu- BinFlowPt distribution",
//                  200,0,20);

  CreatePairHisto(physics,triggerClassName,"PtRecVsSim","#mu+#mu- Pt distribution rec vs sim",
                  200,0,20,200,0,20);

  TIter next(bins);
  AliAnalysisMuMuBinning::Range* r;
  while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(next()) ) )
  {
    TString hname(GetMinvHistoName(*r));    
    
    if ( IsHistogramDisabled(hname.Data()) ) continue;
    
    AliDebug(1,Form("histoname = %s",hname.Data()));
    
    CreatePairHisto(physics,triggerClassName,hname.Data(),
                    Form("#mu+#mu- inv. mass %s",r->AsString().Data()),
                    nMinvBins,minvMin,minvMax);

    TH1* h = fHistogramCollection->Histo("/INPUT/ALL",hname.Data());
    if (!h)
    {
      h = new TH1F(hname.Data(),Form("MC #mu+#mu- inv. mass %s",r->AsString().Data()),
                   nMCMinvBins,minvMin,minvMax);

      fHistogramCollection->Adopt("/INPUT/ALL",h);

      fHistogramCollection->Adopt("/INPUT/INYRANGE",static_cast<TH1*>(h->Clone()));
    }
  }
  
  
  delete bins;
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
  
  if ( IsHistogramDisabled(hname) ) return;
  
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

  if ( IsHistogramDisabled(hname) ) return;

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
  
  if ( IsHistogramDisabled(hname) ) return;

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
    
    fHistogramCollection->Adopt(Form("/%s/%s/%s/",physics,triggerClassName,cent->String().Data()),h);
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
  
  if ( IsHistogramDisabled(hname) ) return;

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
    
      fHistogramCollection->Adopt(Form("/%s/%s/%s/%s",physics,triggerClassName,cent->String().Data(),tcn->String().Data()),h);
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
//      // default values
//      fCentralityLimits.push_back(10.0);
//      fCentralityLimits.push_back(30.0);
//      fCentralityLimits.push_back(50.0);
//      fCentralityLimits.push_back(80.0);
    }
    else
    {
      for ( Int_t i = 0; i < centralities->GetSize(); ++i ) 
      {
//        fCentralityLimits.push_back(centralities->At(i));
      }
    }
  }
  
//  for ( std::vector<double>::size_type i = 0; i < fCentralityLimits.size(); ++i )
//  {
//    Double_t limit = fCentralityLimits[i];
//    fCentralityNames->Add(new TObjString(CentralityName(limit)));
//  }
  
  fCentralityNames->Add(new TObjString(DefaultCentralityName()));
  fCentralityNames->SetOwner(kTRUE);
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::DefineDefaultBinning()
{
  fBinning = new AliAnalysisMuMuBinning("BIN");
  fBinning->AddBin("psi","pt vs y");
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::DisableHistograms(const char* pattern)
{
  /// Disable the histogramming of all the histograms matching the pattern
  
  TString spattern(pattern);
  if (spattern=="*")
  {
    delete fHistogramToDisable;
    fHistogramToDisable = 0x0;
  }
  
  if (!fHistogramToDisable)
  {
    fHistogramToDisable = new TList;
    fHistogramToDisable->SetOwner(kTRUE);
  }
    
  fHistogramToDisable->Add(new TObjString(spattern));
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::EAComputeTrackMasks()
{
  // compute the track masks for the event
  
  fPrecomputedTrackMasks.Reset();
  Int_t n = EAGetNumberOfMuonTracks();
  fPrecomputedTrackMasks.Set(n);
  
  for ( Int_t i = 0; i < n; ++i )
  {
    AliVParticle* track = AliAnalysisMuonUtility::GetTrack(i,Event());
    ComputeTrackMask(*track,i);
  }
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskMuMu::EAGetNumberOfMuonTracks() const
{
  // Get the number of muon tracks *that are not ghosts*
  
  Int_t ntracks = AliAnalysisMuonUtility::GetNTracks(Event());
  
  for ( Int_t i = 0; i < ntracks; ++i )
  {
    AliVParticle* track = AliAnalysisMuonUtility::GetTrack(i,Event());
    if (AliAnalysisMuonUtility::IsMuonGhost(track)) --ntracks;
  }
  
  return ntracks;
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskMuMu::EAGetTrackDCA(const AliVParticle& track) const
{
  // Get track DCA
  
  Double_t xdca = AliAnalysisMuonUtility::GetXatDCA(&track);
  Double_t ydca = AliAnalysisMuonUtility::GetYatDCA(&track);
  
  return TMath::Sqrt(xdca*xdca+ydca*ydca);
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::EAGetTZEROFlags(Bool_t& backgroundFlag, Bool_t& pileupFlag, Bool_t& satelliteFlag) const
{
  // get the TZERO decisions
  // return false if there's no tzero information in this event
  
  Bool_t rv(kFALSE);
  
  if ( Event()->IsA() == AliESDEvent::Class() )
  {
    const AliESDTZERO* tzero = static_cast<AliESDEvent*>(const_cast<AliVEvent*>(Event()))->GetESDTZERO();
    if ( tzero )
    {
      backgroundFlag = tzero->GetBackgroundFlag();
      pileupFlag = tzero->GetPileupFlag();
      satelliteFlag = tzero->GetSatellite();
      rv = kTRUE;
    }    
  }
  else if ( Event()->IsA() == AliAODEvent::Class() )
  {
    AliAODTZERO* tzero = static_cast<const AliAODEvent*>(Event())->GetTZEROData();
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
    AliError(Form("Unknown class for the event = %s",Event()->ClassName()));
  }  
  
  return rv;

}

//_____________________________________________________________________________
AliVEvent*
AliAnalysisTaskMuMu::Event() const
{
  // some const-dirty-dancing
  return const_cast<AliAnalysisTaskMuMu*>(this)->InputEvent();
}

//_____________________________________________________________________________
AliMuonEventCuts*
AliAnalysisTaskMuMu::EventCuts() const
{
  /// Return the single instance of AliMuonEventCuts object we're using
  
  if (!fEventCuts)
  {
    fEventCuts = new AliMuonEventCuts("EventCut","");
  }
  return fEventCuts;
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
void AliAnalysisTaskMuMu::FillMC()
{
  // Fill the input MC histograms
  
  if (!HasMC()) return;

  // Specific things for MC
  if (!Histo("INPUT","ALL","Pt"))
  {
    fHistogramCollection->Adopt("/INPUT/ALL",new TH1F("Pt","Pt",200,0,20));
    fHistogramCollection->Adopt("/INPUT/INYRANGE",new TH1F("Pt","Pt",200,0,20));
    
    Double_t rapidityMin = -5;
    Double_t rapidityMax = -2;
    Int_t nbinsRapidity = GetNbins(rapidityMin,rapidityMax,0.05);
    
    fHistogramCollection->Adopt("/INPUT/ALL",new TH1F("Y","Y",nbinsRapidity,rapidityMin,rapidityMax));
    fHistogramCollection->Adopt("/INPUT/INYRANGE",new TH1F("Y","Y",nbinsRapidity,rapidityMin,rapidityMax));
    
    Double_t etaMin = -5;
    Double_t etaMax = -2;
    Int_t nbinsEta = GetNbins(etaMin,etaMax,0.05);
    
    fHistogramCollection->Adopt("/INPUT/ALL",new TH1F("Eta","Eta",nbinsEta,etaMin,etaMax));
    fHistogramCollection->Adopt("/INPUT/INYRANGE",new TH1F("Eta","Eta",nbinsEta,etaMin,etaMax));
  }

  Int_t nMCTracks = MCEvent()->GetNumberOfTracks();

  if (!fBinArray)
  {
    fBinArray = fBinning->CreateBinObjArray("psi","pt vs y,pt,y");
  }

  TIter nextBin(fBinArray);
  AliAnalysisMuMuBinning::Range* r;
  
  for ( Int_t i = 0; i < nMCTracks; ++i )
  {
    AliVParticle* part = MCEvent()->GetTrack(i);
    
//    std::cout << "part " << i << " isprimary=" << AliAnalysisMuonUtility::IsPrimary(part,MCEvent()) << " motherindex=" << AliAnalysisMuonUtility::GetMotherIndex(part) << std::endl;
//    
//    part->Print();
    
    if  (AliAnalysisMuonUtility::IsPrimary(part,MCEvent()) &&
         AliAnalysisMuonUtility::GetMotherIndex(part)==-1)
    {
      
      Histo("INPUT","ALL","Pt")->Fill(part->Pt());
      Histo("INPUT","ALL","Y")->Fill(part->Y());
      Histo("INPUT","ALL","Eta")->Fill(part->Eta());
      
      if ( part->Y() < -2.5 && part->Y() > -4.0 )
      {
        Histo("INPUT","INYRANGE","Pt")->Fill(part->Pt());
        Histo("INPUT","INYRANGE","Y")->Fill(part->Y());
        Histo("INPUT","INYRANGE","Eta")->Fill(part->Eta());
      }
      
      nextBin.Reset();
      
      while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
      {
        Bool_t ok(kFALSE);
        
        if ( r->Is2D() || r->IsNullObject() )
        {
          ok = r->IsInRange(part->Y(),part->Pt());
        }
        else
        {
          if ( r->Type() == "PT" )
          {
            ok = r->IsInRange(part->Pt());
          }
          else if ( r->Type() == "Y" )
          {
            ok = r->IsInRange(part->Y());
          }
        }
        
        if ( ok )
        {
          TString hname = GetMinvHistoName(*r);
          
          if (!IsHistogramDisabled(hname.Data()))
          {
            TH1* h = Histo("INPUT","ALL",hname.Data());
            
            if (!h)
            {
              AliError(Form("Could not get ALL %s",hname.Data()));
              continue;
            }
            
            h->Fill(part->M());
            
            if ( part->Y() < -2.5 && part->Y() > -4.0 )
            {
              h = Histo("INPUT","INYRANGE",hname.Data());
              if (!h)
              {
                AliError(Form("Could not get INYRANGE %s",hname.Data()));
                continue;
              }
              h->Fill(part->M());
            }
            
          }

        }
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::Fill(const char* eventtype,
                           TObjString* tname,
                           const char* centrality,
                           float fcent)
{
  // Fill one set of histograms
  
  TString seventtype(eventtype);
  seventtype.ToLower();
  
  fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", seventtype.Data(), tname->GetName(), fCurrentRunNumber));

  if ( !IsHistogrammingDisabled() )
  {
    AssertHistogramCollection(eventtype,tname->String().Data());
    FillHistos(eventtype,tname->String().Data(),centrality);
    if (!IsHistogramDisabled("Centrality"))
    {
      fHistogramCollection->Histo(Form("/%s/%s/Centrality",eventtype,tname->String().Data()))->Fill(fcent);
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillEventHistos(const char* physics, const char* triggerClassName,
                                          const char* centrality)
{
  // Fill event-wise histograms
  
  if (!IsHistogramDisabled("BCX"))
  {
    Histo(physics,triggerClassName,centrality,"BCX")->Fill(1.0*Event()->GetBunchCrossNumber());
  }
  
  if (!IsHistogramDisabled("Nevents"))
  {
    Histo(physics,triggerClassName,centrality,"Nevents")->Fill(1.0);
  }
  
  const AliVVertex* vertex = Event()->GetPrimaryVertex();
  
  if ( vertex )
  {
    if (!IsHistogramDisabled("Xvertex"))
    {
      Histo(physics,triggerClassName,centrality,"Xvertex")->Fill(vertex->GetX());
    }
    if (!IsHistogramDisabled("Yvertex"))
    {
      Histo(physics,triggerClassName,centrality,"Yvertex")->Fill(vertex->GetY());
    }
    if (!IsHistogramDisabled("Zvertex"))
    {
      Histo(physics,triggerClassName,centrality,"Zvertex")->Fill(vertex->GetZ());
    }
  }
  
  if ( !fIsFromESD )
  {
    const AliAODTZERO* tzero = static_cast<const AliAODEvent*>(Event())->GetTZEROData();
    
    if (tzero && !IsHistogramDisabled("T0Zvertex"))
    {
      Histo(physics,triggerClassName,centrality,"T0Zvertex")->Fill(tzero->GetT0VertexRaw());
    }
  }
  else
  {
    const AliESDTZERO* tzero = static_cast<const AliESDEvent*>(Event())->GetESDTZERO();
    
    if (tzero && !IsHistogramDisabled("T0Zvertex"))
    {
      Histo(physics,triggerClassName,centrality,"T0Zvertex")->Fill(tzero->GetT0zVertex());
    }
  }
  
  AliVVZERO* vzero = Event()->GetVZEROData();
  
  if (vzero)
  {
    Float_t v0a = vzero->GetV0ATime();
    Float_t v0c = vzero->GetV0CTime();
    
    Float_t x = v0a-v0c;
    Float_t y = v0a+v0c;
    
    if (!IsHistogramDisabled("V02D"))
    {
      Histo(physics,triggerClassName,centrality,"V02D")->Fill(x,y);
    }
    
    Bool_t background,pileup,satellite;
    
    Bool_t tzero = EAGetTZEROFlags(background,pileup,satellite);
    
    if (tzero)
    {
      if ( background )
      {
        if (!IsHistogramDisabled("V02DwT0BG"))
        {
          Histo(physics,triggerClassName,centrality,"V02DwT0BG")->Fill(x,y);
        }
      }
      
      if ( pileup )
      {
        if (!IsHistogramDisabled("V02DwT0PU"))
        {
          Histo(physics,triggerClassName,centrality,"V02DwT0PU")->Fill(x,y);
        }
      }
      
      if ( satellite )
      {
        if (!IsHistogramDisabled("V02DwT0SAT"))
        {
          Histo(physics,triggerClassName,centrality,"V02DwT0SAT")->Fill(x,y);
        }
      }
      
      if ( !background && !pileup && !satellite )
      {
        if (!IsHistogramDisabled("V02DwT0BB"))
        {
          Histo(physics,triggerClassName,centrality,"V02DwT0BB")->Fill(x,y);
        }
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
  
  Double_t phiMin = -TMath::Pi();
  Double_t phiMax = TMath::Pi();
  Int_t nbinsPhi = GetNbins(phiMin,phiMax,0.05);
  
  CreateSingleHisto(physics,triggerClassName,"Chi2MatchTrigger","Chi2 Match Trigger",72,0,72);

  CreateSingleHisto(physics,triggerClassName,"EtaRapidityMu", "Eta distribution vs Rapidity for #mu", nbinsRapidity,rapidityMin,rapidityMax,nbinsEta,etaMin,etaMax, fShouldSeparatePlusAndMinus);

  CreateSingleHisto(physics,triggerClassName,"PtEtaMu", "P_{T} distribution vs Eta for #mu", nbinsEta,etaMin,etaMax, nbinsPt,ptMin,ptMax,fShouldSeparatePlusAndMinus);

  CreateSingleHisto(physics,triggerClassName,"PtRapidityMu", "P_{T} distribution vs Rapidity for #mu", nbinsRapidity,rapidityMin,rapidityMax, nbinsPt,ptMin,ptMax,fShouldSeparatePlusAndMinus);

  CreateSingleHisto(physics,triggerClassName,"PtPhiMu", "P_{T} distribution vs phi for #mu", nbinsPhi,phiMin,phiMax, nbinsPt,ptMin,ptMax,fShouldSeparatePlusAndMinus);
  

  CreateSingleHisto(physics,triggerClassName,"PEtaMu", "P distribution for #mu",nbinsEta,etaMin,etaMax,nbinsP,pMin,pMax,fShouldSeparatePlusAndMinus);

  Double_t chi2min = 0;
  Double_t chi2max = 20;
  Int_t nbinchi2 = GetNbins(chi2min,chi2max,0.05);
  
  CreateSingleHisto(physics, triggerClassName, "Chi2Mu", "chisquare per NDF #mu", nbinchi2, chi2min, chi2max,fShouldSeparatePlusAndMinus);
  
  CreateMinvHistograms(physics,triggerClassName);
  
  CreatePairHisto(physics,triggerClassName,"Chi12","Chi2MatchTrigger of muon 1 vs muon 2",72,0,72,72,0,72);
  CreatePairHisto(physics,triggerClassName,"Rabs12","Rabs of muon 1 vs muon ",100,0,100,100,0,100);

  Double_t xmin = -40;
  Double_t xmax = +40;
  Int_t nbins = GetNbins(xmin,xmax,0.5);
  
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

  
  if ( !IsHistogramDisabled("Centrality") )
  {
    TH1* h = new TH1F("Centrality","Centrality",12,-10,110);

    fHistogramCollection->Adopt(Form("/%s/%s",physics,triggerClassName),h);
  }
  
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
  
  Double_t theta = AliAnalysisMuonUtility::GetThetaAbsDeg(&track);
  
  TIter next(fSingleTrackCutNames);
  TObjString* str;
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    Bool_t test = ( ( str->GetUniqueID() & mask ) == str->GetUniqueID() );
    
    if ( test ) 
    {
      if (!IsHistogramDisabled("Chi2MatchTrigger"))
      {
        Histo(physics,triggerClassName,centrality,str->String().Data(),"Chi2MatchTrigger")->Fill(AliAnalysisMuonUtility::GetChi2MatchTrigger(&track));
      }
      
      if (!IsHistogramDisabled("EtaRapidityMu*"))
      {
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("EtaRapidityMu%s",charge.Data()))->Fill(p.Rapidity(),p.Eta());
      }
      
      if (!IsHistogramDisabled("PtEtaMu*"))
      {
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PtEtaMu%s",charge.Data()))->Fill(p.Eta(),p.Pt());
      }
      
      if (!IsHistogramDisabled("PtRapidityMu*"))
      {
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PtRapidityMu%s",charge.Data()))->Fill(p.Rapidity(),p.Pt());
      }
      
      if (!IsHistogramDisabled("PEtaMu*"))
      {
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PEtaMu%s",charge.Data()))->Fill(p.Eta(),p.P());
      }
      
      if (!IsHistogramDisabled("PtPhiMu*"))
      {
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("PtPhiMu%s",charge.Data()))->Fill(p.Phi(),p.Pt());
      }
      
      if (!IsHistogramDisabled("Chi2Mu*"))
      {
        Histo(physics,triggerClassName,centrality,str->String().Data(),Form("Chi2Mu%s",charge.Data()))->Fill(AliAnalysisMuonUtility::GetChi2perNDFtracker(&track));
      }
      
      if ( theta >= 2.0 && theta < 3.0 )
      {
        
        if (!IsHistogramDisabled("dcaP23Mu*"))
        {
          Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaP23Mu%s",charge.Data()))->Fill(p.P(),dca);
        }

        if ( p.Pt() > 2 )
        {
          if (!IsHistogramDisabled("dcaPwPtCut23Mu*"))
          {
            Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaPwPtCut23Mu%s",charge.Data()))->Fill(p.P(),dca);
          }
        }
      }
      else if ( theta >= 3.0 && theta < 10.0 )
      {
        if (!IsHistogramDisabled("dcaP310Mu*"))
        {
          Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaP310Mu%s",charge.Data()))->Fill(p.P(),dca);
        }
        if ( p.Pt() > 2 )
        {
          if (!IsHistogramDisabled("dcaPwPtCut310Mu*"))
          {
            Histo(physics,triggerClassName,centrality,str->String().Data(),Form("dcaPwPtCut310Mu%s",charge.Data()))->Fill(p.P(),dca);
          }
        }
      }
    }
  }
  
  
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::FillHistos(const char* physics, const char* triggerClassName, 
                                      const char* centrality)
{
  /// Fill histograms for /physics/triggerClassName/centrality
  
  FillEventHistos(physics,triggerClassName,centrality);
  
  // Track loop
  
  if (!fBinArray)
  {
    fBinArray = fBinning->CreateBinObjArray("psi","pt vs y,pt,y,phi");
  }
  
  Int_t nMuonTracks = AliAnalysisMuonUtility::GetNTracks(Event());
  
  for (Int_t i = 0; i < nMuonTracks; ++i) 
  {
    AliVParticle* tracki = AliAnalysisMuonUtility::GetTrack(i,Event());
    
    if (!tracki) continue;
    
    FillHistosForTrack(physics,triggerClassName,centrality,*tracki,i);
    
    TLorentzVector pi(tracki->Px(),tracki->Py(),tracki->Pz(),
                      TMath::Sqrt(MuonMass2()+tracki->P()*tracki->P()));
    
    for (Int_t j = i+1; j < nMuonTracks; ++j) 
    {
      AliVParticle* trackj = AliAnalysisMuonUtility::GetTrack(j,Event());
      
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
          if (!IsHistogramDisabled("Chi12"))
          {
            Histo(physics,triggerClassName,centrality,str->GetName(),"Chi12")
            ->Fill(
                   AliAnalysisMuonUtility::GetChi2perNDFtracker(tracki),
                   AliAnalysisMuonUtility::GetChi2perNDFtracker(trackj));
          }
          
          if (!IsHistogramDisabled("Rabs12"))
          {
            Histo(physics,triggerClassName,centrality,str->GetName(),"Rabs12")
            ->Fill(AliAnalysisMuonUtility::GetRabs(tracki),
                   AliAnalysisMuonUtility::GetRabs(trackj));
          }
          
          if ( ( tracki->Charge() != trackj->Charge() ) )
          {
            TIter nextBin(fBinArray);
            AliAnalysisMuMuBinning::Range* r;
            
            Histo(physics,triggerClassName,centrality,str->GetName(),"Pt")->Fill(pj.Pt());

            if ( HasMC() )
            {
              Int_t labeli = tracki->GetLabel();
              Int_t labelj = trackj->GetLabel();
            
              if ( labeli < 0 || labelj < 0 )
              {
                AliError("Got negative labels!");
              }
              else
              {
                AliVParticle* mcTracki = MCEvent()->GetTrack(labeli);
                AliVParticle* mcTrackj = MCEvent()->GetTrack(labelj);
                
                TLorentzVector mcpi(mcTracki->Px(),mcTracki->Py(),mcTracki->Pz(),
                                  TMath::Sqrt(MuonMass2()+mcTracki->P()*mcTracki->P()));
                TLorentzVector mcpj(mcTrackj->Px(),mcTrackj->Py(),mcTrackj->Pz(),
                                  TMath::Sqrt(MuonMass2()+mcTrackj->P()*mcTrackj->P()));

                mcpj += mcpi;
                
                Histo(physics,triggerClassName,centrality,str->GetName(),"PtRecVsSim")->Fill(mcpj.Pt(),pj.Pt());
                
              }
            }
            
            while ( ( r = static_cast<AliAnalysisMuMuBinning::Range*>(nextBin()) ) )
            {
              AliDebug(2,Form("bin %s pt %e rad %e = %d",r->AsString().Data(),
                              pj.Pt(),pj.Rapidity(),r->IsInRange(pj.Pt(),pj.Rapidity())));
              
              Bool_t ok(kFALSE);
              
              if ( r->Is2D() || r->IsNullObject() )
              {
                ok = r->IsInRange(pj.Rapidity(),pj.Pt());
              }
              else
              {
                if ( r->Type() == "PT" )
                {
                  ok = r->IsInRange(pj.Pt());
                }
                else if ( r->Type() == "Y" )
                {
                  ok = r->IsInRange(pj.Rapidity());
                }
                else if ( r->Type() == "PHI" )
                {
                  ok = r->IsInRange(pj.Phi());
                }
              }
              
              if ( ok )
              {
                TString hname = GetMinvHistoName(*r);
              
                if (!IsHistogramDisabled(hname.Data()))
                {
                  TH1* h = Histo(physics,triggerClassName,centrality,str->GetName(),hname.Data());
                
                  if (!h)
                  {
                    AliError(Form("Could not get %s",hname.Data()));
                  }
                  else
                  {
                    AliDebug(1,Form("filling %s",hname.Data()));
                  }
                  h->Fill(pj.M());
                }
              }
            }
              
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
    fHistogramCollection->PruneEmptyObjects();
  }
}

//_____________________________________________________________________________
UInt_t AliAnalysisTaskMuMu::GetEventMask() const
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
  kEventNOPILEUP = events with the T0 pile-up flag not present
   
   */
  
  UInt_t eMask = EventCuts()->GetSelectionMask(fInputHandler);

  UInt_t m(AliAnalysisTaskMuMu::kEventAll);

  if ( eMask & AliMuonEventCuts::kPhysicsSelected ) m |= AliAnalysisTaskMuMu::kEventPS;
  
  UInt_t trigger = AliAnalysisMuonUtility::GetL0TriggerInputs(Event());
  
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
  
  const AliVVertex* vertex = Event()->GetPrimaryVertex();
  
  if ( vertex->IsA() == AliAODVertex::Class() )
  {
    AliAODVertex* spdVertex = static_cast<const AliAODEvent*>(Event())->GetPrimaryVertexSPD();

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
  
  AliVVZERO* vzero = Event()->GetVZEROData();
  
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
  
  Bool_t backgroundFlag(kFALSE);
  Bool_t pileupFlag(kFALSE);
  Bool_t satelliteFlag(kFALSE);
  
  EAGetTZEROFlags(backgroundFlag,pileupFlag,satelliteFlag);
  
  if ( !pileupFlag )
  {
    m |= AliAnalysisTaskMuMu::kEventNOTZEROPILEUP;
  }
  
  int nmu = EAGetNumberOfMuonTracks();

  if ( nmu >=1 )
  {
    m |= AliAnalysisTaskMuMu::kEventOFFLINEMUL1;
  }

  if ( nmu >=2 )
  {
    m |= AliAnalysisTaskMuMu::kEventOFFLINEMUL2;
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

    /*
    if ( sInputName == "0MSL" ) return (1<<12);

//    if ( sInputName == "0MSH" ) return (1<<13);
    if ( sInputName == "0MSH" ) return (1<<8);
    
    if ( sInputName == "0MUL" ) return (1<<14);
    if ( sInputName == "0MLL" ) return (1<<15);
     */
    
    if ( sInputName == "0MSL") return (1<<12);
    if ( sInputName == "0MSH") return (1<<13);
    if ( sInputName == "0MUL") return (1<<14);
    if ( sInputName == "0MLL") return (1<<15);
    
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
TH1* AliAnalysisTaskMuMu::Histo(const char* physics, const char* triggerClassName, const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s",physics,triggerClassName,histoname)) : 0x0;
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
  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s",physics,triggerClassName,what),histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisTaskMuMu::Histo(const char* physics,
                            const char* triggerClassName, 
                            const char* cent,
                            const char* what,
                            const char* histoname)
{
  /// Get one histo back

  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s/%s",physics,triggerClassName,cent,what),histoname) : 0x0;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::IsHistogramDisabled(const char* hname) const
{
  /// Whether or not a given histogram (identified by its name)
  /// is disabled or not
  if ( !fHistogramToDisable )
  {
    return kFALSE;
  }
  TString shname(hname);
  TIter next(fHistogramToDisable);
  TObjString* str(0x0);
  
  while ( ( str = static_cast<TObjString*>(next()) ) )
  {
    if ( shname.Contains(TRegexp(str->String()) ) )
    {
      return kTRUE;
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::IsHistogrammingDisabled() const
{
  /// Whether or not *all* histograms are disabled
  
  if ( fHistogramToDisable && fHistogramToDisable->GetEntries()==1 )
  {
    TObjString* r = static_cast<TObjString*>(fHistogramToDisable->First());
    if ( r->String() == "*" )
    {
      return kTRUE;
    }
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMuMu::IsPP() const
{
  // whether we're dealing with proton proton collisions
  return fBeamYear.Contains("pp");
}

//_____________________________________________________________________________
//void AliAnalysisTaskMuMu::MergeCentralities(AliHistogramCollection* histogramCollection)
//{
  /// FIXME : Reimplement using AliMergeableCollection::GetSum ?
  
  /// Merge CENT10 + CENT20 + ... into CENTMB
  
//  TList* listA = histogramCollection->CreateListOfKeysA();
//  TList* listB = histogramCollection->CreateListOfKeysB();
//  TList* listC = histogramCollection->CreateListOfKeysC();
//  TList* listD = histogramCollection->CreateListOfKeysD();
//  
//  if (!listA) 
//  {
//    AliErrorClass("listA=0x0");
//    return;
//  }
//  
//  if (!listB) 
//  {
//    AliErrorClass("listB=0x0");
//    return;
//  }
//  
//  if (!listC) 
//  {
//    AliErrorClass("listC=0x0");
//    return;
//  }
//  
//  if (!listD) 
//  {
//    AliErrorClass("listD=0x0");
//    return;
//  }
//  
//  for ( Int_t id = 0; id <= listD->GetLast(); ++id ) 
//  {
//    TString keyD = static_cast<TObjString*>(listD->At(id))->String();
//    
//    for ( Int_t ia = 0; ia <= listA->GetLast(); ++ia )
//    {
//      TString keyA = static_cast<TObjString*>(listA->At(ia))->String();
//      
//      for ( Int_t ib = 0; ib <= listB->GetLast(); ++ib ) 
//      {
//        TString keyB = static_cast<TObjString*>(listB->At(ib))->String();
//        
//        TList* list = new TList;
//        list->SetOwner(kTRUE);
//        
//        AliHistogramCollection* hmerge(0x0);
//        
//        for ( Int_t ic = 0; ic <= listC->GetLast(); ++ic ) 
//        {
//          TString keyC = static_cast<TObjString*>(listC->At(ic))->String();
//          
//          if ( keyC != "CENTX" && keyC != "CENTMB" )
//          {
//            AliHistogramCollection* hc = histogramCollection->Project(keyA.Data(),keyB.Data(),keyC.Data(),keyD.Data());
//            if (!hmerge) 
//            {
//              hmerge = hc;
//            }
//            else
//            {
//              list->Add(hc);
//            }
//          }
//        }
//        if (hmerge)
//        {
//          hmerge->Merge(list);
//          TIter next(hmerge->CreateIterator());
//          TH1* h;
//          while ( ( h = static_cast<TH1*>(next()) ) )
//          {
//            histogramCollection->Adopt(keyA.Data(),keyB.Data(),"CENTMB",keyD.Data(),static_cast<TH1*>(h->Clone()));
//          }
//        }
//        delete list;
//        delete hmerge;
//      }      
//    }
//  }
//  
//  delete listA;
//  delete listB;
//  delete listC;
//  delete listD;
//}

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
      
    fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
    
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
  
  return  ( y < -2.5 && y > -4.0 );
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
  
//  if ( !fTriggerClasses || !fTriggerClasses->First() ) 
//  {
//    cout << "No trigger classes defined yet" << endl;
//  }
//  else
//  {
//    cout << "Trigger classes that will be considered:" << endl;
//    TIter next(fTriggerClasses);
//    TObjString* s;
//    
//    while ( ( s = static_cast<TObjString*>(next()) ) )
//    {
//      cout << s->String().Data() << endl;
//    }
//  }
//  
  if ( fBinning )
  {
    cout << "Binning for Minv plots" << endl;
    fBinning->Print();
  }
}

//_____________________________________________________________________________
void
AliAnalysisTaskMuMu::Terminate(Option_t *)
{
  /// Called once at the end of the query
  /// Just a simple printout of the stat we analyse and how many histograms
  /// we got
  
  fHistogramCollection = dynamic_cast<AliMergeableCollection*>(GetOutputData(1));

  if (!fHistogramCollection)
  {
    AliError("Could not find back histogram collection in output...");
    return;
  }
  
  fHistogramCollection->PruneEmptyObjects();

  UInt_t size2 = fHistogramCollection->EstimateSize();

  AliInfo(Form("size after prune histograms = %5.1f MB",size2/1024.0/1024.0));
  
//  if ( !IsPP() && fCentralityLimits.size() > 1 ) 
//  {
//    MergeCentralities(fHistogramCollection);    
//  }
  
  BeautifyHistos();

  fHistogramCollection->Print("^MinvUS$");
  
  fEventCounters = dynamic_cast<AliCounterCollection*>(GetOutputData(2));
  
  if (!fEventCounters)
  {
    AliError("Could not find back counters in output...");
    return;
  }
  
  fEventCounters->Print("trigger");
  
  // post param container(s)
  PostData(3,fBinning);
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
  
  fHasMC = (MCEvent()!=0x0);

  TString firedTriggerClasses(AliAnalysisMuonUtility::GetFiredTriggerClasses(Event()));

  TString centrality(DefaultCentralityName());
  
  Float_t fcent(EventCuts()->GetCentrality(Event()));
  
  Double_t cent(-100);
  
  if ( fcent > 0 )
  {
//    for ( std::vector<double>::size_type i = 0 ; i < fCentralityLimits.size() && cent < 0 ; ++i )
//    {
//      if ( fcent < fCentralityLimits[i] ) 
//      {
//        cent = fCentralityLimits[i];
//      }
//    }
  }
  
  if ( cent > -1 ) 
  {
    centrality = CentralityName(cent);
  }
  
  int nmu = EAGetNumberOfMuonTracks();
  
  EAComputeTrackMasks();

  // first loop to count things not associated to a specific trigger
  TIter nextEventCut(fEventCutNames);
  TObjString* et;
  
  UInt_t mask = GetEventMask();

//  TString eventType;
//  eventType.Form("EVENTTYPE%d",event->GetEventType());

  while ( ( et = static_cast<TObjString*>(nextEventCut()) ) )
  {
    Bool_t test = ( ( et->GetUniqueID() & mask ) == et->GetUniqueID() );
    
    if ( test )
    {
      fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "EVERYTHING", fCurrentRunNumber));

      if ( HasMC() )
      {
        fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "HASMC", fCurrentRunNumber));
      }
      
//      fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), eventType.Data(), fCurrentRunNumber));

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

      if ( AtLeastOneEmcalTrigger(firedTriggerClasses) )
      {
        fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "ATLEASTONEEMCALTRIGGER", fCurrentRunNumber));

        if ( AtLeastOneMBTrigger(firedTriggerClasses) )
        {
          fEventCounters->Count(Form("event:%s/trigger:%s/run:%d", et->String().Data(), "ATLEASTONEEMCALORMBTRIGGER", fCurrentRunNumber));
        }
      }

    }
  }

  // second loop to count only the triggers we're interested in
  
  TIter next(EventCuts()->GetSelectedTrigClassesInEvent(Event()));
  TObjString* tname;
  
  while ( ( tname = static_cast<TObjString*>(next()) ) )
  {
    nextEventCut.Reset();
    
    while ( ( et = static_cast<TObjString*>(nextEventCut()) ) )
    {
      Bool_t test = ( ( et->GetUniqueID() & mask ) == et->GetUniqueID() );
      
      if ( test ) 
      {
        Fill(et->String().Data(),tname,centrality,fcent);
      }
    }
  }

  FillMC();
  
  // Post output data.
  PostData(1, fHistogramCollection);
  PostData(2, fEventCounters);
  PostData(3, fBinning);
}

//_____________________________________________________________________________
void AliAnalysisTaskMuMu::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  OpenFile(1);
  
  fHistogramCollection = new AliMergeableCollection("MC");
  
  // initialize event counters
  fEventCounters = new AliCounterCollection("CC");

  TIter nextEventCutName(fEventCutNames);
  TObjString* str;
  TString eventRubric;
  while ( ( str = static_cast<TObjString*>(nextEventCutName()) ) )
  {
    if ( eventRubric.Length() > 0 ) eventRubric += "/"; 
    eventRubric += str->String();
  }
  
  fEventCounters->AddRubric("event", eventRubric.Data());
  
  fEventCounters->AddRubric("trigger", 100);
  
  fEventCounters->AddRubric("run", 1000000);
  
  fEventCounters->Init();
  
  // Post output data.
  PostData(1,fHistogramCollection);
  PostData(2,fEventCounters);
  PostData(3,fBinning);
}
