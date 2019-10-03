#include "AliAnalysisMuMuBase.h"

/**
 *
 * \ingroup pwg-muon-mumu
 *
 * \class AliAnalysisMuMuBase
 *
 * Defines the interface of a sub-analysis to be steered by AliAnalysisTaskMuMu
 *
 * Daugther class has to implement one method :
 *
 * - \ref DefineHistogramCollection
 *
 * It may implement one or several of the FillHistosForXXX methods :
 *
 * - \ref FillHistosForEvent to fill histograms for one event
 * - \ref FillHistosForMCEvent to fill histograms for one MC event
 * - \ref FillHistosForTrack to fill histograms for one track
 * - \ref FillHistosForPair to fill histograms for one muon track pair
 *
 * More rarely it may also implement :
 *
 * - \ref SetRun to do some changes when run number changes
 * - \ref Terminate to finalize before ending
 * - \ref SetEvent to e.g. append information to VEvent
 *
 * Daugther class can use the following methods :
 *
 * - \ref HasMC to know if MC information is available in the analyzed data
 * - \ref Event to access the current event
 * - \ref MCEvent to access to current MC event (if available)
 *
 * A few trivial cut methods (\ref AlwaysTrue and \ref AlwaysFalse) are defined as well and
 * can be used to register some control cut combinations (see \ref AliAnalysisMuMuCutCombination)
 *
 */

#include "AliMergeableCollection.h"
#include "AliCounterCollection.h"
#include "TList.h"
#include "TObjString.h"
#include "TMath.h"
#include "TObjArray.h"
#include "AliLog.h"
#include "AliAnalysisMuMuBinning.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TRegexp.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliAnalysisMuMuCutCombination.h"
#include "AliAnalysisMuMuCutRegistry.h"

ClassImp(AliAnalysisMuMuBase)

//_____________________________________________________________________________
AliAnalysisMuMuBase::AliAnalysisMuMuBase()
:
TObject(),
fEventCounters(0x0),
fHistogramCollection(0x0),
fBinning(0x0),
fCutRegistry(0x0),
fEvent(0x0),
fMCEvent(0x0),
fHistogramToDisable(0x0),
fHasMC(kFALSE)
{
 /// default ctor
}

//_____________________________________________________________________________
TString AliAnalysisMuMuBase::BuildPath(const char* eventSelection, const char* triggerClassName,
                                       const char* centrality, const char* cut) const
{
  TString path;

  path.Form("/%s/%s/%s/",eventSelection,triggerClassName,centrality);
  if ( strlen(cut) > 0 )
  {
    path += cut;
    path += "/";
  }

  return path;
}

//_____________________________________________________________________________
TString AliAnalysisMuMuBase::BuildMCPath(const char* eventSelection, const char* triggerClassName,
                                          const char* centrality, const char* cut) const
{
  TString path = BuildPath(eventSelection,triggerClassName,centrality,cut);

  TString mcPath;

  mcPath.Form("/%s%s",MCInputPrefix(),path.Data());

  return mcPath;
}


//_____________________________________________________________________________
void
AliAnalysisMuMuBase::CreateEventHistos(UInt_t dataType,
                                       const char* what,
                                       const char* hname, const char* htitle,
                                       Int_t nbinsx, Double_t xmin, Double_t xmax,
                                       Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  /** Append histograms at the event level. Depending on dataType the created histograms
   * are for real data only, MC data only, or both
   */

  if ( IsHistogramDisabled(hname) ) return;

  TObjArray pathNames;
  pathNames.SetOwner(kTRUE);

  if ( dataType & kHistoForData )
  {
    pathNames.Add(new TObjString(Form("/%s",what)));
  }
  if ( ( dataType & kHistoForMCInput ) && HasMC() )
  {
    pathNames.Add(new TObjString(Form("/%s/%s",MCInputPrefix(),what)));
  }

  CreateHistos(pathNames,hname,htitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
}

//_____________________________________________________________________________
void
AliAnalysisMuMuBase::CreateEventHistos(UInt_t dataType,
                                       const char* eventSelection,
                                       const char* triggerClassName,
                                       const char* centrality,
                                       const char* hname, const char* htitle,
                                       Int_t nbinsx, Double_t xmin, Double_t xmax,
                                       Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  /** Append histograms at the event level. Depending on dataType the created histograms
   * are for real data only, MC data only, or both
   */

  if ( IsHistogramDisabled(hname) ) return;

  TObjArray pathNames;
  pathNames.SetOwner(kTRUE);

  if ( dataType & kHistoForData )
  {
    pathNames.Add(new TObjString(Form("/%s/%s/%s",eventSelection,triggerClassName,centrality)));
  }
  if ( ( dataType & kHistoForMCInput ) && HasMC() )
  {
    pathNames.Add(new TObjString(Form("/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,centrality)));
  }

  CreateHistos(pathNames,hname,htitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
}

//_____________________________________________________________________________
void
AliAnalysisMuMuBase::CreateHistos(const TObjArray& paths,
                                  const char* hname, const char* htitle,
                                  Int_t nbinsx, Double_t xmin, Double_t xmax,
                                  Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  /// Create multiple copies of histogram hname, one per path

  if ( IsHistogramDisabled(hname) ) return;

  TIter next(&paths);
  TObjString* pathName;

  while ( ( pathName = static_cast<TObjString*>(next()) ) )
  {
    TH1* h(0x0);

    if ( nbinsy > 0 )
    {
      h = new TH2F(hname,htitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
    }
    else if ( nbinsy == 0 )
    {
      h = new TProfile(hname,htitle,nbinsx,xmin,xmax,ymin,ymax);
      h->Sumw2();
      static_cast<TProfile*>(h)->Approximate();
    }
    else
    {
      h = new TH1F(hname,htitle,nbinsx,xmin,xmax);

      if ( nbinsy < -1 )
      {
        h->Sumw2();
      }
    }

    HistogramCollection()->Adopt(pathName->String().Data(),h);
  }
}

//_____________________________________________________________________________
void AliAnalysisMuMuBase::CreateSemaphoreHistogram(const char* eventSelection,
                                                   const char* triggerClassName,
                                                   const char* centrality)
{
  /// Create the semaphore histogram, a dummy histogram used
  /// to check whether or not the DefineHistogramCollection has already been done `
  /// for a tuple (eventSelection,triggerClassName,centrality)

  CreateEventHistos(kHistoForData | kHistoForMCInput,eventSelection,triggerClassName,centrality,ClassName(),ClassName(),1,0.0,1.0);
}

//_____________________________________________________________________________
void
AliAnalysisMuMuBase::CreateTrackHistos(UInt_t dataType,
                                       const char* eventSelection,
                                       const char* triggerClassName,
                                       const char* centrality,
                                       const char* hname, const char* htitle,
                                       Int_t nbinsx, Double_t xmin, Double_t xmax,
                                       Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  /// Create n copies of an histogram (with name=hname, title=htitle, etc..)
  /// where n = # of track cut combinations defined
  /// see also CreateHistos


  if ( IsHistogramDisabled(hname) ) return;

  TObjArray pathNames;
  pathNames.SetOwner(kTRUE);

  TIter nextCutCombination(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kTrack));
  AliAnalysisMuMuCutCombination* cutCombination;

  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextCutCombination())) )
  {
    if ( dataType & kHistoForData )
    {
      pathNames.Add(new TObjString(Form("/%s/%s/%s/%s",eventSelection,triggerClassName,centrality,cutCombination->GetName())));
    }
    if ( ( dataType & kHistoForMCInput ) && HasMC() )
    {
      pathNames.Add(new TObjString(Form("/%s/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,centrality,cutCombination->GetName())));
    }

  }

  CreateHistos(pathNames,hname,htitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
}

//_____________________________________________________________________________
void
AliAnalysisMuMuBase::CreatePairTHnSparse(UInt_t dataType,
                                      const char* eventSelection,
                                      const char* triggerClassName,
                                      const char* centrality,
                                      const char* hname, const char* htitle,
                                      Int_t nDim, Int_t* nbinsx, Double_t* xmin, Double_t* xmax) const
{
  /// Create n copies of an histogram (with name=hname, title=htitle, etc..)
  /// where n = # of track *pair* cut combinations defined
  /// see also CreateHistos

  if ( IsHistogramDisabled(hname) ) return;

  TObjArray pathNames;
  pathNames.SetOwner(kTRUE);

  TIter nextCutCombination(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kTrackPair));
  AliAnalysisMuMuCutCombination* cutCombination;

  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextCutCombination())) )
  {
    if ( dataType & kHistoForData )
    {
      pathNames.Add(new TObjString(Form("/%s/%s/%s/%s",eventSelection,triggerClassName,centrality,cutCombination->GetName())));
    }
    if ( ( dataType & kHistoForMCInput ) && HasMC() )
    {
      pathNames.Add(new TObjString(Form("/%s/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,centrality,cutCombination->GetName())));
    }
  }

  TIter next(&pathNames);
  TObjString* pathName;

  while ( ( pathName = static_cast<TObjString*>(next()) ) )
  {
    THnSparse* h(0x0);
    h = new THnSparseT<TArrayF>(hname,htitle,nDim,nbinsx,xmin,xmax);
    h->Sumw2();

    if( HistogramCollection()->Adopt(pathName->String().Data(),h))
    printf("%s/%s adopted\n",pathName->String().Data(),h->GetName() );
  }
}

//_____________________________________________________________________________
void
AliAnalysisMuMuBase::CreateTrackTHnSparse(UInt_t dataType,
                                      const char* eventSelection,
                                      const char* triggerClassName,
                                      const char* centrality,
                                      const char* hname, const char* htitle,
                                      Int_t nDim, Int_t* nbinsx, Double_t* xmin, Double_t* xmax) const
{
  /// Create n copies of an histogram (with name=hname, title=htitle, etc..)
  /// where n = # of track  cut combinations defined
  /// see also CreateHistos

  if ( IsHistogramDisabled(hname) ) return;

  TObjArray pathNames;
  pathNames.SetOwner(kTRUE);

  TIter nextCutCombination(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kTrack));
  AliAnalysisMuMuCutCombination* cutCombination;

  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextCutCombination())) )
  {
    if ( dataType & kHistoForData )
    {
      pathNames.Add(new TObjString(Form("/%s/%s/%s/%s",eventSelection,triggerClassName,centrality,cutCombination->GetName())));
    }
    if ( ( dataType & kHistoForMCInput ) && HasMC() )
    {
      pathNames.Add(new TObjString(Form("/%s/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,centrality,cutCombination->GetName())));
    }
  }

  TIter next(&pathNames);
  TObjString* pathName;

  while ( ( pathName = static_cast<TObjString*>(next()) ) )
  {
    THnSparse* h(0x0);
    h = new THnSparseT<TArrayF>(hname,htitle,nDim,nbinsx,xmin,xmax);
    h->Sumw2();

    if( HistogramCollection()->Adopt(pathName->String().Data(),h))
    printf("%s/%s adopted\n",pathName->String().Data(),h->GetName() );
  }
}

//_____________________________________________________________________________
void
AliAnalysisMuMuBase::CreatePairHistos(UInt_t dataType,
                                      const char* eventSelection,
                                      const char* triggerClassName,
                                      const char* centrality,
                                      const char* hname, const char* htitle,
                                      Int_t nbinsx, Double_t xmin, Double_t xmax,
                                      Int_t nbinsy, Double_t ymin, Double_t ymax) const
{
  /// Create n copies of an histogram (with name=hname, title=htitle, etc..)
  /// where n = # of track *pair* cut combinations defined
  /// see also CreateHistos

  if ( IsHistogramDisabled(hname) ) return;

  TObjArray pathNames;
  pathNames.SetOwner(kTRUE);

  TIter nextCutCombination(CutRegistry()->GetCutCombinations(AliAnalysisMuMuCutElement::kTrackPair));
  AliAnalysisMuMuCutCombination* cutCombination;

  while ( ( cutCombination = static_cast<AliAnalysisMuMuCutCombination*>(nextCutCombination())) )
  {
    if ( dataType & kHistoForData )
    {
      pathNames.Add(new TObjString(Form("/%s/%s/%s/%s",eventSelection,triggerClassName,centrality,cutCombination->GetName())));
    }
    if ( ( dataType & kHistoForMCInput ) && HasMC() )
    {
      pathNames.Add(new TObjString(Form("/%s/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,centrality,cutCombination->GetName())));
    }
  }

  CreateHistos(pathNames,hname,htitle,nbinsx,xmin,xmax,nbinsy,ymin,ymax);
}

//_____________________________________________________________________________
void AliAnalysisMuMuBase::DisableHistograms(const char* pattern)
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
Bool_t AliAnalysisMuMuBase::ExistSemaphoreHistogram(const char* eventSelection,
                                                    const char* triggerClassName,
                                                    const char* centrality) const
{
  /// Test for the existence of the semaphore histogram
  /// @see CreateSemaphoreHistogram

  return ( HistogramCollection()->Histo(Form("/%s/%s/%s/%s",eventSelection,triggerClassName,centrality,ClassName())) != 0x0 );
}

//_____________________________________________________________________________
Int_t AliAnalysisMuMuBase::GetNbins(Double_t xmin, Double_t xmax, Double_t xstep)
{
  /// Compute number of bins to get a given step between each values in the xmin,xmax range
  if ( TMath::AreEqualRel(xstep,0.0,1E-9) ) return 1;

  return TMath::Nint(TMath::Abs((xmax-xmin)/xstep));
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuBase::Histo(const char* eventSelection, const char* triggerClassName, const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s",eventSelection,triggerClassName,histoname)) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuBase::Histo(const char* eventSelection, const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(eventSelection,histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuBase::Histo(const char* eventSelection,
                                const char* triggerClassName,
                                const char* cent,
                                const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s",eventSelection,triggerClassName,cent),histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuBase::Histo(const char* eventSelection,
                                const char* triggerClassName,
                                const char* cent,
                                const char* what,
                                const char* histoname)
{
  /// Get one histo back

  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s/%s",eventSelection,triggerClassName,cent,what),histoname) : 0x0;
}

//_____________________________________________________________________________
TProfile* AliAnalysisMuMuBase::Prof(const char* eventSelection,
                                    const char* histoname)
{
	/// Get one histo profile back

	return fHistogramCollection ? static_cast<TProfile*>(fHistogramCollection->GetObject(Form("/%s",eventSelection),histoname)) : 0x0;
}

//_____________________________________________________________________________
TProfile* AliAnalysisMuMuBase::Prof(const char* eventSelection,
                                    const char* triggerClassName,
                                    const char* histoname)
{
	/// Get one histo profile back

	return fHistogramCollection ? static_cast<TProfile*>(fHistogramCollection->GetObject(Form("/%s/%s",eventSelection,triggerClassName),histoname)) : 0x0;
}

//_____________________________________________________________________________
TProfile* AliAnalysisMuMuBase::Prof(const char* eventSelection,
                                    const char* triggerClassName,
                                    const char* cent,
                                    const char* histoname)
{
	/// Get one histo profile back

	return fHistogramCollection ? static_cast<TProfile*>(fHistogramCollection->GetObject(Form("/%s/%s/%s",eventSelection,triggerClassName,cent),histoname)) : 0x0;
}

//_____________________________________________________________________________
TProfile* AliAnalysisMuMuBase::Prof(const char* eventSelection,
                                    const char* triggerClassName,
                                    const char* cent,
                                    const char* what,
                                    const char* histoname)
{
	/// Get one histo profile back

	return fHistogramCollection ? static_cast<TProfile*>(fHistogramCollection->GetObject(Form("/%s/%s/%s/%s",eventSelection,triggerClassName,cent,what),histoname)) : 0x0;
}

//_____________________________________________________________________________
void AliAnalysisMuMuBase::Init(AliCounterCollection& cc,
                               AliMergeableCollection& hc,
                               const AliAnalysisMuMuBinning& binning,
                               const AliAnalysisMuMuCutRegistry& registry)
{
  /// Set the internal references
  fEventCounters       = &cc;
  fHistogramCollection = &hc;
  fBinning             = &binning;
  fCutRegistry         = &registry;
}

//_____________________________________________________________________________
Bool_t AliAnalysisMuMuBase::IsHistogramDisabled(const char* hname) const
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
Bool_t AliAnalysisMuMuBase::IsHistogrammingDisabled() const
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
TH1* AliAnalysisMuMuBase::MCHisto(const char* eventSelection, const char* triggerClassName, const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,histoname)) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuBase::MCHisto(const char* eventSelection, const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s",MCInputPrefix(),eventSelection,histoname)) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuBase::MCHisto(const char* eventSelection,
                                  const char* triggerClassName,
                                  const char* cent,
                                  const char* histoname)
{
  /// Get one histo back
  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,cent),histoname) : 0x0;
}

//_____________________________________________________________________________
TH1* AliAnalysisMuMuBase::MCHisto(const char* eventSelection,
                                  const char* triggerClassName,
                                  const char* cent,
                                  const char* what,
                                  const char* histoname)
{
  /// Get one histo back

  return fHistogramCollection ? fHistogramCollection->Histo(Form("/%s/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,cent,what),histoname) : 0x0;
}

//_____________________________________________________________________________
TProfile* AliAnalysisMuMuBase::MCProf(const char* eventSelection,
                                    const char* histoname)
{
	/// Get one histo profile back

	return fHistogramCollection ? static_cast<TProfile*>(fHistogramCollection->GetObject(Form("/%s/%s",MCInputPrefix(),eventSelection),histoname)) : 0x0;
}

//_____________________________________________________________________________
TProfile* AliAnalysisMuMuBase::MCProf(const char* eventSelection,
                                    const char* triggerClassName,
                                    const char* histoname)
{
	/// Get one histo profile back

	return fHistogramCollection ? static_cast<TProfile*>(fHistogramCollection->GetObject(Form("/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName),histoname)) : 0x0;
}

//_____________________________________________________________________________
TProfile* AliAnalysisMuMuBase::MCProf(const char* eventSelection,
                                    const char* triggerClassName,
                                    const char* cent,
                                    const char* histoname)
{
	/// Get one histo profile back

	return fHistogramCollection ? static_cast<TProfile*>(fHistogramCollection->GetObject(Form("/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,cent),histoname)) : 0x0;
}

//_____________________________________________________________________________
TProfile* AliAnalysisMuMuBase::MCProf(const char* eventSelection,
                                    const char* triggerClassName,
                                    const char* cent,
                                    const char* what,
                                    const char* histoname)
{
	/// Get one histo profile back

	return fHistogramCollection ? static_cast<TProfile*>(fHistogramCollection->GetObject(Form("/%s/%s/%s/%s/%s",MCInputPrefix(),eventSelection,triggerClassName,cent,what),histoname)) : 0x0;
}

//_____________________________________________________________________________
void AliAnalysisMuMuBase::SetEvent(AliVEvent* event, AliMCEvent* mcEvent)
{
  /// Set pointers to event (and mcEvent)
  /// This is the place where a sub analysis might (if really needed)
  /// append things to the event. Please DO NOT alter the original content,
  /// only add to it. Unless of course you want to shoot yourself in the
  /// foot...

  fEvent = event;
  fMCEvent = mcEvent;
}
