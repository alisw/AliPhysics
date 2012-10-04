/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliMuonEventCuts.h"

#include "TMath.h"
#include "THashList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TFile.h"
#include "TParameter.h"
#include "TKey.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TArrayI.h"

#include "AliLog.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliCentrality.h"
#include "AliTimeStamp.h"

#include "AliAnalysisMuonUtility.h"

/// \cond CLASSIMP
ClassImp(AliMuonEventCuts) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliMuonEventCuts::AliMuonEventCuts() :
  AliAnalysisCuts(),
  fPhysicsSelectionMask(0),
  fVertexMinNContributors(0),
  fVertexVzMin(0.),
  fVertexVzMax(0.),
  fDefaultTrigClassPatterns(""),
  fSelectedTrigPattern(0x0),
  fRejectedTrigPattern(0x0),
  fSelectedTrigLevel(0x0),
  fAllSelectedTrigClasses(0x0),
  fCentralityClasses(0x0),
  fTimeStamp(0x0),
  fSelectedTrigClassesInEvent(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliMuonEventCuts::AliMuonEventCuts(const char* name, const char* title ) :
AliAnalysisCuts(name, title),
  fPhysicsSelectionMask(0),
  fVertexMinNContributors(0),
  fVertexVzMin(0.),
  fVertexVzMax(0.),
  fDefaultTrigClassPatterns(""),
  fSelectedTrigPattern(new TObjArray()),
  fRejectedTrigPattern(new TObjArray()),
  fSelectedTrigLevel(new TObjArray()),
  fAllSelectedTrigClasses(new THashList()),
  fCentralityClasses(0x0),
  fTimeStamp(0x0),
  fSelectedTrigClassesInEvent(new TObjArray())
{
  /// Constructor
  SetDefaultParameters();
  SetDefaultFilterMask();
  SetDefaultTrigClassPatterns();
  SetTrigClassLevels();
  SetCentralityClasses();
  fAllSelectedTrigClasses->SetOwner();
  fSelectedTrigClassesInEvent->SetOwner();
}

//________________________________________________________________________
AliMuonEventCuts::AliMuonEventCuts(const AliMuonEventCuts& obj) :
  AliAnalysisCuts(obj),
  fPhysicsSelectionMask(obj.fPhysicsSelectionMask),
  fVertexMinNContributors(obj.fVertexMinNContributors),
  fVertexVzMin(obj.fVertexVzMin),
  fVertexVzMax(obj.fVertexVzMax),
  fDefaultTrigClassPatterns(obj.fDefaultTrigClassPatterns),
  fSelectedTrigPattern(obj.fSelectedTrigPattern),
  fRejectedTrigPattern(obj.fRejectedTrigPattern),
  fSelectedTrigLevel(obj.fSelectedTrigLevel),
  fAllSelectedTrigClasses(obj.fAllSelectedTrigClasses),
  fCentralityClasses(obj.fCentralityClasses),
  fTimeStamp(obj.fTimeStamp),
  fSelectedTrigClassesInEvent(obj.fSelectedTrigClassesInEvent)
{
  /// Copy constructor
}


//________________________________________________________________________
AliMuonEventCuts& AliMuonEventCuts::operator=(const AliMuonEventCuts& obj)
{
  /// Assignment operator
  if ( this != &obj ) { 
    AliAnalysisCuts::operator=(obj);
    fPhysicsSelectionMask = obj.fPhysicsSelectionMask;
    fVertexMinNContributors = obj.fVertexMinNContributors,
    fVertexVzMin = obj.fVertexVzMin;
    fVertexVzMax = obj.fVertexVzMax;
    fDefaultTrigClassPatterns = obj.fDefaultTrigClassPatterns;
    fSelectedTrigPattern = obj.fSelectedTrigPattern;
    fRejectedTrigPattern = obj.fRejectedTrigPattern;
    fSelectedTrigLevel = obj.fSelectedTrigLevel;
    fAllSelectedTrigClasses = obj.fAllSelectedTrigClasses;
    fCentralityClasses = obj.fCentralityClasses;
    fTimeStamp = obj.fTimeStamp;
    fSelectedTrigClassesInEvent = obj.fSelectedTrigClassesInEvent;
  }
  return *this;
}


//________________________________________________________________________
AliMuonEventCuts::~AliMuonEventCuts()
{
  /// Destructor
  delete fSelectedTrigPattern;
  delete fRejectedTrigPattern;
  delete fSelectedTrigLevel;
  delete fAllSelectedTrigClasses;
  delete fSelectedTrigClassesInEvent;
  delete fCentralityClasses;
  delete fTimeStamp;
}

//________________________________________________________________________
Bool_t AliMuonEventCuts::IsSelected( TObject* obj )
{
  /// Track is selected
  UInt_t filterMask = GetFilterMask();
  UInt_t selectionMask = GetSelectionMask(obj);
  
  return ( ( selectionMask & filterMask ) == filterMask );
}


//________________________________________________________________________
UInt_t AliMuonEventCuts::GetSelectionMask( const TObject* obj )
{
  /// Get selection mask
  
  UInt_t selectionMask = 0;
  
  const AliInputEventHandler* inputHandler = static_cast<const AliInputEventHandler*> ( obj );
  
  if ( const_cast<AliInputEventHandler*>(inputHandler)->IsEventSelected() & fPhysicsSelectionMask ) selectionMask |= kPhysicsSelected;
  
  const AliVEvent* event = inputHandler->GetEvent();
  
  Double_t centrality = GetCentrality(event);
  if ( centrality >= fCentralityClasses->GetXmin() && centrality <= fCentralityClasses->GetXmax() ) selectionMask |= kSelectedCentrality;
  
  UpdateEvent(event);
  
  if ( fSelectedTrigClassesInEvent->GetEntries() > 0 ) selectionMask |= kSelectedTrig;
  
  AliVVertex* vertex = AliAnalysisMuonUtility::GetVertexSPD(event);
  if ( vertex->GetNContributors() >= GetVertexMinNContributors() && 
      vertex->GetZ() >= GetVertexVzMin() && vertex->GetZ() <= GetVertexVzMax() ) selectionMask |= kGoodVertex;
  
  AliDebug(1, Form("Selection mask 0x%x\n", selectionMask));

  return selectionMask;
}


//________________________________________________________________________
Bool_t AliMuonEventCuts::IsSelected( TList* /* list */)
{
  /// Not implemented
  AliError("Function not implemented: Use IsSelected(TObject*)");
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliMuonEventCuts::UpdateEvent ( const AliVEvent* event )
{
  /// Update the transient data member per event
  
  AliTimeStamp currTimeStamp(event->GetOrbitNumber(), event->GetPeriodNumber(), event->GetBunchCrossNumber());
  if ( fTimeStamp && fTimeStamp->Compare(&currTimeStamp) == 0 ) return kFALSE;
  
  BuildTriggerClasses(AliAnalysisMuonUtility::GetFiredTriggerClasses(event));
  
  delete fTimeStamp;
  fTimeStamp = new AliTimeStamp(currTimeStamp);
  
  return kTRUE;
}

//________________________________________________________________________
void AliMuonEventCuts::SetDefaultTrigClassPatterns ()
{
  /// Set the default patterns
  /// (done in such a way to get all muon triggers)
  fDefaultTrigClassPatterns = "CINT CMU CMBAC CPBI !-ACE- !-AC- !-E- !WU !EGA !EJE !PHS";
  SetTrigClassPatterns(fDefaultTrigClassPatterns);
}

//________________________________________________________________________
void AliMuonEventCuts::SetTrigClassPatterns ( TString pattern )
{
  /// Set trigger classes
  ///
  /// Classes are filled dynamically according to the pattern
  /// - if name contains ! (without spaces): reject it
  /// - otherwise, keep it
  /// example:
  /// SetTrigClassPatterns("CMBAC !ALLNOTRD")
  /// keeps classes containing CMBAC, and not containing ALLNOTRD.
  ///
  /// CAVEAT: if you use an fCFContainer and you want an axis to contain the trigger classes,
  ///         please be sure that each pattern matches only 1 trigger class, or triggers will be mixed up
  ///         when merging different chuncks.
  
  fSelectedTrigPattern->SetOwner();
  if ( fSelectedTrigPattern->GetEntries() > 0 ) fSelectedTrigPattern->Delete();
  fRejectedTrigPattern->SetOwner();
  if ( fRejectedTrigPattern->GetEntries() > 0 ) fRejectedTrigPattern->Delete();
  
  pattern.ReplaceAll("  "," ");
  pattern.ReplaceAll("! ","!");
  
  TObjArray* fullList = pattern.Tokenize(" ");
  
  for ( Int_t ipat=0; ipat<fullList->GetEntries(); ++ipat ) {
    TString currPattern = fullList->At(ipat)->GetName();
    if ( currPattern.Contains("!") ) {
      currPattern.ReplaceAll("!","");
      fRejectedTrigPattern->AddLast(new TObjString(currPattern));
    }
    else fSelectedTrigPattern->AddLast(new TObjString(currPattern));
  }
  
  delete fullList;
}

//________________________________________________________________________
void AliMuonEventCuts::SetTrigClassLevels ( TString pattern )
{
  /// Set trigger cut level associated to the trigger class
  ///
  /// example:
  /// SetTrigClassLevels("MSL:Lpt MSH:Hpt MUL:LptLpt")
  ///
  /// For the trigger classes defined in SetTrigClassPatterns
  /// it check if they contains the keywords MSL or MSH
  /// Hence, in the analysis, the function
  /// TrackPtCutMatchTrigClass(track, "CPBIMSL") returns true if track match Lpt
  /// TrackPtCutMatchTrigClass(track, "CPBIMSH") returns true if track match Hpt
  /// TrackPtCutMatchTrigClass(track, "CMBAC") always returns true
  
  fSelectedTrigLevel->SetOwner();
  if ( fSelectedTrigLevel->GetEntries() > 0 ) fSelectedTrigLevel->Delete();
  
  pattern.ReplaceAll("  "," ");
  pattern.ReplaceAll(" :",":");
  
  TObjArray* fullList = pattern.Tokenize(" ");
  UInt_t offset = 2;
  for ( Int_t ipat=0; ipat<fullList->GetEntries(); ++ipat ) {
    TString currPattern = fullList->At(ipat)->GetName();
    TObjArray* arr = currPattern.Tokenize(":");
    TObjString* trigClassPattern = new TObjString(arr->At(0)->GetName());
    TString selTrigLevel = arr->At(1)->GetName();
    selTrigLevel.ToUpper();
    UInt_t trigLevel = 0;
    if ( selTrigLevel.Contains("LPT") ) {
      trigLevel = 2;
      if ( selTrigLevel.Contains("LPTLPT") ) trigLevel += 2<<offset;
    }
    else if ( selTrigLevel.Contains("HPT") ) {
      trigLevel = 3;
      if ( selTrigLevel.Contains("HPTHPT") ) trigLevel += 3<<offset;
    }
    trigClassPattern->SetUniqueID(trigLevel);
    fSelectedTrigLevel->AddLast(trigClassPattern);
    delete arr;
  }
  
  delete fullList;
}

//________________________________________________________________________
TArrayI AliMuonEventCuts::GetTrigClassPtCutLevel ( const TString trigClassName ) const
{
  /// Get trigger class pt cut level for tracking/trigger matching
  TObject* obj = fAllSelectedTrigClasses->FindObject(trigClassName.Data());
  if ( ! obj ) {
    AliWarning(Form("Class %s not in the list!", trigClassName.Data()));
    return -1;
  }
  
  TArrayI ptCutLevel(2);
  ptCutLevel.Reset();
  ptCutLevel[0] = obj->GetUniqueID() & 0x3;
  ptCutLevel[1] = ( obj->GetUniqueID() >> 2 ) & 0x3;
  
  AliDebug(1,Form("Class %s ptCutLevel %i %i",trigClassName.Data(),ptCutLevel[0],ptCutLevel[1]));
  
  return ptCutLevel;
}


//________________________________________________________________________
TObjArray* AliMuonEventCuts::GetSelectedTrigClassesInEvent( const AliVEvent* event )
{
  /// Return the selected trigger classes in the current event
  UpdateEvent(event);
  return fSelectedTrigClassesInEvent;
}


//________________________________________________________________________
void AliMuonEventCuts::BuildTriggerClasses ( const TString firedTrigClasses )
{
  //
  /// Return the list of trigger classes to be considered
  /// for current event. Update the global list if necessary
  //
  
  delete fSelectedTrigClassesInEvent;
  fSelectedTrigClassesInEvent = new TObjArray(0);
  fSelectedTrigClassesInEvent->SetOwner();

  TString firedTrigClassesAny = "ANY " + firedTrigClasses;
  TObjArray* firedTrigClassesList = firedTrigClassesAny.Tokenize(" ");
  
  UInt_t trigLevel = 0;
  for ( Int_t itrig=0; itrig<firedTrigClassesList->GetEntries(); ++itrig ) {
    TString trigName = ((TObjString*)firedTrigClassesList->At(itrig))->GetString();
    
    TObject* foundTrig = fAllSelectedTrigClasses->FindObject(trigName.Data());
    if ( foundTrig ) trigLevel = foundTrig->GetUniqueID();
    else {
      Bool_t rejectTrig = kFALSE;
      for ( Int_t ipat=0; ipat<fRejectedTrigPattern->GetEntries(); ++ipat ) {
        if ( trigName.Contains(fRejectedTrigPattern->At(ipat)->GetName() ) ) {
          rejectTrig = kTRUE;
          break;
        }
      } // loop on reject pattern
      if ( rejectTrig ) continue;
      
      rejectTrig = kTRUE;
      for ( Int_t ipat=0; ipat<fSelectedTrigPattern->GetEntries(); ++ipat ) {
        if ( trigName.Contains(fSelectedTrigPattern->At(ipat)->GetName() ) ) {
          rejectTrig = kFALSE;
          break;
        }
      } // loop on keep pattern
      if ( rejectTrig ) continue;
      
      trigLevel = 0;
      for ( Int_t ipat=0; ipat<fSelectedTrigLevel->GetEntries(); ++ipat ) {
        if ( trigName.Contains(fSelectedTrigLevel->At(ipat)->GetName() ) ) {
          trigLevel = fSelectedTrigLevel->At(ipat)->GetUniqueID();
          break;
        }
      } // loop on trig level patterns      
    }
    TObjString* currTrig = new TObjString(trigName);
    currTrig->SetUniqueID(trigLevel);
    fSelectedTrigClassesInEvent->AddLast(currTrig);
    
    if ( foundTrig ) continue;
    TObjString* addTrig = new TObjString(trigName);
    addTrig->SetUniqueID(trigLevel);
    fAllSelectedTrigClasses->Add(addTrig);
    TString trigLevelInfo = Form("trig level %i ", trigLevel & 0x3);
    trigLevelInfo += ( trigLevel > 3 ) ? "di-muon" : "single-muon";
    AliInfo(Form("Adding %s (%s) to considered trigger classes",trigName.Data(),trigLevelInfo.Data()));
  } // loop on trigger classes
  
  delete firedTrigClassesList;
}

//________________________________________________________________________
void AliMuonEventCuts::SetCentralityClasses(Int_t nCentralityBins, Double_t* centralityBins)
{
  //
  /// Set centrality classes
  //
  Double_t* bins = centralityBins;
  Int_t nbins = nCentralityBins;
  
  Double_t defaultCentralityBins[] = {-5., 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 100., 105.};
  if ( ! centralityBins ) {
    bins = defaultCentralityBins;
    nbins = sizeof(defaultCentralityBins)/sizeof(defaultCentralityBins[0])-1;
  }
  
  TString centralityEstimator = "V0M";
  if ( fCentralityClasses ) {
    centralityEstimator = GetCentralityEstimator();
    delete fCentralityClasses;
  }
  fCentralityClasses = new TAxis(nbins, bins);
  TString currClass = "";
  for ( Int_t ibin=1; ibin<=fCentralityClasses->GetNbins(); ++ibin ){
    currClass = Form("%.0f_%.0f",fCentralityClasses->GetBinLowEdge(ibin),fCentralityClasses->GetBinUpEdge(ibin));
    fCentralityClasses->SetBinLabel(ibin, currClass.Data());
  }
  
  SetCentralityEstimator(centralityEstimator);
}

//________________________________________________________________________
void AliMuonEventCuts::SetCentralityEstimator ( const TString centralityEstimator )
{
  /// Set centrality estimator
  fCentralityClasses->SetName(centralityEstimator.Data());
}


//________________________________________________________________________
TString AliMuonEventCuts::GetCentralityEstimator () const
{
  /// Get centrality estimator
  return fCentralityClasses->GetName();
}

//________________________________________________________________________
Double_t AliMuonEventCuts::GetCentrality ( const AliVEvent* event ) const
{
  /// Get centrality
  AliVEvent* evt = const_cast<AliVEvent*>(event);
  return evt->GetCentrality()->GetCentralityPercentile(GetCentralityEstimator());
}

//________________________________________________________________________
void AliMuonEventCuts::SetDefaultParameters ()
{
  /// Standard parameters for muon event
  SetPhysicsSelectionMask(AliVEvent::kAny);
  SetVertexMinNContributors(1);
  SetVertexVzLimits();
}


//________________________________________________________________________
void AliMuonEventCuts::SetDefaultFilterMask ()
{
  /// Standard cuts for muon event
  SetFilterMask ( kPhysicsSelected | kSelectedTrig | kGoodVertex );
}

//________________________________________________________________________
void AliMuonEventCuts::Print(Option_t* option) const
{
  //
  /// Print info
  //
  TString sopt(option);
  sopt.ToLower();
  if ( sopt.IsNull() || sopt.Contains("*") || sopt.Contains("all") ) sopt = "mask param";
  UInt_t filterMask = GetFilterMask();
  if ( sopt.Contains("mask") ) {
    printf(" *** Muon event filter mask: *** \n");
    printf("  0x%x\n", filterMask);
    if ( filterMask & kPhysicsSelected ) printf("  Pass physics selection 0x%x\n", fPhysicsSelectionMask);
    if ( filterMask & kSelectedCentrality ) printf(  "%g < centrality < %g", fCentralityClasses->GetXmin(), fCentralityClasses->GetXmax() );
    if ( filterMask & kSelectedTrig )    printf("  Has selected trigger classes\n");
    if ( filterMask & kGoodVertex )      printf("  SPD vertex with %i contributors && %g < Vz < %g\n", GetVertexMinNContributors(), GetVertexVzMin(), GetVertexVzMax());
    printf(" ******************** \n");
  }
}
