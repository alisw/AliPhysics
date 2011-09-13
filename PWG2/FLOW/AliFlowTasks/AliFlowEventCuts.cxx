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

/* $Id$ */ 

// AliFlowEventCuts:
// An event cut class for the flow framework
//
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#include <limits.h>
#include <float.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TBrowser.h>
#include "TMath.h"
#include "TNamed.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"
#include "AliMCEvent.h"
#include "AliFlowEventCuts.h"
#include "AliFlowTrackCuts.h"
#include "AliTriggerAnalysis.h"

ClassImp(AliFlowEventCuts)

//-----------------------------------------------------------------------
AliFlowEventCuts::AliFlowEventCuts():
  TNamed(),
  fQA(NULL),
  fCutNumberOfTracks(kFALSE),
  fNumberOfTracksMax(INT_MAX),
  fNumberOfTracksMin(INT_MIN),
  fCutRefMult(kFALSE),
  fRefMultMethod(kTPConly),
  fUseAliESDtrackCutsRefMult(kFALSE),
  fRefMultMethodAliESDtrackCuts(AliESDtrackCuts::kTrackletsITSTPC),
  fRefMultMax(INT_MAX),
  fRefMultMin(INT_MIN),
  fRefMultCuts(NULL),
  fMeanPtCuts(NULL),
  fStandardTPCcuts(NULL),
  fStandardGlobalCuts(NULL),
  fCutPrimaryVertexX(kFALSE),
  fPrimaryVertexXmax(INT_MAX),
  fPrimaryVertexXmin(INT_MIN),
  fCutPrimaryVertexY(kFALSE),
  fPrimaryVertexYmax(INT_MAX),
  fPrimaryVertexYmin(INT_MIN),
  fCutPrimaryVertexZ(kFALSE),
  fPrimaryVertexZmax(INT_MAX),
  fPrimaryVertexZmin(INT_MIN),
  fCutNContributors(kFALSE),
  fNContributorsMax(INT_MAX),
  fNContributorsMin(INT_MIN),
  fCutMeanPt(kFALSE),
  fMeanPtMax(-DBL_MAX),
  fMeanPtMin(DBL_MAX),
  fCutSPDvertexerAnomaly(kFALSE),
  fCutTPCmultiplicityOutliers(kFALSE),
  fCutCentralityPercentile(kFALSE),
  fUseCentralityUnchecked(kFALSE),
  fCentralityPercentileMethod(kTPConly),
  fCentralityPercentileMax(100.),
  fCentralityPercentileMin(0.),
  fCutZDCtiming(kFALSE),
  fTrigAna()
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowEventCuts::AliFlowEventCuts(const char* name, const char* title):
  TNamed(name, title),
  fQA(NULL),
  fCutNumberOfTracks(kFALSE),
  fNumberOfTracksMax(INT_MAX),
  fNumberOfTracksMin(INT_MIN),
  fCutRefMult(kFALSE),
  fRefMultMethod(kTPConly),
  fUseAliESDtrackCutsRefMult(kFALSE),
  fRefMultMethodAliESDtrackCuts(AliESDtrackCuts::kTrackletsITSTPC),
  fRefMultMax(INT_MAX),
  fRefMultMin(INT_MIN),
  fRefMultCuts(NULL),
  fMeanPtCuts(NULL),
  fStandardTPCcuts(AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010()),
  fStandardGlobalCuts(AliFlowTrackCuts::GetStandardGlobalTrackCuts2010()),
  fCutPrimaryVertexX(kFALSE),
  fPrimaryVertexXmax(INT_MAX),
  fPrimaryVertexXmin(INT_MIN),
  fCutPrimaryVertexY(kFALSE),
  fPrimaryVertexYmax(INT_MAX),
  fPrimaryVertexYmin(INT_MIN),
  fCutPrimaryVertexZ(kFALSE),
  fPrimaryVertexZmax(INT_MAX),
  fPrimaryVertexZmin(INT_MIN),
  fCutNContributors(kFALSE),
  fNContributorsMax(INT_MAX),
  fNContributorsMin(INT_MIN),
  fCutMeanPt(kFALSE),
  fMeanPtMax(-DBL_MAX),
  fMeanPtMin(DBL_MAX),
  fCutSPDvertexerAnomaly(kFALSE),
  fCutTPCmultiplicityOutliers(kFALSE),
  fCutCentralityPercentile(kFALSE),
  fUseCentralityUnchecked(kFALSE),
  fCentralityPercentileMethod(kTPConly),
  fCentralityPercentileMax(100.),
  fCentralityPercentileMin(0.),
  fCutZDCtiming(kFALSE),
  fTrigAna()
{
  //constructor 
}

////-----------------------------------------------------------------------
AliFlowEventCuts::AliFlowEventCuts(const AliFlowEventCuts& that):
  TNamed(that),
  fQA(NULL),
  fCutNumberOfTracks(that.fCutNumberOfTracks),
  fNumberOfTracksMax(that.fNumberOfTracksMax),
  fNumberOfTracksMin(that.fNumberOfTracksMin),
  fCutRefMult(that.fCutRefMult),
  fRefMultMethod(that.fRefMultMethod),
  fUseAliESDtrackCutsRefMult(that.fUseAliESDtrackCutsRefMult),
  fRefMultMethodAliESDtrackCuts(that.fRefMultMethodAliESDtrackCuts),
  fRefMultMax(that.fRefMultMax),
  fRefMultMin(that.fRefMultMin),
  fRefMultCuts(NULL),
  fMeanPtCuts(NULL),
  fStandardTPCcuts(NULL),
  fStandardGlobalCuts(NULL),
  fCutPrimaryVertexX(that.fCutPrimaryVertexX),
  fPrimaryVertexXmax(that.fPrimaryVertexXmax),
  fPrimaryVertexXmin(that.fPrimaryVertexXmin),
  fCutPrimaryVertexY(that.fCutPrimaryVertexX),
  fPrimaryVertexYmax(that.fPrimaryVertexYmax),
  fPrimaryVertexYmin(that.fPrimaryVertexYmin),
  fCutPrimaryVertexZ(that.fCutPrimaryVertexX),
  fPrimaryVertexZmax(that.fPrimaryVertexZmax),
  fPrimaryVertexZmin(that.fPrimaryVertexZmin),
  fCutNContributors(that.fCutNContributors),
  fNContributorsMax(that.fNContributorsMax),
  fNContributorsMin(that.fNContributorsMin),
  fCutMeanPt(that.fCutMeanPt),
  fMeanPtMax(that.fMeanPtMax),
  fMeanPtMin(that.fMeanPtMin),
  fCutSPDvertexerAnomaly(that.fCutSPDvertexerAnomaly),
  fCutTPCmultiplicityOutliers(that.fCutTPCmultiplicityOutliers),
  fCutCentralityPercentile(that.fCutCentralityPercentile),
  fUseCentralityUnchecked(that.fUseCentralityUnchecked),
  fCentralityPercentileMethod(that.fCentralityPercentileMethod),
  fCentralityPercentileMax(that.fCentralityPercentileMax),
  fCentralityPercentileMin(that.fCentralityPercentileMin),
  fCutZDCtiming(that.fCutZDCtiming),
  fTrigAna()
{
  if (that.fQA) DefineHistograms();
  //copy constructor 
  if (that.fRefMultCuts)
    fRefMultCuts = new AliFlowTrackCuts(*(that.fRefMultCuts));
  if (that.fMeanPtCuts)
    fMeanPtCuts = new AliFlowTrackCuts(*(that.fMeanPtCuts));
  fStandardTPCcuts = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();
  fStandardGlobalCuts = AliFlowTrackCuts::GetStandardGlobalTrackCuts2010();
}

////-----------------------------------------------------------------------
AliFlowEventCuts::~AliFlowEventCuts()
{
  //dtor
  delete fMeanPtCuts;
  delete fRefMultCuts;
  delete fStandardGlobalCuts;
  delete fStandardTPCcuts;
  if (fQA) { fQA->SetOwner(); fQA->Delete(); delete fQA; }
}

////-----------------------------------------------------------------------
AliFlowEventCuts& AliFlowEventCuts::operator=(const AliFlowEventCuts& that)
{
  //assignment
  if (this==&that) return *this;

  if (that.fQA)
  {
    if (fQA)
    {
      fQA->Delete();
      delete fQA;
    }
    fQA = static_cast<TList*>(that.fQA->Clone());
  }
  else
  {
    fQA->Delete();
    delete fQA;
    fQA=NULL;
  }

  fCutNumberOfTracks=that.fCutNumberOfTracks;
  fNumberOfTracksMax=that.fNumberOfTracksMax;
  fNumberOfTracksMin=that.fNumberOfTracksMin;
  fCutRefMult=that.fCutRefMult;
  fRefMultMethod=that.fRefMultMethod;
  fUseAliESDtrackCutsRefMult=that.fUseAliESDtrackCutsRefMult;
  fRefMultMethodAliESDtrackCuts=that.fRefMultMethodAliESDtrackCuts;
  fRefMultMax=that.fRefMultMax;
  fRefMultMin=that.fRefMultMin;
  if (that.fRefMultCuts) *fRefMultCuts=*(that.fRefMultCuts);
  if (that.fMeanPtCuts) *fMeanPtCuts=*(that.fMeanPtCuts);
  fStandardTPCcuts = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();
  fStandardGlobalCuts = AliFlowTrackCuts::GetStandardGlobalTrackCuts2010();
  fCutPrimaryVertexX=that.fCutPrimaryVertexX;
  fPrimaryVertexXmax=that.fPrimaryVertexXmax;
  fPrimaryVertexXmin=that.fPrimaryVertexXmin;
  fCutPrimaryVertexY=that.fCutPrimaryVertexY;
  fPrimaryVertexYmax=that.fPrimaryVertexYmax;
  fPrimaryVertexYmin=that.fPrimaryVertexYmin;
  fCutPrimaryVertexZ=that.fCutPrimaryVertexZ;
  fPrimaryVertexZmax=that.fPrimaryVertexZmax;
  fPrimaryVertexZmin=that.fPrimaryVertexZmin;
  fCutNContributors=that.fCutNContributors;
  fNContributorsMax=that.fNContributorsMax;
  fNContributorsMin=that.fNContributorsMin;
  fCutMeanPt=that.fCutMeanPt;
  fMeanPtMax=that.fMeanPtMax;
  fMeanPtMin=that.fMeanPtMin;
  fCutSPDvertexerAnomaly=that.fCutSPDvertexerAnomaly;
  fCutTPCmultiplicityOutliers=that.fCutTPCmultiplicityOutliers;
  fCutCentralityPercentile=that.fCutCentralityPercentile;
  fUseCentralityUnchecked=that.fUseCentralityUnchecked;
  fCentralityPercentileMethod=that.fCentralityPercentileMethod;
  fCentralityPercentileMax=that.fCentralityPercentileMax;
  fCentralityPercentileMin=that.fCentralityPercentileMin;
  fCutZDCtiming=that.fCutZDCtiming;
  return *this;
}

//----------------------------------------------------------------------- 
Bool_t AliFlowEventCuts::IsSelected(TObject* obj)
{
  //check cuts
  AliVEvent* vevent = dynamic_cast<AliVEvent*>(obj);
  if (vevent) return PassesCuts(vevent);
  return kFALSE;  //when passed wrong type of object
}
//----------------------------------------------------------------------- 
Bool_t AliFlowEventCuts::PassesCuts(AliVEvent *event)
{
  ///check if event passes cuts
  const AliVVertex* pvtx=event->GetPrimaryVertex();
  Double_t pvtxx = pvtx->GetX();
  Double_t pvtxy = pvtx->GetY();
  Double_t pvtxz = pvtx->GetZ();
  Int_t ncontrib = pvtx->GetNContributors();
  Bool_t pass=kTRUE;
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*>(event);
  Int_t multTPC = 0;
  Int_t multGlobal = 0; 
  if (fQA)
  {
    multTPC = fStandardTPCcuts->Count(event);
    multGlobal = fStandardGlobalCuts->Count(event);
    QAbefore(0)->Fill(pvtxz);
    QAbefore(1)->Fill(multGlobal,multTPC);
  }
  if (fCutTPCmultiplicityOutliers)
  {
    //this is pretty slow as we check the event track by track twice
    //this cut will work for 2010 PbPb data and is dependent on
    //TPC and ITS reco efficiency (e.g. geometry, calibration etc)
    if (!fQA)
    {
      multTPC = fStandardTPCcuts->Count(event);
      multGlobal = fStandardGlobalCuts->Count(event);
    }
    if (multTPC > ( 23+1.216*multGlobal)) {pass=kFALSE;}
    if (multTPC < (-20+1.087*multGlobal)) {pass=kFALSE;}
  }
  if (fCutNContributors)
  {
    if (ncontrib < fNContributorsMin || ncontrib >= fNContributorsMax) pass=kFALSE;
  }
  if (fCutPrimaryVertexX)
  {
    if (pvtxx < fPrimaryVertexXmin || pvtxx >= fPrimaryVertexXmax) pass=kFALSE;
  }
  if (fCutPrimaryVertexY)
  {
    if (pvtxy < fPrimaryVertexYmin || pvtxy >= fPrimaryVertexYmax) pass=kFALSE;
  }
  if (fCutPrimaryVertexZ)
  {
    if (pvtxz < fPrimaryVertexZmin || pvtxz >= fPrimaryVertexZmax)
      pass=kFALSE;
  }
  if (fCutCentralityPercentile&&esdevent)
  {
    AliCentrality* centr = esdevent->GetCentrality();
    if (fUseCentralityUnchecked)
    {
      if (!centr->IsEventInCentralityClassUnchecked( fCentralityPercentileMin,
                                                     fCentralityPercentileMax,
                                                     CentrMethName(fCentralityPercentileMethod) ))
      {
        pass=kFALSE;
      }
    }
    else
    {
      if (!centr->IsEventInCentralityClass( fCentralityPercentileMin,
                                            fCentralityPercentileMax,
                                            CentrMethName(fCentralityPercentileMethod) ))
      {
        pass=kFALSE;
      }
    }
  }
  if (fCutSPDvertexerAnomaly&&esdevent)
  {
    const AliESDVertex* sdpvertex = esdevent->GetPrimaryVertexSPD();
    if (sdpvertex->GetNContributors()<1) pass=kFALSE;
    if (sdpvertex->GetDispersion()>0.04) pass=kFALSE;
    if (sdpvertex->GetZRes()>0.25) pass=kFALSE;
    const AliESDVertex* tpcvertex = esdevent->GetPrimaryVertexTPC();
    if (tpcvertex->GetNContributors()<1) pass=kFALSE;
    const AliMultiplicity* tracklets = esdevent->GetMultiplicity();
    if (tpcvertex->GetNContributors()<(-10.0+0.25*tracklets->GetNumberOfITSClusters(0)))
    {
      pass=kFALSE;
    }
  }
  if (fCutZDCtiming&&esdevent)
  {
    if (!fTrigAna.ZDCTimeTrigger(esdevent))
    {
      pass=kFALSE;
    }
  }
  if(fCutNumberOfTracks) {if ( event->GetNumberOfTracks() < fNumberOfTracksMin ||
                               event->GetNumberOfTracks() >= fNumberOfTracksMax ) pass=kFALSE;}
  if(fCutRefMult&&esdevent)
  {
    //reference multiplicity still to be defined
    Double_t refMult = RefMult(event);
    if (refMult < fRefMultMin || refMult >= fRefMultMax )
    {
      pass=kFALSE;
    }
  }
  if (fCutMeanPt)
  {
    Float_t meanpt=0.0;
    Int_t ntracks=event->GetNumberOfTracks();
    Int_t nselected=0;
    for (Int_t i=0; i<ntracks; i++)
    {
      AliVParticle* track = event->GetTrack(i);
      if (!track) continue;
      Bool_t localpass=kTRUE;
      if (fMeanPtCuts) localpass=fMeanPtCuts->IsSelected(track);
      if (localpass) 
      {
        meanpt += track->Pt();
        nselected++;
      }
    }
    meanpt=meanpt/nselected;
    if (meanpt<fMeanPtMin || meanpt >= fMeanPtMax) pass=kFALSE;
  }
  if (fQA&&pass) 
  {
    QAafter(1)->Fill(multGlobal,multTPC);
    QAafter(0)->Fill(pvtxz);
  }
  return pass;
}

//----------------------------------------------------------------------- 
const char* AliFlowEventCuts::CentrMethName(refMultMethod method) const
{
  //get the string for refmultmethod, for use with AliCentrality in
  //the cut on centrality percentile
  switch (method)
  {
    case kSPDtracklets:
      return "TKL";
    case kSPD1clusters:
      return "CL1";
    case kTPConly:
      return "TRK";
    case kV0:
      return "V0M";
    default:
      return "";
  }
}
//----------------------------------------------------------------------- 
AliFlowEventCuts* AliFlowEventCuts::StandardCuts()
{
  //make a set of standard event cuts, caller becomes owner
  AliFlowEventCuts* cuts = new AliFlowEventCuts();
  return cuts;
}

//----------------------------------------------------------------------- 
Int_t AliFlowEventCuts::RefMult(AliVEvent* event)
{
  //calculate the reference multiplicity, if all fails return 0
  AliESDVZERO* vzero = NULL;
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*>(event);

  if (fUseAliESDtrackCutsRefMult && esdevent)
  {
    //use the standard ALICE reference multiplicity with the default eta range
    return AliESDtrackCuts::GetReferenceMultiplicity(esdevent, fRefMultMethodAliESDtrackCuts);
  }

  if (fRefMultMethod==kTPConly && !fRefMultCuts)
  {
    fRefMultCuts = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts();
    fRefMultCuts->SetEtaRange(-0.8,0.8);
    fRefMultCuts->SetPtMin(0.15);
  }
  else if (fRefMultMethod==kSPDtracklets && !fRefMultCuts)
  {
    fRefMultCuts = new AliFlowTrackCuts("tracklet refmult cuts");
    fRefMultCuts->SetParamType(AliFlowTrackCuts::kSPDtracklet);
    fRefMultCuts->SetEtaRange(-0.8,0.8);
  }
  else if (fRefMultMethod==kV0)
  {
    if (!esdevent) return 0;
    vzero=esdevent->GetVZEROData();
    if (!vzero) return 0;
    return TMath::Nint(vzero->GetMTotV0A()+vzero->GetMTotV0C());
  }
  else if (fRefMultMethod==kSPD1clusters)
  {
    if (!esdevent) return 0;
    const AliMultiplicity* mult = esdevent->GetMultiplicity();
    if (!mult) return 0;
    return mult->GetNumberOfITSClusters(1);
  }

  Int_t refmult=0;
  fRefMultCuts->SetEvent(event);
  for (Int_t i=0; i<fRefMultCuts->GetNumberOfInputObjects(); i++)
  {
    if (fRefMultCuts->IsSelected(fRefMultCuts->GetInputObject(i),i))
      refmult++;
  }
  return refmult;
}
//_____________________________________________________________________________
void AliFlowEventCuts::DefineHistograms()
{
  //define QA histos
  if (fQA) return;

  Bool_t adddirstatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fQA = new TList(); fQA->SetOwner();
  fQA->SetName(Form("%s QA",GetName()));
  TList* before = new TList(); before->SetOwner();
  before->SetName("before");
  TList* after = new TList(); after->SetOwner();
  after->SetName("after");
  fQA->Add(before);
  fQA->Add(after);
  before->Add(new TH1F("zvertex",";z;event cout",500,-15.,15.)); //0
  after->Add(new TH1F("zvertex",";z;event cout",500,-15.,15.)); //0
  before->Add(new TH2F("fTPCvsGlobalMult","TPC only vs Global track multiplicity;global;TPC only",500,0,2500,500,0,3500));//1
  after->Add(new TH2F("fTPCvsGlobalMult","TPC only vs Global track multiplicity;global;TPC only",500,0,2500,500,0,3500));//1
  TH1::AddDirectory(adddirstatus);
}

//---------------------------------------------------------------//
void AliFlowEventCuts::Browse(TBrowser* b)
{
  //some browsing capabilities
  if (fQA) b->Add(fQA);
}

//---------------------------------------------------------------//
Long64_t AliFlowEventCuts::Merge(TCollection* list)
{
  //merge
  Int_t number=0;
  AliFlowEventCuts* obj;
  if (!list) return 0;
  if (list->GetEntries()<1) return 0;
  TIter next(list);
  while ( (obj = dynamic_cast<AliFlowEventCuts*>(next())) )
  {
    if (obj==this) continue;
    TList listwrapper;
    listwrapper.Add(obj->GetQA());
    fQA->Merge(&listwrapper);
    number++;
  }
  return number;
}

