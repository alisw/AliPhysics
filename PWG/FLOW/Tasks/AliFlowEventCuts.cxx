 /* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliCentrality.h"
#include "AliMultSelection.h" // available from November 2015
#include "AliESDVZERO.h"
#include "AliMultiplicity.h"
#include "AliMCEvent.h"
#include "AliFlowEventCuts.h"
#include "AliFlowTrackCuts.h"
#include "AliTriggerAnalysis.h"
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include <iostream>
using namespace std;

ClassImp(AliFlowEventCuts)

//-----------------------------------------------------------------------
AliFlowEventCuts::AliFlowEventCuts():
  AliFlowEventSimpleCuts(),
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
  fUtils(NULL),
  fMultSelection(NULL),
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
  fCutSPDTRKVtxZ(kFALSE),
  fCutTPCmultiplicityOutliers(kFALSE),
  fCutTPCmultiplicityOutliersAOD(kFALSE),
  fUseCentralityUnchecked(kFALSE),
  fCentralityPercentileMethod(kTPConly),
  fCutZDCtiming(kFALSE),
  fTrigAna(),
  fCutImpactParameter(kFALSE),
  fImpactParameterMin(0.0),
  fImpactParameterMax(100.0),
  fhistTPCvsGlobalMult(0),
  fData2011(kFALSE),
  fCheckPileUp(kFALSE)
{
  //constructor 
}

//-----------------------------------------------------------------------
AliFlowEventCuts::AliFlowEventCuts(const char* name, const char* title):
  AliFlowEventSimpleCuts(name, title),
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
  fUtils(NULL),
  fMultSelection(NULL),
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
  fCutSPDTRKVtxZ(kFALSE),
  fCutTPCmultiplicityOutliers(kFALSE),
  fCutTPCmultiplicityOutliersAOD(kFALSE),
  fUseCentralityUnchecked(kFALSE),
  fCentralityPercentileMethod(kTPConly),
  fCutZDCtiming(kFALSE),
  fTrigAna(),
  fCutImpactParameter(kFALSE),
  fImpactParameterMin(0.0),
  fImpactParameterMax(100.0),
  fhistTPCvsGlobalMult(0),
  fData2011(kFALSE),
  fCheckPileUp(kFALSE)
{
  //constructor 
}

////-----------------------------------------------------------------------
AliFlowEventCuts::AliFlowEventCuts(const AliFlowEventCuts& that):
  AliFlowEventSimpleCuts(that),
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
  fUtils(NULL),
  fMultSelection(NULL),
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
  fCutSPDTRKVtxZ(that.fCutSPDTRKVtxZ),
  fCutTPCmultiplicityOutliers(that.fCutTPCmultiplicityOutliers),
  fCutTPCmultiplicityOutliersAOD(that.fCutTPCmultiplicityOutliersAOD),
  fUseCentralityUnchecked(that.fUseCentralityUnchecked),
  fCentralityPercentileMethod(that.fCentralityPercentileMethod),
  fCutZDCtiming(that.fCutZDCtiming),
  fTrigAna(),
  fCutImpactParameter(that.fCutImpactParameter),
  fImpactParameterMin(that.fImpactParameterMin),
  fImpactParameterMax(that.fImpactParameterMax),
  fhistTPCvsGlobalMult(that.fhistTPCvsGlobalMult),
  fData2011(that.fData2011),
  fCheckPileUp(that.fCheckPileUp)
{
  if (that.fQA) DefineHistograms();
  //copy constructor
  if (that.fRefMultCuts)
    fRefMultCuts = new AliFlowTrackCuts(*(that.fRefMultCuts));
  if (that.fMeanPtCuts)
    fMeanPtCuts = new AliFlowTrackCuts(*(that.fMeanPtCuts));
  fStandardTPCcuts = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();
  fStandardGlobalCuts = AliFlowTrackCuts::GetStandardGlobalTrackCuts2010();
    
  if (that.fUtils) {
      fUtils = new AliAnalysisUtils();
      fUtils->SetUseMVPlpSelection(kTRUE);
      fUtils->SetUseOutOfBunchPileUp(kTRUE);
  }
  if (that.fMultSelection) {
    fMultSelection = new AliMultSelection();
  }
}

////-----------------------------------------------------------------------
AliFlowEventCuts::~AliFlowEventCuts()
{
  //dtor
  delete fMeanPtCuts;
  delete fRefMultCuts;
  delete fStandardGlobalCuts;
  delete fStandardTPCcuts;
  if(fUtils) {
      delete fUtils;
      fUtils = NULL;
  }
  if (fMultSelection) delete fMultSelection;
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
  if (that.fUtils) {
      fUtils = new AliAnalysisUtils();
      fUtils->SetUseMVPlpSelection(kTRUE);
      fUtils->SetUseOutOfBunchPileUp(kTRUE);
  }
  if (that.fMultSelection) {
    fMultSelection = new AliMultSelection();
  }
  
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
  fCutSPDTRKVtxZ=that.fCutSPDTRKVtxZ;
  fCutTPCmultiplicityOutliers=that.fCutTPCmultiplicityOutliers;
  fCutTPCmultiplicityOutliersAOD=that.fCutTPCmultiplicityOutliersAOD;
  fUseCentralityUnchecked=that.fUseCentralityUnchecked;
  fCentralityPercentileMethod=that.fCentralityPercentileMethod;
  fCutZDCtiming=that.fCutZDCtiming;
  fhistTPCvsGlobalMult=that.fhistTPCvsGlobalMult;
  fData2011=that.fData2011;
  fCheckPileUp=that.fCheckPileUp;
  return *this;
}

//----------------------------------------------------------------------- 
Bool_t AliFlowEventCuts::IsSelected(TObject* obj, TObject* objmc)
{
  //check cuts
  AliVEvent* vevent = dynamic_cast<AliVEvent*>(obj);
  AliMCEvent* mcevent = dynamic_cast<AliMCEvent*>(objmc);
  if (vevent) return PassesCuts(vevent,mcevent);;
  return kFALSE;  //when passed wrong type of object
}
//----------------------------------------------------------------------- 
Bool_t AliFlowEventCuts::PassesCuts(AliVEvent *event, AliMCEvent *mcevent)
{
  ///check if event passes cuts
  const AliVVertex* pvtx=event->GetPrimaryVertex();
    Double_t pvtxx = 0.;
    Double_t pvtxy = 0.;
    Double_t pvtxz = 0.;
    Int_t ncontrib = 0;
    
    if(pvtx){
        pvtxx = pvtx->GetX();
        pvtxy = pvtx->GetY();
        pvtxz = pvtx->GetZ();
        ncontrib = pvtx->GetNContributors();
    }
    
  Bool_t pass=kTRUE;
  AliESDEvent* esdevent = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent* aodevent = dynamic_cast<AliAODEvent*>(event);
    
    
  if(fCheckPileUp){
      // pile-up rejection
      if(!fUtils) {
          fUtils = new AliAnalysisUtils();
          fUtils->SetUseMVPlpSelection(kTRUE);
          fUtils->SetUseOutOfBunchPileUp(kTRUE);
      }
      if(fUtils->IsPileUpEvent(aodevent)) pass = kFALSE;
   }
    

  Int_t multTPC = 0;
  Int_t multGlobal = 0; 

  // to remove multiplicity outliers, an explicit cut on the correlation 
  // between global and tpc only tracks can be made by counting the two
  // sets. as counting is expensive, only count when qa is requested or cut is enabeled
  // the cut criteria are different for different data takign periods
  // and (of course) cut criteria, specific routines exist for 10h, 11h data
  // and esds (by calling AliFlowTrackCuts) or aods (evaluated here explicitely)
  if(esdevent && (fQA || fCutTPCmultiplicityOutliers)) 
  {
    //this is pretty slow as we check the event track by track twice
    //this cut will work for 2010 PbPb data and is dependent on
    //TPC and ITS reco efficiency (e.g. geometry, calibration etc)
    multTPC = fStandardTPCcuts->Count(event);
    multGlobal = fStandardGlobalCuts->Count(event);
    if(fCutTPCmultiplicityOutliers) {
      if (multTPC > ( 23+1.216*multGlobal)) pass = kFALSE;
      if (multTPC < (-20+1.087*multGlobal)) pass = kFALSE;
    }
  }
  if(aodevent && (fQA || fCutTPCmultiplicityOutliersAOD)) 
  {
    //similar (slow) cut for aod's. will work for both 2010 and 2010 pbpb data
    //but the user is responsible that this object is configured
    //correctly to select the dataset
    //FIXME data could dynamically be determined by this class via the
    //runnumber
    Int_t nTracks(aodevent->GetNumberOfTracks());
    for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) { 
      AliAODTrack* track = dynamic_cast<AliAODTrack*>(aodevent->GetTrack(iTracks));
      if(!track) continue;
      if (!track || track->Pt() < .2 || track->Pt() > 5.0 || TMath::Abs(track->Eta()) > .8 || track->GetTPCNcls() < 70 || !track->GetDetPid() || track->GetDetPid()->GetTPCsignal() < 10.0)  continue;  // general quality cut
      if (track->TestFilterBit(1) && track->Chi2perNDF() > 0.2) multTPC++;
      if (!track->TestFilterBit(16) || track->Chi2perNDF() < 0.1) continue;
      Double_t b[2] = {-99., -99.};
      Double_t bCov[3] = {-99., -99., -99.};
      AliAODTrack copy(*track);
      if (copy.PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov) && TMath::Abs(b[0]) < 0.3 && TMath::Abs(b[1]) < 0.3) multGlobal++;
    }
    if(fCutTPCmultiplicityOutliersAOD) 
    {
      if(event->GetRunNumber() < 209122) {
        if(!fData2011 && (multTPC < (-40.3+1.22*multGlobal) || multTPC > (32.1+1.59*multGlobal))) pass = kFALSE;
        else if(fData2011  && (multTPC < (-36.73 + 1.48*multGlobal) || multTPC > (62.87 + 1.78*multGlobal))) pass = kFALSE;
      } else {
        if(multTPC < (-15.5+1.16*multGlobal) || multTPC > (15.8+1.28*multGlobal)) pass = kFALSE;
      }
    }
  }

  if (fQA)
  {
    QAbefore(0)->Fill(pvtxz);
    QAbefore(1)->Fill(multGlobal,multTPC);
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
   if(!fUseNewCentralityFramework)
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
   } // if(!fUseNewCentralityFramework)
   else
   {
    fMultSelection = (AliMultSelection *) esdevent->FindListObject("MultSelection");
    Float_t lPercentile = fMultSelection->GetMultiplicityPercentile(CentrMethName(fCentralityPercentileMethod));
    if(!(fCentralityPercentileMin <= lPercentile && lPercentile < fCentralityPercentileMax))
    {
      pass=kFALSE;
    }
   } // else
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
  if((fCutRefMult&&mcevent)||(fCutRefMult&&esdevent))
  {
    //reference multiplicity still to be defined
    Double_t refMult = RefMult(event,mcevent);
    if (refMult < fRefMultMin || refMult >= fRefMultMax )
    {
      pass=kFALSE;
    }
  }

  // Handles AOD event
  if(aodevent) {
    if(fCutSPDTRKVtxZ) {
      Double_t tVtxZ = aodevent->GetPrimaryVertex()->GetZ();
      Double_t tSPDVtxZ = aodevent->GetPrimaryVertexSPD()->GetZ();
      if( TMath::Abs(tVtxZ-tSPDVtxZ) > 0.5 ) pass = kFALSE;
    }
    if(!fUseNewCentralityFramework) {
      AliCentrality* centr = ((AliVAODHeader*)aodevent->GetHeader())->GetCentralityP();
      if(fCutTPCmultiplicityOutliers || fCutTPCmultiplicityOutliersAOD){
        Double_t v0Centr  = centr->GetCentralityPercentile("V0M");
        Double_t trkCentr = centr->GetCentralityPercentile("TRK");
        if(TMath::Abs(v0Centr-trkCentr) > 5) pass = kFALSE;
      }
      if (fCutCentralityPercentile) {
        if (fUseCentralityUnchecked) {
          if (!centr->IsEventInCentralityClassUnchecked( fCentralityPercentileMin,
                                                        fCentralityPercentileMax,
                                                        CentrMethName(fCentralityPercentileMethod) )) {
            pass = kFALSE;
          }
        } else {
          if (!centr->IsEventInCentralityClass( fCentralityPercentileMin,
                                               fCentralityPercentileMax,
                                               CentrMethName(fCentralityPercentileMethod) )) {
            pass = kFALSE;
          }
        }
      }
    }
    else {
      if (fCutCentralityPercentile) {
    	fMultSelection = (AliMultSelection *) aodevent->FindListObject("MultSelection");
    	Float_t lPercentile = fMultSelection->GetMultiplicityPercentile(CentrMethName(fCentralityPercentileMethod));
    	if(!fMultSelection){
    	  AliWarning("AliMultSelection not found, did you Run AliMultSelectionTask? \n");
    	}
        if(!(fCentralityPercentileMin <= lPercentile && lPercentile < fCentralityPercentileMax))
        {
          pass=kFALSE;
        }
      }
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
    if (nselected) meanpt=meanpt/nselected;
    if (meanpt<fMeanPtMin || meanpt >= fMeanPtMax) pass=kFALSE;
  }

  //impact parameter cut
  if(fCutImpactParameter) {
    Double_t gImpactParameter = 0.;
    if(mcevent) {
      AliCollisionGeometry* headerH = dynamic_cast<AliCollisionGeometry*>(dynamic_cast<AliMCEvent*>(mcevent)->GenEventHeader());
      if(headerH)
	gImpactParameter = headerH->ImpactParameter();
    }
    if ((gImpactParameter < fImpactParameterMin) || (gImpactParameter >= fImpactParameterMax ))
      pass=kFALSE;
  }

  if (fQA&&pass) 
  {
    QAafter(1)->Fill(multGlobal,multTPC);
    QAafter(0)->Fill(pvtxz);
  }
  return pass;
}

//----------------------------------------------------------------------- 
Float_t AliFlowEventCuts::GetCentrality(AliVEvent* event, AliMCEvent* /*mcEvent*/)
{
  //get the centrality percentile of the event
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(event);
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(event);

  Float_t centrality=-1.;

  if(!fUseNewCentralityFramework)
  {
  AliCentrality* centr = NULL;
  if (esdEvent)
    centr = esdEvent->GetCentrality();
  if (aodEvent) 
    centr = ((AliVAODHeader*)aodEvent->GetHeader())->GetCentralityP();
  
  if (!centr) return -1.;

  if (fUseCentralityUnchecked) 
    centrality=centr->GetCentralityPercentileUnchecked(CentrMethName(fCentralityPercentileMethod));
  else 
    centrality=centr->GetCentralityPercentile(CentrMethName(fCentralityPercentileMethod));
  }
  else // if(!fUseNewCentralityFramework)
  {
   if (esdEvent)
   {
    fMultSelection = (AliMultSelection *) esdEvent->FindListObject("MultSelection");
    Float_t lPercentile = fMultSelection->GetMultiplicityPercentile(CentrMethName(fCentralityPercentileMethod));
    centrality = lPercentile;
   }
   else if (aodEvent)
   {
     fMultSelection = (AliMultSelection *) aodEvent->FindListObject("MultSelection");
     Float_t lPercentile = fMultSelection->GetMultiplicityPercentile(CentrMethName(fCentralityPercentileMethod));
     if(!fMultSelection){
       AliWarning("AliMultSelection not found, did you Run AliMultSelectionTask? \n");
     }
     centrality = lPercentile;
   } // else if (aodEvent)
  } // else // if(!fUseNewCentralityFramework)

  return centrality;
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
    case kVZERO:
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
Int_t AliFlowEventCuts::RefMult(AliVEvent* event, AliMCEvent *mcEvent)
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
  else if (fRefMultMethod==kVZERO)
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
  fRefMultCuts->SetEvent(event,mcEvent);
  Int_t numberOfInputObjects = fRefMultCuts->GetNumberOfInputObjects();
  for (Int_t i=0; i<numberOfInputObjects; i++) {
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
  before->Add(new TH2F("fTPCvsGlobalMult","TPC only vs Global track multiplicity;global;TPC only",600,0,3000,600,0,4000));//1
  after->Add(new TH2F("fTPCvsGlobalMult","TPC only vs Global track multiplicity;global;TPC only",600,0,3000,600,0,4000));//1
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


