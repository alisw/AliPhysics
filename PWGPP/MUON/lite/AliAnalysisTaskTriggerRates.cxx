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

#include <TString.h>

#include "AliLog.h"
#include "AliCounterCollection.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliVParticle.h"
#include "AliESDMuonTrack.h"
#include "AliInputEventHandler.h"
#include "AliCentrality.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisMuonUtility.h"

#include "AliAnalysisTaskTriggerRates.h"

ClassImp(AliAnalysisTaskTriggerRates)

//________________________________________________________________________
AliAnalysisTaskTriggerRates::AliAnalysisTaskTriggerRates()
: AliAnalysisTaskSE(),
fTriggerCounters(0),
fnCent(0),
fDeltaDev(0),
fMULPattern("none"),
fMLLPattern("none"),
fMSLPattern("none"),
fMSHPattern("none"),
fPrinfCounts(kFALSE)
{
  // Dummy constructor ALWAYS needed for I/O.
  
  InitCentralityBins();
  
}

//________________________________________________________________________
AliAnalysisTaskTriggerRates::AliAnalysisTaskTriggerRates(const char *name)
: AliAnalysisTaskSE(name),
fTriggerCounters(0),
fnCent(0),
fDeltaDev(0),
fMULPattern("none"),
fMLLPattern("none"),
fMSLPattern("none"),
fMSHPattern("none"),
fPrinfCounts(kFALSE)
{
  // Constructor
  
  InitCentralityBins();
  
  DefineOutput(1, AliCounterCollection::Class());
}

//________________________________________________________________________
AliAnalysisTaskTriggerRates::~AliAnalysisTaskTriggerRates()
{
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fTriggerCounters;
}

//________________________________________________________________________
void AliAnalysisTaskTriggerRates::UserCreateOutputObjects()
{
  // Create outputs
  
  TString centKeys = "any";
  for (Int_t iCent = 0; iCent < fnCent; ++iCent) centKeys += Form("/%s", fCentBinName[iCent].Data());
  
  fTriggerCounters = new AliCounterCollection(GetOutputSlot(1)->GetContainer()->GetName());
  fTriggerCounters->AddRubric("run", 1000);
  fTriggerCounters->AddRubric("ps", "yes/no");
  fTriggerCounters->AddRubric("cent", centKeys.Data());
  fTriggerCounters->AddRubric("mul", "yes/no");
  fTriggerCounters->AddRubric("mll", "yes/no");
  fTriggerCounters->AddRubric("msl", "yes/no");
  fTriggerCounters->AddRubric("msh", "yes/no");
  fTriggerCounters->AddRubric("muloff", "yes/no");
  fTriggerCounters->AddRubric("mlloff", "yes/no");
  fTriggerCounters->AddRubric("msloff", "yes/no");
  fTriggerCounters->AddRubric("mshoff", "yes/no");
  fTriggerCounters->Init();
  
  PostData(1, fTriggerCounters);
}

//________________________________________________________________________
void AliAnalysisTaskTriggerRates::UserExec(Option_t *)
{
  // One event processing
  
  AliVEvent *evt = InputEvent();
  if (!dynamic_cast<AliESDEvent*>(evt)) {
    AliFatal("can only run on ESDs");
    return;
  }
  
  // select physics events
  if (evt->GetEventType() != 7) return;
  
  // Get Physics selection
  TString psKey = (fInputHandler->IsEventSelected()) ? "ps:yes" : "ps:no";
  
  // Get Centrality
  Float_t cent = evt->GetCentrality()->GetCentralityPercentile("V0M");
  TString centKey = "";
  for (Int_t iCent = 0; iCent < fnCent; iCent++) {
    if (cent > fCentBinRange[iCent][0] && cent <= fCentBinRange[iCent][1]) {
      centKey = Form("cent:%s",fCentBinName[iCent].Data());
      break;
    }
  }
  
  // Get trigger class
  TString trigger = evt->GetFiredTriggerClasses();
  TString onlineMULKey = (trigger.Contains(fMULPattern.Data())) ? "mul:yes" : "mul:no";
  TString onlineMLLKey = (trigger.Contains(fMLLPattern.Data())) ? "mll:yes" : "mll:no";
  TString onlineMSLKey = (trigger.Contains(fMSLPattern.Data())) ? "msl:yes" : "msl:no";
  TString onlineMSHKey = (trigger.Contains(fMSHPattern.Data())) ? "msh:yes" : "msh:no";
  
  // Get trigger class from trigger tracks
  TString offlineMULKey = "muloff:no";
  TString offlineMLLKey = "mlloff:no";
  TString offlineMSLKey = "msloff:no";
  TString offlineMSHKey = "mshoff:no";
  Int_t ntracks = AliAnalysisMuonUtility::GetNTracks(evt);
  Int_t *loCircuit = new Int_t[ntracks];
  Int_t nTrgTracks = 0;
  for (Int_t i = 0; i < ntracks; ++i) {
    
    AliVParticle *tracki = AliAnalysisMuonUtility::GetTrack(i, evt);
    if (!AliAnalysisMuonUtility::MatchLpt(tracki)) continue;
    
    // make sure this track is not already accounted for
    Bool_t trackExist = kFALSE;
    Int_t loCircuiti = AliAnalysisMuonUtility::GetLoCircuit(tracki);
    for (Int_t k = 0; k < nTrgTracks; ++k) {
      if (loCircuiti == loCircuit[k]) {
        trackExist = kTRUE;
        break;
      }
    }
    if (trackExist) continue;
    loCircuit[nTrgTracks++] = loCircuiti;
    
    // fire single muon trigger(s)
    offlineMSLKey = "msloff:yes";
    if (AliAnalysisMuonUtility::MatchHpt(tracki)) offlineMSHKey = "mshoff:yes";
    
    Int_t trgDevSigni = TriggerDevSign(tracki, fDeltaDev);
    
    for (Int_t j = i+1; j < ntracks; ++j) {
      
      AliVParticle *trackj = AliAnalysisMuonUtility::GetTrack(j, evt);
      if (!AliAnalysisMuonUtility::MatchLpt(trackj)) continue;
      
      // make sure this track is not already accounted for
      trackExist = kFALSE;
      Int_t loCircuitj = AliAnalysisMuonUtility::GetLoCircuit(trackj);
      for (Int_t k = 0; k < nTrgTracks; ++k) {
        if (loCircuitj == loCircuit[k]) {
          trackExist = kTRUE;
          break;
        }
      }
      if (trackExist) continue;
      
      Int_t trgDevSignj = TriggerDevSign(trackj, fDeltaDev);
      
      // fire dimuon trigger(s)
      Int_t trgSign = trgDevSigni * trgDevSignj;
      if (trgSign <= 0) offlineMULKey = "muloff:yes";
      if (trgSign >= 0) offlineMLLKey = "mlloff:yes";
      
    }
  }
  delete[] loCircuit;
  
  // Count
  fTriggerCounters->Count(Form("run:%d/%s/cent:any/%s/%s/%s/%s/%s/%s/%s/%s",
                               fCurrentRunNumber, psKey.Data(),
                               onlineMULKey.Data(), onlineMLLKey.Data(), onlineMSLKey.Data(), onlineMSHKey.Data(),
                               offlineMULKey.Data(), offlineMLLKey.Data(), offlineMSLKey.Data(), offlineMSHKey.Data()));
  if (!centKey.IsNull())
    fTriggerCounters->Count(Form("run:%d/%s/%s/%s/%s/%s/%s/%s/%s/%s/%s",
                                 fCurrentRunNumber, psKey.Data(), centKey.Data(),
                                 onlineMULKey.Data(), onlineMLLKey.Data(), onlineMSLKey.Data(), onlineMSHKey.Data(),
                                 offlineMULKey.Data(), offlineMLLKey.Data(), offlineMSLKey.Data(), offlineMSHKey.Data()));
  
  PostData(1, fTriggerCounters);
}

//_____________________________________________________________________________
void AliAnalysisTaskTriggerRates::Terminate(Option_t *)
{
  // Post-processing
  
  fTriggerCounters = dynamic_cast<AliCounterCollection*>(GetOutputData(1));
  if (!fTriggerCounters) return;
  
  printf("\nwithout Physics Selection:\n\n");
  PrintRates("","");
  
  printf("\nwith Physics Selection:\n\n");
  PrintRates("ps:yes","");
  
  printf("\nwith Physics Selection in 0-90%%:\n\n");
  PrintRates("ps:yes","cent:010,1020,2030,3040,4050,5060,6070,7080,8090");
  
}

//_____________________________________________________________________________
void AliAnalysisTaskTriggerRates::InitCentralityBins()
{
  /// Set the default centrality bins
  
  for (Int_t iCent = 0; iCent < 10; ++iCent) {
    fCentBinRange[iCent][0] = 10.*iCent;
    fCentBinRange[iCent][1] = 10.*(iCent+1);
    fCentBinName[iCent] = Form("%d%d", 10*iCent, 10*(iCent+1));
  }
  
}

//_____________________________________________________________________________
Int_t AliAnalysisTaskTriggerRates::TriggerDevSign(AliVParticle *track, UInt_t deltaDev) const
{
  /// get the sign (±1) of track deviation in the trigger (0 = unknown)
  
  AliESDMuonTrack *esdTrack = dynamic_cast<AliESDMuonTrack*>(track);
  if (!esdTrack) return 0;
  
  Int_t deviation = esdTrack->LoDev();
  if (deviation > 15+((Int_t)deltaDev)) return 1;
  else if (deviation < 15-((Int_t)deltaDev)) return -1;
  else return 0;
  
}

//_____________________________________________________________________________
void AliAnalysisTaskTriggerRates::PrintRates(TString ps, TString cent) const
{
  /// print trigger online versus trigger offline
  
  if (!ps.IsNull()) ps += "/";
  if (!cent.IsNull()) cent += "/";
  
  const char* sformat = "%12s";
  const char* dformat = "%12ld";
  const char* fformat = "%11.1f%%";
  
  // define online trigger selection (format = selectKey!rejectKey)
  const Int_t ni = 8;
  TString titlei[ni] = {"MUL", "MLL", "MUL|MLL", "MSL", "MSLonly", "MSH", "MUL|MSH", "MUL|MLL|MSH"};
  TString selecti[ni] = {"mul:yes", "mll:yes", "!mul:no/mll:no", "msl:yes", "msl:yes/mul:no/mll:no/msh:no", "msh:yes", "!mul:no/msh:no", "!mul:no/mll:no/msh:no"};
  
  // define offline trigger selection (format = selectKey!rejectKey)
  const Int_t nj = 12;
  TString titlej[nj] = {"any", "MULonly", "MLLonly", "MUL&MLL", "MUL", "MLL", "MUL|MLL", "MSL", "MSLonly", "MSH", "MUL|MSH", "MUL|MLL|MSH"};
  TString selectj[nj] = {"", "muloff:yes/mlloff:no", "muloff:no/mlloff:yes", "muloff:yes/mlloff:yes", "muloff:yes", "mlloff:yes", "!muloff:no/mlloff:no", "msloff:yes", "msloff:yes/muloff:no/mlloff:no/mshoff:no", "mshoff:yes", "!muloff:no/mshoff:no", "!muloff:no/mlloff:no/mshoff:no"};
  
  Long_t valij[ni][nj], rvali[ni][nj], rvalj[ni][nj], rvalij[ni][nj], val[ni][nj];
  
  printf(sformat,"");
  for (Int_t j = 0; j < nj; ++j) printf(sformat,titlej[j].Data());
  printf("\n");
  for (Int_t i = 0; i < ni; ++i) {
    
    printf(sformat, titlei[i].Data());
    
    // build online selectKey and selectKey/rejectKey
    TString si = selecti[i];
    TString ri = "";
    if (si.Contains("!")) {
      ri = si;
      ri.ReplaceAll("!","/");
      ri.Remove(TString::kLeading,'/');
      si.Remove(si.Index("!"));
    }
    
    for (Int_t j = 0; j < nj; ++j) {
      
      // build offline selectKey and selectKey/rejectKey
      TString sj = selectj[j];
      TString rj = "";
      if (sj.Contains("!")) {
        rj = sj;
        rj.ReplaceAll("!","/");
        rj.Remove(TString::kLeading,'/');
        sj.Remove(sj.Index("!"));
      }
      
      // val = selectKeyi/selectKeyj/(any-rejectKeyi)/(any-rejectKeyj)
      valij[i][j] = (Long_t)fTriggerCounters->GetSum(Form("%s%s%s/%s",ps.Data(),cent.Data(),si.Data(),sj.Data()));
      rvali[i][j] = (!ri.IsNull()) ? (Long_t)fTriggerCounters->GetSum(Form("%s%s%s/%s",ps.Data(),cent.Data(),ri.Data(),sj.Data())) : 0;
      rvalj[i][j] = (!rj.IsNull()) ? (Long_t)fTriggerCounters->GetSum(Form("%s%s%s/%s",ps.Data(),cent.Data(),si.Data(),rj.Data())) : 0;
      rvalij[i][j] = (!ri.IsNull() && !rj.IsNull()) ? (Long_t)fTriggerCounters->GetSum(Form("%s%s%s/%s",ps.Data(),cent.Data(),ri.Data(),rj.Data())) : 0;
      val[i][j] = valij[i][j]-rvali[i][j]-rvalj[i][j]+rvalij[i][j];
      
      if (fPrinfCounts) printf(dformat, val[i][j]);
      else printf(fformat, ((Float_t)val[i][j])/((Float_t)val[i][0])*100.);
      
    }
    printf("\n");
    
  }
  printf("\n");
  
}
