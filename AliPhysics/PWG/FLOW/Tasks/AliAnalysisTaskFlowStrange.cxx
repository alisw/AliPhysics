/*************************************************************************
* Copyright(c) 1998-2008,ALICE Experiment at CERN,All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use,copy,modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee,provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/////////////////////////////////////////////////////
// AliAnalysisTaskFlowStrange:
// Analysis task to select K0/Lambda candidates for flow analysis.
// Authors: Cristian Ivan (civan@cern.ch)
//          Carlos Perez  (cperez@cern.ch)
//          Pawel Debski  (pdebski@cern.ch)
//////////////////////////////////////////////////////

#include "TChain.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TFile.h"

#include "TRandom3.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliVVertex.h"
#include "AliVVZERO.h"
#include "AliStack.h"
#include "AliMCEvent.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"
#include "AliAODTracklets.h"
#include "AliAODHeader.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "TMath.h"
#include "TObjArray.h"
#include "AliFlowCandidateTrack.h"

#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliFlowEvent.h"
#include "AliFlowBayesianPID.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowVector.h"

#include "AliAnalysisTaskFlowStrange.h"

ClassImp(AliAnalysisTaskFlowStrange)

//=======================================================================
AliAnalysisTaskFlowStrange::AliAnalysisTaskFlowStrange() :
  AliAnalysisTaskSE(),
  fPIDResponse(NULL),
  fFB1(NULL),
  fFB1024(NULL),
  fTPCevent(NULL),
  fVZEevent(NULL),
  fCandidates(NULL),
  fList(NULL),
  fRunNumber(-1),
  fDebug(0),
  fQAlevel(0),
  fReadESD(kFALSE),
  fReadMC(kFALSE),
  fAddPiToMCReactionPlane(kTRUE),
  fPostMatched(0),
  fAvoidExec(kFALSE),
  fSkipSelection(kFALSE),
  fSkipVn(kFALSE),
  fUseFP(kFALSE),
  fRunOnpA(kFALSE),
  fRunOnpp(kFALSE),
  fExtraEventRejection(kFALSE),
  fSkipCentralitySelection(kFALSE),
  fCentMethod("V0MTRK"),
  fCentPerMin(0),
  fCentPerMax(100),
  fThisCent(-1.0),
  fV0M(0.0),
  fTRK(0.0),
  fPriVtxZ(0.0),
  fSPDVtxZ(0.0),
  fSPDtracklets(0),
  fVZETotM(0.0),
  fRefMultTPC(0),
  fRefMultHyb(0),
  fVertexZcut(10.0),
  fExcludeTPCEdges(kFALSE),
  fSpecie(0),
  fOnline(kFALSE),
  fHomemade(kFALSE),
  fWhichPsi(1),
  fVZEsave(kFALSE),
  fVZEload(NULL),
  fVZEResponse(NULL),
  fVZEmb(kFALSE),
  fVZEByDisk(kTRUE),
  fVZECa(0),
  fVZECb(3),
  fVZEAa(0),
  fVZEAb(3),
  fVZEQA(NULL),
  fHarmonic(2),
  fPsi2(0.0),
  fMCEP(0.0),
  fQVZEACos(0.0),
  fQVZEASin(0.0),
  fQVZECCos(0.0),
  fQVZECSin(0.0),
  fQVZEA(0.0),
  fQVZEC(0.0),
  fVZEWarning(kFALSE),
  fQTPCACos(0.0),
  fQTPCASin(0.0),
  fQTPCCCos(0.0),
  fQTPCCSin(0.0),
  fQTPC2hCos(0.0),
  fQTPC2hSin(0.0),
  fQTPCA(0.0),
  fQTPCC(0.0),
  fQTPCA_nTracks(0),
  fQTPCC_nTracks(0),
  fSkipTerminate(kTRUE),
  fMassBins(0),
  fMinMass(0.0),
  fMaxMass(0.0),
  fPtBins(0),
  fRFPFilterBit(1),
  fRFPminPt(0.2),
  fRFPmaxPt(5.0),
  fRFPAminEta(0.0),
  fRFPAmaxEta(+0.8),
  fRFPCminEta(-0.8),
  fRFPCmaxEta(0.0),
  fRFPTPCsignal(10.0),
  fRFPmaxIPxy(2.4),
  fRFPmaxIPz(3.2),
  fRFPTPCncls(70),
  fDecayMass(0.0),
  fDecayPhi(0.0),
  fDecayEta(0.0),
  fDecayPt(0.0),
  fDecayDCAdaughters(0.0),
  fDecayCosinePointingAngleXY(0.0),
  fDecayRadXY(0.0),
  fDecayDecayLength(0.0),
  fDecayDecayLengthLab(0.0),
  fDecayQt(0.0),
  fDecayAlpha(0.0),
  fDecayRapidity(0.0),
  fDecayProductIPXY(0.0),
  fDecayIPneg(0.0),
  fDecayIPpos(0.0),
  fDecayXneg(0.0),
  fDecayXpos(0.0),
  fDecayIDneg(-1),
  fDecayIDpos(-1),
  fDecayID(-1),
  fDecayMatchOrigin(0.0),
  fDecayMatchPhi(0.0),
  fDecayMatchEta(0.0),
  fDecayMatchPt(0.0),
  fDecayMatchRadXY(0.0),
  fDecayMinEta(0.0),
  fDecayMaxEta(0.0),
  fDecayMinPt(0.0),
  fDecayMaxDCAdaughters(0.0),
  fDecayMinCosinePointingAngleXY(0.0),
  fDecayMinQt(0.0),
  fDecayAPCutPie(kTRUE),
  fDecayStopPIDAtPt(3.0),
  fDecayMinRadXY(0.0),
  fDecayMaxDecayLength(0.0),
  fDecayMaxProductIPXY(0.0),
  fDecayMaxRapidity(0.0),
  fDaughterPhi(0.0),
  fDaughterEta(0.0),
  fDaughterPt(0.0),
  fDaughterNClsTPC(0),
  fDaughterNClsITS(0),
  fDaughterCharge(0),
  fDaughterNFClsTPC(0),
  fDaughterNSClsTPC(0),
  fDaughterChi2PerNClsTPC(0.0),
  fDaughterXRows(0.0),
  fDaughterImpactParameterXY(0.0),
  fDaughterImpactParameterZ(0.0),
  fDaughterStatus(0),
  fDaughterITScm(0),
  fDaughterNSigmaPID(0.0),
  fDaughterKinkIndex(0),
  fDaughterAtSecPhi(0.0),
  fDaughterAtSecEta(0.0),
  fDaughterAtSecPt(0.0),
  fDaughterMatchPhi(0.0),
  fDaughterMatchEta(0.0),
  fDaughterMatchPt(0.0),
  fDaughterMatchImpactParameterXY(0.0),
  fDaughterMatchImpactParameterZ(0.0),
  fDaughterUnTag(kTRUE),
  fDaughterMinEta(0.0),
  fDaughterMaxEta(0.0),
  fDaughterMinPt(0.0),
  fDaughterMinNClsTPC(0),
  fDaughterMinNClsITS(-1),
  fDaughterMinXRows(0),
  fDaughterMaxChi2PerNClsTPC(0.0),
  fDaughterMinXRowsOverNClsFTPC(0.0),
  fDaughterMinImpactParameterXY(0.0),
  fDaughterMaxNSigmaPID(0.0),
  fDaughterSPDRequireAny(kFALSE),
  fDaughterITSrefit(kFALSE) {
  //ctor
  for(Int_t i=0; i!=100; ++i) fPtBinEdge[i]=0;
  for(Int_t i=0; i!=6; ++i) fDaughterITSConfig[i]=-1;
  for(Int_t i=0; i!=2000; ++i) fQTPCA_fID[i]=-1;
  for(Int_t i=0; i!=2000; ++i) fQTPCC_fID[i]=-1;
  for(Int_t i=0; i!=64; ++i) fVZEextW[i]=1;
}
//=======================================================================
AliAnalysisTaskFlowStrange::AliAnalysisTaskFlowStrange(const char *name) :
  AliAnalysisTaskSE(name),
  fPIDResponse(NULL),
  fFB1(NULL),
  fFB1024(NULL),
  fTPCevent(NULL),
  fVZEevent(NULL),
  fCandidates(NULL),
  fList(NULL),
  fRunNumber(-1),
  fDebug(0),
  fQAlevel(0),
  fReadESD(kFALSE),
  fReadMC(kFALSE),
  fAddPiToMCReactionPlane(kTRUE),
  fPostMatched(0),
  fAvoidExec(kFALSE),
  fSkipSelection(kFALSE),
  fSkipVn(kFALSE),
  fUseFP(kFALSE),
  fRunOnpA(kFALSE),
  fRunOnpp(kFALSE),
  fExtraEventRejection(kFALSE),
  fSkipCentralitySelection(kFALSE),
  fCentMethod("V0MTRK"),
  fCentPerMin(0),
  fCentPerMax(100),
  fThisCent(-1.0),
  fV0M(0.0),
  fTRK(0.0),
  fPriVtxZ(0.0),
  fSPDVtxZ(0.0),
  fSPDtracklets(0),
  fVZETotM(0.0),
  fRefMultTPC(0),
  fRefMultHyb(0),
  fVertexZcut(10.0),
  fExcludeTPCEdges(kFALSE),
  fSpecie(0),
  fOnline(kFALSE),
  fHomemade(kFALSE),
  fWhichPsi(1),
  fVZEsave(kFALSE),
  fVZEload(NULL),
  fVZEResponse(NULL),
  fVZEmb(kFALSE),
  fVZEByDisk(kTRUE),
  fVZECa(0),
  fVZECb(3),
  fVZEAa(0),
  fVZEAb(3),
  fVZEQA(NULL),
  fHarmonic(2),
  fPsi2(0.0),
  fMCEP(0.0),
  fQVZEACos(0.0),
  fQVZEASin(0.0),
  fQVZECCos(0.0),
  fQVZECSin(0.0),
  fQVZEA(0.0),
  fQVZEC(0.0),
  fVZEWarning(kFALSE),
  fQTPCACos(0.0),
  fQTPCASin(0.0),
  fQTPCCCos(0.0),
  fQTPCCSin(0.0),
  fQTPC2hCos(0.0),
  fQTPC2hSin(0.0),
  fQTPCA(0.0),
  fQTPCC(0.0),
  fQTPCA_nTracks(0),
  fQTPCC_nTracks(0),
  fSkipTerminate(kTRUE),
  fMassBins(0),
  fMinMass(0.0),
  fMaxMass(0.0),
  fPtBins(0),
  fRFPFilterBit(1),
  fRFPminPt(0.2),
  fRFPmaxPt(5.0),
  fRFPAminEta(0.0),
  fRFPAmaxEta(+0.8),
  fRFPCminEta(-0.8),
  fRFPCmaxEta(0.0),
  fRFPTPCsignal(10.0),
  fRFPmaxIPxy(2.4),
  fRFPmaxIPz(3.2),
  fRFPTPCncls(70),
  fDecayMass(0.0),
  fDecayPhi(0.0),
  fDecayEta(0.0),
  fDecayPt(0.0),
  fDecayDCAdaughters(0.0),
  fDecayCosinePointingAngleXY(0.0),
  fDecayRadXY(0.0),
  fDecayDecayLength(0.0),
  fDecayDecayLengthLab(0.0),
  fDecayQt(0.0),
  fDecayAlpha(0.0),
  fDecayRapidity(0.0),
  fDecayProductIPXY(0.0),
  fDecayIPneg(0.0),
  fDecayIPpos(0.0),
  fDecayXneg(0.0),
  fDecayXpos(0.0),
  fDecayIDneg(-1),
  fDecayIDpos(-1),
  fDecayID(-1),
  fDecayMatchOrigin(0.0),
  fDecayMatchPhi(0.0),
  fDecayMatchEta(0.0),
  fDecayMatchPt(0.0),
  fDecayMatchRadXY(0.0),
  fDecayMinEta(0.0),
  fDecayMaxEta(0.0),
  fDecayMinPt(0.0),
  fDecayMaxDCAdaughters(0.0),
  fDecayMinCosinePointingAngleXY(0.0),
  fDecayMinQt(0.0),
  fDecayAPCutPie(kTRUE),
  fDecayStopPIDAtPt(3.0),
  fDecayMinRadXY(0.0),
  fDecayMaxDecayLength(0.0),
  fDecayMaxProductIPXY(0.0),
  fDecayMaxRapidity(0.0),
  fDaughterPhi(0.0),
  fDaughterEta(0.0),
  fDaughterPt(0.0),
  fDaughterNClsTPC(0),
  fDaughterNClsITS(0),
  fDaughterCharge(0),
  fDaughterNFClsTPC(0),
  fDaughterNSClsTPC(0),
  fDaughterChi2PerNClsTPC(0.0),
  fDaughterXRows(0.0),
  fDaughterImpactParameterXY(0.0),
  fDaughterImpactParameterZ(0.0),
  fDaughterStatus(0),
  fDaughterITScm(0),
  fDaughterNSigmaPID(0.0),
  fDaughterKinkIndex(0),
  fDaughterAtSecPhi(0.0),
  fDaughterAtSecEta(0.0),
  fDaughterAtSecPt(0.0),
  fDaughterMatchPhi(0.0),
  fDaughterMatchEta(0.0),
  fDaughterMatchPt(0.0),
  fDaughterMatchImpactParameterXY(0.0),
  fDaughterMatchImpactParameterZ(0.0),
  fDaughterUnTag(kTRUE),
  fDaughterMinEta(0.0),
  fDaughterMaxEta(0.0),
  fDaughterMinPt(0.0),
  fDaughterMinNClsTPC(0),
  fDaughterMinNClsITS(-1),
  fDaughterMinXRows(0),
  fDaughterMaxChi2PerNClsTPC(0.0),
  fDaughterMinXRowsOverNClsFTPC(0.0),
  fDaughterMinImpactParameterXY(0.0),
  fDaughterMaxNSigmaPID(0.0),
  fDaughterSPDRequireAny(kFALSE),
  fDaughterITSrefit(kFALSE) {
  //ctor
  for(Int_t i=0; i!=100; ++i) fPtBinEdge[i]=0;
  for(Int_t i=0; i!=6; ++i) fDaughterITSConfig[i]=-1;
  for(Int_t i=0; i!=2000; ++i) fQTPCA_fID[i]=-1;
  for(Int_t i=0; i!=2000; ++i) fQTPCC_fID[i]=-1;
  for(Int_t i=0; i!=64; ++i) fVZEextW[i]=1;
  DefineInput( 0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,AliFlowEventSimple::Class()); // TPC object
  DefineOutput(3,AliFlowEventSimple::Class()); // VZE object
}
//=======================================================================
AliAnalysisTaskFlowStrange::~AliAnalysisTaskFlowStrange() {
  //dtor
  if (fCandidates) delete fCandidates;
  if (fTPCevent)   delete fTPCevent;
  if (fVZEevent)   delete fVZEevent;
  if (fList)       delete fList;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::SetPtEdges(Int_t n, Double_t *p) {
  fPtBins = n;
  for(int i=0;i!=n+1;++i) fPtBinEdge[i] = p[i];
}
//=======================================================================
TList* AliAnalysisTaskFlowStrange::RunTerminateAgain(TList *lst) {
  if(!lst) return NULL;
  fList = lst;
  fSpecie = Int_t( ((TProfile*)((TList*)fList->FindObject("Event"))->FindObject("Configuration"))->GetBinContent(kSpecie) );
  fSkipSelection = ((TProfile*)((TList*)fList->FindObject("Event"))->FindObject("Configuration"))->GetBinContent(kSkipSelection);
  fReadMC = ((TProfile*)((TList*)fList->FindObject("Event"))->FindObject("Configuration"))->GetBinContent(kReadMC);
  Terminate(NULL);
  return fList;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::PrintConfig() {
  //DUMP for main task
  printf("******************************\n");
  printf("<TASK Configuration> %s\n",GetName());
  printf("  fDebug %d\n",fDebug);
  printf("  fQAlevel %d\n",fQAlevel);
  printf("  fExtraEventRejection %s\n",fExtraEventRejection?"kTRUE":"kFALSE");
  printf("  fCentMethod %s\n",fCentMethod.Data());
  printf("    fCentPerMin %d\n",fCentPerMin);
  printf("    fCentPerMax %d\n",fCentPerMax);
  printf("  fVextexZcut %f\n",fVertexZcut);
  printf("  fRunOnpA %s\n",fRunOnpA?"kTRUE":"kFALSE");
  printf("  fRunOnpp %s\n",fRunOnpp?"kTRUE":"kFALSE");
  printf("  fReadESD %s\n",fReadESD?"kTRUE":"kFALSE");
  printf("  fReadMC %s\n",fReadMC?"kTRUE":"kFALSE");
  if(fReadMC) {
    printf("    fAddPiToMCReactionPlane %s\n",fAddPiToMCReactionPlane?"kTRUE":"kFALSE");
    printf("    fPostMatched %d\n",fPostMatched);
    printf("    fAvoidExec %s\n",fAvoidExec?"kTRUE":"kFALSE");
    printf("    fSkipCentralitySelection %s\n",fSkipCentralitySelection?"kTRUE":"kFALSE");
  }
  printf("  fVZEsave %s\n",fVZEsave?"kTRUE":"kFALSE");
  if(fVZEload) {
    printf("  fVZEload %d runs\n",fVZEload->GetEntries());
    printf("    fVZEmb %s\n",fVZEmb?"kTRUE":"kFALSE");
    printf("    fVZEByDisk %s\n",fVZEByDisk?"kTRUE":"kFALSE");
  }
  printf("  fHarmonic %d\n",fHarmonic);
  printf("    fWhichPsi %d\n",fWhichPsi);
  printf("    fVZECa %d\n",fVZECa);
  printf("    fVZECb %d\n",fVZECb);
  printf("    fVZEAa %d\n",fVZEAa);
  printf("    fVZEAb %d\n",fVZEAb);
  printf("    fRFPFilterBit %d\n",fRFPFilterBit);
  printf("    fRFPminPt %f\n",fRFPminPt);
  printf("    fRFPmaxPt %f\n",fRFPmaxPt);
  printf("    fRFPAminEta %f\n",fRFPAminEta);
  printf("    fRFPAmaxEta %f\n",fRFPAmaxEta);
  printf("    fRFPCminEta %f\n",fRFPCminEta);
  printf("    fRFPCmaxEta %f\n",fRFPCmaxEta);
  printf("    fRFPmaxIPxy %f\n",fRFPmaxIPxy);
  printf("    fRFPmaxIPz %f\n",fRFPmaxIPz);
  printf("    fRFPTPCsignal %f\n",fRFPTPCsignal);
  printf("    fRFPTPCncls %d\n",fRFPTPCncls);
  printf("    fExcludeTPCEdges %s\n",fExcludeTPCEdges?"kTRUE":"kFALSE");
  printf("  fSkipSelection %s\n",fSkipSelection?"kTRUE":"kFALSE");
  if(!fSkipSelection) {
    printf("    fSpecie %d\n",fSpecie);
    printf("    fPtBins %d\n      |",fPtBins);
    for(int i=0; i!=fPtBins+1; ++i) printf("%f|",fPtBinEdge[i]); printf("\n");
    if(fSpecie<90) {
      printf("    fMassBins %d\n",fMassBins);
      printf("    fMinMass %f\n",fMinMass);
      printf("    fMaxMass %f\n",fMaxMass);
    }
  }
  printf("  fSkipVn %s\n",fSkipVn?"kTRUE":"kFALSE");
  if(!fSkipVn) {
    printf("    fUseFP %s\n",fUseFP?"kTRUE":"kFALSE");
  }
  MyPrintConfig();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MyPrintConfig() {
  // Dump for derived task
  printf("==================================\n");
  printf("<FlowStrange> \n");
  if(!fSkipSelection) {
    if(fReadESD) {
      printf("  fOnline %s\n",fOnline?"kTRUE":"kFALSE");
      printf("  fHomemade %s\n",fHomemade?"kTRUE":"kFALSE");
    }
    printf("  fDecayMinEta %f\n",fDecayMinEta);
    printf("  fDecayMaxEta %f\n",fDecayMaxEta);
    printf("  fDecayMinPt %f\n",fDecayMinPt);
    printf("  fDecayMaxDCAdaughters %f\n",fDecayMaxDCAdaughters);
    printf("  fDecayMinCosinePointingAngleXY %f\n",fDecayMinCosinePointingAngleXY);
    printf("  fDecayMinQt %f\n",fDecayMinQt);
    printf("  fDecayAPCutPie %s\n",fDecayAPCutPie?"kTRUE":"kFALSE");
    printf("  fDecayStopPIDAtPt %f\n",fDecayStopPIDAtPt);
    printf("  fDecayMinRadXY %f\n",fDecayMinRadXY);
    printf("  fDecayMaxDecayLength %f\n",fDecayMaxDecayLength);
    printf("  fDecayMaxProductIPXY %f\n",fDecayMaxProductIPXY);
    printf("  fDecayMaxRapidity %f\n",fDecayMaxRapidity);
  }
  printf("  fDaughterUnTag %s\n",fDaughterUnTag?"kTRUE":"kFALSE");
  printf("  fDaughterMinEta %f\n",fDaughterMinEta);
  printf("  fDaughterMaxEta %f\n",fDaughterMaxEta);
  printf("  fDaughterMinPt %f\n",fDaughterMinPt);
  printf("  fDaughterMinNClsTPC %d\n",fDaughterMinNClsTPC);
  printf("  fDaughterMinXRows %d\n",fDaughterMinXRows);
  printf("  fDaughterMaxChi2PerNClsTPC %f\n",fDaughterMaxChi2PerNClsTPC);
  printf("  fDaughterMinXRowsOverNClsFTPC %f\n",fDaughterMinXRowsOverNClsFTPC);
  printf("  fDaughterMinImpactParameterXY %f\n",fDaughterMinImpactParameterXY);
  printf("  fDaughterMaxNSigmaPID %f\n",fDaughterMaxNSigmaPID);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::UserCreateOutputObjects() {
  //UserCreateOutputObjects
  if(fDebug) PrintConfig();
  fList=new TList();
  fList->SetOwner();
  AddQAEvents();
  AddQACandidates();
  if(fReadESD) MakeFilterBits();

  AliFlowCommonConstants *cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(3000); cc->SetMultMin(0);   cc->SetMultMax(30000);
  cc->SetNbinsPt(100); cc->SetPtMin(0.0);   cc->SetPtMax(20.0);
  cc->SetNbinsPhi(100);  cc->SetPhiMin(0.0);  cc->SetPhiMax(TMath::TwoPi());
  cc->SetNbinsEta(100);  cc->SetEtaMin(-5.0); cc->SetEtaMax(+5.0);
  cc->SetNbinsQ(100);    cc->SetQMin(0.0);    cc->SetQMax(3.0);
  cc->SetNbinsMass(fMassBins);
  cc->SetMassMin(fMinMass);
  cc->SetMassMax(fMaxMass);

  //loading pid response
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  if(fUseFP) {
    fTPCevent = new AliFlowEvent(100);
    fVZEevent = new AliFlowEvent(100);
    //array of candidates
    fCandidates = new TObjArray(100);
    fCandidates->SetOwner();
  }
  PostData(1,fList);
  if(fUseFP) { // for connection to the flow package
    PostData(2,fTPCevent);
    PostData(3,fVZEevent);
  }

  gRandom->SetSeed();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MyUserCreateOutputObjects() {
  TList *tList;
  TH1D *tH1D;
  TH2D *tH2D;

  //reconstruction
  if(fReadESD) {
    tList=new TList(); tList->SetName("ESD_TrkAll"); tList->SetOwner(); AddTrackSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("ESD_TrkSel"); tList->SetOwner(); AddTrackSpy(tList); fList->Add(tList);
    tH2D = new TH2D("NPAIR", "NPAIR;NPOS;NNEG",1000,0,5000,1000,0,5000); tList->Add(tH2D);
    tH2D = new TH2D("PtIPXY","PtIPXY;Pt;IPxy", 100,0,10,200,-10,+10); tList->Add(tH2D);
  }
  //aod prefilter candidates
  tList=new TList(); tList->SetName("V0SAll"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
  tH2D = new TH2D("V0SADC","V0S AFTER DAUGHTER CUTS;V0ALL;V0IMW",100,0,1000,100,0,1000); tList->Add(tH2D);
  tList=new TList(); tList->SetName("AllDau"); tList->SetOwner(); AddTrackSpy(tList); fList->Add(tList);
  //candidates
  tList=new TList(); tList->SetName("V0SSel"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
  tList=new TList(); tList->SetName("SelDau"); tList->SetOwner(); AddTrackSpy(tList); fList->Add(tList);
  //flow
  if(!fSkipVn) {
    tList=new TList(); tList->SetName("V0SAllVn"); tList->SetOwner(); AddDecayVn(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("V0SSelVn"); tList->SetOwner(); AddDecayVn(tList); fList->Add(tList);
  }
  // IN-OUT
  if(fQAlevel>1) {
    tList=new TList(); tList->SetName("V0SAllIP"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("V0SAllOP"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("V0SSelIP"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("V0SSelOP"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
  }
  //match
  if(fReadMC) {
    tList=new TList(); tList->SetName("STATMC"); tList->SetOwner(); fList->Add(tList);
    tH1D = new TH1D("Events", "Events",5,0.5,5.5); tList->Add(tH1D);
    tH1D->GetXaxis()->SetBinLabel(1,"Selected events");
    tH1D->GetXaxis()->SetBinLabel(2,"Stack found");
    tH1D->GetXaxis()->SetBinLabel(3,"Daughters in stack");
    tH1D->GetXaxis()->SetBinLabel(4,"Correspond to decay");
    tH1D->GetXaxis()->SetBinLabel(5,"Decay has mother");
    tList=new TList(); tList->SetName("Mth"); tList->SetOwner(); AddCandidatesSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("MthDau"); tList->SetOwner(); AddTrackSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("MthPosPos"); tList->SetOwner(); AddCandidatesSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("MthNegNeg"); tList->SetOwner(); AddCandidatesSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("MthPosNeg"); tList->SetOwner(); AddCandidatesSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("MthNegDau"); tList->SetOwner(); AddTrackSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("MthPosDau"); tList->SetOwner(); AddTrackSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("MthFeedDown"); tList->SetOwner(); AddCandidatesSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("UnMth"); tList->SetOwner(); AddCandidatesSpy(tList,false); fList->Add(tList);
    tList=new TList(); tList->SetName("UnMthDau"); tList->SetOwner(); AddTrackSpy(tList,false); fList->Add(tList);
    if(!fSkipVn) {
      tList=new TList(); tList->SetName("V0SMthVn"); tList->SetOwner(); AddDecayVn(tList); fList->Add(tList);
      tList=new TList(); tList->SetName("V0SMthPosPosVn"); tList->SetOwner(); AddDecayVn(tList); fList->Add(tList);
      tList=new TList(); tList->SetName("V0SMthNegNegVn"); tList->SetOwner(); AddDecayVn(tList); fList->Add(tList);
      tList=new TList(); tList->SetName("V0SMthPosNegVn"); tList->SetOwner(); AddDecayVn(tList); fList->Add(tList);
      tList=new TList(); tList->SetName("V0SUnMthVn"); tList->SetOwner(); AddDecayVn(tList); fList->Add(tList);
    }
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddQAEvents() {
  // function to add event qa
  TH1D *tH1D;
  TProfile *tProfile;
  TList *tQAEvents=new TList();
  tQAEvents->SetName("Event");
  tQAEvents->SetOwner();
  tH1D = new TH1D("Events","Number of Events",6,0,6); tQAEvents->Add(tH1D);
  tH1D->GetXaxis()->SetBinLabel(1,"exec");
  tH1D->GetXaxis()->SetBinLabel(2,"userexec");
  tH1D->GetXaxis()->SetBinLabel(3,"reached");
  tH1D->GetXaxis()->SetBinLabel(4,"selected");
  tH1D->GetXaxis()->SetBinLabel(5,"rejectedByLowQw");
  tH1D->GetXaxis()->SetBinLabel(6,"rejectedByErrorLoadVZEcal");
  tProfile = new TProfile("Configuration","Configuration",10,0.5,10.5); tQAEvents->Add(tProfile);
  tProfile->Fill(kSpecie,fSpecie,1);
  tProfile->GetXaxis()->SetBinLabel(kSpecie,"fSpecie");
  tProfile->Fill(kHarmonic,fHarmonic,1);
  tProfile->GetXaxis()->SetBinLabel(kHarmonic,"fHarmonic");
  tProfile->Fill(kReadMC,fReadMC,1);
  tProfile->GetXaxis()->SetBinLabel(kReadMC,"fReadMC");
  tProfile->Fill(kSkipSelection,fSkipSelection,1);
  tProfile->GetXaxis()->SetBinLabel(kSkipSelection,"fSkipSelection");
  tH1D = new TH1D("POI","POIs;multiplicity",800,0,800);         tQAEvents->Add(tH1D);
  tH1D = new TH1D("UNTAG","UNTAG;Untagged Daughters",800,0,800);tQAEvents->Add(tH1D);
  tH1D = new TH1D("RealTime","RealTime;LogT sec",2000,-10,+10); tQAEvents->Add(tH1D);
  fList->Add(tQAEvents);
  AddEventSpy("EventsRaw");
  AddEventSpy("EventsReached");
  AddEventSpy("EventsSelected");
  AddEventSpy("EventsAnalyzed");
  AddMakeQSpy();
  AddVZEQA();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddEventSpy(TString name) {
  TH1D *tH1D;
  TH2D *tH2D;
  TList *tList=new TList();
  tList->SetName(name.Data());
  tList->SetOwner();
  tH2D = new TH2D("VTXZ","VTXZ;PriVtxZ;SPDVtxZ",60,-25,+25,60,-25,+25); tList->Add( tH2D );
  tH2D = new TH2D("CCCC","CCCC;V0M;TRK",60,-10,110,60,-10,110);         tList->Add( tH2D );
  tH2D = new TH2D("HYBTPC","HYBTPC;TPC ONLY;HYBRID",100,0,3000,100,0,3000); tList->Add( tH2D );
  tH1D = new TH1D("HYBTPCRat","HYBTPCRat;TPC/HYB",120,0.2,2.2); tList->Add( tH1D );
  tH2D = new TH2D("SPDVZE","SPDVZE;SPD Tracklets;Total Multiplicity in VZERO",100,0,3500,100,0,25000); tList->Add( tH2D );
  tH1D = new TH1D("SPDVZERat","SPDVZERat;TotalMultiplicityVZERO/SPDTracklets",120,2,+12); tList->Add( tH1D );
  if(fReadMC) {
    tH1D = new TH1D("MCEP","MCEP;MCEP",100,-TMath::TwoPi(),TMath::TwoPi()); tList->Add( tH1D );
  }
  fList->Add(tList);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillEventSpy(TString name) {
  ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject("VTXZ"))->Fill( fPriVtxZ, fSPDVtxZ );
  ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject("CCCC"))->Fill( fV0M, fTRK );
  ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject("HYBTPC"))->Fill( fRefMultTPC, fRefMultHyb );
  if(fRefMultHyb>0)
    ((TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject("HYBTPCRat"))->Fill( double(fRefMultTPC)/double(fRefMultHyb) );
  ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject("SPDVZE"))->Fill( fSPDtracklets, fVZETotM );
  if(fSPDtracklets>0)
    ((TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject("SPDVZERat"))->Fill( fVZETotM/fSPDtracklets );
  if(fReadMC) {
    ((TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject("MCEP"))->Fill( fMCEP );
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddMakeQSpy() {
  TH1D *tH1D;
  TH2D *tH2D;
  TProfile *tPF1;
  TList *tList=new TList();
  tList->SetName("MakeQSpy");
  tList->SetOwner();
  fList->Add(tList);
  tH1D = new TH1D("RFPTPC","TPC Refrence Multiplicity;multiplicity",3000,0,3000);     tList->Add( tH1D );
  tH1D = new TH1D("RFPVZE","VZERO Reference Multiplicity;multiplicity",3000,0,30000); tList->Add( tH1D );
  tH1D = new TH1D("QmTPC","TPC Normalized Q vector;|Q|/#sqrt{M}",360,0,7);   tList->Add( tH1D );
  tH1D = new TH1D("QmVZEA","VZEROA Normalized Q vector;|Q|/#sqrt{W}",360,0,7); tList->Add( tH1D );
  tH1D = new TH1D("QmVZEC","VZEROC Normalized Q vector;|Q|/#sqrt{W}",360,0,7); tList->Add( tH1D );
  tH2D = new TH2D("TPCAllPhiEta","TPCall;Phi;Eta",180,0,TMath::TwoPi(),80,-0.9,+0.9); tList->Add( tH2D );
  tH2D = new TH2D("VZEAllPhiEta","VZEall;Phi;Eta",20,0,TMath::TwoPi(),40,-4.0,+6.0);  tList->Add( tH2D );
  tH1D = new TH1D("TPCPSI","TPCPSI;PSI",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("TPCPSIA","TPCPSIA;PSIA",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("TPCPSIC","TPCPSIC;PSIC",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("VZEPSI","VZEPSI;PSI",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("VZEPSIA","VZEPSIA;PSIA",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("VZEPSIC","VZEPSIC;PSIC",72,0,TMath::Pi()); tList->Add( tH1D );
  tH2D = new TH2D("PSI_TPCAVZEC","PSI_TPCAVZEC",72,0,TMath::Pi(),72,0,TMath::Pi()); tList->Add( tH2D );
  tH2D = new TH2D("PSI_TPCCVZEA","PSI_TPCAVZEC",72,0,TMath::Pi(),72,0,TMath::Pi()); tList->Add( tH2D );
  tH2D = new TH2D("PSI_TPCVZE","PSI_TPCVZE",72,0,TMath::Pi(),72,0,TMath::Pi()); tList->Add( tH2D );
  tPF1 = new TProfile("TPCQm","TPCQm",6,0.5,6.5); tList->Add( tPF1 );
  tPF1->GetXaxis()->SetBinLabel(1,"Qcy"); tPF1->GetXaxis()->SetBinLabel(2,"Qcx");
  tPF1->GetXaxis()->SetBinLabel(3,"Qay"); tPF1->GetXaxis()->SetBinLabel(4,"Qax");
  tPF1->GetXaxis()->SetBinLabel(5,"Qy");  tPF1->GetXaxis()->SetBinLabel(6,"Qx");
  tPF1 = new TProfile("VZEQm","VZEQm",6,0.5,6.5); tList->Add( tPF1 );
  tPF1->GetXaxis()->SetBinLabel(1,"Qcy"); tPF1->GetXaxis()->SetBinLabel(2,"Qcx");
  tPF1->GetXaxis()->SetBinLabel(3,"Qay"); tPF1->GetXaxis()->SetBinLabel(4,"Qax");
  tPF1->GetXaxis()->SetBinLabel(5,"Qy");  tPF1->GetXaxis()->SetBinLabel(6,"Qx");
  tPF1 = new TProfile("QmVZEAQmVZEC","QmVZEAQmVZEC",1,0.5,1.5,"s"); tList->Add( tPF1 );
  tPF1 = new TProfile("QmVZEASQUARED","QmVZEASQUARED",1,0.5,1.5,"s"); tList->Add( tPF1 );
  tPF1 = new TProfile("QmVZECSQUARED","QmVZECSQUARED",1,0.5,1.5,"s"); tList->Add( tPF1 );
  tPF1 = new TProfile("QmTPCQmVZEA","QmTPCQmVZEA",1,0.5,1.5,"s"); tList->Add( tPF1 );
  tPF1 = new TProfile("QmTPCQmVZEC","QmTPCQmVZEC",1,0.5,1.5,"s"); tList->Add( tPF1 );
  tH1D = new TH1D("ChiSquaredVZEA","ChiSquaredVZEC",1,0.5,1.5); tList->Add( tH1D );
  tH1D = new TH1D("ChiSquaredVZEC","ChiSquaredVZEC",1,0.5,1.5); tList->Add( tH1D );
  if(fReadMC) {
    tH1D = new TH1D("PSIMCDIFFTPC","PSIMCDIFFTPC;MC-TPC",72,-TMath::TwoPi(),TMath::TwoPi()); tList->Add( tH1D );
    tH1D = new TH1D("PSIMCDIFFTPCA","PSIMCDIFFTPCA;MC-TPCA",72,-TMath::TwoPi(),TMath::TwoPi()); tList->Add( tH1D );
    tH1D = new TH1D("PSIMCDIFFTPCC","PSIMCDIFFTPCC;MC-TPCC",72,-TMath::TwoPi(),TMath::TwoPi()); tList->Add( tH1D );
    tH1D = new TH1D("PSIMCDIFFVZE","PSIMCDIFFVZE;MC-VZE",72,-TMath::TwoPi(),TMath::TwoPi()); tList->Add( tH1D );
    tH1D = new TH1D("PSIMCDIFFVZEA","PSIMCDIFFVZEA;MC-VZEA",72,-TMath::TwoPi(),TMath::TwoPi()); tList->Add( tH1D );
    tH1D = new TH1D("PSIMCDIFFVZEC","PSIMCDIFFVZEC;MC-VZEC",72,-TMath::TwoPi(),TMath::TwoPi()); tList->Add( tH1D );
  }
  tList=new TList(); tList->SetName("TPCRFPall"); tList->SetOwner(); AddTPCRFPSpy(tList); fList->Add(tList);
  tList=new TList(); tList->SetName("TPCRFPsel"); tList->SetOwner(); AddTPCRFPSpy(tList); fList->Add(tList);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddQACandidates() {
  // function to add histogramming for candidates
  if(fSkipSelection) return;
  TList *tList;
  TH1D *tH1D;

  //charge particles (benchmark)
  if(fSpecie>=90) {
    tList=new TList(); tList->SetName("TrkAll"); tList->SetOwner(); AddTrackSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("TrkSel"); tList->SetOwner(); AddTrackSpy(tList); fList->Add(tList);
    if(!fSkipVn) {
      tList=new TList(); tList->SetName("TrkAllVn"); tList->SetOwner(); AddTrackVn(tList); fList->Add(tList);
      tList=new TList(); tList->SetName("TrkSelVn"); tList->SetOwner(); AddTrackVn(tList); fList->Add(tList);
    }
    //match
    if(fReadMC) {
      tList=new TList(); tList->SetName("STATMC"); tList->SetOwner(); fList->Add(tList);
      tH1D = new TH1D("Events", "Events",3,0.5,3.5); tList->Add(tH1D);
      tH1D->GetXaxis()->SetBinLabel(1,"Selected events");
      tH1D->GetXaxis()->SetBinLabel(2,"Stack found");
      tH1D->GetXaxis()->SetBinLabel(3,"Track in stack");
      tList=new TList(); tList->SetName("Mth"); tList->SetOwner(); AddTrackSpy(tList); fList->Add(tList);
      tList=new TList(); tList->SetName("MthPos"); tList->SetOwner(); AddTrackSpy(tList); fList->Add(tList);
      tList=new TList(); tList->SetName("MthNeg"); tList->SetOwner(); AddTrackSpy(tList); fList->Add(tList);
      if(!fSkipVn) {
	tList=new TList(); tList->SetName("MthVn"); tList->SetOwner(); AddTrackVn(tList); fList->Add(tList);
	tList=new TList(); tList->SetName("MthPosVn"); tList->SetOwner(); AddTrackVn(tList); fList->Add(tList);
	tList=new TList(); tList->SetName("MthNegVn"); tList->SetOwner(); AddTrackVn(tList); fList->Add(tList);
      }
    }
  }
  //stack
  if(fReadMC) {
    tList=new TList(); tList->SetName("MCTPionGenAcc"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTKaonGenAcc"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTK0sGenAcc"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTProtonGenAcc"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTLdaGenAcc"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTPhiGenAcc"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTXiGenAcc");  tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTOmegaGenAcc");  tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTPion"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTKaon"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTK0s"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTLda"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTProton"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
  }
  MyUserCreateOutputObjects();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::Exec(Option_t* option) {
  // bypassing ::exec (needed because of AMPT)
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(0);
  if(fAvoidExec) {
    AliAnalysisTaskFlowStrange::UserExec(option);
  } else {
    AliAnalysisTaskSE::Exec(option);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::UserExec(Option_t *option) {
  // bridge
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(1);
  AliAnalysisTaskFlowStrange::MyUserExec(option);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MyNotifyRun() {
  if(fVZEsave) AddVZEROResponse();
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::CalibrateEvent() {
  if(fVZEsave) SaveVZEROResponse();
  Bool_t okay=kTRUE;
  if(fVZEload) {
    LoadVZEROResponse();
    if(!fVZEResponse) okay = kFALSE;
  }
  return okay;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptAAEvent(AliESDEvent*) {
  // ESD reading discontinued: TO BE UPDATED
  /*
  Double_t acceptEvent=kTRUE;
  Double_t tTPCVtxZ = tESD->GetPrimaryVertexTPC()->GetZ();
  if(tESD->GetPrimaryVertexTPC()->GetNContributors()<=0) return kFALSE;
  Double_t tSPDVtxZ = tESD->GetPrimaryVertexSPD()->GetZ();
  if(tESD->GetPrimaryVertexSPD()->GetNContributors()<=0) return kFALSE;
  // EventCuts
  AliCentrality *cent = tESD->GetCentrality();
  Double_t cc1, cc2;
  cc1 = cent->GetCentralityPercentile("V0M");
  cc2 = cent->GetCentralityPercentile("TRK");
  TString mycent = fCentMethod;
  if(fCentMethod.Contains("V0MTRK")) {
    acceptEvent = TMath::Abs(cc1-cc2)>5.0?kFALSE:acceptEvent; // a la Alex
    mycent = "V0M";
  }
  fThisCent = cent->GetCentralityPercentile( mycent );
  acceptEvent = (fThisCent<fCentPerMin||fThisCent>fCentPerMax)?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tTPCVtxZ-tSPDVtxZ)>0.5?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tTPCVtxZ)>fVertexZcut?kFALSE:acceptEvent;
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("VTXZ"))->Fill( tTPCVtxZ, tSPDVtxZ );
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  // EndOfCuts
  if(acceptEvent) {
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("VTXZ"))->Fill( tTPCVtxZ, tSPDVtxZ );
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  }
  return acceptEvent;
  */
  return kFALSE;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::MinimumRequirementsAA(AliAODEvent *tAOD) {
  fRunNumber = tAOD->GetRunNumber();
  AliCentrality *cent = ((AliVAODHeader*)tAOD->GetHeader())->GetCentralityP();
  fV0M = cent->GetCentralityPercentile("V0M");
  fTRK = cent->GetCentralityPercentile("TRK");
  TString mycent = fCentMethod;
  if(fCentMethod.Contains("V0MTRK")) {
    mycent = "V0M";
  }
  fThisCent = cent->GetCentralityPercentile( mycent );
  fPriVtxZ = tAOD->GetPrimaryVertex()->GetZ();
  fSPDVtxZ = tAOD->GetPrimaryVertexSPD()->GetZ();
  fSPDtracklets = tAOD->GetTracklets()->GetNumberOfTracklets();
  fVZETotM = tAOD->GetVZEROData()->GetMTotV0A() + tAOD->GetVZEROData()->GetMTotV0C();
  int hyb_fb = 272; // for 2010h::AOD086
  if(fRunNumber>=166529&&fRunNumber<=170593) {
    hyb_fb = 768; // for 2011h::AOD145
  }
  fRefMultTPC = RefMult(tAOD,128);
  fRefMultHyb = RefMult(tAOD,hyb_fb);
  if(fReadMC) {
    fMCEP = -999;
    AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(tAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if(mcHeader) {
      fMCEP = mcHeader->GetReactionPlaneAngle();
      if(fAddPiToMCReactionPlane) fMCEP += (gRandom->Rndm()>0.5)*TMath::Pi();
    }
  }
  // centrality selection health
  // cut in Vtx 10 & NContributors
  if(!fSkipCentralitySelection) if(fThisCent<0||fThisCent>100) return kFALSE;
  // vtx z position compatibility within 5 mm
  if(TMath::Abs(fPriVtxZ-fSPDVtxZ)>0.5) return kFALSE;
  if(fExtraEventRejection) {
    // specific cuts for 2010h (AOD086)
    if(fRunNumber>=136851&&fRunNumber<=139517) {
      if(fRefMultTPC>1.118*fRefMultHyb+100) return kFALSE;
      if(fRefMultTPC<1.118*fRefMultHyb-100) return kFALSE;
    }
    // specific cuts for 2011h (AOD145)
    if(fRunNumber>=166529&&fRunNumber<=170593) {
      if(fRefMultTPC>1.205*fRefMultHyb+100) return kFALSE;
      if(fRefMultTPC<1.205*fRefMultHyb-100) return kFALSE;
    }
  }
  return kTRUE;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptAAEvent(AliAODEvent *tAOD) {
  Bool_t minimum = MinimumRequirementsAA(tAOD);
  FillEventSpy("EventsRaw");
  if(!minimum) return kFALSE;

  Double_t acceptEvent=kTRUE;
  TString mycent = fCentMethod;
  if(fCentMethod.Contains("V0MTRK")) {
    acceptEvent = TMath::Abs(fV0M-fTRK)>5.0?kFALSE:acceptEvent;
    mycent = "V0M";
  }
  if(!fSkipCentralitySelection) acceptEvent = (fThisCent<fCentPerMin||fThisCent>fCentPerMax)?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(fPriVtxZ)>fVertexZcut?kFALSE:acceptEvent;
  // HISTOGRAMMING
  FillEventSpy("EventsReached");
  if(acceptEvent) FillEventSpy("EventsSelected");
  return acceptEvent;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptPPEvent(AliAODEvent*) {
  // PP reading discontinued: TO BE UPDATED
  /*
  Double_t acceptEvent=kTRUE;
  Double_t tVtxZ = tAOD->GetPrimaryVertex()->GetZ();
  if(tAOD->GetPrimaryVertex()->GetNContributors()<=0) return kFALSE;
  Double_t tSPDVtxZ = tAOD->GetPrimaryVertexSPD()->GetZ();
  if(tAOD->GetPrimaryVertexSPD()->GetNContributors()<=0) return kFALSE;
  // EventCuts
  AliCentrality *cent = tAOD->GetHeader()->GetCentralityP();
  Double_t cc1, cc2;
  cc1 = cent->GetCentralityPercentile("V0M");
  cc2 = cent->GetCentralityPercentile("TRK");
  fThisCent = GetReferenceMultiplicity();
  //for pp i use fCentPerXXX to select on multiplicity
  acceptEvent = (fThisCent<fCentPerMin||fThisCent>fCentPerMax)?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tVtxZ-tSPDVtxZ)>0.5?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tVtxZ)>fVertexZcut?kFALSE:acceptEvent;
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("VTXZ"))->Fill( tVtxZ, tSPDVtxZ );
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  // EndOfCuts
  if(acceptEvent) {
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("VTXZ"))->Fill( tVtxZ, tSPDVtxZ );
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  }
  return acceptEvent;
  */
  return kFALSE;
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::GetReferenceMultiplicity() { //toberefined
  AliAODEvent *tAOD = (AliAODEvent *) InputEvent();
  if(!tAOD) return -1;
  AliAODTrack *track;
  Int_t rawN = tAOD->GetNumberOfTracks();
  Int_t ref=0;
  for(Int_t id=0; id!=rawN; ++id) {
    track = dynamic_cast<AliAODTrack*>(tAOD->GetTrack(id));
    if(!track) {
        AliFatal("Not a standard AOD");
        continue;
    }
    if(!track->TestFilterBit(fRFPFilterBit)) continue;
    ++ref;
  }
  return ref;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptPAEvent(AliAODEvent*) {
  // PA reading discontinued: TO BE UPDATED
  /*
  //if(aod->GetHeader()->GetEventNumberESDFile() == 0) return; //rejecting first chunk NOT NEEDED ANYMORE
  Int_t bc2 = ((AliVAODHeader*)tAOD->GetHeader())->GetIRInt2ClosestInteractionMap();
  if(bc2!=0) return kFALSE;
  Int_t bc1 = ((AliVAODHeader*)tAOD->GetHeader())->GetIRInt1ClosestInteractionMap();
  if(bc1!=0) return kFALSE;
  Short_t isPileup = tAOD->IsPileupFromSPD(5);
  if(isPileup!=0) return kFALSE;
  if(tAOD->GetHeader()->GetRefMultiplicityComb08()<0) return kFALSE;

  const AliAODVertex* spdVtx = tAOD->GetPrimaryVertexSPD();
  if(!spdVtx) return kFALSE;
  if(spdVtx->GetNContributors()<=0) return kFALSE;

  const AliAODVertex* tpcVtx=NULL;
  Int_t nVertices = tAOD->GetNumberOfVertices();
  for(Int_t iVertices = 0; iVertices < nVertices; iVertices++){
    const AliAODVertex* vertex = tAOD->GetVertex(iVertices);
    if (vertex->GetType() != AliAODVertex::kMainTPC) continue;
    tpcVtx = vertex;
  }
  if(!tpcVtx) return kFALSE;
  if(tpcVtx->GetNContributors()<=0) return kFALSE;
  Double_t tTPCVtxZ = tpcVtx->GetZ();
  Double_t tSPDVtxZ = spdVtx->GetZ();
  if (TMath::Abs(tSPDVtxZ - tTPCVtxZ)>2.0) return kFALSE;
  if(plpMV(tAOD)) return kFALSE;

  Double_t acceptEvent=kTRUE;
  // EventCuts
  AliCentrality *cent = tAOD->GetHeader()->GetCentralityP();
  Double_t cc1, cc2;
  cc1 = cent->GetCentralityPercentile("V0M");
  cc2 = cent->GetCentralityPercentile("TRK");
  if(fCentMethod.Contains("V0MTRK")) fCentMethod = "V0M";
  fThisCent = cent->GetCentralityPercentile( fCentMethod );
  acceptEvent = (fThisCent<fCentPerMin||fThisCent>fCentPerMax)?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tTPCVtxZ)>fVertexZcut?kFALSE:acceptEvent;
  // EndOfCuts
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("VTXZ"))->Fill( tTPCVtxZ, tSPDVtxZ );
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  if(acceptEvent) {
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("VTXZ"))->Fill( tTPCVtxZ, tSPDVtxZ );
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("CCCC"))->Fill( cc1, cc2 );    
  }
  return acceptEvent;
  */
  return kFALSE;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MyUserExec(Option_t *) {
  // MAIN ROUTINE
  TStopwatch tTime;
  tTime.Start();
  if(fDebug) {
    printf("****************\n");
    printf("****************\n");
    printf("**::MyUserExec()\n");
  }
  if(fUseFP) fCandidates->SetLast(-1);
  AliESDEvent *tESD=dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent *tAOD=dynamic_cast<AliAODEvent*>(InputEvent());
  Int_t prevRun = fRunNumber;
  //=>check event
  Bool_t acceptEvent=kFALSE;
  if(fReadESD) {
    if(!tESD) {ResetContainers(); Publish(); return;}
    acceptEvent = fRunOnpp?kFALSE:fRunOnpA?kFALSE:AcceptAAEvent(tESD);
  } else {
    if(!tAOD) {ResetContainers(); Publish(); return;}
    acceptEvent = fRunOnpp?AcceptPPEvent(tAOD):fRunOnpA?AcceptPAEvent(tAOD):AcceptAAEvent(tAOD);
  }
  if(prevRun!=fRunNumber) {
    MyNotifyRun();
  }
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(2);
  //=>does the event clear?
  if(!acceptEvent) {ResetContainers(); Publish(); return;}
  // healthy event incomming
  if( !CalibrateEvent() ) { // saves/retrieves/qas VZEROCAL
    ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(5);
    ResetContainers(); Publish(); return; // issue retrieving callibration
  }
  // loads Q vectors
  MakeQVectors();
  if(fPsi2<-0.1) {
    ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(4);
    ResetContainers(); Publish(); return;
  }
  //}
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(3);
  //=>great, lets do our stuff!
  FillEventSpy("EventsAnalyzed");
  FillVZEQA();
  //=>load candidates
  if(!fSkipSelection) {
    if(fReadESD) {
      ReadFromESD(tESD);
    } else {
      if(fSpecie<10) ReadFromAODv0(tAOD);
      else ChargeParticles(tAOD);
    }
    if(fUseFP) AddCandidates();
    //=>flow
    //=>done
  }
  tTime.Stop();
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("RealTime"))->Fill( TMath::Log( tTime.RealTime() ) );
  Publish();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::Publish() {
  PostData(1,fList);
  if(fUseFP) {
    PostData(2,fTPCevent);
    PostData(3,fVZEevent);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ReadFromESD(AliESDEvent *tESD) {
  AliStack *stack=NULL;
  if(fReadMC) {
    AliMCEvent *mcevent=NULL;
    mcevent = MCEvent();
    if(mcevent) stack = mcevent->Stack();
  }

  Int_t num = tESD->GetNumberOfTracks();
  AliESDtrack *myTrack;
  Int_t plist[3000], nlist[3000], np=0, nn=0;
  Double_t pd0[3000], nd0[3000];
  for (Int_t i=0; i!=num; ++i) {
    myTrack = (AliESDtrack*) tESD->GetTrack(i);
    if(!myTrack) continue;
    LoadTrack(myTrack);
    FillTrackSpy("ESD_TrkAll");
    if(!AcceptDaughter()) continue;
    FillTrackSpy("ESD_TrkSel");
    ((TH2D*)((TList*)fList->FindObject("ESD_TrkSel"))->FindObject("PtIPXY" ))->Fill( myTrack->Pt(), fDaughterImpactParameterXY );
    if( myTrack->Charge()>0 ) {
      pd0[np] = fDaughterImpactParameterXY;
      plist[np++] = i;
    } else {
      nd0[nn] = fDaughterImpactParameterXY;
      nlist[nn++] = i;
    }
  }
  ((TH1D*)((TList*)fList->FindObject("ESD_TrkSel"))->FindObject("NPAIR" ))->Fill( np,nn );
  const AliESDVertex *vtx = tESD->GetPrimaryVertex();
  AliESDtrack *pT, *nT;
  for(int p=0; p!=np; ++p) {
    pT = (AliESDtrack*) tESD->GetTrack( plist[p] );
    for(int n=0; n!=nn; ++n) {
      nT = (AliESDtrack*) tESD->GetTrack( nlist[n] );
      fDecayProductIPXY = pd0[p]*nd0[n];
      AliExternalTrackParam pETP(*pT), nETP(*nT);
      Double_t xa, xb;
      pETP.GetDCA(&nETP,tESD->GetMagneticField(),xa,xb);
      fDecayDCAdaughters = pETP.PropagateToDCA(&nETP,tESD->GetMagneticField());
      AliESDv0 vertex( nETP,nlist[n], pETP,plist[p] );
      fDecayCosinePointingAngleXY = CosThetaPointXY( &vertex, vtx );
      fDecayRadXY = DecayLengthXY( &vertex, vtx );
      fDecayPt = vertex.Pt();
      fDecayPhi = vertex.Phi();
      fDecayEta = vertex.Eta();
      Double_t pmx, pmy, pmz, nmx, nmy, nmz;
      vertex.GetNPxPyPz(nmx,nmy,nmz);
      vertex.GetPPxPyPz(pmx,pmy,pmz);
      TVector3 mom1(pmx,pmy,pmz), mom2(nmx,nmy,nmz), mom(vertex.Px(),vertex.Py(),vertex.Pz());
      Double_t qlpos = mom1.Dot(mom)/mom.Mag();
      Double_t qlneg = mom2.Dot(mom)/mom.Mag();
      fDecayQt = mom1.Perp(mom);
      fDecayAlpha = (qlpos-qlneg)/(qlpos+qlneg);
      Double_t mpi = 0.13957018;
      if(fSpecie==0) {
        Double_t eppi = TMath::Sqrt( mpi*mpi + pmx*pmx + pmy*pmy + pmz*pmz );
        Double_t enpi = TMath::Sqrt( mpi*mpi + nmx*nmx + nmy*nmy + nmz*nmz );
        fDecayMass = TMath::Sqrt( mpi*mpi + mpi*mpi + 2*(eppi*enpi - pmx*nmx - pmy*nmy - pmz*nmz ) );
        fDecayRapidity = vertex.RapK0Short();
      } else {
        Double_t mpr = 0.938272013;
        Double_t epi, epr;
        if(fDecayAlpha>0) {
          epr = TMath::Sqrt( mpr*mpr + pmx*pmx + pmy*pmy + pmz*pmz );
          epi = TMath::Sqrt( mpi*mpi + nmx*nmx + nmy*nmy + nmz*nmz );
        } else {
          epi = TMath::Sqrt( mpi*mpi + pmx*pmx + pmy*pmy + pmz*pmz );
          epr = TMath::Sqrt( mpr*mpr + nmx*nmx + nmy*nmy + nmz*nmz );
        }
        fDecayMass = TMath::Sqrt( mpi*mpi + mpr*mpr + 2*(epi*epr - pmx*nmx - pmy*nmy - pmz*nmz ) );
        fDecayRapidity = vertex.RapLambda();
      }
      Double_t energy = TMath::Sqrt( fDecayMass*fDecayMass + vertex.Px()*vertex.Px() + vertex.Py()*vertex.Py() + vertex.Pz()*vertex.Pz() );
      Double_t gamma = energy/fDecayMass;
      fDecayDecayLength = DecayLength( &vertex, vtx )/gamma;
      fDecayDecayLengthLab = DecayLength( &vertex, vtx );
      Double_t dPHI = fDecayPhi;
      Double_t dDPHI = dPHI - fPsi2;
      if( dDPHI < 0 ) dDPHI += TMath::TwoPi();
      if( dDPHI > TMath::Pi() ) dDPHI = TMath::TwoPi()-dDPHI;
      if(fQAlevel>1) {
        if( (dDPHI>TMath::PiOver4()) && (dDPHI<3*TMath::PiOver4()) ) FillCandidateSpy("V0SAllOP");
        else FillCandidateSpy("V0SAllIP");
      }
      FillCandidateSpy("V0SAll");
      ((TH2D*)((TList*)fList->FindObject("V0SAll"))->FindObject("D0PD0N"))->Fill( pd0[p],nd0[n] );
      ((TH2D*)((TList*)fList->FindObject("V0SAll"))->FindObject("XPOSXNEG"))->Fill( xa, xb );
      if(!AcceptCandidate()) continue;
      if(fDecayMass<fMinMass) continue;
      if(fDecayMass>fMaxMass) continue;
      // PID missing
      if(fQAlevel>1) {
        if( (dDPHI>TMath::PiOver4()) && (dDPHI<3*TMath::PiOver4()) ) FillCandidateSpy("V0SSelOP");
        else FillCandidateSpy("V0SSelIP");
      }
      FillCandidateSpy("V0SSel");
      ((TH2D*)((TList*)fList->FindObject("V0SSel"))->FindObject("D0PD0N"))->Fill( pd0[p],nd0[n] );
      ((TH2D*)((TList*)fList->FindObject("V0SSel"))->FindObject("XPOSXNEG"))->Fill( xa, xb );

      fDecayIDneg = nT->GetID();
      fDecayIDpos = pT->GetID();
      if(fUseFP) MakeTrack();
      LoadTrack(pT); FillTrackSpy("SelDau");
      LoadTrack(nT); FillTrackSpy("SelDau");

      //===== BEGIN OF MCMATCH
      if(stack) {
        bool matched = false;
        Int_t labelpos = pT->GetLabel();
        Int_t labelneg = nT->GetLabel();
        Double_t rOri=-1;
        if( labelpos>0 && labelneg>0 ) {
          TParticle *mcpos = stack->Particle( labelpos );
          TParticle *mcneg = stack->Particle( labelneg );
          Int_t pdgRecPos = mcpos->GetPdgCode();
          Int_t pdgRecNeg = mcneg->GetPdgCode();
          if( pdgRecPos==211&&pdgRecNeg==-211 ) if(mcpos->GetMother(0)>0) {
            if( mcpos->GetMother(0)==mcneg->GetMother(0) ) {
              TParticle *mcmot = stack->Particle( mcpos->GetMother(0) );
              rOri = TMath::Sqrt( mcmot->Vx()*mcmot->Vx() + mcmot->Vy()*mcmot->Vy() );
              if( TMath::Abs(mcmot->GetPdgCode())==310) {
                if(mcmot->GetNDaughters()==2) matched=true;
              }
            }
          }
        }
        if(matched) {
          FillCandidateSpy("Mth");
          ((TH2D*)((TList*)fList->FindObject("Mth"))->FindObject("D0PD0N"))->Fill( pd0[p],nd0[n] );
          ((TH2D*)((TList*)fList->FindObject("Mth"))->FindObject("XPOSXNEG"))->Fill( xa, xb );
          ((TH1D*)((TList*)fList->FindObject("Mth"))->FindObject("MCOrigin"))->Fill( rOri );
          LoadTrack(pT); FillTrackSpy("MthDau");
          LoadTrack(nT); FillTrackSpy("MthDau");
        }
      }
      //===== END OF MCMATCH
    }
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ReadStack(TClonesArray* mcArray) {
  if(!mcArray) return;
  AliAODMCParticle *myMCTrack;//, *iMCDau, *jMCDau;
  for(int i=0; i!=mcArray->GetEntriesFast(); ++i) {
    myMCTrack = dynamic_cast<AliAODMCParticle*>(mcArray->At( i ));
    if(!myMCTrack) continue;
    /*
    int tPDG=310;
    if(fSpecie>0) tPDG = 3122;
    if( TMath::Abs(myMCTrack->PdgCode())==tPDG )
      if( myMCTrack->GetNDaughters() == 2 ) {
        Int_t iDau = myMCTrack->GetDaughter(0);
        Int_t jDau = myMCTrack->GetDaughter(1);
        AliAODMCParticle *posDau=NULL;
        AliAODMCParticle *negDau=NULL;
        if(iDau>0&&jDau>0) {
          iMCDau = dynamic_cast<AliAODMCParticle*>(mcArray->At( iDau ));
          jMCDau = dynamic_cast<AliAODMCParticle*>(mcArray->At( jDau ));
          if(iMCDau) {
            if(iMCDau->Charge()>0) posDau=iMCDau;
            else negDau=iMCDau;
          }
          if(jMCDau) {
            if(jMCDau->Charge()>0) posDau=jMCDau;
            else negDau=jMCDau;
          }
        } //got two daughters
        if(posDau&&negDau) {
          Double_t dx = myMCTrack->Xv() - posDau->Xv();
          Double_t dy = myMCTrack->Yv() - posDau->Yv();
          Double_t dz = myMCTrack->Zv() - posDau->Zv();
          fDecayRadXY = TMath::Sqrt( dx*dx + dy*dy );
          TVector3 momPos(posDau->Px(),posDau->Py(),posDau->Pz());
          TVector3 momNeg(negDau->Px(),negDau->Py(),negDau->Pz());
          TVector3 momTot(myMCTrack->Px(),myMCTrack->Py(),myMCTrack->Pz());
          Double_t qlpos = momPos.Dot(momTot)/momTot.Mag();
          Double_t qlneg = momNeg.Dot(momTot)/momTot.Mag();
          fDecayQt = momPos.Perp(momTot);
          fDecayAlpha = 1.-2./(1.+qlpos/qlneg);
          fDecayMass = myMCTrack->GetCalcMass();
          Double_t energy = myMCTrack->E();
          Double_t gamma = energy/fDecayMass;
          fDecayDecayLength = TMath::Sqrt(dx*dx+dy*dy+dz*dz)/gamma;
          fDecayPt = myMCTrack->Pt();
          fDecayPhi = myMCTrack->Phi();
          fDecayEta = myMCTrack->Eta();
          fDecayRapidity = myMCTrack->Y();
          fDecayDCAdaughters = 0;
          fDecayCosinePointingAngleXY = 1;
          fDecayProductIPXY = -1;
          if(AcceptCandidate()) FillCandidateSpy("GenTru");
        }
      } // k0/lda with two daughters
    */
    //==== BEGIN TRACK CUTS
    if(myMCTrack->Eta()<-0.8) continue;
    if(myMCTrack->Eta()>+0.8) continue;
    if(myMCTrack->Y()<-0.5) continue;
    if(myMCTrack->Y()>+0.5) continue;
    //==== END TRACK CUTS
    switch( TMath::Abs(myMCTrack->PdgCode()) ) {
    case (211): //pi
      FillMCParticleSpy( "MCTPion", myMCTrack );
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTPionGenAcc", myMCTrack );
      break;
    case (321): //kaon
      FillMCParticleSpy( "MCTKaon", myMCTrack );
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTKaonGenAcc", myMCTrack );
      break;
    case (310): //k0s
      FillMCParticleSpy( "MCTK0s", myMCTrack );
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTK0sGenAcc", myMCTrack );
      break;
    case (2212): //proton
      FillMCParticleSpy( "MCTProton", myMCTrack );
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTProtonGenAcc", myMCTrack );
      break;
    case (3122): //lda
      FillMCParticleSpy( "MCTLda", myMCTrack );
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTLdaGenAcc", myMCTrack );
      break;
    case (333): //phi
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTPhiGenAcc", myMCTrack );
      break;
    case (3312): //xi
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTXiGenAcc", myMCTrack );
      break;
    case (3334): //omega
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTOmegaGenAcc", myMCTrack );
      break;
    }
  }
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::CosThetaPointXY(AliESDv0 *me, const AliVVertex *vtx) {
  TVector3 mom( me->Px(), me->Py(), 0 );
  TVector3 fli( me->Xv()-vtx->GetX(), me->Yv()-vtx->GetY(), 0 );
  Double_t ctp = mom.Dot(fli) / mom.Mag() / fli.Mag();
  return ctp;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::CosThetaPointXY(AliAODv0 *me, const AliVVertex *vtx) {
  TVector3 mom( me->Px(), me->Py(), 0 );
  TVector3 fli( me->Xv()-vtx->GetX(), me->Yv()-vtx->GetY(), 0 );
  Double_t ctp = mom.Dot(fli) / mom.Mag() / fli.Mag();
  return ctp;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::DecayLengthXY(AliESDv0 *me, const AliVVertex *vtx) {
  Double_t dx = me->Xv()-vtx->GetX();
  Double_t dy = me->Yv()-vtx->GetY();
  Double_t dxy = TMath::Sqrt( dx*dx + dy*dy );
  return dxy;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::DecayLengthXY(AliAODv0 *me, const AliVVertex *vtx) {
  Double_t dx = me->Xv()-vtx->GetX();
  Double_t dy = me->Yv()-vtx->GetY();
  Double_t dxy = TMath::Sqrt( dx*dx + dy*dy );
  return dxy;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::DecayLength(AliESDv0 *me, const AliVVertex *vtx) {
  Double_t dx = me->Xv()-vtx->GetX();
  Double_t dy = me->Yv()-vtx->GetY();
  Double_t dz = me->Zv()-vtx->GetZ();
  Double_t dxy = TMath::Sqrt( dx*dx + dy*dy + dz*dz );
  return dxy;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::DecayLength(AliAODv0 *me, const AliVVertex *vtx) {
  Double_t dx = me->Xv()-vtx->GetX();
  Double_t dy = me->Yv()-vtx->GetY();
  Double_t dz = me->Zv()-vtx->GetZ();
  Double_t dxy = TMath::Sqrt( dx*dx + dy*dy + dz*dz );
  return dxy;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ReadFromAODv0(AliAODEvent *tAOD) {
  TClonesArray* mcArray=NULL;
  if(fReadMC) {
    mcArray = dynamic_cast<TClonesArray*>(tAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    ReadStack(mcArray);
  }

  Int_t nV0s = tAOD->GetNumberOfV0s();
  AliAODv0 *myV0;
  Int_t v0all=0, v0imw=0;
  for (Int_t i=0; i!=nV0s; ++i) {
    myV0 = (AliAODv0*) tAOD->GetV0(i);
    if(!myV0) continue;
    if(!fOnline) if(myV0->GetOnFlyStatus() ) continue;
    if(fOnline) if(!myV0->GetOnFlyStatus() ) continue;

    fDecayPt = myV0->Pt();
    fDecayPhi = myV0->Phi();
    fDecayEta = myV0->Eta();

    AliAODTrack *iT, *jT;
    AliAODVertex *vtx = tAOD->GetPrimaryVertex();
    Double_t pos[3],cov[6];
    vtx->GetXYZ(pos);
    vtx->GetCovarianceMatrix(cov);
    const AliESDVertex vESD(pos,cov,100.,100);
    // TESTING CHARGE
    int iPos, iNeg;
    iT=(AliAODTrack*) myV0->GetDaughter(0);
    if(iT->Charge()>0) {
      iPos = 0; iNeg = 1;
    } else {
      iPos = 1; iNeg = 0;
    }
    // END OF TEST

    iT=(AliAODTrack*) myV0->GetDaughter(iPos); // positive
    AliESDtrack ieT( iT );
    ieT.SetTPCClusterMap( iT->GetTPCClusterMap() );
    ieT.SetTPCSharedMap( iT->GetTPCSharedMap() );
    ieT.SetTPCPointsF( iT->GetTPCNclsF() );
    ieT.PropagateToDCA(&vESD, tAOD->GetMagneticField(), 100);
    LoadTrack(&ieT,iT->Chi2perNDF());
    Float_t ip[2];
    ieT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];
    fDecayIPpos = fDaughterImpactParameterXY; //ieT.GetD(pos[0], pos[1], tAOD->GetMagneticField());
    FillTrackSpy("AllDau");
    if(!AcceptDaughter(fDecayPt<2.0?kTRUE:kFALSE)) continue;

    jT=(AliAODTrack*) myV0->GetDaughter(iNeg); // negative
    AliESDtrack jeT( jT );
    jeT.SetTPCClusterMap( jT->GetTPCClusterMap() );
    jeT.SetTPCSharedMap( jT->GetTPCSharedMap() );
    jeT.SetTPCPointsF( jT->GetTPCNclsF() );
    jeT.PropagateToDCA(&vESD, tAOD->GetMagneticField(), 100);
    LoadTrack(&jeT,jT->Chi2perNDF());
    jeT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];
    fDecayIPneg = fDaughterImpactParameterXY; //jeT.GetD(pos[0], pos[1], tAOD->GetMagneticField());
    FillTrackSpy("AllDau");
    if(!AcceptDaughter(fDecayPt<2.0?kTRUE:kFALSE)) continue;

    if( fExcludeTPCEdges ) {
      if( IsAtTPCEdge(iT->Phi(),iT->Pt(),+1,tAOD->GetMagneticField()) ) continue;
      if( IsAtTPCEdge(jT->Phi(),jT->Pt(),-1,tAOD->GetMagneticField()) ) continue;
    }
    ieT.GetDCA(&jeT,tAOD->GetMagneticField(),fDecayXpos,fDecayXneg);
    /*
    // cutting out population close to TPC edges :: strange excess saw in 2010
    if( fExcludeTPCEdges ) {
    Double_t phimod = myV0->Phi();
    int sectors[6] = {5,6,9,10,11,12};
    for(int ii=0; ii!=6; ++ii)
    if( (phimod<(sectors[ii]+1)*TMath::Pi()/9.0) && (phimod>sectors[ii]*TMath::Pi()/9.0) )
    return 0;
    }
    */
    if(fSpecie==0)
      fDecayRapidity = myV0->RapK0Short();
    else
      fDecayRapidity = myV0->RapLambda();
    fDecayDCAdaughters = myV0->DcaV0Daughters();
    fDecayCosinePointingAngleXY = CosThetaPointXY( myV0, vtx );
    fDecayRadXY = DecayLengthXY( myV0, vtx );
    fDecayProductIPXY = fDecayIPpos*fDecayIPneg;
    fDecayQt = myV0->PtArmV0();
    fDecayAlpha = myV0->AlphaV0(); // AlphaV0 -> AODRecoDecat::Alpha -> return 1.-2./(1.+QlProng(0)/QlProng(1));
    if(myV0->ChargeProng(iPos)<0) fDecayAlpha = -fDecayAlpha; // protects for a change in convention
    fDecayPt = myV0->Pt();
    fDecayEta = myV0->Eta();
    if( fSpecie==0 ) {
      fDecayMass = myV0->MassK0Short();
    } else {
      if(fDecayAlpha>0) fDecayMass = myV0->MassLambda();
      else fDecayMass = myV0->MassAntiLambda();
    }
    v0all++;
    if(fDecayMass<fMinMass) continue;
    if(fDecayMass>fMaxMass) continue;
    v0imw++;
    Double_t energy = TMath::Sqrt( fDecayMass*fDecayMass + myV0->Px()*myV0->Px() + myV0->Py()*myV0->Py() + myV0->Pz()*myV0->Pz() );
    Double_t gamma = energy/fDecayMass;
    fDecayDecayLength = DecayLength( myV0, vtx )/gamma;
    fDecayDecayLengthLab = DecayLength( myV0, vtx );
    Double_t dPHI = fDecayPhi;
    Double_t dDPHI = dPHI - fPsi2;
    if( dDPHI < 0 ) dDPHI += TMath::TwoPi();
    if( dDPHI > TMath::Pi() ) dDPHI = TMath::TwoPi()-dDPHI;
    if(fQAlevel>1) {
      if( (dDPHI>TMath::PiOver4()) && (dDPHI<3*TMath::PiOver4()) ) FillCandidateSpy("V0SAllOP");
      else FillCandidateSpy("V0SAllIP");
    }
    FillCandidateSpy("V0SAll");
    if(!fSkipVn)
      FillDecayVn("V0SAllVn",fDecayMass,fDecayPt,fDecayPhi,fDecayEta,fDecayIDpos,fDecayIDneg);
    
    if(!AcceptCandidate()) continue;

    if(fDecayPt<fDecayStopPIDAtPt) {
      if( fSpecie==0 ) {//PID for kzero::pion+pion
        if( !PassesPIDCuts(&ieT,AliPID::kPion) ) continue; //positive track
        if( !PassesPIDCuts(&jeT,AliPID::kPion) ) continue; //negative track
      } else { //PID for lambda::proton+pion
        if(fDecayAlpha>0) {
          if( !PassesPIDCuts(&ieT,AliPID::kProton) ) continue; //positive track
	  if( !PassesPIDCuts(&jeT,AliPID::kPion) ) continue; //negative track
        } else {
          if( !PassesPIDCuts(&jeT,AliPID::kProton) ) continue; //negative track
	  if( !PassesPIDCuts(&ieT,AliPID::kPion) ) continue; //positive track
        }
      }
    }
    if(fQAlevel>1) {
      if( (dDPHI>TMath::PiOver4()) && (dDPHI<3*TMath::PiOver4()) ) FillCandidateSpy("V0SSelOP");
      else FillCandidateSpy("V0SSelIP");
    }
    FillCandidateSpy("V0SSel");
    if(!fSkipVn)
      FillDecayVn("V0SSelVn",fDecayMass,fDecayPt,fDecayPhi,fDecayEta,fDecayIDpos,fDecayIDneg);
    // ============================
    // Posting for FlowAnalysis
    if(!fPostMatched) {
      fDecayIDneg = iT->GetID();
      fDecayIDpos = jT->GetID();
      if(fUseFP) MakeTrack();
    }
    // ============================
    LoadTrack(&ieT,iT->Chi2perNDF());
    ieT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];
    FillTrackSpy("SelDau");
    LoadTrack(&jeT,jT->Chi2perNDF()); 
    jeT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];
    FillTrackSpy("SelDau");
    //===== BEGIN OF MCMATCH
    if(fReadMC) ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 1 ); // Selected event
    if(mcArray) {
      ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 2 ); // Stack found
      bool matched = false;
      bool feeddown = false;
      Int_t labelpos = iT->GetLabel();
      Int_t labelneg = jT->GetLabel();
      AliAODMCParticle *mcpos = (AliAODMCParticle*) mcArray->At( TMath::Abs(labelpos) );
      AliAODMCParticle *mcneg = (AliAODMCParticle*) mcArray->At( TMath::Abs(labelneg) );
      if( mcpos && mcneg ) {
        ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 3 ); // Daughters in stack
        Int_t pdgRecPos = mcpos->GetPdgCode();
        Int_t pdgRecNeg = mcneg->GetPdgCode();
        int pospdg=211, negpdg=211;
        int mompdg=310, fdwpdg=333;
        if(fSpecie>0) {
          mompdg=3122;
          fdwpdg=3312;
          if(fDecayAlpha>0) {
            pospdg=2212; negpdg=211;
          } else {
            negpdg=2212; pospdg=211;
          }
        }
        if( TMath::Abs(pdgRecPos)==pospdg&&TMath::Abs(pdgRecNeg)==negpdg )
          if(mcpos->GetMother()>-1)
            if( mcpos->GetMother()==mcneg->GetMother() ) {
              AliAODMCParticle *mcmot = (AliAODMCParticle*) mcArray->At( mcpos->GetMother() );
              fDecayMatchOrigin = TMath::Sqrt( mcmot->Xv()*mcmot->Xv() + mcmot->Yv()*mcmot->Yv() );
              fDecayMatchPt = mcmot->Pt();
              fDecayMatchEta = mcmot->Eta();
              fDecayMatchPhi = mcmot->Phi();
              if( TMath::Abs(mcmot->GetPdgCode())==mompdg) {
                if(mcmot->GetNDaughters()==2) {
		  ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 4 ); // Correspond to decay
                  matched=true;
                  Double_t dx = mcmot->Xv() - mcpos->Xv();
                  Double_t dy = mcmot->Yv() - mcpos->Yv();
                  fDecayMatchRadXY = TMath::Sqrt( dx*dx + dy*dy );
                }
                if(mcmot->GetMother()>-1) {
		  ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 5 ); // Decay has mother
                  AliAODMCParticle *mcfdw = (AliAODMCParticle*) mcArray->At( mcmot->GetMother() );
                  if( TMath::Abs(mcfdw->GetPdgCode())==fdwpdg)
                    feeddown=true;
                } // k0/lda have mother
              } // mother matches k0/lda
            } // both have same mother
      }
      if(matched) {
        FillCandidateSpy("Mth",true);
	if(!fSkipVn)
	  FillDecayVn("V0SMthVn",fDecayMass,fDecayPt,fDecayPhi,fDecayEta,fDecayIDpos,fDecayIDneg);
	if(fPostMatched>0) {
	  fDecayIDneg = iT->GetID();
	  fDecayIDpos = jT->GetID();
	  if(fUseFP) MakeTrack();
	}
	if(labelpos<0&&labelneg<0) {
	  FillCandidateSpy("MthNegNeg",true);
	  if(!fSkipVn)
	    FillDecayVn("V0SMthNegNegVn",fDecayMass,fDecayPt,fDecayPhi,fDecayEta,fDecayIDpos,fDecayIDneg);
	} else if(labelpos>0&&labelneg>0) {
	  if(!fSkipVn)
	    FillDecayVn("V0SMthPosPosVn",fDecayMass,fDecayPt,fDecayPhi,fDecayEta,fDecayIDpos,fDecayIDneg);
	} else if(labelpos*labelneg<0) {
	  FillCandidateSpy("MthPosNeg",true);
	  if(!fSkipVn)
	    FillDecayVn("V0SMthPosNegVn",fDecayMass,fDecayPt,fDecayPhi,fDecayEta,fDecayIDpos,fDecayIDneg);
	}
	AliAODVertex *secvtx = myV0->GetSecondaryVtx();
	Double_t possec[3],covsec[6];
	secvtx->GetXYZ(possec);
	secvtx->GetCovarianceMatrix(covsec);
	const AliESDVertex vSecVtx(possec,covsec,100.,100);
	AliESDtrack trackAtSecI( iT );
	trackAtSecI.SetTPCClusterMap( iT->GetTPCClusterMap() );
	trackAtSecI.SetTPCSharedMap( iT->GetTPCSharedMap() );
	trackAtSecI.SetTPCPointsF( iT->GetTPCNclsF() );
	trackAtSecI.PropagateToDCA(&vSecVtx, tAOD->GetMagneticField(), 100);
	fDaughterAtSecPhi = trackAtSecI.Phi();
	fDaughterAtSecEta = trackAtSecI.Eta();
	fDaughterAtSecPt = trackAtSecI.Pt();
        LoadTrack(&ieT,iT->Chi2perNDF());
	fDaughterMatchPhi=mcpos->Phi();
	fDaughterMatchEta=mcpos->Eta();
	fDaughterMatchPt=mcpos->Pt();
        ieT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
        fDaughterImpactParameterXY = ip[0];
        fDaughterImpactParameterZ = ip[1];
        FillTrackSpy("MthDau",true);
	if(labelpos<0||labelneg<0) FillTrackSpy("MthNegDau",true);
	else FillTrackSpy("MthPosDau",true);
	AliESDtrack trackAtSecJ( jT );
	trackAtSecJ.SetTPCClusterMap( jT->GetTPCClusterMap() );
	trackAtSecJ.SetTPCSharedMap( jT->GetTPCSharedMap() );
	trackAtSecJ.SetTPCPointsF( jT->GetTPCNclsF() );
	trackAtSecJ.PropagateToDCA(&vSecVtx, tAOD->GetMagneticField(), 100);
	fDaughterAtSecPhi = trackAtSecJ.Phi();
	fDaughterAtSecEta = trackAtSecJ.Eta();
	fDaughterAtSecPt = trackAtSecJ.Pt();
        LoadTrack(&jeT,jT->Chi2perNDF());
	fDaughterMatchPhi=mcneg->Phi();
	fDaughterMatchEta=mcneg->Eta();
	fDaughterMatchPt=mcneg->Pt();
        jeT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
        fDaughterImpactParameterXY = ip[0];
        fDaughterImpactParameterZ = ip[1];
        FillTrackSpy("MthDau",true);
	if(labelpos<0||labelneg<0) FillTrackSpy("MthNegDau",true);
	else FillTrackSpy("MthPosDau",true);
      } else {
        FillCandidateSpy("UnMth",false);
	if(!fSkipVn)
	  FillDecayVn("V0SUnMthVn",fDecayMass,fDecayPt,fDecayPhi,fDecayEta,fDecayIDpos,fDecayIDneg);
	if(fPostMatched<0) {
	  fDecayIDneg = iT->GetID();
	  fDecayIDpos = jT->GetID();
	  if(fUseFP) MakeTrack();
	}
        LoadTrack(&ieT,iT->Chi2perNDF());
        ieT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
        fDaughterImpactParameterXY = ip[0];
        fDaughterImpactParameterZ = ip[1];
        FillTrackSpy("UnMthDau",false);
        LoadTrack(&jeT,jT->Chi2perNDF());
        jeT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
        fDaughterImpactParameterXY = ip[0];
        fDaughterImpactParameterZ = ip[1];
        FillTrackSpy("UnMthDau",false);
      }
      if(feeddown) {
        FillCandidateSpy("MthFeedDown",true);
      }
    }
    //===== END OF MCMATCH
  }
  ((TH2D*)((TList*)fList->FindObject("V0SAll"))->FindObject("V0SADC"))->Fill( v0all,v0imw );
  if(!fSkipVn) {
    QCStoreDecayVn("V0SAllVn");
    QCStoreDecayVn("V0SSelVn");
    if(fReadMC) {
      QCStoreDecayVn("V0SMthVn");
      QCStoreDecayVn("V0SMthNegNegVn");
      QCStoreDecayVn("V0SMthPosPosVn");
      QCStoreDecayVn("V0SMthPosNegVn");
    }
  }
  return;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::PassesPIDCuts(AliESDtrack *myTrack, AliPID::EParticleType pid) {
  Bool_t pass=kTRUE;
  if(fPIDResponse) {
    fDaughterNSigmaPID = fPIDResponse->NumberOfSigmasTPC(myTrack,pid);
    if( TMath::Abs(fDaughterNSigmaPID) > fDaughterMaxNSigmaPID )
      pass = kFALSE;
  }
  return pass;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ChargeParticles(AliAODEvent *tAOD) {
  //benchmark purposes
  if(!tAOD) return;
  TClonesArray* mcArray=NULL;
  if(fReadMC) {
    mcArray = dynamic_cast<TClonesArray*>(tAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    ReadStack(mcArray);
  }
  for(int i=0; i!=tAOD->GetNumberOfTracks(); ++i) {
    AliAODTrack *t = dynamic_cast<AliAODTrack*>(tAOD->GetTrack( i ));
    if(!t) continue;
    if( !t->TestFilterBit(1) ) continue;
    fDecayMass=0.0; // using mass as nsigmas control plot
    if(fPIDResponse) { // PID
      switch(fSpecie) { // TPC PID only
      case(kPION):
        fDecayMass = fPIDResponse->NumberOfSigmasTPC(t,AliPID::kPion);
        break;
      case(kKAON):
        fDecayMass = fPIDResponse->NumberOfSigmasTPC(t,AliPID::kKaon);
        break;
      case(kPROTON):
        fDecayMass = fPIDResponse->NumberOfSigmasTPC(t,AliPID::kProton);
        break;
      }
    }
    Bool_t pass = kTRUE;
    if( TMath::Abs(fDecayMass) > 3.0 ) pass=kFALSE;
    if( t->Eta()<-0.5 || t->Eta()>+0.5 ) pass=kFALSE;
    if( t->Pt()<0.2 || t->Pt()>20.0 ) pass=kFALSE;
    AliESDtrack et( t );
    et.SetTPCClusterMap( t->GetTPCClusterMap() );
    et.SetTPCSharedMap( t->GetTPCSharedMap() );
    et.SetTPCPointsF( t->GetTPCNclsF() );
    Float_t ip[2];
    LoadTrack(&et,t->Chi2perNDF()); 
    AliAODVertex *vtx = tAOD->GetPrimaryVertex();
    Double_t pos[3];
    vtx->GetXYZ(pos);
    et.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];

    FillTrackSpy("TrkAll");
    if(!fSkipVn)
      FillTrackVn("TrkAllVn",t->Pt(),t->Phi(),t->Eta(),t->GetID());
    if(!pass) continue;
    FillTrackSpy("TrkSel");
    if(!fSkipVn)
      FillTrackVn("TrkSelVn",t->Pt(),t->Phi(),t->Eta(),t->GetID());
    if(fReadMC) {
      ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 1 ); // Selected event 
      if(mcArray) {
	((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 2 ); // Stack found
	bool matched = false;
	Int_t label = t->GetLabel();
	AliAODMCParticle *mcpar = (AliAODMCParticle*) mcArray->At( TMath::Abs(label) );
	if( mcpar ) {
	  ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 3 ); // Particle in stack
	  Int_t pdgmcpar = TMath::Abs(mcpar->GetPdgCode());
	  switch(fSpecie) {
	  case(kPION):
	    if(pdgmcpar==211) matched = true;
	    break;
	  case(kKAON):
	    if(pdgmcpar==211) matched = true;
	    break;
	  case(kPROTON):
	    if(pdgmcpar==2212) matched = true;
	    break;
	  }
	  if(!mcpar->IsPrimary()) matched = false;
	}
	if(matched) {
	  FillTrackSpy("Mth");
	  if(!fSkipVn)
	    FillTrackVn("MthVn",t->Pt(),t->Phi(),t->Eta(),t->GetID());
	  if(label<0) {
	    FillTrackSpy("MthNeg");
	    if(!fSkipVn)
	      FillTrackVn("MthNegVn",t->Pt(),t->Phi(),t->Eta(),t->GetID());
	  } else {
	    FillTrackSpy("MthPos");
	    if(!fSkipVn)
	      FillTrackVn("MthPosVn",t->Pt(),t->Phi(),t->Eta(),t->GetID());
	  }
	}
      }
    }
    if(fUseFP) {
      fDecayPt=t->Pt();
      fDecayPhi=t->Phi();
      fDecayEta=t->Eta();
      fDecayID=t->GetID();
      MakeTrack();
    }
  }
  if(!fSkipVn) {
    QCStoreTrackVn("TrkAllVn");
    QCStoreTrackVn("TrkSelVn");
    if(fReadMC) {
      QCStoreTrackVn("MthVn");
      QCStoreTrackVn("MthNegVn");
      QCStoreTrackVn("MthPosVn");
    }
  }
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ComputeChi2VZERO() {
  Double_t MeanQaQc = ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmVZEAQmVZEC"))->GetBinContent( 1 );
  Double_t MeanQaQa = ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmVZEASQUARED"))->GetBinContent( 1 );
  Double_t MeanQcQc = ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmVZECSQUARED"))->GetBinContent( 1 );
  Double_t MeanQaQt = ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmTPCQmVZEA"))->GetBinContent( 1 );
  Double_t MeanQcQt = ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmTPCQmVZEC"))->GetBinContent( 1 );
  if(!TMath::AreEqualAbs(MeanQaQt,0,1e-10)&&!TMath::AreEqualAbs(MeanQcQt,0,1e-10)&&!TMath::AreEqualAbs(MeanQaQc,0,1e-10)) {
    Double_t OneOverChiSquaredVZEA = MeanQaQa*MeanQcQt/MeanQaQc/MeanQaQt-1;
    Double_t OneOverChiSquaredVZEC = MeanQcQc*MeanQaQt/MeanQaQc/MeanQcQt-1;
    if(!TMath::AreEqualAbs(OneOverChiSquaredVZEA,0,1e-10)&&!TMath::AreEqualAbs(OneOverChiSquaredVZEC,0,1e-10)) {
      ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("ChiSquaredVZEA"))->SetBinContent( 1, 1/OneOverChiSquaredVZEA );
      ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("ChiSquaredVZEC"))->SetBinContent( 1, 1/OneOverChiSquaredVZEC );
    }
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::Terminate(Option_t *) {
  //terminate
  if(fSkipTerminate) return;
  ComputeChi2VZERO();
  if(fSkipSelection) return;
  if(fSkipVn) return;
  if(fSpecie<10) {
    ComputeDecayVn("V0SAllVn");
    ComputeDecayVn("V0SSelVn");
    if(fReadMC) {
      ComputeDecayVn("V0SMthVn");
      ComputeDecayVn("V0SMthPosPosVn");
      ComputeDecayVn("V0SMthNegNegVn");
      ComputeDecayVn("V0SMthPosNegVn");
      ComputeDecayVn("V0SUnMthVn");
    }
  } else {
    ComputeTrackVn("TrkAllVn");
    ComputeTrackVn("TrkSelVn");
    if(fReadMC) {
      ComputeTrackVn("MthVn");
      ComputeTrackVn("MthPosVn");
      ComputeTrackVn("MthNegVn");
    }
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeTrack() {
  // create track for flow tasks
  if(fCandidates->GetLast()+5>fCandidates->GetSize()) {
    fCandidates->Expand( fCandidates->GetSize()+20 );
  }
  Bool_t overwrite = kTRUE;
  AliFlowCandidateTrack *oTrack = (static_cast<AliFlowCandidateTrack*> (fCandidates->At( fCandidates->GetLast()+1 )));
  if( !oTrack ) { // creates new
    oTrack = new AliFlowCandidateTrack();
    overwrite = kFALSE;
  } else { // overwrites
    oTrack->ClearMe();
  }
  oTrack->SetMass(fDecayMass);
  oTrack->SetPt(fDecayPt);
  oTrack->SetPhi(fDecayPhi);
  oTrack->SetEta(fDecayEta);
  if(fSpecie<10) {
    oTrack->AddDaughter(fDecayIDpos);
    oTrack->AddDaughter(fDecayIDneg);
  } else {
    oTrack->SetID( fDecayID );
  }
  oTrack->SetForPOISelection(kTRUE);
  oTrack->SetForRPSelection(kFALSE);
  if(overwrite) {
    fCandidates->SetLast( fCandidates->GetLast()+1 );
  } else {
    fCandidates->AddLast(oTrack);
  }
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddCandidates() {
  // adds candidates to flow events (untaging if necessary)
  if(fDebug) printf("FlowEventTPC %d tracks | %d RFP | %d POI\n",fTPCevent->NumberOfTracks(),fTPCevent->GetNumberOfRPs(),fTPCevent->GetNumberOfPOIs());
  if(fDebug) printf("FlowEventVZE %d tracks | %d RFP | %d POI\n",fVZEevent->NumberOfTracks(),fVZEevent->GetNumberOfRPs(),fVZEevent->GetNumberOfPOIs());
  if(fDebug) printf("I received %d candidates\n",fCandidates->GetEntriesFast());
  Int_t untagged=0;
  Int_t poi=0;
  for(int iCand=0; iCand!=fCandidates->GetEntriesFast(); ++iCand ) {
    AliFlowCandidateTrack *cand = static_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
    if(!cand) continue;
    cand->SetForPOISelection(kTRUE);
    cand->SetForRPSelection(kFALSE);
    poi++;
    //if(fDebug) printf(" >Checking at candidate %d with %d daughters: mass %f\n",iCand,cand->GetNDaughters(),cand->Mass());
    if(fSpecie<10) { // DECAYS
      // untagging ===>
      if(fDaughterUnTag) {
	for(int iDau=0; iDau!=cand->GetNDaughters(); ++iDau) {
	  if(fDebug) printf("  >Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau));
	  for(int iRPs=0; iRPs!=fTPCevent->NumberOfTracks(); ++iRPs ) {
	    AliFlowTrack *iRP = static_cast<AliFlowTrack*>(fTPCevent->GetTrack( iRPs ));
	    if(!iRP) continue;
	    if(!iRP->InRPSelection()) continue;
	    if(cand->GetIDDaughter(iDau) == iRP->GetID()) {
	      if(fDebug) printf(" was in RP set");
	      ++untagged;
	      iRP->SetForRPSelection(kFALSE);
	      fTPCevent->SetNumberOfRPs( fTPCevent->GetNumberOfRPs() -1 );
	    }
	  }
	  if(fDebug) printf("\n");
	}
      }
      // <=== untagging 
      fTPCevent->InsertTrack( ((AliFlowTrack*) cand) );
    } else {  // CHARGED
      // adding only new tracks and tagging accordingly ===>
      Bool_t found=kFALSE;
      for(int iRPs=0; iRPs!=fTPCevent->NumberOfTracks(); ++iRPs ) {
        AliFlowTrack *iRP = static_cast<AliFlowTrack*>(fTPCevent->GetTrack( iRPs ));
        if(!iRP) continue;
        if(!iRP->InRPSelection()) continue;
        if(cand->GetID() == iRP->GetID()) {
          if(fDebug) printf("  >charged track (%d) was also found in RP set (adding poi tag)\n",cand->GetID());
	  iRP->SetMass( cand->Mass() );
          iRP->SetForPOISelection(kTRUE);
          found = kTRUE;
        }
      }
      if(!found) // not found adding track
        fTPCevent->InsertTrack( ((AliFlowTrack*) cand) );
    }
    fVZEevent->InsertTrack( ((AliFlowTrack*) cand) );
  } //END OF LOOP
  fTPCevent->SetNumberOfPOIs( poi );
  fVZEevent->SetNumberOfPOIs( poi );
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("POI"))->Fill( poi );
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("UNTAG"))->Fill( untagged );
  if(fDebug) printf("FlowEventTPC %d tracks | %d RFP | %d POI\n",fTPCevent->NumberOfTracks(),fTPCevent->GetNumberOfRPs(),fTPCevent->GetNumberOfPOIs());
  if(fDebug) printf("FlowEventVZE %d tracks | %d RFP | %d POI\n",fVZEevent->NumberOfTracks(),fVZEevent->GetNumberOfRPs(),fVZEevent->GetNumberOfPOIs());
}
//=======================================================================
void AliAnalysisTaskFlowStrange::PushBackFlowTrack(AliFlowEvent *flowevent, Double_t pt, Double_t phi, Double_t eta, Double_t w, Int_t id) {
  AliFlowTrack rfp;
  rfp.SetPt(pt);
  rfp.SetPhi(phi);
  rfp.SetEta(eta);
  rfp.SetWeight(w);
  rfp.SetForRPSelection(kTRUE);
  rfp.SetForPOISelection(kFALSE);
  rfp.SetMass(-999);
  rfp.SetID(id);
  flowevent->InsertTrack( &rfp );
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::IsAtTPCEdge(Double_t phi,Double_t pt,Int_t charge,Double_t b) {
  // Origin: Alex Dobrin
  // Implemented by Carlos Perez
  TF1 cutLo("cutLo", "-0.01/x+pi/18.0-0.015", 0, 100);
  TF1 cutHi("cutHi", "0.55/x/x+pi/18.0+0.03", 0, 100);
  Double_t phimod = phi;
  if(b<0) phimod = TMath::TwoPi()-phimod;  //for negatve polarity field
  if(charge<0) phimod = TMath::TwoPi()-phimod; //for negatve charge
  phimod += TMath::Pi()/18.0;
  phimod = fmod(phimod, TMath::Pi()/9.0);
  if( phimod<cutHi.Eval(pt) && phimod>cutLo.Eval(pt) )
    return kTRUE;

  return kFALSE;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQVectors() {
  //computes event plane and updates fPsi2
  //if there is a problem fPsi->-1
  fPsi2=-1;
  fVZEWarning=kFALSE;
  //=>loading event
  MakeQVZE(InputEvent());
  MakeQTPC(InputEvent());
  if(fUseFP&&fReadMC) {
    fVZEevent->SetMCReactionPlaneAngle( fMCEP );    
    fTPCevent->SetMCReactionPlaneAngle( fMCEP );    
  }
  if(fDebug) {
    printf("**::MakeQVectors()");
    printf("  fQVZEACos %.16f | fQVZEASin %.16f || fQVZEA %.3f | fQVZEC %.3f \n",fQVZEACos, fQVZEASin, fQVZEA, fQVZEC);
    printf("  nQTPA_nTracks %d | fQTPC_nTracks %d || fQTPCA %.3f | fQTPCC %.3f \n",fQTPCA_nTracks, fQTPCC_nTracks, fQTPCA, fQTPCC);
    printf("  fQTPCACos %.16f | fQTPCASin %.16f || fQTPC2hCos %.16f | fQTPC2hSin %.16f \n",fQTPCACos, fQTPCASin, fQTPC2hCos, fQTPC2hSin);
   }
  FillMakeQSpy();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillMakeQSpy() {
  //=>computing psi
  //VZERO
  Double_t qvzecos,qvzesin,psivzea,psivzec,psivze,qvze, qvzea, qvzec;
  psivzea = ( TMath::Pi()+TMath::ATan2(-fQVZEASin,-fQVZEACos) )/fHarmonic;
  psivzec = ( TMath::Pi()+TMath::ATan2(-fQVZECSin,-fQVZECCos) )/fHarmonic;
  qvzecos = fQVZEACos + fQVZECCos;
  qvzesin = fQVZEASin + fQVZECSin;
  qvzea = fQVZEA;
  qvzec = fQVZEC;
  qvze = fQVZEA + fQVZEC;
  psivze = ( TMath::Pi()+TMath::ATan2(-qvzesin,-qvzecos) )/fHarmonic;
  //TPC
  Double_t qtpccos,qtpcsin,psitpca,psitpcc,psitpc,qtpc;
  psitpca = ( TMath::Pi()+TMath::ATan2(-fQTPCASin,-fQTPCACos) )/fHarmonic;
  psitpcc = ( TMath::Pi()+TMath::ATan2(-fQTPCCSin,-fQTPCCCos) )/fHarmonic;
  qtpccos = fQTPCACos + fQTPCCCos;
  qtpcsin = fQTPCASin + fQTPCCSin;
  qtpc = fQTPCA + fQTPCC;
  psitpc = ( TMath::Pi()+TMath::ATan2(-qtpcsin,-qtpccos) )/fHarmonic;
  //=>does the event clear?
  switch(fWhichPsi) {
  case(1): //VZERO
    if(fVZEWarning) return;
    fPsi2 = psivze;
    break;
  case(2): //TPC
    if(fQTPCA<2||fQTPCC<2) return;
    fPsi2 = psitpc;
    break;
  }
  //computing physical Qm vectors
  Double_t vzec_qmcos = fQVZECCos/fQVZEC;
  Double_t vzec_qmsin = fQVZECSin/fQVZEC;
  Double_t vzec_qmnor = TMath::Sqrt( vzec_qmcos*vzec_qmcos + vzec_qmsin*vzec_qmsin );
  Double_t vzea_qmcos = fQVZEACos/fQVZEA;
  Double_t vzea_qmsin = fQVZEASin/fQVZEA;
  Double_t vzea_qmnor = TMath::Sqrt( vzea_qmcos*vzea_qmcos + vzea_qmsin*vzea_qmsin );
  Double_t vze_qmcos = qvzecos/qvze;
  Double_t vze_qmsin = qvzesin/qvze;
  Double_t vze_qmnor = TMath::Sqrt( vze_qmcos*vze_qmcos + vze_qmsin*vze_qmsin );
  Double_t tpcc_qmcos = fQTPCCCos/fQTPCC;
  Double_t tpcc_qmsin = fQTPCCSin/fQTPCC;
  Double_t tpcc_qmnor = TMath::Sqrt( tpcc_qmcos*tpcc_qmcos + tpcc_qmsin*tpcc_qmsin );
  Double_t tpca_qmcos = fQTPCACos/fQTPCA;
  Double_t tpca_qmsin = fQTPCASin/fQTPCA;
  Double_t tpca_qmnor = TMath::Sqrt( tpca_qmcos*tpca_qmcos + tpca_qmsin*tpca_qmsin );
  Double_t tpc_qmcos = qtpccos/qtpc;
  Double_t tpc_qmsin = qtpcsin/qtpc;
  Double_t tpc_qmnor = TMath::Sqrt( tpc_qmcos*tpc_qmcos + tpc_qmsin*tpc_qmsin );
  //=>great! recording
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEPSI"))->Fill( psivze );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEPSIA"))->Fill( psivzea );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEPSIC"))->Fill( psivzec );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("RFPVZE"))->Fill( qvze );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmVZEA"))->Fill( vzea_qmnor*TMath::Sqrt(qvzea) );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmVZEC"))->Fill( vzec_qmnor*TMath::Sqrt(qvzec) );
  //------
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCPSI"))->Fill( psitpc );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCPSIA"))->Fill( psitpca );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCPSIC"))->Fill( psitpcc );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("RFPTPC"))->Fill( qtpc );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmTPC"))->Fill( tpc_qmnor*TMath::Sqrt(qtpc) );
  //------
  ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("PSI_TPCAVZEC"))->Fill( psitpca, psivzec );
  ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("PSI_TPCCVZEA"))->Fill( psitpcc, psivzea );
  ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("PSI_TPCVZE"))->Fill( psitpc, psivze );

  if(fReadMC) {
    ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("PSIMCDIFFTPC"))->Fill( psitpc-fMCEP );
    ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("PSIMCDIFFTPCA"))->Fill( psitpca-fMCEP );
    ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("PSIMCDIFFTPCC"))->Fill( psitpcc-fMCEP );
    ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("PSIMCDIFFVZE"))->Fill( psivze-fMCEP );
    ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("PSIMCDIFFVZEA"))->Fill( psivzea-fMCEP );
    ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("PSIMCDIFFVZEC"))->Fill( psivzec-fMCEP );
  }
  //------
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQm"))->Fill( 1., tpcc_qmsin, tpcc_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQm"))->Fill( 2., tpcc_qmcos, tpcc_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQm"))->Fill( 3., tpca_qmsin, tpca_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQm"))->Fill( 4., tpca_qmcos, tpca_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQm"))->Fill( 5., tpc_qmsin, tpc_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQm"))->Fill( 6., tpc_qmcos, tpc_qmnor );
  //------
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQm"))->Fill( 1., vzec_qmsin, vzec_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQm"))->Fill( 2., vzec_qmcos, vzec_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQm"))->Fill( 3., vzea_qmsin, vzea_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQm"))->Fill( 4., vzea_qmcos, vzea_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQm"))->Fill( 5., vze_qmsin, vze_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQm"))->Fill( 6., vze_qmcos, vze_qmnor );
  //------
  Double_t vzeqaqc = vzec_qmcos*vzea_qmcos + vzec_qmsin*vzea_qmsin;
  Double_t vzeqatpcq = vzea_qmcos*tpc_qmcos + vzea_qmsin*tpc_qmsin;
  Double_t vzeqctpcq = vzec_qmcos*tpc_qmcos + vzec_qmsin*tpc_qmsin;
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmVZEAQmVZEC"))->Fill( 1., vzeqaqc, vze_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmVZEASQUARED"))->Fill( 1., vzea_qmnor*vzea_qmnor, vze_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmVZECSQUARED"))->Fill( 1., vzec_qmnor*vzec_qmnor, vze_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmTPCQmVZEA"))->Fill( 1., vzeqatpcq, vze_qmnor );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmTPCQmVZEC"))->Fill( 1., vzeqctpcq, vze_qmnor );
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQVZE(AliVEvent *tevent) {
  //=>cleaning
  if(fUseFP) fVZEevent->ClearFast(); // flowpackage
  //=>computing
  fQVZEACos=fQVZEASin=fQVZEA=fQVZECCos=fQVZECSin=fQVZEC=0;
  Int_t rfp=0;
  Double_t eta, phi, w;
  //v0c -> qa
  for(int id=fVZECa*8;id!=8+fVZECb*8;++id) {
    eta = -3.45+0.5*(id/8);
    phi = TMath::PiOver4()*(0.5+id%8);
    w = tevent->GetVZEROEqMultiplicity(id);
    if(w<3) fVZEWarning=kTRUE;
    w *= fVZEextW[id];
    fQVZECCos += w*TMath::Cos(fHarmonic*phi);
    fQVZECSin += w*TMath::Sin(fHarmonic*phi);
    fQVZEC += w;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEAllPhiEta"))->Fill( phi, eta, w );
    rfp++;
    if(fUseFP) PushBackFlowTrack(fVZEevent,0,phi,eta,w,0); // flowpackage
  }
  //v0a -> qb
  for(int id=32+fVZEAa*8;id!=40+fVZEAb*8;++id) {
    eta = +4.8-0.6*((id/8)-4);
    phi = TMath::PiOver4()*(0.5+id%8);
    w = tevent->GetVZEROEqMultiplicity(id);
    if(w<3) fVZEWarning=kTRUE;
    w *= fVZEextW[id];
    fQVZEACos += w*TMath::Cos(fHarmonic*phi);
    fQVZEASin += w*TMath::Sin(fHarmonic*phi);
    fQVZEA += w;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEAllPhiEta"))->Fill( phi, eta, w );
    rfp++;
    if(fUseFP) PushBackFlowTrack(fVZEevent,0,phi,eta,w,0); // flowpackage
  }
  if(fUseFP) { // flowpackage
    fVZEevent->SetNumberOfRPs(rfp);
    if(fDebug>0) printf("Inserted tracks in FlowEventVZE %d ==> %.1f\n",fVZEevent->NumberOfTracks(),fQVZEC+fQVZEA);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddTPCRFPSpy(TList *me) {
  TH1D *tH1D;
  tH1D = new TH1D("PT",   "PT",           50,0,5);     me->Add(tH1D);
  tH1D = new TH1D("PHI",  "PHI", 90,0,TMath::TwoPi()); me->Add(tH1D);
  tH1D = new TH1D("ETA",  "ETA",          40,-1,+1);   me->Add(tH1D);
  tH1D = new TH1D("TPCS", "TPC Signal",   100,0,500);  me->Add(tH1D);
  tH1D = new TH1D("IPXY", "IPXY",         100,-2,+2);  me->Add(tH1D);
  tH1D = new TH1D("IPZ",  "IPZ",          100,-2,+2);  me->Add(tH1D);
  // TPC
  tH1D = new TH1D("TPCNCLS", "NCLS", 170,-0.5,+169.5);   me->Add(tH1D);
  tH1D = new TH1D("TPCSHCL", "NSCLS / NCLS", 100,0,1);   me->Add(tH1D);
  tH1D = new TH1D("TPCFICL", "NCLS1I / NCLS",100,0,1);   me->Add(tH1D);
  tH1D = new TH1D("TPCXRNF", "XROW / NFCLS", 100,0,1.5); me->Add(tH1D);
  tH1D = new TH1D("TPCRCHI", "CHI2 / NCLS",  50,0,5);    me->Add(tH1D);
  // ITS
  tH1D = new TH1D("ITSNCLS", "NCLS",   7,-0.5,+6.5); me->Add(tH1D);
  tH1D = new TH1D("ITSRCHI", "CHI2 / NCLS", 50,0,5); me->Add(tH1D);

}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::PassesRFPTPCCuts(AliESDtrack *track, Double_t aodchi2cls, Float_t aodipxy, Float_t aodipz) {
  if(track->GetKinkIndex(0)>0) return kFALSE;
  if( (track->GetStatus()&AliESDtrack::kTPCrefit)==0 ) return kFALSE;
  Double_t pt = track->Pt();
  Double_t phi = track->Phi();
  Double_t eta = track->Eta();
  Double_t tpcs = track->GetTPCsignal();
  Float_t ipxy, ipz;
  track->GetImpactParameters(ipxy,ipz);
  Int_t cls = track->GetTPCclusters(0);
  Double_t xrows, findcls, chi2;
  findcls = track->GetTPCNclsF();
  xrows = track->GetTPCCrossedRows();
  chi2 = track->GetTPCchi2();
  Double_t rchi2 = chi2/cls;
  if(!fReadESD) {
    rchi2 = aodchi2cls;
    ipxy = aodipxy;
    ipz = aodipz;
  }
  Double_t xrnfcls = xrows/findcls;
  Double_t scls, cls1i, itschi2;
  Int_t itscls;
  cls1i = track->GetTPCNclsIter1();
  scls = track->GetTPCnclsS();
  itscls = track->GetITSclusters(0);
  itschi2 = track->GetITSchi2();
  Double_t shcl = scls/cls;
  Double_t ficl = cls1i/cls;
  Double_t itsrchi2 = itscls/itschi2;
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("PT"))->Fill( pt );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("PHI"))->Fill( phi );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("ETA"))->Fill( eta );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCS"))->Fill( tpcs );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("IPXY"))->Fill( ipxy );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("IPZ"))->Fill( ipz );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCNCLS"))->Fill( cls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCSHCL"))->Fill( shcl );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCFICL"))->Fill( ficl );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCXRNF"))->Fill( xrnfcls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCRCHI"))->Fill( rchi2 );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("ITSNCLS"))->Fill( itscls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("ITSRCHI"))->Fill( itsrchi2 );
  if(pt<fRFPminPt) return kFALSE; //0.2
  if(pt>fRFPmaxPt) return kFALSE; //5.0
  if(eta<fRFPCminEta||(eta>fRFPCmaxEta&&eta<fRFPAminEta)||eta>fRFPAmaxEta) return kFALSE; // -0.8 0.0 0.0 +0.8
  if(tpcs<fRFPTPCsignal) return kFALSE; //10.0
  if( TMath::Sqrt(ipxy*ipxy/fRFPmaxIPxy/fRFPmaxIPxy+ipz*ipz/fRFPmaxIPz/fRFPmaxIPz)>1 ) return kFALSE; //2.4 3.2
  if(cls<fRFPTPCncls) return kFALSE; //70
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("PT"))->Fill( pt );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("PHI"))->Fill( phi );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("ETA"))->Fill( eta );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCS"))->Fill( tpcs );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("IPXY"))->Fill( ipxy );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("IPZ"))->Fill( ipz );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCNCLS"))->Fill( cls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCSHCL"))->Fill( shcl );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCFICL"))->Fill( ficl );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCXRNF"))->Fill( xrnfcls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCRCHI"))->Fill( rchi2 );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("ITSNCLS"))->Fill( itscls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("ITSRCHI"))->Fill( itsrchi2 );
  return kTRUE;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQTPC(AliVEvent *tevent) {
  AliESDEvent *tESD = (AliESDEvent*) (tevent);
  AliAODEvent *tAOD = (AliAODEvent*) (tevent);
  if(fReadESD) {
    if(!tESD) return;
    MakeQTPC(tESD);
  } else {
    if(!tAOD) return;
    MakeQTPC(tAOD);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQTPC(AliAODEvent *tAOD) {
  //=>cleaning
  if(fUseFP) fTPCevent->ClearFast(); // flowpackage
  fQTPCACos=fQTPCASin=fQTPCA=fQTPC2hCos=0;
  fQTPCCCos=fQTPCCSin=fQTPCC=fQTPC2hSin=0;
  fQTPCA_nTracks = 0;
  fQTPCC_nTracks = 0;
  Int_t rfp=0;
  Double_t eta, phi, w;
  //=>aod stuff
  AliAODVertex *vtx = tAOD->GetPrimaryVertex();
  Double_t pos[3],cov[6];
  vtx->GetXYZ(pos);
  vtx->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  AliAODTrack *track;
  //=>looping
  Int_t rawN = tAOD->GetNumberOfTracks();
  for(Int_t id=0; id!=rawN; ++id) {
    track = dynamic_cast<AliAODTrack*>(tAOD->GetTrack(id));
    if(!track) {
        AliFatal("Not a standard AOD");
        return; // shut up coverity
    }
    //=>cuts
    if(!track->TestFilterBit(fRFPFilterBit)) continue;
    if( fExcludeTPCEdges )
      if( IsAtTPCEdge( track->Phi(), track->Pt(), track->Charge(), tAOD->GetMagneticField() ) )	continue;
    AliESDtrack etrack( track );
    etrack.SetTPCClusterMap( track->GetTPCClusterMap() );
    etrack.SetTPCSharedMap( track->GetTPCSharedMap() );
    etrack.SetTPCPointsF( track->GetTPCNclsF() );
    Float_t ip[2];
    etrack.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    if(!PassesRFPTPCCuts(&etrack,track->Chi2perNDF(),ip[0],ip[1])) continue;
    //=>collecting info
    phi = track->Phi();
    eta = track->Eta();
    w = 1;
    if(eta<0) {
      fQTPCCCos += w*TMath::Cos(fHarmonic*phi);
      fQTPCCSin += w*TMath::Sin(fHarmonic*phi);
      fQTPCC += w;
      fQTPCC_fID[fQTPCC_nTracks++] = track->GetID();
    } else {
      fQTPCACos += w*TMath::Cos(fHarmonic*phi);
      fQTPCASin += w*TMath::Sin(fHarmonic*phi);
      fQTPCA += w;
      fQTPCA_fID[fQTPCA_nTracks++] = track->GetID();
    }
    fQTPC2hCos += w*TMath::Cos(2.0*fHarmonic*phi);
    fQTPC2hSin += w*TMath::Sin(2.0*fHarmonic*phi);
    rfp++;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCAllPhiEta"))->Fill( phi, eta, w );
    if(fUseFP) PushBackFlowTrack(fTPCevent,track->Pt(),phi,eta,w,track->GetID()); // flow package
  }
  if(fUseFP) {
    fTPCevent->SetNumberOfRPs(rfp);
    if(fDebug) printf("Inserted tracks in FlowEventTPC %d ==> %.1f\n",fTPCevent->NumberOfTracks(),fQTPCA+fQTPCC);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQTPC(AliESDEvent *tESD) {
  //=>cleaning
  if(fUseFP) fTPCevent->ClearFast(); // flow package
  fQTPCACos=fQTPCASin=fQTPCA=0;
  fQTPCCCos=fQTPCCSin=fQTPCC=0;
  fQTPCA_nTracks = 0;
  fQTPCC_nTracks = 0;
  Int_t rfp=0;
  Double_t eta, phi, w;
  //=>looping
  AliESDtrack *track;
  Int_t rawN = tESD->GetNumberOfTracks();
  for(Int_t id=0; id!=rawN; ++id) {
    track = tESD->GetTrack(id);
    //=>cuts
    if( fExcludeTPCEdges )
      if( IsAtTPCEdge( track->Phi(), track->Pt(), track->Charge(), tESD->GetMagneticField() ) )	continue;
    if(!PassesFilterBit(track)) continue;
    if(!PassesRFPTPCCuts(track)) continue;
    //=>collecting info
    phi = track->Phi();
    eta = track->Eta();
    w = 1;
    if(eta<0) {
      fQTPCCCos += w*TMath::Cos(fHarmonic*phi);
      fQTPCCSin += w*TMath::Sin(fHarmonic*phi);
      fQTPCC += w;
      fQTPCC_fID[fQTPCC_nTracks++] = track->GetID();
    } else {
      fQTPCACos += w*TMath::Cos(fHarmonic*phi);
      fQTPCASin += w*TMath::Sin(fHarmonic*phi);
      fQTPCA += w;
      fQTPCA_fID[fQTPCA_nTracks++] = track->GetID();
    }
    rfp++;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCAllPhiEta"))->Fill( phi, eta, w );
    if(fUseFP) PushBackFlowTrack(fTPCevent,track->Pt(),phi,eta,w,track->GetID()); // flowpackage
  }
  if(fUseFP) {
    fTPCevent->SetNumberOfRPs(rfp);
    if(fDebug) printf("Inserted tracks in FlowEventTPC %d ==> %.1f\n",fTPCevent->NumberOfTracks(),fQTPCA+fQTPCC);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddMCParticleSpy(TList *me) {
  TH1D *tH1D;
  TH2D *tH2D;
  TProfile *tPro;
  tH1D = new TH1D("Pt",   "Pt",   fPtBins,fPtBinEdge);  me->Add(tH1D);
  tH1D = new TH1D("Phi",  "Phi",  100,0,TMath::TwoPi()); me->Add(tH1D);
  tH1D = new TH1D("Eta",  "Eta",  100,-1,+1);   me->Add(tH1D);
  tH1D = new TH1D("Y",    "Y",    100,-1,+1);   me->Add(tH1D);
  tH1D = new TH1D("Rad2", "Rad2", 1000,0,+100); me->Add(tH1D);
  tH2D = new TH2D("Dphi", "phi-MCEP;pt;dphi",fPtBins,fPtBinEdge, 72,0,TMath::Pi()); me->Add(tH2D);
  tPro = new TProfile("Cos2dphi","Cos2dphi",fPtBins,fPtBinEdge); me->Add(tPro);
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillMCParticleSpy(TString listName, AliAODMCParticle *p) {
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Pt" ))->Fill( p->Pt() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Eta" ))->Fill( p->Eta() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Y" ))->Fill( p->Y() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Phi" ))->Fill( p->Phi() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Rad2" ))->Fill( TMath::Sqrt( p->Xv()*p->Xv() +
												 p->Yv()*p->Yv() ) );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Dphi" ))->Fill( p->Pt(), GetMCDPHI(p->Phi()) );
  ((TProfile*)((TList*)fList->FindObject(listName.Data()))->FindObject("Cos2dphi" ))->Fill( p->Pt(), TMath::Cos( 2*GetMCDPHI(p->Phi()) ), 1 );
  return;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::GetMCDPHI(Double_t phi) {
  Double_t dDPHI = phi - fMCEP;
  //if( dDPHI < 0 ) dDPHI += TMath::TwoPi();
  //if( dDPHI > TMath::Pi() ) dDPHI = TMath::TwoPi()-dDPHI;
  return dDPHI;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillMCParticleSpy(TString listName, TParticle *p) {
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Pt" ))->Fill( p->Pt() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Eta" ))->Fill( p->Eta() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Phi" ))->Fill( p->Phi() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Rad2" ))->Fill( TMath::Sqrt( p->Vx()*p->Vx() +
												 p->Vy()*p->Vy() ) );
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddCandidatesSpy(TList *me,Bool_t res) {
  TH1D *tH1D;
  TH2D *tH2D;
  TProfile *tPro;
  TProfile2D *tPro2;
  tH2D = new TH2D("PhiEta",  "PhiEta;Phi;Eta", 100,0,TMath::TwoPi(),100,-1,+1);  me->Add(tH2D);
  tH2D = new TH2D("PtRAP",   "PtRAP;Pt;Y",     fPtBins,fPtBinEdge,100,-1,+1);  me->Add(tH2D);
  tH2D = new TH2D("PtDCA",   "PtDCA;Pt;DCA",   fPtBins,fPtBinEdge,100,0,1.5);  me->Add(tH2D);
  tH2D = new TH2D("PtCTP",   "PtCTP;Pt;CTP",   fPtBins,fPtBinEdge,100,0.95,+1);me->Add(tH2D);
  tH2D = new TH2D("PtD0D0",  "PtD0D0;Pt;D0D0", fPtBins,fPtBinEdge,100,-5,+5);  me->Add(tH2D);
  tH2D = new TH2D("PtRad2",  "PtRad2;Pt;RadXY",fPtBins,fPtBinEdge,100,0,+50);  me->Add(tH2D);
  tH2D = new TH2D("PtDL",    "PtDL;Pt;DL",     fPtBins,fPtBinEdge,100,0,+50);  me->Add(tH2D);
  tH2D = new TH2D("PtDLlab", "PtDL;Pt;DLlab",  fPtBins,fPtBinEdge,100,0,+100); me->Add(tH2D);
  tH2D = new TH2D("PtMASS",  "PtMASS;Pt;MASS", fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); me->Add(tH2D);
  tH2D = new TH2D("APPOS",   "APPOS;alphaPOS;QtPOS",100,-2,+2,100,0,0.3);     me->Add(tH2D);
  tH2D = new TH2D("D0PD0N",  "D0PD0N;D0P;D0N",      200,-10,+10,200,-10,+10); me->Add(tH2D);
  tH2D = new TH2D("XPOSXNEG","XPOSXNEG;XPOS;XNEG",  200,-50,+50,200,-50,+50); me->Add(tH2D);
  if(fReadMC) {
    if(res) {
      tH1D = new TH1D("MCOrigin","MCOrigin;Rad2",   1000,0,50); me->Add(tH1D);
      tH2D = new TH2D("PHIRes","PHIRes;PHI;MC-DAT", 72,0,TMath::TwoPi(),100,-0.12,+0.12); me->Add(tH2D);
      tH2D = new TH2D("ETARes","ETARes;ETA;MC-DAT", 16,-0.8,+0.8,100,-0.2,+0.2); me->Add(tH2D);
      tH2D = new TH2D("PTRes", "PTRes;Pt;MC-DAT",   fPtBins,fPtBinEdge,100,-0.4,+0.4); me->Add(tH2D);
      tH2D = new TH2D("RXYRes","RXYRes;RXY;MC-DAT", 100,0,50,100,-4.0,+4.0); me->Add(tH2D);
    }
    tH2D = new TH2D("PTDPHIMC","PtDPHIMC;Pt;PHI-MCEP",fPtBins,fPtBinEdge,72,0,TMath::Pi()); me->Add(tH2D);
    tPro = new TProfile("Cos2dphiMC",  "Cos2dphiMC",fPtBins,fPtBinEdge); me->Add(tPro);
    tPro2=new TProfile2D("C2DPHIMCMASS","C2DPHIMCMASS",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); me->Add(tPro2);
  }
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillCandidateSpy(TString listName, Bool_t fillRes) {
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PhiEta"))->Fill( fDecayPhi, fDecayEta );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtRAP" ))->Fill( fDecayPt, fDecayRapidity );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtDCA" ))->Fill( fDecayPt, fDecayDCAdaughters );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtCTP" ))->Fill( fDecayPt, fDecayCosinePointingAngleXY );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtD0D0"))->Fill( fDecayPt, fDecayProductIPXY );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtRad2"))->Fill( fDecayPt, fDecayRadXY );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtDL"  ))->Fill( fDecayPt, fDecayDecayLength );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtDLlab"))->Fill( fDecayPt, fDecayDecayLengthLab );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtMASS"))->Fill( fDecayPt, fDecayMass );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("APPOS" ))->Fill( fDecayAlpha, fDecayQt );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("D0PD0N"))->Fill( fDecayIPpos, fDecayIPneg );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("XPOSXNEG"))->Fill( fDecayXpos, fDecayXneg );
  if(fReadMC) {
    ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PTDPHIMC" ))->Fill( fDecayPt, GetMCDPHI( fDecayPhi ) );
    ((TProfile*)((TList*)fList->FindObject(listName.Data()))->FindObject("Cos2dphiMC" ))->Fill( fDecayPt, TMath::Cos( 2*GetMCDPHI(fDecayPhi) ), 1 );
    ((TProfile2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("C2DPHIMCMASS" ))->Fill(fDecayPt,fDecayMass,TMath::Cos(2*GetMCDPHI(fDecayPhi)), 1 );
    if(fillRes) {
      ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("MCOrigin"))->Fill( fDecayMatchOrigin );
      ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PHIRes"))->Fill( fDecayPhi, fDecayMatchPhi-fDecayPhi );
      ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("ETARes"))->Fill( fDecayEta, fDecayMatchEta-fDecayEta );
      ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PTRes"))->Fill( fDecayPt, fDecayMatchPt-fDecayPt );
      ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("RXYRes"))->Fill( fDecayRadXY, fDecayMatchRadXY-fDecayRadXY );
    }
  }
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptCandidate() {
  if(fDecayEta<fDecayMinEta) return kFALSE;
  if(fDecayEta>fDecayMaxEta) return kFALSE;
  if(fDecayPt<fDecayMinPt) return kFALSE;
  if(fDecayProductIPXY>fDecayMaxProductIPXY) return kFALSE;
  if(fDecayDCAdaughters>fDecayMaxDCAdaughters) return kFALSE;
  if(fDecayCosinePointingAngleXY<fDecayMinCosinePointingAngleXY) return kFALSE;
  if(fDecayRadXY<fDecayMinRadXY) return kFALSE;
  if(TMath::Abs(fDecayRapidity)>fDecayMaxRapidity) return kFALSE;
  if(fSpecie==0) {
    if(fDecayAPCutPie) {
      if(fDecayQt/TMath::Abs(fDecayAlpha)<fDecayMinQt) return kFALSE;
    } else {
      if(fDecayQt<fDecayMinQt) return kFALSE;
    }
    if(fDecayDecayLength>fDecayMaxDecayLength*2.6842) return kFALSE;
  } else {
    if(fDecayDecayLength>fDecayMaxDecayLength*7.89) return kFALSE;
  }
  if(fSpecie==1) if(fDecayAlpha>0) return kFALSE;
  if(fSpecie==2) if(fDecayAlpha<0) return kFALSE;
  return kTRUE;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddTrackSpy(TList *me,Bool_t res) {
  TH2D *tH2D;
  tH2D = new TH2D("PHIETA",       "PHIETA;PHI;ETA",       100,0,TMath::TwoPi(),100,-2,2); me->Add(tH2D);
  tH2D = new TH2D("PTTRACKDECAY", "PTTRACKDECAY;PT;PT",   100,0,10,fPtBins,fPtBinEdge); me->Add(tH2D);
  tH2D = new TH2D("IPXYIPZ",      "IPXYIPZ;IPXY;IPZ",     100,-10,+10,100,-10,+10); me->Add(tH2D);
  tH2D = new TH2D("PTTPCNCLS",    "PTTPCNCLS;PT;NCLS",    fPtBins,fPtBinEdge,170,0,170);  me->Add(tH2D);
  tH2D = new TH2D("PTITSNCLS",    "PTITSNCLS;PT;NCLS",    fPtBins,fPtBinEdge,7,-0.5,6.5); me->Add(tH2D);
  tH2D = new TH2D("PTITSLAY",     "PTITSLAY;PT;ITSLAYER", fPtBins,fPtBinEdge,6,-0.5,+5.5);me->Add(tH2D);
  tH2D = new TH2D("PTITSTPCrefit","PTITSTPCrefit;PT",     fPtBins,fPtBinEdge,2,-0.5,+1.5);me->Add(tH2D);
  tH2D->GetYaxis()->SetBinLabel(1,"ITS refit");
  tH2D->GetYaxis()->SetBinLabel(2,"TPC refit");
  tH2D = new TH2D("POSTPCNCLCHI2","POSTPCNCLCHI2;NCLS;CHI2/NCLS", 170,0,170,100,0,8);   me->Add(tH2D);
  tH2D = new TH2D("POSTPCNFCLNXR","POSTPCNFCLNXR;NFCLS;NXR",      170,0,170,170,0,170); me->Add(tH2D);
  tH2D = new TH2D("POSTPCNCLNFCL","POSTPCNCLNFCL;NCLS;NFCLS",     170,0,170,170,0,170); me->Add(tH2D);
  tH2D = new TH2D("POSTPCNCLNSCL","POSTPCNCLNSCL;NCLS;NSCLS",     170,0,170,170,0,170); me->Add(tH2D);
  tH2D = new TH2D("NEGTPCNCLCHI2","NEGTPCNCLCHI2;NCLS;CHI2/NCLS", 170,0,170,100,0,8);   me->Add(tH2D);
  tH2D = new TH2D("NEGTPCNFCLNXR","NEGTPCNFCLNXR;NFCLS;NXR",      170,0,170,170,0,170); me->Add(tH2D);
  tH2D = new TH2D("NEGTPCNCLNFCL","NEGTPCNCLNFCL;NCLS;NFCLS",     170,0,170,170,0,170); me->Add(tH2D);
  tH2D = new TH2D("NEGTPCNCLNSCL","NEGTPCNCLNSCL;NCLS;NSCLS",     170,0,170,170,0,170); me->Add(tH2D);
  if(fReadMC) {
    TProfile *tPro;
    tPro = new TProfile("COSNDPHIMC","COSNDPHIMC",fPtBins,fPtBinEdge); me->Add(tPro);
  }
  if(res) {
    tH2D = new TH2D("PHIRes", "PHIRes;PHI;MC-DAT", 72,0,TMath::TwoPi(),100,-0.12,+0.12); me->Add(tH2D);
    tH2D = new TH2D("ETARes", "ETARes;ETA;MC-DAT", 16,-0.8,+0.8,100,-0.2,+0.2); me->Add(tH2D);
    tH2D = new TH2D("PTRes",  "PTRes;Pt;MC-DAT", fPtBins,fPtBinEdge,100,-0.4,+0.4); me->Add(tH2D);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillTrackSpy(TString listName,Bool_t res) {
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PHIETA" ))->Fill( fDaughterPhi, fDaughterEta );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTTRACKDECAY" ))->Fill( fDaughterPt, fDecayPt );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "IPXYIPZ" ))->Fill( fDaughterImpactParameterXY, fDaughterImpactParameterZ );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTTPCNCLS" ))->Fill( fDaughterPt, fDaughterNClsTPC );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSNCLS" ))->Fill( fDaughterPt, fDaughterNClsITS );
  if( TESTBIT(fDaughterITScm,0) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 0 );
  if( TESTBIT(fDaughterITScm,1) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 1 );
  if( TESTBIT(fDaughterITScm,2) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 2 );
  if( TESTBIT(fDaughterITScm,3) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 3 );
  if( TESTBIT(fDaughterITScm,4) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 4 );
  if( TESTBIT(fDaughterITScm,5) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 5 );
  if( (fDaughterStatus&AliESDtrack::kITSrefit) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSTPCrefit" ))->Fill( fDaughterPt, 0 );
  if( (fDaughterStatus&AliESDtrack::kTPCrefit) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSTPCrefit" ))->Fill( fDaughterPt, 1 );
  TString ch="NEG";
  if(fDaughterCharge>0) ch="POS";
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( Form("%sTPCNCLCHI2",ch.Data()) ))->Fill( fDaughterNClsTPC, fDaughterChi2PerNClsTPC );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( Form("%sTPCNFCLNXR",ch.Data()) ))->Fill( fDaughterNFClsTPC, fDaughterXRows );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( Form("%sTPCNCLNFCL",ch.Data()) ))->Fill( fDaughterNClsTPC, fDaughterNFClsTPC );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( Form("%sTPCNCLNSCL",ch.Data()) ))->Fill( fDaughterNClsTPC, fDaughterNSClsTPC );
  if(fReadMC) {
    Double_t cosn = TMath::Cos( fHarmonic*GetMCDPHI(fDaughterPhi) );
    ((TProfile*)((TList*)fList->FindObject(listName.Data()))->FindObject("COSNDPHIMC" ))->Fill( fDaughterPt, cosn, 1 );
  }
  if(res) {
    ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PHIRes"))->Fill( fDaughterPhi, fDaughterMatchPhi-fDaughterAtSecPhi );
    ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("ETARes"))->Fill( fDaughterEta, fDaughterMatchEta-fDaughterAtSecEta );
    ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PTRes"))->Fill( fDaughterPt, fDaughterMatchPt-fDaughterAtSecPt );
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::LoadTrack(AliESDtrack *myTrack, Double_t aodChi2NDF) {
  fDaughterCharge = myTrack->Charge();
  fDaughterXRows = myTrack->GetTPCCrossedRows();
  fDaughterNFClsTPC = myTrack->GetTPCNclsF();
  fDaughterNSClsTPC = myTrack->GetTPCnclsS();
  fDaughterNClsTPC = myTrack->GetTPCclusters(0);
  if(fReadESD) {
    if(fDaughterNClsTPC>0) fDaughterChi2PerNClsTPC = myTrack->GetTPCchi2()/fDaughterNClsTPC;
  } else {
    fDaughterChi2PerNClsTPC = aodChi2NDF;
  }
  myTrack->GetImpactParameters(fDaughterImpactParameterXY,fDaughterImpactParameterZ);
  fDaughterStatus = myTrack->GetStatus();
  fDaughterITScm = myTrack->GetITSClusterMap();
  fDaughterPhi = myTrack->Phi();
  fDaughterEta = myTrack->Eta();
  fDaughterPt = myTrack->Pt();
  fDaughterKinkIndex = myTrack->GetKinkIndex(0);
  fDaughterNClsITS=0;
  for(Int_t lay=0; lay!=6; ++lay)
    if(TESTBIT(fDaughterITScm,lay)) fDaughterNClsITS++;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptDaughter(Bool_t strongITS) {
  if(fDaughterKinkIndex>0) return kFALSE;
  if( (fDaughterStatus&AliESDtrack::kTPCrefit)==0 ) return kFALSE;
  if(fDaughterNFClsTPC<1) return kFALSE;
  if(fDaughterPt<fDaughterMinPt) return kFALSE;
  if(fDaughterEta<fDaughterMinEta) return kFALSE;
  if(fDaughterEta>fDaughterMaxEta) return kFALSE;
  if(fDaughterNClsTPC<fDaughterMinNClsTPC) return kFALSE;
  if(fDaughterXRows<fDaughterMinXRows) return kFALSE;
  if(fDaughterChi2PerNClsTPC>fDaughterMaxChi2PerNClsTPC) return kFALSE;
  if(TMath::Abs(fDaughterImpactParameterXY)<fDaughterMinImpactParameterXY) return kFALSE;
  if(fDaughterXRows<fDaughterMinXRowsOverNClsFTPC*fDaughterNFClsTPC) return kFALSE;
  if(strongITS) {
    if( (fDaughterITSrefit) & ((fDaughterStatus&AliESDtrack::kITSrefit)==0) ) return kFALSE;
    for(Int_t lay=0; lay!=6; ++lay)
      if(fDaughterITSConfig[lay]>-0.5) {
	if(fDaughterITSConfig[lay]) {
	  if(!TESTBIT(fDaughterITScm,lay)) return kFALSE;
	} else {
	  if(TESTBIT(fDaughterITScm,lay)) return kFALSE;
	}
      }
    if(fDaughterNClsITS<fDaughterMinNClsITS) return kFALSE;
    if(fDaughterSPDRequireAny) {
      if( !TESTBIT(fDaughterITScm,0)&&!TESTBIT(fDaughterITScm,1)) return kFALSE;
    }
  }
  return kTRUE;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::GetWDist(const AliVVertex* v0, const AliVVertex* v1) {
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    printf("One of vertices is not valid\n");
    return 0;
  }
  static TMatrixDSym vVb(3);
  double dist = -1;
  double dx = v0->GetX()-v1->GetX();
  double dy = v0->GetY()-v1->GetY();
  double dz = v0->GetZ()-v1->GetZ();
  double cov0[6],cov1[6];
  v0->GetCovarianceMatrix(cov0);
  v1->GetCovarianceMatrix(cov1);
  vVb(0,0) = cov0[0]+cov1[0];
  vVb(1,1) = cov0[2]+cov1[2];
  vVb(2,2) = cov0[5]+cov1[5];
  vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
  vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
  vVb.InvertFast();
  if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
    +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::plpMV(const AliVEvent *event) {
  // check for multi-vertexer pile-up
  const AliAODEvent *aod = (const AliAODEvent*)event;
  const AliESDEvent *esd = (const AliESDEvent*)event;
  //
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2 = 5.0;
  const double kMinWDist = 15;
  //
  if (!aod && !esd) {
    printf("Event is neither of AOD nor ESD\n");
    exit(1);
  }
  //
  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;
  int nPlp = 0;
  //
  if (aod) {
    if ( !(nPlp=aod->GetNumberOfPileupVerticesTracks()) ) return kFALSE;
    vtPrm = aod->GetPrimaryVertex();
    if (vtPrm == aod->GetPrimaryVertexSPD()) return kTRUE; // there are pile-up vertices but no primary
  }
  else {
    if ( !(nPlp=esd->GetNumberOfPileupVerticesTracks())) return kFALSE;
    vtPrm = esd->GetPrimaryVertexTracks();
    if (((AliESDVertex*)vtPrm)->GetStatus()!=1) return kTRUE; // there are pile-up vertices but no primary
  }
  //int bcPrim = vtPrm->GetBC();
  //
  for (int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = aod ? (const AliVVertex*)aod->GetPileupVertexTracks(ipl) : (const AliVVertex*)esd->GetPileupVertexTracks(ipl);
    //
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF() > kMaxPlpChi2) continue;
    //  int bcPlp = vtPlp->GetBC();
    //  if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other BC
    //
    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist) continue;
    //
    return kTRUE; // pile-up: well separated vertices
  }
  //
  return kFALSE;
  //
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeFilterBits() {
  //FilterBit 1
  fFB1 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  //FilterBit1024
  fFB1024 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  fFB1024->SetMinNCrossedRowsTPC(120);
  fFB1024->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fFB1024->SetMaxChi2PerClusterITS(36);
  fFB1024->SetMaxFractionSharedTPCClusters(0.4);
  fFB1024->SetMaxChi2TPCConstrainedGlobal(36);
  fFB1024->SetEtaRange(-0.9,0.9);
  fFB1024->SetPtRange(0.15, 1e10);
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::PassesFilterBit(AliESDtrack *track) {
  Bool_t ret=kFALSE;
  switch(fRFPFilterBit) {
    case(1024):
      ret = fFB1024->AcceptTrack(track);
      break;
    default:
      ret = fFB1->AcceptTrack(track);
  }
  return ret;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::LoadVZEROResponse() {
  if(fVZEResponse) {
    TString run = fVZEResponse->GetTitle();
    if( run.Atoi() == fRunNumber ) return;
    fVZEResponse = NULL;
  }
  //==>loading
  fVZEResponse = dynamic_cast<TH2D*> (fVZEload->FindObject( Form("%d",fRunNumber) ));
  if(fVZEResponse) {
    printf("New VZE calibration: run %d || %s -> Entries %.0f\n",fRunNumber, fVZEResponse->GetTitle(),fVZEResponse->GetEntries());
    //=>external weights
    for(int i=0;i!=64;++i) fVZEextW[i]=1;
    if(!fVZEsave) {
      Double_t minC = fCentPerMin, maxC = fCentPerMax;
      if(fVZEmb) {
	minC = 0;
	maxC = 80;
      }
      Int_t ybinmin = fVZEResponse->GetYaxis()->FindBin(minC+1e-6);
      Int_t ybinmax = fVZEResponse->GetYaxis()->FindBin(maxC-1e-6);
      if(fSkipCentralitySelection) {
	ybinmin=-1;
	ybinmax=-1;
      }
      for(int i=0;i!=64;++i) fVZEextW[i] = fVZEResponse->Integral(i+1,i+1,ybinmin,ybinmax)/(maxC-minC);
      //ring-wise normalization
      Double_t ring[8];
      for(int j=0; j!=8; ++j) {
	ring[j]=0;
	for(int i=0;i!=8;++i) ring[j] += fVZEextW[j*8+i]/8;
      }
      //disk-wise normalization
      Double_t disk[2];
      int xbinmin, xbinmax;
      xbinmin = 1+8*fVZECa;
      xbinmax = 8+8*fVZECb;
      disk[0] = fVZEResponse->Integral(xbinmin,xbinmax,ybinmin,ybinmax)/(maxC-minC)/(xbinmax-xbinmin+1);
      xbinmin = 33+8*fVZEAa;
      xbinmax = 40+8*fVZEAb;
      disk[1] = fVZEResponse->Integral(xbinmin,xbinmax,ybinmin,ybinmax)/(maxC-minC)/(xbinmax-xbinmin+1);
      //for(int i=0;i!=64;++i) printf("CELL %d -> W = %f ||",i,fVZEextW[i]);
      if(fVZEByDisk) {
	for(int i=0;i!=64;++i) fVZEextW[i] = disk[i/32]/fVZEextW[i];
      } else {
	for(int i=0;i!=64;++i) fVZEextW[i] = ring[i/8]/fVZEextW[i];
      }
      //for(int i=0;i!=64;++i) printf(" W = %f \n",fVZEextW[i]);
    }
  } else {
    printf("VZE calibration: requested but not found!!!\n");
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddVZEQA() {
  TProfile2D *prof;
  TH2D *tH2D;
  TList *tList = new TList();
  tList->SetName( "VZEQA" );
  tList->SetOwner();
  fList->Add( tList );
  tH2D = new TH2D("EQU","EQU;VZEeqmult-VZEmult;cell",100,-5,+5,64,0,64); tList->Add( tH2D );
  prof = new TProfile2D("LINbefCAL","LINbef;VZEcell;VZEeqmult;SPDtrkl", 64,0,64,350,0,700,0,10000); tList->Add( prof );
  prof = new TProfile2D("LINaftCAL","LINaft;VZEcell;VZEeqmult;SPDtrkl", 64,0,64,350,0,700,0,10000); tList->Add( prof );
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillVZEQA() {
  AliESDEvent *tESD = (AliESDEvent*) (InputEvent());
  AliAODEvent *tAOD = (AliAODEvent*) (InputEvent());
  if(fReadESD) {
    if(!tESD) return;
    //FillVZEQA(tESD);
  } else {
    if(!tAOD) return;
    FillVZEQA(tAOD);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillVZEQA(AliAODEvent *tAOD) {
  AliVVZERO *vzero = tAOD->GetVZEROData();
  AliAODTracklets *tracklets = tAOD->GetTracklets();
  if(!vzero) return;
  if(!tracklets) return;
  Double_t mult, eqmult;
  Int_t trkl = tracklets->GetNumberOfTracklets();
  for(int id=0; id!=64; ++id) {
    mult = vzero->GetMultiplicity(id);
    eqmult = tAOD->GetVZEROEqMultiplicity(id);
    ((TH2D*) ((TList*) fList->FindObject("VZEQA"))->FindObject("EQU"))->Fill(eqmult-mult,id);
    ((TProfile2D*) ((TList*) fList->FindObject("VZEQA"))->FindObject( "LINbefCAL" ))->Fill(id,eqmult,trkl,1);
    ((TProfile2D*) ((TList*) fList->FindObject("VZEQA"))->FindObject( "LINaftCAL" ))->Fill(id,eqmult*fVZEextW[id],trkl,1);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddVZEROResponse() {
  fVZEResponse = NULL;
  AliVEvent *event = InputEvent();
  if(!event) return;
  Int_t thisrun = event->GetRunNumber();
  fVZEResponse = new TH2D( Form("%d",thisrun), Form("%d;cell;CC",thisrun), 64,0,64, 110, -10, 100);
  fList->Add(fVZEResponse);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::SaveVZEROResponse() {
  if(!fVZEResponse) return;
  AliVEvent *event = InputEvent();
  if(!event) return;
  Double_t w;
  // reject event with low ocupancy in VZERO
  // for centralities below 60% this should not happen anyway
  Double_t rejectEvent = kFALSE;
  for(int id=0; id!=64; ++id) {
    w = event->GetVZEROEqMultiplicity(id);
    if(w<3) rejectEvent = kTRUE;
  }
  if(rejectEvent) return;
  // saves weights
  for(int id=0; id!=64; ++id) {
    w = event->GetVZEROEqMultiplicity(id);
    fVZEResponse->Fill(id,fThisCent,w);
  }
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::RefMult(AliAODEvent *tAOD, Int_t fb) {
  Int_t found = 0;
  for(int i=0; i!=tAOD->GetNumberOfTracks(); ++i) {
    AliAODTrack *t = dynamic_cast<AliAODTrack*>(tAOD->GetTrack( i ));
    if(!t) continue;
    if( !t->TestFilterBit(fb) ) continue;
    if( t->Eta()<-0.8 || t->Eta()>+0.8 ) continue;
    if( t->Pt()<0.2 || t->Pt()>5.0 ) continue;
    if( t->GetTPCNcls()<70 ) continue;
    //if( t->GetTPCsignal()<10.0 ) continue;
    if( t->Chi2perNDF()<0.2 ) continue;
    ++found;
  }
  return found;
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::RefMultTPC() {
  AliAODEvent *ev = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!ev) return -1;
  Int_t found = 0;
  for(int i=0; i!=ev->GetNumberOfTracks(); ++i) {
    AliAODTrack *t = dynamic_cast<AliAODTrack*>(ev->GetTrack( i ));
    if(!t) continue;
    if( !t->TestFilterBit(1) ) continue;
    if( t->Eta()<-0.8 || t->Eta()>+0.8 ) continue;
    if( t->Pt()<0.2 || t->Pt()>5.0 ) continue;
    if( t->GetTPCNcls()<70 ) continue;
    if( t->GetTPCsignal()<10.0 ) continue;
    if( t->Chi2perNDF()<0.2 ) continue;    
    ++found;
  }
  return found;
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::RefMultGlobal() {
  AliAODEvent *ev = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!ev) return -1;
  Int_t found = 0;
  for(int i=0; i!=ev->GetNumberOfTracks(); ++i) {
    AliAODTrack *t = dynamic_cast<AliAODTrack*>(ev->GetTrack( i ));
    if(!t) continue;
    if( !t->TestFilterBit(16) ) continue;
    if( t->Eta()<-0.8 || t->Eta()>+0.8 ) continue;
    if( t->Pt()<0.2 || t->Pt()>5.0 ) continue;
    if( t->GetTPCNcls()<70 ) continue;
    if( t->GetTPCsignal()<10.0 ) continue;
    if( t->Chi2perNDF()<0.1 ) continue;    
    Double_t b[3], bcov[3];
    if( !t->PropagateToDCA(ev->GetPrimaryVertex(),ev->GetMagneticField(),100,b,bcov) ) continue;
    if( b[0]>+0.3 || b[0]<-0.3 || b[1]>+0.3 || b[1]<-0.3) continue;
    ++found;
  }
  return found;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ResetContainers() {
  if(fUseFP) {
    fTPCevent->ClearFast();
    fVZEevent->ClearFast();
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddTrackVn(TList *tList) {
  TProfile *tProfile;
  TH1D *tH1D;
  // vze
  tProfile = new TProfile("SP_uVZEA","u x Q_{VZEA}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("SP_uVZEC","u x Q_{VZEC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("SP_VZEAVZEC","Q_{VZEA} x Q_{VZEC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("SP_VZEATPC","Q_{VZEA} x Q_{TPC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("SP_VZECTPC","Q_{VZEC} x Q_{TPC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // error
  tProfile = new TProfile("SP_uVZEAuVZEC","u x Q_{VZEA} . u x Q_{VZEC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("SP_uVZEAVZEAVZEC","u x Q_{VZEA} . Q_{VZEA} x Q_{VZEC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("SP_uVZECVZEAVZEC","u x Q_{VZEC} . Q_{VZEA} x Q_{VZEC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // tpc
  tProfile = new TProfile("SP_uTPCA","u x Q_{TPCA}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("SP_uTPCC","u x Q_{TPCC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("SP_TPCATPCC","Q_{TPCA} x Q_{TPCC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // error
  tProfile = new TProfile("SP_uTPCATPCATPCC","u x Q_{TPCA} . Q_{TPCA} x Q_{TPCC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("SP_uTPCCTPCATPCC","u x Q_{TPCC} . Q_{TPCA} x Q_{TPCC}",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // control
  tH1D = new TH1D("QC_HistPt_P","HistPt_P",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("QC_HistPt_Q","HistPt_Q",fPtBins,fPtBinEdge); tList->Add( tH1D );
  // qc
  tProfile = new TProfile("QC_C2","QC_C2",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_C4","QC_C4",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_DC2","QC_DC2",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_DC4","QC_DC4",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_C2C4","QC_C2C4",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_C2DC2","QC_C2DC2",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_C2DC4","QC_C2DC4",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_C4DC2","QC_C4DC2",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_C4DC4","QC_C4DC4",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_DC2DC4","QC_DC2DC4",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // qc transient
  tProfile = new TProfile("QC_pCos","QC_pCos",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_pSin","QC_pSin",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_qCos","QC_qCos",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_qSin","QC_qSin",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_q2hCos","QC_q2hCos",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile("QC_q2hSin","QC_q2hSin",fPtBins,fPtBinEdge,-3,+3,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // measurements
  tH1D = new TH1D("SP_vnVZEA","SP_vnVZEA",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("SP_vnVZEC","SP_vnVZEC",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("SP_vnTPCA","SP_vnTPCA",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("SP_vnTPCC","SP_vnTPCC",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("QC_Cum2","QC_Cum2",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("QC_Cum4","QC_Cum4",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("QC_DCum2","QC_DCum2",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("QC_DCum4","QC_DCum4",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("SP_vnVZEGA","SP_vnVZEGA",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("SP_vnVZEWA","SP_vnVZEWA",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("SP_vnTPCAA","SP_vnTPCAA",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("QC_vn2","QC_vn2",fPtBins,fPtBinEdge); tList->Add( tH1D );
  tH1D = new TH1D("QC_vn4","QC_vn4",fPtBins,fPtBinEdge); tList->Add( tH1D );
  if(fReadMC) {
    TH2D *tH2D;
    tProfile = new TProfile("MC_COSNDPHI","MC_COSNDPHI",fPtBins,fPtBinEdge,-3,+3); tList->Add( tProfile );
    tH2D = new TH2D("MC_COSNDPHI_uQVZEA","MC_COSNDPHI_uQVZEA",100,-1,+1,100,-0.3,+0.3); tList->Add( tH2D );
    tH2D = new TH2D("MC_COSNDPHI_uQVZEC","MC_COSNDPHI_uQVZEC",100,-1,+1,100,-0.3,+0.3); tList->Add( tH2D );
    tH2D = new TH2D("MC_COSNDPHI_uQTPCA","MC_COSNDPHI_uQTPCA",100,-1,+1,100,-0.3,+0.3); tList->Add( tH2D );
    tH2D = new TH2D("MC_COSNDPHI_uQTPCC","MC_COSNDPHI_uQTPCC",100,-1,+1,100,-0.3,+0.3); tList->Add( tH2D );
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddDecayVn(TList *tList) {
  TProfile2D *tProfile;
  TH2D *tH2D;
  // decay
  tH2D = new TH2D("DecayYield_PtMass","Decay_PtMass",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tProfile = new TProfile2D("DecayAvgPt_PtMass","Decay_PtMass",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tProfile ); tProfile->Sumw2();
  // vze
  tProfile = new TProfile2D("SP_uVZEA","u x Q_{VZEA}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("SP_uVZEC","u x Q_{VZEC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("SP_VZEAVZEC","Q_{VZEA} x Q_{VZEC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("SP_VZEATPC","Q_{VZEA} x Q_{TPC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("SP_VZECTPC","Q_{VZEC} x Q_{TPC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // error
  tProfile = new TProfile2D("SP_uVZEAuVZEC","u x Q_{VZEA} . u x Q_{VZEC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("SP_uVZEAVZEAVZEC","u x Q_{VZEA} . Q_{VZEA} x Q_{VZEC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("SP_uVZECVZEAVZEC","u x Q_{VZEC} . Q_{VZEA} x Q_{VZEC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // tpc
  tProfile = new TProfile2D("SP_uTPCA","u x Q_{TPCA}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("SP_uTPCC","u x Q_{TPCC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("SP_TPCATPCC","Q_{TPCA} x Q_{TPCC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // error
  tProfile = new TProfile2D("SP_uTPCATPCATPCC","u x Q_{TPCA} . Q_{TPCA} x Q_{TPCC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("SP_uTPCCTPCATPCC","u x Q_{TPCC} . Q_{TPCA} x Q_{TPCC}",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // control
  tH2D = new TH2D("QC_HistPt_P","HistPt_P",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("QC_HistPt_Q","HistPt_Q",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  // qc
  tProfile = new TProfile2D("QC_C2","QC_C2",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_C4","QC_C4",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_DC2","QC_DC2",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_DC4","QC_DC4",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_C2C4","QC_C2C4",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_C2DC2","QC_C2DC2",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_C2DC4","QC_C2DC4",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_C4DC2","QC_C4DC2",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_C4DC4","QC_C4DC4",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_DC2DC4","QC_DC2DC4",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // qc transient
  tProfile = new TProfile2D("QC_pCos","QC_pCos",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_pSin","QC_pSin",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_qCos","QC_qCos",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_qSin","QC_qSin",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_q2hCos","QC_q2hCos",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  tProfile = new TProfile2D("QC_q2hSin","QC_q2hSin",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass,"s"); tList->Add( tProfile ); tProfile->Sumw2();
  // measurements
  tH2D = new TH2D("SP_vnVZEA","SP_vnVZEA",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("SP_vnVZEC","SP_vnVZEC",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("SP_vnTPCA","SP_vnTPCA",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("SP_vnTPCC","SP_vnTPCC",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("QC_Cum2","QC_Cum2",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("QC_Cum4","QC_Cum4",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("QC_DCum2","QC_DCum2",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("QC_DCum4","QC_DCum4",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("SP_vnVZEGA","SP_vnVZEGA",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("SP_vnVZEWA","SP_vnVZEWA",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("SP_vnTPCAA","SP_vnTPCAA",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("QC_vn2","QC_vn2",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  tH2D = new TH2D("QC_vn4","QC_vn4",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tH2D );
  if(fReadMC) {
    tProfile = new TProfile2D("MC_COSNDPHI","MC_COSNDPHI",fPtBins,fPtBinEdge,fMassBins,fMinMass,fMaxMass); tList->Add( tProfile );
  }
}
//=======================================================================
TList* AliAnalysisTaskFlowStrange::RebinDecayVn(Int_t nbins, Int_t *bins) {
  TList *out = new TList();
  out->SetOwner();
  out->SetName( fList->GetName() );

  TList *ori, *end;

  ori = (TList*) fList->FindObject("Event");
  out->Add( ori );

  ori = (TList*) fList->FindObject("MakeQSpy");
  out->Add( ori );

  ori = (TList*) fList->FindObject("V0SAllVn");
  end = RebinDecayVn(ori,nbins,bins);
  end->SetName( ori->GetName() ); out->Add( end );

  ori = (TList*) fList->FindObject("V0SSelVn");
  end = RebinDecayVn(ori,nbins,bins);
  end->SetName( ori->GetName() ); out->Add( end );

  if(fReadMC) {
    ori = (TList*) fList->FindObject("V0SMthVn");
    end = RebinDecayVn(ori,nbins,bins);
    end->SetName( ori->GetName() ); out->Add( end );

    ori = (TList*) fList->FindObject("V0SMthPosPosVn");
    end = RebinDecayVn(ori,nbins,bins);
    end->SetName( ori->GetName() ); out->Add( end );

    ori = (TList*) fList->FindObject("V0SMthNegNegVn");
    end = RebinDecayVn(ori,nbins,bins);
    end->SetName( ori->GetName() ); out->Add( end );

    ori = (TList*) fList->FindObject("V0SMthPosNegVn");
    end = RebinDecayVn(ori,nbins,bins);
    end->SetName( ori->GetName() ); out->Add( end );

    ori = (TList*) fList->FindObject("V0SUnMthVn");
    end = RebinDecayVn(ori,nbins,bins);
    end->SetName( ori->GetName() ); out->Add( end );
  }

  return out;
}
//=======================================================================
TList* AliAnalysisTaskFlowStrange::RebinDecayVn(TList *tList,Int_t nbins, Int_t *bins) {
  // getting expected number of mass bins
  int sum=0;
  for(int i=0; i!=nbins; ++i) sum += bins[i];

  TList *list = new TList();
  list->SetOwner();

  Int_t npt;
  Double_t pt[200], mass[200];
  TH2D *tH2D;

  tH2D = ((TH2D*)tList->FindObject( "DecayYield_PtMass" ));
  list->Add(tH2D); //keeping it as it is
  int nmassbins = tH2D->GetNbinsY();
  //consistency check
  if( nmassbins!=sum  ) {
    printf("Error: incompatible binning %d vs %d\nBYE\n",nmassbins,sum);
    return NULL;
  }
  //reading pts
  npt = tH2D->GetNbinsX();
  for(int i=0; i!=npt+1; ++i) pt[i] = tH2D->GetXaxis()->GetBinLowEdge(i+1);
  //making mass bins
  for(int i=0,j=0; i!=nbins+1; ++i) {
    mass[i] = tH2D->GetYaxis()->GetBinLowEdge(j+1);
    j += bins[i];
  }
  //TProfile2D migrating info
  TProfile2D *tProfileOld, *tProfileNew;
  TString tprofiles[31] = {"DecayAvgPt_PtMass","SP_uVZEA","SP_uVZEC","SP_VZEAVZEC","SP_VZEATPC",
			   "SP_VZECTPC","SP_uVZEAuVZEC","SP_uVZEAVZEAVZEC","SP_uVZECVZEAVZEC","SP_uTPCA",
			   "SP_uTPCC","SP_TPCATPCC","SP_uTPCATPCATPCC","SP_uTPCCTPCATPCC","QC_C2",
			   "QC_C4","QC_DC2","QC_DC4","QC_C2C4","QC_C2DC2",
			   "QC_C2DC4","QC_C4DC2","QC_C4DC4","QC_DC2DC4","QC_pCos",
			   "QC_pSin","QC_qCos","QC_qSin","QC_q2hCos","QC_q2hSin",
			   "MC_COSNDPHI"};
  for(int i=0; i!=31; ++i) {
    tProfileOld = (TProfile2D*) tList->FindObject( tprofiles[i].Data() );
    if(!tProfileOld) continue;
    TArrayD *oldsw2 = tProfileOld->GetSumw2();
    TArrayD *oldbsw2 = tProfileOld->GetBinSumw2();
    tProfileNew = new TProfile2D( tprofiles[i].Data(), tprofiles[i].Data(),
				  npt,pt,nbins,mass,"s");
    tProfileNew->Sumw2();
    list->Add( tProfileNew );
    TArrayD *newsw2 = tProfileNew->GetSumw2();
    TArrayD *newbsw2 = tProfileNew->GetBinSumw2();
    for(int x=1; x!=tProfileOld->GetNbinsX()+1; ++x) { //pt
      for(int y=1; y!=tProfileNew->GetNbinsY()+1; ++y) { //mass
	Double_t sContent=0;
	Double_t sEntries=0;
	Double_t sSqWeigh=0;
	Double_t sSqBWeig=0;
	Int_t binnew = tProfileNew->GetBin(x,y);
	Double_t minmass = tProfileNew->GetYaxis()->GetBinLowEdge( y );
	Double_t maxmass = tProfileNew->GetYaxis()->GetBinLowEdge( y+1 );
	Int_t minbin = tProfileOld->GetYaxis()->FindBin( minmass+1e-10 );
	Int_t maxbin = tProfileOld->GetYaxis()->FindBin( maxmass-1e-10 );
	for(int k=minbin; k!=maxbin+1; ++k) {
	  Int_t binold = tProfileOld->GetBin(x,k);
	  Double_t wk = tProfileOld->GetBinEntries( binold );
	  Double_t yk = tProfileOld->GetBinContent( binold );
	  sEntries += wk;
	  sContent += wk*yk;
	  sSqWeigh += oldsw2->At(binold);
	  if(oldbsw2->GetSize()) sSqBWeig += oldbsw2->At(binold);
	}
	tProfileNew->SetBinEntries( binnew, sEntries );
	tProfileNew->SetBinContent( binnew, sContent );
	newsw2->SetAt(sSqWeigh, binnew);
	if(oldbsw2->GetSize()) newbsw2->SetAt(sSqBWeig, binnew);
      }
    }
    tProfileOld = NULL;
  }
  //TH2D all new blank
  TString th2ds[15] = {"QC_HisPt_P","QC_HistPt_Q","SP_vnVZEA","SP_vnVZEC","SP_vnTPCA",
		       "SP_vnTPCC","QC_Cum2","QC_Cum4","QC_DCum2","QC_DCum4",
		       "SP_vnVZEGA","SP_vnVZEWA","SP_vnTPCAA","QC_vn2","QC_vn4"};
  for(int i=0; i!=15; ++i) {
    tH2D = new TH2D( th2ds[i].Data(),th2ds[i].Data(),
		     npt,pt,nbins,mass);
    list->Add(tH2D);
  }
  return list;

  //int nmass = tH2Dold->GetNbinsY();

  /*
  // decay
  TString nam[38] = {
  for(int i=0; i!=38; ++i) {
    tProfile = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( nam[i].Data() ));
    if(i==o) { // binning
      nxs = tProfile->GetXaxis()->GetNbins();
      nys = tProfile->GetYaxis()->GetNbins();
      if(pt) {
	nxs = nbins;
	*nx = *bins;
	for(int y=0; y!=nys; ++y)
	  ny[y] = tProfile->GetYaxis()->GetBinLowEdge(y+1);
      } else {
	nys = nbins;
	*ny = *bins;
	for(int x=0; x!=nxs; ++x)
	  nx[x] = tProfile->GetXaxis()->GetBinLowEdge(x+1);
      }
    }
    tProfileNew = new TProfile2D(tProfile->GetName(),tProfile->GetTitle(),nxs,nx,nys,ny,"s"); list->Add( tProfileNew );
    if(pt) {
      for(int y=0; y!=nys; ++y) {
	
	for(int x=0; x!=nxs; ++x) {
	}
      }
    } else {

    }
  }
  */
}
//=======================================================================
void AliAnalysisTaskFlowStrange::QCStoreTrackVn(TString name) {
  // getting transients
  TProfile *pCos = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pCos" ));
  TProfile *pSin = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pSin" ));
  TProfile *qCos = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qCos" ));
  TProfile *qSin = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qSin" ));
  TProfile *q2hCos = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hCos" ));
  TProfile *q2hSin = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hSin" ));
  if(fDebug) {
    printf("**QCStoreTrackVn( %s )\n",name.Data());
    printf(" Read %.0f entries in p set and %.0f entries in q set\n",pCos->GetEntries(),qCos->GetEntries());
    printf(" pCos[5] %.16f | pSin[5] %.16f \n", pCos->GetBinContent(5), pSin->GetBinContent(5));
    printf(" qCos[5] %.16f | qSin[5] %.16f \n", qCos->GetBinContent(5), qSin->GetBinContent(5));
    printf(" q2hCos[5] %.16f | q2hSin[5] %.16f \n", q2hCos->GetBinContent(5), q2hSin->GetBinContent(5));
  }
  // filling {2} and {4} correlator
  TProfile *c2 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2" ));
  TProfile *c4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4" ));
  TProfile *dc2 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC2" ));
  TProfile *dc4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC4" ));
  TProfile *c2c4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2C4" ));
  TProfile *c2dc2 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2DC2" ));
  TProfile *c2dc4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2DC4" ));
  TProfile *c4dc2 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4DC2" ));
  TProfile *c4dc4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4DC4" ));
  TProfile *dc2dc4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC2DC4" ));
  double tpc_qcos = fQTPCACos + fQTPCCCos;
  double tpc_qsin = fQTPCASin + fQTPCCSin;
  double tpc_qsqr = tpc_qcos*tpc_qcos + tpc_qsin*tpc_qsin;
  double tpc_q2hsqr = fQTPC2hCos*fQTPC2hCos + fQTPC2hSin*fQTPC2hSin;
  double tpc_qmul = fQTPCA + fQTPCC;
  int n = c2->GetNbinsX();
  for(int i=1; i!=n+1; ++i) {
    double mp = pCos->GetBinEntries( i );
    if( mp<1 ) { if(fDebug>2) printf(" bin %d:: mp (%.16f) < 1!\n",i,mp); continue; }
    double mm1 = tpc_qmul*(tpc_qmul-1);
    if( mm1<1e-100 ) { if(fDebug>2) printf(" bin %d:: mm1<1e-100!\n",i); continue; }
    double mq = qCos->GetBinEntries( i );
    double mpmmq = mp*tpc_qmul-mq;
    if( mpmmq<1e-100 )  { if(fDebug>2) printf(" bin %d:: mpmmq<1e-100!\n",i); continue; }
    double pt = c2->GetBinCenter( i );
    double pcos = pCos->GetBinContent( i )*mp;
    double psin = pSin->GetBinContent( i )*mp;
    double qcos = qCos->GetBinContent( i )*mq;
    double qsin = qSin->GetBinContent( i )*mq;
    double q2hcos = q2hCos->GetBinContent( i )*mq;
    double q2hsin = q2hSin->GetBinContent( i )*mq;
    double pQ = pcos*tpc_qcos+psin*tpc_qsin;
    double q2nQnQn = (qcos*tpc_qcos + qsin*tpc_qsin)*tpc_qcos + (qsin*tpc_qcos-qcos*tpc_qsin)*tpc_qsin;
    double pnQnQ2n = (pcos*tpc_qcos - psin*tpc_qsin)*fQTPC2hCos + (psin*tpc_qcos+pcos*tpc_qsin)*fQTPC2hSin;
    double tC2 = (tpc_qsqr-tpc_qmul)/mm1;
    double tDC2 = (pQ-mq)/mpmmq;
    c2->Fill( pt, tC2, mm1 );
    dc2->Fill( pt, tDC2, mpmmq );
    c2dc2->Fill( pt, tC2*tDC2, mm1*mpmmq );
    double mm1m2m3 = tpc_qmul*(tpc_qmul-1)*(tpc_qmul-2)*(tpc_qmul-3);
    if(mm1m2m3<1e-100) continue;
    double mpm3mqm1m2 = (mp*tpc_qmul-3*mq)*(tpc_qmul-1)*(tpc_qmul-2);
    if(mpm3mqm1m2<1e-100) continue;
    double req2hqnqn = fQTPC2hCos*(tpc_qcos*tpc_qcos-tpc_qsin*tpc_qsin)+2*fQTPC2hSin*tpc_qcos*tpc_qsin;
    double tC4 = (tpc_qsqr*tpc_qsqr + tpc_q2hsqr - 2*req2hqnqn - 2*(2*tpc_qsqr*(tpc_qmul-2)-tpc_qmul*(tpc_qmul-3)))/mm1m2m3;
    double tDC4 = pQ*tpc_qsqr -q2nQnQn -pnQnQ2n -2*(tpc_qmul-1)*pQ -2*mq*(tpc_qsqr-tpc_qmul+3) +6*(qcos*tpc_qcos+qsin*tpc_qsin) +(q2hcos*fQTPC2hCos+q2hsin*fQTPC2hSin);
    tDC4 /= mpm3mqm1m2;
    c4->Fill( pt, tC4, mm1m2m3 );
    dc4->Fill( pt, tDC4, mpm3mqm1m2 );
    c2c4->Fill( pt, tC2*tC4, mm1*mm1m2m3 );
    c2dc4->Fill( pt, tC2*tDC4, mm1*mpm3mqm1m2 );
    c4dc2->Fill( pt, tC4*tDC2, mm1m2m3*mpmmq );
    c4dc4->Fill( pt, tC4*tDC4, mm1m2m3*mpm3mqm1m2 );
    dc2dc4->Fill( pt, tDC2*tDC4, mpmmq*mpm3mqm1m2 );
  }
  // clean for next event
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qCos" ))->Reset();
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qSin" ))->Reset();
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hCos" ))->Reset();
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hSin" ))->Reset();
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pCos" ))->Reset();
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pSin" ))->Reset();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::QCStoreDecayVn(TString name) {
  // getting transients
  TProfile2D *pCos = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pCos" ));
  TProfile2D *pSin = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pSin" ));
  TProfile2D *qCos = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qCos" ));
  TProfile2D *qSin = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qSin" ));
  TProfile2D *q2hCos = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hCos" ));
  TProfile2D *q2hSin = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hSin" ));
  if(fDebug) {
    printf("**QCStoreTrackVn( %s )\n",name.Data());
    printf(" Read %.0f entries in p set and %.0f entries in q set\n",pCos->GetEntries(),qCos->GetEntries());
    printf(" pCos[5][5] %.16f | pSin[5][5] %.16f \n", pCos->GetBinContent(5,5), pSin->GetBinContent(5,5));
    printf(" qCos[5][5] %.16f | qSin[5][5] %.16f \n", qCos->GetBinContent(5,5), qSin->GetBinContent(5,5));
    printf(" q2hCos[5][5] %.16f | q2hSin[5][5] %.16f \n", q2hCos->GetBinContent(5,5), q2hSin->GetBinContent(5,5));
  }
  // filling {2} and {4} correlator
  TProfile2D *c2 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2" ));
  TProfile2D *c4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4" ));
  TProfile2D *dc2 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC2" ));
  TProfile2D *dc4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC4" ));
  TProfile2D *c2c4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2C4" ));
  TProfile2D *c2dc2 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2DC2" ));
  TProfile2D *c2dc4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2DC4" ));
  TProfile2D *c4dc2 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4DC2" ));
  TProfile2D *c4dc4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4DC4" ));
  TProfile2D *dc2dc4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC2DC4" ));
  double tpc_qcos = fQTPCACos + fQTPCCCos;
  double tpc_qsin = fQTPCASin + fQTPCCSin;
  double tpc_qsqr = tpc_qcos*tpc_qcos + tpc_qsin*tpc_qsin;
  double tpc_q2hsqr = fQTPC2hCos*fQTPC2hCos + fQTPC2hSin*fQTPC2hSin;
  double tpc_qmul = fQTPCA + fQTPCC;
  int n = c2->GetNbinsX();
  int m = c2->GetNbinsY();
  for(int i=1; i!=n+1; ++i) {
    double pt = c2->GetXaxis()->GetBinCenter( i );
    for(int j=1; j!=m+1; ++j) {
      double ms = c2->GetYaxis()->GetBinCenter( j );
      int k = pCos->GetBin(i,j);
      double mp = pCos->GetBinEntries( k );
      if( mp<1 ) { if(fDebug>2) printf(" bin %d,%d:: mp (%.16f) < 1!\n",i,j,mp); continue; }
      double mm1 = tpc_qmul*(tpc_qmul-1);
      if( mm1<1e-100 ) { if(fDebug>2) printf(" bin %d,%d:: mm1<1e-100!\n",i,j); continue; }
      double mq = qCos->GetBinEntries( k );
      double mpmmq = mp*tpc_qmul-mq;
      if( mpmmq<1e-100 )  { if(fDebug>2) printf(" bin %d,%d:: mpmmq<1e-100!\n",i,j); continue; }
      double pcos = pCos->GetBinContent( i,j )*mp;
      double psin = pSin->GetBinContent( i,j )*mp;
      double qcos = qCos->GetBinContent( i,j )*mq;
      double qsin = qSin->GetBinContent( i,j )*mq;
      double q2hcos = q2hCos->GetBinContent( i,j )*mq;
      double q2hsin = q2hSin->GetBinContent( i,j )*mq;
      double pQ = pcos*tpc_qcos+psin*tpc_qsin;
      double q2nQnQn = (qcos*tpc_qcos + qsin*tpc_qsin)*tpc_qcos + (qsin*tpc_qcos-qcos*tpc_qsin)*tpc_qsin;
      double pnQnQ2n = (pcos*tpc_qcos - psin*tpc_qsin)*fQTPC2hCos + (psin*tpc_qcos+pcos*tpc_qsin)*fQTPC2hSin;
      double tC2 = (tpc_qsqr-tpc_qmul)/mm1;
      double tDC2 = (pQ-mq)/mpmmq;
      c2->Fill( pt, ms, tC2, mm1 );
      dc2->Fill( pt, ms, tDC2, mpmmq );
      c2dc2->Fill( pt, ms, tC2*tDC2, mm1*mpmmq );
      double mm1m2m3 = tpc_qmul*(tpc_qmul-1)*(tpc_qmul-2)*(tpc_qmul-3);
      if(mm1m2m3<1e-100) continue;
      double mpm3mqm1m2 = (mp*tpc_qmul-3*mq)*(tpc_qmul-1)*(tpc_qmul-2);
      if(mpm3mqm1m2<1e-100) continue;
      double req2hqnqn = fQTPC2hCos*(tpc_qcos*tpc_qcos-tpc_qsin*tpc_qsin)+2*fQTPC2hSin*tpc_qcos*tpc_qsin;
      double tC4 = (tpc_qsqr*tpc_qsqr + tpc_q2hsqr - 2*req2hqnqn - 2*(2*tpc_qsqr*(tpc_qmul-2)-tpc_qmul*(tpc_qmul-3)))/mm1m2m3;
      double tDC4 = pQ*tpc_qsqr -q2nQnQn -pnQnQ2n -2*(tpc_qmul-1)*pQ -2*mq*(tpc_qsqr-tpc_qmul+3) +6*(qcos*tpc_qcos+qsin*tpc_qsin) +(q2hcos*fQTPC2hCos+q2hsin*fQTPC2hSin);
      tDC4 /= mpm3mqm1m2;
      c4->Fill( pt, ms, tC4, mm1m2m3 );
      dc4->Fill( pt, ms, tDC4, mpm3mqm1m2 );
      c2c4->Fill( pt, ms, tC2*tC4, mm1*mm1m2m3 );
      c2dc4->Fill( pt, ms, tC2*tDC4, mm1*mpm3mqm1m2 );
      c4dc2->Fill( pt, ms, tC4*tDC2, mm1m2m3*mpmmq );
      c4dc4->Fill( pt, ms, tC4*tDC4, mm1m2m3*mpm3mqm1m2 );
      dc2dc4->Fill( pt, ms, tDC2*tDC4, mpmmq*mpm3mqm1m2 );
    }
  }
  // clean for next event
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qCos" ))->Reset();
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qSin" ))->Reset();
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hCos" ))->Reset();
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hSin" ))->Reset();
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pCos" ))->Reset();
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pSin" ))->Reset();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillTrackVn(TString name,Double_t pt,Double_t phi,Double_t eta,Int_t fid) {
  // reading vze qm
  Double_t vzec_qmcos = fQVZECCos/fQVZEC;
  Double_t vzec_qmsin = fQVZECSin/fQVZEC;
  Double_t vzea_qmcos = fQVZEACos/fQVZEA;
  Double_t vzea_qmsin = fQVZEASin/fQVZEA;
  // reading tpc qm
  Double_t tpcc_qmcos = fQTPCCCos/fQTPCC;
  Double_t tpcc_qmsin = fQTPCCSin/fQTPCC;
  Double_t tpca_qmcos = fQTPCACos/fQTPCA;
  Double_t tpca_qmsin = fQTPCASin/fQTPCA;
  Double_t qtpc = fQTPCA+fQTPCC;
  Double_t tpc_qmcos = (fQTPCACos+fQTPCCCos)/qtpc;
  Double_t tpc_qmsin = (fQTPCASin+fQTPCCSin)/qtpc;
  // computing u
  Double_t cosn = TMath::Cos(fHarmonic*phi);
  Double_t sinn = TMath::Sin(fHarmonic*phi);
  Double_t cos2n = TMath::Cos(2.0*fHarmonic*phi);
  Double_t sin2n = TMath::Sin(2.0*fHarmonic*phi);
  // Scalar Product
  Double_t uQ, uQa, uQc, qaqc;
  // filling flow with vze
  qaqc = (vzea_qmcos*vzec_qmcos + vzea_qmsin*vzec_qmsin);
  uQa = (cosn*vzea_qmcos + sinn*vzea_qmsin);
  uQc = (cosn*vzec_qmcos + sinn*vzec_qmsin);
  Double_t cosmc = TMath::Cos( fHarmonic*GetMCDPHI(phi) );
  if(fReadMC) {
    ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "MC_COSNDPHI_uQVZEA" ))->Fill( cosmc,uQa );
    ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "MC_COSNDPHI_uQVZEC" ))->Fill( cosmc,uQc );
    ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "MC_COSNDPHI" ))->Fill( pt,cosmc );
  }
  Double_t qaqt = (vzea_qmcos*tpc_qmcos + vzea_qmsin*tpc_qmsin);
  Double_t qcqt = (vzec_qmcos*tpc_qmcos + vzec_qmsin*tpc_qmsin);
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEA" ))->Fill( pt,uQa,fQVZEA );
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEC" ))->Fill( pt,uQc,fQVZEC );
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZEAVZEC" ))->Fill( pt,qaqc,fQVZEA*fQVZEC );
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZEATPC" ))->Fill( pt,qaqt,fQVZEA*qtpc );
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZECTPC" ))->Fill( pt,qcqt,fQVZEC*qtpc );
  // error vze
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEAuVZEC" ))->Fill( pt,uQa*uQc,fQVZEA*fQVZEC );
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEAVZEAVZEC" ))->Fill( pt,uQa*qaqc,fQVZEA*fQVZEA*fQVZEC );
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZECVZEAVZEC" ))->Fill( pt,uQc*qaqc,fQVZEC*fQVZEA*fQVZEC );
  // filling flow with tpc
  qaqc = (tpca_qmcos*tpcc_qmcos + tpca_qmsin*tpcc_qmsin);
  if(eta<0) {
    uQ = (cosn*tpca_qmcos + sinn*tpca_qmsin);
    ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCA" ))->Fill( pt,uQ,fQTPCA );
    ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCATPCATPCC" ))->Fill( pt,uQ*qaqc,fQTPCA*fQTPCA*fQTPCC );
    if(fReadMC) {
      ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "MC_COSNDPHI_uQTPCA" ))->Fill( cosmc,uQ );
    }
  } else {
    uQ = (cosn*tpcc_qmcos + sinn*tpcc_qmsin);
    ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCC" ))->Fill( pt,uQ,fQTPCC );
    ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCCTPCATPCC" ))->Fill( pt,uQ*qaqc,fQTPCC*fQTPCA*fQTPCC );
    if(fReadMC) {
      ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "MC_COSNDPHI_uQTPCC" ))->Fill( cosmc,uQ );
    }
  }
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_TPCATPCC" ))->Fill( pt,qaqc,fQTPCA*fQTPCC );
  // QC
  ((TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_HistPt_P" ))->Fill( pt );
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pCos" ))->Fill( pt, cosn, 1.0 );
  ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pSin" ))->Fill( pt, sinn, 1.0 );
  if(InQTPC(fid)) {
    ((TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_HistPt_Q" ))->Fill( pt );
    ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qCos" ))->Fill( pt, cosn, 1.0 );
    ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qSin" ))->Fill( pt, sinn, 1.0 );
    ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hCos" ))->Fill( pt, cos2n, 1.0 );
    ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hSin" ))->Fill( pt, sin2n, 1.0 );
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillDecayVn(TString name,Double_t ms,Double_t pt,Double_t phi,Double_t eta,Int_t fid1,Int_t fid2) {
  ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "DecayYield_PtMass" ))->Fill( pt,ms );
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "DecayAvgPt_PtMass" ))->Fill( pt,ms,pt );

  // reading vze qm
  Double_t vzec_qmcos = fQVZECCos/fQVZEC;
  Double_t vzec_qmsin = fQVZECSin/fQVZEC;
  Double_t vzea_qmcos = fQVZEACos/fQVZEA;
  Double_t vzea_qmsin = fQVZEASin/fQVZEA;
  // reading tpc qm
  Double_t tpcc_qmcos = fQTPCCCos/fQTPCC;
  Double_t tpcc_qmsin = fQTPCCSin/fQTPCC;
  Double_t tpca_qmcos = fQTPCACos/fQTPCA;
  Double_t tpca_qmsin = fQTPCASin/fQTPCA;
  Double_t qtpc = fQTPCA+fQTPCC;
  Double_t tpc_qmcos = (fQTPCACos+fQTPCCCos)/qtpc;
  Double_t tpc_qmsin = (fQTPCASin+fQTPCCSin)/qtpc;
  // computing u
  Double_t cosn = TMath::Cos(fHarmonic*phi);
  Double_t sinn = TMath::Sin(fHarmonic*phi);
  Double_t cos2n = TMath::Cos(2.0*fHarmonic*phi);
  Double_t sin2n = TMath::Sin(2.0*fHarmonic*phi);
  // Scalar Product
  Double_t uQ, uQa, uQc, qaqc;
  // filling flow with vze
  qaqc = (vzea_qmcos*vzec_qmcos + vzea_qmsin*vzec_qmsin);
  uQa = (cosn*vzea_qmcos + sinn*vzea_qmsin);
  uQc = (cosn*vzec_qmcos + sinn*vzec_qmsin);
  Double_t cosmc = TMath::Cos( fHarmonic*GetMCDPHI(phi) );
  if(fReadMC) {
    ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "MC_COSNDPHI" ))->Fill( pt,ms,cosmc );
  }
  Double_t qaqt = (vzea_qmcos*tpc_qmcos + vzea_qmsin*tpc_qmsin);
  Double_t qcqt = (vzec_qmcos*tpc_qmcos + vzec_qmsin*tpc_qmsin);
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEA" ))->Fill( pt,ms,uQa,fQVZEA );
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEC" ))->Fill( pt,ms,uQc,fQVZEC );
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZEAVZEC" ))->Fill( pt,ms,qaqc,fQVZEA*fQVZEC );
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZEATPC" ))->Fill( pt,ms,qaqt,fQVZEA*qtpc );
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZECTPC" ))->Fill( pt,ms,qcqt,fQVZEC*qtpc );
  // error vze
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEAuVZEC" ))->Fill( pt,ms,uQa*uQc,fQVZEA*fQVZEC );
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEAVZEAVZEC" ))->Fill( pt,ms,uQa*qaqc,fQVZEA*fQVZEA*fQVZEC );
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZECVZEAVZEC" ))->Fill( pt,ms,uQc*qaqc,fQVZEC*fQVZEA*fQVZEC );
  // filling flow with tpc
  qaqc = (tpca_qmcos*tpcc_qmcos + tpca_qmsin*tpcc_qmsin);
  if(eta<0) {
    uQ = (cosn*tpca_qmcos + sinn*tpca_qmsin);
    ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCA" ))->Fill( pt,ms,uQ,fQTPCA );
    ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCATPCATPCC" ))->Fill( pt,ms,uQ*qaqc,fQTPCA*fQTPCA*fQTPCC );
  } else {
    uQ = (cosn*tpcc_qmcos + sinn*tpcc_qmsin);
    ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCC" ))->Fill( pt,ms,uQ,fQTPCC );
    ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCCTPCATPCC" ))->Fill( pt,ms,uQ*qaqc,fQTPCC*fQTPCA*fQTPCC );
  }
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_TPCATPCC" ))->Fill( pt,ms,qaqc,fQTPCA*fQTPCC );
  // QC
  ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_HistPt_P" ))->Fill( pt,ms );
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pCos" ))->Fill( pt, ms, cosn, 1.0 );
  ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_pSin" ))->Fill( pt, ms, sinn, 1.0 );
  if(InQTPC(fid1)||InQTPC(fid2)) {
    ((TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_HistPt_Q" ))->Fill( pt,ms );
    ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qCos" ))->Fill( pt, ms, cosn, 1.0 );
    ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_qSin" ))->Fill( pt, ms, sinn, 1.0 );
    ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hCos" ))->Fill( pt, ms, cos2n, 1.0 );
    ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_q2hSin" ))->Fill( pt, ms, sin2n, 1.0 );
  }
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::InQTPC(Int_t id) {
  Bool_t ret = kFALSE;
  for(int i=0; i!=fQTPCA_nTracks; ++i)
    if(fQTPCA_fID[i]==id) {
      ret=kTRUE;
      break;
    }
  if(ret) return kTRUE;
  for(int i=0; i!=fQTPCC_nTracks; ++i)
    if(fQTPCC_fID[i]==id) {
      ret=kTRUE;
      break;
    }
  return ret;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ComputeTrackVn(TString name) {
  TProfile *uQa, *uQc, *qaqc, *uQaqaqc, *uQcqaqc;
  TArrayD *pasww, *pbsww, *pcsww;
  //ScalarProducr TPC
  printf("<<%s>> SP TPC\n",name.Data());
  uQa  = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCA" );
  uQc  = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCC" );
  qaqc = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_TPCATPCC" );
  uQaqaqc = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCATPCATPCC" );
  uQcqaqc = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCCTPCATPCC" );
  pasww = uQa->GetBinSumw2();
  pbsww = uQc->GetBinSumw2();
  pcsww = qaqc->GetBinSumw2();
  //
  TH1D *sptpca = (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnTPCA" );
  TH1D *sptpcc = (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnTPCC" );
  TH1D *sptpcaa = (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnTPCAA" );
  //
  for(Int_t i=1; i!=sptpcaa->GetNbinsX()+1; ++i) {
    sptpcaa->SetBinContent(i,0);
    sptpcaa->SetBinError(i,0);
    sptpca->SetBinContent(i,0);
    sptpca->SetBinError(i,0);
    sptpcc->SetBinContent(i,0);
    sptpcc->SetBinError(i,0);
    double a = uQa->GetBinContent(i);
    double b = uQc->GetBinContent(i);
    double c = qaqc->GetBinContent(i);
    //if(TMath::AreEqualAbs(a,0,1e-100)) continue;
    //if(TMath::AreEqualAbs(b,0,1e-100)) continue;
    if(c<1e-100) continue;
    // nominal sptpca
    double vna = a/TMath::Sqrt(c);
    sptpca->SetBinContent(i,vna);
    // nominal sptpcc
    double vnc = b/TMath::Sqrt(c);
    sptpcc->SetBinContent(i,vnc);
    // nominal sptpc
    double vn = (vna + vnc)/2.0;
    sptpcaa->SetBinContent(i,vn);
    // errors
    double asw = uQa->GetBinEntries(i);
    double bsw = uQc->GetBinEntries(i);
    double csw = qaqc->GetBinEntries(i);
    if(asw<1e-100||bsw<1e-100||csw<1e-100) continue;
    double asww = pasww->At(i);
    double bsww = pbsww->At(i);
    double csww = pcsww->At(i);
    if(asww<1e-100||bsww<1e-100||csww<1e-100) continue;
    if((1<1e-100+asww/asw/asw)||(1<1e-100+bsww/bsw/bsw)||(1<1e-100+csww/csw/csw)) continue;
    if(TMath::AreEqualAbs(asww,asw*asw,1e-100)||
       TMath::AreEqualAbs(bsww,bsw*bsw,1e-100)||
       TMath::AreEqualAbs(csww,csw*csw,1e-100)) continue;
    double ac = uQaqaqc->GetBinContent(i);
    double bc = uQcqaqc->GetBinContent(i);
    double acsw = uQaqaqc->GetBinEntries(i);
    double bcsw = uQcqaqc->GetBinEntries(i);
    double ea = uQa->GetBinError(i)*TMath::Sqrt(asww)/asw/TMath::Sqrt(1-asww/asw/asw);
    double eb = uQc->GetBinError(i)*TMath::Sqrt(bsww)/bsw/TMath::Sqrt(1-bsww/bsw/bsw);
    double ec = qaqc->GetBinError(i)*TMath::Sqrt(csww)/csw/TMath::Sqrt(1-csww/csw/csw);
    //printf("%d >> ea^2 %.16f |||| asww %.16f | asw %.16f | 1-asww/asw/asw %.16f \n", i,ea, asww, asw, 1-asww/asw/asw);
    //printf("%d >> eb^2 %.16f |||| bsww %.16f | bsw %.16f | 1-bsww/bsw/bsw %.16f \n", i,eb, bsww, bsw, 1-bsww/bsw/bsw);
    //printf("%d >> ec^2 %.16f |||| csww %.16f | csw %.16f | 1-csww/csw/csw %.16f \n", i,ec, csww, csw, 1-csww/csw/csw);
    double ebc = (bc-b*c)/(1-bcsw/bsw/csw)*bcsw/bsw/csw;
    double eac = (ac-a*c)/(1-acsw/asw/csw)*acsw/asw/csw;
    double evna = 1.0/TMath::Abs(c) * ( ea*ea + vna*vna/TMath::Abs(c)/4.0*ec*ec - vna/TMath::Sqrt(c)*eac );
    double evnc = 1.0/TMath::Abs(c) * ( eb*eb + vnc*vnc/TMath::Abs(c)/4.0*ec*ec - vnc/TMath::Sqrt(c)*ebc );
    //printf("%d >> evna^2 %.16f |||| ea %.16f | ec %.16f | eac %.16f | c %.16f\n", i,evna, ea, ec, eac,c);
    //printf("%d >> evnc^2 %.16f |||| eb %.16f | ec %.16f | ebc %.16f | c %.16f\n", i,evnc, eb, ec, ebc,c);
    if(evna>1e-100) evna = TMath::Sqrt(evna); else evna=0;
    if(evnc>1e-100) evnc = TMath::Sqrt(evnc); else evnc=0;
    sptpca->SetBinError(i,evna);
    sptpcc->SetBinError(i,evnc);
    sptpcaa->SetBinError(i,TMath::Sqrt(evna*evna+evnc*evnc)/2.0);
  }
  //ScalarProduct VZE
  printf("<<%s>> SP VZE\n",name.Data());
  double cvzea2 = ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("ChiSquaredVZEA"))->GetBinContent( 1 );
  double cvzec2 = ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("ChiSquaredVZEC"))->GetBinContent( 1 );
  if( TMath::AreEqualAbs(cvzea2+cvzec2,0,1e-100) ) return;
  uQa  = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEA" );
  uQc  = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEC" );
  qaqc = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZEAVZEC" );
  uQaqaqc = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEAVZEAVZEC" );
  uQcqaqc = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZECVZEAVZEC" );
  pasww = uQa->GetBinSumw2();
  pbsww = uQc->GetBinSumw2();
  pcsww = qaqc->GetBinSumw2();
  //
  TProfile *qaqt = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZEATPC" );
  TProfile *qcqt = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZECTPC" );
  TProfile *uQauQc  = (TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEAuVZEC" );
  //
  TH1D *spvzea = (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnVZEA" );
  TH1D *spvzec = (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnVZEC" );
  TH1D *spvzega = (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnVZEGA" );
  TH1D *spvzewa = (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnVZEWA" );
  for(Int_t i=1; i!=spvzewa->GetNbinsX()+1; ++i) {
    spvzega->SetBinContent(i,0);
    spvzega->SetBinError(i,0);
    spvzewa->SetBinContent(i,0);
    spvzewa->SetBinError(i,0);
    spvzea->SetBinContent(i,0);
    spvzea->SetBinError(i,0);
    spvzec->SetBinContent(i,0);
    spvzec->SetBinError(i,0);
    double asw = uQa->GetBinEntries(i);
    double bsw = uQc->GetBinEntries(i);
    double csw = qaqc->GetBinEntries(i);
    if(asw<1e-1||bsw<1e-1||csw<1e-1) continue;
    double asww = pasww->At(i);
    double bsww = pbsww->At(i);
    double csww = pcsww->At(i);
    if(asww<1e-1||bsww<1e-1||csww<1e-1) continue;
    if((1<asww/asw/asw)||(1<bsww/bsw/bsw)||(1<csww/csw/csw)) continue;
    double a = uQa->GetBinContent(i);
    double b = uQc->GetBinContent(i);
    double c = qaqc->GetBinContent(i);
    double at = qaqt->GetBinContent(i);
    double bt = qcqt->GetBinContent(i);
    if(TMath::AreEqualAbs(a,0,1e-10)) continue;
    if(TMath::AreEqualAbs(b,0,1e-10)) continue;
    if(TMath::AreEqualAbs(c,0,1e-10)) continue;
    if(TMath::AreEqualAbs(at,0,1e-10)) continue;
    if(TMath::AreEqualAbs(bt,0,1e-10)) continue;
    // nominal spvzea
    double aa = c*at/bt;
    if(aa<1e-100) continue;
    double vna = a/TMath::Sqrt(aa);
    spvzea->SetBinContent(i,vna);
    // nominal spvzec
    double bb = c*bt/at;
    if(bb<1e-100) continue;
    double vnc = b/TMath::Sqrt(bb);
    spvzec->SetBinContent(i,vnc);
    //nominal spvzewa
    double vnwa = (cvzea2*vna + cvzec2*vnc) / (cvzea2+cvzec2);
    spvzewa->SetBinContent(i,vnwa);
    // nominal spvzega
    double vnga = a*b/c;
    if(vnga<1e-100) continue;
    vnga = TMath::Sqrt(vnga);
    spvzega->SetBinContent(i,vnga);
    // errors
    double ab = uQauQc->GetBinContent(i);
    double ac = uQaqaqc->GetBinContent(i);
    double bc = uQcqaqc->GetBinContent(i);
    double absw = uQauQc->GetBinEntries(i);
    double acsw = uQaqaqc->GetBinEntries(i);
    double bcsw = uQcqaqc->GetBinEntries(i);
    double ea = uQa->GetBinError(i)*TMath::Sqrt(asww)/asw/TMath::Sqrt(1-asww/asw/asw);
    double eb = uQc->GetBinError(i)*TMath::Sqrt(bsww)/bsw/TMath::Sqrt(1-bsww/bsw/bsw);
    double ec = qaqc->GetBinError(i)*TMath::Sqrt(csww)/csw/TMath::Sqrt(1-csww/csw/csw);
    if(TMath::AreEqualAbs(1,absw/asw/bsw,1e-100)||TMath::AreEqualAbs(1,bcsw/bsw/csw,1e-100)||TMath::AreEqualAbs(1,acsw/asw/csw,1e-100)) continue;
    double eab = (ab-a*b)/(1-absw/asw/bsw)*absw/asw/bsw;
    double ebc = (bc-b*c)/(1-bcsw/bsw/csw)*bcsw/bsw/csw;
    double eac = (ac-a*c)/(1-acsw/asw/csw)*acsw/asw/csw;
    double nc, nec, neac, nebc;
    nc = c*at/bt;
    nec = ec*at/bt;
    neac = eac*at/bt;
    nebc = ebc*at/bt;
    double evna = 1.0/TMath::Abs(nc*nc*nc) * ( nc*nc*ea*ea + a*a/4.0*nec*nec - a*TMath::Abs(nc)*neac*neac );
    nc = c*bt/at;
    nec = ec*bt/at;
    neac = eac*bt/at;
    nebc = ebc*bt/at;
    double evnc = 1.0/TMath::Abs(nc*nc*nc) * ( nc*nc*eb*eb + b*b/4.0*nec*nec - b*TMath::Abs(nc)*nebc*nebc );
    if(evna>1e-100) evna = TMath::Sqrt(evna); else evna=0;
    if(evnc>1e-100) evnc = TMath::Sqrt(evnc); else evnc=0;
    spvzea->SetBinError(i,evna);
    spvzec->SetBinError(i,evnc);
    double evnwa = TMath::Sqrt( cvzea2*cvzea2*evna*evna + cvzec2*cvzec2*evnc*evnc )/(cvzea2+cvzec2);
    spvzewa->SetBinError(i,evnwa);
    double evnga = 0.25/c/c * ( TMath::Abs(b*c/a)*ea*ea + TMath::Abs(a*c/b)*eb*eb + TMath::Abs(a*b/c)*ec*ec + 2*c*eab - 2*a*ebc - 2*b*eac );
    if(evnga>1e-100) evnga = TMath::Sqrt(evnga); else evnga=0;
    spvzega->SetBinError(i,evnga);
  }
  printf("<<%s>> QC TPC\n",name.Data());
  //Qcumulants
  TH1D *resC2 = (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_Cum2" );
  TH1D *resC4 = (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_Cum4" );
  TH1D *resDC2= (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DCum2" );
  TH1D *resDC4= (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DCum4" );
  TH1D *resvn2= (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_vn2" );
  TH1D *resvn4= (TH1D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_vn4" );
  //correlators
  TProfile *c2 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2" ));
  TProfile *c4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4" ));
  TProfile *dc2 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC2" ));
  TProfile *dc4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC4" ));
  TProfile *c2c4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2C4" ));
  TProfile *c2dc2 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2DC2" ));
  TProfile *c2dc4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2DC4" ));
  TProfile *c4dc2 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4DC2" ));
  TProfile *c4dc4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4DC4" ));
  TProfile *dc2dc4 = ((TProfile*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC2DC4" ));
  TArrayD *c2sww = c2->GetBinSumw2();
  TArrayD *c4sww = c4->GetBinSumw2();
  TArrayD *dc2sww= dc2->GetBinSumw2();
  TArrayD *dc4sww= dc4->GetBinSumw2();
  for(Int_t i=1; i!=resvn2->GetNbinsX()+1; ++i) {
    // cn{2}
    double v_c2sw = c2->GetBinEntries(i);
    if(v_c2sw<1e-100) continue;
    double v_c2sww = c2sww->At(i);
    double v_c2 = c2->GetBinContent(i);
    double e_c2 = TMath::Sqrt(v_c2sww)/v_c2sw*c2->GetBinError(i);
    double cum2 = v_c2;
    double ecum2= e_c2;
    resC2->SetBinContent(i, cum2 );
    resC2->SetBinError(i, ecum2 );
    // cn{4}
    double v_c4sw = c4->GetBinEntries(i);
    if(v_c4sw<1e-100) continue;
    double v_c4sww = c4sww->At(i);
    double v_c4 = c4->GetBinContent(i);
    double e_c4 = TMath::Sqrt(v_c4sww)/v_c4sw*c4->GetBinError(i);
    double v_c2c4 = c2c4->GetBinContent(i);
    double v_c2c4sw = c2c4->GetBinEntries(i);
    if(TMath::AreEqualAbs(v_c2c4sw/v_c2sw/v_c4sw,1,1e-100)) continue;
    double covc2c4 = v_c2c4sw/v_c2sw/v_c4sw*(v_c2c4 - v_c2*v_c4)/(1-v_c2c4sw/v_c2sw/v_c4sw);
    double cum4 = v_c4 - 2*v_c2*v_c2;
    double ecum4= 16.0*v_c2*v_c2*e_c2*e_c2 + e_c4*e_c4 - 8.0*v_c2*covc2c4;
    if(ecum4<1e-100) continue;
    ecum4 = TMath::Sqrt( ecum4 );
    resC4->SetBinContent(i, cum4 );
    resC4->SetBinError(i, ecum4 );
    // dn{2}
    double v_dc2sw = dc2->GetBinEntries(i);
    if(v_dc2sw<1) continue;
    double v_dc2 = dc2->GetBinContent(i);
    double v_dc2sww = dc2sww->At(i);
    double e_dc2 = TMath::Sqrt(v_dc2sww)/v_dc2sw*dc2->GetBinError(i);
    double dcum2 = v_dc2;
    double edcum2= e_dc2;
    resDC2->SetBinContent(i, dcum2 );
    resDC2->SetBinError(i, edcum2 );
    // v2{2}
    if(v_c2<1e-100) continue;
    double dv22 = v_dc2/TMath::Sqrt(v_c2);
    double v_c2dc2 = c2dc2->GetBinContent(i);
    double v_c2dc2sw = c2dc2->GetBinEntries(i);
    if(TMath::AreEqualAbs(v_c2dc2sw/v_c2sw/v_dc2sw,1,1e-100)) continue;
    double covc2dc2 = v_c2dc2sw/v_c2sw/v_dc2sw*(v_c2dc2 - v_c2*v_dc2)/(1-v_c2dc2sw/v_c2sw/v_dc2sw);
    double edv22 = 0.25/v_c2/v_c2/v_c2*(v_dc2*v_dc2*e_c2*e_c2 + 4*v_c2*v_c2*e_dc2*e_dc2 - 4*v_c2*v_dc2*covc2dc2);
    //printf("%d >> dv22 %.16f || edv22^2 %.16f |||| v_c2dc2 %.16f | v_c2dc2sw %.16f | covc2dc2 %.16f \n", i,dv22,edv22,v_c2dc2,v_c2dc2sw,covc2dc2);
    if(edv22<1e-100) continue;
    edv22 = TMath::Sqrt(edv22);
    resvn2->SetBinContent(i,dv22);
    resvn2->SetBinError(i,edv22);
    // dn{4}
    double v_dc4sw = dc4->GetBinEntries(i);
    if(v_dc4sw<1) continue;
    double v_dc4 = dc4->GetBinContent(i);
    double v_dc4sww = dc4sww->At(i);
    double e_dc4 = TMath::Sqrt(v_dc4sww)/v_dc4sw*dc4->GetBinError(i);
    double dcum4 = v_dc4 - 2*v_c2*v_dc2;
    double v_c2dc4 = c2dc4->GetBinContent(i);
    double v_c2dc4sw = c2dc4->GetBinEntries(i);
    if(TMath::AreEqualAbs(v_c2dc4sw/v_c2sw/v_dc4sw,1,1e-100)) continue;
    double covc2dc4 = v_c2dc4sw/v_c2sw/v_dc4sw*(v_c2dc4 - v_c2*v_dc4)/(1-v_c2dc4sw/v_c2sw/v_dc4sw);
    double v_dc2dc4 = dc2dc4->GetBinContent(i);
    double v_dc2dc4sw = dc2dc4->GetBinEntries(i);
    if(TMath::AreEqualAbs(v_dc2dc4sw/v_dc2sw/v_dc4sw,1,1e-100)) continue;
    double covdc2dc4 = v_dc2dc4sw/v_dc2sw/v_dc4sw*(v_dc2dc4 - v_dc2*v_dc4)/(1-v_dc2dc4sw/v_dc2sw/v_dc4sw);
    double edcum4= ( +4.0*v_dc2*v_dc2*e_c2*e_c2
		     +4.0*v_c2*v_c2*e_dc2*e_dc2
		     +e_dc4*e_dc4
		     +8.0*v_c2*v_dc2*covc2dc2
		     -4.0*v_dc2*covc2dc4
		     -4.0*v_c2*covdc2dc4 );
    if(edcum4<1e-100) continue;
    edcum4 = TMath::Sqrt(edcum4);
    resDC4->SetBinContent(i, dcum4 );
    resDC4->SetBinError(i, edcum4 );
    // v2{4}
    if(cum4>1e-100) continue;
    double dv24 = -dcum4/TMath::Power(-cum4,0.75);
    double dterm1 = 2*v_c2*v_c2*v_dc2 - 3*v_c2*v_dc4 + 2*v_c4*v_dc2;
    double dterm2 = 9.0/16.0*dcum4*dcum4;
    double dterm3 = 4.0*v_c2*v_c2*cum4*cum4;
    double dterm4 = cum4*cum4;
    double dterm5 = -3.0/2.0*dcum4*dterm1;
    double dterm6 = -4.0*v_c2*cum4*dterm1;
    double dterm7 = -2.0*cum4*dterm1;
    double dterm8 = 3.0*v_c2*cum4*dcum4;
    double dterm9 = 3.0/2.0*cum4*dcum4;
    double dterm10= 4*v_c2*cum4*cum4;
    double v_c4dc2 = c4dc2->GetBinContent(i);
    double v_c4dc2sw = c4dc2->GetBinEntries(i);
    if(TMath::AreEqualAbs(v_c4dc2sw/v_c4sw/v_dc2sw,1,1e-100)) continue;
    double covc4dc2 = v_c4dc2sw/v_c4sw/v_dc2sw*(v_c4dc2 - v_c4*v_dc2)/(1-v_c4dc2sw/v_c4sw/v_dc2sw);
    double v_c4dc4 = c4dc4->GetBinContent(i);
    double v_c4dc4sw = c4dc4->GetBinEntries(i);
    if(TMath::AreEqualAbs(v_c4dc4sw/v_c4sw/v_dc4sw,1,1e-100)) continue;
    double covc4dc4 = v_c4dc4sw/v_c4sw/v_dc4sw*(v_c4dc4 - v_c4*v_dc4)/(1-v_c4dc4sw/v_c4sw/v_dc4sw);
    double edv24= 1.0/TMath::Power(-cum4,3.5)*(+dterm1*dterm1*e_c2*e_c2
					       +dterm2*e_c4*e_c4
					       +dterm3*e_dc2*e_dc2
					       +dterm4*e_dc4*e_dc4
					       -dterm5*covc2c4
					       -dterm6*covc2dc2
					       +dterm7*covc2dc4
					       +dterm8*covc4dc2
					       -dterm9*covc4dc4
					       -dterm10*covdc2dc4);
    if(edv24<1e-100) continue;
    edv24 = TMath::Sqrt(edv24);
    resvn4->SetBinContent(i,dv24);
    resvn4->SetBinError(i,edv24);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ComputeDecayVn(TString name) {
  TProfile2D *uQa, *uQc, *qaqc, *uQaqaqc, *uQcqaqc;
  TArrayD *pasww, *pbsww, *pcsww;
  //ScalarProducr TPC
  printf("<<%s>> SP TPC\n",name.Data());
  uQa  = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCA" );
  uQc  = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCC" );
  qaqc = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_TPCATPCC" );
  uQaqaqc = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCATPCATPCC" );
  uQcqaqc = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uTPCCTPCATPCC" );
  pasww = uQa->GetBinSumw2();
  pbsww = uQc->GetBinSumw2();
  pcsww = qaqc->GetBinSumw2();
  //
  TH2D *sptpca = (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnTPCA" );
  TH2D *sptpcc = (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnTPCC" );
  TH2D *sptpcaa = (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnTPCAA" );
  //
  for(Int_t i=1; i!=sptpcaa->GetNbinsX()+1; ++i) {
    for(Int_t j=1; j!=sptpcaa->GetNbinsY()+1; ++j) {
      sptpcaa->SetBinContent(i,j,0);
      sptpcaa->SetBinError(i,j,0);
      sptpca->SetBinContent(i,j,0);
      sptpca->SetBinError(i,j,0);
      sptpcc->SetBinContent(i,j,0);
      sptpcc->SetBinError(i,j,0);
      double a = uQa->GetBinContent(i,j);
      double b = uQc->GetBinContent(i,j);
      double c = qaqc->GetBinContent(i,j);
      //if(TMath::AreEqualAbs(a,0,1e-100)) continue;
      //if(TMath::AreEqualAbs(b,0,1e-100)) continue;
      if(c<1e-100) {printf("skipping i=%d, j=%d due to c=%.16f\n",i,j,c); continue;}
      // nominal sptpca
      double vna = a/TMath::Sqrt(c);
      sptpca->SetBinContent(i,j,vna);
      // nominal sptpcc
      double vnc = b/TMath::Sqrt(c);
      sptpcc->SetBinContent(i,j,vnc);
      // nominal sptpc
      double vn = (vna + vnc)/2.0;
      sptpcaa->SetBinContent(i,j,vn);
      // errors
      int k = sptpcaa->GetBin(i,j);
      double asw = uQa->GetBinEntries(k);
      double bsw = uQc->GetBinEntries(k);
      double csw = qaqc->GetBinEntries(k);
      if(asw<1e-100||bsw<1e-100||csw<1e-100) {printf("skipping i=%d, j=%d due to asw=%f or bsw=%f or csw=%f\n",i,j,asw,bsw,csw); continue;}
      double asww = pasww->At(k);
      double bsww = pbsww->At(k);
      double csww = pcsww->At(k);
      if(asww<1e-100||bsww<1e-100||csww<1e-100) {printf("skipping i=%d, j=%d due to asww=%f or bsww=%f or csww=%f\n",i,j,asww,bsww,csww); continue;}
      if((1<1e-100+asww/asw/asw)||(1<1e-100+bsww/bsw/bsw)||(1<1e-100+csww/csw/csw)) {printf("skipping i=%d, j=%d due to COVa=%f or COVb=%f or COVc=%f\n",i,j,asww/asw/asw,bsww/bsw/bsw,csww/csw/csw); continue;}
      if(TMath::AreEqualAbs(asww,asw*asw,1e-100)||
	 TMath::AreEqualAbs(bsww,bsw*bsw,1e-100)||
	 TMath::AreEqualAbs(csww,csw*csw,1e-100)) {printf("skipping i=%d, j=%d due to funny coincidence\n",i,j); continue;}
      double ac = uQaqaqc->GetBinContent(i,j);
      double bc = uQcqaqc->GetBinContent(i,j);
      double acsw = uQaqaqc->GetBinEntries(k);
      double bcsw = uQcqaqc->GetBinEntries(k);
      double ea = uQa->GetBinError(i,j)*TMath::Sqrt(asww)/asw/TMath::Sqrt(1-asww/asw/asw);
      double eb = uQc->GetBinError(i,j)*TMath::Sqrt(bsww)/bsw/TMath::Sqrt(1-bsww/bsw/bsw);
      double ec = qaqc->GetBinError(i,j)*TMath::Sqrt(csww)/csw/TMath::Sqrt(1-csww/csw/csw);
      //printf("%d >> ea^2 %.16f |||| asww %.16f | asw %.16f | 1-asww/asw/asw %.16f \n", i,ea, asww, asw, 1-asww/asw/asw);
      //printf("%d >> eb^2 %.16f |||| bsww %.16f | bsw %.16f | 1-bsww/bsw/bsw %.16f \n", i,eb, bsww, bsw, 1-bsww/bsw/bsw);
      //printf("%d >> ec^2 %.16f |||| csww %.16f | csw %.16f | 1-csww/csw/csw %.16f \n", i,ec, csww, csw, 1-csww/csw/csw);
      double ebc = (bc-b*c)/(1-bcsw/bsw/csw)*bcsw/bsw/csw;
      double eac = (ac-a*c)/(1-acsw/asw/csw)*acsw/asw/csw;
      double evna = 1.0/TMath::Abs(c) * ( ea*ea + vna*vna/TMath::Abs(c)/4.0*ec*ec - vna/TMath::Sqrt(c)*eac );
      double evnc = 1.0/TMath::Abs(c) * ( eb*eb + vnc*vnc/TMath::Abs(c)/4.0*ec*ec - vnc/TMath::Sqrt(c)*ebc );
      //printf("%d >> evna^2 %.16f |||| ea %.16f | ec %.16f | eac %.16f | c %.16f\n", i,evna, ea, ec, eac,c);
      //printf("%d >> evnc^2 %.16f |||| eb %.16f | ec %.16f | ebc %.16f | c %.16f\n", i,evnc, eb, ec, ebc,c);
      if(evna>1e-100) evna = TMath::Sqrt(evna); else evna=0;
      if(evnc>1e-100) evnc = TMath::Sqrt(evnc); else evnc=0;
      sptpca->SetBinError(i,j,evna);
      sptpcc->SetBinError(i,j,evnc);
      sptpcaa->SetBinError(i,j,TMath::Sqrt(evna*evna+evnc*evnc)/2.0);
    }
  }
  //ScalarProduct VZE
  printf("<<%s>> SP VZE\n",name.Data());
  double cvzea2 = ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("ChiSquaredVZEA"))->GetBinContent( 1 );
  double cvzec2 = ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("ChiSquaredVZEC"))->GetBinContent( 1 );
  if( TMath::AreEqualAbs(cvzea2+cvzec2,0,1e-100) ) return;
  uQa  = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEA" );
  uQc  = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEC" );
  qaqc = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZEAVZEC" );
  uQaqaqc = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEAVZEAVZEC" );
  uQcqaqc = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZECVZEAVZEC" );
  pasww = uQa->GetBinSumw2();
  pbsww = uQc->GetBinSumw2();
  pcsww = qaqc->GetBinSumw2();
  //
  TProfile2D *qaqt = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZEATPC" );
  TProfile2D *qcqt = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_VZECTPC" );
  TProfile2D *uQauQc  = (TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_uVZEAuVZEC" );
  //
  TH2D *spvzea = (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnVZEA" );
  TH2D *spvzec = (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnVZEC" );
  TH2D *spvzega = (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnVZEGA" );
  TH2D *spvzewa = (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "SP_vnVZEWA" );
  for(Int_t i=1; i!=spvzewa->GetNbinsX()+1; ++i) {
    for(Int_t j=1; j!=spvzewa->GetNbinsY()+1; ++j) {
      spvzega->SetBinContent(i,j,0);
      spvzega->SetBinError(i,j,0);
      spvzewa->SetBinContent(i,j,0);
      spvzewa->SetBinError(i,j,0);
      spvzea->SetBinContent(i,j,0);
      spvzea->SetBinError(i,j,0);
      spvzec->SetBinContent(i,j,0);
      spvzec->SetBinError(i,j,0);
      double a = uQa->GetBinContent(i,j);
      double b = uQc->GetBinContent(i,j);
      double c = qaqc->GetBinContent(i,j);
      double at = qaqt->GetBinContent(i,j);
      double bt = qcqt->GetBinContent(i,j);
      if(TMath::AreEqualAbs(a,0,1e-100)) {printf("skipping A\n"); continue;}
      if(TMath::AreEqualAbs(b,0,1e-100)) {printf("skipping B\n"); continue;}
      if(TMath::AreEqualAbs(c,0,1e-100)) {printf("skipping C\n"); continue;}
      if(TMath::AreEqualAbs(at,0,1e-100)) {printf("skipping AT\n"); continue;}
      if(TMath::AreEqualAbs(bt,0,1e-100)) {printf("skipping CT\n"); continue;}
      // nominal spvzea
      double aa = c*at/bt;
      if(aa<1e-100) {printf("AA\n"); continue;}
      double vna = a/TMath::Sqrt(aa);
      spvzea->SetBinContent(i,j,vna);
      // nominal spvzec
      double bb = c*bt/at;
      if(bb<1e-100) {printf("BB\n"); continue;}
      double vnc = b/TMath::Sqrt(bb);
      spvzec->SetBinContent(i,j,vnc);
      //nominal spvzewa
      double vnwa = (cvzea2*vna + cvzec2*vnc) / (cvzea2+cvzec2);
      spvzewa->SetBinContent(i,j,vnwa);
      // nominal spvzega
      double vnga = a*b/c;
      if(vnga<1e-100) continue;
      vnga = TMath::Sqrt(vnga);
      spvzega->SetBinContent(i,j,vnga);
      // errors
      int k = spvzea->GetBin(i,j);
      double asw = uQa->GetBinEntries(k);
      double bsw = uQc->GetBinEntries(k);
      double csw = qaqc->GetBinEntries(k);
      if(asw<1e-100||bsw<1e-100||csw<1e-100) continue;
      double asww = pasww->At(k);
      double bsww = pbsww->At(k);
      double csww = pcsww->At(k);
      if(asww<1e-100||bsww<1e-100||csww<1e-100) continue;
      if((1<asww/asw/asw)||(1<bsww/bsw/bsw)||(1<csww/csw/csw)) continue;
      double ab = uQauQc->GetBinContent(i,j);
      double ac = uQaqaqc->GetBinContent(i,j);
      double bc = uQcqaqc->GetBinContent(i,j);
      double absw = uQauQc->GetBinEntries(k);
      double acsw = uQaqaqc->GetBinEntries(k);
      double bcsw = uQcqaqc->GetBinEntries(k);
      if(TMath::AreEqualAbs(1,absw/asw/bsw,1e-100)||TMath::AreEqualAbs(1,bcsw/bsw/csw,1e-100)||TMath::AreEqualAbs(1,acsw/asw/csw,1e-100)) continue;
      double ea = uQa->GetBinError(i,j)*TMath::Sqrt(asww)/asw/TMath::Sqrt(1-asww/asw/asw);
      double eb = uQc->GetBinError(i,j)*TMath::Sqrt(bsww)/bsw/TMath::Sqrt(1-bsww/bsw/bsw);
      double ec = qaqc->GetBinError(i,j)*TMath::Sqrt(csww)/csw/TMath::Sqrt(1-csww/csw/csw);
      double eab = (ab-a*b)/(1-absw/asw/bsw)*absw/asw/bsw;
      double ebc = (bc-b*c)/(1-bcsw/bsw/csw)*bcsw/bsw/csw;
      double eac = (ac-a*c)/(1-acsw/asw/csw)*acsw/asw/csw;
      double nc, nec, neac, nebc;
      nc = c*at/bt;
      nec = ec*at/bt;
      neac = eac*at/bt;
      nebc = ebc*at/bt;
      double evna = 1.0/TMath::Abs(nc*nc*nc) * ( nc*nc*ea*ea + a*a/4.0*nec*nec - a*TMath::Abs(nc)*neac*neac );
      nc = c*bt/at;
      nec = ec*bt/at;
      neac = eac*bt/at;
      nebc = ebc*bt/at;
      double evnc = 1.0/TMath::Abs(nc*nc*nc) * ( nc*nc*eb*eb + b*b/4.0*nec*nec - b*TMath::Abs(nc)*nebc*nebc );
      if(evna>1e-100) evna = TMath::Sqrt(evna); else evna=0;
      if(evnc>1e-100) evnc = TMath::Sqrt(evnc); else evnc=0;
      spvzea->SetBinError(i,j,evna);
      spvzec->SetBinError(i,j,evnc);
      double evnwa = TMath::Sqrt( cvzea2*cvzea2*evna*evna + cvzec2*cvzec2*evnc*evnc )/(cvzea2+cvzec2);
      spvzewa->SetBinError(i,j,evnwa);
      double evnga = 0.25/c/c * ( TMath::Abs(b*c/a)*ea*ea + TMath::Abs(a*c/b)*eb*eb + TMath::Abs(a*b/c)*ec*ec + 2*c*eab - 2*a*ebc - 2*b*eac );
      if(evnga>1e-100) evnga = TMath::Sqrt(evnga); else evnga=0;
      spvzega->SetBinError(i,j,evnga);
    }
  }
  printf("<<%s>> QC TPC\n",name.Data());
  //Qcumulants
  TH2D *resC2 = (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_Cum2" );
  TH2D *resC4 = (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_Cum4" );
  TH2D *resDC2= (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DCum2" );
  TH2D *resDC4= (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DCum4" );
  TH2D *resvn2= (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_vn2" );
  TH2D *resvn4= (TH2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_vn4" );
  //correlators
  TProfile2D *c2 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2" ));
  TProfile2D *c4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4" ));
  TProfile2D *dc2 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC2" ));
  TProfile2D *dc4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC4" ));
  TProfile2D *c2c4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2C4" ));
  TProfile2D *c2dc2 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2DC2" ));
  TProfile2D *c2dc4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C2DC4" ));
  TProfile2D *c4dc2 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4DC2" ));
  TProfile2D *c4dc4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_C4DC4" ));
  TProfile2D *dc2dc4 = ((TProfile2D*)((TList*)fList->FindObject(name.Data()))->FindObject( "QC_DC2DC4" ));
  TArrayD *c2sww = c2->GetBinSumw2();
  TArrayD *c4sww = c4->GetBinSumw2();
  TArrayD *dc2sww= dc2->GetBinSumw2();
  TArrayD *dc4sww= dc4->GetBinSumw2();
  for(Int_t i=1; i!=resvn2->GetNbinsX()+1; ++i) {
    for(Int_t j=1; j!=resvn2->GetNbinsY()+1; ++j) {
      // cn{2}
      int k = c2->GetBin(i,j);
      double v_c2sw = c2->GetBinEntries(k);
      if(v_c2sw<1e-100) continue;
      double v_c2sww = c2sww->At(k);
      double v_c2 = c2->GetBinContent(i,j);
      double e_c2 = TMath::Sqrt(v_c2sww)/v_c2sw*c2->GetBinError(i,j);
      double cum2 = v_c2;
      double ecum2= e_c2;
      resC2->SetBinContent(i,j, cum2 );
      resC2->SetBinError(i,j, ecum2 );
      // cn{4}
      double v_c4sw = c4->GetBinEntries(k);
      if(v_c4sw<1e-100) continue;
      double v_c4sww = c4sww->At(k);
      double v_c4 = c4->GetBinContent(i,j);
      double e_c4 = TMath::Sqrt(v_c4sww)/v_c4sw*c4->GetBinError(i,j);
      double v_c2c4 = c2c4->GetBinContent(i,j);
      double v_c2c4sw = c2c4->GetBinEntries(k);
      if(TMath::AreEqualAbs(v_c2c4sw/v_c2sw/v_c4sw,1,1e-100)) continue;
      double covc2c4 = v_c2c4sw/v_c2sw/v_c4sw*(v_c2c4 - v_c2*v_c4)/(1-v_c2c4sw/v_c2sw/v_c4sw);
      double cum4 = v_c4 - 2*v_c2*v_c2;
      double ecum4= 16.0*v_c2*v_c2*e_c2*e_c2 + e_c4*e_c4 - 8.0*v_c2*covc2c4;
      if(ecum4<1e-100) continue;
      ecum4 = TMath::Sqrt( ecum4 );
      resC4->SetBinContent(i,j, cum4 );
      resC4->SetBinError(i,j, ecum4 );
      // dn{2}
      double v_dc2sw = dc2->GetBinEntries(k);
      if(v_dc2sw<1) continue;
      double v_dc2 = dc2->GetBinContent(i,j);
      double v_dc2sww = dc2sww->At(k);
      double e_dc2 = TMath::Sqrt(v_dc2sww)/v_dc2sw*dc2->GetBinError(i,j);
      double dcum2 = v_dc2;
      double edcum2= e_dc2;
      resDC2->SetBinContent(i,j, dcum2 );
      resDC2->SetBinError(i,j, edcum2 );
      // v2{2}
      if(v_c2<1e-100) continue;
      double dv22 = v_dc2/TMath::Sqrt(v_c2);
      double v_c2dc2 = c2dc2->GetBinContent(i,j);
      double v_c2dc2sw = c2dc2->GetBinEntries(k);
      if(TMath::AreEqualAbs(v_c2dc2sw/v_c2sw/v_dc2sw,1,1e-100)) continue;
      double covc2dc2 = v_c2dc2sw/v_c2sw/v_dc2sw*(v_c2dc2 - v_c2*v_dc2)/(1-v_c2dc2sw/v_c2sw/v_dc2sw);
      double edv22 = 0.25/v_c2/v_c2/v_c2*(v_dc2*v_dc2*e_c2*e_c2 + 4*v_c2*v_c2*e_dc2*e_dc2 - 4*v_c2*v_dc2*covc2dc2);
      //printf("%d >> dv22 %.16f || edv22^2 %.16f |||| v_c2dc2 %.16f | v_c2dc2sw %.16f | covc2dc2 %.16f \n", i,dv22,edv22,v_c2dc2,v_c2dc2sw,covc2dc2);
      if(edv22<1e-100) continue;
      edv22 = TMath::Sqrt(edv22);
      resvn2->SetBinContent(i,j,dv22);
      resvn2->SetBinError(i,j,edv22);
      // dn{4}
      double v_dc4sw = dc4->GetBinEntries(k);
      if(v_dc4sw<1) continue;
      double v_dc4 = dc4->GetBinContent(i,j);
      double v_dc4sww = dc4sww->At(k);
      double e_dc4 = TMath::Sqrt(v_dc4sww)/v_dc4sw*dc4->GetBinError(i,j);
      double dcum4 = v_dc4 - 2*v_c2*v_dc2;
      double v_c2dc4 = c2dc4->GetBinContent(i,j);
      double v_c2dc4sw = c2dc4->GetBinEntries(k);
      if(TMath::AreEqualAbs(v_c2dc4sw/v_c2sw/v_dc4sw,1,1e-100)) continue;
      double covc2dc4 = v_c2dc4sw/v_c2sw/v_dc4sw*(v_c2dc4 - v_c2*v_dc4)/(1-v_c2dc4sw/v_c2sw/v_dc4sw);
      double v_dc2dc4 = dc2dc4->GetBinContent(i,j);
      double v_dc2dc4sw = dc2dc4->GetBinEntries(k);
      if(TMath::AreEqualAbs(v_dc2dc4sw/v_dc2sw/v_dc4sw,1,1e-100)) continue;
      double covdc2dc4 = v_dc2dc4sw/v_dc2sw/v_dc4sw*(v_dc2dc4 - v_dc2*v_dc4)/(1-v_dc2dc4sw/v_dc2sw/v_dc4sw);
      double edcum4= ( +4.0*v_dc2*v_dc2*e_c2*e_c2
		       +4.0*v_c2*v_c2*e_dc2*e_dc2
		       +e_dc4*e_dc4
		       +8.0*v_c2*v_dc2*covc2dc2
		       -4.0*v_dc2*covc2dc4
		       -4.0*v_c2*covdc2dc4 );
      if(edcum4<1e-100) continue;
      edcum4 = TMath::Sqrt(edcum4);
      resDC4->SetBinContent(i,j, dcum4 );
      resDC4->SetBinError(i,j, edcum4 );
      // v2{4}
      if(cum4>1e-100) continue;
      double dv24 = -dcum4/TMath::Power(-cum4,0.75);
      double dterm1 = 2*v_c2*v_c2*v_dc2 - 3*v_c2*v_dc4 + 2*v_c4*v_dc2;
      double dterm2 = 9.0/16.0*dcum4*dcum4;
      double dterm3 = 4.0*v_c2*v_c2*cum4*cum4;
      double dterm4 = cum4*cum4;
      double dterm5 = -3.0/2.0*dcum4*dterm1;
      double dterm6 = -4.0*v_c2*cum4*dterm1;
      double dterm7 = -2.0*cum4*dterm1;
      double dterm8 = 3.0*v_c2*cum4*dcum4;
      double dterm9 = 3.0/2.0*cum4*dcum4;
      double dterm10= 4*v_c2*cum4*cum4;
      double v_c4dc2 = c4dc2->GetBinContent(i,j);
      double v_c4dc2sw = c4dc2->GetBinEntries(k);
      if(TMath::AreEqualAbs(v_c4dc2sw/v_c4sw/v_dc2sw,1,1e-100)) continue;
      double covc4dc2 = v_c4dc2sw/v_c4sw/v_dc2sw*(v_c4dc2 - v_c4*v_dc2)/(1-v_c4dc2sw/v_c4sw/v_dc2sw);
      double v_c4dc4 = c4dc4->GetBinContent(i,j);
      double v_c4dc4sw = c4dc4->GetBinEntries(k);
      if(TMath::AreEqualAbs(v_c4dc4sw/v_c4sw/v_dc4sw,1,1e-100)) continue;
      double covc4dc4 = v_c4dc4sw/v_c4sw/v_dc4sw*(v_c4dc4 - v_c4*v_dc4)/(1-v_c4dc4sw/v_c4sw/v_dc4sw);
      double edv24= 1.0/TMath::Power(-cum4,3.5)*(+dterm1*dterm1*e_c2*e_c2
						 +dterm2*e_c4*e_c4
						 +dterm3*e_dc2*e_dc2
						 +dterm4*e_dc4*e_dc4
						 -dterm5*covc2c4
						 -dterm6*covc2dc2
						 +dterm7*covc2dc4
						 +dterm8*covc4dc2
						 -dterm9*covc4dc4
						 -dterm10*covdc2dc4);
      if(edv24<1e-100) continue;
      edv24 = TMath::Sqrt(edv24);
      resvn4->SetBinContent(i,j,dv24);
      resvn4->SetBinError(i,j,edv24);
    }
  }
}

//=======================================================================
void AliAnalysisTaskFlowStrange::OpenToyModel() {
  fList = new TList();
  fList->SetOwner();

  TList *tList;
  tList=new TList(); tList->SetName("ToyVn"); tList->SetOwner(); AddDecayVn(tList); fList->Add(tList);
  AddMakeQSpy();

  fRFPAminEta=-0.9;
  fRFPAmaxEta=0.0;
  fRFPCminEta=0.0;
  fRFPCmaxEta=+0.9;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::CloseToyModel() {
  ComputeChi2VZERO();
  ComputeDecayVn("ToyVn");
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeToyEvent(  Int_t seed, Int_t m_decay,Double_t v_decay,
						Double_t mass_decay_mu,Double_t mass_decay_sg,
						Int_t m_bgr,Double_t v_bgr,
						Int_t mtpc_a,Double_t v_tpca,Int_t mtpc_c,Double_t v_tpcc,
						Int_t mvze_a,Double_t v_vzea,Int_t mvze_c,Double_t v_vzec ) {
  gRandom->SetSeed( seed );
  // QVectors
  fMCEP = gRandom->Rndm()*TMath::Pi();
  TF1 tf1_tpca( "dphitpca", Form("1+2*%f*TMath::Cos(2*x)",v_tpca),0,TMath::TwoPi() );
  TF1 tf1_tpcc( "dphitpcc", Form("1+2*%f*TMath::Cos(2*x)",v_tpcc),0,TMath::TwoPi() );
  TF1 tf1_vzea( "dphivzea", Form("1+2*%f*TMath::Cos(2*x)",v_vzea),0,TMath::TwoPi() );
  TF1 tf1_vzec( "dphivzec", Form("1+2*%f*TMath::Cos(2*x)",v_vzec),0,TMath::TwoPi() );
  TF1 tf1_decay( "dphidecay", Form("1+2*%f*TMath::Cos(2*x)",v_decay),0,TMath::TwoPi() );
  TF1 tf1_bgr( "dphibgr", Form("1+2*%f*TMath::Cos(2*x)",v_bgr),0,TMath::TwoPi() );
  Double_t phi, eta;
  fQTPCACos=fQTPCASin=fQTPCA=0;
  fQTPCCCos=fQTPCCSin=fQTPCC=0;
  fQTPC2hCos=fQTPC2hSin=0;
  fQVZEACos=fQVZEASin=fQVZEA=0;
  fQVZECCos=fQVZECSin=fQVZEC=0;
  for(int m=0; m!=mtpc_a; ++m) {
    phi = tf1_tpca.GetRandom() + fMCEP;
    if(phi>TMath::TwoPi()) phi -= TMath::TwoPi();
    eta = gRandom->Rndm()*(fRFPAmaxEta-fRFPAminEta)+fRFPAminEta;
    fQTPCACos += TMath::Cos(fHarmonic*phi);
    fQTPCASin += TMath::Sin(fHarmonic*phi);
    fQTPCA += 1;
    fQTPC2hCos += TMath::Cos(2*fHarmonic*phi);
    fQTPC2hSin += TMath::Sin(2*fHarmonic*phi);
    fQTPCA_fID[m] = -99;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCAllPhiEta"))->Fill( phi, eta );
  }
  for(int m=0; m!=mtpc_c; ++m) {
    phi = tf1_tpcc.GetRandom() + fMCEP;
    if(phi>TMath::TwoPi()) phi -= TMath::TwoPi();
    eta = gRandom->Rndm()*(fRFPCmaxEta-fRFPCminEta)+fRFPCminEta;
    fQTPCCCos += TMath::Cos(fHarmonic*phi);
    fQTPCCSin += TMath::Sin(fHarmonic*phi);
    fQTPCC += 1;
    fQTPC2hCos += TMath::Cos(2*fHarmonic*phi);
    fQTPC2hSin += TMath::Sin(2*fHarmonic*phi); 
    fQTPCC_fID[m] = -99;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCAllPhiEta"))->Fill( phi, eta );
  }
  for(int m=0; m!=mvze_a; ++m) {
    phi = tf1_vzea.GetRandom() + fMCEP;
    if(phi>TMath::TwoPi()) phi -= TMath::TwoPi();
    eta = gRandom->Rndm()*2-3.5;
    fQVZEACos += TMath::Cos(fHarmonic*phi);
    fQVZEASin += TMath::Sin(fHarmonic*phi);
    fQVZEA += 1;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEAllPhiEta"))->Fill( phi, eta );
  }
  for(int m=0; m!=mvze_c; ++m) {
    phi = tf1_vzec.GetRandom() + fMCEP;
    if(phi>TMath::TwoPi()) phi -= TMath::TwoPi();
    eta = gRandom->Rndm()*2+2.5;
    fQVZECCos += TMath::Cos(fHarmonic*phi);
    fQVZECSin += TMath::Sin(fHarmonic*phi);
    fQVZEC += 1;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEAllPhiEta"))->Fill( phi, eta );
  }
  fQTPCA_nTracks = mtpc_a;
  fQTPCC_nTracks = mtpc_c;
  FillMakeQSpy();

  //decays
  double ptrange = fPtBinEdge[fPtBins] - fPtBinEdge[0];
  double pt, mass;
  for(int m=0; m!=m_decay; ++m) {
    phi = tf1_decay.GetRandom() + fMCEP;
    if(phi>TMath::TwoPi()) phi -= TMath::TwoPi();
    eta = gRandom->Rndm()*1.6-0.8;
    pt = gRandom->Rndm()*ptrange + fPtBinEdge[0];
    mass = gRandom->Gaus(mass_decay_mu,mass_decay_sg);
    FillDecayVn("ToyVn",mass,pt,phi,eta,+999,+999);
  }
  for(int m=0; m!=m_bgr; ++m) {
    phi = tf1_bgr.GetRandom() + fMCEP;
    if(phi>TMath::TwoPi()) phi -= TMath::TwoPi();
    eta = gRandom->Rndm()*1.6-0.8;
    pt = gRandom->Rndm()*ptrange + fPtBinEdge[0];
    mass = gRandom->Rndm()*(fMaxMass-fMinMass)+fMinMass;
    FillDecayVn("ToyVn",mass,pt,phi,eta,+999,+999);
  }
  QCStoreDecayVn("ToyVn");

}
