#include "AliAnalysisTaskTPCTOFCascade.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
//#include "AliMCParticle.h"
//#include "AliStack.h"
#include "AliPhysicsSelection.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"
#include "AliTOFcalib.h"
#include "AliTOFT0maker.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliESDtrack.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "AliAnalysisPIDCascadeTrack.h"
#include "AliAnalysisPIDCascadeParticle.h"
#include "AliAnalysisPIDCascadeEvent.h"
#include "TClonesArray.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "TTree.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliPIDResponse.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "TRandom.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliAnalysisPIDCascadeV0.h"
#include "AliAnalysisPIDCascade.h"
#include "AliAODVertex.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliVVZERO.h"
#include "AliVCluster.h"
#include "TMath.h"
#include "AliNeutralTrackParam.h"

ClassImp(AliAnalysisTaskTPCTOFCascade)
  
//_______________________________________________________
  
AliAnalysisTaskTPCTOFCascade::AliAnalysisTaskTPCTOFCascade() :
  AliAnalysisTaskSE("AnalysisResults"),
  fInitFlag(kFALSE),
  fMCFlag(kFALSE),
  fMCTuneFlag(kFALSE),
  fPbPbFlag(kFALSE),
  fVertexSelectionFlag(kFALSE),
  fPrimaryDCASelectionFlag(kFALSE),
  fPIDTree(0),
  fEvHist(0),
  fPIDResponse(0),
  fAnUtils(0),
  fRunNumber(0),
  fStartTime(0),
  fEndTime(0),
  fESDEvent(NULL),
  fMCEvent(NULL),
//fMCStack(NULL),
  fTrackCutsV0(NULL),
  fTrackCuts2010(NULL),
  fTrackCuts2011(NULL),
  fTrackCutsTPCRefit(NULL),
  fTrackCuts2011Sys(NULL),
  fESDpid(new AliESDpid()),
  fIsCollisionCandidate(kFALSE),
  fIsEventSelected(0),
  fIsPileupFromSPD(kFALSE),
  fHasVertex(kFALSE),
  fVertexZ(0.),
  fMCTimeZero(0.),
  fCentrality(NULL),
  fAnalysisEvent(new AliAnalysisPIDCascadeEvent()),
  fAnalysisTrackArray(new TClonesArray("AliAnalysisPIDCascadeTrack")),
  fAnalysisTrack(new AliAnalysisPIDCascadeTrack()),
  fAnalysisParticleArray(new TClonesArray("AliAnalysisPIDCascadeParticle")),
  fAnalysisParticle(new AliAnalysisPIDCascadeParticle()),
  fAnalysisV0TrackArray(new TClonesArray("AliAnalysisPIDCascadeV0")),
  fAnalysisV0Track(new AliAnalysisPIDCascadeV0()),
  fAnalysisCascadeTrackArray(new TClonesArray("AliAnalysisPIDCascade")),
  fAnalysisCascadeTrack(new AliAnalysisPIDCascade()),
  fTOFcalib(new AliTOFcalib()),
  fTOFT0maker(new AliTOFT0maker(fESDpid)),
  fTimeResolution(80.),
  fVertexCut(10.0),
  fRapidityCut(1.0),
  V0MBinCount(0),
  fHistoList(new TList()),
  fMCHistoList(new TList())
{
  /* 
   * default constructor 
   */
  fTrackCuts2010 = new AliESDtrackCuts("AliESDtrackCuts2010","AliESDtrackCuts2010");
  fTrackCuts2010 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE); //If not running, set to kFALSE;
  fTrackCuts2010->SetEtaRange(-0.8,0.8);
  fTrackCuts2011 = new AliESDtrackCuts("AliESDtrackCuts2011","AliESDtrackCuts2011");
  fTrackCuts2011 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); //If not running, set to kFALSE;
  fTrackCuts2011->SetEtaRange(-0.8,0.8);
  fTrackCutsTPCRefit = new AliESDtrackCuts("AliESDtrackCutsTPCRefit","AliESDtrackCutsTPCRefit");
  fTrackCutsTPCRefit = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); //If not running, set to kFALSE;
  fTrackCutsTPCRefit->SetRequireTPCRefit(kTRUE);
  fTrackCutsTPCRefit->SetEtaRange(-0.8,0.8);
  //Following is TC for systematics estimation. To compensate, once should probably reduce gamma DeltaM even further :)
  fTrackCuts2011Sys = new AliESDtrackCuts("AliESDtrackCuts2011Sys","AliESDtrackCuts2011Sys");
  fTrackCuts2011Sys = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); //If not running, set to kFALSE;                                                                                            
  fTrackCuts2011Sys->SetEtaRange(-0.8,0.8);
  fTrackCuts2011Sys->SetMinNCrossedRowsTPC(60);
  fTrackCuts2011Sys->SetMaxChi2PerClusterTPC(5);
  fTrackCuts2011Sys->SetMaxDCAToVertexZ(3);
  fTrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0", "AliESDtrackCutsV0");
  fTrackCutsV0 = AliESDtrackCuts::GetStandardV0DaughterCuts();
  fTrackCutsV0->SetEtaRange(-0.8,0.8);

}




AliAnalysisTaskTPCTOFCascade::AliAnalysisTaskTPCTOFCascade(Bool_t isMC) :
  AliAnalysisTaskSE("AnalysisResults"),
  fInitFlag(kFALSE),
  fMCFlag(kFALSE),
  fMCTuneFlag(kFALSE),
  fPbPbFlag(kFALSE),
  fVertexSelectionFlag(kFALSE),
  fPrimaryDCASelectionFlag(kFALSE),
  fPIDTree(0),
  fEvHist(0),
  fPIDResponse(0),
  fAnUtils(0),
  fRunNumber(0),
  fStartTime(0),
  fEndTime(0),
  fESDEvent(NULL),
  fMCEvent(NULL),
  //fMCStack(NULL),
  fTrackCutsV0(NULL),
  fTrackCuts2010(NULL),
  fTrackCuts2011(NULL),
  fTrackCutsTPCRefit(NULL),
  fTrackCuts2011Sys(NULL),
  fESDpid(new AliESDpid()),
  fIsCollisionCandidate(kFALSE),
  fIsEventSelected(0),
  fIsPileupFromSPD(kFALSE),
  fHasVertex(kFALSE),
  fVertexZ(0.),
  fMCTimeZero(0.),
  fCentrality(NULL),
  fAnalysisEvent(new AliAnalysisPIDCascadeEvent()),
  fAnalysisTrackArray(new TClonesArray("AliAnalysisPIDCascadeTrack")),
  fAnalysisTrack(new AliAnalysisPIDCascadeTrack()),
  fAnalysisParticleArray(new TClonesArray("AliAnalysisPIDCascadeParticle")),
  fAnalysisParticle(new AliAnalysisPIDCascadeParticle()),
  fAnalysisV0TrackArray(new TClonesArray("AliAnalysisPIDCascadeV0")),
  fAnalysisV0Track(new AliAnalysisPIDCascadeV0()),
  fAnalysisCascadeTrackArray(new TClonesArray("AliAnalysisPIDCascade")),
  fAnalysisCascadeTrack(new AliAnalysisPIDCascade()),
  fTOFcalib(new AliTOFcalib()),
  fTOFT0maker(new AliTOFT0maker(fESDpid)),
  fTimeResolution(80.),
  fVertexCut(10.0),
  fRapidityCut(1.0),
  V0MBinCount(0),
  fHistoList(new TList()),
  fMCHistoList(new TList())
{
  /* 
   * default constructor 
   */
  fTrackCuts2010 = new AliESDtrackCuts("AliESDtrackCuts2010","AliESDtrackCuts2010");
  fTrackCuts2010 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kFALSE); //If not running, set to kFALSE;
  fTrackCuts2010->SetEtaRange(-0.8,0.8);
  fTrackCuts2011 = new AliESDtrackCuts("AliESDtrackCuts2011","AliESDtrackCuts2011");
  fTrackCuts2011 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); //If not running, set to kFALSE;
  fTrackCuts2011->SetEtaRange(-0.8,0.8);
  fTrackCutsTPCRefit = new AliESDtrackCuts("AliESDtrackCutsTPCRefit","AliESDtrackCutsTPCRefit");
  fTrackCutsTPCRefit = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts(); //If not running, set to kFALSE;
  fTrackCutsTPCRefit->SetRequireTPCRefit(kTRUE);
  fTrackCutsTPCRefit->SetEtaRange(-0.8,0.8);
  //Following is TC for systematics estimation. To compensate, once should probably reduce gamma DeltaM even further :)                                                                                    
  fTrackCuts2011Sys = new AliESDtrackCuts("AliESDtrackCuts2011Sys","AliESDtrackCuts2011Sys");
  fTrackCuts2011Sys = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); //If not running, set to kFALSE;                                                                                            
  fTrackCuts2011Sys->SetEtaRange(-0.8,0.8);
  fTrackCuts2011Sys->SetMinNCrossedRowsTPC(60);
  fTrackCuts2011Sys->SetMaxChi2PerClusterTPC(5);
  fTrackCuts2011Sys->SetMaxDCAToVertexZ(3);
  fTrackCutsV0 = new AliESDtrackCuts("AliESDtrackCutsV0", "AliESDtrackCutsV0");
  fTrackCutsV0 = AliESDtrackCuts::GetStandardV0DaughterCuts();
  fTrackCutsV0->SetEtaRange(-0.8,0.8);

  fMCFlag = isMC;
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TH1D::Class());
}













//_______________________________________________________

AliAnalysisTaskTPCTOFCascade::~AliAnalysisTaskTPCTOFCascade()
{
  /*
   * default destructor
   */
  if (fTrackCutsV0) delete fTrackCutsV0;
  if (fTrackCuts2010) delete fTrackCuts2010;
  if (fTrackCuts2011) delete fTrackCuts2011;
  if (fTrackCutsTPCRefit) delete fTrackCutsTPCRefit;
  if (fTrackCuts2011Sys) delete fTrackCuts2011Sys;
  delete fESDpid;
  delete fTOFcalib;
  delete fTOFT0maker;
  delete fHistoList;
  delete fMCHistoList;
}


//________________________________________________________________________

void
AliAnalysisTaskTPCTOFCascade::UserCreateOutputObjects()
{
  /*
   * user create output objects
   */
  OpenFile(1);
  OpenFile(2);
  /* output tree */
  fPIDTree = new TTree("PIDTree","PIDTree");
  fPIDTree->Branch("AnalysisEvent", "AliAnalysisPIDCascadeEvent", &fAnalysisEvent);  
  fPIDTree->Branch("AnalysisTrack", "TClonesArray", &fAnalysisTrackArray); 
  fPIDTree->Branch("AnalysisV0Track","TClonesArray",&fAnalysisV0TrackArray);
  fPIDTree->Branch("AnalysisCascadeTrack","TClonesArray",&fAnalysisCascadeTrackArray);
  if (fMCFlag)
    fPIDTree->Branch("AnalysisParticle", "TClonesArray", &fAnalysisParticleArray);


  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse(); 
  fAnUtils = new AliAnalysisUtils();
  Double_t V0MbinsDefault[13] = {0, 0.01, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
  V0MBinCount = new TH1F("V0MBinCount","V0MBinCount",12,V0MbinsDefault);

  fEvHist = new TH1D("StatHist","StatHist",10,-0.5,9.5);
  PostData(1,fPIDTree);
  PostData(2,fEvHist);
}

//_______________________________________________________





Bool_t
AliAnalysisTaskTPCTOFCascade::InitRun()
{
  /*
   * init run. 
   */

  /* get ESD event */
  fESDEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESDEvent) {
    AliError("cannot get ESD event");
    return kFALSE;
  }

  /* get run number */
  Int_t runNb = fESDEvent->GetRunNumber();
  /* check run already initialized */
  if (fInitFlag && fRunNumber == runNb) return kTRUE;
  /* init cdb */
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(runNb);
  /* init TOF calib */
  if (!fTOFcalib->Init()) {
    AliError("cannot init TOF calib");
    return kFALSE;
  }
  AliInfo(Form("initialized for run %d", runNb));
  fInitFlag = kTRUE;
  fRunNumber = runNb;
  return kTRUE;
}

//_______________________________________________________
void AliAnalysisTaskTPCTOFCascade::FillHist(Double_t myflag) {
  fEvHist->Fill(myflag);
};
Bool_t AliAnalysisTaskTPCTOFCascade::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
};
Bool_t AliAnalysisTaskTPCTOFCascade::SelectVertex2015pp(AliESDEvent *esd,  Bool_t checkSPDres, Bool_t *SPDandTrkExists, Bool_t *PassProximityCut) 
{
  if (!esd) return kFALSE;
  const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
  const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
  Bool_t hasSPD = spdVertex->GetStatus();
  Bool_t hasTrk = trkVertex->GetStatus();
  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  //MOD: do not reject if SPD&Trk vtx. not there, but store it to the variable, if requested:
  if(SPDandTrkExists)
    (*SPDandTrkExists) = hasSPD&&hasTrk;
  //Set initial value for proximity check. Only checking the proximity if SPDandTrkExists,
  //the default value should be 1, tracking vtx is not required
  if(PassProximityCut)
    (*PassProximityCut) = 1;
  
  //reject events if none between the SPD or track verteces are available
  //if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (!hasSPD) return kFALSE;
    //on demand check the spd vertex resolution and reject if not satisfied
    if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
  } else {
    if (hasSPD) {
      //if enabled check the spd vertex resolution and reject if not satisfied
      //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
      if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
      if(PassProximityCut)
	(*PassProximityCut) = TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())<=0.5;
  //if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE; 
    }
  }

  //Not needed here, done separately
  //Cut on the vertex z position
  //const AliESDVertex * vertex = esd->GetPrimaryVertex();
  //if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
  return kTRUE;
};


Bool_t
AliAnalysisTaskTPCTOFCascade::InitEvent()
{
  /*
   * init event
   */

  /* get ESD event */
  fESDEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  FillHist(0);
  if (!fESDEvent) return kFALSE;
  /* get MC event */
  FillHist(1);
  if (fMCFlag) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCEvent) return kFALSE;
  }
  /* get stack */
  //Stack is gone. Farewll, stack, you've served us well.
  /*  if (fMCFlag) {
    fMCStack = fMCEvent->Stack();
    if (!fMCStack) return kFALSE;
    }*/
  FillHist(2);
  /* event selection */
  fIsCollisionCandidate = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny);
  fIsEventSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  fIsPileupFromSPD = fESDEvent->IsPileupFromSPDInMultBins();
  FillHist(3);  
  if(fESDEvent->IsIncompleteDAQ()) return kFALSE;
  if(fAnUtils->IsSPDClusterVsTrackletBG(fESDEvent)) return kFALSE;
  /* vertex selection */
  const AliESDVertex *vertex = fESDEvent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1) {
    vertex = fESDEvent->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1) fHasVertex = kFALSE;
    else fHasVertex = kTRUE;
    TString vtxTyp = vertex->GetTitle();
    Double_t cov[6]={0};
    vertex->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) fHasVertex = kFALSE;
  }
  else fHasVertex = kTRUE;

  
  fVertexZ = vertex->GetZ();
  /* centrality in PbPb */
  if (fPbPbFlag) {
    fCentrality = fESDEvent->GetCentrality();
  }
  /* calibrate ESD (also in MC to correct TExp) */
  fTOFcalib->CalibrateESD(fESDEvent);

  /* check MC flags */
  if (fMCFlag && fMCTuneFlag)
    fMCTimeZero = fTOFcalib->TuneForMC(fESDEvent, fTimeResolution);

#if 0 /* NEW METHOD */
  /* compute TOF-T0, fill ESD and make PID */
  fTOFT0maker->ComputeT0TOF(fESDEvent);
  fTOFT0maker->WriteInESD(fESDEvent);
  fESDpid->SetTOFResponse(fESDEvent, AliESDpid::kTOF_T0);
  fESDpid->MakePID(fESDEvent, kFALSE, 0.);
#endif

#if 1 /* OLD METHOD */
  /* compute and apply TOF-T0 */
  fTOFT0maker->ComputeT0TOF(fESDEvent);
  fESDpid->MakePID(fESDEvent, kFALSE, 0.);
#endif
  
  return kTRUE;
}

//_______________________________________________________

Bool_t
AliAnalysisTaskTPCTOFCascade::HasPrimaryDCA(AliESDtrack *track)
{
  /*
   * has primary DCA
   */
  
  // cut on transverse impact parameter
  Float_t d0z0[2],covd0z0[3];
  track->GetImpactParameters(d0z0, covd0z0);
  Float_t sigma= 0.0050 + 0.0060 / TMath::Power(track->Pt(), 0.9);
  Float_t d0max = 7. * sigma;
  //
  Float_t sigmaZ = 0.0146 + 0.0070 / TMath::Power(track->Pt(), 1.114758);
  if (track->Pt() > 1.) sigmaZ = 0.0216;
  Float_t d0maxZ = 5. * sigmaZ;
  //
  if(TMath::Abs(d0z0[0]) > d0max || TMath::Abs(d0z0[1]) > d0maxZ)
    return kFALSE;
  
  /* primary DCA ok */
  return kTRUE;
}
Int_t AliAnalysisTaskTPCTOFCascade::GetTrackCutsFlag(AliESDtrack *LocalTrack) {
  Int_t ReturnFlag = 0;
  if(fTrackCuts2010->AcceptTrack(LocalTrack)) ReturnFlag+=1;
  if(fTrackCuts2011->AcceptTrack(LocalTrack)) ReturnFlag+=2;
  if(fTrackCutsTPCRefit->AcceptTrack(LocalTrack)) ReturnFlag+=4;
  if(fTrackCuts2011Sys->AcceptTrack(LocalTrack)) ReturnFlag+=8;
  if(fTrackCutsV0->AcceptTrack(LocalTrack)) ReturnFlag+=16;

  return ReturnFlag;
};
//_______________________________________________________
void AliAnalysisTaskTPCTOFCascade::ProcessV0s() {
  Int_t NV0s = fESDEvent->GetNumberOfV0s();
  if(NV0s<1) return;
  const AliESDVertex *BestPrimaryVertex = fESDEvent->GetPrimaryVertex();
  if(!BestPrimaryVertex) return;
  if(!(BestPrimaryVertex->GetStatus())) return;
  Double_t IPrimaryVtxPosition[3];
  BestPrimaryVertex->GetXYZ(IPrimaryVtxPosition);
  Double_t IPrimaryVtxCov[6];
  BestPrimaryVertex->GetCovMatrix(IPrimaryVtxCov);
  Double_t IPrimaryVtxChi2 = BestPrimaryVertex->GetChi2toNDF();
  AliAODVertex *PrimaryVertex = new AliAODVertex(IPrimaryVtxPosition,IPrimaryVtxCov,IPrimaryVtxChi2,NULL,-1,AliAODVertex::kPrimary);

  /*Calculate DCA to prim. vertex
    Taken from AliAnalysisVertexingHF*/
  /*AliAODVertex *tmpVtx = new AliAODVertex(IPrimaryVtxPosition,BestPrimaryVertex->GetChi2V0(),AliAODVertex::kV0, 2);
  Double_t xyz[3], pxpypz[3];
  BestPrimaryVert*/

  /*Done w/ calculation*/


  Double_t InvMasses[4];
  for(Int_t iV0=0;iV0<NV0s;iV0++) {

    AliESDv0 *V0Vertex = fESDEvent->GetV0(iV0);
    if(!V0Vertex) 
      continue;
    if(V0Vertex->GetV0CosineOfPointingAngle()<0.95)
      continue;
    if(V0Vertex->GetDcaV0Daughters()>2.0)
      continue;
    if(TMath::Abs(V0Vertex->Eta())>0.8)
      continue;
    if(V0Vertex->GetOnFlyStatus()) 
      continue;
    Double_t IV0Position[3];
    V0Vertex->GetXYZ(IV0Position[0],IV0Position[1],IV0Position[2]);
    const Double_t IV0Radius = TMath::Sqrt(IV0Position[0]*IV0Position[0] + IV0Position[1]*IV0Position[1]);
    if(IV0Radius>200||IV0Radius<0.1) 
      continue;

    //fAnalysisTrack->Update(track, fMCStack, fMCEvent,fPIDResponse, fTrackCuts->AcceptTrack(track));
    AliESDtrack *pEsdTrack = fESDEvent->GetTrack((UInt_t)TMath::Abs(V0Vertex->GetPindex()));
    AliESDtrack *nEsdTrack = fESDEvent->GetTrack((UInt_t)TMath::Abs(V0Vertex->GetNindex()));
    
    if(!pEsdTrack||!nEsdTrack) 
      continue;

    if(!(pEsdTrack->IsOn(AliESDtrack::kTPCrefit)))
      continue;
    if(!(nEsdTrack->IsOn(AliESDtrack::kTPCrefit)))
      continue;
    
    if(TMath::Abs(pEsdTrack->Eta())>0.8||(TMath::Abs(nEsdTrack->Eta())>0.8))
      continue;
    if(pEsdTrack->Pt()<0.15||nEsdTrack->Pt()<0.15)
      continue;

    if(pEsdTrack->GetSign()==nEsdTrack->GetSign()) continue; //Remove like-sign
    //if(TMath::Abs(pTrack->GetEta())>0.8 || TMath::Abs(nTrack->GetEta())>0.8) continue; //Eta cut
    //if(pTrack->GetPt()<2) continue; //pT cut on decay product
    //if(nTrack->GetPt()<2) continue; //pT cut on decay product
    Bool_t ChargesSwitched=kFALSE;
    if(pEsdTrack->GetSign()<0) {
      AliESDtrack *ttr = nEsdTrack;
      nEsdTrack = pEsdTrack;//fESDEvent->GetTrack((UInt_t)TMath::Abs(V0Vertex->GetPindex()));
      pEsdTrack = ttr;//nTrack;//fESDEvent->GetTrack((UInt_t)TMath::Abs(V0Vertex->GetNindex()));
      ChargesSwitched=kTRUE;
    };
    //Double_t alpha = V0Vertex->AlphaV0(); //Probably save these
    //Double_t ptarm = V0Vertex->PtArmV0(); //Probably save these
    AliKFVertex PrimaryVtxKF(*PrimaryVertex);
    AliKFParticle::SetField(fESDEvent->GetMagneticField());


    Int_t indecies[2] = {211, 2212}; //pi,p
    Double_t myMasses[] = {0.497614, 1.11568, 1.11568}; //K0s, lambda, anti-lambda
    AliKFParticle *negKF[2] = {0,0}; //-pi, -p
    AliKFParticle *posKF[2] = {0,0}; // pi,  p
    if(ChargesSwitched)
      for(Int_t i=0;i<2;i++) {
	negKF[i] = new AliKFParticle(*(V0Vertex->GetParamP()), -indecies[i]);
	posKF[i] = new AliKFParticle(*(V0Vertex->GetParamN()), indecies[i]);
      } 
    else
      for(Int_t i=0;i<2;i++) {
	negKF[i] = new AliKFParticle(*(V0Vertex->GetParamN()), -indecies[i]);
	posKF[i] = new AliKFParticle(*(V0Vertex->GetParamP()), indecies[i]);
      };
    AliKFParticle V0KFs[3];
    Bool_t TrashTracks=kTRUE; //Trash tracks if all below inv. mass
    for(Int_t i=0;i<3;i++) { 
      V0KFs[i]+=(*posKF[(i==1)?1:0]);
      V0KFs[i]+=(*negKF[(i==2)?1:0]);
      V0KFs[i].SetProductionVertex(PrimaryVtxKF);
      InvMasses[i] = V0KFs[i].GetMass()-myMasses[i];
      if(i==0)
	TrashTracks = TrashTracks&&(TMath::Abs(InvMasses[i])>0.1);
      else
	TrashTracks = TrashTracks&&(TMath::Abs(InvMasses[i])>0.04);
    };
    for(Int_t i=0;i<2;i++) {

      delete negKF[i];
      delete posKF[i];
    }
    

    if(TrashTracks) 
      continue;

    const Double_t lpT = V0Vertex->Pt();
    const Double_t lEta = V0Vertex->Eta();
    
    /*Calculate DCA to prim. vertex
      Taken from AliAnalysisVertexingHF*/
    Double_t lDCAtoPrim = -999;
    Double_t xyz[3], pxpypz[3];
    V0Vertex->XvYvZv(xyz);
    V0Vertex->PxPyPz(pxpypz);
    Double_t cv[21]; 
    for(Int_t i=0;i<21;i++) 
      cv[i] = 0;

    AliNeutralTrackParam *trackesdV0 = new AliNeutralTrackParam(xyz,pxpypz,cv,0);
    if(!trackesdV0) 
      lDCAtoPrim=-999; 
    else {
      Double_t d0[2], covd0[3];
      trackesdV0->PropagateToDCA(PrimaryVertex,fESDEvent->GetMagneticField(),kVeryBig,d0,covd0);
      lDCAtoPrim = TMath::Sqrt(covd0[0]);
    };
    delete trackesdV0;

    /*AliAODVertex *tmpVtx = new AliAODVertex(IPrimaryVtxPosition,BestPrimaryVertex->GetChi2V0(),AliAODVertex::kV0, 2);
      Double_t xyz[3], pxpypz[3];
      BestPrimaryVert*/

    /*Done w/ calculation*/
    AliAnalysisPIDCascadeV0* v0Tree = new ((*fAnalysisV0TrackArray)[fAnalysisV0TrackArray->GetEntries()]) AliAnalysisPIDCascadeV0();
    
    AliAnalysisPIDCascadeTrack *pTrack = new AliAnalysisPIDCascadeTrack();
    AliAnalysisPIDCascadeTrack *nTrack = new AliAnalysisPIDCascadeTrack();
    pTrack->Update(pEsdTrack, fMCEvent, fPIDResponse, GetTrackCutsFlag(pEsdTrack));
    nTrack->Update(nEsdTrack, fMCEvent, fPIDResponse, GetTrackCutsFlag(nEsdTrack));

    v0Tree->Update(pTrack, nTrack, InvMasses, IV0Radius, V0Vertex->GetDcaV0Daughters(), V0Vertex->GetV0CosineOfPointingAngle(), lpT, lEta, lDCAtoPrim);
  };

};
//_______________________________________________________
void AliAnalysisTaskTPCTOFCascade::ProcessCascades() {
  Int_t NCascades = fESDEvent->GetNumberOfCascades();
  if(NCascades<1) return;
  const AliESDVertex *BestPrimaryVertex = fESDEvent->GetPrimaryVertex();
  if(!BestPrimaryVertex) return;
  if(!BestPrimaryVertex->GetStatus()) return;
  Double_t lPrimaryVtxPosition[3] = {-100.0 , -100.0, -100.0};
  BestPrimaryVertex->GetXYZ(lPrimaryVtxPosition);

  Double_t InvMasses[2] = {0};
  Double_t InvV0Masses[4] = {0};
  Double_t ArrayCounter = 0;
  for(Int_t iCasc = 0; iCasc<NCascades; iCasc++) {
    AliESDcascade* casc = (AliESDcascade*)fESDEvent->GetCascade(iCasc);
    if(!casc)
      continue;
/////////////////////////////////////////////////////////////////////////////
    // Begin to evaluate Dca, Pointing Angle & Spatial position of Cascade
    Double_t lDcaXiDaughters          = casc->GetDcaXiDaughters();
    Double_t lXiCosineOfPointingAngle = casc->GetCascadeCosineOfPointingAngle(lPrimaryVtxPosition[0], lPrimaryVtxPosition[1], lPrimaryVtxPosition[2]);   
    //Double_t lDcaXiToPrimVertex = casc->GetD(lPrimaryVtxPosition[0], lPrimaryVtxPosition[1], lPrimaryVtxPosition[2]);

    Double_t lPosXi[3] = {-1000.0, -1000.0, -1000.0};
    casc->GetXYZcascade(lPosXi[0], lPosXi[1], lPosXi[2]);
/////////////////////////////////////////////////////////////////////////////
    // Begin to evaluate Dca, Pointing Angle & Spatial position of V0
    Double_t lDcaV0DaughtersXi = casc->GetDcaV0Daughters();
    Double_t lV0toXiCosineOfPointingAngle = casc->GetV0CosineOfPointingAngle(lPosXi[0], lPosXi[1], lPosXi[2]);
    Double_t lDcaXiToPrimVertex = casc->GetD(lPosXi[0], lPosXi[1], lPosXi[2]);
    // Double_t lDcaV0ToPrimVertexXi = casc->GetD(lPosXi[0], lPosXi[1], lPosXi[2]);
    Double_t lPosV0[3] = {-1000.0, -1000.0, -1000.0};
    casc->GetXYZ(lPosV0[0], lPosV0[1], lPosV0[2]);
    Double_t lDcaV0ToPrimVertexXi = casc->GetD(lPosV0[0], lPosV0[1], lPosV0[2]);
    Double_t lDcaV0ToPrimVertex = casc->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
/////////////////////////////////////////////////////////////////////////////
    //Calculate Radius for both V0 and Cascade
    Double_t lXiRadius = TMath::Sqrt(lPosXi[0]*lPosXi[0] + lPosXi[1]*lPosXi[1]); 
    Double_t lV0Radius = TMath::Sqrt(lPosV0[0]*lPosV0[0] + lPosV0[1]*lPosV0[1]);
    /* ALL DCA & VERTEX CUTS ARE TO BE DONE IN POST-PROCESSING  (NOPE; WE HAVE SOME CUTS HERE)*/
    
/////////////////////////////////////////////////////////////////////////////
    //Defining Mass Variables
    Double_t lInvMassLambdaAsCascDghter = casc->GetEffMass(); //Lambda Mass
    InvV0Masses[1] = lInvMassLambdaAsCascDghter;
    Double_t lInvMassXi    = 0;
    Double_t lInvMassOmega = 0;
    Double_t lV0quality = 0.;
/////////////////////////////////////////////////////////////////////////////
    //Reconstructing & Check bachelor & V0 daughter tracks
    AliESDtrack* temp_b = fESDEvent->GetTrack((UInt_t)TMath::Abs(casc->GetBindex()));
    AliESDtrack* temp_p = fESDEvent->GetTrack((UInt_t)TMath::Abs(casc->GetPindex()));
    AliESDtrack* temp_n = fESDEvent->GetTrack((UInt_t)TMath::Abs(casc->GetNindex()));
    if (!temp_b || !temp_p || !temp_n)
      continue;
    UInt_t lIDtemp_p  = (UInt_t) TMath::Abs(temp_p->GetID());
    UInt_t lIDtemp_n  = (UInt_t) TMath::Abs(temp_n->GetID());
    UInt_t lIDtemp_b   = (UInt_t) TMath::Abs(temp_b->GetID());
    if(lIDtemp_b == lIDtemp_n)
      continue;
    if(lIDtemp_b == lIDtemp_p)
      continue;
    //Cuts

    if (!(temp_p->IsOn(AliESDtrack::kTPCrefit)))
      continue;
    if (!(temp_n->IsOn(AliESDtrack::kTPCrefit)))
      continue;
    if (!(temp_b->IsOn(AliESDtrack::kTPCrefit)))
      continue;
    if (temp_p->GetTPCNcls()<70||temp_n->GetTPCNcls()<70||temp_b->GetTPCNcls()<70)
      continue;
    
    if (TMath::Abs(temp_p->Eta())>0.8||(TMath::Abs(temp_n->Eta())>0.8||TMath::Abs(temp_b->Eta())>0.8))
      continue;
    if (temp_p->Pt()<0.15||temp_n->Pt()<0.15||temp_b->Pt()<0.15)
      continue;
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
    if(temp_b->GetSign()<0){
      casc->ChangeMassHypothesis(lV0quality, 3312); //Xi-
      lInvMassXi    = casc->GetEffMassXi();
      casc->ChangeMassHypothesis(lV0quality, 3334); //Omega-
      lInvMassOmega = casc->GetEffMassXi();
    }
    if(temp_b->GetSign()>0)
      {
	casc->ChangeMassHypothesis(lV0quality, -3312); //Xi+
	lInvMassXi    = casc->GetEffMassXi();
	casc->ChangeMassHypothesis(lV0quality, -3334);// Omega+
	lInvMassOmega = casc->GetEffMassXi();
      }

    if(TMath::Abs(lInvMassXi - 1.32171) > 0.1 && TMath::Abs(lInvMassOmega - 1.6725) > 0.1 )
      continue;

    InvMasses[1]  = lInvMassOmega;
    InvMasses[0]  = lInvMassXi;
      
  
    Double_t V0Px = casc->AliESDv0::Px();
    Double_t V0Py = casc->AliESDv0::Py();
    Double_t V0Pz = casc->AliESDv0::Pz();
    Double_t V0P  = casc->AliESDv0::P();
    Double_t CascPx = casc->AliESDcascade::Px();
    Double_t CascPy = casc->AliESDcascade::Py();
    Double_t CascPz = casc->AliESDcascade::Pz();
    Double_t CascP  = casc->AliESDcascade::P();
    Double_t V0P_Pz = V0P - V0Pz + 1.e-13;
    Double_t CascP_Pz = CascP - CascPz + 1.e-13;
    Double_t V0P_plus_Pz =  V0P + V0Pz;
    Double_t CascP_plus_Pz = CascP + CascPz;
    if(V0Pz>V0P)
      continue;
    if(CascPz>CascP)
      continue;
    if(V0P_plus_Pz/V0P_Pz <= 0)
      continue;
    if(CascP_plus_Pz/CascP_Pz <= 0)
      continue;
    Double_t V0Pt = TMath::Sqrt(V0Px*V0Px + V0Py*V0Py);
    Double_t V0Eta = 0.5*TMath::Log(TMath::Abs((V0P_plus_Pz)/(V0P_Pz)));
    if(TMath::Abs(V0Eta)>0.8)
      continue;
    Double_t CascPt = TMath::Sqrt(CascPx*CascPx + CascPy*CascPy);
    Double_t CascEta = 0.5*TMath::Log(TMath::Abs((CascP_plus_Pz)/(CascP_Pz)));
    if(TMath::Abs(CascEta)>0.8)
      continue;
/////////////////////////////////////////////////////////////////////////////
    //Cascade is good! Bag & Tag.
     AliAnalysisPIDCascade* Cascade = new ((*fAnalysisCascadeTrackArray)[ArrayCounter]) AliAnalysisPIDCascade();
     ArrayCounter++;
        
     AliAnalysisPIDCascadeTrack* bTrack = (AliAnalysisPIDCascadeTrack*)Cascade->GetBachAnalysisTrack();
     AliAnalysisPIDCascadeV0* V0 = (AliAnalysisPIDCascadeV0*)Cascade->GetV0();
     
     AliAnalysisPIDCascadeTrack* pTrack = new AliAnalysisPIDCascadeTrack();
     AliAnalysisPIDCascadeTrack* nTrack = new AliAnalysisPIDCascadeTrack();

     bTrack->Update(temp_b,fMCEvent,fPIDResponse, GetTrackCutsFlag(temp_b));
     pTrack->Update(temp_p,fMCEvent,fPIDResponse, GetTrackCutsFlag(temp_p));
     nTrack->Update(temp_n,fMCEvent,fPIDResponse, GetTrackCutsFlag(temp_n));
     
     
     Int_t cascade_pdg = 0;      // 0 would mean not matched to a unqie
     // cascade (or data)
     Bool_t cascade_primary = kFALSE;  // 1 = primary, 0 = secondary (or data)

     if(fMCEvent){
       const Int_t p_label = TMath::Abs(pTrack->GetLabel());
       const Int_t n_label = TMath::Abs(nTrack->GetLabel());
       // We check if the mother label matches in the next method.  If they do
       // not match then we checks if a pion might show up as a muon
       // (Xi->Lambda->pi->mu or Xi->pi->mu) and then one should go one step
       // further back.
       const Int_t v0_label = FindCommonMother(p_label, n_label);
       if(v0_label) {
      
	 // In principle we could check here that the V0 is the correct one but
	 // we think this is not necessary as we next check that the V0 and the
	 // bachelor comes from the same mother
      
	 const Int_t bach_label    = casc->GetBindex();
	 const Int_t cascade_label = FindCommonMother(v0_label, bach_label); 
	 if(cascade_label) {
	
	   TParticle* cascade_mcTrack = fMCEvent->Particle(cascade_label);	    
	   if(cascade_mcTrack) {
	  
	     if(fMCEvent->IsPhysicalPrimary(cascade_label))
	       cascade_primary = kTRUE;
	  
	     cascade_pdg = cascade_mcTrack->GetPdgCode();
	   }
	 }
       }
     }
     V0->Update(pTrack, nTrack, InvV0Masses, lV0Radius, lDcaV0DaughtersXi, lV0toXiCosineOfPointingAngle, V0Pt, V0Eta, lDcaV0ToPrimVertexXi);
     Cascade->Update(V0, bTrack, InvMasses,  lXiRadius, lDcaXiDaughters, lXiCosineOfPointingAngle, CascPt, CascEta, lDcaXiToPrimVertex, temp_b->GetSign(), lDcaV0ToPrimVertex, cascade_pdg, cascade_primary);

  };   

};




void
AliAnalysisTaskTPCTOFCascade::UserExec(Option_t *option)
{
  /*
   * user exec
   */

  /*** INITIALIZATION ***/

  /* init run */
  if (!InitRun()) return;
  FillHist(4);
  /* init event */
  if (!InitEvent()) return;
  FillHist(5);
  fAnalysisEvent->Reset();
  Int_t EventSelectionFlag = 0;
  Float_t V0MPercentile = -1000;
  AliMultSelection *ams = (AliMultSelection*)fESDEvent->FindListObject("MultSelection");
  if(!ams)
    V0MPercentile = -999;
  else {
    V0MPercentile = ams->GetMultiplicityPercentile("V0M");
    if(ams->GetThisEventIsNotPileup()) EventSelectionFlag += AliAnalysisPIDCascadeEvent::kNotPileupInSPD;
    if(ams->GetThisEventIsNotPileupMV()) EventSelectionFlag += AliAnalysisPIDCascadeEvent::kNotPileupInMV;
    if(ams->GetThisEventIsNotPileupInMultBins()) EventSelectionFlag += AliAnalysisPIDCascadeEvent::kNotPileupInMB;
    if(ams->GetThisEventINELgtZERO()) EventSelectionFlag+=AliAnalysisPIDCascadeEvent::kINELgtZERO;
    if(ams->GetThisEventHasNoInconsistentVertices()) EventSelectionFlag+=AliAnalysisPIDCascadeEvent::kNoInconsistentVtx;
    if(ams->GetThisEventIsNotAsymmetricInVZERO()) EventSelectionFlag+=AliAnalysisPIDCascadeEvent::kNoV0Asym;
  };
  Bool_t lSPDandTrkVtxExists=kFALSE;
  Bool_t lPassProximityCut=kTRUE;
  if(SelectVertex2015pp(fESDEvent,kTRUE,&lSPDandTrkVtxExists,&lPassProximityCut)) EventSelectionFlag+=AliAnalysisPIDCascadeEvent::kVertexSelected2015pp;
  if(lSPDandTrkVtxExists) EventSelectionFlag+=AliAnalysisPIDCascadeEvent::kSPDandTrkVtxExists;
  if(lPassProximityCut) EventSelectionFlag+=AliAnalysisPIDCascadeEvent::kPassProximityCut;
  fAnalysisEvent->SetV0Mmultiplicity(V0MPercentile);
  fAnalysisEvent->SetEventFlags(EventSelectionFlag);
  //  AliVVZERO *v0 = fESDEvent->GetVZEROData();
  // for(Int_t i=0;i<32;i++) fAnalysisEvent->SetV0CellAmplitude(i,v0->GetMultiplicityV0A(i));
  // for(Int_t i=32;i<64;i++) fAnalysisEvent->SetV0CellAmplitude(i,v0->GetMultiplicityV0C(i-32));
  
  Float_t RefMult08 = -1000;
  Float_t RefMult05 = -1000;
  Float_t SPDTracklets = -1000;
  if(!ams){ 
   RefMult08 = -999;
   RefMult05 = -999;
   SPDTracklets = -999;
  }
  if(ams){
    RefMult08 = ams->GetMultiplicityPercentile("RefMult08");
    RefMult05 = ams->GetMultiplicityPercentile("RefMult05");
    SPDTracklets = ams->GetMultiplicityPercentile("SPDTracklets");
  }

  fAnalysisEvent->SetRefMult08(RefMult08);
  fAnalysisEvent->SetRefMult05(RefMult05);
  fAnalysisEvent->SetSPDTracklets(SPDTracklets);
  /*** MC PRIMARY PARTICLES ***/

  Int_t mcmulti = 0;
  if (fMCFlag) {

    /* reset track array */
    fAnalysisParticleArray->Clear();
    
    /* loop over primary particles */
    Int_t nPrimaries = fMCEvent->GetNumberOfPrimaries();//fMCStack->GetNprimary();
    TParticle *particle;
    TParticlePDG *particlePDG;
    /* loop over primary particles */
    for (Int_t ipart = 0; ipart < nPrimaries; ipart++) {
      /* get particle */
      //particle = fMCEvent->Particle(ipart);//((AliMCParticle*)fMCEvent->GetTrack(ipart))->Particle();//fMCStack->Particle(ipart);
      particle = ((AliMCParticle*)fMCEvent->GetTrack(ipart))->Particle();
      if (!particle) continue;
      /* get particlePDG */
      particlePDG = particle->GetPDG();
      if (!particlePDG) continue;

      /* check primary */
      if (!fMCEvent->IsPhysicalPrimary(ipart)) 
	continue;

      /* check charged */
      //if ((particlePDG->Charge()==0.)&&(!OWSave)) continue;
      mcmulti++;
      /* check rapidity and pt cuts */
       if ( (TMath::Abs(particle->Energy() - particle->Pz()) < 1e-6)  ) continue;
       if ( ((particle->Energy() + particle->Pz())/(particle->Energy() - particle->Pz())) < 0) continue;
       Double_t PRap = particle->Y();
       if(std::isnan(PRap))
	 continue;
	
       if (TMath::Abs(particle->Y()) > fRapidityCut) continue;
       //if (particle->Pt() < 0.15) continue; //Maybe remove to properly correct for feeddown?
      //Get mother PDG code. In principle, can be optimized by only doing if for OWSace, as the rest of the particles are physical primaries
      Int_t indexMother = particle->GetFirstMother();
      Int_t lMotherPDG=0; //Just to be safe
      if(indexMother>=0) {
	TParticle *MotherParticle = fMCEvent->Particle(indexMother);
	lMotherPDG = MotherParticle->GetPdgCode();
      }; 

      /* update and add analysis particle */
      fAnalysisParticle->Update(particle, ipart, lMotherPDG);
      new ((*fAnalysisParticleArray)[fAnalysisParticleArray->GetEntries()]) AliAnalysisPIDCascadeParticle(*fAnalysisParticle);
    } /* end of loop over primary particles */

    
  }

  /*** GLOBAL EVENT INFORMATION ***/

  //  fAnalysisEvent->Reset(); // Moved up
  /* update global event info */
  fAnalysisEvent->SetIsCollisionCandidate(fIsCollisionCandidate);
  fAnalysisEvent->SetIsEventSelected(fIsEventSelected);
  fAnalysisEvent->SetIsPileupFromSPD(fIsPileupFromSPD);
  fAnalysisEvent->SetHasVertex(fHasVertex);
  fAnalysisEvent->SetVertexZ(fVertexZ);
  // fAnalysisEvent->SetMCTimeZero(fMCTimeZero);
  fAnalysisEvent->SetRunNumber(fRunNumber);
  fAnalysisEvent->SetMagneticField(fESDEvent->GetMagneticField());

  Int_t refmulti;
  refmulti = AliESDtrackCuts::GetReferenceMultiplicity(fESDEvent, AliESDtrackCuts::kTrackletsITSTPC,0.8);  
  fAnalysisEvent->SetReferenceMultiplicity(refmulti);
  fAnalysisEvent->SetMCMultiplicity(mcmulti);
  /*** RECONSTRUCTED TRACKS ***/

  /* reset track array */
  fAnalysisTrackArray->Clear();
  // fAnalysisV0TrackArray->Clear();
  /* loop over ESD tracks */
  Int_t nTracks = fESDEvent->GetNumberOfTracks();
  AliESDtrack *track;
  for (Int_t itrk = 0; itrk < nTracks; itrk++) {
    /* get track */
    track = fESDEvent->GetTrack(itrk);
    if (!track) continue;
    /* check accept track */
    Int_t trflag = GetTrackCutsFlag(track);
    if(!trflag) continue;
    
    /* update and add analysis track */
    fAnalysisTrack->Update(track, fMCEvent,fPIDResponse, trflag);
    const AliESDVertex *vtx = fESDEvent->GetPrimaryVertexTracks();
    if(!vtx || !vtx->GetStatus())
      vtx = fESDEvent->GetPrimaryVertexSPD();
    // if(vtx) {
    //   if(vtx->GetStatus()) {
    // 	Double_t ChiConstrained = track->GetChi2TPCConstrainedVsGlobal(vtx);
    // 	fAnalysisTrack->SetChi2TPCConstrainedVsGlobal(ChiConstrained);
    //   } else
    // 	fAnalysisTrack->SetChi2TPCConstrainedVsGlobal(-8);
    // };
    // if(track->IsEMCAL()) {
    //   AliVCluster *lvcl = fESDEvent->GetCaloCluster(track->GetEMCALcluster());
    //   if(lvcl)
    // 	fAnalysisTrack->SetEMCalPars(lvcl->E(),track->GetTrackPOnEMCal());
    // };
    new ((*fAnalysisTrackArray)[fAnalysisTrackArray->GetEntries()]) AliAnalysisPIDCascadeTrack(*fAnalysisTrack);

    

  } /* end of loop over ESD tracks */
  ProcessV0s();
  ProcessCascades();
  fPIDTree->Fill();

  PostData(1,fPIDTree);
  PostData(2,fEvHist);

  fAnalysisV0TrackArray->Clear();
  fAnalysisCascadeTrackArray->Clear();
}

void AliAnalysisTaskTPCTOFCascade::Terminate(Option_t *) {
  printf("Terminate!\n");
}

Int_t AliAnalysisTaskTPCTOFCascade::FindCommonMother(Int_t label_1, Int_t label_2)
{

  if (!fMCEvent) 
    return 0;
  
  TParticle* mcTrack_1 = fMCEvent->Particle(label_1);	    
  TParticle* mcTrack_2 = fMCEvent->Particle(label_2);	    
  if (!mcTrack_1 || !mcTrack_2)
    return 0;
  
  Int_t mother_label_1 = TMath::Abs(mcTrack_1->GetMother(0));
  Int_t mother_label_2 = TMath::Abs(mcTrack_2->GetMother(0));
  
  if(mother_label_1 == mother_label_2)
    return mother_label_1;
  
  // label_1 or label_2 was primary tracks
  if(mother_label_1 < 0 || mother_label_2 < 0)
    return 0;
  
  // we just check that one track could be a muon from pi->muon
  TParticle* mcTrack_mother_1 = fMCEvent->Particle(mother_label_1);	    
  TParticle* mcTrack_mother_2 = fMCEvent->Particle(mother_label_2);	    
  if(!mcTrack_mother_1 || !mcTrack_mother_2)
    return 0;
  
  if(TMath::Abs(mcTrack_mother_1->GetPdgCode()) == 13)
    mother_label_1 = TMath::Abs(mcTrack_mother_1->GetMother(0));
  if(TMath::Abs(mcTrack_mother_2->GetPdgCode()) == 13)
    mother_label_2 = TMath::Abs(mcTrack_mother_2->GetMother(0));
  
  if(mother_label_1 == mother_label_2)
    return mother_label_1;

  return 0;
}
    
