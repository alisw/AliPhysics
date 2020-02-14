#include "AliAnalysisTaskGFWFlow.h"
#include "AliEventCuts.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliVVertex.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"
#include "TList.h"
#include "TProfile.h"
#include "TH3D.h"
#include "AliEventCuts.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "TParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliGFWWeights.h"
#include "AliGFWCuts.h"
#include "AliGFWFlowContainer.h"
#include "TObjArray.h"
#include "TNamed.h"
#include "AliGFW.h"
#include "TRandom.h"
#include <vector>

ClassImp(AliAnalysisTaskGFWFlow);

AliAnalysisTaskGFWFlow::AliAnalysisTaskGFWFlow():
  AliAnalysisTaskSE(),
  debugpar(0),
  fTriggerType(AliVEvent::kINT7),
  fProduceWeights(kTRUE),
  fSelections(0),
  fWeightList(0),
  fWeights(0),
  fExtraWeights(0),
  fFC(0),
  fOutputTree(0),
  fMCEvent(0),
  fIsMC(kFALSE),
  fPtAxis(0),
  fPOIpTMin(0.2),
  fPOIpTMax(20),
  fRFpTMin(0.2),
  fRFpTMax(3.0),
  fWeightPath(""),
  fWeightDir(""),
  fTotFlags(15), //Total number of flags: 1 (nominal) + fTotTrackFlags + N_Event flags
  fTotTrackFlags(8), //Total number of track flags (without nominal)
  fRunNo(-1),
  fCurrSystFlag(0),
  fAddQA(kFALSE),
  fQAList(0),
  fBypassCalculations(kFALSE)
{
};
AliAnalysisTaskGFWFlow::AliAnalysisTaskGFWFlow(const char *name, Bool_t ProduceWeights, Bool_t IsMC, Bool_t AddQA):
  AliAnalysisTaskSE(name),
  debugpar(0),
  fTriggerType(AliVEvent::kINT7),
  fProduceWeights(ProduceWeights),
  fSelections(0),
  fWeightList(0),
  fWeights(0),
  fExtraWeights(0),
  fFC(0),
  fOutputTree(0),
  fIsMC(IsMC),
  fPtAxis(new TAxis()),
  fPOIpTMin(0.2),
  fPOIpTMax(20),
  fRFpTMin(0.2),
  fRFpTMax(3.0),
  fWeightPath(""),
  fWeightDir(""),
  fTotFlags(15),
  fTotTrackFlags(8),
  fRunNo(-1),
  fCurrSystFlag(0),
  fAddQA(AddQA),
  fQAList(0),
  fBypassCalculations(kFALSE)
{
  if(!fProduceWeights) DefineInput(1,TList::Class());
  DefineOutput(1,(fProduceWeights?TList::Class():AliGFWFlowContainer::Class()));
  if(fAddQA)
    DefineOutput(2,TList::Class());
};
AliAnalysisTaskGFWFlow::~AliAnalysisTaskGFWFlow() {
};
void AliAnalysisTaskGFWFlow::UserCreateOutputObjects(){
  printf("**************************************\n");
  printf("************  AliGFW v2  *************\n");
  printf("**************************************\n");
  OpenFile(1);
  fTotTrackFlags = AliGFWCuts::fNTrackFlags;
  fTotFlags = fTotTrackFlags+AliGFWCuts::fNEventFlags+1;
  fSelections = new AliGFWCuts*[fTotFlags]; //0 for normal, the rest for systematics
  for(Int_t i=0; i<fTotFlags;i++) {
    fSelections[i] = new AliGFWCuts();
    fSelections[i]->SetupCuts(i);
  };
  if(fProduceWeights) {
    //Initialize selection objects
    fWeightList = new TList();
    fWeightList->SetName("WeightList");
    fWeightList->SetOwner(kTRUE);
    for(Int_t i=0;i<fTotFlags;i++) {
      if(!fSelections[i]->NeedsExtraWeight()) continue;
      fWeightList->Add(new AliGFWWeights());
      fWeights = (AliGFWWeights*)fWeightList->Last();
      fWeights->SetName(Form("weights%s",fSelections[i]->GetSystPF()));
      fWeights->Init(!fIsMC,fIsMC); // AddData = !fIsMC; AddMC = fIsMC
    }
  } else {
    //Setup the structure of FC:
    TObjArray *OAforPt=new TObjArray();
    //No gap:
    OAforPt->Add(new TNamed("MidV22","MidV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV22_pt_%i",i+1),"MidV22_pTDiff"));
    OAforPt->Add(new TNamed("MidV24","MidV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV24_pt_%i",i+1),"MidV24_pTDiff"));
    OAforPt->Add(new TNamed("MidV26","MidV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV26_pt_%i",i+1),"MidV26_pTDiff"));
    OAforPt->Add(new TNamed("MidV28","MidV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV28_pt_%i",i+1),"MidV28_pTDiff"));
    OAforPt->Add(new TNamed("MidV32","MidV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV32_pt_%i",i+1),"MidV32_pTDiff"));
    OAforPt->Add(new TNamed("MidV34","MidV34"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV34_pt_%i",i+1),"MidV34_pTDiff"));
    OAforPt->Add(new TNamed("MidV42","MidV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV42_pt_%i",i+1),"MidV42_pTDiff"));
    OAforPt->Add(new TNamed("MidV52","MidV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV52_pt_%i",i+1),"MidV52_pTDiff"));

    //2SENeg:
    OAforPt->Add(new TNamed("Mid2SENV22","Mid2SENV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV22_pt_%i",i+1),"Mid2SENV22_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV24","Mid2SENV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV24_pt_%i",i+1),"Mid2SENV24_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV26","Mid2SENV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV26_pt_%i",i+1),"Mid2SENV26_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV28","Mid2SENV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV28_pt_%i",i+1),"Mid2SENV28_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV32","Mid2SENV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV32_pt_%i",i+1),"Mid2SENV32_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV42","Mid2SENV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV42_pt_%i",i+1),"Mid2SENV42_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SENV52","Mid2SENV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SENV52_pt_%i",i+1),"Mid2SENV52_pTDiff"));

    //2SEPos:
    OAforPt->Add(new TNamed("Mid2SEPV22","Mid2SEPV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV22_pt_%i",i+1),"Mid2SEPV22_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV24","Mid2SEPV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV24_pt_%i",i+1),"Mid2SEPV24_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV26","Mid2SEPV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV26_pt_%i",i+1),"Mid2SEPV26_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV28","Mid2SEPV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV28_pt_%i",i+1),"Mid2SEPV28_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV32","Mid2SEPV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV32_pt_%i",i+1),"Mid2SEPV32_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV42","Mid2SEPV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV42_pt_%i",i+1),"Mid2SEPV42_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEPV52","Mid2SEPV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEPV52_pt_%i",i+1),"Mid2SEPV52_pTDiff"));

    //|eta|>0.5
    OAforPt->Add(new TNamed("MidGapNV22","MidGapNV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV22_pt_%i",i+1),"MidGapNV22_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV24","MidGapNV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV24_pt_%i",i+1),"MidGapNV24_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV26","MidGapNV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV26_pt_%i",i+1),"MidGapNV26_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV28","MidGapNV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV28_pt_%i",i+1),"MidGapNV28_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV32","MidGapNV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV32_pt_%i",i+1),"MidGapNV32_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV42","MidGapNV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV42_pt_%i",i+1),"MidGapNV42_pTDiff"));
    OAforPt->Add(new TNamed("MidGapNV52","MidGapNV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapNV52_pt_%i",i+1),"MidGapNV52_pTDiff"));

    //|eta|>0.5
    OAforPt->Add(new TNamed("MidGapPV22","MidGapPV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV22_pt_%i",i+1),"MidGapPV22_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV24","MidGapPV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV24_pt_%i",i+1),"MidGapPV24_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV26","MidGapPV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV26_pt_%i",i+1),"MidGapPV26_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV28","MidGapPV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV28_pt_%i",i+1),"MidGapPV28_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV32","MidGapPV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV32_pt_%i",i+1),"MidGapPV32_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV42","MidGapPV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV42_pt_%i",i+1),"MidGapPV42_pTDiff"));
    OAforPt->Add(new TNamed("MidGapPV52","MidGapPV52"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapPV52_pt_%i",i+1),"MidGapPV52_pTDiff"));


    //Multi bins:
    Double_t multibins[] = {5,10,20,30,40,50,60,70};
    fFC = new AliGFWFlowContainer();
    fFC->SetName(Form("FC%s",fSelections[fCurrSystFlag]->GetSystPF()));
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(OAforPt,7,multibins,10); //Statistics only required for nominal profiles, so do not create randomized profiles for systematics
    //Powers per harmonic:
    Int_t NoGap[] = {9,0,8,4,7,2,6,0,5};
    Int_t WithGap[] = {5,0,2,2,3,2,4,0,5};
    fGFW = new AliGFW();
    //Full regions
    fGFW->AddRegion("poiMid",9,NoGap,-0.8,0.8,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refMid",9,NoGap,-0.8,0.8,1,2);
    //2 subevets:
    fGFW->AddRegion("poiSENeg",9,WithGap,-0.8,0.,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refSENeg",9,WithGap,-0.8,0.,1,2);
    fGFW->AddRegion("poiSEPos",9,WithGap,0.,0.8,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refSEPos",9,WithGap,0.,0.8,1,2);
    //With gap
    fGFW->AddRegion("poiGapNeg",9,WithGap,-0.8,-0.5,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refGapNeg",9,WithGap,-0.8,-0.5,1,2);
    fGFW->AddRegion("poiGapPos",9,WithGap,0.5,0.8,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refGapPos",9,WithGap,0.5,0.8,1,2);
  };
  if(fProduceWeights) PostData(1,fWeightList);
  else PostData(1,fFC);
  if(fAddQA) {
    fQAList = new TList();
    fQAList->SetOwner(kTRUE);
    fEventCuts.AddQAplotsToList(fQAList);
    PostData(2,fQAList);
  };
  if(!fProduceWeights) {
    fWeightList = (TList*) GetInputData(1);
    if(!fWeightList) { AliFatal("Could not retrieve weight list!\n"); return; };
    CreateCorrConfigs();
  };
  // printf("\n******************\nStarting the watch\n*****************\n");
  // mywatchFill.Reset();
  // mywatchStore.Reset();
  // mywatch.Start(kTRUE);
};
void AliAnalysisTaskGFWFlow::UserExec(Option_t*) {
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t cent = lMultSel->GetMultiplicityPercentile("V0M");
  if(!CheckTriggerVsCentrality(cent)) return;
  if(!fProduceWeights)
    if(!InitRun()) return;
  if(!AcceptEvent()) return;
  if(!AcceptAODVertex(fAOD)) return;
  if(fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if(!fMCEvent)
      return;
  };
  if(fCurrSystFlag==fTotTrackFlags+4) cent = lMultSel->GetMultiplicityPercentile("CL1"); //CL1 flag is EvFlag 4 = N_TrackFlags + 4
  if(fCurrSystFlag==fTotTrackFlags+5) cent = lMultSel->GetMultiplicityPercentile("CL0"); //CL0 flag is EvFlag 5 = N_TrackFlags + 5
  if(cent<5) return; //Do not consider 0-5%
  if(cent>70) return; //Also, peripheral cutoff
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  Int_t vtxb = GetVtxBit(fAOD);
  if(!vtxb) return; //If no vertex pass, then do not consider further
  if(fBypassCalculations) return;
  //fFlowEvent->SetRunNumber(fAOD->GetRunNumber());
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(fAOD->GetPrimaryVertex());
  Double_t POSvtx[] = {0.,0.,0.};
  vtx->GetXYZ(POSvtx);
  TClonesArray *tca = 0;
  if(fProduceWeights) {
    if(fIsMC) {
      tca = (TClonesArray*)fInputEvent->FindListObject("mcparticles");
      for(Int_t i=0;i<tca->GetEntries();i++) {
      	AliVParticle *lPart;
      	lPart = (AliAODMCParticle*)tca->At(i);
      	if(!AcceptParticle(lPart)) continue;
      	Int_t partbit = GetParticleBit(lPart);
      	if(!partbit) continue; //if no particle bit is set, no need to continue
      	Int_t combinedbit = CombineBits(vtxb,partbit);
      	Int_t count=0;
      	for(Int_t i=0;i<fTotFlags;++i)
      	  if(fSelections[i]->NeedsExtraWeight()) {
      	    if(combinedbit&(1<<i))
      	      ((AliGFWWeights*)fWeightList->At(count))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),cent,2);
      	    count++;
      	  };
      };
    };
    AliAODTrack *lTrack;
    for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
      lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
      Double_t POStrk[] = {0.,0.,0.};
      lTrack->GetXYZ(POStrk);
      Double_t DCA[] = {0.,0.,0.};
      for(Int_t i=0;i<3;i++) DCA[i] = POSvtx[i]-POStrk[i];
      Double_t dcaxy = TMath::Sqrt(DCA[0]*DCA[0]+DCA[1]*DCA[1]);
      Double_t lDCA[] = {TMath::Abs(DCA[2]),dcaxy};
      Int_t trackbit = GetTrackBit(lTrack,lDCA);
      if(!trackbit) continue; //if no track bit is set, no need to continue
      Int_t combinedbit = CombineBits(vtxb,trackbit);
      if(fIsMC) {
        AliVParticle *lPart = (AliAODMCParticle*)tca->At(TMath::Abs(lTrack->GetLabel()));
        Int_t count=0;
        for(Int_t i=0;i<fTotFlags;++i)
          if(fSelections[i]->NeedsExtraWeight()) {
            if(combinedbit&(1<<i))
              ((AliGFWWeights*)fWeightList->At(count))->Fill(lPart->Phi(),lPart->Eta(),vz,lPart->Pt(),cent,1);
            count++;
          };
        } else {
          Int_t count=0;
          for(Int_t i=0;i<fTotFlags;++i)
            if(fSelections[i]->NeedsExtraWeight()) {
              if((combinedbit&(1<<i)))
                ((AliGFWWeights*)fWeightList->At(count))->Fill(lTrack->Phi(),lTrack->Eta(),vz,lTrack->Pt(),cent,0);
              count++;
            };
        };
    };
    PostData(1,fWeightList);
    if(fAddQA) PostData(2,fQAList);
    return;
  } else {
    fGFW->Clear();
    AliAODTrack *lTrack;
    if(!fSelections[fCurrSystFlag]->AcceptVertex(fAOD,1)) return;
    // mywatchFill.Start(kFALSE);
    for(Int_t lTr=0;lTr<fAOD->GetNumberOfTracks();lTr++) {
      lTrack = (AliAODTrack*)fAOD->GetTrack(lTr);
      //if(!AcceptAODTrack(lTrack,tca)) continue;
      Double_t POStrk[] = {0.,0.,0.};
      lTrack->GetXYZ(POStrk);
      Double_t DCA[] = {0.,0.,0.};
      for(Int_t i=0;i<3;i++) DCA[i] = POSvtx[i]-POStrk[i];
      Double_t dcaxy = TMath::Sqrt(DCA[0]*DCA[0]+DCA[1]*DCA[1]);
      Double_t lDCA[] = {TMath::Abs(DCA[2]),dcaxy};
      if(!fSelections[fCurrSystFlag]->AcceptTrack(lTrack,lDCA) &&
      	 !fSelections[9]->AcceptTrack(lTrack,lDCA)) continue;

      Double_t l_pT=lTrack->Pt();
      Bool_t WithinPtPOI = (fPOIpTMin<l_pT) && (l_pT<fPOIpTMax); //within POI pT range
      Bool_t WithinPtRF  = (fRFpTMin <l_pT) && (l_pT<fRFpTMax);  //within RF pT range
      if(!WithinPtPOI && !WithinPtRF) continue; //if the track is not within any pT range, then continue
      if(!fWeights) printf("Weights do not exist!\n");
      Double_t nua = fWeights->GetNUA(lTrack->Phi(),lTrack->Eta(),vz);
      Double_t nue = 1; //Since doing pT-diff., we can set this to one for speed up.
      //To speed up, call getter for NUA directly
      //Double_t nua = fWeights->GetWeight(lTrack->Phi(),lTrack->Eta(),vz,l_pT,cent,0);
      //Double_t nuaITS = fExtraWeights->GetWeight(lTrack->Phi(),lTrack->Eta(),vz,lTrack->Pt(),cent,0);
      //Double_t nue = fPtAxis->GetNbins()>1?1:fWeights->GetWeight(lTrack->Phi(),lTrack->Eta(),vz,cent,l_pT,1);
      if(fSelections[fCurrSystFlag]->AcceptTrack(lTrack, lDCA)) {
      	if(WithinPtPOI) fGFW->Fill(lTrack->Eta(),fPtAxis->FindBin(l_pT)-1,lTrack->Phi(),nua*nue,1); //Fill POI (mask = 1)
        if(WithinPtRF)  fGFW->Fill(lTrack->Eta(),fPtAxis->FindBin(l_pT)-1,lTrack->Phi(),nua*nue,2); //Fit RF (mask = 2)
      }
      /*if(fSelections[9]->AcceptTrack(lTrack, lDCA)) //No ITS for now
	fGFW->Fill(lTrack->Eta(),fPtAxis->FindBin(lTrack->Pt())-1,lTrack->Phi(),nuaITS*nue,2);*/
    };
    // mywatchFill.Stop();
    TRandom rndm(0);
    Double_t rndmn=rndm.Rndm();
    //Calculate & fill profiles:
    //V_2{n}, full acceptance
    // mywatchStore.Start(kFALSE);
    Bool_t filled;
    for(Int_t l_ind=0; l_ind<corrconfigs.size(); l_ind++) {
      //Bool_t DisableOL=kFALSE;
      //if(l_ind<14) DisableOL = (l_ind%2); //Only for 1, 3, 5 ... 13
      filled = FillFCs(corrconfigs.at(l_ind),cent,rndmn);//,DisableOL);
    };
    // mywatchStore.Stop();
    PostData(1,fFC);
    if(fAddQA) PostData(2,fQAList);
  };
};
void AliAnalysisTaskGFWFlow::Terminate(Option_t*) {
  // printf("\n********* Time: %f\n**********",mywatch.RealTime());
  // printf("Filling time: %f\n",mywatchFill.RealTime());
  // printf("Storing time: %f\n",mywatchStore.RealTime());
};
void AliAnalysisTaskGFWFlow::SetPtBins(Int_t nBins, Double_t *bins, Double_t RFpTMin, Double_t RFpTMax) {
  fPtAxis->Set(nBins, bins);
  //Set pT range given by the defined pT axis
  fPOIpTMin=bins[0];
  fPOIpTMax=bins[nBins];
  //if RF pT range is not defined, use the same as for POI
  fRFpTMin=(RFpTMin<0)?fPOIpTMin:RFpTMin;
  fRFpTMax=(RFpTMax<0)?fPOIpTMax:RFpTMax;
}

Bool_t AliAnalysisTaskGFWFlow::AcceptEvent() {
  if(!fEventCuts.AcceptEvent(fInputEvent)) return 0;
  if(fCurrSystFlag==15) if(!fEventCutsForPU.AcceptEvent(fInputEvent)) return 0; //For a tight PU cut, an additional selection
  return kTRUE;
};
Bool_t AliAnalysisTaskGFWFlow::CheckTriggerVsCentrality(Double_t l_cent) {
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  //if(!(fSelMask&fTriggerType)) return kFALSE;
  Bool_t centTrigger = (fSelMask&(AliVEvent::kCentral)) && l_cent>0 && l_cent<10; //fTriggerType& removed
  Bool_t semiCentTri = (fSelMask&(AliVEvent::kSemiCentral)) && l_cent>30 && l_cent<50; //fTriggerType& removed
  Bool_t MBTrigger   = (fSelMask&(AliVEvent::kMB+AliVEvent::kINT7)); //fTriggerType& removed
  if(centTrigger || semiCentTri || MBTrigger) return kTRUE;
  return kFALSE;

}
Bool_t AliAnalysisTaskGFWFlow::AcceptAODVertex(AliAODEvent *inEv) {
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;

  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertexSPD());
  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;

  const Double_t aodVtxZ = vtx->GetZ();
  if(TMath::Abs(aodVtxZ) > 10)
    return kFALSE;

  return kTRUE;
};

Bool_t AliAnalysisTaskGFWFlow::AcceptParticle(AliVParticle *mpa) {
  if(!mpa->IsPhysicalPrimary()) return kFALSE;
  if(mpa->Charge()==0) return kFALSE;
  if(TMath::Abs(mpa->Eta())>0.8) return kFALSE;
  if(mpa->Pt()<0.3) return kFALSE;
  if(mpa->Pt()>20) return kFALSE;
  return kTRUE;
};
Int_t AliAnalysisTaskGFWFlow::GetVtxBit(AliAODEvent *mev) {
  Int_t retbit=0;
  for(Int_t vtxst=fTotTrackFlags+1;vtxst<fTotFlags; vtxst++) //vtx flags are 9-14
    retbit+=fSelections[vtxst]->AcceptVertex(mev,vtxst);

  //PU unc. bit has been set with a wide cut. Override with a check w/ smaller cut:
  Int_t PUSystFlag = 1<<(fTotTrackFlags+6);//PU syst. flag is 6-th in ev/vtx selection
  //Only override if the event has passed the initial selection:
  if((retbit&PUSystFlag)==PUSystFlag) //if event accepted by nominal cuts
    if(!fEventCutsForPU.AcceptEvent(mev)) //Check if it gets rejected by the PU cut
    retbit-=PUSystFlag; //and if so, set the PU flag to zero
  retbit+=fSelections[0]->AcceptVertex(mev); //the standard one
  return retbit;
};
Int_t AliAnalysisTaskGFWFlow::GetParticleBit(AliVParticle *mpa) {
  Int_t retbit=0;
  for(Int_t i=1;i<=fTotTrackFlags;i++)
    retbit+=fSelections[i]->AcceptParticle(mpa,i);
  retbit+=fSelections[0]->AcceptParticle(mpa);
  return retbit;
};
Int_t AliAnalysisTaskGFWFlow::GetTrackBit(AliAODTrack* mtr, Double_t *lDCA) {
  Int_t retbit=0;
  for(Int_t i=1;i<=fTotTrackFlags; i++) //track flags are 1-8
    retbit+=fSelections[i]->AcceptTrack(mtr,lDCA,i);
  //Nominal (TPC) tracks:
  retbit+=fSelections[0]->AcceptTrack(mtr,lDCA);
  //Also, add the ITS tracks:
  //retbit+=fSelections[12]->AcceptTrack(mtr,lDCA,12); //
  return retbit;
};
Int_t AliAnalysisTaskGFWFlow::CombineBits(Int_t VtxBit, Int_t TrkBit) {
  Int_t retbit=((VtxBit&1)*TrkBit | (TrkBit&1)*VtxBit);
  //Also add ITS track (bit 12), if the nominal vtx bit is set:
  //retbit=retbit|((VtxBit&1)*(TrkBit&(1<<13)));
  return retbit;
};
Bool_t AliAnalysisTaskGFWFlow::SetInputWeightList(TList *inlist) {
  if(!inlist) {
    return kFALSE;
  };
  fWeightList = inlist;
  return kTRUE;
};
Bool_t AliAnalysisTaskGFWFlow::LoadWeights(Int_t runno) { //Cannot be used when running on the trains
  if(fWeightList) {
    fWeights = (AliGFWWeights*)fWeightList->FindObject(Form("w%i%s",runno,fSelections[fCurrSystFlag]->NeedsExtraWeight()?
							    fSelections[fCurrSystFlag]->GetSystPF():""));
    if(!fWeights) {
      AliFatal("Weights could not be found in the list!\n");
      return kFALSE;
    };
    fWeights->CreateNUA();
    fWeights->CreateNUE();
    return kTRUE;
  } else {
    AliFatal("Weight list (for some reason) not set!\n");
    return kFALSE;
  };
  printf("You should not be here!\n");
  return kFALSE;
//If weights not set, attempting to fetch them from pre-set directory. This will definitely fail if running on train
  fWeightPath.Clear();
  fWeightPath.Append(fWeightDir.Data());
  fWeightPath.Append(Form("%i.root",runno));
  TFile *tfWeights = TFile::Open(fWeightPath.Data());
  if(!tfWeights) {
    printf("Could not open %s!\n",fWeightPath.Data());
    return kFALSE;
  };
  TList *l_List = (TList*)tfWeights->Get("OutputList");
  if(!l_List) printf("\n\n\n\n\n\n**************************\n\n\n\n\n\n Could not fetch OutputList!\n\n\n\n");
  TString l_weightname(fSelections[fSelections[fCurrSystFlag]->NeedsExtraWeight()?fCurrSystFlag:0]->GetSystPF());
  fWeights = (AliGFWWeights*)l_List->FindObject(Form("weights%s",l_weightname.Data()));
  if(!fWeights) {
    printf("Could not fetch weights%s from %s!\n",fSelections[fCurrSystFlag]->GetSystPF(),fWeightPath.Data());
    return kFALSE;
  };
  fExtraWeights =  (AliGFWWeights*)tfWeights->Get(Form("weights%s",fSelections[9]->GetSystPF()));//Used for when several weights are need (e.g. ITS): (Weights*)tfWeights->Get(Form("weights%s",fSelections[12]->GetSystPF()));
  tfWeights->Close();
  if(!fWeights) return kFALSE;
  fWeights->CreateNUA();
  fWeights->CreateNUE();
  if(fExtraWeights) {
    fExtraWeights->CreateNUA();
    fExtraWeights->CreateNUE();
  };

  return kTRUE;
};
Bool_t AliAnalysisTaskGFWFlow::InitRun() {
  if(!fInputEvent) return kFALSE;
  Int_t runno = fInputEvent->GetRunNumber();
  if(fRunNo!=runno) {
    if(!LoadWeights(runno))
      return kFALSE;
    else
      fRunNo = runno;
    //Run has changed; need to re-override the PU cut in AliEventCuts:
    Bool_t dump1 = fEventCuts.AcceptEvent(fInputEvent);//To setup all the event cuts (automatically)
    dump1 = fEventCutsForPU.AcceptEvent(fInputEvent);//To setup all the event cuts (automatically)
    fEventCuts.fUseVariablesCorrelationCuts = kTRUE; //By default this is not used (for some reason?)
    //Also for the PU systematics:
    fEventCutsForPU.fUseVariablesCorrelationCuts = kTRUE;
    fEventCutsForPU.fESDvsTPConlyLinearCut[0] = 1500; //Cut for systematic (nominal is 15 000)
  };
  return kTRUE;
};
Bool_t AliAnalysisTaskGFWFlow::FillFCs(TString head, TString hn, Double_t cent, Bool_t diff, Double_t rndmn) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(hn,kTRUE).Re();
  if(dnx==0) return kFALSE;
  if(!diff) {
    val = fGFW->Calculate(hn).Re()/dnx;
    if(TMath::Abs(val)<1)
      fFC->FillProfile(head.Data(),cent,val,dnx,rndmn);
    return kTRUE;
  };
  for(Int_t i=1;i<=fPtAxis->GetNbins();i++) {
    TString tss(hn);
    tss.Prepend(Form("(%i) ",i-1));
    dnx = fGFW->Calculate(tss,kTRUE).Re();
    if(dnx==0) continue;
    val = fGFW->Calculate(tss).Re()/dnx;
    if(TMath::Abs(val)<1)
      fFC->FillProfile(Form("%s_pt_%i",head.Data(),i),cent,val,dnx,rndmn);
  };
  return kTRUE;
};
Bool_t AliAnalysisTaskGFWFlow::FillFCs(AliGFW::CorrConfig corconf, Double_t cent, Double_t rndmn, Bool_t DisableOverlap) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).Re();
  if(dnx==0) return kFALSE;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).Re()/dnx;
    if(TMath::Abs(val)<1)
      fFC->FillProfile(corconf.Head.Data(),cent,val,dnx,rndmn);
    return kTRUE;
  };
  /*Int_t binDisableOLFrom = fPtAxis->GetNbins()+1;
  if(DisableOverlap)
    binDisableOLFrom = fPtAxis->FindBin(fRFpTMax); //To stay in the right bin*/
  Bool_t NeedToDisable=kFALSE;
  for(Int_t i=1;i<=fPtAxis->GetNbins();i++) {
    //if(DisableOverlap) NeedToDisable=(i>=binDisableOLFrom);
    dnx = fGFW->Calculate(corconf,i-1,kTRUE,NeedToDisable).Re();
    if(dnx==0) continue;
    val = fGFW->Calculate(corconf,i-1,kFALSE,NeedToDisable).Re()/dnx;
    if(TMath::Abs(val)<1)
      fFC->FillProfile(Form("%s_pt_%i",corconf.Head.Data(),i),cent,val,dnx,rndmn);
  };
  return kTRUE;
};
void AliAnalysisTaskGFWFlow::CreateCorrConfigs() {
//  corrconfigs = new AliGFW::CorrConfig[90];
  corrconfigs.push_back(GetConf("MidV22","refMid {2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV22","poiMid refMid {2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidV24","refMid {2 2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV24","poiMid refMid {2 2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidV26","refMid {2 2 2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV26","poiMid refMid {2 2 2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidV28","refMid {2 2 2 2 -2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidV28","poiMid refMid {2 2 2 2 -2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidV32","refMid {3 -3}", kFALSE));
  corrconfigs.push_back(GetConf("MidV32","poiMid refMid {3 -3}", kTRUE));
  corrconfigs.push_back(GetConf("MidV34","refMid {3 3 -3 -3}", kFALSE));
  corrconfigs.push_back(GetConf("MidV34","poiMid refMid {3 3 -3 -3}", kTRUE));
  corrconfigs.push_back(GetConf("MidV42","refMid {4 -4}", kFALSE));
  corrconfigs.push_back(GetConf("MidV42","poiMid refMid {4 -4}", kTRUE));
  corrconfigs.push_back(GetConf("MidV52","poiMid refMid {5 -5}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV22","refSENeg {2} refSEPos {-2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV24","refSENeg {2 2} refSEPos {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV26","refSENeg {2 2 2} refSEPos {-2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV28","refSENeg {2 2 2 2} refSEPos {-2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV22","poiSENeg refSENeg {2} refSEPos {-2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV24","poiSENeg refSENeg {2 2} refSEPos {-2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV26","poiSENeg refSENeg {2 2 2} refSEPos {-2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV28","poiSENeg refSENeg {2 2 2 2} refSEPos {-2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV22","refSEPos {2} refSENeg {-2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV24","refSEPos {2 2} refSENeg {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV26","refSEPos {2 2 2} refSENeg {-2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV28","refSEPos {2 2 2 2} refSENeg {-2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV22","poiSEPos refSEPos {2} refSENeg {-2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV24","poiSEPos refSEPos {2 2} refSENeg {-2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV26","poiSEPos refSEPos {2 2 2} refSENeg {-2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV28","poiSEPos refSEPos {2 2 2 2} refSENeg {-2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV32","refSENeg {3} refSEPos {-3}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV32","poiSENeg refSENeg {3} refSEPos {-3}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV32","refSEPos {3} refSENeg {-3}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV32","poiSEPos refSEPos {3} refSENeg {-3}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV42","refSENeg {4} refSEPos {-4}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV42","poiSENeg refSENeg {4} refSEPos {-4}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV42","refSEPos {4} refSENeg {-4}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV42","poiSEPos refSEPos {4} refSENeg {-4}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SENV52","refSENeg {5} refSEPos {-5}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SENV52","poiSENeg refSENeg {5} refSEPos {-5}", kTRUE));
  corrconfigs.push_back(GetConf("Mid2SEPV52","refSEPos {5} refSENeg {-5}", kFALSE));
  corrconfigs.push_back(GetConf("Mid2SEPV52","poiSEPos refSEPos {5} refSENeg {-5}", kTRUE));

  corrconfigs.push_back(GetConf("MidGapNV22","refGapNeg {2} refGapPos {-2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV24","refGapNeg {2 2} refGapPos {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV26","refGapNeg {2 2 2} refGapPos {-2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV28","refGapNeg {2 2 2 2} refGapPos {-2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV22","poiGapNeg refGapNeg {2} refGapPos {-2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV24","poiGapNeg refGapNeg {2 2} refGapPos {-2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV26","poiGapNeg refGapNeg {2 2 2} refGapPos {-2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV28","poiGapNeg refGapNeg {2 2 2 2} refGapPos {-2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV22","refGapPos {2} refGapNeg {-2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV24","refGapPos {2 2} refGapNeg {-2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV26","refGapPos {2 2 2} refGapNeg {-2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV28","refGapPos {2 2 2 2} refGapNeg {-2 -2 -2 -2}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV22","poiGapPos refGapPos {2} refGapNeg {-2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV24","poiGapPos refGapPos {2 2} refGapNeg {-2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV26","poiGapPos refGapPos {2 2 2} refGapNeg {-2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV28","poiGapPos refGapPos {2 2 2 2} refGapNeg {-2 -2 -2 -2}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV32","refGapNeg {3} refGapPos {-3}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV32","poiGapNeg refGapNeg {3} refGapPos {-3}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV32","refGapPos {3} refGapNeg {-3}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV32","poiGapPos refGapPos {3} refGapNeg {-3}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV42","refGapNeg {4} refGapPos {-4}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV42","poiGapNeg refGapNeg {4} refGapPos {-4}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV42","refGapPos {4} refGapNeg {-4}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV42","poiGapPos refGapPos {4} refGapNeg {-4}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapNV52","refGapNeg {5} refGapPos {-5}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapNV52","poiGapNeg refGapNeg {5} refGapPos {-5}", kTRUE));
  corrconfigs.push_back(GetConf("MidGapPV52","refGapPos {5} refGapNeg {-5}", kFALSE));
  corrconfigs.push_back(GetConf("MidGapPV52","poiGapPos refGapPos {5} refGapNeg {-5}", kTRUE));

}
