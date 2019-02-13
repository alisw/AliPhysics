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

ClassImp(AliAnalysisTaskGFWFlow);

AliAnalysisTaskGFWFlow::AliAnalysisTaskGFWFlow():
  AliAnalysisTaskSE(),
  debugpar(0),
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
  fWeightPath(""),
  fWeightDir(""),
  fTotFlags(15), //Total number of flags: 1 (nominal) + fTotTrackFlags + N_Event flags
  fTotTrackFlags(8), //Total number of track flags (without nominal)
  fRunNo(-1),
  fCurrSystFlag(0),
  fAddQA(kFALSE),
  fQAList(0)
{
};
AliAnalysisTaskGFWFlow::AliAnalysisTaskGFWFlow(const char *name, Bool_t ProduceWeights, Bool_t IsMC, Bool_t AddQA):
  AliAnalysisTaskSE(name),
  debugpar(0),
  fProduceWeights(ProduceWeights),
  fSelections(0),
  fWeightList(0),
  fWeights(0),
  fExtraWeights(0),
  fFC(0),
  fOutputTree(0),
  fIsMC(IsMC),
  fPtAxis(new TAxis()),
  fWeightPath(""),
  fWeightDir(""),
  fTotFlags(15),
  fTotTrackFlags(8),
  fRunNo(-1),
  fCurrSystFlag(0),
  fAddQA(kFALSE),
  fQAList(0)
{
  DefineOutput(1,(fProduceWeights?TList::Class():AliGFWFlowContainer::Class()));
  if(fAddQA)
    DefineOutput(2,TList::Class());
};
AliAnalysisTaskGFWFlow::~AliAnalysisTaskGFWFlow() {
};
void AliAnalysisTaskGFWFlow::UserCreateOutputObjects(){
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
    OAforPt->Add(new TNamed("MidV36","MidV36"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV36_pt_%i",i+1),"MidV36_pTDiff"));
    OAforPt->Add(new TNamed("MidV42","MidV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV42_pt_%i",i+1),"MidV42_pTDiff"));
    OAforPt->Add(new TNamed("MidV44","MidV44"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidV44_pt_%i",i+1),"MidV44_pTDiff"));

    //2SE:
    OAforPt->Add(new TNamed("Mid2SEV22","Mid2SEV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEV22_pt_%i",i+1),"Mid2SEV22_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEV24","Mid2SEV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEV24_pt_%i",i+1),"Mid2SEV24_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEV26","Mid2SEV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEV26_pt_%i",i+1),"Mid2SEV26_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEV28","Mid2SEV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEV28_pt_%i",i+1),"Mid2SEV28_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEV32","Mid2SEV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEV32_pt_%i",i+1),"Mid2SEV32_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEV34","Mid2SEV34"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEV34_pt_%i",i+1),"Mid2SEV34_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEV36","Mid2SEV36"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEV36_pt_%i",i+1),"Mid2SEV36_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEV42","Mid2SEV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEV42_pt_%i",i+1),"Mid2SEV42_pTDiff"));
    OAforPt->Add(new TNamed("Mid2SEV44","Mid2SEV44"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("Mid2SEV44_pt_%i",i+1),"Mid2SEV44_pTDiff"));

    //|eta|>0.5
    OAforPt->Add(new TNamed("MidGapV22","MidGapV22"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapV22_pt_%i",i+1),"MidGapV22_pTDiff"));
    OAforPt->Add(new TNamed("MidGapV24","MidGapV24"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapV24_pt_%i",i+1),"MidGapV24_pTDiff"));
    OAforPt->Add(new TNamed("MidGapV26","MidGapV26"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapV26_pt_%i",i+1),"MidGapV26_pTDiff"));
    OAforPt->Add(new TNamed("MidGapV28","MidGapV28"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapV28_pt_%i",i+1),"MidGapV28_pTDiff"));
    OAforPt->Add(new TNamed("MidGapV32","MidGapV32"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapV32_pt_%i",i+1),"MidGapV32_pTDiff"));
    OAforPt->Add(new TNamed("MidGapV34","MidGapV34"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapV34_pt_%i",i+1),"MidGapV34_pTDiff"));
    OAforPt->Add(new TNamed("MidGapV36","MidGapV36"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapV36_pt_%i",i+1),"MidGapV36_pTDiff"));
    OAforPt->Add(new TNamed("MidGapV42","MidGapV42"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapV42_pt_%i",i+1),"MidGapV42_pTDiff"));
    OAforPt->Add(new TNamed("MidGapV44","MidGapV44"));
    for(Int_t i=0;i<fPtAxis->GetNbins();i++)
      OAforPt->Add(new TNamed(Form("MidGapV44_pt_%i",i+1),"MidGapV44_pTDiff"));


    //Multi bins:
    Double_t multibins[] = {0,5,10,15,20,30,40,50,70};
    fFC = new AliGFWFlowContainer();
    fFC->SetName(Form("FC%s",fSelections[fCurrSystFlag]->GetSystPF()));
    fFC->Initialize(OAforPt,8,multibins,fCurrSystFlag?1:10); //Statistics only required for nominal profiles, so do not create randomized profiles for systematics
    fGFW = new AliGFW();
    //Full regions
    fGFW->AddRegion("poiMid",10,10,-0.8,0.8,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refMid",10,10,-0.8,0.8,1,1);
    //2 subevets:
    fGFW->AddRegion("poiSENeg",10,10,-0.8,0.,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refSENeg",10,10,-0.8,0.8,1,1);
    fGFW->AddRegion("refSEPos",10,10,0.,0.8,1,1);
    //With gap
    fGFW->AddRegion("poiGapNeg",10,10,-0.8,-0.5,1+fPtAxis->GetNbins(),1);
    fGFW->AddRegion("refGapNeg",10,10,-0.8,-0.5,1,1);
    fGFW->AddRegion("refGapPos",10,10,0.5,0.8,1,1);


  };
  if(fProduceWeights) PostData(1,fWeightList);
  else PostData(1,fFC);
  if(fAddQA) {
    fQAList = new TList();
    fQAList->SetOwner(kTRUE);
    fEventCuts.AddQAplotsToList(fQAList);
    PostData(2,fQAList);
  };
};
void AliAnalysisTaskGFWFlow::UserExec(Option_t*) {
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  if(!fProduceWeights)
    if(!InitRun()) return;
  if(!AcceptEvent()) return;
  if(!AcceptAODVertex(fAOD)) return;
  if(fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if(!fMCEvent)
      return;
  };
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t cent = lMultSel->GetMultiplicityPercentile("V0M");
  if(fCurrSystFlag==fTotTrackFlags+4) cent = lMultSel->GetMultiplicityPercentile("CL1"); //CL1 flag is EvFlag 4 = N_TrackFlags + 4
  if(fCurrSystFlag==fTotTrackFlags+5) cent = lMultSel->GetMultiplicityPercentile("CL0"); //CL0 flag is EvFlag 5 = N_TrackFlags + 5
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  Double_t vtxb = GetVtxBit(fAOD);
  if(!vtxb) return; //If no vertex pass, then do not consider further
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
      if(!fWeights) printf("Weights do not exist!\n");
      Double_t nua = fWeights->GetWeight(lTrack->Phi(),lTrack->Eta(),vz,lTrack->Pt(),cent,0);
      //Double_t nuaITS = fExtraWeights->GetWeight(lTrack->Phi(),lTrack->Eta(),vz,lTrack->Pt(),cent,0);
      Double_t nue = fPtAxis->GetNbins()>1?1:fWeights->GetWeight(lTrack->Phi(),lTrack->Eta(),vz,cent,lTrack->Pt(),1);
      if(fSelections[fCurrSystFlag]->AcceptTrack(lTrack, lDCA))
	fGFW->Fill(lTrack->Eta(),fPtAxis->FindBin(lTrack->Pt())-1,lTrack->Phi(),nua*nue,1);
      /*if(fSelections[9]->AcceptTrack(lTrack, lDCA)) //No ITS for now
	fGFW->Fill(lTrack->Eta(),fPtAxis->FindBin(lTrack->Pt())-1,lTrack->Phi(),nuaITS*nue,2);*/ 
    };
    TRandom rndm(0);
    Double_t rndmn=rndm.Rndm();
    //Calculate & fill profiles:
    //V_2{n}, full acceptance
    Bool_t filled = FillFCs("MidV22","refMid {2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("MidV22","poiMid refMid {2 -2}", cent, kTRUE,rndmn);
    filled = FillFCs("MidV24","refMid {2 2 -2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("MidV24","poiMid refMid {2 2 -2 -2}", cent, kTRUE,rndmn);
    filled = FillFCs("MidV26","refMid {2 2 2 -2 -2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("MidV26","poiMid refMid {2 2 2 -2 -2 -2}", cent, kTRUE,rndmn);
    filled = FillFCs("MidV28","refMid {2 2 2 2 -2 -2 -2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("MidV28","poiMid refMid {2 2 2 2 -2 -2 -2 -2}", cent, kTRUE,rndmn);
    //V_3{n}, full acceptance:
    filled = FillFCs("MidV32","refMid {3 -3}", cent, kFALSE,rndmn);
    filled = FillFCs("MidV32","poiMid refMid {3 -3}", cent, kTRUE,rndmn);
    filled = FillFCs("MidV34","refMid {3 3 -3 -3}", cent, kFALSE,rndmn);
    filled = FillFCs("MidV34","poiMid refMid {3 3 -3 -3}", cent, kTRUE,rndmn);
    filled = FillFCs("MidV36","refMid {2 2 -2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("MidV36","poiMid refMid {2 2 -2 -2}", cent, kTRUE,rndmn);
    //V_4{n}, full acceptance:
    filled = FillFCs("MidV42","refMid {4 -4}", cent, kFALSE,rndmn);
    filled = FillFCs("MidV42","poiMid refMid {4 -4}", cent, kTRUE,rndmn);
    filled = FillFCs("MidV44","refMid {4 4 -4 -4}", cent, kFALSE,rndmn);
    filled = FillFCs("MidV44","poiMid refMid {4 4 -4 -4}", cent, kTRUE,rndmn);
    //V_2{n}, 2 subevents:
    filled = FillFCs("Mid2SEV22","refSENeg {2} refSEPos {-2}", cent, kFALSE,rndmn);
    filled = FillFCs("Mid2SEV22","poiSENeg refSENeg {2} refSEPos {-2}", cent, kTRUE,rndmn);
    filled = FillFCs("Mid2SEV24","refSENeg {2 2} refSEPos {-2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("Mid2SEV24","poiSENeg refSENeg {2 2} refSEPos {-2 -2}", cent, kTRUE,rndmn);
    filled = FillFCs("Mid2SEV26","refSENeg {2 2 2} refSEPos {-2 -2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("Mid2SEV26","poiSENeg refSENeg {2 2 2} refSEPos {-2 -2 -2}", cent, kTRUE,rndmn);
    filled = FillFCs("Mid2SEV28","refSENeg {2 2 2 2} refSEPos {-2 -2 -2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("Mid2SEV28","poiSENeg refSENeg {2 2 2 2} refSEPos {-2 -2 -2 -2}", cent, kTRUE,rndmn);
    //V_3{n}, 2 subevents:
    filled = FillFCs("Mid2SEV32","refSENeg {3} refSEPos {-3}", cent, kFALSE,rndmn);
    filled = FillFCs("Mid2SEV32","poiSENeg refSENeg {3} refSEPos {-3}", cent, kTRUE,rndmn);
    filled = FillFCs("Mid2SEV34","refSENeg {3 3} refSEPos {-3 -3}", cent, kFALSE,rndmn);
    filled = FillFCs("Mid2SEV34","poiSENeg refSENeg {3 3} refSEPos {-3 -3}", cent, kTRUE,rndmn);
    filled = FillFCs("Mid2SEV36","refSENeg {3 3 3} refSEPos {-3 -3 -3}", cent, kFALSE,rndmn);
    filled = FillFCs("Mid2SEV36","poiSENeg refSENeg {3 3 3} refSEPos {-3 -3 -3}", cent, kTRUE,rndmn);
    //V_4{n}, 2 subevents:
    filled = FillFCs("Mid2SEV42","refSENeg {4} refSEPos {-4}", cent, kFALSE,rndmn);
    filled = FillFCs("Mid2SEV42","poiSENeg refSENeg {4} refSEPos {-4}", cent, kTRUE,rndmn);
    filled = FillFCs("Mid2SEV44","refSENeg {4 4} refSEPos {-4 -4}", cent, kFALSE,rndmn);
    filled = FillFCs("Mid2SEV44","poiSENeg refSENeg {4 4} refSEPos {-4 -4}", cent, kTRUE,rndmn);

    //V_2{n}, eta gap 1:
    filled = FillFCs("MidGapV22","refGapNeg {2} refGapPos {-2}", cent, kFALSE,rndmn);
    filled = FillFCs("MidGapV22","poiGapNeg refGapNeg {2} refGapPos {-2}", cent, kTRUE,rndmn);
    filled = FillFCs("MidGapV24","refGapNeg {2 2} refGapPos {-2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("MidGapV24","poiGapNeg refGapNeg {2 2} refGapPos {-2 -2}", cent, kTRUE,rndmn);
    filled = FillFCs("MidGapV26","refGapNeg {2 2 2} refGapPos {-2 -2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("MidGapV26","poiGapNeg refGapNeg {2 2 2} refGapPos {-2 -2 -2}", cent, kTRUE,rndmn);
    filled = FillFCs("MidGapV28","refGapNeg {2 2 2 2} refGapPos {-2 -2 -2 -2}", cent, kFALSE,rndmn);
    filled = FillFCs("MidGapV28","poiGapNeg refGapNeg {2 2 2 2} refGapPos {-2 -2 -2 -2}", cent, kTRUE,rndmn);
    //V_3{n}, eta gap 1:
    filled = FillFCs("MidGapV32","refGapNeg {3} refGapPos {-3}", cent, kFALSE,rndmn);
    filled = FillFCs("MidGapV32","poiGapNeg refGapNeg {3} refGapPos {-3}", cent, kTRUE,rndmn);
    filled = FillFCs("MidGapV34","refGapNeg {3 3} refGapPos {-3 -3}", cent, kFALSE,rndmn);
    filled = FillFCs("MidGapV34","poiGapNeg refGapNeg {3 3} refGapPos {-3 -3}", cent, kTRUE,rndmn);
    filled = FillFCs("MidGapV36","refGapNeg {3 3 3} refGapPos {-3 -3 -3}", cent, kFALSE,rndmn);
    filled = FillFCs("MidGapV36","poiGapNeg refGapNeg {3 3 3} refGapPos {-3 -3 -3}", cent, kTRUE,rndmn);
    //V_4{n}, eta gap 1:
    filled = FillFCs("MidGapV42","refGapNeg {4} refGapPos {-4}", cent, kFALSE,rndmn);
    filled = FillFCs("MidGapV42","poiGapNeg refGapNeg {4} refGapPos {-4}", cent, kTRUE,rndmn);
    filled = FillFCs("MidGapV44","refGapNeg {4 4} refGapPos {-4 -4}", cent, kFALSE,rndmn);
    filled = FillFCs("MidGapV44","poiGapNeg refGapNeg {4 4} refGapPos {-4 -4}", cent, kTRUE,rndmn);

    PostData(1,fFC);
    if(fAddQA) PostData(2,fQAList);
  }; 
};
void AliAnalysisTaskGFWFlow::Terminate(Option_t*) {
};
Bool_t AliAnalysisTaskGFWFlow::AcceptEvent() {
  if(!fEventCuts.AcceptEvent(fInputEvent)) return 0;
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(!(fSelMask&AliVEvent::kINT7)) return 0;
  return kTRUE;
};

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
Bool_t AliAnalysisTaskGFWFlow::AcceptAODTrack(AliAODTrack *mtr, TClonesArray *tca) {
  if(TMath::Abs(mtr->Eta())>1.6) return kFALSE;
  if(mtr->Pt()<0.2) return kFALSE;
  if(mtr->Pt()>20) return kFALSE;
  if(!mtr->TestFilterBit(2+96+768)) return kFALSE;
  //if(mtr->GetTPCNclsF()<70) return kFALSE; //This will be done in the Selection object
  if(tca) {
    Int_t mcind = TMath::Abs(mtr->GetLabel());
    AliAODMCParticle *lp = (AliAODMCParticle*)tca->At(mcind);
    if(!lp->IsPhysicalPrimary()) return kFALSE;
    if(lp->Charge()==0) return kFALSE;
  };
  return kTRUE;
};
Bool_t AliAnalysisTaskGFWFlow::AcceptParticle(AliVParticle *mpa) {
  if(!mpa->IsPhysicalPrimary()) return kFALSE;
  if(mpa->Charge()==0) return kFALSE;
  if(TMath::Abs(mpa->Eta())>0.8) return kFALSE;
  if(mpa->Pt()<0.2) return kFALSE;
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
    fWeights = (AliGFWWeights*)fWeightList->FindObject(Form("%i",runno));
    if(!fWeights) {
      return kFALSE;
    };
    return kTRUE;
  }; //If weights not set, attempting to fetch them from pre-set directory. This will definitely fail if running on train
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
    Bool_t dump1 = AcceptEvent(); //To setup all the event cuts (automatically)
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
    val = fGFW->Calculate(hn).Re();
    fFC->FillProfile(head.Data(),cent,val/dnx,dnx,rndmn);
    return kTRUE;
  };
  for(Int_t i=1;i<=fPtAxis->GetNbins();i++) {
    TString tss(hn);
    tss.Prepend(Form("(%i) ",i-1));
    dnx = fGFW->Calculate(tss,kTRUE).Re();
    if(dnx==0) continue;
    val = fGFW->Calculate(tss).Re();
    fFC->FillProfile(Form("%s_pt_%i",head.Data(),i),cent,val/dnx,dnx,rndmn);
  };
  return kTRUE;
};

