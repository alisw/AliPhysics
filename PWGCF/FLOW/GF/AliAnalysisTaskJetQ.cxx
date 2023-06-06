#include "AliAnalysisTaskJetQ.h"
#include "AliAODForwardMult.h"
AliAnalysisTaskJetQ::AliAnalysisTaskJetQ():
 AliAnalysisTaskSE(),
 fAOD(0),
 fMCEvent(0),
 fPoolMgr(0),
 fPoolTrackArray(0),
 fTriggerType(AliVEvent::kINT7),
 fOutList(0),
 fCentAxis(0),
 fVzAxis(0),
 fPtAxis(0),
 fNormCounter(0),
 fCorrPlot(0),
 fMixCorrPlot(0),
 fNtriggers(0),
 fHMaxPt(0),
 fPtAssocMin(2.),
 fPtAssocMax(2.5),
 fPtTriggMin(6.),
 fPtTriggMax(8.),
 fCalculateFlow(kFALSE),
 fRndmGen(0),
 fFCIncl(0),
 fFCTrig(0),
 fGFW(0),
 fFMDHist(0)
{}
//_____________________________________________________________________________
AliAnalysisTaskJetQ::AliAnalysisTaskJetQ(const char* name, Bool_t lCalcFlow):
 AliAnalysisTaskSE(name),
 fAOD(0),
 fMCEvent(0),
 fPoolMgr(0),
 fPoolTrackArray(0),
 fTriggerType(AliVEvent::kINT7),
 fOutList(0),
 fCentAxis(0),
 fVzAxis(0),
 fPtAxis(0),
 fNormCounter(0),
 fCorrPlot(0),
 fMixCorrPlot(0),
 fNtriggers(0),
 fHMaxPt(0),
 fPtAssocMin(2.),
 fPtAssocMax(2.5),
 fPtTriggMin(6.),
 fPtTriggMax(8.),
 fCalculateFlow(lCalcFlow),
 fRndmGen(0),
 fFCIncl(0),
 fFCTrig(0),
 fGFW(0),
 fFMDHist(0)
{
    DefineOutput(1, TList::Class());
    DefineOutput(2, TH1D::Class());
    if(fCalculateFlow) {
      DefineOutput(3, AliGFWFlowContainer::Class());
      DefineOutput(4, AliGFWFlowContainer::Class());
    }
}
//_____________________________________________________________________________
AliAnalysisTaskJetQ::~AliAnalysisTaskJetQ()
{
  delete fPoolTrackArray;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetQ::UserCreateOutputObjects()
{
    OpenFile(1);
    //Setting up variables
    if(fCentBins.size()<2) {
      Double_t tempCentBins[] = {0,100};
      SetCentralityBins(1,tempCentBins);
    };
    if(fVzBins.size()<2) {
      Double_t vzs[]   = {-10.,10.};
      SetVtxZBins(1,vzs);
    };
    if(fPtBins.size()<2) {
      Double_t ptbs[] = {2.,2.5};
      SetPtBins(1,ptbs);
    };


    fPtAssocMin = fPtAxis->GetBinLowEdge(1);
    fPtAssocMax = fPtAxis->GetBinUpEdge(fPtAxis->GetNbins());
    fPtDif = fPtBins.size()>2; //if only one pT bin, then assoc is not pt-dif. and use TH2 instead of TH3
    if(!fEvMixPars[0]) SetEventMixingCapacity(30,5000,10,10);
    //Setting up the pool manager
    // Double_t cents[] = {0.,50.,70.};
    fPoolMgr = new AliEventPoolManager(fEvMixPars[0], fEvMixPars[1], fCentBins.size()-1, fCentBins.data(), fVzBins.size()-1, fVzBins.data());
      //fNCentBins,fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    if (!fPoolMgr) { AliError("Event Pool manager not created!"); return; }
    fPoolMgr->SetTargetValues(fEvMixPars[1], fEvMixPars[2]*0.01, fEvMixPars[3]);
    fPoolTrackArray = new TObjArray();
    fPoolTrackArray->SetOwner(kTRUE);
    // creating output lists
    fOutList = new TList();
    fOutList->SetOwner(kTRUE);
    TH1 *l_ptBins = new TH1D("ptBins","ptBins",fPtBins.size()-1,fPtBins.data());
    TH1 *l_centBins = new TH1D("centBins","centBins",fCentBins.size()-1,fCentBins.data());
    TH1D *vtzBefore = new TH1D("vtzBefore","vtzBefore",20,-10,10);
    TH1D *vtzAfter  = new TH1D("vtzAfter","vtzAfter",20,-10,10);
    fNormCounter    = new TH2D("NormCounter","NormCounter; multi/cent; index",fCentBins.size()-1, fCentBins.data(), 4, -0.5, 3.5);
    fOutList->Add(vtzBefore);
    fOutList->Add(vtzAfter);
    fOutList->Add(fNormCounter);
    fOutList->Add(l_ptBins);
    fOutList->Add(l_centBins);
    //Setting up correlation plots
    fCorrPlot    = new TH1**[fCentAxis->GetNbins()];
    fMixCorrPlot = new TH1**[fCentAxis->GetNbins()];
    fNtriggers   = new TH1*[fCentAxis->GetNbins()];
    for(Int_t iCent=0;iCent<fCentAxis->GetNbins();iCent++) {
      fCorrPlot[iCent] = new TH1*[fPtAxis->GetNbins()];
      fMixCorrPlot[iCent] = new TH1*[fPtAxis->GetNbins()];
      for(Int_t iPt=0;iPt<fPtAxis->GetNbins();iPt++) {
        fCorrPlot[iCent][iPt]       = new TH3D(Form("fCorr_Cent%i_Pt%i",iCent,iPt),Form("fCorr_Cent%i_Pt%i; #Delta#phi (rad); #Delta#eta; z_{vtx}",iCent,iPt),100,-C_PI_HALF,C_PI_TH,200, -1.6, 1.6,fVzBins.size()-1,-10,10);
        fMixCorrPlot[iCent][iPt]    = new TH3D(Form("fMixCorr_Cent%i_Pt%i",iCent,iPt),Form("fMixCorr_Cent%i_Pt%i; #Delta#phi (rad); #Delta#eta; z_{vtx}",iCent,iPt),100,-C_PI_HALF,C_PI_TH,200, -1.6, 1.6,fVzBins.size()-1,-10,10);
        fCorrPlot[iCent][iPt]->GetZaxis()->Set(fVzBins.size()-1,fVzBins.data());
        fMixCorrPlot[iCent][iPt]->GetZaxis()->Set(fVzBins.size()-1,fVzBins.data());
        fOutList->Add(fCorrPlot[iCent][iPt]);
        fOutList->Add(fMixCorrPlot[iCent][iPt]);
      };
      fNtriggers[iCent]    = new TH1D(Form("fNtriggers_Cent%i",iCent),Form("fNtriggers_Cent%i; v_{z}; dN_{trig}/dv_{z}",iCent),fVzBins.size()-1,-10,10);
      fNtriggers[iCent]->GetZaxis()->Set(fVzBins.size()-1,fVzBins.data());
      fOutList->Add(fNtriggers[iCent]);

    };
    PostData(1, fOutList);
    //Prepare a histogram to keep highest pT track. Should be outside the fOutList, because the list only gets posted for events with trigger track
    Int_t lnTrigBins = TMath::Nint(fPtTriggMax/0.1)+1; //One extra bin for "overflow"
    Double_t *lPtBinsForTrig = new Double_t[lnTrigBins+1]; //And one extra value for upper bin of overflow
    for(Int_t i=0;i<=lnTrigBins;i++) lPtBinsForTrig[i] = i*0.1;
    fHMaxPt = new TH1D("MaxPt","MaxPt; #it{p}_{T, max}; d#it{N}/d#it{p}_{T, max}",lnTrigBins,lPtBinsForTrig);
    PostData(2, fHMaxPt);
    if(fCalculateFlow) {
      SetupFlowOutput();
      PostData(3,fFCIncl);
      PostData(4,fFCTrig);
    };
}
//_____________________________________________________________________________
void AliAnalysisTaskJetQ::UserExec(Option_t *)
{
  Bool_t fIsMC = kFALSE;
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  if(fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCEvent) return;
  }
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t l_Cent = lMultSel->GetMultiplicityPercentile("V0M");
  if(!CheckTrigger(l_Cent)) return;
  Double_t vtxXYZ[] = {0.,0.,0.};
  if(!AcceptAOD(fAOD, vtxXYZ)) return;
  // Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  ((TH1D*)fOutList->At(0))->Fill(vz);
  Int_t ind = FindGivenPt(fPtTriggMin,fPtTriggMax);
  Int_t i_Cent = fCentAxis->FindBin(l_Cent);
  Int_t i_vz   = fVzAxis->FindBin(vz);
  if(!i_Cent || i_Cent>fCentAxis->GetNbins()) return; //out of centrality/vz range
  if(!i_vz || i_vz>fVzAxis->GetNbins()) return; //out of centrality/vz range
  //Need to fill flow containers here, b/c I also want to calculate it in events where no trigger track is present
  if(fCalculateFlow) {
    Double_t rndm = fRndmGen->Rndm();
    FillFCs(l_Cent,rndm,ind>=0); //Both FCs are filled to save CPU time, but only inclusive is posted for now
    PostData(3,fFCIncl);
  };
  if(ind<0) return;
  //Test for FMD
  AliAODForwardMult* aodForward=static_cast<AliAODForwardMult*>(fAOD->FindListObject("Forward"));
  if(!aodForward) { printf("\n\n\n\n\n\n**************************\nFMD stuff not found!\n\n\n\n\n\n"); return; }
  else {
    const TH2D& d2Ndetadphi = aodForward->GetHistogram();
    if(!fFMDHist) fFMDHist = (TH2*)d2Ndetadphi.Clone("FMDDist");
    else fFMDHist->Add(&d2Ndetadphi);
  }
  //End of test
  ((TH1D*)fOutList->At(1))->Fill(vz);
  fNormCounter->Fill(l_Cent,0); //Number of triggers
  Int_t nPairs = FillCorrelations(ind,i_Cent,vz);
  // PostData(4, fFCTrig);
  AliEventPool *pool = fPoolMgr->GetEventPool(i_Cent-1, i_vz-1);
  if(!pool) { printf("Could not find the event pool!\n"); return; };
  // printf("Current numbe of events in pool: %i\n",pool->GetCurrentNEvents());
  if(pool->IsReady()) { Int_t nMixPairs = FillMixedEvent(ind, pool, i_Cent, vz);  };
  pool->UpdatePool((TObjArray*)fPoolTrackArray->Clone());
  fPoolTrackArray->Clear();
  PostData(1, fOutList);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetQ::Terminate(Option_t *)
{
  fFMDHist->Draw("colz");
}
//_____________________________________________________________________________

Bool_t AliAnalysisTaskJetQ::CheckTrigger(Double_t l_Cent) {
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  //Apparently, MB trigger can also mark special triggers, leaving depleted regions in multi. To avoid this, pass true, if MB has been triggered.
  //This would fail if spec. triggers would also flag MB trigger, which seems to NOT be the case.
  if(!(fTriggerType&fSelMask)) { return kFALSE; }; //printf("Returning from the generic check\n");
  if(fSelMask&(fTriggerType&(AliVEvent::kINT7+AliVEvent::kMB))) {return kTRUE; }; //printf("Passed by MB trigger!\n");
  if((fSelMask&fTriggerType&AliVEvent::kCentral) && l_Cent>10) {return kFALSE; }; //printf("Returnning from kCent case\n");
  if((fSelMask&fTriggerType&AliVEvent::kSemiCentral) && (l_Cent<30 || l_Cent>50)) {return kFALSE; }; //printf("Returning from kSC case\n");
  return kTRUE;
}
Bool_t AliAnalysisTaskJetQ::AcceptAOD(AliAODEvent *inEv, Double_t *lvtxXYZ) {
  if(!fEventCuts.AcceptEvent(inEv)) return 0;
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
  vtx->GetXYZ(lvtxXYZ);
  return kTRUE;
};
Int_t AliAnalysisTaskJetQ::FindGivenPt(const Double_t &ptMin, const Double_t &ptMax) {
  Int_t ind=-1;
  AliAODTrack *lTrack;
  Double_t lPtMax=0;
  Double_t lPtCur=0;
  Double_t lAbsMaxPt=0;
  if(fCalculateFlow) fGFW->Clear(); //Need to clear up the GFW before filling it in
  Int_t ptBinMax = fPtAxis->GetNbins();
  for(Int_t i=0;i<fAOD->GetNumberOfTracks();i++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(i);
    if(!lTrack->TestFilterBit(768)) continue;
    lPtCur = lTrack->Pt();
    if(lAbsMaxPt < lPtCur) lAbsMaxPt = lPtCur;
    if(lPtCur > ptMin && lPtCur < ptMax && lPtCur > lPtMax) { lPtMax=lPtCur; ind=i;  };
    //Fill in the GFW here:
    if(fCalculateFlow) {
      Int_t ptind = fPtAxis->FindBin(lPtCur);
      if(!ptind) continue; //Beyond the lower end of our pt: throw away
      if(lPtCur>ptMax) continue; //If above our cutoff, throw away, too
      if(ptind>ptBinMax && lPtCur<ptMin) continue; //Don't consider particles that are inbetween pT regime & trigger
      fGFW->Fill(lTrack->Eta(),ptind-1,lTrack->Phi(),1,1); //Use weights 1 here
      //For debugging purposes
      // if(lPtCur > 8 && lPtCur < 15) fGFW->Fill(lTrack->Eta(),0,lTrack->Phi(),1,4);
      // if(lPtCur > 0.5 && lPtCur < 2.) fGFW->Fill(lTrack->Eta(),0,lTrack->Phi(),1,2);
    }
  };
  fHMaxPt->Fill((lAbsMaxPt>fPtTriggMax)?(fPtTriggMax+0.01):lAbsMaxPt); //In case it's overflow, then store it in appropriate bin
  PostData(2,fHMaxPt);
  return ind;
}
Int_t AliAnalysisTaskJetQ::FillCorrelations(Int_t &triggerIndex, Int_t &centVal, Double_t &vzValue) {
  Int_t nPairs=0;
  AliAODTrack *lTrack;
  AliAODTrack *lTriggerTrack = (AliAODTrack*)fAOD->GetTrack(triggerIndex);
  Double_t l_TrEta = lTriggerTrack->Eta();
  Double_t l_TrPhi = lTriggerTrack->Phi();
  Int_t centBin = centVal-1;
  fNtriggers[centBin]->Fill(vzValue);
  //Fetch pT trigger & make sure we are not in the same pT bin below
  lTrack = (AliAODTrack*)fAOD->GetTrack(triggerIndex);
  Int_t trigBin = fPtAxis->FindBin(lTrack->Pt());
  for(Int_t i=0;i<fAOD->GetNumberOfTracks();i++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(i);
    if(!lTrack->TestFilterBit(768)) continue;
    Int_t ptInd = fPtAxis->FindBin(lTrack->Pt());
    if(!ptInd || ptInd>fPtAxis->GetNbins() || ptInd == trigBin) continue;
    ptInd--;
    Double_t d_Eta = l_TrEta-lTrack->Eta();
    Double_t d_Phi = l_TrPhi-lTrack->Phi();
    fixPhi(d_Phi);
    // printf("Attempting to fill a track with %f %f %f (%s)\n",d_Phi,d_Eta,lPt,fPtDif);
    fill3DHist(fCorrPlot[centBin][ptInd],d_Phi,d_Eta,vzValue);
    nPairs++;
    fPoolTrackArray->Add(new AliBasicParticle(lTrack->Eta(),lTrack->Phi(),lTrack->Pt(), lTrack->Charge()));
  };
  return nPairs;
};
Int_t AliAnalysisTaskJetQ::FillMixedEvent(Int_t &triggerIndex, AliEventPool *l_pool, Int_t &centVal, Double_t &vzValue) {
  Int_t nPairs=0;
  AliAODTrack *lTrack;
  AliAODTrack *lTriggerTrack = (AliAODTrack*)fAOD->GetTrack(triggerIndex);
  Double_t l_TrEta = lTriggerTrack->Eta();
  Double_t l_TrPhi = lTriggerTrack->Phi();
  Double_t l_AsEta;
  Double_t l_AsPhi;
  Double_t l_AsPt;
  TObjArray *mixTracks;
  Int_t centInd = centVal-1;
  //Fetch pT trigger & make sure we are not in the same pT bin below
  lTrack = (AliAODTrack*)fAOD->GetTrack(triggerIndex);
  Int_t trigBin = fPtAxis->FindBin(lTrack->Pt());
  for(Int_t i=0;i<l_pool->GetCurrentNEvents();i++) {
    mixTracks = l_pool->GetEvent(i);
    TObjArrayIter *myIt = new TObjArrayIter(mixTracks);
    while(AliBasicParticle *rTr = (AliBasicParticle*)myIt->Next()) {
      Int_t ptInd = fPtAxis->FindBin(rTr->Pt());
      if(!ptInd || ptInd>fPtAxis->GetNbins() || ptInd==trigBin) continue;
      ptInd-=1;
      l_AsEta = l_TrEta - rTr->Eta();
      l_AsPhi = l_TrPhi - rTr->Phi();
      fixPhi(l_AsPhi);
      fill3DHist(fMixCorrPlot[centInd][ptInd],l_AsPhi,l_AsEta,vzValue);
      nPairs++;
    };
  };
  return nPairs;
};
void AliAnalysisTaskJetQ::SetupFlowOutput() {
  if(!fGFW) fGFW = new AliGFW();
  //Setting up regions
  Int_t powers[] = {3,2,2,2};
  fGFW->AddRegion("pos",4,powers,0.5,0.8,1+fPtAxis->GetNbins(),1); //In principle, should be enough with fPtAxis->GetNbins()?
  fGFW->AddRegion("neg",4,powers,-0.8,-0.5,1+fPtAxis->GetNbins(),1); //In principle, should be enough with fPtAxis->GetNbins()?
  //Debugging purposes
  // fGFW->AddRegion("testN",4,powers,-0.8,-0.5,1,2); //In principle, should be enough with fPtAxis->GetNbins()?
  // fGFW->AddRegion("testP",4,powers,0.5,0.8,1,4); //In principle, should be enough with fPtAxis->GetNbins()?
  //Preparing AliGFWFlowContainer:
  //First, the structure:
  TObjArray *oba = new TObjArray();
  for(Int_t hr=1;hr<=3;hr++)
    for(Int_t i=0;i<=fPtAxis->GetNbins();i++) {
      //Also ordering so that same configurations go next to each other, so that postprocessing is easier
      oba->Add(new TNamed(Form("GapPV%iRef%i",hr,i),Form("GapPV%iRef%i",hr,i)));
      for(Int_t j=i+1;j<=fPtAxis->GetNbins();j++)
        oba->Add(new TNamed(Form("GapPV%iRef%i_pt_%i",hr,i,j),Form("GapPV%iRef%i_pt_%i",hr,i,j)));

      oba->Add(new TNamed(Form("GapNV%iRef%i",hr,i),Form("GapNV%iRef%i",hr,i)));
      for(Int_t j=i+1;j<=fPtAxis->GetNbins();j++)
        oba->Add(new TNamed(Form("GapNV%iRef%i_pt_%i",hr,i,j),Form("GapNV%iRef%i_pt_%i",hr,i,j)));

    };
  //Debugging purposes
  // oba->Add(new TNamed("TestNP","TestNP"));
  //Then also calculate centrality bins
  Double_t *lCentBins = new Double_t[fCentAxis->GetNbins()+1];
  for(Int_t i=0;i<fCentAxis->GetNbins();i++) lCentBins[i] = fCentAxis->GetBinLowEdge(i+1);
  lCentBins[fCentAxis->GetNbins()] = fCentAxis->GetBinUpEdge(fCentAxis->GetNbins());
  //Initialize flow containers
  //First for inclusive (despite the trigger track presence)
  fFCIncl = new AliGFWFlowContainer();
  fFCIncl->SetName("FCInclusive");
  fFCIncl->SetXAxis(fPtAxis);
  fFCIncl->Initialize(oba,fCentAxis->GetNbins(),lCentBins,10); //Also probably let's do the bootstrapping
  //Also, for events with trigger track
  fFCTrig = new AliGFWFlowContainer();
  fFCTrig->SetName("FCTrigger");
  fFCTrig->SetXAxis(fPtAxis);
  fFCTrig->Initialize(oba,fCentAxis->GetNbins(),lCentBins,10); //Also probably let's do the bootstrapping

    //Create correlator configurations
  for(Int_t hr=1;hr<=3;hr++)
    for(Int_t i=0;i<=fPtAxis->GetNbins();i++) {
      corrconfigs.push_back(GetConf(Form("GapPV%iRef%i",hr,i),Form("pos (%i) {%i} neg (%i) {%i}",i,hr,i,-hr), kFALSE));
      corrconfigs.push_back(GetConf(Form("GapPV%iRef%i",hr,i),Form("pos {%i} neg (%i) {%i}",hr,i,-hr), kTRUE));
      corrconfigs.push_back(GetConf(Form("GapNV%iRef%i",hr,i),Form("neg (%i) {%i} pos (%i) {%i}",i,hr,i,-hr), kFALSE));
      corrconfigs.push_back(GetConf(Form("GapNV%iRef%i",hr,i),Form("neg {%i} pos (%i) {%i}",hr,i,-hr), kTRUE));
    };
  //Also make random number gen. for bootstrapping
  fRndmGen = new TRandom(0);
  //Debugging purposes
  // corrconfigs.push_back(GetConf("TestNP","testN {3} testP {-3}",kFALSE));
}
//Have to be very careful here to fill the correct pt bins!
void AliAnalysisTaskJetQ::FillFCs(AliGFW::CorrConfig corconf, Double_t cent, Double_t rndmn, const Bool_t &TrFound) {
  Double_t dnx, val;
  dnx = fGFW->Calculate(corconf,0,kTRUE).real();
  if(dnx==0) return;
  if(!corconf.pTDif) {
    val = fGFW->Calculate(corconf,0,kFALSE).real()/dnx;
    if(TMath::Abs(val)<1) {
      fFCIncl->FillProfile(corconf.Head.c_str(),cent,val,dnx,rndmn);
      if(TrFound) fFCTrig->FillProfile(corconf.Head.c_str(),cent,val,dnx,rndmn);
    };
    return;
  };
  //It seems that here filling I end up filling non-existing bins. Try to check in GFWFlowContainer what is (attempted) to access
  for(Int_t i=corconf.ptInd.back()+1;i<=fPtAxis->GetNbins();i++) { //It's kind of a triangle. If ref flow starts at ptInd 1, then (1,1) is calculated as integrated, so we start from (1,2)
    dnx = fGFW->Calculate(corconf,i,kTRUE).real();
    if(dnx==0) continue;
    val = fGFW->Calculate(corconf,i,kFALSE).real()/dnx;
    if(TMath::Abs(val)<1) {
      fFCIncl->FillProfile(Form("%s_pt_%i",corconf.Head.c_str(),i),cent,val,dnx,rndmn);
      if(TrFound) fFCTrig->FillProfile(Form("%s_pt_%i",corconf.Head.c_str(),i),cent,val,dnx,rndmn);
    };
  };
};
