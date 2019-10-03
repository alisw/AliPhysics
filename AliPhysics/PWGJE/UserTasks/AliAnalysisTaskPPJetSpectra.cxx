#include "AliAnalysisTaskSE.h"
#include "AliAnalysisHelperJetTasks.h"
#include "THnSparse.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "TAxis.h"
#include "AliAODMCParticle.h"


#include "AliAnalysisTaskPPJetSpectra.h"
#include "AliAnalysisHelperJetTasks.h"
#include <iostream>

ClassImp(AliAnalysisTaskPPJetSpectra)

//---------------------------------------------------------------------------------------------------
AliAnalysisTaskPPJetSpectra::AliAnalysisTaskPPJetSpectra()
 :AliAnalysisTaskSE()
  ,fOutputList(0)
  ,fESD(0)
  ,fAOD(0)
  ,fAODIn(0)
  ,fAODOut(0)
  ,fAODExt(0)
  ,fNonStdFile("")
  ,fDebug(0)
  ,fUseMC(kFALSE)
  ,fEvtSelectionMask(0)
  ,fEventClass(-1)
  ,nVtxContCut(2)
  ,fVtxZcut(10)
  ,fVtxRcut(1)
  ,nTrackFilter(272)
  ,trackPtMin(0.15)
  ,trackPtMax(100.)
  ,trackEtaAbsMax(0.9)
  ,jetPtMin(0.)
  ,jetEtaCut(0.)
  ,jetZmax(0.)
  ,fhnEvent(0)
  ,fhnTracks(0)
  ,fhnMC(0)
  ,fhnMC2(0)
  ,fhnRecJetsNoCut(0)
  ,fhnGenJetsNoCut(0)
  ,fhnRecJetsCut(0)
  ,fhnGenJetsCut(0)
  ,fhnRecBckg(0)
  ,fhnGenBckg(0)
  ,fhnRecJetsTrackUEcor(0)
  ,fhnGenJetsTrackUEcor(0)
  ,fhnRecJetsBckgUEcor(0)
  ,fhnGenJetsBckgUEcor(0)
  ,fhnTrackUE(0)
  ,fhnParticleUE(0)
  ,fhnBckgRecUE(0)
  ,fhnBckgGenUE(0)
  ,fRecJetBranch("")
  ,fGenJetBranch("")
  ,fRecBckgBranch("")
  ,fGenBckgBranch("")
  ,fRecJetR(0.)
  ,fGenJetR(0.)
  ,fRecBckgR(0.)
  ,fGenBckgR(0.)
  ,fTrackType(0)
  ,fParticleType(0)
  ,fhnMatching(0)
  ,fhnTrackCorrMatching(0)
  ,fhnBckgCorrMatching(0)
  ,fhnTrackUEanal(0)
  ,kDoUEanalysis(kFALSE)
  ,fRejectPileUp(1)
{
  if(fDebug) printf("%s: %d  Constructor\n",(char*)__FILE__, __LINE__);
  //DefineOutput(1, TList::Class());
}

//---------------------------------------------------------------------------------------------------
AliAnalysisTaskPPJetSpectra::AliAnalysisTaskPPJetSpectra(const char* name):AliAnalysisTaskSE(name)
  ,fOutputList(0)
  ,fESD(0)
  ,fAOD(0)
  ,fAODIn(0)
  ,fAODOut(0)
  ,fAODExt(0)
  ,fNonStdFile("")
  ,fDebug(0)
  ,fUseMC(kFALSE)
  ,fEvtSelectionMask(0)
  ,fEventClass(-1)
  ,nVtxContCut(2)
  ,fVtxZcut(10)
  ,fVtxRcut(1)
  ,nTrackFilter(272)
  ,trackPtMin(0.15)
  ,trackPtMax(100.)
  ,trackEtaAbsMax(0.9)
  ,jetPtMin(0.)
  ,jetEtaCut(0.)
  ,jetZmax(0.)
  ,fhnEvent(0)
  ,fhnTracks(0)
  ,fhnMC(0)
  ,fhnMC2(0)
  ,fhnRecJetsNoCut(0)
  ,fhnGenJetsNoCut(0)
  ,fhnRecJetsCut(0)
  ,fhnGenJetsCut(0)
  ,fhnRecBckg(0)
  ,fhnGenBckg(0)
  ,fhnRecJetsTrackUEcor(0)
  ,fhnGenJetsTrackUEcor(0)
  ,fhnRecJetsBckgUEcor(0)
  ,fhnGenJetsBckgUEcor(0)
  ,fhnTrackUE(0)
  ,fhnParticleUE(0)
  ,fhnBckgRecUE(0)
  ,fhnBckgGenUE(0)
  ,fRecJetBranch("")
  ,fGenJetBranch("")
  ,fRecBckgBranch("")
  ,fGenBckgBranch("")
  ,fRecJetR(0.)
  ,fGenJetR(0.)
  ,fRecBckgR(0.)
  ,fGenBckgR(0.)
  ,fTrackType(0)
  ,fParticleType(0)
  ,fhnMatching(0)
  ,fhnTrackCorrMatching(0)
  ,fhnBckgCorrMatching(0)
  ,fhnTrackUEanal(0)
  ,kDoUEanalysis(kFALSE)
  ,fRejectPileUp(1)
{
  if(fDebug) printf("%s: %d  Constructor\n",(char*)__FILE__, __LINE__);
  DefineOutput(1, TList::Class());
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::Init()
{
  if(fDebug) printf("%s: %d  Initialize\n", (char*)__FILE__, __LINE__);
  return;
}

//---------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskPPJetSpectra::Notify()
{

  if(fDebug) printf("%s: %d  Notify\n", (char*)__FILE__, __LINE__);
  return kTRUE;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::UserCreateOutputObjects()
{
  fOutputList = new TList;
  PostData(1,fOutputList);
  fOutputList->SetOwner(1);

  const Double_t binsPt[] = {0.,2.,4.,6.,8.,10.,12., 14.,16., 18.,20.,24.,28.,32.,38.,44.,50.,58.,66.,76.,86.,100.,114.,128.,150.,250.};
  const Int_t nBinsPt = 25;

  const Double_t trackPtBins[42] = {0,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,8.,10.,12.,14.,16.,18.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100};
  const Int_t nTrackPtBins = 41;

  if(fDebug) printf("%s: %d  Creating output objects\n",(char*)__FILE__,__LINE__);
  // triggered - centrality - z-vtxs - r-vtx - ncontrib - ntracks - IsPileUp?- njetsrec - njetsgen - accepted
  const Int_t evtsDim = 7;
  Int_t     evtsBins[evtsDim] = {   2,   22,   40,  20,  20,    2,  2};
  Double_t  evtsMins[evtsDim] = {-0.5,  -5., -20.,  0.,  0., -0.5, -0.5};
  Double_t  evtsMaxs[evtsDim] = { 1.5, 105.,  20.,  2., 40.,  1.5, 1.5};
  // pt - eta - phi - is272 - charge
  const Int_t trksDim = 7;
  Int_t     trksBins[trksDim] = {  nTrackPtBins,  10,             10,    2,   2,   2,    3};
  Double_t  trksMins[trksDim] = {            0.,  -1.,            0., -0.5,-0.5,-0.5, -1.5};
  Double_t  trksMaxs[trksDim] = {            0.,   1., 2*TMath::Pi(),  1.5, 1.5, 1.5,  1.5};
  // pt - eta - phi - AreaCh - AreaN - ptLeadingTrack/ptJet - nTracks
  const Int_t jetsDim = 7;
  Int_t     jetsBins[jetsDim] = { nBinsPt,  20,            10,3, 5, 20, 30};
  Double_t  jetsMins[jetsDim] = {      0., -1.,            0., 0., 0., 0., 0.};
  Double_t  jetsMaxs[jetsDim] = {      0.,  1., 2*TMath::Pi(), 3., 5., 1.,30.};
  // pt - eta - phi - AreaCH - AreaN - UEdensity - ptLeadingTrack/ptJet - ptOrig 
  const Int_t jetsUEcorDim = 8;
  Int_t     jetsUEcorBins[jetsUEcorDim] = { nBinsPt,  10,            10, 3,  5, 100,  20, nBinsPt};
  Double_t  jetsUEcorMins[jetsUEcorDim] = {  	 0., -1.,            0., 0., 0.,  0.,  0.,  0.};
  Double_t  jetsUEcorMaxs[jetsUEcorDim] = {  	 0.,  1., 2*TMath::Pi(), 3., 5., 10.,  1.,  0.};
  // leading jet pt - eta - sumAllPt1 - sumAllPt2 - Rpar
  const Int_t ueDim = 5;
  Int_t       ueBins[ueDim] = {  nBinsPt, 20, 100, 100, 100};
  Double_t    ueMins[ueDim] = {       0., -1,  0.,  0., 0.};
  Double_t    ueMaxs[ueDim] = {       0.,  1, 10., 10., 10.};
  // ptGen - ptRec - ptFraction - maxDist - bothWays - (1 - rec/gen) - dR - dEta - dPhi - isSamePtBin
  const Int_t matchingDim = 10;
  Int_t	      matchingBins[matchingDim] = {nBinsPt, nBinsPt,  20,    5,    2, 100, 10,  10, 10, 2  };
  Double_t    matchingMins[matchingDim] = {	0., 	 0.,  0., 0.05, -0.5,   0,  0,   0,  0,-0.5};
  Double_t    matchingMaxs[matchingDim] = {	0., 	 0.,  1., 0.55,  1.5,   1,0.5, 0.5,0.5, 1.5};

  fhnEvent    = new THnSparseF("fhnEvent", "fhnEvent;triggered;centrality;zVtx;rVtx;Ncontrib;IsPileUp;isAccepted", evtsDim, evtsBins, evtsMins, evtsMaxs);
  fhnTracks   = new THnSparseF("fhnTracks", "fhnTraks;pT;eta;phi;isFilterMask272;isMask16;isMask256;charge;", trksDim, trksBins, trksMins, trksMaxs);
  fhnMC       = new THnSparseF("fhnMC",     "fhnTraks;pT;eta;phi;isFilterMask272;isMask16;isMask256;charge;", trksDim, trksBins, trksMins, trksMaxs);
  fhnMC2      = new THnSparseF("fhnMC2",    "fhnTraks;pT;eta;phi;isFilterMask272;isMask16;isMask256;charge;", trksDim, trksBins, trksMins, trksMaxs);

  fhnRecJetsNoCut  = new THnSparseF("fhnRecJetsNoCut", "fhnRecJetsNoCut;pT;eta;phi;AreaCH;AreaN;z;nTracks;", jetsDim, jetsBins, jetsMins, jetsMaxs);
  fhnGenJetsNoCut  = new THnSparseF("fhnGenJetsNoCut", "fhnGenJetsNoCut;pT;eta;phi;AreaCH;AreaN;z;nTracks;", jetsDim, jetsBins, jetsMins, jetsMaxs);
  fhnRecJetsCut    = new THnSparseF("fhnRecJetsCut",     "fhnRecJetsCut;pT;eta;phi;AreaCH;AreaN;z;nTracks;", jetsDim, jetsBins, jetsMins, jetsMaxs);
  fhnGenJetsCut    = new THnSparseF("fhnGenJetsCut",     "fhnGenJetsCut;pT;eta;phi;AreaCH;AreaN;z;nTracks;", jetsDim, jetsBins, jetsMins, jetsMaxs);
  //UE sparses
  fhnTrackUE     = new THnSparseF("fhnTrackUE",      "fhnTrackUE;pT-leading jet;eta;sum1;sum2;avgSum;", ueDim, ueBins, ueMins, ueMaxs);
  fhnParticleUE  = new THnSparseF("fhnParticleUE","fhnParticleUE;pT-leading jet;eta;sum1;sum2;avgSum;", ueDim, ueBins, ueMins, ueMaxs);
  // pt - eta - phi - AreaCH - AreaN - UEdensity - ptLeadingTrack/ptJet - ptOrig 
  fhnRecJetsTrackUEcor  = new THnSparseF("fhnRecJetsTrackUEcor", "fhnRecJetsTrackUEcor;pT;eta;phi;AreaCH;AreaN;UE;z;pTorig;", jetsUEcorDim, jetsUEcorBins, jetsUEcorMins, jetsUEcorMaxs);
  fhnGenJetsTrackUEcor  = new THnSparseF("fhnGenJetsTrackUEcor", "fhnGenJetsTrackUEcor;pT;eta;phi;AreaCH;AreaN;UE;z;pTorig;", jetsUEcorDim, jetsUEcorBins, jetsUEcorMins, jetsUEcorMaxs);
  // ptGen - ptRec - ptFraction - maxDist - bothWays
  fhnMatching = new THnSparseF("fhnMatching","fhnMatching;pTgen;pTrec;ptFraction;maxDist;bothWays;1-p_{T,rec}/p_{T,gen}; dR; dEta; dPhi;;isSamePtBin;", matchingDim, matchingBins, matchingMins, matchingMaxs);
  fhnTrackCorrMatching = new THnSparseF("fhnTrackCorrMatching","fhnTrackCorrMatching;pTgen;pTrec;ptFraction;maxDist;bothWays;1-p_{T,rec}/p_{T,gen}; dR; dEta; dPhi;isSamePtBin;", matchingDim, matchingBins, matchingMins, matchingMaxs);

  fhnTracks->SetBinEdges(0, trackPtBins);
  fhnTrackUE->SetBinEdges(0, trackPtBins);
  fhnParticleUE->SetBinEdges(0, trackPtBins);

  fhnRecJetsNoCut->SetBinEdges(0, binsPt);
  fhnGenJetsNoCut->SetBinEdges(0, binsPt);
  fhnRecJetsCut->SetBinEdges(0, binsPt);
  fhnGenJetsCut->SetBinEdges(0, binsPt);

  fhnRecJetsTrackUEcor->SetBinEdges(0,binsPt);
  fhnRecJetsTrackUEcor->SetBinEdges(7,binsPt);
  fhnGenJetsTrackUEcor->SetBinEdges(0,binsPt);
  fhnGenJetsTrackUEcor->SetBinEdges(7,binsPt);

  fhnMatching->SetBinEdges(0,binsPt);
  fhnMatching->SetBinEdges(1,binsPt);
  fhnTrackCorrMatching->SetBinEdges(0,binsPt);
  fhnTrackCorrMatching->SetBinEdges(1,binsPt);

  fhnEvent->Sumw2();
  fhnTracks->Sumw2();
  fhnMC->Sumw2();
  fhnMC2->Sumw2();
  fhnRecJetsNoCut->Sumw2();
  fhnGenJetsNoCut->Sumw2();
  fhnRecJetsCut->Sumw2();
  fhnGenJetsCut->Sumw2();
  fhnRecJetsTrackUEcor->Sumw2();
  fhnGenJetsTrackUEcor->Sumw2();
  fhnTrackUE->Sumw2();
  fhnParticleUE->Sumw2();
  fhnMatching->Sumw2();
  fhnTrackCorrMatching->Sumw2();

  fOutputList->Add(fhnEvent);
  fOutputList->Add(fhnTracks);
  fOutputList->Add(fhnMC);
  fOutputList->Add(fhnMC2);
  fOutputList->Add(fhnRecJetsNoCut);
  fOutputList->Add(fhnGenJetsNoCut);
  fOutputList->Add(fhnRecJetsCut);
  fOutputList->Add(fhnGenJetsCut);
  fOutputList->Add(fhnRecJetsTrackUEcor);
  fOutputList->Add(fhnGenJetsTrackUEcor);
  fOutputList->Add(fhnTrackUE);
  fOutputList->Add(fhnParticleUE);
  fOutputList->Add(fhnMatching);
  fOutputList->Add(fhnTrackCorrMatching);

  const Int_t trackUEdim = 11;
  //max track pt - ch_part_dens - ch_pt_dens - avg_pt_track - zero particles -tMAX_part_dens-tMIN_part_dens-tDIF_part_dens-tMAX_pt_dens-tMIN_pt_dens-tDIF_pt_dens
  Int_t trackUEbins[trackUEdim] = { 100, 50, 50, 50,   2, 50, 50, 50, 50, 50, 50};
  Double_t trackUEmins[trackUEdim] = {  0., 0., 0.,  0.,-0.5, 0., 0., 0., 0., 0., 0.};
  Double_t trackUEmaxs[trackUEdim] = {100., 5., 5., 5., 1.5, 5., 5., 5., 5., 5., 5.};

  fhnTrackUEanal = new THnSparseF("fhnTrackUEanal","fhnTrackUEanal;PTmax;#rho^{ch.part}_{T};#rho^{pT}_{T};avg-pt/track;zero_particles;#rho^{ch.part}_{transMAX};#rho^{ch.part}_{transMIN};#rho^{ch.part}_{transDIF};#rho^{pT}_{transMAX};#rho^{pT}_{transMIN};#rho^{pT}_{transDIF}", trackUEdim, trackUEbins, trackUEmins, trackUEmaxs);
  fhnTrackUEanal->Sumw2();
  fOutputList->Add(fhnTrackUEanal);

}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::UserExec(Option_t */*option*/)
{
  Double_t evtContainer[7];
  Bool_t isEventSelected = EventSelection(evtContainer);

  TList trackList; 
  TList particleList;
  Int_t nTracks         = GetListOfTracks(fTrackType, &trackList);
  Int_t nParticles      = GetListOfTracks(fParticleType, &particleList);

 if(fDebug > 10 ) 
    std::cout<< nTracks<< " " << nParticles<<std::endl;  

  TClonesArray *tcaRecJets = 0;
  TClonesArray *tcaGenJets = 0;
  TClonesArray *tcaRecBckg = 0;
  TClonesArray *tcaGenBckg = 0;

  Int_t rJ, gJ, rB, gB;
  if(fAODOut && !tcaRecJets && fRecJetBranch.Length() > 0) tcaRecJets = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fRecJetBranch.Data()));
  if(fAODOut && !tcaGenJets && fGenJetBranch.Length() > 0) tcaGenJets = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fGenJetBranch.Data()));
  if(fAODOut && !tcaRecBckg && fRecBckgBranch.Length() > 0) tcaRecBckg = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fRecBckgBranch.Data()));
  if(fAODOut && !tcaGenBckg && fGenBckgBranch.Length() > 0) tcaGenBckg = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fGenBckgBranch.Data()));

  if(fAODExt) {
    if(!tcaRecJets && fRecJetBranch.Length() > 0) tcaRecJets = dynamic_cast<TClonesArray*>(fAODExt->GetAOD()->FindListObject(fRecJetBranch.Data()));
    if(!tcaGenJets && fGenJetBranch.Length() > 0) tcaGenJets = dynamic_cast<TClonesArray*>(fAODExt->GetAOD()->FindListObject(fGenJetBranch.Data()));
    if(!tcaRecBckg && fRecBckgBranch.Length() > 0) tcaRecBckg = dynamic_cast<TClonesArray*>(fAODExt->GetAOD()->FindListObject(fRecBckgBranch.Data()));
    if(!tcaGenBckg && fGenBckgBranch.Length() > 0) tcaGenBckg = dynamic_cast<TClonesArray*>(fAODExt->GetAOD()->FindListObject(fGenBckgBranch.Data()));
  }

  if(fAODIn) {
    if(!tcaRecJets && fRecJetBranch.Length() > 0) tcaRecJets = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fRecJetBranch.Data()));
    if(!tcaGenJets && fGenJetBranch.Length() > 0) tcaGenJets = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fGenJetBranch.Data()));
    if(!tcaRecBckg && fRecBckgBranch.Length() > 0) tcaRecBckg = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fRecBckgBranch.Data()));
    if(!tcaGenBckg && fGenBckgBranch.Length() > 0) tcaGenBckg = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fGenBckgBranch.Data()));
  }

  if(!tcaRecJets) rJ = -1; 
  else rJ = tcaRecJets->GetEntries();
  if(!tcaGenJets) gJ = -1; 
  else gJ = tcaGenJets->GetEntries();
  if(!tcaRecBckg) rB = -1; 
  else rB = tcaRecBckg->GetEntries();
  if(!tcaGenBckg) gB = -1; 
  else gB = tcaGenBckg->GetEntries();

  TList jetRecListNoCut; 
  TList jetGenListNoCut; 
  Int_t nRecJetsNoCut   = (rJ>0?GetListOfJets(tcaRecJets, &jetRecListNoCut, kFALSE):0);
  Int_t nGenJetsNoCut   = (gJ>0?GetListOfJets(tcaGenJets, &jetGenListNoCut, kFALSE):0);

  TList jetRecListCut;   
  TList jetGenListCut;   
  Int_t nRecJetsCut     = (rJ>0?GetListOfJets(tcaRecJets, &jetRecListCut, kTRUE):0);
  Int_t nGenJetsCut     = (gJ>0?GetListOfJets(tcaGenJets, &jetGenListCut, kTRUE):0);

  TList bckgRecListCut;  
  TList bckgGenListCut;  
  Int_t nRecBckgCut     = (rB>0?GetListOfJets(tcaRecBckg, &bckgRecListCut, kTRUE):0);
  Int_t nGenBckgCut     = (gB>0?GetListOfJets(tcaGenBckg, &bckgGenListCut, kTRUE):0);

  if(fDebug > 10) std::cout<<nRecJetsNoCut<<" "<<nGenJetsNoCut<<" "<<nRecJetsCut<< " "<< nGenJetsCut<<" "<<nRecBckgCut<<" "<<nGenBckgCut<<std::endl;

  evtContainer[6] = (isEventSelected?1:0);
  fhnEvent->Fill(evtContainer);

  if(!isEventSelected) {
    PostData(1,fOutputList);
    return;
  }

  if(kDoUEanalysis) DoUEAnalysis(&trackList, (Double_t)0.15, (Double_t)0.9);

  Double_t  recUEdensity = GetUE(&jetRecListCut, &trackList, fRecJetR, fhnTrackUE);
  Double_t  genUEdensity = GetUE(&jetGenListCut, &particleList, fGenJetR, fhnParticleUE);

  TList  trackUEcorrJetRecListCut;
  TList  trackUEcorrJetGenListCut;

  CorrectForUE(&jetRecListCut, recUEdensity, &trackUEcorrJetRecListCut, fhnRecJetsTrackUEcor);
  CorrectForUE(&jetGenListCut, genUEdensity, &trackUEcorrJetGenListCut, fhnGenJetsTrackUEcor);

  FillJetContainer(&jetRecListNoCut, fhnRecJetsNoCut);
  FillJetContainer(&jetGenListNoCut, fhnGenJetsNoCut);
  FillJetContainer(&jetRecListCut,   fhnRecJetsCut);
  FillJetContainer(&jetGenListCut,   fhnGenJetsCut);

    if(trackUEcorrJetRecListCut.GetEntries() != 0  && trackUEcorrJetGenListCut.GetEntries() != 0 ) {
    	MatchJets(kTRUE,  &trackUEcorrJetRecListCut,	&trackUEcorrJetGenListCut, 0.1,  fhnTrackCorrMatching);
    	MatchJets(kTRUE,  &trackUEcorrJetRecListCut,	&trackUEcorrJetGenListCut, 0.3,  fhnTrackCorrMatching);
    	MatchJets(kTRUE,  &trackUEcorrJetRecListCut,	&trackUEcorrJetGenListCut, 0.5,  fhnTrackCorrMatching);
    }
    if(jetRecListCut.GetEntries() != 0  && jetGenListCut.GetEntries() != 0 ) {
    	MatchJets(kTRUE,  &jetRecListCut,		&jetGenListCut,		  0.1,  fhnMatching);
    	MatchJets(kTRUE,  &jetRecListCut,		&jetGenListCut,		  0.3,  fhnMatching);
    	MatchJets(kTRUE,  &jetRecListCut,		&jetGenListCut,		  0.5,  fhnMatching);
    }


  PostData(1,fOutputList);
  return;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::Terminate(Option_t *)
{
  if(fDebug) printf("%s: %d Terminating\n",(char*)__FILE__, __LINE__);
  return;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::SetRecJetBranch(TString s) {
  fRecJetBranch = s;
  if(s.Contains("KT02") || s.Contains("SISCONE02") || s.Contains("UA102")) fRecJetR = 0.2;
  else if(s.Contains("KT03") || s.Contains("SISCONE03") || s.Contains("UA103")) fRecJetR = 0.3;
  else if(s.Contains("KT04") || s.Contains("SISCONE04") || s.Contains("UA104")) fRecJetR = 0.4;
  else if(s.Contains("KT05") || s.Contains("SISCONE05") || s.Contains("UA105")) fRecJetR = 0.5;
  else if(s.Contains("KT06") || s.Contains("SISCONE06") || s.Contains("UA106")) fRecJetR = 0.6;
  else fRecJetR = 0;
  return;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::SetGenJetBranch(TString s) {
  fGenJetBranch = s;
  if(s.Contains("KT02") || s.Contains("SISCONE02") || s.Contains("UA102")) fGenJetR = 0.2;
  else if(s.Contains("KT03") || s.Contains("SISCONE03") || s.Contains("UA103")) fGenJetR = 0.3;
  else if(s.Contains("KT04") || s.Contains("SISCONE04") || s.Contains("UA104")) fGenJetR = 0.4;
  else if(s.Contains("KT05") || s.Contains("SISCONE05") || s.Contains("UA105")) fGenJetR = 0.5;
  else if(s.Contains("KT06") || s.Contains("SISCONE06") || s.Contains("UA106")) fGenJetR = 0.6;
  else fGenJetR = 0;
  return;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::SetRecBckgBranch(TString s) {
  fRecBckgBranch = s;
  if(s.Contains("KT02") || s.Contains("SISCONE02") || s.Contains("UA102")) fRecBckgR = 0.2;
  else if(s.Contains("KT03") || s.Contains("SISCONE03") || s.Contains("UA103")) fRecBckgR = 0.3;
  else if(s.Contains("KT04") || s.Contains("SISCONE04") || s.Contains("UA104")) fRecBckgR = 0.4;
  else if(s.Contains("KT05") || s.Contains("SISCONE05") || s.Contains("UA105")) fRecBckgR = 0.5;
  else if(s.Contains("KT06") || s.Contains("SISCONE06") || s.Contains("UA106")) fRecBckgR = 0.6;
  else fRecBckgR = 0;
  return;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::SetGenBckgBranch(TString s) {
  fGenBckgBranch = s;
  if(s.Contains("KT02") || s.Contains("SISCONE02") || s.Contains("UA102")) fGenBckgR = 0.2;
  else if(s.Contains("KT03") || s.Contains("SISCONE03") || s.Contains("UA103")) fGenBckgR = 0.3;
  else if(s.Contains("KT04") || s.Contains("SISCONE04") || s.Contains("UA104")) fGenBckgR = 0.4;
  else if(s.Contains("KT05") || s.Contains("SISCONE05") || s.Contains("UA105")) fGenBckgR = 0.5;
  else if(s.Contains("KT06") || s.Contains("SISCONE06") || s.Contains("UA106")) fGenBckgR = 0.6;
  else fGenBckgR = 0;
  return;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::SetVertexCuts(Int_t nCont, Double_t Zcut, Double_t Rcut) {
  if(fDebug)
    printf("%s: %d Setting Vertex Cuts\n",(char*)__FILE__,__LINE__);
  nVtxContCut = nCont;
  fVtxZcut = Zcut;
  fVtxRcut = Rcut;
  return;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::SetTrackCuts(Double_t ptMin, Double_t ptMax, Double_t etaMax) {
  if(fDebug > 1) 
    printf("%s: %d Track cuts set\n",(char*)__FILE__,__LINE__);
  trackPtMin = ptMin;
  trackPtMax = ptMax;
  trackEtaAbsMax = etaMax;
  return;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::SetJetCuts(Double_t ptMin, Double_t eta, Double_t z) {
  jetPtMin = ptMin;
  jetEtaCut = eta;
  jetZmax = z;
  return;
}

//---------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskPPJetSpectra::EventSelection(Double_t evtContainer[6]) {
  if(fDebug > 1) 
    printf("%s: %d UserExec analysing event %d.\n", (char*)__FILE__, __LINE__, (Int_t)fEntry);

  // Trigger selection
  AliInputEventHandler* inputHandler = (AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(fDebug > 2 ) 
    printf("IsEventSelected: %d eventSelectionMask: %d\n", (Int_t)inputHandler->IsEventSelected(), (Int_t)fEvtSelectionMask);

  Bool_t isTriggerSelected = inputHandler->IsEventSelected() & fEvtSelectionMask;
  evtContainer[0] = (Int_t)isTriggerSelected;

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if(!fESD && fDebug > 2) 
    printf("%s: %d No ESD event found\n",(char*)__FILE__, __LINE__);

  if(!fESD) fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
  fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());

  if(!fESD)fAOD = fAODIn;
  else fAOD = fAODOut;

  if(fNonStdFile.Length()!=0) {
    AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
    fAODExt = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
  }

  TObject* handler = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  Float_t centrality = -1;
  if(fEventClass > 0)
  {
    if(handler->InheritsFrom("AliAODHandler")) centrality = ((AliVAODHeader*)fAODIn->GetHeader())->GetCentrality();
    else if(fESD) centrality = fESD->GetCentrality()->GetCentralityPercentile("V0M");
    else centrality = AliAnalysisHelperJetTasks::EventClass();
  }
  evtContainer[1] = centrality;

  // Event selection
  AliAODVertex* primVtx = fAOD->GetPrimaryVertex();
  if(!primVtx) {
    evtContainer[2]  =  1e10;
    evtContainer[3]  =  1e10;
    evtContainer[4]  =  1e10;
    return kFALSE;
  } else {  
    primVtx->Print();
    evtContainer[2] = primVtx->GetZ();
    evtContainer[3] = TMath::Sqrt(primVtx->GetX()*primVtx->GetX() + primVtx->GetY()*primVtx->GetY() );
    evtContainer[4] = primVtx->GetNContributors();
  }

  if(fRejectPileUp && AliAnalysisHelperJetTasks::IsPileUp()){
    if (fDebug > 1) Printf("%s:%d SPD pileup: event REJECTED...", (char*)__FILE__,__LINE__);
    evtContainer[5] = 0.;
    return kFALSE;
  }
  else 
    evtContainer[5] = 1.;

  if(evtContainer[0] == 0) 
  {
    if(fDebug > 2) 
      printf("%s: %d Event rejected: Trigger\n",(char*)__FILE__, __LINE__);
    return kFALSE;
  }
  if(evtContainer[1] >= 0 && evtContainer[1] < 10*((Int_t)fEventClass/10) && evtContainer[1] > 10*((Int_t)fEventClass/10 + 1) )
  {
    if(fDebug > 2) 
      printf("%s: %d Event rejected: Centrality\n",(char*)__FILE__, __LINE__);
    return kFALSE;
  }
  if(TMath::Abs(evtContainer[2]) > fVtxZcut)
  {
    if(fDebug > 2) 
      printf("%s: %d Event rejected: Vertex position (z)\n",(char*)__FILE__, __LINE__);
    return kFALSE;
  }
  if(evtContainer[3] > fVtxRcut)
  {
    if(fDebug > 2)
      printf("%s: %d Event rejected: Vertex position (x-y)\n",(char*)__FILE__,__LINE__);
    return kFALSE;
  }
  if(evtContainer[4] < nVtxContCut)
  {
    if(fDebug > 2) 
      printf("%s: %d Event rejected: Vertex contributors \n",(char*)__FILE__, __LINE__);
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskPPJetSpectra::GetListOfTracks(Int_t trackType, TList *trackList) {
  if(fDebug) 
    printf("%s: %d Getting Tracks\n",(char*)__FILE__,__LINE__);
  if(trackType == AliAnalysisTaskPPJetSpectra::kNone) {
    if(fDebug) printf("%s: %d No track type selected", (char*)__FILE__,__LINE__);
    return 0;
  }

  if(1 == trackType) {
    for(Int_t i = 0; i < fAOD->GetNumberOfTracks(); i++) {
      Bool_t isOK = kTRUE;
      Bool_t isOK16 = kTRUE;
      Bool_t isOK256 = kTRUE;
      AliAODTrack* track = (AliAODTrack*)fAOD->GetTrack(i);
      if(!track) continue;
      if(nTrackFilter>0 && !track->TestFilterBit(nTrackFilter) ) isOK = kFALSE;
      if(nTrackFilter>0 && !track->TestFilterBit(16) ) isOK16 = kFALSE;
      if(nTrackFilter>0 && !track->TestFilterBit(256) ) isOK256 = kFALSE;
      if(track->Pt() < trackPtMin || track->Pt() > trackPtMax) isOK = kFALSE;
      if(trackEtaAbsMax < TMath::Abs(track->Eta())) isOK = kFALSE;
      Double_t tmpContainer[] = {track->Pt(), track->Eta(), track->Phi(), (Double_t)isOK, (Double_t)isOK16, (Double_t)isOK256, (Double_t)track->Charge()/3};
      if(isOK) trackList->Add(track);
      fhnTracks->Fill(tmpContainer);
    }
  }
  else {
    TClonesArray *tca=dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca) {
      if(fDebug > 2) printf("no branch %s\n", AliAODMCParticle::StdBranchName());
    }
    for(Int_t i = 0; i < tca->GetEntries(); i++) {
      AliAODMCParticle *particle = (AliAODMCParticle*)tca->At(i);
      if(!particle) continue;
      if(!particle->IsPhysicalPrimary()) continue;
      Double_t charge = particle->Charge();
      if( charge > 0 ) charge = 1;
      if( charge < 0 ) charge = -1;
      if( (trackType == 3) && (charge == 0) ) continue;
      trackList->Add(particle);
      Double_t tmpContainer[] = {particle->Pt(), particle->Eta(), particle->Phi(), -999, -999, -999, charge };
      if(trackType == 2) fhnMC->Fill(tmpContainer);
      else fhnMC2->Fill(tmpContainer);
    }
  }
  trackList->Sort();
  return trackList->GetEntries();
}

//---------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskPPJetSpectra::GetListOfJets(TClonesArray* jca, TList* jetList, Bool_t useCuts) {
  if(fDebug) 
    printf("%s: %d Getting list of jets\n",(char*)__FILE__,__LINE__);

  Double_t ptCut = (useCuts   &&  jetPtMin>0 ?jetPtMin :0.1);
  Double_t etaCut = (useCuts  &&  jetEtaCut>0?jetEtaCut:0.9);
  Double_t zCut  = (useCuts   &&  jetZmax>0  ?jetZmax  :1);

  if(!jca){
    return -2;
  } 

  for(Int_t i =0 ; i < jca->GetEntries(); i++) {
    AliAODJet* jet = (AliAODJet*)jca->At(i);
    if(!jet) continue;
    Double_t zLeading = -1;
//    if(jet->GetPtLeading() <= 0) {
//      TRefArray* jetRef = (TRefArray*)jet->GetRefTracks();
//      for(int j = 0; j < jetRef->GetEntries(); j++) {
//	AliVTrack *track = (AliVTrack*)jet->GetTrack(j);
//        if(!track) continue;
//        cout<<track->Pt()<<" "<<zLeading<<endl;
//	if(track->Pt() > zLeading) zLeading = track->Pt();
//        cout<<__LINE__<<endl;
//      }
//      if(jet->Pt() > 0) zLeading /= jet->Pt();
//      else zLeading = -1;
//    }
    //cuty
    if(ptCut > 0  && ptCut > jet->Pt()                ) continue;
    if(etaCut > 0 && etaCut < TMath::Abs(jet->Eta())  ) continue;
    if(zCut > 0   && zCut < zLeading  && zLeading > 0 ) continue;

    jetList->Add(jet);
  }

  jetList->Sort();
  return jetList->GetEntries();
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::FillJetContainer(TList* jets, THnSparseF* container) {
  if(!jets || !container) return;

  for(Int_t i =0 ; i < jets->GetEntries(); i++) {
    AliAODJet* jet = (AliAODJet*)jets->At(i);
    if(!jet) continue;
  // pt - eta - phi - AreaCH - AreaN -  ptLeadingTrack/ptJet - ptOrig 
    Double_t tmpContainer[] = {jet->Pt(), jet->Eta(), jet->Phi(), jet->EffectiveAreaCharged(), jet->EffectiveAreaNeutral(), jet->GetPtLeading()/jet->Pt(), (Double_t)((TRefArray*)jet->GetRefTracks())->GetEntries() };
    container->Fill(tmpContainer);
  }
  return;
}

//---------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskPPJetSpectra::GetUE(TList* jets, TList* tracks, Double_t Rpar,THnSparseF* thnCont) {
  if(!jets || !tracks) {
    if(fDebug) printf("%s: %d No jets or tracks\n",(char*)__FILE__,__LINE__);
    return 0;
  }

  //find leading jet
  AliAODJet *jet = 0;
  for(Int_t iTmp = 0; iTmp < jets->GetEntries(); iTmp++) {
    AliAODJet* tmp = (AliAODJet*)jets->At(iTmp);
    if(!tmp) continue;
    if(!jet) jet = tmp;
    else 
      if( jet->Pt() < tmp->Pt() ) jet = tmp;
  }
  if(!jet) return 0;
  
  //get track in perpendicular cone
  //TVector3 for jet and perpendicular cones with same eta axis
  TVector3 j(jet->Px(), jet->Py(), jet->Pz());
  TVector3 p1(j);
  TVector3 p2(j);

  p1.RotateZ(TMath::Pi()/2.);
  p2.RotateZ(-TMath::Pi()/2.);

  Double_t sumAllPt1 = 0;
  Double_t sumAllPt2 = 0;

  for(int i = 0; i < tracks->GetEntries(); i++)  {
    AliVParticle* tr = (AliVParticle*)tracks->At(i);
    if(!tr) continue;
    TVector3 v(tr->Px(), tr->Py(), tr->Pz());
    Double_t dR1 = v.DrEtaPhi(p1);
    Double_t dR2 = v.DrEtaPhi(p2);

    if(dR1 < Rpar) sumAllPt1+=v.Pt();
    if(dR2 < Rpar) sumAllPt2+=v.Pt();
  }

  if(Rpar != 0) {
    sumAllPt1/=(TMath::Pi()*Rpar*Rpar);
    sumAllPt2/=(TMath::Pi()*Rpar*Rpar);
  }

  Double_t container[] = {j.Pt(), j.Eta(), sumAllPt1, sumAllPt2, (sumAllPt1+sumAllPt2)/2.};
  thnCont->Fill(container);

  return (sumAllPt1+sumAllPt2)/2.;

}

//---------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskPPJetSpectra::GetBckgUE(TList* jets, Double_t Rpar, Bool_t isSkipped, THnSparseF* cont) {
  if(!jets) return 0;

  Int_t leading01 = -1;
  Int_t leading02 = -1;
  Double_t leading01pt = 0;
  Double_t leading02pt = 0;

  if(!isSkipped) {
    for(Int_t i = 0; i < jets->GetEntries(); i++) {
      AliAODJet* tmp = (AliAODJet*)jets->At(i);
      if(!tmp) continue;
      if(tmp->Pt() > leading01pt) {
        leading02 = leading01;
        leading01 = i;
        leading02pt = leading01pt;
        leading01pt = tmp->Pt();
      } else if(tmp->Pt() > leading02) {
        leading02 = i;
        leading02pt = tmp->Pt();
      }
    }
  }

  Int_t subLeading01 = -1;
  Int_t subLeading02 = -1;
  Double_t subLeading01pt = 0;
  Double_t subLeading02pt = 0;
  Double_t subLeading01area = 0;
  Double_t subLeading02area = 0;

  Double_t avgAll = 0;
  Double_t areaAll = 0;
  //search for subleading jets assuming leading 2 (skip02) jets are skipped
  if(jets->GetEntries() > 0) {
    for(Int_t i = 0; i < jets->GetEntries(); i++) {
      if(i == leading01 || i == leading02) continue;
      AliAODJet* tmp = (AliAODJet*)jets->At(i);
      if(!tmp) continue;
      if( subLeading01pt < tmp->Pt() ) {
        subLeading02 = subLeading01;
        subLeading01 = i;
        subLeading02pt = subLeading01pt;
        subLeading01pt = tmp->Pt();
        subLeading02area = subLeading01area;
        subLeading01area = tmp->EffectiveAreaCharged();
      } else if( subLeading02pt < tmp->Pt() ) {
        subLeading02 = i;
        subLeading02pt = tmp->Pt();
        subLeading02area = tmp->EffectiveAreaCharged();
      }
    }

    //do avg jet pt
    Int_t nJets = 0;
    for(int i = 0; i < jets->GetEntries();i++) {
      if(i == leading01 || i == leading02) continue;
      AliAODJet* tmp = (AliAODJet*)jets->At(i);
      if(!tmp) continue;
      nJets++;
      avgAll+= tmp->Pt();
      areaAll += tmp->EffectiveAreaCharged();
    }
  }


  if(avgAll!=0) avgAll/=areaAll;
  if(subLeading01pt!=0) subLeading01pt/=subLeading01area;
  if(subLeading02pt!=0)subLeading02pt/=subLeading02area;
  Double_t tmpContainer[]={subLeading01pt,subLeading02pt,avgAll,Rpar*10};

  cont->Fill(tmpContainer);

  return avgAll; 
}

//---------------------------------------------------------------------------------------------------
Int_t  AliAnalysisTaskPPJetSpectra::CorrectForUE(TList* jets, Double_t UE, TList* newList, THnSparseF* container) {
  if(!jets) return 0;

  for(Int_t i = 0; i < jets->GetEntries(); i++) {
    AliAODJet* tmpJet = (AliAODJet*)jets->At(i);
    if(!tmpJet) continue;

    AliAODJet* jet = (AliAODJet*)tmpJet->Clone(Form("%s%s",tmpJet->GetName(),"_clone"));
    Double_t corrUE = UE*(jet->EffectiveAreaCharged() + jet->EffectiveAreaNeutral() );

    if(jet->Pt() > corrUE) jet->SetPtEtaPhiM(jet->Pt()- corrUE, jet->Eta(), jet->Phi(), jet->M());
    else jet->SetPtEtaPhiM(0, jet->Eta(), jet->Phi(), jet->M());
    newList->Add(jet);

  // pt - eta - phi - AreaCH - AreaN - UEdensity - ptLeadingTrack/ptJet - ptOrig 
    Double_t cont[8]={jet->Pt(), jet->Eta(), jet->Phi(), jet->EffectiveAreaCharged(),jet->EffectiveAreaNeutral(), corrUE, (jet->Pt() > 0?jet->GetPtLeading() / jet->Pt():-1),tmpJet->Pt()};
    container->Fill(cont);
  }
  newList->Sort();
  return newList->GetEntries();
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::MatchJets(Bool_t doBothWay, TList* recList, TList* genList, Float_t maxDist, THnSparseF* container) {
  if(!genList || !recList) return;
  Int_t mode = 1;

  const Int_t kGenJets = genList->GetEntries();
  const Int_t kRecJets = recList->GetEntries();

  TArrayI iMatchIndex(kGenJets);
  TArrayF fPtFraction(kGenJets);

  TArrayI aGenIndex(kRecJets);
  TArrayI aRecIndex(kGenJets);

  Double_t genPt = -1;
  Double_t recPt = -1;
  Double_t ptFraction = -1;

  //-- Closest jets to generated and checked for E-fraction
  if(!doBothWay) {
    AliAnalysisHelperJetTasks::GetJetMatching(genList, kGenJets, recList, kRecJets, iMatchIndex, fPtFraction, fDebug, maxDist, mode);
    for(Int_t ig = 0; ig < kGenJets; ig++) {
      //get gen-jet
      AliAODJet* genJet = (AliAODJet*)genList->At(ig);
      if(!genJet) continue;
      genPt = genJet->Pt();
      
      //get rec-jet
      Int_t ir = iMatchIndex[ig];
      if(ir<0||ir>=recList->GetEntries()) continue;
      AliAODJet* recJet = (AliAODJet*)recList->At(ir);
      if(!recJet) continue;
      recPt = recJet->Pt();

      ptFraction = fPtFraction[ig];

      Double_t dR = recJet->DeltaR(genJet);
      Double_t dEta = recJet->Eta() - genJet->Eta();
      Double_t dPhi2 = dR*dR - dEta*dEta;
      
	Bool_t isSame = kFALSE;
	if(CheckPtBin(genPt) == CheckPtBin(recPt)) isSame = kTRUE;
      // ptGen - ptRec - ptFraction - maxDist - bothWays - (1 - rec/gen) - dR - dEta  - dPhi
      Double_t cont[] = {genPt,recPt,ptFraction,maxDist,(Float_t)doBothWay, (genPt != 0 ? (1. - recPt/genPt):1), dR, dEta, (dPhi2>=0?TMath::Sqrt(dPhi2):-1), (Double_t)isSame};
      container->Fill(cont);
    }
  }
  //-- Closest rec jet to gen and vice-versa & check for E-fraction
  else {
    AliAnalysisHelperJetTasks::GetClosestJets(genList, kGenJets, recList, kRecJets, aGenIndex, aRecIndex,fDebug, maxDist);
    for(Int_t ig = 0; ig < kGenJets; ig++) {
      //get gen-jet
      AliAODJet* genJet = (AliAODJet*)genList->At(ig);
      if(!genJet) continue;
      genPt = genJet->Pt();

      //get rec-jet
      Int_t ir = aRecIndex[ig];
      if(ir>=0&&ir<recList->GetEntries()) {
	AliAODJet *recJet = (AliAODJet*)recList->At(ir);
	if(!recJet) continue;
	recPt = recJet->Pt();
	Double_t dR = recJet->DeltaR(genJet);
	Double_t dEta = recJet->Eta() - genJet->Eta();
	Double_t dPhi2 = dR*dR - dEta*dEta;

	ptFraction = AliAnalysisHelperJetTasks::GetFractionOfJet(recJet,genJet,mode);

	Bool_t isSame = kFALSE;
	if(CheckPtBin(genPt) == CheckPtBin(recPt)) isSame = kTRUE;
	// ptGen - ptRec - ptFraction - maxDist - bothWays - (1 - rec/gen) - dR - dEta  - dPhi
	Double_t cont[] = {genPt,recPt,ptFraction,maxDist,(Float_t)doBothWay, (genPt != 0 ? (1. - recPt/genPt):1), dR, dEta, (dPhi2>=0?TMath::Sqrt(dPhi2):-1), (Double_t)isSame  };
        container->Fill(cont);
      }
    }
  }

  return;
}

//---------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskPPJetSpectra::CheckPtBin(Double_t pT) {
  Double_t bins[] 	= {0.,2.,4.,6.,8.,10.,12., 14.,16., 18.,20.,24.,28.,32.,38.,44.,50.,58.,66.,76.,86.,100.,114.,128.,150.,250.};
  Double_t nBins 	= 25;

  for(Int_t i = 0; i <= nBins; i++) {
    if(pT < bins[i] ) return i;
  }
  return nBins+1;
}

//---------------------------------------------------------------------------------------------------
void AliAnalysisTaskPPJetSpectra::DoUEAnalysis(TList* tracks, Double_t PTcut, Double_t ETAcut) {
	Double_t PTmax = 0;
	int iMax = -1;
	for(int i = 0; i < tracks->GetEntries(); i++) {
		AliVParticle *part = (AliVParticle*)tracks->At(i);
		if(!part) continue;

		if(part->Pt() < PTcut) continue;
		if(TMath::Abs(part->Eta()) > ETAcut ) continue;

		if(PTmax < part->Pt()) {
			PTmax = part->Pt();
			iMax = i;
		}
	}

	if(PTmax < trackPtMin) return;

	Double_t ch_part_density = 0;
	Double_t ch_part_pt_density = 0;
	Double_t avg_track_pt = 0;
	Double_t zero_ch_part = 0;
	Double_t transMAX_part_density = 0;
	Double_t transMIN_part_density = 0;
	Double_t transMAX_pt_density = 0;
	Double_t transMIN_pt_density = 0;
	Double_t transDIF_pt_density = 0;
	Double_t transDIF_part_density = 0;

	Double_t A1_part_density = 0;
	Double_t A2_part_density = 0;
	Double_t A1_pt_density = 0;
	Double_t A2_pt_density = 0;

	AliVParticle *maxParticle = (AliVParticle*)tracks->At(iMax);

	Double_t PHImax = maxParticle->Phi();

	for(int i = 0; i < tracks->GetEntries(); i++) {
		AliVParticle* part = (AliVParticle*)tracks->At(i);
		if(!part) continue;
		if(part->Pt() < PTcut) continue;

		Bool_t isArea1 = kFALSE;
		Bool_t isArea2 = kFALSE;

		Double_t dPhi = part->Phi() - PHImax;
		if(dPhi > TMath::Pi()) dPhi -= TMath::Pi();
		if(dPhi < -TMath::Pi()) dPhi +=TMath::Pi();

		if(TMath::Abs(dPhi) > TMath::Pi()/3 && TMath::Abs(dPhi) < 2*TMath::Pi()/3 ) {
			if(dPhi>0) isArea1 = kTRUE;
			else isArea2 = kTRUE;
		}

		if(!isArea1 && !isArea2) continue;

		ch_part_density += 1;
		ch_part_pt_density += part->Pt();
		avg_track_pt += part->Pt();

		if(isArea1) {
			A1_part_density += 1;
			A1_pt_density += part->Pt();
		} else if(isArea2) {
			A2_part_density += 1;
			A2_pt_density += part->Pt();
		}
	}


	if(ch_part_density > 0) avg_track_pt/=ch_part_density;
	else zero_ch_part = 1;
	ch_part_density*=3./(4.*TMath::Pi()*ETAcut);
	ch_part_pt_density*= 3./(4.*TMath::Pi()*ETAcut);

	A1_part_density*=3./(4.*TMath::Pi()*ETAcut);
	A1_pt_density*= 3./(4.*TMath::Pi()*ETAcut);
	A2_part_density*=3./(4.*TMath::Pi()*ETAcut);
	A2_pt_density*= 3./(4.*TMath::Pi()*ETAcut);

	if(A1_part_density < A2_part_density) {
		transMAX_part_density = A2_part_density;
		transMIN_part_density = A1_part_density;
	} else {
		transMAX_part_density = A1_part_density;
		transMIN_part_density = A2_part_density;
	}

	if(A1_pt_density < A2_pt_density) {
		transMAX_pt_density = A2_pt_density;
		transMIN_pt_density = A1_pt_density;
	} else {
		transMAX_pt_density = A1_pt_density;
		transMIN_pt_density = A2_pt_density;
	}

	transDIF_pt_density 	= transMAX_pt_density - transMIN_pt_density;
	transDIF_part_density 	= transMAX_part_density - transMIN_part_density;

	Double_t container[] = {PTmax, ch_part_density, ch_part_pt_density, avg_track_pt, zero_ch_part,
							transMAX_part_density, transMIN_part_density, transDIF_part_density,
							transMAX_pt_density, transMIN_pt_density, transDIF_pt_density};

	fhnTrackUEanal->Fill(container);
}
