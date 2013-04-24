/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////////
// Class AliTaskITSUPerf                                                 //
// Analysis task to produce data and MC histos needed for tracklets      //
// dNdEta extraction in multiple bins in one go                          //
// Author:  ruben.shahoyan@cern.ch                                       //
///////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h" 
#include "TList.h"
#include "TObjArray.h"
#include "TGeoGlobalMagField.h"

#include "AliAnalysisManager.h"

#include "AliESDEvent.h"  
#include "AliESDInputHandler.h"
#include "AliESDInputHandlerRP.h"
#include "AliCDBPath.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"
#include "AliMagF.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "AliPhysicsSelection.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliESDtrackCuts.h"

#include "AliTaskITSUPerf.h"
#include "../ITS/UPGRADE/AliITSURecoDet.h"
#include "../ITS/UPGRADE/AliITSUGeomTGeo.h"
#include "../ITS/UPGRADE/AliITSUClusterPix.h"
#include "../ITS/UPGRADE/AliITSUTrackCond.h"
#include <TGeoManager.h>


ClassImp(AliTaskITSUPerf)


const char* AliTaskITSUPerf::fgkLabelTypes[AliTaskITSUPerf::kNLabelTypes] = {
  "ITSokTPCok"
  ,"ITSokTPCfk"
  ,"ITSfkTPCok"
  ,"ITSfkTPCfk"
  ,"ITSTPCmismatch"
  ,"ITSTPCnoMatch"
};

//________________________________________________________________________
/*//Default constructor
AliTaskITSUPerf::AliTaskITSUPerf(const char *name)
  : AliAnalysisTaskSE(name),
*/  
//________________________________________________________________________
AliTaskITSUPerf::AliTaskITSUPerf(const char *name) 
  : AliAnalysisTaskSE(name)
    //
  ,fOutput(0)
  ,fHistosCentMCLb(0)
  ,fHistosCent(0)
  ,fRPTree(0)
  ,fStack(0)
  ,fMCEvent(0)
  ,fVtxSPD(0)
  ,fVtxTrc(0)
  ,fESDEvent(0)
  ,fUseSpecialOutput(kTRUE)
  ,fUseMC(kTRUE)
  ,fTrackingCond(0)
    //
  ,fGeom(0)
  ,fITS(0)
  ,fNSelTracksMC(0)
  ,fMCStatus(0)
  ,fNPtBins(20)
  ,fNResBins(100)
  ,fMinTPCclusters(70)
  ,fPtMin(0)
  ,fPtMax(10.)
  ,fEtaMin(-3.)
  ,fEtaMax( 3.)
  ,fZVertexMin(-15.)
  ,fZVertexMax(15.)
    //
  ,fCurrCentBin(0)
  ,fNCentBins(1)
  ,fUseCentralityVar(0)
  ,fTPCCut(0)
{
  // Constructor

  DefineOutput(1, TList::Class());
  //
}

//________________________________________________________________________
AliTaskITSUPerf::~AliTaskITSUPerf()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {  //RRR
    printf("Deleteing output\n");
    delete fOutput;
    fOutput = 0;
  }
  //
  fHistosCentMCLb.Clear();
  fHistosCent.Clear();
  //
  delete fITS;
  delete fGeom;
  delete fTPCCut;
  //
}

//________________________________________________________________________
void AliTaskITSUPerf::UserCreateOutputObjects() 
{
  //
  if (fUseSpecialOutput) OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner(); 
  //
  // init geometry
  if (!gGeoManager) {
    TString geom = "geometry.root";
    AliGeomManager::LoadGeometry(geom.Data());
    if (!gGeoManager) AliFatal("Failed to load geometry");
  }
  //
  // create ITS interface
  fGeom = new AliITSUGeomTGeo(kTRUE,kTRUE);
  AliITSUClusterPix::SetGeom(fGeom);
  fITS = new AliITSURecoDet(fGeom,"ITSURecoInterface");
  fITS->CreateClusterArrays();
  //
  // Create histograms
  for (int ib=0;ib<fNCentBins;ib++) {
    BookHistos(ib);
  }
  //
  int nhist = fOutput->GetEntries();
  for (int i=0;i<nhist;i++) {
    TObject* hst = fOutput->At(i);
    if (!hst || !(hst->InheritsFrom(TH1::Class()))) continue;
    ((TH1*)hst)->Sumw2();
  }
  //
  fTPCCut = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fTPCCut->SetEtaRange(fEtaMin,fEtaMax);
  fTPCCut->SetPtRange(fPtMin,fPtMax);
  fTPCCut->SetPtRange(fPtMin,fPtMax);
  fTPCCut->SetMinNClustersTPC(fMinTPCclusters);
  //
  PostData(1, fOutput);
  //
}

//________________________________________________________________________
void AliTaskITSUPerf::UserExec(Option_t *) 
{
  // Main loop
  static int evID = 0;
  AliInfo(Form("Event: %d",evID++));
  //
  AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
  fRPTree = 0;
  AliESDInputHandlerRP *handRP = (AliESDInputHandlerRP*)anMan->GetInputEventHandler();
  if (!handRP) { AliFatal("No RP handler"); return; }
  //
  fESDEvent  = handRP->GetEvent();
  if (!fESDEvent) { AliFatal("No AliESDEvent"); return; }
  //
  fRPTree = handRP->GetTreeR("ITS");
  if (!fRPTree) { AliFatal("Invalid ITS cluster tree"); return; }
  //
  AliMCEventHandler* eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
  if (!eventHandler) { AliFatal("Could not retrieve MC event handler"); return; }
  fMCEvent = eventHandler->MCEvent();
  if (!fMCEvent) { AliFatal("Could not retrieve MC event"); return; }
  fStack = fMCEvent->Stack();
  if (!fStack) { AliFatal("Stack is not available"); return; }
  //
  // MC Generator info
  AliGenEventHeader* mcGenH = 0;
  mcGenH = fMCEvent->GenEventHeader();
  TArrayF vtmc(3);
  mcGenH->PrimaryVertex(vtmc);
  for (int i=3;i--;) fVtxMC[i] = vtmc[i];
  //
  printf("MagField: %f\n",fESDEvent->GetMagneticField());
  printf("VtxMC : %+e %+e %+e\n",fVtxMC[0],fVtxMC[1],fVtxMC[2]);
  //
  fVtxSPD = fESDEvent->GetPrimaryVertexSPD();
  //
  if (fVtxSPD->GetStatus()) printf("VtxSPD: %+e %+e %+e | Nc=%d\n",fVtxSPD->GetX(),fVtxSPD->GetY(),fVtxSPD->GetZ(),fVtxSPD->GetNContributors());
  else printf("VtxSPD: N/A\n");
  //
  fVtxTrc = fESDEvent->GetPrimaryVertexTracks();
  if (fVtxTrc->GetStatus()) printf("VtxTRC: %+e %+e %+e | Nc=%d\n",fVtxTrc->GetX(),fVtxTrc->GetY(),fVtxTrc->GetZ(),fVtxTrc->GetNContributors());
  else printf("VtxTRC: N/A\n");
  //
  BuildMCInfo();
  CheckTracks();
  //
}      


//________________________________________________________________________
void AliTaskITSUPerf::Terminate(Option_t *) 
{
  Printf("Terminating...");
}

//_________________________________________________________________________
void AliTaskITSUPerf::BookStandardHistosCentMCLb(Int_t bin, Int_t mcLb)
{
  // book standard set of histos for specific centrality bin and MC status type, 
  // adding them to output list and management array
  //
  const double kMaxDPt = 0.2;
  const double kMaxDCA = 0.1;

  TH2* h2=0;
  // pt resolution
  h2 = new TH2F(Form("B%d_%s_%s",bin,fgkLabelTypes[mcLb],"dPT2PTvsPT"),
		Form("B%d %s %s",bin,fgkLabelTypes[mcLb],"#Deltap_{T}/p_{T} vs p_{T}"),
		fNPtBins,fPtMin,fPtMax,fNResBins,-kMaxDPt,kMaxDPt);
  AddHisto(&fHistosCentMCLb,h2, GetHistoID(kHResPTvsPTMC,mcLb,bin) );
  //
  // transverse DCA resolution
  h2 = new TH2F(Form("B%d_%s_%s",bin,fgkLabelTypes[mcLb],"DCARvsPT"),
		Form("B%d %s %s",bin,fgkLabelTypes[mcLb],"DCA R vs p_{T}"),
		fNPtBins,fPtMin,fPtMax,fNResBins,-kMaxDCA,kMaxDCA);
  AddHisto(&fHistosCentMCLb,h2, GetHistoID(kHResDCARvsPTMC,mcLb,bin) );
  //
  // Z DCA resolution
  h2 = new TH2F(Form("B%d_%s_%s",bin,fgkLabelTypes[mcLb],"DCAZvsPT"),
		Form("B%d %s %s",bin,fgkLabelTypes[mcLb],"DCA Z vs p_{T}"),
		fNPtBins,fPtMin,fPtMax,fNResBins,-kMaxDCA,kMaxDCA);
  AddHisto(&fHistosCentMCLb,h2, GetHistoID(kHResDCAZvsPTMC,mcLb,bin) );
  //
}

//_________________________________________________________________________
void AliTaskITSUPerf::BookStandardHistosCent(Int_t bin)
{
  // book standard set of histos for specific centrality bin
  //
  TH2* h2=0;
  //
  // MC labels combination vs pt
  h2 = new TH2F(Form("B%d_%s",bin,"MCLabvsPT"),
		Form("B%d %s",bin,"MCLabComb vs p_{T}"),
		fNPtBins,fPtMin,fPtMax,kNLabelTypes,-0.5,kNLabelTypes-0.5);
  AddHisto(&fHistosCent,h2, GetHistoID(kHMatchStatus,-1,bin) );
  for (int i=0;i<kNLabelTypes;i++) h2->GetYaxis()->SetBinLabel(i+1,fgkLabelTypes[i]);
  //  

}

//_________________________________________________________________________
void AliTaskITSUPerf::BookHistos(Int_t bin)
{
  // book standard set of histos, adding them to output list and management array
  //
  // standard histos
  for (int ilb=0;ilb<kNLabelTypes;ilb++) {  // create similar set of histos for all kind of MC labels
    BookStandardHistosCentMCLb(bin,ilb);
  }
  BookStandardHistosCent(bin);              // st.histos regardless of MClabels comb.
}

//_________________________________________________________________________
void AliTaskITSUPerf::AddHisto(TObjArray* array,TObject* h, Int_t at)
{
  // add single histo to the set
  if (at>=0) array->AddAtAndExpand(h,at);
  else       array->Add(h);
  fOutput->Add(h);
}

//_________________________________________________________________________
Int_t AliTaskITSUPerf::GetHistoID(Int_t htype, Int_t mcStat, Int_t centBin) const
{
  // generate ID for the histo. htype must be <100, mcStat<kNLabelTypes
  //
  if (mcStat>=kNLabelTypes) AliFatal(Form("MCStatus ID=%d > max.allowed %d",mcStat,kNLabelTypes));
  if (centBin>=fNCentBins)  AliFatal(Form("CentBin  ID=%d > max.allowed %d",centBin,fNCentBins));
  //
  if (centBin>=0) {
    if (mcStat>=0) { // MCtruthness dependent histos
      if (htype>=kHNStdHistosCentMC) AliFatal(Form("StdCentMC histo ID=%d > max.allowed %d",htype,kHNStdHistosCentMC));
      return htype + kHNStdHistosCentMC*(mcStat + kNLabelTypes*centBin);
    }
    else {
      if (htype>=kHNStdHistosCent) AliFatal(Form("StdCent histo ID=%d > max.allowed %d",htype,kHNStdHistosCent));
      return htype + kHNStdHistosCent*centBin;
    }
  }
  // custom histo
  return htype;
  //
}

//_________________________________________________________________________
void AliTaskITSUPerf::CheckTracks()
{
  // build mc truth info for tracks 
  //
  int ntr = fESDEvent->GetNumberOfTracks();
  int ntrMC = fStack->GetNtrack();
  //
  double fldBz = fESDEvent->GetMagneticField();
  float dcaRZ[2];
  //
  for (int itr=0;itr<ntr;itr++) {
    AliESDtrack* trc = fESDEvent->GetTrack(itr);
    //
    // at the moment we consider only TPC/ITS tracks
    //if (!trc->IsOn(AliESDtrack::kTPCin)) continue;
    if (!fTPCCut->IsSelected(trc)) continue;

    Int_t labMC    = trc->GetLabel();
    Int_t labMCabs = TMath::Abs(labMC);
    Int_t labMCTPC = trc->GetTPCLabel();
    Int_t labMCITS = trc->GetITSLabel();
    Int_t nClTPC   = trc->GetNcls(1);
    Int_t nClITS   = trc->GetNcls(0);
    //
    UInt_t& mcStatus = (UInt_t &)fMCStatus[labMCabs];
    Bool_t reject = mcStatus & BIT(kTrCondFail);
    if (reject) continue;
    //
    int mcLabType = GetMCLabType(labMCTPC,labMCITS,nClTPC,nClITS);
    //
    double pt  = trc->Pt();
    double eta = trc->Eta();
    //
    //    if (eta>fEtaMax || eta<fEtaMin) continue;
    //    if (nClTPC<fMinTPCclusters) continue;
    //    printf("#%3d pt:%5.2f eta:%+5.2f | Lb:%+6d LbTPC:%+6d LbITS:%+6d MCTp:%d | Ntpc:%3d Nits:%3d\n",itr,pt,eta,labMC,labMCTPC,labMCITS,mcLabType, nClTPC,nClITS);
    //
    GetHisto(&fHistosCent,kHMatchStatus)->Fill(pt,mcLabType);  // matching status
    //
    AliMCParticle *part  = 0;
    //
    if (labMCabs<ntrMC) part = (AliMCParticle*)fMCEvent->GetTrack(labMCabs);
    //
    if (part) {
      //
      double ptMC = part->Pt();
      double etaMC = part->Eta();
      //
      if ( (mcStatus&BIT(kMCPrimBit)) ) {
	printf("#%4d Pt:%5.2f Eta:%5.2f |",itr,ptMC,etaMC);
	for (int k=0;k<32;k++) printf("%d", (mcStatus&(0x1<<k)) ? 1:0);
	printf("| %+5d %+5d %+5d |%3d %d -> %d\n",labMC,labMCTPC,labMCITS,nClTPC,nClITS,mcLabType);

	// compare MC vs reco track params
	if (ptMC>0) GetHisto(&fHistosCentMCLb,kHResPTvsPTMC, mcLabType)->Fill(ptMC, (ptMC-pt)/ptMC);
	//
	// DCA resolution
	trc->GetDZ(fVtxMC[0],fVtxMC[1],fVtxMC[2], fldBz, dcaRZ);
	GetHisto(&fHistosCentMCLb,kHResDCARvsPTMC, mcLabType)->Fill(ptMC, dcaRZ[0]);
	GetHisto(&fHistosCentMCLb,kHResDCAZvsPTMC, mcLabType)->Fill(ptMC, dcaRZ[1]);
	//
      }
      //
    }
    //


  }
  //

}

//_________________________________________________________________________
void AliTaskITSUPerf::BuildMCInfo()
{
  // build mc truth info for tracks 
  //
  int ntrMC = fStack->GetNtrack();
  //
  if (fMCStatus.GetSize()<ntrMC) fMCStatus.Set(ntrMC);
  fMCStatus.Reset();
  //
  for (int itr=0;itr<ntrMC;itr++) {
    //
    AliMCParticle *part  = (AliMCParticle*)fMCEvent->GetTrack(itr);
    if (!part->Charge()) continue;
    //
    UInt_t& mcStatus = (UInt_t &)fMCStatus[itr];
    mcStatus = 0;
    // 
    if (fStack->IsPhysicalPrimary(itr)) mcStatus |= BIT(kMCPrimBit);
    //
    /*
    TParticle* tp = part->Particle();
    printf("mc%4d IsPri=%d",itr,fStack->IsPhysicalPrimary(itr)); tp->Print();
    */
    //
  }
  //
  int nTrCond = fTrackingCond ? fTrackingCond->GetEntriesFast() : 0;
  // fill MC clusters info (if available)
  if (fRPTree) {
    // load clusters
    fITS->LoadClusters(fRPTree);
    //
    for (int ilr=fITS->GetNLayersActive();ilr--;) {
      AliITSURecoLayer* lr = fITS->GetLayerActive(ilr);
      for (int icl=lr->GetNClusters();icl--;) {
	AliCluster* cl = lr->GetCluster(icl); 
	// check labels
	for (int ilb=0;ilb<3;ilb++) {
	  Int_t lab = cl->GetLabel(ilb);
	  if (lab<0||lab>=ntrMC) continue;
	  UInt_t& mcStatus = (UInt_t &)fMCStatus[lab];
	  mcStatus |= 0x1<<(ilr+kITSHitBits);
	}
      }
    }
    //
    if (nTrCond) {
      for (int itr=0;itr<ntrMC;itr++) {
	//
	UInt_t& mcStatus = (UInt_t &)fMCStatus[itr];
	Bool_t fail = kTRUE;
	UShort_t patt = (mcStatus&0xffff);
	for (int ic=0;ic<nTrCond;ic++) {
	  AliITSUTrackCond* cnd = (AliITSUTrackCond*) fTrackingCond->At(ic);
	  if (cnd->CheckPattern(patt)) {fail=kFALSE; break;}
	}
	if (fail) mcStatus |= BIT(kTrCondFail);
      }
    }
    //
  }
  //
}


//_________________________________________________________________________
Int_t AliTaskITSUPerf::GetMCLabType(Int_t labMCTPC,Int_t labMCITS, Int_t nClTPC, Int_t nClITS)
{
  // build type of track (0: both ITS and TPC correct, 1: ITS corr. TPC fake, 2: ITS fake TPC corr., 3: ITS fake TPC fake, 4: mismatch
  //
  int tp = 0;
  if (nClTPC>0 && nClITS<1) return kITSTPCNoMatch;
  if (TMath::Abs(labMCTPC)!=TMath::Abs(labMCITS)) tp = kITSTPCMismatch;
  else {
    if   (labMCITS<0) tp = labMCTPC<0 ? kITS0TPC0 : kITS0TPC1;
    else              tp = labMCTPC<0 ? kITS1TPC0 : kITS1TPC1;
  }
  return tp;
  //
}

//_________________________________________________________________________
Int_t AliTaskITSUPerf::GetCentralityBin() const
{
  // calculate centrality bin
  return 0;
}

