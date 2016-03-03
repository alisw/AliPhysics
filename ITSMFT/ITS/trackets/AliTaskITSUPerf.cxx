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
#include "../ITS/UPGRADE/AliITSMFTAux.h"
#include <TGeoManager.h>
#include <TTree.h>


ClassImp(AliTaskITSUPerf)


const char* AliTaskITSUPerf::fgkLabelTypes[AliTaskITSUPerf::kNLabelTypes] = {
  "ITSokTPCok"
  ,"ITSokTPCfk"
  ,"ITSfkTPCok"
  ,"ITSfkTPCfk"
  ,"ITSTPCmismatch"
  ,"ITSTPCnoMatch"
};


typedef struct {
  Bool_t  prim;  
  Bool_t  rcbl;
  Char_t  nClITS;
  Char_t  nClTPC;
  Char_t  nClITSMC;
  Char_t  mcCl[7];
  Char_t  rcCl[7];
  Char_t  fcCl[7];
  Char_t  charge;
  Float_t zV;
  Float_t ptMC;
  Float_t etaMC;
  Float_t pt;
  Float_t eta;
  Float_t ptTPC;
  Float_t etaTPC;
  Float_t dcaR;
  Float_t dcaZ;
  Float_t dcaRE;
  Float_t dcaZE;
  Float_t ptE;
  Float_t ptTPCE;
  Float_t phi;
  Int_t   pdg;
  Int_t   lbl;
  Int_t   lblTPC;
  Int_t   spl;
  Int_t   evID;
} trInfo_t;

trInfo_t trackInfo;

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
    //
#ifdef _CLUS_LIST_
  ,fClLists(0)
#endif
    //
  ,fNPtBins(20)
  ,fNResBins(100)
  ,fMinTPCclusters(60)
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
  ,fTree(0)
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
  fTree = new TTree("trinfo","trinfo");
  fTree->Branch("prim",&trackInfo.prim,"prim/O");
  fTree->Branch("rcbl",&trackInfo.rcbl,"rcbl/O");
  fTree->Branch("nClITS",&trackInfo.nClITS,"nClITS/b");
  fTree->Branch("nClTPC",&trackInfo.nClTPC,"nClTPC/b");
  fTree->Branch("nClITSMC",&trackInfo.nClITSMC,"nClITSMC/b");
  fTree->Branch("mcCl",&trackInfo.mcCl,"mcCl[7]/b");
  fTree->Branch("rcCl",&trackInfo.rcCl,"rcCl[7]/b");
  fTree->Branch("fcCl",&trackInfo.fcCl,"fcCl[7]/b");
  fTree->Branch("charge",&trackInfo.charge,"charge/B");
  fTree->Branch("zV",    &trackInfo.zV,"zV/F");
  fTree->Branch("ptMC",  &trackInfo.ptMC,"ptMC/F");
  fTree->Branch("etaMC", &trackInfo.etaMC,"etaMC/F");
  fTree->Branch("pt",  &trackInfo.pt,"pt/F");
  fTree->Branch("eta", &trackInfo.eta,"eta/F");
  fTree->Branch("ptTPC",  &trackInfo.ptTPC,"ptTPC/F");
  fTree->Branch("etaTPC", &trackInfo.etaTPC,"etaTPC/F");
  fTree->Branch("dcaR", &trackInfo.dcaR,"dcaR/F");
  fTree->Branch("dcaZ", &trackInfo.dcaZ,"dcaZ/F");
  fTree->Branch("dcaRE", &trackInfo.dcaRE,"dcaRE/F");
  fTree->Branch("dcaZE", &trackInfo.dcaZE,"dcaZE/F");
  fTree->Branch("ptE",  &trackInfo.ptE,"ptE/F");
  fTree->Branch("ptTPCE",  &trackInfo.ptTPCE,"ptTPCE/F");
  fTree->Branch("phi",  &trackInfo.phi,"phi/F");
  fTree->Branch("pdg",  &trackInfo.pdg,"pdg/I");
  fTree->Branch("lbl",  &trackInfo.lbl,"lbl/I");
  fTree->Branch("lblTPC",  &trackInfo.lblTPC,"lblTPC/I");
  fTree->Branch("spl",  &trackInfo.spl,"spl/I");
  fTree->Branch("evID",  &trackInfo.evID,"evID/I");
  fOutput->Add(fTree);
  //
  fTPCCut = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  fTPCCut->SetEtaRange(fEtaMin,fEtaMax);
  fTPCCut->SetPtRange(fPtMin,fPtMax);
  fTPCCut->SetPtRange(fPtMin,fPtMax);
  fTPCCut->SetMinNClustersTPC(fMinTPCclusters);
  printf("Cuts: Eta: %f %f  | Pt: %f %f\n",fEtaMin,fEtaMax,fPtMin,fPtMax);
  //
#ifdef _CLUS_LIST_
  fClLists.SetOwner();
#endif
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
  if (!fRPTree) { AliWarning("Invalid ITS cluster tree"); }
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
#ifdef _CLUS_LIST_
  fClLists.Delete();
#endif
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
  h2 = new TH2F(Form("B%d_%s_RCBL",bin,"MCLabvsPT"),
		Form("B%d %s Reconstructable",bin,"MCLabComb vs p_{T}"),
		fNPtBins,fPtMin,fPtMax,kNLabelTypes,-0.5,kNLabelTypes-0.5);
  AddHisto(&fHistosCent,h2, GetHistoID(kHMatchStatusRcbl,-1,bin) );
  for (int i=0;i<kNLabelTypes;i++) h2->GetYaxis()->SetBinLabel(i+1,fgkLabelTypes[i]);
  //  
  h2 = new TH2F(Form("B%d_%s_NotRCBL_prim",bin,"MCLabvsPT"),
		Form("B%d %s Not Reconstructable, Prim.",bin,"MCLabComb vs p_{T}"),
		fNPtBins,fPtMin,fPtMax,kNLabelTypes,-0.5,kNLabelTypes-0.5);
  AddHisto(&fHistosCent,h2, GetHistoID(kHMatchStatusNRcblPrim,-1,bin) );
  for (int i=0;i<kNLabelTypes;i++) h2->GetYaxis()->SetBinLabel(i+1,fgkLabelTypes[i]);
  //
  h2 = new TH2F(Form("B%d_%s_NotRCBL_sec",bin,"MCLabvsPT"),
		Form("B%d %s Not Reconstructable, Sec.",bin,"MCLabComb vs p_{T}"),
		fNPtBins,fPtMin,fPtMax,kNLabelTypes,-0.5,kNLabelTypes-0.5);
  AddHisto(&fHistosCent,h2, GetHistoID(kHMatchStatusNRcblSec,-1,bin) );
  for (int i=0;i<kNLabelTypes;i++) h2->GetYaxis()->SetBinLabel(i+1,fgkLabelTypes[i]);
  //
  h2 = new TH2F(Form("B%d_MCLr_Prim",bin),
		Form("B%d MCLr Prim vs pt",bin),
		fNPtBins,fPtMin,fPtMax, 7,-0.5,6.5);
  AddHisto(&fHistosCent,h2, GetHistoID(kHMCLrPresPrim,-1,bin) );
  //
  h2 = new TH2F(Form("B%d_MCLr_Sec",bin),
		Form("B%d MCLr Sec vs pt",bin),
		fNPtBins,fPtMin,fPtMax, 7,-0.5,6.5);
  AddHisto(&fHistosCent,h2, GetHistoID(kHMCLrPresSec,-1,bin) );
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
  //  if (mcStat>=kNLabelTypes) AliFatal(Form("MCStatus ID=%d > max.allowed %d",mcStat,kNLabelTypes));
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
  Bool_t dump = kFALSE;
  Bool_t dump1 = kFALSE;
  static int evid = 0;
  trackInfo.evID = evid++;
  trackInfo.zV   = fESDEvent->GetPrimaryVertex()->GetZ();
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
    //
    trc->GetDZ(fVtxMC[0],fVtxMC[1],fVtxMC[2], fldBz, dcaRZ);
    //
    Int_t labMC    = trc->GetLabel();
    Int_t labMCabs = TMath::Abs(labMC);
    Int_t labMCTPC = trc->GetTPCLabel();
    Int_t labMCITS = trc->GetITSLabel();
    Int_t nClTPC   = trc->GetNcls(1);
    Int_t nClITS   = trc->GetNcls(0);
    //
    UInt_t& mcStatus = (UInt_t &)fMCStatus[labMCabs];
    double pt  = trc->Pt();
    double eta = trc->Eta();
    //
    TH2* hlr = (TH2*)( (mcStatus&BIT(kMCPrimBit)) ? GetHisto(&fHistosCent,kHMCLrPresPrim) : GetHisto(&fHistosCent,kHMCLrPresSec));
    if (hlr) {
      for (int il=0;il<7;il++) {
	if (mcStatus&(0x1<<(il+kITSHitBits))) hlr->Fill(pt,il);
      }      
      hlr->Fill(pt,-99); // count tracks
    }
    // 

    Bool_t reject = mcStatus & BIT(kTrCondFail);
    //
    int mcLabType = GetMCLabType(labMCTPC,labMCITS,nClTPC,nClITS);
    //
    //    if (eta>fEtaMax || eta<fEtaMin) continue;
    //    if (nClTPC<fMinTPCclusters) continue;
    //    printf("#%3d pt:%5.2f eta:%+5.2f | Lb:%+6d LbTPC:%+6d LbITS:%+6d MCTp:%d | Ntpc:%3d Nits:%3d | Rej:%d\n",
    //	   itr,pt,eta,labMC,labMCTPC,labMCITS,mcLabType, nClTPC,nClITS,reject);
    AliMCParticle *part  = 0;
    if (labMCabs>=ntrMC) continue;
    part = (AliMCParticle*)fMCEvent->GetTrack(labMCabs);
    double ptMC = part->Pt();
    double etaMC = part->Eta();
    //
    if (fTree) {
      trackInfo.lbl = labMC;
      trackInfo.lblTPC = labMCTPC;
      trackInfo.prim = (mcStatus&BIT(kMCPrimBit)) != 0;
      trackInfo.rcbl = !reject;
      trackInfo.nClITS = trc->GetNcls(0);
      trackInfo.nClTPC = trc->GetNcls(1);
      trackInfo.ptMC = ptMC;
      trackInfo.etaMC = etaMC;
      trackInfo.charge = trc->Charge();
      trackInfo.pt = trc->Pt();
      trackInfo.eta = trc->Eta();
      //
      
      //
      trackInfo.spl = trc->GetITSModuleIndex(10);
      trackInfo.phi = part->Phi();    
      trackInfo.pdg = part->PdgCode();
      trackInfo.dcaR = dcaRZ[0];
      trackInfo.dcaZ = dcaRZ[1];
      trackInfo.dcaRE = TMath::Sqrt(trc->GetSigmaY2());
      trackInfo.dcaZE = TMath::Sqrt(trc->GetSigmaZ2());
      trackInfo.ptE   = TMath::Sqrt(trc->GetSigma1Pt2());

      const AliExternalTrackParam* ttpc = trc->GetTPCInnerParam();
      if (ttpc) {
	trackInfo.ptTPC = ttpc->Pt();
	trackInfo.etaTPC = ttpc->Eta();
	trackInfo.ptTPCE = TMath::Sqrt(ttpc->GetSigma1Pt2());//*ttpc->Pt()*ttpc->Pt();
      }
      else {
	trackInfo.ptTPC = -1;
	trackInfo.etaTPC = 0;
	trackInfo.ptTPCE = -1;
      }
      //
      trackInfo.nClITSMC = 0;
      for (int il=0;il<7;il++) {
	trackInfo.mcCl[il] = (mcStatus & (0x1<<(il+kITSHitBits))) != 0;
	if (trackInfo.mcCl[il]) trackInfo.nClITSMC++;
	trackInfo.rcCl[il] = trc->HasPointOnITSLayer(il);
	trackInfo.fcCl[il] = trc->HasSharedPointOnITSLayer(il);
      }
      fTree->Fill();
    }
    if (reject) {
      GetHisto(&fHistosCent,(mcStatus & BIT(kMCPrimBit)) ? kHMatchStatusNRcblPrim:kHMatchStatusNRcblSec)->Fill(pt,mcLabType);  // matching status
      if (dump) {
	printf("--#%4d Pt:%5.2f Eta:%5.2f |",itr,pt,eta);
	for (int k=0;k<32;k++) printf("%d", (mcStatus&(0x1<<k)) ? 1:0);
	printf("| %+5d %+5d %+5d |%3d %d -> %d\n",labMC,labMCTPC,labMCITS,nClTPC,nClITS,mcLabType);
      }
      continue;
    }
    GetHisto(&fHistosCent,kHMatchStatusRcbl)->Fill(pt,mcLabType);  // matching status
    //
    TObjArray* clList = 0;
    //
    if (part) {
      /*
      int pdg = part->PdgCode();
      if (TMath::Abs(pdg)!=211) continue;
      if (nClITS<7) continue;
      */
#ifdef _CLUS_LIST_
      //
      clList = (TObjArray*)fClLists.At(labMCabs);
      if (clList) {
	printf("Clusters for track %d MCLabel:%d\n",itr,labMC);
	for (int ic=0;ic<clList->GetEntriesFast();ic++) {
	  AliITSUClusterPix* cl = (AliITSUClusterPix*) clList->At(ic);
	  cl->Print();
	}
      }
#endif
    }

    if (part) {
      //
      //
      if (dump) {
	printf("  #%4d Pt:%5.2f Eta:%5.2f |",itr,ptMC,etaMC);
	for (int k=0;k<32;k++) printf("%d", (mcStatus&(0x1<<k)) ? 1:0);
	printf("| %+5d %+5d %+5d |%3d %d -> %d\n",labMC,labMCTPC,labMCITS,nClTPC,nClITS,mcLabType);
      }
      //
      if ( (mcStatus&BIT(kMCPrimBit)) ) {
	if (dump1 && TMath::Abs(ptMC-0.5)<0.1) {
	  printf("#%4d Pt:%5.2f Eta:%5.2f | ",itr,ptMC,etaMC);
	  for (int k=0;k<7;k++) printf("(%d/%d)", (mcStatus&(0x1<<k)) ? 1:0, trc->HasPointOnITSLayer(k));
	  printf("| %+5d %+5d %+5d |%3d %d -> %d\n",labMC,labMCTPC,labMCITS,nClTPC,nClITS,mcLabType);
	}
	// compare MC vs reco track params
	if (ptMC>0) GetHisto(&fHistosCentMCLb,kHResPTvsPTMC, mcLabType)->Fill(ptMC, (ptMC-pt)/ptMC);
	//
	// DCA resolution
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
#ifdef _CLUS_LIST_
  fClLists.Delete();
  if (fClLists.GetSize()<ntrMC) fClLists.Expand(ntrMC);
  //
  // create empty lists for reconstructed tracks
  int ntr = fESDEvent->GetNumberOfTracks();
  for (int itr=0;itr<ntr;itr++) {
    AliESDtrack* trc = fESDEvent->GetTrack(itr);
    int lb = TMath::Abs(trc->GetLabel());
    if (lb<ntrMC && !fClLists.At(itr)) {
      fClLists.AddAt(new TObjArray(7),lb);
    }
  }
#endif
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
    fITS->ProcessClusters(); // order in the same way as in reconstruction
    //
    for (int ilr=0;ilr<fITS->GetNLayersActive();ilr++) {
      AliITSURecoLayer* lr = fITS->GetLayerActive(ilr);
      for (int icl=lr->GetNClusters();icl--;) {
	AliCluster* cl = lr->GetCluster(icl); 
	// check labels
	for (int ilb=0;ilb<3;ilb++) {
	  Int_t lab = cl->GetLabel(ilb);
	  if (lab<0||lab>=ntrMC) continue;
	  UInt_t& mcStatus = (UInt_t &)fMCStatus[lab];
	  mcStatus |= 0x1<<(ilr+kITSHitBits);
#ifdef _CLUS_LIST_
	  TObjArray *clList = (TObjArray*)fClLists.At(lab);
	  if (clList) clList->AddLast(cl);  // create cl.list only for reconstructed tracks
#endif
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
	//
	/*
	int nlr = AliITSMFTAux::NumberOfBitsSet(patt);	
	AliMCParticle *part  = (AliMCParticle*)fMCEvent->GetTrack(itr);
	double etaMC = part->Eta();
	if (etaMC>=fEtaMin && etaMC<=fEtaMax && nlr) {
	  TH2* hlr = (TH2*)( (mcStatus&BIT(kMCPrimBit)) ? GetHisto(&fHistosCent,kHMCLrPresPrim) : GetHisto(&fHistosCent,kHMCLrPresSec));
	  if (hlr) {
	    double ptMC = part->Pt();
	    for (int il=0;il<7;il++) {
	      if (mcStatus&(0x1<<(il+kITSHitBits))) hlr->Fill(ptMC,il);
	    }      
	    hlr->Fill(ptMC,-99); // count tracks
	  }
	}
	*/
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

