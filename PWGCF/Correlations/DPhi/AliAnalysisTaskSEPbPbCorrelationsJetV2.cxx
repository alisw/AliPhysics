#include "AliLog.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TAxis.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TROOT.h"

#include "AliAnalysisTaskSEPbPbCorrelationsJetV2.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"

#include "AliCodeTimer.h"
#include "AliOADBContainer.h"

#include "AliEventPoolManager.h"
#include "AliTHn.h"


ClassImp(AliAnalysisTaskSEPbPbCorrelationsJetV2)
ClassImp(AliBasicParticleST)


 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::cont = NULL;
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2am[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2am[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2as[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2as[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2cm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2cm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2cs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2cs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2trm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2trm[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQx2trs[14] = { NULL };
 AliOADBContainer *AliAnalysisTaskSEPbPbCorrelationsJetV2::contQy2trs[14] = { NULL };


//====================================================================================================================================================

AliAnalysisTaskSEPbPbCorrelationsJetV2::AliAnalysisTaskSEPbPbCorrelationsJetV2() : 
  AliAnalysisTaskSE(), 
  fAOD(0x0),
  fPoolMgr(0x0),
  ftrk(0x0),
  fNbinsCent(1), 
  fNbinsPt(1),
  fNbinsPtTrig(1),
  fNbinsAssocPt(1),
  fNbinsZvtx(1),
  fCentAxis(0x0), 
  fPtAxis(0x0),
  fPtTrigAxis(0x0),
  fPtAssocAxis(0x0),
  fEtaAxis(0x0),
  fZvtxAxis(0x0),
  fRemovePileup(kFALSE),
  fRemovePileup2(kFALSE),
  fRemovePileup3(kFALSE),
  fPtOrder(kTRUE),
  fSameSign(kTRUE),
  fUseRes(kTRUE),
  fForceCL1(kFALSE),
  fN1(0),
  fN2(-1),

  fHistACv2(0x0),
  fHistATv2(0x0),
  fHistCTv2(0x0),
  fHistQxT2(0x0),
  fHistQyT2(0x0),
  fHistQxT3(0x0),
  fHistQyT3(0x0),

  fHistdPhidEtaPt(0x0), 
  fHistdPhidEtaPt_Mixed(0x0), 
  fHistdPhidEtaPt_SS(0x0), 
  fHistdPhidEtaPt_Mixed_SS(0x0), 

  fHistV0Multiplicity(0x0), 
  fHistITSMultiplicity(0x0),
  fHistCentrality(0x0),
  fHistEvStat(0x0),
  fHistNtrVsCent(0x0),
  fHistV0CVsCent(0x0),
  fHistESDvsTPC(0x0),

  fHistCentVsZ(0x0),
  fHistCentVsZMixed(0x0),
  
  fCentMethod("V0M"),
  fMode("ref_V0A"),
  flist_Res(0x0),
  flist_contQ(0x0),
  fOutputList(0x0),
  fOutputList1(0x0),
  fRunN(-1),
  fv0mult(-1.),
  fv0multonline(-1.),
  fitsmult(-1.),
  nITSTrkls(-1),
  fv0mpercentile(-1.),
  fv0meqpercentile(-1.),
  fv0apercentile(-1.),
  fv0cpercentile(-1.),
  fcl0percentile(-1.),
  fcl1percentile(-1.),
  fSPDpercentile(-1.),
  fzvtx(0.),
  multtrkcut(-1),
  Qx2trkcut(0.),
  Qy2trkcut(0.),

  sumMa(-1.),
  sumMc(-1.),
  Qxa2(0.),
  Qya2(0.),
  Qxc2(0.),
  Qyc2(0.),
   
  Qya2Cor(0.),
  Qxa2Cor(0.),
  Qyc2Cor(0.),
  Qxc2Cor(0.),
  
  Qytr2Cor(0.),
  Qxtr2Cor(0.),

  fResACv2(0x0),
  fResATv2(0x0),
  fResCTv2(0x0),
  
  fLowCenCut(0x0),
  fHighCenCut(0x0),
  fV0MultOfOnCut(0x0),
  fFilterBit(32),
  fTPCNcls(70),
  fMultV0(0x0)
  
{

  // Default constructor
  for (Int_t iCent=0; iCent<fNMaxBinsCentrality; iCent++) {
    for(Int_t iZvtx = 0; iZvtx<fNMaxBinsZvtx; ++iZvtx) {
      for(Int_t index = 0; index < 3; ++index) {
      fHistSP2A[iCent][iZvtx][index] = NULL;
      
     }
      for(Int_t ipt = 0; ipt < fNMaxBinsPt; ++ipt) {
	for(Int_t jpt = 0; jpt < fNMaxBinsAssocPt; ++jpt) {
	  fHistSP2AdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2AdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2CdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2TdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	}
      }
    }
  }

  for(Int_t i = 0; i < 14; ++i) {
    fQx2mV0A[i] = fQy2mV0A[i] = fQx2sV0A[i] = fQy2sV0A[i] = NULL;
    fQx2mV0C[i] = fQy2mV0C[i] = fQx2sV0C[i] = fQy2sV0C[i] = NULL;
    fQx2mTrk[i] = fQy2mTrk[i] = fQx2sTrk[i] = fQy2sTrk[i] = NULL;
  }

 //flist_Res = new TList();
 //flist_contQ = new TList();  
 
}


//====================================================================================================================================================

AliAnalysisTaskSEPbPbCorrelationsJetV2::AliAnalysisTaskSEPbPbCorrelationsJetV2(const char *name) : 
  AliAnalysisTaskSE(name), 
  fAOD(0x0),
  fPoolMgr(0x0),
  ftrk(0x0),
  fNbinsCent(1), 
  fNbinsPt(1),
  fNbinsPtTrig(1),
  fNbinsAssocPt(1),
  fNbinsZvtx(1),
  fCentAxis(0x0), 
  fPtAxis(0x0),
  fPtTrigAxis(0x0),
  fPtAssocAxis(0x0),
  fEtaAxis(0x0),
  fZvtxAxis(0x0),
  fRemovePileup(kFALSE),
  fRemovePileup2(kFALSE),
  fRemovePileup3(kFALSE),
  fPtOrder(kTRUE),
  fSameSign(kTRUE),
  fUseRes(kTRUE),
  fForceCL1(kFALSE),
  fN1(0),
  fN2(-1),

  fHistACv2(0x0),
  fHistATv2(0x0),
  fHistCTv2(0x0),
  fHistQxT2(0x0),
  fHistQyT2(0x0),
  fHistQxT3(0x0),
  fHistQyT3(0x0),
 
  fHistdPhidEtaPt(0x0),
  fHistdPhidEtaPt_Mixed(0x0),
  fHistdPhidEtaPt_SS(0x0),
  fHistdPhidEtaPt_Mixed_SS(0x0),
 
  fHistV0Multiplicity(0x0), 
  fHistITSMultiplicity(0x0),
  fHistCentrality(0x0),
  fHistEvStat(0x0),
  fHistNtrVsCent(0x0),
  fHistV0CVsCent(0x0),
  fHistESDvsTPC(0x0),

  fHistCentVsZ(0x0),
  fHistCentVsZMixed(0x0),
  
  fCentMethod("V0M"),
  fMode("ref_V0A"),
  flist_Res(0x0),
  flist_contQ(0x0),
  fOutputList(0x0),
  fOutputList1(0x0),
  fRunN(-1),
  fv0mult(-1.),
  fv0multonline(-1.),
  fitsmult(-1.),
  nITSTrkls(-1),
  fv0mpercentile(-1.),
  fv0meqpercentile(-1.),
  fv0apercentile(-1.),
  fv0cpercentile(-1.),
  fcl0percentile(-1.),
  fcl1percentile(-1.),
  fSPDpercentile(-1.),
  fzvtx(0.),
  multtrkcut(-1),
  Qx2trkcut(0.),
  Qy2trkcut(0.),

  sumMa(-1.),
  sumMc(-1.),
  Qxa2(0.),
  Qya2(0.),
  Qxc2(0.),
  Qyc2(0.),
   
  Qya2Cor(0.),
  Qxa2Cor(0.),
  Qyc2Cor(0.),
  Qxc2Cor(0.),
  
  Qytr2Cor(0.),
  Qxtr2Cor(0.),

  fResACv2(0x0),
  fResATv2(0x0),
  fResCTv2(0x0),
  
  fLowCenCut(0x0),
  fHighCenCut(0x0),
  fV0MultOfOnCut(0x0),
  fFilterBit(32),
  fTPCNcls(70),
  
  fMultV0(0x0)

{

  // Constructor
  for (Int_t iCent=0; iCent<fNMaxBinsCentrality; iCent++) {
    for(Int_t iZvtx = 0; iZvtx<fNMaxBinsZvtx; ++iZvtx) {
      for(Int_t index = 0; index < 3; ++index) {
      fHistSP2A[iCent][iZvtx][index] = NULL;
     }
      for(Int_t ipt = 0; ipt < fNMaxBinsPt; ++ipt) {
	for(Int_t jpt = 0; jpt < fNMaxBinsAssocPt; ++jpt) {
	  fHistSP2AdPhidEta[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2AdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2CdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	  fHistSP2TdPhidEtaSS[iCent][iZvtx][ipt][jpt] = NULL;
	}
      }
    }
  }

  for(Int_t i = 0; i < 14; ++i) {
    fQx2mV0A[i] = fQy2mV0A[i] = fQx2sV0A[i] = fQy2sV0A[i] = NULL;
    fQx2mV0C[i] = fQy2mV0C[i] = fQx2sV0C[i] = fQy2sV0C[i] = NULL;
    fQx2mTrk[i] = fQy2mTrk[i] = fQx2sTrk[i] = fQy2sTrk[i] = NULL;
  }
  //flist_Res = new TList();
  //flist_contQ = new TList();
  
  // Define input and output slots here
  DefineInput(1, TList::Class());
  DefineInput(2, TList::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//====================================================================================================================================================

AliAnalysisTaskSEPbPbCorrelationsJetV2::~AliAnalysisTaskSEPbPbCorrelationsJetV2() {
  
  if (fOutputList  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) 
    delete fOutputList;

  if (fOutputList1  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    delete fOutputList1;

  if (flist_contQ  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    delete flist_contQ;

  if (flist_Res  && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
    delete flist_Res;

}

//___________________________________________________________________________
void AliAnalysisTaskSEPbPbCorrelationsJetV2::NotifyRun()
{
  /// Notify run
  printf("Loading calibration for run %d\n",fInputHandler->GetEvent()->GetRunNumber());
  OpenInfoCalbration(fInputHandler->GetEvent()->GetRunNumber());
}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::UserCreateOutputObjects() {
  
  //=====================================================================
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fOutputList1 = new TList();
  fOutputList1->SetOwner(kTRUE);

  fOutputList->Add(new TH2D("fHist_Qx2_V0A","Q_{x2}(V0A) vs CL1",100,0.,100.,3000,-1500,1500));
  fOutputList->Add(new TH2D("Corr_fHist_Qx2_V0A","Corr Q_{x2}(V0A) vs CL1",100,0.,100.,100,-50,50));
 
  fOutputList->Add(new TH2D("fHist_Qy2_V0A","Q_{y2}(V0A) vs CL1",100,0.,100.,3000,-1500,1500));
  fOutputList->Add(new TH2D("Corr_fHist_Qy2_V0A","Corr Q_{y2}(V0A) vs CL1",100,0.,100.,100,-50,50));
 
  fOutputList->Add(new TH2D("fHist_Qx2_V0C","Q_{x2}(V0C) vs CL1",100,0.,100.,3000,-1500,1500));
  fOutputList->Add(new TH2D("Corr_fHist_Qx2_V0C","Corr Q_{x2}(V0C) vs CL1",100,0.,100.,100,-50,50));

  fOutputList->Add(new TH2D("fHist_Qy2_V0C","Q_{y2}(V0C) vs CL1",100,0.,100.,3000,-1500,1500));
  fOutputList->Add(new TH2D("Corr_fHist_Qy2_V0C","Corr Q_{y2}(V0C) vs CL1",100,0.,100.,100,-50,50));


    for (Int_t iCent=0; iCent<fNbinsCent; iCent++) {
      for(Int_t iZvtx = 0; iZvtx < fNbinsZvtx; ++iZvtx) {
      for(Int_t index = 0; index < 3; ++index) {
	fHistSP2A[iCent][iZvtx][index] = new TProfile(Form("fHist%sSP2A_Cent%02d_Z%02d_%d",fUseRes?"Res":"",iCent,iZvtx,index), 
						      Form("%d-%d %%, %.1f<z<%.1f",
							   Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
							   Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
							   fZvtxAxis->GetBinLowEdge(iZvtx+1),
							   fZvtxAxis->GetBinUpEdge(iZvtx+1)),
						      fNbinsPt,fPtAxis->GetXbins()->GetArray());
	fOutputList -> Add(fHistSP2A[iCent][iZvtx][index]);
      }
      for(Int_t iPt = 0; iPt < fNbinsPtTrig; ++iPt) {
	for(Int_t jPt = 0; jPt < fNbinsAssocPt; ++jPt) {
	  fHistSP2AdPhidEta[iCent][iZvtx][iPt][jPt] = new TProfile2D(Form("fHist%sSP2AdPhidEta_Cent%02d_Z%02d_Pt%02d_%02d",fUseRes?"Res":"",iCent,iZvtx,iPt,jPt), 
								Form("%d-%d %%, %.1f<z<%.1f",
								     Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								     Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								     fZvtxAxis->GetBinLowEdge(iZvtx+1),
								     fZvtxAxis->GetBinUpEdge(iZvtx+1)),
								//24,-0.5*TMath::Pi(),1.5*TMath::Pi(),24,-1.6,1.6);
								36,-0.5*TMath::Pi(),1.5*TMath::Pi(),32,-1.6,1.6);
	  fOutputList1 -> Add(fHistSP2AdPhidEta[iCent][iZvtx][iPt][jPt]);
	  
	  // same-sign track pairs
	  fHistSP2AdPhidEtaSS[iCent][iZvtx][iPt][jPt] = new TProfile2D(Form("fHist%sSP2AdPhidEtaSS_Cent%02d_Z%02d_Pt%02d_%02d",fUseRes?"Res":"",iCent,iZvtx,iPt,jPt), 
								  Form("%d-%d %%, %.1f<z<%.1f",
								       Int_t(fCentAxis->GetBinLowEdge(iCent+1)),
								       Int_t(fCentAxis->GetBinUpEdge(iCent+1)),
								       fZvtxAxis->GetBinLowEdge(iZvtx+1),
								       fZvtxAxis->GetBinUpEdge(iZvtx+1)),
								  36,-0.5*TMath::Pi(),1.5*TMath::Pi(),32,-1.6,1.6);
	  fOutputList1 -> Add(fHistSP2AdPhidEtaSS[iCent][iZvtx][iPt][jPt]);
	}
      }
     }
    }

 const Int_t nbins_dPhidEtaPt[] = {72, 32, 20, fNbinsPtTrig, fNbinsAssocPt, fNbinsCent};
 const Int_t nVar = sizeof(nbins_dPhidEtaPt) / sizeof(Int_t); 

 const TArrayD *aBin_Trig(fPtTrigAxis->GetXbins());
 const Int_t nBin_Trig(aBin_Trig->GetSize() - 1);
 Double_t dBin_Trig[nBin_Trig+1];

 const TArrayD *aBin_Asso(fPtAssocAxis->GetXbins());
 const Int_t nBin_Asso(aBin_Asso->GetSize() - 1);
 Double_t dBin_Asso[nBin_Asso+1];

 const TArrayD *aBin_Cen(fCentAxis->GetXbins());
 const Int_t nBin_Cen(aBin_Cen->GetSize() - 1);
 Double_t dBin_Cen[nBin_Cen+1];

 for (Int_t i=0; i<=nBin_Trig; ++i) dBin_Trig[i] = (*aBin_Trig)[i];
 for (Int_t i=0; i<=nBin_Asso; ++i) dBin_Asso[i] = (*aBin_Asso)[i];
 for (Int_t i=0; i<=nBin_Cen;  ++i) dBin_Cen[i]  = (*aBin_Cen)[i];


 fHistdPhidEtaPt = new AliTHn("fHistdPhidEtaPt", "fHistdPhidEtaPt", 1, nVar, nbins_dPhidEtaPt);

 fHistdPhidEtaPt->SetBinLimits(0, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistdPhidEtaPt->SetBinLimits(1, -1.6, 1.6);
 fHistdPhidEtaPt->SetBinLimits(2, -10, 10);
 fHistdPhidEtaPt->SetBinLimits(3, dBin_Trig);
 fHistdPhidEtaPt->SetBinLimits(4, dBin_Asso);
 fHistdPhidEtaPt->SetBinLimits(5, dBin_Cen);
 
 fHistdPhidEtaPt->SetVarTitle(0, "#Delta#phi");
 fHistdPhidEtaPt->SetVarTitle(1, "#Delta#eta");
 fHistdPhidEtaPt->SetVarTitle(2, "Vz");
 fHistdPhidEtaPt->SetVarTitle(3, "leading p_{T} GeV/c");
 fHistdPhidEtaPt->SetVarTitle(4, "asso p_{T} GeV/c");
 fHistdPhidEtaPt->SetVarTitle(5, "centrality");

 fHistdPhidEtaPt_SS = new AliTHn("fHistdPhidEtaPt_SS", "fHistdPhidEtaPt_SS", 1, nVar, nbins_dPhidEtaPt);

 fHistdPhidEtaPt_SS->SetBinLimits(0, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistdPhidEtaPt_SS->SetBinLimits(1, -1.6, 1.6);
 fHistdPhidEtaPt_SS->SetBinLimits(2, -10, 10);
 fHistdPhidEtaPt_SS->SetBinLimits(3, dBin_Trig);
 fHistdPhidEtaPt_SS->SetBinLimits(4, dBin_Asso);
 fHistdPhidEtaPt_SS->SetBinLimits(5, dBin_Cen);

 fHistdPhidEtaPt_SS->SetVarTitle(0, "#Delta#phi");
 fHistdPhidEtaPt_SS->SetVarTitle(1, "#Delta#eta");
 fHistdPhidEtaPt_SS->SetVarTitle(2, "Vz");
 fHistdPhidEtaPt_SS->SetVarTitle(3, "leading p_{T} GeV/c");
 fHistdPhidEtaPt_SS->SetVarTitle(4, "asso p_{T} GeV/c");
 fHistdPhidEtaPt_SS->SetVarTitle(5, "centrality");
 
 fHistdPhidEtaPt_Mixed = new AliTHn("fHistdPhidEtaPt_Mixed", "fHistdPhidEtaPt_Mixed", 1, nVar, nbins_dPhidEtaPt);

 fHistdPhidEtaPt_Mixed->SetBinLimits(0, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistdPhidEtaPt_Mixed->SetBinLimits(1, -1.6, 1.6);
 fHistdPhidEtaPt_Mixed->SetBinLimits(2, -10, 10);
 fHistdPhidEtaPt_Mixed->SetBinLimits(3, dBin_Trig);
 fHistdPhidEtaPt_Mixed->SetBinLimits(4, dBin_Asso);
 fHistdPhidEtaPt_Mixed->SetBinLimits(5, dBin_Cen);

 fHistdPhidEtaPt_Mixed->SetVarTitle(0, "#Delta#phi");
 fHistdPhidEtaPt_Mixed->SetVarTitle(1, "#Delta#eta");
 fHistdPhidEtaPt_Mixed->SetVarTitle(2, "Vz");
 fHistdPhidEtaPt_Mixed->SetVarTitle(3, "leading p_{T} GeV/c");
 fHistdPhidEtaPt_Mixed->SetVarTitle(4, "asso p_{T} GeV/c"); 
 fHistdPhidEtaPt_Mixed->SetVarTitle(5, "centrality"); 

 fHistdPhidEtaPt_Mixed_SS = new AliTHn("fHistdPhidEtaPt_Mixed_SS", "fHistdPhidEtaPt_Mixed_SS", 1, nVar, nbins_dPhidEtaPt);

 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(0, -0.5*TMath::Pi(), 1.5*TMath::Pi());
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(1, -1.6, 1.6);
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(2, -10, 10);
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(3, dBin_Trig);
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(4, dBin_Asso);
 fHistdPhidEtaPt_Mixed_SS->SetBinLimits(5, dBin_Cen);

 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(0, "#Delta#phi");
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(1, "#Delta#eta");
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(2, "Vz");
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(3, "leading p_{T} GeV/c");
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(4, "asso p_{T} GeV/c"); 
 fHistdPhidEtaPt_Mixed_SS->SetVarTitle(5, "centrality"); 

 
 const Int_t nbins_dTrig[] = {fNbinsPtTrig, 20, fNbinsCent};
 const Int_t nVar_Trig = sizeof(nbins_dTrig) / sizeof(Int_t);

 fHistTrig = new AliTHn("fHistTrig", "fHistTrig", 1, nVar_Trig, nbins_dTrig);
 fHistTrig->SetBinLimits(0, dBin_Trig);
 fHistTrig->SetBinLimits(1, -10, 10);
 fHistTrig->SetBinLimits(2, dBin_Cen);
 fHistTrig->SetVarTitle(0, "leading p_{T} GeV/c");
 fHistTrig->SetVarTitle(1, "Vz");
 fHistTrig->SetVarTitle(2, "centrality");


 fOutputList1->Add(fHistdPhidEtaPt); 
 fOutputList1->Add(fHistdPhidEtaPt_SS); 
 fOutputList1->Add(fHistdPhidEtaPt_Mixed); 
 fOutputList1->Add(fHistdPhidEtaPt_Mixed_SS); 
 fOutputList1->Add(fHistTrig);


  fHistACv2 = new TProfile("fHistACv2", "; centrality percentile; <Q^{2}_{V0A}.Q*^{2}_{V0C}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistACv2);
  fHistATv2 = new TProfile("fHistATv2", "; centrality percentile; <Q^{2}_{V0A}.Q*^{2}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistATv2);
  fHistCTv2 = new TProfile("fHistCTv2", "; centrality percentile; <Q^{2}_{V0C}.Q*^{2}_{Trk}>", fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray());
  fOutputList->Add(fHistCTv2);

  fHistV0Multiplicity = new TH2D("fHistV0Multiplicity", "V0 Multiplicity (Online vs Offline)",
				 500, 0, 50000,
				 500, 0, 50000);
  fHistV0Multiplicity -> SetXTitle("Multiplicity Offline");
  fHistV0Multiplicity -> SetYTitle("Multiplicity Online");
  fHistV0Multiplicity -> Sumw2();

  fHistITSMultiplicity = new TH1D("fHistITSMultiplicity", "ITS Multiplicity", 5000, 0, 25000);
  fHistITSMultiplicity -> SetXTitle("N_{Clusters}");
  fHistITSMultiplicity -> Sumw2();

  fHistCentrality = new TH1D("fHistCentrality", Form("%s_Centrality",fCentMethod.Data()), 400, -100, 300);
  fHistCentrality -> SetXTitle("Centrality  [%]");
  fHistCentrality -> Sumw2();


  fPileup1_Before = new TH2D("fPileup1_Before", "Pile_up_before (Online vs Offline)",
                                 500, 0, 50000,
                                 500, 0, 50000);
  fPileup1_Before -> SetXTitle("Multiplicity Online");
  fPileup1_Before -> SetYTitle("Multiplicity Offline");
  fPileup1_Before -> Sumw2();

  fPileup1_After = new TH2D("fPileup1_After", "Pile_up_After (Online vs Offline)",
                                 500, 0, 50000,
                                 500, 0, 50000);
  fPileup1_After -> SetXTitle("Multiplicity Online");
  fPileup1_After -> SetYTitle("Multiplicity Offline");
  fPileup1_After -> Sumw2();

  fPileup2_Before = new TH2D("fPileup2_Before", "ITS Multiplicity vs Num Tracklets", 5000, 0, 25000, 2000, 0 , 10000);
  fPileup2_Before -> SetXTitle("N_{Clusters}");
  fPileup2_Before -> SetYTitle("N_{Tracklets}");
  fPileup2_Before -> Sumw2();
  
  fPileup2_After = new TH2D("fPileup2_After", "ITS Multiplicity vs Num Tracklets", 5000, 0, 25000, 2000, 0 , 10000);
  fPileup2_After -> SetXTitle("N_{Clusters}");
  fPileup2_After -> SetYTitle("N_{Tracklets}");
  fPileup2_After -> Sumw2();

  fPileup3_Before_Low = new TH2D("fPileup3_Before_Low", "CL0 vs V0M", 100, 0, 100, 100, 0 , 100);
  fPileup3_Before_Low -> SetXTitle("Centrality CL0 [%]");
  fPileup3_Before_Low -> SetYTitle("Centrality V0M [%]");
  fPileup3_Before_Low -> Sumw2();

  fPileup3_After_Low = new TH2D("fPileup3_After_Low", "CL0 vs V0M", 100, 0, 100, 100, 0 , 100);
  fPileup3_After_Low -> SetXTitle("Centrality CL0 [%]");
  fPileup3_After_Low -> SetYTitle("Centrality V0M [%]");
  fPileup3_After_Low -> Sumw2();

  fPileup3_Before_High = new TH2D("fPileup3_Before_High", "CL0 vs V0M", 100, 0, 100, 100, 0 , 100);
  fPileup3_Before_High -> SetXTitle("Centrality CL0 [%]");
  fPileup3_Before_High -> SetYTitle("Centrality V0M [%]");
  fPileup3_Before_High -> Sumw2();

  fPileup3_After_High = new TH2D("fPileup3_After_High", "CL0 vs V0M", 100, 0, 100, 100, 0 , 100);
  fPileup3_After_High -> SetXTitle("Centrality CL0 [%]");
  fPileup3_After_High -> SetYTitle("Centrality V0M [%]");
  fPileup3_After_High -> Sumw2();


  fOutputList -> Add(fHistV0Multiplicity);
  fOutputList -> Add(fHistITSMultiplicity);
  fOutputList -> Add(fHistCentrality);

  fOutputList -> Add(fPileup1_Before);
  fOutputList -> Add(fPileup1_After);

  fOutputList -> Add(fPileup2_Before);
  fOutputList -> Add(fPileup2_After);

  fOutputList -> Add(fPileup3_Before_Low);
  fOutputList -> Add(fPileup3_After_Low);

  fOutputList -> Add(fPileup3_Before_High);
  fOutputList -> Add(fPileup3_After_High);



  fHistEvStat = new TH1D("fHistEvStat","Event cuts statistics",25,-0.5,24.5);
  fHistEvStat->SetXTitle("Cut index");
  fOutputList->Add(fHistEvStat);

  fHistNtrVsCent = new TH2D("fHistNtrVsCent","Number of tracks vs centrality",
			    fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(),
			    2000,0.,20000.);
  fOutputList->Add(fHistNtrVsCent);

  fHistV0CVsCent = new TH2D("fHistV0CVsCent","V0C vs centrality",
			    fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(),
			    500,0.,50000.);
  fOutputList->Add(fHistV0CVsCent);

  fHistESDvsTPC = new TH2D("fHistESDvsTPC","ESD vs TPC multiplicity",
			   500,0.,20000.,
			   500,0.,20000.);
  fOutputList->Add(fHistESDvsTPC);

  fHistCentVsZ = new TH2D("fHistCentVsZ","",20,-10,10,20,0,100);
  fOutputList->Add(fHistCentVsZ);
  fHistCentVsZMixed = new TH2D("fHistCentVsZMixed","Mixed events",20,-10,10,20,0,100);
  fOutputList->Add(fHistCentVsZMixed);
 
  PostData(1, fOutputList); 
  PostData(2, fOutputList1); 

  //================================================================
  flist_Res = dynamic_cast<TList*>(GetInputData(2));   
 
  if (!flist_Res) {
      printf("list_Res cannot be opened \n");
      return; 
  }

  fResACv2 = (TF1*)flist_Res->FindObject("fResACv2");
  fResATv2 = (TF1*)flist_Res->FindObject("fResATv2");
  fResCTv2 = (TF1*)flist_Res->FindObject("fResCTv2");
   
  
  Double_t parCent[6] = {0.0252185, 0.975157, 0.677622, 0.0290711, -0.000545769, 5.83219e-06};
    
  fHighCenCut = new TF1("fHighCenCut", "[0]+[1]*x + 6*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fHighCenCut->SetParameters(parCent);
    
  fLowCenCut = new TF1("fLowCenCut", "[0]+[1]*x - 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fLowCenCut->SetParameters(parCent);
    
    
  Double_t parV0[8] = {24.819, 0.952938, 0.0775265, 184.972, 6.9353, 0.0433712, -0.000192522, 7.67571e-07};
  fV0MultOfOnCut = new TF1("fV0MultOfOnCut", "[0]+[1]*x - 5.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0MultOfOnCut->SetParameters(parV0);

  ftrk = new TClonesArray("AliBasicParticleST",500);
  
  const Int_t nzvtx = 10;
  Double_t zvtxbins[nzvtx+1] = {-10.,-8.,-6.,-4.,-2.,0.,2.,4.,6.,8.,10.};
 
  const Int_t nCentralityBins  = 14;
  Double_t centBins[] = {0.,1.,2.,3.,4.,5., 10.,20.,30.,40.,50.,60.,70.,80.,90.};

  fPoolMgr = new AliEventPoolManager(1000, 50000, nCentralityBins, centBins, nzvtx, zvtxbins);
  //fPoolMgr = new AliEventPoolManager(500, 30000, fNbinsCent, (Double_t*)fCentAxis->GetXbins()->GetArray(), nzvtx, zvtxbins);
  if (!fPoolMgr) return;
  fPoolMgr->SetTargetValues(50000, 0.1, 5);

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::UserExec(Option_t *) {

  Int_t cutIndex = 0;
  fAOD = dynamic_cast<AliAODEvent *>(InputEvent());  
  if (!fAOD) return;  

  fRunN = fAOD->GetRunNumber();

  //Exclude some runs which are not in calibrated list
  //if(fRunN == 246928 || fRunN == 246810) return;
  //if(fRunN == 246766 || fRunN == 246271) return;
  //if(fRunN == 246185 || fRunN == 246180) return;
  //if(fRunN == 246052 || fRunN == 245923) return;
  //if(fRunN == 245702) return;

  fHistEvStat->Fill(cutIndex++); // Index = 0

  // Trigger selection
  TString trigStr(fAOD->GetFiredTriggerClasses());

  if (trigStr.Contains("CINT7-B-NOPF")) fHistEvStat->Fill(cutIndex); //Index = 1
  cutIndex++;

  AliAODVZERO* aodV0 = fAOD->GetVZEROData();
  if (trigStr.Contains("CINT7-B-NOPF") &&
      aodV0->GetV0ADecision()==1 &&
      aodV0->GetV0CDecision()==1) fHistEvStat->Fill(cutIndex); // Index = 2
  cutIndex++;
  
  if (!(IsTriggerFired())) return;
  fHistEvStat->Fill(cutIndex++);     // Index = 3
  
  fv0mult = GetV0Multiplicity();
  fv0multonline = aodV0->GetTriggerChargeA()+aodV0->GetTriggerChargeC();
  fHistV0Multiplicity  -> Fill(fv0mult,fv0multonline);
  fitsmult = GetITSMultiplicity();
  fHistITSMultiplicity -> Fill(fitsmult);

  Double_t percentile;
  AliMultSelection *multSelection = (AliMultSelection *)fAOD->FindListObject("MultSelection");
  if(!multSelection) return;

  if(fCentMethod == "") {
   if(fForceCL1) fCentMethod = "CL1";
   else fCentMethod = "V0M";
   fHistEvStat->Fill(cutIndex);   // Index = 4
  }
  //percentile = multSelection->GetMultiplicityPercentile(fCentMethod.Data());
  percentile = multSelection->GetMultiplicityPercentile(fCentMethod);
  fHistCentrality->Fill(percentile);

  fv0mpercentile = multSelection->GetMultiplicityPercentile("V0M");
  fv0meqpercentile = multSelection->GetMultiplicityPercentile("V0MEq");
  fv0apercentile = multSelection->GetMultiplicityPercentile("V0A");
  fv0cpercentile = multSelection->GetMultiplicityPercentile("V0C");
  fcl0percentile = multSelection->GetMultiplicityPercentile("CL0");
  fcl1percentile = multSelection->GetMultiplicityPercentile("CL1");
  fSPDpercentile = multSelection->GetMultiplicityPercentile("SPDTracklets");

  ++cutIndex;             

  Int_t centBin = GetCentBin(percentile);
  if (centBin<0) return;
  fHistEvStat->Fill(cutIndex++); // Index = 5

  if (fAOD->IsIncompleteDAQ()) return;
  fHistEvStat->Fill(cutIndex++); // Index = 6

  // Vertex selection
  const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
  if (trkVtx->GetNContributors() < 2) return;
  fHistEvStat->Fill(cutIndex++); // Index = 7
  Double_t covTrk[6] = {0};
  trkVtx->GetCovarianceMatrix(covTrk);
  Double_t zvtxTrk = trkVtx->GetZ();
  
  const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
  if (spdVtx->GetNContributors()<=0) return;
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (spdVtx->IsFromVertexerZ() && (zRes>0.25)) return;
  fHistEvStat->Fill(cutIndex++); // Index = 8

  fzvtx = spdVtx->GetZ();
  if (fzvtx >= fZvtxAxis->GetXmax()) return;
  if (fzvtx <= fZvtxAxis->GetXmin()) return;
  fHistEvStat->Fill(cutIndex++); // Index = 9

  Double_t dz = zvtxTrk - fzvtx;
  Double_t errTot = TMath::Sqrt(covTrk[5]+cov[5]);
  Double_t errTrk = TMath::Sqrt(covTrk[5]);
  Double_t nsigTot = dz/errTot;
  Double_t nsigTrk = dz/errTrk;
  if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrk)>20)
    return; // bad vertexing
  fHistEvStat->Fill(cutIndex++); // Index = 10
  
  //fPileup1_Before->Fill(fv0multonline,fV0MultOfOnCut->Eval(fv0mult));
  fPileup1_Before->Fill(fv0multonline,fv0mult);

  if (fRemovePileup) {
    if (fv0multonline < fV0MultOfOnCut->Eval(fv0mult))
      return;
    fPileup1_After->Fill(fv0multonline,fv0mult);
  }
  fHistEvStat->Fill(cutIndex++); // Index = 11

  AliAODTracklets* aodTrkl = (AliAODTracklets*)fAOD->GetTracklets();
  nITSTrkls = aodTrkl->GetNumberOfTracklets();

  fPileup2_Before->Fill(fitsmult,nITSTrkls);
  if (fRemovePileup2) {
    if (fitsmult > (400.+4.*Double_t(nITSTrkls)))
      return;
    fPileup2_After->Fill(fitsmult,nITSTrkls);
  }
  fHistEvStat->Fill(cutIndex++); // Index = 12

  fPileup3_Before_Low->Fill(fcl0percentile, fv0mpercentile);
  fPileup3_Before_High->Fill(fcl0percentile,fv0mpercentile);
  if (fRemovePileup3) {
    if (fcl0percentile < fLowCenCut->Eval(fv0mpercentile))
      return;
    fPileup3_After_Low->Fill(fcl0percentile, fv0mpercentile);
    if (fcl0percentile > fHighCenCut->Eval(fv0mpercentile))
      return;
    fPileup3_After_High->Fill(fcl0percentile,fv0mpercentile);
  }

  fHistEvStat->Fill(cutIndex++); // Index = 13

  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multEsd = ((AliAODHeader*)fAOD->GetHeader())->GetNumberOfESDTracks();
  Int_t multTpc = 0;
  for (Int_t it = 0; it < nTracks; it++) {
    AliAODTrack* aodTrk = (AliAODTrack*)fAOD->GetTrack(it);
    if (!aodTrk) continue;
    if (aodTrk->TestFilterBit(128)) multTpc++;
  }
  fHistESDvsTPC->Fill(multTpc,multEsd);
  Double_t multESDTPCDif = Double_t(multEsd) - Double_t(multTpc)*3.38;
  if (multESDTPCDif > 15000.) return;
  fHistEvStat->Fill(cutIndex++); // Index = 14
  
  if (!ComputeQ(fAOD,fzvtx)) return;
  
  fHistEvStat->Fill(cutIndex++); // Index = 15

  Int_t zvtxBin = GetZvtxBin(fzvtx);

  fHistV0CVsCent->Fill(percentile,fAOD->GetVZEROData()->GetMTotV0C());
  
  //ProcessEvent(percentile, centBin, zvtxBin);

  fHistNtrVsCent->Fill(percentile,nTracks);
 
  // Calculate the q vect resolutions
  Double_t resA2=1.,resC2=1.,resT2=1.;
  Double_t resA3=1.,resC3=1.,resT3=1.;
  if (fUseRes) CalcResolutions(percentile, resA2, resC2, resT2);

//=================
  TObjArray *selectedTrackArray = new TObjArray; selectedTrackArray->SetOwner(kTRUE);
  selectedTrackArray = GetAcceptedTracks(fAOD,selectedTrackArray); 
 
 
  FillHistogramsdPhidEta(selectedTrackArray,centBin,percentile,zvtxBin,
                             resA2, resC2, resT2,
                             resA3, resC3, resT3); 
  
  FillHistogramsdPhidEtaMixed(selectedTrackArray, percentile, zvtxBin);
 
  for (Int_t iTr=0; iTr<nTracks; iTr++) {
    AliAODTrack *track = (AliAODTrack*) fAOD->GetTrack(iTr);
    if (!track) continue;
    if (track->TestFilterBit(fFilterBit) && TMath::Abs(track->Eta()) < 0.8 && track->GetTPCNcls() >= fTPCNcls) {
      FillHistogramsV2(track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
		       resA2, resC2, resT2, 0);
      FillHistogramsV2(track->Pt(),track->Eta(),track->Phi(),centBin,percentile,zvtxBin,
		       resA2, resC2, resT2, 1);
    }
  }

  Double_t Qv0aQv0c2  = Qxa2Cor*Qxc2Cor  + Qya2Cor*Qyc2Cor;

  fHistACv2->Fill(percentile, Qv0aQv0c2);

  ((TH2D*)(fOutputList->FindObject("fHist_Qx2_V0A")))->Fill(percentile,Qxa2);
  ((TH2D*)(fOutputList->FindObject("Corr_fHist_Qx2_V0A")))->Fill(percentile,Qxa2Cor);

  ((TH2D*)(fOutputList->FindObject("fHist_Qy2_V0A")))->Fill(percentile,Qya2);
  ((TH2D*)(fOutputList->FindObject("Corr_fHist_Qy2_V0A")))->Fill(percentile,Qya2Cor);

  ((TH2D*)(fOutputList->FindObject("fHist_Qx2_V0C")))->Fill(percentile,Qxc2);
  ((TH2D*)(fOutputList->FindObject("Corr_fHist_Qx2_V0C")))->Fill(percentile,Qxc2Cor);

  ((TH2D*)(fOutputList->FindObject("fHist_Qy2_V0C")))->Fill(percentile,Qyc2);
  ((TH2D*)(fOutputList->FindObject("Corr_fHist_Qy2_V0C")))->Fill(percentile,Qyc2Cor);

  selectedTrackArray->Clear();
  delete selectedTrackArray;
  PostData(1, fOutputList); 
  PostData(2, fOutputList1); 
}


void AliAnalysisTaskSEPbPbCorrelationsJetV2::FillHistogramsV2(Double_t pt,Double_t eta,Double_t phi,Int_t centrality,Double_t percentile,Int_t zvtxBin,
					     Double_t resA2, Double_t resC2, Double_t resT2,
					     Int_t index)
{
 Double_t u2x = TMath::Cos(2.*phi);
 Double_t u2y = TMath::Sin(2.*phi);
 if(index == 0) fHistSP2A[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxa2Cor+u2y*Qya2Cor)/resA2); 
 if(index == 1) fHistSP2A[centrality][zvtxBin][index]->Fill(pt,(u2x*Qxc2Cor+u2y*Qyc2Cor)/resC2); 
}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::FillHistogramsdPhidEta(TObjArray *selectedArray, Int_t centrality,Double_t percentile,Int_t zvtxBin,
						   Double_t resA2, Double_t resC2, Double_t resT2,
						   Double_t resA3, Double_t resC3, Double_t resT3)
{
 Double_t binscont_trig[3];
 Double_t binscont[6];

 for(Int_t i = 0; i < selectedArray->GetEntriesFast(); i++) {
  AliBasicParticleST *trigger = (AliBasicParticleST*)selectedArray->At(i);  
  if (!trigger)    continue;
  Int_t trigID = trigger->GetID();
  Double_t triggerPt = trigger->Pt();
  Double_t triggerEta = trigger->Eta();
  Double_t triggerPhi = trigger->Phi();

  Int_t ptBin = fPtTrigAxis->FindBin(triggerPt);

  if (ptBin<1 || ptBin>fNbinsPtTrig) continue;

  binscont_trig[0] = triggerPt;
  binscont_trig[1] = fzvtx;
  binscont_trig[2] = percentile;
  fHistTrig->Fill(binscont_trig,0);

  Double_t u2x = TMath::Cos(2.*triggerPhi);
  Double_t u2y = TMath::Sin(2.*triggerPhi);

  for (Int_t j = 0; j < selectedArray->GetEntriesFast(); j++) {
   AliBasicParticleST *associate = (AliBasicParticleST*)selectedArray->At(j); 
   if (!associate) continue;
   if (trigID == associate->GetID())  continue;   
   Int_t assocPtBin = fPtAssocAxis->FindBin(associate->Pt());

   if(fPtOrder && triggerPt < associate->Pt()) continue; 

   if (assocPtBin<1 || assocPtBin>fNbinsAssocPt) continue;  
   Double_t dphi = triggerPhi - associate->Phi();
   if (dphi >  1.5*TMath::Pi()) dphi -= TMath::TwoPi();
   if (dphi < -0.5*TMath::Pi()) dphi += TMath::TwoPi();
   Double_t deta = triggerEta - associate->Eta();
   if(!fSameSign) fHistSP2AdPhidEta[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u2x*Qxa2Cor+u2y*Qya2Cor)/resA2); 
   binscont[0] = dphi;
   binscont[1] = deta;
   binscont[2] = fzvtx;
   binscont[3] = triggerPt;
   binscont[4] = associate->Pt();
   binscont[5] = percentile;
   if(!fSameSign) fHistdPhidEtaPt->Fill(binscont,0);
   if (trigger->Charge()*associate->Charge()>0) {
     if(fMode == "ref_V0A") fHistSP2AdPhidEtaSS[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u2x*Qxa2Cor+u2y*Qya2Cor)/resA2);
     if(fMode == "ref_V0C") fHistSP2AdPhidEtaSS[centrality][zvtxBin][ptBin-1][assocPtBin-1]->Fill(dphi,deta,(u2x*Qxc2Cor+u2y*Qyc2Cor)/resC2);
     fHistdPhidEtaPt_SS->Fill(binscont,0);
   }
  }
 }
}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::FillHistogramsdPhidEtaMixed(TObjArray *selectedArray, Double_t percentile,Int_t zvtxBin)
{
  // check if mixed event pool is ready
  AliEventPool* pool = fPoolMgr->GetEventPool(percentile,fzvtx);
  if (!pool){
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", percentile,
                  fzvtx));
  }
  
  Bool_t poolReady = (pool->IsReady() || pool->NTracksInPool() > 5000 || pool->GetCurrentNEvents() > 5);
    
  Double_t binscont[6];
  if(poolReady)
  {
   fHistCentVsZMixed->Fill(fzvtx,percentile,pool->GetCurrentNEvents());
   for(Int_t j = 0; j < selectedArray->GetEntriesFast(); ++j) {
    AliBasicParticleST *trigger = (AliBasicParticleST*)selectedArray->At(j);
    Double_t triggerPt = trigger->Pt();
    Double_t triggerEta = trigger->Eta();
    Double_t triggerPhi = trigger->Phi();
 
    Int_t ptBin = fPtTrigAxis->FindBin(triggerPt);
    if (ptBin<1 || ptBin>fNbinsPtTrig) continue;
    
    for(Int_t jMix=0; jMix<pool->GetCurrentNEvents(); jMix++) {
     TObjArray *mixEvents = pool->GetEvent(jMix);
     for (Int_t jTrk=0; jTrk<mixEvents->GetEntriesFast(); jTrk++) {       
      AliBasicParticleST* associate = (AliBasicParticleST*)mixEvents->At(jTrk);

      Int_t assocPtBin = fPtAssocAxis->FindBin(associate->Pt());
      if (assocPtBin<1 || assocPtBin>fNbinsAssocPt) continue;

      Double_t dphi = triggerPhi - associate->Phi();
      if (dphi >  1.5*TMath::Pi()) dphi -= TMath::TwoPi();
      if (dphi < -0.5*TMath::Pi()) dphi += TMath::TwoPi();   
         
      binscont[0] = dphi;
      binscont[1] = triggerEta - associate->Eta();
      binscont[2] = fzvtx;
      binscont[3] = triggerPt;
      binscont[4] = associate->Pt();
      binscont[5] = percentile;
      if(!fSameSign) fHistdPhidEtaPt_Mixed->Fill(binscont,0);
   
      if (trigger->Charge()*associate->Charge()>0) {
       fHistdPhidEtaPt_Mixed_SS->Fill(binscont,0);
      }
     }
    }
   }
  }
   
  TObjArray* tracksClone=CloneTrack(selectedArray);
  pool->UpdatePool(tracksClone);

}

//====================================================================================================================================================

Bool_t AliAnalysisTaskSEPbPbCorrelationsJetV2::IsTriggerFired() {
  
  Bool_t isSelected = kFALSE;
  isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7); 
  return isSelected;
}

//====================================================================================================================================================

Float_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetV0Multiplicity() {
  
  Float_t multiplicity=0;
  for (Int_t iChannel=0; iChannel<64; iChannel++) multiplicity += fAOD->GetVZEROData()->GetMultiplicity(iChannel);
  return multiplicity;

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetCentBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsCentrality) {
    AliInfo(Form("WARNING : only %d centrality bins (out of the %d proposed) will be considered",fNMaxBinsCentrality,nBins));
    nBins = fNMaxBinsCentrality;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one centrality bin must be considered");
    nBins = 1;
  }
  
  fNbinsCent = nBins;
  fCentAxis  = new TAxis(fNbinsCent, limits);

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetZvtxBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsZvtx) {
    AliInfo(Form("WARNING : only %d Zvtx bins (out of the %d proposed) will be considered",fNMaxBinsZvtx,nBins));
    nBins = fNMaxBinsZvtx;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one Zvtx bin must be considered");
    nBins = 1;
  }
  
  fNbinsZvtx = nBins;
  fZvtxAxis  = new TAxis(fNbinsZvtx, limits);

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetPtBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsPt) {
    AliInfo(Form("WARNING : only %d pt bins (out of the %d proposed) will be considered",fNMaxBinsPt,nBins));
    nBins = fNMaxBinsPt;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one pt bin must be considered");
    nBins = 1;
  }
  
  fNbinsPt = nBins;
  fPtAxis  = new TAxis(fNbinsPt, limits);

}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetTrigPtBinning(Int_t nBins, Double_t *limits) {

  fNbinsPtTrig = nBins;
  fPtTrigAxis  = new TAxis(fNbinsPtTrig, limits);

}



//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetAssocPtBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsAssocPt) {
    AliInfo(Form("WARNING : only %d assoc pt bins (out of the %d proposed) will be considered",fNMaxBinsAssocPt,nBins));
    nBins = fNMaxBinsAssocPt;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one assoc pt bin must be considered");
    nBins = 1;
  }
  
  fNbinsAssocPt = nBins;
  fPtAssocAxis  = new TAxis(fNbinsAssocPt, limits);

}

//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::SetEtaBinning(Int_t nBins, Double_t *limits) {

  if (nBins>fNMaxBinsEta) {
    AliInfo(Form("WARNING : only %d pt bins (out of the %d proposed) will be considered",fNMaxBinsEta,nBins));
    nBins = fNMaxBinsEta;
  }
  if (nBins<=0) {
    AliInfo("WARNING : at least one pt bin must be considered");
    nBins = 1;
  }
  
  fEtaAxis  = new TAxis(nBins, limits);

}

//====================================================================================================================================================

Int_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetCentBin() {

  Double_t percentile;
    AliMultSelection *multSelection = (AliMultSelection *)fAOD->FindListObject("MultSelection");
    percentile = multSelection->GetMultiplicityPercentile(fCentMethod.Data());

  Int_t bin = fCentAxis->FindBin(percentile) - 1;
  if (bin >= fNbinsCent) bin = -1;
  return bin;
  
}

//====================================================================================================================================================

Int_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetCentBin(Double_t percentile) {

  Int_t bin = fCentAxis->FindBin(percentile) - 1;
  if (bin >= fNbinsCent) bin = -1;
  return bin;
  
}

//====================================================================================================================================================

Int_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetZvtxBin(Double_t zvtx) {

  Int_t bin = fZvtxAxis->FindBin(zvtx) - 1;
  if (bin >= fNbinsZvtx) bin = -1;
  return bin;
  
}

//====================================================================================================================================================

Double_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetITSMultiplicity() {

  Double_t multiplicity = ((AliVAODHeader*)fAOD->GetHeader())->GetNumberOfITSClusters(0)+
    ((AliVAODHeader*)fAOD->GetHeader())->GetNumberOfITSClusters(1);

  return multiplicity;

}

Double_t AliAnalysisTaskSEPbPbCorrelationsJetV2::CalcCorrectedPhi(Double_t phi, Double_t dPhi) const{
  phi += 39./34.*dPhi;
  if (phi < 0.)  phi+=2*TMath::Pi();
  if (phi > 2.*TMath::Pi()) phi-=2*TMath::Pi();
  return phi;
}


//====================================================================================================================================================

void AliAnalysisTaskSEPbPbCorrelationsJetV2::Terminate(Option_t *) {

  if (fPoolMgr)    delete fPoolMgr;

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }  

  fOutputList1 = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputList1) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//__________________________________________________________
Short_t AliAnalysisTaskSEPbPbCorrelationsJetV2::GetVertexZ(Double_t vtxZ) const
{
  // this method gives the z vtx bin used to access flow vector calibrations
    Short_t zvtx = -10;
    
    if (vtxZ >= -10. && vtxZ < -8.)
        zvtx = 0;
    if (vtxZ >= -8. && vtxZ < -6.)
        zvtx = 1;
    if (vtxZ >= -6. && vtxZ < -4.)
        zvtx = 2;
    if (vtxZ >= -4. && vtxZ < -3.)
        zvtx = 3;
    if (vtxZ >= -3. && vtxZ < -2.)
        zvtx = 4;
    if (vtxZ >= -2. && vtxZ < -1.)
        zvtx = 5;
    if (vtxZ >= -1. && vtxZ < 0)
        zvtx = 6;
    if (vtxZ >= 0 && vtxZ < 1.)
        zvtx = 7;
    if (vtxZ >= 1. && vtxZ < 2.)
        zvtx = 8;
    if (vtxZ >= 2. && vtxZ < 3.)
        zvtx = 9;
    if (vtxZ >= 3. && vtxZ < 4.)
        zvtx = 10;
    if (vtxZ >= 4. && vtxZ < 6.)
        zvtx = 11;
    if (vtxZ >= 6. && vtxZ < 8.)
        zvtx = 12;
    if (vtxZ >= 8. && vtxZ <= 10.)
        zvtx = 13;
    
    return zvtx;
    
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEPbPbCorrelationsJetV2::ComputeQ(AliAODEvent* aod, Double_t Zvtx)
{  

    AliAODTracklets* aodTrkl = (AliAODTracklets*)aod->GetTracklets();
    Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
    
    AliAODVZERO* aodV0 = (AliAODVZERO*)aod->GetVZEROData();

    Int_t iCentV0 = Int_t(fv0mpercentile);
    if (iCentV0 >= 90)
      return kFALSE;
    
    Int_t iCentSPD = Int_t(fcl1percentile);
    if (iCentSPD >= 90)
      return kFALSE;
    
    
    Short_t zvt = GetVertexZ(Zvtx);
    if (zvt < 0)
      return kFALSE;
   
    zvt = 0;   
 
    //Tracklets
    Qx2trkcut = Qy2trkcut = 0;
    multtrkcut = 0;
    for (Int_t it = 0; it < nITSTrkls; it++){

      // We do not use these cuts anymore
      //        if (TMath::Abs(aodTrkl->GetEta(it)) > 0.5 || TMath::Abs(aodTrkl->GetDeltaPhi(it)) >= 0.005)
      //            continue;
        
        Double_t corPhi = CalcCorrectedPhi(aodTrkl->GetPhi(it), aodTrkl->GetDeltaPhi(it));
        
        Qx2trkcut += TMath::Cos(2.*corPhi);
        Qy2trkcut += TMath::Sin(2.*corPhi);
        
        multtrkcut++;
    }
    
    //V0 info
    Qxa2 = Qya2 = 0;
    Qxc2 = Qyc2 = 0;
    
    sumMa = sumMc = 0;
    
    for (Int_t iV0 = 0; iV0 < 64; iV0++) {
        
        Double_t phiV0 = TMath::PiOver4()*(0.5 + iV0 % 8);
        Float_t multv0 = aodV0->GetMultiplicity(iV0);
        
        if (iV0 < 32){
            
            Double_t multCorC = -10;
            
            if (iV0 < 8)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(1);
            else if (iV0 >= 8 && iV0 < 16)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(9);
            else if (iV0 >= 16 && iV0 < 24)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(17);
            else if (iV0 >= 24 && iV0 < 32)
                multCorC = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(25);
            
            if (multCorC < 0){
              //  cout<<"Problem with multiplicity in V0C"<<endl;
                continue;
            }
            
            Qxc2 += TMath::Cos(2.*phiV0) * multCorC;
            Qyc2 += TMath::Sin(2.*phiV0) * multCorC;
 
            sumMc = sumMc + multCorC;
            
        } else {
            
            Double_t multCorA = -10;
            
            if (iV0 >= 32 && iV0 < 40)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(33);
            else if (iV0 >= 40 && iV0 < 48)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(41);
            else if (iV0 >= 48 && iV0 < 56)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(49);
            else if (iV0 >= 56 && iV0 < 64)
                multCorA = multv0/fMultV0->GetBinContent(iV0+1)*fMultV0->GetBinContent(57);
            
            if (multCorA < 0){
            //    cout<<"Problem with multiplicity in V0A"<<endl;
                continue;
            }
            
            Qxa2 += TMath::Cos(2.*phiV0) * multCorA;
            Qya2 += TMath::Sin(2.*phiV0) * multCorA;

            sumMa = sumMa + multCorA;
            
        }
    }

    if (sumMa <= 0 || sumMc <= 0 || multtrkcut <= 0)
      return kFALSE;
    
    
    //might need to check the sigma values != 0
    if (fQy2sV0A[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQy2sV0C[zvt]->GetBinContent(iCentSPD+1)<1e-8)
      return kFALSE;
    if (fQx2sV0A[zvt]->GetBinContent(iCentSPD+1)<1e-8 ||
	fQx2sV0C[zvt]->GetBinContent(iCentSPD+1)<1e-8)
      return kFALSE;
    
    Qya2Cor = (Qya2 - fQy2mV0A[zvt]->GetBinContent(iCentSPD+1))/fQy2sV0A[zvt]->GetBinContent(iCentSPD+1);
    Qxa2Cor = (Qxa2 - fQx2mV0A[zvt]->GetBinContent(iCentSPD+1))/fQx2sV0A[zvt]->GetBinContent(iCentSPD+1);
    
    Qyc2Cor = (Qyc2 - fQy2mV0C[zvt]->GetBinContent(iCentSPD+1))/fQy2sV0C[zvt]->GetBinContent(iCentSPD+1);
    Qxc2Cor = (Qxc2 - fQx2mV0C[zvt]->GetBinContent(iCentSPD+1))/fQx2sV0C[zvt]->GetBinContent(iCentSPD+1);

    return kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskSEPbPbCorrelationsJetV2::OpenInfoCalbration(Int_t run)
{
  flist_contQ = dynamic_cast<TList*>(GetInputData(1));
          
  if(!flist_contQ){
        printf("OADB V0-Trkl calibration list cannot be opened\n");
        return;
  }
    
  if (!cont) cont = (AliOADBContainer*) flist_contQ->FindObject("hMultV0BefCorPfpx");
  if(!cont){
        printf("OADB object hMultV0BefCorPfpx is not available in the file\n");
        return;
  }
  if(!(cont->GetObject(run))){
        printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
        return;
  }
  fMultV0 = ((TH1D*) cont->GetObject(run));
       
  for (Int_t k = 0; k < 1; k++){
        
        if (!contQx2am[k]) contQx2am[k] = (AliOADBContainer*) flist_contQ->FindObject("fqxa2m");
        if(!contQx2am[k]){
            printf("OADB object fqxa2m is not available in the file\n");
            return;
        }
        if(!(contQx2am[k]->GetObject(run))){
            printf("OADB object fqxa2m is not available for run %i\n", run);
            return;
        }
        fQx2mV0A[k]= ((TH1D*) contQx2am[k]->GetObject(run));
        
        if (!contQy2am[k]) contQy2am[k] = (AliOADBContainer*) flist_contQ->FindObject("fqya2m");
        if(!contQy2am[k]){
            printf("OADB object fqya2m is not available in the file\n");
            return;
        }
        if(!(contQy2am[k]->GetObject(run))){
            printf("OADB object fqya2m is not available for run %i\n", run);
            return;
        }
        fQy2mV0A[k]= ((TH1D*) contQy2am[k]->GetObject(run));
        
        if (!contQx2as[k]) contQx2as[k] = (AliOADBContainer*) flist_contQ->FindObject("fqxa2s");
        if(!contQx2as[k]){
            printf("OADB object fqxa2s is not available in the file\n");
            return;
        }
        if(!(contQx2as[k]->GetObject(run))){
            printf("OADB object fqxa2s is not available for run %i\n", run);
            return;
        }
        fQx2sV0A[k]= ((TH1D*) contQx2as[k]->GetObject(run));
        
        
        if (!contQy2as[k]) contQy2as[k] = (AliOADBContainer*) flist_contQ->FindObject("fqya2s");
        if(!contQy2as[k]){
            printf("OADB object fqya2s is not available in the file\n");
            return;
        }
        if(!(contQy2as[k]->GetObject(run))){
            printf("OADB object fqya2s is not available for run %i\n", run);
            return;
        }
        fQy2sV0A[k]= ((TH1D*) contQy2as[k]->GetObject(run));
        
        
        if (!contQx2cm[k]) contQx2cm[k] = (AliOADBContainer*) flist_contQ->FindObject("fqxc2m");
        if(!contQx2cm[k]){
            printf("OADB object fqxc2m is not available in the file\n");
            return;
        }
        if(!(contQx2cm[k]->GetObject(run))){
            printf("OADB object fqxc2m is not available for run %i\n", run);
            return;
        }
        fQx2mV0C[k]= ((TH1D*) contQx2cm[k]->GetObject(run));
        
        
        if (!contQy2cm[k]) contQy2cm[k] = (AliOADBContainer*) flist_contQ->FindObject("fqyc2m");
        if(!contQy2cm[k]){
            printf("OADB object fqyc2m is not available in the file\n");
            return;
        }
        if(!(contQy2cm[k]->GetObject(run))){
            printf("OADB object fqyc2m is not available for run %i\n", run);
            return;
        }
        fQy2mV0C[k]= ((TH1D*) contQy2cm[k]->GetObject(run));
        
        
        if (!contQx2cs[k]) contQx2cs[k] = (AliOADBContainer*) flist_contQ->FindObject("fqxc2s");
        if(!contQx2cs[k]){
            printf("OADB object fqxc2s is not available in the file\n");
            return;
        }
        if(!(contQx2cs[k]->GetObject(run))){
            printf("OADB object fqxc2s is not available for run %i\n", run);
            return;
        }
        fQx2sV0C[k]= ((TH1D*) contQx2cs[k]->GetObject(run));
        
        
        if (!contQy2cs[k]) contQy2cs[k] = (AliOADBContainer*) flist_contQ->FindObject("fqyc2s");
        if(!contQy2cs[k]){
            printf("OADB object fqyc2s is not available in the file\n");
            return;
        }
        if(!(contQy2cs[k]->GetObject(run))){
            printf("OADB object fqyc2s is not available for run %i\n", run);
            return;
        }
        fQy2sV0C[k]= ((TH1D*) contQy2cs[k]->GetObject(run));
    }

    
}

void AliAnalysisTaskSEPbPbCorrelationsJetV2::CalcResolutions(Double_t percentile,
					      Double_t &resA2, Double_t &resC2, Double_t &resT2)
{
    // average resolution
    // v2
      Double_t ac2 = fResACv2->Eval(percentile);
      Double_t at2 = fResATv2->Eval(percentile);
      Double_t ct2 = fResCTv2->Eval(percentile);
      if (ac2 <= 0 || at2 <= 0 || ct2 <= 0) {
	resA2 = resC2 = resT2 = 1e6;
      }
      else {
	resA2 = TMath::Sqrt(ac2*at2/ct2);
	resC2 = TMath::Sqrt(ac2*ct2/at2);
	resT2 = TMath::Sqrt(at2*ct2/ac2);
      }
}

TObjArray *AliAnalysisTaskSEPbPbCorrelationsJetV2::GetAcceptedTracks(AliAODEvent *fAOD, TObjArray*tracks)
{
 Int_t nTracks = fAOD->GetNumberOfTracks();
 for (Int_t i = 0; i < nTracks; i++) {
  AliAODTrack *aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i));
  if (!aodTrack)      continue;  
  if (aodTrack->Charge() == 0)      continue;
  //if (aodTrack->TestFilterBit(768) && TMath::Abs(aodTrack->Eta()) < 0.8 && aodTrack->GetTPCNcls() >= 70) {
  if (aodTrack->TestFilterBit(fFilterBit) && TMath::Abs(aodTrack->Eta()) < 0.8 && aodTrack->GetTPCNcls() >= fTPCNcls) {
   tracks->Add(new AliBasicParticleST(aodTrack->Charge(), aodTrack->Eta(), aodTrack->Phi(), aodTrack->Pt(), aodTrack->GetID(), -999, -999, 0, 1));
  }
 }
 return tracks;
}

TObjArray* AliAnalysisTaskSEPbPbCorrelationsJetV2::CloneTrack(TObjArray*selectedTrackArray){
  TObjArray *tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  for (Int_t i = 0; i < selectedTrackArray->GetEntriesFast(); i++) {
    AliBasicParticleST *particle =  (AliBasicParticleST *)selectedTrackArray->At(i);
    tracksClone->Add(new AliBasicParticleST(particle->Charge(), particle->Eta(), particle->Phi(), particle->Pt(),
                                              particle->GetID(), particle->GetIDFirstDaughter(),
                                              particle->GetIDSecondDaughter(), particle->WhichCandidate(),
                                              particle->Multiplicity()));
  }

  return tracksClone;
}

