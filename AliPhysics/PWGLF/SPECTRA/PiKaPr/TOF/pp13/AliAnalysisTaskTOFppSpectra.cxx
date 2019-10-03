/**********************************************************************************
 * Date   :: 13/03/2018                                                           *
 * Author :: Pranjal Sarma & B. Bhattacharjee, Gauhati University, India          *
 * Analysis Task for Pi,Ka & Pr vs Multiplicity in pp @ 13 TeV with TOF        *
 **********************************************************************************/
#define LOG_NO_INFO
#define LOG_NO_DEBUG

#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TFile.h"
#include "AliPIDResponse.h"
#include <AliPID.h>
#include <TSpline.h>
#include "AliESDtrack.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliVTrack.h"
#include "AliAnalysisTaskTOFppSpectra.h"
#include "AliVEvent.h"

#include "AliCentrality.h"
//#include "TProfile.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliESDVZERO.h"
#include "AliESDtrackCuts.h"
#include "AliTrackerBase.h"
#include "fstream"


#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"

#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"

#include "TRandom.h"
#include "TGrid.h"


//#endif
Double_t V0mpc;
Int_t value_Sigma, value_Slope;
// Authors:Pranjal Sarma  (22/03/17)

ClassImp(AliAnalysisTaskTOFppSpectra)
//________________________________________________________________________
AliAnalysisTaskTOFppSpectra::AliAnalysisTaskTOFppSpectra()
  : AliAnalysisTaskSE(), fESD(0), fOutputList(0), fPIDResponse(0), fesdTrackCuts(0x0),fesdTrackCuts_no_dca(0x0),
fTrigSel(AliVEvent::kINT7),fMultSelection(0x0),
fEventCounter(0), fEventPS(0), fEventVtx(0),fEventVtx10(0),fZVertex(0),

fTOFTimeV0MPtPosPi(0),fTOFTimeV0MPtNegPi(0),

fTOFExpTimeV0MPtPosPi_El(0),fTOFExpTimeV0MPtNegPi_El(0),
fTOFExpTimeV0MPtPosPi_Mu(0),fTOFExpTimeV0MPtNegPi_Mu(0),
fTOFExpTimeV0MPtPosPi_Pi(0),fTOFExpTimeV0MPtNegPi_Pi(0),
fTOFExpTimeV0MPtPosPi_K(0),fTOFExpTimeV0MPtNegPi_K(0),
fTOFExpTimeV0MPtPosPi_P(0),fTOFExpTimeV0MPtNegPi_P(0),
fTOFExpTimeV0MPtPosPi_D(0),fTOFExpTimeV0MPtNegPi_D(0),
fTOFMismatchTimeV0MPtPosPi(0),fTOFMismatchTimeV0MPtNegPi(0),

fTOFDCAxyV0MPtPosPi(0),fTOFDCAxyV0MPtNegPi(0),
fEventV0M(0),
fT0Resolution(0x0),fTimeOfFlightRes(0x0),fTimeOfFlightTOFRes(0x0),fTimeOfFlightGoodRes(0x0),
fTPCdEdxP(0x0),fTPCdEdxPt(0x0),
fBetaP(0x0),fBetaPNoMismatch(0x0),fBetaPNoMismatchEtaCut(0x0),
fBetaPt(0x0),fBetaPtNoMismatch(0x0),fBetaPtNoMismatchEtaCut(0x0),
fTPCTOFnSigmaPi(0),fTOFChannelVsTime(0),

fGausTime(0),fTOFGausTime(0),


fTOFTimeV0MPtPosK(0),fTOFTimeV0MPtNegK(0),

fTOFExpTimeV0MPtPosK_El(0),fTOFExpTimeV0MPtNegK_El(0),
fTOFExpTimeV0MPtPosK_Mu(0),fTOFExpTimeV0MPtNegK_Mu(0),
fTOFExpTimeV0MPtPosK_Pi(0),fTOFExpTimeV0MPtNegK_Pi(0),
fTOFExpTimeV0MPtPosK_K(0),fTOFExpTimeV0MPtNegK_K(0),
fTOFExpTimeV0MPtPosK_P(0),fTOFExpTimeV0MPtNegK_P(0),
fTOFExpTimeV0MPtPosK_D(0),fTOFExpTimeV0MPtNegK_D(0),
fTOFMismatchTimeV0MPtPosK(0),fTOFMismatchTimeV0MPtNegK(0),

fTOFDCAxyV0MPtPosK(0),fTOFDCAxyV0MPtNegK(0),
fTPCTOFnSigmaK(0),


fTOFTimeV0MPtPosP(0),fTOFTimeV0MPtNegP(0),

fTOFExpTimeV0MPtPosP_El(0),fTOFExpTimeV0MPtNegP_El(0),
fTOFExpTimeV0MPtPosP_Mu(0),fTOFExpTimeV0MPtNegP_Mu(0),
fTOFExpTimeV0MPtPosP_Pi(0),fTOFExpTimeV0MPtNegP_Pi(0),
fTOFExpTimeV0MPtPosP_K(0),fTOFExpTimeV0MPtNegP_K(0),
fTOFExpTimeV0MPtPosP_P(0),fTOFExpTimeV0MPtNegP_P(0),
fTOFExpTimeV0MPtPosP_D(0),fTOFExpTimeV0MPtNegP_D(0),
fTOFMismatchTimeV0MPtPosP(0),fTOFMismatchTimeV0MPtNegP(0),

fTOFDCAxyV0MPtPosP(0),fTOFDCAxyV0MPtNegP(0),
fTPCTOFnSigmaP(0),


fEventV0MPS(0),fEventV0MVtx(0), fSec(0),fSecondary(0),fSec_p(0),fSec_k(0),
fV0MPC(0),  ftail(0),fV0MPC_vertexcut(0),ftail_Random(0),
fPtTPC_AllP(0),fPtTOF_AllP(0), fPtTPC_AllN(0),fPtTOF_AllN(0),
fTPC_CR(0), fChi2TPCcluster(0), fDCAZ(0),fDCAxy(0),
fMinTPCcr(0),fMaxChi2PerTPC(0),fMaxDCAz(0),fMaxDCAxy(0), fSigma(value_Sigma), fSlope(value_Slope)


{}

//________________________________________________________________________
AliAnalysisTaskTOFppSpectra::AliAnalysisTaskTOFppSpectra(const char *Periodname, Int_t nTPC_CR, Int_t Chi2_TPCcluser, Int_t DCAz, Int_t DCAxy, Int_t value_Sigma, Int_t value_Slope)
  : AliAnalysisTaskSE("name"), fESD(0), fOutputList(0), fPIDResponse(0), fesdTrackCuts(0x0),fesdTrackCuts_no_dca(0x0),
fTrigSel(AliVEvent::kINT7),fMultSelection(0x0), 
fEventCounter(0), fEventPS(0), fEventVtx(0),fEventVtx10(0),fZVertex(0),

fTOFTimeV0MPtPosPi(0),fTOFTimeV0MPtNegPi(0),

fTOFExpTimeV0MPtPosPi_El(0),fTOFExpTimeV0MPtNegPi_El(0),
fTOFExpTimeV0MPtPosPi_Mu(0),fTOFExpTimeV0MPtNegPi_Mu(0),
fTOFExpTimeV0MPtPosPi_Pi(0),fTOFExpTimeV0MPtNegPi_Pi(0),
fTOFExpTimeV0MPtPosPi_K(0),fTOFExpTimeV0MPtNegPi_K(0),
fTOFExpTimeV0MPtPosPi_P(0),fTOFExpTimeV0MPtNegPi_P(0),
fTOFExpTimeV0MPtPosPi_D(0),fTOFExpTimeV0MPtNegPi_D(0),
fTOFMismatchTimeV0MPtPosPi(0),fTOFMismatchTimeV0MPtNegPi(0),

fTOFDCAxyV0MPtPosPi(0),fTOFDCAxyV0MPtNegPi(0),
fEventV0M(0),
fT0Resolution(0x0),fTimeOfFlightRes(0x0),fTimeOfFlightTOFRes(0x0),fTimeOfFlightGoodRes(0x0),
fTPCdEdxP(0x0),fTPCdEdxPt(0x0),
fBetaP(0x0),fBetaPNoMismatch(0x0),fBetaPNoMismatchEtaCut(0x0),
fBetaPt(0x0),fBetaPtNoMismatch(0x0),fBetaPtNoMismatchEtaCut(0x0),
fTPCTOFnSigmaPi(0),fTOFChannelVsTime(0),

fGausTime(0),fTOFGausTime(0),


fTOFTimeV0MPtPosK(0),fTOFTimeV0MPtNegK(0),

fTOFExpTimeV0MPtPosK_El(0),fTOFExpTimeV0MPtNegK_El(0),
fTOFExpTimeV0MPtPosK_Mu(0),fTOFExpTimeV0MPtNegK_Mu(0),
fTOFExpTimeV0MPtPosK_Pi(0),fTOFExpTimeV0MPtNegK_Pi(0),
fTOFExpTimeV0MPtPosK_K(0),fTOFExpTimeV0MPtNegK_K(0),
fTOFExpTimeV0MPtPosK_P(0),fTOFExpTimeV0MPtNegK_P(0),
fTOFExpTimeV0MPtPosK_D(0),fTOFExpTimeV0MPtNegK_D(0),
fTOFMismatchTimeV0MPtPosK(0),fTOFMismatchTimeV0MPtNegK(0),

fTOFDCAxyV0MPtPosK(0),fTOFDCAxyV0MPtNegK(0),
fTPCTOFnSigmaK(0),


fTOFTimeV0MPtPosP(0),fTOFTimeV0MPtNegP(0),

fTOFExpTimeV0MPtPosP_El(0),fTOFExpTimeV0MPtNegP_El(0),
fTOFExpTimeV0MPtPosP_Mu(0),fTOFExpTimeV0MPtNegP_Mu(0),
fTOFExpTimeV0MPtPosP_Pi(0),fTOFExpTimeV0MPtNegP_Pi(0),
fTOFExpTimeV0MPtPosP_K(0),fTOFExpTimeV0MPtNegP_K(0),
fTOFExpTimeV0MPtPosP_P(0),fTOFExpTimeV0MPtNegP_P(0),
fTOFExpTimeV0MPtPosP_D(0),fTOFExpTimeV0MPtNegP_D(0),
fTOFMismatchTimeV0MPtPosP(0),fTOFMismatchTimeV0MPtNegP(0),

fTOFDCAxyV0MPtPosP(0),fTOFDCAxyV0MPtNegP(0),
fTPCTOFnSigmaP(0),


fEventV0MPS(0),fEventV0MVtx(0), fSec(0),fSecondary(0),fSec_p(0),fSec_k(0),
fV0MPC(0),  ftail(0),fV0MPC_vertexcut(0),ftail_Random(0),
fPtTPC_AllP(0),fPtTOF_AllP(0), fPtTPC_AllN(0),fPtTOF_AllN(0),
fTPC_CR(0), fChi2TPCcluster(0), fDCAZ(0),fDCAxy(0),
fMinTPCcr(0),fMaxChi2PerTPC(0),fMaxDCAz(0),fMaxDCAxy(0), fSigma(value_Sigma), fSlope(value_Slope)



{  
// Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());

}
//________________________________________________________________________
void AliAnalysisTaskTOFppSpectra::Exit(const char *msg) {

  AliInfo(msg);
  PostData(1, fOutputList);
  return;
}
//================================================================
Double_t FuncTOFsignal(const Double_t *x, const Double_t *par)
{
  Double_t norm = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t tof_tail = par[3];
  Double_t slope = par[4];
  
	if (x[0] <= (tof_tail + mean))
    return  norm * TMath::Gaus(x[0], mean, sigma);
  else
    return norm * TMath::Gaus(tof_tail + mean, mean, sigma) * TMath::Exp(-slope * (x[0] - tof_tail - mean));//for pp
//    return  Norm * TMath::Gaus(tail + mean, mean, sigma) * TMath::Exp(-tail * (x[0] - tail - mean)/(sigma * sigma));//for pbpb
}
//________________________________________________________________________
void AliAnalysisTaskTOFppSpectra::UserCreateOutputObjects()
{

//fPPVsMultUtils=new AliPPVsMultUtils();
//AnalysisUtils = new AliAnalysisUtils();
  //pid response object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

inputHandler->SetNeedField();

  // Create histograms
  // Called once
  fOutputList = new TList();


	Double_t nPtbins=59;
	Double_t Ptbins[]={0.01, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};//PbPb bining

	Double_t nV0Mbins=11;
	Double_t V0Mbins[12] ={0.0, 0.1, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0};


	fEventCounter = new TH1F( "fEventCounter", ";Evt. Sel. Step;Count",9,0,9);
        fEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fEventCounter->GetXaxis()->SetBinLabel(2, "Pass Physics selection and trigger");
        fEventCounter->GetXaxis()->SetBinLabel(3, "INEL>0");
        fEventCounter->GetXaxis()->SetBinLabel(4, "Incomplete DAQ");
        fEventCounter->GetXaxis()->SetBinLabel(5, "SPD cluster vs tracklets cut");
        fEventCounter->GetXaxis()->SetBinLabel(6, "PileupSPD in mult bin");
        fEventCounter->GetXaxis()->SetBinLabel(7, "Reconstructed Vertex");
        fEventCounter->GetXaxis()->SetBinLabel(8, "Vertexcut<10");
        fEventCounter->GetXaxis()->SetBinLabel(9, "Selected by Analysis");
	fOutputList->Add(fEventCounter);

	fEventPS = new TH1F("fEventPS","Event after Physics selection",3,0,3);
	fOutputList->Add(fEventPS);

	fEventVtx = new TH1F("fEventVtx","Event after Phy sel & w/o vertex cut",3,0,3);
	fOutputList->Add(fEventVtx);

	fEventVtx10 = new TH1F("fEventVtx10","Event after Phy sel & w vertex cut",3,0,3);
	fOutputList->Add(fEventVtx10);

	fZVertex = new TH1F("fZVertex","Z vertex dist;Vtx_{z};Counts",40,-20,20);
	fOutputList->Add(fZVertex);

	Int_t nTimebins=2000;
//	Int_t nSigmabins=1000;
//	Int_t nResobins=1000;
	Double_t Timebins[2001];
//	Double_t Sigmabins[1001];
//	Double_t Resobins[1001];
	
	Timebins[0]=-20000.0;
//	Sigmabins[0]=-300.0;
//	Resobins[0]=0;
        for (Int_t i=1;i<2001;i++){
	Timebins[i]=Timebins[i-1]+20.0;
//	Sigmabins[i]=Sigmabins[i-1]+0.6;
//	Resobins[i]=Resobins[i-1]+1.0;
}
	
	fTOFTimeV0MPtPosPi=new TH3F("fTOFTimeV0MPtPosPi","TOF Time vs pT #pi^{+};p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFTimeV0MPtPosPi);
	fTOFTimeV0MPtNegPi=new TH3F("fTOFTimeV0MPtNegPi","TOF Time vs pT #pi^{-};p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFTimeV0MPtNegPi);



	fTOFExpTimeV0MPtPosPi_El=new TH3F("fTOFExpTimeV0MPtPosPi_El","Exp TOF Time vs pT #pi^{+} Hypo. El;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosPi_El);
	fTOFExpTimeV0MPtNegPi_El=new TH3F("fTOFExpTimeV0MPtNegPi_El","Exp TOF Time vs pT #pi^{-} Hypo. El;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegPi_El);


	fTOFExpTimeV0MPtPosPi_Mu=new TH3F("fTOFExpTimeV0MPtPosPi_Mu","Exp TOF Time vs pT #pi^{+} Hypo. Mu;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosPi_Mu);
	fTOFExpTimeV0MPtNegPi_Mu=new TH3F("fTOFExpTimeV0MPtNegPi_Mu","Exp TOF Time vs pT #pi^{-} Hypo. Mu;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegPi_Mu);


	fTOFExpTimeV0MPtPosPi_Pi=new TH3F("fTOFExpTimeV0MPtPosPi_Pi","Exp TOF Time vs pT #pi^{+} Hypo. Pi;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosPi_Pi);
	fTOFExpTimeV0MPtNegPi_Pi=new TH3F("fTOFExpTimeV0MPtNegPi_Pi","Exp TOF Time vs pT #pi^{-} Hypo. Pi;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegPi_Pi);


	fTOFExpTimeV0MPtPosPi_K=new TH3F("fTOFExpTimeV0MPtPosPi_K","Exp TOF Time vs pT #pi^{+} Hypo. K;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosPi_K);
	fTOFExpTimeV0MPtNegPi_K=new TH3F("fTOFExpTimeV0MPtNegPi_K","Exp TOF Time vs pT #pi^{-} Hypo. K;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegPi_K);


	fTOFExpTimeV0MPtPosPi_P=new TH3F("fTOFExpTimeV0MPtPosPi_P","Exp TOF Time vs pT #pi^{+} Hypo. P;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosPi_P);
	fTOFExpTimeV0MPtNegPi_P=new TH3F("fTOFExpTimeV0MPtNegPi_P","Exp TOF Time vs pT #pi^{-} Hypo. P;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegPi_P);


	fTOFExpTimeV0MPtPosPi_D=new TH3F("fTOFExpTimeV0MPtPosPi_D","Exp TOF Time vs pT #pi^{+} Hypo. D;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosPi_D);
	fTOFExpTimeV0MPtNegPi_D=new TH3F("fTOFExpTimeV0MPtNegPi_D","Exp TOF Time vs pT #pi^{-} Hypo. D;p_{T} (GeV/c);T-T_{0}-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegPi_D);



	//Mismatch time and sigma	
	fTOFMismatchTimeV0MPtPosPi=new TH3F("fTOFMismatchTimeV0MPtPosPi","TOF Mismatch Time vs pT #pi^{+};p_{T} (GeV/c);T-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFMismatchTimeV0MPtPosPi);
	fTOFMismatchTimeV0MPtNegPi=new TH3F("fTOFMismatchTimeV0MPtNegPi","TOF Mismatch Time vs pT #pi^{-};p_{T} (GeV/c);T-T_{exp #pi} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFMismatchTimeV0MPtNegPi);


	Double_t DCAbin[1201];
	DCAbin[0]=-3.00;
        for (Int_t i=1;i<1201;i++) DCAbin[i]=DCAbin[i-1]+0.005;	


	fTOFDCAxyV0MPtPosPi = new TH3F("fTOFDCAxyV0MPtPosPi","Pt vs V0M vs DCAxy #pi^{+};p_{T} (GeV/c);V0M p.c.;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fTOFDCAxyV0MPtPosPi);
	fTOFDCAxyV0MPtNegPi = new TH3F("fTOFDCAxyV0MPtNegPi","Pt vs V0M vs DCAxy #pi^{-};p_{T} (GeV/c);V0M p.c.;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fTOFDCAxyV0MPtNegPi);



	Double_t eventbins[]={0,1,2,3};
        fEventV0M = new TH2F("fEventV0M","Event vs V0M;Events;V0M percentile",3,eventbins,nV0Mbins,V0Mbins);
	fOutputList->Add(fEventV0M);



	//Quality check
	fT0Resolution = new TH1F("fT0Resolution", "T0Resolution;T0 #sigma", 500, -250, 250);
    	fOutputList->Add(fT0Resolution);
	fTimeOfFlightRes = new TH1F("fTimeOfFlightRes", "TOF Resolution in pt [0.9, 1.1];t_{TOF}-t_{0}-t_{exp #pi} (ps)", 100, -500, 500);
    	fOutputList->Add(fTimeOfFlightRes);
	fTimeOfFlightTOFRes = new TH1F("fTimeOfFlightTOFRes", "TOF Resolution with TOF T0 in pt [0.9, 1.1];t_{TOF}-t_{0}-t_{exp #pi} (ps)", 100, -500, 500);
    	fOutputList->Add(fTimeOfFlightTOFRes);
     fTimeOfFlightGoodRes = new TH1F("fTimeOfFlightGoodRes", "TOF Resolution Good in pt [0.9, 1.1];t_{TOF}-t_{0}-t_{exp #pi} (ps)",100,-500,500);
  	fTimeOfFlightGoodRes->Sumw2();
	fOutputList->Add(fTimeOfFlightGoodRes);

	
	//performance plot
	fBetaP = new TH2F("fBetaP","Distribution of the beta;p (GeV/c);TOF #beta",400,0.,10.,400,0.,1.5);
	fOutputList->Add(fBetaP);
	fBetaPNoMismatch = new TH2F("fBetaPNoMismatch","Distribution of the beta w/o Mismatch;p (GeV/c);TOF #beta",400,0.,10.,400,0.,1.5);
	fOutputList->Add(fBetaPNoMismatch);
	fBetaPNoMismatchEtaCut=new TH2F("fBetaPNoMismatchEtaCut","Distribution of beta w/o Mismatch and a |#eta|<0.5;p (GeV/c);TOF #beta",400,0.,10.,400,0.,1.5);
	fOutputList->Add(fBetaPNoMismatchEtaCut);

	fBetaPt = new TH2F("fBetaPt","Distribution of the beta;p_{T} (GeV/c);TOF #beta",400,0.,10.,400,0.,1.5);
	fOutputList->Add(fBetaPt);
	fBetaPtNoMismatch = new TH2F("fBetaPtNoMismatch","Distribution of the beta w/o Mismatch;p_{T} (GeV/c);TOF #beta",400,0.,10.,400,0.,1.5);
	fOutputList->Add(fBetaPtNoMismatch);
	fBetaPtNoMismatchEtaCut=new TH2F("fBetaPtNoMismatchEtaCut","Distribution of beta w/o Mismatch and a |#eta|<0.5;p_{T} (GeV/c);TOF #beta",400,0.,10.,400,0.,1.5);
	fOutputList->Add(fBetaPtNoMismatchEtaCut);

	fTPCdEdxP=new TH2F("fTPCdEdxP","Distribution of the TPC dE/dx;p (GeV/c);d#it{E}/d#it{x} in TPC (a. u.)",1000,0.25,30.,1000,0.,1000);
	fOutputList->Add(fTPCdEdxP);
	fTPCdEdxPt=new TH2F("fTPCdEdxPt","Distribution of the TPC dE/dx;p_{T} (GeV/c);d#it{E}/d#it{x} in TPC (a. u.)",1000,0.25,30.,1000,0.,1000);
        fOutputList->Add(fTPCdEdxPt);


	fTPCTOFnSigmaPi = new TH2F("fTPCTOFnSigmaPi","TPC TOF PID separation for pion;TOF (n#sigma);TPC (n#sigma)",200,-20.10,19.90,200,-20.10,19.90);
        fOutputList->Add(fTPCTOFnSigmaPi);

fTOFChannelVsTime = new TH2F("fTOFChannelVsTime", "TOF Raw time Vs TOF Channel;TOF Channel;TOF raw Times(ps)",500,0.,170000,7000,10000,80000);
	fOutputList->Add(fTOFChannelVsTime);

	fGausTime = new TH1F("fGausTime", "Gaussian;TOF raw Times(ps);Counts",200,-1000,1000);
	fOutputList->Add(fGausTime);
	fTOFGausTime = new TH1F("fTOFGausTime", "TOF Gaussian;TOF raw Times(ps);Counts",200,-1000,1000);
	fOutputList->Add(fTOFGausTime);


	//==================================================================Kaon=========================================

	fTOFTimeV0MPtPosK=new TH3F("fTOFTimeV0MPtPosK","TOF Time vs pT K^{+};p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFTimeV0MPtPosK);
	fTOFTimeV0MPtNegK=new TH3F("fTOFTimeV0MPtNegK","TOF Time vs pT K^{-};p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFTimeV0MPtNegK);



	fTOFExpTimeV0MPtPosK_El=new TH3F("fTOFExpTimeV0MPtPosK_El","Exp TOF Time vs pT K^{+} Hypo. El;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosK_El);
	fTOFExpTimeV0MPtNegK_El=new TH3F("fTOFExpTimeV0MPtNegK_El","Exp TOF Time vs pT K^{-} Hypo. El;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegK_El);


	fTOFExpTimeV0MPtPosK_Mu=new TH3F("fTOFExpTimeV0MPtPosK_Mu","Exp TOF Time vs pT K^{+} Hypo. Mu;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosK_Mu);
	fTOFExpTimeV0MPtNegK_Mu=new TH3F("fTOFExpTimeV0MPtNegK_Mu","Exp TOF Time vs pT K^{-} Hypo. Mu;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegK_Mu);



	fTOFExpTimeV0MPtPosK_Pi=new TH3F("fTOFExpTimeV0MPtPosK_Pi","Exp TOF Time vs pT K^{+} Hypo. Pi;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosK_Pi);
	fTOFExpTimeV0MPtNegK_Pi=new TH3F("fTOFExpTimeV0MPtNegK_Pi","Exp TOF Time vs pT K^{-} Hypo. Pi;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegK_Pi);



	fTOFExpTimeV0MPtPosK_K=new TH3F("fTOFExpTimeV0MPtPosK_K","Exp TOF Time vs pT K^{+} Hypo. K;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosK_K);
	fTOFExpTimeV0MPtNegK_K=new TH3F("fTOFExpTimeV0MPtNegK_K","Exp TOF Time vs pT K^{-} Hypo. K;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegK_K);


	fTOFExpTimeV0MPtPosK_P=new TH3F("fTOFExpTimeV0MPtPosK_P","Exp TOF Time vs pT K^{+} Hypo. P;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosK_P);
	fTOFExpTimeV0MPtNegK_P=new TH3F("fTOFExpTimeV0MPtNegK_P","Exp TOF Time vs pT K^{-} Hypo. P;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegK_P);


	fTOFExpTimeV0MPtPosK_D=new TH3F("fTOFExpTimeV0MPtPosK_D","Exp TOF Time vs pT K^{+} Hypo. D;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosK_D);
	fTOFExpTimeV0MPtNegK_D=new TH3F("fTOFExpTimeV0MPtNegK_D","Exp TOF Time vs pT K^{-} Hypo. D;p_{T} (GeV/c);T-T_{0}-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegK_D);



	//Mismatch time and sigma	
	fTOFMismatchTimeV0MPtPosK=new TH3F("fTOFMismatchTimeV0MPtPosK","TOF Mismatch Time vs pT K^{+};p_{T} (GeV/c);T-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFMismatchTimeV0MPtPosK);
	fTOFMismatchTimeV0MPtNegK=new TH3F("fTOFMismatchTimeV0MPtNegK","TOF Mismatch Time vs pT K^{-};p_{T} (GeV/c);T-T_{exp K} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFMismatchTimeV0MPtNegK);


	fTOFDCAxyV0MPtPosK = new TH3F("fTOFDCAxyV0MPtPosK","Pt vs V0M vs DCAxy K^{+};p_{T} (GeV/c);V0M p.c.;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fTOFDCAxyV0MPtPosK);
	fTOFDCAxyV0MPtNegK = new TH3F("fTOFDCAxyV0MPtNegK","Pt vs V0M vs DCAxy K^{-};p_{T} (GeV/c);V0M p.c.;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fTOFDCAxyV0MPtNegK);


	fTPCTOFnSigmaK = new TH2F("fTPCTOFnSigmaK","TPC TOF separation for Kaon;TOF (n#sigma);TPC (n#sigma)",200,-20.10,19.90,200,-20.10,19.90);
        fOutputList->Add(fTPCTOFnSigmaK);

//==========================================================Proton====================================================================


	fTOFTimeV0MPtPosP=new TH3F("fTOFTimeV0MPtPosP","TOF Time vs pT P^{+};p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFTimeV0MPtPosP);
	fTOFTimeV0MPtNegP=new TH3F("fTOFTimeV0MPtNegP","TOF Time vs pT P^{-};p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFTimeV0MPtNegP);



	fTOFExpTimeV0MPtPosP_El=new TH3F("fTOFExpTimeV0MPtPosP_El","Exp TOF Time vs pT P^{+} Hypo. El;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosP_El);
	fTOFExpTimeV0MPtNegP_El=new TH3F("fTOFExpTimeV0MPtNegP_El","Exp TOF Time vs pT P^{-} Hypo. El;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegP_El);


	fTOFExpTimeV0MPtPosP_Mu=new TH3F("fTOFExpTimeV0MPtPosP_Mu","Exp TOF Time vs pT P^{+} Hypo. Mu;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosP_Mu);
	fTOFExpTimeV0MPtNegP_Mu=new TH3F("fTOFExpTimeV0MPtNegP_Mu","Exp TOF Time vs pT P^{-} Hypo. Mu;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegP_Mu);



	fTOFExpTimeV0MPtPosP_Pi=new TH3F("fTOFExpTimeV0MPtPosP_Pi","Exp TOF Time vs pT P^{+} Hypo. P;p_{T} (GeV/c);T-T_{0}-T_{exp #P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosP_Pi);
	fTOFExpTimeV0MPtNegP_Pi=new TH3F("fTOFExpTimeV0MPtNegP_Pi","Exp TOF Time vs pT P^{-} Hypo. P;p_{T} (GeV/c);T-T_{0}-T_{exp #P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegP_Pi);


	fTOFExpTimeV0MPtPosP_K=new TH3F("fTOFExpTimeV0MPtPosP_K","Exp TOF Time vs pT P^{+} Hypo. K;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosP_K);
	fTOFExpTimeV0MPtNegP_K=new TH3F("fTOFExpTimeV0MPtNegP_K","Exp TOF Time vs pT P^{-} Hypo. K;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegP_K);


	fTOFExpTimeV0MPtPosP_P=new TH3F("fTOFExpTimeV0MPtPosP_P","Exp TOF Time vs pT P^{+} Hypo. P;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosP_P);
	fTOFExpTimeV0MPtNegP_P=new TH3F("fTOFExpTimeV0MPtNegP_P","Exp TOF Time vs pT P^{-} Hypo. P;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegP_P);


	fTOFExpTimeV0MPtPosP_D=new TH3F("fTOFExpTimeV0MPtPosP_D","Exp TOF Time vs pT P^{+} Hypo. D;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtPosP_D);
	fTOFExpTimeV0MPtNegP_D=new TH3F("fTOFExpTimeV0MPtNegP_D","Exp TOF Time vs pT P^{-} Hypo. D;p_{T} (GeV/c);T-T_{0}-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
	fOutputList->Add(fTOFExpTimeV0MPtNegP_D);



	//Mismatch time and sigma	
	fTOFMismatchTimeV0MPtPosP=new TH3F("fTOFMismatchTimeV0MPtPosP","TOF Mismatch Time vs pT P^{+};p_{T} (GeV/c);T-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFMismatchTimeV0MPtPosP);
	fTOFMismatchTimeV0MPtNegP=new TH3F("fTOFMismatchTimeV0MPtNegP","TOF Mismatch Time vs pT P^{-};p_{T} (GeV/c);T-T_{exp P} (ps);V0M PC",nPtbins,Ptbins,nTimebins,Timebins,nV0Mbins,V0Mbins);
        fOutputList->Add(fTOFMismatchTimeV0MPtNegP);


	fTOFDCAxyV0MPtPosP = new TH3F("fTOFDCAxyV0MPtPosP","Pt vs V0M vs DCAxy P^{+};p_{T} (GeV/c);V0M p.c.;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fTOFDCAxyV0MPtPosP);
	fTOFDCAxyV0MPtNegP = new TH3F("fTOFDCAxyV0MPtNegP","Pt vs V0M vs DCAxy P^{-};p_{T} (GeV/c);V0M p.c.;DCA_{xy}",nPtbins,Ptbins,nV0Mbins,V0Mbins,1200,DCAbin);
        fOutputList->Add(fTOFDCAxyV0MPtNegP);


	fTPCTOFnSigmaP = new TH2F("fTPCTOFnSigmaP","TPC TOF separation for Proton;TOF (n#sigma);TPC (n#sigma)",200,-20.10,19.90,200,-20.10,19.90);
        fOutputList->Add(fTPCTOFnSigmaP);


        fEventV0MPS = new TH2F("fEventV0MPS","Event vs V0M after PS;Events;V0M percentile",3,eventbins,nV0Mbins,V0Mbins);
	fOutputList->Add(fEventV0MPS);
        fEventV0MVtx = new TH2F("fEventV0MVtx","Event vs V0M;Events;V0M percentile",3,eventbins,nV0Mbins,V0Mbins);
	fOutputList->Add(fEventV0MVtx);


	fV0MPC=new TH1F("fV0MPC","V0M percentile;V0M PC",120,0,120);
	fOutputList->Add(fV0MPC);
	
	ftail=new TH1F("ftail","Seco tail",300, -2000,10000);
	fOutputList->Add(ftail);
	
	fV0MPC_vertexcut=new TH1F("fV0MPC_vertexcut","V0M percentile;V0M PC",120,0,120);
	//fV0mpc=new TH1F("fV0mpc","V0M percentile;V0M PC",nV0Mbins,V0Mbins);
	fOutputList->Add(fV0MPC_vertexcut);
	
	fPtTPC_AllP=new TH1F("fPtTPC_AllP","TPC pt distribution;p_{T} (GeV/c)",nPtbins,Ptbins);
	fOutputList->Add(fPtTPC_AllP);
	fPtTOF_AllP=new TH1F("fPtTOF_AllP","TOF pt distribution;p_{T} (GeV/c)",nPtbins,Ptbins);
	fOutputList->Add(fPtTOF_AllP);
	fPtTPC_AllN=new TH1F("fPtTPC_AllN","TPC pt distribution;p_{T} (GeV/c)",nPtbins,Ptbins);
	fOutputList->Add(fPtTPC_AllN);
	fPtTOF_AllN=new TH1F("fPtTOF_AllN","TOF pt distribution;p_{T} (GeV/c)",nPtbins,Ptbins);
	fOutputList->Add(fPtTOF_AllN);

	ftail_Random=new TH1F("ftail_Random","Seco tail random",300, -2000,10000);
	//ftail_Random=new TH1F("ftail_Random","Seco tail random",600, -2000,10000);
	fOutputList->Add(ftail_Random);

	fTPC_CR=new TH1F("fTPC_CR","TPC cr distribution;# of CR",200,0,200);
        fOutputList->Add(fTPC_CR);
        fChi2TPCcluster=new TH1F("fChi2TPCcluster","Chi2 /TPC cluster distribution;# of CR",100,0,10);
        fOutputList->Add(fChi2TPCcluster);
        fDCAZ=new TH1F("fDCAZ","DCA z distribution;# of CR",100,0,10);
        fOutputList->Add(fDCAZ);
        fDCAxy=new TH1F("fDCAxy","DCA xy distribution;# of CR",1200,-3,3);
        fOutputList->Add(fDCAxy);


	//ESD track cut
        fesdTrackCuts =  AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
	fesdTrackCuts->SetMinNCrossedRowsTPC(fMinTPCcr);
        fesdTrackCuts->SetMaxChi2PerClusterTPC(fMaxChi2PerTPC);
        fesdTrackCuts->SetMaxDCAToVertexZ(fMaxDCAz);
        fesdTrackCuts->SetMaxDCAToVertexXYPtDep(fMaxDCAxy);

        //no DCAxy cut for secondaries
        fesdTrackCuts_no_dca =  AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);// no DCA xy cut
        fesdTrackCuts_no_dca->SetMinNCrossedRowsTPC(fMinTPCcr);
        fesdTrackCuts_no_dca->SetMaxChi2PerClusterTPC(fMaxChi2PerTPC);
        fesdTrackCuts_no_dca->SetMaxDCAToVertexZ(fMaxDCAz);

 PostData(1, fOutputList);
}
//________________________________________________________________________
void AliAnalysisTaskTOFppSpectra::UserExec(Option_t *)
{


	Double_t y_value[301]={0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 2, 0, 3, 2, 2, 6, 4, 4, 9, 12, 21, 40, 69, 103, 188, 296, 410, 547, 556, 605, 739, 1015, 1277, 1677, 2581, 4393, 6739, 13724, 55179, 271841, 410963, 121659, 17730, 6138, 3521, 2488, 1898, 1439, 1238, 1048, 911, 786, 775, 680, 649, 567, 507, 481, 409, 425, 372, 355, 318, 280, 269, 221, 227, 207, 184, 191, 179, 134, 167, 143, 146, 142, 138, 115, 104, 112, 101, 115, 119, 95, 102, 110, 74, 82, 85, 70, 78, 76, 73, 71, 60, 53, 67, 54, 58, 58, 47, 52, 38, 51, 46, 40, 47, 48, 36, 32, 38, 23, 33, 28, 31, 30, 34, 34, 32, 30, 21, 20, 26, 23, 23, 30, 20, 14, 24, 18, 17, 13, 17, 19, 13, 28, 16, 16, 16, 10, 12, 7, 17, 9, 13, 6, 12, 20, 12, 15, 15, 14, 12, 11, 15, 15, 6, 10, 6, 6, 8, 5, 11, 4, 10, 9, 6, 9, 6, 7, 10, 6, 6, 4, 5, 5, 2, 3, 7, 6, 7, 5, 5, 5, 8, 4, 6, 7, 4, 4, 6, 5, 4, 8, 7, 5, 2, 4, 6, 1, 4, 2, 4, 4, 3, 1, 6, 3, 8, 3, 4, 5, 4, 5, 2, 3, 5, 3, 3, 2, 2, 8, 5, 2, 2, 6, 4, 3, 2, 4, 1, 5, 2, 0, 1, 6, 3, 7, 2, 3, 4, 1, 1, 4, 3, 2, 1, 2, 3, 0, 3, 2, 2, 4, 2, 2, 2, 3, 1, 0, 3, 3, 3, 3, 2, 4, 3, 6, 2, 3, 2, 1, 2, 1, 2, 5, 2, 1, 2, 4, 1, 2, 3, 0, 1, 1, 2, 0, 4, 1};

 
        for (Int_t i=1;i<=300;i++) {
        //for (Int_t i=1;i<=600;i++) {
       ftail->SetBinContent(i,y_value[i]);
}

// Main loop
  // Called for each event
  // Post output data.
        AliESDEvent *fESD = 0x0;
          fESD = dynamic_cast<AliESDEvent*>(InputEvent());
          if (!fESD) {
            printf("ERROR: fESD not available\n");
            return;
          }


        AliESDVZERO *esdV0=fESD->GetVZEROData();

        fEventCounter->Fill(0.5);

//================================================physics selection=================================
	UInt_t maskPhysSel = ((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  	maskPhysSel &= fTrigSel;
  	TString firedTriggerClasses = fESD->GetFiredTriggerClasses();
	if (maskPhysSel != fTrigSel) 
	  return Exit(Form("Event doesn't pass physics evt. sel. for trigger %d", fTrigSel));
	fEventCounter->Fill(1.5);
	
/*	if(!firedTriggerClasses.Contains("ALLNOTRD"))
	  return;
*/
	//Bool_t INELgt0 = (AliMultSelection*) fESD->GetThisEventINELgtZERO();
	//INEL >0 cut
	Double_t INELgt0=(AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0) >= 1);
	if(!INELgt0)
	   return;
	fEventCounter->Fill(2.5);

	//incomplete DAQ
	if (fESD->IsIncompleteDAQ())
	  return;
	fEventCounter->Fill(3.5);

	//tracklet vs cluster cut
	AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
	Double_t IsCluVstrk = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD);
	if(IsCluVstrk)
	  return;	
	fEventCounter->Fill(4.5);
      
	//pileup (type-1)
        Bool_t ispileup = fESD->IsPileupFromSPDInMultBins();
        if(ispileup) 
	  return Exit(Form("Pile up does't pass"));
	fEventCounter->Fill(5.5);
	
	
	fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
        if (!fMultSelection)
           //cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
           AliWarning("No AliMultSelection Object Found --------");
        else
          V0mpc = fMultSelection->GetMultiplicityPercentile("V0M",kFALSE);

	fV0MPC->Fill(V0mpc);
	fEventPS->Fill(1);
	fEventV0MPS->Fill(1,V0mpc);

        //z vertex cut<10
        Bool_t IsVertex = selectVertex2015pp(fESD,kTRUE,kFALSE,kTRUE);
 	if(!IsVertex)
	  return Exit(Form("Not in vertex cut"));
	
	
	fEventCounter->Fill(6.5);
	fEventVtx->Fill(1);
	fEventV0MVtx->Fill(1,V0mpc);
	fV0MPC_vertexcut->Fill(V0mpc);
	
	const AliESDVertex * fVertex = fESD->GetPrimaryVertex();
	if (TMath::Abs(fVertex->GetZ())>10) return;
	
	fEventCounter->Fill(7.5);
	fEventVtx10->Fill(1);
	fEventV0M->Fill(1,V0mpc);

	fZVertex->Fill(fVertex->GetZ());


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++track loop++++++++++++++++++++++++++++++
	Double_t N_ch=0;


	fPIDResponse->SetTOFResponse(fESD,AliPIDResponse::kBest_T0);

	 for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
        AliESDtrack *track= fESD->GetTrack(iTracks);
                if (!track) {
                 printf("ERROR: Could not receive track %d\n", iTracks);
         continue;
}

	Bool_t TPCPIDStatus = TPCPID(track);
        Bool_t TOFPIDStatus = TOFPID(track);
	
	Float_t dcaxy[2];
        Float_t dcaz[3];
        track->GetImpactParameters(dcaxy,dcaz);


        if(fesdTrackCuts->AcceptTrack(track)){

	fTPC_CR->Fill(track->GetTPCCrossedRows());
        fChi2TPCcluster->Fill(track->GetTPCchi2()/track->GetTPCNcls());
        fDCAZ->Fill(dcaxy[1]);
        fDCAxy->Fill(dcaxy[0]);


	
	if(TMath::Abs(track->Eta())<0.8){
	if(TPCPIDStatus){
	if(track->Charge()>0.) fPtTPC_AllP->Fill(track->Pt());
	if(track->Charge()<0.) fPtTPC_AllN->Fill(track->Pt());
}	
	if(TOFPIDStatus){
	if(track->Charge()>0.) fPtTOF_AllP->Fill(track->Pt());
	if(track->Charge()<0.) fPtTOF_AllN->Fill(track->Pt());
}	
}

	if(TOFPIDStatus){

	Double_t pt=track->Pt();
        Double_t eta=track->Eta();


	//performance plot
	Double_t nTOFClusters = track->GetNTOFclusters();
	Double_t trkLength = track->GetIntegratedLength();
	const Double_t C_Value = TMath::C() * 1.e2 / 1.e12; // cm/ps

	 const Double_t beta = trkLength / ((track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(track->P())) * C_Value);
        fBetaP->Fill(track->P(), beta);
        fBetaPt->Fill(track->Pt(), beta);
        if(nTOFClusters < 2){
        fBetaPNoMismatch->Fill(track->P(),beta);
        if(TMath::Abs(eta) < 0.5) fBetaPNoMismatchEtaCut->Fill(track->P(),beta);
          
	fBetaPtNoMismatch->Fill(track->Pt(),beta);
        if(TMath::Abs(eta) < 0.5) fBetaPtNoMismatchEtaCut->Fill(track->Pt(),beta);
        }

	//TPC
        fTPCdEdxP->Fill(track->P(),track->GetTPCsignal());
        fTPCdEdxPt->Fill(track->Pt(),track->GetTPCsignal());
	
	//TOF channel
	 Double_t fTOFchan = track->GetTOFCalChannel();    // Channel Index of the TOF Signal
	 fTOFChannelVsTime->Fill(fTOFchan,track->GetTOFsignal());



	if(TMath::Abs(eta)<0.8){


	Double_t Y_Pi=Rapidity(track,AliPID::ParticleMass(AliPID::kPion));
	Double_t Y_K=Rapidity(track,AliPID::ParticleMass(AliPID::kKaon));
	Double_t Y_P=Rapidity(track,AliPID::ParticleMass(AliPID::kProton));


	
	Double_t TOFImpactDx = track->GetTOFsignalDx();  // local x  of track's impact on the TOF pad
	Double_t TOFImpactDz = track->GetTOFsignalDz();  // local z  of track's impact on the TOF pad
    	Double_t T0TrkSigma = fPIDResponse->GetTOFResponse().GetStartTimeRes(track->P());  // T0best resolution time
	fT0Resolution->Fill(T0TrkSigma);


	//PID with the TOF
  	Double_t fTOFTime = track->GetTOFsignal(); //Gets the TOF signal
        Double_t fT0TrkTime = fPIDResponse->GetTOFResponse().GetStartTime(track->P()); // T0best time

	Double_t pidTime[6]; track->GetIntegratedTimes(pidTime,6);
        //Double_t TExpTimeEl=fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron); or
        Double_t fTOFExpTimeEl = pidTime[AliPID::kElectron];
        Double_t fTOFExpTimeMu = pidTime[AliPID::kMuon];
        Double_t fTOFExpTimePi = pidTime[AliPID::kPion];
        Double_t fTOFExpTimeK = pidTime[AliPID::kKaon];
        Double_t fTOFExpTimeP = pidTime[AliPID::kProton];
        Double_t fTOFExpTimeD = pidTime[AliPID::kDeuteron];
        //cout<<TTOFExpTimeK<<"            "<<pidTime[3]<<endl;




	//TOF sigma Expected
	Double_t fTOFExpSigmaEl = fPIDResponse->GetTOFResponse().GetExpectedSigma(track->P(), pidTime[AliPID::kElectron], AliPID::kElectron);
	Double_t fTOFExpSigmaMu = fPIDResponse->GetTOFResponse().GetExpectedSigma(track->P(), pidTime[AliPID::kMuon], AliPID::kMuon);
	Double_t fTOFExpSigmaPi = fPIDResponse->GetTOFResponse().GetExpectedSigma(track->P(), pidTime[AliPID::kPion], AliPID::kPion);
	Double_t fTOFExpSigmaK = fPIDResponse->GetTOFResponse().GetExpectedSigma(track->P(), pidTime[AliPID::kKaon], AliPID::kKaon);
	Double_t fTOFExpSigmaP = fPIDResponse->GetTOFResponse().GetExpectedSigma(track->P(), pidTime[AliPID::kProton], AliPID::kProton);
	Double_t fTOFExpSigmaD = fPIDResponse->GetTOFResponse().GetExpectedSigma(track->P(), pidTime[AliPID::kDeuteron], AliPID::kDeuteron);


	//TOF nsigma
	Double_t fTOFSigmaEl = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
	Double_t fTOFSigmaMu = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kMuon);
	Double_t fTOFSigmaPi = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
	Double_t fTOFSigmaK = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
	Double_t fTOFSigmaP = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
	Double_t fTOFSigmaD = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);

	//TPC nsigma
	Double_t fTPCSigmaEl = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
	Double_t fTPCSigmaMu = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kMuon);
	Double_t fTPCSigmaPi = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
	Double_t fTPCSigmaK = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
	Double_t fTPCSigmaP = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
	Double_t fTPCSigmaD = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
	

	//TPC TOF Combined
	fTPCTOFnSigmaPi->Fill(fTOFSigmaPi,fTPCSigmaPi);
	fTPCTOFnSigmaK->Fill(fTOFSigmaK,fTPCSigmaK);
	fTPCTOFnSigmaP->Fill(fTOFSigmaP,fTPCSigmaP);


	if((track->Pt()>1.0) && (track->Pt()<1.1)){//P range selected for TOF resolution measurements//previously P()
        Float_t deltatime = fTOFTime-fT0TrkTime-fTOFExpTimePi;
        fTimeOfFlightRes->Fill(deltatime);
        if(fT0TrkTime != 0) fTimeOfFlightTOFRes->Fill(deltatime);
        if((TMath::Abs(TOFImpactDx)<1.25) && (TMath::Abs(TOFImpactDz)<1.75)) fTimeOfFlightGoodRes->Fill(deltatime);
      }
	
	const Double_t TdiffPi = fTOFTime-fT0TrkTime-fTOFExpTimePi;
	const Double_t TdiffK = fTOFTime-fT0TrkTime-fTOFExpTimeK;
	const Double_t TdiffP = fTOFTime-fT0TrkTime-fTOFExpTimeP;
  	

	TF1 *fTOFsignal = new TF1("fTOFsignal", FuncTOFsignal, -3000., 3000., 5);//changes to 5 parameter on 10 apr
	fTOFsignal->SetParameter(0, 1);
	fTOFsignal->SetParameter(1, 0);
	fTOFsignal->SetParameter(2, fSigma);
	fTOFsignal->SetParameter(3, 95);

	if (fSlope==0)	fTOFsignal->SetParameter(4, 0.0125);
	else if (fSlope==+1)	fTOFsignal->SetParameter(4, 0.01375);
	else if (fSlope==-1)	fTOFsignal->SetParameter(4, 0.01125);

	Double_t tof_sig=fTOFsignal->GetRandom();
	Double_t sec_tail=ftail->GetRandom()-20;
	ftail_Random->Fill(sec_tail);

	fTOFGausTime->Fill(tof_sig);
	
	Double_t sigma_el=TMath::Sqrt(fTOFExpSigmaEl*fTOFExpSigmaEl-80*80);
	Double_t sigma_mu=TMath::Sqrt(fTOFExpSigmaMu*fTOFExpSigmaMu-80*80);
	Double_t sigma_pi=TMath::Sqrt(fTOFExpSigmaPi*fTOFExpSigmaPi-80*80);
	Double_t sigma_k=TMath::Sqrt(fTOFExpSigmaK*fTOFExpSigmaK-80*80);
	Double_t sigma_p=TMath::Sqrt(fTOFExpSigmaP*fTOFExpSigmaP-80*80);
	Double_t sigma_d=TMath::Sqrt(fTOFExpSigmaD*fTOFExpSigmaD-80*80);
	Double_t extra_sm_el= gRandom->Gaus(0, sigma_el);
	Double_t extra_sm_mu= gRandom->Gaus(0, sigma_mu);
	Double_t extra_sm_pi= gRandom->Gaus(0, sigma_pi);
	Double_t extra_sm_k= gRandom->Gaus(0, sigma_k);
	Double_t extra_sm_p= gRandom->Gaus(0, sigma_p);
	Double_t extra_sm_d= gRandom->Gaus(0, sigma_d);

       // for pion 
        if(TMath::Abs(Y_Pi)<0.5){
	
	//TOF Time
	if(track->Charge()>0)fTOFTimeV0MPtPosPi->Fill(track->Pt(),TdiffPi,V0mpc);
	if(track->Charge()<0)fTOFTimeV0MPtNegPi->Fill(track->Pt(),TdiffPi,V0mpc);

/*	const Double_t TdiffSigmaPi = TdiffPi/fTOFExpSigmaPi;
	//TOF Sigma
	fTOFSigmaV0MPtPi->Fill(track->Pt(),TdiffSigmaPi,V0mpc);
	if(track->Charge()>0)fTOFSigmaV0MPtPosPi->Fill(track->Pt(),TdiffSigmaPi,V0mpc);
	if(track->Charge()<0)fTOFSigmaV0MPtNegPi->Fill(track->Pt(),TdiffSigmaPi,V0mpc);
*/


//	Double_t gaus= gRandom->Gaus(0., 95.);
//	fGausTime->Fill(gaus);
	

	Double_t pi_expTdiffEl=0,pi_expTdiffMu=0,pi_expTdiffPi=0,pi_expTdiffK=0,pi_expTdiffP=0,pi_expTdiffD=0;
	
	pi_expTdiffEl=fTOFExpTimeEl-fTOFExpTimePi+tof_sig+extra_sm_el;
	pi_expTdiffMu=fTOFExpTimeMu-fTOFExpTimePi+tof_sig+extra_sm_mu;
	pi_expTdiffPi=fTOFExpTimePi-fTOFExpTimePi+tof_sig+extra_sm_pi+sec_tail;
	pi_expTdiffK=fTOFExpTimeK-fTOFExpTimePi+tof_sig+extra_sm_k+sec_tail;
	pi_expTdiffP=fTOFExpTimeP-fTOFExpTimePi+tof_sig+extra_sm_p+sec_tail;
	pi_expTdiffD=fTOFExpTimeD-fTOFExpTimePi+tof_sig+extra_sm_d;

	
	//Expected time with different particles hypothesis
	//electron
	//fTOFExpTimeV0MPtPi_El->Fill(track->Pt(),(fTOFExpTimePi-fT0TrkTime-fTOFExpTimeEl),V0mpc);//test1. not matching

	if(track->Charge()>0) fTOFExpTimeV0MPtPosPi_El->Fill(track->Pt(),pi_expTdiffEl,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegPi_El->Fill(track->Pt(),pi_expTdiffEl,V0mpc);


	//muon
	if(track->Charge()>0) fTOFExpTimeV0MPtPosPi_Mu->Fill(track->Pt(),pi_expTdiffMu,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegPi_Mu->Fill(track->Pt(),pi_expTdiffMu,V0mpc);

	//Pion
	if(track->Charge()>0) fTOFExpTimeV0MPtPosPi_Pi->Fill(track->Pt(),(pi_expTdiffPi),V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegPi_Pi->Fill(track->Pt(),(pi_expTdiffPi),V0mpc);
	
	//Kaon
	if(track->Charge()>0) fTOFExpTimeV0MPtPosPi_K->Fill(track->Pt(),pi_expTdiffK,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegPi_K->Fill(track->Pt(),pi_expTdiffK,V0mpc);

	//Proton
	if(track->Charge()>0) fTOFExpTimeV0MPtPosPi_P->Fill(track->Pt(),pi_expTdiffP,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegPi_P->Fill(track->Pt(),pi_expTdiffP,V0mpc);

	//Deuteron
	if(track->Charge()>0) fTOFExpTimeV0MPtPosPi_D->Fill(track->Pt(),pi_expTdiffD,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegPi_D->Fill(track->Pt(),pi_expTdiffD,V0mpc);


	Double_t mismtch_Pi=AliTOFPIDResponse::GetMismatchRandomValue(track->Eta());

	//TOF Mismatch Time
	if(track->Charge()>0) fTOFMismatchTimeV0MPtPosPi->Fill(track->Pt(),(mismtch_Pi-fTOFExpTimePi),V0mpc);
	if(track->Charge()<0) fTOFMismatchTimeV0MPtNegPi->Fill(track->Pt(),(mismtch_Pi-fTOFExpTimePi),V0mpc);

}//rapidity pion

//=============================================================================Kaon=========================================================

	// for Kaon 
        if(TMath::Abs(Y_K)<0.5){


	//TOF Time
	if(track->Charge()>0)fTOFTimeV0MPtPosK->Fill(track->Pt(),TdiffK,V0mpc);
	if(track->Charge()<0)fTOFTimeV0MPtNegK->Fill(track->Pt(),TdiffK,V0mpc);

/*	const Double_t TdiffSigmaK = TdiffK/fTOFExpSigmaK;
	//TOF Sigma
	fTOFSigmaV0MPtK->Fill(track->Pt(),TdiffSigmaK,V0mpc);
	if(track->Charge()>0)fTOFSigmaV0MPtPosK->Fill(track->Pt(),TdiffSigmaK,V0mpc);
	if(track->Charge()<0)fTOFSigmaV0MPtNegK->Fill(track->Pt(),TdiffSigmaK,V0mpc);
*/


/*	Double_t K_sigma_el=TMath::Sqrt(fTOFExpSigmaEl*fTOFExpSigmaEl-80*80);
	Double_t K_sigma_mu=TMath::Sqrt(fTOFExpSigmaMu*fTOFExpSigmaMu-80*80);
	Double_t K_sigma_pi=TMath::Sqrt(fTOFExpSigmaPi*fTOFExpSigmaPi-80*80);
	Double_t K_sigma_k=TMath::Sqrt(fTOFExpSigmaK*fTOFExpSigmaK-80*80);
	Double_t K_sigma_p=TMath::Sqrt(fTOFExpSigmaP*fTOFExpSigmaP-80*80);
	Double_t K_sigma_d=TMath::Sqrt(fTOFExpSigmaD*fTOFExpSigmaD-80*80);
	Double_t K_extra_sm_el= gRandom->Gaus(0, K_sigma_el);
	Double_t K_extra_sm_mu= gRandom->Gaus(0, K_sigma_mu);
	Double_t K_extra_sm_pi= gRandom->Gaus(0, K_sigma_pi);
	Double_t K_extra_sm_k= gRandom->Gaus(0, K_sigma_k);
	Double_t K_extra_sm_p= gRandom->Gaus(0, K_sigma_p);
	Double_t K_extra_sm_d= gRandom->Gaus(0, K_sigma_d);
*/

	Double_t K_expTdiffEl=0,K_expTdiffMu=0,K_expTdiffPi=0,K_expTdiffK=0,K_expTdiffP=0,K_expTdiffD=0;
	
	K_expTdiffEl=fTOFExpTimeEl-fTOFExpTimeK+tof_sig+extra_sm_el;
	K_expTdiffMu=fTOFExpTimeMu-fTOFExpTimeK+tof_sig+extra_sm_mu;
	K_expTdiffPi=fTOFExpTimePi-fTOFExpTimeK+tof_sig+extra_sm_pi+sec_tail;
	K_expTdiffK=fTOFExpTimeK-fTOFExpTimeK+tof_sig+extra_sm_k+sec_tail;
	K_expTdiffP=fTOFExpTimeP-fTOFExpTimeK+tof_sig+extra_sm_p+sec_tail;
	K_expTdiffD=fTOFExpTimeD-fTOFExpTimeK+tof_sig+extra_sm_d;

	
	//Expected time with different particles hypothesis
	//electron
	if(track->Charge()>0) fTOFExpTimeV0MPtPosK_El->Fill(track->Pt(),K_expTdiffEl,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegK_El->Fill(track->Pt(),K_expTdiffEl,V0mpc);

	//muon
	if(track->Charge()>0) fTOFExpTimeV0MPtPosK_Mu->Fill(track->Pt(),K_expTdiffMu,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegK_Mu->Fill(track->Pt(),K_expTdiffMu,V0mpc);

	//Pion
	if(track->Charge()>0) fTOFExpTimeV0MPtPosK_Pi->Fill(track->Pt(),(K_expTdiffPi),V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegK_Pi->Fill(track->Pt(),(K_expTdiffPi),V0mpc);
	
	//Kaon
	if(track->Charge()>0) fTOFExpTimeV0MPtPosK_K->Fill(track->Pt(),K_expTdiffK,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegK_K->Fill(track->Pt(),K_expTdiffK,V0mpc);

	//Proton
	if(track->Charge()>0) fTOFExpTimeV0MPtPosK_P->Fill(track->Pt(),K_expTdiffP,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegK_P->Fill(track->Pt(),K_expTdiffP,V0mpc);

	//Deuteron
	if(track->Charge()>0) fTOFExpTimeV0MPtPosK_D->Fill(track->Pt(),K_expTdiffD,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegK_D->Fill(track->Pt(),K_expTdiffD,V0mpc);

	Double_t mismtch_K=AliTOFPIDResponse::GetMismatchRandomValue(track->Eta());

	//TOF Mismatch Time
	if(track->Charge()>0) fTOFMismatchTimeV0MPtPosK->Fill(track->Pt(),(mismtch_K-fTOFExpTimeK),V0mpc);
	if(track->Charge()<0) fTOFMismatchTimeV0MPtNegK->Fill(track->Pt(),(mismtch_K-fTOFExpTimeK),V0mpc);

}//rapidity Kaon
//==============================================================Proton=======================================================


	// for Proton 
        if(TMath::Abs(Y_P)<0.5){
	
	//TOF Time
	if(track->Charge()>0)fTOFTimeV0MPtPosP->Fill(track->Pt(),TdiffP,V0mpc);
	if(track->Charge()<0)fTOFTimeV0MPtNegP->Fill(track->Pt(),TdiffP,V0mpc);

/*	const Double_t TdiffSigmaP = TdiffP/fTOFExpSigmaP;
	//TOF Sigma
	fTOFSigmaV0MPtP->Fill(track->Pt(),TdiffSigmaP,V0mpc);
	if(track->Charge()>0)fTOFSigmaV0MPtPosP->Fill(track->Pt(),TdiffSigmaP,V0mpc);
	if(track->Charge()<0)fTOFSigmaV0MPtNegP->Fill(track->Pt(),TdiffSigmaP,V0mpc);
*/

/*	Double_t P_sigma_el=TMath::Sqrt(fTOFExpSigmaEl*fTOFExpSigmaEl-80*80);
	Double_t P_sigma_mu=TMath::Sqrt(fTOFExpSigmaMu*fTOFExpSigmaMu-80*80);
	Double_t P_sigma_pi=TMath::Sqrt(fTOFExpSigmaPi*fTOFExpSigmaPi-80*80);
	Double_t P_sigma_k=TMath::Sqrt(fTOFExpSigmaK*fTOFExpSigmaK-80*80);
	Double_t P_sigma_p=TMath::Sqrt(fTOFExpSigmaP*fTOFExpSigmaP-80*80);
	Double_t P_sigma_d=TMath::Sqrt(fTOFExpSigmaD*fTOFExpSigmaD-80*80);
	Double_t P_extra_sm_el= gRandom->Gaus(0, P_sigma_el);
	Double_t P_extra_sm_mu= gRandom->Gaus(0, P_sigma_mu);
	Double_t P_extra_sm_pi= gRandom->Gaus(0, P_sigma_pi);
	Double_t P_extra_sm_k= gRandom->Gaus(0, P_sigma_k);
	Double_t P_extra_sm_p= gRandom->Gaus(0, P_sigma_p);
	Double_t P_extra_sm_d= gRandom->Gaus(0, P_sigma_d);
*/

	Double_t P_expTdiffEl=0,P_expTdiffMu=0,P_expTdiffPi=0,P_expTdiffK=0,P_expTdiffP=0,P_expTdiffD=0;
	
	P_expTdiffEl=fTOFExpTimeEl-fTOFExpTimeP+tof_sig+extra_sm_el;
	P_expTdiffMu=fTOFExpTimeMu-fTOFExpTimeP+tof_sig+extra_sm_mu;
	P_expTdiffPi=fTOFExpTimePi-fTOFExpTimeP+tof_sig+extra_sm_pi+sec_tail;
	P_expTdiffK=fTOFExpTimeK-fTOFExpTimeP+tof_sig+extra_sm_k+sec_tail;
	P_expTdiffP=fTOFExpTimeP-fTOFExpTimeP+tof_sig+extra_sm_p+sec_tail;
	P_expTdiffD=fTOFExpTimeD-fTOFExpTimeP+tof_sig+extra_sm_d;

	
	//Expected time with different particles hypothesis
	//electron
	if(track->Charge()>0) fTOFExpTimeV0MPtPosP_El->Fill(track->Pt(),P_expTdiffEl,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegP_El->Fill(track->Pt(),P_expTdiffEl,V0mpc);

	//muon
	if(track->Charge()>0) fTOFExpTimeV0MPtPosP_Mu->Fill(track->Pt(),P_expTdiffMu,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegP_Mu->Fill(track->Pt(),P_expTdiffMu,V0mpc);

	//Pion
	if(track->Charge()>0) fTOFExpTimeV0MPtPosP_Pi->Fill(track->Pt(),(P_expTdiffPi),V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegP_Pi->Fill(track->Pt(),(P_expTdiffPi),V0mpc);
	
	//Kaon
	if(track->Charge()>0) fTOFExpTimeV0MPtPosP_K->Fill(track->Pt(),P_expTdiffK,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegP_K->Fill(track->Pt(),P_expTdiffK,V0mpc);

	//Proton
	if(track->Charge()>0) fTOFExpTimeV0MPtPosP_P->Fill(track->Pt(),P_expTdiffP,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegP_P->Fill(track->Pt(),P_expTdiffP,V0mpc);

	//Deuteron
	if(track->Charge()>0) fTOFExpTimeV0MPtPosP_D->Fill(track->Pt(),P_expTdiffD,V0mpc);
	if(track->Charge()<0) fTOFExpTimeV0MPtNegP_D->Fill(track->Pt(),P_expTdiffD,V0mpc);

	Double_t mismtch_P=AliTOFPIDResponse::GetMismatchRandomValue(track->Eta());

	//TOF Mismatch Time
	if(track->Charge()>0) fTOFMismatchTimeV0MPtPosP->Fill(track->Pt(),(mismtch_P-fTOFExpTimeP),V0mpc);
	if(track->Charge()<0) fTOFMismatchTimeV0MPtNegP->Fill(track->Pt(),(mismtch_P-fTOFExpTimeP),V0mpc);

}//rapidity Proton
}//eta cut
}//TOFpid	
}//track cut with dcaxy


        if(fesdTrackCuts_no_dca->AcceptTrack(track)){
	if(TMath::Abs(track->Eta())<0.8){
	//dca xy information

	Double_t Y_Pi=Rapidity(track,AliPID::ParticleMass(AliPID::kPion));
	Double_t Y_K=Rapidity(track,AliPID::ParticleMass(AliPID::kKaon));
	Double_t Y_P=Rapidity(track,AliPID::ParticleMass(AliPID::kProton));


	Double_t fTOFSigmaPi = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        Double_t fTOFSigmaK = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
        Double_t fTOFSigmaP = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

        //TPC nsigma
        Double_t fTPCSigmaPi = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        Double_t fTPCSigmaK = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        Double_t fTPCSigmaP = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);


	if(TMath::Abs(Y_Pi)<0.5){

	//for DCA xy TPCTOF Combined is used
	Double_t fTPCTOFSigmaPi = TMath::Sqrt(fTPCSigmaPi*fTPCSigmaPi + fTOFSigmaPi*fTOFSigmaPi);
	if(TMath::Abs(fTPCTOFSigmaPi)<2.){
	if(track->Charge()>0) fTOFDCAxyV0MPtPosPi->Fill(track->Pt(),V0mpc,dcaxy[0]);
	if(track->Charge()<0) fTOFDCAxyV0MPtNegPi->Fill(track->Pt(),V0mpc,dcaxy[0]);
}
}

	if(TMath::Abs(Y_K)<0.5){

	//for DCA xy TPCTOF Combined is used
	Double_t fTPCTOFSigmaK = TMath::Sqrt(fTPCSigmaK*fTPCSigmaK + fTOFSigmaK*fTOFSigmaK);
	if(TMath::Abs(fTPCTOFSigmaK)<2.){
	if(track->Charge()>0) fTOFDCAxyV0MPtPosK->Fill(track->Pt(),V0mpc,dcaxy[0]);
	if(track->Charge()<0) fTOFDCAxyV0MPtNegK->Fill(track->Pt(),V0mpc,dcaxy[0]);
}
}
		
	if(TMath::Abs(Y_P)<0.5){

	//for DCA xy TPCTOF Combined is used
        Double_t fTPCTOFSigmaP = TMath::Sqrt(fTPCSigmaP*fTPCSigmaP + fTOFSigmaP*fTOFSigmaP);
        if(TMath::Abs(fTPCTOFSigmaP)<2.){
        if(track->Charge()>0) fTOFDCAxyV0MPtPosP->Fill(track->Pt(),V0mpc,dcaxy[0]);
        if(track->Charge()<0) fTOFDCAxyV0MPtNegP->Fill(track->Pt(),V0mpc,dcaxy[0]);
}
}

}//eta cut
}// no dca xy cut for contamination estimation
}//track loop end


        fEventCounter->Fill(8.5);
PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskTOFppSpectra::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
	fEventCounter = dynamic_cast<TH1F*> (fOutputList->At(0));
	fEventPS = dynamic_cast<TH1F*> (fOutputList->At(1));
	fEventVtx = dynamic_cast<TH1F*> (fOutputList->At(2));
	fEventVtx10 = dynamic_cast<TH1F*> (fOutputList->At(3));
	fZVertex = dynamic_cast<TH1F*> (fOutputList->At(4));

	fTOFTimeV0MPtPosPi= dynamic_cast<TH3F*> (fOutputList->At(5));
	fTOFTimeV0MPtNegPi= dynamic_cast<TH3F*> (fOutputList->At(6));

	fTOFExpTimeV0MPtPosPi_El= dynamic_cast<TH3F*> (fOutputList->At(7));
	fTOFExpTimeV0MPtNegPi_El= dynamic_cast<TH3F*> (fOutputList->At(8));

	fTOFExpTimeV0MPtPosPi_Mu= dynamic_cast<TH3F*> (fOutputList->At(9));
	fTOFExpTimeV0MPtNegPi_Mu= dynamic_cast<TH3F*> (fOutputList->At(10));

	fTOFExpTimeV0MPtPosPi_Pi= dynamic_cast<TH3F*> (fOutputList->At(11));
	fTOFExpTimeV0MPtNegPi_Pi= dynamic_cast<TH3F*> (fOutputList->At(12));

	fTOFExpTimeV0MPtPosPi_K= dynamic_cast<TH3F*> (fOutputList->At(13));
	fTOFExpTimeV0MPtNegPi_K= dynamic_cast<TH3F*> (fOutputList->At(14));

	fTOFExpTimeV0MPtPosPi_P= dynamic_cast<TH3F*> (fOutputList->At(15));
	fTOFExpTimeV0MPtNegPi_P= dynamic_cast<TH3F*> (fOutputList->At(16));

	fTOFExpTimeV0MPtPosPi_D= dynamic_cast<TH3F*> (fOutputList->At(17));
	fTOFExpTimeV0MPtNegPi_D= dynamic_cast<TH3F*> (fOutputList->At(18));

	fTOFMismatchTimeV0MPtPosPi= dynamic_cast<TH3F*> (fOutputList->At(19));
	fTOFMismatchTimeV0MPtNegPi= dynamic_cast<TH3F*> (fOutputList->At(20));


	fTOFDCAxyV0MPtPosPi= dynamic_cast<TH3F*> (fOutputList->At(21));
	fTOFDCAxyV0MPtNegPi= dynamic_cast<TH3F*> (fOutputList->At(22));

	fEventV0M= dynamic_cast<TH2F*> (fOutputList->At(23));

	fT0Resolution= dynamic_cast<TH1F*> (fOutputList->At(24));
	fTimeOfFlightRes= dynamic_cast<TH1F*> (fOutputList->At(25));
	fTimeOfFlightTOFRes= dynamic_cast<TH1F*> (fOutputList->At(26));
	fTimeOfFlightGoodRes= dynamic_cast<TH1F*> (fOutputList->At(27));

	fBetaP= dynamic_cast<TH2F*> (fOutputList->At(28));
	fBetaPNoMismatch= dynamic_cast<TH2F*> (fOutputList->At(29));
	fBetaPNoMismatchEtaCut= dynamic_cast<TH2F*> (fOutputList->At(30));
	fBetaPt= dynamic_cast<TH2F*> (fOutputList->At(31));
	fBetaPtNoMismatch= dynamic_cast<TH2F*> (fOutputList->At(32));
	fBetaPtNoMismatchEtaCut= dynamic_cast<TH2F*> (fOutputList->At(33));

	fTPCdEdxP= dynamic_cast<TH2F*> (fOutputList->At(34));
	fTPCdEdxPt= dynamic_cast<TH2F*> (fOutputList->At(35));

	fTPCTOFnSigmaPi= dynamic_cast<TH2F*> (fOutputList->At(36));
	fTOFChannelVsTime= dynamic_cast<TH2F*> (fOutputList->At(37));

	fGausTime= dynamic_cast<TH1F*> (fOutputList->At(38));
	fTOFGausTime= dynamic_cast<TH1F*> (fOutputList->At(39));


	fTOFTimeV0MPtPosK= dynamic_cast<TH3F*> (fOutputList->At(40));
	fTOFTimeV0MPtNegK= dynamic_cast<TH3F*> (fOutputList->At(41));


	fTOFExpTimeV0MPtPosK_El= dynamic_cast<TH3F*> (fOutputList->At(42));
	fTOFExpTimeV0MPtNegK_El= dynamic_cast<TH3F*> (fOutputList->At(43));

	fTOFExpTimeV0MPtPosK_Mu= dynamic_cast<TH3F*> (fOutputList->At(44));
	fTOFExpTimeV0MPtNegK_Mu= dynamic_cast<TH3F*> (fOutputList->At(45));

	fTOFExpTimeV0MPtPosK_Pi= dynamic_cast<TH3F*> (fOutputList->At(46));
	fTOFExpTimeV0MPtNegK_Pi= dynamic_cast<TH3F*> (fOutputList->At(47));

	fTOFExpTimeV0MPtPosK_K= dynamic_cast<TH3F*> (fOutputList->At(48));
	fTOFExpTimeV0MPtNegK_K= dynamic_cast<TH3F*> (fOutputList->At(49));

	fTOFExpTimeV0MPtPosK_P= dynamic_cast<TH3F*> (fOutputList->At(50));
	fTOFExpTimeV0MPtNegK_P= dynamic_cast<TH3F*> (fOutputList->At(51));

	fTOFExpTimeV0MPtPosK_D= dynamic_cast<TH3F*> (fOutputList->At(52));
	fTOFExpTimeV0MPtNegK_D= dynamic_cast<TH3F*> (fOutputList->At(53));

	fTOFMismatchTimeV0MPtPosK= dynamic_cast<TH3F*> (fOutputList->At(54));
	fTOFMismatchTimeV0MPtNegK= dynamic_cast<TH3F*> (fOutputList->At(55));

	fTOFDCAxyV0MPtPosK= dynamic_cast<TH3F*> (fOutputList->At(56));
	fTOFDCAxyV0MPtNegK= dynamic_cast<TH3F*> (fOutputList->At(57));

	fTPCTOFnSigmaK= dynamic_cast<TH2F*> (fOutputList->At(58));
	

	fTOFTimeV0MPtPosP= dynamic_cast<TH3F*> (fOutputList->At(59));
	fTOFTimeV0MPtNegP= dynamic_cast<TH3F*> (fOutputList->At(60));

	fTOFExpTimeV0MPtPosP_El= dynamic_cast<TH3F*> (fOutputList->At(61));
	fTOFExpTimeV0MPtNegP_El= dynamic_cast<TH3F*> (fOutputList->At(62));

	fTOFExpTimeV0MPtPosP_Mu= dynamic_cast<TH3F*> (fOutputList->At(63));
	fTOFExpTimeV0MPtNegP_Mu= dynamic_cast<TH3F*> (fOutputList->At(64));

	fTOFExpTimeV0MPtPosP_Pi= dynamic_cast<TH3F*> (fOutputList->At(65));
	fTOFExpTimeV0MPtNegP_Pi= dynamic_cast<TH3F*> (fOutputList->At(66));

	fTOFExpTimeV0MPtPosP_K= dynamic_cast<TH3F*> (fOutputList->At(67));
	fTOFExpTimeV0MPtNegP_K= dynamic_cast<TH3F*> (fOutputList->At(68));

	fTOFExpTimeV0MPtPosP_P= dynamic_cast<TH3F*> (fOutputList->At(69));
	fTOFExpTimeV0MPtNegP_P= dynamic_cast<TH3F*> (fOutputList->At(70));

	fTOFExpTimeV0MPtPosP_D= dynamic_cast<TH3F*> (fOutputList->At(71));
	fTOFExpTimeV0MPtNegP_D= dynamic_cast<TH3F*> (fOutputList->At(72));

	fTOFMismatchTimeV0MPtPosP= dynamic_cast<TH3F*> (fOutputList->At(73));
	fTOFMismatchTimeV0MPtNegP= dynamic_cast<TH3F*> (fOutputList->At(74));

	fTOFDCAxyV0MPtPosP= dynamic_cast<TH3F*> (fOutputList->At(75));
	fTOFDCAxyV0MPtNegP= dynamic_cast<TH3F*> (fOutputList->At(76));

	fTPCTOFnSigmaP= dynamic_cast<TH2F*> (fOutputList->At(77));
	fEventV0MPS= dynamic_cast<TH2F*> (fOutputList->At(78));
        fEventV0MVtx= dynamic_cast<TH2F*> (fOutputList->At(79));
        fV0MPC= dynamic_cast<TH1F*> (fOutputList->At(80));
	ftail= dynamic_cast<TH1F*> (fOutputList->At(81));
	fV0MPC_vertexcut= dynamic_cast<TH1F*> (fOutputList->At(82));
	
	fPtTPC_AllP= dynamic_cast<TH1F*> (fOutputList->At(83));
	fPtTPC_AllN= dynamic_cast<TH1F*> (fOutputList->At(84));
	fPtTOF_AllP= dynamic_cast<TH1F*> (fOutputList->At(85));
	fPtTOF_AllN= dynamic_cast<TH1F*> (fOutputList->At(86));

	ftail_Random= dynamic_cast<TH1F*> (fOutputList->At(87));
	fTPC_CR= dynamic_cast<TH1F*> (fOutputList->At(88));
        fChi2TPCcluster= dynamic_cast<TH1F*> (fOutputList->At(89));
        fDCAZ= dynamic_cast<TH1F*> (fOutputList->At(90));
        fDCAxy= dynamic_cast<TH1F*> (fOutputList->At(91));

/*
TFile *f=new TFile("result/12dec/ALICE_final_TOF_output_run1.root","recreate");//runlist1
f->cd();

fTPC_CR->Write();
fChi2TPCcluster->Write();
fDCAZ->Write();
fDCAxy->Write();

ftail->Write();
ftail_Random->Write();
fEventCounter->Write();
fEventPS->Write();
fEventVtx->Write();
fV0MPC->Write();
fV0MPC_vertexcut->Write();
fEventVtx10->Write();
fZVertex->Write();
fT0Resolution->Write();
fTimeOfFlightRes->Write();
fTimeOfFlightTOFRes->Write();
fTimeOfFlightGoodRes->Write();
fTPCdEdxP->Write();
fTPCdEdxPt->Write();
fBetaP->Write();
fBetaPNoMismatch->Write();
fBetaPNoMismatchEtaCut->Write();
fBetaPt->Write();
fBetaPtNoMismatch->Write();
fBetaPtNoMismatchEtaCut->Write();
fTPCTOFnSigmaPi->Write();
fTOFChannelVsTime->Write();
fGausTime->Write();
fTOFGausTime->Write();
fGausTime_K->Write();
fTOFGausTime_K->Write();
fGausTime_P->Write();
fTOFGausTime_P->Write();

fTOFTimeV0MPtPi->Write();
fTOFTimeV0MPtPosPi->Write();
fTOFTimeV0MPtNegPi->Write();
fTOFSigmaV0MPtPi->Write();
fTOFSigmaV0MPtPosPi->Write();
fTOFSigmaV0MPtNegPi->Write();

fTOFNoMismatchTimeV0MPtPi->Write();
fTOFNoMismatchTimeV0MPtPosPi->Write();
fTOFNoMismatchTimeV0MPtNegPi->Write();
fTOFNoMismatchSigmaV0MPtPi->Write();
fTOFNoMismatchSigmaV0MPtPosPi->Write();
fTOFNoMismatchSigmaV0MPtNegPi->Write();

fTOFResolutionV0MPtPi->Write();
fTOFResolutionV0MPtPosPi->Write();
fTOFResolutionV0MPtNegPi->Write();

fTOFExpTimeV0MPtPi_El->Write();
fTOFExpTimeV0MPtPosPi_El->Write();
fTOFExpTimeV0MPtNegPi_El->Write();
fTOFExpSigmaV0MPtPi_El->Write();
fTOFExpSigmaV0MPtPosPi_El->Write();
fTOFExpSigmaV0MPtNegPi_El->Write();

fTOFExpTimeV0MPtPi_Mu->Write();
fTOFExpTimeV0MPtPosPi_Mu->Write();
fTOFExpTimeV0MPtNegPi_Mu->Write();
fTOFExpSigmaV0MPtPi_Mu->Write();
fTOFExpSigmaV0MPtPosPi_Mu->Write();
fTOFExpSigmaV0MPtNegPi_Mu->Write();

fTOFExpTimeV0MPtPi_Pi->Write();
fTOFExpTimeV0MPtPosPi_Pi->Write();
fTOFExpTimeV0MPtNegPi_Pi->Write();
fTOFExpSigmaV0MPtPi_Pi->Write();
fTOFExpSigmaV0MPtPosPi_Pi->Write();
fTOFExpSigmaV0MPtNegPi_Pi->Write();

fTOFExpTimeV0MPtPi_K->Write();
fTOFExpTimeV0MPtPosPi_K->Write();
fTOFExpTimeV0MPtNegPi_K->Write();
fTOFExpSigmaV0MPtPi_K->Write();
fTOFExpSigmaV0MPtPosPi_K->Write();
fTOFExpSigmaV0MPtNegPi_K->Write();

fTOFExpTimeV0MPtPi_P->Write();
fTOFExpTimeV0MPtPosPi_P->Write();
fTOFExpTimeV0MPtNegPi_P->Write();
fTOFExpSigmaV0MPtPi_P->Write();
fTOFExpSigmaV0MPtPosPi_P->Write();
fTOFExpSigmaV0MPtNegPi_P->Write();

fTOFExpTimeV0MPtPi_D->Write();
fTOFExpTimeV0MPtPosPi_D->Write();
fTOFExpTimeV0MPtNegPi_D->Write();
fTOFExpSigmaV0MPtPi_D->Write();
fTOFExpSigmaV0MPtPosPi_D->Write();
fTOFExpSigmaV0MPtNegPi_D->Write();

fTOFMismatchTimeV0MPtPi->Write();
fTOFMismatchTimeV0MPtPosPi->Write();
fTOFMismatchTimeV0MPtNegPi->Write();

fTOFMismatchSigmaV0MPtPi->Write();
fTOFMismatchSigmaV0MPtPosPi->Write();
fTOFMismatchSigmaV0MPtNegPi->Write();

fTOFDCAxyV0MPtPi->Write();
fTOFDCAxyV0MPtPosPi->Write();
fTOFDCAxyV0MPtNegPi->Write();
fEventV0MPS->Write();
fEventV0MVtx->Write();
fEventV0M->Write();




fTPCTOFnSigmaK->Write();

fTOFTimeV0MPtK->Write();
fTOFTimeV0MPtPosK->Write();
fTOFTimeV0MPtNegK->Write();

fTOFSigmaV0MPtK->Write();
fTOFSigmaV0MPtPosK->Write();
fTOFSigmaV0MPtNegK->Write();

fTOFNoMismatchTimeV0MPtK->Write();
fTOFNoMismatchTimeV0MPtPosK->Write();
fTOFNoMismatchTimeV0MPtNegK->Write();
fTOFNoMismatchSigmaV0MPtK->Write();
fTOFNoMismatchSigmaV0MPtPosK->Write();
fTOFNoMismatchSigmaV0MPtNegK->Write();


fTOFResolutionV0MPtK->Write();
fTOFResolutionV0MPtPosK->Write();
fTOFResolutionV0MPtNegK->Write();

fTOFExpTimeV0MPtK_El->Write();
fTOFExpTimeV0MPtPosK_El->Write();
fTOFExpTimeV0MPtNegK_El->Write();
fTOFExpSigmaV0MPtK_El->Write();
fTOFExpSigmaV0MPtPosK_El->Write();
fTOFExpSigmaV0MPtNegK_El->Write();

fTOFExpTimeV0MPtK_Mu->Write();
fTOFExpTimeV0MPtPosK_Mu->Write();
fTOFExpTimeV0MPtNegK_Mu->Write();
fTOFExpSigmaV0MPtK_Mu->Write();
fTOFExpSigmaV0MPtPosK_Mu->Write();
fTOFExpSigmaV0MPtNegK_Mu->Write();

fTOFExpTimeV0MPtK_Pi->Write();
fTOFExpTimeV0MPtPosK_Pi->Write();
fTOFExpTimeV0MPtNegK_Pi->Write();
fTOFExpSigmaV0MPtK_Pi->Write();
fTOFExpSigmaV0MPtPosK_Pi->Write();
fTOFExpSigmaV0MPtNegK_Pi->Write();

fTOFExpTimeV0MPtK_K->Write();
fTOFExpTimeV0MPtPosK_K->Write();
fTOFExpTimeV0MPtNegK_K->Write();
fTOFExpSigmaV0MPtK_K->Write();
fTOFExpSigmaV0MPtPosK_K->Write();
fTOFExpSigmaV0MPtNegK_K->Write();

fTOFExpTimeV0MPtK_P->Write();
fTOFExpTimeV0MPtPosK_P->Write();
fTOFExpTimeV0MPtNegK_P->Write();
fTOFExpSigmaV0MPtK_P->Write();
fTOFExpSigmaV0MPtPosK_P->Write();
fTOFExpSigmaV0MPtNegK_P->Write();

fTOFExpTimeV0MPtK_D->Write();
fTOFExpTimeV0MPtPosK_D->Write();
fTOFExpTimeV0MPtNegK_D->Write();
fTOFExpSigmaV0MPtK_D->Write();
fTOFExpSigmaV0MPtPosK_D->Write();
fTOFExpSigmaV0MPtNegK_D->Write();

fTOFMismatchTimeV0MPtK->Write();
fTOFMismatchTimeV0MPtPosK->Write();
fTOFMismatchTimeV0MPtNegK->Write();

fTOFMismatchSigmaV0MPtK->Write();
fTOFMismatchSigmaV0MPtPosK->Write();
fTOFMismatchSigmaV0MPtNegK->Write();

fTOFDCAxyV0MPtK->Write();
fTOFDCAxyV0MPtPosK->Write();
fTOFDCAxyV0MPtNegK->Write();


fTPCTOFnSigmaP->Write();

fTOFTimeV0MPtP->Write();
fTOFTimeV0MPtPosP->Write();
fTOFTimeV0MPtNegP->Write();

fTOFSigmaV0MPtP->Write();
fTOFSigmaV0MPtPosP->Write();
fTOFSigmaV0MPtNegP->Write();


fTOFNoMismatchTimeV0MPtP->Write();
fTOFNoMismatchTimeV0MPtPosP->Write();
fTOFNoMismatchTimeV0MPtNegP->Write();
fTOFNoMismatchSigmaV0MPtP->Write();
fTOFNoMismatchSigmaV0MPtPosP->Write();
fTOFNoMismatchSigmaV0MPtNegP->Write();


fTOFResolutionV0MPtP->Write();
fTOFResolutionV0MPtPosP->Write();
fTOFResolutionV0MPtNegP->Write();

fTOFExpTimeV0MPtP_El->Write();
fTOFExpTimeV0MPtPosP_El->Write();
fTOFExpTimeV0MPtNegP_El->Write();
fTOFExpSigmaV0MPtP_El->Write();
fTOFExpSigmaV0MPtPosP_El->Write();
fTOFExpSigmaV0MPtNegP_El->Write();

fTOFExpTimeV0MPtP_Mu->Write();
fTOFExpTimeV0MPtPosP_Mu->Write();
fTOFExpTimeV0MPtNegP_Mu->Write();
fTOFExpSigmaV0MPtP_Mu->Write();
fTOFExpSigmaV0MPtPosP_Mu->Write();
fTOFExpSigmaV0MPtNegP_Mu->Write();

fTOFExpTimeV0MPtP_Pi->Write();
fTOFExpTimeV0MPtPosP_Pi->Write();
fTOFExpTimeV0MPtNegP_Pi->Write();
fTOFExpSigmaV0MPtP_Pi->Write();
fTOFExpSigmaV0MPtPosP_Pi->Write();
fTOFExpSigmaV0MPtNegP_Pi->Write();

fTOFExpTimeV0MPtP_K->Write();
fTOFExpTimeV0MPtPosP_K->Write();
fTOFExpTimeV0MPtNegP_K->Write();
fTOFExpSigmaV0MPtP_K->Write();
fTOFExpSigmaV0MPtPosP_K->Write();
fTOFExpSigmaV0MPtNegP_K->Write();

fTOFExpTimeV0MPtP_P->Write();
fTOFExpTimeV0MPtPosP_P->Write();
fTOFExpTimeV0MPtNegP_P->Write();
fTOFExpSigmaV0MPtP_P->Write();
fTOFExpSigmaV0MPtPosP_P->Write();
fTOFExpSigmaV0MPtNegP_P->Write();

fTOFExpTimeV0MPtP_D->Write();
fTOFExpTimeV0MPtPosP_D->Write();
fTOFExpTimeV0MPtNegP_D->Write();
fTOFExpSigmaV0MPtP_D->Write();
fTOFExpSigmaV0MPtPosP_D->Write();
fTOFExpSigmaV0MPtNegP_D->Write();

fTOFMismatchTimeV0MPtP->Write();
fTOFMismatchTimeV0MPtPosP->Write();
fTOFMismatchTimeV0MPtNegP->Write();

fTOFMismatchSigmaV0MPtP->Write();
fTOFMismatchSigmaV0MPtPosP->Write();
fTOFMismatchSigmaV0MPtNegP->Write();

fTOFDCAxyV0MPtP->Write();
fTOFDCAxyV0MPtPosP->Write();
fTOFDCAxyV0MPtNegP->Write();

fPtTPC_AllP->Write();
fPtTPC_AllN->Write();
fPtTOF_AllP->Write();
fPtTOF_AllN->Write();
*/
}
/*
//_________________________________________________________________________
Float_t AliAnalysisTaskTOFppSpectra::GetVertex(AliESDEvent* esd) const
{
  Float_t zvtx = -999;
  //const AliAODVertex* vtxAOD = aod->GetPrimaryVertex();
 const AliESDVertex* vtxESD = esd->GetPrimaryVertex();
  if (!vtxESD)
    return zvtx;
  if(vtxESD->GetNContributors()>0)
    zvtx = vtxESD->GetZ();
  return zvtx;
}
*/
//-----------------------------------------------------------------
Bool_t AliAnalysisTaskTOFppSpectra::selectVertex2015pp(AliESDEvent *esd,
			  Bool_t checkSPDres, //enable check on vtx resolution
			  Bool_t requireSPDandTrk, //ask for both trk and SPD vertex 
			  Bool_t checkProximity) //apply cut on relative position of spd and trk verteces 
{

  if (!esd) return kFALSE;
  
  const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
  const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
  Bool_t hasSPD = spdVertex->GetStatus();
  Bool_t hasTrk = trkVertex->GetStatus();
 
  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;
  
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
      if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE; 
    }
  }

/*
  //Cut on the vertex z position
  const AliESDVertex *vertex = esd->GetPrimaryVertex();
  if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
  return kTRUE;
*/
}
//_________________________________________________________________________________________________
Bool_t AliAnalysisTaskTOFppSpectra::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}
//_______________________________________________________________________________________
Bool_t AliAnalysisTaskTOFppSpectra::TPCPID(AliESDtrack *track)
{
    if ((track->GetStatus() & AliESDtrack::kTPCin   ) == 0) return kFALSE;
    if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) return kFALSE;
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) return kFALSE;
    return kTRUE;
}

// =============== For TOF PID CHECK =================

Bool_t AliAnalysisTaskTOFppSpectra:: TOFPID(AliESDtrack * track)
{
    Double_t TOFPtRange =0.6;         // Controling Agent 1

// *************** Check if the particle has TOF Matching ***************
    UInt_t status;
    status=track->GetStatus();

    if ((status & AliESDtrack::kITSrefit) == 0) return kFALSE;
    if ((status & AliESDtrack::kTPCrefit) == 0) return kFALSE;

//    if((status&AliESDtrack::kTOFout)==0) return kFALSE;
//    if((status&AliESDtrack::kTIME)==0) return kFALSE; 
//    if((status&AliESDtrack::kTRDout)==0) return kFALSE;

// same condition as above   
 if((status&AliESDtrack::kTOFout)==0 && (status&AliESDtrack::kTIME)==0)// || (status&AliESDtrack::kTRDout)==0)
	     return kFALSE;

 // TPC TOF mismatch is be implemented
    Float_t length = track->GetIntegratedLength();
    if (length < 350.)      
        return kFALSE;

// -------  in addition to KTOFout and kTIME we look at the pt  ------
//    if(track->Pt()<TOFPtRange) return kFALSE;
    return kTRUE;
}
//================================================================
Double_t AliAnalysisTaskTOFppSpectra::Rapidity(AliESDtrack *track , Double_t mass)
{
    Double_t E,rap,pz,pt;
    pt=track->Pt();
    pz = track->Pz();
    E = TMath::Sqrt(pt*pt+pz*pz+mass*mass);
    rap = 0.5 * TMath::Log ((E+pz)/(E-pz));
    return rap;
}
