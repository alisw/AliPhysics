/*   This macro produces: Delta phi correclation, leading particle spectra in different multiplicity classes and RT classess. in addition pT spectra for near side, away side and transverse side in diffreent Multiplicity classes.
     Aditya Nath Mishra ICN-UNAM
     Please report bugs to: amishra@cern.ch / aditya.mishra@correo.nucleares.unam.mx 
     First version: 25/07/2019
 */

#include "AliAnalysisTaskUeMultDep.h"

// ROOT includes
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2D.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TFile.h>


// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliEventCuts.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include <TTreeStream.h>
#include <AliHeader.h>
#include <AliAnalysisUtils.h>
#include <AliAODInputHandler.h> 
#include <AliAODHandler.h> 
#include <AliAODVertex.h>
#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliDataFile.h>

#include "AliPPVsMultUtils.h"
#include "AliAnalysisUtils.h"

#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliMultSelectionTask.h"

#include <iostream>
using namespace std;

const Char_t * legTrackRT[15] = {"#it{R}_{T}>0", "0<#it{R}_{T}<1", "1<#it{R}_{T}<2", "2<#it{R}_{T}<3", "3<#it{R}_{T}<4", "4<#it{R}_{T}<5", "5<#it{R}_{T}<6", "6<#it{R}_{T}<7", "7<#it{R}_{T}<8", "8<#it{R}_{T}<9", "9 < #it{R}_{T} < 10", "10 < #it{R}_{T} < 11", "11 < #it{R}_{T} < 12", "12 < #it{R}_{T} < 13", "13 < #it{R}_{T} < 14"};

const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;

ClassImp(AliAnalysisTaskUeMultDep)

//_____________________________________________________________________________
AliAnalysisTaskUeMultDep::AliAnalysisTaskUeMultDep():
AliAnalysisTaskSE(),
  fESD(0x0),
  fEventCuts(0x0),
  fTrackFilter(0x0),
  fAnalysisType("ESD"),
  ftrigBit(0x0),
  fNcl(70),
  fListOfObjects(0x0),
  fHistEventCounter(0x0),
  fTriggeredEventMB(-999),
  fVtxCut(10.0),
  fEtaCut(0.8),
  fisPS(kFALSE),
  fisTracklet(kFALSE),
  isINEL0Rec(kFALSE),
  fVtxBeforeCuts(0x0), 
  fVtxAfterCuts(0x0),
  fPileUpRej(kFALSE),
  fUeRTUtils(0x0),
  fMultSelection(0x0),
  hPhi(0x0),
  hpT(0x0),
  hPtL(0x0),
  hPhiL(0x0),
  hEtaL(0x0),
  hDphi(0x0),
  hDphiNS(0x0),
  hDphiAS(0x0),
  hDphiTS(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  ProfpTLvsNch(0x0),
  ProfpTLvsNchNS(0x0),
  ProfpTLvsNchAS(0x0),
  ProfpTLvsNchTS(0x0),
  hpTvsRefMult08(0x0),
  hpTLvsRefMult08(0x0),
  hpTvsV0Mmult(0x0),
  hpTLvsV0Mmult(0x0),
  hRefMultvsV0Mmult(0x0),
  hpTLvsRefMult08vsDphi(0x0),
  hpTLvsV0MmultvsDphi(0x0),
  hpTvspTLvsRefMult08(0x0),
  hpTvspTLvsRefMult08NS(0x0),
  hpTvspTLvsRefMult08AS(0x0),
  hpTvspTLvsRefMult08TS(0x0),
  hpTvspTLvsV0Mmult(0x0),
  hpTvspTLvsV0MmultNS(0x0),
  hpTvspTLvsV0MmultAS(0x0),
  hpTvspTLvsV0MmultTS(0x0),
  ftrackmult08(-999),   
  fv0mpercentile(-999),
  hSumptVsRefMult08(0x0),
  hINEL0(0x0),
  hTrackRT(0x0),
  hTrackRTvsRefMult08(0x0),
  hMultTSvsRefMult08(0x0),
  hSumptVsTrackRT(0x0),
  hpTLvsRefMult08vsRT(0x0),
  hDeltaphiVspTLVsRT(0x0),
  hDeltaphiVsMultVsRT(0x0),
  hINEL0pTL5(0x0),
  hTrackRTpTL5(0x0),
  hTrackRTvsRefMult08pTL5(0x0),
  hMultTSvsRefMult08pTL5(0x0),
  hSumptVsTrackRTpTL5(0x0),
  hpTLvsRefMult08vsRTpTL5(0x0),
  hDeltaphiVspTLVsRTpTL5(0x0),
  hDeltaphiVsMultVsRTpTL5(0x0)
{
  for ( int i = 0 ; i < 15 ; i++ )
    {
     hpTLRT[i] =  0;
     hpTRT[i] =  0;
     hpTNSRT[i] = 0;
     hpTASRT[i] = 0;
     hpTTSRT[i] = 0;
     hDEtaDPhiRT[i] = 0;
     hDeltaphiVspTLRT[i] = 0;

     hpTLRTpTL5[i] =  0;
     hpTRTpTL5[i] =  0;
     hpTNSRTpTL5[i] = 0;
     hpTASRTpTL5[i] = 0;
     hpTTSRTpTL5[i] = 0;
     hDEtaDPhiRTpTL5[i] = 0;
     hDeltaphiVspTLRTpTL5[i] = 0;
    }
  
  // Default constructor (should not be used)
}

//______________________________________________________________________________
AliAnalysisTaskUeMultDep::AliAnalysisTaskUeMultDep(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fEventCuts(0x0),
  fTrackFilter(0x0),
  fAnalysisType("ESD"),
  ftrigBit(0x0),
  fNcl(70),
  fListOfObjects(0x0),
  fHistEventCounter(0x0),
  fTriggeredEventMB(-999),
  fisPS(kFALSE),
  fisTracklet(kFALSE),
  isINEL0Rec(kFALSE),
  fVtxBeforeCuts(0x0), 
  fVtxAfterCuts(0x0),
  fPileUpRej(kFALSE),
  fUeRTUtils(0x0),
  fMultSelection(0x0),
  fVtxCut(10.0),
  fEtaCut(0.8),
  hPhi(0x0),
  hpT(0x0),
  hPtL(0x0),
  hPhiL(0x0),
  hEtaL(0x0),
  hDphi(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  ProfpTLvsNch(0x0),
  ProfpTLvsNchNS(0x0),
  ProfpTLvsNchAS(0x0),
  ProfpTLvsNchTS(0x0),
  hpTvsRefMult08(0x0),
  hpTLvsRefMult08(0x0),
  hpTvsV0Mmult(0x0),
  hpTLvsV0Mmult(0x0),
  hRefMultvsV0Mmult(0x0),
  hpTLvsRefMult08vsDphi(0x0),
  hpTLvsV0MmultvsDphi(0x0),
  hpTvspTLvsRefMult08(0x0),
  hpTvspTLvsRefMult08NS(0x0),
  hpTvspTLvsRefMult08AS(0x0),
  hpTvspTLvsRefMult08TS(0x0),
  hpTvspTLvsV0Mmult(0x0),
  hpTvspTLvsV0MmultNS(0x0),
  hpTvspTLvsV0MmultAS(0x0),
  hpTvspTLvsV0MmultTS(0x0),
  ftrackmult08(-999),   
  fv0mpercentile(-999),
  hSumptVsRefMult08(0x0),
  hINEL0(0x0),
  hTrackRT(0x0),
  hTrackRTvsRefMult08(0x0),
  hMultTSvsRefMult08(0x0),
  hSumptVsTrackRT(0x0),
  hpTLvsRefMult08vsRT(0x0),
  hDeltaphiVspTLVsRT(0x0),
  hDeltaphiVsMultVsRT(0x0),
  hINEL0pTL5(0x0),
  hTrackRTpTL5(0x0),
  hTrackRTvsRefMult08pTL5(0x0),
  hMultTSvsRefMult08pTL5(0x0),
  hSumptVsTrackRTpTL5(0x0),
  hpTLvsRefMult08vsRTpTL5(0x0),
  hDeltaphiVspTLVsRTpTL5(0x0),
  hDeltaphiVsMultVsRTpTL5(0x0)
{
  for ( int i = 0 ; i < 15 ; i++ )
    {
     hpTLRT[i] =  0;
     hpTRT[i] =  0;
     hpTNSRT[i] = 0;
     hpTASRT[i] = 0;
     hpTTSRT[i] = 0;
     hDEtaDPhiRT[i] = 0;
     hDeltaphiVspTLRT[i] = 0;

     hpTLRTpTL5[i] =  0;
     hpTRTpTL5[i] =  0;
     hpTNSRTpTL5[i] = 0;
     hpTASRTpTL5[i] = 0;
     hpTTSRTpTL5[i] = 0;
     hDEtaDPhiRTpTL5[i] = 0;
     hDeltaphiVspTLRTpTL5[i] = 0;
    }
  
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskUeMultDep::Exit(const char *msg) {

  Printf("%s", msg);
  return;
}


//_____________________________________________________________________________
AliAnalysisTaskUeMultDep::~AliAnalysisTaskUeMultDep()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fListOfObjects;
    fListOfObjects = 0x0;
  }

   if (fMultSelection)
  {
    delete fMultSelection;
    fMultSelection = 0x0;
  }

   if (fUeRTUtils)
  {
    delete fUeRTUtils;
    fUeRTUtils = 0x0;
  }

}

//______________________________________________________________________________
void AliAnalysisTaskUeMultDep::UserCreateOutputObjects()
{
  const Int_t nNchBins = 200;
  Double_t NchBins[nNchBins+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
				21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,
				39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,
				57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,
				75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,
				93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,
				109,110,111,112,113,114,115,116,117,118,119,120,121,122,
				123,124,125,126,127,128,129,130,131,132,133,134,135,136,
				137,138,139,140,141,142,143,144,145,146,147,148,149,150,
				151,152,153,154,155,156,157,158,159,160,161,162,163,164,
				165,166,167,168,169,170,171,172,173,174,175,176,177,178,
				179,180,181,182,183,184,185,186,187,188,189,190,191,192,
				193,194,195,196,197,198,199,200};

  const Int_t nPtBins      = 79; 
  Double_t PtBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
				0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
				1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
				2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
				4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
				11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
				26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0, 55.0,
				60.0, 65.0, 70.0, 80.0, 90.0, 100.0,120.0,140.0,170.0,200.0 };

 const Int_t nDeltabins = 64;
Double_t Deltabins[nDeltabins+1]={-1.0472, -0.957204, -0.867211, -0.777217, -0.687223, -0.59723, -0.507236, -0.417243, -0.327249, -0.237256, -0.147262, -0.0572686, 0.0327249, 0.122718, 0.212712, 0.302706, 0.392699, 0.482693, 0.572686, 0.66268, 0.752673, 0.842667, 0.93266, 1.02265, 1.11265, 1.20264, 1.29263, 1.38263, 1.47262, 1.56262, 1.65261, 1.7426, 1.8326, 1.92259, 2.01258, 2.10258, 2.19257, 2.28256, 2.37256, 2.46255, 2.55254, 2.64254, 2.73253, 2.82252, 2.91252, 3.00251, 3.09251, 3.1825, 3.27249, 3.36249, 3.45248, 3.54247, 3.63247, 3.72246, 3.81245, 3.90245, 3.99244, 4.08243, 4.17243, 4.26242, 4.35241, 4.44241, 4.5324, 4.6224, 4.71239};

 const Int_t nRTbins = 15;
 Double_t RTbins[nRTbins+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
 const Double_t RTmin[nRTbins+1]={0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,  10.0, 11.0, 12.0, 13.0};
 const Double_t RTmax[nRTbins+1]={100.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0};
 
  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested 

  // Definition of trackcuts
  if(!fTrackFilter){	
    fTrackFilter = new AliAnalysisFilter("trackFilter2015");
    SetTrackCuts(fTrackFilter);
  }

  OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  //
  // Histograms
  //  
  if(! fHistEventCounter ){
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",10,0,10);

    fHistEventCounter->GetXaxis()->SetBinLabel(1, "All events");// 
    fHistEventCounter->GetXaxis()->SetBinLabel(2, "MB Trigger "); //fTriggeredEventMB
    fHistEventCounter->GetXaxis()->SetBinLabel(3, "+ Incomplete DAQ"); //fTriggeredEventMB && !IncompleteDAQ 
    fHistEventCounter->GetXaxis()->SetBinLabel(4, "+ BG rejection"); //fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG
    fHistEventCounter->GetXaxis()->SetBinLabel(5, "+ Pile-Up rejection"); //fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp
    fHistEventCounter->GetXaxis()->SetBinLabel(6, "+ |#eta| < 1.0"); //fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && |#eta| < 1.0
    fHistEventCounter->GetXaxis()->SetBinLabel(7, "+ |v_{z}| < 10"); //fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && |#eta| < 1.0 && isVtxGood 
    fHistEventCounter->GetXaxis()->SetBinLabel(8, "Analyzed event"); //isINEL0Rec  = kTRUE;
    fListOfObjects->Add(fHistEventCounter);
  }

  fVtxBeforeCuts = new TH1D("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 12, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);
  
  fVtxAfterCuts = new TH1D("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 12, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);
  
  hPhi = new TH1D("hPhi", ";#phi (rad); count", 64,0,2.0*TMath::Pi());
  fListOfObjects->Add(hPhi);

  hpT = new TH1D("hpT",";#it{p}_{T} (GeV/#it{c});counts",nPtBins,PtBins);
  fListOfObjects->Add(hpT);

  hPtL = 0;
  hPtL = new TH1D("hPtL",";#it{p}_{T}^{leading} (GeV/#it{c});counts",nPtBins,PtBins);
  fListOfObjects->Add(hPtL);
	
  hEtaL = 0;
  hEtaL = new TH1D("hEtaL","; #eta^{leading};counts",20,-1,1);
  fListOfObjects->Add(hEtaL);
	
  hPhiL = 0;
  hPhiL = new TH1D("hPhiL","; #phi^{leading} (rad);counts",64,0,2.0*TMath::Pi());
  fListOfObjects->Add(hPhiL);

  hDphi = 0;
  hDphi = new TH1D("hDphi",";#Delta#phi (rad); count",64,-2*TMath::Pi(),2*TMath::Pi());
  fListOfObjects->Add(hDphi);

   hDphiNS = 0;
  hDphiNS = new TH1D("hDphiNS","Near Side",64,-2*TMath::Pi(),2*TMath::Pi());
  fListOfObjects->Add(hDphiNS);

  hDphiAS = 0;
  hDphiAS = new TH1D("hDphiAS","Away Side",64,-2*TMath::Pi(),2*TMath::Pi());
  fListOfObjects->Add(hDphiAS);

  hDphiTS = 0;
  hDphiTS = new TH1D("hDphiTS","Transverse Side",64,-2*TMath::Pi(),2*TMath::Pi());
  fListOfObjects->Add(hDphiTS);

  hRefMult08 = 0;
  hRefMult08 = new TH1D("hRefMult08","Multiplicity (-0.8 < #eta < 0.8);N_{ch};count",nNchBins,NchBins);   
  fListOfObjects->Add(hRefMult08);
	
  hV0Mmult = 0;
  hV0Mmult = new TH1D("hV0Mmult","V0M ;V0M percentile;count",nNchBins,NchBins);
  fListOfObjects->Add(hV0Mmult);
		
  hpTvsRefMult08 = 0;
  hpTvsRefMult08 = new TH2D("hpTvsRefMult08","p_{T} vs RefMult08; #it{p}_{T} (GeV/#it{c}); RefMult08 (|#eta| < 0.8 & p_{T} > 0.15 GeV/c)",nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvsRefMult08);
	
  hpTLvsRefMult08 = 0;
  hpTLvsRefMult08 = new TH2D("hpTLvsRefMult08","p_{T}^{leading} vs RefMult08; #it{p}_{T}^{leading} (GeV/#it{c}); RefMult08 (|#eta| < 0.8 & p_{T} > 0.15 GeV/c)",nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTLvsRefMult08);
	
  hpTvsV0Mmult = 0;
  hpTvsV0Mmult = new TH2D("hpTvsV0Mmult","p_{T} vs V0Mmult; #it{p}_{T} (GeV/#it{c}); V0Mmult (|#eta| < 0.8 & p_{T} > 0.15 GeV/c)",nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvsV0Mmult);

   ProfpTLvsNch = 0;
  ProfpTLvsNch = new TProfile("ProfpTLvsNch","Full #phi; #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.15 GeV/#it{c})",nPtBins,PtBins,0,20);
  fListOfObjects->Add(ProfpTLvsNch);

   ProfpTLvsNchNS = 0;
  ProfpTLvsNchNS = new TProfile("ProfpTLvsNchNS","Near side; #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.15 GeV/#it{c})",nPtBins,PtBins,0,20);
  fListOfObjects->Add(ProfpTLvsNchNS);

  ProfpTLvsNchAS = 0;
  ProfpTLvsNchAS = new TProfile("ProfpTLvsNchAS","Away side; #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.15 GeV/#it{c})",nPtBins,PtBins,0,20);
  fListOfObjects->Add(ProfpTLvsNchAS);

  ProfpTLvsNchTS = 0;
  ProfpTLvsNchTS = new TProfile("ProfpTLvsNchTS","Transverse side; #it{p}^{leading} (GeV/#it{c}); #it{N}_{ch} (#it{p}_{T}>0.15 GeV/#it{c})",nPtBins,PtBins,0,20);
  fListOfObjects->Add(ProfpTLvsNchTS);
	
  hpTLvsV0Mmult = 0;
  hpTLvsV0Mmult = new TH2D("hpTLvsV0Mmult","p_{T}^{leading} vs V0Mmult; #it{p}_{T}^{leading} (GeV/#it{c}); V0Mmult (|#eta| < 0.8 & p_{T} > 0.15 GeV/c)",nPtBins,PtBins,200,0,200);
  fListOfObjects->Add(hpTLvsV0Mmult);
	
  hRefMultvsV0Mmult = 0;
  hRefMultvsV0Mmult = new TH2D("hRefMultvsV0Mmult","N_{ch} vs V0M percentile;RefMult08 (|#eta| < 0.8 & p_{T} > 0.15 GeV/c); v0M percentile",nNchBins,NchBins,nNchBins,NchBins);
  fListOfObjects->Add(hRefMultvsV0Mmult);

  //for pTLeading vs RefMult vs Dphi...
  hpTLvsRefMult08vsDphi = 0;
  // hpTLvsRefMult08vsDphi = new TH3D("hpTLvsRefMult08vsDphi","p_{T}^{Leading} vs RefMult08 vs Dphi;#it{p}_{T}^{Leading} (GeV/c);RefMult08; #Delta#phi (rad)",nPtBins,PtBins,nNchBins,NchBins,nDeltabins,Deltabins);
  hpTLvsRefMult08vsDphi = new TH3D("hpTLvsRefMult08vsDphi","p_{T}^{Leading} vs RefMult08 vs Dphi;#it{p}_{T}^{Leading} (GeV/c);RefMult08; #Delta#phi (rad)",400,0,200,200,0,200,64,0,2.0*TMath::Pi());
  fListOfObjects->Add(hpTLvsRefMult08vsDphi);
  
  //for pTLeading vs V0M mult vs Dphi...
  hpTLvsV0MmultvsDphi = 0;
  hpTLvsV0MmultvsDphi = new TH3D("hpTLvsV0MmultvsDphi","p_{T}^{Leading} vs V0Mmult vs Dphi;#it{p}_{T}^{Leading} (GeV/c);V0Mmult; #Delta#phi (rad)",400,0,200,200,0,200,64,0,2.0*TMath::Pi());
  fListOfObjects->Add(hpTLvsV0MmultvsDphi);
  
  //for pT...
  hpTvspTLvsRefMult08 = 0;
  hpTvspTLvsRefMult08 = new TH3D("hpTvspTLvsRefMult08","p_{T} vs p_{T}^{Leading} vs RefMult08;#it{p}_{T} (GeV/c);#it{p}_{T}^{Leading} (GeV/c);RefMult08",nPtBins,PtBins,nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvspTLvsRefMult08);

  hpTvspTLvsRefMult08NS = 0;
  hpTvspTLvsRefMult08NS = new TH3D("hpTvspTLvsRefMult08NS","Near Side;#it{p}_{T} (GeV/c);#it{p}_{T}^{Leading} (GeV/c);RefMult08;",nPtBins,PtBins,nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvspTLvsRefMult08NS);
	
  hpTvspTLvsRefMult08AS = 0;
  hpTvspTLvsRefMult08AS = new TH3D("hpTvspTLvsRefMult08AS","Away Side;#it{p}_{T} (GeV/c);#it{p}_{T}^{Leading} (GeV/c);RefMult08",nPtBins,PtBins,nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvspTLvsRefMult08AS);

  hpTvspTLvsRefMult08TS = 0;
  hpTvspTLvsRefMult08TS = new TH3D("hpTvspTLvsRefMult08TS","Transverse Side;#it{p}_{T} (GeV/c);#it{p}_{T}^{Leading} (GeV/c);RefMult08",nPtBins,PtBins,nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvspTLvsRefMult08TS);

  hpTvspTLvsV0Mmult = 0;
  hpTvspTLvsV0Mmult = new TH3D("hpTvspTLvsV0Mmult","p_{T} vs p_{T}^{Leading} vs V0Mmult;#it{p}_{T} (GeV/c);#it{p}_{T}^{Leading} (GeV/c);V0Mmult",nPtBins,PtBins,nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvspTLvsV0Mmult);

  hpTvspTLvsV0MmultNS = 0;
  hpTvspTLvsV0MmultNS = new TH3D("hpTvspTLvsV0MmultNS","Near Side;#it{p}_{T} (GeV/c);#it{p}_{T}^{Leading} (GeV/c);V0Mmult",nPtBins,PtBins,nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvspTLvsV0MmultNS);

  hpTvspTLvsV0MmultAS = 0;
  hpTvspTLvsV0MmultAS = new TH3D("hpTvspTLvsV0MmultAS","Away Side;#it{p}_{T} (GeV/c);#it{p}_{T}^{Leading} (GeV/c);V0Mmult",nPtBins,PtBins,nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvspTLvsV0MmultAS);

  hpTvspTLvsV0MmultTS = 0;
  hpTvspTLvsV0MmultTS = new TH3D("hpTvspTLvsV0MmultTS","Transverse Side;#it{p}_{T} (GeV/c);#it{p}_{T}^{Leading} (GeV/c);V0Mmult",nPtBins,PtBins,nPtBins,PtBins,nNchBins,NchBins);
  fListOfObjects->Add(hpTvspTLvsV0MmultTS);
	
  hSumptVsRefMult08 = 0;
  hSumptVsRefMult08= new TH2D("hSumptVsRefMult08","#Sigma p_{T} vs RefMult08;RefMult08;#Sigma p_{T} (GeV/c)",nNchBins,NchBins,1000,0,500);
  fListOfObjects->Add(hSumptVsRefMult08);

    //hists for RT....
  hINEL0 = 0;
  hINEL0 = new TH1D( "hINEL0", ";INEL Evt. Step;Count",20,0,20);
  fListOfObjects->Add(hINEL0);
  
  hTrackRT = 0;
  hTrackRT = new TH1D( "hTrackRT", "#it{R}_{T};#it{R}_{T};count", 20,0,20);
  fListOfObjects->Add(hTrackRT);

  hTrackRTvsRefMult08 = 0;
  hTrackRTvsRefMult08 = new TH2D("hTrackRTvsRefMult08","#it{R}_{T} vs N_{ch};#it{R}_{T};N_{ch} (|#eta| < 0.8 & p_{T} > 0.15 GeV/c)",200,0,20,nNchBins,NchBins);
  fListOfObjects->Add(hTrackRTvsRefMult08);

  hMultTSvsRefMult08  = 0;
  hMultTSvsRefMult08  = new TH2D("hMultTSvsRefMult08","N_{ch} vs N_{ch}^{TS} ;N_{ch}^{TS};N_{ch} (|#eta| < 0.8 & p_{T} > 0.15 GeV/c)",nNchBins,NchBins,nNchBins,NchBins);
   fListOfObjects->Add(hMultTSvsRefMult08);

   hSumptVsTrackRT = 0;
   hSumptVsTrackRT = new TH2D("hSumptVsTrackRT","#Sigma p_{T} vs #it{R}_{T}; #it{R}_{T};#Sigma p_{T} (GeV/c)",200,0,20,5000,0,500);
   fListOfObjects->Add(hSumptVsTrackRT);

   hpTLvsRefMult08vsRT = 0;
   hpTLvsRefMult08vsRT = new TH3D("hpTLvsRefMult08vsRT","p_{T}^{Leading} vs RefMult08 vs #it{R}_{T};#it{p}_{T}^{Leading} (GeV/c);RefMult08; #Delta#phi (rad)",nPtBins,PtBins,nNchBins,NchBins,nRTbins,RTbins);
   fListOfObjects->Add(hpTLvsRefMult08vsRT);

   hDeltaphiVspTLVsRT = 0;
   hDeltaphiVspTLVsRT  =new TH3D("hDeltaphiVspTLVsRT",";#Delta#phi(#phi_{Leading}-#phi) (rad);#it{p}_{T}^{leading} (GeV/c);#it{R}_{T}",64,-pi/2.0,3.0*pi/2.0,400,0,200,20,0,20);
   fListOfObjects->Add(hDeltaphiVspTLVsRT);
   
   hDeltaphiVsMultVsRT = 0;
   hDeltaphiVsMultVsRT  =new TH3D("hDeltaphiVsMultVsRT",";#Delta#phi(#phi_{Leading}-#phi) (rad);N_{ch} (|#eta| < 0.8 & p_{T} > 0.15 GeV/c) ;#it{R}_{T}",64,-pi/2.0,3.0*pi/2.0,200,0,200,20,0,20);
   fListOfObjects->Add(hDeltaphiVsMultVsRT);

  //for ptl>5 GeV/c...
   hINEL0pTL5 = 0;
   hINEL0pTL5 = new TH1D( "hINEL0pTL5", ";INEL Evt. Step;Count",20,0,20);
   fListOfObjects->Add(hINEL0pTL5);

   hTrackRTpTL5 = 0;
   hTrackRTpTL5 = new TH1D( "hTrackRTpTL5", "#it{R}_{T};#it{R}_{T};count", 20,0,20);
  fListOfObjects->Add(hTrackRTpTL5);

  hTrackRTvsRefMult08pTL5 = 0;
  hTrackRTvsRefMult08pTL5 = new TH2D("hTrackRTvsRefMult08pTL5","#it{R}_{T} vs N_{ch};#it{R}_{T};N_{ch} (|#eta| < 0.8 & p_{T} > 0.15 GeV/c)",200,0,20,nNchBins,NchBins);
  fListOfObjects->Add(hTrackRTvsRefMult08pTL5);

  hMultTSvsRefMult08pTL5  = 0;
  hMultTSvsRefMult08pTL5  = new TH2D("hMultTSvsRefMult08pTL5","N_{ch} vs N_{ch}^{TS} ;N_{ch}^{TS};N_{ch} (|#eta| < 0.8 & p_{T} > 0.15 GeV/c)",nNchBins,NchBins,nNchBins,NchBins);
   fListOfObjects->Add(hMultTSvsRefMult08pTL5);

   hSumptVsTrackRTpTL5 = 0;
   hSumptVsTrackRTpTL5 = new TH2D("hSumptVsTrackRTpTL5","#Sigma p_{T} vs #it{R}_{T}; #it{R}_{T};#Sigma p_{T} (GeV/c)",200,0,20,5000,0,500);
   fListOfObjects->Add(hSumptVsTrackRTpTL5);
   
   hpTLvsRefMult08vsRTpTL5 = 0;
   hpTLvsRefMult08vsRTpTL5 = new TH3D("hpTLvsRefMult08vsRTpTL5","p_{T}^{Leading} vs RefMult08 vs #it{R}_{T};#it{p}_{T}^{Leading} (GeV/c);RefMult08; #Delta#phi (rad)",nPtBins,PtBins,nNchBins,NchBins,nRTbins,RTbins);
   fListOfObjects->Add(hpTLvsRefMult08vsRTpTL5);

   hDeltaphiVspTLVsRTpTL5 = 0;
   hDeltaphiVspTLVsRTpTL5  =new TH3D("hDeltaphiVspTLVsRTpTL5",";#Delta#phi(#phi_{Leading}-#phi) (rad);#it{p}_{T}^{leading} (GeV/c);#it{R}_{T}",64,-pi/2.0,3.0*pi/2.0,400,0,200,20,0,20);
   fListOfObjects->Add(hDeltaphiVspTLVsRTpTL5);
   
   hDeltaphiVsMultVsRTpTL5 = 0;
   hDeltaphiVsMultVsRTpTL5  =new TH3D("hDeltaphiVsMultVsRTpTL5",";#Delta#phi(#phi_{Leading}-#phi) (rad);N_{ch} (|#eta| < 0.8 & p_{T} > 0.15 GeV/c) ;#it{R}_{T}",64,-pi/2.0,3.0*pi/2.0,200,0,200,20,0,20);
   fListOfObjects->Add(hDeltaphiVsMultVsRTpTL5);

  for ( int i = 0 ; i < 10 ; i++ ){
    hpTLRT[i] = 0;
    hpTLRT[i] = new TH1D(Form("hpTLRT%d",i),Form("#it{p}_{T}^{leading} (%s);#it{p}_{T}^{leading} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTLRT[i]);

    hpTRT[i] = 0;
    hpTRT[i] = new TH1D(Form("hpTRT%d",i),Form("#it{p}_{T} (%s);#it{p}_{T} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTRT[i]);

    hpTNSRT[i] = 0;
    hpTNSRT[i] = new TH1D(Form("hpTNSRT%d",i),Form("#it{p}_{T} for NS (%s);#it{p}_{T} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTNSRT[i]);

    hpTASRT[i] = 0;
    hpTASRT[i] = new TH1D(Form("hpTASRT%d",i),Form("#it{p}_{T} for AS (%s);#it{p}_{T} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTASRT[i]);

    hpTTSRT[i] = 0;
    hpTTSRT[i] = new TH1D(Form("hpTTSRT%d",i),Form("#it{p}_{T} for TS (%s);#it{p}_{T} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTTSRT[i]);
    
    hDEtaDPhiRT[i] = 0;
    hDEtaDPhiRT[i] = new TH2D(Form("hDEtaDPhiRT%d", i),Form("%s;#Delta#phi(rad);#Delta#eta",legTrackRT[i]),64,-pi/2.0,3.0*pi/2.0,20,-1,1);
    fListOfObjects->Add(hDEtaDPhiRT[i]);
     
    hDeltaphiVspTLRT[i] = 0;
    hDeltaphiVspTLRT[i]  =new TH2D(Form("hDeltaphiVspTLRT%d", i),Form("%s;#Delta#phi(#phi_{Leading}-#phi) (rad);#it{p}_{T}^{leading} (GeV/c)",legTrackRT[i]),64,-pi/2.0,3.0*pi/2.0,400,0,200);
    fListOfObjects->Add(hDeltaphiVspTLRT[i]);

    //for pTL > 5 GeV/c...
    hpTLRTpTL5[i] = 0;
    hpTLRTpTL5[i] = new TH1D(Form("hpTLRTpTL5%d",i),Form("#it{p}_{T}^{leading} (%s);#it{p}_{T}^{leading} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTLRTpTL5[i]);
    
    hpTRTpTL5[i] = 0;
    hpTRTpTL5[i] = new TH1D(Form("hpTRTpTL5%d",i),Form("#it{p}_{T} (%s);#it{p}_{T} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTRTpTL5[i]);
    
    hpTNSRTpTL5[i] = 0;
    hpTNSRTpTL5[i] = new TH1D(Form("hpTNSRTpTL5%d",i),Form("#it{p}_{T} for NS (%s);#it{p}_{T} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTNSRTpTL5[i]);

    hpTASRTpTL5[i] = 0;
    hpTASRTpTL5[i] = new TH1D(Form("hpTASRTpTL5%d",i),Form("#it{p}_{T} for AS (%s);#it{p}_{T} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTASRTpTL5[i]);

    hpTTSRTpTL5[i] = 0;
    hpTTSRTpTL5[i] = new TH1D(Form("hpTTSRTpTL5%d",i),Form("#it{p}_{T} for TS (%s);#it{p}_{T} (GeV/#it{c});count",legTrackRT[i]),nPtBins,PtBins);
    fListOfObjects->Add(hpTTSRTpTL5[i]);

    hDEtaDPhiRTpTL5[i] = 0;
    hDEtaDPhiRTpTL5[i] = new TH2D(Form("hDEtaDPhiRTpTL5%d", i),Form("%s;#Delta#phi(rad);#Delta#eta",legTrackRT[i]),64,-pi/2.0,3.0*pi/2.0,20,-1,1);
    fListOfObjects->Add(hDEtaDPhiRTpTL5[i]);
    
    hDeltaphiVspTLRTpTL5[i] = 0;
    hDeltaphiVspTLRTpTL5[i]  =new TH2D(Form("hDeltaphiVspTLRTpTL5%d", i),Form("%s;#Delta#phi(#phi_{Leading}-#phi) (rad);#it{p}_{T}^{leading} (GeV/c)",legTrackRT[i]),64,-pi/2.0,3.0*pi/2.0,400,0,200);
    fListOfObjects->Add(hDeltaphiVspTLRTpTL5[i]);
  }
  
  fEventCuts.AddQAplotsToList(fListOfObjects);

  PostData(1, fListOfObjects);

}

//______________________________________________________________________________
void AliAnalysisTaskUeMultDep::UserExec(Option_t *)
{

  // -----------------------------------------------------
  //			 InputEvent
  // -----------------------------------------------------

  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  // -----------------------------------------------------
  //			 E S D
  // -----------------------------------------------------

  if (fAnalysisType == "ESD"){
    fESD = dynamic_cast<AliESDEvent*>(event);

    if(!fESD){
      Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }
  }

  /************ BEGGINING OF EVENT SELECTION ************************************************************************************/
  
  // Get trigger decision
  fTriggeredEventMB = 0; //init    
  if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ftrigBit )
    {
      fTriggeredEventMB = 1;  //event triggered as minimum bias
    }
  
  // check for spd vs clusters background 
  Bool_t SPDvsClustersBG = kFALSE;
  
  AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
  if (!AnalysisUtils)
    {
      cout<<"------- No AnalysisUtils Object Found --------"<<AnalysisUtils<<endl;
      return;
    }
  else
    SPDvsClustersBG = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD); // We want NO BG
  
  // is pile up ? We want isNotPileUp
  Bool_t isNotPileUp = AliPPVsMultUtils::IsNotPileupSPDInMultBins( fESD ); // my analysis 
  //Bool_t isNotPileUp = !fESD->IsPileupFromSPD(5,0.8);                    // GSI analysis
  //Bool_t isNotPileUp = !fESD->IsPileupFromSPD();                         // Gyula's note 
  
  Bool_t IncompleteDAQ = fESD->IsIncompleteDAQ(); // we want is not incomplete DAQ
  // -------------------------------------- multiplcity estimators section ------------------------------------------ //

  ftrackmult08 = -999;
  fv0mpercentile = -999;

  //ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);     //reference
  ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, fEtaCut);     //tracklets
  //ftrackmult08 = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD); //Combined estimator

  hRefMult08->Fill(ftrackmult08);

  fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection"); // Esto es para 13 TeV
  if (!fMultSelection)
    cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
  else
    fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  hV0Mmult->Fill(fv0mpercentile);

  cout<<"------- V0M mult ==  "<<fv0mpercentile<<"--------"<<endl;
	
  hRefMultvsV0Mmult->Fill(ftrackmult08,fv0mpercentile);

  // ------------------------------------------ end of mult estimators -------------------------------------------------//

  // vertex bussines 
  //const AliESDVertex * spdVertex =    fESD->GetPrimaryVertexSPD();
  const AliESDVertex * vertex    =    fESD->GetPrimaryVertex(); // tracks vertex, if not -> spd vertex, if not TPC vertex
  
  Bool_t isVtxGood = vertex->GetStatus() &&
    selectVertex2015pp( fESD ,kTRUE,kFALSE,kTRUE); // requires Tracks and spd vertex, and Zconsistency of 5mm
  
  double vertex_z = vertex->GetZ();
  Bool_t isVtxInZCut = (TMath::Abs(vertex_z)   <= fVtxCut); // Zvtx in +- 10
  
  // physics selection		     
  fisTracklet = (AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 1.0) >= 1);


  	/****************** IS PHYSICS SELECTION FLAG *************************************************/
	fisPS = fTriggeredEventMB && 
			!IncompleteDAQ && 
			!SPDvsClustersBG && 
			isNotPileUp &&
			fisTracklet; // *tracklet now included in the physics selection !!!
		

	// recontructed INEL > 0 is PS + vtx + Zvtx inside +-10 ------
	if ( fisPS && isVtxGood && isVtxInZCut)
		isINEL0Rec = kTRUE;
	else
		isINEL0Rec = kFALSE;

	if (fisPS)
		fVtxBeforeCuts->Fill(vertex_z);
	if (isINEL0Rec)
		fVtxAfterCuts->Fill(vertex_z);

	fHistEventCounter->Fill(0.5); 	//	All events

	if (!fEventCuts.AcceptEvent(event)) {
	  PostData(1, fListOfObjects);
	  return;
	}

	if(fTriggeredEventMB)  fHistEventCounter->Fill(1.5); // triggered events
	
	if ( fTriggeredEventMB && !IncompleteDAQ ) fHistEventCounter->Fill(2.5); // trigger + IsIncompleteDAQ
	
	if (fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG) fHistEventCounter->Fill(3.5); // trigger + IsIncompleteDAQ + BG rejection

	if(fPileUpRej)
	{
		if(fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp) fHistEventCounter->Fill(4.5); // trigger + IsIncompleteDAQ + BG rejection + PileUp
	}
	if (fisPS) fHistEventCounter->Fill(5.5); //PS: trigger + IsIncompleteDAQ + BG rejection + PileUp + 1 Tracklet in eta +-1
	if (fisPS && isVtxGood)	fHistEventCounter->Fill(6.5); //PS + GetPrimaryVertex
	if (isINEL0Rec) fHistEventCounter->Fill(7.5); //PS + GetPrimaryVertex + isVtxInZCut

	if (isINEL0Rec) DataAnaMult(fEtaCut);

	
	// Post output data.
	PostData(1, fListOfObjects);
	
}
//________________________________________________________________________
void AliAnalysisTaskUeMultDep::DataAnaMult( Double_t etaCut ){


  // selection on leading particle
  Double_t pt_leading    = 0;
  Double_t p_leading    = 0;
  Double_t eta_leading    = 0;
  Double_t phi_leading    = 0;

  Int_t    i_leading = 0;
 
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

    AliESDtrack* esdTrack = fESD->GetTrack(i);

    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t momentum = esdTrack->P();
    Double_t pt       = esdTrack->Pt();

    hPhi->Fill(phi);
	
    if(TMath::Abs(eta) > etaCut) continue;

    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack)) continue;

    if(pt<0.15) continue;

    Short_t ncl = esdTrack->GetTPCsignalN();
    if(ncl<fNcl) continue;
	
    if(pt>pt_leading){
      pt_leading      = pt;
      p_leading       = momentum;
      eta_leading     = eta;
      phi_leading     = phi;
      i_leading = i;
    }

    hpT->Fill(pt);

  }// end loop over tracks

  if(pt_leading<0.15) return;

  hPtL->Fill(pt_leading);
  hEtaL->Fill(eta_leading);
  hPhiL->Fill(phi_leading);

  // Next step: pTL vs Number density (NS, AS, TS)
  Double_t mult_ns = 0;
  Double_t mult_nns = 0;
  Double_t mult_as = 0;
  Double_t mult_ts = 0;
  Double_t mult    = 0;

  Double_t TrackRT=0;
  Double_t sumpt = 0;
  
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

    // exclude the auto-correlation
   if(i==i_leading) continue;

    AliESDtrack* esdTrack = fESD->GetTrack(i);

    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t pt       = esdTrack->Pt();

    if(TMath::Abs(eta) > etaCut)
      continue;

    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack))
      continue;

    if(pt<0.15)// only above 500 GeV/c
      continue;
    mult++;
    sumpt+= pt;
    
    Double_t DPhi = DeltaPhi( phi, phi_leading );
    Double_t DEta = TMath::Abs( eta -  eta_leading );
    Double_t R = TMath::Sqrt(DPhi*DPhi+DEta*DEta);
    Double_t DeltaEta = eta -  eta_leading;

    hDphi->Fill(DPhi);
    hpTvsRefMult08->Fill(pt,ftrackmult08);
    hpTvsV0Mmult->Fill(pt,fv0mpercentile);
    hpTLvsRefMult08vsDphi->Fill(pt_leading,ftrackmult08,DPhi);
    hpTLvsV0MmultvsDphi->Fill(pt_leading,fv0mpercentile,DPhi);
    hpTvspTLvsRefMult08->Fill(pt,pt_leading,ftrackmult08);
    hpTvspTLvsV0Mmult->Fill(pt,pt_leading,fv0mpercentile);
	
    // near side
    if(TMath::Abs(DPhi)<pi/3.0){ mult_ns++;
      hpTvspTLvsRefMult08NS->Fill(pt,pt_leading,ftrackmult08);
      hpTvspTLvsV0MmultNS->Fill(pt,pt_leading,fv0mpercentile);
    }
    else if(TMath::Abs(DPhi-pi)<pi/3.0){ mult_as++;
      hpTvspTLvsRefMult08AS->Fill(pt,pt_leading,ftrackmult08);
      hpTvspTLvsV0MmultAS->Fill(pt,pt_leading,fv0mpercentile);
    }
    else{ mult_ts++;
      hpTvspTLvsRefMult08TS->Fill(pt,pt_leading,ftrackmult08);
      hpTvspTLvsV0MmultTS->Fill(pt,pt_leading,fv0mpercentile);
    }  
  }// end loop over tracks
  hpTLvsRefMult08->Fill(pt_leading,ftrackmult08);
  hpTLvsV0Mmult->Fill(pt_leading,fv0mpercentile);
 
  hSumptVsRefMult08->Fill(ftrackmult08,sumpt);
        		
  // areas
  Double_t total_a = 2*TMath::Pi()*1.6;
  Double_t ns_a = 2*(TMath::Pi()/3.0)*1.6;
  Double_t as_a = 2*(TMath::Pi()/3.0)*1.6;
  Double_t ts_a = total_a-(ns_a+as_a);

  ProfpTLvsNch->Fill(pt_leading,mult/total_a);
  ProfpTLvsNchNS->Fill(pt_leading,mult_ns/ns_a);
  ProfpTLvsNchAS->Fill(pt_leading,mult_as/as_a);
  ProfpTLvsNchTS->Fill(pt_leading,mult_ts/ts_a);

  //*********************************
  // Calculations for RT
  //*********************************
  TrackRT = mult_ts/4.939;
  //  cout<<"TrackRT == "<<TrackRT<<endl;
  
  const Int_t nRTbins = 15;
  Double_t RTbins[nRTbins+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  const Double_t RTmin[nRTbins+1]={0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,  10.0, 11.0, 12.0, 13.0};
  const Double_t RTmax[nRTbins+1]={100.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0};

  Double_t sumptRT=0;
  Int_t nchRT = 0;
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
    
    AliESDtrack* esdTrack = fESD->GetTrack(i);
    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t pt       = esdTrack->Pt();
    
    if(TMath::Abs(eta) > etaCut) continue;
    
    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack)) continue;
    if(pt<0.15) continue;
    if(i==i_leading) continue;
    sumptRT+= pt;
    nchRT++;
  }
  
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
    
    AliESDtrack* esdTrack = fESD->GetTrack(i);
    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t pt       = esdTrack->Pt();
    
    if(TMath::Abs(eta) > etaCut) continue;
    
    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack)) continue;
    if(pt<0.15) continue;
    if(i==i_leading) continue;
    
    Double_t DPhi = DeltaPhi( phi, phi_leading );
    Double_t DEta = TMath::Abs( eta -  eta_leading );
    Double_t R = TMath::Sqrt(DPhi*DPhi+DEta*DEta);
    Double_t DeltaEta = eta -  eta_leading;
    
    hDeltaphiVspTLVsRT->Fill(DPhi,pt_leading,TrackRT);
    hDeltaphiVsMultVsRT->Fill(DPhi,nchRT,TrackRT);
    
    for(Int_t i=0;i<nRTbins;i++){// Loop for MPI STARTS
      if(TrackRT >= RTmin[i] && TrackRT < RTmax[i]) {
	hpTRT[i]->Fill(pt);
	hDEtaDPhiRT[i]->Fill(DPhi,DeltaEta);
	hDeltaphiVspTLRT[i]->Fill(DPhi,pt_leading);
      }
    }
    
    if(TMath::Abs(DPhi)<pi/3.0){
      for(Int_t i=0;i<nRTbins;i++){// Loop for MPI STARTS
	if(TrackRT >= RTmin[i] && TrackRT < RTmax[i]) {
	  hpTNSRT[i]->Fill(pt);
	}
      }
    }
    else if(TMath::Abs(DPhi-pi)<pi/3.0){
      for(Int_t i=0;i<nRTbins;i++){// Loop for MPI STARTS
	if(TrackRT >= RTmin[i] && TrackRT < RTmax[i]) {
	  hpTASRT[i]->Fill(pt);
	}
      } 
    }
    else{
      for(Int_t i=0;i<nRTbins;i++){// Loop for MPI STARTS
	if(TrackRT >= RTmin[i] && TrackRT < RTmax[i]) {
	  hpTTSRT[i]->Fill(pt);
	}
      }
    }
  }// end loop over tracks
 
  hTrackRT->Fill(TrackRT);
  hTrackRTvsRefMult08->Fill(TrackRT,ftrackmult08);
  hMultTSvsRefMult08->Fill(mult_ts,ftrackmult08);
  hSumptVsTrackRT->Fill(TrackRT,sumptRT);
  hpTLvsRefMult08vsRT->Fill(pt_leading,ftrackmult08,TrackRT);

  for(Int_t i=0;i<nRTbins;i++){
    if(TrackRT >= RTmin[i] && TrackRT< RTmax[i]){
      hINEL0->Fill(i);
      hpTLRT[i]->Fill(pt_leading);
    }
  }

  //
  //*********************************
  // Calculations for RT For pTL > 5 GeV/c......
  //*********************************
  if(pt_leading<5.0) return;
  //  cout<<"TrackRT == "<<TrackRT<<endl;
  Double_t mult_ts5 = 0;

  Double_t sumptRTpTL5=0;
  Int_t nchRTpTL5 = 0;
 for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

    AliESDtrack* esdTrack = fESD->GetTrack(i);
    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t pt       = esdTrack->Pt();

    if(TMath::Abs(eta) > etaCut) continue;

    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack)) continue;
    if(pt<0.15) continue;
    if(i==i_leading) continue;
    nchRTpTL5++;    
    sumptRTpTL5+= pt;
    
    Double_t DPhi = DeltaPhi( phi, phi_leading );
    Double_t DEta = TMath::Abs( eta -  eta_leading );
    Double_t R = TMath::Sqrt(DPhi*DPhi+DEta*DEta);
    Double_t DeltaEta = eta -  eta_leading;

   
    if(TMath::Abs(DPhi)<pi/3.0){}
    else if(TMath::Abs(DPhi-pi)<pi/3.0){}
    else mult_ts5++;
 }
 //RT calculations for PTL>5 GeV/c....
 Double_t TrackRTpTL5=0;
 TrackRTpTL5 = mult_ts5/4.939;
  
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

    AliESDtrack* esdTrack = fESD->GetTrack(i);
    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t pt       = esdTrack->Pt();

    if(TMath::Abs(eta) > etaCut) continue;
    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack)) continue;
    if(pt<0.15) continue;
		
    Double_t DPhi = DeltaPhi( phi, phi_leading );
    Double_t DEta = TMath::Abs( eta -  eta_leading );
    Double_t R = TMath::Sqrt(DPhi*DPhi+DEta*DEta);
    Double_t DeltaEta = eta -  eta_leading;

    hDeltaphiVspTLVsRTpTL5->Fill(DPhi,pt_leading,TrackRTpTL5);
    hDeltaphiVsMultVsRTpTL5->Fill(DPhi,nchRTpTL5,TrackRTpTL5);

    for(Int_t i=0;i<nRTbins;i++){// Loop for MPI STARTS
      if(TrackRTpTL5 >= RTmin[i] && TrackRTpTL5 < RTmax[i]) {
	hpTRTpTL5[i]->Fill(pt);
	hDEtaDPhiRTpTL5[i]->Fill(DPhi,DeltaEta);
	hDeltaphiVspTLRTpTL5[i]->Fill(DPhi,pt_leading);
      }
    }
    
    if(TMath::Abs(DPhi)<pi/3.0){
      for(Int_t i=0;i<nRTbins;i++){// Loop for MPI STARTS
	if(TrackRTpTL5 >= RTmin[i] && TrackRTpTL5 < RTmax[i]) {
	  hpTNSRTpTL5[i]->Fill(pt);
	}
      }
    }
    else if(TMath::Abs(DPhi-pi)<pi/3.0){
      for(Int_t i=0;i<nRTbins;i++){// Loop for MPI STARTS
	if(TrackRTpTL5 >= RTmin[i] && TrackRTpTL5 < RTmax[i]) {
	  hpTASRTpTL5[i]->Fill(pt);
	}
      } 
    }
    else{
      for(Int_t i=0;i<nRTbins;i++){// Loop for MPI STARTS
	if(TrackRTpTL5 >= RTmin[i] && TrackRTpTL5 < RTmax[i]) {
	  hpTTSRTpTL5[i]->Fill(pt);
	}
      }
    }
  }// end loop over tracks

  for(Int_t i=0;i<nRTbins;i++){
    if(TrackRTpTL5 >= RTmin[i] && TrackRTpTL5< RTmax[i]){
      hINEL0pTL5->Fill(i);
      hpTLRTpTL5[i]->Fill(pt_leading);
    }
  }
  hTrackRTpTL5->Fill(TrackRTpTL5);
  hTrackRTvsRefMult08pTL5->Fill(TrackRTpTL5,ftrackmult08);
  hMultTSvsRefMult08pTL5->Fill(mult_ts5,ftrackmult08);
  hSumptVsTrackRTpTL5->Fill(TrackRTpTL5,sumptRTpTL5);
  hpTLvsRefMult08vsRTpTL5->Fill(pt_leading,ftrackmult08,TrackRTpTL5);
  
}
//____________________________________________________________
void AliAnalysisTaskUeMultDep::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  // TPC
  esdTrackCuts->SetMinNCrossedRowsTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(fEtaCut);

  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  // 7*(0.0015+0.0050/pt^1.1)
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");

  esdTrackCuts->SetMaxDCAToVertexZ(2);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);

  esdTrackCuts->SetMaxChi2PerClusterITS(36);

  /*
    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);//
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
    esdTrackCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);//
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);//
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);//
    esdTrackCuts->SetRequireTPCRefit(kTRUE);//
    esdTrackCuts->SetRequireITSRefit(kTRUE);//
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
    AliESDtrackCuts::kAny);//
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//
    esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);//
    esdTrackCuts->SetMaxDCAToVertexZ(2);//
    esdTrackCuts->SetDCAToVertex2D(kFALSE);//
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);//
    esdTrackCuts->SetMaxChi2PerClusterITS(36);//
  */
  fTrackFilter->AddCuts(esdTrackCuts);

}

Double_t AliAnalysisTaskUeMultDep::DeltaPhi(Double_t phia, Double_t phib,
					    Double_t rangeMin, Double_t rangeMax)
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();

  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += 2*pi;
  else if (dphi > rangeMax) dphi -= 2*pi;

  return dphi;
}

//______________________________________________________________________
Bool_t AliAnalysisTaskUeMultDep::selectVertex2015pp(AliESDEvent *esd,
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
return kTRUE;
}

//______________________________________________________________________
Bool_t AliAnalysisTaskUeMultDep::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}
//______________________________________________________________________
