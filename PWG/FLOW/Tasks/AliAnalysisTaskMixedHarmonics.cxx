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

/**************************************
 * analysis task for mixed harmomics  * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Raimond Snellings         *
 *           (snelling@nikhef.nl)     * 
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/
 
class TFile;
class TString;
class TList;
class AliAnalysisTaskSE; 

#include "Riostream.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskMixedHarmonics.h"
#include "AliFlowAnalysisWithMixedHarmonics.h"

#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliVVertex.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "TMatrixDSym.h"


using std::cout;
using std::endl;
ClassImp(AliAnalysisTaskMixedHarmonics)

//================================================================================================================

AliAnalysisTaskMixedHarmonics::AliAnalysisTaskMixedHarmonics(const char *name, Bool_t useParticleWeights): 
AliAnalysisTaskSE(name), 
fMultSelection(NULL),
fAnalysisUtil(NULL),
fEvent(NULL),
fMH(NULL), 
fListHistos(NULL),
fHarmonic(1),
fNoOfMultipicityBins(100),
fMultipicityBinWidth(1.),
fMinMultiplicity(3.),
fOppositeChargesPOI(kFALSE),
fEvaluateDifferential3pCorrelator(kFALSE),
fCorrectForDetectorEffects(kFALSE),
fPrintOnTheScreen(kTRUE),
fCalculateVsM(kFALSE),
fShowBinLabelsVsM(kFALSE),
fUseParticleWeights(useParticleWeights),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fRejectPileUp(kFALSE),
fRejectPileUpTight(kFALSE),
fFillQAHistograms(kFALSE),
fTPCvsGlobalTrkBefore(NULL),
fTPCvsGlobalTrkAfter(NULL),
fTPCvsESDTrk(NULL),
fWeightsList(NULL),
fCachedRunNum(0),
fListV0MCorr(NULL),
fHCorrectV0ChWeghts(NULL),
fHCorrectQNxV0C(NULL),
fHCorrectQNyV0C(NULL),
fHCorrectQNxV0A(NULL),
fHCorrectQNyV0A(NULL),
fHCorrectQ3xV0C(NULL),
fHCorrectQ3yV0C(NULL),
fHCorrectQ3xV0A(NULL),
fHCorrectQ3yV0A(NULL),   
fListZDCCorr(NULL),
fHZDCCparameters(NULL),
fHZDCAparameters(NULL),
fIs2018Data(kFALSE),
fV0CutPU(NULL),
fSPDCutPU(NULL),
fMultCutPU(NULL),
fCenCutLowPU(NULL),
fCenCutHighPU(NULL),
fHistTPConlyVsCL1Before(NULL),
fHistTPConlyVsCL1After(NULL), 
fHistTPConlyVsV0MBefore(NULL),
fHistTPConlyVsV0MAfter(NULL), 
fHistCentCL0VsV0MBefore(NULL),
fHistCentCL0VsV0MAfter(NULL), 
fHistTPCVsESDTrkBefore(NULL), 
fHistTPCVsESDTrkAfter(NULL)
{
 // constructor
 cout<<"AliAnalysisTaskMixedHarmonics::AliAnalysisTaskMixedHarmonics(const char *name, Bool_t useParticleWeights)"<<endl;
 
 // Define input and output slots here
 // Input slot #0 works with an AliFlowEventSimple
 DefineInput(0, AliFlowEventSimple::Class());  
 // Input slot #1 is needed for the weights input file:
 if(useParticleWeights)
 {
  DefineInput(1, TList::Class());   
 }  
 // Output slot #0 is reserved              
 // Output slot #1 writes into a TList container
 DefineOutput(1, TList::Class());  
}

AliAnalysisTaskMixedHarmonics::AliAnalysisTaskMixedHarmonics(): 
AliAnalysisTaskSE(),
fMultSelection(NULL),
fAnalysisUtil(NULL),
fEvent(NULL),
fMH(NULL),
fListHistos(NULL),
fHarmonic(0),
fNoOfMultipicityBins(0),
fMultipicityBinWidth(0),
fMinMultiplicity(0),
fOppositeChargesPOI(kFALSE),
fEvaluateDifferential3pCorrelator(kFALSE),
fCorrectForDetectorEffects(kFALSE),
fPrintOnTheScreen(kFALSE),
fCalculateVsM(kFALSE),
fShowBinLabelsVsM(kFALSE),
fUseParticleWeights(kFALSE),
fUsePhiWeights(kFALSE),
fUsePtWeights(kFALSE),
fUseEtaWeights(kFALSE),
fRejectPileUp(kFALSE),
fRejectPileUpTight(kFALSE),
fFillQAHistograms(kFALSE),
fTPCvsGlobalTrkBefore(NULL),
fTPCvsGlobalTrkAfter(NULL),
fTPCvsESDTrk(NULL),
fWeightsList(NULL),
fCachedRunNum(0),
fListV0MCorr(NULL),
fHCorrectV0ChWeghts(NULL),
fHCorrectQNxV0C(NULL),
fHCorrectQNyV0C(NULL),
fHCorrectQNxV0A(NULL),
fHCorrectQNyV0A(NULL),
fHCorrectQ3xV0C(NULL),
fHCorrectQ3yV0C(NULL),
fHCorrectQ3xV0A(NULL),
fHCorrectQ3yV0A(NULL),   
fListZDCCorr(NULL),
fHZDCCparameters(NULL),
fHZDCAparameters(NULL),
fIs2018Data(kFALSE),
fV0CutPU(NULL),
fSPDCutPU(NULL),
fMultCutPU(NULL),
fCenCutLowPU(NULL),
fCenCutHighPU(NULL),
fHistTPConlyVsCL1Before(NULL),
fHistTPConlyVsCL1After(NULL), 
fHistTPConlyVsV0MBefore(NULL),
fHistTPConlyVsV0MAfter(NULL), 
fHistCentCL0VsV0MBefore(NULL),
fHistCentCL0VsV0MAfter(NULL), 
fHistTPCVsESDTrkBefore(NULL), 
fHistTPCVsESDTrkAfter(NULL)
{
 // Dummy constructor
 cout<<"AliAnalysisTaskMixedHarmonics::AliAnalysisTaskMixedHarmonics()"<<endl;
}

//================================================================================================================

void AliAnalysisTaskMixedHarmonics::UserCreateOutputObjects() 
{

 fAnalysisUtil = new AliAnalysisUtils();
 fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
 fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);


 // Called at every worker node to initialize
 cout<<"AliAnalysisTaskMixedHarmonics::UserCreateOutputObjects()"<<endl;

 // Analyser:
 fMH = new AliFlowAnalysisWithMixedHarmonics();
  
 // Common:
 fMH->SetHarmonic(fHarmonic);
 fMH->SetNoOfMultipicityBins(fNoOfMultipicityBins);
 fMH->SetMultipicityBinWidth(fMultipicityBinWidth);
 fMH->SetMinMultiplicity(fMinMultiplicity);
 fMH->SetOppositeChargesPOI(fOppositeChargesPOI);
 fMH->SetEvaluateDifferential3pCorrelator(fEvaluateDifferential3pCorrelator); 
 fMH->SetCorrectForDetectorEffects(fCorrectForDetectorEffects);
 fMH->SetPrintOnTheScreen(fPrintOnTheScreen); 
 fMH->SetCalculateVsM(fCalculateVsM); 
 fMH->SetShowBinLabelsVsM(fShowBinLabelsVsM);
 if(fUseParticleWeights)
 {
  // Pass the flags to class:
  if(fUsePhiWeights) fMH->SetUsePhiWeights(fUsePhiWeights);
  if(fUsePtWeights) fMH->SetUsePtWeights(fUsePtWeights);
  if(fUseEtaWeights) fMH->SetUseEtaWeights(fUseEtaWeights);
  // Get data from input slot #1 which is used for weights:
  if(GetNinputs()==2) 
  {                   
   fWeightsList = (TList*)GetInputData(1); 
  }
  // Pass the list with weights to class:
  if(fWeightsList) fMH->SetWeightsList(fWeightsList);
 }
 
 fMH->Init();
 
 if(fMH->GetHistList()) 
 {
  fListHistos = fMH->GetHistList();
  // fListHistos->Print();
 } else 
   {
    Printf("ERROR: Could not retrieve histogram list (MH, Task::UserCreateOutputObjects()) !!!!"); 
   }
 

 //Add QA histograms:
  fTPCvsGlobalTrkBefore = new TH2F("fTPCvsGlobalTrkBefore","Global(Fb32) vs TPC(FB128)",250,0,5000,250,0,5000);
  fListHistos->Add(fTPCvsGlobalTrkBefore);
  fTPCvsGlobalTrkAfter = new TH2F("fTPCvsGlobalTrkAfter","Global(Fb32) vs TPC(FB128)",250,0,5000,250,0,5000);
  fListHistos->Add(fTPCvsGlobalTrkAfter);


  fTPCvsESDTrk = new TH2F("fTPCvsESDTrk","ESDTrk vs TPC(FB128)",1000,0,20000,250,0,5000);
  fListHistos->Add(fTPCvsESDTrk);

  // Centrality Correlations and PileUp QA
  fHistTPConlyVsCL1Before = new TH2F("fHistTPConlyVsCL1Before","Before;Cent(CL1); TPC(FB128)",100,0,100,250,0,5000);
  fListHistos->Add(fHistTPConlyVsCL1Before);
  fHistTPConlyVsCL1After  = new TH2F("fHistTPConlyVsCL1After","After; Cent(CL1); TPC(FB128) ",100,0,100,250,0,5000);
  fListHistos->Add(fHistTPConlyVsCL1After);

  fHistTPConlyVsV0MBefore = new TH2F("fHistTPConlyVsV0MBefore","Before;Cent(V0M); TPC(FB128)",100,0,100,250,0,5000);
  fListHistos->Add(fHistTPConlyVsV0MBefore);
  fHistTPConlyVsV0MAfter  = new TH2F("fHistTPConlyVsV0MAfter","After; Cent(V0M); TPC(FB128) ",100,0,100,250,0,5000);
  fListHistos->Add(fHistTPConlyVsV0MAfter);

  fHistCentCL0VsV0MBefore = new TH2F("fHistCentCL0VsV0MBefore","Before;Cent(V0M); Cent(CL0)",100,0,100,100,0,100);
  fListHistos->Add(fHistCentCL0VsV0MBefore);
  fHistCentCL0VsV0MAfter  = new TH2F("fHistCentCL0VsV0MAfter"," After; Cent(V0M); Cent(CL0)",100,0,100,100,0,100);
  fListHistos->Add(fHistCentCL0VsV0MAfter);

  fHistTPCVsESDTrkBefore = new TH2F("fHistTPCVsESDTrkBefore","Before; TPC1; ESD trk",100,0,5000,200,0,20000);
  fListHistos->Add(fHistTPCVsESDTrkBefore);
  fHistTPCVsESDTrkAfter  = new TH2F("fHistTPCVsESDTrkAfter"," After;  TPC1; ESD trk",100,0,5000,200,0,20000);
  fListHistos->Add(fHistTPCVsESDTrkAfter);

 PostData(1,fListHistos);
  
} // end of void AliAnalysisTaskMixedHarmonics::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskMixedHarmonics::UserExec(Option_t *) 
{


 AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());

 Bool_t kPileupEvent = kFALSE;

 //cout<<" Run = "<<aodEvent->GetRunNumber();

 if(fRejectPileUp || fFillQAHistograms){
   if (!fIs2018Data) 
     kPileupEvent = CheckEventIsPileUp(aodEvent);
   else if (fIs2018Data)
     kPileupEvent = CheckEventIsPileUp2018(aodEvent);
 }

 if(fRejectPileUp && kPileupEvent)  return;

 // main loop (called for each event)
 fEvent = dynamic_cast<AliFlowEventSimple*>(GetInputData(0));
 
 //cout<<" tracks = "<<fEvent->NumberOfTracks()<<endl;
 
 // ================================== V0 and ZDC Q-vector ============================
 // Read current run number
 Int_t runNumber = aodEvent->GetRunNumber();

 //-------------- Vtx ---------------
 const AliVVertex *pointVtx = aodEvent->GetPrimaryVertex();
 Double_t pVtxZ = -999, pVtxX = -999, pVtxY=-999;
 pVtxZ  = pointVtx->GetZ();
 pVtxX  = pointVtx->GetX();
 pVtxY  = pointVtx->GetY();
 
 //-------------- OrbitNum -----------
 UInt_t period = aodEvent->GetPeriodNumber();
 UInt_t orbit24 = aodEvent->GetOrbitNumber();
  
 if (period > 255) { // 8 bits
  cout<<"invalid period number"<<endl;
  period = 255;
  orbit24 = (1<<24)-1;
 }
    
 if (orbit24 >= (1<<24)) { // 24 bits
  cout<<"invalid orbit number"<<endl;
  period = 255;
  orbit24 = (1<<24)-1;
 }
  
 UInt_t orbit = period * (1<<24) + orbit24;

 Double_t fOrbitNumber = static_cast<double>(orbit)/1000000.; // scale down by 10^6 to fit to the scale used in least square fit. In least square fit, orbit number is scaled down by 10^6 to weight down the contribution to it. Maybe not necessary to scale down, but it should make the fit more stable in principle
  
  
 // Set Up Correction Map for this run:
 if(runNumber!=fCachedRunNum) {
  if(fListV0MCorr){
   GetV0MCorrectionHist(runNumber);
  }
  if(fListZDCCorr){
   GetZDCCorrectionHist(runNumber);
  }
  fCachedRunNum = runNumber;
 }
 //-------------- centrality --------
 Double_t centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
 Double_t centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");

 // Get V0 Q-vector
 Double_t fQ2xV0C=0, fQ2yV0C=0, fSumM2V0C=0, fQ2xV0A=0, fQ2yV0A=0, fSumM2V0A=0; 
 Bool_t kPassV0 = GetGainCorrectedV0Qvector(aodEvent, pVtxZ, 2, fQ2xV0C, fQ2yV0C, fSumM2V0C, fQ2xV0A, fQ2yV0A, fSumM2V0A); 

 Double_t fQ3xV0C=0, fQ3yV0C=0, fSumM3V0C=0, fQ3xV0A=0, fQ3yV0A=0, fSumM3V0A=0;
 kPassV0 = GetGainCorrectedV0Qvector(aodEvent, pVtxZ, 3, fQ3xV0C, fQ3yV0C, fSumM3V0C, fQ3xV0A, fQ3yV0A, fSumM3V0A);

 if(!kPassV0) return;           // V0 does not have signal for this event.  

 ApplyV0XqVectRecenter(centrCL1, 2, fQ2xV0C, fQ2yV0C, fQ2xV0A, fQ2yV0A);
 ApplyV0XqVectRecenter(centrCL1, 3, fQ3xV0C, fQ3yV0C, fQ3xV0A, fQ3yV0A);

 fEvent->SetV02Qsub(fQ2xV0C,fQ2yV0C,fSumM2V0C,fQ2xV0A,fQ2yV0A,fSumM2V0A,2); // 2nd order cos(2Psi_Vo) & sin(2Psi_Vo)
 fEvent->SetV02Qsub(fQ3xV0C,fQ3yV0C,fSumM3V0C,fQ3xV0A,fQ3yV0A,fSumM3V0A,3); // 3rd order cos(3Psi_Vo) & sin(3Psi_Vo)
 
 // Get ZDC Q-vector
 Double_t fQxZNCC=0, fQyZNCC=0, fQxZNCA=0, fQyZNCA=0; 
 Double_t denZNC=0, denZNA=0;
 
 Bool_t kPassZNC = GetGainCorrectedZNCQvector(aodEvent, fQxZNCC, fQyZNCC, denZNC, fQxZNCA, fQyZNCA, denZNA);

 if(!kPassZNC) return;           // ZNC does not have signal for this event.  

 ApplyZNCqVectRecenter(centrV0M, pVtxX, pVtxY, pVtxZ, fOrbitNumber, fQxZNCC, fQyZNCC, fQxZNCA, fQyZNCA);
 Double_t xyZNC[2], xyZNA[2];
 xyZNC[0] = fQxZNCC;
 xyZNC[1] = fQyZNCC;
 xyZNA[0] = fQxZNCA;
 xyZNA[1] = fQyZNCA;

 fEvent->SetZDC2Qsub(xyZNC,denZNC,xyZNA,denZNA); // 1st order ZNC qvec cos(Psi_ZNC) & sin(Psi_ZNC)
 
 // ================================== end of V0 and ZDC Q-vector ============================

 // Mixed Harmonics:
 if(fEvent) 
 {

  fMH->Make(fEvent);

 }
 else 
 {
  cout<<"WARNING: No input data (MH, Task::UserExec()) !!!!"<<endl;
  cout<<endl;
 }
  
 PostData(1,fListHistos);

}







//================================================================================================================

void AliAnalysisTaskMixedHarmonics::Terminate(Option_t *) 
{
 //accessing the merged output list: 
 fListHistos = (TList*)GetOutputData(1);
 
 fMH = new AliFlowAnalysisWithMixedHarmonics(); 
 
 if(fListHistos) 
 {
  fMH->GetOutputHistograms(fListHistos);
  fMH->Finish();
  PostData(1,fListHistos);
 } else
   {
    cout<<" WARNING: histogram list pointer is empty (MH, Task::Terminate()) !!!!"<<endl;
    cout<<endl;
   }
    
} // end of void AliAnalysisTaskMixedHarmonics::Terminate(Option_t *)




//----- PileUp removal function -------

Bool_t AliAnalysisTaskMixedHarmonics::CheckEventIsPileUp(AliAODEvent *faod) {

  Bool_t BisPileup=kFALSE;

  Double_t centrV0M=300;
  Double_t centrCL1=300;
  Double_t centrCL0=300;
  Double_t centrTRK=300;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
 
  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(1);
  }
  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
  centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");

  //-- pile-up a la Dobrin for LHC15o -----
  if(PileUpMultiVertex(faod)) {
    //fPileUpCount->Fill(0.5);
    BisPileup=kTRUE;
  }
  Int_t isPileup = faod->IsPileupFromSPD(3);
  if(isPileup != 0) {
    //fPileUpCount->Fill(1.5);
    BisPileup=kTRUE;          
  }
  if(((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
    //fPileUpCount->Fill(2.5);
    BisPileup=kTRUE;
  }
  if(faod->IsIncompleteDAQ())  {
    //fPileUpCount->Fill(3.5);
    BisPileup=kTRUE;
  }
  if(fabs(centrV0M-centrCL1)> 5.0)  {//default: 7.5
    //fPileUpCount->Fill(4.5);
    BisPileup=kTRUE;
  }

  // check vertex consistency
  const AliAODVertex* vtTrc = faod->GetPrimaryVertex();
  const AliAODVertex* vtSPD = faod->GetPrimaryVertexSPD();

  if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
    //fPileUpCount->Fill(5.5);
    BisPileup=kTRUE;
  }

  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);

  double dz = vtTrc->GetZ() - vtSPD->GetZ();

  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = dz/errTot;
  double nsigTrc = dz/errTrc;

  if(TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
    //fPileUpCount->Fill(6.5);
    BisPileup=kTRUE;
  }

  //cuts on tracks
  //Int_t multTrk = 0;
  //Int_t multTrkBefC = 0;
  //Int_t multTrkTOFBefC = 0;

  Int_t multTPC = 0;
  Int_t multITSfb96 = 0;
  Int_t multITSfb32 = 0;

  Int_t multTPCFE = 0;
  Int_t multGlobal = 0;
  Int_t multTPCuncut = 0;

  Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks();

  const Int_t nTracks = faod->GetNumberOfTracks();

  for(Int_t iTracks = 0; iTracks < nTracks; iTracks++) {
    //AliNanoAODTrack* track = dynamic_cast<AliNanoAODTrack*>(faod->GetTrack(iTracks));
    AliAODTrack* track = (AliAODTrack*)faod->GetTrack(iTracks);
    if(!track)  continue;
    //---------- old method -----------
    if(track->TestFilterBit(128))
      multTPC++;
    if(track->TestFilterBit(96))
      multITSfb96++;
    //----------------------------------
    if(track->TestFilterBit(1))  multTPCuncut++;
    if(track->TestFilterBit(32)) multITSfb32++;


    if(track->Pt()<0.2 || track->Pt()>5.0 || TMath::Abs(track->Eta())>0.8 || track->GetTPCNcls()<70 || track->GetTPCsignal()<10.0)
      continue;
    if(track->TestFilterBit(1) && track->Chi2perNDF()>0.2)  multTPCFE++;
    if(!track->TestFilterBit(16) || track->Chi2perNDF()<0.1)   continue;
                
    Double_t b[2]    = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};
                
    AliAODTrack copy(*track);
    Double_t magField = faod->GetMagneticField();
                
    if(magField!=0){     
      if(track->PropagateToDCA(faod->GetPrimaryVertex(), magField, 100., b, bCov) && TMath::Abs(b[0]) < 0.3 && TMath::Abs(b[1]) < 0.3) multGlobal++;    
    }
  }

  Double_t multTPCn      = multTPC;
  Double_t multEsdn      = multEsd;

  //fixed for test:
  Double_t fPileUpSlopeParm = 3.55;
  Double_t fPileUpConstParm = 90;

  Double_t multESDTPCDif  = multEsdn  - fPileUpSlopeParm*multTPCn;
  //Double_t multTPCGlobDif = multTPCFE - fPileUpSlopeParm*multGlobal;

  if(fFillQAHistograms){
    //fGlobalTracks->Fill(multGlobal);
    fTPCvsGlobalTrkBefore->Fill(multITSfb32,multTPC);
    fTPCvsESDTrk->Fill(multEsd,multTPC);
  }

  /*if(multESDTPCDif > (fRejectPileUpTight?700.:15000.)) {
  //fPileUpCount->Fill(7.5);
  BisPileup=kTRUE;
  }*/
  /*if(multESDTPCDif > 15000.) { //default: 15000
  //fPileUpCount->Fill(7.5);
  BisPileup=kTRUE;
  }*/


     
  if(fRejectPileUp) {
    if(multESDTPCDif > fPileUpConstParm) { 
      //fPileUpCount->Fill(7.5);
      BisPileup=kTRUE;
    }
    if(BisPileup==kFALSE) {
      if(!fMultSelection->GetThisEventIsNotPileup()) BisPileup=kTRUE;
      if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
      if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
      if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
      if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
      if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
      if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
      //if(BisPileup)     fPileUpCount->Fill(9.5);
    }  
  }
  
  if(!BisPileup){
    fTPCvsGlobalTrkAfter->Fill(multITSfb32,multTPC);
  }


  return BisPileup; 
} //-------pile up function ------








 Bool_t AliAnalysisTaskMixedHarmonics::PileUpMultiVertex(const AliAODEvent* faod)
 {  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if(!(nPlp=faod->GetNumberOfPileupVerticesTracks()))
  return kFALSE;

  vtPrm = faod->GetPrimaryVertex();
  if(vtPrm == faod->GetPrimaryVertexSPD())
  return kTRUE;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for(int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)faod->GetPileupVertexTracks(ipl);
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF()    > kMaxPlpChi2)    continue;
    //int bcPlp = vtPlp->GetBC();
    //if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2)
    // return kTRUE; // pile-up from other BC

    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist)        continue;

    return kTRUE; // pile-up: well separated vertices
    }
  return kFALSE;
}

void AliAnalysisTaskMixedHarmonics::SetupPileUpRemovalFunctions18qPass3() { //@Shi for 2018 period Pass3 data
	// 18q pass3
	fSPDCutPU = new TF1("fSPDCutPU", "480. + 3.95*x", 0, 50000);
   
    Double_t parV0[8] = {41.3226, 0.822835, 0.0880984, 206.961, 3.56337, 0.0965816, -0.00076483, 2.11591e-06};
    fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);
   
    Double_t parV0CL0[6] = {0.362458, 0.962768, 0.995134, 0.0331353, -0.000692428, 6.59962e-06};
    fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);
   
    Double_t parFB32[9] = {-812.555, 6.38397, 5379.01, -0.394814, 0.0296228, -26.1633, 317.365, -0.842175, 0.0165651};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCutPU->SetParameters(parFB32);
	
}

void AliAnalysisTaskMixedHarmonics::SetupPileUpRemovalFunctions18rPass3() { //@Shi for 2018 period Pass3 data
    // 18r pass3
    fSPDCutPU = new TF1("fSPDCutPU", "480. + 3.95*x", 0, 50000);
   
    Double_t parV0[8] = {42.4921, 0.823255, 0.0824939, 139.826, 7.27032, 0.0488425, -0.00045769, 1.40891e-06};
    fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
    fV0CutPU->SetParameters(parV0);
   
    Double_t parV0CL0[6] = {0.317973, 0.961823, 1.02383, 0.0330231, -0.000721551, 6.92564e-06};
    fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutLowPU->SetParameters(parV0CL0);
    fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fCenCutHighPU->SetParameters(parV0CL0);
   
    Double_t parFB32[9] = {-817.169, 6.40836, 5380.3, -0.394358, 0.0295209, -25.9573, 316.586, -0.843951, 0.0165442};
    fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 6.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCutPU->SetParameters(parFB32);
	
}

Bool_t AliAnalysisTaskMixedHarmonics::CheckEventIsPileUp2018(AliAODEvent *faod) {
  // check pile up for 2018 
  Bool_t BisPileup=kFALSE;

  Double_t centrV0M=-99.0;
  Double_t centrCL1=-99.0;
  Double_t centrCL0=-99.0;

  fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");

  if(!fMultSelection) {
    printf("\n\n **WARNING** ::UserExec() AliMultSelection object not found.\n\n");
    exit(111);
  }

  centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
  centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
  centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");

  Int_t nITSClsLy0 = faod->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = faod->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

  AliAODTracklets* aodTrkl = (AliAODTracklets*)faod->GetTracklets();
  Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

  const Int_t nTracks = faod->GetNumberOfTracks();

  Int_t multTrk = 0;

  for (Int_t it = 0; it < nTracks; it++) {
    
    AliAODTrack* aodTrk = (AliAODTrack*)faod->GetTrack(it);

    if (!aodTrk){
      delete aodTrk;
      continue;
    }

    if (aodTrk->TestFilterBit(32)){
      if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2))
        multTrk++;
    }
  }

  AliAODVZERO* aodV0 = faod->GetVZEROData();
  Float_t multV0a = aodV0->GetMTotV0A();
  Float_t multV0c = aodV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = aodV0->GetTriggerChargeA();
  UShort_t multV0cOn = aodV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;

  Int_t tpcClsTot = faod->GetNumberOfTPCClusters();
  Float_t nclsDif = Float_t(tpcClsTot) - (60932.9 + 69.2897*multV0Tot - 0.000217837*multV0Tot*multV0Tot);

  if (centrCL0 < fCenCutLowPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (centrCL0 > fCenCutHighPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (Float_t(nITSCls) > fSPDCutPU->Eval(nITSTrkls)) {
    BisPileup=kTRUE;
  }     
  if (multV0On < fV0CutPU->Eval(multV0Tot)) {
    BisPileup=kTRUE;
  }
  if (Float_t(multTrk) < fMultCutPU->Eval(centrV0M)) {
    BisPileup=kTRUE;
  }
  if (((AliAODHeader*)faod->GetHeader())->GetRefMultiplicityComb08() < 0) {
    BisPileup=kTRUE;
  }
  if (faod->IsIncompleteDAQ()) {
    BisPileup=kTRUE;
  }    
  //if (nclsDif > 200000)//can be increased to 200000
  // BisPileup=kTRUE;

  Int_t multEsd = ((AliAODHeader*)faod->GetHeader())->GetNumberOfESDTracks();

  if(fFillQAHistograms){
    fHistCentCL0VsV0MBefore->Fill(centrV0M,centrCL0);
    fHistTPCVsESDTrkBefore->Fill(multTrk,multEsd);  
    fHistTPConlyVsCL1Before->Fill(centrCL1,multTrk);
    fHistTPConlyVsV0MBefore->Fill(centrV0M,multTrk);
    
    
    if (!BisPileup) { 
      fHistCentCL0VsV0MAfter->Fill(centrV0M,centrCL0);
      fHistTPCVsESDTrkAfter->Fill(multTrk,multEsd);  
      fHistTPConlyVsCL1After->Fill(centrCL1,multTrk);
      fHistTPConlyVsV0MAfter->Fill(centrV0M,multTrk);      
    }
  }

  return BisPileup; 
}


double AliAnalysisTaskMixedHarmonics::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    AliDebug(2,"\n\n ::GetWDist => One of vertices is not valid\n\n");
    return 0;
  }
  static TMatrixDSym vVb(3);
  double dist = -1;
  double dx = v0->GetX()-v1->GetX();
  double dy = v0->GetY()-v1->GetY();
  double dz = v0->GetZ()-v1->GetZ();
  double cov0[6],cov1[6];
  v0->GetCovarianceMatrix(cov0);
  v1->GetCovarianceMatrix(cov1);
  vVb(0,0) = cov0[0]+cov1[0];
  vVb(1,1) = cov0[2]+cov1[2];
  vVb(2,2) = cov0[5]+cov1[5];
  vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
  vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
  vVb.InvertFast();
  if (!vVb.IsValid()) {
    AliDebug(2,"Singular Matrix\n");
    return dist;
  }
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
  +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;
}

void AliAnalysisTaskMixedHarmonics::GetV0MCorrectionHist(Int_t run){ 

  if(fListV0MCorr){
    //V0 Channel Gains:
    fHCorrectV0ChWeghts = (TH2F *) fListV0MCorr->FindObject(Form("hWgtV0ChannelsvsVzRun%d",run));
    if(fHCorrectV0ChWeghts){
      printf("\n ===========> Info:: V0 Channel Weights Found for Run %d \n ",run);
    }
    //Get V0A, V0C <Q> Vectors:
    fHCorrectQNxV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNxvsCentV0CRun%d",run));
    fHCorrectQNyV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNyvsCentV0CRun%d",run));    
    fHCorrectQNxV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNxvsCentV0ARun%d",run));
    fHCorrectQNyV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQNyvsCentV0ARun%d",run));
	
    fHCorrectQ3xV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3xvsCentV0CRun%d",run));
    fHCorrectQ3yV0C = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3yvsCentV0CRun%d",run));    
    fHCorrectQ3xV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3xvsCentV0ARun%d",run));
    fHCorrectQ3yV0A = (TH1D *) fListV0MCorr->FindObject(Form("fHisAvgQ3yvsCentV0ARun%d",run));    
    if(fHCorrectQNxV0C && fHCorrectQNyV0C && fHCorrectQNxV0A && fHCorrectQNyV0A){
      printf(" ===========> Info:: V0A,V0C <Q> Found for Run %d \n ",run);
    }    
  }
  else{
    fHCorrectV0ChWeghts=NULL;
  } 
}

void AliAnalysisTaskMixedHarmonics::GetZDCCorrectionHist(Int_t run){ 
  if(fListZDCCorr){
	fHZDCCparameters = (TH1D*)(fListZDCCorr->FindObject(Form("Run %d", run))->FindObject(Form("fZDCCparameters[%d]",run)));
	fHZDCAparameters = (TH1D*)(fListZDCCorr->FindObject(Form("Run %d", run))->FindObject(Form("fZDCAparameters[%d]",run)));
	if(fHZDCCparameters && fHZDCAparameters){
      printf("\n ===========> Info:: ZDC Channel Weights Found for Run %d \n ",run);
    }
  }
  else{
	fHZDCCparameters=NULL;
	fHZDCAparameters=NULL;
	printf("\n ===========> Info:: ZDC Channel Weights NOT Found for Run %d \n ",run);
  }
}

Bool_t AliAnalysisTaskMixedHarmonics::GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &sumV0C,Double_t &qnxV0A,Double_t &qnyV0A,Double_t &sumV0A){

  const AliAODVZERO *fAODV0 = (AliAODVZERO *) faod->GetVZEROData();
  Float_t fMultV0 = 0.;
  Float_t fPhiV0  = 0.;
  Float_t fV0chGain = 1.0;

  Double_t fQxV0CHarmN=0,fQyV0CHarmN=0,fQxV0AHarmN=0,fQyV0AHarmN=0;

  Double_t fSumMV0A = 0;
  Double_t fSumMV0C = 0;
  Int_t ibinV0=0;

  for(int iV0 = 0; iV0 < 64; iV0++) { //0-31 is V0C, 32-63 VOA

    fMultV0 = fAODV0->GetMultiplicity(iV0);

    /// V0 Channel Gain Correction:
    if(fHCorrectV0ChWeghts){ 
      ibinV0    = fHCorrectV0ChWeghts->FindBin(fVtxZ,iV0);
      fV0chGain = fHCorrectV0ChWeghts->GetBinContent(ibinV0); 
    }
    
    fMultV0 = fMultV0*fV0chGain;   //Corrected Multiplicity
    
    fPhiV0  = TMath::PiOver4()*(0.5 + iV0 % 8);

    if(iV0 < 32){
      qnxV0C   += TMath::Cos(gPsiN*fPhiV0) * fMultV0;
      qnyV0C   += TMath::Sin(gPsiN*fPhiV0) * fMultV0;
      fSumMV0C += fMultV0;
    }
    else if(iV0 >= 32){
      qnxV0A   += TMath::Cos(gPsiN*fPhiV0) * fMultV0;
      qnyV0A   += TMath::Sin(gPsiN*fPhiV0) * fMultV0;
      fSumMV0A += fMultV0;
    } 
  }///V0 Channel loop

  /// Now the q vectors:
  if(fSumMV0A<=1e-4 || fSumMV0C<=1e-4){
    qnxV0C = 0;
    qnyV0C = 0;
    sumV0C = 0;
    qnxV0A = 0;
    qnyV0A = 0;   
    sumV0A = 0; 
    return kFALSE;       
  }
  else{
    qnxV0C = qnxV0C/fSumMV0C;
    qnyV0C = qnyV0C/fSumMV0C;
    sumV0C = fSumMV0C;
    qnxV0A = qnxV0A/fSumMV0A;
    qnyV0A = qnyV0A/fSumMV0A;
    sumV0A = fSumMV0A;
    return kTRUE;  
  }
  
}

void AliAnalysisTaskMixedHarmonics::ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A){

  Int_t icentbin = 0;
  Double_t avgqx=0,avgqy=0; 
  //cout<<" => Before qnxV0C "<<qnxV0C<<"\tqnyV0C "<<qnyV0C<<"\tqnxV0A "<<qnxV0A<<"\tqnyV0A "<<qnyV0A<<endl;
  if(gPsiN==3){  ///<Q> correction for Psi3:  
    if(fHCorrectQ3xV0C && fHCorrectQ3yV0C){
      icentbin = fHCorrectQ3xV0C->FindBin(fCent);
      avgqx = fHCorrectQ3xV0C->GetBinContent(icentbin);
      avgqy = fHCorrectQ3yV0C->GetBinContent(icentbin);
      qnxV0C -= avgqx;
      qnyV0C -= avgqy;	
      //cout<<" V0C PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
    }
    if(fHCorrectQ3xV0A && fHCorrectQ3yV0A){
      icentbin = fHCorrectQ3xV0A->FindBin(fCent);
      avgqx = fHCorrectQ3xV0A->GetBinContent(icentbin);
      avgqy = fHCorrectQ3yV0A->GetBinContent(icentbin);
      qnxV0A -= avgqx;
      qnyV0A -= avgqy;
      //cout<<" V0A PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
    }
    //cout<<" => After qnxV0C "<<qnxV0C<<"\tqnyV0C "<<qnyV0C<<" qnxV0A "<<qnxV0A<<"\tqnyV0A"<<qnyV0A<<endl;
  }
  else{ /// Proper File which contain <q> for harmonic 'N' should be set in AddTask!! 
    if(fHCorrectQNxV0C && fHCorrectQNyV0C){
      icentbin = fHCorrectQNxV0C->FindBin(fCent);
      avgqx = fHCorrectQNxV0C->GetBinContent(icentbin);
      avgqy = fHCorrectQNyV0C->GetBinContent(icentbin);
      qnxV0C -= avgqx;
      qnyV0C -= avgqy;      
      //cout<<" V0C PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
    }
    if(fHCorrectQNxV0A && fHCorrectQNyV0A){
      icentbin = fHCorrectQNxV0A->FindBin(fCent);
      avgqx = fHCorrectQNxV0A->GetBinContent(icentbin);
      avgqy = fHCorrectQNyV0A->GetBinContent(icentbin);
      qnxV0A -= avgqx;
      qnyV0A -= avgqy;           
      //cout<<" V0A PsiN: "<<gPsiN<<" Cent: "<<fCent<<"\t <qx> "<<avgqx<<"\t <qy> "<<avgqy<<endl;
    }
    //cout<<" => After qnxV0C "<<qnxV0C<<"\tqnyV0C "<<qnyV0C<<" qnxV0A "<<qnxV0A<<"\tqnyV0A "<<qnyV0A<<endl;
  }

 
  return;

}

Bool_t AliAnalysisTaskMixedHarmonics::GetGainCorrectedZNCQvector(AliAODEvent *faod,Double_t &qnxZNCC,Double_t &qnyZNCC,Double_t &fdenZNC,Double_t &qnxZNCA,Double_t &qnyZNCA,Double_t &fdenZNA){
 AliAODZDC *aodZDC = faod->GetZDCData();
	
 if(!aodZDC) {
  printf("\n ********* Error: could not find ZDC data ************ \n ");
  return kFALSE;
 }
 else {
  const Double_t *fZNATowerRawAOD = aodZDC->GetZNATowerEnergy();
  const Double_t *fZNCTowerRawAOD = aodZDC->GetZNCTowerEnergy();
  if((fZNATowerRawAOD[0]<0) || (fZNATowerRawAOD[1]<0) || (fZNATowerRawAOD[2]<0) || (fZNATowerRawAOD[3]<0) || (fZNATowerRawAOD[4] < 0)) {
	return kFALSE;
  }
	
  if((fZNCTowerRawAOD[0]<0) || (fZNCTowerRawAOD[1]<0) || (fZNCTowerRawAOD[2]<0) || (fZNCTowerRawAOD[3]<0) || (fZNCTowerRawAOD[4] < 0)) {
    return kFALSE;
  }
	
  Double_t towZNCraw1GainEq = 0, towZNCraw2GainEq = 0, towZNCraw3GainEq = 0, towZNCraw4GainEq = 0;
  towZNCraw1GainEq = fZNCTowerRawAOD[1]*fHZDCCparameters->GetBinContent(1);
  towZNCraw2GainEq = fZNCTowerRawAOD[2]*fHZDCCparameters->GetBinContent(2);
  towZNCraw3GainEq = fZNCTowerRawAOD[3]*fHZDCCparameters->GetBinContent(3);
  towZNCraw4GainEq = fZNCTowerRawAOD[4]*fHZDCCparameters->GetBinContent(4);

  Double_t towZNAraw1GainEq = 0, towZNAraw2GainEq = 0, towZNAraw3GainEq = 0, towZNAraw4GainEq = 0;
  towZNAraw1GainEq = fZNATowerRawAOD[1]*fHZDCAparameters->GetBinContent(1);
  towZNAraw2GainEq = fZNATowerRawAOD[2]*fHZDCAparameters->GetBinContent(2);
  towZNAraw3GainEq = fZNATowerRawAOD[3]*fHZDCAparameters->GetBinContent(3);
  towZNAraw4GainEq = fZNATowerRawAOD[4]*fHZDCAparameters->GetBinContent(4);
  
  const Double_t xZDCC[4] = {-1, 1, -1, 1}; // directional vector
  const Double_t yZDCC[4] = {-1, -1, 1, 1};
  const Double_t xZDCA[4] = {1, -1, 1, -1};
  const Double_t yZDCA[4] = {-1, -1, 1, 1};
    
  Double_t towZNC[5] = {fZNCTowerRawAOD[0], towZNCraw1GainEq, towZNCraw2GainEq, towZNCraw3GainEq, towZNCraw4GainEq};
  Double_t towZNA[5] = {fZNATowerRawAOD[0], towZNAraw1GainEq, towZNAraw2GainEq, towZNAraw3GainEq, towZNAraw4GainEq};
    
  Double_t EZNC = 0, wZNC = 0, denZNC = 0, numXZNC = 0, numYZNC = 0;
  Double_t EZNA = 0, wZNA = 0, denZNA = 0, numXZNA = 0, numYZNA = 0; 

  for(Int_t i=0; i<4; i++){
   // ZNC part
   // get energy
   EZNC = towZNC[i+1];
       
   // build ZDCC centroid
   wZNC = TMath::Max(0., 4.0 + TMath::Log(towZNC[i+1]/fZNCTowerRawAOD[0]));
   numXZNC += xZDCC[i]*wZNC;
   numYZNC += yZDCC[i]*wZNC;
   denZNC += wZNC;
   
   // ZNA part
   // get energy
   EZNA = towZNA[i+1];

   // build ZDCA centroid
   wZNA = TMath::Max(0., 4.0 + TMath::Log(towZNA[i+1]/fZNATowerRawAOD[0]));
   numXZNA += xZDCA[i]*wZNA;
   numYZNA += yZDCA[i]*wZNA;
   denZNA += wZNA;
  }
  
  if (denZNC==0) {return kFALSE;}
  if (denZNA==0) {return kFALSE;}
  fdenZNC = denZNC;
  qnxZNCC = numXZNC/denZNC;
  qnyZNCC = numYZNC/denZNC;
  fdenZNA = denZNA;
  qnxZNCA = numXZNA/denZNA;
  qnyZNCA = numYZNA/denZNA;
  
  return kTRUE;
 }
 
}

void AliAnalysisTaskMixedHarmonics::ApplyZNCqVectRecenter(Float_t centrality,Double_t pVtxX,Double_t pVtxY,Double_t pVtxZ,Double_t fOrbitNumber,Double_t &qnxZNCC,Double_t &qnyZNCC,Double_t &qnxZNCA,Double_t &qnyZNCA){
	
 Double_t ZDCCAvgxPosFromVtxFit = 0;
 Double_t ZDCCAvgyPosFromVtxFit = 0;
	
 Double_t ZDCAAvgxPosFromVtxFit = 0;
 Double_t ZDCAAvgyPosFromVtxFit = 0;

 ZDCCAvgxPosFromVtxFit = fHZDCCparameters->GetBinContent(6)*centrality + fHZDCCparameters->GetBinContent(7)*pow(centrality,2) + fHZDCCparameters->GetBinContent(8)*pow(centrality,3) + fHZDCCparameters->GetBinContent(9)*pVtxX + fHZDCCparameters->GetBinContent(10)*pVtxY + fHZDCCparameters->GetBinContent(11)*pVtxZ + fHZDCCparameters->GetBinContent(12)*fOrbitNumber + fHZDCCparameters->GetBinContent(13);
 ZDCCAvgyPosFromVtxFit = fHZDCCparameters->GetBinContent(14)*centrality + fHZDCCparameters->GetBinContent(15)*pow(centrality,2) + fHZDCCparameters->GetBinContent(16)*pow(centrality,3) + fHZDCCparameters->GetBinContent(17)*pVtxX + fHZDCCparameters->GetBinContent(18)*pVtxY + fHZDCCparameters->GetBinContent(19)*pVtxZ + fHZDCCparameters->GetBinContent(20)*fOrbitNumber + fHZDCCparameters->GetBinContent(21);
	
 ZDCAAvgxPosFromVtxFit = fHZDCAparameters->GetBinContent(6)*centrality + fHZDCAparameters->GetBinContent(7)*pow(centrality,2) + fHZDCAparameters->GetBinContent(8)*pow(centrality,3) + fHZDCAparameters->GetBinContent(9)*pVtxX + fHZDCAparameters->GetBinContent(10)*pVtxY + fHZDCAparameters->GetBinContent(11)*pVtxZ + fHZDCAparameters->GetBinContent(12)*fOrbitNumber + fHZDCAparameters->GetBinContent(13);
 ZDCAAvgyPosFromVtxFit = fHZDCAparameters->GetBinContent(14)*centrality + fHZDCAparameters->GetBinContent(15)*pow(centrality,2) + fHZDCAparameters->GetBinContent(16)*pow(centrality,3) + fHZDCAparameters->GetBinContent(17)*pVtxX + fHZDCAparameters->GetBinContent(18)*pVtxY + fHZDCAparameters->GetBinContent(19)*pVtxZ + fHZDCAparameters->GetBinContent(20)*fOrbitNumber + fHZDCAparameters->GetBinContent(21);


 qnxZNCC -= ZDCCAvgxPosFromVtxFit;
 qnyZNCC -= ZDCCAvgyPosFromVtxFit;

 qnxZNCA -= ZDCAAvgxPosFromVtxFit;
 qnyZNCA -= ZDCAAvgyPosFromVtxFit;
 
 return;
}










