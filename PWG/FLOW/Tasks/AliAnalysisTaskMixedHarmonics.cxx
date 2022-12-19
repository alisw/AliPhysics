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
fCalculateVsZDC(kFALSE),
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
fCalculateVsZDC(kFALSE),
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
 fMH->SetCalculateVsZDC(fCalculateVsZDC); 
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


















