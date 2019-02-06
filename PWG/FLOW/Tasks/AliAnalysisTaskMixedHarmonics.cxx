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
fWeightsList(NULL)
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
fWeightsList(NULL)
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


 PostData(1,fListHistos);
  
} // end of void AliAnalysisTaskMixedHarmonics::UserCreateOutputObjects() 

//================================================================================================================

void AliAnalysisTaskMixedHarmonics::UserExec(Option_t *) 
{



 AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());

 Bool_t kPileupEvent = kFALSE;

 //cout<<" Run = "<<aodEvent->GetRunNumber();

 if(fRejectPileUp || fFillQAHistograms){
   kPileupEvent = CheckEventIsPileUp(aodEvent);
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



















