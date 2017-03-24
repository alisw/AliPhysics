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

///////////////////////////////////////////////
// AliAnalysisTaskVnZDC:
// Author: Rihan Haque (mhaque@cern.ch)
///////////////////////////////////////////////


#include "Riostream.h" //needed as include

class TFile;
class TList;
class AliAnalysisTaskSE;

#include "TProfile.h"  //needed as include
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "AliMultSelection.h"
#include "AliVVertex.h"

#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskVnZDC.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliLog.h"

using std::endl;
using std::cout;


ClassImp(AliAnalysisTaskVnZDC)

//________________________________________________________________________
AliAnalysisTaskVnZDC::AliAnalysisTaskVnZDC(const char *name, Bool_t usePhiWeights) :
AliAnalysisTaskSE(name),
fEvent(NULL),
fListHistos(NULL),
fMinimalBook(kFALSE),
fUsePhiWeights(usePhiWeights),
fListWeights(NULL),
fRelDiffMsub(1.0),
fApplyCorrectionForNUA(kFALSE),
fHarmonic(2),
fNormalizationType(1),
fNCentBins(9),
fTotalQvector(NULL),
fievent(0),
frunflag(0),
fHist_Vertex_XY(NULL),
fHist_Vx_vs_runnum(NULL),
fHist_Vy_vs_runnum(NULL),
fHist_Vz_vs_runnum(NULL),
fHist_tracks_vs_runnum(NULL),
fHist_Vx_ArrayFinder(NULL),
fHist_Vy_ArrayFinder(NULL),
fHist_Vz_ArrayFinder(NULL),
fHist_Vertex_Vz(NULL),
fHist_Event_count(NULL),
fHist_Cent_count1(NULL),
fHist_Cent_count2(NULL),
fHist_EventperRun(NULL),
fHist_Psi1_zdnA(NULL),
fHist_Psi1_zdnC(NULL),
fHist_Psi1_zdnA_after1(NULL),
fHist_Psi1_zdnC_after1(NULL),
fHist_Psi1_zdnA_after2(NULL),
fHist_Psi1_zdnC_after2(NULL),
fHist_vBincount(NULL),
fZDCESEList(NULL),
fZListDummy(NULL),
fFBEffiList1(NULL),
fPileUpMultSelCount(NULL),
fPileUpCount(NULL),
fMultSelection(NULL),
fAnalysisUtil(NULL),
fHist_v2xV1_ZDN_Norm_cosXX(NULL),
fHist_v2xV1_ZDN_Refm_cosXX(NULL),
fHist_v2xV1_ZDN_Cent_cosXX(NULL),
fHist_v2xV1_ZDN_Norm_cosYY(NULL),
fHist_v2xV1_ZDN_Refm_cosYY(NULL),
fHist_v2xV1_ZDN_Cent_cosYY(NULL),
fHist_v2xV1_ZDN_Norm_sinXY(NULL),
fHist_v2xV1_ZDN_Refm_sinXY(NULL),
fHist_v2xV1_ZDN_Cent_sinXY(NULL),
fHist_v2xV1_ZDN_Norm_sinYX(NULL),
fHist_v2xV1_ZDN_Refm_sinYX(NULL),
fHist_v2xV1_ZDN_Cent_sinYX(NULL),
fHist_v2xV1_ZDN_Norm_All(NULL),
fHist_v2xV1_ZDN_Refm_All(NULL),
fHist_v2xV1_ZDN_Cent_All(NULL),
fHist_ZDN_resol_Norm_All(NULL),     
fHist_ZDN_resol_Cent_All(NULL),
fHist_ZDN_resol_Refm_All(NULL),        
fHist_ZDN_resol_Norm_cos(NULL),
fHist_ZDN_resol_Refm_cos(NULL),
fHist_ZDN_resol_Cent_cos(NULL),
fHist_ZDN_resol_Norm_sin(NULL),
fHist_ZDN_resol_Refm_sin(NULL),
fHist_ZDN_resol_Cent_sin(NULL), 
fHist_ZDCn_A_XYvsRun(NULL),
fHist_ZDCn_C_XYvsRun(NULL),
fAvg_Cent_vs_Vz_Cent_woCut(NULL),
fAvg_Cent_vs_Vz_Peri_woCut(NULL),
fAvg_Cent_vs_Vx_Cent_woCut(NULL),
fAvg_Cent_vs_Vx_Peri_woCut(NULL),
fAvg_Cent_vs_Vy_Cent_woCut(NULL),
fAvg_Cent_vs_Vy_Peri_woCut(NULL),
fAvg_Cent_vs_Vz_Cent_wCuts(NULL),
fAvg_Cent_vs_Vz_Peri_wCuts(NULL),
fAvg_Cent_vs_Vx_Cent_wCuts(NULL),
fAvg_Cent_vs_Vx_Peri_wCuts(NULL),
fAvg_Cent_vs_Vy_Cent_wCuts(NULL),
fAvg_Cent_vs_Vy_Peri_wCuts(NULL),
fTPCV0M_CentDiff_vs_Vz(NULL),
fTPCV0M_CentDiff_vs_Vx(NULL),
fTPCV0M_CentDiff_vs_Vy(NULL),
fWeight_Cent(NULL),
fCent_fromDATA(NULL),
fRejectPileUpTight(kFALSE),
fRejectPileUp(kFALSE),
checkOnce(0),
vxBin(0),
vyBin(0),
vzBin(0)
{
 for(int i=0;i<90;i++)
     runNums[i] = 0;

 for(int i=0;i<90;i++){
    for(int j=0;j<10;j++){
        fHist_znCx_V0_VxVy[i][j] = NULL;
        fHist_znCy_V0_VxVy[i][j] = NULL;
        fHist_znAx_V0_VxVy[i][j] = NULL;
        fHist_znAy_V0_VxVy[i][j] = NULL;
     }
    fHist_ZDCA_En_Run[i]  = NULL;
    fHist_ZDCC_En_Run[i]  = NULL;
 }

 for(int i=0;i<2;i++)
 {
      VxCut[i] = 0;
      VyCut[i] = 0;
      VzCut[i] = 0;
 }

 for(int i=0;i<4;i++){
  for(int j=0;j<5;j++){
      fHist_X_vs_Obs_before[i][j] = NULL;
      fHist_X_vs_Obs_after1[i][j] = NULL;
      fHist_XX_vs_Obs_before[i][j] = NULL;
      fHist_XX_vs_Obs_after1[i][j] = NULL;
   }
 }

 for(int i=0;i<10;i++){
     fFB_Efficiency_Cent[i] = NULL;
 }


     DefineInput(1, AliFlowEventSimple::Class()); // Input slot #0 works with an AliFlowEventSimple
    DefineOutput(1,TList::Class());
    DefineOutput(2,TList::Class());

   fDataSet="2010";
   fAnalysisSet="recenter1";

   fTotalQvector = new TString("QaQb");         // Total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb"
   //AliDebug(2,
   cout<<"AliAnalysisTaskVnZDC::Constructor is called..!! fAnalysisSet = "<<fAnalysisSet.Data()<<endl;

}//-------------constructor-----------------

//________________________________________________
AliAnalysisTaskVnZDC::AliAnalysisTaskVnZDC() :
AliAnalysisTaskSE(),
fEvent(NULL),
fListHistos(NULL),
fMinimalBook(kFALSE),
fUsePhiWeights(kFALSE),
fListWeights(NULL),
fRelDiffMsub(1.0),
fApplyCorrectionForNUA(kFALSE),
fHarmonic(2),
fNormalizationType(1),
fNCentBins(9),
fTotalQvector(NULL),
fievent(0),
frunflag(0),
fHist_Vertex_XY(NULL),
fHist_Vx_vs_runnum(NULL),
fHist_Vy_vs_runnum(NULL),
fHist_Vz_vs_runnum(NULL),
fHist_tracks_vs_runnum(NULL),
fHist_Vx_ArrayFinder(NULL),
fHist_Vy_ArrayFinder(NULL),
fHist_Vz_ArrayFinder(NULL),
fHist_Vertex_Vz(NULL),
fHist_Event_count(NULL),
fHist_Cent_count1(NULL),
fHist_Cent_count2(NULL),
fHist_EventperRun(NULL),
fHist_Psi1_zdnA(NULL),
fHist_Psi1_zdnC(NULL),
fHist_Psi1_zdnA_after1(NULL),
fHist_Psi1_zdnC_after1(NULL),
fHist_Psi1_zdnA_after2(NULL),
fHist_Psi1_zdnC_after2(NULL),
fHist_vBincount(NULL),
fZDCESEList(NULL),
fZListDummy(NULL),
fFBEffiList1(NULL),
fPileUpMultSelCount(NULL),
fPileUpCount(NULL),
fMultSelection(NULL),
fAnalysisUtil(NULL),
fHist_v2xV1_ZDN_Norm_cosXX(NULL),
fHist_v2xV1_ZDN_Refm_cosXX(NULL),
fHist_v2xV1_ZDN_Cent_cosXX(NULL),
fHist_v2xV1_ZDN_Norm_cosYY(NULL),
fHist_v2xV1_ZDN_Refm_cosYY(NULL),
fHist_v2xV1_ZDN_Cent_cosYY(NULL),
fHist_v2xV1_ZDN_Norm_sinXY(NULL),
fHist_v2xV1_ZDN_Refm_sinXY(NULL),
fHist_v2xV1_ZDN_Cent_sinXY(NULL),
fHist_v2xV1_ZDN_Norm_sinYX(NULL),
fHist_v2xV1_ZDN_Refm_sinYX(NULL),
fHist_v2xV1_ZDN_Cent_sinYX(NULL),
fHist_ZDN_resol_Norm_cos(NULL),
fHist_ZDN_resol_Refm_cos(NULL),
fHist_ZDN_resol_Cent_cos(NULL),
fHist_ZDN_resol_Norm_sin(NULL),
fHist_ZDN_resol_Refm_sin(NULL),
fHist_ZDN_resol_Cent_sin(NULL), 
fHist_v2xV1_ZDN_Norm_All(NULL),
fHist_v2xV1_ZDN_Refm_All(NULL),
fHist_v2xV1_ZDN_Cent_All(NULL),
fHist_ZDN_resol_Norm_All(NULL),     
fHist_ZDN_resol_Cent_All(NULL),
fHist_ZDN_resol_Refm_All(NULL),        
fHist_ZDCn_A_XYvsRun(NULL),
fHist_ZDCn_C_XYvsRun(NULL),
fAvg_Cent_vs_Vz_Cent_woCut(NULL),
fAvg_Cent_vs_Vz_Peri_woCut(NULL),
fAvg_Cent_vs_Vx_Cent_woCut(NULL),
fAvg_Cent_vs_Vx_Peri_woCut(NULL),
fAvg_Cent_vs_Vy_Cent_woCut(NULL),
fAvg_Cent_vs_Vy_Peri_woCut(NULL),
fAvg_Cent_vs_Vz_Cent_wCuts(NULL),
fAvg_Cent_vs_Vz_Peri_wCuts(NULL),
fAvg_Cent_vs_Vx_Cent_wCuts(NULL),
fAvg_Cent_vs_Vx_Peri_wCuts(NULL),
fAvg_Cent_vs_Vy_Cent_wCuts(NULL),
fAvg_Cent_vs_Vy_Peri_wCuts(NULL),
fTPCV0M_CentDiff_vs_Vz(NULL),
fTPCV0M_CentDiff_vs_Vx(NULL),
fTPCV0M_CentDiff_vs_Vy(NULL),
fWeight_Cent(NULL),
fCent_fromDATA(NULL),
fRejectPileUpTight(kFALSE),
fRejectPileUp(kFALSE),
checkOnce(0),
vxBin(0),
vyBin(0),
vzBin(0)
{
 for(int i=0;i<90;i++)
     runNums[i] = 0;

 for(int i=0;i<90;i++){
    for(int j=0;j<10;j++){
        fHist_znCx_V0_VxVy[i][j] = NULL;
        fHist_znCy_V0_VxVy[i][j] = NULL;
        fHist_znAx_V0_VxVy[i][j] = NULL;
        fHist_znAy_V0_VxVy[i][j] = NULL;
     }
    fHist_ZDCA_En_Run[i]  = NULL;
    fHist_ZDCC_En_Run[i]  = NULL;
 }

 for(int i=0;i<2;i++)
 {
      VxCut[i] = 0;
      VyCut[i] = 0;
      VzCut[i] = 0;
 }

 for(int i=0;i<4;i++){
  for(int j=0;j<5;j++){
      fHist_X_vs_Obs_before[i][j] = NULL;
      fHist_X_vs_Obs_after1[i][j] = NULL;
      fHist_XX_vs_Obs_before[i][j] = NULL;
      fHist_XX_vs_Obs_after1[i][j] = NULL;
   }
 }

 for(int i=0;i<10;i++){
     fFB_Efficiency_Cent[i] = NULL;
 }


   fDataSet="2010";
   fAnalysisSet="recenter1";

	// AliDebug(2,"AliAnalysisTaskVnZDC::AliAnalysisTaskVnZDC()");
}





//________________________________________________________________________
AliAnalysisTaskVnZDC::~AliAnalysisTaskVnZDC()
{
 //delete [] fSP; //is this not deleting all fSP[i] ?
 delete fListHistos;
 delete  fHist_Vx_ArrayFinder;
 delete  fHist_Vy_ArrayFinder;
 delete  fHist_Vz_ArrayFinder;
 delete        fMultSelection;
 delete         fAnalysisUtil;
 delete                fEvent;
 delete        fCent_fromDATA;
 delete          fWeight_Cent;
 //delete  fHist_ZDCn_A_XYvsRun;
 //delete  fHist_ZDCn_C_XYvsRun;

 if(fAnalysisSet=="recenter2"){
  for(int i=0;i<90;i++){
    for(int j=0;j<10;j++){
       delete fHist_znCx_V0_VxVy[i][j];
       delete fHist_znCy_V0_VxVy[i][j];
       delete fHist_znAx_V0_VxVy[i][j];
       delete fHist_znAy_V0_VxVy[i][j];
      }
    delete fHist_ZDCA_En_Run[i];
    delete fHist_ZDCC_En_Run[i];
  }
 }


 for(int i=0;i<10;i++){
    delete fFB_Efficiency_Cent[i];
  }

 delete           fListHistos;
 delete          fFBEffiList1;
 delete          fListWeights;
 delete           fZDCESEList;
 delete           fZListDummy;

 printf("\n\n ~AliAnalysisTaskVnZDC::Destructor is called..!!\n\n");

}

//________________________________________________________________________
void AliAnalysisTaskVnZDC::UserCreateOutputObjects()
{
  //printf("\n\n... AliAnalysisTaskVnZDC::CreateOutputObjects() beginning 0 ...\n");
  //AliDebug(2,"AliAnalysisTaskVnZDC::CreateOutputObjects() is called....");

 Int_t runArray_2010[89] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137243, 137236, 137235, 137232, 137231, 137162, 137161};

 Int_t runArray_2011[68] = {167915, 168115, 168460, 169035, 169238, 169859, 170228, 167920, 168310, 168464, 169091, 169411, 169923, 170230, 167985, 168311, 168467, 169094, 169415, 170027, 170268, 167987, 168322, 168511, 169138, 169417, 170081, 170269, 167988, 168325, 168512, 169144, 169835, 170155, 170270, 168069, 168341, 168514, 169145, 169837, 170159, 170306, 168076, 168342, 168777, 169148, 169838, 170163, 170308, 168105, 168361, 168826, 169156, 169846, 170193, 170309, 168107, 168362, 168988, 169160, 169855, 170203, 168108, 168458, 168992, 169167, 169858, 170204};

 Int_t runArray_2015[90] = {246994, 246991, 246989, 246984, 246982, 246980, 246948, 246945, 246928, 246871, 246870, 246867, 246865, 246864, 246859, 246858, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246676, 246675, 246540, 246495, 246493, 246488, 246487, 246434, 246431, 246428, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245700, 245692, 245683};

 //Int_t runArray_2015[77] = {246994, 246991, 246989, 246984, 246982, 246948, 246945, 246928, 246851, 246847, 246846, 246845, 246844, 246810, 246809, 246808, 246807, 246805, 246804, 246766, 246765, 246763, 246760, 246759, 246758, 246757, 246751, 246750, 246495, 246493, 246488, 246487, 246434, 246431, 246424, 246276, 246275, 246272, 246271, 246225, 246222, 246217, 246185, 246182, 246181, 246180, 246178, 246153, 246152, 246151, 246148, 246115, 246113, 246089, 246087, 246053, 246052, 246049, 246048, 246042, 246037, 246036, 246012, 246003, 246001, 245963, 245954, 245952, 245949, 245923, 245833, 245831, 245829, 245705, 245702, 245692, 245683};


 if(fDataSet=="2010"){
    frunflag = 89;
    for(int i=0;i<frunflag;i++)
       runNums[i] = runArray_2010[i];
 }
 if(fDataSet=="2011"){
    frunflag = 10;                 //<--------- 2011 runnumbers have to be checked and put...
    for(int i=0;i<frunflag;i++)
       runNums[i] = runArray_2011[i];
 }
 if(fDataSet=="2015"){
    frunflag = 90;
    for(int i=0;i<frunflag;i++)
       runNums[i] = runArray_2015[i];
 }


  Int_t totalQ = 0;
  if(fTotalQvector->Contains("Qa")) totalQ += 1;
  if(fTotalQvector->Contains("Qb")) totalQ += 2;

  fListHistos = new TList();
  fListHistos->SetOwner(kTRUE);


  fHist_Event_count    = new TH1F("fHist_Event_count"," ",20,0,20);
  fListHistos->Add(fHist_Event_count);
  fHist_vBincount      = new TH1F("fHist_vBincount"," ",500,0,500);
  fListHistos->Add(fHist_vBincount);
  fHist_EventperRun    = new TH1F("fHist_EventperRun","", frunflag,0,frunflag);
  fListHistos->Add(fHist_EventperRun);
  fHist_Psi1_zdnA      = new TH1F("fHist_Psi1_zdnA","",  200, 0,2.*TMath::Pi());
  fListHistos->Add(fHist_Psi1_zdnA);
  fHist_Psi1_zdnC      = new TH1F("fHist_Psi1_zdnC","",  200, 0,2.*TMath::Pi());
  fListHistos->Add(fHist_Psi1_zdnC);
  fHist_Psi1_zdnA_after1      = new TH1F("fHist_Psi1_zdnA_after1","",  200, 0,2.*TMath::Pi());
  fListHistos->Add(fHist_Psi1_zdnA_after1);
  fHist_Psi1_zdnC_after1      = new TH1F("fHist_Psi1_zdnC_after1","",  200, 0,2.*TMath::Pi());
  fListHistos->Add(fHist_Psi1_zdnC_after1);
  fHist_Psi1_zdnA_after2      = new TH1F("fHist_Psi1_zdnA_after2","",  200, 0,2.*TMath::Pi());
  fListHistos->Add(fHist_Psi1_zdnA_after2);
  fHist_Psi1_zdnC_after2      = new TH1F("fHist_Psi1_zdnC_after2","",  200, 0,2.*TMath::Pi());
  fListHistos->Add(fHist_Psi1_zdnC_after2);


  fHist_Vertex_Vz    = new TH1F("fHist_Vertex_Vz_dist","Vertex Z Distribution",200,  -15,    15);
  fListHistos->Add(fHist_Vertex_Vz);
  fHist_Vertex_XY    = new TH2F("fHist_Vertex_XY_dist","Vertex XY Distributn.",100,-0.2,0.2,100,0.0,0.4);
  fListHistos->Add(fHist_Vertex_XY);

  // Vx,Vy,Vz vs runnumber:
  fHist_Vx_vs_runnum = new TProfile("fHist_Vx_vs_runnum","<Vx>_vs_runnum",frunflag,0,frunflag);
  fListHistos->Add(fHist_Vx_vs_runnum);
  fHist_Vy_vs_runnum = new TProfile("fHist_Vy_vs_runnum","<Vy>_vs_runnum",frunflag,0,frunflag);
  fListHistos->Add(fHist_Vy_vs_runnum);
  fHist_Vz_vs_runnum = new TProfile("fHist_Vz_vs_runnum","<Vz>_vs_runnum",frunflag,0,frunflag);
  fListHistos->Add(fHist_Vz_vs_runnum);
  fHist_tracks_vs_runnum = new TProfile("fHist_tracks_vs_runnum","<nTracks>_vs_runnum",frunflag,0,frunflag);
  fListHistos->Add(fHist_tracks_vs_runnum);

  fHist_ZDCn_A_XYvsRun = new TH3F("fHist_ZDCn_A_XYvsRun","ZDC XY vs run",80,-2.0,2.0,80,-2.0,2.0,89,0,89);
  fListHistos->Add(fHist_ZDCn_A_XYvsRun);
  fHist_ZDCn_C_XYvsRun = new TH3F("fHist_ZDCn_C_XYvsRun","ZDC XY vs run",80,-2.0,2.0,80,-2.0,2.0,89,0,89);
  fListHistos->Add(fHist_ZDCn_C_XYvsRun);


  Double_t centRange[12]    = {0,5,10,20,30,40,50,60,70,80,90,100};

  fPileUpCount = new TH1F("fPileUpCount", "fPileUpCount", 9, 0., 9.);
  fPileUpCount->GetXaxis()->SetBinLabel(1,"plpMV");
  fPileUpCount->GetXaxis()->SetBinLabel(2,"fromSPD");
  fPileUpCount->GetXaxis()->SetBinLabel(3,"RefMultiplicityComb08");
  fPileUpCount->GetXaxis()->SetBinLabel(4,"IncompleteDAQ");
  fPileUpCount->GetXaxis()->SetBinLabel(5,"abs(V0M-CL1)>7.5");
  fPileUpCount->GetXaxis()->SetBinLabel(6,"missingVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(7,"inconsistentVtx");
  fPileUpCount->GetXaxis()->SetBinLabel(8,"multESDTPCDif");
  fPileUpCount->GetXaxis()->SetBinLabel(9,"extraPileUpMultSel");
  fListHistos->Add(fPileUpCount);

  fPileUpMultSelCount = new TH1F("fPileUpMultSelCount", "fPileUpMultSelCount", 8, 0., 8.);
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(1,"IsNotPileup");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(2,"IsNotPileupMV");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(3,"IsNotPileupInMultBins");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(4,"InconsistentVertices");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(5,"TrackletVsCluster");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(6,"AsymmetricInVZERO");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(7,"IncompleteDAQ");
  fPileUpMultSelCount->GetXaxis()->SetBinLabel(8,"GoodVertex2016");
  fListHistos->Add(fPileUpMultSelCount);




 //-------- define and add the recenter histograms ----------------
if(fDataSet=="2010") {
   vxBin = 8;
   vyBin = 8;
 //vxBin = 20;
 //vyBin = 20;
   vzBin = 10;   // j<10

   VxCut[0] =  -0.035;
   VxCut[1] =   0.020;
   VyCut[0] =   0.144;
   VyCut[1] =   0.214;
   VzCut[0] =   -10.0;
   VzCut[1] =    10.0;
 }

 if(fDataSet=="2015"){
   vxBin = 8;
   vyBin = 8;
 //vxBin = 20;
 //vyBin = 20;

   vzBin = 10;   // j<10

   VxCut[0] =   0.060;
   VxCut[1] =   0.086;
   VyCut[0] =   0.321;
   VyCut[1] =   0.345;
   VzCut[0] =   -10.0;
   VzCut[1] =    10.0;
  }
 if(fDataSet=="2011"){
    AliDebug(2,"\n\n !!** WARNING ***!!  \n cuts not defined for DATASET: %s\n ...EXIT...\n\n)");
    exit(0);
 }


  const int NbinVt =  (vyBin-1)*vxBin + vxBin;

  fHist_Vx_ArrayFinder = new TH1F("fHist_Vx_ArrayFinder","",vxBin,VxCut[0],VxCut[1]);
  fHist_Vy_ArrayFinder = new TH1F("fHist_Vy_ArrayFinder","",vyBin,VyCut[0],VyCut[1]);
  fHist_Vz_ArrayFinder = new TH1F("fHist_Vz_ArrayFinder","",vzBin,VzCut[0],VzCut[1]);





 if(fAnalysisSet=="recenter2") {

  TString     sNameQn[4] = {"Xa","Xc","Ya","Yc"};
  TString    sNameQn2[4] = {"XaXc","YaYc","XcYa","YcXa"};
  TString     sNameVs[5] = {"Cent","Mult","Vx","Vy","Vz"};       // ********** to add vs cent ****** !!!!!!!!!!!!
  Int_t     nBinNumVs[5] = {100, 500,  40, 40,400};  // number of bins for "Mult","Vx","Vy","Vz"
  Double_t  lBinLowVs[5] = {  0,   0,-0.1,0.1,-10}; // bin low value: "Mult","Vx","Vy","Vz"
  Double_t lBinHighVs[5] = {100,2000, 0.1,0.3, 10};// bin high value: "Mult","Vx","Vy","Vz"

  char name[100];


    for(int i=0;i<4;i++) {
      for(int j=0;j<5;j++) {
       //store: X,Y position for recenter:
         sprintf(name,"fHist_%s_vs_%s_before",static_cast<const char*>(sNameQn[i]),static_cast<const char*>(sNameVs[j]));
         fHist_X_vs_Obs_before[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
         fListHistos->Add(fHist_X_vs_Obs_before[i][j]);

         sprintf(name,"fHist_%s_vs_%s_after1",static_cast<const char*>(sNameQn[i]),static_cast<const char*>(sNameVs[j]));
         fHist_X_vs_Obs_after1[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
         fListHistos->Add(fHist_X_vs_Obs_after1[i][j]);

         sprintf(name,"fHist_%s_vs_%s_before",static_cast<const char*>(sNameQn2[i]),static_cast<const char*>(sNameVs[j]));
         fHist_XX_vs_Obs_before[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
         fListHistos->Add(fHist_XX_vs_Obs_before[i][j]);

         sprintf(name,"fHist_%s_vs_%s_after1",static_cast<const char*>(sNameQn2[i]),static_cast<const char*>(sNameVs[j]));
         fHist_XX_vs_Obs_after1[i][j] = new TProfile(name,"", nBinNumVs[j], lBinLowVs[j], lBinHighVs[j],"");
         fListHistos->Add(fHist_XX_vs_Obs_after1[i][j]);
       }
    }

  Double_t centRange[12] = {0,5,10,20,30,40,50,60,70,80,90,100};

  fHist_v2xV1_ZDN_Norm_cosXX  = new TProfile("fHist_v2xV1_ZDN_Norm_cosXX","v2 X V1^2 (ZDC) cos(2p)XX",11,centRange);
  fHist_v2xV1_ZDN_Norm_cosXX->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Norm_cosXX);
  fHist_v2xV1_ZDN_Refm_cosXX  = new TProfile("fHist_v2xV1_ZDN_Refm_cosXX","v2 X V1^2 (ZDC) cos(2p)XX", 40,0, 4000);
  fHist_v2xV1_ZDN_Refm_cosXX->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Refm_cosXX);
  fHist_v2xV1_ZDN_Cent_cosXX  = new TProfile("fHist_v2xV1_ZDN_Cent_cosXX","v2 X V1^2 (ZDC) cos(2p)XX",100,0,  100);
  fHist_v2xV1_ZDN_Cent_cosXX->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Cent_cosXX);

  fHist_v2xV1_ZDN_Norm_cosYY  = new TProfile("fHist_v2xV1_ZDN_Norm_cosYY","v2 X V1^2 (ZDC) cos(2p)YY",11,centRange);
  fHist_v2xV1_ZDN_Norm_cosYY->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Norm_cosYY);
  fHist_v2xV1_ZDN_Refm_cosYY  = new TProfile("fHist_v2xV1_ZDN_Refm_cosYY","v2 X V1^2 (ZDC) cos(2p)YY", 40,0, 4000);
  fHist_v2xV1_ZDN_Refm_cosYY->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Refm_cosYY);
  fHist_v2xV1_ZDN_Cent_cosYY  = new TProfile("fHist_v2xV1_ZDN_Cent_cosYY","v2 X V1^2 (ZDC) cos(2p)YY",100,0,  100);
  fHist_v2xV1_ZDN_Cent_cosYY->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Cent_cosYY);

  fHist_v2xV1_ZDN_Norm_sinXY  = new TProfile("fHist_v2xV1_ZDN_Norm_sinXY","v2  X V1^2 (ZDC) sin(2p)XY",11,centRange);
  fHist_v2xV1_ZDN_Norm_sinXY->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Norm_sinXY);
  fHist_v2xV1_ZDN_Refm_sinXY  = new TProfile("fHist_v2xV1_ZDN_Refm_sinXY"," v2 X V1^2 (ZDC) sin(2p)XY", 40,0, 4000);
  fHist_v2xV1_ZDN_Refm_sinXY->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Refm_sinXY);
  fHist_v2xV1_ZDN_Cent_sinXY  = new TProfile("fHist_v2xV1_ZDN_Cent_sinXY"," v2 X V1^2 (ZDC) sin(2p)XY",100,0,  100);
  fHist_v2xV1_ZDN_Cent_sinXY->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Cent_sinXY);  

  fHist_v2xV1_ZDN_Norm_sinYX  = new TProfile("fHist_v2xV1_ZDN_Norm_sinYX","v2  X V1^2 (ZDC) sin(2p)YX",11,centRange);
  fHist_v2xV1_ZDN_Norm_sinYX->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Norm_sinYX);
  fHist_v2xV1_ZDN_Refm_sinYX  = new TProfile("fHist_v2xV1_ZDN_Refm_sinYX"," v2 X V1^2 (ZDC) sin(2p)YX", 40,0, 4000);
  fHist_v2xV1_ZDN_Refm_sinYX->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Refm_sinYX);
  fHist_v2xV1_ZDN_Cent_sinYX  = new TProfile("fHist_v2xV1_ZDN_Cent_sinYX"," v2 X V1^2 (ZDC) sin(2p)YX",100,0,  100);
  fHist_v2xV1_ZDN_Cent_sinYX->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Cent_sinYX);  

  fHist_ZDN_resol_Norm_cos  = new TProfile("fHist_ZDN_resol_Norm_XX","Resol. V1^2(ZDC) XX",11,centRange);
  fHist_ZDN_resol_Norm_cos->Sumw2();
  fListHistos->Add(fHist_ZDN_resol_Norm_cos);
  fHist_ZDN_resol_Refm_cos  = new TProfile("fHist_ZDN_resol_Refm_XX","Resol. V1^2(ZDC) XX", 40,0, 4000);
  fHist_ZDN_resol_Refm_cos->Sumw2();
  fListHistos->Add(fHist_ZDN_resol_Refm_cos);
  fHist_ZDN_resol_Cent_cos  = new TProfile("fHist_ZDN_resol_Cent_XX","Resol. V1^2(ZDC) XX",100,0,  100);
  fHist_ZDN_resol_Cent_cos->Sumw2();
  fListHistos->Add(fHist_ZDN_resol_Cent_cos); 

  fHist_ZDN_resol_Norm_sin  = new TProfile("fHist_ZDN_resol_Norm_YY","Resol. V1^2(ZDC) YY",11,centRange);
  fHist_ZDN_resol_Norm_sin->Sumw2();
  fListHistos->Add(fHist_ZDN_resol_Norm_sin);
  fHist_ZDN_resol_Refm_sin  = new TProfile("fHist_ZDN_resol_Refm_YY","Resol. V1^2(ZDC) YY",40,0, 4000);
  fHist_ZDN_resol_Refm_sin->Sumw2();
  fListHistos->Add(fHist_ZDN_resol_Refm_sin);
  fHist_ZDN_resol_Cent_sin  = new TProfile("fHist_ZDN_resol_Cent_YY","Resol. V1^2(ZDC) YY",100,0,  100);
  fHist_ZDN_resol_Cent_sin->Sumw2();
  fListHistos->Add(fHist_ZDN_resol_Cent_sin);  

  //all terms:

  fHist_v2xV1_ZDN_Norm_All  = new TProfile("fHist_v2xV1_ZDN_Norm_All","v2 X V1^2 (ZDC) all terms",11,centRange);
  fHist_v2xV1_ZDN_Norm_All->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Norm_All);
  fHist_v2xV1_ZDN_Refm_All  = new TProfile("fHist_v2xV1_ZDN_Refm_All","v2 X V1^2 (ZDC) all terms", 40,0, 4000);
  fHist_v2xV1_ZDN_Refm_All->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Refm_All);
  fHist_v2xV1_ZDN_Cent_All  = new TProfile("fHist_v2xV1_ZDN_Cent_All","v2 X V1^2 (ZDC) all terms",100,0,  100);
  fHist_v2xV1_ZDN_Cent_All->Sumw2();
  fListHistos->Add(fHist_v2xV1_ZDN_Cent_All);

  fHist_ZDN_resol_Norm_All  = new TProfile("fHist_ZDN_resol_Norm_All","Resol. V1^2(ZDC) All",11,centRange);
  fHist_ZDN_resol_Norm_All->Sumw2();
  fListHistos->Add(fHist_ZDN_resol_Norm_All);
  fHist_ZDN_resol_Refm_All  = new TProfile("fHist_ZDN_resol_Refm_All","Resol. V1^2(ZDC) All",40,0, 4000);
  fHist_ZDN_resol_Refm_All->Sumw2();
  fListHistos->Add(fHist_ZDN_resol_Refm_All);
  fHist_ZDN_resol_Cent_All  = new TProfile("fHist_ZDN_resol_Cent_All","Resol. V1^2(ZDC) All",100,0,  100);
  fHist_ZDN_resol_Cent_All->Sumw2();
  fListHistos->Add(fHist_ZDN_resol_Cent_All);  


  //------------ calculate centrality weight: ------------
   if(fZDCESEList) {
      fCent_fromDATA = (TH1F *)  fZDCESEList->FindObject("fHist_Cent_afterr_ZDCcut");
    }
    else{
      AliDebug(2,"\n\n ******** Running without Centrality weight !!******** \n\n");
      fCent_fromDATA = new TH1F("fCent_fromDATA","Centrality distribution",100,0,100);
       for(int i=1;i<=100;i++){
          fCent_fromDATA->SetBinContent(i,100);
       }
    }

   Double_t maxCount = -1;
   Int_t    maxCent  = -1;
   Double_t Content,weight,error;

   for(int i=1;i<=fCent_fromDATA->GetNbinsX();i++){
     Content = fCent_fromDATA->GetBinContent(i);
      if(maxCount < Content){
         maxCount = Content;
         maxCent = i;
      }
   }

   fWeight_Cent = new TH1F("fWeight_Cent","Weight for centrality",100,0,100);

   for(int i=1;i<=fCent_fromDATA->GetNbinsX();i++)
   {
     Content = fCent_fromDATA->GetBinContent(i);
     error   = fCent_fromDATA->GetBinError(i);
     weight  = maxCount/Content;
     if(weight>1e9) continue;
     fWeight_Cent->SetBinContent(i,weight);
     fWeight_Cent->SetBinError(i,weight/sqrt(Content));
     //cout<<"cent "<<i-1<<"-"<<i<<" \t wgt = "<<fWeight_Cent->GetBinContent(i)<<" error = "<<fWeight_Cent->GetBinError(i)<<endl;
   }
 
   fListHistos->Add(fWeight_Cent);
   //------------------------------------------------------



 //---------Filter bit efficiency----------
   if(!fFBEffiList1) {
     printf("\n\n !!**  Warning ***!!  \n TList FilterBit efficiency not found !!\n\n");
   }

 TString name2;

 if(fFBEffiList1) {
   for(int i=0;i<10;i++) {
      fFB_Efficiency_Cent[i] = (TH1D *) fFBEffiList1->FindObject(Form("eff_unbiased_%d",i));
      name2 = fFB_Efficiency_Cent[i]->GetName();
    //cout<<i<<" name = "<<name2<<endl;
   }
 }
 else{ // if MC efficiency not found then use weight = 1. Define histograms so that code doesn't crash
    for(int i=0;i<10;i++){
       fFB_Efficiency_Cent[i] = new TH1D(Form("eff_unbiased_%d",i),"",100,0,20); 
        for(int j=1;j<=100;j++){
           fFB_Efficiency_Cent[i]->SetBinContent(j,1.000);
     }
   } 
 }
 //---------Filter bit efficiency----------


 if(!fZDCESEList) {
    printf("\n\n !!** ERROR ***!!  \n \n TList fZDCESEList not found !!\n .......EXIT...... \n\n)");
    exit(1);
  }
  if(fZDCESEList) {
   for(int i=0;i<frunflag;i++) { 
     for(int j=0;j<10;j++){
         fHist_znCx_V0_VxVy[i][j] = (TProfile2D *) fZDCESEList->FindObject(Form("fHist_znCx_V0_Run%d_Vz%d",runNums[i],j+1));
         fHist_znCy_V0_VxVy[i][j] = (TProfile2D *) fZDCESEList->FindObject(Form("fHist_znCy_V0_Run%d_Vz%d",runNums[i],j+1));
         fHist_znAx_V0_VxVy[i][j] = (TProfile2D *) fZDCESEList->FindObject(Form("fHist_znAx_V0_Run%d_Vz%d",runNums[i],j+1));
         fHist_znAy_V0_VxVy[i][j] = (TProfile2D *) fZDCESEList->FindObject(Form("fHist_znAy_V0_Run%d_Vz%d",runNums[i],j+1));
          if(!fHist_znCx_V0_VxVy[i][j] || !fHist_znCy_V0_VxVy[i][j] || !fHist_znAx_V0_VxVy[i][j] || !fHist_znAy_V0_VxVy[i][j]) {
  	     printf("\n\n !!** WARNING *** One/more Recenter1 histograms NOT FOUND\n ...........!! \n\n");
	   //exit(0);
      }
     }
    }
   }
  else{ // if recenter file not found then define empty profile histograms so that code doesn't crash
    //AliDebug(2,"\n\n ******** Running without Recentering 1 !!******** \n\n");
    printf("\n\n ******** Running without Recentering 1 !!******** \n\n");
    for(int i=0;i<frunflag;i++){
     for(int j=0;j<10;j++){
         fHist_znCx_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znCx_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90,"");
         fHist_znCy_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znCy_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90,"");
         fHist_znAx_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znAx_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90,"");
         fHist_znAy_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znAy_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90,"");
     }
    }
   }
 }//--------------"recenter2"------------------


    fHist_Cent_count1         = new TH1F("fHist_Cent_before_ZDCcut"," ",100,0,100);
    fHist_Cent_count2         = new TH1F("fHist_Cent_afterr_ZDCcut"," ",100,0,100);


 if(fAnalysisSet=="recenter1") {

     fZDCESEList = new TList();
     fZDCESEList->SetOwner(kTRUE);
  
     fAvg_Cent_vs_Vz_Cent_woCut  = new TProfile("fAvg_Cent_vs_Vz_Cent_woCut","<cent> vs Vz",200,-10,10);
     fListHistos->Add(fAvg_Cent_vs_Vz_Cent_woCut);
     fAvg_Cent_vs_Vz_Peri_woCut  = new TProfile("fAvg_Cent_vs_Vz_Peri_woCut","<peri> vs Vz",200,-10,10);
     fListHistos->Add(fAvg_Cent_vs_Vz_Peri_woCut); 
     fAvg_Cent_vs_Vx_Cent_woCut  = new TProfile("fAvg_Cent_vs_Vx_Cent_woCut","<cent> vs Vx",100,VxCut[0], VxCut[1]);
     fListHistos->Add(fAvg_Cent_vs_Vx_Cent_woCut);
     fAvg_Cent_vs_Vx_Peri_woCut  = new TProfile("fAvg_Cent_vs_Vx_Peri_woCut","<peri> vs Vx",100,VxCut[0], VxCut[1]);
     fListHistos->Add(fAvg_Cent_vs_Vx_Peri_woCut); 
     fAvg_Cent_vs_Vy_Cent_woCut  = new TProfile("fAvg_Cent_vs_Vy_Cent_woCut","<cent> vs Vy",100,VyCut[0], VyCut[1]);
     fListHistos->Add(fAvg_Cent_vs_Vy_Cent_woCut);
     fAvg_Cent_vs_Vy_Peri_woCut  = new TProfile("fAvg_Cent_vs_Vy_Peri_woCut","<peri> vs Vy",100,VyCut[0], VyCut[1]);
     fListHistos->Add(fAvg_Cent_vs_Vy_Peri_woCut); 

     fAvg_Cent_vs_Vz_Cent_wCuts  = new TProfile("fAvg_Cent_vs_Vz_Cent_wCuts","<cent> vs Vz",200,-10,10);
     fListHistos->Add(fAvg_Cent_vs_Vz_Cent_wCuts);
     fAvg_Cent_vs_Vz_Peri_wCuts  = new TProfile("fAvg_Cent_vs_Vz_Peri_wCuts","<peri> vs Vz",200,-10,10);
     fListHistos->Add(fAvg_Cent_vs_Vz_Peri_wCuts); 
     fAvg_Cent_vs_Vx_Cent_wCuts  = new TProfile("fAvg_Cent_vs_Vx_Cent_wCuts","<cent> vs Vx",100,VxCut[0], VxCut[1]);
     fListHistos->Add(fAvg_Cent_vs_Vx_Cent_wCuts);
     fAvg_Cent_vs_Vx_Peri_wCuts  = new TProfile("fAvg_Cent_vs_Vx_Peri_wCuts","<peri> vs Vx",100,VxCut[0], VxCut[1]);
     fListHistos->Add(fAvg_Cent_vs_Vx_Peri_wCuts); 
     fAvg_Cent_vs_Vy_Cent_wCuts  = new TProfile("fAvg_Cent_vs_Vy_Cent_wCuts","<cent> vs Vy",100,VyCut[0], VyCut[1]);
     fListHistos->Add(fAvg_Cent_vs_Vy_Cent_wCuts);
     fAvg_Cent_vs_Vy_Peri_wCuts  = new TProfile("fAvg_Cent_vs_Vy_Peri_wCuts","<peri> vs Vy",100,VyCut[0], VyCut[1]);
     fListHistos->Add(fAvg_Cent_vs_Vy_Peri_wCuts); 


        fTPCV0M_CentDiff_vs_Vz   = new TProfile("fTPCV0M_CentDiff_vs_Vz","<TPC-V0M> vs Vz",200,-10,10);
        fListHistos->Add(fTPCV0M_CentDiff_vs_Vz);
        fTPCV0M_CentDiff_vs_Vx   = new TProfile("fTPCV0M_CentDiff_vs_Vx","<TPC-V0M> vs Vx",100,VxCut[0], VxCut[1]);
        fListHistos->Add(fTPCV0M_CentDiff_vs_Vx);
        fTPCV0M_CentDiff_vs_Vy   = new TProfile("fTPCV0M_CentDiff_vs_Vy","<TPC-V0M> vs Vy",100,VyCut[0], VyCut[1]);
        fListHistos->Add(fTPCV0M_CentDiff_vs_Vy);

   //std::cout<<"\n\n Rihan: AliAnalysisTaskVnZDC::CreateOutputObjects() is called....\n\n"<<std::endl;

    for(int i=0;i<frunflag;i++){
    //store: ZDC energy for gain calibration:
        fHist_ZDCA_En_Run[i]  = new TProfile2D(Form("fHist_ZDCA_En_Run%d",runNums[i]),"",100,0,100,5,0,5,"");
        fHist_ZDCC_En_Run[i]  = new TProfile2D(Form("fHist_ZDCC_En_Run%d",runNums[i]),"",100,0,100,5,0,5,"");
    //store: <X>,<Y> run by run for recenter://90 centrality bins:
        for(int j=0;j<10;j++){
            fHist_znCx_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znCx_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90,"");
            fHist_znCy_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znCy_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90,"");
            fHist_znAx_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znAx_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90,"");
            fHist_znAy_V0_VxVy[i][j] = new TProfile2D(Form("fHist_znAy_V0_Run%d_Vz%d",runNums[i],j+1),"",NbinVt,0,NbinVt,90,0,90,"");
        }
     }

    for(int i=0;i<frunflag;i++){
        fListHistos->Add(fHist_ZDCA_En_Run[i]);
        fListHistos->Add(fHist_ZDCC_En_Run[i]);
        for(int j=0;j<10;j++) {
            fZDCESEList->Add(fHist_znCx_V0_VxVy[i][j]);
            fZDCESEList->Add(fHist_znCy_V0_VxVy[i][j]);
            fZDCESEList->Add(fHist_znAx_V0_VxVy[i][j]);
            fZDCESEList->Add(fHist_znAy_V0_VxVy[i][j]);
        }
    }

    fZDCESEList->Add(fHist_Cent_count1);
    fZDCESEList->Add(fHist_Cent_count2);

 }//---------- recenter1 ------------


  fAnalysisUtil = new AliAnalysisUtils();
  fAnalysisUtil->SetUseMVPlpSelection(kTRUE);
  fAnalysisUtil->SetUseOutOfBunchPileUp(kTRUE);


  PostData(1,fListHistos); //posting in slot #1,

  if(fAnalysisSet=="recenter1") {
     PostData(2,fZDCESEList); 
  }
  else{
     fZListDummy = new TList();
     fZListDummy->SetOwner(kTRUE);
     fZListDummy->Add(fHist_Cent_count1);
     fZListDummy->Add(fHist_Cent_count2);
     PostData(2,fZListDummy); 
  }

  //AliDebug(2,
  printf("\n\n::UserCreateOutPutObject(). NbinVt= %d, frunflag= %d, dataset: %s, Analysis= %s\n\n",NbinVt,frunflag,fDataSet.Data(),fAnalysisSet.Data());

}

//________________________________________________________________________
void AliAnalysisTaskVnZDC::UserExec(Option_t *)
{

  int stepCount = 0;
  //printf("\n ... ::UserExec() is being called. 1 Step %d...  \n",stepCount);

  fHist_Event_count->Fill(stepCount);
  stepCount++;

  /*if(!checkOnce){
  checkOnce++;
  }*/

  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
  fEvent           = dynamic_cast<AliFlowEventSimple*>(GetInputData(1));


  if(!aod){
     printf("\n ... ::UserExec no aod found.....  \n");
     return;
  }

  fHist_Event_count->Fill(stepCount);
  stepCount++;


  AliAODVertex *pVertex = aod->GetPrimaryVertex();
  Double_t Vxyz[3]    =         {0,0,0};
  Vxyz[0]             = pVertex->GetX();
  Vxyz[1]             = pVertex->GetY();
  Vxyz[2]             = pVertex->GetZ();


  //------- Apply Necessary event cuts ---------
  if(Vxyz[2] >= VzCut[1]  || Vxyz[2] <= VzCut[0])  return;
  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(Vxyz[0] >= VxCut[1]  || Vxyz[0] <= VxCut[0])  return;
  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(Vxyz[1] >= VyCut[1]  || Vxyz[1] <= VyCut[0])  return;
  fHist_Event_count->Fill(stepCount);
  stepCount++;



  Double_t EvtCent = fEvent->GetCentrality();
  fHist_Cent_count1->         Fill(EvtCent);

  //---------- get runindex: --------------
  Int_t runNumber = aod->GetRunNumber();
  Int_t runindex = -111;

 for(int i=0;i<frunflag;i++){
   if(runNumber==runNums[i])
    {
      runindex = i;
      break;
    }
  }
 if(runindex<0) {
    AliDebug(2,"\n ::UserExec Runnumber is not listed .....  \n");
    exit(1);
  }
 //-----------------------------------------


 //--------- starting pileup rejection work: --------
   Double_t centrV0M=300; Double_t centrCL1=300;
   Double_t centrCL0=300; Double_t centrTRK=300;

 if(fDataSet=="2010"||fDataSet=="2011"){
    centrV0M = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
    centrCL1 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");
    centrCL0 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL0");
    centrTRK = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("TRK");
  }// 2010/2011

 else{
     fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      if(!fMultSelection) {
        //If you get this warning, please check if AliMultSelectionTask actually ran (before your task)
        AliDebug(2,Form("\n **WARNING** ::UserExec() AliMultSelection object not found. Step# %d\n",stepCount));
        exit(1);
    }
      centrV0M = fMultSelection->GetMultiplicityPercentile("V0M");
      centrCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
      centrCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
      centrTRK = fMultSelection->GetMultiplicityPercentile("TRK");
  }// 2015


     Bool_t BisPileup=kFALSE;

  if(fRejectPileUp && InputEvent()) {
    //if(!fCutsEvent->IsSelected(InputEvent(),MCEvent())) return;
    //if(fDataSet!=k2015 && fDataSet!=k2015v6) { //jacopo
    if(fDataSet!="2015") {
          if(plpMV(aod)) {
            fPileUpCount->Fill(0.5);
            BisPileup=kTRUE;
          }
          Int_t isPileup = aod->IsPileupFromSPD(3);
          if(isPileup != 0) {
            fPileUpCount->Fill(1.5);
            BisPileup=kTRUE;
          }
          if(((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
            fPileUpCount->Fill(2.5);
            BisPileup=kTRUE;
          }
          if(aod->IsIncompleteDAQ())  {
            fPileUpCount->Fill(3.5);
            BisPileup=kTRUE;
          }

    //check vertex consistency
     const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
     const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();

     if(vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
            fPileUpCount->Fill(5.5);
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
            fPileUpCount->Fill(6.5);
            BisPileup=kTRUE;
          }
          if(fAnalysisUtil->IsPileUpEvent(InputEvent())) {
            fPileUpCount->Fill(7.5);
            BisPileup=kTRUE;
          }
        }

        else {
          // pileup from AliMultSelection
          if(!fMultSelection->GetThisEventIsNotPileup())
             fPileUpMultSelCount->Fill(0.5);
          if(!fMultSelection->GetThisEventIsNotPileupMV())
             fPileUpMultSelCount->Fill(1.5);
          if(!fMultSelection->GetThisEventIsNotPileupInMultBins())
             fPileUpMultSelCount->Fill(2.5);
          if(!fMultSelection->GetThisEventHasNoInconsistentVertices())
             fPileUpMultSelCount->Fill(3.5);
          if(!fMultSelection->GetThisEventPassesTrackletVsCluster())
             fPileUpMultSelCount->Fill(4.5);
          if(!fMultSelection->GetThisEventIsNotAsymmetricInVZERO())
             fPileUpMultSelCount->Fill(5.5);
          if(!fMultSelection->GetThisEventIsNotIncompleteDAQ())
             fPileUpMultSelCount->Fill(6.5);
          if(!fMultSelection->GetThisEventHasGoodVertex2016())
             fPileUpMultSelCount->Fill(7.5);

             BisPileup=kFALSE;

      // pile-up a la Dobrin for LHC15o
          if(plpMV(aod)) {
            fPileUpCount->Fill(0.5);
            BisPileup=kTRUE;
          }

          Int_t isPileup = aod->IsPileupFromSPD(3);
          if(isPileup != 0) {
            fPileUpCount->Fill(1.5);
            BisPileup=kTRUE;
          }

          if(((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
            fPileUpCount->Fill(2.5);
            BisPileup=kTRUE;
          }

          if(aod->IsIncompleteDAQ())  {
            fPileUpCount->Fill(3.5);
            BisPileup=kTRUE;
          }

          if(fabs(centrV0M-centrCL1)>7.5)  {
            fPileUpCount->Fill(4.5);
            BisPileup=kTRUE;
          }

     // check vertex consistency
     const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
     const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();

     if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1) {
            fPileUpCount->Fill(5.5);
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
            fPileUpCount->Fill(6.5);
            BisPileup=kTRUE;
          }

       // cuts on tracks
          const Int_t nTracks = aod->GetNumberOfTracks();
          Int_t multEsd = ((AliAODHeader*)aod->GetHeader())->GetNumberOfESDTracks();

          Int_t multTrk = 0;
          Int_t multTrkBefC = 0;
          Int_t multTrkTOFBefC = 0;
          Int_t multTPC = 0;

    for(Int_t it = 0; it < nTracks; it++) {
       AliAODTrack* aodTrk = (AliAODTrack*)aod->GetTrack(it);
        if(!aodTrk) {
          delete aodTrk;
          continue;
         }
//      if(aodTrk->TestFilterBit(32)){
//         multTrkBefC++;
//       if(TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
//         multTrkTOFBefC++;
//         if((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2) && (aodTrk->Pt() < 20.))
//            multTrk++;
//       }
         if(aodTrk->TestFilterBit(128))
              multTPC++;
     } // end of for (Int_t it = 0; it < nTracks; it++)

      Double_t multTPCn = multTPC;
      Double_t multEsdn = multEsd;
      Double_t multESDTPCDif = multEsdn - multTPCn*3.38;

     if(multESDTPCDif > (fRejectPileUpTight?700.:15000.)) {
        fPileUpCount->Fill(7.5);
        BisPileup=kTRUE;
      }

     if(fRejectPileUpTight) {
       if(BisPileup==kFALSE) {
              if(!fMultSelection->GetThisEventIsNotPileup()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
              if(!fMultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
              if(BisPileup) fPileUpCount->Fill(8.5);
            }
          }
	}//-> 2015
      }

 //----------- pile up rejection done ---------

  //printf("\n ... ::UserExec() is being called. 4 Step %d...  \n",stepCount);

  if(BisPileup){
    AliDebug(2,Form("\n ::UserExec Pileup event found, skipping.\n Dataset: %s \n",fDataSet.Data()));
    return;
  }

    fHist_Event_count->Fill(stepCount);
    stepCount++;

    fTPCV0M_CentDiff_vs_Vx->Fill(Vxyz[0],centrTRK-centrV0M);
    fTPCV0M_CentDiff_vs_Vy->Fill(Vxyz[1],centrTRK-centrV0M);
    fTPCV0M_CentDiff_vs_Vz->Fill(Vxyz[2],centrTRK-centrV0M);

  if(EvtCent<20){
     fAvg_Cent_vs_Vx_Cent_woCut->Fill(Vxyz[0],EvtCent);
     fAvg_Cent_vs_Vy_Cent_woCut->Fill(Vxyz[1],EvtCent);
     fAvg_Cent_vs_Vz_Cent_woCut->Fill(Vxyz[2],EvtCent);
   }
  if(EvtCent>50){
     fAvg_Cent_vs_Vx_Peri_woCut->Fill(Vxyz[0],EvtCent);
     fAvg_Cent_vs_Vy_Peri_woCut->Fill(Vxyz[1],EvtCent);
     fAvg_Cent_vs_Vz_Peri_woCut->Fill(Vxyz[2],EvtCent);
   }


 // ********** fZDCgain alpha = 0.50 instead of 0.35 *********
  AliAODZDC *aodZDC = aod->GetZDCData();
  Float_t                            fZDCGainAlpha = 0.500;
  Float_t energyZNC  = (Float_t)  (aodZDC->GetZNCEnergy());
  Float_t energyZPC  = (Float_t)  (aodZDC->GetZPCEnergy());
  Float_t energyZNA  = (Float_t)  (aodZDC->GetZNAEnergy());
  Float_t energyZPA  = (Float_t)  (aodZDC->GetZPAEnergy());
  Float_t energyZEM1 = (Float_t) (aodZDC->GetZEM1Energy());
  Float_t energyZEM2 = (Float_t) (aodZDC->GetZEM2Energy());

  const Double_t * towZNC   = aodZDC->GetZNCTowerEnergy();
  const Double_t * towZPC   = aodZDC->GetZPCTowerEnergy();
  const Double_t * towZNA   = aodZDC->GetZNATowerEnergy();
  const Double_t * towZPA   = aodZDC->GetZPATowerEnergy();

  const Double_t * towZNClg = aodZDC->GetZNCTowerEnergyLR(); // Low gain something, should not be used.
  const Double_t * towZNAlg = aodZDC->GetZNATowerEnergyLR();

  Double_t towZPClg[5] = {0.,};
  Double_t towZPAlg[5] = {0.,};

  for(Int_t it=0; it<5; it++) {
    towZPClg[it] = 8*towZPC[it];
    towZPAlg[it] = 8*towZNA[it];
  }

  //sanity: remove if any of ZDC_C_A has negetive Energy:
  if(towZNC[1]<0 || towZNC[2]<0 || towZNC[3]<0 || towZNC[4]<0)  return;
  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(towZNA[1]<0 || towZNA[2]<0 || towZNA[3]<0 || towZNA[4]<0)  return;
  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(fAnalysisSet=="recenter1"){
   for(Int_t it=0; it<5; it++) {
       fHist_ZDCA_En_Run[runindex]->Fill(EvtCent,it,towZNA[it]);
       fHist_ZDCC_En_Run[runindex]->Fill(EvtCent,it,towZNC[it]);
    }
  }

//********** Get centroid from ZDCs **************

  Double_t xyZNC[2]={999.,999.};
  Double_t xyZNA[2]={999.,999.};

  Float_t zncEnergy=0., znaEnergy=0.;

  for(Int_t i=0; i<5; i++){
     zncEnergy += towZNC[i];
     znaEnergy += towZNA[i];
  }


  Double_t AvTowerGain[8] = {1., 1., 1., 1., 1., 1., 1., 1.};

/*---------------------------------------------------------------
  Int_t CenBin = GetCenBin(centrperc);
  Double_t zvtxpos[3]={0.,0.,0.};
  fFlowEvent->GetVertexPosition(zvtxpos);
  Int_t RunNum=fFlowEvent->GetRun();
  if(fTowerEqList) {
   if(RunNum!=fCachedRunNum) {
      for(Int_t i=0; i<8; i++) {
      fTowerGainEq[i] = (TH1D*)(fTowerEqList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fhnTowerGainEqFactor[%d][%d]",RunNum,i)));
     }
   }
 }
  Bool_t fUseMCCen = kFALSE; //rihan:hardcoded
  if (fUseMCCen) {
   if(aod->GetRunNumber() < 209122)
    aodZDC->GetZNCentroidInPbPb(1380., xyZNC, xyZNA);
    else
    aodZDC->GetZNCentroidInPbPb(2510., xyZNC, xyZNA);
  }
  else {
  //set tower gain equalization, if available
  if(fTowerEqList) {
   for(Int_t i=0; i<8; i++)
   {
    if(fTowerGainEq[i])
    AvTowerGain[i] = fTowerGainEq[i]->GetBinContent(fTowerGainEq[i]->FindBin(centrperc));
   }
 }//---------------------------------------------------------------  */



  const Float_t x[4] = {-1.75,  1.75,-1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};

  Float_t numXZNC=0., numYZNC=0., denZNC=0., wZNC;
  Float_t numXZNA=0., numYZNA=0., denZNA=0., wZNA;

  for(Int_t i=0; i<4; i++)
   {
    if(towZNC[i+1]>0.)
      {
       wZNC = TMath::Power(towZNC[i+1], fZDCGainAlpha)*AvTowerGain[i];
       numXZNC += x[i]*wZNC;
       numYZNC += y[i]*wZNC;
       denZNC  += wZNC;
       }

    if(towZNA[i+1]>0.) {
       wZNA = TMath::Power(towZNA[i+1], fZDCGainAlpha)*AvTowerGain[i+4];
       numXZNA += x[i]*wZNA;
       numYZNA += y[i]*wZNA;
       denZNA  += wZNA;
       }
    }

  if(denZNC!=0) {
     xyZNC[0] = numXZNC/denZNC;
     xyZNC[1] = numYZNC/denZNC;
  }
  else{
     xyZNC[0]  = 999.;
     xyZNC[1]  = 999.;
     zncEnergy =   0.;
    }
  if(denZNA!=0) {
     xyZNA[0] = numXZNA/denZNA;
     xyZNA[1] = numYZNA/denZNA;
    }
   else{
     xyZNA[0]  = 999.;
     xyZNA[1]  = 999.;
     znaEnergy =   0.;
   }


  //----- Important: zdcA_X = -zdcA_X ---------
  xyZNA[0] = -1.*xyZNA[0];
  //===========================================


  if(sqrt(xyZNC[0]*xyZNC[0] + xyZNC[1]*xyZNC[1])>1.5)  return;
  fHist_Event_count->Fill(stepCount);
  stepCount++;

  if(sqrt(xyZNA[0]*xyZNA[0] + xyZNA[1]*xyZNA[1])>1.5)  return;
  fHist_Event_count->Fill(stepCount);
  stepCount++;

 // ----------------  ZDC multiplicities posted at the ende of this code-----------------------

  fHist_EventperRun->Fill(runindex);


 //----- calculate RefMult for event ---------
    AliFlowTrackSimple*   pTrack = NULL;
    Int_t iTracks = fEvent->NumberOfTracks();

    Double_t    Qnx_TPC[3]  = {0,};
    Double_t    Qny_TPC[3]  = {0,};
    Double_t    psi2,dPhi,dPt,dEta;
    Int_t                   ipTBin;
    Double_t                 pTwgt;
    Int_t             nRefMult = 0;
    Double_t          npoiMult = 0;

    Int_t cIndex = -1;
    if(EvtCent<5.0){ cIndex  = 0;}
    else if(EvtCent>=5.0 && EvtCent<10){
          cIndex  = 1;
    }
    else {
          cIndex = abs(EvtCent/10.0) + 1;
    }


    Double_t dUx = 0;
    Double_t dUy = 0;

 if(fAnalysisSet=="recenter2") {
    for(Int_t i=0; i<iTracks; i++)
     {
      pTrack     =  fEvent->GetTrack(i);
      if (!pTrack)             continue;
      dPhi      =         pTrack->Phi();
      dPt       =         pTrack-> Pt();
      dEta      =         pTrack->Eta();
      //dCharge   =      pTrack->Charge();
      if(fabs(dEta)>0.8)       continue;
      if(dPt<0.15 || dPt>10.0) continue;
      nRefMult++;
      if(dPt<0.20)             continue;

      //fHistQA_etaphi->Fill(dPhi,dEta);
      ipTBin = fFB_Efficiency_Cent[cIndex]->FindBin(dPt);
      pTwgt  = 1.0/fFB_Efficiency_Cent[cIndex]->GetBinContent(ipTBin);

      Qnx_TPC[0] += TMath::Cos(2.*dPhi)*pTwgt;
      Qny_TPC[0] += TMath::Sin(2.*dPhi)*pTwgt;

      npoiMult   += pTwgt;
      //cout<<"c = "<<EvtCent<<" indC = "<<cIndex<<"  "<<npoiMult<<" pt = "<<dPt<<" weight = "<<pTwgt<<endl;
    }//track loop ends
  } //----- call track loop only for recenter2 -------


  if(npoiMult>0){
     dUx = Qnx_TPC[0]/npoiMult;
     dUy = Qny_TPC[0]/npoiMult;
  }
  else{
     dUx = 0;
     dUy = 0;
  }

   //double psi2 = 0.5*TMath::ATan2(Qny_TPC[0],Qnx_TPC[0]);
   //if(psi2<0) psi2 += TMath::Pi();
   //fHist_EventPlane2->Fill(psi2);
   //fill q vectors for TPC recentering
   //fHist_QnxRecent->Fill(EvtCent,runindex,(Qnx_TPC[0]/(double)npoiMult));
   //fHist_QnyRecent->Fill(EvtCent,runindex,(Qny_TPC[0]/(double)npoiMult));
 
  fHist_Cent_count2->        Fill(EvtCent);



  if(EvtCent<20){
     fAvg_Cent_vs_Vx_Cent_wCuts->Fill(Vxyz[0],EvtCent);
     fAvg_Cent_vs_Vy_Cent_wCuts->Fill(Vxyz[1],EvtCent);
     fAvg_Cent_vs_Vz_Cent_wCuts->Fill(Vxyz[2],EvtCent);
   }
  if(EvtCent>50){
     fAvg_Cent_vs_Vx_Peri_wCuts->Fill(Vxyz[0],EvtCent);
     fAvg_Cent_vs_Vy_Peri_wCuts->Fill(Vxyz[1],EvtCent);
     fAvg_Cent_vs_Vz_Peri_wCuts->Fill(Vxyz[2],EvtCent);
   }

  
  Int_t nTracks = aod->GetNumberOfTracks();      //number of AOD tracks

  fHist_Vertex_Vz->Fill(Vxyz[2]);
  fHist_Vertex_XY->Fill(Vxyz[0],Vxyz[1]);

  fHist_Vx_vs_runnum    ->Fill(runindex,Vxyz[0]);
  fHist_Vy_vs_runnum    ->Fill(runindex,Vxyz[1]);
  fHist_Vz_vs_runnum    ->Fill(runindex,Vxyz[2]);
  fHist_tracks_vs_runnum->Fill(runindex,nTracks);
 



  psi2 = TMath::ATan2(xyZNC[1],xyZNC[0]);
  if(psi2<0) psi2 += 2.*TMath::Pi();
  fHist_Psi1_zdnC->Fill(psi2);

  psi2 = TMath::ATan2(xyZNA[1],xyZNA[0]);
  if(psi2<0) psi2 += 2.*TMath::Pi();
  fHist_Psi1_zdnA->Fill(psi2);



  Int_t indexVx = fHist_Vx_ArrayFinder->FindBin(Vxyz[0]);
  Int_t indexVy = fHist_Vy_ArrayFinder->FindBin(Vxyz[1]);
  Int_t indexVz = fHist_Vz_ArrayFinder->FindBin(Vxyz[2]);


 /*if(fAnalysisSet=="recenter1"){
  Double_t tVertexBin1 = (Double_t) (indexVy-1)*vxBin + (Double_t)indexVx - 0.5 ; // 
  fHist_vBincount->Fill(tVertexBin1);
  fHist_znCx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[0]); //EvtCent
  fHist_znCy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNC[1]);
  fHist_znAx_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[0]);
  fHist_znAy_V0_VxVy[runindex][indexVz-1]->Fill(tVertexBin1,EvtCent,xyZNA[1]);
  } */

  Int_t    tVertexBin2 =  (indexVy-1)*vxBin + indexVx;         // 'Int_t' because bin starts with 1.
  Int_t    tCentBin    =  abs(EvtCent) + 1; 
  Double_t CentWgt     =               1.0;


  //printf("\n ... ::UserExec() is being called. Step 6 %d...  \n",stepCount);

  Double_t meanCx[2] = {0,};
  Double_t meanCy[2] = {0,};
  Double_t meanAx[2] = {0,};
  Double_t meanAy[2] = {0,};

  if(fAnalysisSet=="recenter2") {

   Double_t tVertexBinD =  (Double_t)tVertexBin2 - 0.5;
   fHist_vBincount->Fill(tVertexBinD);

   Double_t  FillVsWith[5] = {EvtCent,static_cast<Double_t>(nRefMult), Vxyz[0], Vxyz[1], Vxyz[2]};
   Double_t  FillValue1[4] = {xyZNA[0],xyZNC[0],xyZNA[1],xyZNC[1]};
   Double_t FillValue11[4] = {xyZNA[0]*xyZNC[0],xyZNA[1]*xyZNC[1],xyZNC[0]*xyZNA[1],xyZNC[1]*xyZNA[0]}; //XaXc,YaYc,XcYa,YcXa

  //fill the uncorrected Qns first:
   for(int i=0;i<4;i++){
     for(int j=0;j<5;j++){
        fHist_X_vs_Obs_before[i][j]->Fill(FillVsWith[j],FillValue1[i]);
        fHist_XX_vs_Obs_before[i][j]->Fill(FillVsWith[j],FillValue11[i]);
     }
   }

   meanCx[0] = fHist_znCx_V0_VxVy[runindex][indexVz-1]->GetBinContent(tVertexBin2,tCentBin); 
   meanCy[0] = fHist_znCy_V0_VxVy[runindex][indexVz-1]->GetBinContent(tVertexBin2,tCentBin);
   meanAx[0] = fHist_znAx_V0_VxVy[runindex][indexVz-1]->GetBinContent(tVertexBin2,tCentBin);
   meanAy[0] = fHist_znAy_V0_VxVy[runindex][indexVz-1]->GetBinContent(tVertexBin2,tCentBin); 





   xyZNC[0] = xyZNC[0] - meanCx[0];
   xyZNC[1] = xyZNC[1] - meanCy[0];
   xyZNA[0] = xyZNA[0] - meanAx[0];
   xyZNA[1] = xyZNA[1] - meanAy[0];

   Double_t  FillValue2[4] = {xyZNA[0],xyZNC[0],xyZNA[1],xyZNC[1]};
   Double_t FillValue21[4] = {xyZNA[0]*xyZNC[0],xyZNA[1]*xyZNC[1],xyZNC[0]*xyZNA[1],xyZNC[1]*xyZNA[0]}; //XaXc,YaYc,XcYa,YcXa

   //fill QA after recenter 1:
   for(int i=0;i<4;i++){
     for(int j=0;j<5;j++){
        fHist_X_vs_Obs_after1[i][j]->Fill(FillVsWith[j],FillValue2[i]);
       fHist_XX_vs_Obs_after1[i][j]->Fill(FillVsWith[j],FillValue21[i]);
     }
   }

    psi2 = TMath::ATan2(xyZNC[1],xyZNC[0]);
    if(psi2<0) psi2 += 2.*TMath::Pi();
    fHist_Psi1_zdnC_after1->Fill(psi2);

    psi2 = TMath::ATan2(xyZNA[1],xyZNA[0]);
    if(psi2<0) psi2 += 2.*TMath::Pi();
    fHist_Psi1_zdnA_after1->Fill(psi2);

    //fill results histograms:
    fHist_v2xV1_ZDN_Norm_cosXX ->Fill(EvtCent, dUx*xyZNA[0]*xyZNC[0],CentWgt);
    fHist_v2xV1_ZDN_Cent_cosXX ->Fill(EvtCent, dUx*xyZNA[0]*xyZNC[0],CentWgt);
    fHist_v2xV1_ZDN_Refm_cosXX ->Fill(nRefMult,dUx*xyZNA[0]*xyZNC[0],CentWgt);

    fHist_v2xV1_ZDN_Norm_cosYY ->Fill(EvtCent, dUx*xyZNA[1]*xyZNC[1],CentWgt);
    fHist_v2xV1_ZDN_Cent_cosYY ->Fill(EvtCent, dUx*xyZNA[1]*xyZNC[1],CentWgt);
    fHist_v2xV1_ZDN_Refm_cosYY ->Fill(nRefMult,dUx*xyZNA[1]*xyZNC[1],CentWgt);

    fHist_v2xV1_ZDN_Norm_sinXY ->Fill(EvtCent, dUy*xyZNA[0]*xyZNC[1],CentWgt);
    fHist_v2xV1_ZDN_Cent_sinXY ->Fill(EvtCent, dUy*xyZNA[0]*xyZNC[1],CentWgt);
    fHist_v2xV1_ZDN_Refm_sinXY ->Fill(nRefMult,dUy*xyZNA[0]*xyZNC[1],CentWgt);

    fHist_v2xV1_ZDN_Norm_sinYX ->Fill(EvtCent, dUy*xyZNA[1]*xyZNC[0],CentWgt);
    fHist_v2xV1_ZDN_Cent_sinYX ->Fill(EvtCent, dUy*xyZNA[1]*xyZNC[0],CentWgt);
    fHist_v2xV1_ZDN_Refm_sinYX ->Fill(nRefMult,dUy*xyZNA[1]*xyZNC[0],CentWgt);

    fHist_ZDN_resol_Norm_cos ->Fill(EvtCent, xyZNA[0]*xyZNC[0],CentWgt);
    fHist_ZDN_resol_Cent_cos ->Fill(EvtCent, xyZNA[0]*xyZNC[0],CentWgt);
    fHist_ZDN_resol_Refm_cos ->Fill(nRefMult,xyZNA[0]*xyZNC[0],CentWgt);

    fHist_ZDN_resol_Norm_sin ->Fill(EvtCent, xyZNA[1]*xyZNC[1],CentWgt);
    fHist_ZDN_resol_Cent_sin ->Fill(EvtCent, xyZNA[1]*xyZNC[1],CentWgt);
    fHist_ZDN_resol_Refm_sin ->Fill(nRefMult,xyZNA[1]*xyZNC[1],CentWgt);

    Double_t fullTerm = dUx*xyZNA[0]*xyZNC[0]-dUx*xyZNA[1]*xyZNC[1]+dUy*xyZNA[0]*xyZNC[1]+dUy*xyZNA[1]*xyZNC[0];
    Double_t fullReso = xyZNA[0]*xyZNC[0]+xyZNA[1]*xyZNC[1];

    fHist_v2xV1_ZDN_Norm_All ->Fill(EvtCent,fullTerm,CentWgt);
    fHist_v2xV1_ZDN_Cent_All ->Fill(EvtCent,fullTerm,CentWgt);
    fHist_v2xV1_ZDN_Refm_All ->Fill(nRefMult,fullTerm,CentWgt);
    fHist_ZDN_resol_Norm_All ->Fill(EvtCent,fullReso,CentWgt);
    fHist_ZDN_resol_Cent_All ->Fill(EvtCent,fullReso,CentWgt);    
    fHist_ZDN_resol_Refm_All ->Fill(nRefMult,fullReso,CentWgt);  
  }



  //fHist_ZDCn_A_XYvsRun->Fill(xyZNA[0],xyZNA[1],runindex);
  //fHist_ZDCn_C_XYvsRun->Fill(xyZNC[0],xyZNC[1],runindex);


  
  //if(fievent%10==0){
  //std::cout<<fievent<<" cTPC= "<<EvtCent<<" cV0M = "<<centrV0M<<" cWgt = "<<CentWgt<<" vz= "<<Vxyz[2]<<"\tCx= "<<meanCx[0]
  //<<"\tCy= "<<meanCy[0]<<"\tAx= "<<meanAx[0]<<"\tAy= "<<meanAy[0]<<"\tRefm= "<<nRefMult<<std::endl;  } 
 

  PostData(1,fListHistos);

  if(fAnalysisSet=="recenter2"){
    PostData(2,fZListDummy); 
  }
  else{
    PostData(2,fZDCESEList); 
  }



  fievent++;

}




void AliAnalysisTaskVnZDC::Terminate(Option_t *)
{
  // Called once at the end of the query
       /* AliFlowAnalysisIDCSP *fSPTerm = new AliFlowAnalysisIDCSP();
	fListHistos = (TList*) GetOutputData(1);
        if(fListHistos)
	  {
	   fSPTerm->GetOutputHistograms(fListHistos);
	   fSPTerm->Finish();
	   PostData(1,fListHistos);
	  }
        else
	{
          std::cout << "histgram list pointer is empty in Scalar Product" << endl;
         } */
      AliDebug(2,"\n ... AliAnalysisTaskVnZDC::Terminate() is being called ...  \n");
}


double AliAnalysisTaskVnZDC::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
  // calculate sqrt of weighted distance to other vertex
  if(!v0 || !v1) {
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
  return dist;}
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
  +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;
}

 Bool_t AliAnalysisTaskVnZDC::plpMV(const AliAODEvent* aod)
 {  // check for multi-vertexer pile-up
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2    = 5.0;
  const double kMinWDist      = 15;

  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;

  int nPlp = 0;

  if(!(nPlp=aod->GetNumberOfPileupVerticesTracks()))
  return kFALSE;

  vtPrm = aod->GetPrimaryVertex();
  if(vtPrm == aod->GetPrimaryVertexSPD())
  return kTRUE;  // there are pile-up vertices but no primary

  //int bcPrim = vtPrm->GetBC();

  for(int ipl=0;ipl<nPlp;ipl++) {
     vtPlp = (const AliVVertex*)aod->GetPileupVertexTracks(ipl);
      if(vtPlp->GetNContributors() < kMinPlpContrib) continue;
      if(vtPlp->GetChi2perNDF()    > kMaxPlpChi2)    continue;
    //int bcPlp = vtPlp->GetBC();
    //if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2)
    // return kTRUE; // pile-up from other BC

    double wDst = GetWDist(vtPrm,vtPlp);
      if(wDst<kMinWDist)        continue;

    return kTRUE; // pile-up: well separated vertices
    }

   return kFALSE;
 }




