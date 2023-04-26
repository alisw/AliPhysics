#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "TString.h"
#include "TParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"

#include "AliAODpidUtil.h"
#include "AliAODHeader.h"

#include "AliAnalysisTaskParticleEffDYM3.h"




ClassImp(AliAnalysisTaskParticleEffDYM3)
//ClassImp(AliAnalysisTaskParticleEff)

double fV1M3[3];


void AliAnalysisTaskParticleEffDYM3::SetFB(int fb)
{
  fFB = fb;
}

void AliAnalysisTaskParticleEffDYM3::SetPidMethod(PidMethod method)
{
  fPidMethod = method;
}

int AliAnalysisTaskParticleEffDYM3::GetPidMethod()
{
  return (int)fPidMethod;
}

void AliAnalysisTaskParticleEffDYM3::SetPidMethod(int method)
{
  switch(method){
  case 0: fPidMethod=kNSigma;
    break;
  case 1: fPidMethod=kNSigmaNoDoubleCounting;
    break;
  case 2: fPidMethod=kExclusivePID;
    break;
  case 3: fPidMethod=kExclusivePIDDiffRejection;
    break;
  }
}

//_______________________________________________________

AliAnalysisTaskParticleEffDYM3::AliAnalysisTaskParticleEffDYM3(TString name, int pidMethod, int filterbit) :
  AliAnalysisTaskSE(name), centrality(0), fHistoList(0),  fMassInvLambdaPass(0),fMassInvAntiLambdaPass(0), fMassInvLambdaFail(0), fMassInvAntiLambdaFail(0),fYLambda(0),fPtLambda(0), fYAntiLambda(0),fPtAntiLambda(0), fCutsLambda(0), fCutsAntiLambda(0), fTruePtLambdaMC(0), fRecPtLambdaMC(0), fTruePtAntiLambdaMC(0),fRecPtAntiLambdaMC(0), fMassInvXimPass(0),fMassInvXipPass(0), fMassInvXimFail(0), fMassInvXipFail(0),fYXim(0),fPtXim(0), fYXip(0),fPtXip(0), fCutsXim(0), fCutsXip(0), recoParticleArrayXi(0), fTruePtXimMC(0), fRecPtXimMC(0), fTruePtXipMC(0), fRecPtXipMC(0), fDCAtoPrimVtx(0), fIfAliEventCuts(kFALSE), fFB(128), fPidMethod(kExclusivePIDDiffRejection),  fEstEventMult(kV0M),fIfXiAnalysis(kFALSE), fpidResponse(0), fAODpidUtil(0), fEventCuts(0)

{

  for(Int_t i = 0; i < MULTBINS*PARTTYPES; i++)  {
    for(Int_t chg=0;chg<2;chg++){
      fGeneratedMCPrimaries[i][chg] = NULL;
      fMCPrimariesThatAreReconstructed[i][chg] = NULL;
      fMCPrimariesThatAreReconstructedNoNsigma[i][chg] = NULL;
      fReconstructedAfterCuts[i][chg] = NULL;
      fReconstructedNotPrimaries[i][chg] = NULL;
      fReconstructedPrimaries[i][chg] = NULL;
      fContamination[i][chg] = NULL;

      fPrimVsDCA[i][chg] = NULL;
      fSecWeakVsDCA[i][chg] = NULL;
      fSecMatVsDCA[i][chg] = NULL;
      fFakeVsDCA[i][chg] = NULL;
      
      fPrimVsCosPointingAngle[i][chg] = NULL;
      fSecWeakVsCosPointingAngle[i][chg] = NULL;
      fSecMatVsCosPointingAngle[i][chg] = NULL;
      fFakeVsCosPointingAngle[i][chg] = NULL;

      fMCPrimariesThatAreReconstructed4D[i][chg] = NULL;
      fGeneratedMCPrimaries4D[i][chg] = NULL;
    }
  }
 for ( Int_t i = 0; i < 11; i++) {
    if(i<4) fHistEv[i] = NULL;
    fHistQA[i] = NULL;
    if(i<3) fHistQA2D[i] = NULL;
    if(i<4) fHistP[i]=NULL;
  }

  /* init track cuts */
  //fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
  /*fTrackCuts =  AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    if( !fTrackCuts ) return;
    fTrackCuts->SetMinNClustersTPC(70);*/
  
  


  //DefineInput(0, TChain::Class());
  //DefineOutput(0, TTree::Class()); 
  DefineOutput(1, TList::Class());
}

//_______________________________________________________

AliAnalysisTaskParticleEffDYM3::~AliAnalysisTaskParticleEffDYM3()
{
  /* if(centrality) delete centrality;
s128     if(fHistoList) delete fHistoList;
     if(vertex) delete vertex;
     if(vtxSPD) delete vtxSPD;*/

  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis)
    delete fHistoList;
}

//_______________________________________________________

void AliAnalysisTaskParticleEffDYM3::UserCreateOutputObjects()
{

  /* create output */
  fHistoList = new TList();
  fHistoList->SetOwner(kTRUE);
   
  TString hname1, hname2, hname3, hname4, hname5, hname6, hname7, hname8, hname9, hname10, hname11,hname12,hname13,hname14,hname15,hname16,hname17,hname18,hname19,hname20,hname21,hname22,hname23,hname24,hname25,hname26 ;
  
  TString htitle1, htitle2, htitle3, htitle4,htitle6,htitle5,htitle7,htitle8,htitle9,htitle10,htitle11,htitle12,htitle13,htitle14,htitle15,htitle16,htitle17,htitle18,htitle19,htitle20,htitle21,htitle22,htitle23,htitle24,htitle25,htitle26;
  
  TString hname1M, hname2M, hname3M, hname4M, hname5M, hname;
  
  TString htitle1M, htitle2M, htitle3M, htitle4M, htitle5M, htitle;

  TString parttypename = "None";

  for(Int_t j = 0; j < PARTTYPES; j++)  {
    if (j==0) parttypename="All";
    else if (j==1) parttypename="Pion";
    else if (j==2) parttypename="Kaon";
    else if (j==3) parttypename="Proton";
    else if (j==4) parttypename="Lambda";

    for(Int_t i = 0; i < MULTBINS; i++)  {
      hname1  = "hGeneratedMCPrimariesEffM"; hname1+=i; hname1+=parttypename;
      htitle1 = "Kinematic level Y_pT (prim only) M"; htitle1+=i; htitle1+=parttypename;
      fGeneratedMCPrimaries[i*PARTTYPES+j][0] = new TH2F(hname1.Data(),htitle1.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname1+="Minus";htitle1+="Minus";
      fGeneratedMCPrimaries[i*PARTTYPES+j][1] = new TH2F(hname1.Data(),htitle1.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname3  = "hMCPrimariesThatAreReconstructedM"; hname3+=i; hname3+=parttypename;
      htitle3 = "Reconstructed level Y_pT (prim only) M"; htitle3+=i; htitle3+=parttypename;
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][0] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname3+="Minus";htitle3+="Minus";
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][1] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname3  = "hMCPrimariesThatAreReconstructedNoNsigmaM"; hname3+=i; hname3+=parttypename;
      htitle3 = "Reconstructed level Y_pT (prim only) no Nsigma cut only PDG M"; htitle3+=i; htitle3+=parttypename;
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][0] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname3+="Minus";htitle3+="Minus";
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][1] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname2  = "hHistoReconstructedAfterCutsM"; hname2+=i; hname2+=parttypename;
      htitle2 = "Total Reconstructed tracks M "; htitle2+=i; htitle2+=parttypename;
      fReconstructedAfterCuts[i*PARTTYPES+j][0] = new TH2F(hname2.Data(),htitle2.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname2+="Minus";htitle2+="Minus";
      fReconstructedAfterCuts[i*PARTTYPES+j][1] = new TH2F(hname2.Data(),htitle2.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname4  = "hHistoReconstructedNotPrimariesM"; hname4+=i; hname4+=parttypename;
      htitle4 = "Reconstructed level Y_pT (not primaries) M"; htitle4+=i; htitle4+=parttypename;
      fReconstructedNotPrimaries[i*PARTTYPES+j][0] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname4+="Minus";htitle4+="Minus";
      fReconstructedNotPrimaries[i*PARTTYPES+j][1] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname4  = "hHistoReconstructedPrimariesM"; hname4+=i; hname4+=parttypename;
      htitle4 = "Reconstructed level Y_pT (primaries) M"; htitle4+=i; htitle4+=parttypename;
      fReconstructedPrimaries[i*PARTTYPES+j][0] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname4+="Minus";htitle4+="Minus";
      fReconstructedPrimaries[i*PARTTYPES+j][1] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname5  = "hContaminationM"; hname5+=i; hname5+=parttypename;
      htitle5 = "Contamination M"; htitle5+=i; htitle5+=parttypename;
      fContamination[i*PARTTYPES+j][0] = new TH2F(hname5.Data(),htitle5.Data(),6000, -3000, 3000.,100,0.,10.0); //50
      hname5+="Minus";htitle5+="Minus";
      fContamination[i*PARTTYPES+j][1] = new TH2F(hname5.Data(),htitle5.Data(),6000, -3000, 3000.,100,0.,10.0); //50

      fReconstructedAfterCuts[i*PARTTYPES+j][0]->Sumw2();
      fReconstructedNotPrimaries[i*PARTTYPES+j][0]->Sumw2();
      fReconstructedPrimaries[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][0]->Sumw2();
      fGeneratedMCPrimaries[i*PARTTYPES+j][0]->Sumw2();
      fContamination[i*PARTTYPES+j][0]->Sumw2();
      fContamination[i*PARTTYPES+j][1]->Sumw2();

     
      hname6 = "hPrimVsDCAM"; hname6+=i; hname6+=parttypename;
      htitle6 = "Primaries vs DCA M"; htitle6+=i; htitle6+=parttypename;      
      fPrimVsDCA[i*PARTTYPES+j][0] = new TH2F(hname6.Data(),htitle6.Data(),500,0,5.0,7,0.5,4);
      hname6+="Minus"; htitle6+="Minus";
      fPrimVsDCA[i*PARTTYPES+j][1] = new TH2F(hname6.Data(),htitle6.Data(),500,0,5.0,7,0.5,4);
      hname7 = "hSecWeakVsDCAM"; hname7+=i; hname7+=parttypename;
      htitle7 = "Sec. weak decay vs DCA M"; htitle7+=i; htitle7+=parttypename;  
      fSecWeakVsDCA[i*PARTTYPES+j][0] = new TH2F(hname7.Data(),htitle7.Data(),500,0,5.0,7,0.5,4);
      hname7+="Minus"; htitle7+="Minus";
      fSecWeakVsDCA[i*PARTTYPES+j][1] = new TH2F(hname7.Data(),htitle7.Data(),500,0,5.0,7,0.5,4);
      hname8 = "hSecMatVsDCAM"; hname8+=i; hname8+=parttypename;
      htitle8 = "Sec. material vs DCA M"; htitle8+=i; htitle8+=parttypename;
      fSecMatVsDCA[i*PARTTYPES+j][0] = new TH2F(hname8.Data(),htitle8.Data(),500,0,5.0,7,0.5,4);
      hname8+="Minus"; htitle8+="Minus";
      fSecMatVsDCA[i*PARTTYPES+j][1] = new TH2F(hname8.Data(),htitle8.Data(),500,0,5.0,7,0.5,4);
      hname9 = "hFakeVsDCAM"; hname9+=i; hname9+=parttypename;
      htitle9 = "Fake vs DCA M"; htitle9+=i; htitle9+=parttypename;
      fFakeVsDCA[i*PARTTYPES+j][0] = new TH2F(hname9.Data(),htitle9.Data(),500,0,5.0,7,0.5,4);
      hname9+="Minus"; htitle9+="Minus";
      fFakeVsDCA[i*PARTTYPES+j][1] = new TH2F(hname9.Data(),htitle9.Data(),500,0,5.0,7,0.5,4);


      hname10 = "hPrimVsCosPointingAngle"; hname10+=i; hname10+=parttypename;
      htitle10 = "Primaries vs CosPointingAngle M"; htitle10+=i; htitle10+=parttypename;      
      fPrimVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname10.Data(),htitle10.Data(),200,0.95,1.0,7,0.5,4);
      hname10+="Minus"; htitle10+="Minus";
      fPrimVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname10.Data(),htitle10.Data(),200,0.95,1.0,7,0.5,4);
      hname11 = "hSecWeakVsCosPointingAngleM"; hname11+=i; hname11+=parttypename;
      htitle11 = "Sec. weak decay vs CosPointingAngle M"; htitle11+=i; htitle11+=parttypename;  
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname11.Data(),htitle11.Data(),200,0.95,1.0,7,0.5,4);
      hname11+="Minus"; htitle11+="Minus";
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname11.Data(),htitle11.Data(),200,0.95,1.0,7,0.5,4);
      hname12 = "hSecMatVsCosPointingAngleM"; hname12+=i; hname12+=parttypename;
      htitle12 = "Sec. material vs CosPointingAngle M"; htitle12+=i; htitle12+=parttypename;
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname12.Data(),htitle12.Data(),200,0.95,1.0,7,0.5,4);
      hname12+="Minus"; htitle12+="Minus";
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname12.Data(),htitle12.Data(),200,0.95,1.0,7,0.5,4);
      hname13 = "hFakeVsCosPointingAngleM"; hname13+=i; hname13+=parttypename;
      htitle13 = "Fake vs CosPointingAngle M"; htitle13+=i; htitle13+=parttypename;
      fFakeVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname13.Data(),htitle13.Data(),200,0.95,1.0,7,0.5,4);
      hname13+="Minus"; htitle13+="Minus";
      fFakeVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname13.Data(),htitle13.Data(),200,0.95,1.0,7,0.5,4);


      hname14 = "hPrimVsDecayRadius"; hname14+=i; hname14+=parttypename;
      htitle14 = "Primaries vs DecayRadius M"; htitle14+=i; htitle14+=parttypename;      
      fPrimVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname14.Data(),htitle14.Data(),500,0,5.0,7,0.5,4);
      hname14+="Minus"; htitle14+="Minus";
      fPrimVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname14.Data(),htitle14.Data(),500,0,5.0,7,0.5,4);
      hname15 = "hSecWeakVsDecayRadiusM"; hname15+=i; hname15+=parttypename;
      htitle15 = "Sec. weak decay vs DecayRadius M"; htitle15+=i; htitle15+=parttypename;  
      fSecWeakVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname15.Data(),htitle15.Data(),500,0,5.0,7,0.5,4);
      hname15+="Minus"; htitle15+="Minus";
      fSecWeakVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname15.Data(),htitle15.Data(),500,0,5.0,7,0.5,4);
      hname16 = "hSecMatVsDecayRadiusM"; hname16+=i; hname16+=parttypename;
      htitle16 = "Sec. material vs DecayRadius M"; htitle16+=i; htitle16+=parttypename;
      fSecMatVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname16.Data(),htitle16.Data(),500,0,5.0,7,0.5,4);
      hname16+="Minus"; htitle16+="Minus";
      fSecMatVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname16.Data(),htitle16.Data(),500,0,5.0,7,0.5,4);
      hname17 = "hFakeVsDecayRadiusM"; hname17+=i; hname17+=parttypename;
      htitle17 = "Fake vs DecayRadius M"; htitle17+=i; htitle17+=parttypename;
      fFakeVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);
      hname17+="Minus"; htitle17+="Minus";
      fFakeVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);
      
      hname17 = "hAllVsDecayRadiusM"; hname17+=i; hname17+=parttypename;
      htitle17 = "All vs DecayRadius M"; htitle17+=i; htitle17+=parttypename;
      fAllVsDecayRadius[i*PARTTYPES+j][0] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);
      hname17+="Minus"; htitle17+="Minus";
      fAllVsDecayRadius[i*PARTTYPES+j][1] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);

      hname17 = "hAllVsCosPointingAngleM"; hname17+=i; hname17+=parttypename;
      htitle17 = "All vs CosPointingAngle M"; htitle17+=i; htitle17+=parttypename;
      fAllVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname17.Data(),htitle17.Data(),200,0.95,1.0,7,0.5,4);
      hname17+="Minus"; htitle17+="Minus";
      fAllVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname17.Data(),htitle17.Data(),200,0.95,1.0,7,0.5,4);

      hname17 = "hAllVsDCAM";hname17+=i; hname17+=parttypename;
      htitle17 = "All vs DCA M";htitle17+=i; htitle17+=parttypename;
      fAllVsDCA[i*PARTTYPES+j][0] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);
      hname17+="Minus"; htitle17+="Minus";
      fAllVsDCA[i*PARTTYPES+j][1] = new TH2F(hname17.Data(),htitle17.Data(),500,0,5.0,7,0.5,4);



      Double_t min[]={-1.5,0,-10,0};
      Double_t max[]={1.5,10,10,2*TMath::Pi()};
      Int_t nbins[]={20,100,10,20};
      
      hname1  = "hGeneratedMCPrimariesEff4DM"; hname1+=i; hname1+=parttypename;
      htitle1 = "Kinematic level Y_pT (prim only) M"; htitle1+=i; htitle1+=parttypename;
      fGeneratedMCPrimaries4D[i*PARTTYPES+j][0] = new THnSparseF(hname1.Data(),htitle1.Data(),4,nbins,min,max);
      hname1+="Minus";htitle1+="Minus";
      fGeneratedMCPrimaries4D[i*PARTTYPES+j][1] = new THnSparseF(hname1.Data(),htitle1.Data(),4,nbins,min,max);

      hname3  = "hMCPrimariesThatAreReconstructed4DM"; hname3+=i; hname3+=parttypename;
      htitle3 = "Reconstructed level Y_pT (prim only) M"; htitle3+=i; htitle3+=parttypename;
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][0] = new THnSparseF(hname3.Data(),htitle3.Data(),4,nbins,min,max);
      hname3+="Minus";htitle3+="Minus";
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][1] = new THnSparseF(hname3.Data(),htitle3.Data(),4,nbins,min,max);
      
      fGeneratedMCPrimaries4D[i*PARTTYPES+j][0]->Sumw2();
      fGeneratedMCPrimaries4D[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][1]->Sumw2();
      
      fPrimVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fSecWeakVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fSecMatVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fFakeVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fPrimVsDCA[i*PARTTYPES+j][1]->Sumw2();
      fSecWeakVsDCA[i*PARTTYPES+j][1]->Sumw2();
      fSecMatVsDCA[i*PARTTYPES+j][1]->Sumw2();
      fFakeVsDCA[i*PARTTYPES+j][1]->Sumw2();
      
      fPrimVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fFakeVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fPrimVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();
      fFakeVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();

      fPrimVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fSecWeakVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fSecMatVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fFakeVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fPrimVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();
      fSecWeakVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();
      fSecMatVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();
      fFakeVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();

      fAllVsDCA[i*PARTTYPES+j][0]->Sumw2();
      fAllVsDCA[i*PARTTYPES+j][1]->Sumw2();
      fAllVsCosPointingAngle[i*PARTTYPES+j][0]->Sumw2();
      fAllVsCosPointingAngle[i*PARTTYPES+j][1]->Sumw2();
      fAllVsDecayRadius[i*PARTTYPES+j][0]->Sumw2();
      fAllVsDecayRadius[i*PARTTYPES+j][1]->Sumw2();


    }
    
    hname  = "pidTPCdEdx";  hname+=parttypename;
    htitle = parttypename + " TPC dEdx vs. momentum";
    fHistQAPID[0][j][0] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[0][j][1] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    hname  = "pidTOFTime";  hname+=parttypename;
    htitle = parttypename + " TOF Time vs. momentum";
    fHistQAPID[1][j][0] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[1][j][1] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    hname  = "pidTOFNSigma";  hname+=parttypename;
    htitle = parttypename + " TOF NSigma vs. momentum";
    fHistQAPID[2][j][0]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[2][j][1]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCNSigma";  hname+=parttypename;
    htitle = parttypename + " TPC NSigma vs. momentum";
    fHistQAPID[3][j][0] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[3][j][1] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCTOFNSigma";  hname+=parttypename;
    htitle = parttypename + " TPC vs TOF NSigma";
    fHistQAPID[4][j][0] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[4][j][1] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);

    hname  = "pidTPCdEdxFail";  hname+=parttypename;
    htitle = parttypename + " TPC dEdx vs. momentum Fail";
    fHistQAPIDFail[0][j][0] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[0][j][1] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    hname  = "pidTOFTimeFail";  hname+=parttypename;
    htitle = parttypename + " TOF Time vs. momentum Fail";
    fHistQAPIDFail[1][j][0] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[1][j][1] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    hname  = "pidTOFNSigmaFail";  hname+=parttypename;
    htitle = parttypename + " TOF NSigma vs. momentum Fail";
    fHistQAPIDFail[2][j][0]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[2][j][1]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCNSigmaFail";  hname+=parttypename;
    htitle = parttypename + " TPC NSigma vs. momentum Fail";
    fHistQAPIDFail[3][j][0] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[3][j][1] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCTOFNSigmaFail";  hname+=parttypename;
    htitle = parttypename + " TPC vs TOF NSigma Fail";
    fHistQAPIDFail[4][j][0] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[4][j][1] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);
  }

  fHistEv[0] = new TH1F("fHistEv", "Multiplicity", 100, 0, 5000);
  fHistEv[1] = new TH1F("fHistEvFB16", "Multiplicity FB16", 100, 0, 200);
  fHistEv[2] = new TH1F("fHistEvFB96", "Multiplicity FB96", 100, 0, 200);
  fHistEv[3] = new TH1F("fHistEvFB128", "Multiplicity FB128", 100, 0, 200);
  for(Int_t i = 0; i < 4; i++)
    fHistoList->Add(fHistEv[i]);


  for(Int_t i = 0; i < 4; i++)  {
    hname = "fHistEventCutsM";
    hname+= i;
    

    fHistEvCuts[i] = new TH1F(hname,Form("Event Cuts M%d",i) , 10, 0, 11);
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(1,"All");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(2,"NoVertex");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(3,"PileUp");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(4,"z-vertex>10");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(5,"eventspions");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(6,"eventskaons");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(7,"eventsprotons");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(8,"eventspions T");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(9,"eventskaons T");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(10,"eventsprotons T");
    fHistoList->Add(fHistEvCuts[i]);

    for(Int_t chg=0;chg<2;chg++){
      hname  = "hMisidentificationM"; hname+=i; if(chg==0) hname+="Plus"; else hname+="Minus"; 
      htitle = "Misidentification Fraction M"; htitle+=i; if(chg==0) htitle+="Plus"; else htitle+="Minus";
      fMisidentification[i][chg] = new TH2F(hname.Data(),htitle.Data(), 3, 0.5, 3.5, 4 , 0, 4);
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(1,"Pions, MC");
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(2,"Kaons, MC");
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(3,"Protons, MC");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(1,"Pions, Data");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(2,"Kaons, Data");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(3,"Protons, Data");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(4,"Other, Data");
      fHistoList->Add(fMisidentification[i][chg]);
    }
  }

  fHistQA[0] = new TH1F("fHistVtx", "Z vertex distribution", 100, -15., 15.);
  fHistQA[1] = new TH1F("fHistnTpcCluster", "n TPC Cluster", 100, 0., 200.);
  fHistQA[2] = new TH1F("fHistnTpcClusterF", "n TPC Cluster findable", 100, 0., 200.);
  fHistQA[3] = new TH1F("dcaHistDcaXY1D", "DCA XY", 210, -2.1, 2.1);
  fHistQA[4] = new TH1F("dcaHistDcaZ1D", "DCA Z", 210, -2.1, 2.1);
  fHistQA[5] = new TH1F("fHistChi2Tpc", "Chi2 TPC", 100, 0., 8.);
  fHistQA[6] = new TH1F("fHistpT", "pT distribution",1000,0.,10.0);
  fHistQA[7] = new TH1F("fHistPhi", "Phi distribution" , 100, -TMath::Pi(), TMath::Pi());
  fHistQA[8] = new TH1F("fHistY", "Y distribution" , 100, -2, 2);
 
  fHistQA[9] = new TH1F("fHistEventCuts", "Event Cuts" , 10, 0, 11);
  fHistQA[9]->GetXaxis()->SetBinLabel(1,"All");
  fHistQA[9]->GetXaxis()->SetBinLabel(2,"NoVertex");
  fHistQA[9]->GetXaxis()->SetBinLabel(3,"PileUp");
  fHistQA[9]->GetXaxis()->SetBinLabel(4,"z-vertex>10");
  fHistQA[9]->GetXaxis()->SetBinLabel(5,"eventspions");
  fHistQA[9]->GetXaxis()->SetBinLabel(6,"eventskaons");
  fHistQA[9]->GetXaxis()->SetBinLabel(7,"eventsprotons");
  fHistQA[9]->GetXaxis()->SetBinLabel(8,"eventspions T");
  fHistQA[9]->GetXaxis()->SetBinLabel(9,"eventskaons T");
  fHistQA[9]->GetXaxis()->SetBinLabel(10,"eventsprotons T");


  fHistQA[10] = new TH1F("fHistTrackCuts", "Track Cuts" , 7, 0.5, 7.5);
  fHistQA[10]->GetXaxis()->SetBinLabel(1,"AllTracksInEvents");
  fHistQA[10]->GetXaxis()->SetBinLabel(2,"GetTrack");
  fHistQA[10]->GetXaxis()->SetBinLabel(3,"Filter bit");
  fHistQA[10]->GetXaxis()->SetBinLabel(4,"Y");
  fHistQA[10]->GetXaxis()->SetBinLabel(5,"Pt");
  fHistQA[10]->GetXaxis()->SetBinLabel(6,"DCA");
  fHistQA[10]->GetXaxis()->SetBinLabel(7,"Electron Rejection");

  for(Int_t i = 0; i < 4; i++)  {
   hname = "fHistParticleCounterM";
   hname+= i;
    
  fHistP[i] = new TH1F(hname,Form("fHistParticleCounterM%d",i), 7, 0, 8);
  fHistP[i]->GetXaxis()->SetBinLabel(1,"pions R");
  fHistP[i]->GetXaxis()->SetBinLabel(2,"kaons R");
  fHistP[i]->GetXaxis()->SetBinLabel(3,"protons R");
  fHistP[i]->GetXaxis()->SetBinLabel(4,"pions T");
  fHistP[i]->GetXaxis()->SetBinLabel(5,"kaons T");
  fHistP[i]->GetXaxis()->SetBinLabel(6,"protons T");
  fHistoList->Add(fHistP[i]);
}

  fHistQA2D[0] = new TH2F("dcaHistDcaXY","DCA XY",50, 0, 5,210, -2.1, 2.1);
  fHistQA2D[1] = new TH2F("dcaHistDcaZ","DCA Z", 50, 0, 5, 210, -2.1, 2.1);
  fHistQA2D[2] = new TH2F("fPhiY","Y-Phi",100, -2, 2, 100, -TMath::Pi(), TMath::Pi());



  fHistQALambdas[0] = new TH2F("fHistQALambda", "V0 Details" , 8, 0.5, 8.5,100,0,20);
  fHistQALambdas[1] = new TH2F("fHistQALambdaMinus", "V0 Details" , 8, 0.5, 8.5,100,0,20);
  for(int i=0;i<2;i++){
    fHistQALambdas[i]->GetXaxis()->SetBinLabel(1,"AllV0s");
    fHistQALambdas[i]->GetXaxis()->SetBinLabel(2,"DaugtersNotPrimary");
    fHistQALambdas[i]->GetXaxis()->SetBinLabel(3,"Mother Pos=Neg");
    fHistQALambdas[i]->GetXaxis()->SetBinLabel(4,"MotherFoundInMCarray");
    fHistQALambdas[i]->GetXaxis()->SetBinLabel(5,"PDG of Lambda");
    fHistQALambdas[i]->GetXaxis()->SetBinLabel(6,"IsNotFrom Weak");
    fHistQALambdas[i]->GetXaxis()->SetBinLabel(7,"IsNotFrom Material");
    fHistQALambdas[i]->GetXaxis()->SetBinLabel(8,"Is Primary");
    fHistoList->Add(fHistQALambdas[i]);
  }

  TString originLambdas[]={"PosDaughterPDG","NegDaughterPDG","PosDaugterMotherPDG","NegDaugterMotherPDG","V0PDG"};
  for(int j=0;j<5;j++){
    fOriginLambdas[j][0] = new TH2F("fOriginLambdas"+originLambdas[j],"Origin of particles "+originLambdas[j] , 4000, -4000.0, 4000.0 ,100,0,20);
    fOriginLambdas[j][1] = new TH2F("fOriginLambdas" + originLambdas[j]+"Minus", "Origin of particles "+originLambdas[j]+"Minus" , 400, -4000.0, 4000.0 ,100,0,20);
    fHistoList->Add(fOriginLambdas[j][0]);
    fHistoList->Add(fOriginLambdas[j][1]);
  }

  
  for ( Int_t i = 0; i < 11; i++)
    {
      fHistoList->Add(fHistQA[i]);
      if(i<4) fHistoList->Add(fHistP[i]);
      if(i<3) fHistoList->Add(fHistQA2D[i]);
      
      if(i<5) {
	for(Int_t j = 0 ; j<PARTTYPES; j++)
	  for(int chg=0;chg<2;chg++)
	    {
		fHistoList->Add(fHistQAPID[i][j][chg]);
		fHistoList->Add(fHistQAPIDFail[i][j][chg]);
	    }
      }
    }
    
  for (Int_t i = 0; i < MULTBINS*PARTTYPES; i++){
    for(Int_t chg=0;chg<2;chg++){
      fHistoList->Add(fGeneratedMCPrimaries[i][chg]);
      fHistoList->Add(fMCPrimariesThatAreReconstructed[i][chg]);
      fHistoList->Add(fGeneratedMCPrimaries4D[i][chg]);
      fHistoList->Add(fMCPrimariesThatAreReconstructed4D[i][chg]);
      fHistoList->Add(fMCPrimariesThatAreReconstructedNoNsigma[i][chg]);
      fHistoList->Add(fReconstructedAfterCuts[i][chg]);
      fHistoList->Add(fReconstructedNotPrimaries[i][chg]);
      fHistoList->Add(fReconstructedPrimaries[i][chg]);
      fHistoList->Add(fContamination[i][chg]);


      if(i==4){ //works only for MULTBINS == 1!!!!
	fHistoList->Add(fPrimVsDCA[i][chg]);
	fHistoList->Add(fSecWeakVsDCA[i][chg]);
	fHistoList->Add(fSecMatVsDCA[i][chg]);
	fHistoList->Add(fFakeVsDCA[i][chg]);
	fHistoList->Add(fPrimVsCosPointingAngle[i][chg]);
	fHistoList->Add(fSecWeakVsCosPointingAngle[i][chg]);
	fHistoList->Add(fSecMatVsCosPointingAngle[i][chg]);
	fHistoList->Add(fFakeVsCosPointingAngle[i][chg]);
      
	fHistoList->Add(fPrimVsDecayRadius[i][chg]);
	fHistoList->Add(fSecWeakVsDecayRadius[i][chg]);
	fHistoList->Add(fSecMatVsDecayRadius[i][chg]);
	fHistoList->Add(fFakeVsDecayRadius[i][chg]);
	fHistoList->Add(fAllVsDCA[i][chg]);
	fHistoList->Add(fAllVsCosPointingAngle[i][chg]);
	fHistoList->Add(fAllVsDecayRadius[i][chg]);
      }
    }

  }

  fMassInvLambdaPass = new TH1D("fMassfInvLambdaPass","Mass Assuming Lambda Hypothesis Pass", 10000, 0, 5);
  fHistoList->Add(fMassInvLambdaPass);

  fMassInvAntiLambdaPass = new TH1D("fMassfInvAntiLambdaPass","Mass Assuming AntiLambda Hypothesis Pass", 10000, 0, 5);
  fHistoList->Add(fMassInvAntiLambdaPass);

  fMassInvLambdaFail = new TH1D("fMassfInvLambdaFail","Mass Assuming Lambda Hypothesis Fail", 10000, 0, 5);
  fHistoList->Add(fMassInvLambdaFail);

  fMassInvAntiLambdaFail = new TH1D("fMassfInvAntiLambdaFail","Mass Assuming AntiLambda Hypothesis Fail", 10000, 0, 5);
  fHistoList->Add(fMassInvAntiLambdaFail);

  fYLambda = new TH1D("fYLambda", "|Y| distribution of Lambda", 500, 0.0, 8.);
  fPtLambda = new TH1D("fPtLambda", "Pt distribution of Lambda", 500, 0.0, 8.);
  fYAntiLambda = new TH1D("fYAntiLambda", "|Y| distribution of AntiLambda", 500, 0.0, 8.);
  fPtAntiLambda = new TH1D("fPtAntiLambda", "Pt distribution of AntiLambda", 500, 0.0, 8.);
  fHistoList->Add(fYLambda);
  fHistoList->Add(fPtLambda);
  fHistoList->Add(fYAntiLambda);
  fHistoList->Add(fPtAntiLambda);
  
  fCutsLambda = new TH1D("fCutsLambda","Cuts Lambda", 20, 0.5, 20.5);
  fCutsAntiLambda = new TH1D("fCutsAntiLambda","Cuts AntiLambda", 20, 0.5, 20.5);

  fHistoList->Add(fCutsLambda);
  fHistoList->Add(fCutsAntiLambda);

  fTruePtLambdaMC = new TH2D("fTruePtLambdaMC","True pT of Lambdas MC",10000, -5000, 5000.,50,0.,10.0);
  fTruePtAntiLambdaMC = new TH2D("fTruePtAntiLambdaMC","True pT of AntiLambdas MC",10000, -5000, 5000.,50,0.,10.0);

  fRecPtLambdaMC = new TH2D("fRecPtLambdaMC","Rec pT of Lambdas MC",10000, -5000, 5000.,50,0.,10.0);
  fRecPtAntiLambdaMC = new TH2D("fRecPtAntiLambdaMC","Rec pT of AntiLambdas MC",10000, -5000, 5000.,50,0.,10.0);

  fHistoList->Add(fTruePtLambdaMC);
  fHistoList->Add(fTruePtAntiLambdaMC);
  fHistoList->Add(fRecPtLambdaMC);
  fHistoList->Add(fRecPtAntiLambdaMC);

  
  
//********** PID ****************

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fpidResponse = inputHandler->GetPIDResponse();
  std::cout<<"*******"<< fpidResponse<<std::endl;
  
  // ************************
  PostData(1, fHistoList);
}


//_____________________________________________________________________

bool IsPionNSigmaM3(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime)
{

    if (mom > 0.5) {
        if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 2)
            return true;
	}
    else {
        if (TMath::Abs(nsigmaTPCPi) < 2)
            return true;
    }

  return false;
}

bool IsPionNSigmaV0M3(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{

  if (TMath::Abs(nsigmaTPCPi) < 3.0) return true;

  return false;
}

bool IsPionNSigmaV0TPC5M3(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{

  if (TMath::Abs(nsigmaTPCPi) < 5.0) return true;

  return false;
}


bool IsProtonNSigmaV0TPC5M3(float mom, float nsigmaTPCP, float nsigmaTOFP)
{

  if (TMath::Abs(nsigmaTPCP) < 5.0) return true;

  return false;
}

bool IsPionNSigma3M3(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime)
{

    if (mom > 0.5) {
        if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 3)
            return true;
	}
    else {
        if (TMath::Abs(nsigmaTPCPi) < 3)
            return true;
    }

  return false;
}

bool IsKaonNSigmaM3(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
    if (mom > 0.5) {
      //rejection of unwanted contamination
      if(mom>1 && TOFtime<-400)
	return false;

      if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < 2)
	return true;
	}
    else {
        if (TMath::Abs(nsigmaTPCK) < 2)
            return true;
    }


  return false;
}

bool IsKaonNSigma3M3(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
  if (mom > 0.5) {
    //rejection of unwanted contamination
   // if(mom>1 && TOFtime<-400)
     // return false;
 
    if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < 3)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCK) < 3)
      return true;
  }


  return false;
}

bool IsProtonNSigmaM3(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{

    if (mom > 0.5) {
      if(mom>1.8 && TOFtime<-300)
	return false;



        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 2)
            return true;
	}
    else {
        if (TMath::Abs(nsigmaTPCP) < 2)
            return true;
    }



  return false;
}

bool IsProtonNSigmaV0M3(float mom, float nsigmaTPCP, float nsigmaTOFP)
{

  if (mom < 0.8) {
    if (TMath::Abs(nsigmaTPCP) < 3.0) return true;
  } else {
    if (nsigmaTOFP < -999.) {
      if (TMath::Abs(nsigmaTPCP) < 3.0) return true;
    } else {
      if (TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 3.0) return true;
    }
    }

  return false;
}

bool IsProtonNSigma3M3(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{
  if (mom > 0.5) {
   
   // if(mom>1.8 && TOFtime<-300)
     
     // return false;


    if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 3)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCP) < 3)
      return true;
  }
  

  return false;
}


bool IsElectronM3(float nsigmaTPCe, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
{
  if(TMath::Abs(nsigmaTPCe)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3)
      return true;
   else
     return false;
}

//_______________________________________________________

void AliAnalysisTaskParticleEffDYM3::UserExec(Option_t *)
{


  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliAODEvent *fAOD = aodH->GetEvent();
  fAODpidUtil = aodH->GetAODpidUtil();
  

  /***Get Event****/
  //AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
 /* AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aodEvent) return;
  AliAODHeader *fAODheader = (AliAODHeader*)aodEvent->GetHeader();
  Double_t mult = fAODheader->GetRefMultiplicity();
// AliCentrality* alicent= aodEvent->GetCentrality(); //in PbPb and pPb
//  Double_t mult = alicent->GetCentralityPercentile("V0A"); //in pPb
//  Double_t mult = alicent->GetCentralityPercentile("V0A"); //in PbPb
  fHistEv->Fill(mult); */
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aodEvent) return;
  fHistEvCuts[0]->Fill(0);
  
  Double_t mult;
  
    if(fEstEventMult == kV0M)
    {
      AliCentrality* alicent= aodEvent->GetCentrality(); //in PbPb and pPb
      AliMultSelection *mult_selection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");
      mult = mult_selection->GetMultiplicityPercentile("V0M"); 
    }
    else if(fEstEventMult == kRefMult)
    {
      AliAODHeader *fAODheader = (AliAODHeader*)aodEvent->GetHeader();
      mult = fAODheader->GetRefMultiplicity(); 

   }
   else if(fEstEventMult == kV0A)
    {
      AliCentrality* alicent= aodEvent->GetCentrality(); //in PbPb and pPb
      mult = alicent->GetCentralityPercentile("V0A"); 

    }
  //if(mult < 0.001 || mult > 100000) return;
  //fHistEv[0]->Fill(mult);


  if(fIfAliEventCuts){
    //******* Ali Event Cuts - applied on AOD event - standard cuts for Run2 as prepared by DPG group ************
    if (!fEventCuts->AcceptEvent(aodEvent)) {
      return;
    }
    //******************************************
  }

  //******************
  // load MC array
  // arrayMC =  (TClonesArray*)aodEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  //  if(!arrayMC) {
  //  printf("AliAnalysisTaskParticleEficiency::UserExec: MC particles branch not found!\n");
  //return;
  //}
    
  // load MC header
  //mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  //if(!mcHeader) {
  //printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC header branch not found!\n");
  //return;
  // }
  //*********************   
 


  // EVENT SELECTION ********************

  fHistQA[9]->Fill(1);

  // collision candidate 
  // if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return;


  //****** Multiplicity selection *********
  Int_t fcent = -999;  
  //if(mult >= 0 && mult <=20)  fcent = 0;
  //else if(mult >= 20 && mult <=39) fcent = 1;
  //else if(mult >= 40 && mult <=59) fcent = 2;
  //else if(mult >= 60 && mult <=90) fcent = 3;
  //else if(mult >= 99990 && mult <=99936) fcent = 4;
  //else if(mult >= 999937 && mult <=99944) fcent = 5;
  //else if(mult >= 999945 && mult <=99957) fcent = 6;
  //else if(mult >= 999958 && mult <=99149) fcent = 6;
  //else fcent = 7;
  //if (fcent == 7) return;

fcent=0;
int fcent2=0;
  if(mult >= 0 && mult <=20)  fcent2 = 10;
  else if(mult >= 20 && mult <=40) fcent2 = 1;
  else if(mult >= 40 && mult <=70) fcent2 = 2;
  else if(mult >= 70 && mult <=100) fcent2 = 3;
  if(fcent2!=3) return; // check this
 // else return;

  if(fcent2==10)fHistEvCuts[0]->Fill(1);
  else if(fcent2==1)fHistEvCuts[1]->Fill(1);
  else if(fcent2==2)fHistEvCuts[2]->Fill(1);
  else if(fcent2==3)fHistEvCuts[3]->Fill(1);

  //"ESDs/pass2/AOD049/*AliAOD.root");
  const AliAODVertex* vertex =(AliAODVertex*) aodEvent->GetPrimaryVertex();
  vertex->GetPosition(fV1M3);
  if (!vertex || vertex->GetNContributors()<=0) return;

  fHistQA[9]->Fill(2);
 if(fcent2==10)fHistEvCuts[0]->Fill(2);
  else if(fcent2==1)fHistEvCuts[1]->Fill(2);
  else if(fcent2==2)fHistEvCuts[2]->Fill(2);
  else if(fcent2==3)fHistEvCuts[3]->Fill(2);


//********* Pile-up removal*******************
  //check this: https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsPileup
  AliAnalysisUtils *anaUtil=new AliAnalysisUtils();
    
  Bool_t fpA2013 = kFALSE;
  Bool_t fMVPlp = kFALSE;
  Bool_t fisPileUp = kFALSE;
  Int_t fMinPlpContribMV = 0;
  Int_t fMinPlpContribSPD = 3;

  if(fpA2013)
  if(anaUtil->IsVertexSelected2013pA(aodEvent)==kFALSE) return;
 
  if(fMVPlp) anaUtil->SetUseMVPlpSelection(kTRUE);
  else anaUtil->SetUseMVPlpSelection(kFALSE);
 
  if(fMinPlpContribMV) anaUtil->SetMinPlpContribMV(fMinPlpContribMV);
  if(fMinPlpContribSPD) anaUtil->SetMinPlpContribSPD(fMinPlpContribSPD);

  if(fisPileUp)
  if(anaUtil->IsPileUpEvent(aodEvent)) return;

  delete anaUtil;   

  fHistQA[9]->Fill(3);
  if(fcent2==10)fHistEvCuts[0]->Fill(3);
  else if(fcent2==1)fHistEvCuts[1]->Fill(3);
  else if(fcent2==2)fHistEvCuts[2]->Fill(3);
  else if(fcent2==3)fHistEvCuts[3]->Fill(3);

  //TString vtxTtl = vertex->GetTitle();
  //if (!vtxTtl.Contains("VertexerTracks")) return;
  Float_t zvtx = vertex->GetZ();
  if (TMath::Abs(zvtx) > 10) return;
  fHistQA[0]->Fill(zvtx);
  fHistQA[9]->Fill(4);
  if(fcent2==10)fHistEvCuts[0]->Fill(4);
  else if(fcent2==1)fHistEvCuts[1]->Fill(4);
  else if(fcent2==2)fHistEvCuts[2]->Fill(4);
  else if(fcent2==3)fHistEvCuts[3]->Fill(4);
  
 //**** getting MC array ******
  TClonesArray  *arrayMC;

  arrayMC = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));


  //get stack 
  //AliStack *mcStack = mcEvent->Stack();
  //if (!mcStack) return;
  //***********************


  // old vertex selection 
  /*const AliESDVertex *vertex = esdEvent->GetPrimaryVertex();
    if (vertex->GetNContributors() < 1) return;

    //z-vertex cut
    if (TMath::Abs(vertex->GetZ()) > 10.) return;
  

    const AliESDVertex *vtxSPD = esdEvent->GetPrimaryVertexSPD();*/
  // Double_t zVertex = vtxSPD->GetZ();

  //std::cout << "Event  Z vtx ==========> " << vertex->GetZ() <<std::endl;
  // centrality selection 
  //AliCentrality *centrality = aodEvent->GetHeader()->GetCentralityP();
  //if (centrality->GetQuality() != 0) return;
  //Double_t cent = centrality->GetCentralityPercentileUnchecked("V0M");
  //if(cent < 0 || cent > 100.) return;


//copying pid information for FB 128
  int labels[20000];
  for (int il=0; il<20000; il++) labels[il] = -1;

  // looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
  for (int i=0;i<aodEvent->GetNumberOfTracks();i++) {
  
    const AliAODTrack *aodtrack=(AliAODTrack*)aodEvent->GetTrack(i);
    if (!aodtrack->TestFilterBit(128)) {
      if(aodtrack->GetID() < 0) continue;
      labels[aodtrack->GetID()] = i;
    }
  }



  //RECONSTRUCTED TRACKS 

  TObjArray recoParticleArray[PARTTYPES];

  bool evpass=true;  	
  bool collect[3] = {false,false,false};  	
 
  int hmPionsR=0, hmKaonsR=0, hmProtonsR=0;


  fHistQA[10]->Fill(1,aodEvent->GetNumberOfTracks());
  //loop over AOD tracks 
  for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++) {
  	AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iTracks); 
	if (!track)continue;
	if(track->Y() < -0.5 || track->Y() > 0.5)
      		continue; 
      	
	
	UInt_t filterBit = fFB;
	if(!track->TestFilterBit(filterBit))continue;		
        
	bool isPionNsigma = 0;
	bool isKaonNsigma = 0;
	bool isProtonNsigma  = 0;

    	if (isPionNsigma){
		if (track->Pt() > 0.2 || track->Pt() < 2.5){
		continue;
		}
	}
	if (isKaonNsigma){
     		if (track->Pt() > 0.5 || track->Pt() < 2.5){
     		continue;
     		}
        }
     	if (isProtonNsigma){
     		if (track->Pt() > 0.5 || track->Pt() < 2.5){
     		continue;
     		}
        }
        
	AliAODTrack* aodtrackpid;
		
	if(filterBit==(1 << (7)))
	aodtrackpid =(AliAODTrack*)aodEvent->GetTrack(labels[-1-aodEvent->GetTrack(iTracks)->GetID()]);
	else
	aodtrackpid = track;

	float nSigmaTPCPi = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kPion);
	float nSigmaTPCK = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kKaon);
	float nSigmaTPCP = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kProton);
    
	double nSigmaTOFPi = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kPion);
	double nSigmaTOFK = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kKaon);
	double nSigmaTOFP = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kProton);


	float tTofSig = aodtrackpid->GetTOFsignal();
	double pidTime[5]; aodtrackpid->GetIntegratedTimes(pidTime);

	isPionNsigma = (IsPionNSigmaM3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]) && !IsKaonNSigma3M3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigma3M3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
	isKaonNsigma = (!IsPionNSigma3M3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && IsKaonNSigmaM3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigma3M3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
	isProtonNsigma = (!IsPionNSigma3M3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && !IsKaonNSigma3M3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && IsProtonNSigmaM3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
   
   
   
   	if (isPionNsigma){
	collect[0]=true; 
	if(fcent2==10)fHistP[0]->Fill(1);
	else if(fcent2==1)fHistP[1]->Fill(1);
	else if(fcent2==2)fHistP[2]->Fill(1);
	else if(fcent2==3)fHistP[3]->Fill(1);
	}
  	if (isKaonNsigma){
	collect[1]=true; 
	if(fcent2==10)fHistP[0]->Fill(2);
	else if(fcent2==1)fHistP[1]->Fill(2);
	else if(fcent2==2)fHistP[2]->Fill(2);
	else if(fcent2==3)fHistP[3]->Fill(2);
	}
	if (isProtonNsigma){
	collect[2]=true;
	if(fcent2==10)fHistP[0]->Fill(3);
	else if(fcent2==1)fHistP[1]->Fill(3);
	else if(fcent2==2)fHistP[2]->Fill(3);
	else if(fcent2==3)fHistP[3]->Fill(3);
	}	
  	}



if(collect[0]==true){
	fHistQA[9]->Fill(5);
	if(fcent2==10)fHistEvCuts[0]->Fill(5);
	else if(fcent2==1)fHistEvCuts[1]->Fill(5);
	else if(fcent2==2)fHistEvCuts[2]->Fill(5);
	else if(fcent2==3)fHistEvCuts[3]->Fill(5);
}

if(collect[1]==true){
	fHistQA[9]->Fill(6);
	if(fcent2==10)fHistEvCuts[0]->Fill(6);
	else if(fcent2==1)fHistEvCuts[1]->Fill(6);
	else if(fcent2==2)fHistEvCuts[2]->Fill(6);
	else if(fcent2==3)fHistEvCuts[3]->Fill(6);
}


if(collect[2]==true){
	fHistQA[9]->Fill(7);
	if(fcent2==10)fHistEvCuts[0]->Fill(7);
	else if(fcent2==1)fHistEvCuts[1]->Fill(7);
	else if(fcent2==2)fHistEvCuts[2]->Fill(7);
	else if(fcent2==3)fHistEvCuts[3]->Fill(7);
}




	//if(!evpass) return;
	
  for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++) {
    //get track 
    

    AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iTracks); 
    if (!track)continue;
    fHistQA[10]->Fill(2);

      UInt_t filterBit = fFB;
    if(!track->TestFilterBit(filterBit))continue;		

    Int_t charge = 0;
    if(track->Charge() > 0 ) charge=0;
    else if (track->Charge() < 0 ) charge=1; 
  
    fHistQA[10]->Fill(3);
     
    if(track->Y() < -0.5 || track->Y() > 0.5)
      continue; 
    fHistQA[10]->Fill(4);

    
    //DCA
    
    Double_t DCAXY;
    Double_t DCAZ;

    DCAXY = -TMath::Abs(track->DCA());
    DCAZ = -TMath::Abs(track->ZAtDCA());
 
      if(!(DCAXY==-999 || DCAZ==-999)){
	//if(TMath::Abs(DCAXY) > 0.0182 + 0.035*TMath::Power(track->Pt(), -1.01)) continue; //XY, Pt dep
	//no DCA cut
	//if(TMath::Abs(DCAXY) > 1000.0) {continue;} //XY
	//if(TMath::Abs(DCAZ) > 1000.0) {continue;} //Z
      }
    else {
      // code from Michael and Prabhat from AliAnalysisTaskDptDptCorrelations
      // const AliAODVertex* vertex = (AliAODVertex*) aodEvent->GetPrimaryVertex(); (already defined above)
      float vertexX  = -999.;
      float vertexY  = -999.;
      float vertexZ  = -999.;

      if(vertex) {
	Double32_t fCov[6];
	vertex->GetCovarianceMatrix(fCov);
	if(vertex->GetNContributors() > 0) {
	  if(fCov[5] != 0) {
	    vertexX = vertex->GetX();
	    vertexY = vertex->GetY();
	    vertexZ = vertex->GetZ();

	  }
	}
      }

      Double_t pos[3];
      track->GetXYZ(pos);

      Double_t DCAX = pos[0] - vertexX;
      Double_t DCAY = pos[1] - vertexY;
      DCAZ = pos[2] - vertexZ;
      DCAXY = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));
    }

    fHistQA[10]->Fill(6);

    AliAODTrack* aodtrackpid;

    //for FB 128 - tpc only tracks
    if(filterBit==(1 << (7)))
      aodtrackpid =(AliAODTrack*)aodEvent->GetTrack(labels[-1-aodEvent->GetTrack(iTracks)->GetID()]);
    else
      aodtrackpid = track;

    //Electron rejection
    float nSigmaTPCPi = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kPion);
    float nSigmaTPCK = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kKaon);
    float nSigmaTPCP = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kProton);
    float nSigmaTPCe = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kElectron);
    if(IsElectronM3(nSigmaTPCe,nSigmaTPCPi,nSigmaTPCK,nSigmaTPCP))
      continue;
    fHistQA[10]->Fill(7);
    
    fHistQA[1]->Fill(track->GetTPCClusterInfo(2,1)); 
    //fHistQA[2]->Fill(track->GetTPCNclsF());
    fHistQA[3]->Fill(DCAXY);
    fHistQA[4]->Fill(DCAZ);
    Float_t chi2Tpc = track->Chi2perNDF();
    fHistQA[5]->Fill(chi2Tpc);
    fHistQA[6]->Fill(track->Pt());

    float px=track->Px(); float py=track->Py();  float ph=atan2(py,px); //track->Phi()
    float tPt = track->Pt();

    fHistQA[7]->Fill(ph);
    fHistQA[8]->Fill(track->Y());
    fHistQA2D[2]->Fill(track->Y(),ph);

    fHistQA2D[0]->Fill(tPt,DCAXY);
    fHistQA2D[1]->Fill(tPt,DCAZ);

    //PID monitors
    double nSigmaTOFPi = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kPion);
    double nSigmaTOFK = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kKaon);
    double nSigmaTOFP = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kProton);

    float tdEdx = aodtrackpid->GetTPCsignal();
    float tTofSig = aodtrackpid->GetTOFsignal();
    double pidTime[5]; aodtrackpid->GetIntegratedTimes(pidTime);


    fHistQAPID[0][0][charge]->Fill(tPt,tdEdx);
    fHistQAPID[1][0][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
    fHistQAPID[2][0][charge]->Fill(tPt,nSigmaTOFPi);
    fHistQAPID[3][0][charge]->Fill(tPt,nSigmaTPCPi);
    fHistQAPID[4][0][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);

    fHistQAPIDFail[0][0][charge]->Fill(tPt,tdEdx);
    fHistQAPIDFail[1][0][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
    fHistQAPIDFail[2][0][charge]->Fill(tPt,nSigmaTOFPi);
    fHistQAPIDFail[3][0][charge]->Fill(tPt,nSigmaTPCPi);
    fHistQAPIDFail[4][0][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);

    bool isPionNsigma = 0;
    bool isKaonNsigma = 0;
    bool isProtonNsigma  = 0;

    if(fPidMethod==kNSigma){
    //******** With double counting *******************
      isPionNsigma = (IsPionNSigmaM3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]));
      isKaonNsigma = (IsKaonNSigmaM3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]));
      isProtonNsigma = (IsProtonNSigmaM3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
    }
    else if(fPidMethod==kNSigmaNoDoubleCounting){
      //******** Without double counting *******************
      double nSigmaPIDPi = 0, nSigmaPIDK = 0, nSigmaPIDP = 0;

      if(track->Pt()<0.5){
	nSigmaPIDPi = abs(nSigmaTPCPi);
	nSigmaPIDK  = abs(nSigmaTPCK);
	nSigmaPIDP  = abs(nSigmaTPCP);
      }
      else{
	nSigmaPIDPi = TMath::Hypot(nSigmaTPCPi,nSigmaTOFPi);
	nSigmaPIDK= TMath::Hypot(nSigmaTPCK,nSigmaTOFK);
	nSigmaPIDP= TMath::Hypot(nSigmaTPCP,nSigmaTOFP);
      }

      if(nSigmaPIDPi<nSigmaPIDK && nSigmaPIDPi<nSigmaPIDP){
	isPionNsigma = (IsPionNSigmaM3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]));
      }
      else if(nSigmaPIDK<nSigmaPIDPi && nSigmaPIDK<nSigmaPIDP){
	isKaonNsigma = (IsKaonNSigmaM3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]));
      }
      else if(nSigmaPIDP<nSigmaPIDPi && nSigmaPIDP<nSigmaPIDK){
	isProtonNsigma = (IsProtonNSigmaM3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      }
    }
    else if(fPidMethod==kExclusivePID){
      //******** Exclusive PID ********************
      isPionNsigma = (IsPionNSigmaM3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]) && !IsKaonNSigmaM3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigmaM3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isKaonNsigma = (!IsPionNSigmaM3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && IsKaonNSigmaM3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigmaM3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isProtonNsigma = (!IsPionNSigmaM3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && !IsKaonNSigmaM3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && IsProtonNSigmaM3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
    }
    else if(fPidMethod==kExclusivePIDDiffRejection){
      //******** Exclusive PID, different rejection  ********************
      isPionNsigma = (IsPionNSigmaM3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]) && !IsKaonNSigma3M3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigma3M3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isKaonNsigma = (!IsPionNSigma3M3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && IsKaonNSigmaM3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigma3M3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isProtonNsigma = (!IsPionNSigma3M3(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && !IsKaonNSigma3M3(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && IsProtonNSigmaM3(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
    }
    if (isPionNsigma){
      fHistQAPID[0][1][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][1][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
      fHistQAPID[2][1][charge]->Fill(tPt,nSigmaTOFPi);
      fHistQAPID[3][1][charge]->Fill(tPt,nSigmaTPCPi);
      fHistQAPID[4][1][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);
    }
    else
      {
	fHistQAPIDFail[0][1][charge]->Fill(tPt,tdEdx);
	fHistQAPIDFail[1][1][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
	fHistQAPIDFail[2][1][charge]->Fill(tPt,nSigmaTOFPi);
	fHistQAPIDFail[3][1][charge]->Fill(tPt,nSigmaTPCPi);
	fHistQAPIDFail[4][1][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);
      }
    if (isKaonNsigma){
      fHistQAPID[0][2][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][2][charge]->Fill(tPt,tTofSig-pidTime[3]);//kaon
      fHistQAPID[2][2][charge]->Fill(tPt,nSigmaTOFK);
      fHistQAPID[3][2][charge]->Fill(tPt,nSigmaTPCK);
      fHistQAPID[4][2][charge]->Fill(nSigmaTPCK,nSigmaTOFK);
    }
    else
      {
	fHistQAPIDFail[0][2][charge]->Fill(tPt,tdEdx);
	fHistQAPIDFail[1][2][charge]->Fill(tPt,tTofSig-pidTime[3]);//kaon
	fHistQAPIDFail[2][2][charge]->Fill(tPt,nSigmaTOFK);
	fHistQAPIDFail[3][2][charge]->Fill(tPt,nSigmaTPCK);
	fHistQAPIDFail[4][2][charge]->Fill(nSigmaTPCK,nSigmaTOFK);
      }
    if (isProtonNsigma){
      fHistQAPID[0][3][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][3][charge]->Fill(tPt,tTofSig-pidTime[4]);//proton
      fHistQAPID[2][3][charge]->Fill(tPt,nSigmaTOFP);
      fHistQAPID[3][3][charge]->Fill(tPt,nSigmaTPCP);
      fHistQAPID[4][3][charge]->Fill(nSigmaTPCP,nSigmaTOFP);
    }
    else
      {
	fHistQAPIDFail[0][3][charge]->Fill(tPt,tdEdx);
	fHistQAPIDFail[1][3][charge]->Fill(tPt,tTofSig-pidTime[4]);//proton
	fHistQAPIDFail[2][3][charge]->Fill(tPt,nSigmaTOFP);
	fHistQAPIDFail[3][3][charge]->Fill(tPt,nSigmaTPCP);
	fHistQAPIDFail[4][3][charge]->Fill(nSigmaTPCP,nSigmaTOFP);
      }

	fReconstructedAfterCuts[PARTTYPES*fcent][charge]->Fill(track->Y(), track->Pt());//Fills hist. for all reconstructed particles after cuts
    
    if(!arrayMC){
      continue;
    }
    //get coresponding MC particle 
    Int_t label = TMath::Abs(track->GetLabel());
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(label);

   //getting no. of tracks for each particle species after all the cuts:

    //********* PID - pions ********
     if (isPionNsigma){
	if (track->Pt() > 0.2 || track->Pt() < 2.5)
       fReconstructedAfterCuts[PARTTYPES*fcent+1][charge]->Fill(track->Y(), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[1].Add(MCtrk);
       }
       //Fills for all identified pions found after cuts (reconstructed) - numerator for Efficiency

     //********* PID - kaons ********
     if (isKaonNsigma){
     if (track->Pt() > 0.3 || track->Pt() < 2.5)
       fReconstructedAfterCuts[PARTTYPES*fcent+2][charge]->Fill(track->Y(), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[2].Add(MCtrk);
       }
       //Fills for all identified kaons found after cuts (reconstructed) - numerator for Efficiency

    //********* PID - protons ********
     if (isProtonNsigma){
     if (track->Pt() > 0.5 || track->Pt() < 2.5)
       fReconstructedAfterCuts[PARTTYPES*fcent+3][charge]->Fill(track->Y(), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[3].Add(MCtrk);
       }

      //Fills for all identified protos found after cuts (reconstructed) - numerator for Efficiency
   //******************************

     //get coresponding MC particle 
     // Int_t label = TMath::Abs(track->GetLabel()); //moved up
     // if(!label) std::cout<<"no label"<<std::endl;
     //if(label) std::cout<<"label = "<<label<<std::endl;
       
    //AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(label); //moved up
    if (!MCtrk) continue;
    if(MCtrk->Charge()==0){std::cout<<"!!!"<<std::endl; continue;}
    recoParticleArray[0].Add(MCtrk);


    //Fills histogram for particles that are contamination from secondaries:
    if (!MCtrk->IsPhysicalPrimary()) {
      fReconstructedNotPrimaries[PARTTYPES*fcent][charge]->Fill(track->Y(), track->Pt());
    }
    else{
      fReconstructedPrimaries[PARTTYPES*fcent][charge]->Fill(track->Y(), track->Pt());
    }

 
    int PDGcode = MCtrk->GetPdgCode();

   //And secondaries for different particle species:
    if (!MCtrk->IsPhysicalPrimary() && (isPionNsigma && abs(PDGcode)==211)) { //secondaries in pions
      fReconstructedNotPrimaries[PARTTYPES*fcent+1][charge]->Fill(track->Y(), track->Pt());
    }
    else if(MCtrk->IsPhysicalPrimary() && (isPionNsigma && abs(PDGcode)==211)) {
      fReconstructedPrimaries[PARTTYPES*fcent+1][charge]->Fill(track->Y(), track->Pt());
    }

    if (!MCtrk->IsPhysicalPrimary() && (isKaonNsigma && abs(PDGcode)==321)) { //secondaries in kaons
      fReconstructedNotPrimaries[PARTTYPES*fcent+2][charge]->Fill(track->Y(), track->Pt());
    }
    else if(MCtrk->IsPhysicalPrimary() && (isKaonNsigma && abs(PDGcode)==321)) {
      fReconstructedPrimaries[PARTTYPES*fcent+2][charge]->Fill(track->Y(), track->Pt());
    }

    if (!MCtrk->IsPhysicalPrimary() && (isProtonNsigma && abs(PDGcode)==2212)) { //secondaries in protons
      fReconstructedNotPrimaries[PARTTYPES*fcent+3][charge]->Fill(track->Y(), track->Pt());
    } 
    else if(MCtrk->IsPhysicalPrimary() && (isProtonNsigma && abs(PDGcode)==2212)) {
      fReconstructedPrimaries[PARTTYPES*fcent+3][charge]->Fill(track->Y(), track->Pt());
    } 



    //step 1, TOF Matching
    UInt_t status;
    status=track->GetStatus();
    if((status&AliVTrack::kTOFout)==0 || (status&AliVTrack::kTIME)==0)
      status=0;
    if(track->Pt()<0.5) status = 1;

    //Misidentification fraction
    if(abs(PDGcode)==211)
      {
	if(isPionNsigma)
	  fMisidentification[fcent][charge]-> Fill(1,0.5);
	if(isKaonNsigma)
	  fMisidentification[fcent][charge]-> Fill(1,1.5);
	if(isProtonNsigma)
	  fMisidentification[fcent][charge]-> Fill(1,2.5);
	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma)
	  if(status)
	    fMisidentification[fcent][charge]-> Fill(1,3.5);


      }
    else if(abs(PDGcode)==321)
      {
	if(isPionNsigma)
	  fMisidentification[fcent][charge]-> Fill(2,0.5);
	if(isKaonNsigma)
	  fMisidentification[fcent][charge]-> Fill(2,1.5);
	if(isProtonNsigma)
	  fMisidentification[fcent][charge]-> Fill(2,2.5);
	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma)
	  if(status)
	    fMisidentification[fcent][charge]-> Fill(2,3.5);


      }
    else if(abs(PDGcode) == 2212)
      {
	if(isPionNsigma)
	  fMisidentification[fcent][charge]-> Fill(3,0.5);
	if(isKaonNsigma)
	  fMisidentification[fcent][charge]-> Fill(3,1.5);
	if(isProtonNsigma)
	  {
	  fMisidentification[fcent][charge]-> Fill(3,2.5);
	  }
	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma)
	  if(status)
	    fMisidentification[fcent][charge]-> Fill(3,3.5);
      }


      fContamination[PARTTYPES*fcent][charge]-> Fill(PDGcode,track->Pt()); 
    //Contaminations: "how many pions are in the kaons sample"? etc.
    //Do not use for corrections: using those values will be dependant on i.e. Pi/K ratio in MC
    //Use misidentification fraction instead
    if(isPionNsigma)
      {
	fContamination[PARTTYPES*fcent+1][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for pions
      }
    if(isKaonNsigma)
      {
	fContamination[PARTTYPES*fcent+2][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for kaons
      }
    if(isProtonNsigma)
      {
	fContamination[PARTTYPES*fcent+3][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for protons
      }
  
  }

  //loop over V0s

  for (Int_t i = 0; i < aodEvent->GetNumberOfV0s(); i++) 
    {
      double LambdaMass = 1.115683;
      double K0sMass = 0.497613;
      
      int cutLam = 1;
      int cutALam = 1;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);
      AliAODv0 *aodv0 = aodEvent->GetV0(i);
      if (!aodv0) continue;
      if (aodv0->GetNDaughters() > 2) continue;
      if (aodv0->GetNProngs() > 2) continue;
      if (aodv0->GetCharge() != 0) continue;
      if (aodv0->ChargeProng(0) == aodv0->ChargeProng(1)) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);
      
      if (aodv0->CosPointingAngle(fV1M3) < 0.99) continue;
      //if (aodv0->CosPointingAngle(fV1M3) < 0.95) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);
      
      AliAODTrack *daughterTrackPos = (AliAODTrack *)aodv0->GetDaughter(0); //getting positive daughter track
      AliAODTrack *daughterTrackNeg = (AliAODTrack *)aodv0->GetDaughter(1); //getting negative daughter track
      if (!daughterTrackPos) continue; //daughter tracks must exist
      if (!daughterTrackNeg) continue;
      if (daughterTrackNeg->Charge() == daughterTrackPos->Charge()) continue; //and have different charge
      
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);
      fMassInvLambdaFail->Fill(aodv0->MassLambda());
      fMassInvAntiLambdaFail->Fill(aodv0->MassAntiLambda());



      if(aodv0->Pt() < 0.5 || aodv0->Pt() > 4) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);
      
      if(TMath::Abs(aodv0->Y()) > 0.8) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);
 
      if(aodv0->GetOnFlyStatus() == kTRUE) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);

      if(aodv0->DcaV0Daughters() > 1.0) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);

      if(aodv0->DcaPosToPrimVertex() < 0.06) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);
      
      if(aodv0->DcaNegToPrimVertex() < 0.06) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);
      
      fDCAtoPrimVtx = 0.6; // aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx
      //fCutsLambda->Fill(cutLam++);
      //fCutsAntiLambda->Fill(cutALam++);
      
      if(aodv0->DecayLength(fV1M3) > 60) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);
 
      Double_t radius = aodv0->RadiusV0();
      	if ( radius < 0.5 ) {
      	AliDebugClass(2, "Failed fiducial volume");
      	continue;
      }

           
      int negid = aodv0->GetNegID();
      int posid = aodv0->GetPosID();
      
      AliAODTrack *trackpos = (AliAODTrack*)aodv0->GetDaughter(0);
      AliAODTrack *trackneg = (AliAODTrack*)aodv0->GetDaughter(1);

      if((!trackpos) || (!trackneg)) continue;
      
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);

      if(TMath::Abs(trackpos->Y()) > 0.8) continue;
      if(TMath::Abs(trackneg->Y()) > 0.8) continue;
      if(trackpos->GetTPCNcls() < 70) continue;
      if(trackneg->GetTPCNcls() < 70) continue;
      if(trackpos->Chi2perNDF() > 4.0) continue;
      if(trackneg->Chi2perNDF() > 4.0) continue;
	  
      if(!(trackpos->GetStatus() & (AliESDtrack::kTPCrefit))) continue;
      if(!(trackneg->GetStatus() & (AliESDtrack::kTPCrefit))) continue;
      fCutsLambda->Fill(cutLam++);
      fCutsAntiLambda->Fill(cutALam++);

	  
      //PID
      //double nSigmaTPCPiPos = fpidResponse->NumberOfSigmasTPC(trackpos,AliPID::kPion);
      //double nSigmaTPCPPos = fpidResponse->NumberOfSigmasTPC(trackpos,AliPID::kProton);
      //double nSigmaTPCPiNeg = fpidResponse->NumberOfSigmasTPC(trackneg,AliPID::kPion);
      //double nSigmaTPCPNeg = fpidResponse->NumberOfSigmasTPC(trackneg,AliPID::kProton);

      //double nSigmaTOFPiPos = fpidResponse->NumberOfSigmasTOF(trackpos,AliPID::kPion);
      //double nSigmaTOFPPos = fpidResponse->NumberOfSigmasTOF(trackpos,AliPID::kProton);

      //double nSigmaTOFPiNeg = fpidResponse->NumberOfSigmasTOF(trackneg,AliPID::kPion);
      //double nSigmaTOFPNeg = fpidResponse->NumberOfSigmasTOF(trackneg,AliPID::kProton);

      double nSigmaTPCPiPos = fAODpidUtil->NumberOfSigmasTPC(trackpos,AliPID::kPion);
      double nSigmaTPCPPos = fAODpidUtil->NumberOfSigmasTPC(trackpos,AliPID::kProton);
      double nSigmaTPCPiNeg = fAODpidUtil->NumberOfSigmasTPC(trackneg,AliPID::kPion);
      double nSigmaTPCPNeg = fAODpidUtil->NumberOfSigmasTPC(trackneg,AliPID::kProton);


      //TOF time
      float tdEdxPos = trackpos->GetTPCsignal();
      float tTofSigPos = trackpos->GetTOFsignal();
      double pidTimePos[5]; trackpos->GetIntegratedTimes(pidTimePos);

      float tdEdxNeg = trackneg->GetTPCsignal();
      float tTofSigNeg = trackneg->GetTOFsignal();
      double pidTimeNeg[5]; trackneg->GetIntegratedTimes(pidTimeNeg);

      Float_t probMisPos = 1.0;
      Float_t probMisNeg = 1.0;


      double nSigmaTOFPiPos = 0;
      double nSigmaTOFPPos = 0;
	  
      double nSigmaTOFPiNeg = 0;
      double nSigmaTOFPNeg = 0;
	  
      if(((trackpos->GetStatus() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((trackpos->GetStatus() & AliVTrack::kTIME) == AliVTrack::kTIME))
	{	      
	  probMisPos = fAODpidUtil->GetTOFMismatchProbability(trackpos);
	}
	  
      if(((trackneg->GetStatus() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((trackneg->GetStatus() & AliVTrack::kTIME) == AliVTrack::kTIME))
	{	      
	  probMisNeg = fAODpidUtil->GetTOFMismatchProbability(trackneg);
	}

	  
	  
      if(!(((trackpos->GetStatus() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((trackpos->GetStatus() & AliVTrack::kTIME) == AliVTrack::kTIME)) || probMisPos > 0.01)
	{

	  if(!(((trackneg->GetStatus() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((trackpos->GetStatus() & AliVTrack::kTIME) == AliVTrack::kTIME)) || probMisNeg > 0.01)
	    {
	      nSigmaTOFPiPos = -1000;
	      nSigmaTOFPPos = -1000;
	      nSigmaTOFPiNeg = -1000;
	      nSigmaTOFPNeg = -1000;
		  	      
	      //tFemtoV0->SetTOFProtonTimePos(-1000);
	      //tFemtoV0->SetTOFPionTimePos(-1000);
	      //tFemtoV0->SetTOFKaonTimePos(-1000);
	      //tFemtoV0->SetTOFProtonTimeNeg(-1000);
	      //tFemtoV0->SetTOFPionTimeNeg(-1000);
	      //tFemtoV0->SetTOFKaonTimeNeg(-1000);
	    }
	}
      else
	{
	  if(((trackpos->GetStatus() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((trackpos->GetStatus() & AliVTrack::kTIME) == AliVTrack::kTIME) && probMisPos < 0.01)
	    {

	      nSigmaTOFPiPos = fAODpidUtil->NumberOfSigmasTOF(trackpos, AliPID::kPion);
	      nSigmaTOFPPos = fAODpidUtil->NumberOfSigmasTOF(trackpos, AliPID::kProton);
		  

	    }
	  if(((trackneg->GetStatus() & AliVTrack::kTOFout) == AliVTrack::kTOFout) && ((trackneg->GetStatus() & AliVTrack::kTIME) == AliVTrack::kTIME) && probMisNeg < 0.01)
	    {

	      nSigmaTOFPiNeg = fAODpidUtil->NumberOfSigmasTOF(trackneg, AliPID::kPion);
	      nSigmaTOFPNeg = fAODpidUtil->NumberOfSigmasTOF(trackneg, AliPID::kProton);

	    }
	    

	}

      bool isPionNsigmaPos = 0;
      bool isProtonNsigmaPos  = 0;
      bool isPionNsigmaNeg = 0;
      bool isProtonNsigmaNeg  = 0;

      isPionNsigmaPos = IsPionNSigmaV0TPC5M3(trackpos->Pt(),nSigmaTPCPiPos,nSigmaTOFPiPos);
      isProtonNsigmaPos = IsProtonNSigmaV0TPC5M3(trackpos->Pt(),nSigmaTPCPPos,nSigmaTOFPPos);
      isPionNsigmaNeg = IsPionNSigmaV0TPC5M3(trackneg->Pt(),nSigmaTPCPiNeg,nSigmaTOFPiNeg);
      isProtonNsigmaNeg = IsProtonNSigmaV0TPC5M3(trackneg->Pt(),nSigmaTPCPNeg,nSigmaTOFPNeg);

      bool Lambda = false;
      bool AntiLambda = false;

      //antilambda
      if(isProtonNsigmaNeg && isPionNsigmaPos)
	{
	  fCutsAntiLambda->Fill(cutALam++);  
	  if(trackpos->Pt() < 0.16 || trackpos->Pt() > 4.0) continue; //pions plus
	  fCutsAntiLambda->Fill(cutALam++);
	  if(trackneg->Pt() < 0.3 || trackneg->Pt() > 4.0)  continue; //antiproton
	  fCutsAntiLambda->Fill(cutALam++);
	  if(aodv0->MassAntiLambda() < (LambdaMass-0.0038) || aodv0->MassAntiLambda() > (LambdaMass+0.0038)) continue;
	  if(aodv0->MassK0Short() > 0.48 && aodv0->MassK0Short() < 0.515) continue;
	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	    fCutsAntiLambda->Fill(cutALam++);
	    fMassInvAntiLambdaPass->Fill(aodv0->MassAntiLambda());
	    fPtAntiLambda->Fill(aodv0->Pt());
	    fYAntiLambda->Fill(TMath::Abs(aodv0->Y()));
	    fHistQALambdas[1]->Fill(1,aodv0->Pt());
	    fAllVsCosPointingAngle[PARTTYPES*fcent+4][1]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt());
	    fAllVsDecayRadius[PARTTYPES*fcent+4][1]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	  }
	  fAllVsDCA[PARTTYPES*fcent+4][1]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt());
	  AntiLambda = true;
	 
	}
	  
      //lambda
      if(isProtonNsigmaPos && isPionNsigmaNeg)
	{
	  fCutsLambda->Fill(cutLam++);  
	  if(trackpos->Pt() < 0.3 || trackpos->Pt() > 4.0) continue; //proton
	  fCutsLambda->Fill(cutLam++);
	  if(trackneg->Pt() < 0.16 || trackneg->Pt() > 4.0) continue; //pion minus
	  fCutsLambda->Fill(cutLam++);
	  if(aodv0->MassLambda() < (LambdaMass-0.0038) || aodv0->MassLambda() > (LambdaMass+0.0038)) continue;
	  if(aodv0->MassK0Short() > 0.48 && aodv0->MassK0Short() < 0.515) continue;
	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	    fCutsLambda->Fill(cutLam++);
	    fMassInvLambdaPass->Fill(aodv0->MassLambda());
	    fPtLambda->Fill(aodv0->Pt());
	    fYLambda->Fill(TMath::Abs(aodv0->Y()));
	    fHistQALambdas[0]->Fill(1,aodv0->Pt());
	    fAllVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt());
	    fAllVsDecayRadius[PARTTYPES*fcent+4][0]->Fill(aodv0->RadiusV0(),aodv0->Pt());

	  }
	  fAllVsDCA[PARTTYPES*fcent+4][0]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt());
	  Lambda = true;
	}


      if(!arrayMC) continue;
	      
      //get coresponding MC particles
      Int_t labelPos = TMath::Abs(trackpos->GetLabel());
      Int_t labelNeg = TMath::Abs(trackneg->GetLabel());
      AliAODMCParticle *MCtrkPos = (AliAODMCParticle*)arrayMC->At(labelPos);
      AliAODMCParticle *MCtrkNeg = (AliAODMCParticle*)arrayMC->At(labelNeg);
      
      Int_t motherPos = MCtrkPos->GetMother();
      Int_t motherNeg = MCtrkNeg->GetMother();

      if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	//********* PID - lambdas, antilambda ********
	if(Lambda)
	  {
	    int charge = 0;
	    fReconstructedAfterCuts[PARTTYPES*fcent+4][charge]->Fill(aodv0->Y(), aodv0->Pt());
	    fOriginLambdas[0][0]->Fill(MCtrkPos->GetPdgCode(), aodv0->Pt());
	    fOriginLambdas[1][0]->Fill(MCtrkNeg->GetPdgCode(), aodv0->Pt());
	  }
	if(AntiLambda)
	  {
	    int charge = 1;
	    fReconstructedAfterCuts[PARTTYPES*fcent+4][charge]->Fill(aodv0->Y(), aodv0->Pt());
	    fOriginLambdas[0][1]->Fill(MCtrkPos->GetPdgCode(), aodv0->Pt());
	    fOriginLambdas[1][1]->Fill(MCtrkNeg->GetPdgCode(), aodv0->Pt());
	  }
	//********************************
      }

      if(MCtrkPos->IsPhysicalPrimary() || MCtrkNeg->IsPhysicalPrimary())
	continue;

      if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	if(Lambda) {
	  fHistQALambdas[0]->Fill(2,aodv0->Pt());
	  fOriginLambdas[2][0]->Fill(((AliAODMCParticle*)arrayMC->At(motherPos))->GetPdgCode(), aodv0->Pt());
	  fOriginLambdas[3][0]->Fill(((AliAODMCParticle*)arrayMC->At(motherNeg))->GetPdgCode(), aodv0->Pt());
	}
	if(AntiLambda) {
	  fHistQALambdas[1]->Fill(2,aodv0->Pt());
	  fOriginLambdas[2][1]->Fill(((AliAODMCParticle*)arrayMC->At(motherPos))->GetPdgCode(), aodv0->Pt());
	  fOriginLambdas[3][1]->Fill(((AliAODMCParticle*)arrayMC->At(motherNeg))->GetPdgCode(), aodv0->Pt());
	}
      }
     
      if(motherPos != motherNeg) continue;
 
      if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	if(Lambda) fHistQALambdas[0]->Fill(3,aodv0->Pt());
	if(AntiLambda) fHistQALambdas[1]->Fill(3,aodv0->Pt());
      }

      AliAODMCParticle *MCtrkMother = (AliAODMCParticle*)arrayMC->At(motherPos);
      if(!MCtrkMother) continue;

      if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	if(Lambda) fHistQALambdas[0]->Fill(4,aodv0->Pt());
	if(AntiLambda) fHistQALambdas[1]->Fill(4,aodv0->Pt());
      }

      int pdgMother = MCtrkMother->GetPdgCode();


      if (MCtrkMother->IsPhysicalPrimary() && Lambda && pdgMother==3122)
	{
	  fPrimVsDCA[PARTTYPES*fcent+4][0]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt()); 
	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt()); 
	    fPrimVsDecayRadius[PARTTYPES*fcent+4][0]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	  }
	}
      else if (MCtrkMother->IsPhysicalPrimary() && AntiLambda && pdgMother==-3122)
	{
	  fPrimVsDCA[PARTTYPES*fcent+4][1]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt()); 
	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+4][1]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt()); 
	    fPrimVsDecayRadius[PARTTYPES*fcent+4][1]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	  }
	}
      else if(MCtrkMother->IsSecondaryFromWeakDecay() && Lambda && pdgMother==3122)
	{
	  fSecWeakVsDCA[PARTTYPES*fcent+4][0]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt()); 
	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt()); 
	    fSecWeakVsDecayRadius[PARTTYPES*fcent+4][0]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	  }
	}
      else if(MCtrkMother->IsSecondaryFromWeakDecay() && AntiLambda && pdgMother==-3122)
	{
	  fSecWeakVsDCA[PARTTYPES*fcent+4][1]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt()); 
	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+4][1]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt()); 
	    fSecWeakVsDecayRadius[PARTTYPES*fcent+4][1]->Fill(aodv0->RadiusV0(),aodv0->Pt()); 
	  }
	}
      else if(MCtrkMother->IsSecondaryFromMaterial() && Lambda && pdgMother==3122)
	{
	fSecMatVsDCA[PARTTYPES*fcent+4][0]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt());
	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt()); 
	    fSecMatVsDecayRadius[PARTTYPES*fcent+4][0]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	  }
	}
      else if(MCtrkMother->IsSecondaryFromMaterial() && AntiLambda && pdgMother==-3122)
	{
	fSecMatVsDCA[PARTTYPES*fcent+4][1]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt());
	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){ 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+4][1]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt()); 
	    fSecMatVsDecayRadius[PARTTYPES*fcent+4][1]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	  }
	}
      else if(Lambda)
	{

	  fFakeVsDCA[PARTTYPES*fcent+4][0]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt()); 

	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	    fFakeVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt()); 
	    fFakeVsDecayRadius[PARTTYPES*fcent+4][0]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	  }
	
	}
      else if(AntiLambda)
	{
 
	  fFakeVsDCA[PARTTYPES*fcent+4][1]->Fill(aodv0->DcaV0ToPrimVertex(),aodv0->Pt()); 
	  if( aodv0->DcaV0ToPrimVertex() < fDCAtoPrimVtx){
	    fFakeVsCosPointingAngle[PARTTYPES*fcent+4][1]->Fill(aodv0->CosPointingAngle(fV1M3),aodv0->Pt()); 
	    fFakeVsDecayRadius[PARTTYPES*fcent+4][1]->Fill(aodv0->RadiusV0(),aodv0->Pt());
	  }	  
	}
      if( aodv0->DcaV0ToPrimVertex() > fDCAtoPrimVtx) continue; // we do not longer need full DCA V0 to prim vertex sample

      //contamination from secondaries
      if (!MCtrkMother->IsPhysicalPrimary() && Lambda && pdgMother==3122) { //secondaries in lambdas
	fReconstructedNotPrimaries[PARTTYPES*fcent+4][0]->Fill(aodv0->Y(), aodv0->Pt());
      } 
      else if(MCtrkMother->IsPhysicalPrimary() && Lambda && pdgMother==3122) {
	fReconstructedPrimaries[PARTTYPES*fcent+4][0]->Fill(aodv0->Y(), aodv0->Pt());
      } 
      else if (!MCtrkMother->IsPhysicalPrimary() && AntiLambda && pdgMother==-3122) { //secondaries in lambdas
	fReconstructedNotPrimaries[PARTTYPES*fcent+4][1]->Fill(aodv0->Y(), aodv0->Pt());
      } 
      else if(MCtrkMother->IsPhysicalPrimary() && AntiLambda && pdgMother==-3122) {
	fReconstructedPrimaries[PARTTYPES*fcent+4][1]->Fill(aodv0->Y(), aodv0->Pt());
      }
      


      if(Lambda)
	{
	  fTruePtLambdaMC->Fill(pdgMother,MCtrkMother->Pt());
	  fRecPtLambdaMC->Fill(pdgMother,aodv0->Pt());
	  recoParticleArray[4].Add(MCtrkMother);
	  if(pdgMother==3122)  {fHistQALambdas[0]->Fill(5,aodv0->Pt());} 
	  else continue;
	  if(!MCtrkMother->IsSecondaryFromWeakDecay()) {fHistQALambdas[0]->Fill(6,aodv0->Pt());  } 
	  else continue;
	  if(!MCtrkMother->IsSecondaryFromMaterial()) {fHistQALambdas[0]->Fill(7,aodv0->Pt());  }
	  else continue;
	  if(MCtrkMother->IsPhysicalPrimary()) {fHistQALambdas[0]->Fill(8,aodv0->Pt());  }

	  fOriginLambdas[4][0]->Fill(pdgMother, aodv0->Pt());
	}
      if(AntiLambda)
	{
	  fTruePtAntiLambdaMC->Fill(pdgMother,MCtrkMother->Pt());
	  fRecPtAntiLambdaMC->Fill(pdgMother,aodv0->Pt());
	  recoParticleArray[4].Add(MCtrkMother);
	  if(pdgMother==-3122)  {fHistQALambdas[1]->Fill(5,aodv0->Pt());  } 
	  else continue;
	  if(!MCtrkMother->IsSecondaryFromWeakDecay()) {fHistQALambdas[1]->Fill(6,aodv0->Pt());}
	  else continue;
	  if(!MCtrkMother->IsSecondaryFromMaterial()) {fHistQALambdas[1]->Fill(7,aodv0->Pt()); }
	  else continue;
	  if(MCtrkMother->IsPhysicalPrimary()) {fHistQALambdas[1]->Fill(8,aodv0->Pt()); }
	  fOriginLambdas[4][1]->Fill(pdgMother, aodv0->Pt());
	}


     

    }


  // MONTECARLO PARTICLES 
  if(!arrayMC){
    AliError("Array of MC particles not found");
    return;
  }

  bool evpassT=true;  	
  bool collectT[3] = {false,false,false};  	
  	
  // loop over MC stack 
  int hmPionsT=0, hmKaonsT=0, hmProtonsT=0;
  for (Int_t ipart = 0; ipart < arrayMC->GetEntries(); ipart++) {
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(ipart);

    if (!MCtrk) continue;
    
	if(MCtrk->Y() < -0.5 || MCtrk->Y() > 0.5){
	continue; }
	
      if(MCtrk->GetPdgCode() == 211){
      	if (MCtrk->Pt() < 0.2 || MCtrk->Pt() > 2.5){
	continue;
	}
      }
      if(MCtrk->GetPdgCode() == 321){
      	if (MCtrk->Pt() < 0.5 || MCtrk->Pt() > 2.5){
	continue;
	}
       }
       if(MCtrk->GetPdgCode() == 3122){
      	if (MCtrk->Pt() < 0.5 || MCtrk->Pt() > 2.5){
	continue;
	}
       }
      
      // check physical primary 

      if(MCtrk->IsPhysicalPrimary()) // Not from weak decay!
	{
    
    Int_t PDGcode = TMath::Abs(MCtrk->GetPdgCode()); 

   	if (PDGcode==211){
	collectT[0]=true; 
	if(fcent2==10)fHistP[0]->Fill(4);
	else if(fcent2==1)fHistP[1]->Fill(4);
	else if(fcent2==2)fHistP[2]->Fill(4);
	else if(fcent2==3)fHistP[3]->Fill(4);
	}
  	if (PDGcode==321){
	collectT[1]=true; 
	if(fcent2==10)fHistP[0]->Fill(5);
	else if(fcent2==1)fHistP[1]->Fill(5);
	else if(fcent2==2)fHistP[2]->Fill(5);
	else if(fcent2==3)fHistP[3]->Fill(5);
	}
	if (PDGcode==2212){
	collectT[2]=true;
	if(fcent2==10)fHistP[0]->Fill(6);
	else if(fcent2==1)fHistP[1]->Fill(6);
	else if(fcent2==2)fHistP[2]->Fill(6);
	else if(fcent2==3)fHistP[3]->Fill(6);
	}
	}
   }
  	

        if(collectT[0]==true){
        	fHistQA[9]->Fill(8);
		if(fcent2==10)fHistEvCuts[0]->Fill(8);
		else if(fcent2==1)fHistEvCuts[1]->Fill(8);
		else if(fcent2==2)fHistEvCuts[2]->Fill(8);
		else if(fcent2==3)fHistEvCuts[3]->Fill(8);
        }

	if(collectT[1]==true){
		fHistQA[9]->Fill(9);
		if(fcent2==10)fHistEvCuts[0]->Fill(9);
		else if(fcent2==1)fHistEvCuts[1]->Fill(9);
		else if(fcent2==2)fHistEvCuts[2]->Fill(9);
		else if(fcent2==3)fHistEvCuts[3]->Fill(9);
	}


	if(collectT[2]==true){
		fHistQA[9]->Fill(10);
		if(fcent2==10)fHistEvCuts[0]->Fill(10);
		else if(fcent2==1)fHistEvCuts[1]->Fill(10);
		else if(fcent2==2)fHistEvCuts[2]->Fill(10);
		else if(fcent2==3)fHistEvCuts[3]->Fill(10);

	}
  
  
  
  
      for (Int_t ipart = 0; ipart < arrayMC->GetEntries(); ipart++) {
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(ipart);

    if (!MCtrk) continue;
    //std::cout<<"particle obtained"<<std::endl;
    
    Int_t PDGcode = TMath::Abs(MCtrk->GetPdgCode()); 

    
    //if(MCtrk->Charge() == 0) continue;
    Int_t charge=0;

    if(MCtrk->Charge() < 0) charge=1;
    else if(MCtrk->Charge() > 0) charge=0;


    if(MCtrk->Charge() == 0)
      {	
	if(MCtrk->GetPdgCode() == 3122) charge = 0;
	else if(MCtrk->GetPdgCode() == -3122) charge = 1;
      }


     
    //*** PID - check if pion ***
    //if(PDGcode!=211) continue; //(PDGcode==11 || PDGcode==321 || PDGcode==2212 || PDGcode==13)

      if(MCtrk->Y() < -0.5 || MCtrk->Y() > 0.5){
	continue; }
	
      if(MCtrk->GetPdgCode() == 211){
      	if (MCtrk->Pt() < 0.2 || MCtrk->Pt() > 2.5){
	continue;
	}
      }
      if(MCtrk->GetPdgCode() == 321){
      	if (MCtrk->Pt() < 0.5 || MCtrk->Pt() > 2.5){
	continue;
	}
       }
       if(MCtrk->GetPdgCode() == 3122){
      	if (MCtrk->Pt() < 0.5 || MCtrk->Pt() > 2.5){
	continue;
	}
       }


      
      // check physical primary 

      if(MCtrk->IsPhysicalPrimary()) // Not from weak decay!
	{

	// Filling histograms for MC truth particles
	fGeneratedMCPrimaries[fcent*PARTTYPES][charge]->Fill(MCtrk->Y(), MCtrk->Pt());

	 Double_t val[] = {MCtrk->Y(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	 fGeneratedMCPrimaries4D[fcent*PARTTYPES][charge]->Fill(val);

	 if(PDGcode==211){
	  fGeneratedMCPrimaries[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	  fGeneratedMCPrimaries4D[fcent*PARTTYPES+1][charge]->Fill(val);}
	 else if(PDGcode==321){
	  fGeneratedMCPrimaries[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	  fGeneratedMCPrimaries4D[fcent*PARTTYPES+2][charge]->Fill(val);}
	 else if(PDGcode==2212){
	  fGeneratedMCPrimaries[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	  fGeneratedMCPrimaries4D[fcent*PARTTYPES+3][charge]->Fill(val);}
	 else if(PDGcode==3122){
	  fGeneratedMCPrimaries[fcent*PARTTYPES+4][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	  fGeneratedMCPrimaries4D[fcent*PARTTYPES+4][charge]->Fill(val);}

	  //Filling data from MC truth particles only for particles that were reconstruced
	if (recoParticleArray[0].Contains(MCtrk)){ //All
	  fMCPrimariesThatAreReconstructed[fcent*PARTTYPES][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	  Double_t val[] = {MCtrk->Y(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	  fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES][charge]->Fill(val);
	  
	  fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	  if(PDGcode==211)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	  if(PDGcode==321)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	  if(PDGcode==2212)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	  if(PDGcode==3212)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+4][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	}
	if (recoParticleArray[1].Contains(MCtrk)){ //Pions
	
	  if(PDGcode==211){
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	    Double_t val[] = {MCtrk->Y(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+1][charge]->Fill(val);
	  }
	}
	if (recoParticleArray[2].Contains(MCtrk)){ //Kaons
	  if(PDGcode==321){
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	    Double_t val[] = {MCtrk->Y(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+2][charge]->Fill(val);
	  }
	}
	if (recoParticleArray[3].Contains(MCtrk)){ //Protons
	  if(PDGcode==2212){
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	    Double_t val[] = {MCtrk->Y(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+3][charge]->Fill(val);
	  }
	}
	if (recoParticleArray[4].Contains(MCtrk)){ //Lambdas
	  if(PDGcode==3122){
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+4][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
	    Double_t val[] = {MCtrk->Y(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+4][charge]->Fill(val);
	  }
	}

      }

  }
  PostData(1, fHistoList);
}
//-----------------------------------------------------------------

//void AliAnalysisTaskParticleEff::Terminate(Option_t *) 
//{}


