#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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
#include "AliAODpidUtil.h"
#include "AliAODHeader.h"
#include "AliFemtoEventReaderAODMultSelection.h"
#include "AliMultSelection.h"

#include "AliAnalysisTaskParticleEffpbpb5020.h"


ClassImp(AliAnalysisTaskParticleEffpbpb5020)


double fVer1[3];

//_______________________________________________________

AliAnalysisTaskParticleEffpbpb5020::AliAnalysisTaskParticleEffpbpb5020(const Char_t *partName) :
  AliAnalysisTaskSE(partName), centrality(0), fHistoList(0),  fHistEv(0), fpidResponse(0), fAODpidUtil(0)
{
  for(Int_t i = 0; i < MULTBINS*PARTTYPES; i++)  {
    for(Int_t chg=0;chg<2;chg++){
      //part 1
      fGeneratedMCPrimaries[i][chg] = NULL;
      fMCPrimariesThatAreReconstructed[i][chg] = NULL;
      fMCPrimariesThatAreReconstructedNoNsigma[i][chg] = NULL;
      fReconstructedAfterCuts[i][chg] = NULL;
      fReconstructedNotPrimaries[i][chg] = NULL;
      fReconstructedPrimaries[i][chg] = NULL;
      fContamination[i][chg] = NULL;
      //part 2
      fPrimVsCosPointingAngle[i][chg] = NULL;
      fSecWeakVsCosPointingAngle[i][chg] = NULL;
      fSecMatVsCosPointingAngle[i][chg] = NULL;
      fFakeVsCosPointingAngle[i][chg] = NULL;
      //part 3
      fPrim_DCAxy_Pt[i][chg] = NULL;
      fSecMat_DCAxy_Pt[i][chg] = NULL;
      fSecWeak_DCAxy_Pt[i][chg] = NULL;
      fPrim_DCAz_Pt[i][chg] = NULL;
      fSecMat_DCAz_Pt[i][chg] = NULL;
      fSecWeak_DCAz_Pt[i][chg] = NULL;
      //part 4
      fMCPrimariesThatAreReconstructed4D[i][chg] = NULL;
      fGeneratedMCPrimaries4D[i][chg] = NULL;
    }
  }
  for ( Int_t i = 0; i < 11; i++) { 
    fHistQA[i] = NULL;
    if(i<3) fHistQA2D[i] = NULL;
  }

  DefineOutput(1, TList::Class());
}

//_______________________________________________________

AliAnalysisTaskParticleEffpbpb5020::~AliAnalysisTaskParticleEffpbpb5020()
{
  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis)
    delete fHistoList;
}

//_______________________________________________________

void AliAnalysisTaskParticleEffpbpb5020::UserCreateOutputObjects()
{
  std::cout<<"Create Output Objects"<<std::endl;

  fHistoList = new TList();
  fHistoList->SetOwner(kTRUE);
   
  TString hna, hname[20], htit, htitle[20];
  
  TString parttypename = "None";

  for(Int_t j = 0; j < PARTTYPES; j++)  {

    if (j==0) parttypename="All";
    else if (j==1) parttypename="Pion";
    else if (j==2) parttypename="Kaon";
    else if (j==3) parttypename="Proton";
    else if (j==4) parttypename="Deuteron";
    
    for(Int_t i = 0; i < MULTBINS; i++)  {
      //---------------part 1----------------------
      hname[0] = "hGeneratedMCPrimariesEffM"; hname[0]+=i; hname[0]+=parttypename;
      htitle[0] = "Kinematic level eta_pT (prim only) M"; htitle[0]+=i; htitle[0]+=parttypename;
      fGeneratedMCPrimaries[i*PARTTYPES+j][0] = new TH2F(hname[0].Data(),htitle[0].Data(),50, -1.5, 1.5,500,0.,5.0);
      hname[0]+="Minus";htitle[0]+="Minus";
      fGeneratedMCPrimaries[i*PARTTYPES+j][1] = new TH2F(hname[0].Data(),htitle[0].Data(),50, -1.5, 1.5,500,0.,5.0);

      hname[1]  = "hMCPrimariesThatAreReconstructedM"; hname[1]+=i; hname[1]+=parttypename;
      htitle[1] = "Reconstructed level eta_pT (prim only) M"; htitle[1]+=i; htitle[1]+=parttypename;
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][0] = new TH2F(hname[1].Data(),htitle[1].Data(),50, -1.5, 1.5,500,0.,5.0);
      hname[1]+="Minus";htitle[1]+="Minus";
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][1] = new TH2F(hname[1].Data(),htitle[1].Data(),50, -1.5, 1.5,500,0.,5.0);

      hname[2]  = "hMCPrimariesThatAreReconstructedNoNsigmaM"; hname[2]+=i; hname[2]+=parttypename;
      htitle[2] = "Reconstructed level eta_pT (prim only) no Nsigma cut only PDG M"; htitle[2]+=i; htitle[2]+=parttypename;
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][0] = new TH2F(hname[2].Data(),htitle[2].Data(),50, -1.5, 1.5,500,0.,5.0);
      hname[2]+="Minus";htitle[2]+="Minus";
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][1] = new TH2F(hname[2].Data(),htitle[2].Data(),50, -1.5, 1.5,500,0.,5.0);

      hname[3] = "hHistoReconstructedAfterCutsM"; hname[3]+=i; hname[3]+=parttypename;
      htitle[3] = "Total Reconstructed tracks M "; htitle[3]+=i; htitle[3]+=parttypename;
      fReconstructedAfterCuts[i*PARTTYPES+j][0] = new TH2F(hname[3].Data(),htitle[3].Data(),50, -1.5, 1.5,500,0.,5.0);
      hname[3]+="Minus";htitle[3]+="Minus";
      fReconstructedAfterCuts[i*PARTTYPES+j][1] = new TH2F(hname[3].Data(),htitle[3].Data(),50, -1.5, 1.5,500,0.,5.0);

      hname[4]  = "hHistoReconstructedNotPrimariesM"; hname[4]+=i; hname[4]+=parttypename;
      htitle[4] = "Reconstructed level eta_pT (not primaries) M"; htitle[4]+=i; htitle[4]+=parttypename;
      fReconstructedNotPrimaries[i*PARTTYPES+j][0] = new TH2F(hname[4].Data(),htitle[4].Data(),50, -1.5, 1.5,500,0.,5.0);
      hname[4]+="Minus";htitle[4]+="Minus";
      fReconstructedNotPrimaries[i*PARTTYPES+j][1] = new TH2F(hname[4].Data(),htitle[4].Data(),50, -1.5, 1.5,500,0.,5.0);

      hname[5]  = "hHistoReconstructedPrimariesM"; hname[5]+=i; hname[5]+=parttypename;
      htitle[5] = "Reconstructed level eta_pT (primaries) M"; htitle[5]+=i; htitle[5]+=parttypename;
      fReconstructedPrimaries[i*PARTTYPES+j][0] = new TH2F(hname[5].Data(),htitle[5].Data(),50, -1.5, 1.5,500,0.,5.0);
      hname[5]+="Minus";htitle[5]+="Minus";
      fReconstructedPrimaries[i*PARTTYPES+j][1] = new TH2F(hname[5].Data(),htitle[5].Data(),50, -1.5, 1.5,500,0.,5.0);

      hname[6]  = "hContaminationM"; hname[6]+=i; hname[6]+=parttypename;
      htitle[6] = "Contamination M"; htitle[6]+=i; htitle[6]+=parttypename;
      fContamination[i*PARTTYPES+j][0] = new TH2F(hname[6].Data(),htitle[6].Data(),6000,-3000 ,3000,500,0.,5.0); //50
      hname[6]+="Minus";htitle[6]+="Minus";
      fContamination[i*PARTTYPES+j][1] = new TH2F(hname[6].Data(),htitle[6].Data(),6000, -3000, 3000,500,0.,5.0); //50

      
      //------------------part 2---------------------------
      hname[7] = "hPrimVsCosPointingAngle"; hname[7]+=i; hname[7]+=parttypename;
      htitle[7] = "Primaries vs CosPointingAngle M"; htitle[7]+=i; htitle[7]+=parttypename;      
      fPrimVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname[7].Data(),htitle[7].Data(),200,-1.0,1.0,50,0.0,4.0);//0.95->-1.0
      hname[7]+="Minus"; htitle[7]+="Minus";
      fPrimVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname[7].Data(),htitle[7].Data(),200,-1.0,1.0,50,0.0,4.0);
      
      hname[8] = "hSecWeakVsCosPointingAngleM"; hname[8]+=i; hname[8]+=parttypename;
      htitle[8] = "Sec. weak decay vs CosPointingAngle M"; htitle[8]+=i; htitle[8]+=parttypename;  
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname[8].Data(),htitle[8].Data(),200,-1.0,1.0,50,0.0,4.0);
      hname[8]+="Minus"; htitle[8]+="Minus";
      fSecWeakVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname[8].Data(),htitle[8].Data(),200,-1.0,1.0,50,0.0,4.0);
      
      hname[9] = "hSecMatVsCosPointingAngleM"; hname[9]+=i; hname[9]+=parttypename;
      htitle[9] = "Sec. material vs CosPointingAngle M"; htitle[9]+=i; htitle[9]+=parttypename;
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname[9].Data(),htitle[9].Data(),200,-1.0,1.0,50,0.0,4.0);
      hname[9]+="Minus"; htitle[9]+="Minus";
      fSecMatVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname[9].Data(),htitle[9].Data(),200,-1.0,1.0,50,0.0,4.0);
      
      hname[10] = "hFakeVsCosPointingAngleM"; hname[10]+=i; hname[10]+=parttypename;
      htitle[10] = "Fake vs CosPointingAngle M"; htitle[10]+=i; htitle[10]+=parttypename;
      fFakeVsCosPointingAngle[i*PARTTYPES+j][0] = new TH2F(hname[10].Data(),htitle[10].Data(),200,-1.0,1.0,7,0.0,4.0);
      hname[10]+="Minus"; htitle[10]+="Minus";
      fFakeVsCosPointingAngle[i*PARTTYPES+j][1] = new TH2F(hname[10].Data(),htitle[10].Data(),200,-1.0,1.0,7,0.0,4.0);
      
      
      //---------------part 3-------------------------------
      hname[11] = "Prim_DCAxy_Pt"; hname[11]+=i; hname[11]+=parttypename;
      htitle[11] = "Prim_DCAxy_Pt M"; htitle[11]+=i; htitle[11]+=parttypename;     
      fPrim_DCAxy_Pt[i*PARTTYPES+j][0] = new TH2D(Form("Prim_DCAxy_Pt%s",hname[11].Data()),Form("Prim_DCAxy_Pt%s",hname[11].Data()),480,-2.4,2.4,50,0,4.0);
      hname[11]+="Minus"; htitle[11]+="Minus";
      fPrim_DCAxy_Pt[i*PARTTYPES+j][1] = new TH2D(Form("Prim_DCAxy_Pt1%s",htitle[11].Data()),Form("Prim_DCAxy_Pt1%s",htitle[11].Data()),480,-2.4,2.4,50,0,4.0);
      
      hname[12] = "Prim_DCAz_Pt"; hname[12]+=i; hname[12]+=parttypename;
      htitle[12] = "Prim_DCAz_Pt M"; htitle[12]+=i; htitle[12]+=parttypename;   
      fPrim_DCAz_Pt[i*PARTTYPES+j][0] = new TH2D(Form("Prim_DCAz_Pt%s",hname[12].Data()),Form("Prim_DCAz_Pt%s",htitle[12].Data()),320,-3.2,3.2,50,0,4.0);
      hname[12]+="Minus"; htitle[12]+="Minus";
      fPrim_DCAz_Pt[i*PARTTYPES+j][1] = new TH2D(Form("Prim_DCAz_Pt1%s",hname[12].Data()),Form("Prim_DCAz_Pt1%s",htitle[12].Data()),320,-3.2,3.2,50,0,4.0);

      hname[13] = "secweak_DCAxy_Pt M"; hname[13]+=i; hname[13]+=parttypename;      
      htitle[13] = "SecWeak_DCAxy_Pt"; htitle[13]+=i; htitle[13]+=parttypename;
      fSecWeak_DCAxy_Pt[i*PARTTYPES+j][0] = new TH2D(Form("SecWeak_DCAxy_Pt%s",hname[13].Data()),Form("SecWeak_DCAxy_Pt%s",htitle[13].Data()),480,-2.4,2.4,50,0,4.0);
      hname[13]+="Minus"; htitle[13]+="Minus";
      fSecWeak_DCAxy_Pt[i*PARTTYPES+j][1] = new TH2D(Form("SecWeak_DCAxy_Pt1%s",hname[13].Data()),Form("SecWeak_DCAxy_Pt1%s",htitle[13].Data()),480,-2.4,2.4,50,0,4.0);
      
      hname[14] = "secweak_DCAz_Pt M"; hname[14]+=i; hname[14]+=parttypename;      
      htitle[14] = "SecWeak_DCAz_Pt"; htitle[14]+=i; htitle[14]+=parttypename;
      fSecWeak_DCAz_Pt[i*PARTTYPES+j][0] = new TH2D(Form("SecWeak_DCAz_Pt%s",hname[14].Data()),Form("SecWeak_DCAz_Pt%s",htitle[14].Data()),320,-3.2,3.2,50,0,4.0);
      hname[14]+="Minus"; htitle[14]+="Minus";
      fSecWeak_DCAz_Pt[i*PARTTYPES+j][1] = new TH2D(Form("SecWeak_DCAz_Pt1%s",hname[14].Data()),Form("SecWeak_DCAz_Pt1%s",htitle[14].Data()),320,-3.2,3.2,50,0,4.0);

      hname[15] = "secmat_DCAxy_Pt M"; hname[15]+=i; hname[15]+=parttypename;      
      htitle[15] = "SecMat_DCAxy_Pt"; htitle[15]+=i; htitle[15]+=parttypename;
      fSecMat_DCAxy_Pt[i*PARTTYPES+j][0] = new TH2D(Form("SecMat_DCAxy_Pt%s",hname[15].Data()),Form("SecMat_DCAxy_Pt%s",htitle[15].Data()),480,-2.4,2.4,50,0,4.0);
      hname[15]+="Minus"; htitle[15]+="Minus";
      fSecMat_DCAxy_Pt[i*PARTTYPES+j][1] = new TH2D(Form("SecMat_DCAxy_Pt1%s",hname[15].Data()),Form("SecMat_DCAxy_Pt1%s",htitle[15].Data()),480,-2.4,2.4,50,0,4.0);
      
      hname[16] = "secmat_DCAz_Pt M"; hname[16]+=i; hname[16]+=parttypename;      
      htitle[16] = "SecMat_DCAz_Pt"; htitle[16]+=i; htitle[16]+=parttypename;
      fSecMat_DCAz_Pt[i*PARTTYPES+j][0] = new TH2D(Form("SecMat_DCAz_Pt%s",hname[16].Data()),Form("SecMat_DCAz_Pt%s",htitle[16].Data()),320,-3.2,3.2,50,0,4.0);
      hname[16]+="Minus"; htitle[16]+="Minus";
      fSecMat_DCAz_Pt[i*PARTTYPES+j][1] = new TH2D(Form("SecMat_DCAz_Pt1%s",hname[16].Data()),Form("SecMat_DCAz_Pt1%s",htitle[16].Data()),320,-3.2,3.2,50,0,4.0);


      //---------------part 4--------------------------------
      Double_t min[]={-1.0,0,-10,0};
      Double_t max[]={1.0,10,10,2*TMath::Pi()};
      Int_t nbins[]={10,100,10,10};
      
      hname[17]  = "hGeneratedMCPrimariesEff4DM"; hname[17]+=i; hname[17]+=parttypename;
      htitle[17] = "Kinematic level eta_pT (prim only) M"; htitle[17]+=i; htitle[17]+=parttypename;
      fGeneratedMCPrimaries4D[i*PARTTYPES+j][0] = new THnSparseF(hname[17].Data(),htitle[17].Data(),4,nbins,min,max);
      hname[17]+="Minus";htitle[17]+="Minus";
      fGeneratedMCPrimaries4D[i*PARTTYPES+j][1] = new THnSparseF(hname[17].Data(),htitle[17].Data(),4,nbins,min,max);

      hname[18]  = "hMCPrimariesThatAreReconstructed4DM"; hname[18]+=i; hname[18]+=parttypename;
      htitle[18] = "Reconstructed level eta_pT (prim only) M"; htitle[18]+=i; htitle[18]+=parttypename;
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][0] = new THnSparseF(hname[18].Data(),htitle[18].Data(),4,nbins,min,max);
      hname[18]+="Minus";htitle[18]+="Minus";
      fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][1] = new THnSparseF(hname[18].Data(),htitle[18].Data(),4,nbins,min,max);
      


      for(int iter = 0;iter < 2; iter++){
         //part 1
         fReconstructedAfterCuts[i*PARTTYPES+j][iter]->Sumw2();
         fReconstructedNotPrimaries[i*PARTTYPES+j][iter]->Sumw2();
         fReconstructedPrimaries[i*PARTTYPES+j][iter]->Sumw2();
         fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][iter]->Sumw2();
         fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][iter]->Sumw2();
         fGeneratedMCPrimaries[i*PARTTYPES+j][iter]->Sumw2();
         fContamination[i*PARTTYPES+j][iter]->Sumw2();
         //part 3
         fPrim_DCAxy_Pt[i*PARTTYPES+j][iter]->Sumw2();
         fSecMat_DCAxy_Pt[i*PARTTYPES+j][iter]->Sumw2();
         fSecWeak_DCAxy_Pt[i*PARTTYPES+j][iter]->Sumw2();
         fPrim_DCAz_Pt[i*PARTTYPES+j][iter]->Sumw2();
         fSecMat_DCAz_Pt[i*PARTTYPES+j][iter]->Sumw2();
         fSecWeak_DCAz_Pt[i*PARTTYPES+j][iter]->Sumw2();
         //part 2

         fPrimVsCosPointingAngle[i*PARTTYPES+j][iter]->Sumw2();
         fSecWeakVsCosPointingAngle[i*PARTTYPES+j][iter]->Sumw2();
         fSecMatVsCosPointingAngle[i*PARTTYPES+j][iter]->Sumw2();
         fFakeVsCosPointingAngle[i*PARTTYPES+j][iter]->Sumw2();    

         //part 4
         fGeneratedMCPrimaries4D[i*PARTTYPES+j][iter]->Sumw2();
         fMCPrimariesThatAreReconstructed4D[i*PARTTYPES+j][iter]->Sumw2(); 


      }

    }

    hna = "pidTPCdEdx";  hna+=parttypename;
    htit = parttypename + " TPC dEdx vs. momentum";
    fHistQAPID[0][j][0] = new TH2F(hna, htit, 100, 0.0, 5.0, 250, 0.0, 500.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPID[0][j][1] = new TH2F(hna, htit, 100, 0.0, 5.0, 250, 0.0, 500.0);
    hna = "pidTOFTime";  hna+=parttypename;
    htit = parttypename + " TOF Time vs. momentum";
    fHistQAPID[1][j][0] = new TH2F(hna, htit, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPID[1][j][1] = new TH2F(hna, htit, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    hna = "pidTOFNSigmaNew";  hna+=parttypename;
    htit = parttypename + " TOF NSigmaNew vs. momentum";
    fHistQAPID[2][j][0]= new TH2F(hna,htit, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPID[2][j][1]= new TH2F(hna,htit, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hna = "pidTPCNSigmaNew";  hna+=parttypename;
    htit = parttypename + " TPC NSigmaNew vs. momentum";
    fHistQAPID[3][j][0] = new TH2F(hna,htit, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPID[3][j][1] = new TH2F(hna,htit, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hna = "pidTPCTOFNSigmaNew";  hna+=parttypename;
    htit = parttypename + " TPC vs TOF NSigmaNew";
    fHistQAPID[4][j][0] = new TH2F(hna,htit, 200, -10.0, 10.0, 200, -10.0, 10.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPID[4][j][1] = new TH2F(hna,htit, 200, -10.0, 10.0, 200, -10.0, 10.0);
    hna = "pidTPCdEdxFail";  hna+=parttypename;
    htit = parttypename + " TPC dEdx vs. momentum Fail";
    fHistQAPIDFail[0][j][0] = new TH2F(hna, htit, 100, 0.0, 5.0, 250, 0.0, 500.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPIDFail[0][j][1] = new TH2F(hna, htit, 100, 0.0, 5.0, 250, 0.0, 500.0);
    hna = "pidTOFTimeFail";  hna+=parttypename;
    htit = parttypename + " TOF Time vs. momentum Fail";
    fHistQAPIDFail[1][j][0] = new TH2F(hna, htit, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPIDFail[1][j][1] = new TH2F(hna, htit, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    hna = "pidTOFNSigmaNewFail";  hna+=parttypename;
    htit = parttypename + " TOF NSigmaNew vs. momentum Fail";
    fHistQAPIDFail[2][j][0]= new TH2F(hna,htit, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPIDFail[2][j][1]= new TH2F(hna,htit, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hna = "pidTPCNSigmaNewFail";  hna+=parttypename;
    htit = parttypename + " TPC NSigmaNew vs. momentum Fail";
    fHistQAPIDFail[3][j][0] = new TH2F(hna,htit, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPIDFail[3][j][1] = new TH2F(hna,htit, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hna = "pidTPCTOFNSigmaNewFail";  hna+=parttypename;
    htit = parttypename + " TPC vs TOF NSigmaNew Fail";
    fHistQAPIDFail[4][j][0] = new TH2F(hna,htit, 200, -10.0, 10.0, 200, -10.0, 10.0);
    htit+="Minus"; hna+="Minus";
    fHistQAPIDFail[4][j][1] = new TH2F(hna,htit, 200, -10.0, 10.0, 200, -10.0, 10.0);
  }
    
  fHistEv = new TH1F("fHistEv", "Multiplicity", 100, 0, 100);
  fHistoList->Add(fHistEv);
    std::cout<<"Dupaaaa 3"<<std::endl;
  for(Int_t i = 0; i < MULTBINS; i++)  {
    hna= "fHistEventCutsM";
    hna+= i;
    
    fHistEvCuts[i] = new TH1F(hna,Form("Event Cuts M%d",i) , 4, 0, 5);
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(1,"All");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(2,"NoVertex");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(3,"PileUp");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(4,"z-vertex>10");
    fHistoList->Add(fHistEvCuts[i]);

    for(Int_t chg=0;chg<2;chg++){
      hna = "hMisidentificationM"; hna+=i; if(chg==0) hna+="Plus"; else hna+="Minus"; 
      htit = "Misidentification Fraction M"; htit+=i; if(chg==0) htit+="Plus"; else htit+="Minus";
      fMisidentification[i][chg] = new TH3F(hna.Data(),htit.Data(), 4, 0.5, 3.5, 4 , 0, 5, 60, 0.1, 2.5);
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(1,"Pions, MC");
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(2,"Kaons, MC");
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(3,"Protons, MC");
      fMisidentification[i][chg]->GetXaxis()->SetBinLabel(4,"Deuterons, MC");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(1,"Pions, Data");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(2,"Kaons, Data");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(3,"Protons, Data");
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(4,"Deuterons, Data");      
      fMisidentification[i][chg]->GetYaxis()->SetBinLabel(5,"Other, Data");
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
  fHistQA[8] = new TH1F("fHistEta", "Eta distribution" , 100, -2, 2);
 
  fHistQA[9] = new TH1F("fHistEventCuts", "Event Cuts" , 4, 0, 5);
  fHistQA[9]->GetXaxis()->SetBinLabel(1,"All");
  fHistQA[9]->GetXaxis()->SetBinLabel(2,"NoVertex");
  fHistQA[9]->GetXaxis()->SetBinLabel(3,"PileUp");
  fHistQA[9]->GetXaxis()->SetBinLabel(4,"z-vertex>10");


  fHistQA[10] = new TH1F("fHistTrackCuts", "Track Cuts" , 7, 0.5, 7.5);
  fHistQA[10]->GetXaxis()->SetBinLabel(1,"AllTracksInEvents");
  fHistQA[10]->GetXaxis()->SetBinLabel(2,"GetTrack");
  fHistQA[10]->GetXaxis()->SetBinLabel(3,"Filter bit");
  fHistQA[10]->GetXaxis()->SetBinLabel(4,"Eta");
  fHistQA[10]->GetXaxis()->SetBinLabel(5,"Pt");
  fHistQA[10]->GetXaxis()->SetBinLabel(6,"DCA");
  fHistQA[10]->GetXaxis()->SetBinLabel(7,"Electron Rejection");

  fHistQA2D[0] = new TH2F("dcaHistDcaXY","DCA XY",50, 0, 5,210, -2.1, 2.1);
  fHistQA2D[1] = new TH2F("dcaHistDcaZ","DCA Z", 50, 0, 5, 210, -2.1, 2.1);
  fHistQA2D[2] = new TH2F("fPhiEta","Eta-Phi",100, -2, 2, 100, -TMath::Pi(), TMath::Pi());

  for ( Int_t i = 0; i < 11; i++)
    {
      fHistoList->Add(fHistQA[i]);
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

      fHistoList->Add(fPrim_DCAxy_Pt[i][chg]);
      fHistoList->Add(fSecWeak_DCAxy_Pt[i][chg]);
      fHistoList->Add(fSecMat_DCAxy_Pt[i][chg]);
      fHistoList->Add(fPrim_DCAz_Pt[i][chg]);
      fHistoList->Add(fSecWeak_DCAz_Pt[i][chg]);
      fHistoList->Add(fSecMat_DCAz_Pt[i][chg]);
      
      //if(i==4){ //works only for MULTBINS == 1!!!!
	fHistoList->Add(fPrimVsCosPointingAngle[i][chg]);
	fHistoList->Add(fSecWeakVsCosPointingAngle[i][chg]);
	fHistoList->Add(fSecMatVsCosPointingAngle[i][chg]);
	fHistoList->Add(fFakeVsCosPointingAngle[i][chg]);
      
     // }
    }

  }


  
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  fAODpidUtil = aodH->GetAODpidUtil();


  PostData(1, fHistoList);
}
    
    
     
bool IsPionNSigmaNew(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime)
{

    if (mom > 0.5 && mom<1.5) {
        if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 3*1.41)
            return true;
	}
    else  {
        if (TMath::Abs(nsigmaTPCPi) < 3)
            return true;
    }

  return false;
}
    

bool IsKaonNSigmaNew(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
  if (mom > 1) {
    if ((TMath::Abs(nsigmaTOFK) < 1) && (TMath::Abs(nsigmaTPCK) < 3))
      return true;
  }
  else if (mom <1 && mom>0.8) {
    if ((TMath::Abs(nsigmaTOFK) < 1.5) && (TMath::Abs(nsigmaTPCK) < 3))
      return true;
  }
  else if (mom <0.8 && mom>0.5) {
    if ((TMath::Abs(nsigmaTOFK) < 2) && (TMath::Abs(nsigmaTPCK) < 3))
      return true;
  }  
  else if (mom <0.5 && mom>0.45) {
    if ((TMath::Abs(nsigmaTOFK) < 2) && (TMath::Abs(nsigmaTPCK) < 3))
      return true;
  }
  else if (mom <0.45 && mom>0.4) {
    if (TMath::Abs(nsigmaTPCK) < 1)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCK) < 2)
      return true;
  }


  return false;
}

bool IsProtonNSigmaNew(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{
  if (mom > 0.5 && mom< 3.0) {
    if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 3)
      return true;
  }

  return false;
}


bool IsDeuteronNSigmaNew(float mom, float nsigmaTPCD, float nsigmaTOFD)
{
     if((mom >= 0.8) && (mom < 1.3)){
       if(nsigmaTPCD < 2)
          return true;
     }
     else if(mom>=1.3 && mom<2.2){
       if((nsigmaTPCD < 2) && (nsigmaTOFD < 2))    
          return true; 
     }

  return false;
}


bool IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP, float nsigmaTPCD)
{
  if(TMath::Abs(nsigmaTPCE)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3 &&TMath::Abs(nsigmaTPCD)>3)
      return true;
   else
     return false;
}


bool IsDeuteronTPCdEdx(float mom, float dEdx, float maxmom)
{

  double a1 = -250.0,  b1 = 400.0;
  double a2 = -135.0,  b2 = 270.0;
  double a3 = -80,   b3 = 190.0;
  double a4 = 0.0,   b4 = 20.0;

  double a5 = 125.0,   b5 = -100.0; 

  if (mom < 1.1) {
    if (dEdx < a1*mom+b1) return false;
  }
  else if (mom < 1.4) {
    if (dEdx < a2*mom+b2) return false;
  }
  else if (mom < 2) {
    if (dEdx < a3*mom+b3) return false;
  }
  else if (mom >= 2) {
    if (dEdx < a4*mom+b4) return false;
  }

  if (mom > maxmom) return false;
  
  return true;
}

//_______________________________________________________

void AliAnalysisTaskParticleEffpbpb5020::UserExec(Option_t *)
{

  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliAODEvent *fAOD = aodH->GetEvent();
  fAODpidUtil = aodH->GetAODpidUtil();
  

  /***Get Event****/
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aodEvent) return;
  AliAODHeader *fAODheader = (AliAODHeader*)aodEvent->GetHeader();
  AliCentrality* alicent= aodEvent->GetCentrality(); //in PbPb and pPb
  AliMultSelection *MultSelection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");

  Double_t mult = 0.0;
  if(!MultSelection) {
    cout << "AliMultSelection object not found!" << endl;
    return;
  }

  else mult = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
  if(mult < 0 || mult >= 80) return;

  cout<<"########fcentrality#######   "<<mult<<endl;
  fHistEv->Fill(mult);



  // EVENT SELECTION ********************                                                                                                                                                                  

  fHistQA[9]->Fill(1);


  //****** Multiplicity selection *********                                                                                                                                                                
  Int_t fcent = -999;

  if (mult >= 0 && mult < 5) fcent = 0;
  else if (mult >= 5 && mult < 10) fcent = 1;
  else if (mult >= 10 && mult < 20) fcent = 2;
  else if (mult >= 20 && mult < 30) fcent = 3;
  //else if (mult >= 30 && mult < 40) fcent = 4;
  //else if (mult >= 40 && mult < 50) fcent = 5;
  //else if (mult >= 50 && mult < 60) fcent = 6;
  //else if (mult >= 60 && mult < 70) fcent = 7;
  //else if (mult >= 70 && mult < 80) fcent = 8;
  else return;


  // EVENT SELECTION ********************
  fHistQA[9]->Fill(1);

  fHistEvCuts[0]->Fill(1);

  const AliAODVertex* vertex =(AliAODVertex*) aodEvent->GetPrimaryVertex();
  vertex->GetPosition(fVer1);
  if (!vertex || vertex->GetNContributors()<=0) return;

  fHistQA[9]->Fill(2);
  fHistEvCuts[0]->Fill(2);
 
  AliAnalysisUtils *anaUtil=new AliAnalysisUtils();
    
  Bool_t fpA2013 = kFALSE;
  Bool_t fMVPlp = kFALSE;
  Bool_t fisPileUp = kTRUE;
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
  fHistEvCuts[0]->Fill(3);
  Float_t zvtx = vertex->GetZ();
  if (TMath::Abs(zvtx) > 7) return;
  fHistQA[0]->Fill(zvtx);
  fHistQA[9]->Fill(4);
  fHistEvCuts[0]->Fill(4);

  //**** getting MC array ******
  TClonesArray  *arrayMC;
  arrayMC = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));

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
  fHistQA[10]->Fill(1,aodEvent->GetNumberOfTracks());
  
      cout<<"  "<<endl<<endl;
  //loop over AOD tracks
  for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++) {
    //get track 
    Int_t label;
    AliAODMCParticle *MCtrk;
    int PDGcode;
    //AliESDtrack* track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(esdEvent),iTracks);
    AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iTracks); 
    if (!track)continue;
    fHistQA[10]->Fill(2);

    UInt_t filterBit = (1 << (0));
  //  UInt_t filterBit = (1 << (4));
  //  UInt_t filterBit = 1;//(1 << 0);
    if(!track->TestFilterBit(filterBit))continue;	

    //charge
    Int_t charge = 0;
    if(track->Charge() > 0 ) charge=0;
    else if (track->Charge() < 0 ) charge=1; 

    fHistQA[10]->Fill(3);
    if(track->Eta() < -0.8 || track->Eta() > 0.8)
      continue; 
      
    fHistQA[10]->Fill(4);
    if (track->Pt() < 0.1 || track->Pt() > 2.5)
      continue;
    fHistQA[10]->Fill(5);

    //single track cuts
     if(track->GetTPCNcls() < 70) continue;

    //DCA
    float vertexX  = -999.;
    float vertexY  = -999.;
    float vertexZ  = -999.;

    if(vertex) {
       Double32_t fCov[6] = {0.0};
       vertex->GetCovarianceMatrix(fCov);
       if(fCov[5] != 0) {
         vertexX = vertex->GetX();
         vertexY = vertex->GetY();
         vertexZ = vertex->GetZ();
       }
   }

   Double_t pos[3];
   track->GetXYZ(pos);

   Double_t DCAX = pos[0] - vertexX;
   Double_t DCAY = pos[1] - vertexY;
   Double_t DCAZ = pos[2] - vertexZ;
   Double_t DCAXY = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));

   label = TMath::Abs(track->GetLabel());
   MCtrk = (AliAODMCParticle*)arrayMC->At(label);
   PDGcode = MCtrk->GetPdgCode();

   if( (DCAZ > 0.3) || (DCAZ < -0.3) || (DCAXY > 0.3) || (DCAXY < -0.3) ) continue;

   fHistQA[10]->Fill(6);
   AliAODTrack* aodtrackpid;
   //if(filterBit==(1 << (7)))
   //for FB 128 - tpc only tracks
    if(filterBit==(1 << (7)))
      aodtrackpid =(AliAODTrack*)aodEvent->GetTrack(labels[-1-aodEvent->GetTrack(iTracks)->GetID()]);
    else
      aodtrackpid = track;
      
    //Electron rejection
    double nSigmaTPCPi = fAODpidUtil->NumberOfSigmasTPC(aodtrackpid,AliPID::kPion);
    double nSigmaTPCK = fAODpidUtil->NumberOfSigmasTPC(aodtrackpid,AliPID::kKaon);
    double nSigmaTPCP = fAODpidUtil->NumberOfSigmasTPC(aodtrackpid,AliPID::kProton);
    double nSigmaTPCe = fAODpidUtil->NumberOfSigmasTPC(aodtrackpid,AliPID::kElectron);
    double nSigmaTPCD = fAODpidUtil->NumberOfSigmasTPC(aodtrackpid,AliPID::kDeuteron);
    
    if(IsElectron(nSigmaTPCe,nSigmaTPCPi,nSigmaTPCK,nSigmaTPCP,nSigmaTPCD))
      continue;
      
    fHistQA[10]->Fill(7);     
    fHistQA[1]->Fill(track->GetTPCClusterInfo(2,1)); 
    fHistQA[3]->Fill(DCAXY);
    fHistQA[4]->Fill(DCAZ);
    Float_t chi2Tpc = track->Chi2perNDF();
    fHistQA[5]->Fill(chi2Tpc);
    fHistQA[6]->Fill(track->Pt());

    float px=track->Px(); float py=track->Py();  float ph=atan2(py,px); //track->Phi()
    float tPt = track->Pt();
    float tP = track->P();
    fHistQA[7]->Fill(ph);
    fHistQA[8]->Fill(track->Eta());
    fHistQA2D[2]->Fill(track->Eta(),ph);
    fHistQA2D[0]->Fill(tPt,DCAXY);
    fHistQA2D[1]->Fill(tPt,DCAZ);

    //PID monitors
    double nSigmaTOFPi = -1000;
    double nSigmaTOFK = -1000;
    double nSigmaTOFP = -1000;
    double nSigmaTOFD = -1000;
     
    ULong_t status = aodtrackpid->GetStatus();
    Float_t probMis;
    if (((status & AliVTrack::kTOFout) == AliVTrack::kTOFout)
	&& ((status & AliVTrack::kTIME) == AliVTrack::kTIME)) {
      
      probMis = fAODpidUtil->GetTOFMismatchProbability(aodtrackpid);
      if(probMis < 0.01){
	nSigmaTOFPi = fAODpidUtil->NumberOfSigmasTOF(aodtrackpid,AliPID::kPion);
	nSigmaTOFK = fAODpidUtil->NumberOfSigmasTOF(aodtrackpid,AliPID::kKaon);
	nSigmaTOFP = fAODpidUtil->NumberOfSigmasTOF(aodtrackpid,AliPID::kProton);
        nSigmaTOFD = fAODpidUtil->NumberOfSigmasTOF(aodtrackpid,AliPID::kDeuteron);
      }
    }
    
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
    bool isDeuteronNsigma  = 0;

    isPionNsigma = (IsPionNSigmaNew(track->P(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]) && !IsKaonNSigmaNew(track->P(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigmaNew(track->P(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]) && !IsDeuteronNSigmaNew(track->P(), nSigmaTPCD,nSigmaTOFD));
    isKaonNsigma = (!IsPionNSigmaNew(track->P(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && IsKaonNSigmaNew(track->P(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigmaNew(track->P(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4])  && !IsDeuteronNSigmaNew(track->P(), nSigmaTPCD,nSigmaTOFD));
    isProtonNsigma = (!IsPionNSigmaNew(track->P(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && !IsKaonNSigmaNew(track->P(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && IsProtonNSigmaNew(track->P(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4])  && !IsDeuteronNSigmaNew(track->P(), nSigmaTPCD,nSigmaTOFD));
   // isDeuteronNsigma = (track->Pt() < 2.2 && !IsPionNSigmaNew(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && !IsKaonNSigmaNewReal(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigmaNew(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4])  && IsDeuteronNSigmaNew(track->Pt(),track->P(),nSigmaTPCD,nSigmaTOFD,1.5));
    isDeuteronNsigma = (track->Pt() < 2.0 && IsDeuteronNSigmaNew(track->P(), nSigmaTPCD, nSigmaTOFD));
     if(isDeuteronNsigma)
      if(track->P() < 2.2)
    	if(!(IsDeuteronTPCdEdx(track->P(), tdEdx, 2.2)))
    	  isDeuteronNsigma = false;


    if (isPionNsigma){
      fHistQAPID[0][1][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][1][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
      fHistQAPID[2][1][charge]->Fill(tPt,nSigmaTOFPi);
      fHistQAPID[3][1][charge]->Fill(tPt,nSigmaTPCPi);
      fHistQAPID[4][1][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);
    }
    else {
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
    else {
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
    else {
	fHistQAPIDFail[0][3][charge]->Fill(tPt,tdEdx);
	fHistQAPIDFail[1][3][charge]->Fill(tPt,tTofSig-pidTime[4]);//proton
	fHistQAPIDFail[2][3][charge]->Fill(tPt,nSigmaTOFP);
	fHistQAPIDFail[3][3][charge]->Fill(tPt,nSigmaTPCP);
	fHistQAPIDFail[4][3][charge]->Fill(nSigmaTPCP,nSigmaTOFP);
      }
    if (isDeuteronNsigma){
      fHistQAPID[0][4][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][4][charge]->Fill(tPt,tTofSig);//deuteron
      fHistQAPID[2][4][charge]->Fill(tPt,nSigmaTOFD);
      fHistQAPID[3][4][charge]->Fill(tP,nSigmaTPCD);
      fHistQAPID[4][4][charge]->Fill(nSigmaTPCD,nSigmaTOFD);
    }
    else {
	fHistQAPIDFail[0][4][charge]->Fill(tPt,tdEdx);
	fHistQAPIDFail[1][4][charge]->Fill(tPt,tTofSig);//deuteron
	fHistQAPIDFail[2][4][charge]->Fill(tPt,nSigmaTOFD);
	fHistQAPIDFail[3][4][charge]->Fill(tP,nSigmaTPCD);
	fHistQAPIDFail[4][4][charge]->Fill(nSigmaTPCD,nSigmaTOFD);
      }
    //int fcent=0;
    fReconstructedAfterCuts[PARTTYPES*fcent][charge]->Fill(track->Eta(), track->Pt());//Fills hist. for all reconstructed particles after cuts
 
      
    if(!arrayMC){
      continue;
    }
   
     if (MCtrk->IsPhysicalPrimary() && (isProtonNsigma && PDGcode==2212))
	{
	    fPrim_DCAxy_Pt[PARTTYPES*fcent+3][0]->Fill(DCAXY,track->Pt());
	    fPrim_DCAz_Pt[PARTTYPES*fcent+3][0]->Fill(DCAZ,track->Pt());  
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+3][0]->Fill(track->Eta(),track->Pt()); 
 	    fPrimVsCosPointingAngle[PARTTYPES*fcent+3][0]->GetXaxis()->SetRangeUser(-1,1);
	}
      else if (MCtrk->IsPhysicalPrimary() && (isProtonNsigma && PDGcode==-2212))
	{
	    fPrim_DCAxy_Pt[PARTTYPES*fcent+3][1]->Fill(DCAXY,track->Pt());
	    fPrim_DCAz_Pt[PARTTYPES*fcent+3][1]->Fill(DCAZ,track->Pt());  
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+3][1]->Fill(track->Eta(),track->Pt()); 
 	    fPrimVsCosPointingAngle[PARTTYPES*fcent+3][1]->GetXaxis()->SetRangeUser(-1,1);

	}
      else if(MCtrk->IsSecondaryFromWeakDecay() && (isProtonNsigma && PDGcode==2212) )
	{
	    fSecWeak_DCAxy_Pt[PARTTYPES*fcent+3][0]->Fill(DCAXY,track->Pt());
	    fSecWeak_DCAz_Pt[PARTTYPES*fcent+3][0]->Fill(DCAZ,track->Pt());    
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+3][0]->Fill(track->Eta(),track->Pt()); 
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+3][0]->GetXaxis()->SetRangeUser(-1,1);
	}
      else if(MCtrk->IsSecondaryFromWeakDecay() && (isProtonNsigma && PDGcode==-2212) )
	{
	    fSecWeak_DCAxy_Pt[PARTTYPES*fcent+3][1]->Fill(DCAXY,track->Pt());
	    fSecWeak_DCAz_Pt[PARTTYPES*fcent+3][1]->Fill(DCAZ,track->Pt());    
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+3][1]->Fill(track->Eta(),track->Pt()); 
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+3][1]->GetXaxis()->SetRangeUser(-1,1);
	 
	}
      else if(MCtrk->IsSecondaryFromMaterial() && (isProtonNsigma && PDGcode==2212) )
	{

	    fSecMat_DCAxy_Pt[PARTTYPES*fcent+3][0]->Fill(DCAXY,track->Pt());
	    fSecMat_DCAz_Pt[PARTTYPES*fcent+3][0]->Fill(DCAZ,track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+3][0]->Fill(track->Eta(),track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+3][0]->GetXaxis()->SetRangeUser(-1,1);

	}
      else if(MCtrk->IsSecondaryFromMaterial() && (isProtonNsigma && PDGcode==-2212) )
	{
	    fSecMat_DCAxy_Pt[PARTTYPES*fcent+3][1]->Fill(DCAXY,track->Pt());
	    fSecMat_DCAz_Pt[PARTTYPES*fcent+3][1]->Fill(DCAZ,track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+3][1]->Fill(track->Eta(),track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+3][1]->GetXaxis()->SetRangeUser(-1,1);
	}         
   
   
      if (MCtrk->IsPhysicalPrimary() && (isKaonNsigma && PDGcode==321))
	{
	    fPrim_DCAxy_Pt[PARTTYPES*fcent+2][0]->Fill(DCAXY,track->Pt());
	    fPrim_DCAz_Pt[PARTTYPES*fcent+2][0]->Fill(DCAZ,track->Pt());  
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+3][0]->Fill(track->Eta(),track->Pt()); 
 	    fPrimVsCosPointingAngle[PARTTYPES*fcent+3][0]->GetXaxis()->SetRangeUser(-1,1);
	}
      else if (MCtrk->IsPhysicalPrimary() && (isKaonNsigma && PDGcode==-321))
	{
	    fPrim_DCAxy_Pt[PARTTYPES*fcent+2][1]->Fill(DCAXY,track->Pt());
	    fPrim_DCAz_Pt[PARTTYPES*fcent+2][1]->Fill(DCAZ,track->Pt());  
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+2][1]->Fill(track->Eta(),track->Pt()); 
 	    fPrimVsCosPointingAngle[PARTTYPES*fcent+2][1]->GetXaxis()->SetRangeUser(-1,1);

	}
      else if(MCtrk->IsSecondaryFromWeakDecay() && (isKaonNsigma && PDGcode==321) )
	{
	    fSecWeak_DCAxy_Pt[PARTTYPES*fcent+2][0]->Fill(DCAXY,track->Pt());
	    fSecWeak_DCAz_Pt[PARTTYPES*fcent+2][0]->Fill(DCAZ,track->Pt());    
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+2][0]->Fill(track->Eta(),track->Pt()); 
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+2][0]->GetXaxis()->SetRangeUser(-1,1);
	}
      else if(MCtrk->IsSecondaryFromWeakDecay() && (isKaonNsigma && PDGcode==-321) )
	{
	    fSecWeak_DCAxy_Pt[PARTTYPES*fcent+2][1]->Fill(DCAXY,track->Pt());
	    fSecWeak_DCAz_Pt[PARTTYPES*fcent+2][1]->Fill(DCAZ,track->Pt());    
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+2][1]->Fill(track->Eta(),track->Pt()); 
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+2][1]->GetXaxis()->SetRangeUser(-1,1);
	 
	}
      else if(MCtrk->IsSecondaryFromMaterial() && (isKaonNsigma && PDGcode==321) )
	{

	    fSecMat_DCAxy_Pt[PARTTYPES*fcent+2][0]->Fill(DCAXY,track->Pt());
	    fSecMat_DCAz_Pt[PARTTYPES*fcent+2][0]->Fill(DCAZ,track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+2][0]->Fill(track->Eta(),track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+2][0]->GetXaxis()->SetRangeUser(-1,1);

	}
      else if(MCtrk->IsSecondaryFromMaterial() && (isKaonNsigma && PDGcode==-321) )
	{
	    fSecMat_DCAxy_Pt[PARTTYPES*fcent+2][1]->Fill(DCAXY,track->Pt());
	    fSecMat_DCAz_Pt[PARTTYPES*fcent+2][1]->Fill(DCAZ,track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+2][1]->Fill(track->Eta(),track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+2][1]->GetXaxis()->SetRangeUser(-1,1);
	}        
        
        
     if (MCtrk->IsPhysicalPrimary() && (isPionNsigma && PDGcode==211))
	{
	    fPrim_DCAxy_Pt[PARTTYPES*fcent+1][0]->Fill(DCAXY,track->Pt());
	    fPrim_DCAz_Pt[PARTTYPES*fcent+1][0]->Fill(DCAZ,track->Pt());  
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+1][0]->Fill(track->Eta(),track->Pt()); 
 	    fPrimVsCosPointingAngle[PARTTYPES*fcent+1][0]->GetXaxis()->SetRangeUser(-1,1);
	}
      else if (MCtrk->IsPhysicalPrimary() && (isPionNsigma && PDGcode==-211))
	{
	    fPrim_DCAxy_Pt[PARTTYPES*fcent+1][1]->Fill(DCAXY,track->Pt());
	    fPrim_DCAz_Pt[PARTTYPES*fcent+1][1]->Fill(DCAZ,track->Pt());  
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+1][1]->Fill(track->Eta(),track->Pt()); 
 	    fPrimVsCosPointingAngle[PARTTYPES*fcent+1][1]->GetXaxis()->SetRangeUser(-1,1);

	}
      else if(MCtrk->IsSecondaryFromWeakDecay() && (isPionNsigma && PDGcode==211) )
	{
	    fSecWeak_DCAxy_Pt[PARTTYPES*fcent+1][0]->Fill(DCAXY,track->Pt());
	    fSecWeak_DCAz_Pt[PARTTYPES*fcent+1][0]->Fill(DCAZ,track->Pt());    
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+1][0]->Fill(track->Eta(),track->Pt()); 
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+1][0]->GetXaxis()->SetRangeUser(-1,1);
	}
      else if(MCtrk->IsSecondaryFromWeakDecay() && (isPionNsigma && PDGcode==-211) )
	{
	    fSecWeak_DCAxy_Pt[PARTTYPES*fcent+1][1]->Fill(DCAXY,track->Pt());
	    fSecWeak_DCAz_Pt[PARTTYPES*fcent+1][1]->Fill(DCAZ,track->Pt());    
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+1][1]->Fill(track->Eta(),track->Pt()); 
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+1][1]->GetXaxis()->SetRangeUser(-1,1);
	 
	}
      else if(MCtrk->IsSecondaryFromMaterial() && (isPionNsigma && PDGcode==211) )
	{

	    fSecMat_DCAxy_Pt[PARTTYPES*fcent+1][0]->Fill(DCAXY,track->Pt());
	    fSecMat_DCAz_Pt[PARTTYPES*fcent+1][0]->Fill(DCAZ,track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+1][0]->Fill(track->Eta(),track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+1][0]->GetXaxis()->SetRangeUser(-1,1);

	}
      else if(MCtrk->IsSecondaryFromMaterial() && (isPionNsigma && PDGcode==-211) )
	{
	    fSecMat_DCAxy_Pt[PARTTYPES*fcent+1][1]->Fill(DCAXY,track->Pt());
	    fSecMat_DCAz_Pt[PARTTYPES*fcent+1][1]->Fill(DCAZ,track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+1][1]->Fill(track->Eta(),track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+1][1]->GetXaxis()->SetRangeUser(-1,1);
	}
        
         

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //step 1, TOF Matching
    UInt_t statusTOF;
    statusTOF=track->GetStatus();
    if((statusTOF&AliVTrack::kTOFout)==0 || (statusTOF&AliVTrack::kTIME)==0)
      statusTOF=0;
    if(track->Pt()<0.5) statusTOF = 1;

    //Misidentification fraction
    if(abs(PDGcode)==211)
      {
	if(isPionNsigma)
	  fMisidentification[fcent][charge]-> Fill(1,0.5,track->Pt());
	if(isKaonNsigma)
	  fMisidentification[fcent][charge]-> Fill(1,1.5,track->Pt());
	if(isProtonNsigma)
	  fMisidentification[fcent][charge]-> Fill(1,2.5,track->Pt());
	if(isDeuteronNsigma)
	  fMisidentification[fcent][charge]-> Fill(1,3.5,track->Pt());
	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma && !isDeuteronNsigma)
	  if(statusTOF)
	    fMisidentification[fcent][charge]-> Fill(1,4.5,track->Pt());
      }
    else if(abs(PDGcode)==321)
      {
	if(isPionNsigma)
	  fMisidentification[fcent][charge]-> Fill(2,0.5,track->Pt());
	if(isKaonNsigma)
	  fMisidentification[fcent][charge]-> Fill(2,1.5,track->Pt());
	if(isProtonNsigma)
	  fMisidentification[fcent][charge]-> Fill(2,2.5,track->Pt());
	if(isDeuteronNsigma)
	  fMisidentification[fcent][charge]-> Fill(2,3.5,track->Pt());
	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma && !isDeuteronNsigma)
	  if(statusTOF)
	    fMisidentification[fcent][charge]-> Fill(2,4.5,track->Pt());
      }
    else if(abs(PDGcode) == 2212)
      {
	if(isPionNsigma)
	  fMisidentification[fcent][charge]-> Fill(3,0.5,track->Pt());
	if(isKaonNsigma)
	  fMisidentification[fcent][charge]-> Fill(3,1.5,track->Pt());
	if(isProtonNsigma)
	  fMisidentification[fcent][charge]-> Fill(3,2.5,track->Pt());
	if(isDeuteronNsigma)
	  fMisidentification[fcent][charge]-> Fill(3,3.5,track->Pt());  
	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma && !isDeuteronNsigma)
	  if(statusTOF)
	    fMisidentification[fcent][charge]-> Fill(3,4.5,track->Pt());
      }
    else if(abs(PDGcode) == 1000010020)
      {
	if(isPionNsigma)
	  fMisidentification[fcent][charge]-> Fill(4,0.5,track->Pt());
	if(isKaonNsigma)
	  fMisidentification[fcent][charge]-> Fill(4,1.5,track->Pt());
	if(isProtonNsigma)
	  fMisidentification[fcent][charge]-> Fill(4,2.5,track->Pt());
	if(isDeuteronNsigma)
	  fMisidentification[fcent][charge]-> Fill(4,3.5,track->Pt());  
	if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma && !isDeuteronNsigma)
	  if(statusTOF)
	    fMisidentification[fcent][charge]-> Fill(4,4.5,track->Pt());
      }

    if(PDGcode==1000010020)
      PDGcode=777;
    else if(PDGcode==-1000010020)
      PDGcode=-777;
    
      fContamination[PARTTYPES*fcent][charge]-> Fill(PDGcode,track->Pt());
      
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
    if(isDeuteronNsigma)
      {
	fContamination[PARTTYPES*fcent+4][charge]-> Fill(PDGcode,track->Pt());// filling contamination histogram for deuterons
      }
  
    if(PDGcode==777)
      PDGcode=1000010020;
    else if(PDGcode==-777)
      PDGcode=-1000010020;
      
  

    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

     if (isDeuteronNsigma){
       // cout<<"deuteron"<<endl;
       if (!MCtrk) continue;
       recoParticleArray[4].Add(MCtrk);
       }

      //Fills for all identified deuterons found after cuts (reconstructed) - numerator for Efficiency
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
      fReconstructedNotPrimaries[PARTTYPES*fcent][charge]->Fill(track->Eta(), track->Pt());
    }
    else{
      fReconstructedPrimaries[PARTTYPES*fcent][charge]->Fill(track->Eta(), track->Pt());
    }

 
//    int PDGcode = MCtrk->GetPdgCode();
   //And secondaries for different particle species:



///////////////////////////////////////////////////////////////////////////////////////////////////////////

     if (MCtrk->IsPhysicalPrimary() && (isDeuteronNsigma && PDGcode==1000010020))
	{
	    fPrim_DCAxy_Pt[PARTTYPES*fcent+4][0]->Fill(DCAXY,track->Pt());
	    fPrim_DCAz_Pt[PARTTYPES*fcent+4][0]->Fill(DCAZ,track->Pt());  
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(track->Eta(),track->Pt()); 
 	    fPrimVsCosPointingAngle[PARTTYPES*fcent+4][0]->GetXaxis()->SetRangeUser(-1,1);
	}
      else if (MCtrk->IsPhysicalPrimary() && (isDeuteronNsigma && PDGcode==-1000010020))
	{
	    fPrim_DCAxy_Pt[PARTTYPES*fcent+4][1]->Fill(DCAXY,track->Pt());
	    fPrim_DCAz_Pt[PARTTYPES*fcent+4][1]->Fill(DCAZ,track->Pt());  
	    fPrimVsCosPointingAngle[PARTTYPES*fcent+4][1]->Fill(track->Eta(),track->Pt()); 
 	    fPrimVsCosPointingAngle[PARTTYPES*fcent+4][1]->GetXaxis()->SetRangeUser(-1,1);

	}
      else if(MCtrk->IsSecondaryFromWeakDecay() && (isDeuteronNsigma && PDGcode==1000010020) )
	{
	    fSecWeak_DCAxy_Pt[PARTTYPES*fcent+4][0]->Fill(DCAXY,track->Pt());
	    fSecWeak_DCAz_Pt[PARTTYPES*fcent+4][0]->Fill(DCAZ,track->Pt());    
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(track->Eta(),track->Pt()); 
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+4][0]->GetXaxis()->SetRangeUser(-1,1);
	}
      else if(MCtrk->IsSecondaryFromWeakDecay() && (isDeuteronNsigma && PDGcode==-1000010020) )
	{
	    fSecWeak_DCAxy_Pt[PARTTYPES*fcent+4][1]->Fill(DCAXY,track->Pt());
	    fSecWeak_DCAz_Pt[PARTTYPES*fcent+4][1]->Fill(DCAZ,track->Pt());    
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+4][1]->Fill(track->Eta(),track->Pt()); 
	    fSecWeakVsCosPointingAngle[PARTTYPES*fcent+4][1]->GetXaxis()->SetRangeUser(-1,1);
	 
	}
      else if(MCtrk->IsSecondaryFromMaterial() && (isDeuteronNsigma && PDGcode==1000010020) )
	{

	    fSecMat_DCAxy_Pt[PARTTYPES*fcent+4][0]->Fill(DCAXY,track->Pt());
	    fSecMat_DCAz_Pt[PARTTYPES*fcent+4][0]->Fill(DCAZ,track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+4][0]->Fill(track->Eta(),track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+4][0]->GetXaxis()->SetRangeUser(-1,1);

	}
      else if(MCtrk->IsSecondaryFromMaterial() && (isDeuteronNsigma && PDGcode==-1000010020) )
	{
	    fSecMat_DCAxy_Pt[PARTTYPES*fcent+4][1]->Fill(DCAXY,track->Pt());
	    fSecMat_DCAz_Pt[PARTTYPES*fcent+4][1]->Fill(DCAZ,track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+4][1]->Fill(track->Eta(),track->Pt()); 
	    fSecMatVsCosPointingAngle[PARTTYPES*fcent+4][1]->GetXaxis()->SetRangeUser(-1,1);
	}

}



  // MONTECARLO PARTICLES 
  if(!arrayMC){
    AliError("Array of MC particles not found");
    return;
  }
  // loop over MC stack 
  for (Int_t ipart = 0; ipart < arrayMC->GetEntries(); ipart++) {
    //std::cout<<"Entered MC loop"<<std::endl;
    
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(ipart);

    if (!MCtrk) continue;
    //std::cout<<"particle obtained"<<std::endl;
    Int_t PDGcode = TMath::Abs(MCtrk->GetPdgCode()); 
    Int_t charge=0;

    if(PDGcode < 0) charge=1;
    else if(PDGcode > 0) charge=0;
    PDGcode = TMath::Abs(MCtrk->GetPdgCode()); 



      //*** PID - check if pion ***

      if(MCtrk->Eta() < -0.8 || MCtrk->Eta() > 0.8){
	continue; }
	
      if(MCtrk->Pt() < 0.1 || MCtrk->Pt() > 2.5){
	continue;}


      // check physical primary 
      if(MCtrk->IsPhysicalPrimary() || (PDGcode==1000010020 && !MCtrk->IsSecondaryFromMaterial())) // Not from weak decay!
	{
	    int fcent=0;
	// Filling histograms for MC truth particles
	fGeneratedMCPrimaries[fcent*PARTTYPES][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());

	 Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	 fGeneratedMCPrimaries4D[fcent*PARTTYPES][charge]->Fill(val);
	 
	 if(PDGcode==211){
	  fGeneratedMCPrimaries[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  fGeneratedMCPrimaries4D[fcent*PARTTYPES+1][charge]->Fill(val);}
	 else if(PDGcode==321){
	  fGeneratedMCPrimaries[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  fGeneratedMCPrimaries4D[fcent*PARTTYPES+2][charge]->Fill(val);}
	 else if(PDGcode==2212){
	  fGeneratedMCPrimaries[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  fGeneratedMCPrimaries4D[fcent*PARTTYPES+3][charge]->Fill(val);}
	 else if(PDGcode==1000010020){
	  fGeneratedMCPrimaries[fcent*PARTTYPES+4][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  fGeneratedMCPrimaries4D[fcent*PARTTYPES+4][charge]->Fill(val);}
	  //Filling data from MC truth particles only for particles that were reconstruced
	if (recoParticleArray[0].Contains(MCtrk)){ //All
	  fMCPrimariesThatAreReconstructed[fcent*PARTTYPES][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  
	  Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	  fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES][charge]->Fill(val);
	  
	  fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  if(PDGcode==211)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  if(PDGcode==321)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  if(PDGcode==2212)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  if(PDGcode==1000010020)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+4][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	}
	if (recoParticleArray[1].Contains(MCtrk)){ //Pions
	  if(PDGcode==211){
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	    Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+1][charge]->Fill(val);
	  }
	}
	if (recoParticleArray[2].Contains(MCtrk)){ //Kaons
	  if(PDGcode==321){
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	    Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+2][charge]->Fill(val);
	  }
	}
	if (recoParticleArray[3].Contains(MCtrk)){ //Protons
	  if(PDGcode==2212){
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	    Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+3][charge]->Fill(val);
	  }
	}
       if (recoParticleArray[4].Contains(MCtrk)){ //Deuterons
	  if(PDGcode==1000010020){
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+4][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	    Double_t val[] = {MCtrk->Eta(), MCtrk->Pt(), MCtrk->Zv() ,MCtrk->Phi()};
	    fMCPrimariesThatAreReconstructed4D[fcent*PARTTYPES+4][charge]->Fill(val);
	  }
	}


      }
    
   }
    
  PostData(1, fHistoList);
}






