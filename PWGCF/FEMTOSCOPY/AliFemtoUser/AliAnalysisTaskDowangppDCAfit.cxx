#include "AliAnalysisTaskDowangppDCAfit.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TMatrixDSym.h>
#include <TF1.h>
#include <TRandom.h>
#include <TVector3.h>
#include <THnSparse.h>
#include <TProfile2D.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

// AliRoot includes
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODHeader.h"
#include "AliAODVZERO.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliOADBContainer.h"
#include "AliAODv0.h"
#include "AliMultSelection.h"
#include "AliAODcascade.h"
#include "AliAODTracklets.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"

// STL includes
#include <iostream>
#include <ctime>
#include <sys/time.h>
#include<vector> 
using namespace std;
using std::cout;
using std::endl;
const float MassP = 0.938272013;

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623
#define fphiL -1.4137167	//default 20bin for phi&eta
#define fphiT 4.8694686
static float TPCradii[9] = { 0.85, 1.05, 1.25, 1.45, 1.65, 1.85, 2.05, 2.25, 2.45 };


ClassImp(AliAnalysisTaskDowangppDCAfit)
//_____________________________________________________________________________
AliAnalysisTaskDowangppDCAfit::AliAnalysisTaskDowangppDCAfit():
  AliAnalysisTaskSE(),
  fAOD(0),
  mcEvent(0),
  fMCarray(0),
  fMCheader(0),
  fAnalysisUtil(0),
  fAODpidUtil(0),
  fPIDResponse(0),
  aodH(0),
  fVtxCut(10.0),  
  fFilterbit(128),
  fEtaCut(0.8),
  fNoClus(70),
  fMinPt(0.2),
  fMaxPt(5.0),
  fListOfObjects(0)
{
    
    EventDis = nullptr;
    /*
    for(int ic=0;ic<9;ic++){
if(ic!=0) continue;
    for(int i=0;i<8;i++){
        for(int j=0;j<3;j++){

            DCAxy_noCut[i][j][ic] = nullptr;
            DCAz_noCut[i][j][ic] = nullptr;

            DCAxy_check_method[i][j][ic] = nullptr;
            DCAz_check_method[i][j][ic] = nullptr;

            DCAxy_allCut[i][j][ic] = nullptr;
            DCAz_allCut[i][j][ic] = nullptr;
            dEdxVspT[i][j][ic] = nullptr;
            m2TOFVspT[i][j][ic] = nullptr;
            pTdis[i][j][ic] = nullptr;
            Etadis[i][j][ic] = nullptr;

            // TH2F * dEdxVspT_noCut[8][3][1];
            //TH2F * dEdxVspT_allCut[8][3][1];
            dEdxVspT_noCut[i][j][ic] = nullptr;
            dEdxVspT_allCut[i][j][ic] = nullptr;


        }
    }
    }

    fHistQA_Event = nullptr;
    fHistQA_Track = nullptr;
   
    for(int j=0;j<3;j++){
        PileUPTracl_QA[j] = nullptr;
        fHistQA_pInEvent[j]  = nullptr;
    }

    for(int i=0;i<8;i++){
        MCTrackQA[i] = nullptr;
    }
    */
    // ddd
    for(int ic=0;ic<5;ic++){
       for(int ichg=0;ichg<2;ichg++){
            for(int i=0;i<8;i++){
                DCAxy_pDetailFrac[ichg][i][ic] = nullptr;
            }
        }
    }
    //ddd
    for(int ic=0;ic<5;ic++){
        for(int ichg=0;ichg<2;ichg++){
            MomSmearing[ichg][ic] = nullptr;
            MomSmearing_mix[ichg][ic] = nullptr;
        }
    }

    // ddd
    for(int ic=0;ic<5;ic++){
        for(int ichg=0;ichg<2;ichg++){
            pTdisReco[ichg][ic] = nullptr;
            //pTdisTrue[ichg][ic] = nullptr;
            pTdisRecoWhenPrimary[ichg][ic] = nullptr;
            pTdisRecoIDAsP[ichg][ic] = nullptr;
        }
    }

    for(int ic=0;ic<5;ic++){
        for(int ichg=0;ichg<2;ichg++){
            fNumDPhiDEtaAvgQA[ichg][ic] = nullptr;
            fDumDPhiDEtaAvgQA[ichg][ic] = nullptr;

            fNumDPhiDEtaAvgQA_afterPairCut[ichg][ic] = nullptr;
            fDumDPhiDEtaAvgQA_afterPairCut[ichg][ic] = nullptr;

        }
    }

    for(int ic=0;ic<5;ic++){
        for(int ichg=0;ichg<2;ichg++){
            kStarVskT2DinMixPID[ichg][ic] = nullptr;
            kStarVskT2DinMixTrue[ichg][ic] = nullptr;
            kStarVsmT2DinMixPID[ichg][ic] = nullptr;
            kStarVsmT2DinMixTrue[ichg][ic] = nullptr;

        }
    }



    // mix pool
    pool_bin empty_pool;

    for (int ip=0;ip<numOfVertexZBins*numOfCent*numOfMultBins;ip++) {
        All_Event_pool.push_back(empty_pool);
    }

}
//______________________________________________________________________________
AliAnalysisTaskDowangppDCAfit::AliAnalysisTaskDowangppDCAfit(const char *name):
AliAnalysisTaskSE(name),
fAOD(0),
mcEvent(0),
fMCarray(0),
fMCheader(0),
fAnalysisUtil(0),
fAODpidUtil(0),
fPIDResponse(0),
aodH(0),
fVtxCut(10.0),
fFilterbit(128),
fEtaCut(0.8),
fNoClus(70),
fMinPt(0.2),
fMaxPt(5.0),
fListOfObjects(0)
{

    // Output slot #1 writes into a TTree
    EventDis = nullptr;
    /*
    for(int ic=0;ic<9;ic++){
    if(ic!=0) continue;
    for(int i=0;i<8;i++){
        for(int j=0;j<3;j++){

            DCAxy_noCut[i][j][ic] = nullptr;
            DCAz_noCut[i][j][ic] = nullptr;

            DCAxy_check_method[i][j][ic] = nullptr;
            DCAz_check_method[i][j][ic] = nullptr;

            DCAxy_allCut[i][j][ic] = nullptr;
            DCAz_allCut[i][j][ic] = nullptr;
            dEdxVspT[i][j][ic] = nullptr;
            m2TOFVspT[i][j][ic] = nullptr;
            pTdis[i][j][ic] = nullptr;
            Etadis[i][j][ic] = nullptr;

            // TH2F * dEdxVspT_noCut[8][3][1];
            //TH2F * dEdxVspT_allCut[8][3][1];
            dEdxVspT_noCut[i][j][ic] = nullptr;
            dEdxVspT_allCut[i][j][ic] = nullptr;


        }
    }
    }

    fHistQA_Event = nullptr;
    fHistQA_Track = nullptr;
   
    for(int j=0;j<3;j++){
        PileUPTracl_QA[j] = nullptr;
        fHistQA_pInEvent[j]  = nullptr;
    }

    for(int i=0;i<8;i++){
        MCTrackQA[i] = nullptr;
    }
    */
   
    // ddd
    for(int ic=0;ic<5;ic++){
       for(int ichg=0;ichg<2;ichg++){
            for(int i=0;i<8;i++){
                DCAxy_pDetailFrac[ichg][i][ic] = nullptr;
            }
        }
    }
    //ddd
    for(int ic=0;ic<5;ic++){
       //if(ic > 0) continue;
        for(int ichg=0;ichg<2;ichg++){
            MomSmearing[ichg][ic] = nullptr;
            MomSmearing_mix[ichg][ic] = nullptr;
        }
    }

    // ddd
    for(int ic=0;ic<5;ic++){
        //if(ic > 0) continue;
        for(int ichg=0;ichg<2;ichg++){
            pTdisReco[ichg][ic] = nullptr;
            //pTdisTrue[ichg][ic] = nullptr;
            pTdisRecoWhenPrimary[ichg][ic] = nullptr;
            pTdisRecoIDAsP[ichg][ic] = nullptr;
        }
    }

        for(int ic=0;ic<5;ic++){
        //if(ic > 0) continue;
        for(int ichg=0;ichg<2;ichg++){
            fNumDPhiDEtaAvgQA[ichg][ic] = nullptr;
            fDumDPhiDEtaAvgQA[ichg][ic] = nullptr;

            fNumDPhiDEtaAvgQA_afterPairCut[ichg][ic] = nullptr;
            fDumDPhiDEtaAvgQA_afterPairCut[ichg][ic] = nullptr;

        }
    }
    for(int ic=0;ic<5;ic++){
        //if(ic > 0) continue;
        for(int ichg=0;ichg<2;ichg++){
            kStarVskT2DinMixPID[ichg][ic] = nullptr;
            kStarVskT2DinMixTrue[ichg][ic] = nullptr;
            kStarVsmT2DinMixPID[ichg][ic] = nullptr;
            kStarVsmT2DinMixTrue[ichg][ic] = nullptr;
        }
    }

    // 默认每个事件都至少有2个p,2anti-p
    pool_bin empty_pool;
    for (int ip=0;ip<numOfVertexZBins*numOfCent*numOfMultBins;ip++) {
        All_Event_pool.push_back(empty_pool);
    }

        


    DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnalysisTaskDowangppDCAfit::~AliAnalysisTaskDowangppDCAfit()
{
    // Destructor
    if (fListOfObjects)
        delete fListOfObjects;
  
    if (fAnalysisUtil)
        delete fAnalysisUtil;
}

//______________________________________________________________________________
void AliAnalysisTaskDowangppDCAfit::UserCreateOutputObjects()
{ 

    fAnalysisUtil = new AliAnalysisUtils();

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();


    OpenFile(1);
    fListOfObjects = new TList();
    fListOfObjects->SetOwner();

    EventDis =new TH2F("fHistEventDis"," ",100,0,100,10,0,10);
    fListOfObjects->Add(EventDis);


    TString filename;
    /*
         for(int ic=0;ic<9;ic++){
        if(ic!=0) continue;
    for(int i=0;i<8;i++){
        for(int j=0;j<3;j++){

            filename.Form("DCAxy_noCut_%d_%d_%d",i,j,ic);
            DCAxy_noCut[i][j][ic] = new TH2F(filename," ",400, -3.0, 3.0, 100,0.0,5.0);
            filename.Form("DCAz_noCut_%d_%d_%d",i,j,ic);
            DCAz_noCut[i][j][ic] = new TH2F(filename," ",800, -4.0, 4.0, 100,0.0,5.0);

            filename.Form("DCAxy_check_method_%d_%d_%d",i,j,ic);
            DCAxy_check_method[i][j][ic] = new TH2F(filename," ",400, -3.0, 3.0, 100,0.0,5.0);
            filename.Form("DCAz_check_method_%d_%d_%d",i,j,ic);
            DCAz_check_method[i][j][ic] = new TH2F(filename," ",800, -4.0, 4.0, 100,0.0,5.0);
            //fListOfObjects->Add(DCAxy_check_method[i][j][ic]);
            //fListOfObjects->Add(DCAz_check_method[i][j][ic]);


            filename.Form("DCAxy_allCut_%d_%d_%d",i,j,ic);
            DCAxy_allCut[i][j][ic] = new TH2F(filename," ",400, -3.0, 3.0, 100,0.0,5.0);
            filename.Form("DCAz_allCut_%d_%d_%d",i,j,ic);
            DCAz_allCut[i][j][ic] = new TH2F(filename," ",800, -4.0, 4.0, 100,0.0,5.0);
          

            filename.Form("dEdxVspT_%d_%d_%d",i,j,ic);
            dEdxVspT[i][j][ic] = new TH2F(filename," ",100,0.0,5.0,250,0.0,500.0);

            filename.Form("m2TOFVspT_%d_%d_%d",i,j,ic);
            m2TOFVspT[i][j][ic] = new TH2F(filename," ",100,0.0,5.0,550,-1,10);

            filename.Form("pTdis_%d_%d_%d",i,j,ic);
            pTdis[i][j][ic] = new TH1F(filename," ",200,0.0,20.0);

            filename.Form("Etadis_%d_%d_%d",i,j,ic);
            Etadis[i][j][ic] = new TH1F(filename," ",100,-1.,1.0);
        
            filename.Form("dEdxVspT_noCut_%d_%d_%d",i,j,ic);
            dEdxVspT_noCut[i][j][ic] = new TH2F(filename," ",100,0.0,5.0,250,0.0,500.0);
            filename.Form("dEdxVspT_allCut_%d_%d_%d",i,j,ic);
            dEdxVspT_allCut[i][j][ic] = new TH2F(filename," ",100,0.0,5.0,250,0.0,500.0);


            // fListOfObjects->Add(DCAxy_allCut[i][j][ic]);
            // fListOfObjects->Add(DCAz_allCut[i][j][ic]);

            // fListOfObjects->Add(dEdxVspT[i][j][ic]);
            // fListOfObjects->Add(dEdxVspT_noCut[i][j][ic]);
            // fListOfObjects->Add(dEdxVspT_allCut[i][j][ic]);

            // fListOfObjects->Add(m2TOFVspT[i][j][ic]);
            // fListOfObjects->Add(pTdis[i][j][ic]);
            // fListOfObjects->Add(Etadis[i][j][ic]);


        }
    }   
    }

    fHistQA_Event = new TH1F("fHistEventCuts", "Event Cuts" , 4, 0, 5);
    fHistQA_Event->GetXaxis()->SetBinLabel(1,"All");
    fHistQA_Event->GetXaxis()->SetBinLabel(2,"NoVertex");
    fHistQA_Event->GetXaxis()->SetBinLabel(3,"PileUp");
    fHistQA_Event->GetXaxis()->SetBinLabel(4,"z-vertex>10");
    fListOfObjects->Add(fHistQA_Event);

    fHistQA_Track = new TH1F("fHistTrackCuts", "Track Cuts" , 7, 0.5, 7.5);
    fHistQA_Track->GetXaxis()->SetBinLabel(1,"AllTracksInEvents");
    fHistQA_Track->GetXaxis()->SetBinLabel(2,"GetTrack");
    fHistQA_Track->GetXaxis()->SetBinLabel(3,"Filter bit");
    fHistQA_Track->GetXaxis()->SetBinLabel(4,"Eta");
    fHistQA_Track->GetXaxis()->SetBinLabel(5,"Pt");
    fHistQA_Track->GetXaxis()->SetBinLabel(6,"DCA");
    fHistQA_Track->GetXaxis()->SetBinLabel(7,"Electron Rejection");
    fListOfObjects->Add(fHistQA_Track);

    for(int j=0;j<3;j++){
        filename.Form("fHistQA_pInEvent%d",j);
        fHistQA_pInEvent[j] = new TH1F(filename, "particle Cuts" , 8, 0.5, 8.5);
        fHistQA_pInEvent[j]->GetXaxis()->SetBinLabel(1,"p");
        fHistQA_pInEvent[j]->GetXaxis()->SetBinLabel(2,"anti-p");
        fHistQA_pInEvent[j]->GetXaxis()->SetBinLabel(3,"d");
        fHistQA_pInEvent[j]->GetXaxis()->SetBinLabel(4,"anti-d");
        fHistQA_pInEvent[j]->GetXaxis()->SetBinLabel(5,"t");
        fHistQA_pInEvent[j]->GetXaxis()->SetBinLabel(6,"anti-t");
        fHistQA_pInEvent[j]->GetXaxis()->SetBinLabel(7,"He3");
        fHistQA_pInEvent[j]->GetXaxis()->SetBinLabel(8,"anti-He3");
        fListOfObjects->Add(fHistQA_pInEvent[j]);
    
    }

    for(int j=0;j<3;j++){
        filename.Form("PileUPTracl_QA_%d",j);
        PileUPTracl_QA[j] = new TH1F(filename, " " , 2, 0, 2);
        fListOfObjects->Add(PileUPTracl_QA[j]);

       

    }

    for(int i=0;i<8;i++){
        filename.Form("MCTrackQA%d",i);
        MCTrackQA[i] = new TH1F(filename, " " , 8, 0.5, 8.5);
        fListOfObjects->Add(MCTrackQA[i]);
    }

    */
    for(int ic=0;ic<5;ic++){
        for(int ichg=0;ichg<2;ichg++){
            for(int i=0;i<8;i++){
                filename.Form("DCAxy_pDetailFrac%d%d%d",ichg,i,ic);
                DCAxy_pDetailFrac[ichg][i][ic] = new TH2F(filename," ",400, -3.0, 3.0, 100,0.0,5.0);
                fListOfObjects->Add(DCAxy_pDetailFrac[ichg][i][ic]);
            }
        }
    }
    for(int ic=0;ic<5;ic++){
        
        for(int ichg=0;ichg<2;ichg++){
            filename.Form("MomSmearing%d%d",ichg,ic);
            MomSmearing[ichg][ic] = new TH2F(filename," ",1000, 0.0, 1.0, 1000,0.0,1.0);
            fListOfObjects->Add(MomSmearing[ichg][ic]);

            filename.Form("MomSmearing_mix%d%d",ichg,ic);
            MomSmearing_mix[ichg][ic] = new TH2F(filename," ",1000, 0.0, 1.0, 1000,0.0,1.0);
            fListOfObjects->Add(MomSmearing_mix[ichg][ic]);


        }
    }
    for(int ic=0;ic<5;ic++){
       
        for(int ichg=0;ichg<2;ichg++){
            filename.Form("pTdisReco%d%d",ichg,ic);
            pTdisReco[ichg][ic] = new TH1F(filename," ",500,0.0,5.0);
            fListOfObjects->Add(pTdisReco[ichg][ic]);

            // filename.Form("pTdisTrue%d%d",ichg,ic);
            // pTdisTrue[ichg][ic] = new TH1F(filename," ",500,0.0,5.0);
            // fListOfObjects->Add(pTdisTrue[ichg][ic]);

            filename.Form("pTdisRecoWhenPrimary%d%d",ichg,ic);
            pTdisRecoWhenPrimary[ichg][ic] = new TH1F(filename," ",500,0.0,5.0);
            fListOfObjects->Add(pTdisRecoWhenPrimary[ichg][ic]);

            filename.Form("pTdisRecoIDAsP%d%d",ichg,ic);
            pTdisRecoIDAsP[ichg][ic] = new TH1F(filename," ",500,0.0,5.0);
            fListOfObjects->Add(pTdisRecoIDAsP[ichg][ic]);

            
        }
    }

    for(int ic=0;ic<5;ic++){
       
        for(int ichg=0;ichg<2;ichg++){
            filename.Form("fNumDPhiDEtaAvgQA%d%d",ichg,ic);
            fNumDPhiDEtaAvgQA[ichg][ic] = new TH2F(filename,"",300, -0.15, 0.15, 400, -0.2, 0.2);
            filename.Form("fDumDPhiDEtaAvgQA%d%d",ichg,ic);
            fDumDPhiDEtaAvgQA[ichg][ic] = new TH2F(filename,"",300, -0.15, 0.15, 400, -0.2, 0.2);

            fListOfObjects->Add(fNumDPhiDEtaAvgQA[ichg][ic]);
            fListOfObjects->Add(fDumDPhiDEtaAvgQA[ichg][ic]);

            filename.Form("fNumDPhiDEtaAvgQA_afterPairCut%d%d",ichg,ic);
            fNumDPhiDEtaAvgQA_afterPairCut[ichg][ic] = new TH2F(filename,"",300, -0.15, 0.15, 400, -0.2, 0.2);
            filename.Form("fDumDPhiDEtaAvgQA_afterPairCut%d%d",ichg,ic);
            fDumDPhiDEtaAvgQA_afterPairCut[ichg][ic] = new TH2F(filename,"",300, -0.15, 0.15, 400, -0.2, 0.2);

            fListOfObjects->Add(fNumDPhiDEtaAvgQA_afterPairCut[ichg][ic]);
            fListOfObjects->Add(fDumDPhiDEtaAvgQA_afterPairCut[ichg][ic]);


            
        }
    }
    for(int ic=0;ic<5;ic++){
      
        for(int ichg=0;ichg<2;ichg++){
            filename.Form("kStarVskT2DinMixPID%d%d",ichg,ic);
            kStarVskT2DinMixPID[ichg][ic] = new TH2F(filename," ",1000, 0.0, 1.0,100,0,5);
            filename.Form("kStarVskT2DinMixTrue%d%d",ichg,ic);
            kStarVskT2DinMixTrue[ichg][ic] = new TH2F(filename," ",1000, 0.0, 1.0,100,0,5);

            fListOfObjects->Add(kStarVskT2DinMixPID[ichg][ic]);
            fListOfObjects->Add(kStarVskT2DinMixTrue[ichg][ic]);

            filename.Form("kStarVsmT2DinMixPID%d%d",ichg,ic);
            kStarVsmT2DinMixPID[ichg][ic] = new TH2F(filename," ",1000, 0.0, 1.0,25,1,3.5);
            filename.Form("kStarVsmT2DinMixTrue%d%d",ichg,ic);
            kStarVsmT2DinMixTrue[ichg][ic] = new TH2F(filename," ",1000, 0.0, 1.0,25,1,3.5);

            fListOfObjects->Add(kStarVsmT2DinMixPID[ichg][ic]);
            fListOfObjects->Add(kStarVsmT2DinMixTrue[ichg][ic]);


        }
    }


    
    aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    fAODpidUtil = aodH->GetAODpidUtil();

    // Post output data.
    PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskDowangppDCAfit::UserExec(Option_t *) 
{
    float CentMax = 50.;
    float CentMin = -1.;
    //cout<<" Main loop"<<endl; 
    // Main loop
    // Called for each event
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD || !fAOD->GetHeader()){
        Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
        this->Dump();
        return;
    }
    
    

    mcEvent = MCEvent();
    //mcEvent = dynamic_cast<AliMCEvent*>(InputEvent());

    if (!mcEvent)
    {
        return;
    }
    else
    {
        //cout<<"AliMCEvent found!"<<endl;
    }

    
    fMCarray = (TClonesArray*)fAOD->FindListObject("mcparticles");
    if(!fMCarray){
        Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
        this->Dump();
        return;
    }
    
    fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject("mcHeader"));
    if(!fMCheader){
        Printf("%s:%d AOD MC header not found in Input Manager",(char*)__FILE__,__LINE__);
        this->Dump();
        return;
    }

  
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    fAODpidUtil = aodH->GetAODpidUtil();

    AliCentrality* alicent= fAOD->GetCentrality(); //in PbPb and pPb
    AliMultSelection *mult_selection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    Float_t v0Centr    = -100.;
    v0Centr = mult_selection->GetMultiplicityPercentile("V0M");
    if( v0Centr < CentMin || v0Centr >  CentMax) return;
	//if( v0Centr < 0 ||v0Centr >50) return;
    // fHistQA_Event->Fill(1);

    Float_t zvtx = AODEventCut(fAOD);
    if (TMath::Abs(zvtx) > 10.) return;
    //   fHistQA_Event->Fill(4);

    int centrCode = -10;
    centrCode = ReCentrCode(v0Centr);
    if (centrCode < 0)
        return;
	
    // if(centrCode!=0) return;

    EventDis->Fill(v0Centr,centrCode);
	Analyze(fAOD,v0Centr,zvtx);
    PostData(1, fListOfObjects);
    // Post output data.
    
}
float AliAnalysisTaskDowangppDCAfit::AODEventCut(AliAODEvent* fAOD){
    int iway = 2;
    float vtxz = -999;

    if(iway==2){
        const AliAODVertex* vertex =(AliAODVertex*) fAOD->GetPrimaryVertex();
        if (!vertex || vertex->GetNContributors()<=0)
            return vtxz;
        vtxz = vertex->GetZ();
        return vtxz;
    }
    return vtxz;
}
//________________________________________________________________________
void AliAnalysisTaskDowangppDCAfit::Analyze(AliAODEvent* aod,Float_t v0Centr,float zvtx)
{


    float Vz_pos = zvtx;
    int centrCode = -10;
    centrCode = ReCentrCode(v0Centr);
    int iCent = centrCode;
    float mag = aod->GetMagneticField();


    // loose DCAxy
    float DCAxyMax[8] = {
                    0.1,0.1,
                    0.2,100,
                    0.5,100,
                    0.5,100,
    };
    // strict DCAz
    float DCAzMax[8] = {
                    1.0,1.0,
                    1.0,3.2,
                    1.0,3.2,
                    1.0,3.2
    };

    // float p2LowPt_v[8] = {
    //             0.5,0.5,
    //             0.8,0.8,
    //             0.1,0.1,
    //             0.8,0.8
	// 			};
    // float p2UpPt_v[8] = {
    //             4.0,4.0,
    //             2.2,2.2,
    //             100.,100.,
    //             100.,100.
	// 			};




 // cout<<"vvv  "<<vertex->GetZ()<<endl;
    //====== getting MC array =====
    fMCarray = (TClonesArray*)fAOD->FindListObject("mcparticles");
    if(!fMCarray){
        Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
        this->Dump();
        return;
    }
   
    fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->FindListObject("mcHeader"));
    if(!fMCheader){
        Printf("%s:%d AOD MC header not found in Input Manager",(char*)__FILE__,__LINE__);
        this->Dump();
        return;
    }

    //cout<<"name "<<AliAODMCParticle::StdBranchName()<<endl;
    TClonesArray  *arrayMC;
    arrayMC = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));

    Int_t nMCtracks = mcEvent->GetNumberOfTracks();

    int nMCprimaries = mcEvent->GetNumberOfPrimaries();
    //===== MC track loop ====
    //int IsThisEventHaveMCP[8] = {0};
    //int IsThisEventHaveRealP[8] = {0};
//cout<<"howmany "<<nMCtracks<<endl;

    //vector<TLorentzVector> RecoTrackLabelV[2];
    //vector<TLorentzVector> TureTrackLabelV[2];

    vector<vector<p_info>> track_info(numOfChTypes);
    //vector<vector<TLorentzVector>> True_info(numOfChTypes);


    for (int itrack = 0; itrack < fAOD->GetNumberOfTracks(); itrack++){
    //            cout<<"sss "<<endl;
        float pTPC = 0.f;
        AliAODTrack *track = (AliAODTrack *)fAOD->GetTrack(itrack);

        if (!track)
        {
        cout << "Error: could not find track " << itrack << endl;
        continue;
        }

        // track cut
        if (!(track->TestFilterBit(1))) continue;
        
        if(track->GetTPCNcls() < 80) continue;

        int CrossRows = track->GetTPCCrossedRows();
        if(CrossRows < 70) continue;

        int NclsF = track->GetTPCNclsF();
        if(NclsF==0) continue;
        float tmp_ratio = (float)CrossRows/(float)NclsF;
        if(tmp_ratio<0.83) continue;

        int ITS_layer[6] = {0};
        for (int ii = 0; ii < 6; ii++) {
            if(track->HasPointOnITSLayer(ii)) ITS_layer[ii] = 1;
        }
        int SPD_hit = ITS_layer[0]+ITS_layer[1];
        if(SPD_hit==0) continue;

        float trackEta = track->Eta();
        if (TMath::Abs(trackEta) > 0.8) continue;

        Float_t chi2Tpc = track->Chi2perNDF();
        if(chi2Tpc > 4.) continue;

        float trackPt = track->Pt();
        if(trackPt < 0.5 || trackPt > 3.0) continue;

        float trackP = track->P();
        float TPCsignal = track->GetTPCsignal();

        double DCAxy = -99.;
        double DCAz = -99.;
        double dcaVals[2] = {-99., -99.};
        double covar[3] = {0., 0., 0.};
        AliAODTrack copy(*track);
        const AliVVertex *AODeventVtx = copy.GetAODEvent()->GetPrimaryVertex();
        Double_t AODeventMagneticF = copy.GetAODEvent()->GetMagneticField();
        if (copy.PropagateToDCA(AODeventVtx, AODeventMagneticF, 10, dcaVals, covar))
        {
            DCAxy = dcaVals[0];
            DCAz = dcaVals[1];
        }
        else
        {
            DCAxy = -99.; // track->DCA();
            DCAz = -99.;  // track->ZAtDCA();
        }
        //if (abs(DCAz) > DCAzMax[pLabel])  continue;
        if (abs(DCAz) > 1.0) continue;

        double nSigmaTPCpion    = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        double nsigmaTPCP       = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        // double nSigmaTPCD       = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kDeuteron);
        // double nSigmaTPCT       = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);
        // double nSigmaTPCH       = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);

        double nSigmaTOFpion    = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
        double nsigmaTOFP       = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
        // double nSigmaTOFD       = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);
        // double nSigmaTOFT       = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kTriton);
        // double nSigmaTOFH       = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kHe3);

        double TOFm2 = GetMass2TOF(GetTOFBeta(track), track);
        int imost = 0;
        
        if(TMath::Hypot(nsigmaTPCP,nsigmaTOFP) < 3.){
            imost = 1;
            if(WiolaRejectPion(trackP,nSigmaTPCpion, nSigmaTOFpion) ) imost = 0;
        }

        // 如果不是鉴别为p/pbar 就过
        if(imost==0) continue;
        // int pLabel = -999;
        
        // int pLabel = ReParticleLabel(pdgCode);
        // if(pLabel==-999) continue;
        // which fill label [10]!
        // (+ -) * (primary, week from Lambda, week from Sigma+, marterial)

        Int_t pidCode       = 0; // pid code
        Int_t primaryFlag   = 0; // 0 = not primary mc track
        Int_t pileupFlag    = 0; // 0 = pileup, 1 = no pileup

        int ChargeLabel = -999;
        int ThisTrackPDGCode = -999;
        if(track->Charge() > 0){
            ThisTrackPDGCode = 2212;
            ChargeLabel = 0;
        }
        if(track->Charge() < 0){
            ThisTrackPDGCode = -2212;
            ChargeLabel = 1;
        }
        // 筛除 track->Charge() == 0
        // 但是应该没有

        if(ChargeLabel==-999 || ThisTrackPDGCode==-999) continue;
        // cout<<"dddd "<<endl;  
        float tmp_DCAxyCut = 0.0105 + 0.0350 * TMath::Power(trackPt,-1.1);

       
        //vector<int> RecoTrackLabelV[2];
        //vector<int> TureTrackLabelV[2];
        const Int_t label = TMath::Abs(track->GetLabel());
        // AliAODMCParticle *mcTrack = (AliAODMCParticle *)mcEvent->GetTrack(iMCtrack);   
        //AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCarray->At(label));
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCarray->At(label));

        Int_t pdgCode = 0;
        // 只要有对应的MC track : primaryFlag 就会 > 0
        // 如果没有对应, primaryFlag = 0, ghost track
        // primaryFlat == 1: Primary + 是p/pbar
        // primaryFlat == 2: sec Lambda + 是p/pbar
        // primaryFlat == 3: sec Sigma+ + 是p/pbar
        // primaryFlat == 4: 其他 sec + 是p/pbar
        // primaryFlat == 5: 不存在
        // primaryFlat == 6: 不存在
        // primaryFlat == 7: 来自材料
        // primaryFlat == 8: 不是 p/pbar
        if (mcTrack){

            Int_t MClabel = mcTrack->Label();
            pdgCode = mcTrack->GetPdgCode();
            
            if(ThisTrackPDGCode == pdgCode){
                if (mcTrack->IsPhysicalPrimary() && MClabel>=0 &&  MClabel<=nMCprimaries){
                    primaryFlag = 1;
                }
                else if (mcTrack->IsSecondaryFromWeakDecay()){

                    if (mcTrack->GetMother() > -1) {
                        AliAODMCParticle *mother = static_cast<AliAODMCParticle *>(fMCarray->At(mcTrack->GetMother()));
                        // Check if this is the same particle stored twice on the stack
                        if (mother) {
                            if (mother->GetPdgCode() ==  track->Charge() * 3122){
                                primaryFlag = 2;
                            }
                            else if(mother->GetPdgCode() ==  track->Charge() * 3222){
                                primaryFlag = 3;
                            }
                            else{
                                primaryFlag = 4;
                            }
                        }
                        else{
                            primaryFlag = 5;
                        }
                    }
                    else{
                        primaryFlag = 6;
                    }
                }
                else if (mcTrack->IsSecondaryFromMaterial()){
                    primaryFlag = 7;
                }
            }
            else{
                // fake !
                primaryFlag = 8;
            } 
           
        }// switch  primaryFlag
        if(primaryFlag==0) continue;
        // DCA cut for get same track as inn real data
        if(abs(DCAxy) < tmp_DCAxyCut){

                pTdisReco[ChargeLabel][iCent]->Fill(track->Pt());
                TLorentzVector ParticleReco;
                ParticleReco.SetPtEtaPhiM(track->Pt(),track->Eta(),track->Phi(),MassP);
                TLorentzVector ParticleTure;
                ParticleTure.SetPtEtaPhiM(mcTrack->Pt(),mcTrack->Eta(),mcTrack->Phi(),MassP);
                // 只要鉴别为p or anti-p 都存了
                p_info tmp_info;
                //Save_p_info(tmp_info,track,ParticleReco,AODeventMagneticF); //same!
                Save_p_info(tmp_info,track,ParticleReco,ParticleTure,mag,primaryFlag);
                track_info[ChargeLabel].push_back(tmp_info);

                // 不是ghost 而且不是fake
                if(primaryFlag!=0 && primaryFlag!=8) pTdisRecoIDAsP[ChargeLabel][iCent]->Fill(trackPt);
        }

        // 只看纯P的部分
        if(primaryFlag > 0){
            //cout<<"here "<<primaryFlag<<endl;
            DCAxy_pDetailFrac[ChargeLabel][primaryFlag-1][iCent]->Fill(DCAxy,trackPt);
            if(primaryFlag==1 && abs(DCAxy) < tmp_DCAxyCut && pdgCode==ThisTrackPDGCode){
                pTdisRecoWhenPrimary[ChargeLabel][iCent]->Fill(trackPt);
            }
 
        }

    }// end event loop

 // same event loop
    int fir_index = 0;
    int sec_index = 0;
    bool pair_exsit[numOfpairTypes] = {false};


    for(int ichg=0;ichg<2;ichg++){
        int HowmanyTrack = track_info[ichg].size();

        for(int ip1=0;ip1<HowmanyTrack;ip1++){
            TLorentzVector Particle1_Reco = track_info[ichg][ip1].kin_info;
            TLorentzVector Particle1_True = track_info[ichg][ip1].MCture_info;
            for(int ip2=ip1+1;ip2<HowmanyTrack;ip2++){
                TLorentzVector Particle2_Reco = track_info[ichg][ip2].kin_info;
                TLorentzVector Particle2_True = track_info[ichg][ip2].MCture_info;

                float AvgDPhi = ReAvgDphi(track_info[ichg][ip1],track_info[ichg][ip2]);
	            double deta = Particle1_Reco.Eta() - Particle2_Reco.Eta();
                fNumDPhiDEtaAvgQA[ichg][iCent]->Fill(deta,AvgDPhi);

                if(Pass(track_info[ichg][ip1],track_info[ichg][ip2],aodH)){
                    float kStar_Reco = re_kstar(Particle1_Reco,Particle2_Reco);
                    float kStar_True = re_kstar(Particle1_True,Particle2_True);
                    //cout<<"ddd "<<kStar_Reco<<" "<<kStar_True<<endl;
                    MomSmearing[ichg][iCent]->Fill(kStar_True,kStar_Reco);
                    fNumDPhiDEtaAvgQA_afterPairCut[ichg][iCent]->Fill(deta,AvgDPhi);
                 
                }
            }
        }
    }

    // save mix event
    int event_inpool_index = 0;
    
    // input vz, cent index
    event_inpool_index = where_pool(Vz_pos,iCent,v0Centr);
    //cout<<"event_inpool_index "<<event_inpool_index<<" "<<zvtx<<" "<<iCent<<" "<<v0Centr<<endl;
    //aodH
    int first_p_index[numOfpairTypes] = {0,1};
    int second_p_index[numOfpairTypes] = {0,1};


    //\ begin mix
    int this_pair_pool_index = All_Event_pool[event_inpool_index].pool_index;
    if(this_pair_pool_index != 0){
        // p1 from mix and p2 current event!
        int HowManyEventInPool = this_pair_pool_index;
        for(int iE=0;iE<HowManyEventInPool;iE++){
            for(int ichg=0;ichg<2;ichg++){
                int p1MixSize = All_Event_pool[event_inpool_index].Event_in_poolBin[iE].alltrack[ichg].size();
                

                // PID track mix
                for(int ip1=0;ip1<p1MixSize;ip1++){
                    p_info mix_p1 = All_Event_pool[event_inpool_index].Event_in_poolBin[iE].alltrack[ichg][ip1];
                    TLorentzVector Particle1_Reco = mix_p1.kin_info;
                    TLorentzVector Particle1_True = mix_p1.MCture_info;
                    int Particle1_MCTruthTrackLabel = mix_p1.MCTruthTrackLabel;

                    for(int ip2=0;ip2<track_info[ichg].size();ip2++){

                            TLorentzVector Particle2_Reco = track_info[ichg][ip2].kin_info;
                            TLorentzVector Particle2_True = track_info[ichg][ip2].MCture_info;
                            int Particle2_MCTruthTrackLabel = track_info[ichg][ip2].MCTruthTrackLabel;

                            float kStar_Reco = re_kstar(Particle1_Reco,Particle2_Reco);
                            float kStar_True = re_kstar(Particle1_True,Particle2_True);

                            TLorentzVector tmp_pair = Particle1_Reco + Particle2_Reco;
                            float kT_tmp = 0.5 * tmp_pair.Pt();
                            float mT_tmp = 0.5 * tmp_pair.Mt();

                            float AvgDPhi = ReAvgDphi(mix_p1,track_info[ichg][ip2]);
                            double deta = Particle1_Reco.Eta() - Particle2_Reco.Eta();
                            fDumDPhiDEtaAvgQA[ichg][iCent]->Fill(deta,AvgDPhi);
                            if(Pass(mix_p1,track_info[ichg][ip2],aodH)){
                                //cout<<"ddd "<<kStar_Reco<<" "<<kStar_True<<endl;
                                MomSmearing_mix[ichg][iCent]->Fill(kStar_True,kStar_Reco);
                                fDumDPhiDEtaAvgQA_afterPairCut[ichg][iCent]->Fill(deta,AvgDPhi);
                                kStarVskT2DinMixPID[ichg][iCent]->Fill(kStar_Reco,kT_tmp);
                                kStarVsmT2DinMixPID[ichg][iCent]->Fill(kStar_Reco,mT_tmp);
                                // 这两条track 都不能是fake
                                if(Particle1_MCTruthTrackLabel!=8 && Particle2_MCTruthTrackLabel!=8){
                                    kStarVskT2DinMixTrue[ichg][iCent]->Fill(kStar_Reco,kT_tmp);
                                    kStarVsmT2DinMixTrue[ichg][iCent]->Fill(kStar_Reco,mT_tmp);
                                }
                                
                            }

                    }
                }// end PID track mix


            }
        }
    }
    // begin fill pool
    event_info this_event_info;
    this_event_info.vz = Vz_pos;
    this_event_info.cent = iCent;
    this_event_info.mag = mag;
    //this_event_info.alltrack = track_info;
    for(int is=0;is<numOfChTypes;is++){
        this_event_info.alltrack.push_back(track_info[is]);
        //this_event_info.PIDandTrueTrack.push_back(True_info[is]);
    }
    if (All_Event_pool[event_inpool_index].pool_index==maxNumEventsToMix) {
        // cout<<"mixing buffer is full "<<" "<<ipair<<endl;
        vector<event_info>::iterator pos;
        pos = All_Event_pool[event_inpool_index].Event_in_poolBin.begin();
        All_Event_pool[event_inpool_index].Event_in_poolBin.erase(pos);
        All_Event_pool[event_inpool_index].Event_in_poolBin.push_back(this_event_info);
    }
    else{
        All_Event_pool[event_inpool_index].Event_in_poolBin.push_back(this_event_info);
        All_Event_pool[event_inpool_index].pool_index++;
    }


/*
    // MC ture paricle loop
    for (int iMCtrack = 0; iMCtrack < nMCtracks; iMCtrack++)
    {
        AliAODMCParticle *part = (AliAODMCParticle *)mcEvent->GetTrack(iMCtrack);
        if (!part) continue;
        //
        Int_t pdgCode = part->PdgCode();
        if (not(pdgCode == 2212 or pdgCode == -2212)) continue;

        int ChargeLabel = -999;
        if(part->Charge() > 0){
            ChargeLabel = 0;
        }
        if(part->Charge() < 0){
            ChargeLabel = 1;
        }
        
        if(ChargeLabel==-999) continue;
        pTdisTrue[ChargeLabel][iCent]->Fill(part->Pt());

    }
*/
}


int AliAnalysisTaskDowangppDCAfit::ReParticleLabel(Int_t pdgCode){
    int pLabel = -999;
        switch (pdgCode)
        {
            case 2212:
                pLabel = 0;
                break;
            case -2212:
                pLabel = 1;
                break;
            case 1000010020:
                pLabel = 2;
                break;
            case -1000010020:
                pLabel = 3;
                break;
            case 1000010030:
                pLabel = 4;
                break;
            case -1000010030:
                pLabel = 5;
                break;
            case 1000020030:
                pLabel = 6;
                break;
            case -1000020030:
                pLabel = 7;
                break;
        }
    return pLabel;
}
Float_t AliAnalysisTaskDowangppDCAfit::GetTOFBeta(AliAODTrack *track)
{
  float beta = -999;
  double integratedTimes[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};

  track->GetIntegratedTimes(integratedTimes);

  const float c = 2.99792457999999984e-02;
  float p = track->P();
  float l = integratedTimes[0] * c;

  float trackT0 = fPIDResponse->GetTOFResponse().GetStartTime(p);

  float timeTOF = track->GetTOFsignal() - trackT0;
  if (timeTOF > 0)
  {
    beta = l / timeTOF / c;
  }
  return beta;
}

Float_t AliAnalysisTaskDowangppDCAfit::GetMass2TOF(Float_t beta, AliAODTrack *track)
{
  Float_t p = track->P();
  Float_t mass2sq = -999;
  if (!(beta == 0))
  {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}
   

void AliAnalysisTaskDowangppDCAfit::GetDCADis(AliAODEvent* fEvent,AliAODTrack *tAodTrack,float *DCA_re,int fDCAglobalTrack){

    if (fDCAglobalTrack == 0) {

        DCA_re[0] = tAodTrack->DCA();
        DCA_re[1] = tAodTrack->ZAtDCA();
    }
    if(fDCAglobalTrack==1){
        auto *vertex = static_cast<const AliAODVertex *>(fEvent->GetPrimaryVertex());
        float vertexX = -999.;
        float vertexY = -999.;
        float vertexZ = -999.;

        if (vertex) {
        Double_t fCov[6] = {0.0};
        vertex->GetCovarianceMatrix(fCov);
            if (fCov[5] != 0.0) {
                vertexX = vertex->GetX();
                vertexY = vertex->GetY();
                vertexZ = vertex->GetZ();
            }
        }

        Double_t pos[3];
        tAodTrack->GetXYZ(pos);

        Double_t DCAX = pos[0] - vertexX;
        Double_t DCAY = pos[1] - vertexY;
        Double_t DCAZ = pos[2] - vertexZ;

        Double_t DCAXY = TMath::Sqrt((DCAX * DCAX) + (DCAY * DCAY));
        DCA_re[0] = DCAXY;
        DCA_re[1] = DCAZ;
    }
if(fDCAglobalTrack==3){
float DCAXY = -999;
        float DCAZ = -999;
tAodTrack->GetImpactParameters(DCAXY,DCAZ);
  DCA_re[0] = DCAXY;
        DCA_re[1] = DCAZ;
}
    if(fDCAglobalTrack==2){
        Double_t DCAXY = -999;
        Double_t DCAZ = -999;

        //DCA for TPC only - from PropagateToDCA method
        AliExternalTrackParam aliextparam;
        aliextparam.CopyFromVTrack(tAodTrack);
        
        if (aliextparam.GetX() > 3.0) {
        DCAXY = -999;
        DCAZ = -999;
        } else {
        Double_t covar[3] = {0, 0, 0};
        Double_t DCA[2] = {0, 0};
        if (!aliextparam.PropagateToDCA(fEvent->GetPrimaryVertex(), fEvent->GetMagneticField(), 99999.0, DCA, covar)) {
            DCAXY = -999;
            DCAZ = -999;
        } else {
            DCAXY = DCA[0];
            DCAZ = DCA[1];
        }
        }
        DCA_re[0] = DCAXY;
        DCA_re[1] = DCAZ;
    }

}

//_____________________________________________________________________________
Float_t AliAnalysisTaskDowangppDCAfit::GetVertex(AliAODEvent* aod) const
{

    Float_t vtxz = -999;

    const AliAODVertex* trkVtx = aod->GetPrimaryVertex();
    if (!trkVtx || trkVtx->GetNContributors()<=0) 
        return vtxz;

    // const AliAODVertex* trkVtx = aod->GetPrimaryVertex();
    // if (!trkVtx || trkVtx->GetNContributors() < 2)
    //     return vtxz;
    
    // const AliAODVertex* spdVtx = aod->GetPrimaryVertexSPD();
    // if (!spdVtx || spdVtx->GetNContributors() < 1)
    //     return vtxz;
    // TString vtxTyp = spdVtx->GetTitle();
    // Double_t cov[6]={0};
    // spdVtx->GetCovarianceMatrix(cov);
    // Double_t zRes = TMath::Sqrt(cov[5]);
    // if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25) && (spdVtx->GetNContributors() < 20))
    //     return vtxz;
    
    
    vtxz = trkVtx->GetZ();
    
    return vtxz;
    
}
int AliAnalysisTaskDowangppDCAfit::GetImost(Short_t pidCode,int TrackCharge){
    int appImost = 0;
    int reImost = -999;
    if(TrackCharge<0) appImost = 1;
    switch (pidCode) {
        case 3:
            reImost = 0; // p
            break;
        case 13:
            reImost = 2; // d
            break;
        case 14:
            reImost = 4; // t
            break;
        case 15:
            reImost = 6; // He3
            break;
	default:
	    reImost = -999;
    }
    return reImost + appImost;
}
int AliAnalysisTaskDowangppDCAfit::GetImostByPDGCode(Int_t pdgCode){
    int reImost = -999;
    switch (pdgCode) {
        case 2212:
            reImost = 0;
            break;
        case -2212:
            reImost = 1;
            break;
        case 1000010020:    //  d
            reImost = 2;
            break; 
        case -1000010020:
            reImost = 3; 
            break;
        case 1000010030:    // t
            reImost = 4;
            break; 
        case -1000010030:  
            reImost = 5; 
            break;
        case 1000020030:    // He3
            reImost = 6;
            break; 
        case -1000020030:  
            reImost = 7; 
            break;
        default:
	        reImost = -999;
    };
    return reImost;

}
//_____________________________________________________________________________
Short_t AliAnalysisTaskDowangppDCAfit::GetPidCode(Int_t pdgCode) const
{
    // return our internal code for pions, kaons, protons. electrons, muons
    
    Short_t pidCode = 0;
    
    switch (pdgCode) {
        case 211:
            pidCode = 1; // pion
            break;
        case 321:
            pidCode = 2; // kaon
            break;
        case 2212:
            pidCode = 3; // proton
            break;
        case 11:
            pidCode = 4; // electron
            break;
        case 13:
            pidCode = 5; // muon
            break;
        case 1000010020:
            pidCode = 13;  // d
            break;
        case 1000010030:
            pidCode = 14;  // t
            break;
        case 1000020030:
            pidCode = 15;  // He3
            break;
        default:
            pidCode = 6;  // something else?
    };
    
    return pidCode;
}


//_____________________________________________________________________________
void AliAnalysisTaskDowangppDCAfit::Terminate(Option_t *)
{
  // Terminate loop
  Printf("Terminate()");
}

bool AliAnalysisTaskDowangppDCAfit::IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP, float nsigmaTPCD)
{
  if(TMath::Abs(nsigmaTPCE)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3 &&TMath::Abs(nsigmaTPCD)>3)
      return true;
   else
     return false;
}

bool AliAnalysisTaskDowangppDCAfit::IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime)
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
    

bool AliAnalysisTaskDowangppDCAfit::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
  if (mom > 0.5) {
  
    if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < 3)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCK) < 3)
      return true;
  }


  return false;
}

bool AliAnalysisTaskDowangppDCAfit::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{
  if (mom > 0.7) {
    if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 3)
        return true;
  }
  else{
    if (TMath::Abs(nsigmaTPCP) < 3) 
        return true;
  }



//   if (mom > 0.5 && mom< 3.0) {
//     if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 3)
//       return true;
//   }


  return false;
}

bool AliAnalysisTaskDowangppDCAfit::IsDeuteronNSigma(float mom, float nsigmaTPCD, float nsigmaTOFD,float swithch_mom)
{
    if (mom > swithch_mom) {
        if (TMath::Hypot( nsigmaTOFD, nsigmaTPCD  ) < 3)
            return true;
    }
    else{
        if (TMath::Abs(nsigmaTPCD) < 3) 
            return true;
    }
    // float sigmaMass = -1;
    // float maxmom = 1.8;
    // if (TMath::Abs(nsigmaTPCD) < 3){
    //      if(sigmaMass>0){
    //        return true;
    //      }  
    //      else{
    //        if( (mom > 0.8) && (nsigmaTPCD > (6*mom-9.3)) && (mom < maxmom) )
    //          return true;
    //      }
    //  }
  

    return false;
}
bool AliAnalysisTaskDowangppDCAfit::IsTritonNSigma(float mom, float nsigmaTPCT, float nsigmaTOFT)
{
    if (mom > 2.) {
        if (TMath::Hypot( nsigmaTPCT, nsigmaTOFT  ) < 3.)
            return true;
    }
    else{
        if (TMath::Abs(nsigmaTPCT) < 3.) 
            return true;
    }
   
    return false;
}
bool AliAnalysisTaskDowangppDCAfit::IsHe3NSigma(float mom, float nsigmaTPCH, float nsigmaTOFH)
{

    if (TMath::Abs(nsigmaTPCH) < 3.) 
        return true;
    
   
    return false;
}

bool AliAnalysisTaskDowangppDCAfit::IsKaonNSigmaReal(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
  if (mom > 1) {
    if ((TMath::Abs(nsigmaTOFK) < 1) && (TMath::Abs(nsigmaTPCK) < 3))
      return true;
  }
  else if (mom <1 && mom>0.8) {
    if ((TMath::Abs(nsigmaTOFK) < 1.5) && (TMath::Abs(nsigmaTPCK) < 3))
      return true;
  }
  else if (mom <0.8 && mom>0.45) {
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
bool AliAnalysisTaskDowangppDCAfit::RejectFakeP(double *NSigmaList, float mom){

	bool rejected = true;
	
    
	double k_NsigmaCombine = 0.;
	double pi_NsigmaCombine = 0.;
    double p_NsigmaCombine = 0.;
    double e_NsigmaCombine = 0.;

	if (mom > 0.7) {
        
		pi_NsigmaCombine    = TMath::Hypot(NSigmaList[0], NSigmaList[4]);
		k_NsigmaCombine     = TMath::Hypot(NSigmaList[1], NSigmaList[5]);
		p_NsigmaCombine     = TMath::Hypot(NSigmaList[2], NSigmaList[6]);
        e_NsigmaCombine     = TMath::Hypot(NSigmaList[3], NSigmaList[7]); 

	    if ( ( pi_NsigmaCombine < p_NsigmaCombine) || 
			( k_NsigmaCombine < p_NsigmaCombine)  ||
            ( e_NsigmaCombine < p_NsigmaCombine) )
	            return rejected;	
	}
    else {
        pi_NsigmaCombine    = abs(NSigmaList[0]);
        k_NsigmaCombine     = abs(NSigmaList[1]);
		p_NsigmaCombine     = abs(NSigmaList[2]);
        e_NsigmaCombine     = abs(NSigmaList[3]);
	        if ( (pi_NsigmaCombine < p_NsigmaCombine) ||
                (k_NsigmaCombine < p_NsigmaCombine) || 
               (e_NsigmaCombine < p_NsigmaCombine) )
	            return rejected;
    }
	return false;

}
int AliAnalysisTaskDowangppDCAfit::ReCentrCode(Float_t v0Centr){

    int centrCode = 999;
        if ((v0Centr >= 0) && (v0Centr < 10.))
        centrCode = 0;
    else if ((v0Centr >= 10.) && (v0Centr < 20.))
        centrCode = 1;
    else if ((v0Centr >= 20.) && (v0Centr < 30.))
        centrCode = 2;
    else if ((v0Centr >= 30.) && (v0Centr < 40.))
        centrCode = 3;
    else if ((v0Centr >= 40.) && (v0Centr < 50.))
        centrCode = 4;
    else if ((v0Centr >= 50.) && (v0Centr < 60.))
        centrCode = 5;
    else if ((v0Centr >= 60.) && (v0Centr < 70.))
        centrCode = 6;
    else if ((v0Centr >= 70.) && (v0Centr < 80.))
        centrCode = 7;
    else if ((v0Centr >= 80.) && (v0Centr < 90.))
        centrCode = 8;

    return centrCode;

}

bool AliAnalysisTaskDowangppDCAfit::WiolaDCut(float mom, float nsigmaTPCD, float nsigmaTOFD){

	if(mom < 1.3){
		if(abs(nsigmaTPCD) < 2)	return true;
	}
	else{
		if((abs(nsigmaTPCD) < 2) && (abs(nsigmaTOFD) < 2)) return true;
	}
	return false;
}
bool AliAnalysisTaskDowangppDCAfit::IsDeuteronTPCdEdx(float mom, float dEdx){
	double a1 = -250.0,  b1 = 400.0;
  	double a2 = -80,   b2 = 190.0;
  	if (mom < 1.1) {
    		if (dEdx < a1 * mom+b1) return false;
  	}
  	else if (mom < 2) {
    		if (dEdx < a2 *mom+b2) return false;
  	}
    	return true;

}
bool AliAnalysisTaskDowangppDCAfit::WiolaRejectPion(float mom,float nsigmaTPCpi,float nsigmaTOFpi){
	if (mom > 0.5) {
	        if (TMath::Hypot( nsigmaTOFpi, nsigmaTPCpi ) < 2.) return true;	
	}
    	else {
        	if (TMath::Abs(nsigmaTPCpi) < 2.) return true;
    	}
    	return false;

}


double AliAnalysisTaskDowangppDCAfit::re_kstar(TLorentzVector &Particle1, TLorentzVector &Particle2){
    double fKStarCalc = 999.;
    //if(Particle2.P()>1) return fKStarCalc;
    double pE1 = Particle1.E();
    double pE2 = Particle2.E();
    
    double px1 = Particle1.Px();
    double px2 = Particle2.Px();
    
    double py1 = Particle1.Py();
    double py2 = Particle2.Py();
    
    double pz1 = Particle1.Pz();
    double pz2 = Particle2.Pz();
    
//    double p1 = pt1/;
//    double p2 = 0.;
    double mass1 = Particle1.M();
    double mass2 = Particle2.M();
    
    double tPx = px1+px2;
    double tPy = py1+py2;
    double tPz = pz1+pz2;
    double tPE = pE1+pE2;
    
    double tPtrans = tPx*tPx + tPy*tPy;
    double tMtrans = tPE*tPE - tPz*tPz;
    
    
    double tPinv =   sqrt(tMtrans - tPtrans);
    double tQinvL = (pE1-pE2)*(pE1-pE2) - (px1-px2)*(px1-px2) - (py1-py2)*(py1-py2) - (pz1-pz2)*(pz1-pz2);
    
    double tQ = (mass1*mass1 - mass2*mass2)/tPinv;
    tQ = sqrt ( tQ*tQ - tQinvL);
    fKStarCalc = tQ/2;
    
    return fKStarCalc;
}
void AliAnalysisTaskDowangppDCAfit::Save_p_info(p_info &tmp_info,AliAODTrack *track, TLorentzVector RecoInfo, TLorentzVector ParticleTure, float MagneticField, int MCtureLabel){
    tmp_info.kin_info = RecoInfo;
    tmp_info.MCture_info = ParticleTure;
    float t_mag = MagneticField;
    tmp_info.mag = t_mag;

    tmp_info.MCTruthTrackLabel = MCtureLabel;

    tmp_info.charge = track->Charge();


    float globalPositionsAtRadii[9][3];
    GetGlobalPositionAtGlobalRadiiThroughTPC(track, t_mag, globalPositionsAtRadii);
    for (int i=0;i<9;i++) {
        tmp_info.TPC_Entrance[i][0] = globalPositionsAtRadii[i][0];
        tmp_info.TPC_Entrance[i][1] = globalPositionsAtRadii[i][1];
        tmp_info.TPC_Entrance[i][2] = globalPositionsAtRadii[i][2];
    }

    tmp_info.fClusters = track->GetTPCClusterMap();
    tmp_info.fShared = track->GetTPCSharedMap();

}
float AliAnalysisTaskDowangppDCAfit::re_mass(int PDGCode){
    switch(abs(PDGCode)){
        case 2212:
        return 0.9382720;
    
        case 1000010020:
        return 1.8756;
    
        case 1000010030:
        return 2.8089;
    
        case 1000020030:
        return 2.8084;
    }
}

int AliAnalysisTaskDowangppDCAfit::where_pool(float &vz,int iCent,float Multi){
    
    int ix = 0;
    int iy = iCent * numOfMultBins;
    
    float tmp_thisCent = Multi - float(iCent) * 10.;    //
    float CentBinWidth = 10./float(numOfMultBins);// 固定是2.5%
    int app_iy = (int)floor(tmp_thisCent/CentBinWidth);
    //cout<<"app_iy "<<app_iy<<endl;
    iy += app_iy;

    ix = (int)floor( (vz - Vz_low)/Vz_step );
    int bin = ix + numOfVertexZBins * iy;
    
    return bin;
}
bool AliAnalysisTaskDowangppDCAfit::Pass(p_info &first_p_info, p_info &sec_p_info,AliAODInputHandler *aodH){

    bool pair_pass = true;
    // AliFemtoPairCutMergedFraction::Pass
    // Prepare variables:
    pair_pass = CheckMergedFraction(fMagSign,first_p_info,sec_p_info,aodH);
    return pair_pass;
    
}
bool AliAnalysisTaskDowangppDCAfit::CheckMergedFraction(int &fMagSign,p_info &first_p_info, p_info &sec_p_info,AliAODInputHandler *aodH){
    //AliFemtoPairCutMergedFraction.cxx
    float phi1 = first_p_info.kin_info.Phi();
    float phi2 = sec_p_info.kin_info.Phi();
    float chg1 = first_p_info.charge;
    float chg2 = sec_p_info.charge;
    float pt1 = first_p_info.kin_info.Pt();
    float pt2 = sec_p_info.kin_info.Pt();
    float eta1 = first_p_info.kin_info.Eta();
    float eta2 = sec_p_info.kin_info.Eta();
    
    
    double magsign = 0.;
    
    if (!aodH) {
        return false;
    }
    else {
        AliAODEvent *fAOD;
        fAOD = aodH->GetEvent();
        magsign = fAOD->GetMagneticField();
    }
    
    //float magsign = this_event_mag;//first_p_info.mag;
    if (magsign > 1)
        fMagSign = 1;
    else if ( magsign < 1)
        fMagSign = -1;
    else
        fMagSign = magsign;
    
   // cout<<"CheckMergedFraction "<<fMagSign<<endl;
    
    float deta = eta2 - eta1;
    //else return false;
    if (TMath::Abs(deta) > TMath::Abs(fDEtaMax)) return true;
    //if(first_p_info.track_id==-4755 && sec_p_info.track_id==-170) cout<<"enter pass TTTTTTTTTTTTTT"<<endl;
    bool pair_pass = true;
    float badpoints = 0.;
    float allpoints = 0.;

    for (float irad = fRadiusMin; irad < fRadiusMax; irad += 0.01) {
        
        // Calculate radius:
        float rad = irad;
        
        // Calculate dPhiStar:
        double afsi0b = -0.15*abs(magsign)*chg1*fMagSign*rad/pt1;
        double afsi1b = -0.15*abs(magsign)*chg2*fMagSign*rad/pt2;
        Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
        dphistar = TVector2::Phi_mpi_pi(dphistar); // returns phi angle in the interval [-PI,PI)
        
        // Calculate distance:
        //double distance = TMath::Sqrt(rad * rad * (2 - 2 * TMath::Cos(dphistar)));
        double distance = 2 * TMath::Sin(TMath::Abs(dphistar) * 0.5) * rad;
        
        // Check if pair parameters meet the requirements:
        if (distance < fDistanceMax) {
            badpoints += 1.0;
        }
        allpoints += 1.0;
    }
     //cout<<chg1<<" "<<chg2<<" "<<badpoints<<" "<<allpoints<<endl;
    if (allpoints != 0.0) {
        // Calculate fraction:
        Double_t fraction = badpoints / allpoints;
        // Remove pair if the fraction is above limit:
        if(fraction > fMergedFractionLimit) {
            pair_pass = false;
        }
    }
    else {
        pair_pass = true;
    }
    if (pair_pass){
        pair_pass = CheckAntiGamma(first_p_info,sec_p_info);
    }
    return pair_pass;

}
bool AliAnalysisTaskDowangppDCAfit::CheckAntiGamma(p_info &first_p_info, p_info &sec_p_info){
    //AliFemtoPairCutAntiGamma.cxx
    bool temp = true;
    float chg1 = first_p_info.charge;
    float chg2 = sec_p_info.charge;
    if (chg1 * chg2 <0.) {
        float theta1 = first_p_info.kin_info.Theta();
        float theta2 = sec_p_info.kin_info.Theta();
        float dtheta = TMath::Abs(theta1 - theta2);

        float p1 = first_p_info.kin_info.P();
        float p2 = sec_p_info.kin_info.P();

        float px1 = first_p_info.kin_info.Px();
        float px2 = sec_p_info.kin_info.Px();
        float py1 = first_p_info.kin_info.Py();
        float py2 = sec_p_info.kin_info.Py();
        float pz1 = first_p_info.kin_info.Pz();
        float pz2 = sec_p_info.kin_info.Pz();

        float E1 = TMath::Sqrt(Mass_e*Mass_e + p1*p1);
        float E2 = TMath::Sqrt(Mass_e*Mass_e + p2*p2);
    
        float minv = 2.* Mass_e*Mass_e + 2.*(E1 * E2 - px1*px2 - py1*py2 - pz1*pz2);
        if ((minv < fMaxEEMinv) && (dtheta < fMaxDTheta)) temp = false;
    }
    
    // check separation at TPC entrance
    bool tempTPCEntrance = true;
    
    float tpc_x1 = first_p_info.TPC_Entrance[0][0];
    float tpc_x2 = sec_p_info.TPC_Entrance[0][0];
    
    float tpc_y1 = first_p_info.TPC_Entrance[0][1];
    float tpc_y2 = sec_p_info.TPC_Entrance[0][1];
    
    float tpc_z1 = first_p_info.TPC_Entrance[0][2];
    float tpc_z2 = sec_p_info.TPC_Entrance[0][2];
    
    float dist = sqrt(TMath::Power(tpc_x1-tpc_x2,2) + TMath::Power(tpc_y1-tpc_y2,2) + TMath::Power(tpc_z1-tpc_z2,2));
    tempTPCEntrance = dist > fDTPCMin;
    
    // check average separation
    bool avgsepCheck = true;
    float avgSep = 0.;
    int p_count = 0;
    for (int i = 0; i < 8; i++) {
        if (TpcPointIsUnset(first_p_info,i) || TpcPointIsUnset(sec_p_info,i)) break;
        avgSep += TpcPointSep(first_p_info,sec_p_info,i);
        p_count++;
    }
    
    // this maybe wrong?
    avgSep /= p_count;
    avgsepCheck = avgSep > fMinAvgsep;
    
    if (temp && tempTPCEntrance && avgsepCheck) {
        temp = CheckShareQuality(first_p_info,sec_p_info);
        return temp;
        
    }
    else return false;
}
bool AliAnalysisTaskDowangppDCAfit::CheckShareQuality(p_info &first_p_info, p_info &sec_p_info){
    // AliFemtoShareQualityPairCut.cxx
    bool share_pass = true;
    if (share_pass && (fShareFractionMax < 1.0 || fShareQualityMax < 1.0)) {
        
        int nh = 0;
        int an = 0;
        int ns = 0;
        
        unsigned int n_bits = first_p_info.fClusters.GetNbits();
        
        auto &tpc_clusters_1 = first_p_info.fClusters,
        &tpc_clusters_2 = sec_p_info.fClusters;
        
        auto &tpc_sharing_1 = first_p_info.fShared,
        &tpc_sharing_2 = sec_p_info.fShared;
        
        for (unsigned int imap = 0; imap < n_bits; imap++) {
                const bool  cluster_bit_1 = tpc_clusters_1.TestBitNumber(imap),
                            cluster_bit_2 = tpc_clusters_2.TestBitNumber(imap);
                // If both have clusters in the same row
                if (cluster_bit_1 && cluster_bit_2) {
                    // Do they share it ?
                    if (tpc_sharing_1.TestBitNumber(imap) && tpc_sharing_2.TestBitNumber(imap)){
                        an++;
                        nh+=2;
                        ns+=2;
                    }
                    // Different hits on the same padrow
                    else {
                        an--;
                        nh+=2;
                    }
                }
                else if (cluster_bit_1 || cluster_bit_2) {
                    // One track has a hit, the other does not
                    an++;
                    nh++;
                }
        }
        float hsmval = 0.;
        float hsfval = 0.;
        if (nh > 0) {
            hsmval = an*1.0/nh;
            hsfval = ns*1.0/nh;
        }
        
        
        if (fShareQualityMax < 1.0) {
            share_pass &= (hsmval < fShareQualityMax);
        }
        if (fShareFractionMax < 1.0) {
            share_pass &= (hsfval < fShareFractionMax);
        }
        
    }
    return share_pass;
    
}
bool AliAnalysisTaskDowangppDCAfit::TpcPointIsUnset(p_info & info,int i){
    
 return info.TPC_Entrance[i][0] < -9000. || info.TPC_Entrance[i][1] < -9000. || info.TPC_Entrance[i][2] < -9000.;
}

float AliAnalysisTaskDowangppDCAfit::TpcPointSep(p_info &first_p_info, p_info &sec_p_info,int i){
    
    float tpc_sepx1 = first_p_info.TPC_Entrance[i][0];
    float tpc_sepx2 = sec_p_info.TPC_Entrance[i][0];
    
    float tpc_sepy1 = first_p_info.TPC_Entrance[i][1];
    float tpc_sepy2 = sec_p_info.TPC_Entrance[i][1];
    
    float tpc_sepz1 = first_p_info.TPC_Entrance[i][2];
    float tpc_sepz2 = sec_p_info.TPC_Entrance[i][2];
    
    float sep = 0.;
    
    sep = sqrt( TMath::Power(tpc_sepx1-tpc_sepx2,2) + TMath::Power(tpc_sepy1-tpc_sepy2,2) + TMath::Power(tpc_sepz1-tpc_sepz2,2) );
    
    return sep;
}
void AliAnalysisTaskDowangppDCAfit::GetGlobalPositionAtGlobalRadiiThroughTPC(AliAODTrack *track, float bfield, float globalPositionsAtRadii[9][3]){
    // Gets the global position of the track at nine different radii in the TPC
    // params:
    //   track - the track to propagate
    //   bfield - magnetic field of event
    //   globalPositionsAtRadii - Output array of global positions in the radii and xyz
    const Float_t DEFAULT_VALUE = -9999.0;
    
    // The radii at which we get the global positions
    // IROC (OROC) from 84.1 cm to 132.1 cm (134.6 cm to 246.6 cm)
    const Float_t Rwanted[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};
    
    // Make a copy of the track to not change parameters of the track
    AliExternalTrackParam etp;
    etp.CopyFromVTrack(track);
    
    // index of global position we are filling
    //  - first we use AliExternalTrackParam, then just default value
    Int_t radius_index = 0;
    
    // loop over the array of radii
    for (; radius_index < 9; radius_index++) {
        
        // extracted radius
        const Float_t radius = Rwanted[radius_index];
        // buffer to store position
        Double_t pos_buffer[3] = {0.,0.,0.};
        
        //AliFemtoThreeVector(pos_buffer).Perp()
        
        // get the global position of the track at this radial location
        bool good = etp.GetXYZatR(radius, bfield, pos_buffer, NULL);
        float t_r = sqrt(pos_buffer[0]*pos_buffer[0] + pos_buffer[1]*pos_buffer[1]);
        // if value is not good, break loading loop
        if (!good || fabs(t_r - radius) > 0.5) {
            radius_index--; // decrement to fill current location with default value
            break;
        }
        
        // store the global position
        globalPositionsAtRadii[radius_index][0] = pos_buffer[0];
        globalPositionsAtRadii[radius_index][1] = pos_buffer[1];
        globalPositionsAtRadii[radius_index][2] = pos_buffer[2];
    }
    
    // Fill any remaining positions with the default value
    for (; radius_index < 9; radius_index++) {
        std::fill_n(globalPositionsAtRadii[radius_index], 3, DEFAULT_VALUE);
    }
}
float AliAnalysisTaskDowangppDCAfit::ReAvgDphi(p_info &first_p_info, p_info &sec_p_info){
	
    float phi1 = first_p_info.kin_info.Phi();
    float phi2 = sec_p_info.kin_info.Phi();
    float chg1 = first_p_info.charge;
    float chg2 = sec_p_info.charge;
    float pt1 = first_p_info.kin_info.Pt();
    float pt2 = sec_p_info.kin_info.Pt();
    float eta1 = first_p_info.kin_info.Eta();
    float eta2 = sec_p_info.kin_info.Eta();
    
    
    double magsign = 0.;
    
    if (!aodH) {
        return false;
    }
    else {
        AliAODEvent *fAOD;
        fAOD = aodH->GetEvent();
        magsign = fAOD->GetMagneticField();
    }
    
    //float magsign = this_event_mag;//first_p_info.mag;
    if (magsign > 1)
        fMagSign = 1;
    else if ( magsign < 1)
        fMagSign = -1;
    else
        fMagSign = magsign;
    
   // cout<<"CheckMergedFraction "<<fMagSign<<endl;

    float magval = 0.5;
    float dphiAvg = 0.;

    for(int i=0;i<9;i++){
        Double_t rad = TPCradii[i];
        // Calculate dPhiStar:
        //double afsi0b = -0.07510020733*chg1*fMagSign*rad/pt1;
        //double afsi1b = -0.07510020733*chg2*fMagSign*rad/pt2;
        double afsi0b = -0.15*magval*chg1*fMagSign*rad/pt1;
        double afsi1b = -0.15*magval*chg2*fMagSign*rad/pt2;
        Double_t dphistar =  phi2 - phi1 + TMath::ASin(afsi1b) - TMath::ASin(afsi0b);
        dphistar = TVector2::Phi_mpi_pi(dphistar); // returns phi angle in the interval [-PI,PI)
        dphiAvg += dphistar;
    }
    dphiAvg = dphiAvg/9.;
    return dphiAvg; 
}




