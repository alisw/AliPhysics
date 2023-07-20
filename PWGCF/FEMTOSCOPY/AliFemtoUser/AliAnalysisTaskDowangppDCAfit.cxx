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
       //if(ic > 0) continue;
        for(int ichg=0;ichg<2;ichg++){
            MomSmearing[ichg][ic] = nullptr;
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

    EventDis =new TH1F("fHistEventDis"," ",100,0,100);
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

    
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    fAODpidUtil = aodH->GetAODpidUtil();

    // Post output data.
    PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnalysisTaskDowangppDCAfit::UserExec(Option_t *) 
{
    float CentMax = 50.;
    float CentMin = 10.;
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
    Analyze(fAOD,v0Centr);

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
void AliAnalysisTaskDowangppDCAfit::Analyze(AliAODEvent* aod,Float_t v0Centr)
{


    int centrCode = -10;
    centrCode = ReCentrCode(v0Centr);
    int iCent = centrCode;

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
    int IsThisEventHaveMCP[8] = {0};
    int IsThisEventHaveRealP[8] = {0};
//cout<<"howmany "<<nMCtracks<<endl;

    vector<TLorentzVector> RecoTrackLabelV[2];
    vector<TLorentzVector> TureTrackLabelV[2];
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
        
        if(ChargeLabel==-999 || ThisTrackPDGCode==-999) continue;
  // cout<<"dddd "<<endl;  
      float tmp_DCAxyCut = 0.0105 + 0.0350 * TMath::Power(trackPt,-1.1);

        if(abs(DCAxy) < tmp_DCAxyCut){
            pTdisReco[ChargeLabel][iCent]->Fill(track->Pt());
        }
        //vector<int> RecoTrackLabelV[2];
        //vector<int> TureTrackLabelV[2];
        const Int_t label = TMath::Abs(track->GetLabel());
        // AliAODMCParticle *mcTrack = (AliAODMCParticle *)mcEvent->GetTrack(iMCtrack);   
        //AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCarray->At(label));
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCarray->At(label));
        Int_t pdgCode = 0;
        if (mcTrack){

            Int_t MClabel = mcTrack->Label();

            if(abs(DCAxy) < tmp_DCAxyCut){
                TLorentzVector ParticleReco;
                ParticleReco.SetPtEtaPhiM(track->Pt(),track->Eta(),track->Phi(),MassP);
                TLorentzVector ParticleTure;
                ParticleTure.SetPtEtaPhiM(mcTrack->Pt(),mcTrack->Eta(),mcTrack->Phi(),MassP);
                RecoTrackLabelV[ChargeLabel].push_back(ParticleReco);
                TureTrackLabelV[ChargeLabel].push_back(ParticleTure);
            }
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
        if(primaryFlag > 0){
            //cout<<"here "<<primaryFlag<<endl;
            DCAxy_pDetailFrac[ChargeLabel][primaryFlag-1][iCent]->Fill(DCAxy,trackPt);
            if(primaryFlag==1 && abs(DCAxy) < tmp_DCAxyCut && pdgCode==ThisTrackPDGCode){
                pTdisRecoWhenPrimary[ChargeLabel][iCent]->Fill(trackPt);
            }
            if(abs(DCAxy) < tmp_DCAxyCut && pdgCode==ThisTrackPDGCode){
                pTdisRecoIDAsP[ChargeLabel][iCent]->Fill(trackPt);
            }
        }

    }// end event loop

    for(int ichg=0;ichg<2;ichg++){
        for(int ip1=0;ip1<RecoTrackLabelV[ichg].size();ip1++){
            TLorentzVector Particle1_Reco = RecoTrackLabelV[ichg][ip1];
            TLorentzVector Particle1_True = TureTrackLabelV[ichg][ip1];

            for(int ip2=ip1+1;ip2<RecoTrackLabelV[ichg].size();ip2++){
                TLorentzVector Particle2_Reco = RecoTrackLabelV[ichg][ip2];
                TLorentzVector Particle2_True = TureTrackLabelV[ichg][ip2];

                float kStar_Reco = re_kstar(Particle1_Reco,Particle2_Reco);
                float kStar_True = re_kstar(Particle1_True,Particle2_True);
                //cout<<"ddd "<<kStar_Reco<<" "<<kStar_True<<endl;
                MomSmearing[ichg][iCent]->Fill(kStar_True,kStar_Reco);
            }
        }
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

