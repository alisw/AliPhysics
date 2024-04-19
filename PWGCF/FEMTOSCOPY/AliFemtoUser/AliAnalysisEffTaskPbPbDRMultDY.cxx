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

#include "AliAnalysisEffTaskPbPbDRMultDY.h"

ClassImp(AliAnalysisEffTaskPbPbDRMultDY)
//MDR -> 9 multiplicity in PbPb in DR task 
double fV1MDR[3];


void AliAnalysisEffTaskPbPbDRMultDY::SetFB(int fb)
{
  fFB = fb;
}

void AliAnalysisEffTaskPbPbDRMultDY::SetPidMethod(PidMethod method)
{
  fPidMethod = method;
}

int AliAnalysisEffTaskPbPbDRMultDY::GetPidMethod()
{
  return (int)fPidMethod;
}

void AliAnalysisEffTaskPbPbDRMultDY::SetPidMethod(int method)
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

AliAnalysisEffTaskPbPbDRMultDY::AliAnalysisEffTaskPbPbDRMultDY() :
  AliAnalysisTaskSE(), centrality(0), fHistoList(0), fDCAtoPrimVtx(0), fIfAliEventCuts(kFALSE), fFB(128), fPidMethod(kExclusivePIDDiffRejection), fpidResponse(0), fAODpidUtil(0), fEventCuts(0)

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
    }
  }
 for ( Int_t i = 0; i < 11; i++) {
    if(i<4) fHistEv[i] = NULL;
    fHistQA[i] = NULL;
    if(i<3) fHistQA2D[i] = NULL;
  }
}
AliAnalysisEffTaskPbPbDRMultDY::AliAnalysisEffTaskPbPbDRMultDY(TString name, int pidMethod, int filterbit) :
  AliAnalysisTaskSE(name), centrality(0), fHistoList(0), fDCAtoPrimVtx(0), fIfAliEventCuts(kFALSE), fFB(128), fPidMethod(kExclusivePIDDiffRejection), fpidResponse(0), fAODpidUtil(0), fEventCuts(0)

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
    }
  }
 for ( Int_t i = 0; i < 11; i++) {
    if(i<4) fHistEv[i] = NULL;
    fHistQA[i] = NULL;
    if(i<3) fHistQA2D[i] = NULL;
  }

 //DefineInput(0, TChain::Class());
 //DefineOutput(0, TTree::Class());
 DefineOutput(1, TList::Class());
}

//_______________________________________________________

AliAnalysisEffTaskPbPbDRMultDY::~AliAnalysisEffTaskPbPbDRMultDY()
{
  // Destructor
  if(AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis)
    delete fHistoList;
}

//_______________________________________________________

void AliAnalysisEffTaskPbPbDRMultDY::UserCreateOutputObjects()
{

  /* create output */
  fHistoList = new TList();
  fHistoList->SetOwner(kTRUE);
   
  TString hname1, hname2, hname3, hname4, hname5;
  
  TString htitle1, htitle2, htitle3, htitle4,htitle5;
  
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
      fReconstructedAfterCuts[i*PARTTYPES+j][0] = new TH2F(hname2.Data(),htitle2.Data(),50, -1.5, 1.5,100,0.,5.0);
      hname2+="Minus";htitle2+="Minus";
      fReconstructedAfterCuts[i*PARTTYPES+j][1] = new TH2F(hname2.Data(),htitle2.Data(),50, -1.5, 1.5,100,0.,5.0);

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
      fGeneratedMCPrimaries[i*PARTTYPES+j][1]->Sumw2();
      fContamination[i*PARTTYPES+j][0]->Sumw2();
      fContamination[i*PARTTYPES+j][1]->Sumw2();
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


  for(Int_t i = 0; i < 5; i++)  {
    hname = "fHistEventCutsM";
    hname+= i;
    

    fHistEvCuts[i] = new TH1F(hname,Form("Event Cuts M%d",i) , 5, 0, 5);
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(1,"All");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(2,"MultCut");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(3,"NoVertex");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(4,"PileUp");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(5,"z-vertex>10");
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
 
  fHistQA[9] = new TH1F("fHistEventCuts", "Event Cuts" , 4, 0, 5);
  fHistQA[9]->GetXaxis()->SetBinLabel(1,"All");
  fHistQA[9]->GetXaxis()->SetBinLabel(2,"NoVertex");
  fHistQA[9]->GetXaxis()->SetBinLabel(3,"PileUp");
  fHistQA[9]->GetXaxis()->SetBinLabel(4,"z-vertex>10");


  fHistQA[10] = new TH1F("fHistTrackCuts", "Track Cuts" , 7, 0.5, 7.5);
  fHistQA[10]->GetXaxis()->SetBinLabel(1,"AllTracksInEvents");
  fHistQA[10]->GetXaxis()->SetBinLabel(2,"GetTrack");
  fHistQA[10]->GetXaxis()->SetBinLabel(3,"Filter bit");
  fHistQA[10]->GetXaxis()->SetBinLabel(4,"Y");
  fHistQA[10]->GetXaxis()->SetBinLabel(5,"Pt");
  fHistQA[10]->GetXaxis()->SetBinLabel(6,"DCA");
  fHistQA[10]->GetXaxis()->SetBinLabel(7,"Electron Rejection");

  fHistQA2D[0] = new TH2F("dcaHistDcaXY","DCA XY",50, 0, 5,210, -2.1, 2.1);
  fHistQA2D[1] = new TH2F("dcaHistDcaZ","DCA Z", 50, 0, 5, 210, -2.1, 2.1);
  fHistQA2D[2] = new TH2F("fPhiY","Y-Phi",100, -2, 2, 100, -TMath::Pi(), TMath::Pi());
 

 for (Int_t i = 0; i < 11; i++){
   fHistoList->Add(fHistQA[i]);
   if(i<3) fHistoList->Add(fHistQA2D[i]);
   if(i<5) {
   for(Int_t j = 0 ; j<PARTTYPES; j++)
     for(int chg=0;chg<2;chg++){
       fHistoList->Add(fHistQAPID[i][j][chg]);
       fHistoList->Add(fHistQAPIDFail[i][j][chg]);
     }
   }
 }
 for (Int_t i = 0; i < MULTBINS*PARTTYPES; i++){
  for(Int_t chg=0;chg<2;chg++){
      fHistoList->Add(fGeneratedMCPrimaries[i][chg]);
      fHistoList->Add(fMCPrimariesThatAreReconstructed[i][chg]);
      fHistoList->Add(fMCPrimariesThatAreReconstructedNoNsigma[i][chg]);
      fHistoList->Add(fReconstructedAfterCuts[i][chg]);
      fHistoList->Add(fReconstructedNotPrimaries[i][chg]);
      fHistoList->Add(fReconstructedPrimaries[i][chg]);
      fHistoList->Add(fContamination[i][chg]);
    }
  }
//********** PID ****************


  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fpidResponse = inputHandler->GetPIDResponse();
  std::cout<<"*******"<< fpidResponse<<std::endl;  
  // ************************
  PostData(1, fHistoList);
}


//_____________________________________________________________________

bool IsPionNSigmaMDR(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime)
{
   if (mom > 0.5) {
     if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 2) return true;
   }
   else {
     if (TMath::Abs(nsigmaTPCPi) < 2) return true;
   }
   return false;
}

bool IsPionNSigma3MDR(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime)
{
    if(mom > 0.5) {
     if(TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 3) return true;
    }
    else{
     if (TMath::Abs(nsigmaTPCPi) < 3) return true;
    }
    return false;
}

bool IsKaonNSigmaMDR(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
  if (mom > 0.5) {
     //rejection of unwanted contamination
    //  if(mom>1 && TOFtime<-400) return false;
      if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < 2) return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCK) < 2) return true;
  }
  return false;
}

bool IsKaonNSigma3MDR(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime)
{
  if (mom > 0.5) {
    if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < 3) return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCK) < 3) return true;
  }
  return false;
}

bool IsProtonNSigmaMDR(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{
    if (mom > 0.5) {
   //  if(mom>1.8 && TOFtime<-300) return false;
     if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 2) return true;
    }
    else {
      if (TMath::Abs(nsigmaTPCP) < 2) return true;
    }
    return false;
}


bool IsProtonNSigma3MDR(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime)
{
  if (mom > 0.5) {
    if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 3) return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCP) < 3) return true;
  }
  return false;
}


bool IsElectronMDR(float nsigmaTPCe, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
{
  if(TMath::Abs(nsigmaTPCe)<1)
    return true;
  else
    return false;
}

//_______________________________________________________
void AliAnalysisEffTaskPbPbDRMultDY::UserExec(Option_t *)
{
    AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliAODEvent *fAOD = aodH->GetEvent();
 
 // std::cout<<GetPidMethod()<<std::endl;

  /***Get Event****/
  //AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());


  if (!aodEvent) return;
  if(aodEvent->GetNumberOfTracks() <= 0)
   return;

  fHistEvCuts[0]->Fill(0);

  if(fIfAliEventCuts){
    //******* Ali Event Cuts - applied on AOD event - standard cuts for Run2 as prepared by DPG group ************
    if (!fEventCuts->AcceptEvent(aodEvent)) {
      return;
    }
  }


  AliAODHeader *fAODheader = (AliAODHeader*)aodEvent->GetHeader();
  AliCentrality* alicent= aodEvent->GetCentrality(); //in PbPb and pPb
  AliMultSelection *MultSelection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");

  Double_t mult = 0.0;
  if(!MultSelection) {
   // std::cout << "AliMultSelection object not found!" << std::endl;
    return;
  }

  else mult = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
  if(mult <= 0 || mult >= 10000) return;

  //std::cout<<"########fcentrality#######   "<<mult<<std::endl;
  fHistEv[0]->Fill(mult);
  fHistQA[9]->Fill(1);

  //****** Multiplicity selection *********
  Int_t fcent = -999;
  if(mult >= 0 && mult <5)  fcent = 0;
  else if(mult >= 5 && mult <10) fcent = 1;
  else if(mult >= 10 && mult <20) fcent = 2;
  else if(mult >= 20 && mult <30) fcent = 3;
  else if(mult >= 30 && mult <40) fcent = 4;
  else return;

  if(fcent==0)fHistEvCuts[0]->Fill(1);
  else if(fcent==1)fHistEvCuts[1]->Fill(1);
  else if(fcent==2)fHistEvCuts[2]->Fill(1);
  else if(fcent==3)fHistEvCuts[3]->Fill(1);
  else if(fcent==4)fHistEvCuts[4]->Fill(1);

  const AliAODVertex* vertex =(AliAODVertex*) aodEvent->GetPrimaryVertex();
  vertex->GetPosition(fV1MDR);
  if (!vertex || vertex->GetNContributors()<1) return;

  fHistQA[9]->Fill(2);
  if(fcent==0)fHistEvCuts[0]->Fill(2);
  else if(fcent==1)fHistEvCuts[1]->Fill(2);
  else if(fcent==2)fHistEvCuts[2]->Fill(2);
  else if(fcent==3)fHistEvCuts[3]->Fill(2);
  else if(fcent==4)fHistEvCuts[4]->Fill(2);



//********* Pile-up removal*******************
  //check this: https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsPileup
  AliAnalysisUtils *anaUtil=new AliAnalysisUtils();

  Bool_t fpA2013 = kFALSE;
  Bool_t fMVPlp = kFALSE;
  Bool_t fOutOfBunchPlp = kFALSE;
  Bool_t fisPileUp = kFALSE;

  if(fpA2013)
    if(anaUtil->IsVertexSelected2013pA(aodEvent)==kFALSE) return;


  //Multiple vertices with tracks
  if(fMVPlp) anaUtil->SetUseMVPlpSelection(kTRUE);
  else anaUtil->SetUseMVPlpSelection(kFALSE);
  //if this (fMVPlp) is false, than rejection based on Multiple SPD vertices is used

  //out-of-bunch pile-up rejection from AnaUtil
  anaUtil->SetUseOutOfBunchPileUp(fOutOfBunchPlp);
  /*
  Int_t fMinPlpContribMV = 0; //default value: fMinPlpContribMV(5),
  Int_t fMinPlpContribSPD = 3; //default value: fMinPlpContribSPD(5),

  if(fMinPlpContribMV) anaUtil->SetMinPlpContribMV(fMinPlpContribMV);
  if(fMinPlpContribSPD) anaUtil->SetMinPlpContribSPD(fMinPlpContribSPD);
  */


  if(fisPileUp)
  if(anaUtil->IsPileUpEvent(aodEvent)) return;

  delete anaUtil;

  //pileup for LHC20e3a -> Injective Pileup over events 
  /*AliAODMCHeader *mcHeader = 0;
  mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if(!mcHeader) {
    printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
    return;
  }
  */
  fHistQA[9]->Fill(3);
  if(fcent==0)fHistEvCuts[0]->Fill(3);
  else if(fcent==1)fHistEvCuts[1]->Fill(3);
  else if(fcent==2)fHistEvCuts[2]->Fill(3);
  else if(fcent==3)fHistEvCuts[3]->Fill(3);
  else if(fcent==4)fHistEvCuts[4]->Fill(3);


  //TString vtxTtl = vertex->GetTitle();
  //if (!vtxTtl.Contains("VertexerTracks")) return;
  Float_t zvtx = vertex->GetZ();
  if (TMath::Abs(zvtx) > 10) return;

  fHistQA[0]->Fill(zvtx);
  fHistQA[9]->Fill(4);
  if(fcent==0)fHistEvCuts[0]->Fill(4);
  else if(fcent==1)fHistEvCuts[1]->Fill(4);
  else if(fcent==2)fHistEvCuts[2]->Fill(4);
  else if(fcent==3)fHistEvCuts[3]->Fill(4);
  else if(fcent==4)fHistEvCuts[4]->Fill(4);


 //**** getting MC array ******
  TClonesArray  *arrayMC;

  arrayMC = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));


  //copying pid information for FB 128
  int nofTracks=0;
  
  nofTracks = aodEvent->GetNumberOfTracks();
  const int UNDEFINED_LABEL = -1;

  // 'labels' maps a track's id to the track's index in the Event
  // i.e. labels[Event->GetTrack(x)->GetID()] == x
  std::vector<int> labels(nofTracks, UNDEFINED_LABEL);
  // looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
  for (int i=0;i<nofTracks;i++) {
  
   const AliAODTrack *aodtrack=(AliAODTrack*)aodEvent->GetTrack(i);
   if (!aodtrack->TestFilterBit(128)) {
      // Skip TPC-only tracks
      const int id = aodtrack->GetID();
      if (id < 0) {
        continue;
      }

      // Resize labels vector if "id" is larger than mapping allows
      if (static_cast<size_t>(id) >= labels.size()) {
        labels.resize(id + 1024, UNDEFINED_LABEL);
      }
      labels[id] = i;
    }
  }


//*****************Track reconstruction****************************************
 TObjArray recoParticleArray[PARTTYPES];

  fHistQA[10]->Fill(1,aodEvent->GetNumberOfTracks());
  //loop over AOD tracks

  for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++) // MC reconstructed loop
  {

    AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iTracks); 

    if (!track)
      continue;
    fHistQA[10]->Fill(2);
    
    //UInt_t filterBit = (1 << (7));
    UInt_t filterBit = fFB;
    if(filterBit && !track->TestFilterBit(filterBit))
      continue;        
    fHistQA[10]->Fill(3);
    
    Int_t charge = 0;
    if(track->Charge() > 0 ) charge=0;
    else if (track->Charge() < 0 ) charge=1; 
    
    double PionMass = 0.13957;
    double KaonMass = 0.4937;
    double ProtonMass = 0.9382720;
     
   
    fHistQA[10]->Fill(4);
    
    if(track->Pt() < 0.2 || track->Pt() > 2.5)
      continue; 
    fHistQA[10]->Fill(5);
    
   // if(track->GetTPCNcls() < 70) continue;
    
   //DCA
    Double_t DCAXY;
    Double_t DCAZ;

    
    AliAODTrack* aodtrackpid;

    //for FB 128 - tpc only tracks
    if(filterBit==(1 << (7)))
      aodtrackpid =(AliAODTrack*)aodEvent->GetTrack(labels[-1-aodEvent->GetTrack(iTracks)->GetID()]);
    else
      aodtrackpid = track;
    
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
     aodtrackpid->GetXYZ(pos);

     Double_t DCAX = pos[0] - vertexX;
     Double_t DCAY = pos[1] - vertexY;
     DCAZ = pos[2] - vertexZ;
     DCAXY = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));
     fHistQA[10]->Fill(6);
        
    //Electron rejection
    float nSigmaTPCPi = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kPion);
    float nSigmaTPCK = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kKaon);
    float nSigmaTPCP = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kProton);
    float nSigmaTPCe = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kElectron);
    if(IsElectronMDR(nSigmaTPCe,nSigmaTPCPi,nSigmaTPCK,nSigmaTPCP))
      continue;
   
    fHistQA[10]->Fill(7);
    
    fHistQA[1]->Fill(track->GetTPCClusterInfo(2,1)); 
    fHistQA[2]->Fill(track->GetTPCNclsF());
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
    bool isProtonNsigma  =0 ;

    if(fPidMethod==kNSigma){
    //******** With double counting *******************
      isPionNsigma = (IsPionNSigmaMDR(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]));
      isKaonNsigma = (IsKaonNSigmaMDR(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]));
      isProtonNsigma = (IsProtonNSigmaMDR(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
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
       isPionNsigma = (IsPionNSigmaMDR(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]));
      }
      else if(nSigmaPIDK<nSigmaPIDPi && nSigmaPIDK<nSigmaPIDP){
       isKaonNsigma = (IsKaonNSigmaMDR(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]));
      }
      else if(nSigmaPIDP<nSigmaPIDPi && nSigmaPIDP<nSigmaPIDK){
       isProtonNsigma = (IsProtonNSigmaMDR(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      }
    }
    else if(fPidMethod==kExclusivePID){
      //******** Exclusive PID ********************
      isPionNsigma = (IsPionNSigmaMDR(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]) && !IsKaonNSigmaMDR(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigmaMDR(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isKaonNsigma = (!IsPionNSigmaMDR(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && IsKaonNSigmaMDR(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigmaMDR(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isProtonNsigma = (!IsPionNSigmaMDR(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && !IsKaonNSigmaMDR(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && IsProtonNSigmaMDR(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
    }
    else if(fPidMethod==kExclusivePIDDiffRejection){
      //******** Exclusive PID, different rejection  ********************
      isPionNsigma = (IsPionNSigmaMDR(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2]) && !IsKaonNSigma3MDR(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigma3MDR(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isKaonNsigma = (!IsPionNSigma3MDR(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && IsKaonNSigmaMDR(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && !IsProtonNSigma3MDR(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
      isProtonNsigma = (!IsPionNSigma3MDR(track->Pt(),nSigmaTPCPi, nSigmaTOFPi, tTofSig-pidTime[2])  && !IsKaonNSigma3MDR(track->Pt(),nSigmaTPCK, nSigmaTOFK, tTofSig-pidTime[3]) && IsProtonNSigmaMDR(track->Pt(),nSigmaTPCP, nSigmaTOFP, tTofSig-pidTime[4]));
    }
    if (isPionNsigma){
      fHistQAPID[0][1][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][1][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
      fHistQAPID[2][1][charge]->Fill(tPt,nSigmaTOFPi);
      fHistQAPID[3][1][charge]->Fill(tPt,nSigmaTPCPi);
      fHistQAPID[4][1][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);
    }
    else{
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
    else{
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
      if(track->Y(PionMass) < -0.5 || track->Y(PionMass) > 0.5) continue; 
      if (track->Pt() > 0.2 && track->Pt() < 2.5)
        fReconstructedAfterCuts[PARTTYPES*fcent+1][charge]->Fill(track->Y(PionMass), track->Pt());
      if (!MCtrk) continue;
      recoParticleArray[1].Add(MCtrk);
     }
     //Fills for all identified pions found after cuts (reconstructed) - numerator for Efficiency

     //********* PID - kaons ********
     if (isKaonNsigma){
       if(track->Y(KaonMass) < -0.5 || track->Y(KaonMass) > 0.5) continue; 
       if (track->Pt() > 0.5 && track->Pt() < 2.5)
         fReconstructedAfterCuts[PARTTYPES*fcent+2][charge]->Fill(track->Y(KaonMass), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[2].Add(MCtrk);
     }
     //Fills for all identified kaons found after cuts (reconstructed) - numerator for Efficiency

    //********* PID - protons ********
     if (isProtonNsigma){
       if(track->Y(ProtonMass) < -0.5 || track->Y(ProtonMass) > 0.5) continue;
       if (track->Pt() > 0.5 && track->Pt() < 2.5)
         fReconstructedAfterCuts[PARTTYPES*fcent+3][charge]->Fill(track->Y(ProtonMass), track->Pt());
       if (!MCtrk) continue;
       recoParticleArray[3].Add(MCtrk);
     }

     //Fills for all identified protos found after cuts (reconstructed) - numerator for Efficiency
     //******************************

     //get coresponding MC particle 
       
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
    if(abs(PDGcode)==211){
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
    else if(abs(PDGcode)==321){
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
    else if(abs(PDGcode) == 2212){
       if(isPionNsigma)
         fMisidentification[fcent][charge]-> Fill(3,0.5);
       if(isKaonNsigma)
         fMisidentification[fcent][charge]-> Fill(3,1.5);
       if(isProtonNsigma)
         fMisidentification[fcent][charge]-> Fill(3,2.5);
       if(!isPionNsigma && !isKaonNsigma && !isProtonNsigma)
         if(status)
           fMisidentification[fcent][charge]-> Fill(3,3.5);
    }
    fContamination[PARTTYPES*fcent][charge]-> Fill(PDGcode,track->Pt());
    //Contaminations: "how many pions are in the kaons sample"? etc.
    //Do not use for corrections: using those values will be dependant on i.e. Pi/K ratio in MC
    //Use misidentification fraction instead

    if(isPionNsigma){
       fContamination[PARTTYPES*fcent+1][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for pions
    }
    if(isKaonNsigma){
       fContamination[PARTTYPES*fcent+2][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for kaons
    }
    if(isProtonNsigma){
       fContamination[PARTTYPES*fcent+3][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for protons
    }

  }


  // MONTECARLO PARTICLES 
  if(!arrayMC){
    AliError("Array of MC particles not found");
    return;
  }
  // loop over MC stack
  for (Int_t ipart = 0; ipart < arrayMC->GetEntriesFast(); ipart++) {
  AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(ipart);

  if (!MCtrk) continue;
  //std::cout<<"particle obtained"<<std::endl;
    
  Int_t PDGcode = TMath::Abs(MCtrk->GetPdgCode()); 
  Int_t charge=0;
  if(MCtrk->Charge() < 0) charge=1;
  else if(MCtrk->Charge() > 0) charge=0;


  if(MCtrk->Y() < -0.5 || MCtrk->Y() > 0.5) continue; 
 
  if(MCtrk->GetPdgCode() == 211){
    if (MCtrk->Pt() < 0.2 || MCtrk->Pt() > 2.5) continue;
  }
  if(MCtrk->GetPdgCode() == 321){
    if (MCtrk->Pt() < 0.5 || MCtrk->Pt() > 2.5) continue;
  }
  if(MCtrk->GetPdgCode() == 2212){
    if (MCtrk->Pt() < 0.5 || MCtrk->Pt() > 2.5) continue;
  }
  // check physical primary 
  if(MCtrk->IsPhysicalPrimary()){
    // Filling histograms for MC truth particles
    fGeneratedMCPrimaries[fcent*PARTTYPES][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
    if(PDGcode==211)
      fGeneratedMCPrimaries[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
    else if(PDGcode==321)
      fGeneratedMCPrimaries[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
    else if(PDGcode==2212)
      fGeneratedMCPrimaries[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Y(), MCtrk->Pt());

    //Filling data from MC truth particles only for particles that were reconstruced
    if (recoParticleArray[0].Contains(MCtrk)){ //All
      fMCPrimariesThatAreReconstructed[fcent*PARTTYPES][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
      fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
      if(PDGcode==211)
       fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
      if(PDGcode==321)
       fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
      if(PDGcode==2212)
        fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
      }
      if (recoParticleArray[1].Contains(MCtrk)){ //Pions
       if(PDGcode==211){
         fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
       }
      }
      if (recoParticleArray[2].Contains(MCtrk)){ //Kaons
       if(PDGcode==321){
           fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
       }
      }
      if (recoParticleArray[3].Contains(MCtrk)){ //Protons
       if(PDGcode==2212){
           fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Y(), MCtrk->Pt());
       }
      }
    }//end of check primary tracks loop 
  } //end of mc truth loop
  PostData(1, fHistoList);
}
//-----------------------------------------------------------------




