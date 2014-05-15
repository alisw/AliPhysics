/*
  Int_t year=2012;
  Int_t period=5; // 0=a, 1=b, 2=c
  .x $HOME/rootlogon.C   
  .L ~/alice/reconstruction/trunk/QAplots/source/highPt/AnalyzedEdx.C+ 
  //SelectV0sforPID();
  MakeSelectedArrays(4,4000);
  FitSelected(10000,100000,year,period);

 */

/*
  To implement:
  1. Selection of runs:
     a.)  Redo run wise gain calibration. Important for the bad periods 
     b.)  Remove runs with outlier peak position.  
     c.)  Make sector vise calibration. (To eliminate trigger bias)

   

*/

#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH2.h"
#include "TF1.h"
#include "TTreeStream.h"
#include "AliMathBase.h"
#include "TSystem.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TRandom.h"
#include "AliTPCcalibBase.h"
//
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "TMath.h"
#include "AliXRDPROOFtoolkit.h"


TTree * tree  = 0;
Int_t run=0;
Int_t period=0;
TTreeSRedirector *pcstream = 0; //new TTreeSRedirector("trend.root");
TObjArray * arrayPions=0;
TObjArray * arrayProtons=0;
TObjArray * arrayElectrons=0;
TMatrixD *ppullPionsV0=0;
TMatrixD *ppullProtonsV0=0;
TMatrixD *ppullElectronsV0=0;



//void SetPionAliases();
//void FitSlopeROC();
//void AnalyzedEdxFile(const char * fname="");



void SelectV0sforPID( const char * finput="highptAll.list"){
  //
  // Code to select identified V0 for the PID 
  // As an input chain of filter trees is used
  // Parameter:
  //   finput - name of the list file
  // Oputput:
  //   file - V0Selected.root
  TChain * chain  = AliXRDPROOFtoolkit::MakeChainRandom(finput,"V0s",0,1000);
  chain->SetCacheSize(1000000000);
  //
  TDatabasePDG pdg;
  Double_t massLambda = pdg.GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg.GetParticle("K0")->Mass();
  Double_t massPion = pdg.GetParticle("pi+")->Mass();
  Double_t massProton = pdg.GetParticle("proton")->Mass();
  //
  //
  chain->SetAlias("massPion",Form("(%f+0)",massPion));
  chain->SetAlias("massProton",Form("(%f+0)",massProton));
  chain->SetAlias("massK0",Form("(%f+0)",massK0));
  chain->SetAlias("massLambda",Form("(%f+0)",massLambda));
  // delta of mass
  chain->SetAlias("K0Delta","(v0.GetEffMass(2,2)-massK0)");
  chain->SetAlias("LDelta","(v0.GetEffMass(4,2)-massLambda)");
  chain->SetAlias("ALDelta","(v0.GetEffMass(2,4)-massLambda)");
  chain->SetAlias("EDelta","(v0.GetEffMass(0,0))");
  // pull of the mass
  chain->SetAlias("K0Pull","(v0.GetEffMass(2,2)-massK0)/v0.GetKFInfo(2,2,1)");
  chain->SetAlias("LPull","(v0.GetEffMass(4,2)-massLambda)/v0.GetKFInfo(4,2,1)");
  chain->SetAlias("ALPull","(v0.GetEffMass(2,4)-massLambda)/v0.GetKFInfo(2,4,1)");
  chain->SetAlias("EPull","EDelta/v0.GetKFInfo(0,0,1)");
  // effective pull of the mass - (empirical values form fits)
  chain->SetAlias("K0PullEff","K0Delta/sqrt((3.63321e-03)**2+(5.68795e-04*v0.Pt())**2)");
  chain->SetAlias("LPullEff","LDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  chain->SetAlias("ALPullEff","ALDelta/sqrt((1.5e-03)**2+(1.8e-04*v0.Pt())**2)");
  chain->SetAlias("EPullEff","v0.GetEffMass(0,0)/sqrt((5e-03)**2+(1.e-04*v0.Pt())**2)");
  //
  //
  chain->SetAlias("dEdx0DProton","AliMathBase::BetheBlochAleph(track0.fIp.P()/massProton)");
  chain->SetAlias("dEdx1DProton","AliMathBase::BetheBlochAleph(track1.fIp.P()/massProton)");
  chain->SetAlias("dEdx0DPion","AliMathBase::BetheBlochAleph(track0.fIp.P()/massPion)");
  chain->SetAlias("dEdx1DPion","AliMathBase::BetheBlochAleph(track1.fIp.P()/massPion)");
  //
  // V0 - cuts -PID, 
  //	
  chain->SetAlias("cutDist","sqrt((track0.fIp.fP[0]-track1.fIp.fP[0])**2+(track0.fIp.fP[1]-track1.fIp.fP[1])**2)>3");
  chain->SetAlias("cutLong","track0.GetTPCClusterInfo(3,1,0)-5*abs(track0.fP[4])>130&&track1.GetTPCClusterInfo(3,1,0)>130-5*abs(track0.fP[4])");
  chain->SetAlias("cutPID","track0.fTPCsignal>0&&track1.fTPCsignal>0");
  chain->SetAlias("cutResol","sqrt(track0.fC[14]/track0.fP[4])<0.15&&sqrt(track1.fC[14]/track1.fP[4])<0.15");
  chain->SetAlias("cutV0","cutPID&&cutDist&&cutLong&&cutResol");	
  //
  //  
  chain->SetAlias("K0Selected",      "abs(K0Pull)<3. &&abs(K0PullEff)<3.  && abs(LPull)>3  && abs(ALPull)>3  &&v0.PtArmV0()>0.11"); 
  chain->SetAlias("LambdaSelected",  "abs(LPull)<3.  &&abs(LPullEff)<3.   && abs(K0Pull)>3 && abs(EPull)>3 && abs(EDelta)>0.05");  
  chain->SetAlias("ALambdaSelected", "abs(ALPull)<3. &&abs(ALPullEff)<3   && abs(K0Pull)>3 && abs(EPull)>3 &&abs(EDelta)>0.05");
  //
  chain->SetAlias("GammaSelected", "abs(EPull)<3     && abs(K0Pull)>3 && abs(LPull)>3 && abs(ALPull)>3");
  //
  // Gamma PID selection
  //
  

  //
  //
  TFile *fselected = TFile::Open("V0Selected.root","recreate");
  TTree * treeK0     =   chain->CopyTree("type==8&&cutV0&&K0Selected");
  TTree * treeLambda =   chain->CopyTree("type==4&&cutV0&&LambdaSelected");
  TTree * treeALambda =   chain->CopyTree("type==2&&cutV0&&ALambdaSelected");
  TTree * treeGamma =   chain->CopyTree("type==1&&cutV0&&GammaSelected");
  //
  TTree * trees[4]={treeK0,treeLambda, treeGamma,treeALambda};
  TList * aliases = chain->GetListOfAliases();
  Int_t nalias= aliases->GetEntries();

  for (Int_t i=0; i<4; i++){
    for (Int_t ialias=0; ialias<nalias; ialias++){
      TNamed *alias = (TNamed*)aliases->At(ialias);
      trees[i]->SetAlias(alias->GetName(),alias->GetTitle());
    }
  }  
  treeK0->Write("treeK0");
  treeLambda->Write("treeLambda");
  treeALambda->Write("treeALambda");
  treeGamma->Write("treeGamma");
  fselected->Close();
  //
}



Double_t GetPullMass(AliESDv0 * v0, Int_t p0, Int_t p1, Int_t pdgCode){
  //
  // reeturn mass pull
  //  Test values: p0=2; p1=2; pdgCode=321
  //
  TDatabasePDG *pdg= TDatabasePDG::Instance();
  Double_t pdgMass = pdg->GetParticle(pdgCode)->Mass();
  Double_t recMass = v0->GetKFInfo(p0,p1, 0);
  Double_t rmsMass = v0->GetKFInfo(p0,p1, 1);
  if (rmsMass<=0) return -1;
  return   (recMass-pdgMass)/rmsMass;
}



void MakeSelectedArrays(Int_t scaling=20,  const Int_t maxEntries=10000){
  //
  //  Select pions, proton and electron based on the V0 information
  //  Selected tracks are written to the array in parallel with the V0 pull information
  // 
  //
  

  TFile *fselected = TFile::Open("V0Selected.root");
  TTree * treeK0= (TTree*)fselected->Get("treeK0");  
  TTree * treeLambda= (TTree*)fselected->Get("treeLambda");  
  TTree * treeALambda= (TTree*)fselected->Get("treeALambda");  
  TTree * treeGamma= (TTree*)fselected->Get("treeGamma");  
  Int_t entriesK0= treeK0->GetEntries();  
  Int_t entriesLambda= treeLambda->GetEntries();  
  Int_t entriesALambda= treeALambda->GetEntries();  
  Int_t entriesGamma= treeGamma->GetEntries();  
  
  //overestimated size 
  arrayPions=new TObjArray(entriesK0*2+entriesLambda+entriesALambda);
  arrayProtons=new TObjArray(entriesLambda+entriesALambda);
  arrayElectrons=new TObjArray(entriesGamma*2);
  //
  ppullElectronsV0=new TMatrixD(entriesGamma*2,2); 
  ppullPionsV0=new TMatrixD(entriesK0*2+entriesLambda+entriesALambda,2);
  ppullProtonsV0=new TMatrixD(entriesLambda+entriesALambda,2);
  TMatrixD & pullElectronsV0=*ppullElectronsV0;
  TMatrixD & pullPionsV0=*ppullPionsV0;
  TMatrixD & pullProtonsV0=*ppullProtonsV0;
  //
  // Use K0s, Lambda and ALambda to select pions
  //
  AliESDtrack * track0=0, *track1=0;
  AliESDv0 * v0=0;
  treeK0->SetBranchAddress("track0.",&track0);
  treeK0->SetBranchAddress("track1.",&track1);
  treeK0->SetBranchAddress("v0.",&v0);
  //
  treeLambda->SetBranchAddress("track0.",&track0);
  treeLambda->SetBranchAddress("track1.",&track1);
  treeLambda->SetBranchAddress("v0.",&v0);
  //
  treeALambda->SetBranchAddress("track0.",&track1);
  treeALambda->SetBranchAddress("track1.",&track0);
  treeALambda->SetBranchAddress("v0.",&v0);
  //
  treeGamma->SetBranchAddress("track0.",&track1);
  treeGamma->SetBranchAddress("track1.",&track0);
  treeGamma->SetBranchAddress("v0.",&v0);
  //
  Int_t counterPions=0; 
  Int_t counterProtons=0;
  Int_t counterElectrons=0;
 
  TH2F *hisPionPt = new TH2F("hisPionPt","hisPionPt", 20,0.2,20,20,-1,1);
  TH2F *hisProtonPt = new TH2F("hisProtonPt","hisProtonPt", 20,0.2,20,20,-1,1);
  TH2F *hisElectronPt = new TH2F("hisElectronPt","hisElectronPt", 20,0.2,20,20,-1,1);
  //
  AliTPCcalibBase::BinLogX(hisPionPt->GetXaxis());
  AliTPCcalibBase::BinLogX(hisProtonPt->GetXaxis());
  AliTPCcalibBase::BinLogX(hisElectronPt->GetXaxis());
  //
  //
  // Filter K0s
  //
  //
  for (Int_t iv0=0; iv0<entriesK0/scaling; iv0++){
    treeK0->GetEntry(iv0);
    Double_t pullK0= GetPullMass(v0,2,2,kK0);
    // 
    hisPionPt->Fill(track0->Pt(),track0->GetTgl());
    Int_t bin =  hisPionPt->FindBin(track0->Pt(),track0->GetTgl());
    Int_t entriesPt=(bin>0)? hisPionPt->GetBinContent(bin):0;
    if (entriesPt< maxEntries) {
      pullPionsV0(counterPions,0) =pullK0;
      arrayPions->AddAt(track0->Clone(),counterPions); counterPions++;	
    }
    
    hisPionPt->Fill(track1->Pt(), track1->GetTgl());
    bin =  hisPionPt->FindBin(track1->Pt(),track1->GetTgl());
    entriesPt=(bin>0)? hisPionPt->GetBinContent(bin):0;
    if (entriesPt< maxEntries) {
      pullPionsV0(counterPions,0) =pullK0;
      arrayPions->AddAt(track1->Clone(),counterPions);
      counterPions++;
    }    
    if (iv0%200==0) printf("K0s %d\t%d\n",iv0, counterPions);
  }  
  //
  // Filter lambdas
  //
  for (Int_t ilambda=0; ilambda<=1; ilambda++){
    TTree * treeL = (ilambda==0)? treeLambda: treeALambda;
    Int_t entriesTree=treeL->GetEntries();
    for (Int_t iv0=0; iv0<entriesTree/scaling; iv0++){
      treeL->GetEntry(iv0);
      Double_t pullLambda= 0;
      if (ilambda==0) pullLambda=GetPullMass(v0,4,2,kLambda0);
      if (ilambda==0) pullLambda=GetPullMass(v0,2,4,kLambda0);
      // 1 Additional cut on Gammas
      if (TMath::Abs(v0->GetEffMass(0,0))<0.05) continue;
      //
      hisPionPt->Fill(track1->Pt(), track1->GetTgl());
      Int_t bin =  hisPionPt->FindBin(track1->Pt(),track1->GetTgl());
      Int_t entriesPt=(bin>0)? hisPionPt->GetBinContent(bin):0;
      if (entriesPt< maxEntries) {
	arrayPions->AddAt(track1->Clone(),counterPions);
	pullPionsV0(counterPions,0) =pullLambda;
	counterPions++;
      }
      //
      hisProtonPt->Fill(track0->Pt(), track0->GetTgl());
      bin =  hisProtonPt->FindBin(track0->Pt(),track0->GetTgl());
      entriesPt=(bin>0)? hisProtonPt->GetBinContent(bin):0;
      if (entriesPt< maxEntries) {
	arrayProtons->AddAt(track0->Clone(),counterProtons);
	pullProtonsV0(counterProtons,0) =pullLambda;
	counterProtons++;
      }
      if (iv0%200==0) printf("Lambda Protons %d\t%d\t%d\n",ilambda, iv0,counterProtons);
    }
  }
  //
  // Electrons selection
  //

  Int_t nel= treeGamma->Draw("track0.fTPCsignal","abs(track0.fTPCsignal/track1.fTPCsignal-1)<0.05","goff",10000);
  Double_t meanElectron0=0, sigmaElectron0=0, meanElectron1 =0 , sigmaElectron1=0;
  AliMathBase::EvaluateUni(nel, treeGamma->GetV1(),meanElectron0,sigmaElectron0,0.55*nel);
  nel= treeGamma->Draw("track1.fTPCsignal",Form("abs(track0.fTPCsignal-%f)<%f&&abs(track0.fTPCsignal/track1.fTPCsignal-1)<0.15",meanElectron0,sigmaElectron0),"goff",10000);
  AliMathBase::EvaluateUni(nel, treeGamma->GetV1(),meanElectron1,sigmaElectron1,0.55*nel);
  //
  for (Int_t iv0=0; iv0<entriesGamma/scaling; iv0++){    
    treeGamma->GetEntry(iv0); 
    Double_t pullGamma= GetPullMass(v0,0,0,kGamma);
    //
    //
    hisElectronPt->Fill(track0->Pt(), track0->GetTgl());
    Int_t bin =  hisElectronPt->FindBin(track0->Pt(),track0->GetTgl());
    Int_t entriesPt=(bin>0)? hisElectronPt->GetBinContent(bin):0;
    if (entriesPt< maxEntries && TMath::Abs(track1->GetTPCsignal()-meanElectron1)<sigmaElectron1) {
      arrayElectrons->AddAt(track0->Clone(),counterElectrons);
      pullElectronsV0(counterElectrons,0) =pullGamma;
      counterElectrons++;
    }
    hisElectronPt->Fill(track1->Pt(), track1->GetTgl());
    bin =  hisElectronPt->FindBin(track1->Pt(),track1->GetTgl());
    entriesPt=(bin>0)? hisElectronPt->GetBinContent(bin):0;
    if (entriesPt< maxEntries  && TMath::Abs(track0->GetTPCsignal()-meanElectron1)<sigmaElectron1) {
      arrayElectrons->AddAt(track1->Clone(),counterElectrons);
      pullElectronsV0(counterElectrons,0) =pullGamma;
      counterElectrons++;      
    }
    if (iv0%100==0) printf("Electrons %d\t%d\n",iv0, counterElectrons);
  }



}

void SelectHPT(){
  //
  Int_t counter=0;
  AliESDtrack * track0=0;
  TChain * chainHPT  = AliXRDPROOFtoolkit::MakeChainRandom("highptAll.list","highPt",0,1000);
  Int_t entriesHPT = chainHPT->GetEntries();
  chainHPT->SetCacheSize(1000000000);
  chainHPT->SetBranchAddress("esdTrack.",&track0);  
  TBranch * branchFlags = chainHPT->GetBranch("esdTrack.fFlags");
  //  TBranch * branchTOFSignal = chainHPT->GetBranch("esdTrack.fITSsignal");
  TBranch * branchTOFSignalDX = chainHPT->GetBranch("esdTrack.fTOFsignalDx");
  TBranch * branchParam = chainHPT->GetBranch("esdTrack.AliExternalTrackParam.fP[5]");
  TBranch * branchTrack = chainHPT->GetBranch("esdTrack.");
  counter=0;

  for (Int_t i=0; i<entriesHPT; i++){    
    //
    if (i%100000==0) printf("Entrye\t %d\t%f\n",i,counter);
    branchFlags->GetEntry(i);
    if (track0->IsOn(0x4)==0) continue;    
    //
    branchTOFSignalDX->GetEntry(i);
    if (track0->GetTOFsignalDx()>3) continue;
    //
    branchParam->GetEntry(i);
    if (TMath::Abs(track0->GetParameter()[4])<0.25) continue; 
    branchTrack->GetEntry(i);
    // Select pions
    //esdTrack.fITSncls>4&&abs(esdTrack.fTOFsignalDx**2+esdTrack.fTOFsignalDz**2)<3&&esdTrack.fTOFsignal<esdTrack.fTrackTime[3]-400&&abs(esdTrack.fTOFsignal-esdTrack.fTrackTime[2])<300"
    //
    // Select Protons
    // chainHPT->Draw("esdTrack.fTPCsignal:esdTrack.P()>>his(60,0.3,6,50,30,200)","esdTrack.fITSncls>4&&abs(esdTrack.fTOFsignalDx**2+esdTrack.fTOFsignalDz**2)<4*(1+abs(esdTrack.fP[4]))&&esdTrack.fTOFsignal>esdTrack.fTrackTime[3]+500&&abs(esdTrack.fTOFsignal-esdTrack.fTrackTime[4])<4*(100+abs(esdTrack.fP[4])*100)","colz",10000000);
    counter++;
    if (counter%100==0) printf("%d\t%d\t%f\n",i, counter,track0->GetTOFsignal());
  }
}







void FitSelected( const Int_t ngener=1000, const Int_t maxEntries=10000, Int_t year=0, Int_t period=0){
  //
  // 1. fit the roubust means and the 
  //
  Int_t counter=0;
  TVectorD vecPt0(maxEntries);
  TVectorD vecP(maxEntries);
  TVectorD vecTheta(maxEntries);
  TVectorD vecdEdxTheor(maxEntries);
  //
  TVectorD vecdEdx0(maxEntries);
  TVectorD vecdEdx1(maxEntries);
  TVectorD vecdEdx2(maxEntries);
  TVectorD vecdEdxTRD(maxEntries);
  TVectorD vecdEdxOROC(maxEntries);
  TVectorD vecdEdxAll(maxEntries);
  TVectorD * vecArray[6]={&vecdEdx0, &vecdEdx1, &vecdEdx2, &vecdEdxTRD, &vecdEdxOROC, &vecdEdxAll};
  TVectorD vecMeanTrunc60(6);
  TVectorD vecRMSTrunc60(6);
  TVectorD vecMeanTrunc80(6);
  TVectorD vecRMSTrunc80(6);
  TVectorD vecMeanTrunc(6);
  TVectorD vecRMSTrunc(6);
  TVectorD vecMeanTruncAll(6);
  TVectorD vecRMSTruncAll(6);
  TVectorD vecSelected(6);
  //
  //
  //
  TDatabasePDG pdg;
  Double_t massPion = pdg.GetParticle("pi+")->Mass();
  Double_t massProton = pdg.GetParticle("proton")->Mass();
  Double_t massElectron = pdg.GetParticle("e+")->Mass();

  TTreeSRedirector * pcstream = new TTreeSRedirector("pidDEDX.root","recreate");
  //
  //
  //
  for (Int_t isel=0; isel<ngener; isel++){
    Double_t mcp=gRandom->Rndm()*20;
    Double_t mctheta=2*(gRandom->Rndm()-0.5);
    Int_t ipid= TMath::Nint(gRandom->Rndm()*3.);
    Double_t mass=massPion;
    TObjArray *arrayTrack = arrayPions;
    if (ipid==1) {arrayTrack=arrayProtons; mass=massProton;}
    if (ipid==2) {arrayTrack=arrayElectrons; mass=massElectron;}
    Int_t entriesArray= arrayTrack->GetEntries();
    //
    Int_t nselected=0;
    Int_t nselectedArray[6]={0};
    TVectorD vselectedArray(6);
    //
    // 1. First loop - select the tracks which fullfill the pt and theta criteria
    //
    for (Int_t itrack=0; itrack<entriesArray; itrack++){
      AliESDtrack * track = (AliESDtrack*)arrayTrack->At(itrack);
      if (!track) continue;
      if (!track->GetInnerParam()) continue;
      if (TMath::Abs(track->P()/mcp-1)>0.15) continue; // pt selection
      if (TMath::Abs(track->GetTgl()-mctheta)>0.15) continue; // theta selection
      Double_t mom = track->GetInnerParam()->P();
      Double_t dEdxTheor=50*AliMathBase::BetheBlochAleph(mom/mass);
      //
      vecPt0[nselected]=track->Pt();
      vecP[nselected]=mom;
      vecTheta[nselected]=track->GetTgl();
      vecdEdxTheor[nselected]=dEdxTheor;
      //
      if (track->GetTPCdEdxInfo()){
	vecdEdx0[nselectedArray[0]++]=track->GetTPCdEdxInfo()->GetSignal(0)/dEdxTheor;      
	vecdEdx1[nselectedArray[1]++]=track->GetTPCdEdxInfo()->GetSignal(1)/dEdxTheor;      
	vecdEdx2[nselectedArray[2]++]=track->GetTPCdEdxInfo()->GetSignal(2)/dEdxTheor;      
	vecdEdxOROC[nselectedArray[4]++]=track->GetTPCdEdxInfo()->GetSignal(3)/dEdxTheor;      
      }
      if (track->GetTRDsignal()>0 && track->GetTRDncls()>80){
	vecdEdxTRD[nselectedArray[3]++]=50.*track->GetTRDsignal()/dEdxTheor;            
      }
      vecdEdxAll[nselectedArray[5]++]=track->GetTPCsignal()/dEdxTheor;            
      nselected++;
    }    
    for (Int_t idedx=0; idedx<6; idedx++) vselectedArray[idedx]=nselectedArray[idedx];
    if (nselected <20) continue;
    //
    Double_t meanPt=TMath::Mean(nselected,vecPt0.GetMatrixArray());
    Double_t meanP=TMath::Mean(nselected,vecP.GetMatrixArray());
    Double_t meanTheta=TMath::Mean(nselected,vecTheta.GetMatrixArray());
    Double_t meanTheor=TMath::Mean(nselected,vecdEdxTheor.GetMatrixArray());
    Double_t meanRob, rmsRob;   
    //
    for (Int_t idet=0; idet<6; idet++){
      vecMeanTrunc60[idet]=0;
      vecRMSTrunc60[idet]=0;
      vecMeanTrunc80[idet]=0;
      vecRMSTrunc80[idet]=0;
      vecMeanTruncAll[idet]=0;
      vecRMSTruncAll[idet]=0;
      if  (nselectedArray[idet]<20) continue;
      AliMathBase::EvaluateUni( nselectedArray[idet], vecArray[idet]->GetMatrixArray(), meanRob, rmsRob, nselectedArray[idet]*0.6);
      vecMeanTrunc60[idet]=meanRob;
      vecRMSTrunc60[idet]=rmsRob;
      AliMathBase::EvaluateUni( nselectedArray[idet], vecArray[idet]->GetMatrixArray(), meanRob, rmsRob, nselectedArray[idet]*0.8);
      vecMeanTrunc80[idet]=meanRob;
      vecRMSTrunc80[idet]=rmsRob;
      vecMeanTruncAll[idet]= TMath::Mean(nselectedArray[idet],vecArray[idet]->GetMatrixArray());
      vecRMSTruncAll[idet] = TMath::RMS(nselectedArray[idet],vecArray[idet]->GetMatrixArray());
    }
    //
    if (isel%10==0) printf("%d\n",isel);
    Double_t dEdxTheorMeanP = AliMathBase::BetheBlochAleph(meanP/mass);
    Double_t dEdxTheorCenterP = AliMathBase::BetheBlochAleph(mcp/mass);
    (*pcstream)<<"pid"<<
      "ipid="<<ipid<<                     // pid type -0 pion, 1- proton, 2 electron
      "nAll="<<nselected<<                // number of primary points
      "year="<<year<<                     // year
      "period="<<period<<                 // period
      //
      "dEdxTheorMean="<<meanTheor<<                 // initial dEdx hypothesis - mean over given bin
      "dEdxTheorMeanP="<<dEdxTheorMeanP<<           // initial dEdx hypothesis - for mean momenta
      "dEdxTheorCenterP="<<dEdxTheorCenterP<<       // initial dEdx hypothesis - for central momenta
      // bin
      "p="<<mcp<<                                   // center of bin for particle momentum
      "theta="<<mctheta<<                           // particle theta
      //
      "mass="<<mass<<                               // mass of particle
      "meanPt="<<meanPt<<                           // mean pt in bin
      "meanP="<<meanP<<                              
      "meanTheta="<<meanTheta<<
      //
      "dedxMean60.="<<&vecMeanTrunc60<<
      "dedxRMS60.="<<&vecRMSTrunc60<<
      "dedxMean80.="<<&vecMeanTrunc80<<
      "dedxRMS80.="<<&vecRMSTrunc80<<
      "dedxMeanAll.="<<&vecMeanTruncAll<<
      "dedxRMSAll.="<<&vecRMSTruncAll<<
      "\n";
    //
    // Test 
    //
  }
  delete pcstream;

}




/*



void  AnalyzedEdx(){
  //
  pcstream = new TTreeSRedirector("trend.root","recreate");
  TString flist = gSystem->GetFromPipe("cat highpt.list");
  TObjArray * array = flist.Tokenize("\n");
  array->Print();
  Int_t entries = array->GetEntries();
  for (Int_t ifile=0; ifile<entries; ifile++){
    //TString name = array->At(ifile)->GetName();
    printf("\n\n\n");
    printf("Analyzing:\t%d\t%s\n", ifile,array->At(ifile)->GetName());
    AnalyzedEdxFile(Form("%s#FilterEvents_Trees.root",array->At(ifile)->GetName())); 
    printf("\n\n\n");
  }
  delete pcstream;
}

void AnalyzedEdxFile(const char * fname){
  //
  //scomment
    const char * fname="/hera/alice/local/filtered/alice/data/2013/LHC13e/000196201/vpass1/root_archive.zip#FilterEvents_Trees.root";
    const char * fname="/hera/alice/local/filtered/alice/data/2013/LHC13c/000195531/ESDs/pass1/root_archive.zip#FilterEvents_Trees.root";
  //ecomment
  //
  TFile * fin=TFile::Open(fname);
  if (!fin){
    printf("File\t%s  not existing\n",fname);
    return ;
  }
  TTree * treeIn = (TTree*)fin->Get("highPt");  
  if (!fin){
    printf("Tree \t%s  not existing\n",fname);
    return ;
  }
  treeIn->SetCacheSize(1000000000.0);  
  // TFile *fout = TFile::Open("aaa.root","recreate");
  
  tree =treeIn->CopyTree("esdTrack.fITSncls>0||(esdTrack.fTOFsignal<100000-1)","",50000);
  if (!tree) return;
  if (tree->GetEntries()<2000) return;
  SetPionAliases();
  FitSlopeROC();
}

void SetPionAliases(){
  //
  // Select the tracks with higher dEdx - at high momenta it is used in order to get clean sample of Pions
  // IROC selection used to select sample for OROC
  // OROC selection used to select sample for IROC
  //
  Int_t entries = 0;
  tree->SetAlias("pionOROC","(1+0)");
  tree->SetAlias("pionIROC","(1+0)");
  
  Double_t mean=0, sigma=0;
  entries = tree->Draw("runNumber","","goff",1000);
  run = tree->GetV1()[0];
  {for (Int_t iter=0; iter<2; iter++){
    entries = tree->Draw("esdTrack.fTPCdEdxInfo.fTPCsignalRegion[3]/AliMathBase::BetheBlochAleph(esdTrack.P()/0.1)>>his","esdTrack.fTPCdEdxInfo.fTPCsignalRegion[3]>0&&pionIROC", "goff");
    AliMathBase::EvaluateUni(entries, tree->GetV1(), mean,sigma, 0.6*entries);
    tree->SetAlias("pionOROC",Form("esdTrack.fTPCdEdxInfo.fTPCsignalRegion[3]/AliMathBase::BetheBlochAleph(esdTrack.P()/0.1)-%f>%f",mean,0.0*sigma));
    entries = tree->Draw("esdTrack.fTPCdEdxInfo.fTPCsignalRegion[0]/AliMathBase::BetheBlochAleph(esdTrack.P()/0.1)>>his","esdTrack.fTPCdEdxInfo.fTPCsignalRegion[0]>0&&pionOROC", "goff");
    AliMathBase::EvaluateUni(entries, tree->GetV1(), mean,sigma, 0.6*entries);
    tree->SetAlias("pionIROC",Form("esdTrack.fTPCdEdxInfo.fTPCsignalRegion[0]/AliMathBase::BetheBlochAleph(esdTrack.P()/0.1)-%f>%f",mean,0.0*sigma));
    printf("%d\n",iter);
    printf("%s\n",tree->GetAlias("pionIROC"));
    printf("%s\n",tree->GetAlias("pionOROC"));
    }}    
  
}


void FitSlopeROC(){
  //
  // fit dEdx/eta for "identified high pt pions"
  //
  TF1 f1IROC("f1IROC","pol1");
  TF1 f1OROC("f1IROC","pol1");
  TVectorD vectorIROC(2);
  TVectorD vectorOROC(2);
  TH2 * hisOROC=0;
  TH2 * hisIROC=0;
  TObjArray fitIROC(3);
  TObjArray fitOROC(3);
  tree->Draw("esdTrack.fTPCdEdxInfo.fTPCsignalRegion[3]/AliMathBase::BetheBlochAleph(esdTrack.P()/0.1):abs(esdTrack.fP[3])>>hisOROC(5,0,1,100,20,80)","esdTrack.fTPCncls>70&&pionIROC","colz");
  hisOROC = (TH2*)tree->GetHistogram()->Clone();
  tree->Draw("esdTrack.fTPCdEdxInfo.fTPCsignalRegion[0]/AliMathBase::BetheBlochAleph(esdTrack.P()/0.1):abs(esdTrack.fP[3])>>hisIROC(5,0,1,100,20,80)","esdTrack.fTPCncls>70&&pionOROC","colz");
  hisIROC = (TH2*)tree->GetHistogram()->Clone();
  hisOROC->FitSlicesY(0,0,-1,0,"QNR",&fitOROC);
  hisIROC->FitSlicesY(0,0,-1,0,"QNR",&fitIROC);
  TH1 * phisOROC[3];
  TH1 * phisIROC[3];
  for (Int_t ihis=0; ihis<3; ihis++){
    phisOROC[ihis]=(TH1*)fitOROC.At(ihis);
    phisIROC[ihis]=(TH1*)fitIROC.At(ihis);
  }
  //
  phisOROC[1]->SetMarkerColor(2);
  phisOROC[1]->SetMarkerStyle(25); 
  phisIROC[1]->SetMarkerColor(2);
  phisIROC[1]->SetMarkerStyle(25);   
  phisOROC[1]->Fit(&f1OROC);
  f1OROC.GetParameters(vectorOROC.GetMatrixArray());
  phisOROC[1]->Draw();
  phisIROC[1]->Fit(&f1IROC);
  f1IROC.GetParameters(vectorIROC.GetMatrixArray());
  phisIROC[1]->Draw();
  //
  (*pcstream)<<"dedxFit"<<
    "run="<<run<<
    "hisIROC1.="<<phisIROC[1]<<  // IROC - dEdx/dEdxExp - mean as function of theta
    "hisIROC2.="<<phisIROC[2]<<  // IROC - dEdx/dEdxExp -        RMS as function of the theta
    "hisOROC1.="<<phisOROC[1]<<  // OROC - dEdx/dEdxExp - mean as function of theta
    "hisOROC2.="<<phisOROC[2]<<  // OROC - dEdx/dEdxExp -        RMS as function of theta 
    "fIROC.="<<&vectorIROC<<     // IROC vector with fit parameters mean (at 0 theta) and splope IROC
    "fOROC.="<<&vectorOROC<<     // OROC vector 
    "f1IROC.="<<&f1IROC<<          
    "f1OROC.="<<&f1OROC<<
    "\n";

  //scomment
  //resolution at given bin divided 
  // 
  Standard plots:
  TStatToolkit::MakeGraphSparse(chain,"hisOROC2.fArray[1]/hisOROC1.fArray[1]:run","",25,1,0)->Draw("alp");
  TStatToolkit::MakeGraphSparse(chain,"hisIROC2.fArray[1]/hisIROC1.fArray[1]:run","",25,1,0)->Draw("alp")
  //
  TStatToolkit::MakeGraphSparse(chain,"fIROC.fElements[0]+fIROC.fElements[1]*0.5:run","",25,1,0)->Draw("alp");
  
  //ecomment
  
}


void DrawdEdx(){
  //
  //pdg.GetParticle("Lambda0")->Mass
  //
  TDatabasePDG pdg;
  Double_t massLambda = pdg.GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg.GetParticle("K0")->Mass();
  Double_t massPion = pdg.GetParticle("pi+")->Mass();
  Double_t massProton = pdg.GetParticle("proton")->Mass();
 
  chain->SetAlias("massPion",Form("(%f+0)",massPion));
  chain->SetAlias("massProton",Form("(%f+0)",massProton));
  chain->SetAlias("massK0",Form("(%f+0)",massK0));
  chain->SetAlias("massLambda",Form("(%f+0)",massLambda));
  chain->SetAlias("K0Pull","(v0.GetEffMass(2,2)-massK0)/v0.GetKFInfo(2,2,1)");
  chain->SetAlias("LPull","(v0.GetEffMass(4,2)-massLambda)/v0.GetKFInfo(4,2,1)");
  chain->SetAlias("K0sel", "abs(v0.GetEffMass(2,2)-massK0)<min(abs(v0.GetEffMass(4,2)-massLambda), abs(v0.GetEffMass(2,4)-massLambda))&&abs(K0Pull)<3");
  //
  //
  //
  chain->Draw("track0.fTPCsignal/AliMathBase::BetheBlochAleph(track0.fIp.P()/massPion)>>hisP(50,20,80)","K0sel","");
  chain->Draw("track1.fTPCsignal/AliMathBase::BetheBlochAleph(track1.fIp.P()/massPion)>>+hisN(50,20,80)","K0sel","");

  chain->Draw("track0.fTPCsignal/AliMathBase::BetheBlochAleph(track0.fIp.P()/massPion):track0.fIp.P()>>his(20,0,10,100,30,70)","K0sel","colz");


  TH3F * hisBGPion3D = new TH3F("hisBGPion3D","hisBGPion3D",40,0.5,200,4,0,1,100,30,70);
  TH3F * hisBGProton3D = new TH3F("hisBGProton3D","hisBGProton3D",40,0.5,200,4,0,1,100,30,70);
  AliTPCcalibBase::BinLogX(hisBGPion3D->GetXaxis());
  AliTPCcalibBase::BinLogX(hisBGProton3D->GetXaxis());
  
  chain->Draw("track0.fTPCsignal/AliMathBase::BetheBlochAleph(track0.fIp.P()/massPion):abs(track0.fIp.fP[3]):track0.fIp.P()/massPion>>hisBGPion3D","K0sel||(track0.fTOFr[2]+track0.fTOFr[1])>0.6","");
  chain->Draw("track0.fTPCsignal/AliMathBase::BetheBlochAleph(track0.fIp.P()/massProton):abs(track0.fIp.fP[3]):track0.fIp.P()/massProton>>hisBGProton3D","(!K0sel&&abs(LPull)<1.5)||track0.fTOFr[4]>0.6","");

  
  TTreeSRedirector *pcstream = new TTreeSRedirector("bb.root","recreate");
  TVectorD vecPionMean(4),vecPionRMS(4), vecPionEntries(4);
  TVectorD vecProtonMean(4),vecProtonRMS(4), vecProtonEntries(4);
  TVectorD vecTheta(4);
  {for (Int_t ibg=2; ibg<38; ibg++){
    for (Int_t itheta=0; itheta<4; itheta++){
      TH1 * hisPion = hisBGPion3D->ProjectionZ("hisPion",ibg-1, ibg+1, itheta+1, itheta+1);
      TH1 * hisProton = hisBGProton3D->ProjectionZ("hisProton",ibg-1, ibg+1, itheta+1, itheta+1);
      hisPion->Fit("gaus");      
      vecPionMean[itheta]=gaus->GetParameter(1);
      vecPionRMS[itheta]=gaus->GetParameter(2);
      vecPionEntries[itheta]=hisPion->GetEntries();
      hisProton->Fit("gaus");      
      vecProtonMean[itheta]=gaus->GetParameter(1);
      vecProtonRMS[itheta]=gaus->GetParameter(2);
      vecProtonEntries[itheta]=hisProton->GetEntries();
      vecTheta[itheta]=hisBGPion3D->GetXaxis()->GetBinCenter(itheta);
      delete hisPion;
      delete hisProton;
    }
    //
    Double_t bg = hisBGPion3D->GetXaxis()->GetBinCenter(ibg);
    Double_t dEdxRef = AliMathBase::BetheBlochAleph(bg);
    (*pcstream)<<"bb"<<
      "bg="<<bg<<
      "dedxRef="<<dEdxRef<<
      "ibg="<<ibg<<
      "pionMean.="<<&vecPionMean<<
      "pionRMS.="<<&vecPionRMS<<
      "pionEntries.="<<&vecPionEntries<<
      "protonMean.="<<&vecProtonMean<<
      "protonRMS.="<<&vecProtonRMS<<
      "protonEntries.="<<&vecProtonEntries<<
      "\n";
    }}

  delete pcstream;
  
  

  
  hisBGPion->FitSlicesY();
  hisBGProton->FitSlicesY();

  hisBGPion_1->SetMarkerStyle(25);
  hisBGProton_1->SetMarkerStyle(25);
  hisBGPion_1->SetMarkerColor(2);
  hisBGProton_1->SetMarkerColor(4);

  hisBGPion_1->Draw();
  hisBGProton_1->Draw("same");



}



void aaa(){
  
  TChain * chain  = AliXRDPROOFtoolkit::MakeChainRandom("/hera/alice/miranov/highptAll.list","V0s",0,1000);
  chain->SetCacheSize(1000000000);

  TChain * chainHPT  = AliXRDPROOFtoolkit::MakeChainRandom("/hera/alice/miranov/highptAll.list","highPt",0,1000);
  chainHPT->SetCacheSize(1000000000);
  // for Kaons not good referenece data
  TFile ftofSelected("tofSelected.root","recreate");
  TTree * treeTOF = chainHPT->CopyTree("abs(esdTrack.fTOFsignalDz)<5&&abs(esdTrack.fTOFsignalDx)<5&&esdTrack.fITSncls>0&&abs(1/esdTrack.fP[4])<2");
  treeTOF->Write("treeTOF");
  //ftofSelected.Write();
  
}






void DrawdEdxRatio(){
  //
  //pdg.GetParticle("Lambda0")->Mass
  //
  TDatabasePDG pdg;
  Double_t massLambda = pdg.GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg.GetParticle("K0")->Mass();
  Double_t massPion = pdg.GetParticle("pi+")->Mass();
  Double_t massProton = pdg.GetParticle("proton")->Mass();

  chain->SetAlias("massPion",Form("(%f+0)",massPion));
  chain->SetAlias("massProton",Form("(%f+0)",massProton));
  chain->SetAlias("massK0",Form("(%f+0)",massK0));
  chain->SetAlias("massLambda",Form("(%f+0)",massLambda));
  chain->SetAlias("K0Pull","(v0.GetEffMass(2,2)-massK0)/v0.GetKFInfo(2,2,1)");
  chain->SetAlias("LPull","(v0.GetEffMass(4,2)-massLambda)/v0.GetKFInfo(4,2,1)");
  chain->SetAlias("K0sel", "abs(v0.GetEffMass(2,2)-massK0)<min(abs(v0.GetEffMass(4,2)-massLambda), abs(v0.GetEffMass(2,4)-massLambda))&&abs(K0Pull)<3");
 

  //
  //
  //

  TH3F * hisBGPion3DR1 = new TH3F("hisBGPion3DR1","hisBGPion3DR1",40,0.5,200,4,0,1,100,0.5,1.5);
  TH3F * hisBGProton3DR1 = new TH3F("hisBGProton3DR1","hisBGProton3DR1",40,0.5,200,4,0,1,100,0.5,1.5);
  TH3F * hisBGPion3DR2 = new TH3F("hisBGPion3DR2","hisBGPion3DR2",40,0.5,200,4,0,1,100,0.5,1.5);
  TH3F * hisBGProton3DR2 = new TH3F("hisBGProton3DR2","hisBGProton3DR2",40,0.5,200,4,0,1,100,0.5,1.5);
  AliTPCcalibBase::BinLogX(hisBGPion3DR1->GetXaxis());
  AliTPCcalibBase::BinLogX(hisBGProton3DR1->GetXaxis());
  AliTPCcalibBase::BinLogX(hisBGPion3DR2->GetXaxis());
  AliTPCcalibBase::BinLogX(hisBGProton3DR2->GetXaxis());

  chain->Draw("track0.fTPCdEdxInfo.fTPCsignalRegion[1]/track0.fTPCdEdxInfo.fTPCsignalRegion[0]:abs(track0.fIp.fP[3]):track0.fIp.P()/massPion>>hisBGPion3DR1","K0sel||(track0.fTOFr[2]+track0.fTOFr[1])>0.6","");
  chain->Draw("track0.fTPCdEdxInfo.fTPCsignalRegion[2]/track0.fTPCdEdxInfo.fTPCsignalRegion[0]:abs(track0.fIp.fP[3]):track0.fIp.P()/massPion>>hisBGPion3DR2","K0sel||(track0.fTOFr[2]+track0.fTOFr[1])>0.6","");
  //
  //
  //
  chain->Draw("track0.fTPCdEdxInfo.fTPCsignalRegion[1]/track0.fTPCdEdxInfo.fTPCsignalRegion[0]:abs(track0.fIp.fP[3]):track0.fIp.P()/massProton>>hisBGProton3DR1","(!K0sel&&abs(LPull)<1.5)||track0.fTOFr[4]>0.6","");
  chain->Draw("track0.fTPCdEdxInfo.fTPCsignalRegion[2]/track0.fTPCdEdxInfo.fTPCsignalRegion[0]:abs(track0.fIp.fP[3]):track0.fIp.P()/massProton>>hisBGProton3DR2","(!K0sel&&abs(LPull)<1.5)||track0.fTOFr[4]>0.6","");

  TFile fratio("bbRatio.root","recreate");
  hisBGPion3DR1->Write();
  hisBGPion3DR2->Write();
  hisBGProton3DR1->Write();
  hisBGProton3DR2->Write();
  fratio.Close();

  Double_t masses[5]={massPion, massProton,0,0,0};
  const Double_t kB2C=-0.299792458e-3;

  AliTPCParamSR par;
  par.Update();

  TTreeSRedirector *pcstream = new TTreeSRedirector("bbRatio.root","update");
  TVectorD vecMean(4),vecRMS(4), vecEntries(4);
  TVectorD vecTheta(4);
  TVectorD vecP(4);
  TVectorD vecPt(4);
  TVectorD vecPhi0(4);  // tracklet angle at 0 
  TVectorD vecPhi1(4);  // tracklet angle at 1 
  TVectorD vecPhi2(4);  // tracklet angle at 2 
  TVectorD vecTheta0(4);  // tracklet angle at 0 
  TVectorD vecTheta1(4);  // tracklet angle at 1 
  TVectorD vecTheta2(4);  // tracklet angle at 2 
  TVectorD vecL0(4);  // tracklet length at 0 
  TVectorD vecL1(4);  // tracklet letgth at 1 
  TVectorD vecL2(4);  // tracklet length at 2 
  {
    for (Int_t ptype=0; ptype<2; ptype++){
      for (Int_t itype=0; itype<2; itype++){
	for (Int_t ibg=2; ibg<38; ibg++){
	  Double_t bg = hisBGPion3DR1->GetXaxis()->GetBinCenter(ibg);
	  Double_t dEdxRef = AliMathBase::BetheBlochAleph(bg);
	  Double_t length0=0;
	  for (Int_t itheta=0; itheta<4; itheta++){
	    TH1 * hisProjection = 0;
	    if (itype==0){
	      if (ptype==0) hisProjection   = hisBGPion3DR1->ProjectionZ("hisPion",ibg-1, ibg+1, itheta+1, itheta+1);
	      if (ptype==1) hisProjection   = hisBGProton3DR1->ProjectionZ("hisProton",ibg-1, ibg+1, itheta+1, itheta+1);
	    }
	    if (itype==1){
	      if (ptype==0) hisProjection   = hisBGPion3DR2->ProjectionZ("hisPion",ibg-1, ibg+1, itheta+1, itheta+1);
	      if (ptype==1) hisProjection   = hisBGProton3DR2->ProjectionZ("hisProton",ibg-1, ibg+1, itheta+1, itheta+1);
	    }
	    hisProjection->Fit("gaus");      
	    vecMean[itheta]=gaus->GetParameter(1);
	    vecRMS[itheta]=gaus->GetParameter(2);
	    vecEntries[itheta]=hisProjection->GetEntries();
	    vecTheta[itheta]=hisBGPion3DR1->GetYaxis()->GetBinCenter(itheta+1);
	    vecP[itheta]= bg*masses[ptype];
	    vecPt[itheta]=bg*masses[ptype]/TMath::Sqrt(1+vecTheta[itheta]*vecTheta[itheta]);  //Pt
	    //
	    // Get angle and length
	    //
	    Double_t crv= 5*kB2C/vecPt[itheta];   //GetC(b); // bz*kB2C/pt;
	    Double_t angleIROC= TMath::ASin(TMath::Min(TMath::Abs(par.GetPadRowRadii(0, par.GetNRow(0)/2.)*crv)*0.5,1.));
	    Double_t angleOROC0= TMath::ASin(TMath::Min(TMath::Abs(par.GetPadRowRadii(36, par.GetNRowUp1()/2.)*crv)*0.5,1.));
	    Double_t angleOROC1= TMath::ASin(TMath::Min(TMath::Abs(par.GetPadRowRadii(36, par.GetNRowUp1()+par.GetNRowUp2()/2.)*crv)*0.5,1.));
	    vecPhi0[itheta]=angleIROC;
	    vecPhi1[itheta]=angleOROC0;
	    vecPhi2[itheta]=angleOROC1;
	    //
	    vecTheta0[itheta]=TMath::Sqrt(1+TMath::Tan(vecPhi0[itheta])*TMath::Tan(vecPhi0[itheta]))*vecTheta[itheta];
	    vecTheta1[itheta]=TMath::Sqrt(1+TMath::Tan(vecPhi1[itheta])*TMath::Tan(vecPhi1[itheta]))*vecTheta[itheta];
	    vecTheta2[itheta]=TMath::Sqrt(1+TMath::Tan(vecPhi2[itheta])*TMath::Tan(vecPhi2[itheta]))*vecTheta[itheta];

	    //
	    vecL0[itheta]  = 0.75*TMath::Sqrt(1+TMath::Tan(vecPhi0[itheta])*TMath::Tan(vecPhi0[itheta]) + vecTheta0[itheta]*vecTheta0[itheta]);
	    vecL1[itheta]  = 1.0*TMath::Sqrt(1+TMath::Tan(vecPhi1[itheta])*TMath::Tan(vecPhi1[itheta])  +  vecTheta1[itheta]*vecTheta1[itheta]);
	    vecL2[itheta]  = 1.5*TMath::Sqrt(1+TMath::Tan(vecPhi2[itheta])*TMath::Tan(vecPhi2[itheta])  +  vecTheta2[itheta]*vecTheta2[itheta]);
	    delete hisProjection;
	  }
	  //
	  (*pcstream)<<"bbRatio"<<
	    "ptype="<<ptype<<            // particle type
	    "rtype="<<itype<<            // ratio type
	    "mass="<<masses[ptype]<<     // mass of the 
	    "bg="<<bg<<
	    "dedxRef="<<dEdxRef<<
	    "ibg="<<ibg<<
	    "vecTheta.="<<&vecTheta<<
	    "vecTheta0.="<<&vecTheta0<<
	    "vecTheta1.="<<&vecTheta1<<
	    "vecTheta2.="<<&vecTheta2<<
	    "ratioMean.="<<&vecMean<<
	    "ratioRMS.="<<&vecRMS<<
	    "ratioEntries.="<<&vecEntries<<
	    "vecP.="<<&vecP<<             // momentum
	    "vecPt.="<<&vecPt<<             // pt
	    "vecPhi0.="<<&vecPhi0<<       // angle in the IROC
	    "vecPhi1.="<<&vecPhi1<<       // angle in the OROC1
	    "vecPhi2.="<<&vecPhi2<<       // angle in the OROC2
	    "vecL0.="<<&vecL0<<       // angle in the IROC
	    "vecL1.="<<&vecL1<<       // angle in the OROC1
	    "vecL2.="<<&vecL2<<       // angle in the OROC2
	    "\n";
	}}
    }
  }

  delete pcstream;
  
}
 
void Fit(){
  
  TFile ff("bbRatio.root");  
  TTree * tree = ff.Get("bbRatio");
  tree->SetMarkerStyle(25);
  //
  Double_t chi2=0;
  Int_t    npoints=0;
  TVectorD param;
  TMatrixD covar;
  Int_t npointsMax=10000000;

  TCut cutFit="ratioRMS.fElements<0.3&&ratioEntries.fElements>100&&vecL1.fElements<10&&vecL2.fElements<10";
  TString fstringFast="";  
  fstringFast+="(1-(rtype==0)*2)++";
  fstringFast+="(rtype==0)*vecTheta.fElements++";
  fstringFast+="(rtype==1)*vecTheta.fElements++";
  //  fstringFast+="vecTheta.fElements++";
  //fstringFast+="vecPhi0.fElements++";
  fstringFast+="1/(vecL0.fElements*dedxRef)++";
  //fstringFast+="(rtype==0)*vecPhi1.fElements++";
  fstringFast+="(rtype==0)/(vecL1.fElements*dedxRef)++";
  //fstringFast+="(rtype==1)*vecPhi2.fElements++";
  fstringFast+="(rtype==1)/(vecL2.fElements*dedxRef)++";
  //
  TString *strDelta = TStatToolkit::FitPlaneConstrain(tree,"ratioMean.fElements:ratioRMS.fElements", fstringFast.Data(),cutFit, chi2,npoints,param,covar,0.9,0, npointsMax, 1);
  TObjArray* tokArr = strDelta->Tokenize("++");
  tokArr->Print();
  tree->SetAlias("corr",strDelta->Data());
}


void SelectTracks(){
  //
  //
  //
  TFile f("tofSelected.root");
  TChain * chain = f.Get("highPt");
  chain->SetAlias("ITSTOFrPion","sqrt((esdTrack.fTOFr[1]+esdTrack.fTOFr[2]+0.05)*(esdTrack.fITSr[1]+esdTrack.fITSr[2]+0.05))");
  chain->SetAlias("ITSTOFrProton","sqrt((esdTrack.fTOFr[4]+0.05)*(esdTrack.fITSr[4]+0.05))");
  chain->SetAlias("ITSTOFrKaon","sqrt((esdTrack.fTOFr[3]+0.05)*(esdTrack.fITSr[3]+0.05))");
  TDatabasePDG pdg;
  Double_t massLambda = pdg.GetParticle("Lambda0")->Mass();
  Double_t massK0 = pdg.GetParticle("K0")->Mass();
  Double_t massPion = pdg.GetParticle("pi+")->Mass();
  Double_t massProton = pdg.GetParticle("proton")->Mass();

  chain->SetAlias("massPion",Form("(%f+0)",massPion));
  chain->SetAlias("massProton",Form("(%f+0)",massProton));
  chain->SetAlias("massK0",Form("(%f+0)",massK0));

  // 
  chain->SetAlias("isPion", "((ITSTOFrPion>1.5*ITSTOFrProton&&ITSTOFrPion>1.5*ITSTOFrKaon)&&esdTrack.fTRDsignal/AliMathBase::BetheBlochAleph(esdTrack.P()/0.13-1)<0.2)");
  //
  chain->SetAlias("isKaon", "((ITSTOFrKaon>1.5*ITSTOFrProton&&ITSTOFrKaon>1.5*ITSTOFrPion))&&abs(esdTrack.fTRDsignal/AliMathBase::BetheBlochAleph(esdTrack.P()/massK0)-1)<0.3");

  chain->SetAlias("isProton", "((ITSTOFrProton>1.5*ITSTOFrKaon&&ITSTOFrProton>1.5*ITSTOFrPion))");
    
}

*/
