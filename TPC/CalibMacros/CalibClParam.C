/// \file CalibClParam.C
///
/// Procedures to create a cluster parametererization
/// Currently only pseudo code - once will be stable - to be updated as a "normal" macro

void PseudoCode(){
  /// That rough sequence to update  a Cluster param calibration using debug streamers

  //
  //0. Load libraries
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCcalib");
  
  
  //1. Load Parameters to be modified:
  //e.g:
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0) ;
  AliTPCClusterParam * param = AliTPCcalibDB::Instance()->GetClusterParam();  
  //
  //2. Load chain from debug streamers
  //
  //e.g
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  AliXRDPROOFtoolkit tool;  
  TChain * chainres = tool.MakeChain("tracks.txt","ResolCl",0,10200);
  chainres->Lookup();  
  TChain * chainGain = tool.MakeChain("gain.txt","dEdx",0,1000000);
  chainGain->Lookup();
  TChain * chainCosmic =  tool.MakeChain("cosmic.txt","Track0",0,1000000);
  chainCosmic->Lookup();
  //
  // 3. Do fits - pos correction
  //
  AliTPCcalibTracks::MakeQPosNormAll(chainres,param,200000,0) ;

  //
  //4. Do fits gain
  // 
  param->fPosQTnorm[0] = AliTPCcalibTracksGain::MakeQPosNorm(chainGain,0,kFALSE,200000,kTRUE);
  param->fPosQTnorm[1] = AliTPCcalibTracksGain::MakeQPosNorm(chainGain,1,kFALSE,200000,kTRUE);
  param->fPosQTnorm[2] = AliTPCcalibTracksGain::MakeQPosNorm(chainGain,2,kFALSE,200000,kTRUE);
  //
  param->fPosQMnorm[0] = AliTPCcalibTracksGain::MakeQPosNorm(chainGain,0,kTRUE,200000,kTRUE);
  param->fPosQMnorm[1] = AliTPCcalibTracksGain::MakeQPosNorm(chainGain,1,kTRUE,200000,kTRUE);
  param->fPosQMnorm[2] = AliTPCcalibTracksGain::MakeQPosNorm(chainGain,2,kTRUE,200000,kTRUE);
  //
  //5. Store results
  //
  TFile f("paramout.root","recreate");
  param->Write("clusterParam");
  f.Close();


  //6. Verification
  TFile f2("paramout.root");
  AliTPCClusterParam *param2 = (AliTPCClusterParam*)f2.Get("clusterParam");
  param2->SetInstance(param2);
  //
  // position correction
  chainres->Draw("AliTPCClusterParam::SPosCorrection(1,0,Cl.fPad,Cl.fTimeBin,Cl.fZ,Cl.fSigmaY2,Cl.fSigmaZ2,Cl.fMax)","Cl.fDetector<36","",10000); // should be line 
  // gain correction
  TCut cutA("dedge>3&&fraction2<0.7");
  chainGain->Draw("(Cl.fMax/gain/dedxM.fElements[2]):AliTPCClusterParam::SQnormPos(0,1,Cl.fPad,Cl.fTimeBin,Cl.fZ,Cl.fSigmaY2,Cl.fSigmaZ2,Cl.fMax,Cl.fQ)","IPad==0"+cutA,"prof",100000);


}


void UpdateParam(){
  /// Pseudo code -to update cluster params

  .L $ALICE_ROOT/TPC/Cal/AliTPCCreateDummyCDB.C
  TFile f2("paramout.root");
  AliTPCClusterParam *param2 = (AliTPCClusterParam*)f2.Get("clusterParam");  
  AliCDBManager *man = AliCDBManager::Instance();
  gStorLoc = man->GetStorage(gCDBpath);
  metaData = CreateMetaObject("AliTPCClusterParam");
  StoreObject("TPC/Calib/ClusterParam/", param2, metaData);
  
}


void TestChainCosmicDedx(){
  /// pseudo cose

  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCcalib");
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(0) ;
  AliTPCClusterParam * param = AliTPCcalibDB::Instance()->GetClusterParam();  
  param->SetInstance(param);
  //e.g
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  AliXRDPROOFtoolkit tool;  
  TChain * chainCosmic =  tool.MakeChain("cosmic.txt","Track0",0,1000000);
  chainCosmic->Lookup();



  TCut cutT("cutT","abs(Tr1.fP[3]+Tr0.fP[3])<0.01");  // OK
  TCut cutD("cutD","abs(Tr0.fP[0]+Tr1.fP[0])<2");     // OK
  TCut cutP1("cutP1","abs(Tr0.fP[1]-Tr1.fP[1])<3");   // OK
  TCut cutPt("cutPt","abs(Tr1.fP[4]+Tr0.fP[4])<0.1&&abs(Tr0.fP[4])+abs(Tr1.fP[4])<10");
  TCut cutN("cutN","min(Orig0.fTPCncls,Orig1.fTPCncls)>100");
  TCut cutA=cutT+cutD+cutPt+cutN+cutP1;
  
  TCut cutS("cutS","Orig0.fIp.fP[1]*Orig1.fIp.fP[1]>0");

  if (gProof) chainCosmic->SetProof(kTRUE);
  //
  //
  //
  //
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0)/Seed1.CookdEdxNorm(0.01,0.65,1,0,159,0)","min(Orig0.fTPCncls,Orig1.fTPCncls)>120","",10000);
  
  //
  // Q normalization - short - medium - long
  //

  TVectorD paramT(3);
  TVectorD paramM(3);
  TH1F * his =0;
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,0,0,64,0,kTRUE)/Seed0.CookdEdxNorm(0.01,0.65,0,0,159,0,kTRUE)>>his(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>140","",100000); 
  TH1 * his = (TH1*)gROOT->FindObject("his");
  paramT[0]=his->GetMean();
  //
 chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,0,64,129,0,kTRUE)/Seed0.CookdEdxNorm(0.01,0.65,0,0,159,0,kTRUE)>>his(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>140","",100000);  
 TH1 * his = (TH1*)gROOT->FindObject("his");
  paramT[1]=his->GetMean();
  //
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,0,129,159,0,kTRUE)/Seed0.CookdEdxNorm(0.01,0.65,0,0,159,0,kTRUE)>>his(100,0.5,2.0)","min(Orig0.fTPCncls,Orig1.fTPCncls)>140","",100000); 
 TH1 * his = (TH1*)gROOT->FindObject("his");
 paramT[2]=his->GetMean();
  //
  //
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,1,0,64,0,kTRUE)/Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0,kTRUE)>>his(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>140","",100000); 
  TH1 * his = (TH1*)gROOT->FindObject("his");
  paramM[0]=his->GetMean();
  //
 chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,1,64,129,0,kTRUE)/Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0,kTRUE)>>his(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>140","",100000); 
 TH1 * his = (TH1*)gROOT->FindObject("his");
 paramM[1]=his->GetMean();
  //
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,1,129,159,0,kTRUE)/Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0,kTRUE)>>his(100,0.5,2.0)","min(Orig0.fTPCncls,Orig1.fTPCncls)>140","",100000); 
 TH1 * his = (TH1*)gROOT->FindObject("his");
 paramM[2]=his->GetMean();

 param->fQpadTnorm=(TVectorD*)paramT->Clone();
 param->fQpadMnorm=(TVectorD*)paramM->Clone();
 
 TFile f("paramout.root","recreate");
 param->Write("clusterParam");
 f.Close();

}


void dedxDemo(){
  chainCosmic->Draw("(Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0,kTRUE,kTRUE)+Seed1.CookdEdxNorm(0.01,0.65,1,0,159,0,kTRUE,kTRUE))*0.5:log(Tr0.P())>>his2(100,-2,5,500,0,200)","min(Orig0.fTPCncls,Orig1.fTPCncls)>110");   
}
