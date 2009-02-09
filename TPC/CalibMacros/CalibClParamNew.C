/*
// Pseudo code
// 1. Load Calibration from DB

   .x $ALICE_ROOT/TPC/macros/ConfigOCDB.C

//
//2.  Load Calibration components
//
.x ~/UliStyle.C
gSystem->Load("libANALYSIS");
gSystem->Load("libSTAT");
gSystem->Load("libTPCcalib");
TFile fcalib("CalibObjects.root");
TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
AliTPCcalibTracksGain * gain = ( AliTPCcalibTracksGain *)array->FindObject("calibTracksGain");

//
//3.
// 
AliTPCClusterParam * paramCl = AliTPCcalibDB::Instance()->GetClusterParam(); 
gain->UpdateClusterParam(paramCl);

//
//4. Test - visulaize results
//
TF1 fdr("fdr","AliTPCClusterParam::SQnorm(0,0,x,0,0)",0,1);

//
//5. Update the ClusterParam in your OCDB
//
//
Int_t runNumber = 0; //61725
AliCDBMetaData *metaData= new AliCDBMetaData();
metaData->SetObjectClassName("AliTPCClusterParam");
metaData->SetResponsible("Marian Ivanov");
metaData->SetBeamPeriod(1);
metaData->SetAliRootVersion("05-06-00"); //root version
metaData->SetComment("October runs calibration");
AliCDBId id1("TPC/Calib/ClusterParam", runNumber, AliCDBRunRange::Infinity());
gStorage = AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT/OCDB");
gStorage->Put(paramCl, id1, metaData);


//
//  6. dEdx matching - Currently only debug streamers
//  Load Proof with OCDB 
    See pseudo code TestChainCosmicdEdx
*/





void TestChainCosmicDedx(){
  //
  // pseudo cose
  //
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCcalib");
  
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
   
  //
  // Q normalization - short - medium - long
  //

  TVectorD paramT(3);
  TVectorD paramM(3);
  TH1 * hisRatio =0;

  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,0,0,64,0,kFALSE)/Seed0.CookdEdxNorm(0.01,0.65,0,0,159,0,kFALSE)>>hisRatio(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>130",""); 
  hisRatio = (TH1*)gROOT->FindObject("hisRatio");
  hisRatio->Fit("gaus");
  paramT[0]=hisRatio->GetFunction("gaus")->GetParameter(1);
  //
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,0,64,127,0,kFALSE)/Seed0.CookdEdxNorm(0.01,0.65,0,0,159,0,kFALSE)>>hisRatio(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>130",""); 
  hisRatio = (TH1*)gROOT->FindObject("hisRatio");
  hisRatio->Fit("gaus");
  paramT[1]=hisRatio->GetFunction("gaus")->GetParameter(1);
  //
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,0,127,159,0,kFALSE)/Seed0.CookdEdxNorm(0.01,0.65,0,0,159,0,kFALSE)>>hisRatio(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>130",""); 
  hisRatio = (TH1*)gROOT->FindObject("hisRatio");
  hisRatio->Fit("gaus");
  paramT[2]=hisRatio->GetFunction("gaus")->GetParameter(1);
  //
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,1,0,64,0,kFALSE)/Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0,kFALSE)>>hisRatio(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>130",""); 
  hisRatio = (TH1*)gROOT->FindObject("hisRatio");
  hisRatio->Fit("gaus");
  paramM[0]=hisRatio->GetFunction("gaus")->GetParameter(1);
  //
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,1,64,127,0,kFALSE)/Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0,kFALSE)>>hisRatio(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>130",""); 
  hisRatio = (TH1*)gROOT->FindObject("hisRatio");
  hisRatio->Fit("gaus");
  paramM[1]=hisRatio->GetFunction("gaus")->GetParameter(1);
  //
  chainCosmic->Draw("Seed0.CookdEdxNorm(0.01,0.65,1,127,159,0,kFALSE)/Seed0.CookdEdxNorm(0.01,0.65,1,0,159,0,kFALSE)>>hisRatio(100,0.5,1.5)","min(Orig0.fTPCncls,Orig1.fTPCncls)>130",""); 
  hisRatio = (TH1*)gROOT->FindObject("hisRatio");
  hisRatio->Fit("gaus");
  paramM[2]=hisRatio->GetFunction("gaus")->GetParameter(1);
  //




		
 paramCl->fQpadTnorm=(TVectorD*)paramT->Clone();
 paramCl->fQpadMnorm=(TVectorD*)paramM->Clone();
 

 //check
 

}






