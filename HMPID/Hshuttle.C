void Hshuttle(Int_t runTime=1500)
{// this macro is to simulate the functionality of SHUTTLE.
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");
//  AliTestShuttle::SetMainCDB(TString("local://$HOME/CDB"));
  AliTestShuttle::SetMainCDB(TString("local://$HOME"));
  
  TMap *pDcsMap = new TMap;       pDcsMap->SetOwner(1);          //DCS archive map
  
  AliTestShuttle* pShuttle = new AliTestShuttle(0,0,1000000);   
//  pShuttle->SetInputRunType("PHYSICS");
//  pShuttle->SetInputRunType("CALIBRATION");
  pShuttle->SetInputRunType("PHYSICS");
  SimPed();   
  for(Int_t ldc=51;ldc<=52;ldc++) 
  {
   if(ldc==51) {for(Int_t iddl=0;iddl<=7;iddl++)  pShuttle->AddInputFile(AliTestShuttle::kDAQ,"HMP",Form("HmpidPedDdl%02i.txt",iddl),Form("LDC%i",ldc),Form("HmpidPedDdl%02i.txt",iddl));}
   if(ldc==52) {for(Int_t iddl=8;iddl<=13;iddl++) pShuttle->AddInputFile(AliTestShuttle::kDAQ,"HMP",Form("HmpidPedDdl%02i.txt",iddl),Form("LDC%i",ldc),Form("HmpidPedDdl%02i.txt",iddl));}
  }
  SimMap(pDcsMap,runTime); pShuttle->SetDCSInput(pDcsMap);                                    //DCS map
  
  AliPreprocessor* pp = new AliHMPIDPreprocessor(pShuttle); pShuttle->Process();  delete pp;  //here goes preprocessor 

  DrawInput(pDcsMap); DrawOutput();
  gSystem->Exec("rm -rf HmpidPedDdl*.txt");

}//Hshuttle()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimPed()
{
  Int_t iDDLmin=0,iDDLmax=13;
  Int_t nSigmas = 3;                                            // value stored in the ddl files of pedestals
  ofstream out;
  for(Int_t ddl=iDDLmin;ddl<=iDDLmax;ddl++){
    out.open(Form("HmpidPedDdl%02i.txt",ddl));
    out << Form("%8s %2d\n","RunNumber",      999);             //read run number
    out << Form("%8s %2d\n","LdcId" ,         999);             //read LDC Id
    out << Form("%8s %2d\n","TimeStamp",      999);             //read time stamp
    out << Form("%8s %2d\n","TotNumEvt",      999);             //read number of total events processed
    out << Form("%8s %2d\n","TotDDLEvt",      999);             //read number of bad events for DDL # nDDL processed
    out << Form("%8s %2d\n","NumBadEvt",      999);             //read number of bad events for DDL # nDDL processed
    out << Form("%8s %2f\n","NBadE(%)",       999.9);           //read number of bad events (in %) for DDL # nDDL processed
    out << Form("%8s %2.2d\n","SigCut",      nSigmas);          //# of sigma cuts
    for(Int_t row=1;row<=24;row++)
      for(Int_t dil=1;dil<=10;dil++)
        for(Int_t adr=0;adr<=47;adr++){
          Float_t mean  = 150+200*gRandom->Rndm();
          Float_t sigma = 1+0.3*gRandom->Gaus();
          Int_t inhard=((Int_t(mean+nSigmas*sigma))<<9)+Int_t(mean); //right calculation, xchecked with Paolo 8/4/2008
          out << Form("%2i %2i %2i %5.3f %5.3f %4.4x \n",row,dil,adr,mean,sigma,inhard);
        }

    out.close();
  }
  Printf("HMPID - All %i DDL pedestal files created successfully",iDDLmax-iDDLmin+1);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimMap(TMap *pDcsMap,Int_t runTime=1500)
{
  Int_t stepTime=100; //time interval between measurements
  Int_t startTime=0;
  
  TObjArray *pHV[7];
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
    TObjArray *pPCH4=new TObjArray;  pPCH4->SetOwner(1);
    TObjArray *pPenv=new TObjArray;  pPenv->SetOwner(1);
    pHV[iCh]=new TObjArray; pHV[iCh]->SetOwner(1);
    for(Int_t time=0;time<runTime;time+=stepTime) {
       pPCH4->Add(new AliDCSValue((Float_t)4.0,time));               //sample CH4 pressure [mBar] differential respect to atm
       pPenv->Add(new AliDCSValue((Float_t)1000.0 ,time));               //also atm. pressure set to the same value
//       pHV[iCh]->Add(new AliDCSValue((Float_t)(1930+iCh*20),time));   //sample chamber HV [V]
       pHV[iCh]->Add(new AliDCSValue((Float_t)(2050),time));   //sample chamber HV [V]
    }
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_GAS/HMP_MP%i_GAS_PMWPC.actual.value"           ,iCh,iCh,iCh    )),pPCH4);         //CH4 pressure wrt atm
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_ENV/HMP_ENV_PENV.actual.value"                                               )),pPenv);         //atm pressure
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC0/HMP_MP%i_SEC0_HV.actual.vMon",iCh,iCh,iCh,iCh)),pHV[iCh]);      //HV SEC0
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC1/HMP_MP%i_SEC1_HV.actual.vMon",iCh,iCh,iCh,iCh)),pHV[iCh]);      //HV SEC1
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC2/HMP_MP%i_SEC2_HV.actual.vMon",iCh,iCh,iCh,iCh)),pHV[iCh]);      //HV SEC2 
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC3/HMP_MP%i_SEC3_HV.actual.vMon",iCh,iCh,iCh,iCh)),pHV[iCh]);      //HV SEC3 
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC4/HMP_MP%i_SEC4_HV.actual.vMon",iCh,iCh,iCh,iCh)),pHV[iCh]);      //HV SEC4 
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC5/HMP_MP%i_SEC5_HV.actual.vMon",iCh,iCh,iCh,iCh)),pHV[iCh]);      //HV SEC5 

    for(Int_t iRad=0;iRad<3;iRad++){//radiators loop
      TObjArray *pT1=new TObjArray; pT1->SetOwner(1);
      TObjArray *pT2=new TObjArray; pT2->SetOwner(1);
      for (Int_t time=0;time<runTime;time+=stepTime)  pT1->Add(new AliDCSValue((Float_t)(13.0+gRandom->Rndm()),time));  //sample inlet temperature                    Nmean=1.292 @ 13 degrees
      for (Int_t time=0;time<runTime;time+=stepTime)  pT2->Add(new AliDCSValue((Float_t)(18.0+gRandom->Rndm()),time));  //sample outlet temperature
      pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iIn_Temp",iCh,iCh,iRad)) ,pT1);               //Temperature in  Rad 
      pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iOut_Temp",iCh,iCh,iRad)),pT2);               //Temperature out Rad
    }//radiators loop
  }//chambers loop
}//SimMap()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawInput(TMap *pDcsMap)
{
  
  TCanvas *c=new TCanvas("cc","Input data",600,600);    c->Divide(3,3);
  
  AliDCSValue *pVal; Int_t cnt;
  
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
    if(iCh==6) c->cd(1);  if(iCh==5) c->cd(2);                          
    if(iCh==4) c->cd(4);  if(iCh==3) c->cd(5);  if(iCh==2) c->cd(6);
                          if(iCh==1) c->cd(8);  if(iCh==0) c->cd(9); 
                          
    TObjArray *pHV=(TObjArray*)pDcsMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC0/HMP_MP%i_SEC0_HV.actual.vMon",iCh,iCh,iCh,iCh)); //HV
    TObjArray *pP =(TObjArray*)pDcsMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_GAS/HMP_MP%i_GAS_PMWC.actual.value",iCh,iCh,iCh)); //P
    TGraph *pGrHV=new TGraph; pGrHV->SetMarkerStyle(5); TIter nextHV(pHV);
    TGraph *pGrP =new TGraph; pGrP ->SetMarkerStyle(5); TIter nextP (pP );
    
    for(Int_t iRad=0;iRad<3;iRad++){
      TObjArray *pT1=(TObjArray*)pDcsMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iIn_Temp",iCh,iCh,iRad));  TIter nextT1(pT1);
      TObjArray *pT2=(TObjArray*)pDcsMap->GetValue(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iOut_Temp",iCh,iCh,iRad)); TIter nextT2(pT2);
      TGraph *pGrT1=new TGraph; pGrT1->SetMarkerStyle(5); 
      TGraph *pGrT2=new TGraph; pGrT2->SetMarkerStyle(5); 
      pGrT1->SetMaximum(40);
      pGrT2->SetMaximum(40);
      cnt=0; while((pVal=(AliDCSValue*)nextT1())) pGrT1->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());
      cnt=0; while((pVal=(AliDCSValue*)nextT2())) pGrT2->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());
      if(iRad==0) pGrT1->Draw("AP"); else pGrT1->Draw("same");
      pGrT2->Draw("same");
    }//radiators loop
  }//chambers loop  
}//DrawInput()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawOutput()
{
//  AliCDBManager::Instance()->SetDefaultStorage("local://$HOME/CDB"); AliCDBManager::Instance()->SetRun(0);
  AliCDBManager::Instance()->SetDefaultStorage("local://$HOME"); AliCDBManager::Instance()->SetRun(0);
  AliCDBEntry *pQthreEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Qthre");
  AliCDBEntry *pNmeanEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Nmean");
  AliCDBEntry *pDaqSigEnt=AliCDBManager::Instance()->Get("HMPID/Calib/DaqSig");
  
  if(!pQthreEnt || !pNmeanEnt || !pDaqSigEnt) return;
  
  TObjArray *pNmean =(TObjArray*)pNmeanEnt ->GetObject(); 
  TObjArray *pQthre =(TObjArray*)pQthreEnt ->GetObject(); 
  TObjArray *pDaqSig=(TObjArray*)pDaqSigEnt->GetObject();
   
  TF1 *pRad0,*pRad1,*pRad2;  
  TCanvas *c2=new TCanvas("c2","Output"); c2->Divide(3,3);
  
  TH1F *pSig[7];
  
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
    TMatrix *pM=(TMatrix*)pDaqSig->At(iCh);
    
    pSig[iCh]=new TH1F(Form("sig%i",iCh),"Sigma;[QDC]",100,-5,20); //pSig[iCh]->SetLineColor(iCh+kRed);
    for(Int_t padx=0;padx<160;padx++) for(Int_t pady=0;pady<144;pady++) pSig[iCh]->Fill((*pM)(padx,pady));
    
    c2->cd(7);    if(iCh==0) pSig[iCh]->Draw(); else pSig[iCh]->Draw("same");
    
    if(iCh==6) c2->cd(1);  if(iCh==5) c2->cd(2);                          
    if(iCh==4) c2->cd(4);  if(iCh==3) c2->cd(5);  if(iCh==2) c2->cd(6);
                           if(iCh==1) c2->cd(8);  if(iCh==0) c2->cd(9); 
                          
    TF1 *pRad0=(TF1*)pNmean->At(iCh*3+0); pRad0->Draw();  pRad0->GetXaxis()->SetTimeDisplay(kTRUE); pRad0->GetYaxis()->SetRangeUser(1.28,1.3);
    TF1 *pRad1=(TF1*)pNmean->At(iCh*3+1); pRad1->Draw("same");
    TF1 *pRad2=(TF1*)pNmean->At(iCh*3+2); pRad2->Draw("same");
  }//chambers loop    
}//DrawOutput()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
