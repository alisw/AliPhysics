void Hshuttle(Int_t runTime=1500)
{// this macro is to simulate the functionality of SHUTTLE.
  gSystem->Load("libTestShuttle.so");
  AliTestShuttle::SetMainCDB(TString("local://$HOME/CDB"));
  
  TMap *pDcsMap = new TMap;       pDcsMap->SetOwner(1);          //DCS archive map
  
  AliTestShuttle* pShuttle = new AliTestShuttle(0,0,1000000);   
  SimPed();   for(Int_t ldc=1;ldc<=4;ldc++) pShuttle->AddInputFile(AliTestShuttle::kDAQ,"HMP","pedestals",Form("LDC%i",ldc),Form("HmpidPeds%i.tgz",ldc));
  SimMap(pDcsMap,runTime); pShuttle->SetDCSInput(pDcsMap);                                    //DCS map
  
  AliPreprocessor* pp = new AliHMPIDPreprocessor(pShuttle); pShuttle->Process();  delete pp;  //here goes preprocessor 

  DrawInput(pDcsMap); DrawOutput();
}//Hshuttle()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimPed()
{
  ofstream out;
  for(Int_t ddl=0;ddl<=13;ddl++){
    out.open(Form("HmpidPedDdl%02i.txt",ddl));
    out << 3 <<endl;
    for(Int_t row=1;row<=24;row++)
      for(Int_t dil=1;dil<=10;dil++)
        for(Int_t adr=0;adr<=47;adr++){
          Float_t mean  = 150+200*gRandom->Rndm();
          Float_t sigma = ddl+gRandom->Gaus();
          Int_t inhard=((Int_t(mean))<<9)+Int_t(mean+3*sigma);
          out << Form("%2i %2i %2i %5.2f %5.2f %x\n",row,dil,adr,mean,sigma,inhard);
        }

    Printf("file ped %02i created",ddl);
    out.close();
  }
  gSystem->Exec("tar cf HmpidPeds1.tgz HmpidPedDdl00.txt HmpidPedDdl01.txt HmpidPedDdl02.txt HmpidPedDdl03.txt");
  gSystem->Exec("tar cf HmpidPeds2.tgz HmpidPedDdl04.txt HmpidPedDdl05.txt HmpidPedDdl06.txt HmpidPedDdl07.txt");
  gSystem->Exec("tar cf HmpidPeds3.tgz HmpidPedDdl08.txt HmpidPedDdl09.txt HmpidPedDdl10.txt HmpidPedDdl11.txt");
  gSystem->Exec("tar cf HmpidPeds4.tgz HmpidPedDdl12.txt HmpidPedDdl13.txt");
  gSystem->Exec("rm -rf ped*.txt");
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimMap(TMap *pDcsMap,Int_t runTime=1500)
{
  Int_t stepTime=100; //time interval between mesuraments
  Int_t startTime=0;
  
  
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
    TObjArray *pP=new TObjArray;  pP->SetOwner(1);
    TObjArray *pHV=new TObjArray; pHV->SetOwner(1); 
    for(Int_t time=0;time<runTime;time+=stepTime)  pP->Add(new AliDCSValue((Float_t)1005.0 ,time));   //sample CH4 pressure [mBar]
                                                   pHV->Add(new AliDCSValue((Float_t)2010.0,time));   //sample chamber HV [V]
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_GAS/HMP_MP%i_GAS_PMWC.actual.value"           ,iCh,iCh,iCh)),pP); 
    pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_PW/HMP_MP%i_SEC0/HMP_MP%i_SEC0_HV.actual.vMon",iCh,iCh,iCh)),pHV); 
        
    for(Int_t iRad=0;iRad<3;iRad++){//radiators loop
      TObjArray *pT1=new TObjArray; pT1->SetOwner(1); 
      TObjArray *pT2=new TObjArray; pT2->SetOwner(1); 
      for (Int_t time=0;time<runTime;time+=stepTime)  pT1->Add(new AliDCSValue(13,time));  //sample inlet temperature    Nmean=1.292 @ 13 degrees
      for (Int_t time=0;time<runTime;time+=stepTime)  pT2->Add(new AliDCSValue(14,time));  //sample outlet temperature
      pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iIn_Temp",iCh,iCh,iRad)) ,pT1); 
      pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iOut_Temp",iCh,iCh,iRad)),pT2);
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
      cnt=0; while((pVal=(AliDCSValue*)nextT1())) pGrT1->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());
      cnt=0; while((pVal=(AliDCSValue*)nextT2())) pGrT2->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());
      pGrT1->Draw("AP");  pGrT2->Draw("same");
    }//radiators loop
  }//chambers loop  
}//DrawInput()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawOutput()
{
  AliCDBManager::Instance()->SetDefaultStorage("local://$HOME/CDB"); AliCDBManager::Instance()->SetRun(0);
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
