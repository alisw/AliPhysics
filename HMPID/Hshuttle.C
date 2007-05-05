void Hshuttle(Int_t runTime=1500)
{
// this macro is to simulate the functionality of SHUTTLE.
// Here the list of DCS aliases is created and packed in TMap of structure "alias name" - TObjArray of AliDCSValue; AliDCSValue is a pair value-time stamp     
// currently simulated: freon temperature 2 per radiator (inlet,outlet)
//                      methane pressure 1 per chamber   
  gSystem->Load("libTestShuttle.so");
  
  AliTestShuttle::SetMainCDB(TString("local://$HOME"));
  
  TMap        *pDcsMap = new TMap;       pDcsMap->SetOwner(1);          //DCS archive map

  SimPed();
  SimMap(pDcsMap,runTime);
  TestShuttle(pDcsMap);  
  TCanvas *c=new TCanvas("cc","ff",600,600);  
  DrawInput(c,pDcsMap);
  DrawOutput();
}//Hshuttle()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimPed()
{
  ofstream out;
  for(Int_t ddl=0;ddl<=13;ddl++){
    out.open(Form("ped_%02i.txt",ddl));
    out << 3 <<endl;
    for(Int_t row=1;row<=24;row++)
      for(Int_t dil=1;dil<=10;dil++)
        for(Int_t adr=0;adr<=47;adr++){
          Float_t mean  = 150+200*gRandom->Rndm();
          Float_t sigma = 1+0.2*gRandom->Rndm();
          Int_t inhard=((Int_t(mean))<<9)+Int_t(mean+3*sigma);
          out << Form("%2i %2i %2i %2i %5.2f %5.2f %x\n",ddl,row,dil,adr,mean,sigma,inhard);
        }

    Printf("file ped %02i created",ddl);
    out.close();
  }
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
    pDcsMap->Add(new TObjString(Form(AliHMPIDPreprocessor::GetP(),iCh,iCh,iCh)),pP); 
    pDcsMap->Add(new TObjString(Form(AliHMPIDPreprocessor::GetHV(),iCh,iCh,iCh)),pHV); 
        
    for(Int_t iRad=0;iRad<3;iRad++){//radiators loop
      TObjArray *pT1=new TObjArray; pT1->SetOwner(1); 
      TObjArray *pT2=new TObjArray; pT2->SetOwner(1); 
      for (Int_t time=0;time<runTime;time+=stepTime)  pT1->Add(new AliDCSValue(13,time));  //sample inlet temperature    Nmean=1.292 @ 13 degrees
      for (Int_t time=0;time<runTime;time+=stepTime)  pT2->Add(new AliDCSValue(13,time));  //sample outlet temperature
      pDcsMap->Add(new TObjString(Form(AliHMPIDPreprocessor::GetT1(),iCh,iCh,iRad)) ,pT1); 
      pDcsMap->Add(new TObjString(Form(AliHMPIDPreprocessor::GetT2(),iCh,iCh,iRad)),pT2);
    }//radiators loop    
  }//chambers loop
}//SimMap()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void TestShuttle(TMap *pDcsMap)
{
  AliTestShuttle* pShuttle = new AliTestShuttle(0,0,1000000);   
  pShuttle->SetDCSInput(pDcsMap);                                                    //DCS map
  for(Int_t ddl=0;ddl<=13;ddl++) pShuttle->AddInputFile(AliTestShuttle::kDAQ,"HMP","pedestals",Form("DDL%i",ddl),Form("./ped_%02i.txt",ddl));
  AliPreprocessor* pp = new AliHMPIDPreprocessor(pShuttle);                           //actual preprocessor is created here
  pShuttle->Process();                                                               //run SHUTTLE simulator
  delete pp;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawInput(TCanvas *c,TMap *pDcsMap)
{
  c->Divide(3,3);
  
  AliDCSValue *pVal; Int_t cnt=0;
  
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
    if(iCh==6) c->cd(1);  if(iCh==5) c->cd(2);                          
    if(iCh==4) c->cd(4);  if(iCh==3) c->cd(5);  if(iCh==2) c->cd(6);
                          if(iCh==1) c->cd(8);  if(iCh==0) c->cd(9); 
                          
    TObjArray *pHV=(TObjArray*)pDcsMap->GetValue(Form(AliHMPIDPreprocessor::GetHV(),iCh,iCh,iCh,iCh)); //HV
    TObjArray *pP =(TObjArray*)pDcsMap->GetValue(Form(AliHMPIDPreprocessor::GetP(),iCh,iCh,iCh)); //P
    TGraph *pGr=new TGraph; pGr->SetMarkerStyle(5);
    
    TIter nextp(pP); cnt=0; while((pVal=(AliDCSValue*)nextp())){ pGr->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());}//P

    
    pGr->Draw("AP");
  }  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawOutput()
{
  AliCDBManager::Instance()->SetDefaultStorage("local://$HOME");
  AliCDBEntry *pQthreEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Qthre",0);
  AliCDBEntry *pNmeanEnt =AliCDBManager::Instance()->Get("HMPID/Calib/Nmean",0);
  AliCDBEntry *pSigCutEnt=AliCDBManager::Instance()->Get("HMPID/Calib/SigCut",0);
  AliCDBEntry *pDaqSigEnt=AliCDBManager::Instance()->Get("HMPID/Calib/DaqSig",0);
  
  if(!pQthreEnt || ! pNmeanEnt || !pSigCutEnt || !pDaqSigEnt) return;
  
  TObjArray *pNmean =(TObjArray*)pNmeanEnt ->GetObject(); 
  TObjArray *pQthre =(TObjArray*)pQthreEnt ->GetObject(); 
  TObjArray *pSigCut=(TObjArray*)pSigCutEnt->GetObject(); 
  TObjArray *pDaqSig=(TObjArray*)pDaqSigEnt->GetObject();
   
  TF1 *pRad0,*pRad1,*pRad2;  
  TCanvas *c2=new TCanvas("c2","Nmean"); c2->Divide(3,3);
  
  
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
    if(iCh==6) c2->cd(1);  if(iCh==5) c2->cd(2);                          
    if(iCh==4) c2->cd(4);  if(iCh==3) c2->cd(5);  if(iCh==2) c2->cd(6);
                           if(iCh==1) c2->cd(8);  if(iCh==0) c2->cd(9); 
                          
    TF1 *pRad0=(TF1*)pNmean->At(iCh*3+0); pRad0->Draw();  pRad0->GetXaxis()->SetTimeDisplay(kTRUE); pRad0->GetYaxis()->SetRangeUser(1.28,1.3);
    TF1 *pRad1=(TF1*)pNmean->At(iCh*3+1); pRad1->Draw("same");
    TF1 *pRad2=(TF1*)pNmean->At(iCh*3+2); pRad2->Draw("same");
  }//chambers loop  
}
