void Hshuttle(Int_t runTime=1500)
{
// this macro is to simulate the functionality of SHUTTLE.
// Here the list of DCS aliases is created and packed in TMap of structure "alias name" - TObjArray of AliDCSValue; AliDCSValue is a pair value-time stamp     
// currently simulated: freon temperature 2 per radiator (inlet,outlet)
//                      methane pressure 1 per chamber   
  gSystem->Load("libTestShuttle.so");
  
  AliTestShuttle::SetMainCDB(TString("local://$HOME"));
  
  TMap        *pDcsMap = new TMap;       pDcsMap->SetOwner(1);          //DCS archive map

  SimMap(pDcsMap,runTime);
  TestShuttle(pDcsMap);  
  TCanvas *c=new TCanvas("cc","ff",600,600);  
  DrawInput(c,pDcsMap);
  DrawOutput();
}//Hshuttle()
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
      
    pDcsMap->Add(new TObjString(Form(AliHMPIDPreprocessor::fP,iCh,iCh,iCh)),pP); 
    pDcsMap->Add(new TObjString(Form(AliHMPIDPreprocessor::fHV,iCh,iCh,iCh)),pHV); 
        
    for(Int_t iRad=0;iRad<3;iRad++){//radiators loop
      TObjArray *pT1=new TObjArray; pT1->SetOwner(1); 
      TObjArray *pT2=new TObjArray; pT2->SetOwner(1); 
      for (Int_t time=0;time<runTime;time+=stepTime)  pT1->Add(new AliDCSValue(13,time));  //sample inlet temperature    Nmean=1.292 @ 13 degrees
      for (Int_t time=0;time<runTime;time+=stepTime)  pT2->Add(new AliDCSValue(13,time));  //sample outlet temperature
      pDcsMap->Add(new TObjString(Form(AliHMPIDPreprocessor::fT1,iCh,iCh,iRad)) ,pT1); 
      pDcsMap->Add(new TObjString(Form(AliHMPIDPreprocessor::fT2,iCh,iCh,iRad)),pT2);
    }//radiators loop    
  }//chambers loop
}//SimMap()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void TestShuttle(TMap *pDcsMap)
{
  AliTestShuttle* pShuttle = new AliTestShuttle(0,0,1000000);   
  pShuttle->SetDCSInput(pDcsMap);                                                    //DCS map
  AliPreprocessor* pp = new AliHMPIDPreprocessor(pShuttle);                           //actual ipreprocessor is created here
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
                          
    TObjArray *pHV=(TObjArray*)pDcsMap->GetValue(Form(AliHMPIDPreprocessor::fHV,iCh,iCh,iCh,iCh)); //HV
    TObjArray *pP =(TObjArray*)pDcsMap->GetValue(Form(AliHMPIDPreprocessor::fP,iCh,iCh,iCh)); //P
    TGraph *pGr=new TGraph; pGr->SetMarkerStyle(5);
    
    TIter nextp(pP); cnt=0; while((pVal=(AliDCSValue*)nextp())){ pGr->SetPoint(cnt++,pVal->GetTimeStamp(),pVal->GetFloat());}//P

    
    pGr->Draw("AP");
  }  
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DrawOutput()
{
  AliCDBManager::Instance()->SetDefaultStorage("local://$HOME");
  AliCDBEntry *pQthreEnt=AliCDBManager::Instance()->Get("HMPID/Calib/Qthre",0);
  AliCDBEntry *pNmeanEnt=AliCDBManager::Instance()->Get("HMPID/Calib/Nmean",0);
  
  if(!pQthreEnt || ! pNmeanEnt) return;
  
  TObjArray *pNmean=(TObjArray*)pNmeanEnt->GetObject(); 
  TObjArray *pQthre=(TObjArray*)pQthreEnt->GetObject(); 
  
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
