void HMPIDTestShuttle()
{
// this macro is to simulate the functionality of SHUTTLE.
// Here the list of DCS aliases is created and packed in TMap of structure "alias name" - TObjArray of AliDCSValue    
  TMultiGraph *pMG[7]; for(Int_t i=0;i<7;i++) {pMG[i]=new TMultiGraph; pMG[i]->SetTitle("T,grad C;time");}
  TGraph      *pGr[21];for(Int_t i=0;i<21;i++){pGr[i]=new TGraph;      pGr[i]->SetMarkerStyle(i%3+24); pGr[i]->SetMarkerColor(i%3+2); pMG[i/3]->Add(pGr[i]);}    
  TMap        *pDcsMap = new TMap;       pDcsMap->SetOwner(1);          //DCS archive map
  
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
    for(Int_t iRad=0;iRad<3;iRad++){//radiators loop
      TObjArray* pValLst  = new TObjArray;  pValLst->SetOwner(1);                
      Int_t iPoint=0;
      for (Int_t time=0;time<1000;time+=50) {
        AliDCSValue*    pVal = new AliDCSValue(Float_t(iCh*3+iRad+0.1*gRandom->Gaus()), time);         //sample new data point 
        pValLst->Add(pVal);                                                                            //add it to the list
        pGr[3*iCh+iRad]->SetPoint(iPoint++,time,pVal->GetFloat());                                     //and also to the graph  
      }
      pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iIn_Temp",iCh,iCh,iRad)),pValLst); 
      pDcsMap->Add(new TObjString(Form("HMP_DET/HMP_MP%i/HMP_MP%i_LIQ_LOOP.actual.sensors.Rad%iOut_Temp",iCh,iCh,iRad)),pValLst);
    }//radiators loop
  }//chambers loop
  
  
  AliCDBManager::Instance()->SetDefaultStorage("local://$HOME/tstCDB"); // initialize location of CDB
      
  gSystem->Load("libTestShuttle.so"); 
  Int_t iRun=1;   
  AliTestShuttle* pShuttle = new AliTestShuttle(iRun,0,100000);   
  pShuttle->SetDCSInput(pDcsMap);                                                    //DCS map
  AliPreprocessor* pp = new AliHMPIDPreprocessor(pShuttle);                           //actual ipreprocessor is created here
  pShuttle->Process();                                                               //run SHUTTLE simulator
  delete pp;
  
    
  AliCDBEntry *pTempEn=AliCDBManager::Instance()->Get("HMPID/DCS/RadTemp",iRun);
  if(!pTempEn) {Printf("ERROR file is not retrieved!!!");return;}

  TObjArray *pTempLst=(TObjArray*)pTempEn->GetObject(); TF1 *pRad0,*pRad1,*pRad2;  
  TCanvas *pC=new TCanvas; pC->Divide(3,3);
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
    if(iCh==6) pC->cd(1);  if(iCh==5) pC->cd(2);                          //this is just to see the input
    if(iCh==4) pC->cd(4);  if(iCh==3) pC->cd(5);  if(iCh==2) pC->cd(6);
                           if(iCh==1) pC->cd(8);  if(iCh==0) pC->cd(9); 
    pMG[iCh]->Draw("ap"); pMG[iCh]->GetXaxis()->SetTimeDisplay(kTRUE);
    pRad0=(TF1*)pTempLst->At(iCh*3+0); pRad0->Draw("same");
    pRad1=(TF1*)pTempLst->At(iCh*3+1); pRad1->Draw("same");
    pRad2=(TF1*)pTempLst->At(iCh*3+2); pRad2->Draw("same");
  }  
  
  AliCDBEntry *pIdxEn=AliCDBManager::Instance()->Get("HMPID/DCS/MeanIdx",iRun);
  if(!pIdxEn) {Printf("ERROR file is not retrieved!!!");return;}

  TObjArray *pIdxLst=(TObjArray*)pIdxEn->GetObject(); TF1 *pRad0,*pRad1,*pRad2;  
  TCanvas *pC=new TCanvas("c2","Ref Idx"); pC->Divide(3,3);
  for(Int_t iCh=0;iCh<7;iCh++){//chambers loop
    if(iCh==6) pC->cd(1);  if(iCh==5) pC->cd(2);                          //this is just to see the input
    if(iCh==4) pC->cd(4);  if(iCh==3) pC->cd(5);  if(iCh==2) pC->cd(6);
                           if(iCh==1) pC->cd(8);  if(iCh==0) pC->cd(9); 
    pRad0=(TF1*)pIdxLst->At(iCh*3+0); pRad0->Draw();  pRad0->GetXaxis()->SetTimeDisplay(kTRUE); pRad0->GetYaxis()->SetRangeUser(1.28,1.3);
    pRad1=(TF1*)pIdxLst->At(iCh*3+1); pRad1->Draw("same");
    pRad2=(TF1*)pIdxLst->At(iCh*3+2); pRad2->Draw("same");
  }    
}//Test()
