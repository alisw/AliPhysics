void ITSSDDQAChecker(Char_t *filename){

  //-------opens file & histo
  TFile *fileinput = new TFile(filename,"read");
  const Int_t nSDDmodules=260;
  TH2F *ChargeMap;  
  TH2F *CountsMap; 
  TH1F *hModulePattern;
  TH1F *hModuleSidePattern;
  TH1F *mapProjX[2*nSDDmodules];      //260 dx e 260 sx  with A, T, Q
  TH1F *mapProjY[2*nSDDmodules]; 
  TProfile *mapProfX[2*nSDDmodules];
  TProfile *mapProfY[2*nSDDmodules];


  //-------creates module & side histo
  hModulePattern = (TH1F*)fileinput->Get("hModulePattern");
  hModuleSidePattern = (TH1F*)fileinput->Get("hModuleSidePattern");  
  TCanvas *modCanvas = new TCanvas("modCanvas","modCanvas");
  TCanvas *sideCanvas = new TCanvas("sideCanvas","sideCanvas");
  modCanvas->cd();
  hModulePattern->Draw();    
  modCanvas->Update();
  sideCanvas->cd();
  hModuleSidePattern->Draw(); 
  sideCanvas->Update();

  //-------creates Projections & Profiles, one canvas each module
  Char_t *takeChargeMap = new Char_t[50];
  Char_t *takeCountsMap = new Char_t[50];  
  Char_t canvname[100];  
  Char_t canvtitle[100];
 
  TCanvas *moduleCanvas[nSDDmodules];
  gStyle->SetPalette(1);  
  Int_t nActiveModules = 0;
  
  for(Int_t imod=0; imod<nSDDmodules; imod++){   //nSDDmodules
    if(hModulePattern->GetBinContent(imod+1)!=0){
      nActiveModules++;	
      cout<<imod<<" imod" <<endl;
      sprintf(canvtitle,"canvas_module_%d",imod);
      sprintf(canvname,"moduleCanvas[%d]",imod); 
      moduleCanvas[imod]=new TCanvas(canvname,canvtitle);
      moduleCanvas[imod]->Divide(2,6);
      moduleCanvas[imod]->Update();
      
      for(Int_t isid=0;isid<2;isid++){
        Int_t index=2*imod+isid;
	if(hModuleSidePattern->GetBinContent(index+1)!=0){
	  //cout << "Module: " << imod << ", Side: " << isid << ", index: " << index << ", update canvases" << endl;
	  sprintf(takeChargeMap,"chargeMap%d",index);
	  ChargeMap= (TH2F*)fileinput->Get(takeChargeMap);
	  sprintf(takeCountsMap,"countsMap%d",index);
	  CountsMap= (TH2F*)fileinput->Get(takeCountsMap);

	  mapProjX[index] = ChargeMap->ProjectionX();
	  mapProjY[index] = ChargeMap->ProjectionY();
	  mapProjX[index]->GetXaxis()->SetTitle("Time Bin");
	  mapProjX[index]->GetYaxis()->SetTitle("Total Counts");
	  mapProjY[index]->GetXaxis()->SetTitle("Anode");
	  mapProjY[index]->GetYaxis()->SetTitle("Total Counts");
       	
	  mapProfX[index] = ChargeMap->ProfileX();
	  mapProfY[index] = ChargeMap->ProfileY();
	  mapProfX[index]->GetXaxis()->SetTitle("Time Bin");
	  mapProfX[index]->GetYaxis()->SetTitle("Average Counts");
	  mapProfY[index]->GetXaxis()->SetTitle("Anode");
	  mapProfY[index]->GetYaxis()->SetTitle("Average Counts");
	  
	  moduleCanvas[imod] ->cd(1+isid);
	  ChargeMap->Draw("colz"); 
	  moduleCanvas[imod] ->cd(3+isid);
	  CountsMap->Draw("colz");
	  moduleCanvas[imod] ->cd(5+isid);
	  mapProjX[index]->Draw(); 
	  moduleCanvas[imod] ->cd(7+isid);
	  mapProjY[index]->Draw(); 
	  moduleCanvas[imod] ->cd(9+isid);
	  mapProfX[index]->Draw(); 
	  moduleCanvas[imod] ->cd(11+isid);
	  mapProfY[index]->Draw(); 
	  moduleCanvas[imod]->Update();
	  }

      }	
      gSystem->Exec("sleep 0.5");
    }
  } 	

  cout << nActiveModules << " Modules containing Data" << endl;

}

