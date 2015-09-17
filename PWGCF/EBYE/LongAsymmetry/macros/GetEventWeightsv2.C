#include "/Users/rashmiraniwala/ToALICE/ANALYSIS/SetStyle.C" 
TH2F * hVz_Cent[3], *hRatio1, *hRatio2;
Int_t centmin,centmax;
Int_t runmin, runmax;
Int_t nfile;
//Int_t nrun = 2;
Int_t nrun = 60;
Int_t RunNumber[]={
  
  137231, 137232, 137236, 137243, 137432,
  137434, 137440, 137441, 137443, 137541,
  137544, 137686, 137691, 137692, 137704,
  137718, 137722, 137724, 137751, 137752,
  137844, 137848,
  
  138190, 138192, 138197, 138201, 138225,
  138275, 138364, 138396, 138438, 138439,
  138442, 138469, 138534, 138578, 138579,
  138582, 138583,
  
  139028,  139029,139036, 139037, 139038,
  139105, 139107, 139173, 139309, 139310,
  139314, 139328, 139329, 139360, 139437,
  139438, 139465, 139503, 139505, 139507,
  139510
};
Int_t Debug=0;
Float_t zdcasym, fCentPercentile[4],fVertexZ;

void GetEventWeightsv2(Int_t incentmin=15, Int_t incentmax = 20, Int_t inRunMin=0, Int_t inRunMax=59){
  SetStyle();
  gStyle->SetOptTitle(1);
 
  centmin = incentmin;
  centmax = incentmax;
  runmin = inRunMin;
  runmax = inRunMax;
  
  for(Int_t ireg = 0;ireg<3;ireg++){
    Char_t hname[30];
    sprintf(hname,"hVz_Cent_Reg%i",ireg+1);
    hVz_Cent[ireg]=new TH2F(hname,hname,10,-5.0,5.0,centmax-centmin,centmin,centmax);
  }

  Int_t nevent1 = 0, nevent2 =0, nevent3=0;
  for(Int_t irun = runmin;irun<=runmax;irun++){
    nfile = 1;
    cout<<" Runnumber is "<<RunNumber[irun]<<endl;
    Char_t fname[120];
    sprintf(fname,"/Users/rashmiraniwala/ToALICE/ANALYSIS/V0MCentrality/Cent%3.3dto%3.3d/AsymRun%dCent%3.3dto%3.3d.root",centmin,centmax,RunNumber[irun],centmin,centmax);
    cout<<" Reading root file "<<fname<<" for asymTree "<<endl;
    
    TFile * f = new TFile(fname,"READ");
    TTree* asymTree = (TTree*)f->Get("asymTree");
    asymTree->SetBranchAddress("zdcasym",&zdcasym);
    asymTree->SetBranchAddress("fCentPercentile",fCentPercentile);
    asymTree->SetBranchAddress("fVertexZ",&fVertexZ);
    cout<<" Got asymTree with entries "<<asymTree->GetEntries()<<endl;
    for(Int_t iev = 0;iev<asymTree->GetEntries();iev++){
      asymTree->GetEntry(iev);
      //      cout<<" zdcasym, vertex, cent = "<<iev<<" "<<zdcasym<<" "<<fVertexZ<<" "<<fCentPercentile[0]<<endl;
      if(zdcasym>-1.0 && zdcasym<-0.1){
	nevent1++;
	hVz_Cent[0]->Fill(fVertexZ,fCentPercentile[0]);
      }
      if(zdcasym>0.1 && zdcasym<1.0){
	nevent2++;
	hVz_Cent[1]->Fill(fVertexZ,fCentPercentile[0]);
      }
      if(zdcasym>-0.1 && zdcasym<0.1){
	nevent3++;
	hVz_Cent[2]->Fill(fVertexZ,fCentPercentile[0]);
      }
    }
    f->Close();
    cout<<"nevent 1,2,3 = "<<nevent1<<" "<<nevent2<<" "<<nevent3<<endl;
  }
  for(Int_t ireg = 0;ireg<3;ireg++){
    cout<<" Number of events in region "<<ireg+1<<" = "<<hVz_Cent[ireg]->GetEntries()<<endl;
    hVz_Cent[ireg]->Scale(1.0/hVz_Cent[ireg]->GetEntries());
  }

  hRatio1=(TH2F*)hVz_Cent[0]->Clone("hRatio1");
  hRatio1->Divide(hVz_Cent[2]);

  hRatio2=(TH2F*)hVz_Cent[1]->Clone("hRatio2");
  hRatio2->Divide(hVz_Cent[2]);
  
  
  TCanvas * c = new TCanvas("c","c",1000,600);
  c->Divide(3,1);
  c->cd(1);
  hVz_Cent[0]->SetXTitle("Vz");
  hVz_Cent[0]->SetYTitle("Centrality %");
  hVz_Cent[0]->SetZTitle(" Number of Events");
  hVz_Cent[0]->Draw("lego");
  c->cd(2);
  hVz_Cent[1]->SetXTitle("Vz");
  hVz_Cent[1]->SetYTitle("Centrality %");
  hVz_Cent[1]->SetZTitle(" Number of Events");
  hVz_Cent[1]->Draw("lego");
  c->cd(3);
  hVz_Cent[2]->SetXTitle("Vz");
  hVz_Cent[2]->SetYTitle("Centrality %");
  hVz_Cent[2]->SetZTitle(" Number of Events");
  hVz_Cent[2]->Draw("lego");
  
  TCanvas * c1 = new TCanvas("c1","c1",900,500);
  c1->Divide(2,1);

  c1->cd(1);
  hRatio1->SetTitle("Event Weights for correction of Centrality and Vertex-Z distribution: Reg1");
  hRatio1->SetLineColor(4);
  hRatio1->SetXTitle("Vz");
  hRatio1->SetYTitle("Centrality %");
  hRatio1->SetZTitle("Ratio of Reg1/Reg3");
  hRatio1->Draw("lego");

  hRatio2->SetTitle("Event Weights for correction of Centrality and Vertex-Z distribution: Reg2");
  
  c1->cd(2);
  hRatio2->SetLineColor(2);
  hRatio2->SetXTitle("Vz");
  hRatio2->SetYTitle("Centrality %");
  hRatio2->SetZTitle("Ratio of Reg2/Reg3");
  hRatio2->Draw("lego");

  TLatex * l = new TLatex(-10,20,"Vertex-Z and Centrality Correction Weights: Region 2");
  l->SetTextSize(0.7);
  l->Draw();

  
  sprintf(fname,"/Users/rashmiraniwala/ToALICE/ANALYSIS/V0MCentrality/Cent%3.3dto%3.3d/EventWeightsCent%3.3dto%3.3d_Run%ito%i.root",centmin,centmax,centmin,centmax,runmin,runmax);
  TFile * fout = new TFile(fname,"RECREATE");
  fout->cd();
  for(Int_t ireg = 0;ireg<3;ireg++){
    hVz_Cent[ireg]->Write();
  }
  hRatio1->Write();
  hRatio2->Write();

}
