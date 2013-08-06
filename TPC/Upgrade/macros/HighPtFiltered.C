void DrawPTResol(){
  //
  // Example macro to create the 1/pt resolution plot for Marek
  //
  TFile * f = TFile::Open("Filtered.root");
  TTree * treePt= (TTree*)f->Get("highPt");
  //
  //
  // 
  TCut cutNcl = "esdTrack.GetTPCClusterInfo(3,1)>120&&esdTrack.fITSncls>4";
  TH2 *phis1PtPt[3]={0};
  TH1 *his1PtPtRes[3]={0};
  TObjArray  * fitArray = new TObjArray(3);
  //
  //
  treePt->Draw("abs(esdTrack.fP[4])-1/particle.Pt():1/particle.Pt()>>his1Pt1Pt(20,0,0.5,100,-0.01,0.01)",cutNcl,"colz");
  phis1PtPt[0] = (TH2*)treePt->GetHistogram()->Clone();
  phis1PtPt[0]->FitSlicesY(0,0,-1,0,"QNR",fitArray);
  his1PtPtRes[0]=(TH1*)fitArray->At(2)->Clone();
  //
  his1PtPtRes[0]->Draw();
}


void DrawMatchingEffiency(){
  //
  // problem with ITS simulation looks like
  //
  TFile * f = TFile::Open("Filtered.root");  
  TTree * treePt= (TTree*)f->Get("highPt");
  treePt->SetAlias("ITSrefit","(esdTrack.fFlags&0x4)!=0");
  //
  //
  TCut cutNcl = "esdTrack.GetTPCClusterInfo(3,1)>120&&abs(esdTrack.fP[3])<0.9"; 
  TCut cutPileUp = "abs(particle.fVt)<0.000000001&&abs(esdTrack.fP[1])<15"; 
  TCut cutFindable = "particle.R()<0.2&&nrefITS>5"; 
  treePt->Draw("ITSrefit:1/particle.Pt()>>hisMatching(20,0,2)",cutNcl+cutPileUp+cutFindable,"prof");
}




void DrawMatchingEffiency(){
  //
  //
  //  
  TCut cutNcl = "esdTrack.GetTPCClusterInfo(3,1)>120&&abs(esdTrack.fP[3])<0.9"; 
  TCut cutPileUp = "abs(particle.fVt)<0.000000001"; 
  TCut cutFindable = "particle.R()<0.2&&tpcTrackLength"; 
  TCut cutBug="abs(vtxESD.fPosition[2]-particle.fVz)<0.01";  // fix the bug after  
  //
  //
  TChain* chains[20]={0};
  TProfile * hefFindable[20]={0};
  chains[0]=AliXRDPROOFtoolkit::MakeChain("pileup_4.list","MCEffTree",0,1000,0);
  chains[1]=0;
  //
  chains[2]=AliXRDPROOFtoolkit::MakeChain("pileup_6.list","MCEffTree",0,1000,0);
  chains[3]=AliXRDPROOFtoolkit::MakeChain("pileup_6_gem.list","MCEffTree",0,1000,0);
  //
  chains[4]=AliXRDPROOFtoolkit::MakeChain("pileup_8.list","MCEffTree",0,1000,0);
  //
  chains[6]=AliXRDPROOFtoolkit::MakeChain("pileup_10.list","MCEffTree",0,1000,0);
  chains[7]=AliXRDPROOFtoolkit::MakeChain("pileup_10_gem.list","MCEffTree",0,1000,0);
  //
  //  
  TFile *fhisto = TFile::Open("histoEff.root","update");
  for (Int_t ihis=0;ihis<20; ihis++){
    if (!chains[ihis]) continue;
    char hname[1000];
    snprintf(hname,100,"EffFindable_%d_%d",ihis%2,((ihis/2)+2)*2);
    printf("%d\t%s\n",ihis,hname);
    //
    hefFindable[ihis] = (TProfile*)fhisto->Get(hname);
    if (!hefFindable[ihis]){
      chains[ihis]->SetMarkerStyle(25);    
      chains[ihis]->SetCacheSize(1000000000);
      chains[ihis]->Draw("isRec:1/particle.Pt()>>his(10,0,4)",cutPileUp+cutFindable+cutBug,"prof");
      hefFindable[ihis]=(TProfile*)(chains[ihis]->GetHistogram()->Clone());
      hefFindable[ihis]->SetName(hname);
      fhisto->cd();
      hefFindable[ihis]->Write(hname);
    }
  }

  TLegend * legend = new TLegend(0.11,0.11,0.5,0.4,"TPC Efficiency for findable tracks");
  for (Int_t ihis=0; ihis<20; ihis++){
    if (hefFindable[ihis]==0) continue;
    hefFindable[ihis]->GetXaxis()->SetTitle("1/p_{T}");
    hefFindable[ihis]->GetYaxis()->SetTitle("#epsilon");
    hefFindable[ihis]->SetMinimum(0.85);
    hefFindable[ihis]->SetMaximum(1.01);
    if (ihis%2==0) hefFindable[ihis]->SetMarkerStyle(21);
    if (ihis%2==1) hefFindable[ihis]->SetMarkerStyle(25);
    hefFindable[ihis]->SetMarkerColor(1+ihis/2);    
    if (ihis==0) hefFindable[ihis]->Draw();
    hefFindable[ihis]->Draw("same");
    legend->AddEntry(hefFindable[ihis],Form("GEM%d Pileup%d",ihis%2,((ihis/2)+2)*2));
  }
  legend->Draw();

  
}


/*
 for a in `ls -d /hera/alice/mkowalsk/alice/simulations/gas/pileup*/`; do 
   echo $a; 
   dname=`basename $a`;
   ls $a/Filtered/*/F*.root > $dname.list
 done;


  

*/
