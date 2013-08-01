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
  //
  //
  TFile * f = TFile::Open("Filtered.root"); 
  TTree * treeMC = f->Get("MCEffTree");
  reeMC->Draw("isRec:1/particle.Pt()","mult>2500","prof");

  
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
