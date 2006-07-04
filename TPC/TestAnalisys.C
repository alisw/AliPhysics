/*
 .L AliGenInfo.C+
 .L TestAnalisys.C
  AddChains(872);    // AddChains(runNumber);
  AddChains(873);    
  Select();          // make default selection of data
*/


void Select();
void AddChains(Int_t run);
void PRFYZ(TCut cut0, TCut cut1,  char * description);
void ResYZ(TCut cut0, TCut cut1,  char * description);
TProfile * ProfileARow(TCut cut0, Int_t max);




//
// global variables
//
TChain chaincl("Tracks","Tracks");   // tpc tracks and clusters
TChain chainFit("Fit","Fit");        // fitted signals with fit parameters
TChain chainPed("Fit","Fit");        // fitted pedestal with noise
//
AliComparisonDraw comp;
comp.fTree = &chaincl;
AliComparisonDraw compF;
compF.fTree = &chainFit;
AliComparisonDraw compP;
compP.fTree = &chainPed;
//
//
// selection of data for analysis
//
TEventList * listTracks    = new TEventList("listTracks","listTracks");
TEventList * listFitS      = new TEventList("listFitS","listFitS");
TEventList * listFitPed    = new TEventList("listFitPed","listFitPed");


void AddChains(Int_t run){
  //
  // add files to the chains
  //
  ifstream in0;
  ifstream in1;
  ifstream in2;
  TString sfile;
  char strcl[100];
  // TPC tracks
  //
  sprintf(strcl,"ls  run%d*/TPCtracks.root > /tmp/files.txt", run);
  gSystem->Exec(strcl);
  in0.open("/tmp/files.txt");
  for (;in0>>sfile;){
    if (sfile.Length()==0) break;
    printf("%s\n",sfile.Data());
    chaincl.Add(sfile.Data());
  }
  //
  // Fitted signals
  sprintf(strcl,"ls  run%d*/FitSignal.root > /tmp/files.txt", run);
  gSystem->Exec(strcl);
  in1.open("/tmp/files.txt");
  for (;in1>>sfile;){
    if (sfile.Length()==0) break;
    printf("%s\n",sfile.Data());
    chainFit.Add(sfile.Data());
  }
  //
  // Fitted pedestal
  sprintf(strcl,"ls  run%d*/TPCsignal.root > /tmp/files.txt", run);
  gSystem->Exec(strcl);
  in2.open("/tmp/files.txt");
  for (;in2>>sfile;){
    if (sfile.Length()==0) break;
    printf("%s\n",sfile.Data());
    chainPed.Add(sfile.Data());
  }
}

void Select(){
  //
  // base cut on the tracks
  //
  comp.fTree->Draw(">>listTracks","Track.fN>50&&abs(Track.fP4)<0.001");
  comp.fTree->SetEventList(listTracks);
  //
}



void PRFYZ(TCut cut0, TCut cut1,  char * description){
  //
  // plot Pad response function as funtion of drift z
  //
  //
  TF1 * f1 = new TF1("fdiff","sqrt([0]*[0]+(250-x)*[1]*[1])");
  f1->SetParameter(1,0.2);
  f1->SetParameter(0,0.2);
  comp.DrawXY("Cl.fZ","sqrt(Cl.fSigmaY2)","Track.fTrackPoints.GetAngleY()<0.05","Track.fTrackPoints.fTX>0"+cut0,5,10,240,-0,1);
  TH1F * prfInnerY = (TH1F*)comp.fMean->Clone();

  comp.DrawXY("Cl.fZ","sqrt(Cl.fSigmaY2)","Track.fTrackPoints.GetAngleY()<0.05","Track.fTrackPoints.fTX>0"+cut1,5,10,240,-0,1);
  TH1F * prfOuterY = (TH1F*)comp.fMean->Clone();
  //
  //
  prfOuterY->SetMinimum(0);
  prfOuterY->SetMarkerStyle(23);
  prfInnerY->SetMarkerStyle(24);
  prfOuterY->SetXTitle("Z position (cm)");
  prfOuterY->SetYTitle("PRF width (cm)");
  char chouter[100];
  char chinner[100];
  prfOuterY->Fit(f1);
  sprintf(chouter,"Outer sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfInnerY->Fit(f1);
  sprintf(chinner,"Inner sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfOuterY->Draw();
  prfInnerY->Draw("same");
  TString desc = description;
  TLegend *legend = new TLegend(0.25,0.12,0.85,0.35, desc+"\nTPC cluster shape Fit: #sigma = #sqrt{p_{0}^{2}+(z_{d}-z)p_{1}^{2}}");
  legend->SetBorderSize(1);
  legend->AddEntry(prfOuterY,chouter);
  legend->AddEntry(prfInnerY,chinner);
  legend->Draw();
}


void ResYZ(TCut cut0, TCut cut1,  char * description){
  //
  // resolution in y coordinate as function of z
  //
  TF1 * f1 = new TF1("fdiff","sqrt([0]*[0]+(250-x)*[1]*[1])");
  f1->SetParameter(1,0.2);
  f1->SetParameter(0,0.2);
  comp.DrawXY("Cl.fZ","Track.fTrackPoints.GetY()-Cl.GetY()","Track.fTrackPoints.GetAngleY()<0.05","Track.fTrackPoints.fTX>0"+cut0,5,10,240,-0.5,0.5);
  TH1F * prfInnerY = (TH1F*)comp.fRes->Clone();

  comp.DrawXY("Cl.fZ","Track.fTrackPoints.GetY()-Cl.GetY()","Track.fTrackPoints.GetAngleY()<0.05","Track.fTrackPoints.fTX>0"+cut1,5,10,240,-0.5,0.5);
  TH1F * prfOuterY = (TH1F*)comp.fRes->Clone();
  //
  //
  prfOuterY->SetMinimum(0);
  prfOuterY->SetMaximum(0.15);  
  prfOuterY->SetMarkerStyle(23);
  prfInnerY->SetMarkerStyle(24);
  prfOuterY->SetXTitle("Z position (cm)");
  prfOuterY->SetYTitle("Y resolution (cm)");
  char chouter[100];
  char chinner[100];
  prfOuterY->Fit(f1);
  sprintf(chouter,"Outer sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfInnerY->Fit(f1);
  sprintf(chinner,"Inner sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfOuterY->Draw();
  prfInnerY->Draw("same");
  TString desc = description;
  TLegend *legend = new TLegend(0.25,0.12,0.85,0.35, desc+"TPC cluster resolution: #sigma = #sqrt{p_{0}^{2}+(z_{d}-z)p_{1}^{2}}");
  legend->SetBorderSize(1);
  legend->AddEntry(prfOuterY,chouter);
  legend->AddEntry(prfInnerY,chinner);
  legend->Draw();
}


TProfile * ProfileARow(TCut cut0, char *name, Int_t max){ 
  //
  // make profile histrogram of amplitudes
  //
  TProfile *profA = new TProfile(name,name,max,0,max-1);
  char expr[100];
  sprintf(expr,"Cl.fMax:Cl.fRow>>%s",name);
  comp.fTree->Draw(expr,"Cl.fZ>0&&Cl.fMax<500"+cut0,"prof");
  profA->SetXTitle("Pad Row");
  profA->SetYTitle("Amplitude (ADC)");
  return profA;
}
