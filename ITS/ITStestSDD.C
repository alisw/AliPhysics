{

gROOT->Reset();

TFile f("SDD_histos_test.root");

Int_t nbins = 18;
Float_t dmax = 36.;

a = local->ProjectionX();
t = local->ProjectionY();

TF1 *faa = new TF1("faa","gaus",-100.,100.);
a->Fit("faa","R","Q",-100.,100.);
TF1 *fat = new TF1("fat","gaus",-50,50);
t->Fit("fat","R","Q",-50,50);

TH1F *anode_resolution = new TH1F("anoder","Anode resolution vs Drift Path",nbins,0.,dmax);
TH1D *anodes[nbins];
Float_t res_anodes[nbins];
Float_t errres_anodes[nbins];
Float_t dmin = 60.;

for(Int_t i=0;i<nbins;i++){
  if(i>6) dmin = 70.;
  if(i>11) dmin = 90.;
  TString *aa = new TString("aa_");
  Char_t ai[2];
  sprintf(ai,"%d",i+1);
  aa->Append(ai);
  anodes[i] = at->ProjectionY(aa->Data(),i,i+1);
  TF1 *fa = new TF1("fa","gaus",-1.*dmin,dmin);
  anodes[i]->Fit("fa","R","Q",-1.*dmin,dmin);
  res_anodes[i] = fa->GetParameter(2);
  Float_t RMS = anodes[i]->GetRMS();
  //if(res_anodes[i] > RMS) 
  res_anodes[i] = RMS;
  errres_anodes[i] = fa->GetParError(2);
  anode_resolution->Fill(i*dmax/nbins+dmax/(2*nbins),res_anodes[i]);
  anode_resolution->SetBinError(i+1,(Stat_t) errres_anodes[i]);
  anode_resolution->SetMarkerColor(2);
  anode_resolution->SetLineColor(2);
}

f->cd();

TH1F *time_resolution = new TH1F("timer","Time resolution vs Drift Path",nbins,0.,dmax);
TH1D *times[nbins];
Float_t res_times[nbins];
Float_t errres_times[nbins];
for(Int_t i=0;i<nbins;i++){
  TString *ta = new TString("tt_");
  Char_t ti[2];
  sprintf(ti,"%d",i+1);
  ta->Append(ti);
  times[i] = tt->ProjectionY(ta->Data(),i,i+1);
  TF1 *ft = new TF1("ft","gaus",-50,50);
  times[i]->Fit("ft","R","Q",-50,50);
  res_times[i] = ft->GetParameter(2);
  Float_t RMS = times[i]->GetRMS();
  //if(res_times[i] > RMS) 
  res_times[i] = RMS;
  errres_times[i] = ft->GetParError(2);
  time_resolution->Fill(i*dmax/nbins+dmax/(2*nbins),res_times[i]);
  time_resolution->SetBinError(i+1,(Stat_t) errres_times[i]);
  time_resolution->SetMarkerColor(6);
  time_resolution->SetLineColor(6);
}
Float_t x1 = 4000.;
Float_t x2 = x1;
Float_t y1 = 94;
Float_t y2 = 84;

TMarker *m1 = new TMarker(x1,y1,21);
TMarker *m2 = new TMarker(x2,y2,22);
Text_t *text1 = "Anode";
Text_t *text2 = "Time";
TText *t1 = new TText(x1 + 250,y1-2,text1);
TText *t2 = new TText(x2 + 250,y2-2,text2);

anode_resolution->SetMarkerStyle(21);
anode_resolution->SetMaximum(100.);
time_resolution->SetMarkerStyle(23);

c1->Clear();
at->Draw();
c1->SaveAs("ITS_at.ps");

c1->Clear();
tt->Draw();
c1->SaveAs("ITS_tt.ps");

nanodes->ProfileX();
nanodes_pfx->Draw();
nsampls->ProfileX();
nsampls_pfx->Draw();
nanodes_pfx->SetLineColor(2);
nanodes_pfx->Draw("SAME");
c1->SaveAs("ITS_clsize.ps");

f.cd();
anoder->SetMinimum(0.);
anoder->Draw();
anoder->SetXTitle("Drift Path (mm)");
anoder->SetYTitle("Resolution (um)");
TMarker *mk0 = new TMarker(16.,54.,23);
mk0->SetMarkerColor(6);
mk0->Draw();
TMarker *mk3 = new TMarker(16.,46.,21);
mk3->SetMarkerColor(2);
mk3->Draw();
Text_t *text = "294 um pitch detector";
TText *t3 = new TText(2,90,text);
Text_t *text0 = "Spatial Resolution";
TText *t0 = new TText(2,82,text0);
Text_t *text2 = "Simulation (Anode)";
TText *t2 = new TText(18,44,text2);
Text_t *text4 = "Simulation (Time)";
TText *t4 = new TText(18,52,text4);
t4->Draw();
t3->Draw();
t2->Draw();
t0->Draw();
timer->Draw("SAME");
c1->SaveAs("ITS_res_check_294.ps");

f.cd();
nanodes_pfx->SetXTitle("Drift Time (ns)");
nanodes_pfx->SetYTitle("Anodes/Cluster");
nanodes_pfx->SetMaximum(3.);
nanodes_pfx->SetMarkerStyle(21);
nanodes_pfx->SetLineColor(2);
nanodes_pfx->SetMarkerColor(2);
nanodes_pfx->Draw();
TMarker *mk2 = new TMarker(4000.,0.3,21);
mk2->SetMarkerColor(2);
mk2->Draw();
Text_t *text = "294 um pitch detector";
TText *t3 = new TText(200,2.7,text);
Text_t *text0 = "Average Number of Anodes/Cluster";
TText *t0 = new TText(200,2.4,text0);
Text_t *text2 = "Simulation";
TText *t2 = new TText(4200,0.2,text2);
t3->Draw();
t0->Draw();
t2->Draw();
c1->SaveAs("ITS_and_check_294.ps");

f.cd();
nsampls_pfx->SetXTitle("Drift Time (ns)");
nsampls_pfx->SetYTitle("Time bins/Anode/Cluster");
nsampls_pfx->SetMaximum(10.);
nsampls_pfx->SetMarkerStyle(21);
nsampls_pfx->SetLineColor(2);
nsampls_pfx->SetMarkerColor(2);
nsampls_pfx->Draw();
TMarker *mk2 = new TMarker(4000.,1.,21);
mk2->SetMarkerColor(2);
mk2->Draw();
Text_t *text = "294 um pitch detector";
TText *t3 = new TText(200,9.,text);
Text_t *text0 = "Average Number of Time bins/Anode/Cluster";
TText *t0 = new TText(200,8.,text0);
Text_t *text2 = "Simulation";
TText *t2 = new TText(4200,0.8,text2);
t3->Draw();
t0->Draw();
t2->Draw();
c1->SaveAs("ITS_tim_check_294.ps");

f.cd();
amplit->ProfileX();
amplit_pfx->SetLineColor(2);
amplit_pfx->SetMaximum(500.);
amplit_pfx->SetMarkerStyle(21);
amplit_pfx->SetMarkerColor(2);
amplit_pfx->Draw();
TMarker *mk2 = new TMarker(500.,50,21);
mk2->SetMarkerColor(2);
mk2->Draw();
Text_t *text = "294 um pitch detector";
TText *t3 = new TText(1000,450.,text);
Text_t *text0 = "Peak Amplitude";
TText *t0 = new TText(1000,400.,text0);
Text_t *text2 = "Simulation";
TText *t2 = new TText(700,40.,text2);
t3->Draw();
t0->Draw();
t2->Draw();
c1->SaveAs("ITS_amp_check_294.ps");

f.cd();
chp->ProfileX();
chp_pfx->SetLineColor(2);
chp_pfx->SetMaximum(2000.);
chp_pfx->SetMarkerStyle(21);
chp_pfx->SetMarkerColor(2);
chp_pfx->Draw();
TMarker *mk2 = new TMarker(500.,200.,21);
mk2->SetMarkerColor(2);
mk2->Draw();
Text_t *text = "294 um pitch detector";
TText *t3 = new TText(300,1800.,text);
Text_t *text0 = "Total Charge";
TText *t0 = new TText(300,1600.,text0);
Text_t *text2 = "Simulation";
TText *t2 = new TText(700,150.,text2);
t3->Draw();
t0->Draw();
t2->Draw();
c1->SaveAs("ITS_cha_check_294.ps");

TH1F *eff = f.Get("rec_vs_time");
TH1F *efh = f.Get("hit_vs_time");
eff->Divide(efh);

TH1F *effn = new TH1F("effn","Efficiency vs. drift path (mm)",18,0.,36.);
effn->SetMinimum(0.9);

for(Int_t i=1;i<=36;i++) {
  Float_t cont = eff->GetBinContent(i);
  i++;
  cont += eff->GetBinContent(i);
  cont /= 2.;
  effn->SetBinContent(i/2,cont);
}

effn->SetXTitle("Drift Path (mm)");
effn->SetYTitle("Reconstruction Efficiency");
effn->SetMaximum(1.2);
effn->SetMinimum(0.6);
effn->SetMarkerStyle(21);
effn->SetLineColor(2);
effn->SetMarkerColor(2);
effn->Draw("p");
TMarker *mk2 = new TMarker(20.,0.7,21);
mk2->SetMarkerColor(2);
mk2->Draw();
Text_t *text = "294 um pitch detector";
TText *t3 = new TText(2.,1.13,text);
Text_t *text2 = "Simulation";
TText *t2 = new TText(22,0.685,text2);
t3->Draw();
t0->Draw();
t2->Draw();
c1->SaveAs("ITS_eff_check_294.ps");

f.cd();
ntotal->ProfileX();
ntotal->SetMarkerStyle(21);
ntotal->SetLineColor(2);
ntotal->SetMarkerColor(2);
ntotal_pfx->SetMaximum(20);
ntotal_pfx->SetMarkerStyle(21);
ntotal_pfx->SetMarkerColor(2);
ntotal_pfx->SetLineColor(2);
ntotal_pfx->SetXTitle("Drift Time (ns)");
ntotal_pfx->SetYTitle("N_total");
ntotal_pfx->Draw();
/*
TH1F *nsi = new TH1F("nsi","nsi",28,0.,7000.);
for(Int_t i=1;i<=28;i++) {
  f->cd();
  Float_t tscal = nsampls_pfx->GetBinContent(i);
  Float_t cont = nanodes_pfx->GetBinContent(i);
  cont *= tscal;
  nsi->SetMaximum(20.);
  nsi->SetBinContent(i,cont);
}
nsi->SetMarkerStyle(21);
nsi->SetMarkerColor(2);
nsi->SetLineColor(2);
nsi->Draw("p,SAME");
*/
TMarker *mk2 = new TMarker(300.,2.,21);
mk2->SetMarkerColor(2);
mk2->Draw();
Text_t *text = "294 um pitch detector";
TText *t3 = new TText(300.,18.,text);
Text_t *text2 = "Simulation";
TText *t2 = new TText(500,1.8,text2);
t3->Draw();
t2->Draw();

c1->SaveAs("ITS_ncl_check_294.ps");

}
