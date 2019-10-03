void preparepad()
{
  gPad->SetFillColor(0);
  gPad->SetFillStyle(4000);

  gPad->SetTopMargin(0.02);
  //  gPad->SetRightMargin(0.02);
}

void preparehist(TH1D *hist)
{
  hist->SetMarkerSize(0.8);
  hist->SetMarkerStyle(8);
  hist->SetMarkerColor(2);
}

TH2D *plot2ddtpc(TH2D *numdtpc, TH2D *dendtpc, const char *suff)
{
  // Create a 2D CF vs minimum separation distance in the TPC



  int nslices = numdtpc->GetNbinsY();
  int nrange = numdtpc->GetNbinsX();
  char hname[200];
  
  sprintf(hname,"numdtpc%s",suff);
  TH2D *nums = new TH2D(hname, ";q_{inv} [GeV/c];nominal TPC entrance sep - lower cut-off [cm]", 
			numdtpc->GetNbinsX(), numdtpc->GetXaxis()->GetXmin(), numdtpc->GetXaxis()->GetXmax(), 
			numdtpc->GetNbinsY(), numdtpc->GetYaxis()->GetXmin(), numdtpc->GetYaxis()->GetXmax());
  sprintf(hname,"dendtpc%s",suff);
  TH2D *dens = new TH2D(hname, ";q_{inv} [GeV/c];nominal TPC entrance sep - lower cut-off [cm]", 
			numdtpc->GetNbinsX(), numdtpc->GetXaxis()->GetXmin(), numdtpc->GetXaxis()->GetXmax(), 
			numdtpc->GetNbinsY(), numdtpc->GetYaxis()->GetXmin(), numdtpc->GetYaxis()->GetXmax());
  sprintf(hname,"ratdtpc%s",suff);
  TH2D *rats = new TH2D(hname,  ";q_{inv} [GeV/c];nominal TPC entrance sep - lower cut-off [cm]", 
			numdtpc->GetNbinsX(), numdtpc->GetXaxis()->GetXmin(), numdtpc->GetXaxis()->GetXmax(), 
			numdtpc->GetNbinsY(), numdtpc->GetYaxis()->GetXmin(), numdtpc->GetYaxis()->GetXmax());
  
  sprintf(hname,"chi2toone%s",suff);
  TH1D *chi2toone = new TH1D(hname,";nominal TPC entrance sep - lower cut-off [cm];\\chi^2",
			     numdtpc->GetNbinsY(), numdtpc->GetYaxis()->GetXmin(), numdtpc->GetYaxis()->GetXmax());

  char bufname[100];
  for (int iter=0; iter<nslices; iter++) {
    sprintf(bufname, "numbuf%i", iter);
    TH1D *numbuf = numdtpc->ProjectionX(bufname, (iter+1),  nslices, "e");
    sprintf(bufname, "denbuf%i", iter);
    TH1D *denbuf = dendtpc->ProjectionX(bufname, (iter+1),  nslices, "e");
    
    TH1D *ratbuf = new TH1D(*numbuf);
    ratbuf->Divide(denbuf);
    Double_t scale = ratbuf->Integral(nrange-9,nrange)/10.0;

    if (scale > 0.0)
      ratbuf->Scale(1.0/scale);
    
    double chi2 = 0.0;
    for (int ibin = 1; ibin<=ratbuf->GetNbinsX(); ibin++) {
      rats->SetBinContent(ibin, iter+1, ratbuf->GetBinContent(ibin));
      if (ratbuf->GetBinError(ibin) > 0.0) {
	chi2 += (TMath::Power(ratbuf->GetBinContent(ibin)-1.0,2)/
		 TMath::Power(ratbuf->GetBinError(ibin),2));
      }
    }
    chi2toone->SetBinContent(iter+1, chi2);
  }

  rats->GetXaxis()->SetTitleOffset(0.9);
  rats->GetYaxis()->SetTitleOffset(0.95);

  return rats;
}

TH2D * plot2dshare(TH2D *numshare, TH2D *denshare, const char *suff)
{
  // Create a 2D CF vs minimum separation distance in the TPC
  int nrange = numshare->GetNbinsX();
  char hname[200];
    
  sprintf(hname,"numshare%s",suff);
  TH2D *nums = new TH2D(hname, ";q_{inv} [GeV/c];fraction of shared hits - upper cut-off", 
			numshare->GetNbinsX(), numshare->GetXaxis()->GetXmin(), numshare->GetXaxis()->GetXmax(),
			numshare->GetNbinsY(), numshare->GetYaxis()->GetXmin(), numshare->GetYaxis()->GetXmax());

  sprintf(hname,"denshare%s",suff);
  TH2D *dens = new TH2D(hname, ";q_{inv} [GeV/c];fraction of shared hits - upper cut-off", 
			numshare->GetNbinsX(), numshare->GetXaxis()->GetXmin(), numshare->GetXaxis()->GetXmax(),
			numshare->GetNbinsY(), numshare->GetYaxis()->GetXmin(), numshare->GetYaxis()->GetXmax());

  sprintf(hname,"ratshare%s",suff);
  TH2D *rats = new TH2D(hname, ";q_{inv} [GeV/c];fraction of shared hits - upper cut-off", 
			numshare->GetNbinsX(), numshare->GetXaxis()->GetXmin(), numshare->GetXaxis()->GetXmax(),
			numshare->GetNbinsY(), numshare->GetYaxis()->GetXmin(), numshare->GetYaxis()->GetXmax());
  
  char bufname[100];
  for (int iter=0; iter<numshare->GetNbinsY(); iter++) {
    sprintf(bufname, "numbuf%i", iter);
    TH1D *numbuf = numshare->ProjectionX(bufname, 1, iter+1, "e");
    sprintf(bufname, "denbuf%i", iter);
    TH1D *denbuf = denshare->ProjectionX(bufname, 1, iter+1, "e");
    
    ratbuf = new TH1D(*numbuf);
    ratbuf->Divide(denbuf);
    Double_t scale = ratbuf->Integral(nrange-24,nrange)/25.0;
    ratbuf->Scale(1.0/scale);

    for (int ibin = 1; ibin<=ratbuf->GetNbinsX(); ibin++) {
      rats->SetBinContent(ibin, iter+1, ratbuf->GetBinContent(ibin));
    }
  }
  
  rats->GetXaxis()->SetTitleOffset(0.95);
  rats->GetYaxis()->SetTitleOffset(0.9);

  return rats;
}

TH2D *plot2dquality(TH2D *numqual, TH2D *denqual, const char *suff)
{
  // Create a 2D CF vs pair quality 
  int nrange = numqual->GetNbinsX();
  char hname[200];
    
  sprintf(hname,"numqual%s",suff);
  TH2D *nums = new TH2D(hname, ";q_{inv} [GeV/c];pair quality - upper cut-off", 
			numqual->GetNbinsX(), numqual->GetXaxis()->GetXmin(), numqual->GetXaxis()->GetXmax(),
			numqual->GetNbinsY(), numqual->GetYaxis()->GetXmin(), numqual->GetYaxis()->GetXmax());
  sprintf(hname,"denqual%s",suff);
  TH2D *dens = new TH2D(hname, ";q_{inv} [GeV/c];pair quality - upper cut-off", 
			numqual->GetNbinsX(), numqual->GetXaxis()->GetXmin(), numqual->GetXaxis()->GetXmax(),
			numqual->GetNbinsY(), numqual->GetYaxis()->GetXmin(), numqual->GetYaxis()->GetXmax());
  sprintf(hname,"ratqual%s",suff);
  TH2D *rats = new TH2D(hname, ";q_{inv} [GeV/c];pair quality - upper cut-off", 
			numqual->GetNbinsX(), numqual->GetXaxis()->GetXmin(), numqual->GetXaxis()->GetXmax(),
			numqual->GetNbinsY(), numqual->GetYaxis()->GetXmin(), numqual->GetYaxis()->GetXmax());
  
  char bufname[100];
  for (int iter=0; iter<numqual->GetNbinsY(); iter++) {
    sprintf(bufname, "numbuf%i", iter);
    TH1D *numbuf = numqual->ProjectionX(bufname, 1, iter+1, "e");
    sprintf(bufname, "denbuf%i", iter);
    TH1D *denbuf = denqual->ProjectionX(bufname, 1, iter+1, "e");
    
    ratbuf = new TH1D(*numbuf);
    ratbuf->Divide(denbuf);
    Double_t scale = ratbuf->Integral(nrange-24,nrange)/25.0;
    if (scale > 0.0)
      ratbuf->Scale(1.0/scale);

    for (int ibin = 1; ibin<=ratbuf->GetNbinsX(); ibin++) {
      rats->SetBinContent(ibin, iter+1, ratbuf->GetBinContent(ibin));
    }
  }
  
  rats->GetXaxis()->SetTitleOffset(0.95);
  rats->GetYaxis()->SetTitleOffset(0.9);

  return rats;
}

TH1D *getcf(TH1D *numq, TH1D *denq, const char *suff)
{
  char hname[200];

  TH1D *rats = new TH1D(*numq);
  rats->Reset("ICE");
  rats->Divide(numq, denq, 1.0*denq->Integral()/numq->Integral(), 1.0);

  sprintf(hname, "ratqinv%s", suff);
  rats->SetName(hname);
  rats->SetTitle(";q_{inv} [GeV/c];C(q_{inv})");
  preparehist(rats);

  return rats;
}

void FemtoQAPlots(const char *infname)
{
  gStyle->SetStatX(0.8);
  gStyle->SetStatW(0.25);
  gStyle->SetOptStat(11);

  TFile *ifile = new TFile(infname);
  
  gStyle->SetPalette(1,0);

  // Make plot for TPC entrance separation - positive

  TCanvas *candtpcp = new TCanvas("candtpcp","candtpcp",1200,300);
  preparepad();
  candtpcp->Divide(3,1,0.0001,0.0001);
  
  candtpcp->cd(1);  preparepad();
  TH2D *ratdtpcpipstd = plot2ddtpc((TH2D *) gDirectory->Get("NumTPCSepsqqinvcfpipstd"), (TH2D *) gDirectory->Get("DenTPCSepsqqinvcfpipstd"), "pipstd");
  ratdtpcpipstd->Draw("COLZ");

  candtpcp->cd(2);  preparepad();
  TH2D *ratdtpcpipnct = plot2ddtpc((TH2D *) gDirectory->Get("NumTPCSepsqqinvcfpipnct"), (TH2D *) gDirectory->Get("DenTPCSepsqqinvcfpipnct"), "pipnct");
  ratdtpcpipnct->Draw("COLZ");

  candtpcp->cd(3);  preparepad();
  TH2D *ratdtpcpiptpc = plot2ddtpc((TH2D *) gDirectory->Get("NumTPCSepsqqinvcfpiptpc"), (TH2D *) gDirectory->Get("DenTPCSepsqqinvcfpiptpc"), "piptpc");
  ratdtpcpiptpc->Draw("COLZ");

  // Make plot for TPC entrance separation - negative

  TCanvas *candtpcm = new TCanvas("candtpcm","candtpcm",1200,300);
  preparepad();
  candtpcm->Divide(3,1,0.0001,0.0001);
  
  candtpcm->cd(1);  preparepad();
  TH2D *ratdtpcpimstd = plot2ddtpc((TH2D *) gDirectory->Get("NumTPCSepsqqinvcfpimstd"), (TH2D *) gDirectory->Get("DenTPCSepsqqinvcfpimstd"), "pimstd");
  ratdtpcpimstd->Draw("COLZ");

  candtpcm->cd(2);  preparepad();
  TH2D *ratdtpcpimnct = plot2ddtpc((TH2D *) gDirectory->Get("NumTPCSepsqqinvcfpimnct"), (TH2D *) gDirectory->Get("DenTPCSepsqqinvcfpimnct"), "pimnct");
  ratdtpcpimnct->Draw("COLZ");

  candtpcm->cd(3);  preparepad();
  TH2D *ratdtpcpimtpc = plot2ddtpc((TH2D *) gDirectory->Get("NumTPCSepsqqinvcfpimtpc"), (TH2D *) gDirectory->Get("DenTPCSepsqqinvcfpimtpc"), "pimtpc");
  ratdtpcpimtpc->Draw("COLZ");

  // Make plot for pair fraction of shared hits - positive

  TCanvas *cansharep = new TCanvas("cansharep","cansharep",1200,300);
  preparepad();
  cansharep->Divide(3,1,0.0001,0.0001);
  
  cansharep->cd(1);  preparepad();
  TH2D *ratsharepipstd = plot2dshare((TH2D *) gDirectory->Get("NumSharesqqinvcfpipstd"), (TH2D *) gDirectory->Get("DenSharesqqinvcfpipstd"), "pipstd");
  ratsharepipstd->Draw("COLZ");

  cansharep->cd(2);  preparepad();
  TH2D *ratsharepipnct = plot2dshare((TH2D *) gDirectory->Get("NumSharesqqinvcfpipnct"), (TH2D *) gDirectory->Get("DenSharesqqinvcfpipnct"), "pipnct");
  ratsharepipnct->Draw("COLZ");

  cansharep->cd(3);  preparepad();
  TH2D *ratsharepiptpc = plot2dshare((TH2D *) gDirectory->Get("NumSharesqqinvcfpiptpc"), (TH2D *) gDirectory->Get("DenSharesqqinvcfpiptpc"), "piptpc");
  ratsharepiptpc->Draw("COLZ");

  // Make plot for pair fraction of shared hits - negative

  TCanvas *cansharem = new TCanvas("cansharem","cansharem",1200,300);
  preparepad();
  cansharem->Divide(3,1,0.0001,0.0001);
  
  cansharem->cd(1);  preparepad();
  TH2D *ratsharepimstd = plot2dshare((TH2D *) gDirectory->Get("NumSharesqqinvcfpimstd"), (TH2D *) gDirectory->Get("DenSharesqqinvcfpimstd"), "pimstd");
  ratsharepimstd->Draw("COLZ");

  cansharem->cd(2);  preparepad();
  TH2D *ratsharepimnct = plot2dshare((TH2D *) gDirectory->Get("NumSharesqqinvcfpimnct"), (TH2D *) gDirectory->Get("DenSharesqqinvcfpimnct"), "pimnct");
  ratsharepimnct->Draw("COLZ");

  cansharem->cd(3);  preparepad();
  TH2D *ratsharepimtpc = plot2dshare((TH2D *) gDirectory->Get("NumSharesqqinvcfpimtpc"), (TH2D *) gDirectory->Get("DenSharesqqinvcfpimtpc"), "pimtpc");
  ratsharepimtpc->Draw("COLZ");

  // Make signal pair fraction of shared hits - positive

  TCanvas *cansharesignalp = new TCanvas("cansharesignalp","cansharesignalp",1200,300);
  preparepad();
  cansharesignalp->Divide(3,1,0.0001,0.0001);
  
  cansharesignalp->cd(1);  preparepad();
  ((TH2D *) gDirectory->Get("NumSharesqqinvcfpipstd"))->Draw("COLZ");

  cansharesignalp->cd(2);  preparepad();
  ((TH2D *) gDirectory->Get("NumSharesqqinvcfpipnct"))->Draw("COLZ");

  cansharesignalp->cd(3);  preparepad();
  ((TH2D *) gDirectory->Get("NumSharesqqinvcfpiptpc"))->Draw("COLZ");

  // Make signal pair fraction of shared hits - negative

  TCanvas *cansharesignalm = new TCanvas("cansharesignalm","cansharesignalm",1200,300);
  preparepad();
  cansharesignalm->Divide(3,1,0.0001,0.0001);
  
  cansharesignalm->cd(1);  preparepad();
  ((TH2D *) gDirectory->Get("NumSharesqqinvcfpimstd"))->Draw("COLZ");

  cansharesignalm->cd(2);  preparepad();
  ((TH2D *) gDirectory->Get("NumSharesqqinvcfpimnct"))->Draw("COLZ");

  cansharesignalm->cd(3);  preparepad();
  ((TH2D *) gDirectory->Get("NumSharesqqinvcfpimtpc"))->Draw("COLZ");

  // Make background pair fraction of shared hits - positive

  TCanvas *cansharebackp = new TCanvas("cansharebackp","cansharebackp",1200,300);
  preparepad();
  cansharebackp->Divide(3,1,0.0001,0.0001);
  
  cansharebackp->cd(1);  preparepad();
  ((TH2D *) gDirectory->Get("DenSharesqqinvcfpipstd"))->Draw("COLZ");

  cansharebackp->cd(2);  preparepad();
  ((TH2D *) gDirectory->Get("DenSharesqqinvcfpipnct"))->Draw("COLZ");

  cansharebackp->cd(3);  preparepad();
  ((TH2D *) gDirectory->Get("DenSharesqqinvcfpiptpc"))->Draw("COLZ");

  // Make background pair fraction of shared hits - negative

  TCanvas *cansharebackm = new TCanvas("cansharebackm","cansharebackm",1200,300);
  preparepad();
  cansharebackm->Divide(3,1,0.0001,0.0001);
  
  cansharebackm->cd(1);  preparepad();
  ((TH2D *) gDirectory->Get("DenSharesqqinvcfpimstd"))->Draw("COLZ");

  cansharebackm->cd(2);  preparepad();
  ((TH2D *) gDirectory->Get("DenSharesqqinvcfpimnct"))->Draw("COLZ");

  cansharebackm->cd(3);  preparepad();
  ((TH2D *) gDirectory->Get("DenSharesqqinvcfpimtpc"))->Draw("COLZ");

  // Make plot for pair quality - positive

  TCanvas *canqualp = new TCanvas("canqualp","canqualp",1200,300);
  preparepad();
  canqualp->Divide(3,1,0.0001,0.0001);
  
  canqualp->cd(1);  preparepad();
  TH2D *ratqualpipstd = plot2dquality((TH2D *) gDirectory->Get("NumQualitysqqinvcfpipstd"), (TH2D *) gDirectory->Get("DenQualitysqqinvcfpipstd"), "pipstd");
  ratqualpipstd->Draw("COLZ");

  canqualp->cd(2);  preparepad();
  TH2D *ratqualpipnct = plot2dquality((TH2D *) gDirectory->Get("NumQualitysqqinvcfpipnct"), (TH2D *) gDirectory->Get("DenQualitysqqinvcfpipnct"), "pipnct");
  ratqualpipnct->Draw("COLZ");

  canqualp->cd(3);  preparepad();
  TH2D *ratqualpiptpc = plot2dquality((TH2D *) gDirectory->Get("NumQualitysqqinvcfpiptpc"), (TH2D *) gDirectory->Get("DenQualitysqqinvcfpiptpc"), "piptpc");
  ratqualpiptpc->Draw("COLZ");

  // Make plot for pair quality - negative

  TCanvas *canqualm = new TCanvas("canqualm","canqualm",1200,300);
  preparepad();
  canqualm->Divide(3,1,0.0001,0.0001);
  
  canqualm->cd(1);  preparepad();
  TH2D *ratqualpimstd = plot2dquality((TH2D *) gDirectory->Get("NumQualitysqqinvcfpimstd"), (TH2D *) gDirectory->Get("DenQualitysqqinvcfpimstd"), "pimstd");
  ratqualpimstd->Draw("COLZ");

  canqualm->cd(2);  preparepad();
  TH2D *ratqualpimnct = plot2dquality((TH2D *) gDirectory->Get("NumQualitysqqinvcfpimnct"), (TH2D *) gDirectory->Get("DenQualitysqqinvcfpimnct"), "pimnct");
  ratqualpimnct->Draw("COLZ");

  canqualm->cd(3);  preparepad();
  TH2D *ratqualpimtpc = plot2dquality((TH2D *) gDirectory->Get("NumQualitysqqinvcfpimtpc"), (TH2D *) gDirectory->Get("DenQualitysqqinvcfpimtpc"), "pimtpc");
  ratqualpimtpc->Draw("COLZ");

  // Make signal pair quality - positive

  TCanvas *canqualitysignalp = new TCanvas("canqualitysignalp","canqualitysignalp",1200,300);
  preparepad();
  canqualitysignalp->Divide(3,1,0.0001,0.0001);
  
  canqualitysignalp->cd(1);  preparepad();
  ((TH2D *) gDirectory->Get("NumQualitysqqinvcfpipstd"))->Draw("COLZ");

  canqualitysignalp->cd(2);  preparepad();
  ((TH2D *) gDirectory->Get("NumQualitysqqinvcfpipnct"))->Draw("COLZ");

  canqualitysignalp->cd(3);  preparepad();
  ((TH2D *) gDirectory->Get("NumQualitysqqinvcfpiptpc"))->Draw("COLZ");

  // Make signal pair quality - negative

  TCanvas *canqualitysignalm = new TCanvas("canqualitysignalm","canqualitysignalm",1200,300);
  preparepad();
  canqualitysignalm->Divide(3,1,0.0001,0.0001);
  
  canqualitysignalm->cd(1);  preparepad();
  ((TH2D *) gDirectory->Get("NumQualitysqqinvcfpimstd"))->Draw("COLZ");

  canqualitysignalm->cd(2);  preparepad();
  ((TH2D *) gDirectory->Get("NumQualitysqqinvcfpimnct"))->Draw("COLZ");

  canqualitysignalm->cd(3);  preparepad();
  ((TH2D *) gDirectory->Get("NumQualitysqqinvcfpimtpc"))->Draw("COLZ");

  // Make plot for pair quality - positive

  TCanvas *canqinvp = new TCanvas("canqinvp","canqinvp",1200,300);
  preparepad();
  canqinvp->Divide(3,1,0.0001,0.0001);
  
  canqinvp->cd(1); preparepad();
  TH1D *ratqinvpipstd = getcf((TH1D *) gDirectory->Get("Numqinvcfpipstd"), (TH1D *) gDirectory->Get("Denqinvcfpipstd"), "pipstd");
  ratqinvpipstd->Draw();

  canqinvp->cd(2); preparepad();
  TH1D *ratqinvpipnct = getcf((TH1D *) gDirectory->Get("Numqinvcfpipnct"), (TH1D *) gDirectory->Get("Denqinvcfpipnct"), "pipnct");
  ratqinvpipnct->Draw("");

  canqinvp->cd(3); preparepad();
  TH1D *ratqinvpiptpc = getcf((TH1D *) gDirectory->Get("Numqinvcfpiptpc"), (TH1D *) gDirectory->Get("Denqinvcfpiptpc"), "piptpc");
  ratqinvpiptpc->Draw("");

  // Make plot for pair quality - negative

  TCanvas *canqinvm = new TCanvas("canqinvm","canqinvm",1200,300);
  preparepad();
  canqinvm->Divide(3,1,0.0001,0.0001);
  
  canqinvm->cd(1); preparepad();
  TH1D *ratqinvpimstd = getcf((TH1D *) gDirectory->Get("Numqinvcfpimstd"), (TH1D *) gDirectory->Get("Denqinvcfpimstd"), "pimstd");
  ratqinvpimstd->Draw("");

  canqinvm->cd(2); preparepad();
  TH1D *ratqinvpimnct = getcf((TH1D *) gDirectory->Get("Numqinvcfpimnct"), (TH1D *) gDirectory->Get("Denqinvcfpimnct"), "pimnct");
  ratqinvpimnct->Draw("");

  canqinvm->cd(3); preparepad();
  TH1D *ratqinvpimtpc = getcf((TH1D *) gDirectory->Get("Numqinvcfpimtpc"), (TH1D *) gDirectory->Get("Denqinvcfpimtpc"), "pimtpc");
  ratqinvpimtpc->Draw("");

  
}
