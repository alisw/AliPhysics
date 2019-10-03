/***************************************************************
 processFemtoQA.C
 Post Processing of Femto QA task in Analysis QA train
 Author: Maciej Szymanski, maszyman@cern.ch
***************************************************************/

enum collidingSystem {
  PbPb,
  pPb,
  pp
};

Double_t calculateNormalizationFactor(TH1D *num,TH1D *den, Double_t qlo=0.3,Double_t qhi=0.4)
{
  Double_t binlo = num->GetXaxis()->FindFixBin(qlo);
  Double_t binhi = num->GetXaxis()->FindFixBin(qhi);
  Double_t integralNum = num->Integral(binlo, binhi);
  Double_t integralDen = den->Integral(binlo, binhi);
  return integralDen / integralNum;
}

void processFemtoQA(const char *filePath = "AnalysisResults.root",
                    const char *listname = "femtolist",
                    const char *suffix = "png",
                    enum collidingSystem system = PbPb) {

  TFile *_file = new TFile(filePath,"read");
  TList* _femtolist = (TList*)_file->Get(Form("PWG2FEMTO/%s",listname));

  // 1D pion correlation function low kT
  TH1D* numCFlowkT = (TH1D*)_femtolist->FindObject("NumcqinvpimtpcM0kT0");
  TH1D* denCFlowkT = (TH1D*)_femtolist->FindObject("DencqinvpimtpcM0kT0");
  Double_t norm = calculateNormalizationFactor(numCFlowkT,denCFlowkT,0.4,0.45 );
  numCFlowkT->Divide(denCFlowkT);
  numCFlowkT->Scale(norm);
  numCFlowkT->SetXTitle("q_{inv} (GeV/c)");
  numCFlowkT->SetYTitle("C(q_{inv})");
  numCFlowkT->GetXaxis()->SetRangeUser(0,0.5);
  numCFlowkT->GetYaxis()->SetRangeUser(0.5,2.5);

  // 1D pion correlation function high kT
  TH1D* numCFhighkT = (TH1D*)_femtolist->FindObject("NumcqinvpimtpcM0kT3");
  TH1D* denCFhighkt = (TH1D*)_femtolist->FindObject("DencqinvpimtpcM0kT3");
  Double_t norm = calculateNormalizationFactor(numCFhighkT,denCFhighkt,0.4,0.45 );
  numCFhighkT->Divide(denCFhighkt);
  numCFhighkT->Scale(norm);
  numCFhighkT->SetXTitle("q_{inv} (GeV/c)");
  numCFhighkT->SetYTitle("C(q_{inv})");
  numCFhighkT->SetTitle(Form());
  numCFhighkT->GetXaxis()->SetRangeUser(0,0.5);
  numCFhighkT->GetYaxis()->SetRangeUser(0.5,2.5);

  // delta eta - delta phi* low kT
  TH2D* numPhiEtalowkT = (TH2D*)_femtolist->FindObject("NumRadDPhistarEtapimtpcM0kT0");
  TH2D* denPhiEtalowkT = (TH2D*)_femtolist->FindObject("DenRadDPhistarEtapimtpcM0kT0");
  numPhiEtalowkT->Divide(denPhiEtalowkT);
  numPhiEtalowkT->SetXTitle("#Delta #phi*");
  numPhiEtalowkT->SetYTitle("#Delta #eta");
  numPhiEtalowkT->SetTitle(Form());

  // delta eta - delta phi* high kT
  TH2D* numPhiEtahighkT = (TH2D*)_femtolist->FindObject("NumRadDPhistarEtapimtpcM0kT3");
  TH2D* denPhiEtahighkt = (TH2D*)_femtolist->FindObject("DenRadDPhistarEtapimtpcM0kT3");
  numPhiEtahighkT->Divide(denPhiEtahighkt);
  numPhiEtahighkT->SetXTitle("#Delta #phi*");
  numPhiEtahighkT->SetYTitle("#Delta #eta");
  numPhiEtahighkT->SetTitle(Form());

  // qinv vs. separation at TPC entrance low kT
  TH2D* numQinvEtpclowkT = (TH2D*)_femtolist->FindObject("NumDTPCPhistarEtapimtpcM0kT0");
  TH2D* denQinvEtpclowkT = (TH2D*)_femtolist->FindObject("DenDTPCPhistarEtapimtpcM0kT0");
  numQinvEtpclowkT->Divide(denQinvEtpclowkT);
  numQinvEtpclowkT->SetXTitle("q_{inv} (GeV/c)");
  numQinvEtpclowkT->SetYTitle("separation at TPC entrance");
  numQinvEtpclowkT->SetTitle(Form());

  // qinv vs. separation at TPC entrance high kT
  TH2D* numQinvEtpchighkT = (TH2D*)_femtolist->FindObject("NumDTPCPhistarEtapimtpcM0kT3");
  TH2D* denQinvEtpchighkt = (TH2D*)_femtolist->FindObject("DenDTPCPhistarEtapimtpcM0kT3");
  numQinvEtpchighkT->Divide(denQinvEtpchighkt);
  numQinvEtpchighkT->SetXTitle("q_{inv} (GeV/c)");
  numQinvEtpchighkT->SetYTitle("separation at TPC entrance");
  numQinvEtpchighkT->SetTitle(Form());

  if ( system == PbPb ) {
    numCFlowkT->SetTitle("#pi^{-}#pi^{-} 0-10%, 0.2 < k_{T} < 0.3 GeV/c");
    numCFhighkT->SetTitle("#pi^{-}#pi^{-} 0-10%, 0.6 < k_{T} < 0.7 GeV/c");
  }
  else if ( system == pPb ) {
    numCFlowkT->SetTitle("#pi^{-}#pi^{-} 0-20%, 0.2 < k_{T} < 0.3 GeV/c");
    numCFhighkT->SetTitle("#pi^{-}#pi^{-} 0-20%, 0.6 < k_{T} < 0.7 GeV/c");
  }
  else if ( system == pp ) {
    numCFlowkT->SetTitle("#pi^{-}#pi^{-} N_{ch} 50-150, 0.2 < k_{T} < 0.3 GeV/c");
    numCFhighkT->SetTitle("#pi^{-}#pi^{-} N_{ch} 50-150, 0.6 < k_{T} < 0.7 GeV/c");
  }

  gStyle->SetOptStat(0);
  TCanvas* _can = new TCanvas("Femto QA","Femto QA");
  _can->Divide(2,3);
  _can->cd(1);
  numCFlowkT->Draw();
  _can->cd(2);
  numCFhighkT->Draw();
  _can->cd(3);
  numPhiEtalowkT->Draw("colz");
  _can->cd(4);
  numPhiEtahighkT->Draw("colz");
  _can->cd(5);
  numQinvEtpclowkT->Draw("colz");
  _can->cd(6);
  numQinvEtpchighkT->Draw("colz");

  _can->SaveAs(Form("fig_cf_FemtoQA.%s",suffix));

  _file->Close();

  delete _file;
  delete _can;
}
