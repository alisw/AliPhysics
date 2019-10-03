//////////////////////////////////////////////////////////////////////////////////////
// FindOutliers.C (called by AODQAChecks.C)                                         //
//                                                                                  //
// Written by John Groh                                                             //
//////////////////////////////////////////////////////////////////////////////////////

#include <fstream>
using namespace std;

void FindOutliers(Int_t runs[],
		  const Int_t nRuns,
		  Bool_t useMC,
		  Int_t icut,
		  const Float_t nSigmaCut,
		  TFile*& fout,
		  
		  TH1F*& TPCnsigMeanTrendPion,
		  TH1F*& TPCnsigMeanTrendKaon,
		  TH1F*& TPCnsigMeanTrendProton,
		  TH1F*& TPCnsigSigmaTrendPion,
		  TH1F*& TPCnsigSigmaTrendKaon,
		  TH1F*& TPCnsigSigmaTrendProton,
		  TH1F*& TOFnsigMeanTrendPion,
		  TH1F*& TOFnsigMeanTrendKaon,
		  TH1F*& TOFnsigMeanTrendProton,
		  TH1F*& TOFnsigSigmaTrendPion,
		  TH1F*& TOFnsigSigmaTrendKaon,
		  TH1F*& TOFnsigSigmaTrendProton,
		  
		  TH1F*& IntegRawYieldAll,
		  
		  TH1F*& EfficiencyPiPlus,
		  TH1F*& EfficiencyKPlus,
		  TH1F*& EfficiencyProton,
		  TH1F*& EfficiencyPiMinus,
		  TH1F*& EfficiencyKMinus,
		  TH1F*& EfficiencyAntiproton,
		  
		  TH1F*& MatchEffPos,
		  TH1F*& MatchEffNeg)
{
  Printf("\n\n\n--- Inside function FindOutliers() ---\n\n\n");

  // output to a text file
  char name[100];
  if (useMC) sprintf(name, "results/output_MC_icut=%i_nSigmaCut=%.1f.txt",icut,nSigmaCut);
  else sprintf(name, "results/output_DATA_icut=%i_nSigmaCut=%.1f.txt",icut,nSigmaCut);
  ofstream outFile(name);
  outFile << name << endl;
  outFile << "Info about outlers with useMC = " << useMC << ", icut = " << icut << ", and nSigmaCut = " << nSigmaCut << endl;

  // this array stores whether each run is good or not
  // each run is assumed to be good unless its flag is set otherwise
  Bool_t isRunGood[nRuns];
  for (Int_t irun=0; irun<nRuns; irun++)
    isRunGood[irun] = true;

  // for drawing lines at 0, -3, and 3
  TF1 * fLine1 = new TF1("fLine1","0",0,nRuns);
  TF1 * fLine2 = new TF1("fLine2","-3",0,nRuns);
  TF1 * fLine3 = new TF1("fLine3","3",0,nRuns);
  fLine1->SetLineColor(7);
  fLine2->SetLineColor(kBlack);
  fLine3->SetLineColor(kBlack);
  fLine1->SetLineStyle(7);
  fLine2->SetLineStyle(7);
  fLine3->SetLineStyle(7);
  fLine1->SetLineWidth(1);
  fLine2->SetLineWidth(1);
  fLine3->SetLineWidth(1); 

  //-----------------------------------------
  // nsigma projection fit means and sigmas -
  //-----------------------------------------  
  // because the error bars are really small for data and larger for MC, the normal method doesn't work for data
  if (useMC)
    {
      // fit horizontal lines to the histograms and use the value of the y-intercept as the mean
      TF1 * fitFuncTPCnsigMeanPion = new TF1("fitFuncTPCnsigMeanPion","[0]");
      TF1 * fitFuncTPCnsigMeanKaon = new TF1("fitFuncTPCnsigMeanKaon","[0]");
      TF1 * fitFuncTPCnsigMeanProton = new TF1("fitFuncTPCnsigMeanProton","[0]");
      TF1 * fitFuncTPCnsigSigmaPion = new TF1("fitFuncTPCnsigSigmaPion","[0]");
      TF1 * fitFuncTPCnsigSigmaKaon = new TF1("fitFuncTPCnsigSigmaKaon","[0]");
      TF1 * fitFuncTPCnsigSigmaProton = new TF1("fitFuncTPCnsigSigmaProton","[0]");
      TF1 * fitFuncTOFnsigMeanPion = new TF1("fitFuncTOFnsigMeanPion","[0]");
      TF1 * fitFuncTOFnsigMeanKaon = new TF1("fitFuncTOFnsigMeanKaon","[0]");
      TF1 * fitFuncTOFnsigMeanProton = new TF1("fitFuncTOFnsigMeanProton","[0]");
      TF1 * fitFuncTOFnsigSigmaPion = new TF1("fitFuncTOFnsigSigmaPion","[0]");
      TF1 * fitFuncTOFnsigSigmaKaon = new TF1("fitFuncTOFnsigSigmaKaon","[0]");
      TF1 * fitFuncTOFnsigSigmaProton = new TF1("fitFuncTOFnsigSigmaProton","[0]");
      
      TPCnsigMeanTrendPion->Fit(fitFuncTPCnsigMeanPion,"QN");
      TPCnsigMeanTrendKaon->Fit(fitFuncTPCnsigMeanKaon,"QN");
      TPCnsigMeanTrendProton->Fit(fitFuncTPCnsigMeanProton,"QN");
      TPCnsigSigmaTrendPion->Fit(fitFuncTPCnsigSigmaPion,"QN");
      TPCnsigSigmaTrendKaon->Fit(fitFuncTPCnsigSigmaKaon,"QN");
      TPCnsigSigmaTrendProton->Fit(fitFuncTPCnsigSigmaProton,"QN");
      TOFnsigMeanTrendPion->Fit(fitFuncTOFnsigMeanPion,"QN");
      TOFnsigMeanTrendKaon->Fit(fitFuncTOFnsigMeanKaon,"QN");
      TOFnsigMeanTrendProton->Fit(fitFuncTOFnsigMeanProton,"QN");
      TOFnsigSigmaTrendPion->Fit(fitFuncTOFnsigSigmaPion,"QN");
      TOFnsigSigmaTrendKaon->Fit(fitFuncTOFnsigSigmaKaon,"QN");
      TOFnsigSigmaTrendProton->Fit(fitFuncTOFnsigSigmaProton,"QN");
      
      Float_t meanTPCnsigMeanPion = fitFuncTPCnsigMeanPion->GetParameter(0);
      Float_t meanTPCnsigMeanKaon = fitFuncTPCnsigMeanKaon->GetParameter(0);
      Float_t meanTPCnsigMeanProton = fitFuncTPCnsigMeanProton->GetParameter(0);
      Float_t meanTPCnsigSigmaPion = fitFuncTPCnsigSigmaPion->GetParameter(0);
      Float_t meanTPCnsigSigmaKaon = fitFuncTPCnsigSigmaKaon->GetParameter(0);
      Float_t meanTPCnsigSigmaProton = fitFuncTPCnsigSigmaProton->GetParameter(0);
      Float_t meanTOFnsigMeanPion = fitFuncTOFnsigMeanPion->GetParameter(0);
      Float_t meanTOFnsigMeanKaon = fitFuncTOFnsigMeanKaon->GetParameter(0);
      Float_t meanTOFnsigMeanProton = fitFuncTOFnsigMeanProton->GetParameter(0);
      Float_t meanTOFnsigSigmaPion = fitFuncTOFnsigSigmaPion->GetParameter(0);
      Float_t meanTOFnsigSigmaKaon = fitFuncTOFnsigSigmaKaon->GetParameter(0);
      Float_t meanTOFnsigSigmaProton = fitFuncTOFnsigSigmaProton->GetParameter(0);
      
      // report outliers distributed more than nSigmaCut sigma away from the median
      // also fill histograms to plot the deviation in number of sigmas as a function of the run number
      TH1F * hDevTPCnsigMeanPion = new TH1F("hDevTPCnsigMeanPion","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigMeanKaon = new TH1F("hDevTPCnsigMeanKaon","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigMeanProton = new TH1F("hDevTPCnsigMeanProton","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigSigmaPion = new TH1F("hDevTPCnsigSigmaPion","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigSigmaKaon = new TH1F("hDevTPCnsigSigmaKaon","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigSigmaProton = new TH1F("hDevTPCnsigSigmaProton","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigMeanPion = new TH1F("hDevTOFnsigMeanPion","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigMeanKaon = new TH1F("hDevTOFnsigMeanKaon","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigMeanProton = new TH1F("hDevTOFnsigMeanProton","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigSigmaPion = new TH1F("hDevTOFnsigSigmaPion","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigSigmaKaon = new TH1F("hDevTOFnsigSigmaKaon","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigSigmaProton = new TH1F("hDevTOFnsigSigmaProton","",nRuns,0,nRuns);
      
      cout << "\n\n\n!!! RUNS FOR WHICH AN NSIGMA FIT PARAMETER IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      outFile << "\n\n\n!!! RUNS FOR WHICH AN NSIGMA FIT PARAMETER IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      Int_t nBadNSig = 0;
      for (Int_t irun=0; irun<nRuns; irun++)
	{
	  Float_t devTPCnsigMeanPion = (TPCnsigMeanTrendPion->GetBinContent(irun+1) - meanTPCnsigMeanPion) / TPCnsigMeanTrendPion->GetBinError(irun+1);
	  Float_t devTPCnsigMeanKaon = (TPCnsigMeanTrendKaon->GetBinContent(irun+1) - meanTPCnsigMeanKaon) / TPCnsigMeanTrendKaon->GetBinError(irun+1);
	  Float_t devTPCnsigMeanProton = (TPCnsigMeanTrendProton->GetBinContent(irun+1) - meanTPCnsigMeanProton) / TPCnsigMeanTrendProton->GetBinError(irun+1);
	  Float_t devTPCnsigSigmaPion = (TPCnsigSigmaTrendPion->GetBinContent(irun+1) - meanTPCnsigSigmaPion) / TPCnsigSigmaTrendPion->GetBinError(irun+1);
	  Float_t devTPCnsigSigmaKaon = (TPCnsigSigmaTrendKaon->GetBinContent(irun+1) - meanTPCnsigSigmaKaon) / TPCnsigSigmaTrendKaon->GetBinError(irun+1);
	  Float_t devTPCnsigSigmaProton = (TPCnsigSigmaTrendProton->GetBinContent(irun+1) - meanTPCnsigSigmaProton) / TPCnsigSigmaTrendProton->GetBinError(irun+1);
	  Float_t devTOFnsigMeanPion = (TOFnsigMeanTrendPion->GetBinContent(irun+1) - meanTOFnsigMeanPion) / TOFnsigMeanTrendPion->GetBinError(irun+1);
	  Float_t devTOFnsigMeanKaon = (TOFnsigMeanTrendKaon->GetBinContent(irun+1) - meanTOFnsigMeanKaon) / TOFnsigMeanTrendKaon->GetBinError(irun+1);
	  Float_t devTOFnsigMeanProton = (TOFnsigMeanTrendProton->GetBinContent(irun+1) - meanTOFnsigMeanProton) / TOFnsigMeanTrendProton->GetBinError(irun+1);
	  Float_t devTOFnsigSigmaPion = (TOFnsigSigmaTrendPion->GetBinContent(irun+1) - meanTOFnsigSigmaPion) / TOFnsigSigmaTrendPion->GetBinError(irun+1);
	  Float_t devTOFnsigSigmaKaon = (TOFnsigSigmaTrendKaon->GetBinContent(irun+1) - meanTOFnsigSigmaKaon) / TOFnsigSigmaTrendKaon->GetBinError(irun+1);
	  Float_t devTOFnsigSigmaProton = (TOFnsigSigmaTrendProton->GetBinContent(irun+1) - meanTOFnsigSigmaProton) / TOFnsigSigmaTrendProton->GetBinError(irun+1);
	  
	  hDevTPCnsigMeanPion->SetBinContent(irun+1,devTPCnsigMeanPion);
	  hDevTPCnsigMeanKaon->SetBinContent(irun+1,devTPCnsigMeanKaon);
	  hDevTPCnsigMeanProton->SetBinContent(irun+1,devTPCnsigMeanProton);
	  hDevTPCnsigSigmaPion->SetBinContent(irun+1,devTPCnsigSigmaPion);
	  hDevTPCnsigSigmaKaon->SetBinContent(irun+1,devTPCnsigSigmaKaon);
	  hDevTPCnsigSigmaProton->SetBinContent(irun+1,devTPCnsigSigmaProton);
	  hDevTOFnsigMeanPion->SetBinContent(irun+1,devTOFnsigMeanPion);
	  hDevTOFnsigMeanKaon->SetBinContent(irun+1,devTOFnsigMeanKaon);
	  hDevTOFnsigMeanProton->SetBinContent(irun+1,devTOFnsigMeanProton);
	  hDevTOFnsigSigmaPion->SetBinContent(irun+1,devTOFnsigSigmaPion);
	  hDevTOFnsigSigmaKaon->SetBinContent(irun+1,devTOFnsigSigmaKaon);
	  hDevTOFnsigSigmaProton->SetBinContent(irun+1,devTOFnsigSigmaProton);
	  
	  if (TMath::Abs(devTPCnsigMeanPion) > nSigmaCut || TMath::Abs(devTPCnsigMeanKaon) > nSigmaCut || TMath::Abs(devTPCnsigMeanProton) > nSigmaCut ||
	      TMath::Abs(devTPCnsigSigmaPion) > nSigmaCut || TMath::Abs(devTPCnsigSigmaKaon) > nSigmaCut || TMath::Abs(devTPCnsigSigmaProton) > nSigmaCut ||
	      TMath::Abs(devTOFnsigMeanPion) > nSigmaCut || TMath::Abs(devTOFnsigMeanKaon) > nSigmaCut || TMath::Abs(devTOFnsigMeanProton) > nSigmaCut ||
	      TMath::Abs(devTOFnsigSigmaPion) > nSigmaCut || TMath::Abs(devTOFnsigSigmaKaon) > nSigmaCut || TMath::Abs(devTOFnsigSigmaProton) > nSigmaCut)
	    {
	      cout << runs[irun];
	      outFile << runs[irun] << endl;
	      if (TMath::Abs(devTPCnsigMeanPion) > nSigmaCut)
		{
		  cout << " - TPC Pion Mean";
		  outFile << " - TPC Pion Mean";
		}
	      if (TMath::Abs(devTPCnsigMeanKaon) > nSigmaCut)
		{
		  cout << " - TPC Kaon Mean";
		  outFile << " - TPC Kaon Mean";
		}
	      if (TMath::Abs(devTPCnsigMeanProton) > nSigmaCut)
		{
		  cout << " - TPC Proton Mean";
		  outFile << " - TPC Proton Mean";
		}
	      if (TMath::Abs(devTPCnsigSigmaPion) > nSigmaCut)
		{
		  cout << " - TPC Pion Sigma";
		  outFile << " - TPC Pion Sigma";
		}
	      if (TMath::Abs(devTPCnsigSigmaKaon) > nSigmaCut)
		{
		  cout << " - TPC Kaon Sigma";
		  outFile << " - TPC Kaon Sigma";
		}
	      if (TMath::Abs(devTPCnsigSigmaProton) > nSigmaCut)
		{
		  cout << " - TPC Proton Sigma";
		  outFile << " - TPC Proton Sigma";
		}
	      if (TMath::Abs(devTOFnsigMeanPion) > nSigmaCut)
		{
		  cout << " - TOF Pion Mean";
		  outFile << " - TOF Pion Mean";
		}
	      if (TMath::Abs(devTOFnsigMeanKaon) > nSigmaCut)
		{
		  cout << " - TOF Kaon Mean";
		  outFile << " - TOF Kaon Mean";
		}
	      if (TMath::Abs(devTOFnsigMeanProton) > nSigmaCut)
		{
		  cout << " - TOF Proton Mean";
		  outFile << " - TOF Proton Mean";
		}
	      if (TMath::Abs(devTOFnsigSigmaPion) > nSigmaCut)
		{
		  cout << " - TOF Pion Sigma";
		  outFile << " - TOF Pion Sigma";
		}
	      if (TMath::Abs(devTOFnsigSigmaKaon) > nSigmaCut)
		{
		  cout << " - TOF Kaon Sigma";
		  outFile << " - TOF Kaon Sigma";
		}
	      if (TMath::Abs(devTOFnsigSigmaProton) > nSigmaCut)
		{
		  cout << " - TOF Proton Sigma";
		  outFile << " - TOF Proton Sigma";
		}
	      cout << endl;
	      outFile << endl;
	      nBadNSig++;
	      if (isRunGood[irun]) isRunGood[irun] = kFALSE; // flag run as bad
	    }
	}
      cout << "Number of rejected runs: " << nBadNSig << endl;
      outFile << "Number of rejected runs: " << nBadNSig << endl;
    }
  // for data, I have to fit the projection with a gaussian and cut on nSigmaCut using the sigma from the gaussian
  else
    {
      // project the histos onto the y-axis
      TH1F * TPCnsigMeanPionProj = new TH1F("TPCnsigMeanPionProj","",100,-1,1);
      TH1F * TPCnsigMeanKaonProj = new TH1F("TPCnsigMeanKaonProj","",100,-1,1);
      TH1F * TPCnsigMeanProtonProj = new TH1F("TPCnsigMeanProtonProj","",100,-1,1);
      TH1F * TPCnsigSigmaPionProj = new TH1F("TPCnsigSigmaPionProj","",100,0,2);
      TH1F * TPCnsigSigmaKaonProj = new TH1F("TPCnsigSigmaKaonProj","",100,0,2);
      TH1F * TPCnsigSigmaProtonProj = new TH1F("TPCnsigSigmaProtonProj","",100,0,2);
      TH1F * TOFnsigMeanPionProj = new TH1F("TOFnsigMeanPionProj","",100,-1,1);
      TH1F * TOFnsigMeanKaonProj = new TH1F("TOFnsigMeanKaonProj","",100,-1,1);
      TH1F * TOFnsigMeanProtonProj = new TH1F("TOFnsigMeanProtonProj","",100,-1,1);
      TH1F * TOFnsigSigmaPionProj = new TH1F("TOFnsigSigmaPionProj","",100,0,2);
      TH1F * TOFnsigSigmaKaonProj = new TH1F("TOFnsigSigmaKaonProj","",100,0,2);
      TH1F * TOFnsigSigmaProtonProj = new TH1F("TOFnsigSigmaProtonProj","",100,0,2);
      for (Int_t irun=0; irun<nRuns; irun++)
	{
	  TPCnsigMeanPionProj->Fill(TPCnsigMeanTrendPion->GetBinContent(irun+1));
	  TPCnsigMeanKaonProj->Fill(TPCnsigMeanTrendKaon->GetBinContent(irun+1));
	  TPCnsigMeanProtonProj->Fill(TPCnsigMeanTrendProton->GetBinContent(irun+1));
	  TPCnsigSigmaPionProj->Fill(TPCnsigSigmaTrendPion->GetBinContent(irun+1));
	  TPCnsigSigmaKaonProj->Fill(TPCnsigSigmaTrendKaon->GetBinContent(irun+1));
	  TPCnsigSigmaProtonProj->Fill(TPCnsigSigmaTrendProton->GetBinContent(irun+1));
	  TOFnsigMeanPionProj->Fill(TOFnsigMeanTrendPion->GetBinContent(irun+1));
	  TOFnsigMeanKaonProj->Fill(TOFnsigMeanTrendKaon->GetBinContent(irun+1));
	  TOFnsigMeanProtonProj->Fill(TOFnsigMeanTrendProton->GetBinContent(irun+1));
	  TOFnsigSigmaPionProj->Fill(TOFnsigSigmaTrendPion->GetBinContent(irun+1));
	  TOFnsigSigmaKaonProj->Fill(TOFnsigSigmaTrendKaon->GetBinContent(irun+1));
	  TOFnsigSigmaProtonProj->Fill(TOFnsigSigmaTrendProton->GetBinContent(irun+1));
	}

      // fit with gaussians and reject runs that lie more than nSigmaCut from the mean
      TF1 * fTPCnsigMeanPionFit = new TF1("fTPCnsigMeanPionFit","gaus");
      TF1 * fTPCnsigMeanKaonFit = new TF1("fTPCnsigMeanKaonFit","gaus");
      TF1 * fTPCnsigMeanProtonFit = new TF1("fTPCnsigMeanProtonFit","gaus");
      TF1 * fTPCnsigSigmaPionFit = new TF1("fTPCnsigSigmaPionFit","gaus");
      TF1 * fTPCnsigSigmaKaonFit = new TF1("fTPCnsigSigmaKaonFit","gaus");
      TF1 * fTPCnsigSigmaProtonFit = new TF1("fTPCnsigSigmaProtonFit","gaus");
      TF1 * fTOFnsigMeanPionFit = new TF1("fTOFnsigMeanPionFit","gaus");
      TF1 * fTOFnsigMeanKaonFit = new TF1("fTOFnsigMeanKaonFit","gaus");
      TF1 * fTOFnsigMeanProtonFit = new TF1("fTOFnsigMeanProtonFit","gaus");
      TF1 * fTOFnsigSigmaPionFit = new TF1("fTOFnsigSigmaPionFit","gaus");
      TF1 * fTOFnsigSigmaKaonFit = new TF1("fTOFnsigSigmaKaonFit","gaus");
      TF1 * fTOFnsigSigmaProtonFit = new TF1("fTOFnsigSigmaProtonFit","gaus");
      
      TPCnsigMeanPionProj->Fit(fTPCnsigMeanPionFit,"N","",-1,1);
      TPCnsigMeanKaonProj->Fit(fTPCnsigMeanKaonFit,"N","",-1,1);
      TPCnsigMeanProtonProj->Fit(fTPCnsigMeanProtonFit,"N","",-1,1);
      TPCnsigSigmaPionProj->Fit(fTPCnsigSigmaPionFit,"N","",0,2);
      TPCnsigSigmaKaonProj->Fit(fTPCnsigSigmaKaonFit,"N","",0,2);
      TPCnsigSigmaProtonProj->Fit(fTPCnsigSigmaProtonFit,"N","",0,2);
      TOFnsigMeanPionProj->Fit(fTOFnsigMeanPionFit,"N","",-1,1);
      TOFnsigMeanKaonProj->Fit(fTOFnsigMeanKaonFit,"N","",-1,1);
      TOFnsigMeanProtonProj->Fit(fTOFnsigMeanProtonFit,"N","",-1,1);
      TOFnsigSigmaPionProj->Fit(fTOFnsigSigmaPionFit,"N","",0,2);
      TOFnsigSigmaKaonProj->Fit(fTOFnsigSigmaKaonFit,"N","",0,2);
      TOFnsigSigmaProtonProj->Fit(fTOFnsigSigmaProtonFit,"N","",0,2);

      // means and sigmas from the gaussian fits
      Float_t meanTPCnsigMeanPion = fTPCnsigMeanPionFit->GetParameter(1);
      Float_t meanTPCnsigMeanKaon = fTPCnsigMeanKaonFit->GetParameter(1);
      Float_t meanTPCnsigMeanProton = fTPCnsigMeanProtonFit->GetParameter(1);
      Float_t meanTPCnsigSigmaPion = fTPCnsigSigmaPionFit->GetParameter(1);
      Float_t meanTPCnsigSigmaKaon = fTPCnsigSigmaKaonFit->GetParameter(1);
      Float_t meanTPCnsigSigmaProton = fTPCnsigSigmaProtonFit->GetParameter(1);
      Float_t meanTOFnsigMeanPion = fTOFnsigMeanPionFit->GetParameter(1);
      Float_t meanTOFnsigMeanKaon = fTOFnsigMeanKaonFit->GetParameter(1);
      Float_t meanTOFnsigMeanProton = fTOFnsigMeanProtonFit->GetParameter(1);
      Float_t meanTOFnsigSigmaPion = fTOFnsigSigmaPionFit->GetParameter(1);
      Float_t meanTOFnsigSigmaKaon = fTOFnsigSigmaKaonFit->GetParameter(1);
      Float_t meanTOFnsigSigmaProton = fTOFnsigSigmaProtonFit->GetParameter(1);
      Float_t sigTPCnsigMeanPion = fTPCnsigMeanPionFit->GetParameter(2);
      Float_t sigTPCnsigMeanKaon = fTPCnsigMeanKaonFit->GetParameter(2);
      Float_t sigTPCnsigMeanProton = fTPCnsigMeanProtonFit->GetParameter(2);
      Float_t sigTPCnsigSigmaPion = fTPCnsigSigmaPionFit->GetParameter(2);
      Float_t sigTPCnsigSigmaKaon = fTPCnsigSigmaKaonFit->GetParameter(2);
      Float_t sigTPCnsigSigmaProton = fTPCnsigSigmaProtonFit->GetParameter(2);
      Float_t sigTOFnsigMeanPion = fTOFnsigMeanPionFit->GetParameter(2);
      Float_t sigTOFnsigMeanKaon = fTOFnsigMeanKaonFit->GetParameter(2);
      Float_t sigTOFnsigMeanProton = fTOFnsigMeanProtonFit->GetParameter(2);
      Float_t sigTOFnsigSigmaPion = fTOFnsigSigmaPionFit->GetParameter(2);
      Float_t sigTOFnsigSigmaKaon = fTOFnsigSigmaKaonFit->GetParameter(2);
      Float_t sigTOFnsigSigmaProton = fTOFnsigSigmaProtonFit->GetParameter(2);

      // histograms to plot the deviation as a function of the run number
      TH1F * hDevTPCnsigMeanPion = new TH1F("hDevTPCnsigMeanPion","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigMeanKaon = new TH1F("hDevTPCnsigMeanKaon","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigMeanProton = new TH1F("hDevTPCnsigMeanProton","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigSigmaPion = new TH1F("hDevTPCnsigSigmaPion","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigSigmaKaon = new TH1F("hDevTPCnsigSigmaKaon","",nRuns,0,nRuns);
      TH1F * hDevTPCnsigSigmaProton = new TH1F("hDevTPCnsigSigmaProton","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigMeanPion = new TH1F("hDevTOFnsigMeanPion","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigMeanKaon = new TH1F("hDevTOFnsigMeanKaon","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigMeanProton = new TH1F("hDevTOFnsigMeanProton","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigSigmaPion = new TH1F("hDevTOFnsigSigmaPion","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigSigmaKaon = new TH1F("hDevTOFnsigSigmaKaon","",nRuns,0,nRuns);
      TH1F * hDevTOFnsigSigmaProton = new TH1F("hDevTOFnsigSigmaProton","",nRuns,0,nRuns);
 
      cout << "\n\n\n!!! RUNS FOR WHICH AN NSIGMA FIT PARAMETER IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      outFile << "\n\n\n!!! RUNS FOR WHICH AN NSIGMA FIT PARAMETER IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      Int_t nBadNSig = 0;
      for (Int_t irun=0; irun<nRuns; irun++)
	{
	  Float_t devTPCnsigMeanPion = (TPCnsigMeanTrendPion->GetBinContent(irun+1) - meanTPCnsigMeanPion) / sigTPCnsigMeanPion;
	  Float_t devTPCnsigMeanKaon = (TPCnsigMeanTrendKaon->GetBinContent(irun+1) - meanTPCnsigMeanKaon) / sigTPCnsigMeanKaon;
	  Float_t devTPCnsigMeanProton = (TPCnsigMeanTrendProton->GetBinContent(irun+1) - meanTPCnsigMeanProton) / sigTPCnsigMeanProton;
	  Float_t devTPCnsigSigmaPion = (TPCnsigSigmaTrendPion->GetBinContent(irun+1) - meanTPCnsigSigmaPion) / sigTPCnsigSigmaPion;
	  Float_t devTPCnsigSigmaKaon = (TPCnsigSigmaTrendKaon->GetBinContent(irun+1) - meanTPCnsigSigmaKaon) / sigTPCnsigSigmaKaon;
	  Float_t devTPCnsigSigmaProton = (TPCnsigSigmaTrendProton->GetBinContent(irun+1) - meanTPCnsigSigmaProton) / sigTPCnsigSigmaProton;
	  Float_t devTOFnsigMeanPion = (TOFnsigMeanTrendPion->GetBinContent(irun+1) - meanTOFnsigMeanPion) / sigTOFnsigMeanPion;
	  Float_t devTOFnsigMeanKaon = (TOFnsigMeanTrendKaon->GetBinContent(irun+1) - meanTOFnsigMeanKaon) / sigTOFnsigMeanKaon;
	  Float_t devTOFnsigMeanProton = (TOFnsigMeanTrendProton->GetBinContent(irun+1) - meanTOFnsigMeanProton) / sigTOFnsigMeanProton;
	  Float_t devTOFnsigSigmaPion = (TOFnsigSigmaTrendPion->GetBinContent(irun+1) - meanTOFnsigSigmaPion) / sigTOFnsigSigmaPion;
	  Float_t devTOFnsigSigmaKaon = (TOFnsigSigmaTrendKaon->GetBinContent(irun+1) - meanTOFnsigSigmaKaon) / sigTOFnsigSigmaKaon;
	  Float_t devTOFnsigSigmaProton = (TOFnsigSigmaTrendProton->GetBinContent(irun+1) - meanTOFnsigSigmaProton) / sigTOFnsigSigmaProton;
	  
	  hDevTPCnsigMeanPion->SetBinContent(irun+1,devTPCnsigMeanPion);
	  hDevTPCnsigMeanKaon->SetBinContent(irun+1,devTPCnsigMeanKaon);
	  hDevTPCnsigMeanProton->SetBinContent(irun+1,devTPCnsigMeanProton);
	  hDevTPCnsigSigmaPion->SetBinContent(irun+1,devTPCnsigSigmaPion);
	  hDevTPCnsigSigmaKaon->SetBinContent(irun+1,devTPCnsigSigmaKaon);
	  hDevTPCnsigSigmaProton->SetBinContent(irun+1,devTPCnsigSigmaProton);
	  hDevTOFnsigMeanPion->SetBinContent(irun+1,devTOFnsigMeanPion);
	  hDevTOFnsigMeanKaon->SetBinContent(irun+1,devTOFnsigMeanKaon);
	  hDevTOFnsigMeanProton->SetBinContent(irun+1,devTOFnsigMeanProton);
	  hDevTOFnsigSigmaPion->SetBinContent(irun+1,devTOFnsigSigmaPion);
	  hDevTOFnsigSigmaKaon->SetBinContent(irun+1,devTOFnsigSigmaKaon);
	  hDevTOFnsigSigmaProton->SetBinContent(irun+1,devTOFnsigSigmaProton);
	  
	  if (TMath::Abs(devTPCnsigMeanPion) > nSigmaCut || TMath::Abs(devTPCnsigMeanKaon) > nSigmaCut || TMath::Abs(devTPCnsigMeanProton) > nSigmaCut ||
	      TMath::Abs(devTPCnsigSigmaPion) > nSigmaCut || TMath::Abs(devTPCnsigSigmaKaon) > nSigmaCut || TMath::Abs(devTPCnsigSigmaProton) > nSigmaCut ||
	      TMath::Abs(devTOFnsigMeanPion) > nSigmaCut || TMath::Abs(devTOFnsigMeanKaon) > nSigmaCut || TMath::Abs(devTOFnsigMeanProton) > nSigmaCut ||
	      TMath::Abs(devTOFnsigSigmaPion) > nSigmaCut || TMath::Abs(devTOFnsigSigmaKaon) > nSigmaCut || TMath::Abs(devTOFnsigSigmaProton) > nSigmaCut)
	    {
	      cout << runs[irun];
	      outFile << runs[irun];
	      if (TMath::Abs(devTPCnsigMeanPion) > nSigmaCut)
		{
		  cout << " - TPC Pion Mean";
		  outFile << " - TPC Pion Mean";
		}
	      if (TMath::Abs(devTPCnsigMeanKaon) > nSigmaCut)
		{
		  cout << " - TPC Kaon Mean";
		  outFile << " - TPC Kaon Mean";
		}
	      if (TMath::Abs(devTPCnsigMeanProton) > nSigmaCut)
		{
		  cout << " - TPC Proton Mean";
		  outFile << " - TPC Proton Mean";
		}
	      if (TMath::Abs(devTPCnsigSigmaPion) > nSigmaCut)
		{
		  cout << " - TPC Pion Sigma";
		  outFile << " - TPC Pion Sigma";
		}
	      if (TMath::Abs(devTPCnsigSigmaKaon) > nSigmaCut)
		{
		  cout << " - TPC Kaon Sigma";
		  outFile << " - TPC Kaon Sigma";
		}
	      if (TMath::Abs(devTPCnsigSigmaProton) > nSigmaCut)
		{
		  cout << " - TPC Proton Sigma";
		  outFile << " - TPC Proton Sigma";
		}
	      if (TMath::Abs(devTOFnsigMeanPion) > nSigmaCut)
		{
		  cout << " - TOF Pion Mean";
		  outFile << " - TOF Pion Mean";
		}
	      if (TMath::Abs(devTOFnsigMeanKaon) > nSigmaCut)
		{
		  cout << " - TOF Kaon Mean";
		  outFile << " - TOF Kaon Mean";
		}
	      if (TMath::Abs(devTOFnsigMeanProton) > nSigmaCut)
		{
		  cout << " - TOF Proton Mean";
		  outFile << " - TOF Proton Mean";
		}
	      if (TMath::Abs(devTOFnsigSigmaPion) > nSigmaCut)
		{
		  cout << " - TOF Pion Sigma";
		  outFile << " - TOF Pion Sigma";
		}
	      if (TMath::Abs(devTOFnsigSigmaKaon) > nSigmaCut)
		{
		  cout << " - TOF Kaon Sigma";
		  outFile << " - TOF Kaon Sigma";
		}
	      if (TMath::Abs(devTOFnsigSigmaProton) > nSigmaCut)
		{
		  cout << " - TOF Proton Sigma";
		  outFile << " - TOF Proton Sigma";
		}
	      cout << endl;
	      outFile << endl;
	      nBadNSig++;
	      if (isRunGood[irun]) isRunGood[irun] = kFALSE; // flag run as bad
	    }
	}
      cout << "Number of rejected runs: " << nBadNSig << endl;
      outFile << "Number of rejected runs: " << nBadNSig << endl;
    }

  // draw the deviations as a function of the run number
  TCanvas * cDevNSig = new TCanvas("cDevNSig","cDevNSig",300,150,700,500);
  cDevNSig->Divide(2,2);
  // TPC Means
  cDevNSig->cd(1);
  TH2F * hAxesDevTPCnsigMeans = new TH2F("hAxesDevTPCnsigMeans","",nRuns,0,nRuns,100,-10,10);
  hAxesDevTPCnsigMeans->SetTitle("Mean, TPC");
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesDevTPCnsigMeans->GetXaxis()->SetBinLabel(irun+1,Form("%i",runs[irun]));
  hAxesDevTPCnsigMeans->GetYaxis()->SetTitle("# of #sigma from Mean Value (TPC Means)");
  hAxesDevTPCnsigMeans->SetStats(kFALSE);
  hAxesDevTPCnsigMeans->DrawCopy();
  TLegend * lDevTPCnsigMean = new TLegend(.89,.79,.99,.99);
  hDevTPCnsigMeanPion->SetMarkerColor(Color[0]);
  hDevTPCnsigMeanPion->SetMarkerStyle(Marker[0]);
  hDevTPCnsigMeanPion->SetLineColor(Color[0]);
  hDevTPCnsigMeanPion->SetStats(kFALSE);
  hDevTPCnsigMeanPion->DrawCopy("Psame");
  lDevTPCnsigMean->AddEntry(hDevTPCnsigMeanPion,"#pi^{+}, #pi^{-}","p");
  hDevTPCnsigMeanKaon->SetMarkerColor(Color[1]);
  hDevTPCnsigMeanKaon->SetMarkerStyle(Marker[1]);
  hDevTPCnsigMeanKaon->SetLineColor(Color[1]);
  hDevTPCnsigMeanKaon->SetStats(kFALSE);
  hDevTPCnsigMeanKaon->DrawCopy("Psame");
  lDevTPCnsigMean->AddEntry(hDevTPCnsigMeanKaon,"K^{+}, K^{-}","p");
  hDevTPCnsigMeanProton->SetMarkerColor(Color[2]);
  hDevTPCnsigMeanProton->SetMarkerStyle(Marker[2]);
  hDevTPCnsigMeanProton->SetLineColor(Color[2]);
  hDevTPCnsigMeanProton->SetStats(kFALSE);
  hDevTPCnsigMeanProton->DrawCopy("Psame");
  lDevTPCnsigMean->AddEntry(hDevTPCnsigMeanProton,"p, #bar{p}","p");
  lDevTPCnsigMean->SetFillColor(0);
  lDevTPCnsigMean->DrawClone();
  fLine1->DrawCopy("same");
  fLine2->DrawCopy("same");
  fLine3->DrawCopy("same");
  // TPC Sigmas
  cDevNSig->cd(3);
  TH2F * hAxesDevTPCnsigSigmas = new TH2F("hAxesDevTPCnsigSigmas","",nRuns,0,nRuns,100,-10,10);
  hAxesDevTPCnsigSigmas->SetTitle("Sigma, TPC");
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesDevTPCnsigSigmas->GetXaxis()->SetBinLabel(irun+1,Form("%i",runs[irun]));
  hAxesDevTPCnsigSigmas->GetYaxis()->SetTitle("# of #sigma from Sigma Value (TPC Sigmas)");
  hAxesDevTPCnsigSigmas->SetStats(kFALSE);
  hAxesDevTPCnsigSigmas->DrawCopy();
  TLegend * lDevTPCnsigSigma = new TLegend(.89,.79,.99,.99);
  hDevTPCnsigSigmaPion->SetMarkerColor(Color[0]);
  hDevTPCnsigSigmaPion->SetMarkerStyle(Marker[0]);
  hDevTPCnsigSigmaPion->SetLineColor(Color[0]);
  hDevTPCnsigSigmaPion->SetStats(kFALSE);
  hDevTPCnsigSigmaPion->DrawCopy("Psame");
  lDevTPCnsigSigma->AddEntry(hDevTPCnsigSigmaPion,"#pi^{+}, #pi^{-}","p");
  hDevTPCnsigSigmaKaon->SetMarkerColor(Color[1]);
  hDevTPCnsigSigmaKaon->SetMarkerStyle(Marker[1]);
  hDevTPCnsigSigmaKaon->SetLineColor(Color[1]);
  hDevTPCnsigSigmaKaon->SetStats(kFALSE);
  hDevTPCnsigSigmaKaon->DrawCopy("Psame");
  lDevTPCnsigSigma->AddEntry(hDevTPCnsigSigmaKaon,"K^{+}, K^{-}","p");
  hDevTPCnsigSigmaProton->SetMarkerColor(Color[2]);
  hDevTPCnsigSigmaProton->SetMarkerStyle(Marker[2]);
  hDevTPCnsigSigmaProton->SetLineColor(Color[2]);
  hDevTPCnsigSigmaProton->SetStats(kFALSE);
  hDevTPCnsigSigmaProton->DrawCopy("Psame");
  lDevTPCnsigSigma->AddEntry(hDevTPCnsigSigmaProton,"p, #bar{p}","p");
  lDevTPCnsigSigma->SetFillColor(0);
  lDevTPCnsigSigma->DrawClone();
  fLine1->DrawCopy("same");
  fLine2->DrawCopy("same");
  fLine3->DrawCopy("same");
  // TOF Means
  cDevNSig->cd(2);
  TH2F * hAxesDevTOFnsigMeans = new TH2F("hAxesDevTOFnsigMeans","",nRuns,0,nRuns,100,-10,10);
  hAxesDevTOFnsigMeans->SetTitle("Mean, TOF");
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesDevTOFnsigMeans->GetXaxis()->SetBinLabel(irun+1,Form("%i",runs[irun]));
  hAxesDevTOFnsigMeans->GetYaxis()->SetTitle("# of #sigma from Mean Value (TOF Means)");
  hAxesDevTOFnsigMeans->SetStats(kFALSE);
  hAxesDevTOFnsigMeans->DrawCopy();
  TLegend * lDevTOFnsigMean = new TLegend(.89,.79,.99,.99);
  hDevTOFnsigMeanPion->SetMarkerColor(Color[0]);
  hDevTOFnsigMeanPion->SetMarkerStyle(Marker[0]);
  hDevTOFnsigMeanPion->SetLineColor(Color[0]);
  hDevTOFnsigMeanPion->SetStats(kFALSE);
  hDevTOFnsigMeanPion->DrawCopy("Psame");
  lDevTOFnsigMean->AddEntry(hDevTOFnsigMeanPion,"#pi^{+}, #pi^{-}","p");
  hDevTOFnsigMeanKaon->SetMarkerColor(Color[1]);
  hDevTOFnsigMeanKaon->SetMarkerStyle(Marker[1]);
  hDevTOFnsigMeanKaon->SetLineColor(Color[1]);
  hDevTOFnsigMeanKaon->SetStats(kFALSE);
  hDevTOFnsigMeanKaon->DrawCopy("Psame");
  lDevTOFnsigMean->AddEntry(hDevTOFnsigMeanKaon,"K^{+}, K^{-}","p");
  hDevTOFnsigMeanProton->SetMarkerColor(Color[2]);
  hDevTOFnsigMeanProton->SetMarkerStyle(Marker[2]);
  hDevTOFnsigMeanProton->SetLineColor(Color[2]);
  hDevTOFnsigMeanProton->SetStats(kFALSE);
  hDevTOFnsigMeanProton->DrawCopy("Psame");
  lDevTOFnsigMean->AddEntry(hDevTOFnsigMeanProton,"p, #bar{p}","p");
  lDevTOFnsigMean->SetFillColor(0);
  lDevTOFnsigMean->DrawClone();
  fLine1->DrawCopy("same");
  fLine2->DrawCopy("same");
  fLine3->DrawCopy("same");
  // TOF Sigmas
  cDevNSig->cd(4);
  TH2F * hAxesDevTOFnsigSigmas = new TH2F("hAxesDevTOFnsigSigmas","",nRuns,0,nRuns,100,-10,10);
  hAxesDevTOFnsigSigmas->SetTitle("Sigma, TOF");
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesDevTOFnsigSigmas->GetXaxis()->SetBinLabel(irun+1,Form("%i",runs[irun]));
  hAxesDevTOFnsigSigmas->GetYaxis()->SetTitle("# of #sigma from Sigma Value (TOF Sigmas)");
  hAxesDevTOFnsigSigmas->SetStats(kFALSE);
  hAxesDevTOFnsigSigmas->DrawCopy();
  TLegend * lDevTOFnsigSigma = new TLegend(.89,.79,.99,.99);
  hDevTOFnsigSigmaPion->SetMarkerColor(Color[0]);
  hDevTOFnsigSigmaPion->SetMarkerStyle(Marker[0]);
  hDevTOFnsigSigmaPion->SetLineColor(Color[0]);
  hDevTOFnsigSigmaPion->SetStats(kFALSE);
  hDevTOFnsigSigmaPion->DrawCopy("Psame");
  lDevTOFnsigSigma->AddEntry(hDevTOFnsigSigmaPion,"#pi^{+}, #pi^{-}","p");
  hDevTOFnsigSigmaKaon->SetMarkerColor(Color[1]);
  hDevTOFnsigSigmaKaon->SetMarkerStyle(Marker[1]);
  hDevTOFnsigSigmaKaon->SetLineColor(Color[1]);
  hDevTOFnsigSigmaKaon->SetStats(kFALSE);
  hDevTOFnsigSigmaKaon->DrawCopy("Psame");
  lDevTOFnsigSigma->AddEntry(hDevTOFnsigSigmaKaon,"K^{+}, K^{-}","p");
  hDevTOFnsigSigmaProton->SetMarkerColor(Color[2]);
  hDevTOFnsigSigmaProton->SetMarkerStyle(Marker[2]);
  hDevTOFnsigSigmaProton->SetLineColor(Color[2]);
  hDevTOFnsigSigmaProton->SetStats(kFALSE);
  hDevTOFnsigSigmaProton->DrawCopy("Psame");
  lDevTOFnsigSigma->AddEntry(hDevTOFnsigSigmaProton,"p, #bar{p}","p");
  lDevTOFnsigSigma->SetFillColor(0);
  lDevTOFnsigSigma->DrawClone();
  fLine1->DrawCopy("same");
  fLine2->DrawCopy("same");
  fLine3->DrawCopy("same");
  // write to file
  fout->cd();
  cDevNSig->Write();      
  
  //------------------------------------
  // Raw yields integrated over all Pt -
  //------------------------------------
  // fit horizontal lines to the histograms and use the value of the y-intercept as the mean
  TF1 * fitFuncIntegRawYieldAll = new TF1("fitFuncIntegRawYieldAll","[0]");
  TF1 * fitFuncIntegRawYieldPiPlus = new TF1("fitFuncIntegRawYieldPiPlus","[0]");
  TF1 * fitFuncIntegRawYieldKPlus = new TF1("fitFuncIntegRawYieldKPlus","[0]");
  TF1 * fitFuncIntegRawYieldProton = new TF1("fitFuncIntegRawYieldProton","[0]");
  TF1 * fitFuncIntegRawYieldPiMinus = new TF1("fitFuncIntegRawYieldPiMinus","[0]");
  TF1 * fitFuncIntegRawYieldKMinus = new TF1("fitFuncIntegRawYieldKMinus","[0]");
  TF1 * fitFuncIntegRawYieldAntiproton = new TF1("fitFuncIntegRawYieldAntiproton","[0]");

  IntegRawYieldAll->Fit(fitFuncIntegRawYieldAll,"QN");
  IntegRawYieldPiPlus->Fit(fitFuncIntegRawYieldPiPlus,"QN");
  IntegRawYieldKPlus->Fit(fitFuncIntegRawYieldKPlus,"QN");
  IntegRawYieldProton->Fit(fitFuncIntegRawYieldProton,"QN");
  IntegRawYieldPiMinus->Fit(fitFuncIntegRawYieldPiMinus,"QN");
  IntegRawYieldKMinus->Fit(fitFuncIntegRawYieldKMinus,"QN");
  IntegRawYieldAntiproton->Fit(fitFuncIntegRawYieldAntiproton,"QN");

  Float_t meanIntegRawYieldAll = fitFuncIntegRawYieldAll->GetParameter(0);
  Float_t meanIntegRawYieldPiPlus = fitFuncIntegRawYieldPiPlus->GetParameter(0);
  Float_t meanIntegRawYieldKPlus = fitFuncIntegRawYieldKPlus->GetParameter(0);
  Float_t meanIntegRawYieldProton = fitFuncIntegRawYieldProton->GetParameter(0);
  Float_t meanIntegRawYieldPiMinus = fitFuncIntegRawYieldPiMinus->GetParameter(0);
  Float_t meanIntegRawYieldKMinus = fitFuncIntegRawYieldKMinus->GetParameter(0);
  Float_t meanIntegRawYieldAntiproton = fitFuncIntegRawYieldAntiproton->GetParameter(0);

  // report outliers distributed more than nSigmaCut sigma away from the median
  // also fill histograms to plot the deviation in number of sigmas as a function of the run number
  TH1F * hDevIntegRawYieldAll = new TH1F("hDevIntegRawYieldAll","",nRuns,0,nRuns);
  TH1F * hDevIntegRawYieldPiPlus = new TH1F("hDevIntegRawYieldPiPlus","",nRuns,0,nRuns);
  TH1F * hDevIntegRawYieldKPlus = new TH1F("hDevIntegRawYieldKPlus","",nRuns,0,nRuns);
  TH1F * hDevIntegRawYieldProton = new TH1F("hDevIntegRawYieldProton","",nRuns,0,nRuns);
  TH1F * hDevIntegRawYieldPiMinus = new TH1F("hDevIntegRawYieldPiMinus","",nRuns,0,nRuns);
  TH1F * hDevIntegRawYieldKMinus = new TH1F("hDevIntegRawYieldKMinus","",nRuns,0,nRuns);
  TH1F * hDevIntegRawYieldAntiproton = new TH1F("hDevIntegRawYieldAntiproton","",nRuns,0,nRuns);

  cout << "\n\n\n!!! RUNS FOR WHICH THE INTEGRATED RAW YIELD IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
  outFile << "\n\n\n!!! RUNS FOR WHICH THE INTEGRATED RAW YIELD IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
  Int_t nBadIntegRaw = 0;
  for (Int_t irun=0; irun<nRuns; irun++)
    {
      Float_t devIntegRawYieldAll = (IntegRawYieldAll->GetBinContent(irun+1) - meanIntegRawYieldAll) / IntegRawYieldAll->GetBinError(irun+1);
      Float_t devIntegRawYieldPiPlus = (IntegRawYieldPiPlus->GetBinContent(irun+1) - meanIntegRawYieldPiPlus) / IntegRawYieldPiPlus->GetBinError(irun+1);
      Float_t devIntegRawYieldKPlus = (IntegRawYieldKPlus->GetBinContent(irun+1) - meanIntegRawYieldKPlus) / IntegRawYieldKPlus->GetBinError(irun+1);
      Float_t devIntegRawYieldProton = (IntegRawYieldProton->GetBinContent(irun+1) - meanIntegRawYieldProton) / IntegRawYieldProton->GetBinError(irun+1);
      Float_t devIntegRawYieldPiMinus = (IntegRawYieldPiMinus->GetBinContent(irun+1) - meanIntegRawYieldPiMinus) / IntegRawYieldPiMinus->GetBinError(irun+1);
      Float_t devIntegRawYieldKMinus = (IntegRawYieldKMinus->GetBinContent(irun+1) - meanIntegRawYieldKMinus) / IntegRawYieldKMinus->GetBinError(irun+1);
      Float_t devIntegRawYieldAntiproton = (IntegRawYieldAntiproton->GetBinContent(irun+1) - meanIntegRawYieldAntiproton) / IntegRawYieldAntiproton->GetBinError(irun+1);
      
      hDevIntegRawYieldAll->SetBinContent(irun+1,devIntegRawYieldAll);
      hDevIntegRawYieldPiPlus->SetBinContent(irun+1,devIntegRawYieldPiPlus);
      hDevIntegRawYieldKPlus->SetBinContent(irun+1,devIntegRawYieldKPlus);
      hDevIntegRawYieldProton->SetBinContent(irun+1,devIntegRawYieldProton);
      hDevIntegRawYieldPiMinus->SetBinContent(irun+1,devIntegRawYieldPiMinus);
      hDevIntegRawYieldKMinus->SetBinContent(irun+1,devIntegRawYieldKMinus);
      hDevIntegRawYieldAntiproton->SetBinContent(irun+1,devIntegRawYieldAntiproton);      

      if (TMath::Abs(devIntegRawYieldAll) > nSigmaCut || TMath::Abs(devIntegRawYieldPiPlus) > nSigmaCut || TMath::Abs(devIntegRawYieldKPlus) > nSigmaCut ||
	  TMath::Abs(devIntegRawYieldProton) > nSigmaCut || TMath::Abs(devIntegRawYieldPiMinus) > nSigmaCut || TMath::Abs(devIntegRawYieldKMinus) > nSigmaCut ||
	  TMath::Abs(devIntegRawYieldAntiproton) > nSigmaCut)
	{
	  cout << runs[irun];
	  outFile << runs[irun] << endl;
	  if (TMath::Abs(devIntegRawYieldAll) > nSigmaCut)
	    {
	      cout << " - All";
	      outFile << " - All";
	    }
	  if (TMath::Abs(devIntegRawYieldPiPlus) > nSigmaCut)
	    {
	      cout << " - PiPlus";
	      outFile << " - PiPlus";
	    }
	  if (TMath::Abs(devIntegRawYieldKPlus) > nSigmaCut)
	    {
	      cout << " - KPlus";
	      outFile << " - KPlus";
	    }
	  if (TMath::Abs(devIntegRawYieldProton) > nSigmaCut)
	    {
	      cout << " - Proton";
	      outFile << " - Proton";
	    }
	  if (TMath::Abs(devIntegRawYieldPiMinus) > nSigmaCut)
	    {
	      cout << " - PiMinus";
	      outFile << " - PiMinus";
	    }
	  if (TMath::Abs(devIntegRawYieldKMinus) > nSigmaCut)
	    {
	      cout << " - KMinus";
	      outFile << " - KMinus";
	    }
	  if (TMath::Abs(devIntegRawYieldAntiproton) > nSigmaCut)
	    {
	      cout << " - Antiproton";
	      outFile << " - Antiproton";
	    }
	  cout << endl;
	  outFile << endl;
	  nBadIntegRaw++;
	  if (isRunGood[irun]) isRunGood[irun] = kFALSE; // flag run as bad
	}
    }
  cout << "Number of rejected runs: " << nBadIntegRaw << endl;
  outFile << "Number of rejected runs: " << nBadIntegRaw << endl;
  
  // draw the deviations as a function of the run number
  TCanvas * cDevRawYield = new TCanvas("cDevRawYield","cDevRawYield",350,175,700,500);
  TLegend * lDevIntegRawYield = new TLegend(.89,.79,.99,.99);
  lDevIntegRawYield->SetFillColor(0);

  // all particle spectrum
  for (Int_t irun=0; irun<nRuns; irun++)
    hDevIntegRawYieldAll->GetXaxis()->SetBinLabel(irun+1,Form("%i",runs[irun]));
  hDevIntegRawYieldAll->SetTitle("Raw Yield, Integrated over all p_{T};;# of #sigma from Mean");
  hDevIntegRawYieldAll->GetYaxis()->SetTitleOffset(0.6);
  hDevIntegRawYieldAll->GetYaxis()->CenterTitle();
  hDevIntegRawYieldAll->SetStats(kFALSE);
  hDevIntegRawYieldAll->GetYaxis()->SetRangeUser(-7,7);
  hDevIntegRawYieldAll->SetMarkerStyle(34);
  hDevIntegRawYieldAll->SetMarkerColor(kGreen);
  hDevIntegRawYieldAll->DrawCopy("P");
  lDevIntegRawYield->AddEntry(hDevIntegRawYieldAll,"All Particles","p");
  // individual particles
  hDevIntegRawYieldPiPlus->SetMarkerColor(Color[0]);
  hDevIntegRawYieldPiPlus->SetLineColor(Color[0]);
  hDevIntegRawYieldPiPlus->SetMarkerStyle(Marker[0]);
  hDevIntegRawYieldPiPlus->SetStats(kFALSE);
  hDevIntegRawYieldPiPlus->DrawCopy("Psame");
  lDevIntegRawYield->AddEntry(hDevIntegRawYieldPiPlus,"#pi^{+}","p");
  hDevIntegRawYieldKPlus->SetMarkerColor(Color[1]);
  hDevIntegRawYieldKPlus->SetLineColor(Color[1]);
  hDevIntegRawYieldKPlus->SetMarkerStyle(Marker[1]);
  hDevIntegRawYieldKPlus->SetStats(kFALSE);
  hDevIntegRawYieldKPlus->DrawCopy("Psame");
  lDevIntegRawYield->AddEntry(hDevIntegRawYieldKPlus,"K^{+}","p");
  hDevIntegRawYieldProton->SetMarkerColor(Color[2]);
  hDevIntegRawYieldProton->SetLineColor(Color[2]);
  hDevIntegRawYieldProton->SetMarkerStyle(Marker[2]);
  hDevIntegRawYieldProton->SetStats(kFALSE);
  hDevIntegRawYieldProton->DrawCopy("Psame");
  lDevIntegRawYield->AddEntry(hDevIntegRawYieldProton,"p","p");
  hDevIntegRawYieldPiMinus->SetMarkerColor(Color[0]);
  hDevIntegRawYieldPiMinus->SetLineColor(Color[0]);
  hDevIntegRawYieldPiMinus->SetMarkerStyle(Marker[3]);
  hDevIntegRawYieldPiMinus->SetStats(kFALSE);
  hDevIntegRawYieldPiMinus->DrawCopy("Psame");
  lDevIntegRawYield->AddEntry(hDevIntegRawYieldPiMinus,"#pi^{-}","p");
  hDevIntegRawYieldKMinus->SetMarkerColor(Color[1]);
  hDevIntegRawYieldKMinus->SetLineColor(Color[1]);
  hDevIntegRawYieldKMinus->SetMarkerStyle(Marker[4]);
  hDevIntegRawYieldKMinus->SetStats(kFALSE);
  hDevIntegRawYieldKMinus->DrawCopy("Psame");
  lDevIntegRawYield->AddEntry(hDevIntegRawYieldKMinus,"K^{-}","p");
  hDevIntegRawYieldAntiproton->SetMarkerColor(Color[2]);
  hDevIntegRawYieldAntiproton->SetLineColor(Color[2]);
  hDevIntegRawYieldAntiproton->SetMarkerStyle(Marker[5]);
  hDevIntegRawYieldAntiproton->SetStats(kFALSE);
  hDevIntegRawYieldAntiproton->SetTitle("Raw Yield, Integrated over all p_{T};;# of #sigma from Mean");
  hDevIntegRawYieldAntiproton->GetYaxis()->SetTitleOffset(0.6);
  hDevIntegRawYieldAntiproton->GetYaxis()->CenterTitle();
  hDevIntegRawYieldAntiproton->DrawCopy("Psame");
  lDevIntegRawYield->AddEntry(hDevIntegRawYieldAntiproton,"#bar{p}","p");
  lDevIntegRawYield->DrawClone();
  fLine2->DrawCopy("same");
  fLine3->DrawCopy("same");
  // write to file
  fout->cd();
  cDevRawYield->Write();


  //---------------
  // efficiencies -
  //---------------
  if (useMC)
    {
      // fit horizontal lines to the histograms and use the value of the y-intercept as the mean
      TF1 * fitFuncEffPiPlus = new TF1("fitFuncEffPiPlus","[0]");
      TF1 * fitFuncEffKPlus = new TF1("fitFuncEffKPlus","[0]");
      TF1 * fitFuncEffProton = new TF1("fitFuncEffProton","[0]");
      TF1 * fitFuncEffPiMinus = new TF1("fitFuncEffPiMinus","[0]");
      TF1 * fitFuncEffKMinus = new TF1("fitFuncEffKMinus","[0]");
      TF1 * fitFuncEffAntiproton = new TF1("fitFuncEffAntiproton","[0]");

      EfficiencyPiPlus->Fit(fitFuncEffPiPlus,"QN");
      EfficiencyKPlus->Fit(fitFuncEffKPlus,"QN");
      EfficiencyProton->Fit(fitFuncEffProton,"QN");
      EfficiencyPiMinus->Fit(fitFuncEffPiMinus,"QN");
      EfficiencyKMinus->Fit(fitFuncEffKMinus,"QN");
      EfficiencyAntiproton->Fit(fitFuncEffAntiproton,"QN");

      Float_t meanEffPiPlus = fitFuncEffPiPlus->GetParameter(0);
      Float_t meanEffKPlus = fitFuncEffKPlus->GetParameter(0);
      Float_t meanEffProton = fitFuncEffProton->GetParameter(0);
      Float_t meanEffPiMinus = fitFuncEffPiMinus->GetParameter(0);
      Float_t meanEffKMinus = fitFuncEffKMinus->GetParameter(0);
      Float_t meanEffAntiproton = fitFuncEffAntiproton->GetParameter(0);
      
      // report outliers distributed more than nSigmaCut sigma away from the median
      // also fill histograms to plot the deviation in number of sigmas as a function of the run number
      TH1F * hDevEffPiPlus = new TH1F("hDevEffPiPlus","",nRuns,0,nRuns);
      TH1F * hDevEffKPlus = new TH1F("hDevEffKPlus","",nRuns,0,nRuns);
      TH1F * hDevEffProton = new TH1F("hDevEffProton","",nRuns,0,nRuns);
      TH1F * hDevEffPiMinus = new TH1F("hDevEffPiMinus","",nRuns,0,nRuns);
      TH1F * hDevEffKMinus = new TH1F("hDevEffKMinus","",nRuns,0,nRuns);
      TH1F * hDevEffAntiproton = new TH1F("hDevEffAntiproton","",nRuns,0,nRuns);

      cout << "\n\n\n!!! RUNS FOR WHICH THE EFFICIENCY IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      outFile << "\n\n\n!!! RUNS FOR WHICH THE EFFICIENCY IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      Int_t nBadEff = 0;
      for (Int_t irun=0; irun<nRuns; irun++)
	{
	  // I know this is not very robust, but it's only one run that's giving me a problem b/c the proton efficiency has 0 error
	  if (EfficiencyProton->GetBinError(irun+1) == 0)
	    {
	      cout << runs[irun] << " - PROTON EFFICIENCY HAS ZERO ERROR\n";
	      outFile << runs[irun] << " - PROTON EFFICIENCY HAS ZERO ERROR\n";
	      nBadEff++;
	      continue;
	    }

	  Float_t devEffPiPlus = (EfficiencyPiPlus->GetBinContent(irun+1) - meanEffPiPlus) / EfficiencyPiPlus->GetBinError(irun+1);
	  Float_t devEffKPlus = (EfficiencyKPlus->GetBinContent(irun+1) - meanEffKPlus) / EfficiencyKPlus->GetBinError(irun+1);
	  Float_t devEffProton = (EfficiencyProton->GetBinContent(irun+1) - meanEffProton) / EfficiencyProton->GetBinError(irun+1);
	  Float_t devEffPiMinus = (EfficiencyPiMinus->GetBinContent(irun+1) - meanEffPiMinus) / EfficiencyPiMinus->GetBinError(irun+1);
	  Float_t devEffKMinus = (EfficiencyKMinus->GetBinContent(irun+1) - meanEffKMinus) / EfficiencyKMinus->GetBinError(irun+1);
	  Float_t devEffAntiproton = (EfficiencyAntiproton->GetBinContent(irun+1) - meanEffAntiproton) / EfficiencyAntiproton->GetBinError(irun+1);

	  hDevEffPiPlus->SetBinContent(irun+1,devEffPiPlus);
	  hDevEffKPlus->SetBinContent(irun+1,devEffKPlus);
	  hDevEffProton->SetBinContent(irun+1,devEffProton);
	  hDevEffPiMinus->SetBinContent(irun+1,devEffPiMinus);
	  hDevEffKMinus->SetBinContent(irun+1,devEffKMinus);
	  hDevEffAntiproton->SetBinContent(irun+1,devEffAntiproton);

	  if (TMath::Abs(devEffPiPlus) > nSigmaCut || TMath::Abs(devEffKPlus) > nSigmaCut || TMath::Abs(devEffProton) > nSigmaCut ||
	      TMath::Abs(devEffPiMinus) > nSigmaCut || TMath::Abs(devEffKMinus) > nSigmaCut || TMath::Abs(devEffAntiproton) > nSigmaCut)
	    {
	      cout << runs[irun];
	      outFile << runs[irun];
	      if (TMath::Abs(devEffPiPlus) > nSigmaCut)
		{
		  cout << " - PiPlus";
		  outFile << " - PiPlus";
		}
	      if (TMath::Abs(devEffKPlus) > nSigmaCut)
		{
		  cout << " - KPlus";
		  outFile << " - KPlus";
		}
	      if (TMath::Abs(devEffProton) > nSigmaCut)
		{
		  cout << " - Proton";
		  outFile << " - Proton";
		}
	      if (TMath::Abs(devEffPiMinus) > nSigmaCut)
		{
		  cout << " - PiMinus";
		  outFile << " - PiMinus";
		}
	      if (TMath::Abs(devEffKMinus) > nSigmaCut)
		{
		  cout << " - KMinus";
		  outFile << " - KMinus";
		}
	      if (TMath::Abs(devEffAntiproton) > nSigmaCut)
		{
		  cout << " - Antiproton";
		  outFile << " - Antiproton";
		}
	      cout << endl;
	      outFile << endl;
	      nBadEff++;
	      if (isRunGood[irun]) isRunGood[irun] = kFALSE; // flag run as bad
	    }
	}
      cout << "Number of rejected runs: " << nBadEff << endl;
      outFile << "Number of rejected runs: " << nBadEff << endl; 

      // draw the deviations as a function of the run number
      TCanvas * cDevEff = new TCanvas("cDevEff","cDevEff",400,200,700,500);
      TH2F * hAxesDevEff = new TH2F("hAxesDevEff","",nRuns,0,nRuns,100,-15,15);
      hAxesDevEff->SetTitle("MC Correction Factor at p_{T} = 0.9 GeV/c");
      for (Int_t irun=0; irun<nRuns; irun++)
	hAxesDevEff->GetXaxis()->SetBinLabel(irun+1,Form("%i",runs[irun]));
      hAxesDevEff->GetYaxis()->SetTitle("# of #sigma from mean");
      hAxesDevEff->GetYaxis()->CenterTitle();
      hAxesDevEff->SetStats(kFALSE);
      hAxesDevEff->DrawCopy();
      TLegend * lDevEff = new TLegend(.89,.79,.99,.99);
      hDevEffPiPlus->SetMarkerColor(Color[0]);
      hDevEffPiPlus->SetLineColor(Color[0]);
      hDevEffPiPlus->SetMarkerStyle(Marker[0]);
      hDevEffPiPlus->SetStats(kFALSE);
      hDevEffPiPlus->DrawCopy("Psame");
      lDevEff->AddEntry(hDevEffPiPlus,"#pi^{+}","p");
      hDevEffKPlus->SetMarkerColor(Color[1]);
      hDevEffKPlus->SetLineColor(Color[1]);
      hDevEffKPlus->SetMarkerStyle(Marker[1]);
      hDevEffKPlus->SetStats(kFALSE);
      hDevEffKPlus->DrawCopy("Psame");
      lDevEff->AddEntry(hDevEffKPlus,"K^{+}","p");
      hDevEffProton->SetMarkerColor(Color[2]);
      hDevEffProton->SetLineColor(Color[2]);
      hDevEffProton->SetMarkerStyle(Marker[2]);
      hDevEffProton->SetStats(kFALSE);
      hDevEffProton->DrawCopy("Psame");
      lDevEff->AddEntry(hDevEffProton,"p","p");
      hDevEffPiMinus->SetMarkerColor(Color[0]);
      hDevEffPiMinus->SetLineColor(Color[0]);
      hDevEffPiMinus->SetMarkerStyle(Marker[3]);
      hDevEffPiMinus->SetStats(kFALSE);
      hDevEffPiMinus->DrawCopy("Psame");
      lDevEff->AddEntry(hDevEffPiMinus,"#pi^{-}","p");
      hDevEffKMinus->SetMarkerColor(Color[1]);
      hDevEffKMinus->SetLineColor(Color[1]);
      hDevEffKMinus->SetMarkerStyle(Marker[4]);
      hDevEffKMinus->SetStats(kFALSE);
      hDevEffKMinus->DrawCopy("Psame");
      lDevEff->AddEntry(hDevEffKMinus,"K^{-}","p");
      hDevEffAntiproton->SetMarkerColor(Color[2]);
      hDevEffAntiproton->SetLineColor(Color[2]);
      hDevEffAntiproton->SetMarkerStyle(Marker[5]);
      hDevEffAntiproton->SetStats(kFALSE);
      hDevEffAntiproton->DrawCopy("Psame");
      lDevEff->AddEntry(hDevEffAntiproton,"#bar{p}","p");      
      lDevEff->SetFillColor(0);
      lDevEff->DrawClone();
      fLine1->DrawCopy("same");
      fLine2->DrawCopy("same");
      fLine3->DrawCopy("same");
      // write to file
      fout->cd();
      cDevEff->Write();

    } // end if (useMC)

  
  //------------------------
  // matching efficiencies -
  //------------------------
  // because the error bars are so small for data, I can only use the error bars for MC
  if (useMC)
    {
      // fit horizontal lines to the histograms and use the value of the y-intercept as the mean
      TF1 * fitFuncMatchEffPos = new TF1("fitFuncMatchEffPos","[0]");
      TF1 * fitFuncMatchEffNeg = new TF1("fitFuncMatchEffNeg","[0]");
      
      MatchEffPos->Fit(fitFuncMatchEffPos,"QN");  
      MatchEffNeg->Fit(fitFuncMatchEffNeg,"QN");  
      
      Float_t meanMatchEffPos = fitFuncMatchEffPos->GetParameter(0);
      Float_t meanMatchEffNeg = fitFuncMatchEffNeg->GetParameter(0);
      
      // report outliers distributed more than nSigmaCut sigma away from the median
      // also fill histograms to plot the deviation in number of sigmas as a function of the run number
      TH1F * hDevMatchEffPos = new TH1F("hDevMatchEffPos","",nRuns,0,nRuns);
      TH1F * hDevMatchEffNeg = new TH1F("hDevMatchEffNeg","",nRuns,0,nRuns);
      
      cout << "\n\n\n!!! RUNS FOR WHICH THE MATCHING EFFICIENCY IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      outFile << "\n\n\n!!! RUNS FOR WHICH THE MATCHING EFFICIENCY IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      Int_t nBadMatchEff = 0;
      for (Int_t irun=0; irun<nRuns; irun++)
	{
	  Float_t devMatchEffPos = (MatchEffPos->GetBinContent(irun+1) - meanMatchEffPos) / MatchEffPos->GetBinError(irun+1);
	  Float_t devMatchEffNeg = (MatchEffNeg->GetBinContent(irun+1) - meanMatchEffNeg) / MatchEffNeg->GetBinError(irun+1);

	  hDevMatchEffPos->SetBinContent(irun+1,devMatchEffPos);
	  hDevMatchEffNeg->SetBinContent(irun+1,devMatchEffNeg);
	  
	  if (TMath::Abs(devMatchEffPos) > nSigmaCut || TMath::Abs(devMatchEffNeg) > nSigmaCut)
	    {
	      cout << runs[irun];
	      outFile << runs[irun];
	      if (TMath::Abs(devMatchEffPos) > nSigmaCut)
		{
		  cout << " - Pos";
		  outFile << " - Pos";
		}
	      if (TMath::Abs(devMatchEffNeg) > nSigmaCut)
		{
		  cout << " - Neg";
		  outFile << " - Neg";
		}
	      cout << endl;
	      outFile << endl;
	      nBadMatchEff++;
	      if (isRunGood[irun]) isRunGood[irun] = kFALSE; // flag run as bad
	    }
	}
      cout << "Number of rejected runs: " << nBadMatchEff << endl;
      outFile << "Number of rejected runs: " << nBadMatchEff << endl; 
    }
  else
    {
      // project the matching efficiency histograms onto the y-axis
      TH1F * MatchEffPosProj = new TH1F("MatchEffPosProj","",100,0.54,0.68);
      TH1F * MatchEffNegProj = new TH1F("MatchEffNegProj","",100,0.54,0.68);
      for (Int_t irun=0; irun<nRuns; irun++)
	{
	  MatchEffPosProj->Fill(MatchEffPos->GetBinContent(irun+1));
	  MatchEffNegProj->Fill(MatchEffNeg->GetBinContent(irun+1));
	}
      
        // fit with gaussians and reject runs that lie more than 3 sigma from the mean
      TF1 * fMatchEffPosProjFit = new TF1("fMatchEffPosProjFit","gaus");
      TF1 * fMatchEffNegProjFit = new TF1("fMatchEffNegProjFit","gaus");
      
      MatchEffPosProj->Fit(fMatchEffPosProjFit,"N","",.55,.65);
      MatchEffNegProj->Fit(fMatchEffNegProjFit,"N","",.55,.65);
      
      Float_t meanMatchEffPos = fMatchEffPosProjFit->GetParameter(1);
      Float_t sigMatchEffPos = fMatchEffPosProjFit->GetParameter(2);
      Float_t meanMatchEffNeg = fMatchEffNegProjFit->GetParameter(1);
      Float_t sigMatchEffNeg = fMatchEffNegProjFit->GetParameter(2);
      
      // histograms for plotting the deviation against the run number
      TH1F * hDevMatchEffPos = new TH1F("hDevMatchEffPos","",nRuns,0,nRuns);
      TH1F * hDevMatchEffNeg = new TH1F("hDevMatchEffNeg","",nRuns,0,nRuns);

      cout << "\n\n\n!!! RUNS FOR WHICH THE MATCHING EFFICIENCY IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      outFile << "\n\n\n!!! RUNS FOR WHICH THE MATCHING EFFICIENCY IS MORE THAN " << nSigmaCut << " SIGMA FROM THE MEAN:\n";
      Int_t nBadMatchEff = 0;
      for (Int_t irun=0; irun<nRuns; irun++)
	{
	  Float_t devMatchEffPos = (MatchEffPos->GetBinContent(irun+1) - meanMatchEffPos) / sigMatchEffPos;
	  Float_t devMatchEffNeg = (MatchEffNeg->GetBinContent(irun+1) - meanMatchEffNeg) / sigMatchEffNeg;
	  
	  hDevMatchEffPos->SetBinContent(irun+1,devMatchEffPos);
	  hDevMatchEffNeg->SetBinContent(irun+1,devMatchEffNeg);
	  
	  if (TMath::Abs(devMatchEffPos) > nSigmaCut || TMath::Abs(devMatchEffNeg) > nSigmaCut)
	    {
	      cout << runs[irun];
	      outFile << runs[irun];
	      if (TMath::Abs(devMatchEffPos) > nSigmaCut)
		{
		  cout << " - Pos";
		  outFile << " - Pos";
		}
	      if (TMath::Abs(devMatchEffNeg) > nSigmaCut)
		{
		  cout << " - Neg";
		  outFile << " - Neg";
		}
	      cout << endl;
	      outFile << endl;
	      nBadMatchEff++;
	      if (isRunGood[irun]) isRunGood[irun] = kFALSE; // flag run as bad
	    }
	}
      cout << "Number of rejected runs: " << nBadMatchEff << endl;
      outFile << "Number of rejected runs: " << nBadMatchEff << endl; 
    }

  // draw the deviations as a function of the run number
  TCanvas * cDevMatchEff = new TCanvas("cDevMatchEff","cDevMatchEff",450,225,700,500);
  TH2F * hAxesDevMatchEff = new TH2F("hAxesDevMatchEff","",nRuns,0,nRuns,100,-30,30);
  hAxesDevMatchEff->SetTitle("TOF Matching Efficiency at p_{T} = 0.9 GeV/c;;Deviation from the Mean in # of #sigma");
  hAxesDevMatchEff->GetYaxis()->CenterTitle();
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesDevMatchEff->GetXaxis()->SetBinLabel(irun+1,Form("%i",runs[irun]));
  hAxesDevMatchEff->SetStats(kFALSE);
  hAxesDevMatchEff->DrawCopy();
  TLegend * lDevMatchEff = new TLegend(.89,.79,.99,.99);
  hDevMatchEffPos->SetMarkerColor(kBlue);
  hDevMatchEffPos->SetLineColor(kBlue);
  hDevMatchEffPos->SetMarkerStyle(Marker[0]);
  hDevMatchEffPos->SetStats(kFALSE);
  hDevMatchEffPos->DrawCopy("Psame");
  lDevMatchEff->AddEntry(hDevMatchEffPos,"Positive Particles","p");
  hDevMatchEffNeg->SetMarkerColor(kRed);
  hDevMatchEffNeg->SetLineColor(kRed);
  hDevMatchEffNeg->SetMarkerStyle(Marker[1]);
  hDevMatchEffNeg->SetStats(kFALSE);
  hDevMatchEffNeg->DrawCopy("Psame");
  lDevMatchEff->AddEntry(hDevMatchEffNeg,"Negative Particles","p");
  lDevMatchEff->SetFillColor(0);
  lDevMatchEff->DrawClone();
  //  fLine1->DrawCopy("same");
  fLine2->DrawCopy("same");
  fLine3->DrawCopy("same");
  // write to file
  fout->cd();
  cDevMatchEff->Write();

  // also plot value/mean on a new canvas
  TCanvas * cDataOverMeanMatchEff = new TCanvas("cDataOverMeanMatchEff","cDataOverMeanMatchEff");
  TH1F * hDataOverFitMatchEffPos = (TH1F*)MatchEffPos->Clone();
  hDataOverFitMatchEffPos->Scale(1./meanMatchEffPos);
  TH1F * hDataOverFitMatchEffNeg = (TH1F*)MatchEffNeg->Clone();
  hDataOverFitMatchEffNeg->Scale(1./meanMatchEffNeg);
  TH2F * hAxesDataOverMeanMatchEff = new TH2F("hAxesDataOverMeanMatchEff","Value/Mean for TOF Matching Efficiency at p_{T} = 0.9 GeV/c",nRuns,0,nRuns,100,0.5,1.5);
  for (Int_t irun=0; irun<nRuns; irun++)
    hAxesDataOverMeanMatchEff->GetXaxis()->SetBinLabel(irun+1,Form("%i",runs[irun]));
  hAxesDataOverMeanMatchEff->DrawCopy();
  TLegend * lDataOverMeanMatchEff = new TLegend(.69,.69,.99,.99);
  hDataOverFitMatchEffPos->SetMarkerColor(kBlue);
  hDataOverFitMatchEffPos->SetMarkerStyle(Marker[0]);
  hDataOverFitMatchEffPos->SetMarkerSize(.75);
  hDataOverFitMatchEffPos->DrawCopy("histPsame");
  lDataOverMeanMatchEff->AddEntry(hDataOverFitMatchEffPos,"Pos","p");
  hDataOverFitMatchEffNeg->SetMarkerColor(kRed);
  hDataOverFitMatchEffNeg->SetMarkerStyle(Marker[1]);
  hDataOverFitMatchEffNeg->SetMarkerSize(.75);
  hDataOverFitMatchEffNeg->DrawCopy("histPsame");
  lDataOverMeanMatchEff->AddEntry(hDataOverFitMatchEffNeg,"Neg","p");
  lDataOverMeanMatchEff->DrawClone();
  // write to file
  fout->cd();
  cDataOverMeanMatchEff->Write();


  //--------------------
  // All rejected runs -
  //--------------------
  cout << "\n\n\n!!! LIST OF ALL REJECTED RUNS:\n";
  outFile << "\n\n\n!!! LIST OF ALL REJECTED RUNS:\n";
  Int_t nBadTotal=0;
  for (Int_t irun=0; irun<nRuns; irun++)
    {
      if (!isRunGood[irun])
	{
	  cout << runs[irun] << endl;
	  outFile << runs[irun] << endl;
	  nBadTotal++;
	}
    }
  cout << "Total number of rejected runs: " << nBadTotal << endl;
  outFile << "Total number of rejected runs: " << nBadTotal << endl;
  
  // close the output text file
  outFile.close();

  Printf("\n\n\n--- Leaving function FindOutliers() ---\n\n\n");
}






