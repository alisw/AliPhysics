//Oystein Djuvsland
//macro for plotting correction factors for EM et for use with AliAnalysisEt and writes them to a file
#include "TTree.h"
#include "TFile.h"
#include <TList.h>
#include <Rtypes.h>
#include "TCanvas.h"
#include <TH2I.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include "TF1.h"
#include <iostream>
#include <AliAnalysisEtTrackMatchCorrections.h>
#include <AliAnalysisEtRecEffCorrection.h>

TCanvas *c1 = 0;


//creates an empty set of correction factors for debugging purposes
TF1 * generateRecEffFunction(Double_t p0, Double_t p1);
int createDummy(Double_t p0 = 1.0, Double_t p1 = 0.0, char *det = "Phos")
{
  TFile *outfile = TFile::Open("calocorrections.root","RECREATE");
  TF1 fitneutral("fitneutral","0", 0, 100);
  TF1 fitcharged("fitcharged","0", 0, 100);
  TF1 fitgamma("fitgamma","0", 0, 100);
  TF1 fitsecondary("fitsecondary","0", 0, 100);
  AliAnalysisEtTrackMatchCorrections *cor = new AliAnalysisEtTrackMatchCorrections(Form("TmCorrections%s",det), fitcharged, fitneutral, fitgamma, fitsecondary,0,0,0,0);
  TF1 *func = generateRecEffFunction(p0, p1);
  AliAnalysisEtRecEffCorrection *recor = new AliAnalysisEtRecEffCorrection(Form("ReCorrections%s",det), *func, 1000);

  cor->Write();
  recor->Write();
  outfile->Close();


}
//p0 = correction factors for efficiency from track matching
//p1 = correction factors for efficiency from track matching
//determined from fit of efficiency from track matching, 0.366 from simulations, fit as a function of energy
//p0 is the constant
//p1 is linear term
int calculateCorrections(TString filename="rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root", Double_t p0 = 0.366, Double_t p1 = 0.0, char *det = "Emcal", Bool_t isSim = kFALSE)
{
  TString detector = det;
  TString emcal = "Emcal";
  c1 = new TCanvas;
  
  TFile *f = TFile::Open(filename, "READ");
  
  TList *l = (TList*)f->Get("out1");
  
  TTree *primaryTree = (TTree*)l->FindObject(("fPrimaryTree"+detector+"MC").Data());
  TTree *recTree = (TTree*)l->FindObject(("fEventSummaryTree"+detector+"Rec").Data());
  TTree *mcTree = (TTree*)l->FindObject(("fEventSummaryTree"+detector+"MC").Data());
  std::cout << primaryTree << " " << recTree << " " << mcTree << std::endl;
 
  Int_t clusterMult = 0;
  Int_t nChargedNotRemoved = 0;
  Int_t nNeutralNotRemoved = 0;
  Int_t nGammaRemoved = 0;
  Int_t nSecNotRemoved = 0;
  
  Double_t emEtRec = 0;
  Double_t emEtMc = 0;

  recTree->SetBranchAddress("fNeutralMultiplicity", &clusterMult);
  mcTree->SetBranchAddress("fChargedNotRemoved", &nChargedNotRemoved);
  mcTree->SetBranchAddress("fNeutralNotRemoved", &nNeutralNotRemoved);
  mcTree->SetBranchAddress("fGammaRemoved", &nGammaRemoved);
  mcTree->SetBranchAddress("fSecondaryNotRemoved", &nSecNotRemoved);
  
  Float_t maxMult = 99.5;
  Int_t nbins = 100;
//   if(detector==emcal){
//     //100 seems to be sufficient...
//     maxMult = 299.5;
//     nbins = 300;
//  }
  TH2I *hChargedVsClusterMult = new TH2I("hChVsMult", "hChVsMult", nbins, -0.5, maxMult, nbins, -0.5, maxMult);
  TH2I *hNeutralVsClusterMult = new TH2I("hNeutVsMult", "hNeutVsMult", nbins, -0.5, maxMult, nbins, -0.5, maxMult);
  TH2I *hGammaVsClusterMult = new TH2I("hGammaVsMult", "hGammaVsMult", nbins, -0.5, maxMult, nbins, -0.5, maxMult);
  TH2I *hSecVsClusterMult = new TH2I("hSecVsMult", "hSecVsMult", nbins, -0.5, maxMult, nbins, -0.5, maxMult);

  Int_t nEvents = mcTree->GetEntriesFast();
    for(Int_t i = 0; i < nEvents; i++)
    {
      mcTree->GetEvent(i);
      recTree->GetEvent(i);
      hChargedVsClusterMult->Fill(clusterMult, nChargedNotRemoved);
      hNeutralVsClusterMult->Fill(clusterMult, nNeutralNotRemoved);
      hGammaVsClusterMult->Fill(clusterMult, nGammaRemoved);
      hSecVsClusterMult->Fill(clusterMult, nSecNotRemoved);
    }
    
  c1->Divide(2,2);
  c1->cd(1);
  TString title = "Charged particles not removed";
  TString xtitle = "Cluster Multiplicity";
  TString ytitle = "N_{ch}";
  hChargedVsClusterMult->SetTitle(title);
  hChargedVsClusterMult->GetYaxis()->SetTitle(ytitle);
  hChargedVsClusterMult->GetXaxis()->SetTitle(xtitle);
  hChargedVsClusterMult->SetStats(0);
  hChargedVsClusterMult->Draw();
  TProfile *chProf = hChargedVsClusterMult->ProfileX();
  chProf->SetStats(0);
  chProf->Draw("SAME");
  //TF1 fitcharged("fitcharged","([0] + [1]*x)*(0.48/([2] + [3]*[2]))", 0, 100);
  TF1 fitcharged("fitcharged","pol2", 0, 100);//fit of number of charged tracks vs detector multiplicity
  //if straight line track matching roughly not dependent on centrality
//   fitcharged.FixParameter(2, p0);
//   fitcharged.FixParameter(3, p1);
  TFitResultPtr chRes = chProf->Fit(&fitcharged,"S");
  TArrayD ch;
  if(!chRes)
  {
    ch = TArrayD(chRes->NPar(),chRes->GetParams());
  }
  else
  {
    std::cout << "Could not extract charged contribution params" << std::endl;
  }
  c1->cd(2);
  title = "Neutral particles not removed";
  ytitle = "N_{neutral}";
  hNeutralVsClusterMult->SetTitle(title);
  hNeutralVsClusterMult->GetXaxis()->SetTitle(xtitle);
  hNeutralVsClusterMult->GetYaxis()->SetTitle(ytitle);
  hNeutralVsClusterMult->SetStats(0);
  hNeutralVsClusterMult->Draw();
  TProfile *neuProf = hNeutralVsClusterMult->ProfileX();
  neuProf->SetStats(0);
  neuProf->Draw("SAME");
  TF1 fitneutral("fitneutral","pol2", 0, 100);//fit of number of neutral particles in calo that we call background in calo vs detector multiplicity
  //may include K0S
  TFitResultPtr neuRes = neuProf->Fit(&fitneutral,"S");
  TArrayD neu;
  if(!neuRes)
  {
    neu = TArrayD(neuRes->NPar(),neuRes->GetParams());
  }
  else
  {
    std::cout << "Could not extract charged contribution params" << std::endl;
  }
  c1->cd(3);
  
  title = "Gammas removed";
  ytitle = "N_{#gamma}";
  hGammaVsClusterMult->SetTitle(title);
  hGammaVsClusterMult->GetYaxis()->SetTitle(ytitle);
  hGammaVsClusterMult->GetXaxis()->SetTitle(xtitle);
  hGammaVsClusterMult->SetStats(0);
  hGammaVsClusterMult->Draw();
  TProfile *gamProf = hGammaVsClusterMult->ProfileX();
  gamProf->SetStats(0);
  gamProf->Draw("SAME");
  TF1 fitgamma("fitgamma","pol2", 0, 100);//fit of number of gammas removed erroneously vs detector multiplicity
  TFitResultPtr gammaRes = gamProf->Fit(&fitgamma,"S");
  TArrayD gamma;
  if(!gammaRes)
  {
    gamma = TArrayD(gammaRes->NPar(),gammaRes->GetParams());
  }
  else
  {
    std::cout << "Could not extract charged contribution params" << std::endl;
  }
  
  c1->cd(4);
  
  title = "Secondaries not removed";
  ytitle = "N_{sec}";
  hSecVsClusterMult->SetTitle(title);
  hSecVsClusterMult->GetXaxis()->SetTitle(xtitle);
  hSecVsClusterMult->GetYaxis()->SetTitle(ytitle);
  hSecVsClusterMult->SetStats(0);
  hSecVsClusterMult->Draw();
  TProfile *secProf = hSecVsClusterMult->ProfileX();
  secProf->SetStats(0);
  secProf->Draw("SAME");
  TF1 fitsecondary("fitsecondary","pol2", 0, 100);//fit of number of secondary particles that leave deposits in calo vs detector multiplicity
  TFitResultPtr secRes = secProf->Fit(&fitsecondary,"S");
  TArrayD sec;
  if(!secRes)
  {
    sec = TArrayD(secRes->NPar(),secRes->GetParams());
  }
  else
  {
    std::cout << "Could not extract charged contribution params" << std::endl;
  }
  //ugly hack for getting the energy
  //changing number of particles to energy of particles by multiplying by the <energy> and dividing by the efficiency for the average energy
  //average energy from each of these:  in excel file
  Double_t meanCharged = 0.48/(p0 + p1*0.48);
  Double_t meanNeutral = 0.53/(p0 + p1*0.53);
  Double_t meanGamma = 0.51/(p0 + p1*0.51);
  Double_t meanSecondary = meanGamma; 
  
  AliAnalysisEtTrackMatchCorrections *cor = new AliAnalysisEtTrackMatchCorrections(("TmCorrections"+detector).Data(), fitcharged, fitneutral, fitgamma, fitsecondary, 
										   meanCharged, meanNeutral, meanGamma, meanSecondary );
  
  TF1 *func = generateRecEffFunction(p0, p1);
  AliAnalysisEtRecEffCorrection *recor = new AliAnalysisEtRecEffCorrection(("ReCorrections"+detector).Data(), *func, 1000);
  TFile *outfile = TFile::Open("calocorrections.root","RECREATE");
  cor->Write();
  recor->Write();
  return 0;
}


TF1* generateRecEffFunction(Double_t p0, Double_t p1)
{
  TF1 *f = new TF1("receff", "x/([0] + x*[1])", 0, 200);
  
  Double_t params[2] = {p0, p1};
  f->SetParameters(params);
  return f;
}
