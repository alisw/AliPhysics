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

TH2* bayneseffdiv2D(TH2* numerator, TH2* denominator,Char_t* name) ;
TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name) ;

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
Int_t CutSet = 0;//Defaults: 250 MeV for PHOS and 300 MeV for EMCal
int calculateCorrections(TString filename="rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root", Double_t p0 = 0.366, Double_t p1 = 0.0, char *det = "Emcal", Bool_t isSim = kFALSE)
{
  string inline;
  float value = 0;
  float error = 0;
  int i=0;
  TString detector = det;
  TString emcal = "Emcal";
  c1 = new TCanvas;
  
  TFile *f = TFile::Open(filename, "READ");
  
  TList *l = (TList*)f->Get("out1");
  
//   TTree *primaryTree = (TTree*)l->FindObject(("fPrimaryTree"+detector+"MC").Data());
//   TTree *recTree = (TTree*)l->FindObject(("fEventSummaryTree"+detector+"Rec").Data());
//   TTree *mcTree = (TTree*)l->FindObject(("fEventSummaryTree"+detector+"MC").Data());
//   std::cout << primaryTree << " " << recTree << " " << mcTree << std::endl;
 
  Int_t clusterMult = 0;
  Int_t nChargedNotRemoved = 0;
  Int_t nNeutralNotRemoved = 0;
  Int_t nGammaRemoved = 0;
  Int_t nSecNotRemoved = 0;
  
  Double_t emEtRec = 0;
  Double_t emEtMc = 0;

//   recTree->SetBranchAddress("fNeutralMultiplicity", &clusterMult);
//   mcTree->SetBranchAddress("fChargedNotRemoved", &nChargedNotRemoved);
//   mcTree->SetBranchAddress("fNeutralNotRemoved", &nNeutralNotRemoved);
//   mcTree->SetBranchAddress("fGammaRemoved", &nGammaRemoved);
//   mcTree->SetBranchAddress("fSecondaryNotRemoved", &nSecNotRemoved);
  
  Float_t maxMult = 99.5;
  Int_t nbins = 100;
   if(detector==emcal){
//     //100 seems to be sufficient...
     maxMult = 249.5;
     nbins = 250;
  }
  TH2I *hChargedVsClusterMult = new TH2I("hChVsMult", "hChVsMult", nbins, -0.5, maxMult, nbins, -0.5, maxMult);
  TH2I *hNeutralVsClusterMult = new TH2I("hNeutVsMult", "hNeutVsMult", nbins, -0.5, maxMult, nbins, -0.5, maxMult);
  TH2I *hGammaVsClusterMult = new TH2I("hGammaVsMult", "hGammaVsMult", nbins, -0.5, maxMult, nbins, -0.5, maxMult);
  TH2I *hSecVsClusterMult = new TH2I("hSecVsMult", "hSecVsMult", nbins, -0.5, maxMult, nbins, -0.5, maxMult);

//   Int_t nEvents = mcTree->GetEntriesFast();
//     for(Int_t i = 0; i < nEvents; i++)
//     {
//       mcTree->GetEvent(i);
//       recTree->GetEvent(i);
//       hChargedVsClusterMult->Fill(clusterMult, nChargedNotRemoved);
//       hNeutralVsClusterMult->Fill(clusterMult, nNeutralNotRemoved);
//       hGammaVsClusterMult->Fill(clusterMult, nGammaRemoved);
//       hSecVsClusterMult->Fill(clusterMult, nSecNotRemoved);
//     }
    
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
  //fitcharged.FixParameter(2, p0);
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
  
  TH2F  *fHistHadronDepositsAllMult = l->FindObject("fHistHadronDepositsAllCent");
  TH2F  *fHistHadronDepositsRecoMult = l->FindObject("fHistHadronDepositsRecoCent");
  TH2F *eff2D;
//   if(fHistHadronDepositsRecoMult && fHistHadronDepositsAllMult ){
    eff2D = (TH2F*) bayneseffdiv2D(fHistHadronDepositsRecoMult,fHistHadronDepositsAllMult,"eff2D");
//   }
//   else{
//     cerr<<"Warning!  Did not calculate reconstruction efficiency!!"<<endl;
//     eff2D =  new TH2F("eff2D", "eff2D",200, 0, 10,20,-0.5,19.5);
//   }
  //cor->SetReconstructionEfficiency(eff2D);
  
  TH2F  *fHistGammasGeneratedMult = l->FindObject("fHistGammasGeneratedCent");
//   TH2F  *fHistGammasFoundMult = l->FindObject("fHistGammasFoundCent");
//   TH2F  *fHistGammasFoundOutOfAccMult = l->FindObject("fHistGammasFoundOutOfAccCent");
   TH2F  *fHistGammasFoundMult = l->FindObject("fHistGammasFoundRecoEnergyCent");
   TH2F  *fHistGammasFoundOutOfAccMult = l->FindObject("fHistGammasFoundOutOfAccRecoEnergyCent");
  if(fHistGammasFoundOutOfAccMult){
    //cout<<"I have "<<fHistGammasFoundOutOfAccMult->GetEntries()<<" entries"<<endl;
    fHistGammasFoundMult->Add(fHistGammasFoundOutOfAccMult);
    //cout<<"I have "<<fHistGammasFoundMult->GetEntries()<<" entries"<<endl;
  }
  else{cout<<"fHistGammasFoundOutOfAccCent does not exist!"<<endl;}
  //fHistGammasFoundMult->Add(fHistGammasFoundOutOfAccMult);
  TH2F *gammaEff2D;
  //if(fHistGammasGeneratedMult && fHistGammasFoundMult){
  gammaEff2D = (TH2F*) bayneseffdiv2D((TH2F*)fHistGammasFoundMult,(TH2F*)fHistGammasGeneratedMult,"gammaEff2D");
//   }
//   else{
//     cerr<<"Warning!  Did not calculate reconstruction efficiency!!"<<endl;
//     gammaEff2D =  new TH2F("gammaEff2D", "gammaEff2D",200, 0, 10,20,-0.5,19.5);
//   }

  AliAnalysisEtTrackMatchCorrections *cor = new AliAnalysisEtTrackMatchCorrections(("TmCorrections"+detector).Data(), fitcharged, fitneutral, fitgamma, fitsecondary, *eff2D,
										   meanCharged, meanNeutral, meanGamma, meanSecondary );
  
  TString corrfilename = "allcorrections"+detector+".dat";
  cout<<"Reading "<<corrfilename<<endl;
  ifstream myfile (corrfilename.Data());
  Float_t neutroncorr = 0;
  Float_t hadroncorr = 0;
  Float_t kaoncorr = 0;
  Float_t secondarycorr = 0;
  Float_t minetcorr = 0;
  Float_t junk = 0;
  i=0;
  if (myfile.is_open()){
    while ( myfile.good() )
      {
	getline (myfile,inline);
	istringstream tmp(inline);
// 	tmp >> neutroncorr;
// 	tmp >> hadroncorr;
// 	tmp >> kaoncorr;
// 	tmp >> secondarycorr;

	tmp >> junk;
	tmp >> junk;
	tmp >> junk;
	tmp >> junk;
	tmp >> minetcorr;
	tmp >> junk;//phosMinEtError[i];
	tmp >> junk;//phosNonLinError[i];
	tmp >> neutroncorr;//phosNeutronCorr[i];
	tmp >> junk;//phosNeutronError[i];
	tmp >> hadroncorr;//phosHadronCorr[i];
	tmp >> junk;//phosHadronError[i];
	tmp >> kaoncorr;//phosKaonCorr[i];
	tmp >> junk;//phosKaonError[i];
	tmp >> secondarycorr;//phosSecondaryCorr[i];
	tmp >> junk;//phosSecondaryError[i];
	if(i<20){
	  //cout<<" cb "<<i<<" "<<minetcorr<<" "<<neutroncorr<<" "<<hadroncorr<<" "<<kaoncorr<<" "<<secondarycorr<<endl;
	  cor->SetMinEtCorrection(i,minetcorr);
 	  cor->SetNeutronCorrection(i,neutroncorr);
 	  cor->SetHadronCorrection(i,hadroncorr);
 	  cor->SetKaonCorrection(i,kaoncorr);
 	  cor->SetSecondaryCorrection(i,secondarycorr);
	  i++;

	}
      }
    myfile.close();
  }


//   TString cutstring = "";
//   if(CutSet==1) cutstring = "Cut6";
//   if(CutSet==2) cutstring = "Cut7";
//   TString minetInfileNameShort = "MinEt"+detector+"Short"+cutstring+".dat";
//   cout<<"Reading "<<minetInfileNameShort.Data()<<endl;
//   Float_t minEtErrorShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
//   Float_t minEtCorrShort[10] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0};
//   ifstream myminetfileShort (minetInfileNameShort.Data());
//   if (myminetfileShort.is_open()){
//     while ( myminetfileShort.good() )
//       {
// 	getline (myminetfileShort,inline);
// 	istringstream tmp(inline);
// 	tmp >> value;
// 	tmp >> error;
// 	if(i<10){
// 	  minEtCorrShort[i] = value;
// 	  minEtErrorShort[i] = error;
// 	}
// 	cout<<"min et corr cb "<<i<<" "<< value <<" +/- "<< error <<endl;
// 	i++;
//       }
//     myminetfileShort.close();
//   }



  TF1 *func = generateRecEffFunction(p0, p1);
  AliAnalysisEtRecEffCorrection *recor = new AliAnalysisEtRecEffCorrection(("ReCorrections"+detector).Data(), *func,*gammaEff2D, 1000);
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

TH1* bayneseffdiv(TH1* numerator, TH1* denominator,Char_t* name) 
{
    if(!numerator){
      cerr<<"Error:  numerator does not exist!"<<endl;
      return NULL;
    }
    if(!denominator){
      cerr<<"Error:  denominator does not exist!"<<endl;
      return NULL;
    }
    TH1* result = (TH1*)numerator->Clone(name);
    Int_t nbins = numerator->GetNbinsX();
    for (Int_t ibin=0; ibin<= nbins+1; ++ibin) {
      Double_t numeratorVal = numerator->GetBinContent(ibin);
      Double_t denominatorVal = denominator->GetBinContent(ibin);
      // Check if the errors are right or the thing is scaled
      Double_t numeratorValErr = numerator->GetBinError(ibin);
      if (!(numeratorValErr==0. || numeratorVal ==0.) ) {
	Double_t rescale = numeratorValErr*numeratorValErr/numeratorVal;
	numeratorVal /= rescale;
      }
      Double_t denominatorValErr = denominator->GetBinError(ibin);
      if (!(denominatorValErr==0. || denominatorVal==0. )) {
	Double_t rescale = denominatorValErr*denominatorValErr/denominatorVal;
	denominatorVal /= rescale;
      }
      Double_t quotient = 0.;
      if (denominatorVal!=0.) {
	quotient = numeratorVal/denominatorVal;
      }
      Double_t quotientErr=0;
      quotientErr = TMath::Sqrt(
				(numeratorVal+1.0)/(denominatorVal+2.0)*
				((numeratorVal+2.0)/(denominatorVal+3.0)-(numeratorVal+1.0)/(denominatorVal+2.0)));
      result->SetBinContent(ibin,quotient);
      result->SetBinError(ibin,quotientErr);
      //cout<<"Setting bin "<<ibin<<" to "<<quotient<<" "<<numeratorVal<<"/"<<denominatorVal<<endl;
    }
    return result;
}


TH2* bayneseffdiv2D(TH2* numerator, TH2* denominator,Char_t* name) 
{
  if(!numerator){
    cerr<<"Error:  numerator does not exist!"<<endl;
    return NULL;
  }
  if(!denominator){
    cerr<<"Error:  denominator does not exist!"<<endl;
    return NULL;
  }
  TH2* result = (TH2*)numerator->Clone(name);
  Int_t nbinsX = numerator->GetNbinsX();
  Int_t nbinsY = numerator->GetNbinsY();
  for (Int_t ibin=0; ibin<= nbinsX+1; ++ibin) {
    for (Int_t jbin=0; jbin<= nbinsY+1; ++jbin) {
      Double_t numeratorVal = numerator->GetBinContent(ibin,jbin);
      Double_t denominatorVal = denominator->GetBinContent(ibin,jbin);
      // Check if the errors are right or the thing is scaled
      Double_t numeratorValErr = numerator->GetBinError(ibin,jbin);
      if (!(numeratorValErr==0. || numeratorVal ==0.) ) {
	Double_t rescale = numeratorValErr*numeratorValErr/numeratorVal;
	if(rescale != 0.0) numeratorVal /= rescale;
      }
      Double_t denominatorValErr = denominator->GetBinError(ibin,jbin);
      if (!(denominatorValErr==0. || denominatorVal==0. )) {
	Double_t rescale = denominatorValErr*denominatorValErr/denominatorVal;
	if(rescale != 0.0) denominatorVal /= rescale;
      }
      Double_t quotient = 0.;
      if (denominatorVal!=0.) {
	quotient = numeratorVal/denominatorVal;
      }
      Double_t quotientErr=0;
      if(denominatorVal>0){
	//cerr<<(numeratorVal+1.0)<<"/"<<(denominatorVal+2.0)<<"*"<<"("<<(numeratorVal+2.0)<<"/"<<(denominatorVal+3.0)<<"-"<<(numeratorVal+1.0)<<"/"<<(denominatorVal+2.0)<<")"<<endl;
	float val = (numeratorVal+1.0)/(denominatorVal+2.0)*
	  ((numeratorVal+2.0)/(denominatorVal+3.0)-(numeratorVal+1.0)/(denominatorVal+2.0));
	if(val>0) quotientErr = TMath::Sqrt(val);
      }
      result->SetBinContent(ibin,jbin,quotient);
      result->SetBinError(ibin,jbin,quotientErr);
      //cout<<"Setting bin "<<ibin<<" to "<<quotient<<" "<<numeratorVal<<"/"<<denominatorVal<<endl;
    }
  }
  return result;
}
