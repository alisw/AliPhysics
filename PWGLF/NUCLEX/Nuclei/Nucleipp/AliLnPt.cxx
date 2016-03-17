/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// reconstructed pt
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TString.h>
#include <TROOT.h>
#include <TString.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooMsgService.h>
#include <TMath.h>

#include "AliLnPt.h"
#include "B2.h"

ClassImp(AliLnPt)

AliLnPt::AliLnPt(const TString& particle, Double_t trigEff, const TString& inputFilename, const TString& outputFilename, const TString& otag, const TString& corrFilename, const TString& corrtag)
: TObject()
, fParticle(particle)
, fTrigEff(trigEff)
, fInputFilename(inputFilename)
, fOutputFilename(outputFilename)
, fOutputTag(otag)
, fCorrFilename(corrFilename)
, fCorrTag(corrtag)
, fPtMin(0.5)
, fPtMax(3.0)
, fPid(0)
, fSecondaries(1)
, fEfficiency(1)
, fIsOnlyGen(0)
, fPtPid(1.)
, fBkgMin(2.2)
, fBkgMax(5.)
, fIntMin(2.)
, fIntMax(6.)
, fMakeStats(1)
, fMCtoINEL(0)
, fVtxFactor(1)
, fFitFrac(1)
, fFdwnCorr(1)
, fSameFdwn(0)
, fPidProc(kMassSquared)
, fPidEff(1)
, fDebugLevel(0)
{
//
// constructor
//
	TH1::SetDefaultSumw2(); // switch on histogram errors
	
	if(fDebugLevel < 3)
	{
		// disable verbose in RooFit
		RooMsgService::instance().setSilentMode(1);
		RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
	}
}

AliLnPt::~AliLnPt()
{
//
// destructor
//
}

Int_t AliLnPt::Exec()
{
//
// Get the true pt distribution
//
	TFile* finput = new TFile(fInputFilename.Data(), "read");
	if (finput->IsZombie()) exit(1);
	
	TFile* foutput = new TFile(fOutputFilename.Data(), "recreate");
	if (foutput->IsZombie()) exit(1);
	
	TFile* fcorr = 0;
	TFile* fdebug = 0;
	
	TString species = fParticle;
	species.ReplaceAll("Anti","");
	
	TH1D* hStats = FindObj<TH1D>(finput, species + "_Stats");
	
	TH1D* hCorrPt = 0;
	
	if(fIsOnlyGen)
	{
		hCorrPt = dynamic_cast<TH1D*>( (FindObj<TH1D>(finput, fParticle + "_Gen_PhS_Prim_Pt"))->Clone( Form("%s_Cor_Pt", fParticle.Data())) );
	}
	else
	{
		fcorr = new TFile(fCorrFilename.Data(),"read");
		if (fcorr->IsZombie()) exit(1);
		
		fdebug =  new TFile(Form("debug-%s",fOutputFilename.Data()),"recreate");
		if (fdebug->IsZombie()) exit(1);
		
		fdebug->mkdir(fOutputTag.Data());
		
		// pointers to the required histograms
		
		TH1D* hOrigPt  = FindObj<TH1D>(finput, fParticle + "_PID_Pt");
		
		TH2D* hPidPt = 0;
		if(fPidProc == kTimeOfFlight)        hPidPt = FindObj<TH2D>(finput, fParticle + "_PID_DTime_Pt");
		else if(fPidProc == kMassSquared)    hPidPt = FindObj<TH2D>(finput, fParticle + "_PID_M2_Pt");
		else if(fPidProc == kMassDifference) hPidPt = FindObj<TH2D>(finput, fParticle + "_PID_DM2_Pt");
		else fPid = 0;
		
		// corrections
		TH1D* hFracMatPt   = FindObj<TH1D>(fcorr, fCorrTag, fParticle + "_Frac_Mat_Pt");
		TH1D* hFracFdwnPt  = FindObj<TH1D>(fcorr, fCorrTag, fParticle + "_Frac_Fdwn_Pt");
		TF1* fncFracMatPt  = FindObj<TF1>(fcorr, fCorrTag, fParticle + "_Frac_Mat_Fit_Pt");
		TF1* fncFracFdwnPt = FindObj<TF1>(fcorr, fCorrTag, fParticle + "_Frac_Fdwn_Fit_Pt");
		
		if(fSameFdwn)
		{
			hFracFdwnPt = FindObj<TH1D>(fcorr, fCorrTag, species + "_Frac_Fdwn_Pt");
			fncFracFdwnPt = FindObj<TF1>(fcorr, fCorrTag, species + "_Frac_Fdwn_Fit_Pt");
		}
		
		TH1D* hEffVtxPt    = FindObj<TH1D>(fcorr, fCorrTag, fParticle + "_Eff_Vtx_Pt");
		TH1D* hEffAccTrkPt = FindObj<TH1D>(fcorr, fCorrTag, fParticle + "_Eff_AccTrk_Pt");
		
		// apply corrections
		
		fdebug->cd(fOutputTag.Data());
		
		TH1D* hPidCorrPt = 0;
		if(fPid)
		{
			if(fDebugLevel>1) std::cout << "\t\tpid correction" << std::endl;
			
			hPidCorrPt = this->PID(hOrigPt, hPidPt, fPtPid, fPtMax, fParticle + "_PidCorr_Pt");
			hPidCorrPt->Scale(1./fPidEff);
		}
		else
		{
			hPidCorrPt = dynamic_cast<TH1D*>(hOrigPt->Clone(Form("%s_PidCorr_Pt", fParticle.Data())));
		}
		
		TH1D* hSecCorrPt = 0;
		if(fSecondaries)
		{
			if(fDebugLevel>1) std::cout << "\t\tremoval of secondaries" << std::endl;
			
			if(fFitFrac)
			{
				hSecCorrPt = this->Secondaries(hPidCorrPt, fncFracMatPt, fncFracFdwnPt, fParticle + "_SecCorr_Pt");
			}
			else
			{
				hSecCorrPt = this->Secondaries(hPidCorrPt, hFracMatPt, hFracFdwnPt, fParticle + "_SecCorr_Pt");
			}
		}
		else
		{
			hSecCorrPt = dynamic_cast<TH1D*>(hPidCorrPt->Clone(Form("%s_SecCorr_Pt", fParticle.Data())));
		}
		
		TH1D* hEffCorrPt = 0;
		if(fEfficiency)
		{
			if(fDebugLevel>1) std::cout << "\t\tefficiency correction" << std::endl;
			
			if(fMCtoINEL)
			{
				TH1D* hStatsMC = FindObj<TH1D>(fcorr, fCorrTag, species + "_Stats");
				fVtxFactor = this->GetVertexCorrection(hStats, hStatsMC);
			}
			
			hEffCorrPt = this->Efficiency(hSecCorrPt, hEffVtxPt, hEffAccTrkPt, fParticle + "_EffCorr_Pt");
		}
		else
		{
			hEffCorrPt = dynamic_cast<TH1D*>(hSecCorrPt->Clone(Form("%s_EffCorr_Pt", fParticle.Data())));
		}
		
		// corrected pt
		
		hCorrPt = dynamic_cast<TH1D*>(hEffCorrPt->Clone(Form("%s_Corr_Pt", fParticle.Data())));
		
		// debug histograms
		
		fdebug->cd(fOutputTag.Data());
		
		hOrigPt->Write(); // original distribution
		hPidCorrPt->Write();
		hSecCorrPt->Write();
		hEffCorrPt->Write();
		hCorrPt->Write();
		
		// free memory
		delete hPidCorrPt;
		delete hSecCorrPt;
		delete hEffCorrPt;
	}
	
	// save final pt
	
	TH1D* hPt = dynamic_cast<TH1D*>(hCorrPt->Clone(Form("%s_Pt", fParticle.Data())));
	hPt->SetTitle(fParticle.Data());
	hPt->Reset();
	hPt->Sumw2();
	
	Int_t lowbin = hPt->GetXaxis()->FindFixBin(fPtMin);
	Int_t hibin  = hPt->GetXaxis()->FindFixBin(fPtMax);
	
	for(Int_t i=lowbin; i<hibin; ++i)
	{
		hPt->SetBinContent(i, hCorrPt->GetBinContent(i));
		hPt->SetBinError(i, hCorrPt->GetBinError(i));
	}
	
	foutput->cd();
	
	hPt->Write();
	
	delete hPt;
	delete hCorrPt;
	
	delete fdebug;
	delete fcorr;
	
	// stats
	
	if(fMakeStats)
	{
		TH1D* hStatsFix = new TH1D(Form("%s_Stats", fParticle.Data()), "Stats", 6, 0, 6);
		hStatsFix->GetXaxis()->SetBinLabel(1,"Events");
		hStatsFix->GetXaxis()->SetBinLabel(2,"Trig");
		hStatsFix->GetXaxis()->SetBinLabel(3,"AnaT");
		hStatsFix->GetXaxis()->SetBinLabel(4,"Inel");
		hStatsFix->GetXaxis()->SetBinLabel(5,"Vtx");
		hStatsFix->GetXaxis()->SetBinLabel(6,"Vz");
		hStatsFix->SetStats(0);
		
		// Fill
		hStatsFix->Fill("Events",hStats->Integral(1,1));
		hStatsFix->Fill("Trig",hStats->Integral(2,2));
		if(fIsOnlyGen)
		{
			hStatsFix->Fill("Ana",hStats->Integral(2,2));
			hStatsFix->Fill("Inel",hStats->Integral(2,2));
		}
		else
		{
			Double_t nTrig = hStats->Integral(2,2);
			Double_t nAna  = hStats->Integral(3,3);
			Double_t nVtx  = hStats->Integral(4,4);
			Double_t nVz   = hStats->Integral(5,5);
			Double_t nAnaT = nTrig; // all triggering events since we have the MC correction
			Double_t nInel = nTrig/fTrigEff;
			
			if(!fMCtoINEL)
			{
				// Add events without vertex to get the proportional trig. events
				Double_t nNoVtx = (nAna/nVtx)*(nTrig-nVtx);
				nAnaT = nAna + nNoVtx;
				nInel = nAnaT/fTrigEff;
			}
			
			hStatsFix->Fill("AnaT",nAnaT);
			hStatsFix->Fill("Inel",nInel);
			hStatsFix->Fill("Vtx",nVtx);
			hStatsFix->Fill("Vz",nVz);
		}
		
		foutput->cd();
		hStatsFix->Write();
		delete hStatsFix;
	}
	
	delete foutput;
	delete finput;
	
	return 0;
}

Double_t AliLnPt::GetVertexCorrection(const TH1D* hData, const TH1D* hMC) const
{
//
// vertex correction factor
//
	Double_t nAna   = hData->Integral(3,3);
	Double_t nVtx   = hData->Integral(4,4);
	Double_t nAnaMC = hMC->Integral(3,3);
	Double_t nVtxMC = hMC->Integral(4,4);
	
	return (nVtx/nAna)/(nVtxMC/nAnaMC);
}

Double_t AliLnPt::GetM2Width(Double_t pt, Double_t m) const
{
//
// expected M2 width parameterization from LHC12a5b
//
	Double_t c0 = 1.90212e-03;
	Double_t c1 = 1.29814e-02;
	
	TString species = fParticle;
	species.ReplaceAll("Anti","");
	
	if(species == "Deuteron")
	{
		c0 = 1.37785e-02;
		c1 = 6.78951e-02;
	}
	else if(species == "Triton")
	{
		c0 = 3.79471e-02;
		c1 = 1.87478e-01;
	}
	else if(species == "He3")
	{
		c0 = 2.05887e-01;
		c1 = 1.17813e-01;
	}
	else if(species == "Alpha")
	{
		c0 = 4.13329e-01;
		c1 = 2.27076e-01;
	}
	
	return c0/(pt*pt) + c1*(pt*pt/(m*m) + 1.);
}

Double_t AliLnPt::GetExpectedDT(Double_t pt, Double_t m) const
{
//
// estimate of the expected time of flight (ns) for mass hypothesis m
//
	Double_t l = 370.; // cm
	Double_t c = 2.99792458e+1; // cm/ns
	
	return (l/c)*TMath::Sqrt(1.+(m/pt)*(m/pt));
}

RooWorkspace* AliLnPt::GetToFModel(Double_t pt, const TString& name) const
{
//
// model for the time of flight distribution
//
	using namespace RooFit;
	
	RooWorkspace* w = new RooWorkspace(name.Data());
	
	// signal
	w->factory(Form("AliLnGaussianExpTail::Sd(x[%f,%f], mu[0,-0.06,0.06], sigma[0.120,0.05,0.3], tau[0.1,0.08,0.4])", fBkgMin, fBkgMax));
	
	// background
	Double_t expTd   = this->GetExpectedDT(pt, GetMass("deuteron"));
	Double_t expDTp  = this->GetExpectedDT(pt, GetMass("proton")) - expTd;
	Double_t expDTk  = this->GetExpectedDT(pt, GetMass("kaon")) - expTd;
	Double_t expDTpi = this->GetExpectedDT(pt, GetMass("pion")) - expTd;
	
	w->factory(Form("AliLnGaussianExpTail::ProtonBkg(x, muP[%f,%f-0.2,%f+0.2], widthP[0.120,0.05,0.3], tauP[0.08,0.01,0.2])",expDTp, expDTp, expDTp));
	w->factory(Form("AliLnGaussianExpTail::KaonBkg(x, muK[%f,%f-0.2,%f+0.2], widthK[0.120,0.05,0.3], tauP[0.08,0.01,0.2])",expDTk, expDTk, expDTk));
	w->factory(Form("AliLnGaussianExpTail::PionBkg(x, muPi[%f,%f-0.2,%f+0.2], widthPi[0.120,0.05,0.3], tauPi[0.08,0.01,0.2])",expDTpi, expDTpi, expDTpi));
	
	w->factory("SUM::Bkg( Npi[0.3,0,1.e+9]*PionBkg, Nk[0.24,0,1.e+9]*KaonBkg, Np[0.47,0,1.e+9]*ProtonBkg, Np[0.47,0,1.e+9]*ProtonBkg)");
	
	w->factory("SUM::model( Nd[0,1.e+9]*Sd,  Nbkg[0,1.e+9]*Bkg )");
	
	w->var("x")->SetTitle("t - <t> (ns)");
	
	return w;
}

RooWorkspace* AliLnPt::GetM2Model(Double_t pt, Double_t m, const TString& name, Double_t max) const
{
//
// mass squared model for deuterons
//
	using namespace RooFit;
	
	RooWorkspace* w = new RooWorkspace(name.Data());
	
	// background
	
	w->factory(Form("Exponential::PionBkg(x[%f,%f],lambdaPi[-0.3,-10.,0.])", fBkgMin, fBkgMax));
	w->factory("Exponential::KaonBkg(x,lambdaK[-0.4,-10.,0.])");
	w->factory("Exponential::ProtonBkg(x,lambdaP[-0.5,-10.,0.])");
	
	w->factory("SUM::Bkg( Npi[0.3,0.,1.]*PionBkg, Nk[0.24,0.,1.]*KaonBkg, Np[0.47,0.,1.]*ProtonBkg )");
	
	// signal
	
	Double_t sigma = GetM2Width(pt, m);
	Double_t mm = (fPidProc == kMassSquared) ? m*m : 0;
	
	if(pt<2.0)
	{
		w->factory(Form("AliLnM2::Sd(x,mu[%f,%f,%f], sigma[%f,%f,%f],tau[0.1,0.05,0.4],p[%f,%f,%f])",mm, mm-0.15, mm+0.1, sigma, 0.1, sigma+0.1, 1.16*pt,pt,1.33*pt));
	}
	else // neglect non-gaussian tails
	{
		w->factory(Form("Gaussian::Sd(x, mu[%f,%f,%f], sigma[%f,%f,%f])",mm, mm-0.15, mm+0.1, sigma, 0.1, sigma+0.1));
	}
	
	// model
	w->factory(Form("SUM::model( Nd[0.,%f]*Sd, Nbkg[0.,%f]*Bkg )", max, max));
	
	w->var("x")->SetTitle("#it{m}^{2} (GeV^{2}/c^{4})");
	
	return w;
}

void AliLnPt::GetPtFromPid(TH1D* hPt, Double_t ptmin, Double_t ptmax, const TH2D* hPidPt, Double_t pidMin, Double_t pidMax) const
{
//
// Integrate pid distribution
//
	Int_t lowbin = hPt->GetXaxis()->FindFixBin(ptmin);
	Int_t hibin  = hPt->GetXaxis()->FindFixBin(ptmax);
	
	for(Int_t i=lowbin; i<hibin; ++i)
	{
		TH1D* hPid = hPidPt->ProjectionY(Form("%s_PID_%02d",fParticle.Data(),i),i,i);
		
		Int_t pidLowBin = hPid->GetXaxis()->FindFixBin(pidMin);
		Int_t pidHiBin = hPid->GetXaxis()->FindFixBin(pidMax);
		
		Double_t err = 0;
		Double_t n = hPid->IntegralAndError(pidLowBin, pidHiBin, err);
		
		hPt->SetBinContent(i, n);
		hPt->SetBinError(i, err);
		
		delete hPid;
	}
}

TH1D* AliLnPt::PID(const TH1D* hPt, const TH2D* hPidPt2, Double_t ptmin, Double_t ptmax, const TString& name)
{
//
// Remove contamination in the pt distribution
//
	TH1D* hPidCorrPt = dynamic_cast<TH1D*>(hPt->Clone(name.Data()));
	if(hPidCorrPt == 0) return 0;
	
	hPidCorrPt->Reset();
	hPidCorrPt->Sumw2();
	
	// set binning
	const Int_t kNBin = 2;
	
	TH2D* hPidPt = dynamic_cast<TH2D*>(hPidPt2->Clone("_pid_pt_"));
	hPidPt->Rebin2D(1,kNBin);
	
	// integrate around the mean value
	
	this->GetPtFromPid(hPidCorrPt, 0, ptmin, hPidPt, fIntMin, fIntMax);
	
	// slice the m2-pt
	
	using namespace RooFit;
	
	Int_t lowbin = hPt->GetXaxis()->FindFixBin(ptmin);
	Int_t hibin  = hPt->GetXaxis()->FindFixBin(ptmax);
	
	Double_t m = GetMass(fParticle);
	
	for(Int_t i=lowbin; i<hibin; ++i)
	{
		TH1D* hPid = hPidPt->ProjectionY(Form("%s_PID_%02d",fParticle.Data(),i),i,i);
		
		RooWorkspace* w = (fPidProc==kTimeOfFlight) ? GetToFModel( hPt->GetBinCenter(i), Form("%s_M2_%02d",fParticle.Data(),i)) : GetM2Model(hPt->GetBinCenter(i), m, Form("%s_M2_%02d",fParticle.Data(),i),hPid->Integral());
		
		RooRealVar* x = w->var("x");
		RooDataHist data("data", "data to fit", *x, hPid);
		
		w->import(data);
		
		w->pdf("model")->fitTo(data, Minos(kTRUE), PrintLevel(-1), Verbose(kFALSE), Warnings(kFALSE), PrintEvalErrors(-1) );
		
		hPidCorrPt->SetBinContent(i, w->var("Nd")->getVal());
		hPidCorrPt->SetBinError(i, w->var("Nd")->getError());
		
		w->Write();
		
		delete w;
		delete hPid;
	}
	
	delete hPidPt;
	return hPidCorrPt;
}

TH1D* AliLnPt::Secondaries(const TH1D* hPt, const TH1D* hFracMatPt, const TH1D* hFracFdwnPt, const TString& name) const
{
//
// remove contamination of secondaries using fractions
// N_prim = N/(1+fmat+fdwn)
//
	TF1* one = new TF1("one","1",0,10);
	
	TH1D* xfactor = (TH1D*)hFracMatPt->Clone(Form("%s_1_fmat_fdwn_",fParticle.Data()));
	if(fFdwnCorr) xfactor->Add(hFracFdwnPt);
	xfactor->Add(one);
	
	TH1D* hSecCorPt = (TH1D*)hPt->Clone(name.Data());
	
	hSecCorPt->Divide(xfactor);
	
	delete xfactor;
	delete one;
	
	return hSecCorPt;
}

TH1D* AliLnPt::Secondaries(const TH1D* hPt, TF1* fncMatPt, TF1* fncFdwnPt, const TString& name) const
{
//
// remove contamination of secondaries using a function fitted to the fractions
// N_prim = N/(1+fmat+fdwn)
//
	TF1* fnc = new TF1("_secondaries_","1 + [0]*exp(-[1]*x)+[2] + [3]*exp(-[4]*x)+[5]", 0, 10);
	
	Double_t param[6] = {0};
	fncMatPt->GetParameters(param);
	if(fFdwnCorr) fncFdwnPt->GetParameters(&param[3]);
	
	Double_t err[6] = {0};
	for(Int_t i=0; i<3; ++i) err[i] = fncMatPt->GetParError(i);
	if(fFdwnCorr) for(Int_t i=3; i<6; ++i) err[i] = fncFdwnPt->GetParError(i);
	
	fnc->SetParameters(param);
	fnc->SetParErrors(err);
	
	TH1D* hSecCorPt = (TH1D*)hPt->Clone(name.Data());
	hSecCorPt->Divide(fnc);
	
	delete fnc;
	
	return hSecCorPt;
}

TH1D* AliLnPt::Efficiency(const TH1D* hPt, const TH1D* hEffVtxPt, const TH1D* hEffAccTrkPt, const TString& name) const
{
//
// correct by efficiency
//
	TH1D* hEffCorPt = (TH1D*)hPt->Clone(name.Data());
	hEffCorPt->Divide(hEffAccTrkPt);
	
	if(fMCtoINEL)
	{
		hEffCorPt->Divide(hEffVtxPt);
		hEffCorPt->Scale(fVtxFactor);
	}
	
	return hEffCorPt;
}
