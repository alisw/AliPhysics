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

#include "AliLnPt.h"
#include "B2.h"

#ifdef HAVE_ROOUNFOLD
  #include "/opt/RooUnfold/src/RooUnfoldResponse.h"
  #include "/opt/RooUnfold/src/RooUnfoldBayes.h"
#endif

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
, fLowPtBin(1)
, fHiPtBin(19)
, fPidM2(0)
, fUnfolding(0)
, fNIter(4)
, fSecondaries(1)
, fEfficiency(1)
, fIsOnlyGen(0)
, fYMin(-0.5)
, fYMax(0.5)
, fLowM2Bin(9)
, fHiM2Bin(17)
, fMinM2Bkg(2.2)
, fMaxM2Bkg(5.)
, fMinM2tpc(2.)
, fMaxM2tpc(6.)
, fMakeStats(1)
, fMCtoINEL(0)
, fVtxFactor(1)
, fFitFrac(1)
, fFdwnCorr(1)
, fSameFdwn(0)
{
//
// constructor
//
	TH1::SetDefaultSumw2(); // switch on histogram errors
	
	// disable verbose in RooFit
	RooMsgService::instance().setSilentMode(1);
	RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
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
	using namespace std;
	
	TFile* finput = new TFile(fInputFilename.Data(), "read");
	if (finput->IsZombie()) exit(1);
	
	TFile* foutput = new TFile(fOutputFilename.Data(), "recreate");
	if (foutput->IsZombie()) exit(1);
	
	TFile* fcorr = 0;
	TFile* fdebug = 0;
	
	TString species = fParticle;
	species.ReplaceAll("Anti","");
	
	TH1D* hStats = (TH1D*)FindObj(finput, species + "_Stats");
	
	TH1D* hCorPt = 0;
	
	if(fIsOnlyGen)
	{
		hCorPt = (TH1D*)((TH1D*)FindObj(finput, fParticle + "_Gen_PhS_Prim_Pt"))->Clone( Form("%s_Cor_Pt", fParticle.Data()));
	}
	else // correct the reconstructed pt
	{
		fcorr = new TFile(fCorrFilename.Data(),"read");
		if (fcorr->IsZombie()) exit(1);
		
		fdebug =  new TFile(Form("debug-%s",fOutputFilename.Data()),"recreate");
		if (fdebug->IsZombie()) exit(1);
		
		fdebug->mkdir(fOutputTag.Data());
		
		// -------------- get pointers to the required histograms ---------
		
		// pt distributions
		TH1D* hOrigPt = (TH1D*)FindObj(finput, fParticle + "_PID_Pt");
		TH2D* hM2Pt   = (TH2D*)FindObj(finput, fParticle + "_PID_M2_Pt");
		
		// corrections
#ifdef HAVE_ROOUNFOLD
		TH2D* hResponseMtx = (TH2D*)FindObj(fcorr, fCorrTag, fParticle + "_Response_Matrix");
		TH1D* hMeasuredPt  = (TH1D*)FindObj(fcorr, fCorrTag, fParticle + "_Measured_Pt");
		TH1D* hTruePt      = (TH1D*)FindObj(fcorr, fCorrTag, fParticle + "_True_Pt");
#endif
		
		TH1D* hFracMatPt   = (TH1D*)FindObj(fcorr, fCorrTag, fParticle + "_Frac_Mat_Pt");
		TH1D* hFracFdwnPt  = (TH1D*)FindObj(fcorr, fCorrTag, fParticle + "_Frac_Fdwn_Pt");
		TF1* fncFracMatPt  = (TF1*)FindObj(fcorr, fCorrTag, fParticle + "_Frac_Mat_Fit_Pt");
		TF1* fncFracFdwnPt = (TF1*)FindObj(fcorr, fCorrTag, fParticle + "_Frac_Fdwn_Fit_Pt");
		
		if(fSameFdwn)
		{
			hFracFdwnPt = (TH1D*)FindObj(fcorr, fCorrTag, species + "_Frac_Fdwn_Pt");
			fncFracFdwnPt = (TF1*)FindObj(fcorr, fCorrTag, species + "_Frac_Fdwn_Fit_Pt");
		}
		
		TH1D* hEffVtxPt    = (TH1D*)FindObj(fcorr, fCorrTag, fParticle + "_Eff_Vtx_Pt");
		TH1D* hEffAccTrkPt = (TH1D*)FindObj(fcorr, fCorrTag, fParticle + "_Eff_AccTrk_Pt");
		
		// -------------- apply corrections -----------
		
		fdebug->cd(fOutputTag.Data());
		
		TH1D* hPidPt = 0;
		if(fPidM2 && (fHiPtBin > fLowM2Bin))
		{
			Int_t hibin = (fHiPtBin > fHiM2Bin) ? fHiM2Bin : fHiPtBin;
			hPidPt = this->PID(hOrigPt, hM2Pt, fLowM2Bin, hibin, fParticle + "_M2Corr_Pt");
		}
		else
		{
			hPidPt = (TH1D*)hOrigPt->Clone(Form("%s_M2Corr_Pt", fParticle.Data()));
		}
		
		TH1D* hSecCorPt = 0;
		if(fSecondaries)
		{
			if(fFitFrac)
			{
				hSecCorPt = this->Secondaries(hPidPt, fncFracMatPt, fncFracFdwnPt, fParticle + "_SecCor_Pt");
			}
			else
			{
				hSecCorPt = this->Secondaries(hPidPt, hFracMatPt, hFracFdwnPt, fParticle + "_SecCor_Pt");
			}
		}
		else
		{
			hSecCorPt = (TH1D*)hPidPt->Clone(Form("%s_SecCor_Pt", fParticle.Data()));
		}
		
		TH1D* hUnfoldedPt = 0;
#ifdef HAVE_ROOUNFOLD
		if(fUnfolding)
		{
			hUnfoldedPt = this->Unfolding(hSecCorPt, hMeasuredPt, hTruePt, hResponseMtx, fParticle + "_Unfolded_Pt", fNIter);
		}
		else
#endif
		{
			hUnfoldedPt = (TH1D*)hSecCorPt->Clone(Form("%s_Unfolded_Pt", fParticle.Data()));
		}
		
		TH1D* hEffCorPt = 0;
		if(fEfficiency)
		{
			if(fMCtoINEL)
			{
				TH1D* hStatsMC = (TH1D*)FindObj(fcorr, species + "_Stats");
				fVtxFactor = this->GetVertexCorrection(hStats, hStatsMC);
			}
			
			hEffCorPt = this->Efficiency(hUnfoldedPt, hEffVtxPt, hEffAccTrkPt, fParticle + "_EffCor_Pt");
		}
		else
		{
			hEffCorPt = (TH1D*)hUnfoldedPt->Clone(Form("%s_EffCor_Pt", fParticle.Data()));
		}
		
		// --------- corrected pt --------
		
		hCorPt = (TH1D*) hEffCorPt->Clone(Form("%s_Cor_Pt", fParticle.Data()));
		
		// ---------- debug -----
		
		fdebug->cd(fOutputTag.Data());
		
		hOrigPt->Write(); // original distribution
		hPidPt->Write();
		hUnfoldedPt->Write();
		hSecCorPt->Write();
		hEffCorPt->Write();
		hCorPt->Write();
		
		// ------------- clean --------------
		
		delete hPidPt;
		delete hUnfoldedPt;
		delete hSecCorPt;
		delete hEffCorPt;
	}
	
	// ---------- save the final pt ----------
	
	TH1D* hPt = (TH1D*) hCorPt->Clone(Form("%s_Pt", fParticle.Data()));
	hPt->SetTitle(fParticle.Data());
	hPt->Reset();
	hPt->Sumw2();
	
	for(Int_t i=fLowPtBin; i<fHiPtBin; ++i)
	{
		hPt->SetBinContent(i, hCorPt->GetBinContent(i));
		hPt->SetBinError(i, hCorPt->GetBinError(i));
	}
	
	foutput->cd();
	
	hPt->Write();
	
	delete hPt;
	delete hCorPt;
	
	delete fdebug;
	delete fcorr;
	
	// --------- stats -------------------
	
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
// parameterization of the expected M2 width
// as a function of pt from LHC12a5b
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

RooWorkspace* AliLnPt::GetM2Model(Double_t pt, Double_t m, const TString& name, Double_t max) const
{
//
// deuteron m2 model for TPC pid
// gaussian + exponential background of pi, K and p
//
	using namespace RooFit;
	
	RooWorkspace* w = new RooWorkspace(name.Data());
	
	RooRealVar m2("m2", "m^{2} (GeV^{2}/c^{4})", fMinM2Bkg, fMaxM2Bkg);
	w->import(m2);
	
	Double_t sigma = GetM2Width(pt, m);
	
	// signal (m is shifted to the left)
	w->factory(Form("Gaussian::Sd(m2, mean[%f,%f,%f], width[%f,%f,%f])",m*m, m*m-0.3, m*m+0.05, sigma, sigma-0.05, sigma+0.05));
	
	// background
	w->factory("Exponential::PionBkg(m2,lambda1[-0.5,-100.,0.])");
	w->factory("Exponential::KaonBkg(m2,lambda2[-0.5,-100.,0.])");
	w->factory("Exponential::ProtonBkg(m2,lambda3[-0.5,-100.,0.])");
	
	w->factory(Form("SUM::Bkg( Npi[0.3,0,%f]*PionBkg, Nk[0.24,0,%f]*KaonBkg, Np[0.47,0,%f]*ProtonBkg )",max,max,max));
	
	w->factory(Form("SUM::model( Nd[0,%f]*Sd, Nbkg[0,%f]*Bkg )",max,max));
	
	return w;
}

void AliLnPt::GetPtFromM2(TH1D* hPt, Int_t lowbin, Int_t hibin, const TH2D* hM2Pt, Double_t m2min, Double_t m2max) const
{
//
// Fill the pt from lowbin to hibin integrating
// the m2-pt distribution from m2min to m2max
//
	for(Int_t i=lowbin; i<hibin; ++i)
	{
		TH1D* hM2 = hM2Pt->ProjectionY(Form("_M2_%02d",i),i,i);
		
		Int_t m2minbin = hM2->GetXaxis()->FindFixBin(m2min);
		Int_t m2maxbin = hM2->GetXaxis()->FindFixBin(m2max); // m2 has a non-gaussian tail on the right
		
		Double_t err = 0;
		Double_t n = hM2->IntegralAndError(m2minbin, m2maxbin, err);
		
		hPt->SetBinContent(i, n);
		hPt->SetBinError(i, err);
		
		delete hM2;
	}
}

TH1D* AliLnPt::PID(const TH1D* hPt, const TH2D* hM2Pt, Int_t lowbin, Int_t hibin, const TString& name)
{
//
// Remove contamination in the m^2 distribution
//
	const Int_t kNBin = 2;
	Double_t m = GetMass(fParticle);
	
	TH1D* hPidPt = (TH1D*) hPt->Clone(name.Data());
	hPidPt->Reset();
	hPidPt->Sumw2();
	
	// first bins are not contaminated since the identification is clear
	// integrate around the m2 value
	
	this->GetPtFromM2(hPidPt, 0, lowbin, hM2Pt, fMinM2tpc, fMaxM2tpc);
	
	// the remaining bins are contaminated
	// slice the m2-pt and fit each slice with a gaussian + exponentials bkg
	
	using namespace RooFit;
	
	for(Int_t i=lowbin; i<hibin; ++i)
	{
		TH1D* hM2 = hM2Pt->ProjectionY(Form("%s_PID_M2_%02d",fParticle.Data(),i),i,i);
		hM2->Rebin(kNBin);
		
		RooWorkspace* w = GetM2Model(hPt->GetBinCenter(i), m, Form("%s_M2_%02d",fParticle.Data(),i),hM2->Integral());
		
		RooRealVar* m2 = w->var("m2");
		RooDataHist data("data", "data to fit", *m2, hM2);
		
		w->import(data);
		
		w->pdf("model")->fitTo(data, Minos(kTRUE), PrintLevel(-1), Verbose(kFALSE), Warnings(kFALSE), PrintEvalErrors(-1) );
		
		hPidPt->SetBinContent(i, w->var("Nd")->getVal());
		hPidPt->SetBinError(i, w->var("Nd")->getError());
		
		w->Write();
		
		delete w;
		delete hM2;
	}
	
	return hPidPt;
}

#ifdef HAVE_ROOUNFOLD

TH1D* AliLnPt::Unfolding(const TH1D* hPt, const TH1D* hMeasuredPt, const TH1D* hTruePt, const TH2D* hResponseMtx, const TString& name, Int_t iterations) const
{
//
// Bayesian unfolding
//
	RooUnfoldResponse* matrix = new RooUnfoldResponse(hMeasuredPt, hTruePt, hResponseMtx);
	
	RooUnfoldBayes* unfold = new RooUnfoldBayes(matrix, hPt, iterations);
	
	TH1D* hUnfoldedPt = (TH1D*)unfold->Hreco();
	hUnfoldedPt->SetName(name.Data());
	hUnfoldedPt->SetXTitle("p_{T} (GeV/c)");
	
	delete unfold;
	delete matrix;
	
	return hUnfoldedPt;
}

#endif

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
// remove contamination of secondaries using fitted fractions
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
