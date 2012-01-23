/////////////////////////////////////////////////////////
//Macro to read the results of the 3gausfit for spectra
//E. Biolcati, 13-feb-2010
/////////////////////////////////////////////////////////
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TLatex.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TFile.h>
#include <TStyle.h>
#include <TPad.h>
#endif

void Labella(Int_t i);
void Legenda(TH1 *h2, TH1 *h3, TH1 *h4);

//_______________________________________________________
void MakeCorrectedITSsaSpectraMultiBin(Int_t multibin=0){

	TString dirNameSIM, dirNameDATA;
	TFile *fiSIM  = new TFile(Form("../gridmultiplicitybins/LHC10d1_1.5sigma_7DCA_negmag/Spectra_MC_negmag_MultBin%d/SpectraReco.root",multibin));
	TFile *fiDATA = new TFile(Form("../gridmultiplicitybins/data_1.5sigma_7DCA_negmag/Spectra_data_negmag_MultBin%d/SpectraReco.root",multibin));
	TFile *fout   = new TFile(Form("ITSsaSpectraCorr_MultiBin%d.root",multibin),"recreate");

	//dca correction
	TFile *fPanosPosP= new TFile(Form("RootFilesPanosCorr/PanosCorr_1.5sigma_7DCA_PosP_expostrange_MultBin%d.root",multibin));  
	TFile *fPanosNegP= new TFile(Form("RootFilesPanosCorr/PanosCorr_1.5sigma_7DCA_NegP_expostrange_MultBin%d.root",multibin)); 
	TH1F *hPanosPosP= (TH1F*)fPanosPosP->Get("fHistSecTOTCorrDATAMC"); 
	TH1F *hPanosNegP= (TH1F*)fPanosNegP->Get("fHistSecTOTCorrDATAMC"); 
	for(Int_t pbin=0;pbin<=hPanosPosP->GetNbinsX();pbin++)hPanosPosP->SetBinError(pbin,0); 
	for(Int_t pbin=0;pbin<=hPanosNegP->GetNbinsX();pbin++)hPanosNegP->SetBinError(pbin,0);

	//nevts
	TH1F *hstat = (TH1F*)fiDATA->Get("fHistNEvents");
	Double_t NEvts = hstat->GetBinContent(hstat->FindBin(1.));
	cout<<"Event number used for the normalization "<<NEvts<<endl;

	//canvas
	gStyle->SetOptStat(0);
	TCanvas *cs=new TCanvas("cs","cs",1200,700);
	cs->Divide(2,1);
	TCanvas *cmix=new TCanvas("cmix","cmix",800,900);
	cmix->Divide(1,3,0,0);

	TH1F *hPieff[2];
	TH1F *hKeff[2];
	TH1F *hPeff[2];
	TH1F *hPiDATA[2];
	TH1F *hKDATA[2];
	TH1F *hPDATA[2];
	TH1F *hPiDATANorm[2];
	TH1F *hKDATANorm[2];
	TH1F *hPDATANorm[2];
	TH1F *hKsuPi[2]; 
	TH1F *hPsuPi[2]; 
	TH1F *hPsuK[2]; 

	hPieff[0] = (TH1F*)fiSIM->Get("hCorrFacNeg0");
	hKeff[0]  = (TH1F*)fiSIM->Get("hCorrFacNeg1");
	hPeff[0]  = (TH1F*)fiSIM->Get("hCorrFacNeg2");
	hPieff[1] = (TH1F*)fiSIM->Get("hCorrFacPos0");
	hKeff[1]  = (TH1F*)fiSIM->Get("hCorrFacPos1");
	hPeff[1]  = (TH1F*)fiSIM->Get("hCorrFacPos2");
	TCanvas *ceff=new TCanvas("ceff","ceff",800,500);
	ceff->Divide(2,1);
	ceff->cd(1);
	gPad->SetGridy();
	hPieff[0]->Draw();
	hKeff[0]->Draw("same");
	hPeff[0]->Draw("same");
	ceff->cd(2);
	gPad->SetGridy();
	hPieff[1]->Draw();
	hKeff[1]->Draw("same");
	hPeff[1]->Draw("same");

	hPiDATA[0] = (TH1F*)fiDATA->Get("hSpectraNeg0");
	hKDATA[0]  = (TH1F*)fiDATA->Get("hSpectraNeg1");
	hPDATA[0]  = (TH1F*)fiDATA->Get("hSpectraNeg2");
	hPiDATA[1] = (TH1F*)fiDATA->Get("hSpectraPos0");
	hKDATA[1]  = (TH1F*)fiDATA->Get("hSpectraPos1");
	hPDATA[1]  = (TH1F*)fiDATA->Get("hSpectraPos2");

	hPiDATANorm[0] = new TH1F(*hPiDATA[0]);
	hPiDATANorm[0]->SetName("hSpectraPiPlusN");
	hPiDATANorm[1] = new TH1F(*hPiDATA[1]);
	hPiDATANorm[1]->SetName("hSpectraPiMinusN");
	hKDATANorm[0]  = new TH1F(*hKDATA[0]);
	hKDATANorm[0]->SetName("hSpectraKPlusN");
	hKDATANorm[1]  = new TH1F(*hKDATA[1]);
	hKDATANorm[1]->SetName("hSpectraKMinusN");
	hPDATANorm[0]  = new TH1F(*hPDATA[0]);
	hPDATANorm[0]->SetName("hSpectraPPlusN");
	hPDATANorm[1]  = new TH1F(*hPDATA[1]);
	hPDATANorm[1]->SetName("hSpectraPMinusN");

	for(Int_t i=0;i<2;i++){ //0==> negative, 1==>positive
		//line colors
		hPiDATA[i] -> SetLineColor(2);
		hKDATA[i]  -> SetLineColor(3);
		hPDATA[i]  -> SetLineColor(4);
		hPiDATA[i] -> SetMarkerStyle(27);
		hKDATA[i]  -> SetMarkerStyle(23);
		hPDATA[i]  -> SetMarkerStyle(24);
		hPiDATA[i] -> SetMarkerColor(2);
		hKDATA[i]  -> SetMarkerColor(3);
		hPDATA[i]  -> SetMarkerColor(4);
		hPiDATA[i] -> SetMarkerStyle(23);
		hKDATA[i]  -> SetMarkerStyle(23);
		hPDATA[i]  -> SetMarkerStyle(23);

		hPiDATANorm[i] -> SetLineColor(2);
		hKDATANorm[i]  -> SetLineColor(3);
		hPDATANorm[i]  -> SetLineColor(4);
		hPiDATANorm[i] -> SetMarkerStyle(27);
		hKDATANorm[i]  -> SetMarkerStyle(23);
		hPDATANorm[i]  -> SetMarkerStyle(24);
		hPiDATANorm[i] -> SetMarkerColor(2);
		hKDATANorm[i]  -> SetMarkerColor(3);
		hPDATANorm[i]  -> SetMarkerColor(4);
		hPiDATANorm[i] -> SetMarkerStyle(23);
		hKDATANorm[i]  -> SetMarkerStyle(23);
		hPDATANorm[i]  -> SetMarkerStyle(23);

		//division for efficiency
		hPiDATA[i] -> Divide(hPieff[i]);
		hKDATA[i]  -> Divide(hKeff[i]);
		hPDATA[i]  -> Divide(hPeff[i]);
		hPiDATANorm[i] -> Divide(hPieff[i]);
		hKDATANorm[i]  -> Divide(hKeff[i]);
		hPDATANorm[i]  -> Divide(hPeff[i]);

		//normalization number of events
		hPiDATANorm[i] -> Scale(1./NEvts);
		hKDATANorm[i]  -> Scale(1./NEvts);
		hPDATANorm[i]  -> Scale(1./NEvts);

		//correction factor based on fit to DCA distr for P and Pibar 
		hPDATA[0]->Multiply(hPanosPosP); 
		hPDATA[1]->Multiply(hPanosNegP);

		//drawing
		cs->cd(i+1);
		gPad->SetLogy();
		gPad->SetGridy();
		gPad->SetGridx();
		Float_t minim=0.01;
		Float_t maxim=10.;
		hPiDATANorm[i]->Draw("e");
		hPiDATANorm[i]->GetYaxis()->SetRangeUser(minim,maxim);
		hPiDATANorm[i]->GetYaxis()->SetTitle("#frac{1}{N}#frac{d^{2}N}{dp_{t}dy}");
		hPiDATANorm[i]->SetTitle("DATA from ITSsa - corrected");
		hKDATANorm[i]->Draw("esames");
		hPDATANorm[i]->Draw("esames");
		Labella(i);
		Legenda(hPiDATA[i],hKDATA[i],hPDATA[i]);
		cs->Update();
	}

	//fluka correction for pi
	TFile *fGeanFlukaPi= new TFile(Form("RootFilesGeantFlukaCorrection/correctionForCrossSection.211.root"));
	TH1F *hGeantFlukaPiPos=(TH1F*)fGeanFlukaPi->Get("gHistCorrectionForCrossSectionParticles");
	TH1F *hGeantFlukaPiNeg=(TH1F*)fGeanFlukaPi->Get("gHistCorrectionForCrossSectionAntiParticles");
	for(Int_t binPi=0;binPi<=hPiDATA[0]->GetNbinsX();binPi++){
		Float_t FlukaCorrPiPos=hGeantFlukaPiPos->GetBinContent(hGeantFlukaPiPos->FindBin(hPiDATA[0]->GetBinCenter(binPi)));
		Float_t FlukaCorrPiNeg=hGeantFlukaPiNeg->GetBinContent(hGeantFlukaPiNeg->FindBin(hPiDATA[1]->GetBinCenter(binPi)));
		//cout<<"PiPos  "<<FlukaCorrPiPos<<"  "<<hPiDATA[0]->GetBinCenter(binPi)<<"  " <<binPi <<endl;
		//cout<<"PiNeg  "<<FlukaCorrPiNeg<<"  "<<hPiDATA[0]->GetBinCenter(binPi)<<"  " <<binPi <<endl;
		hPiDATA[0]->SetBinContent(binPi,hPiDATA[0]->GetBinContent(binPi)*FlukaCorrPiPos);
		hPiDATA[1]->SetBinContent(binPi,hPiDATA[1]->GetBinContent(binPi)*FlukaCorrPiNeg);
	}
	//fluka correction for pi
	TFile *fGeanFlukaK= new TFile(Form("RootFilesGeantFlukaCorrection/correctionForCrossSection.321.root"));
	TH1F *hGeantFlukaKPos=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionParticles");
	TH1F *hGeantFlukaKNeg=(TH1F*)fGeanFlukaK->Get("gHistCorrectionForCrossSectionAntiParticles");
	for(Int_t binK=0;binK<=hKDATA[0]->GetNbinsX();binK++){
		Float_t FlukaCorrKPos=hGeantFlukaKPos->GetBinContent(hGeantFlukaKPos->FindBin(hKDATA[0]->GetBinCenter(binK)));
		Float_t FlukaCorrKNeg=hGeantFlukaKNeg->GetBinContent(hGeantFlukaKNeg->FindBin(hKDATA[1]->GetBinCenter(binK)));
		//cout<<"KPos :"<<FlukaCorrKPos<<"  "<<hKDATA[0]->GetBinCenter(binK)<<"  " <<binK <<endl;
		//cout<<"KNeg :"<<FlukaCorrKNeg<<"  "<<hKDATA[0]->GetBinCenter(binK)<<"  " <<binK <<endl;
		hKDATA[0]->SetBinContent(binK,hKDATA[0]->GetBinContent(binK)*FlukaCorrKPos);
		hKDATA[1]->SetBinContent(binK,hKDATA[1]->GetBinContent(binK)*FlukaCorrKNeg);
	}
	//fluka correction for P
	//ITS specific file for protons/antiprotons
	Int_t kPos=0;
	Int_t kNeg=1;
	TFile* fITS = new TFile ("RootFilesGeantFlukaCorrection/correctionForCrossSectionITS_20100719.root");
	TH2D * hCorrFlukaITS[2];
	hCorrFlukaITS[kPos] = (TH2D*)fITS->Get("gHistCorrectionForCrossSectionProtons");
	hCorrFlukaITS[kNeg] = (TH2D*)fITS->Get("gHistCorrectionForCrossSectionAntiProtons");
	for(Int_t icharge = 0; icharge < 2; icharge++){
		Int_t nbins = hPDATA[0]->GetNbinsX();
		Int_t nbinsy=hCorrFlukaITS[icharge]->GetNbinsY();
		for(Int_t ibin = 0; ibin < nbins; ibin++){
			Float_t pt = hPDATA[0]->GetBinCenter(ibin);
			Float_t minPtCorrection = hCorrFlukaITS[icharge]->GetYaxis()->GetBinLowEdge(1);
			Float_t maxPtCorrection = hCorrFlukaITS[icharge]->GetYaxis()->GetBinLowEdge(nbinsy+1);
			if (pt < minPtCorrection) pt = minPtCorrection+0.0001;
			if (pt > maxPtCorrection) pt = maxPtCorrection;
			Float_t correction = hCorrFlukaITS[icharge]->GetBinContent(1,hCorrFlukaITS[icharge]->GetYaxis()->FindBin(pt));
			if(icharge==0){
				if (correction != 0) {// If the bin is empty this is a  0
					hPDATA[0]->SetBinContent(ibin,hPDATA[0]->GetBinContent(ibin)*correction);
					hPDATA[0]->SetBinError(ibin,hPDATA[0]->GetBinError  (ibin)*correction);
				}else if (hPDATA[0]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
					cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for protons secondaries, ITS, " << endl;
					cout << " Bin content: " << hPDATA[0]->GetBinContent(ibin) << endl;
				}
			}
			if(icharge==1){
				if (correction != 0) {// If the bin is empty this is a  0
					hPDATA[1]->SetBinContent(ibin,hPDATA[1]->GetBinContent(ibin)*correction);
					hPDATA[1]->SetBinError(ibin,hPDATA[1]->GetBinError  (ibin)*correction);
				}else if (hPDATA[1]->GetBinContent(ibin) > 0) { // If we are skipping a non-empty bin, we notify the user
					cout << "Fluka/GEANT: Not correcting bin "<<ibin << " for Antiprotons secondaries, ITS, " << endl;
					cout << " Bin content: " << hPDATA[1]->GetBinContent(ibin) << endl;
				}
			}
		}
	}


	//mixed particle ratios
	for(Int_t i=0;i<2;i++){
		hKsuPi[i] = (TH1F*)hKDATA[i]->Clone("KsuPi");
		hPsuPi[i] = (TH1F*)hPDATA[i]->Clone("PsuPi");
		hPsuK[i]  = (TH1F*)hPDATA[i]->Clone("PsuK");
		hKsuPi[i] -> Divide(hPiDATA[i]);
		hPsuPi[i] -> Divide(hPiDATA[i]);
		hPsuK[i]  -> Divide(hKDATA[i]);
	}

	//positive/negative ratios
	TH1F *hPiratio = (TH1F*)hPiDATA[1]->Clone("PionsRatio");
	TH1F *hKratio  = (TH1F*)hKDATA[1]->Clone("KaonsRatio");
	TH1F *	hPratio  = (TH1F*)hPDATA[1]->Clone("ProtonsRatio");
	hPiratio -> Divide(hPiDATA[0]);
	hKratio  -> Divide(hKDATA[0]);
	hPratio  -> Divide(hPDATA[0]);


	//drawing positive/negative ratios
	TCanvas *cratio=new TCanvas("cratio","",980,600);
	cratio->Divide(1,3,0,0);
	cratio->SetBottomMargin(0.08);
	cratio->cd(1);
	gPad->SetGridy();
	hPiratio->GetYaxis()->SetTitle("");
	hPiratio->GetYaxis()->SetRangeUser(0.7,1.3);
	hPiratio->SetTitle("");
	hPiratio->Draw("mp");
	TLatex *ll1=new TLatex(0.7,0.7,"#pi^{+}/#pi^{-}");
	ll1->SetNDC();
	ll1->SetTextSize(0.14);
	ll1->Draw();
	cratio->cd(2);
	gPad->SetGridy();
	hKratio->SetTitle("");
	hKratio->GetYaxis()->SetRangeUser(0.7,1.3);
	hKratio->Draw("mp");
	TLatex *ll2=new TLatex(0.7,0.7,"K^{+}/K^{-}");
	ll2->SetNDC();
	ll2->SetTextSize(0.14);
	ll2->Draw();
	cratio->cd(3);
	gPad->SetGridy();
	hPratio->SetTitle("");
	hPratio->GetYaxis()->SetRangeUser(0.7,1.3);
	hPratio->GetXaxis()->SetTitle("p_{t} [GeV/c]");
	hPratio->Draw("mp");
	TLatex *ll3=new TLatex(0.7,0.7,"p/#bar{p}");
	ll3->SetNDC();
	ll3->SetTextSize(0.144);
	ll3->Draw();

	//drawing mixed particle ratios
	gStyle->SetOptTitle(0);
	cmix->cd(1);
	hKsuPi[0]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
	hKsuPi[0]->Draw();
	hKsuPi[1]->Draw("same");
	hKsuPi[0]->SetMinimum(0);
	hKsuPi[0]->SetMarkerStyle(23);
	hKsuPi[1]->SetMarkerStyle(24);
	TLegend *legm1=new TLegend(0.2,0.6,0.39,0.89,NULL,"brNDC");
	legm1->AddEntry(hKsuPi[0],"K^{-}/#pi^{-}","p");
	legm1->AddEntry(hKsuPi[1],"K^{+}/#pi^{+}","p");
	legm1->SetFillColor(0);
	legm1->SetBorderSize(0);
	legm1->Draw();
	hKsuPi[0]->SetMarkerColor(2);
	hKsuPi[1]->SetMarkerColor(4);
	hKsuPi[0]->SetLineColor(2);
	hKsuPi[1]->SetLineColor(4);

	cmix->cd(2);
	hPsuPi[0]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
	hPsuPi[0]->Draw();
	hPsuPi[1]->Draw("same");
	hPsuPi[0]->SetMinimum(0);
	hPsuPi[0]->SetMarkerStyle(23);
	hPsuPi[1]->SetMarkerStyle(24);
	hPsuPi[0]->SetMarkerColor(2);
	hPsuPi[1]->SetMarkerColor(4);
	hPsuPi[0]->SetLineColor(2);
	hPsuPi[1]->SetLineColor(4);
	TLegend *legm2=new TLegend(0.2,0.6,0.39,0.89,NULL,"brNDC");
	legm2->AddEntry(hPsuPi[0],"#bar{p}/#pi^{-}","p");
	legm2->AddEntry(hPsuPi[1],"p/#pi^{+}","p");
	legm2->SetFillColor(0);
	legm2->SetBorderSize(0);
	legm2->Draw();

	cmix->cd(3);
	hPsuK[0]->GetXaxis()->SetTitle("p_{t} [GeV/c]");
	hPsuK[0]->Draw();
	hPsuK[1]->Draw("same");
	hPsuK[0]->SetMinimum(0);
	hPsuK[0]->SetMarkerStyle(23);
	hPsuK[1]->SetMarkerStyle(24);
	hPsuK[0]->SetMarkerColor(2);
	hPsuK[1]->SetMarkerColor(4);
	hPsuK[0]->SetLineColor(2);
	hPsuK[1]->SetLineColor(4);
	TLegend *legm3=new TLegend(0.2,0.6,0.39,0.89,NULL,"brNDC");
	legm3->AddEntry(hPsuPi[0],"#bar{p}/K^{-}","p");
	legm3->AddEntry(hPsuPi[1],"p/K^{+}","p");
	legm3->SetFillColor(0);
	legm3->SetBorderSize(0);
	legm3->Draw();

  //save histograms in the root files
	fout->cd();
	hPiDATA[0]->Write();
	hKDATA[0] ->Write();
	hPDATA[0] ->Write();
	hPiDATA[1]->Write();
	hKDATA[1] ->Write();
	hPDATA[1] ->Write();
	hPiDATANorm[0]->Write();
	hKDATANorm[0] ->Write();
	hPDATANorm[0] ->Write();
	hPiDATANorm[1]->Write();
	hKDATANorm[1] ->Write();
	hPDATANorm[1] ->Write();
	//fout->Close();
	return;

}//end of the main


//_______________________________________________________
void Labella(Int_t i){
	Char_t txt[50];
	if(i==0) sprintf(txt,"negative particles");
	else sprintf(txt,"positive particles");
	TLatex *ltx=new TLatex(0.4,0.3,txt);
	ltx->SetNDC();
	ltx->SetTextColor(6);
	ltx->SetTextFont(22);
	ltx->Draw();
	return;
}

//_______________________________________________________
void Legenda(TH1 *h2, TH1 *h3, TH1 *h4){
	TLegend *leg=new TLegend(0.51,0.11,0.84,0.25,NULL,"brNDC");
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	TLegendEntry *entry2=leg->AddEntry(h2,"pions","p");
	entry2->SetTextColor(2);
	TLegendEntry *entry3=leg->AddEntry(h3,"kaons","p");
	entry3->SetTextColor(3);
	TLegendEntry *entry4=leg->AddEntry(h4,"protons","p");
	entry4->SetTextColor(4);
	leg->Draw("same");
}

//EOF

