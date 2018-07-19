#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TFileMerger.h"
#include "TMap.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TGridResult.h"
#include "TObjString.h"
#endif

void readMCPerform(TString filename="QAresults_AOD.root", Int_t drawOnlyDzerDplus = 1, Int_t runNumber=-1)
{

  const Int_t totTrending=5;
  Float_t vecForTrend[totTrending];
  TString varForTrending[totTrending]={"nDzeroCandperEv","nDplusCandperEv","nDsCandperEv","nLcCandperEv","nDstarCandperEv"};

  

  TTree* trtree=new TTree("trendingHF","tree of trending variables");
  trtree->Branch("nrun",&runNumber,"nrun/I");
  for(Int_t j=0; j<totTrending; j++){
    trtree->Branch(varForTrending[j].Data(),&vecForTrend[j],Form("%s/F",varForTrending[j].Data()));
    vecForTrend[j]=-99.;
  }

  TFile *ff = new TFile(filename.Data());

  Int_t color[5] = {kBlack, kRed, kGreen, kBlue, kOrange};

  TDirectoryFile *dirD2H = (TDirectoryFile *)ff->Get("PWG3_D2H_QA");
  if(!dirD2H){
    printf("Directory PWG3_D2H_QA not found in file %s\n",filename.Data());
    return;
  }
  TList *listD2H = (TList *)dirD2H->Get("nEntriesQA");
  if(!listD2H) listD2H = (TList *)dirD2H->Get("nEntriesQAD0");
  if(!listD2H){
    printf("TList nEntriesQA and nEntriesQAD0 not found in file %s\n",filename.Data());
    return;  
  }
  TH1F *hNentries = (TH1F *)listD2H->FindObject("hNentries");
  TH2F *hHasSelBit = (TH2F *)listD2H->FindObject("HasSelBit");

  TCanvas *cqa = new TCanvas("cqa", "cqa", 800, 500);
  cqa->Divide(2, 1);
  cqa->cd(1);
  hNentries->Draw();
  cqa->cd(2);
  hHasSelBit->Draw("colz");
  cqa->SaveAs("plot_D2HQA.png");

  Double_t nEv=hNentries->GetBinContent(10);
  vecForTrend[0]=hHasSelBit->GetBinContent(1)/nEv;
  vecForTrend[1]=hHasSelBit->GetBinContent(2)/nEv;
  vecForTrend[2]=hHasSelBit->GetBinContent(3)/nEv;
  vecForTrend[3]=hHasSelBit->GetBinContent(4)/nEv;
  vecForTrend[4]=hHasSelBit->GetBinContent(5)/nEv;


  TDirectoryFile *dir = (TDirectoryFile *)ff->Get("PWGHF_D2H_MCPerform");
  TList* list = 0x0;
  if (dir)
  {
    list = (TList *)dir->Get("coutputDperfQA");
    if(list){
    
      TH1F *hn = (TH1F *)list->FindObject("fHistNEvents");
      TH1F *hnGenD = (TH1F *)list->FindObject("fHistNGenD");
      Int_t entries = hn->GetBinContent(3);
      
      TH2F *fHistNCand = (TH2F *)list->FindObject("fHistNCand");
      TH1F *fHistNCandDzero = (TH1F *)fHistNCand->ProjectionY("fHistNCandDzero", 1, 1);
      TH1F *fHistNCandDplus = (TH1F *)fHistNCand->ProjectionY("fHistNCandDplus", 2, 2);
      TH1F *fHistNCandDstar = (TH1F *)fHistNCand->ProjectionY("fHistNCandDstar", 3, 3);
      TH1F *fHistNCandDs = (TH1F *)fHistNCand->ProjectionY("fHistNCandDs", 4, 4);
      TH1F *fHistNCandLc = (TH1F *)fHistNCand->ProjectionY("fHistNCandLc", 5, 5);
      
      TString names[5] = {"Dzero", "Dplus", "Dstar", "Ds", "Lc2pkpi"};
      TString type[2] = {"Prompt", "Feeddown"};
      const Int_t nDecays = 5;
      TH2F *fHistXvtxResVsPt[2 * nDecays];
      TH2F *fHistYvtxResVsPt[2 * nDecays];
      TH2F *fHistZvtxResVsPt[2 * nDecays];
      TH2F *fHistInvMassVsPt[2 * nDecays];
      TH2F *fHistDecLenVsPt[2 * nDecays];
      TH2F *fHistNormDLxyVsPt[2 * nDecays];
      TH2F *fHistCosPointVsPt[2 * nDecays];
      
      TH3F *fHistPtYMultGenDauInAcc[2 * nDecays];
      TH3F *fHistPtYMultRecoFilt[2 * nDecays];
      
      TProfile *fHistXvtxRes[2 * nDecays];
      TProfile *fHistYvtxRes[2 * nDecays];
      TProfile *fHistZvtxRes[2 * nDecays];
      TProfile *fHistXvtxMean[2 * nDecays];
      TProfile *fHistYvtxMean[2 * nDecays];
      TProfile *fHistZvtxMean[2 * nDecays];
      
      TH1F *fHistXvtxRes2[2 * nDecays];
      TH1F *fHistYvtxRes2[2 * nDecays];
      TH1F *fHistZvtxRes2[2 * nDecays];
      
      TProfile *fHistInvMass[2 * nDecays];
      TProfile *fHistDecLen[2 * nDecays];
      TProfile *fHistCosp[2 * nDecays];
      
      TH1F *fHistInvMassRes[2 * nDecays];
      
      TH1F *hEffPt[2 * nDecays];
      TH1F *htemp;
      
      for (Int_t j = 0; j < 5; j++)
	{ //decays
	  for (Int_t i = 0; i < 2; i++)
	    { //prompt and fd
	      Int_t index = j * 2 + i;
	      fHistXvtxResVsPt[index] = (TH2F *)list->FindObject(Form("hXvtxResVsPt%s%s", type[i].Data(), names[j].Data()));
	      fHistYvtxResVsPt[index] = (TH2F *)list->FindObject(Form("hYvtxResVsPt%s%s", type[i].Data(), names[j].Data()));
	      fHistZvtxResVsPt[index] = (TH2F *)list->FindObject(Form("hZvtxResVsPt%s%s", type[i].Data(), names[j].Data()));
	      fHistInvMassVsPt[index] = (TH2F *)list->FindObject(Form("hInvMassVsPt%s%s", type[i].Data(), names[j].Data()));
	      fHistDecLenVsPt[index] = (TH2F *)list->FindObject(Form("hDecLenVsPt%s%s", type[i].Data(), names[j].Data()));
	      fHistCosPointVsPt[index] = (TH2F *)list->FindObject(Form("hCosPointVsPt%s%s", type[i].Data(), names[j].Data()));
	      fHistPtYMultGenDauInAcc[index] = (TH3F *)list->FindObject(Form("hPtYMult%sGenDauInAcc%s", type[i].Data(), names[j].Data()));
	      fHistPtYMultRecoFilt[index] = (TH3F *)list->FindObject(Form("hPtYMult%sRecoFilt%s", type[i].Data(), names[j].Data()));
	      
	      fHistXvtxMean[index] = (TProfile *)fHistXvtxResVsPt[index]->ProfileX(Form("hXvtxMean%s%s", type[i].Data(), names[j].Data()));
	      fHistXvtxMean[index]->SetLineColor(color[j]);
	      fHistXvtxMean[index]->SetLineWidth(2);
	      fHistXvtxMean[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistXvtxMean[index]->GetYaxis()->SetTitle("Xvtx (reco-true) mean (#mum)");
	      fHistXvtxMean[index]->SetTitle("Xvtx residual vs pT");

	      fHistYvtxMean[index] = (TProfile *)fHistYvtxResVsPt[index]->ProfileX(Form("hYvtxMean%s%s", type[i].Data(), names[j].Data()));
	      fHistYvtxMean[index]->SetLineColor(color[j]);
	      fHistYvtxMean[index]->SetLineWidth(2);
	      fHistYvtxMean[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistYvtxMean[index]->GetYaxis()->SetTitle("Yvtx (reco-true) mean (#mum)");
	      fHistYvtxMean[index]->SetTitle("Yvtx residual vs pT");

	      fHistZvtxMean[index] = (TProfile *)fHistZvtxResVsPt[index]->ProfileX(Form("hZvtxMean%s%s", type[i].Data(), names[j].Data()));
	      fHistZvtxMean[index]->SetLineColor(color[j]);
	      fHistZvtxMean[index]->SetLineWidth(2);
	      fHistZvtxMean[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistZvtxMean[index]->GetYaxis()->SetTitle("Zvtx (reco-true) mean (#mum)");
	      fHistZvtxMean[index]->SetTitle("Zvtx residual vs pT");

	      fHistXvtxRes[index] = (TProfile *)fHistXvtxResVsPt[index]->ProfileX(Form("hXvtxRes%s%s", type[i].Data(), names[j].Data()), 1, -1, "s");
	      fHistYvtxRes[index] = (TProfile *)fHistYvtxResVsPt[index]->ProfileX(Form("hYvtxRes%s%s", type[i].Data(), names[j].Data()), 1, -1, "s");
	      fHistZvtxRes[index] = (TProfile *)fHistZvtxResVsPt[index]->ProfileX(Form("hZvtxRes%s%s", type[i].Data(), names[j].Data()), 1, -1, "s");
	      fHistXvtxRes2[index] = (TH1F *)fHistXvtxResVsPt[index]->ProjectionX(Form("hXvtxRes2%s%s", type[i].Data(), names[j].Data()));
	      fHistYvtxRes2[index] = (TH1F *)fHistYvtxResVsPt[index]->ProjectionX(Form("hYvtxRes2%s%s", type[i].Data(), names[j].Data()));
	      fHistZvtxRes2[index] = (TH1F *)fHistZvtxResVsPt[index]->ProjectionX(Form("hZvtxRes2%s%s", type[i].Data(), names[j].Data()));
	      fHistXvtxRes[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistXvtxRes[index]->GetYaxis()->SetTitle("Xvtx (reco-true) RMS (#mum)");
	      fHistXvtxRes[index]->SetTitle("Xvtx resolution vs pT");
	      fHistYvtxRes[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistYvtxRes[index]->GetYaxis()->SetTitle("Yvtx (reco-true) RMS (#mum)");
	      fHistYvtxRes[index]->SetTitle("Yvtx resolution vs pT");
	      fHistZvtxRes[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistZvtxRes[index]->GetYaxis()->SetTitle("Zvtx (reco-true) RMS (#mum)");
	      fHistZvtxRes[index]->SetTitle("Zvtx resolution vs pT");

	      fHistXvtxRes2[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistXvtxRes2[index]->GetYaxis()->SetTitle("Xvtx (reco-true) RMS (#mum)");
	      fHistXvtxRes2[index]->SetTitle("Xvtx resolution vs pT");
	      fHistYvtxRes2[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistYvtxRes2[index]->GetYaxis()->SetTitle("Yvtx (reco-true) RMS (#mum)");
	      fHistYvtxRes2[index]->SetTitle("Yvtx resolution vs pT");
	      fHistZvtxRes2[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistZvtxRes2[index]->GetYaxis()->SetTitle("Zvtx (reco-true) RMS (#mum)");
	      fHistZvtxRes2[index]->SetTitle("Zvtx resolution vs pT");

	      fHistXvtxRes2[index]->SetLineColor(color[j]);
	      fHistYvtxRes2[index]->SetLineColor(color[j]);
	      fHistZvtxRes2[index]->SetLineColor(color[j]);
	      fHistXvtxRes2[index]->SetMarkerColor(color[j]);
	      fHistXvtxRes2[index]->SetMarkerStyle(20);
	      fHistYvtxRes2[index]->SetMarkerColor(color[j]);
	      fHistYvtxRes2[index]->SetMarkerStyle(20);
	      fHistZvtxRes2[index]->SetMarkerColor(color[j]);
	      fHistZvtxRes2[index]->SetMarkerStyle(20);
	      fHistXvtxRes2[index]->SetLineWidth(2);
	      fHistYvtxRes2[index]->SetLineWidth(2);
	      fHistZvtxRes2[index]->SetLineWidth(2);
	      fHistXvtxRes2[index]->Sumw2();
	      fHistYvtxRes2[index]->Sumw2();
	      fHistZvtxRes2[index]->Sumw2();

	      fHistInvMass[index] = (TProfile *)fHistInvMassVsPt[index]->ProfileX(Form("hInvMassVsPt%s%s", type[i].Data(), names[j].Data()));
	      fHistInvMass[index]->SetLineColor(color[j]);
	      fHistInvMass[index]->SetLineWidth(2);
	      fHistInvMass[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistInvMass[index]->GetYaxis()->SetTitle("Inv Mass (GeV/c2)");
	      fHistInvMass[index]->SetTitle("Inv Mass vs pT");

	      fHistDecLen[index] = (TProfile *)fHistDecLenVsPt[index]->ProfileX(Form("hDecLenVsPt%s%s", type[i].Data(), names[j].Data()));
	      fHistDecLen[index]->SetLineColor(color[j]);
	      fHistDecLen[index]->SetLineWidth(2);
	      fHistDecLen[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistDecLen[index]->GetYaxis()->SetTitle("Dec Len (#mum)");
	      fHistDecLen[index]->SetTitle("Prompt Dec Len vs pT");

	      fHistCosp[index] = (TProfile *)fHistCosPointVsPt[index]->ProfileX(Form("hCosPVsPt%s%s", type[i].Data(), names[j].Data()));
	      fHistCosp[index]->SetLineColor(color[j]);
	      fHistCosp[index]->SetLineWidth(2);
	      fHistCosp[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      fHistCosp[index]->GetYaxis()->SetTitle("Cos Point");
	      fHistCosp[index]->SetTitle("Prompt CosPoint vs pT");

	      if (index % 2 == 1)
		fHistDecLen[index]->SetTitle("FeedDown Dec Len vs pT");

	      htemp = (TH1F *)fHistPtYMultGenDauInAcc[index]->ProjectionX(Form("hPtDen%s%s", type[i].Data(), names[j].Data()));
	      hEffPt[index] = (TH1F *)fHistPtYMultRecoFilt[index]->ProjectionX(Form("hPtNum%s%s", type[i].Data(), names[j].Data()));
	      fHistPtYMultGenDauInAcc[index]->Sumw2();
	      fHistPtYMultRecoFilt[index]->Sumw2();
	      hEffPt[index]->Sumw2();
	      hEffPt[index]->Divide(htemp);
	      hEffPt[index]->SetLineColor(color[j]);
	      hEffPt[index]->SetLineWidth(2);
	      hEffPt[index]->GetXaxis()->SetTitle("pT (GeV/c)");
	      hEffPt[index]->GetYaxis()->SetTitle("Prompt Efficiency");
	      hEffPt[index]->SetTitle("Prompt Efficiency");
	      hEffPt[index]->SetStats(0);
	      fHistCosp[index]->SetStats(0);
	      fHistDecLen[index]->SetStats(0);

	      if (index % 2 == 1)
		{
		  hEffPt[index]->GetYaxis()->SetTitle("Feeddown Efficiency");
		  hEffPt[index]->SetTitle("Feeddown Efficiency");
		}

	      fHistInvMassRes[index] = new TH1F(*hEffPt[index]);
	      for (Int_t jj = 1; jj < hEffPt[index]->GetNbinsX() + 1; jj++)
		{
		  TH1F *hTemp = (TH1F *)fHistInvMassVsPt[index]->ProjectionY("htemp", jj, jj);
		  fHistInvMassRes[index]->SetBinContent(jj, hTemp->GetRMS());
		  fHistInvMassRes[index]->SetBinError(jj, hTemp->GetRMSError());
		  fHistInvMassRes[index]->SetLineColor(color[j]);
		  fHistInvMassRes[index]->SetLineWidth(2);
		  fHistInvMassRes[index]->GetXaxis()->SetTitle("pT (GeV/c)");
		  fHistInvMassRes[index]->GetYaxis()->SetTitle("Inv Mass RMS (GeV/c2)");
		  fHistInvMassRes[index]->SetTitle("Inv Mass RMS vs pT");
		  if (index == 0)
		    printf("D0: pt=%f, res=%f \n", fHistInvMassRes[index]->GetBinCenter(jj), fHistInvMassRes[index]->GetBinContent(jj));
		  TH1F *hTempX = (TH1F *)fHistXvtxResVsPt[index]->ProjectionY("htempX", jj, jj);
		  TH1F *hTempY = (TH1F *)fHistYvtxResVsPt[index]->ProjectionY("htempY", jj, jj);
		  TH1F *hTempZ = (TH1F *)fHistZvtxResVsPt[index]->ProjectionY("htempZ", jj, jj);

		  fHistXvtxRes2[index]->SetBinContent(jj, hTempX->GetRMS());
		  fHistXvtxRes2[index]->SetBinError(jj, hTempX->GetRMSError());

		  fHistYvtxRes2[index]->SetBinContent(jj, hTempY->GetRMS());
		  fHistYvtxRes2[index]->SetBinError(jj, hTempY->GetRMSError());

		  fHistZvtxRes2[index]->SetBinContent(jj, hTempZ->GetRMS());
		  fHistZvtxRes2[index]->SetBinError(jj, hTempZ->GetRMSError());
		}
	    }
	}

      fHistNCandDplus->SetLineColor(2);
      fHistNCandDstar->SetLineColor(3);
      fHistNCandDs->SetLineColor(4);
      fHistNCandLc->SetLineColor(kOrange);
      fHistNCandDplus->SetLineWidth(2);
      fHistNCandDstar->SetLineWidth(2);
      fHistNCandDs->SetLineWidth(2);
      fHistNCandLc->SetLineWidth(2);

      fHistNCandDzero->GetXaxis()->SetTitle("pT (GeV/c)");
      fHistNCandDzero->GetYaxis()->SetTitle("counts");
      TLegend *leg = new TLegend(0.6, 0.7, 0.8, 0.9);
      leg->AddEntry(fHistNCandDzero, "Dzero", "l");
      leg->AddEntry(fHistNCandDplus, "Dplus", "l");
      leg->AddEntry(fHistNCandDstar, "Dstar", "l");
      leg->AddEntry(fHistNCandDs, "Ds", "l");
      leg->AddEntry(fHistNCandLc, "Lc", "l");

      TLegend *leg1 = new TLegend(0.5, 0.7, 0.7, 0.9);
      leg1->AddEntry(fHistYvtxRes2[0], "Dzero", "pl");
      leg1->AddEntry(fHistYvtxRes2[2], "Dplus", "pl");
      if (drawOnlyDzerDplus == 0)
	leg1->AddEntry(fHistYvtxRes2[6], "Ds", "pl");
      if (drawOnlyDzerDplus == 0)
	leg1->AddEntry(fHistYvtxRes2[8], "Lc", "pl");

      TLegend *leg2 = new TLegend(0.5, 0.7, 0.7, 0.9);
      leg2->AddEntry(fHistYvtxMean[0], "Dzero", "l");
      leg2->AddEntry(fHistYvtxMean[2], "Dplus", "l");
      if (drawOnlyDzerDplus == 0)
	leg2->AddEntry(fHistYvtxMean[6], "Ds", "l");
      if (drawOnlyDzerDplus == 0)
	leg2->AddEntry(fHistYvtxMean[8], "Lc", "l");

      TLegend *leg3 = new TLegend(0.2, 0.7, 0.4, 0.9);
      leg3->AddEntry(fHistNCandDzero, "Dzero", "l");
      leg3->AddEntry(fHistNCandDplus, "Dplus", "l");
      if (drawOnlyDzerDplus == 0)
	leg3->AddEntry(fHistNCandDstar, "Dstar", "l");
      if (drawOnlyDzerDplus == 0)
	leg3->AddEntry(fHistNCandDs, "Ds", "l");
      if (drawOnlyDzerDplus == 0)
	leg3->AddEntry(fHistNCandLc, "Lc", "l");

      TLegend *leg4 = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg4->AddEntry(fHistNCandDzero, "Dzero", "l");
      leg4->AddEntry(fHistNCandDplus, "Dplus", "l");
      if (drawOnlyDzerDplus == 0)
	leg4->AddEntry(fHistNCandDstar, "Dstar", "l");
      if (drawOnlyDzerDplus == 0)
	leg4->AddEntry(fHistNCandDs, "Ds", "l");
      if (drawOnlyDzerDplus == 0)
	leg4->AddEntry(fHistNCandLc, "Lc", "l");

      TCanvas *c0_1 = new TCanvas("c0_1", "c0_1", 500, 500);
      hnGenD->SetTitle("number of generated D mesons");
      hnGenD->Draw();
      c0_1->SaveAs("plotDgen.png");

      TCanvas *c0_2 = new TCanvas("c0_2", "c0_2", 500, 500);
      c0_2->SetLogy();

      fHistNCandDs->SetTitle("Candidates passing filtering cuts");
      fHistNCandDs->Draw("");
      c0_2->Update();
      TPaveStats *stats = (TPaveStats *)c0_2->GetPrimitive("stats");
      stats->SetName("h1stats");
      stats->SetY1NDC(0.5);
      stats->SetY2NDC(0.35);
      c0_2->Update();

      fHistNCandDplus->Draw("sames");
      c0_2->Update();
      TPaveStats *stats2 = (TPaveStats *)c0_2->GetPrimitive("stats");
      stats2->SetName("h2stats");
      stats2->SetY1NDC(0.8);
      stats2->SetY2NDC(.65);
      c0_2->Update();

      fHistNCandDstar->Draw("sames");
      c0_2->Update();
      TPaveStats *stats3 = (TPaveStats *)c0_2->GetPrimitive("stats");
      stats3->SetName("h3stats");
      stats3->SetY1NDC(0.65);
      stats3->SetY2NDC(.5);
      c0_2->Update();

      fHistNCandDzero->Draw("sames");
      c0_2->Update();
      TPaveStats *stats4 = (TPaveStats *)c0_2->GetPrimitive("stats");
      stats4->SetName("h4stats");
      stats4->SetY1NDC(0.95);
      stats4->SetY2NDC(.8);
      c0_2->Update();

      fHistNCandLc->Draw("sames");
      c0_2->Update();
      TPaveStats *stats5 = (TPaveStats *)c0_2->GetPrimitive("stats");
      stats5->SetName("h1stats");
      stats5->SetY1NDC(0.35);
      stats5->SetY2NDC(.2);
      c0_2->Update();
      leg->Draw();
      c0_2->SaveAs("plotDcandpt.png");

      TCanvas *c0_3 = new TCanvas("c0_3", "c0_3", 500, 500);
      fHistInvMass[0]->SetMinimum(1.6);
      fHistInvMass[0]->SetMaximum(2.4);
      fHistInvMass[0]->Draw();
      fHistInvMass[2]->Draw("sames");
      fHistInvMass[4]->Draw("sames");
      fHistInvMass[6]->Draw("sames");
      fHistInvMass[8]->Draw("sames");
      leg->Draw();
      c0_3->SaveAs("plotDcandInvMass.png");

      TCanvas *c0_4 = new TCanvas("c0_4", "c0_4", 500, 500);
      //fHistInvMassRes[0]->SetMinimum(1.6);
      //fHistInvMassRes[0]->SetMaximum(2.4);

      fHistInvMassRes[0]->GetYaxis()->SetTitleOffset(1.4);
      fHistInvMassRes[0]->SetTitle("D0 Inv Mass RMS vs pT");
      fHistInvMassRes[0]->Draw("");
      //  fHistInvMassRes[2]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistInvMassRes[4]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistInvMassRes[6]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistInvMassRes[8]->Draw("sames");
      //leg->Draw();
      c0_4->SaveAs("plotD0candInvMassWidth.png");

      fHistXvtxMean[0]->SetStats(0);
      fHistYvtxMean[0]->SetStats(0);
      fHistZvtxMean[0]->SetStats(0);
      fHistXvtxMean[2]->SetStats(0);
      fHistYvtxMean[2]->SetStats(0);
      fHistZvtxMean[2]->SetStats(0);
      fHistXvtxRes2[0]->SetStats(0);
      fHistYvtxRes2[0]->SetStats(0);
      fHistZvtxRes2[0]->SetStats(0);
      fHistXvtxRes2[2]->SetStats(0);
      fHistYvtxRes2[2]->SetStats(0);
      fHistZvtxRes2[2]->SetStats(0);

      TCanvas *cc = new TCanvas("cc", "cc", 1200, 500);
      cc->Divide(3, 1);
      cc->cd(1);
      fHistXvtxMean[0]->GetYaxis()->SetTitleOffset(1.4);
      fHistXvtxMean[0]->SetMinimum(-300.);
      fHistXvtxMean[0]->SetMaximum(300.);
      fHistXvtxMean[0]->Draw();

      fHistXvtxMean[2]->Draw("sames");
      leg2->Draw();

      cc->cd(2);

      fHistYvtxMean[0]->GetYaxis()->SetTitleOffset(1.4);
      fHistYvtxMean[0]->SetMinimum(-300.);
      fHistYvtxMean[0]->SetMaximum(300.);
      fHistYvtxMean[0]->Draw();
      fHistYvtxMean[2]->Draw("sames");
      leg2->Draw();

      cc->cd(3);
      fHistZvtxMean[0]->GetYaxis()->SetTitleOffset(1.4);
      fHistZvtxMean[0]->SetMinimum(-300.);
      fHistZvtxMean[0]->SetMaximum(300.);
      fHistZvtxMean[0]->Draw();
      fHistZvtxMean[2]->Draw("sames");
      leg2->Draw();
      cc->SaveAs("plotXYZVtxMean.png");
      /////////

      TCanvas *ccr = new TCanvas("ccr", "ccr", 1200, 500);
      ccr->Divide(3, 1);
      ccr->cd(1);
      fHistXvtxRes2[0]->GetYaxis()->SetTitleOffset(1.4);
      fHistXvtxRes2[0]->SetMinimum(0.);
      fHistXvtxRes2[0]->SetMaximum(500.);
      fHistXvtxRes2[0]->Draw();
      fHistXvtxRes2[2]->Draw("sames");
      leg2->Draw();

      ccr->cd(2);
      fHistYvtxRes2[0]->GetYaxis()->SetTitleOffset(1.4);
      fHistYvtxRes2[0]->SetMinimum(0.);
      fHistYvtxRes2[0]->SetMaximum(500.);
      fHistYvtxRes2[0]->Draw();
      fHistYvtxRes2[2]->Draw("sames");
      leg2->Draw();

      ccr->cd(3);
      fHistZvtxRes2[0]->GetYaxis()->SetTitleOffset(1.4);
      fHistZvtxRes2[0]->SetMinimum(0.);
      fHistZvtxRes2[0]->SetMaximum(500.);
      fHistZvtxRes2[0]->Draw();
      fHistZvtxRes2[2]->Draw("sames");
      leg2->Draw();
      ccr->SaveAs("plotXYZVtxRMS.png");

      TCanvas *ccc = new TCanvas("ccc", "ccc", 1200, 800);
      ccc->Divide(3, 2);
      ccc->cd(1);
      fHistDecLen[0]->GetYaxis()->SetTitleOffset(1.45);
      fHistDecLen[0]->Draw();
      fHistDecLen[2]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistDecLen[4]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistDecLen[6]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistDecLen[8]->Draw("sames");
      leg3->Draw();

      ccc->cd(2);
      fHistCosp[0]->GetYaxis()->SetTitleOffset(1.45);
      fHistCosp[0]->Draw();
      fHistCosp[2]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistCosp[4]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistCosp[6]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistCosp[8]->Draw("sames");
      leg4->Draw();

      ccc->cd(3);
      hEffPt[0]->GetYaxis()->SetTitleOffset(1.45);
      hEffPt[0]->Draw();
      hEffPt[2]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	hEffPt[4]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	hEffPt[6]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	hEffPt[8]->Draw("sames");
      leg3->Draw();

      ccc->cd(4);
      fHistDecLen[1]->GetYaxis()->SetTitleOffset(1.45);
      fHistDecLen[1]->Draw();
      fHistDecLen[3]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistDecLen[5]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistDecLen[7]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistDecLen[9]->Draw("sames");
      leg3->Draw();

      ccc->cd(5);
      fHistCosp[1]->GetYaxis()->SetTitleOffset(1.45);
      fHistCosp[1]->Draw();
      fHistCosp[3]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistCosp[5]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistCosp[7]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	fHistCosp[9]->Draw("sames");
      leg4->Draw();

      ccc->cd(6);
      hEffPt[1]->GetYaxis()->SetTitleOffset(1.45);
      hEffPt[1]->Draw();
      hEffPt[3]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	hEffPt[5]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	hEffPt[7]->Draw("sames");
      if (drawOnlyDzerDplus == 0)
	hEffPt[9]->Draw("sames");
      leg3->Draw();

      ccc->SaveAs("plot_DL_cosp_Eff_prompt_fd.png");
    }
  }

  trtree->Fill();

  if(runNumber>0){
    TFile* foutfile=new TFile("trendingHF.root","recreate");
    trtree->Write();
    TDirectory* outdir=foutfile->mkdir(dirD2H->GetName());
    outdir->cd();
    listD2H->Write(listD2H->GetName(),1);
    foutfile->cd();
    if(dir && list){
      TDirectory* outdir2=foutfile->mkdir(dir->GetName());
      outdir2->cd();
      list->Write(list->GetName(),1);
    }
    foutfile->Close();
    delete foutfile;
  }
}
