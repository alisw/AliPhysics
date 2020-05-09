#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TPaletteAxis.h"

#include <iostream>
#include <fstream>

Bool_t  savePlots = kTRUE;
Bool_t  writeSysErrors = kTRUE;


// --------------------------------------------------

void setHistoStyleColor(TH1* hist, Bool_t isSpectrum = kFALSE, Int_t color=2) {

  hist->SetLineColor(1);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.7);
  hist->SetMarkerColor(color);
  hist->SetTitleSize(0.08);
  
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06); 
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.0); 
	hist->GetYaxis()->SetTitleOffset(isSpectrum ? 0.8 : 0.7); 

  hist->GetXaxis()->SetTitleFont(62); 
  hist->GetYaxis()->SetTitleFont(62); 
  hist->GetXaxis()->SetLabelFont(62); 
  hist->GetYaxis()->SetLabelFont(62); 


  /*
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetTitleSize(0.04); 
  hist->GetYaxis()->SetTitleOffset(1.0); 
  */
}

TLegend* createLegend(Int_t type)
{
  TLegend* leg = 0x0;
  if (type == 0)
    leg = new TLegend(0.55,0.20,0.96,0.41);
  else if (type == 1)
    leg = new TLegend(0.49,0.63,0.83,0.84);
  else
    leg = new TLegend(0.55,0.61,0.96,0.83);
  
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->SetBorderSize(1);
  
  return leg;
}

// ---------------------------------------------------


TCanvas* createCanvas(TString name, Bool_t isSpectrum = kFALSE)
{
  TCanvas* c = new TCanvas(name.Data(),"",760,420); 
  c->Divide(2,2);
  
  for (Int_t i = 1; i <= 4; i++) {
    c->GetPad(i)->SetLeftMargin(isSpectrum ? 0.1 : 0.081);
    c->GetPad(i)->SetRightMargin(0.01);
    c->GetPad(i)->SetBottomMargin(0.13);
    c->GetPad(i)->SetTopMargin(0.12);
  }
  
  return c;
}

// ---------------------------------------------------

Int_t getColor(Int_t col){
  
  Int_t color = kBlack;

  switch(col)
    {
    case 0 : color = kBlue;    break;
    case 1 : color = kCyan;    break;
    case 2 : color = kGreen;   break;
    case 3 : color = kOrange;  break;
    case 4 : color = kRed;     break;
    case 5 : color = kPink;    break;
    case 6 : color = kMagenta; break;
    case 7 : color = kViolet;  break;
    case 8 : color = kBlue-2;  break;
    case 9 : color = kGreen-2;  break;
    case 10: color = kRed+2;  break;
    }

  return color;
  
}

// -------------------------------------------------------------------------
TList* getList(TString name = "jets_noGenJets_trackTypeUndef_jetTypeUndef"){

  TList* list = 0x0; 

  if(gDirectory->Get("obusch_list")){
    list = (TList*) gDirectory->Get("obusch_list");
  }
  else{
    
    TString dirName  = "PWGJE_FragmentationFunction_" + name;
    TString listName = "fracfunc_" + name;
    
    gDirectory->cd(dirName);
    list = (TList*) gDirectory->Get(listName);
        //gDirectory->cd("PWG4_FragmentationFunction_jets_noGenJets_KINEb_KINEb");
    //list = (TList*) gDirectory->Get("fracfunc_jets_noGenJets_KINEb_KINEb");
  }

  return list;
}

// ----------------------------------------------------------------------------

void sysErrorsPythiaFastJet_new(TString fileDir = "files", TString saveDir = "files"){
  
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleH(0.06);//0.054);
  
  const Int_t nModes = 5;
  
  const Bool_t useModes[nModes] = {kTRUE, kTRUE, kTRUE, kTRUE, kTRUE};
  TString modeString[nModes] = {"TrackPt", "Z", "Xi", "R", "jT"};
	const Bool_t useLogX[nModes] = {kTRUE, kFALSE, kFALSE, kFALSE, kTRUE};
  TString xAxeTitles[nModes] = {"#it{p}_{T} (GeV/#it{c})", "#it{z}", "#it{#xi}", "R", "j_{T} (GeV/#it{c})"};

  const Int_t nJetPtBins = 5;
  Double_t jetPtLim[nJetPtBins+1] = {5,10,15,20,30,80}; // nBins+1 entries
  
  TString strSp[] = {"","_pi","_K","_p","_e","_mu"}; 
  const Int_t nSpecies   = 5;

  TString strTitSp[] = {"h^{+} + h^{-}","#pi^{+} + #pi^{-}","K^{+} + K^{-}","p + #bar{p}","e^{+} + e^{-}"};//,"#mu^{+} + #mu^{-}"}; 

  TString strSp_10f6a[] = {"","_pi","_K","_p","_e"};//,"_mu"}; 
  
  TH1F* hSysErrEff[nModes][nSpecies][nJetPtBins];
  
  TH1F* hSysErrRes[nModes][nSpecies][nJetPtBins];

  TH1F* fh1FFGenPrim[nModes][nSpecies][nJetPtBins];
  
  const Int_t nVar = 2;
  const Bool_t useVariations[nVar] = {kTRUE, kTRUE};
  
  TString legendEntry[nVar] = {"Efficiency +/- 5%", "Resolution +/- 20%"};
  TString nameVar[nVar] = {"Eff", "Res"};
  
  TH1F* fh1FFRecPrim_recPt[2*nVar+1][nModes][nSpecies][nJetPtBins];
  
  TH1F* corrFacOriginalMC[nModes][nSpecies][nJetPtBins] = {0x0};
  
  TH1F* corrFac[2*nVar+1][nModes][nSpecies][nJetPtBins];

  TH1F* corrFacSys[nVar][nModes][nSpecies][nJetPtBins];

  //TString strInFile10f6a = "files/outCorrections_10f6a_tpcCut.root";
  TString strOriginalMCResults = fileDir + "/corrections_LHC13b2_efix_p1.root";
  
  TString strResDir = fileDir;

  TString strInFileGen(Form("%s/outCorrections_PythiaFastJet_Eff100_Res100.root",strResDir.Data()));
  
  TString strInFile_100_100(Form("%s/outCorrections_PythiaFastJet_Eff100_Res100.root",strResDir.Data()));
  TString strInFile_095_100(Form("%s/outCorrections_PythiaFastJet_Eff095_Res100.root",strResDir.Data()));
  TString strInFile_105_100(Form("%s/outCorrections_PythiaFastJet_Eff105_Res100.root",strResDir.Data()));
  TString strInFile_100_080(Form("%s/outCorrections_PythiaFastJet_Eff100_Res080.root",strResDir.Data()));
  TString strInFile_100_120(Form("%s/outCorrections_PythiaFastJet_Eff100_Res120.root",strResDir.Data()));
  
  TString strInFile[2*nVar+1] = {strInFile_100_100, strInFile_095_100, strInFile_105_100, strInFile_100_080, strInFile_100_120};
  
  Int_t varPercent[2*(2*nVar+1)] = {100, 100, 95, 100, 105, 100, 100, 80, 100, 120};
 
  // Load original MC results

  TFile f1(strOriginalMCResults,"READ");
  if (f1.IsZombie()) {
    cout << "Could not open full simulation file " << strOriginalMCResults << ", no comparison shown." << endl;
  }
  else {
    for (Int_t sp=0; sp<nSpecies; sp++) {
      for (Int_t i=0; i<nJetPtBins; i++) {
        for (Int_t mode=0;mode<nModes;++mode) {
          if (!useModes[mode])
            continue;
          
          TString strTitle(Form("hBbBCorr%s%s_%02d_%02d",modeString[mode].Data(),strSp_10f6a[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));        
          corrFacOriginalMC[mode][sp][i] = (TH1F*) gDirectory->Get(strTitle);
          corrFacOriginalMC[mode][sp][i]->SetDirectory(0);
        }
      }
    }   
    f1.Close();  
  }
  
  //TODO: catch files not present
  // Load fast MC results particle level

  TFile f2(strInFileGen,"READ");
  for (Int_t sp=0; sp<nSpecies; sp++) {
    gDirectory->cd(Form("Gen%s",strSp[sp].Data()));

    for(Int_t i=0; i<nJetPtBins; i++){
    
      // for genLevel doesn't matter which histos we take - eff == 1 
      for (Int_t mode=0;mode<nModes;++mode) {
        if (!useModes[mode])
          continue;
        
        TString strTitle(Form("fh1FF%sGen%s_%02d_%02d",modeString[mode].Data(),strSp_10f6a[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
        
        fh1FFGenPrim[mode][sp][i] = (TH1F*) gDirectory->Get(strTitle);
        
        fh1FFGenPrim[mode][sp][i]->SetDirectory(0);
      }
    }
    gDirectory->cd("..");
  }
  
  f2.Close();


  // ---
  
  for (Int_t var=0;var<2*nVar+1;var++) {
    if (var != 0 && !useVariations[(var-1)/2]) 
      continue;
          
    TFile f(strInFile[var].Data(),"READ");

    for (Int_t sp=0; sp<nSpecies; sp++) {
      
      gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));
      
      for(Int_t i=0; i<nJetPtBins; i++){
        
        for (Int_t mode=0;mode<nModes;++mode) {
          if (!useModes[mode])
            continue;
          
          TString strTitleRec(Form("fh1FF%sRecCuts%s_%02d_%02d",modeString[mode].Data(),strSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
          
          TH1F* fInput = (TH1F*) gDirectory->Get(strTitleRec);
          
          // -- corr factors
          
          TString strTit(Form("corrFac%s_%03d_%03d_%s_%d",modeString[mode], varPercent[2*var], varPercent[2*var + 1], strSp[sp].Data(), i));
          
          corrFac[var][mode][sp][i] = (TH1F*) fh1FFGenPrim[mode][sp][i]->Clone(strTit);
          corrFac[var][mode][sp][i]->Divide(fh1FFGenPrim[mode][sp][i],fInput,1,1,"B");
          
          /*********/
          
          delete fInput;
          fInput = 0x0;
          
          corrFac[var][mode][sp][i]->SetDirectory(0);  
        }
        
      }
      gDirectory->cd("..");
    }               
    f.Close();
  }


  // --------------------------------------------------------

  // calc syst. error as difference between corr factors 
  // store them as bin error in histos corrFacPt(Xi)SysEff / corrFacPt(Xi)SysRes . 
  
  for(Int_t sp=0; sp<nSpecies; sp++) {
    for(Int_t i=0; i<nJetPtBins; i++) {
      for (Int_t mode=0; mode<nModes; mode++) {
        if (!useModes[mode])
          continue;      
        
        for (Int_t var=0; var<nVar; var++) {
          if (!useVariations[var])
            continue;
          
          TString strTitSys(Form("corrFac%sSys%s%s_%d", modeString[mode].Data(), nameVar[var], strSp[sp].Data(),i));
          
          corrFacSys[var][mode][sp][i] = (TH1F*) corrFac[0][mode][sp][i]->Clone(strTitSys);
          
          for(Int_t bin=1; bin<=corrFac[0][mode][sp][i]->GetNbinsX(); bin++) {
            Double_t err = 0.5*TMath::Abs(corrFac[var*nVar+1][mode][sp][i]->GetBinContent(bin) - corrFac[var*nVar+2][mode][sp][i]->GetBinContent(bin)); 
            corrFacSys[var][mode][sp][i]->SetBinError(bin,err); 
          }
        }
      }
    }
  }

  // --------------------------------------------------------

  gStyle->SetOptStat(0);
  
  //gStyle->SetTitleStyle(1001);      
  //gStyle->SetTitleBorderSize(1);

  Int_t selectBin = 0;
  
  const Int_t nSpeciesPlot = 4;

  // ---
  
  const Int_t nOfCanvases = nVar * nModes;
  
  TCanvas* c[nOfCanvases];
  TLegend* leg[nOfCanvases];
  
  for (Int_t mode=0;mode<nModes;mode++) {
    if (!useModes[mode])
      continue;      
    for (Int_t var=0;var<nVar;var++) {
      if (!useVariations[var])
        continue;

      c[nVar*mode + var] = createCanvas(Form("c%d",(nVar*mode + var)+1));
      leg[nVar*mode + var] = createLegend((nVar*mode + var)/2);
    }
  }
      
  for (Int_t mode=0;mode<nModes;mode++) {
    if (!useModes[mode])
      continue;  
    for(Int_t sp=0; sp<nSpeciesPlot; sp++) {
      for(Int_t i=0; i<nJetPtBins-1; i++) {
        if (i != selectBin)
          continue;
        
        TString strPlotTitle(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}",strTitSp[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
        
        corrFac[0][mode][sp][i]->SetTitle(strPlotTitle.Data());
        corrFac[0][mode][sp][i]->SetXTitle(xAxeTitles[mode].Data());
        corrFac[0][mode][sp][i]->SetYTitle("Correction Factor");
        setHistoStyleColor(corrFac[0][mode][sp][i],4);
				corrFac[0][mode][sp][i]->GetXaxis()->SetRangeUser(0.0, corrFacOriginalMC[mode][sp][i]->GetBinLowEdge(corrFacOriginalMC[mode][sp][i]->FindLastBinAbove() + 1));
				
				Double_t minY = TMath::Min(corrFacOriginalMC[mode][sp][i]->GetBinContent(corrFacOriginalMC[mode][sp][i]->GetMinimumBin()), corrFac[0][mode][sp][i]->GetBinContent(corrFac[0][mode][sp][i]->GetMinimumBin()));
				Double_t maxY = TMath::Max(corrFacOriginalMC[mode][sp][i]->GetBinContent(corrFacOriginalMC[mode][sp][i]->GetMaximumBin()), corrFac[0][mode][sp][i]->GetBinContent(corrFac[0][mode][sp][i]->GetMaximumBin()));
				
				corrFac[0][mode][sp][i]->GetYaxis()->SetRangeUser(0.7*minY, TMath::Min(3.0, 1.3*maxY));
				
//         if (mode == 0) {
//           corrFac[0][mode][sp][i]->GetXaxis()->SetRangeUser(0,jetPtLim[i+1]);
//           corrFac[0][mode][sp][i]->GetYaxis()->SetRangeUser(0.41,2.2);
//         }
//         else if (mode == 1) {
//           corrFac[0][mode][sp][i]->GetXaxis()->SetRangeUser(0,1.0);
//           corrFac[0][mode][sp][i]->GetYaxis()->SetRangeUser(0.35,2.2);
//           if(sp>=2) corrFac[0][mode][sp][i]->GetYaxis()->SetRangeUser(0.35,4.2);          
//         }
//         else if (mode == 2) {
//           corrFac[0][mode][sp][i]->GetXaxis()->SetRangeUser(0,6);
//           corrFac[0][mode][sp][i]->GetYaxis()->SetRangeUser(0.61,2.2);
//         }
        
        corrFac[0][mode][sp][i]->SetLineColor(4); 
        
        for (Int_t var=0;var<nVar;var++) {
          if (!useVariations[var])
            continue;
          
          c[nVar*mode + var]->cd(sp+1)->SetLogx(useLogX[mode]);
          corrFac[0][mode][sp][i]->DrawCopy();
        
          setHistoStyleColor(corrFacSys[var][mode][sp][i],4);
          corrFacSys[var][mode][sp][i]->SetFillColor(7);
          corrFacSys[var][mode][sp][i]->DrawCopy("same,E2");
          
          if (corrFacOriginalMC[mode][sp][i]) {
            setHistoStyleColor(corrFacOriginalMC[mode][sp][i],2);
            corrFacOriginalMC[mode][sp][i]->SetMarkerStyle(24);
            corrFacOriginalMC[mode][sp][i]->DrawCopy("same");  
          }
              
          if(sp==0){
            if (corrFacOriginalMC[mode][sp][i])
              leg[nVar*mode + var]->AddEntry(corrFacOriginalMC[mode][sp][i],"Full simulation","P");
            
            leg[nVar*mode + var]->AddEntry(corrFac[0][mode][sp][i],"Fast simulation","P");
            leg[nVar*mode + var]->AddEntry(corrFacSys[var][mode][sp][i],legendEntry[var].Data(),"F");
            leg[nVar*mode + var]->Draw();
          }
          
          gPad->RedrawAxis("");
          gPad->RedrawAxis("G");
        }
      }
    }
  }
  

  
  if (savePlots) {
    for (Int_t mode=0;mode<nModes;mode++) {
      if (!useModes[mode])
        continue;    
      for (Int_t var=0;var<nVar;var++) {
        if (!useVariations[var])
          continue;
        c[nVar*mode + var]->SaveAs(Form("%scorrFacMod%s_dNd%s_%02d_%02d.pdf", saveDir.Data(), nameVar[var].Data(), modeString[mode].Data(), (int)jetPtLim[selectBin],(int)jetPtLim[selectBin+1]));
      }
    }
  }

  // ----
  // write syst. error to file
  
  TH1F* hSysErr[nVar][nModes][nSpecies][nJetPtBins];

  if (writeSysErrors) {

    for(Int_t sp=0; sp<nSpecies; sp++) {
      for(Int_t i=0; i<nJetPtBins; i++) { // jet slices
        Int_t jetPtLoLim = (int) jetPtLim[i];
        Int_t jetPtUpLim = (int) jetPtLim[i+1];

        for (Int_t var=0;var<nVar;var++) {
          if (!useVariations[var])
            continue;
          for (Int_t mode=0;mode<nModes;mode++) {
            if (!useModes[mode])
              continue;  
            TString strName(Form("hSysErr%s%s_%02d_%02d%s",nameVar[var].Data(), modeString[mode].Data(), jetPtLoLim,jetPtUpLim,strSp[sp].Data()));
            
            hSysErr[var][mode][sp][i] = (TH1F*) corrFac[0][mode][sp][i]->Clone(strName.Data());
            
            hSysErr[var][mode][sp][i]->Reset();
            
            for(Int_t bin=1; bin<=corrFacSys[var][mode][sp][i]->GetNbinsX(); bin++) {
              if (corrFacOriginalMC[mode][sp][i] && corrFacOriginalMC[mode][sp][i]->GetBinContent(bin) == 0) 
                continue; 
              
              Double_t sysErr = 0.0;
              
              Double_t err  = corrFacSys[var][mode][sp][i]->GetBinError(bin);
              Double_t cont = TMath::Abs(corrFacSys[var][mode][sp][i]->GetBinContent(bin));  
              
              if (cont)
                sysErr = err/cont;
              
              hSysErr[var][mode][sp][i]->SetBinContent(bin,sysErr); 
              hSysErr[var][mode][sp][i]->SetBinError(bin,0);      
              
            }
          }
        }
      }
    }
    
    for (Int_t var=0;var<nVar;var++) { 
      if (!useVariations[var])
        continue;  
      TString fname = Form("%soutSysErr_%s.root", saveDir.Data(),nameVar[var].Data());
      TFile fOut(fname,"RECREATE");
      
      cout <<" write to file " << fname << endl;
    
      for (Int_t sp=0; sp<nSpecies; sp++) {
        for (Int_t i=0; i<nJetPtBins; i++) {
          for (Int_t mode=0;mode<nModes;mode++) {  
            if (!useModes[mode])
              continue;  
            hSysErr[var][mode][sp][i]->Write();
          }
        }
      }
    
      fOut.Close();
    }
  } // writeSysErrors
  
	TString strInFile_LowPtEnhancement(Form("%s/outCorrections_PythiaFastJet_LowPtEnhancement.root",strResDir.Data()));
	TString strInFile_LowPtDepletion(Form("%s/outCorrections_PythiaFastJet_LowPtDepletion.root",strResDir.Data()));
  
  const Int_t nFFModes = 2;
	
	TString ffFiles[nFFModes] = {strInFile_LowPtEnhancement, strInFile_LowPtDepletion};
	TString ffModeString[nFFModes] = {"LowPtEnhancement", "LowPtDepletion"};
  
  TH1F* fh1FFGenPrimFFmode[nFFModes][nModes][nSpecies][nJetPtBins];
	TH1F* fh1FFRecPrimFFmode[nFFModes][nModes][nSpecies][nJetPtBins];
	TH1F* corrFacFFmode[nFFModes][nModes][nSpecies][nJetPtBins];
	TH1F* corrFacFFBbB[nModes][nSpecies][nJetPtBins];
	TH1F* hSysErrBbB[nModes][nSpecies][nJetPtBins];
	
	//TODO: catch files not present
  // Load fast MC results particle level for FF modes

	for (Int_t ffMode=0;ffMode<nFFModes;++ffMode) {
		TFile f(ffFiles[ffMode],"READ");
		for (Int_t sp=0; sp<nSpecies; sp++) {
			gDirectory->cd(Form("Gen%s",strSp[sp].Data()));

			for(Int_t i=0; i<nJetPtBins; i++){
			
				// for genLevel doesn't matter which histos we take - eff == 1 
				for (Int_t mode=0;mode<nModes;++mode) {
					if (!useModes[mode])
						continue;
					
					TString strTitle(Form("fh1FF%sGen%s_%02d_%02d",modeString[mode].Data(),strSp_10f6a[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
					
					fh1FFGenPrimFFmode[ffMode][mode][sp][i] = (TH1F*) gDirectory->Get(strTitle);
					
					fh1FFGenPrimFFmode[ffMode][mode][sp][i]->SetDirectory(0);
				}
			}
			gDirectory->cd("..");
		}
		for (Int_t sp=0; sp<nSpecies; sp++) {
			gDirectory->cd(Form("RecCuts%s",strSp[sp].Data()));

			for(Int_t i=0; i<nJetPtBins; i++){
			
				for (Int_t mode=0;mode<nModes;++mode) {
					if (!useModes[mode])
						continue;
					
					TString strTitle(Form("fh1FF%sRecCuts%s_%02d_%02d",modeString[mode].Data(),strSp_10f6a[sp].Data(),(int)jetPtLim[i],(int)jetPtLim[i+1]));
					
					fh1FFRecPrimFFmode[ffMode][mode][sp][i] = (TH1F*) gDirectory->Get(strTitle);
					
					fh1FFRecPrimFFmode[ffMode][mode][sp][i]->SetDirectory(0);
				}
			}
			gDirectory->cd("..");
		}
		f.Close();
		for (Int_t sp=0; sp<nSpecies; sp++) {
			for(Int_t i=0; i<nJetPtBins; i++){
				for (Int_t mode=0;mode<nModes;++mode) {
					if (!useModes[mode])
						
						continue;
					
					TString strTit(Form("corrFac%s_%s_%s_%d", modeString[mode], ffModeString[ffMode], strSp[sp].Data(), i));
					
					corrFacFFmode[ffMode][mode][sp][i] = (TH1F*) fh1FFGenPrimFFmode[ffMode][mode][sp][i]->Clone();
					corrFacFFmode[ffMode][mode][sp][i]->Divide(fh1FFGenPrimFFmode[ffMode][mode][sp][i],fh1FFRecPrimFFmode[ffMode][mode][sp][i],1,1,"B");
				}
			}
		}
	}

	for (Int_t sp=0; sp<nSpecies; sp++) {
		for(Int_t i=0; i<nJetPtBins; i++){
			for (Int_t mode=0;mode<nModes;++mode) {
				if (!useModes[mode])
					continue;
				
				TString strNameCorrSysErrBbB(Form("corrFac%sSysBbB_%02d_%02d%s", modeString[mode].Data(), (int)jetPtLim[i], (int)jetPtLim[i+1], strSp[sp].Data()));
				TString strNamehSysErrBbB(Form("hSysErrBbB%s_%02d_%02d%s", modeString[mode].Data(), (int)jetPtLim[i], (int)jetPtLim[i+1], strSp[sp].Data()));
				
				corrFacFFBbB[mode][sp][i] = (TH1F*) corrFac[0][mode][sp][i]->Clone(strNameCorrSysErrBbB);
				hSysErrBbB[mode][sp][i] =  (TH1F*) corrFac[0][mode][sp][i]->Clone(strNamehSysErrBbB);
				hSysErrBbB[mode][sp][i]->Reset();
				
				for (Int_t bin=1;bin<=corrFac[0][mode][sp][i]->GetNbinsX();++bin) {
					corrFacFFBbB[mode][sp][i]->SetBinError(bin, 0);
					
					if (fh1FFGenPrim[mode][sp][i]->GetBinContent(bin) == 0)
						continue;
					
					Double_t minEl = corrFacFFBbB[mode][sp][i]->GetBinContent(bin);
					Double_t maxEl = minEl;
					
					for (Int_t ffMode=0;ffMode<nFFModes;ffMode++) {
						minEl = TMath::Min(minEl, corrFacFFmode[ffMode][mode][sp][i]->GetBinContent(bin));
						maxEl = TMath::Max(maxEl, corrFacFFmode[ffMode][mode][sp][i]->GetBinContent(bin));
					}
					
					Double_t sysErrBbB = 0; 
					if(minEl) 
						sysErrBbB = 0.5*fabs((maxEl/minEl) - 1);
					
					hSysErrBbB[mode][sp][i]->SetBinContent(bin,sysErrBbB);
					hSysErrBbB[mode][sp][i]->SetBinError(bin,0);
					
					corrFacFFBbB[mode][sp][i]->SetBinError(bin, sysErrBbB);
				}
				
			}
		}
	}
	
	TFile* fOutFFMode = new TFile(Form("%soutSysErr_BbB.root", saveDir.Data()), "RECREATE");
	
	for (Int_t sp=0; sp<nSpecies; sp++) {
		for (Int_t i=0; i<nJetPtBins; i++) {
			for (Int_t mode=0;mode<nModes;mode++) {  
				if (!useModes[mode] || !hSysErrBbB[mode][sp][i])
					continue;  
				
				hSysErrBbB[mode][sp][i]->Write();
			}
		}
	}
	
	fOutFFMode->Close();
	
	gStyle->SetOptStat(0);
  
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleH(0.06);
	
	const Int_t ffModesColor[nFFModes] = {51, 8};
	const TString ffModesNames[nFFModes] = {"Low-pt enhanced", "Low-pt depleted"};
	
	for (Int_t mode=0;mode<nModes;mode++) {
		Int_t selectBin = 0;

		TCanvas* c1 = createCanvas(Form("c%d", mode), kTRUE);
		TLegend* leg1 = createLegend(0);
		
		TCanvas* c2 = createCanvas(Form("c%d", mode + nModes));
		TLegend* leg2 = createLegend(2);
		
		const Int_t nSpeciesPlot = 4;
		
		for(Int_t sp=0; sp<nSpeciesPlot; sp++){
			for(Int_t i=0; i<nJetPtBins; i++){

				if(i != selectBin) 
					continue;
			
				c1->cd(sp+1)->SetLogx(useLogX[mode]);

				TString strPlotTitle(Form("%s, #it{p}_{T}^{jet, ch} = %d-%d GeV/#it{c}", strTitSp[sp].Data(), (int)jetPtLim[i], (int)jetPtLim[i+1]));
				
				fh1FFGenPrim[mode][sp][i]->SetTitle(strPlotTitle);
				fh1FFGenPrim[mode][sp][i]->SetXTitle(xAxeTitles[mode].Data());
				fh1FFGenPrim[mode][sp][i]->SetYTitle(Form("1/#it{N}_{Jets} d#it{N}/d#it{%s}", xAxeTitles[mode].Data()));
				setHistoStyleColor(fh1FFGenPrim[mode][sp][i], kTRUE, 2);
				Int_t lastBinWithContent = fh1FFGenPrim[mode][sp][i]->FindLastBinAbove(0.0);
				fh1FFGenPrim[mode][sp][i]->GetXaxis()->SetRange(1,lastBinWithContent);

				Double_t maxValue = fh1FFGenPrim[mode][sp][i]->GetBinContent(fh1FFGenPrim[mode][sp][i]->GetMaximumBin());

				for (Int_t ffMode=0;ffMode<nFFModes;++ffMode) {
					maxValue = TMath::Max(maxValue, fh1FFGenPrimFFmode[ffMode][mode][sp][i]->GetBinContent(fh1FFGenPrimFFmode[ffMode][mode][sp][i]->GetMaximumBin()));
				}
				
				fh1FFGenPrim[mode][sp][i]->GetYaxis()->SetRangeUser(0,maxValue *1.5);
				fh1FFGenPrim[mode][sp][i]->DrawCopy();  
				for (Int_t ffMode=0;ffMode<nFFModes;++ffMode) {
					setHistoStyleColor(fh1FFGenPrimFFmode[ffMode][mode][sp][i],kTRUE, ffModesColor[ffMode]);
					fh1FFGenPrimFFmode[ffMode][mode][sp][i]->DrawCopy("same");
				}
				fh1FFGenPrim[mode][sp][i]->DrawCopy("same"); //redraw on top   
						
				if(sp==0){
					leg1->AddEntry(fh1FFGenPrim[mode][sp][i],"Reference, gen level","P");
					for (Int_t ffMode=0;ffMode<nFFModes;++ffMode) {
						leg1->AddEntry(fh1FFGenPrimFFmode[ffMode][mode][sp][i],ffModesNames[ffMode].Data(),"P");
					}
					
					leg1->Draw();
				}
				
				gPad->RedrawAxis("");
				gPad->RedrawAxis("G");
				
				c2->cd(sp+1)->SetLogx(useLogX[mode]);

				corrFac[0][mode][sp][i]->SetTitle(strPlotTitle);
				corrFac[0][mode][sp][i]->SetXTitle(xAxeTitles[mode].Data());
				corrFac[0][mode][sp][i]->SetYTitle("Correction Factor");
				setHistoStyleColor(corrFac[0][mode][sp][i],kFALSE,4);
				lastBinWithContent = corrFac[0][mode][sp][i]->FindLastBinAbove(0.0);
				corrFac[0][mode][sp][i]->GetXaxis()->SetRange(1,lastBinWithContent);
// 				corrFac[0][mode][sp][i]->GetXaxis()->SetRangeUser(0,jetPtLim[i+1]);
				corrFac[0][mode][sp][i]->SetLineColor(4); 
				maxValue = corrFac[0][mode][sp][i]->GetBinContent(corrFac[0][mode][sp][i]->GetMaximumBin());
				
				setHistoStyleColor(corrFacFFBbB[mode][sp][i],kFALSE,4);
				corrFacFFBbB[mode][sp][i]->SetFillColor(7);
				maxValue = TMath::Max(maxValue, corrFacFFBbB[mode][sp][i]->GetBinContent(corrFacFFBbB[mode][sp][i]->GetMaximumBin()));
				
				setHistoStyleColor(corrFacOriginalMC[mode][sp][i],kFALSE,2);
				corrFacOriginalMC[mode][sp][i]->SetMarkerStyle(24);
				maxValue = TMath::Max(maxValue, corrFacOriginalMC[mode][sp][i]->GetBinContent(corrFacOriginalMC[mode][sp][i]->GetMaximumBin()));
				
				//Do not use dynamic calculation here because correction can sometimes become very big
// 				corrFac[0][mode][sp][i]->GetYaxis()->SetRangeUser(0.0,maxValue * 1.2);
				corrFac[0][mode][sp][i]->GetYaxis()->SetRangeUser(0.0,3.0);
				corrFac[0][mode][sp][i]->DrawCopy();  
				corrFacFFBbB[mode][sp][i]->DrawCopy("same,E2");
				corrFacOriginalMC[mode][sp][i]->DrawCopy("same");
				
				if(sp==0){
					leg2->AddEntry(corrFacOriginalMC[mode][sp][i],"Full simulation","P");
					leg2->AddEntry(corrFac[0][mode][sp][i],"Fast simulation","P");
					leg2->AddEntry(corrFacFFBbB[mode][sp][i],"Low-pt enhanced/depleted","F");
					leg2->Draw();
				}
				
				gPad->RedrawAxis("");
				gPad->RedrawAxis("G");
				
				if (savePlots) {
					c1->SaveAs(Form("%sspectraModFF_dNd%s_%02d_%02d.pdf", saveDir.Data(), modeString[mode].Data(), (int) jetPtLim[i],(int) jetPtLim[i+1]));
					c2->SaveAs(Form("%scorrFacModFF_dNd%s_%02d_%02d.pdf", saveDir.Data(), modeString[mode].Data(), (int) jetPtLim[i],(int) jetPtLim[i+1]));
				}
					
			}
		}
	}
	
}
