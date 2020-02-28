#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include "Hyp3FindConfig.h"

void SetEfficiencyErrors(TH1D*,TH1D*);
void SetRateErrors(TH1D*,TH1D*);

void efficiency_analysis(char * input_name = "selector_results.root",char * output_name = "efficiency.root",char * folder_name = "plot_folder")
{
  gStyle->SetOptStat(0);
  TFile input_file(input_name);

  const char lAM[3]{"AM"};
  const char *lProj[2]{"pT", "ct"};
  const char *lProjSet[2]{"zx", "zy"};
  const char *lVarName[3]={"NsigmaTPC", "NclusterTPC","NclusterITS"};
  const char *lMeasure[4]={"Efficiency","Resolution","FakeRate","CloneRate"};
  //Gen must be the first of the list
  const char* lProjMeas[5] = {"Gen","Rec","Res","Fake","Clones"};

  gSystem->Exec(Form("mkdir %s",folder_name));
  for(auto VarName : lVarName){
    gSystem->Exec(Form("mkdir %s/%s/",folder_name,VarName));
    for(auto Measure : lMeasure)
      gSystem->Exec(Form("mkdir %s/%s/%s",folder_name,VarName,Measure));
  }

  TFile output_file(output_name,"RECREATE");
  output_file.cd();

  TCanvas cv("","");
  TCanvas Rainbow("","");

  TH3D* fHistVsCuts[5][kNVar][2];
  TH2D* fHistProjRec[5][kNVar][2];
  TH1D* fHist[5][kNVar][2][kNCut];

  //matter or antimatter
  for(int iMatter=0; iMatter<2; iMatter++){
    //projection on pt or ct
    for(int iProj=0; iProj<2; iProj++){
      //variable of the cut
      for(int iVar=0; iVar<kNVar; iVar++){
        //meaurements
        for(int iHist=0; iHist<5; iHist++){
          fHistVsCuts[iHist][iVar][iMatter] = (TH3D*) input_file.Get(Form("fHist%s_%s_%c",lProjMeas[iHist],lVarName[iVar],lAM[iMatter]));
          fHistProjRec[iHist][iVar][iMatter] = (TH2D*) fHistVsCuts[iHist][iVar][iMatter]->Project3D(lProjSet[iProj]);
          //different cuts
          for(int iCut=0; iCut<kNCut; iCut++){
            fHist[iHist][iVar][iMatter][iCut] = (TH1D*) fHistProjRec[iHist][iVar][iMatter]->ProjectionX(Form("fHist%s_%s_%s_cut_%.1f_%c",lProjMeas[iHist],lVarName[iVar],lProj[iProj],kCuts[iVar][iCut][0],lAM[iMatter]),iCut+1,iCut+1);  
              
            if(lProjMeas[iHist]=="Gen") continue;
            if(lProjMeas[iHist]=="Rec"){
              fHist[iHist][iVar][iMatter][iCut]->Divide(fHist[0][iVar][iMatter][iCut]);
              SetEfficiencyErrors(fHist[iHist][iVar][iMatter][iCut],fHist[0][iVar][iMatter][iCut]);
              fHist[iHist][iVar][iMatter][iCut]->GetYaxis()->SetTitle("efficiency");
            }else if(lProjMeas[iHist]=="Fake" || lProjMeas[iHist]=="Clones"){
              fHist[iHist][iVar][iMatter][iCut]->Divide(fHist[0][iVar][iMatter][iCut]);
              SetRateErrors(fHist[iHist][iVar][iMatter][iCut],fHist[0][iVar][iMatter][iCut]); 
              fHist[iHist][iVar][iMatter][iCut]->GetYaxis()->SetTitle(Form("#frac{#%s}{#Gen}",lProjMeas[iHist]));
            }
            else{
              fHist[iHist][iVar][iMatter][iCut]->Scale(1./fHist[iHist][iVar][iMatter][iCut]->GetEntries());                fHist[iHist][iVar][iMatter][iCut]->GetYaxis()->SetTitle("counts/N_{ev}");
            }
            fHist[iHist][iVar][iMatter][iCut]->SetMarkerStyle(8);
            fHist[iHist][iVar][iMatter][iCut]->SetTitle(Form("cut on %s: %.1f %c",lVarName[iVar],kCuts[iVar][iCut][0],lAM[iMatter]));
            cv.cd();
            fHist[iHist][iVar][iMatter][iCut]->Draw();
            cv.SaveAs(Form("%s/%s/%s/%s_%s_%s_cut_%.1f_%c.pdf",folder_name,lVarName[iVar],lMeasure[iHist-1],lMeasure[iHist-1],lVarName[iVar],lProj[iProj],kCuts[iVar][iCut][0],lAM[iMatter]));
              
            output_file.cd();
            fHist[iHist][iVar][iMatter][iCut]->Write();

            Rainbow.cd();
            if(iCut==0)
              fHist[iHist][iVar][iMatter][iCut]->Draw("PMC PLC");
            else
              fHist[iHist][iVar][iMatter][iCut]->Draw("SAME PMC PLC");

          }

          if(iHist!=0){
            Rainbow.SetName(Form("RainbowPlot_%s_%s_%s_%c",lVarName[iVar],lMeasure[iHist-1],lProj[iProj],lAM[iMatter])); 
            Rainbow.BuildLegend();
            Rainbow.SaveAs(Form("%s/%s/RainbowPlot_%s_%s_%s_%c.pdf",folder_name,lVarName[iVar],lVarName[iVar],lMeasure[iHist-1],lProj[iProj],lAM[iMatter])); 
            output_file.cd();
            Rainbow.Write();
          }
        }
      }  
    }
  }
}

void SetEfficiencyErrors(TH1D* HistEff, TH1D* HistGen){
  for(int iBin=1; iBin<=HistEff->GetNbinsX(); iBin++){
    double gen = HistGen->GetBinContent(iBin);
    double eff = HistEff->GetBinContent(iBin);
    if(gen!=0)
      HistEff->SetBinError(iBin,TMath::Sqrt(eff*(1-eff)/gen));
    else{
      HistEff->SetBinError(iBin,0);
      HistEff->SetBinContent(iBin,1);
    }
  }
}

void SetRateErrors(TH1D* HistRate, TH1D* HistGen){
  for(int iBin=1; iBin<=HistRate->GetNbinsX(); iBin++){
    double gen = HistGen->GetBinContent(iBin);
    double rate = HistRate->GetBinContent(iBin);
    if(gen!=0)
      HistRate->SetBinError(iBin,TMath::Sqrt(rate/gen));
    else{
      HistRate->SetBinError(iBin,0);
      HistRate->SetBinContent(iBin,-1);
    }
  }
}

