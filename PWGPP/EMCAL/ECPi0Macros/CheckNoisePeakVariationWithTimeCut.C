///
/// \file CheckNoisePeakVariationWithTimeCut
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief study noise peak and time cut correlation
///
/// Study noise peak variation with time cut.
///
/// \author Astrid Vauthier, <Astrid.Vauthier@cern.ch>, (LPSC-CNRS)
///

#include <stdio.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1.h>
#include "TPostScript.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"
#include "TList.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TF1.h"
#include "TList.h"
#include "TLine.h"

//#include "/Users/vauthier/root/macros/myMacros/MyPaletteColor.c"               //------------->Have a nice "colz" plot
//#include "/Users/vauthier/root/macros/myMacros/ErrorPropagation.c"             //------------->Implemented formula for error propagation
//#include "/Users/vauthier/root/macros/myMacros/SettingsUtils.cpp"              //------------->Settings for canvas and histos


void CheckNoisePeakVariationWithTimeCut()
{
  //Declare needed variables - histos - files
  const int nTimeCut = 2;
  TString TimeCut[nTimeCut] = {"LooseTimeCut_NoNegCoeffs","_NoNegCoeffs"};
  const int nSM = 20;
  const int nCol = 48;
  const int nRow = 24;
  
  int nColForSM = 48;
  int nRowForSM = 24;
  
  TFile *file[nTimeCut];
  TList *list[nTimeCut];
  TH1F *hMgg[nSM][nCol][nRow][nTimeCut];
  
//  TFile *rootFileOut = new TFile("output_calibPi0.root","RECREATE");
  char psfile[100];
  sprintf(psfile,"NoisePeakVariationNoTimeCut.ps");
  const int cWidth=500;
  const int cHeight=(int)(500*(29./21.));
  TCanvas *c1 = new TCanvas("c1","Noise peaks variation as a function of time cuts",cWidth,cHeight);
  TPostScript *ps = new TPostScript(psfile,111);
  
  Color_t HistoColor[nTimeCut] = {kMagenta,kBlack};
  
  
  
  for(int iTimeCut = 0; iTimeCut < nTimeCut; iTimeCut++)
  {
    printf("--\nOPEN FILE %sTimeCutApplied\n",TimeCut[iTimeCut].Data());
    file[iTimeCut] = new TFile(Form("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/input/pass4/AnalysisResults_LHC15sumijMaskTowersByHand24Nov2016_pass4%s.root",TimeCut[iTimeCut].Data()));
    list[iTimeCut] = (TList*) file[iTimeCut]->Get("Pi0Calibration_Trig");
  }//End loop on time cut
  
  
  TFile *fOut = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass4_DCALandEMCALThirds/NoisePeakVariationNoTimeCut.root","RECREATE");
  
  for(int iSM = 0; iSM < nSM; iSM++)
  {
    printf("\nProcessing SM %i ...\n",iSM);
    if((iSM == 10) || (iSM == 11) || (iSM == 18) || (iSM == 19)) nRowForSM = 8;
    else nRowForSM = 24;
    if((iSM > 11) && ((iSM != 18) && (iSM != 19))) nColForSM = 32;
    else nColForSM = 48;
    for(int iCol = 0; iCol < nColForSM; iCol++)
    {
      for(int iRow = 0; iRow < nRowForSM; iRow++)
      {
        if(iRow%(nRow/2) == 0)
        {
          ps->NewPage();
          c1->Clear();
          c1->Divide(3,4);
        }
        c1->cd(iRow%(nRow/2)+1);
        for(int iTimeCut = 0; iTimeCut < nTimeCut; iTimeCut++)
        {
          hMgg[iSM][iCol][iRow][iTimeCut] = (TH1F*) list[iTimeCut]->FindObject(Form("%i_%i_%i",iSM,iCol,iRow));
            hMgg[iSM][iCol][iRow][iTimeCut]->SetLineColor(HistoColor[iTimeCut]);
            hMgg[iSM][iCol][iRow][iTimeCut]->SetMarkerColor(HistoColor[iTimeCut]);
            hMgg[iSM][iCol][iRow][iTimeCut]->SetXTitle("m_{#gamma#gamma} (MeV/c^{2})");
//          SetHisto(hMgg[iSM][iCol][iRow][iTimeCut],HistoColor[iTimeCut],20,"","m_{#gamma#gamma} (MeV/c^{2})","");
          if(iTimeCut == 0)hMgg[iSM][iCol][iRow][iTimeCut]->Draw("");
          else hMgg[iSM][iCol][iRow][iTimeCut]->Draw("same");
          hMgg[iSM][iCol][iRow][iTimeCut]->Write(Form("%i_%i_%i_%s",iSM,iCol,iRow,TimeCut[iTimeCut].Data()));
        }//End loop on time cut
        c1->Update();
      }//End loop on row
    }//End loop on column
    printf("... SM %i done\n",iSM);
  }//End loop on SM
  
  fOut->Write();
  fOut->Close();
  ps->Close();
  
  return;
}
