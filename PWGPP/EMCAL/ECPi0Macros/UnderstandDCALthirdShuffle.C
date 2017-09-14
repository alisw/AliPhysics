

#include <TChain.h>
#include <TNtuple.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1F.h>
#include <TVector.h>
#include <TRefArray.h>
#include <TArrayS.h>
#include "TError.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGraphErrors.h"
#include "TPostScript.h"
#include "TLegend.h"
#include "TH2I.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TFile.h"
#include "TMath.h"
#include "TLeaf.h" 
#include "TBranch.h" 

//#include "/cebaf/faivre/recherche/utilities/defineMyPalette2011.C"

///
/// \file UnderstandDCALthirdShuffle.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief UnderstandDCALthirdShuffle
///
/// How to run :
///
///  .x UnderstandDCALthirdShuffle.C++
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///

///
/// Main
///
void UnderstandDCALthirdShuffle(void)
{int iSM,iCol,iRow,sm,col,row;
 double coeff,HV;
 double tabCalibCoeff[48*8],tabHV1[48*8],tabHV2[48*8],tabHVapplied[48*8],tabHVdesired[48*8];
 
 TCanvas *c1 = new TCanvas("c1","EMCal",550.,500.);
 c1->SetBorderMode(0);
 c1->SetFillColor(10); //Pas kWhite ni 0, sinon fait la bordure rouge dans la frame du plot !
 
  //defineMyPalette2011(30,5);
 
 TH1F *hCoeff = new TH1F("hCoeff","hCoeff",200,0.5,2.);
 hCoeff->SetStats(0);
 hCoeff->SetXTitle("Coeff");
 hCoeff->SetYTitle("Counts");
 TH1F *hHV_A = new TH1F("hHV_A","hHV_A",200,210.,400.);
 hHV_A->SetStats(0);
 hHV_A->SetXTitle("HV (A)");
 hHV_A->SetYTitle("Counts");
 TH1F *hHV_C = new TH1F("hHV_C","hHV_C",200,210.,400.);
 hHV_C->SetStats(0);
 hHV_C->SetXTitle("HV (C)");
 hHV_C->SetYTitle("Counts");
 TH1F *hHVdiff = new TH1F("hHVdiff","hHVdiff",200,-50.,50.);
 hHVdiff->SetStats(0);
 hHVdiff->SetXTitle("HV difference");
 hHVdiff->SetYTitle("Counts");
 TH2F *hCoeffVsHVratio = new TH2F("hCoeffVsHVratio","hCoeffVsHVratio",80,0.8,1.2,80,0.5,2.);
 //TH2F *hCoeffVsHVratio = new TH2F("hCoeffVsHVratio","hCoeffVsHVratio",80,-30.,30.,80,0.5,2.); //To try with HVdiff instead of HVratio.
 hCoeffVsHVratio->SetStats(0);
 hCoeffVsHVratio->SetXTitle("HV_A/HV_C");
 hCoeffVsHVratio->SetYTitle("Coeff");
 hCoeffVsHVratio->SetContour(30);
 
 FILE *fch = new FILE;
 fch = fopen("/cebaf/cebaf/EMCAL/calibPi0_run2/calibPi0_4_with2015data/output/pass2_DCALandThirdSMsVeryHighTowers/multiplyPi0CalibrationFactors_TextToHisto_forIter3_finalFile4HVcalculation.txt","r");
 for (int jm=0;jm<(10*48*24+6*32*24+4*48*8);jm++)
    {fscanf(fch," %d %d %d %lf\n",&sm,&col,&row,&coeff);
     if (sm != 18) continue;
     tabCalibCoeff[8*col+row]=coeff;
     hCoeff->Fill(coeff);
     //printf("(%d,%d) : %f\n",col,row,coeff);
     }
 fclose(fch);

 FILE *fch2 = new FILE;
 fch2 = fopen("/cebaf/faivre/recherche/calibPi0/recalculateDCAL_HV_July2015/output_HVrecalculationDCAL2015July_mergedWithEMCalHV_lowerTemperature_createBias/DCal/FinalBias_CN1.txt","r");
 for (int jm=0;jm<(48*24);jm++)
    {fscanf(fch," %d %d %lf\n",&col,&row,&HV);
     if (row<8)
        {tabHV1[8*col+row]=HV;
	 hHV_A->Fill(HV);
	 }
     if (row>=16)
        {tabHV2[8*col+row-16]=HV;
	 hHV_C->Fill(HV);
	 }
     }
 fclose(fch2);

 for (iCol=0;iCol<48;iCol++)
    {for (iRow=0;iRow<8;iRow++)
        {tabHVapplied[8*iCol+iRow]=tabHV1[8*iCol+iRow];
	 tabHVdesired[8*iCol+iRow]=tabHV2[8*iCol+iRow]; //Correct
	 //tabHVdesired[8*iCol+iRow]=tabHV1[8*(47-iCol)+(7-iRow)]; //Random-1
	 //tabHVdesired[8*iCol+iRow]=tabHV2[8*(47-iCol)+(7-iRow)]; //Random-2
	 if (tabCalibCoeff[8*iCol+iRow] == 1.) continue;
         hCoeffVsHVratio->Fill(tabHVapplied[8*iCol+iRow]/tabHVdesired[8*iCol+iRow],tabCalibCoeff[8*iCol+iRow]);
	   //hCoeffVsHVratio->Fill(tabHVapplied[8*iCol+iRow]-tabHVdesired[8*iCol+iRow],tabCalibCoeff[8*iCol+iRow]); //To look at HVdiff instead of HVratio.
	 hHVdiff->Fill(tabHV1[8*iCol+iRow]-tabHV2[8*iCol+iRow]);
         }
     }
 
 hCoeff->Draw();
 hHV_A->Draw();
 hHV_C->Draw();
 
 TLine *ligneH = new TLine(hCoeffVsHVratio->GetXaxis()->GetXmin(),1.,hCoeffVsHVratio->GetXaxis()->GetXmax(),1.);
 TLine *ligneV = new TLine(1.,hCoeffVsHVratio->GetYaxis()->GetXmin(),1.,hCoeffVsHVratio->GetYaxis()->GetXmax());
 hCoeffVsHVratio->SetMaximum(5.);

 hCoeffVsHVratio->Draw("COLZ");
 hCoeffVsHVratio->SetTitle(0);
 ligneH->Draw();
 ligneV->Draw();
 
 hHVdiff->SetTitle(0);
 hHVdiff->Rebin(2);
 hHVdiff->Draw();

 return;
 }






