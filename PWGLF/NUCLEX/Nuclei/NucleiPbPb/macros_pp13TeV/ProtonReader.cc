#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <fstream>


void ProtonReader(){
  TFile f(Form("%sOutputYields_pp13mult.root",kBaseOutputDir.data()));
  TGraphAsymmErrors* stat = (TGraphAsymmErrors*)f.Get("p_Stat");
  Requires(stat, "ProtonStat");
  TGraphAsymmErrors* syst = (TGraphAsymmErrors*)f.Get("p_Syst");
  Requires(syst, "ProtonSyst");
  TGraphAsymmErrors* syst_uncorr = (TGraphAsymmErrors*)f.Get("p_SystUnc");
  Requires(syst_uncorr, "ProtonSystUnc");

  std::ofstream outtxt(Form("%sproton_yields.txt",kBaseOutputDir.data()));

  int n_points = stat->GetN();
  for(int i=0; i<n_points; i++){
    double px,py;
    stat->GetPoint(i,px,py);
    double errx, staty, systy, systy_uncorr, systy_corr;
    errx = syst->GetErrorX(i);
    staty = stat->GetErrorY(i);
    systy = syst->GetErrorY(i);
    systy_uncorr = syst_uncorr->GetErrorY(i);
    systy_corr = TMath::Sqrt(Sq(systy)-Sq(systy_uncorr));
    outtxt << px << " " << errx << " " << py << " " << staty << " " << systy << " " << systy_uncorr << " " << systy_corr << std::endl;
    printf("******************************************\n");
    printf("x: %f +- %f\n", px, errx);
    printf("y: %f +- %f +- %f (+- %f +-%f)\n", py, staty, systy, systy_uncorr, systy_corr);
    printf("******************************************\n");
  }
  outtxt.close();
}
