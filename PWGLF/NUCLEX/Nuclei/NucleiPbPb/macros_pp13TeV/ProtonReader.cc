#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <fstream>


void ProtonReader(){
  TFile f("../results/ProtonYields.root");
  TGraphAsymmErrors* stat = (TGraphAsymmErrors*)f.Get("ProtonStat");
  Requires(stat, "ProtonStat");
  TGraphAsymmErrors* syst = (TGraphAsymmErrors*)f.Get("ProtonSyst");
  Requires(syst, "ProtonSyst");

  ofstream outtxt("../results/proton_yields.txt");

  int n_points = stat->GetN();
  for(int i=0; i<n_points; i++){
    double px,py;
    stat->GetPoint(i,px,py);
    double errx, staty, systy;
    errx = syst->GetErrorX(i);
    staty = stat->GetErrorY(i);
    systy = syst->GetErrorY(i);
    outtxt << px << " " << errx << " " << py << " " << staty << " " << systy << endl;
    printf("******************************************\n");
    printf("x: %f +- %f\n", px,py);
    printf("y: %f +- %f +- %f\n", py, staty, systy);
    printf("******************************************\n");
  }
  outtxt.close();
}
