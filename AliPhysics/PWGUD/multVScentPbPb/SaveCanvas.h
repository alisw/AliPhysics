/**************************************************************/
/*                                                            */
/* Set of methods to manipulate with canvas                   */
/* its elements etc.                                          */
/* Author:                                                    */
/* ruben.shahoyan@cern.ch                                     */
/*                                                            */
/**************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TString.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TList.h>
#include <TH1.h>
#include <TGraph.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TPaveStats.h>

extern TH1* GetBaseHisto(TPad* pad=0);
extern TFrame* GetFrame(TPad* pad=0);
extern void SetStatPad(TH1* hst,float x1,float x2,float y1,float y2);
extern TPaveStats* GetStatPad(TH1* hst);
extern void SetHStyle(TH1* hst,int col=kRed,int mark=20,float mrsize=0.7);
extern void SetGStyle(TGraph* hst,int col=kRed,int mark=20,float mrsize=0.7);
extern TH1* Cumulate(TH1* histo, Bool_t doErr=kTRUE, const char* copyName=0);
extern TLatex* AddLabel(const char*txt,float x=0.1,float y=0.9,int color=kBlack,float size=0.04);
extern void SaveCanvas(TCanvas* canv,const char* path="canv",const Option_t *option="ecg");
extern void wSum(double v1,double v2, double err1=0,double err2=0, double* wv=0,double *we=0);
extern void wAv(double v1,double v2, double err1=0,double err2=0, double* wv=0,double *we=0);

#endif
