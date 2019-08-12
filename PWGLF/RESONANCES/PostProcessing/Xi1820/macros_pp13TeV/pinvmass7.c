//Use of polynomial to determine systematic error
//also, determine if chi^2 is optimized
//uses c++ style
#include <Riostream.h>
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
#include "TGraph.h"
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
#include <iostream>
#include "overhead.C"
void plot_mass(char* dc);

class AllFunctions{
public:
    double function(double *x, double *par)
    {
        Double_t arg1 = 14.0/22.0; // 2 over pi
        Double_t arg2 = par[3]*par[3]*par[4]*par[4]; //Gamma=par[1]  M=par[2]
        Double_t arg3 = ((x[0]*x[0]) - (par[4]*par[4]))*((x[0]*x[0]) - (par[4]*par[4]));
        Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[3]*par[3])/(par[4]*par[4]));
        return par[2]*arg1*arg2/(arg3 + arg4)+par[0]+par[1]*x[0];
    }
    
    double fline(double *x, double *par)
    {
        if ((x[0] > (1.823-4.0*(0.024/2.35))) && (x[0] < (1.823+4.0*(0.024/2.35)))) {//x[0] > 1.78 &&  x[0] < 1.87
            TF1::RejectPoint();
            return 0;
        }
        return par[0] + par[1]*x[0];
    }
    
    double BWfunction(double *x, double *par)
    {
        Double_t arg1 = 2.0/TMath::Pi(); // 2 over pi
        Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
        Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
        Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
        return par[0]*arg1*arg2/(arg3 + arg4);
    }
    
    double funVoigt(double *x, double *par)
    {
        return par[0]*TMath::Voigt(x[0]-par[1],par[2],par[3],4);
        //return par[0]*TMath::Voigt(x[0]-par[1],par[2],par[3],4)*(1+0.15*par[4]*TMath::Erf((x[0]-par[1])/(0.15*1.1283416)));
    }
    
    double funVoigtMC(double *x, double *par)
    {
        return par[0]*TMath::Voigt(x[0]-par[1],par[2],par[3],4)*(1+0.15*par[4]*TMath::Erf((x[0]-par[1])/(0.15*1.1283416)));
    }
    
    double plainVoigt(double *x, double *par)
    {
        return par[0]*TMath::Voigt(x[0]-par[1],par[2],par[3],4);
    }
    
    double backline(double *x, double *par)
    {
        return par[0] + par[1]*x[0];
    }
    
    double function2nd(double *x, double *par)
    {
        double arg1 = 14.0/22.0; // 2 over pi
        double arg2 = par[4]*par[4]*par[5]*par[5]; //Gamma=par[1]  M=par[2]
        double arg3 = ((x[0]*x[0]) - (par[5]*par[5]))*((x[0]*x[0]) - (par[5]*par[5]));
        double arg4 = x[0]*x[0]*x[0]*x[0]*((par[4]*par[4])/(par[5]*par[5]));
        return par[3]*arg1*arg2/(arg3 + arg4)+par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
    }
    
    double f2nd(double *x, double *par)
    {
        if ((x[0] > (1.823-4.0*(0.024/2.35))) && (x[0] < (1.823+4.0*(0.024/2.35)))) {///*reject &&*/ x[0] > 1.78 && x[0] < 1.87/* && flagXi1820 == 1*/
            TF1::RejectPoint();
            return 0;
        }
        return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
    }
    
    double back2nd(double *x, double *par)
    {
        return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
    }
    
    double function3rd(double *x, double *par)
    {
        Double_t arg1 = 14.0/22.0; // 2 over pi
        Double_t arg2 = par[5]*par[5]*par[6]*par[6]; //Gamma=par[1]  M=par[2]
        Double_t arg3 = ((x[0]*x[0]) - (par[6]*par[6]))*((x[0]*x[0]) - (par[6]*par[6]));
        Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[5]*par[5])/(par[6]*par[6]));
        return par[4]*arg1*arg2/(arg3 + arg4)+par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0];
    }
    
    double f3rd(double *x, double *par)
    {
        if ((x[0] > (1.823-4.0*(0.024/2.35))) && (x[0] < (1.823+4.0*(0.024/2.35)))) {
            TF1::RejectPoint();
            return 0;
        }
        return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
    }
    
    double back3rd(double *x, double *par)
    {
        return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
    }
    
    //4th
    double function4th(double *x, double *par)
    {
        Double_t arg1 = 14.0/22.0; // 2 over pi
        Double_t arg2 = par[6]*par[6]*par[7]*par[7]; //Gamma=par[1]  M=par[2]
        Double_t arg3 = ((x[0]*x[0]) - (par[7]*par[7]))*((x[0]*x[0]) - (par[7]*par[7]));
        Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[6]*par[6])/(par[7]*par[7]));
        return par[5]*arg1*arg2/(arg3 + arg4)+par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
    }
    
    double f4th(double *x, double *par)
    {
        if ((x[0] > (1.823-4.0*(0.024/2.35))) && (x[0] < (1.823+4.0*(0.024/2.35)))) {
            TF1::RejectPoint();
            return 0;
        }
        return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
    }
    
    double back4th(double *x, double *par)
    {
        return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
    }
    //4th
    
    double functionNewK0(double *x, double *par)
    {
        
        Double_t aug5=1.115683+.497611;//mass os Lambda + K0 in GeV/c^2
        Double_t arg1 = 14.0/22.0; // 2 over pi
        Double_t arg2 = par[6]*par[6]*par[7]*par[7]; //Gamma=par[1]  M=par[2]
        Double_t arg3 = ((x[0]*x[0]) - (par[7]*par[7]))*((x[0]*x[0]) - (par[7]*par[7]));
        Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[6]*par[6])/(par[7]*par[7]));
        //return par[4]*arg1*arg2/(arg3 + arg4)+(aug5-x[0])*par[0]*TMath::Exp(par[1]+par[2]*aug5+par[3]*aug5*aug5);
        return par[5]*arg1*arg2/(arg3 + arg4)+(x[0]-aug5)*par[0]*TMath::Exp(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])+par[4];
        
    }
    
    double fNewK0(double *x, double *par)
    {
        if ((x[0] > (1.823-4.0*(0.024/2.35))) && (x[0] < (1.823+4.0*(0.024/2.35)))) {
            TF1::RejectPoint();
            return 0;
        }
        
        Double_t aug5=1.115683+.497611;//mass os Lambda + K0 in GeV/c^2
        //return (aug5-x[0])*par[0]*TMath::Exp(par[1]+par[2]*aug5+par[3]*aug5*aug5);
        return (x[0]-aug5)*par[0]*TMath::Exp(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])+par[4];
        
    }
    
    double backNewK0(double *x, double *par)
    {
        /*if ((x[0] > 1.80) && (x[0] < 1.84)) {
         TF1::RejectPoint();
         return 0;
         }*/
        
        Double_t aug5=1.115683+.497611;//mass os Lambda + K0 in GeV/c^2
        //return (aug5-x[0])*par[0]*TMath::Exp(par[1]+par[2]*aug5+par[3]*aug5*aug5);
        return (x[0]-aug5)*par[0]*TMath::Exp(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])+par[4];
        
        
    }
    
    double functionNewKX(double *x, double *par)
    {
        Double_t aug5=1.115683+0.493677;//mass of Lambda + KX in GeV/c^2
        
        Double_t arg1 = 14.0/22.0; // 2 over pi
        Double_t arg2 = par[6]*par[6]*par[7]*par[7]; //Gamma=par[1]  M=par[2]
        Double_t arg3 = ((x[0]*x[0]) - (par[7]*par[7]))*((x[0]*x[0]) - (par[7]*par[7]));
        Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[6]*par[6])/(par[7]*par[7]));
        //return par[4]*arg1*arg2/(arg3 + arg4)+(aug5-x[0])*par[0]*TMath::Exp(par[1]+par[2]*aug5+par[3]*aug5*aug5);
        return par[5]*arg1*arg2/(arg3 + arg4)+(x[0]-aug5)*par[0]*TMath::Exp(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])+par[4];
        
        
    }
    
    double fNewKX(double *x, double *par)
    {
        if ((x[0] > (1.823-4.0*(0.024/2.35))) && (x[0] < (1.823+4.0*(0.024/2.35)))) {
            TF1::RejectPoint();
            return 0;
        }
        
        Double_t aug5=1.115683+0.493677;//mass of Lambda + KX in GeV/c^2
        //return (aug5-x[0])*par[0]*TMath::Exp(par[1]+par[2]*aug5+par[3]*aug5*aug5);
        return (x[0]-aug5)*par[0]*TMath::Exp(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])+par[4];
        
    }
    
    double backNewKX(double *x, double *par)
    {
        /*if ((x[0] > 1.80) && (x[0] < 1.84)) {
         TF1::RejectPoint();
         return 0;
         }*/
        
        Double_t aug5=1.115683+0.493677;//mass of Lambda + KX in GeV/c^2
        //return (aug5-x[0])*par[0]*TMath::Exp(par[1]+par[2]*aug5+par[3]*aug5*aug5);
        return (x[0]-aug5)*par[0]*TMath::Exp(par[1]+par[2]*x[0]+par[3]*x[0]*x[0])+par[4];
        
    }
    
    double fitMonteCarlo(double *x, double *par)
    {
        return par[0] + par[1]*(x[0]-1.8234);
    }
};

void pinvmass7(int n,char* dc= (char*)"pikx"){
    //gROOT->LoadMacro("~/Desktop/xiTest/overhead.C");
    gStyle->SetOptStat(0);
    if(!n) plot_mass(dc);
    
    return;
}

void plot_mass(char* dc){
    //list of all declarations;
    TCanvas *c0, *c1, *c2, *cA, *cB, *cC, *cD, *cTotalgm, *cTotal0, *cTotal1, *cTotal2, *cTotal3, *cTotalunlike, *cLinearMix, *cLinearGM, *cSubMix, *cSubGM, *cSYSY, *clog;
    TCanvas *cSYSM, *cSYSW, *cSTAM, *cSTAW, *cSYSYGM, *cSYSMGM, *cSYSWGM, *cSTAMGM, *cSTAWGM, *c, *cMass, *cMass2;
    char normal[6];
    char partname[10];
    char partnameComplex[70];
    char othernormal[6];
    char AEword[75];
    int particleloop, hmtflag, writeplots, fitflag, fitsubtract, typemax, flagSB, flag50, newflag, UseProject, UseGMonly, UseSmooth, Use1NR, UseLinear;
    int writetxt, Use1PT, UseNoPoly, Use1Type, UseSigmaRange, UseSigmaRange2, UseVoigt, Fixmean, Nsigmaflag, Nsigmaloop, cutstart, cutend, loopcut, onechange, BellSwitch;
    int Bin1820, rangeloop, backnumber, versionloop, typeloop, removeloop, setloop2;
    double Bin1820Width, fitmin, fitmax, fitmin2, fitmax2;
    double numberofyield, calculateSYSyield, calculateSYSyielderror, calculateSYSyieldBin, calculateSYSyielderrorBin, numberofmean, calculateSYSmean, calculateSYSmeanerror, calculateSTAmeanerror, numberofwidth, calculateSYSwidth, calculateSYSwidtherror, calculateSTAwidtherror;
    TFile *f, *f2, *fevents;
    TH1D *h[4],*a,*aMC,*aEvents;
    TH1D *htemp, *htemp2;
    int typestart, typeend;
    TH1D *AlloftheYield, *AlloftheMean, *AlloftheWidth, *AlloftheChi, *GMAlloftheYield, *GMAlloftheMean, *GMAlloftheWidth, *GMAlloftheChi;
    /*Double_t *function(Double_t , Double_t ), *BWfunction(Double_t , Double_t ), *funVoigt(Double_t , Double_t ), *fline(Double_t , Double_t ), *backline(Double_t , Double_t );
     Double_t *function2nd(Double_t , Double_t ), *f2nd(Double_t , Double_t ), *back2nd(Double_t , Double_t );
     Double_t *function3rd(Double_t , Double_t ), *f3rd(Double_t , Double_t ), *back3rd(Double_t , Double_t );
     Double_t *function4th(Double_t , Double_t ), *f4th(Double_t , Double_t ), *back4th(Double_t , Double_t );
     Double_t *functionNewK0(Double_t , Double_t ), *fNewK0(Double_t , Double_t ), *backNewK0(Double_t , Double_t );
     Double_t *functionNewKX(Double_t , Double_t ), *fNewKX(Double_t , Double_t ), *backNewKX(Double_t , Double_t );*/
    TF1 *fitline, *fitfunction, *fityield, *backfit, *fitConstant, *fitMC, *fitRes;
    TText *tconst;
    TPaveStats *ps1,*psOnePT,*psTemp;
    TAxis *axis;
    double par[9];//8
    int n, m, l, k, i, y, g;
    int polynumber, bmax, bmin;
    double integral, integralbin, binvalue;
    double dy, NumberOfEvents;
    double tempYield, tempMean, tempWidth;
    TF1 *fit1, *simpleline, *simpleline2;
    TMultiGraph *MultiAllYield;
    TMultiGraph *P3MultiAllYield;
    TMultiGraph *P3MultiAllMean;
    TMultiGraph *P3MultiAllWidth;
    TMultiGraph *MultiYieldAverage;
    TMultiGraph *MultiMeanAverage;
    TMultiGraph *MultiWidthAverage;
    TMultiGraph *MultiAE;
    TGraphErrors *gAllYieldP, *gAllMeanP, *gAllWidthP, *OtherData, *AEwitherror;
    TMultiGraph *TestAllWidth;
    ofstream outputFile;
    ofstream outputFileSTAT;
    FILE *Cutoftxt;
    FILE *Percentage;
    AllFunctions * fAll;
    fAll = new AllFunctions();
    int cutloop, hmtloop, cutflag, MCflag, only1flag;
    TLegend *legAE;
    TLegend *legYield;
    TLegend *legMean;
    TLegend *legWidth;
    TLegend *legYieldAverage;
    TLegend *legMeanAverage;
    TLegend *legWidthAverage;
    TLegend *legOnePT;
    TList *fenter;
    //list of all declarations
    
    only1flag=1;//using only1flag will remove the average and lambda+KX mixed event calcualtions
    MCflag=0;
    cutflag=0;//1,2
    hmtflag=1;//if hmtflag==1, activate loop
    writeplots=1;//if writeplots=1, writes mass and sub plots
    //int writeall=0;//if writeall=1, writes allmass plots
    fitflag=1;//if fitflag=1, fits are preformed
    fitsubtract=1;//if fitsubtract=1, fits are subtracted
    typemax=0;
    flagSB=0;
    flag50=1;
    newflag=1;
    UseProject=0;
    UseGMonly=1;
    UseSmooth=0;//if UseSmooth==1, use smooth calculation
    Use1NR=0;//1
    UseLinear=1;
    Nsigmaflag=0; //if Nsigmaflag=1, use Nsigma calculations
    const int versionmax=4;//6+1
    writetxt=1;//if writetxt=1, txt file is written
    Use1PT=1;
    UseNoPoly=1;
    Use1Type=0;//if 1, Only one type is used
    UseSigmaRange=1;//if 1, use sigma as range
    UseSigmaRange2=0;
    UseVoigt=1;//if 1, use voigt
    Fixmean=0;//if 1, fix mean for fit functions
    onechange=1;//if onechange=1, use only the systematic calcualtions that are "one" change from default (free sigma + 3rd, free sigma + 4th etc)
    BellSwitch=0;//if BellSwitch=1, use some of the changes that Dr. Bellwied wanted for some plots.
    
    if (BellSwitch == 1) {//list of caution commands
        if (Use1NR != 1) {
            cout << "ERROR use only 1 Normalization range for BellSwitch" << endl;
            return 0;
        }
    }
    
    Double_t xarray[versionmax], integralarray[2][versionmax], meanarray[2][versionmax], widtharray[2][versionmax];
    Double_t xarrayerror[versionmax], integralarrayerror[2][versionmax], meanarrayerror[2][versionmax], widtharrayerror[2][versionmax];
    Double_t systemxarray[2], systemxarrayerror[2], systemmean[2], systemmeanerror[2], systemwidth[2], systemwidtherror[2], yieldarray[3][versionmax], yieldplotversion[versionmax], yieldplottype[3];
    
    c0=new TCanvas("c0","",1500,1500);
    c1=new TCanvas("c1","",1500,1500);
    c2=new TCanvas("c2","",1500,1500);
    cA=new TCanvas("cA","",1500,1500);
    cB=new TCanvas("cB","",1500,1500);
    cC=new TCanvas("cC","",1500,1500);
    cD=new TCanvas("cD","",1500,1500);
    cTotalgm=new TCanvas("cTotalgm","",2000,2000);
    cTotal0=new TCanvas("cTotal0","",2000,2000);
    cTotal1=new TCanvas("cTotal1","",2000,2000);
    cTotal2=new TCanvas("cTotal2","",2000,2000);
    cTotal3=new TCanvas("cTotal3","",2000,2000);
    cTotalunlike=new TCanvas("cTotalunlike","",2000,2000);
    cLinearMix=new TCanvas("cLinearMix","",2000,2000);
    cLinearGM=new TCanvas("cLinearGM","",2000,2000);
    cSubMix=new TCanvas("cSubMix","",2000,2000);
    cSubGM=new TCanvas("cSubGM","",2000,2000);
    cSYSY=new TCanvas("cSYSY","",2000,2000);
    cSYSM=new TCanvas("cSYSM","",2000,2000);
    cSYSW=new TCanvas("cSYSW","",2000,2000);
    cSTAM=new TCanvas("cSTAM","",2000,2000);
    cSTAW=new TCanvas("cSTAW","",2000,2000);
    cSYSYGM=new TCanvas("cSYSYGM","",2000,2000);
    cSYSMGM=new TCanvas("cSYSMGM","",2000,2000);
    cSYSWGM=new TCanvas("cSYSWGM","",2000,2000);
    cSTAMGM=new TCanvas("cSTAMGM","",2000,2000);
    cSTAWGM=new TCanvas("cSTAWGM","",2000,2000);
    c=new TCanvas("c","",10,10,500,400);
    clog=new TCanvas("clog","",4000,6000);
    cMass=new TCanvas("cMass","",4000,6000);
    cMass2=new TCanvas("cMass2","",4000,6000);
    c->SetFillColor(0);
    c->SetRightMargin(0.02);
    c->Clear();
    c->Update();
    c0->Clear();
    c0->Update();
    c1->Clear();
    c1->Update();
    c2->Clear();
    c2->Update();
    cA->Clear();
    cA->Update();
    cB->Clear();
    cB->Update();
    cC->Clear();
    cC->Update();
    cD->Clear();
    cD->Update();
    cTotal0->Clear();
    cTotal0->Update();
    cTotal1->Clear();
    cTotal1->Update();
    cTotal2->Clear();
    cTotal2->Update();
    cTotal3->Clear();
    cTotal3->Update();
    cTotalgm->Clear();
    cTotalgm->Update();
    cTotalunlike->Clear();
    cTotalunlike->Update();
    cLinearMix->Clear();
    cLinearMix->Update();
    cLinearGM->Clear();
    cLinearGM->Update();
    cSubMix->Clear();
    cSubMix->Update();
    cSubGM->Clear();
    cSubGM->Update();
    cSYSY->Clear();
    cSYSY->Update();
    cSYSM->Clear();
    cSYSM->Update();
    cSYSW->Clear();
    cSYSW->Update();
    cSTAM->Clear();
    cSTAM->Update();
    cSTAW->Clear();
    cSTAW->Update();
    cSYSYGM->Clear();
    cSYSYGM->Update();
    cSYSMGM->Clear();
    cSYSMGM->Update();
    cSYSWGM->Clear();
    cSYSWGM->Update();
    cSTAMGM->Clear();
    cSTAMGM->Update();
    cSTAWGM->Clear();
    cSTAWGM->Update();
    clog->Clear();
    clog->SetLogy();
    clog->Update();
    cMass->Clear();
    cMass->Update();
    cMass2->Clear();
    cMass2->Update();
    int j,scheme;
    int rb;//re-binning factor//5, but 2 is also viable
    int nrd;//For some plots, I only show one normalization region for mixed events; nrd chooses this region.
    int NRstart;
    rb=5;//re-binning factor//5, but 2 is also viable //do not use 3 (number of bin=400)
    nrd=1;//1 //For some plots, I only show one normalization region for mixed events; nrd chooses this region.
    NRstart=0;
    if (Use1NR == 1) {
        NRstart=2;
    }
    versionloop=0;
    rangeloop=0;
    polynumber=0;
    backnumber=0;
    gStyle->SetStatH(0.3);
    j=GetDC(dc); jDC=j; TitlesDC();
    TMultiGraph *AllIntegral = new TMultiGraph();
    TMultiGraph *AllMean = new TMultiGraph();
    TMultiGraph *AllWidth = new TMultiGraph();
    for (hmtloop = 0; hmtloop <= hmtflag; hmtloop++) {//begin hmt loop
        for (particleloop=0; particleloop <= 1; particleloop++) {//begin particle loop
            for (cutloop = 0; cutloop <= cutflag; cutloop++) {//begin cut loop, 0; cutflag
                if (jDC == Klpi) {
                    particleloop=2;
                }
                
                if (particleloop == 0) {
                    jDC=Klkx;
                    if (Nsigmaflag == 1) {
                        cutflag=2;
                    }
                    else{
                        cutflag=0;
                    }
                    //typemax=3;//3
                }
                else if(particleloop == 1){
                    jDC=Klk0;
                    cutflag=0;
                    //typemax=4;//4
                }
                
                if (particleloop == 0) {
                    sprintf(partname,"LambdaKX");
                    if (BellSwitch == 1) {
                        sprintf(partnameComplex,"#Xi^{**#pm}(1820)#rightarrow#Lambda(#bar{#Lambda})K^{#pm}");
                    }
                    else{
                        sprintf(partnameComplex,"#Xi^{**#mp}(1820)#rightarrow#Lambda(#bar{#Lambda})K^{#mp}");
                    }
                }
                else if(particleloop == 1){
                    sprintf(partname,"LambdaK0");
                    sprintf(partnameComplex,"#Xi^{**0}(1820)#rightarrow#Lambda(#bar{#Lambda})K^{0}_{S}");
                }
                else if(particleloop == 2){
                    sprintf(partname,"LambdaPi");
                }
                
                if (hmtloop == 0) {
                    sprintf(normal,"Normal");
                    sprintf(othernormal,"mb");
                }
                else if(hmtloop == 1){
                    sprintf(normal,"hmt");
                    sprintf(othernormal,"hmt");
                }
                
                //TFile* f=TFile::Open(Form("~/Desktop/xitest/plotsAllV5/%s/%s/histogramsAll.root",normal,partname));
                
                if (cutloop == 0) {
                    f=TFile::Open(Form("~/Desktop/xiTest/plotsAllV9/%s/%s/histogramsAll.root",normal,partname));
                    //f=TFile::Open(Form("~/Desktop/xiTest/plotsAllV9/MC/%s/histogramsAll.root",partname));//MC test
                }
                else{
                    f=TFile::Open(Form("~/Desktop/xiTest/plotsAllV9/%s/%s/cut%i/histogramsAll.root",normal,partname,cutloop));
                    //f=TFile::Open(Form("~/Desktop/xiTest/plotsAllV9/MC/%s/cut%i/histogramsAll.root",partname,cutloop));//MC test
                }
                
                if(!f) return;
                
                if (MCflag == 1) {
                    if (cutloop == 0) {
                        f2=TFile::Open(Form("~/Desktop/xiTest/plotsAllV9/MC/%s/histogramsAll.root",partname));
                    }
                    else{
                        f2=TFile::Open(Form("~/Desktop/xiTest/plotsAllV9/MC/%s/cut%i/histogramsAll.root",partname,cutloop));
                    }
                    if(!f2) return;
                }
                
                if (hmtloop == 0) {
                    fevents=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultsAllV9Normal.root"));
                }
                else{
                    fevents=TFile::Open(Form("~/Desktop/xiTest/AnalysisResultsAllV9hmt.root"));
                }
                
                if(!f) return;
                
                char m0[4][100],m1[4][100];//m0 contains the base names of the histograms, m1 contains the titles that could be added to legends (although I do not currently use them)
                sprintf(m0[0],"unlike");
                sprintf(m0[1],"");//will be filled later with "gm" or "sum", depending on the value of scheme
                sprintf(m0[2],"mixunlike");
                sprintf(m0[3],"mixlike");
                if(jDC==Kpikx){
                    scheme=0;
                    sprintf(m1[0],"#pi^{#pm}K^{#mp}");
                    sprintf(m1[1],"Like-Charge");
                    sprintf(m1[2],"Mixed Events");
                    sprintf(m1[3],"");
                }else if(jDC==Kpik0){
                    scheme=1;
                    sprintf(m1[0],"#pi^{#pm}K^{0}_{S}");
                    sprintf(m1[1],"");
                    sprintf(m1[2],"Mixed Events");
                    sprintf(m1[3],"");
                }else if(jDC==Kkxk0){
                    scheme=1;
                    sprintf(m1[0],"K^{#pm}K^{0}_{S}");
                    sprintf(m1[1],"");
                    sprintf(m1[2],"Mixed Events");
                    sprintf(m1[3],"");
                }else if(jDC==Kpkx){
                    scheme=0;
                    sprintf(m1[0],"p(#bar{p})K^{#mp}");
                    sprintf(m1[1],"Like-Charge");
                    sprintf(m1[2],"Mixed Events");
                    sprintf(m1[3],"");
                }else if(jDC==Kpk0){
                    scheme=1;
                    sprintf(m1[0],"p(#bar{p})K^{0}_{S}");
                    sprintf(m1[1],"");
                    sprintf(m1[2],"Mixed Events");
                    sprintf(m1[3],"");
                }else if(jDC==Klpi){
                    //scheme=2;
                    scheme=0;
                    sprintf(DC0,"Lambdapi");
                    sprintf(DC1,"X#rightarrow#Lambda(#bar{#Lambda})#pi^{#pm}");
                    sprintf(m1[0],"#Lambda(#bar{#Lambda})#pi^{#pm}");
                    sprintf(m1[1],"#Lambda(#bar{#Lambda})#pi^{#mp}");
                    sprintf(m1[2],"#Lambda(#bar{#Lambda})#pi^{#pm} Mixed Events");
                    sprintf(m1[3],"#Lambda(#bar{#Lambda})#pi^{#mp} Mixed Events");
                }else if(jDC==Klkx){
                    //scheme=2;
                    scheme=0;
                    sprintf(DC0,"Lambdakx");
                    if (BellSwitch == 1) {
                        sprintf(DC1,"X#rightarrow#Lambda(#bar{#Lambda})K^{#pm}");
                    }
                    else{
                        sprintf(DC1,"X#rightarrow#Lambda(#bar{#Lambda})K^{#mp}");
                    }
                    sprintf(m1[0],"#Lambda(#bar{#Lambda})K^{#pm}");
                    sprintf(m1[1],"#Lambda(#bar{#Lambda})K^{#mp}");
                    sprintf(m1[2],"#Lambda(#bar{#Lambda})K^{#pm} Mixed Events");
                    sprintf(m1[3],"#Lambda(#bar{#Lambda})K^{#mp} Mixed Events");
                }else if(jDC==Klk0){
                    scheme=1;
                    sprintf(DC0,"Lambdak0");
                    sprintf(DC1,"X#rightarrow#Lambda(#bar{#Lambda})K^{0}_{S}");
                    sprintf(m1[0],"#Lambda(#bar{#Lambda})K^{0}_{S}");
                    sprintf(m1[1],"");
                    sprintf(m1[2],"Mixed Events");
                    sprintf(m1[3],"");
                }else if(jDC==Klp){
                    scheme=2;
                    sprintf(m1[0],"#Lambdap(#bar{#Lambda}#bar{p})");
                    sprintf(m1[1],"#Lambda#bar{p}(#bar{#Lambda}p)");
                    sprintf(m1[2],"#Lambdap(#bar{#Lambda}#bar{p}) Mixed Events");
                    sprintf(m1[3],"#Lambda#bar{p}(#bar{#Lambda}p) Mixed Events");
                }else return;
                
                if(!scheme) sprintf(m0[1],"gm");
                else if(scheme==2) sprintf(m0[1],"sum");
                
                ptbins();
                //USe 1 section
                //if (Use1PT == 1) {
                //const int MinPT=0;//3
                //if (Use1PT == 1) {
                const int MinPT=3;//3
                //}
                const int MaxPT=npb;
                
                const int MaxCen=ncb;
                /*}
                 else{
                 const int MinPT=0;
                 const int MaxPT=npb;
                 }
                 //use 1 section
                 if (Use1Type == 1) {
                 const int array8=1;
                 const int array4=1;
                 const int array1=1;
                 }
                 else{*/
                
                /*const int array8=8;
                 const int array4=4;
                 const int array1=1;*/
                //if (cutloop == 0) {
                const int array8=6;//2
                const int array4=6;//2
                const int array1=6;//2
                const int NRBegin=NRstart;
                const int NREnd=nNR;
                //if (Use1NR == 1) {
                    //const int NREnd=nNR-1;
                //}
                //}
                Double_t ALLYield[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array8];//[cen][pt][NR][range][version][type]
                Double_t BinALLYield[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array8];//[cen][pt][NR][range][version][type]
                Double_t ALLMean[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array4];//[pt][NR][range][version][type]
                Double_t ALLWidth[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array4];//[pt][NR][range][version][type]
                Double_t ALLChi2[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array1];//[pt][NR][range][version][type]
                Double_t GMALLYield[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array8];//[pt][NR][range][version][type]
                Double_t BinGMALLYield[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array8];//[pt][NR][range][version][type]
                Double_t GMALLMean[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array4];//[pt][NR][range][version][type]
                Double_t GMALLWidth[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array4];//[pt][NR][range][version][type]
                Double_t GMALLChi2[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array1];//[pt][NR][range][version][type]
                Double_t ALLYieldError[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array8];
                Double_t ALLMeanError[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array4];
                Double_t ALLWidthError[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array4];
                Double_t GMALLYieldError[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array8];
                Double_t GMALLMeanError[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array4];
                Double_t GMALLWidthError[MaxCen][MaxPT-MinPT][NREnd-NRBegin][3][versionmax][array4];
                Double_t AEarray[MaxCen][MaxPT-MinPT];
                Double_t AEarrayerror[MaxCen][MaxPT-MinPT];
                Double_t AEslope[MaxCen][MaxPT-MinPT];
                Double_t AEslopeerror[MaxCen][MaxPT-MinPT];
                Double_t AEXarray[MaxPT-MinPT];
                Double_t AEXarrayerror[MaxPT-MinPT];
                Double_t AEYarray[MaxPT-MinPT];
                Double_t AEYarrayerror[MaxPT-MinPT];
                Double_t MCSigma[MaxCen][MaxPT-MinPT];
                TMultiGraph *MultiAE = new TMultiGraph();
                
                //}
                
                for(cb=0;cb<ncb && ptbins();cb++) for(pb=MinPT;pb<npb;pb++){//loop over multiplicity and pT bins
                    //first loop to look at
                    if (Use1Type != 1) {
                        AlloftheYield = new TH1D(Form("yield_%s_%s_p%i",partname,othernormal,pb),"",8,0,8);
                        AlloftheMean = new TH1D(Form("mean_%s_%s_p%i",partname,othernormal,pb),"",8,0,8);
                        AlloftheMean->SetMinimum(1.800);
                        AlloftheWidth = new TH1D(Form("width_%s_%s_p%i",partname,othernormal,pb),"",8,0,8);
                        AlloftheChi = new TH1D(Form("#Chi^{2}_%s_%s_p%i",partname,othernormal,pb),"",8,0,8);
                        GMAlloftheYield = new TH1D(Form("yield_%s_%s_p%i_GM",partname,othernormal,pb),"",8,0,8);
                        GMAlloftheMean = new TH1D(Form("mean_%s_%s_p%i_GM",partname,othernormal,pb),"",8,0,8);
                        GMAlloftheMean->SetMinimum(1.800);
                        GMAlloftheWidth = new TH1D(Form("width_%s_%s_p%i_GM",partname,othernormal,pb),"",8,0,8);
                        GMAlloftheChi = new TH1D(Form("#Chi^{2}_%s_%s_p%i_GM",partname,othernormal,pb),"",8,0,8);
                    }
                    //possible else
                    
                    if ((hmtloop == 1) && (cb >= 1)) {
                        cout << "one centrality bin for hmt" << endl;
                        continue;
                    }
                    
                    //plot distributions before background subtraction
                    /*for (int pb2=0; pb2 < npb; pb2++) {//terporary loop
                     for(j=0;j<4;j++){
                     h[j]=0;
                     if(scheme==1 && j==1) continue;
                     if(scheme<2 && j==3) continue;
                     
                     //set colors for histograms
                     if(!j) AKColor("black",color);
                     else if(j<3) TitlesBS(j);
                     else AKColor("orange",color);
                     
                     //get histograms
                     if(j<=1) sprintf(s0,"mass_%s_c%i_p%i",m0[j],cb,pb2);//unlike- and like-charge histograms
                     else sprintf(s0,"mass_%sN%s%i_c%i_p%i",m0[j-2],m0[j],nrd,cb,pb2);//mixed-event histograms (normalized, only nrd shown)
                     h[j]=(TH1D*) f->Get(s0);
                     if(!h[j]){cerr<<"missing histogram "<<s0<<endl; return;}
                     cout << "name " << s0 << " j " << j << endl;
                     //if (rangeloop == 0) {
                     h[j]->Rebin(rb);
                     //}
                     h[j]->SetMarkerStyle(1);
                     h[j]->SetMarkerColor(TColor::GetColor(color));
                     h[j]->SetLineColor(TColor::GetColor(color));
                     
                     }
                     //plot histograms
                     h[0]->SetTitle(Form("%s, %1.1f<Centrality<%1.1f, %1.1f<#it{p}_{T}<%1.1f GeV/c",DC1,cbmin[cb],cbmax[cb],ptmin[pb2],ptmax[pb2]));
                     h[0]->SetXTitle("Invariant Mass [GeV/#it{c^{2}}]");
                     h[0]->GetXaxis()->SetTitleOffset(1.3);
                     //h[0]->GetXaxis()->SetRange(0,50);
                     h[0]->GetXaxis()->SetRange(h[0]->GetXaxis()->FindBin(1.60),h[0]->GetXaxis()->FindBin(2.10));
                     h[0]->SetYTitle("Counts/(0.01 GeV/c^{2})");
                     h[0]->GetYaxis()->SetTitleOffset(1.3);
                     c->cd();
                     h[0]->Draw();
                     for(j=1;j<4;j++) {
                     if (j != 2) {
                     if(h[j]) {
                     h[j]->Draw("same");
                     }
                     }
                     }
                     c->Print(Form("plotsCompleteV3/%s/mass_c%i_p%i.pdf",partname,cb,pb2));
                     }//end of temproary loop
                     */
                    // to be reactivated
                    for(j=0;j<4;j++){
                        h[j]=0;
                        if(scheme==1 && j==1) continue;
                        if(scheme<2 && j==3) continue;
                        
                        //set colors for histograms
                        if(!j) AKColor((char *)"black",color);
                        else if(j<3) TitlesBS2(j);
                        else AKColor((char *)"orange",color);
                        
                        //get histograms
                        //if(j<=1) sprintf(s0,"mass_%s_c%i_p%i",m0[j],cb,pb);//unlike- and like-charge histograms
                        if(j == 0) sprintf(s0,"mass_%s_c%i_p%i",m0[j],cb,pb);//unlike- charge histograms
                        else if(j == 1) sprintf(s0,"mass_unlikeN%s%i_c%i_p%i",m0[j],nrd,cb,pb);//gm histograms
                        else sprintf(s0,"mass_%sN%s%i_c%i_p%i",m0[j-2],m0[j],nrd,cb,pb);//mixed-event histograms (normalized, only nrd shown)
                        //cout << "hmt " << hmtloop << " cutloop " << cutloop << " c " << cb << " pb " << pb << endl;
                        h[j]=(TH1D*) f->Get(s0);
                        if(!h[j]){cerr<<"missing histogram "<<s0<<endl; return;}
                        cout << "name " << s0 << " j " << j << endl;
                        //if (rangeloop == 0) {
                        h[j]->Rebin(rb);
                        //}
                        h[j]->SetMarkerStyle(1);
                        h[j]->SetMarkerColor(TColor::GetColor(color));
                        h[j]->SetLineColor(TColor::GetColor(color));
                        
                    }
                    //plot histograms
                    h[0]->SetTitle(Form("%s, %1.1f-%1.1f %% , %1.1f<#it{p}_{T}<%1.1f GeV/c",DC1,cbmin[cb],cbmax[cb],ptmin[pb],ptmax[pb]));
                    h[0]->SetXTitle("Invariant Mass [GeV/#it{c^{2}}]");
                    h[0]->GetXaxis()->SetTitleOffset(1.3);
                    //h[0]->GetXaxis()->SetRange(0,50);
                    h[0]->GetXaxis()->SetRange(h[0]->GetXaxis()->FindBin(1.70),h[0]->GetXaxis()->FindBin(2.10));
                    h[0]->SetYTitle("Counts/(0.01 GeV/c^{2})");
                    h[0]->GetYaxis()->SetTitleOffset(1.3);
                    c->cd();
                    h[0]->Draw();
                    for(j=1;j<4;j++) if(h[j]) h[j]->Draw("same");
                    //for(j=1;j<4;j++) if((h[j]) && ((j != 2) && (particleloop == 0))) h[j]->Draw("same");
                    
                    if (hmtloop == 0) {//Normal
                        if (cutloop == 0) {//no cut
                            c->Print(Form("plotsCompleteV3/%s/mass_c%i_p%i.pdf",partname,cb,pb));
                        }//no cut
                        else{//cut
                            c->Print(Form("plotsCompleteV3/%s/cut%i_mass_c%i_p%i.pdf",partname,cutloop,cb,pb));
                        }//cut
                    }//Normal
                    else{//hmt
                        if (cutloop == 0) {//no cut
                            c->Print(Form("plotsCompleteV3/%s/hmt_mass_c%i_p%i.pdf",partname,cb,pb));
                        }//no cut
                        else{//cut
                            c->Print(Form("plotsCompleteV3/%s/hmtcut%i_mass_c%i_p%i.pdf",partname,cutloop,cb,pb));
                        }//cut
                    }//hmt
                    
                    //MC level
                    if (MCflag == 1) {
                        //resolution
                        sprintf(s0,"pt_sum_res_c%i_p%i",cb,pb);
                        aMC=(TH1D*) f2->Get(s0);
                        aMC->SetTitle(Form("mass resolution for %1.1f-%1.1f %%  %1.1f<pT<%1.1f",cbmin[cb],cbmax[cb],ptmin[pb],ptmax[pb]));
                        aMC->GetXaxis()->SetTitle("Resolution");
                        //fitRes = new TF1("fitRes",fAll,&AllFunctions::plainVoigt,-0.02,0.02,4);
                        //fitRes->FixParameter(0,aMC->GetBinContent(100));
                        //fitRes->SetParameter(1,0);
                        //fitRes->SetParameter(2,0);
                        //fitRes->SetParameter(3,0);
                        //aMC->Fit("fitRes","R","SAME",-0.02,0.02);
                        //aMC->Fit("gaus","B");
                        //fitRes = aMC->GetFunction("gaus");
                        
                        /*
                        fitRes = new TF1("fitRes","gaus");
                        fitRes->FixParameter(0,aMC->GetBinContent(100));
                        aMC->Fit("fitRes");
                        MCSigma[cb][pb-MinPT]=fitRes->GetParameter(2);
                        */
                        
                        //MCSigma[cb][pb-MinPT]=fitRes->GetParameter(3);
                        if (hmtloop == 0) {
                            if ((cb == 0) && (pb == MinPT)) {
                                cMass2->Clear();
                                if (Use1PT == 1) {
                                    //cMass2->Divide(2,2);//to be reinstated
                                    cMass2->Divide(MaxCen/2.0,MaxCen/2.0);
                                }
                                else{
                                    //cMass2->Divide(4,(MaxPT-MinPT));//to be reinstated
                                    cMass2->Divide(MaxCen,(MaxPT-MinPT));
                                }
                                cMass2->Update();
                            }
                            cMass2->cd((cb+1)+(MaxPT-MinPT)*(pb-MinPT));
                        }
                        else{
                            if (pb == MinPT) {
                                cMass2->Clear();
                                if (Use1PT == 1) {
                                    cMass2->Divide(1);
                                }
                                else{
                                    cMass2->Divide((MaxPT-MinPT));
                                }
                                cMass2->Update();
                            }
                            cMass2->cd(pb-MinPT+1);
                        }
                        aMC->DrawCopy();
                        cMass2->Update();
                        if (hmtloop == 0) {//Normal
                            if (cutloop == 0) {//no cut
                                if ((cb == (MaxCen-1)) && (pb == (MaxPT-1))) {
                                    cMass2->Print(Form("plotsCompleteV3/%s/Res.pdf",partname));
                                }
                            }//no cut
                            else{//cut
                                if ((cb == (MaxCen-1)) && (pb == (MaxPT-1))) {
                                    cMass2->Print(Form("plotsCompleteV3/%s/cut%i_Res.pdf",partname,cutloop));
                                }
                            }//cut
                        }//Normal
                        else{//hmt
                            if (cutloop == 0) {//no cut
                                if (pb == (MaxPT-1)) {
                                    cMass2->Print(Form("plotsCompleteV3/%s/hmt_Res.pdf",partname));
                                }
                            }//no cut
                            else{//cut
                                if (pb == (MaxPT-1)) {
                                    cMass2->Print(Form("plotsCompleteV3/%s/hmtcut%i_Res.pdf",partname,cutloop));
                                }
                            }//cut
                        }//hmt
                        //fitRes->ReleaseParameter(0);
                        //resolution
                        
                        //sprintf(s0,"pt_sum_trueMM_gen_c%i_p%i",cb,pb);
                        sprintf(s0,"pt_sum_trueMM_gen_c%i_m%i",cb,pb);
                        aMC=(TH1D*) f2->Get(s0);
                        aMC->SetTitle(Form("Acceptance * Efficiency for %1.1f-%1.1f %%  %1.1f<pT<%1.1f",cbmin[cb],cbmax[cb],ptmin[pb],ptmax[pb]));
                        aMC->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
                        aMC->GetYaxis()->SetTitle("A*E (trueMM/gen)");
                        //rebin
                        /*
                        if (cutloop == 0) {
                            aMC->Rebin(20);
                        }
                         */
                        /*
                        fitConstant = new TF1("fitConstant","pol0",1.720,1.920);
                        aMC->Fit("fitConstant","Q","SAME",1.720,1.920);
                        */
                        par[0]=0;
                        par[1]=0;
                        fitMC = new TF1("fitMC",fAll,&AllFunctions::fitMonteCarlo,1.720,1.920,2);
                        aMC->Fit("fitMC","Q","SAME",1.720,1.920);
                        if (hmtloop == 0) {
                            if ((cb == 0) && (pb == MinPT)) {
                                cMass->Clear();
                                if (Use1PT == 1) {
                                    cMass->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                }
                                else{
                                    cMass->Divide(MaxCen,(MaxPT-MinPT));//4,()
                                }
                                cMass->Update();
                            }
                            cMass->cd((cb+1)+(MaxPT-MinPT)*(pb-MinPT));
                        }
                        else{
                            if (pb == MinPT) {
                                cMass->Clear();
                                if (Use1PT == 1) {
                                    cMass->Divide(1);
                                }
                                else{
                                    cMass->Divide((MaxPT-MinPT));
                                }
                                cMass->Update();
                            }
                            cMass->cd(pb-MinPT+1);
                        }
                        aMC->DrawCopy();
                        //fitConstant->DrawCopy("SAME");
                        fitMC->DrawCopy("SAME");
                        cMass->Update();
                        if (hmtloop == 0) {//Normal
                            if (cutloop == 0) {//no cut
                                if ((cb == (MaxCen-1)) && (pb == (MaxPT-1))) {
                                    cMass->Print(Form("plotsCompleteV3/%s/A*E.pdf",partname));
                                }
                            }//no cut
                            else{//cut
                                if ((cb == (MaxCen-1)) && (pb == (MaxPT-1))) {
                                    cMass->Print(Form("plotsCompleteV3/%s/cut%i_A*E.pdf",partname,cutloop));
                                }
                            }//cut
                        }//Normal
                        else{//hmt
                            if (cutloop == 0) {//no cut
                                if (pb == (MaxPT-1)) {
                                    cMass->Print(Form("plotsCompleteV3/%s/hmt_A*E.pdf",partname));
                                }
                            }//no cut
                            else{//cut
                                if (pb == (MaxPT-1)) {
                                    cMass->Print(Form("plotsCompleteV3/%s/hmtcut%i_A*E.pdf",partname,cutloop));
                                }
                            }//cut
                        }//hmt
                        c->cd();
                        c->Clear();
                        c->Update();
                        //to change
                        /*
                        AEarray[cb][pb-MinPT]=(fitConstant->GetParameter(0))/20.0;
                        AEarrayerror[cb][pb-MinPT]=(fitConstant->GetParError(0))/20.0;
                         */
                        //AEarray[cb][pb-MinPT]=(fitMC->GetParameter(0))/20.0;
                        //AEarrayerror[cb][pb-MinPT]=(fitMC->GetParError(0))/20.0;
                        //AEslope[cb][pb-MinPT]=(fitMC->GetParameter(1))/20.0;
                        //AEslopeerror[cb][pb-MinPT]=(fitMC->GetParError(1))/20.0;
                        if (pb == (npb-1)) {
                            for (i=MinPT; i<npb; i++) {
                                AEXarray[i-MinPT]=(ptmin[i-MinPT]+ptmax[i-MinPT])/2.0;
                                AEXarrayerror[i-MinPT]=(ptmax[i-MinPT]-ptmin[i-MinPT])/2.0;
                                AEYarray[i-MinPT]=AEarray[cb][i-MinPT];
                                AEYarrayerror[i-MinPT]=AEarrayerror[cb][i-MinPT];
                            }
                            AEwitherror = new TGraphErrors(MaxPT-MinPT,AEXarray,AEYarray,AEXarrayerror,AEYarrayerror);
                            AEwitherror->SetLineColor(1+cb);
                            MultiAE->Add(AEwitherror);
                            if (cb == 0) {
                                legAE= new TLegend(0.7,0.7,0.9,0.9);
                            }
                            sprintf(AEword,"%s %1.1f<cen<%1.1f",partname,cbmin[cb],cbmax[cb]);
                            legAE->AddEntry(AEwitherror,AEword,"l");
                            if (cb == (MaxCen-1)) {
                                if (hmtloop == 0) {//Normal
                                    if (cutloop == 0) {//no cut
                                        sprintf(AEword,"PT-specrta of %s",partname);
                                    }//no cut
                                    else{//cut
                                        sprintf(AEword,"PT-specrta of %s cut%i",partname,cutloop);
                                    }//cut
                                }//Normal
                                else{//hmt
                                    if (cutloop == 0) {//no cut
                                        sprintf(AEword,"PT-specrta of %s hmt",partname);
                                    }//no cut
                                    else{//cut
                                        sprintf(AEword,"PT-specrta of %s hmt-cut%i",partname,cutloop);
                                    }//cut
                                }//hmt
                                MultiAE->SetTitle(AEword);
                                MultiAE->DrawClone("apl");
                                legAE->Draw("same");
                                c->Update();
                                if (hmtloop == 0) {//Normal
                                    if (cutloop == 0) {//no cut
                                        c->Print(Form("plotsCompleteV3/%s/A*E_complete.pdf",partname));
                                    }//no cut
                                    else{//cut
                                        c->Print(Form("plotsCompleteV3/%s/cut%i_A*E_complete.pdf",partname,cutloop));
                                    }//cut
                                }//Normal
                                else{//hmt
                                    if (cutloop == 0) {//no cut
                                        c->Print(Form("plotsCompleteV3/%s/hmt_A*E_complete.pdf",partname));
                                    }//no cut
                                    else{//cut
                                        c->Print(Form("plotsCompleteV3/%s/hmtcut%i_A*E_complete.pdf",partname,cutloop));
                                    }//cut
                                }//hmt
                            }
                        }
                    }
                    //MC level
                    
                    if (hmtloop == 0) {//992
                        sprintf(s0,"RsnOut_%s",DC0);//get the TList object corresponding to the desired decay channel
                        /*if (particleloop == 0) {//for pp
                         sprintf(s0,"RsnOut_%s_%s",DC0,DChmt);//get the TList object corresponding to the desired decay channel
                         }//for pp*/
                        
                        if (cutloop >= 1) {
                            sprintf(s0,"RsnOut_%s_cut%i",DC0,cutloop);//get the TList object corresponding to the desired decay channel
                        }
                    }
                    else{//993
                        sprintf(s0,"RsnOut_%s_hmt",DC0);//get the TList object corresponding to the desired decay channel
                        if (cutloop >= 1) {
                            sprintf(s0,"RsnOut_%s_hmt_cut%i",DC0,cutloop);//get the TList object corresponding to the desired decay channel
                        }
                    }
                    fenter=(TList*) fevents->Get(s0);
                    sprintf(s0,"hAEventsVsMulti");
                    aEvents=(TH1D*) fenter->FindObject(s0);
                    NumberOfEvents=aEvents->GetEntries();
                    dy=1.0;
                    
                    for (rangeloop=0; rangeloop <= 2; rangeloop++) {//begin range loop//0,2
                        cout << "rangeloop = " << rangeloop << endl;

                        //int rangeloop=0;
                        fitmin=0.0;
                        fitmax=0.0;
                        fitmin2=0.0;
                        fitmax2=0.0;
                        
                        if ((BellSwitch == 1) && (rangeloop > 0)) {
                            cout << "Bellswitch activate for rangeloop " << rangeloop << endl;
                            continue;
                        }
                        
                        /*if (rangeloop == 1) {
                         fitmin=1.7;
                         fitmax=2.0;
                         }
                         else if(rangeloop == 2){
                         fitmin=1.7;
                         fitmax=1.95;
                         }
                         else if(rangeloop == 3){
                         fitmin=1.73;
                         fitmax=2.0;
                         }
                         else*/
                        /*if(rangeloop == 0){
                         fitmin=1.73;
                         fitmax=1.97;
                         }*/
                        
                        //----- plot background-subtracted distributions
                        backnumber=0;
                        if (Use1Type == 1) {
                            typestart=8;
                            typeend=8;
                        }
                        else{
                            /*int typestart=1;
                             int typeend=8;*/
                            typestart=1;
                            typeend=2;
                        }
                        
                        if(!scheme){//plot unlike-gm
                            ///jNR=2;//set to best 2
                            //Remember that Normalization range 2 is best for unlike-gm for LambdaKx
                            
                            ////rangeloop=1;//set to best
                            ////versionloop=2;//set to best
                            
                            /*if (rangeloop == 1) {
                             fitmin=1.7;
                             fitmax=2.0;
                             }
                             else if(rangeloop == 2){
                             fitmin=1.7;
                             fitmax=1.95;
                             }
                             else if(rangeloop == 3){
                             fitmin=1.73;
                             fitmax=2.0;
                             }
                             else if(rangeloop == 4){
                             fitmin=1.73;
                             fitmax=1.95;
                             }*/
                            
                            if (UseSigmaRange == 1) {//use sigma range
                                //fitmin=(1.823-6.0*(0.024/2.35));//6.0
                                //fitmax=(1.823+6.0*(0.024/2.35));//6.0
                                
                                fitmin2=(1.823-4.0*(0.024/2.35));//4.0
                                fitmax2=(1.823+4.0*(0.024/2.35));//4.0
                                if (rangeloop == 0) {
                                    fitmin=(1.7);//changed from 1.760 to 1.7
                                    fitmax=(2.0);
                                }
                                else if(rangeloop == 1){
                                    fitmin=(1.760);
                                    fitmax=(1.9);
                                }
                                else if(rangeloop == 2){
                                    fitmin=(1.760);//changed from 1.7 to 1.760
                                    fitmax=(2.0);
                                }
                            }//use sigma range
                            
                            for(jNR=NRstart;jNR<nNR;jNR++){//plot unlike-gm
                                
                                TitlesBS2(1);
                                sprintf(s0,"mass_unlikeQgm%i_c%i_p%i",jNR,cb,pb);
                                a=(TH1D*) f->Get(s0);
                                for (versionloop=0; versionloop <= (versionmax-1); versionloop ++) {//begin version loop
                                    cout << "version loop " << versionloop << endl;
                                    for (typeloop=typestart; typeloop <= typeend; typeloop++) {//begin typeloop for unlike
                                        //if (versionloop == 0) {
                                        if (Use1Type == 1) {
                                            if (typeloop == 8) {
                                                a->Rebin(rb);
                                                //a->GetXaxis()->SetRange(0,50);
                                                a->GetXaxis()->SetRange(a->GetXaxis()->FindBin(1.70),a->GetXaxis()->FindBin(2.10));
                                                Bin1820 = a->GetXaxis()->FindBin(1.82);
                                                Bin1820Width = a->GetXaxis()->GetBinWidth(Bin1820);
                                            }
                                        }
                                        else{
                                            if (typeloop == 1) {
                                                if (rangeloop == 0) {
                                                    if (versionloop == 0) {
                                                        //if (jNR == NRstart) {
                                                        a->Rebin(rb);
                                                        a->GetXaxis()->SetRange(a->GetXaxis()->FindBin(1.70),a->GetXaxis()->FindBin(2.10));
                                                        Bin1820 = a->GetXaxis()->FindBin(1.82);
                                                        Bin1820Width = a->GetXaxis()->GetBinWidth(Bin1820);
                                                        //}
                                                    }
                                                }
                                            }
                                        }
                                        a->SetTitle(Form("%s, %1.1f-%1.1f %% , %1.1f<#it{p}_{T}<%1.1f GeV/c",DC1,cbmin[cb],cbmax[cb],ptmin[pb],ptmax[pb]));
                                        a->SetXTitle("Invariant Mass [GeV/#it{c^{2}}]");
                                        a->GetXaxis()->SetTitleOffset(1.3);
                                        a->SetYTitle("Counts/(0.01 GeV/c^{2})");
                                        a->GetYaxis()->SetTitleOffset(1.3);
                                        a->SetMarkerStyle(1);
                                        a->SetMarkerColor(TColor::GetColor(color));
                                        a->SetLineColor(TColor::GetColor(color));
                                        //}
                                        if (Use1Type != 1) {
                                            GMAlloftheYield->GetXaxis()->SetBinLabel((typeloop),Form("peak%i",typeloop));
                                            GMAlloftheMean->GetXaxis()->SetBinLabel((typeloop),Form("peak%i",typeloop));
                                            GMAlloftheWidth->GetXaxis()->SetBinLabel((typeloop),Form("peak%i",typeloop));
                                            xarray[0]=0;
                                            xarrayerror[0]=0;
                                            for (i=0; i <= 1; i++) {
                                                systemxarray[i]=i;
                                                systemxarrayerror[i]=0;
                                            }
                                        }
                                        //possible else
                                        
                                        GMALLYield[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                        GMALLMean[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                        GMALLMeanError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                        GMALLWidth[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                        GMALLWidthError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                        GMALLChi2[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                        backnumber=0;
                                        if (fitsubtract == 1) {
                                            htemp = (TH1D*)a->Clone();//new section
                                            htemp2 = (TH1D*)a->Clone();
                                        }
                                        
                                        if (fitflag == 1) {//fitflag
                                            //double par[8]=0;
                                            if (versionloop == 0) {
                                                polynumber=5;
                                                backnumber=2;
                                            }
                                            else if(versionloop == 1){
                                                polynumber=6;
                                                backnumber=3;
                                            }
                                            else if(versionloop == 2){
                                                polynumber=7;
                                                backnumber=4;
                                            }
                                            else if(versionloop == 3){
                                                polynumber=8;
                                                backnumber=5;
                                            }
                                            
                                            for (int z = 0; z < (polynumber+1); z++) {
                                                par[z]=0;
                                            }
                                            
                                            if (versionloop == 0) {
                                                //fitfunction = new TF1("fitfunction",function,fitmin,fitmax,polynumber);
                                                fitfunction = new TF1("fitfunction",fAll,&AllFunctions::function,fitmin,fitmax,polynumber,"AllFunctions","function");
                                                fitline = new TF1("fitline",fAll,&AllFunctions::fline,fitmin,fitmax,backnumber,"AllFunctions","fline");//new
                                                //fitline = new TF1("fitline",fline,fitmin,fitmax,backnumber);//new
                                            }
                                            else if(versionloop == 1){
                                                fitfunction = new TF1("fitfunction",fAll,&AllFunctions::function2nd,fitmin,fitmax,polynumber);
                                                fitline = new TF1("fitline",fAll,&AllFunctions::f2nd,fitmin,fitmax,backnumber);//new
                                            }
                                            else if(versionloop == 2){
                                                fitfunction = new TF1("fitfunction",fAll,&AllFunctions::function3rd,fitmin,fitmax,polynumber);
                                                fitline = new TF1("fitline",fAll,&AllFunctions::f3rd,fitmin,fitmax,backnumber);//new
                                            }
                                            else if(versionloop == 3){
                                                fitfunction = new TF1("fitfunction",fAll,&AllFunctions::function4th,fitmin,fitmax,polynumber);
                                                fitline = new TF1("fitline",fAll,&AllFunctions::f4th,fitmin,fitmax,backnumber);//new
                                            }
                                            
                                            fitfunction->SetParLimits(backnumber, 0.0, 500000.0);
                                            fitfunction->ReleaseParameter(backnumber+1);
                                            fitfunction->ReleaseParameter(backnumber+2);
                                            fitfunction->SetParLimits(backnumber+1,0.01,0.06);
                                            fitfunction->SetParLimits(backnumber+2,1.8,1.85);
                                            fitfunction->SetParameter(backnumber+1,0.025);
                                            fitfunction->SetParameter(backnumber+2,1.82);
                                            if (Fixmean == 1) {
                                                fitfunction->FixParameter(backnumber+2,1.8218);
                                            }
                                            fitfunction->SetNpx(500);
                                            fitline->SetNpx(500);
                                            a->Fit("fitline","Q","SAME",fitmin,fitmax);
                                            
                                            if ((fitline->GetChisquare()/fitline->GetNDF()) >= 10) {
                                                cout << "***********************ERROR in FITLine****************************************" << endl;
                                                //cout << "data is type " << typeloop << " version " << versionloop << " range " << rangeloop << endl;
                                                a->Fit("fitline","Q0","",fitmin,fitmax);
                                                //////Possible check for bad fits
                                                if ((fitline->GetChisquare()/fitline->GetNDF()) >= 10) {
                                                    //cout << "refit failed, reject" << endl;
                                                    //continue;
                                                }
                                                else{
                                                    cout << "success in refiting" << endl;
                                                }
                                            }
                                            
                                            for (int j=0; j < backnumber; j++) {
                                                if (typeloop <= 7) {
                                                    fitfunction->FixParameter(j,fitline->GetParameter(j));
                                                }
                                                else{
                                                    fitfunction->SetParameter(j,fitline->GetParameter(j));
                                                }
                                            }
                                            
                                            if (Fixmean == 1) {//possible error in code
                                                for (int limiter=0; limiter < backnumber; limiter++) {
                                                    fitfunction->SetParLimits(limiter,(fitline->GetParameter(limiter))/10.0,(fitline->GetParameter(limiter))*10.0);
                                                    cout << "parameter " << limiter << " limits " << (fitline->GetParameter(limiter))/10.0 << " " << (fitline->GetParameter(limiter))*10.0 << endl;
                                                }
                                                //fitfunction->SetParLimits(0,(fitline->GetParameter(0))/10.0,(fitline->GetParameter(0))*10.0);
                                                //cout << "parameter " << 0 << " limits " << (fitline->GetParameter(0))/10.0 << " " << (fitline->GetParameter(0))*10.0 << endl;
                                            }
                                            
                                            a->Fit("fitfunction","QB","SAME",fitmin,fitmax);
                                            if ((fitfunction->GetChisquare()/fitfunction->GetNDF()) >= 10) {
                                                for (int refitloop=0; refitloop <= 10; refitloop++) {
                                                    if ((fitfunction->GetChisquare()/fitfunction->GetNDF()) >= 10) {
                                                        cout << "Begin refit" << endl;
                                                        a->Fit("fitfunction","Q","SAME",fitmin,fitmax);
                                                    }
                                                }
                                                if ((fitfunction->GetChisquare()/fitfunction->GetNDF()) >= 10) {
                                                    cout << "ERROR With Refit procedure ****************************************************" << endl;
                                                    //cout << "reject data" << endl;
                                                    //continue;
                                                }
                                                ////possible error check here
                                            }
                                            
                                            if (typeloop == typestart) {
                                                cTotalgm->Clear();
                                            }
                                            cTotalgm->cd();
                                            gPad->SetRightMargin(0.01);
                                            if (typeloop == typestart) {
                                                a->Draw();
                                                cTotalgm->Update();
                                            }
                                            fitfunction->SetLineColor(typeloop);
                                            fitfunction->SetLineStyle(1);
                                            fitfunction->SetLineWidth(1);
                                            fitfunction->DrawCopy("SAME");
                                            cTotalgm->Update();
                                            
                                            GMALLMean[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fitfunction->GetParameter(backnumber+2);
                                            GMALLMeanError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fitfunction->GetParError(backnumber+2);
                                            GMALLWidth[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fitfunction->GetParameter(backnumber+1);
                                            GMALLWidthError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fitfunction->GetParError(backnumber+1);
                                            GMALLChi2[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=(fitfunction->GetChisquare()/fitfunction->GetNDF());
                                            
                                            if (UseNoPoly != 1) {//old subtraction method
                                                if (fitsubtract == 1) {
                                                    if (versionloop == 1) {
                                                        backfit = new TF1("backfit",fAll,&AllFunctions::backline,fitmin,fitmax,backnumber);
                                                        fityield = new TF1("fityield",fAll,&AllFunctions::function,1.6,2.4,polynumber);
                                                    }
                                                    else if(versionloop == 2){
                                                        backfit = new TF1("backfit",fAll,&AllFunctions::back2nd,fitmin,fitmax,backnumber);
                                                        fityield = new TF1("fityield",fAll,&AllFunctions::function2nd,1.6,2.4,polynumber);
                                                        
                                                    }
                                                    else if(versionloop == 3){
                                                        backfit = new TF1("backfit",fAll,&AllFunctions::back3rd,fitmin,fitmax,backnumber);
                                                        fityield = new TF1("fityield",fAll,&AllFunctions::function3rd,1.6,2.4,polynumber);
                                                    }
                                                    else if(versionloop == 4){
                                                        backfit = new TF1("backfit",fAll,&AllFunctions::back4th,fitmin,fitmax,backnumber);
                                                        fityield = new TF1("fityield",fAll,&AllFunctions::function4th,1.6,2.4,polynumber);
                                                    }
                                                    //htemp->Reset();
                                                    for (int resetloop=1; resetloop <= htemp->GetNbinsX(); resetloop++) {
                                                        htemp->SetBinContent(resetloop,0);
                                                        htemp->SetBinError(resetloop,0);
                                                    }
                                                    for (int setloop=0; setloop < backnumber; setloop++) {
                                                        backfit->SetParameter(setloop,fitfunction->GetParameter(setloop));
                                                        backfit->SetParError(setloop,fitfunction->GetParError(setloop));
                                                    }
                                                    for (int setloop2=0; setloop2 < polynumber; setloop2++) {
                                                        if (setloop2 < backnumber) {
                                                            fityield->SetParameter(setloop2,0);
                                                            fityield->SetParError(setloop2,0);
                                                        }
                                                        else{
                                                            fityield->SetParameter(setloop2,fitfunction->GetParameter(setloop2));
                                                            fityield->SetParError(setloop2,fitfunction->GetParError(setloop2));
                                                        }
                                                    }//end of setloop2
                                                    if (particleloop == 0) {
                                                        yieldarray[typeloop][0]=((fityield->Integral(1.115683+0.493677,2.4))/(Bin1820Width));//KX
                                                        GMALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][(typeloop-1)+(2*cutloop*Nsigmaflag)]=yieldarray[typeloop][0];
                                                    }
                                                    else if(particleloop == 1){
                                                        yieldarray[typeloop][0]=((fityield->Integral(1.115683+.497611,2.4))/(Bin1820Width));//K0
                                                        GMALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][(typeloop-1)+(2*cutloop*Nsigmaflag)]=yieldarray[typeloop][0];
                                                    }
                                                    GMAlloftheYield->SetBinContent((typeloop),yieldarray[typeloop][0]);
                                                    for (int i=1; i <= ((fitmax-fitmin)/Bin1820Width); i++) {
                                                        double xbin = a->GetBinCenter(i+1+((fitmin-1.6)/Bin1820Width));
                                                        double fval = backfit->Eval(xbin);
                                                        double diff = (a->GetBinContent(i+1+((fitmin-1.6)/Bin1820Width)) - fval);
                                                        htemp->Fill(xbin,diff);
                                                        htemp->SetBinError(i+1+((fitmin-1.6)/Bin1820Width),a->GetBinError(i+1+((fitmin-1.6)/Bin1820Width)));
                                                    }
                                                }//end of subtraction
                                            }//old subtraction method
                                            else{//new subtraction method
                                                if (fitsubtract == 1) {
                                                    if (UseSigmaRange == 1) {
                                                        if (UseVoigt == 1) {
                                                            if (MCflag == 1) {
                                                                fityield = new TF1("fityield",fAll,&AllFunctions::funVoigtMC,fitmin,fitmax,5);//No polynomial new
                                                            }
                                                            else{
                                                                fityield = new TF1("fityield",fAll,&AllFunctions::funVoigt,fitmin,fitmax,4);//No polynomial new
                                                            }
                                                        }
                                                        else{
                                                            fityield = new TF1("fityield",fAll,&AllFunctions::BWfunction,fitmin,fitmax,3);//No polynomial
                                                        }
                                                    }
                                                    else{
                                                        fityield = new TF1("fityield",fAll,&AllFunctions::BWfunction,fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),3);//No polynomial
                                                    }
                                                    
                                                    if (versionloop == 0) {
                                                        backfit = new TF1("backfit",fAll,&AllFunctions::backline,fitmin,fitmax,backnumber);
                                                    }
                                                    else if(versionloop == 1){
                                                        backfit = new TF1("backfit",fAll,&AllFunctions::back2nd,fitmin,fitmax,backnumber);
                                                    }
                                                    else if(versionloop == 2){
                                                        backfit = new TF1("backfit",fAll,&AllFunctions::back3rd,fitmin,fitmax,backnumber);
                                                    }
                                                    else if(versionloop == 3){
                                                        backfit = new TF1("backfit",fAll,&AllFunctions::back4th,fitmin,fitmax,backnumber);
                                                    }
                                                    //htemp->Reset();
                                                    for (int resetloop=1; resetloop <= htemp->GetNbinsX(); resetloop++) {
                                                        htemp->SetBinContent(resetloop,0);
                                                        htemp->SetBinError(resetloop,0);
                                                    }
                                                    for (int setloop=0; setloop < backnumber; setloop++) {
                                                        backfit->SetParameter(setloop,fitfunction->GetParameter(setloop));
                                                        backfit->SetParError(setloop,fitfunction->GetParError(setloop));
                                                    }
                                                    if (UseVoigt == 1) {
                                                        fityield->SetParameter(0,fitfunction->GetParameter(backnumber));
                                                        fityield->SetParError(0,fitfunction->GetParError(backnumber));
                                                        fityield->SetParameter(1,fitfunction->GetParameter(2+backnumber));
                                                        fityield->SetParError(1,fitfunction->GetParError(2+backnumber));
                                                        if (typeloop == 1) {
                                                            /*if (MCflag == 1) {
                                                                fityield->FixParameter(2,MCSigma[cb][pb-MinPT]);
                                                            }
                                                            else{*/
                                                                fityield->FixParameter(2,0.002);
                                                            //}
                                                        }
                                                        else{
                                                            fityield->SetParameter(2,0.002);
                                                        }
                                                        fityield->SetParameter(3,fitfunction->GetParameter(1+backnumber));
                                                        fityield->SetParError(3,fitfunction->GetParError(1+backnumber));
                                                        //slope fixed
                                                        if (MCflag == 1) {
                                                            fityield->FixParameter(4,AEslope[cb][pb-MinPT]);
                                                        }
                                                        //slope fixed
                                                    }
                                                    else{
                                                        for (int setloop2=0; setloop2 < 3; setloop2++) {
                                                            fityield->SetParameter(setloop2,fitfunction->GetParameter(setloop2+backnumber));
                                                            fityield->SetParError(setloop2,fitfunction->GetParError(setloop2+backnumber));
                                                        }//end of setloop2
                                                    }
                                                    
                                                    if (UseSigmaRange == 1) {//use sigma range
                                                        if (Use1Type == 1) {//use 1 type and No poly
                                                            
                                                            GMALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][0]=(fityield->Integral(fitmin,fitmax)/(Bin1820Width));//KX
                                                            for (int i=1; i <= ((fitmax-fitmin)/(Bin1820Width)); i++) {//definate change
                                                                double xbin = a->GetBinCenter(i+1+((fitmin-1.6)/Bin1820Width));
                                                                double fval = backfit->Eval(xbin);
                                                                double diff = (a->GetBinContent(i+1+(((fitmin)-1.6)/Bin1820Width)) - fval);
                                                                htemp->Fill(xbin,diff);
                                                                htemp->SetBinError(i+1+(((fitmin)-1.6)/Bin1820Width),a->GetBinError(i+1+(((fitmin)-1.6)/Bin1820Width)));
                                                            }
                                                        }//use 1 type and No poly
                                                        else{//else
                                                            yieldarray[typeloop][0]=(fityield->Integral(fitmin,fitmax)/(Bin1820Width));//KX
                                                            GMALLYield[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=yieldarray[typeloop][0];
                                                            //GMAlloftheYield->SetBinContent((typeloop),yieldarray[typeloop][0]);
                                                            for (int i=1; i <= ((fitmax-fitmin)/(Bin1820Width)); i++) {//definate change
                                                                double xbin = a->GetBinCenter(i+1+((fitmin-1.6)/Bin1820Width));
                                                                double fval = backfit->Eval(xbin);
                                                                double diff = (a->GetBinContent(i+1+(((fitmin)-1.6)/Bin1820Width)) - fval);
                                                                htemp->Fill(xbin,diff);
                                                                htemp->SetBinError(i+1+(((fitmin)-1.6)/Bin1820Width),a->GetBinError(i+1+(((fitmin)-1.6)/Bin1820Width)));
                                                            }
                                                        }//else
                                                    }//use sigma range
                                                    else{
                                                        if (Use1Type == 1) {//use 1 type and No poly
                                                            GMALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][0]=(fityield->Integral(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))/(Bin1820Width));//KX
                                                            for (int i=1; i <= ((fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))-(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))))/(Bin1820Width)); i++) {//definate change
                                                                double xbin = a->GetBinCenter(i+1+((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))-1.6)/Bin1820Width));
                                                                double fval = backfit->Eval(xbin);
                                                                double diff = (a->GetBinContent(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width)) - fval);
                                                                htemp->Fill(xbin,diff);
                                                                htemp->SetBinError(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width),a->GetBinError(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width)));
                                                            }
                                                        }//use 1 type and No poly
                                                        else{//else
                                                            yieldarray[typeloop][0]=(fityield->Integral(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))/(Bin1820Width));//KX
                                                            GMALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][(typeloop-1)+(2*cutloop*Nsigmaflag)]=yieldarray[typeloop][0];
                                                            
                                                            GMAlloftheYield->SetBinContent((typeloop),yieldarray[typeloop][0]);
                                                            for (int i=1; i <= ((fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))-(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))))/(Bin1820Width)); i++) {//definate change
                                                                double xbin = a->GetBinCenter(i+1+((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))-1.6)/Bin1820Width));
                                                                double fval = backfit->Eval(xbin);
                                                                double diff = (a->GetBinContent(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width)) - fval);
                                                                htemp->Fill(xbin,diff);
                                                                htemp->SetBinError(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width),a->GetBinError(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width)));
                                                            }
                                                        }//else
                                                    }
                                                }//end of subtraction
                                            }//new subtraction method
                                            
                                            c->cd();
                                            gStyle->SetOptStat("emr");//adding stat box ()
                                            //gStyle->SetOptTitle(0);
                                            //gStyle->SetStatBorderSize(0);
                                            ///a->DrawCopy();//original to keep
                                            a->Draw();
                                            c->Update();
                                            fitfunction->SetLineColor(1);
                                            fitfunction->Draw("SAME");
                                            c->Update();
                                            backfit->SetLineColor(3);
                                            backfit->Draw("SAME");
                                            c->Update();
                                            
                                            if (c->GetPrimitive("stats") != NULL) {
                                                ps1 = (TPaveStats*)c->GetPrimitive("stats");
                                                if (!ps1) {
                                                    cout << "Error found **111111111111111111" << endl;
                                                }
                                                ps1->SetName("mystats");//mystats
                                                c->Modified();
                                                //ps1->SetOptStat(1110);
                                                TList *listofLines2 = ps1->GetListOfLines();
                                                if (!listofLines2) {
                                                    cout << "Error here **11" << endl;
                                                }
                                                tconst = ps1->GetLineWith("RMS");
                                                listofLines2->Remove(tconst);
                                                tconst = ps1->GetLineWith("Entries");
                                                listofLines2->Remove(tconst);
                                                tconst = ps1->GetLineWith("Mean");
                                                listofLines2->Remove(tconst);
                                                tconst = ps1->GetLineWith("Std Dev");
                                                listofLines2->Remove(tconst);
                                                //new
                                                //tconst = ps1->GetLineWith("mass_");
                                                //listofLines2->Remove(tconst);
                                            }
                                            else{
                                                ps1 = (TPaveStats*)c->GetPrimitive("mystats");
                                                if (!ps1) {
                                                    cout << "Error found **111111111111111111" << endl;
                                                }
                                                ps1->SetName("mystats");//mystats
                                                c->Modified();
                                                TList *listofLines2 = ps1->GetListOfLines();
                                                if (!listofLines2) {
                                                    cout << "Error here **11" << endl;
                                                }
                                                tconst = ps1->GetLineWith("Integral");
                                                listofLines2->Remove(tconst);
                                                tconst = ps1->GetLineWith("Mean");
                                                listofLines2->Remove(tconst);
                                                tconst = ps1->GetLineWith("Width");
                                                listofLines2->Remove(tconst);
                                                tconst = ps1->GetLineWith("chi2/ndf");
                                                listofLines2->Remove(tconst);
                                                //new
                                                //tconst = ps1->GetLineWith("mass_");
                                                //listofLines2->Remove(tconst);
                                            }
                                            a->SetStats(0);
                                            //ps1->AddText(Form("Integral = %5.1f",yieldarray[typeloop][0]));
                                            ps1->AddText(Form("Mean = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+2),fitfunction->GetParError(backnumber+2)));
                                            ps1->AddText(Form("Width = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+1),fitfunction->GetParError(backnumber+1)));
                                            ps1->AddText(Form("chi2/ndf = %5.4f",(fitfunction->GetChisquare())/(fitfunction->GetNDF())));
                                            //new section, moving stat box
                                            if (BellSwitch == 1) {
                                                ps1->SetX1NDC(0.6);
                                                ps1->SetY1NDC(0.6);
                                                ps1->SetY2NDC(0.9);
                                            }
                                            //new section
                                            c->Modified();//keep here
                                            //adding stat box
                                            if (hmtloop == 0) {//Normal
                                                if (cutloop == 0) {//no cut
                                                    c->Print(Form("plotsCompleteV3/%s/Type%i/sub_gm_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cb,pb,rangeloop,jNR,versionloop));
                                                }//no cut
                                                else{//cut
                                                    c->Print(Form("plotsCompleteV3/%s/Type%i/cut%i_sub_gm_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cutloop,cb,pb,rangeloop,jNR,versionloop));
                                                }//cut
                                            }//Normal
                                            else{//hmt
                                                if (cutloop == 0) {//no cut
                                                    c->Print(Form("plotsCompleteV3/%s/Type%i/hmt_sub_gm_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cb,pb,rangeloop,jNR,versionloop));
                                                }//no cut
                                                else{//cut
                                                    c->Print(Form("plotsCompleteV3/%s/Type%i/hmtcut%i_sub_gm_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cutloop,cb,pb,rangeloop,jNR,versionloop));
                                                }//cut
                                            }//hmt
                                            c->Clear();
                                            
                                            if (typeloop == 1) {
                                                if (rangeloop == 0) {
                                                    if (versionloop == 1) {
                                                        if (jNR == 2) {
                                                            if ((cb == 0) && (pb == MinPT)) {
                                                                cSubGM->Clear();
                                                                if (Use1PT == 1) {
                                                                    cSubGM->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                                                }
                                                                else{
                                                                    cSubGM->Divide(MaxCen,(MaxPT-MinPT));//4
                                                                }
                                                                cSubGM->Update();
                                                            }
                                                            cSubGM->cd((cb+1)+(MaxPT-MinPT)*(pb-MinPT));
                                                            a->DrawCopy();
                                                            cSubGM->Update();
                                                            //fitfunction->SetLineColor(1);
                                                            //fitfunction->DrawCopy("SAME");
                                                            //cSubGM->Update();
                                                            backfit->SetLineColor(3);
                                                            backfit->DrawCopy("SAME");
                                                            cSubGM->Update();
                                                            if ((cb == (MaxCen-1)) && (pb == (MaxPT-1))) {
                                                                if (hmtloop == 0) {//Normal
                                                                    if (cutloop == 0) {//no cut
                                                                        cSubGM->Print(Form("plotsCompleteV3/%s/Type%i/Sub_gm_allC_pt.pdf",partname,typeloop));
                                                                    }//no cut
                                                                    else{//cut
                                                                        cSubGM->Print(Form("plotsCompleteV3/%s/Type%i/cut%i_Sub_gm_allC_pt.pdf",partname,typeloop,cutloop));
                                                                    }//cut
                                                                }//Normal
                                                                else{//hmt
                                                                    if (cutloop == 0) {//no cut
                                                                        cSubGM->Print(Form("plotsCompleteV3/%s/Type%i/hmt_Sub_gm_allC_pt.pdf",partname,typeloop));
                                                                    }//no cut
                                                                    else{//cut
                                                                        cSubGM->Print(Form("plotsCompleteV3/%s/Type%i/hmtcut%i_Sub_gm_allC_pt.pdf",partname,typeloop,cutloop));
                                                                    }//cut
                                                                }//hmt
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            c->cd();
                                            if (fitsubtract == 1) {
                                                if (UseNoPoly == 1) {//use No poly
                                                    fityield->SetParLimits(0, 0.0, 500000.0);
                                                    if (UseVoigt == 1) {
                                                        fityield->SetParLimits(3,0.01,0.06);
                                                        //fityield->SetParLimits(2,0.0,0.06);
                                                        fityield->SetParLimits(1,1.8,1.85);
                                                        fityield->SetParameter(3,0.025);
                                                        fityield->SetParLimits(2,0.001,0.003);
                                                        if (typeloop == 1) {
                                                            /*if (MCflag == 1) {
                                                                fityield->FixParameter(2,MCSigma[cb][pb-MinPT]);
                                                            }
                                                            else{*/
                                                                fityield->FixParameter(2,0.002);
                                                            //}
                                                        }
                                                        else{
                                                            fityield->SetParameter(2,0.002);
                                                        }
                                                        fityield->SetParameter(1,1.82);
                                                        if (Fixmean == 1) {
                                                            fityield->FixParameter(1,1.8218);
                                                        }
                                                        
                                                        if (MCflag == 1) {
                                                            fityield->FixParameter(4,AEslope[cb][pb-MinPT]);
                                                        }
                                                    }
                                                    else{
                                                        fityield->SetParLimits(1,0.01,0.06);
                                                        fityield->SetParLimits(2,1.8,1.85);
                                                        fityield->SetParameter(1,0.025);
                                                        fityield->SetParameter(2,1.82);
                                                    }
                                                    
                                                }//use no poly
                                                else{
                                                    fitfunction->SetParLimits(backnumber, 0.0, 500000.0);
                                                    fitfunction->SetParLimits(backnumber+1,0.01,0.06);
                                                    fitfunction->SetParLimits(backnumber+2,1.8,1.85);
                                                    fitfunction->SetParameter(backnumber+1,0.025);
                                                    fitfunction->SetParameter(backnumber+2,1.82);
                                                }
                                            }
                                            //new linear section
                                            if (UseNoPoly == 1) {
                                                for (int removeloop=0; removeloop < backnumber; removeloop++) {
                                                    fitfunction->FixParameter(removeloop,0);
                                                }
                                            }
                                            //new linear section
                                            if (UseNoPoly == 1) {
                                                if (UseSigmaRange == 1) {
                                                    //fityield->SetParameter(0,1000);
                                                    htemp->Fit("fityield","Q","SAME",fitmin2,fitmax2);
                                                    //cout << "usenopoly = 1 usesigmarange = 1 mean = " << fityield->GetParameter(1) << endl;;
                                                }
                                                else{
                                                    htemp->Fit("fityield","Q","SAME",fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)));
                                                }
                                            }
                                            else{
                                                htemp->Fit("fitfunction","Q","SAME",fitmin,fitmax);
                                            }
                                            
                                            if (UseNoPoly == 1) {
                                                for (int setloop=0; setloop < backnumber; setloop++) {
                                                    backfit->FixParameter(setloop,0);
                                                }
                                                backfit->SetLineColor(kBlue);
                                                if (UseSigmaRange == 1) {
                                                    htemp->Fit("backfit","QR+","SAME",fitmin,fitmax);
                                                }
                                                else{
                                                    htemp->Fit("backfit","QR+","SAME",fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)));
                                                }
                                                
                                            } else {
                                                for (int setloop=0; setloop < backnumber; setloop++) {
                                                    backfit->FixParameter(setloop,fitfunction->GetParameter(setloop));
                                                }
                                                backfit->SetLineColor(kBlue);
                                                htemp->Fit("backfit","QR+","SAME",fitmin,fitmax);
                                            }
                                            c->cd();
                                            gStyle->SetOptStat("emr");//()
                                            htemp->Draw();
                                            c->Update();
                                            
                                            if (c->GetPrimitive("stats") != NULL) {
                                                ps1 = (TPaveStats*)c->GetPrimitive("stats");
                                                if (!ps1) {
                                                    cout << "Error found 111111111111111111" << endl;
                                                }
                                                ps1->SetName("mystats");//mystats
                                                TList *listofLines1 = ps1->GetListOfLines();
                                                tconst = ps1->GetLineWith("RMS");
                                                listofLines1->Remove(tconst);
                                                tconst = ps1->GetLineWith("Entries");
                                                listofLines1->Remove(tconst);
                                                tconst = ps1->GetLineWith("Mean");
                                                listofLines1->Remove(tconst);
                                                tconst = ps1->GetLineWith("Std Dev");
                                                listofLines1->Remove(tconst);
                                            }
                                            else{
                                                ps1 = (TPaveStats*)c->GetPrimitive("mystats");
                                                if (c->GetPrimitive("mystats") == NULL) {
                                                    cout << "No primitive" << endl;
                                                }
                                                
                                                if (!ps1) {
                                                    cout << "Error found 111111111111111111" << endl;
                                                }
                                                ps1->SetName("mystats");//mystats
                                                c->Modified();
                                                TList *listofLines1 = ps1->GetListOfLines();
                                                if (!listofLines1) {
                                                    cout << "Error here 11" << endl;
                                                }
                                                tconst = ps1->GetLineWith("Integral");
                                                listofLines1->Remove(tconst);
                                                tconst = ps1->GetLineWith("Mean");
                                                listofLines1->Remove(tconst);
                                                tconst = ps1->GetLineWith("Width");
                                                listofLines1->Remove(tconst);
                                                if (UseVoigt == 1) {
                                                    tconst = ps1->GetLineWith("Sigma");
                                                    listofLines1->Remove(tconst);
                                                }
                                                tconst = ps1->GetLineWith("chi2/ndf");
                                                listofLines1->Remove(tconst);
                                                tconst = ps1->GetLineWith("Yield Bin");
                                                listofLines1->Remove(tconst);
                                            }
                                            
                                            axis = htemp->GetXaxis();
                                            if (UseNoPoly == 1) {
                                                bmax = axis->FindBin(1.82);
                                                //double integral50 = fitfunction->Integral(1.77,1.87)/axis->GetBinWidth(bmax);
                                                //if (flag50 == 1) {
                                                integralbin=0;
                                                if (UseSigmaRange == 1) {//use sigma range
                                                    integral = fityield->Integral(fitmin2,fitmax2)/axis->GetBinWidth(bmax);
                                                    //double integralbin=0;//
                                                    //double integralbackbin=0;
                                                    for (i=1; i <= htemp->GetNbinsX(); i++) {
                                                        binvalue=htemp->GetBinContent(i);
                                                        if (binvalue == 0.0) {
                                                            //cout << "i " << i << " = Null bin value " << binvalue << " bin center " << htemp->GetBinCenter(i) << endl;
                                                        }
                                                        else{
                                                            //cout << "i " << i << " not Null bin value " << binvalue << " bin center " << htemp->GetBinCenter(i) << endl;
                                                            if (((i*Bin1820Width+1.60) >= fitmin2) && ((i*Bin1820Width+1.60) <= fitmax2)) {//range2
                                                                integralbin+=binvalue;
                                                            }
                                                        }
                                                    }
                                                    BinGMALLYield[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=integralbin;
                                                    //cout << "preliminary bin counting: integral bin " << integralbin << " inregral back bin " << integralbackbin << " S/B " << integralbin/integralbackbin << endl;
                                                    //insert bin counting method
                                                }//use sigma range
                                                else{
                                                    integral = fityield->Integral(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))/axis->GetBinWidth(bmax);
                                                }
                                                
                                                //}
                                                htemp->SetStats(0);
                                                ps1->AddText(Form("Integral = %5.1f",integral));
                                                if (UseVoigt == 1) {
                                                    ps1->AddText(Form("Mean = %5.4f +- %5.5f",fityield->GetParameter(1),fityield->GetParError(1)));
                                                    ps1->AddText(Form("Width = %5.4f +- %5.5f",fityield->GetParameter(3),fityield->GetParError(3)));
                                                    ps1->AddText(Form("Sigma = %5.4f +- %5.5f",fityield->GetParameter(2),fityield->GetParError(2)));
                                                }
                                                else{
                                                    ps1->AddText(Form("Mean = %5.4f +- %5.5f",fityield->GetParameter(2),fityield->GetParError(2)));
                                                    ps1->AddText(Form("Width = %5.4f +- %5.5f",fityield->GetParameter(1),fityield->GetParError(1)));
                                                }
                                                
                                                ps1->AddText(Form("chi2/ndf = %5.4f",(fityield->GetChisquare())/(fityield->GetNDF())));
                                                ps1->AddText(Form("Yield Bin = %5.1f",integralbin));
                                                if (BellSwitch == 1) {
                                                    ps1->SetX1NDC(0.6);
                                                    ps1->SetY1NDC(0.6);
                                                    ps1->SetY2NDC(0.9);
                                                }
                                            }
                                            else{
                                                int bmin = axis->FindBin(1.77);
                                                int bmax = axis->FindBin(1.87);
                                                //double integral50 = fitfunction->Integral(1.77,1.87)/axis->GetBinWidth(bmax);
                                                //if (flag50 == 1) {
                                                double integral = fitfunction->Integral(1.77,1.87)/axis->GetBinWidth(bmax);
                                                //}
                                                //double integralback = backfit->Integral(1.77,1.87)/axis->GetBinWidth(bmax);
                                                htemp->SetStats(0);
                                                //if (flag50 == 1) {
                                                ps1->AddText(Form("Integral = %5.1f",integral));
                                                //}
                                                ps1->AddText(Form("Mean = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+2),fitfunction->GetParError(backnumber+2)));
                                                ps1->AddText(Form("Width = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+1),fitfunction->GetParError(backnumber+1)));
                                                ps1->AddText(Form("chi2/ndf = %5.4f",(fitfunction->GetChisquare())/(fitfunction->GetNDF())));
                                                if (BellSwitch == 1) {
                                                    ps1->SetX1NDC(0.6);
                                                    ps1->SetY1NDC(0.6);
                                                    ps1->SetY2NDC(0.9);
                                                }
                                            }
                                            
                                            if (UseLinear == 1) {//use linear
                                                if (UseNoPoly == 1) {//use no poly
                                                    if (UseVoigt == 1) {
                                                        GMALLMean[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fityield->GetParameter(1);
                                                        GMALLMeanError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fityield->GetParError(1);
                                                        GMALLWidth[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fityield->GetParameter(3);
                                                        GMALLWidthError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fityield->GetParError(3);
                                                        GMALLChi2[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=(fityield->GetChisquare()/fityield->GetNDF());
                                                        GMALLYield[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=integral;
                                                    }
                                                    else{
                                                        GMALLMean[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParameter(2);
                                                        GMALLMeanError[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParError(2);
                                                        GMALLWidth[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParameter(1);
                                                        GMALLWidthError[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParError(1);
                                                        GMALLChi2[cb][pb-MinPT][jNR-NRBegin][0][0][0]=(fityield->GetChisquare()/fityield->GetNDF());
                                                        GMALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][0]=integral;
                                                    }
                                                    
                                                }//use no poly
                                                else{
                                                    GMALLMean[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParameter(backnumber+2);
                                                    GMALLMeanError[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParError(backnumber+2);
                                                    GMALLWidth[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParameter(backnumber+1);
                                                    GMALLWidthError[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParError(backnumber+1);
                                                    GMALLChi2[cb][pb-MinPT][jNR-NRBegin][0][0][0]=(fityield->GetChisquare()/fityield->GetNDF());
                                                    GMALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][0]=integral;
                                                }
                                                
                                                //}
                                            }//use linear
                                            
                                            if (Use1Type != 1) {//not 1type
                                                xarray[0]=versionloop;
                                                integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=integral;
                                                meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=fitfunction->GetParameter(backnumber+2);
                                                widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=fitfunction->GetParameter(backnumber+1);
                                                xarrayerror[0]=0;
                                                integralarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=0;
                                                meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=fitfunction->GetParError(backnumber+2);
                                                widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=fitfunction->GetParError(backnumber+1);
                                                
                                                integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                integralarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=0;
                                                yieldarray[typeloop][0]=yieldarray[typeloop][0];
                                                meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                
                                                /*if (versionloop == 1) {
                                                 integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 integralarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=0;
                                                 yieldarray[typeloop][0]=yieldarray[typeloop][0];
                                                 meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]*meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]*widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 }
                                                 else if (versionloop == (versionmax-1)) {
                                                 integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]/(1.0*(versionmax-1));
                                                 meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]/(1.0*(versionmax-1));
                                                 widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]/(1.0*(versionmax-1));
                                                 integralarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=0;
                                                 yieldarray[typeloop][0]+=yieldarray[typeloop][0];
                                                 meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=(meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]*meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]);
                                                 widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=(widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]*widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]);
                                                 meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=sqrt(meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]/((versionmax-1)*1.0));
                                                 widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]=sqrt(widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]/((versionmax-1)*1.0));
                                                 yieldarray[typeloop][0]=yieldarray[typeloop][0]/(1.0*(versionmax-1));
                                                 }
                                                 else{
                                                 integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=integralarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=meanarray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=widtharray[(typeloop-1)+(2*cutloop*Nsigmaflag)][0];
                                                 integralarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=0;
                                                 yieldarray[typeloop][0]+=yieldarray[typeloop][0];
                                                 meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=(meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]*meanarrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]);
                                                 widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]+=(widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]*widtharrayerror[(typeloop-1)+(2*cutloop*Nsigmaflag)][0]);
                                                 }*/
                                                
                                                if (typeloop == 1) {
                                                    yieldarray[0][0]=0;
                                                }
                                                
                                                /*if (typeloop == 8) {
                                                 systemxarray[typeloop-7]=typeloop-7;
                                                 systemxarrayerror[typeloop-7]=0;
                                                 systemwidth[typeloop-7]=fitfunction->GetParameter(backnumber+1);
                                                 systemwidtherror[typeloop-7]=fitfunction->GetParError(backnumber+1);
                                                 systemmean[typeloop-7]=fitfunction->GetParameter(backnumber+2);
                                                 systemmeanerror[typeloop-7]=fitfunction->GetParError(backnumber+2);
                                                 yieldarray[0][0]+=yieldarray[typeloop][0];
                                                 yieldarray[0][0]=yieldarray[0][0]/(1.0*8);
                                                 }
                                                 else if ((typeloop <= 7) && (typeloop >= 5)) {//collect widths for systematic
                                                 systemwidth[typeloop-3]=fitfunction->GetParameter(backnumber+1);
                                                 systemwidtherror[typeloop-3]=fitfunction->GetParError(backnumber+1);
                                                 yieldarray[0][0]+=yieldarray[typeloop][0];
                                                 }//collect widths for systematic
                                                 else if ((typeloop != 1) && (typeloop <= 4)) {//collect mean for systematic
                                                 systemxarray[typeloop]=typeloop;
                                                 systemxarrayerror[typeloop]=0;
                                                 systemmean[typeloop]=fitfunction->GetParameter(backnumber+2);
                                                 systemmeanerror[typeloop]=fitfunction->GetParError(backnumber+2);
                                                 yieldarray[0][0]+=yieldarray[typeloop][0];
                                                 }//collect mean for systematic*/
                                                
                                            }//not 1type
                                            
                                            fitfunction->ReleaseParameter(backnumber+1);
                                            fitfunction->ReleaseParameter(backnumber+2);
                                            for (int setloop=0; setloop < backnumber; setloop++) {
                                                fitfunction->ReleaseParameter(setloop);
                                                fitline->ReleaseParameter(setloop);
                                                backfit->ReleaseParameter(setloop);
                                            }
                                            c->Modified();
                                            if (hmtloop == 0) {//Normal
                                                if (cutloop == 0) {//no cut
                                                    c->Print(Form("plotsCompleteV3/%s/Type%i/linear_sub_gm_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cb,pb,rangeloop,jNR,versionloop));
                                                }//no cut
                                                else{//cut
                                                    c->Print(Form("plotsCompleteV3/%s/Type%i/cut%i_linear_sub_gm_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cutloop,cb,pb,rangeloop,jNR,versionloop));
                                                }//cut
                                            }//Normal
                                            else{//hmt
                                                if (cutloop == 0) {//no cut
                                                    c->Print(Form("plotsCompleteV3/%s/Type%i/hmt_linear_sub_gm_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cb,pb,rangeloop,jNR,versionloop));
                                                }//no cut
                                                else{//cut
                                                    c->Print(Form("plotsCompleteV3/%s/Type%i/hmtcut%i_linear_sub_gm_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cutloop,cb,pb,rangeloop,jNR,versionloop));
                                                }//cut
                                            }//hmt
                                            if (typeloop == 1) {
                                                if (rangeloop == 0) {
                                                    if (versionloop == 1) {
                                                        if (jNR == 2) {
                                                            if ((cb == 0) && (pb == MinPT)) {
                                                                cLinearGM->Clear();
                                                                if (Use1PT == 1) {
                                                                    cLinearGM->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                                                }
                                                                else{
                                                                    cLinearGM->Divide(MaxCen,(MaxPT-MinPT));//4
                                                                }
                                                                cLinearGM->Update();
                                                            }
                                                            cLinearGM->cd((cb+1)+((MaxPT-MinPT)*(pb-MinPT)));
                                                            htemp->DrawCopy();
                                                            /*if ((psTemp = (TPaveStats*)htemp->FindObject("stats")) == NULL) {
                                                                cout << "THIS DID NOT WORK" << endl;
                                                            }
                                                            psTemp = (TPaveStats*)htemp->FindObject("stats");
                                                            psTemp->SetX1NDC(0.6);
                                                            psTemp->SetY1NDC(0.5);
                                                            psTemp->SetY2NDC(0.9);*/
                                                            cLinearGM->Update();
                                                            if ((cb == (MaxCen-1)) && (pb == (MaxPT-1))) {
                                                                if (hmtloop == 0) {//Normal
                                                                    if (cutloop == 0) {//no cut
                                                                        cLinearGM->Print(Form("plotsCompleteV3/%s/Type%i/Linear_gm_allC_pt.pdf",partname,typeloop));
                                                                    }//no cut
                                                                    else{//cut
                                                                        cLinearGM->Print(Form("plotsCompleteV3/%s/Type%i/cut%i_Linear_gm_allC_pt.pdf",partname,typeloop,cutloop));
                                                                    }//cut
                                                                }//Normal
                                                                else{//hmt
                                                                    if (cutloop == 0) {//no cut
                                                                        cLinearGM->Print(Form("plotsCompleteV3/%s/Type%i/hmt_Linear_gm_allC_pt.pdf",partname,typeloop));
                                                                    }//no cut
                                                                    else{//cut
                                                                        cLinearGM->Print(Form("plotsCompleteV3/%s/Type%i/hmtcut%i_Linear_gm_allC_pt.pdf",partname,typeloop,cutloop));
                                                                    }//cut
                                                                }//hmt
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                            
                                            fitfunction->ReleaseParameter(backnumber+1);
                                            fitfunction->ReleaseParameter(backnumber+2);
                                            //new linear section
                                            if (UseNoPoly == 1) {
                                                for (removeloop=0; removeloop <= backnumber; removeloop++) {
                                                    fitfunction->ReleaseParameter(removeloop);
                                                }
                                                for (int setloop2=0; setloop2 < 3; setloop2++) {
                                                    fityield->ReleaseParameter(setloop2);
                                                }
                                                if (UseVoigt == 1) {
                                                    fityield->ReleaseParameter(3);
                                                    //slope
                                                    fityield->ReleaseParameter(4);
                                                    //slope
                                                }
                                            }
                                            //new linear section
                                            //if (versionloop == 1) {
                                            c1->Clear();
                                            //c1->Divide(1,3);
                                            //}
                                            c1->cd();
                                            htemp->DrawCopy();
                                            c1->Update();
                                            c->cd();
                                            for (int setloop=0; setloop < backnumber; setloop++) {
                                                backfit->ReleaseParameter(setloop);
                                            }
                                        }
                                        gStyle->SetOptStat(0);
                                        fitfunction->ReleaseParameter(backnumber+1);
                                        fitfunction->ReleaseParameter(backnumber+2);
                                    }//end of type loop
                                }//end of versionloop
                            }//end of GM NR loop
                        }///gm
                        
                        ////unlikemark
                        TitlesBS2(2);
                        //if ((UseGMonly == 1) && (particleloop == 0)) {
                        //    continue;
                        //}
                        for(jNR=NRstart;jNR<nNR;jNR++){//plot unlike-mixunlike
                            
                            if (particleloop == 0) {
                                //versionloop=2;//to be determined
                                
                                ///jNR=0;//best estamate
                                
                                /*if(pb == 0){
                                 rangeloop=3;
                                 }
                                 else if(pb == 1){
                                 rangeloop=3;
                                 }
                                 else if(pb == 2){
                                 rangeloop=4;
                                 }
                                 else if(pb == 3){
                                 rangeloop=4;
                                 }*/
                            }
                            else if(particleloop == 1){
                                //versionloop=2;//to be determined
                                /*if (pb == 0) {
                                 jNR=2;
                                 rangeloop=4;
                                 }
                                 else if(pb == 1){
                                 jNR=2;
                                 rangeloop=4;
                                 }
                                 else if(pb == 2){
                                 jNR=2;
                                 rangeloop=4;
                                 }
                                 else if(pb == 3){*/
                                
                                //jNR=0;//best estamate
                                
                                //rangeloop=4;
                                //}
                            }
                            //sudden change
                            ///rangeloop=2;
                            //sudden change
                            
                            if (rangeloop == 1) {
                                fitmin=1.7;
                                fitmax=2.0;
                            }
                            else if(rangeloop == 2){
                                fitmin=1.7;
                                fitmax=1.95;
                            }
                            else if(rangeloop == 3){
                                fitmin=1.73;
                                fitmax=2.0;
                            }
                            else if(rangeloop == 4){
                                fitmin=1.73;
                                fitmax=1.95;
                            }
                            
                            if (UseSigmaRange == 1) {//use sigma range
                                //fitmin=(1.823-6.0*(0.024/2.35));
                                //fitmax=(1.823+6.0*(0.024/2.35));
                                
                                fitmin2=(1.823-4.0*(0.024/2.35));//4.0
                                fitmax2=(1.823+4.0*(0.024/2.35));//4.0
                                if (rangeloop == 0) {
                                    fitmin=(1.7);//changed from 1.760 to 1.7
                                    fitmax=(2.0);
                                }
                                else if(rangeloop == 1){
                                    fitmin=(1.760);
                                    fitmax=(1.9);
                                }
                                else if(rangeloop == 2){
                                    fitmin=(1.760);//changed from 1.7 to 1.760
                                    fitmax=(2.0);
                                }
                            }//use sigma range
                            
                            /*if ((jNR == 0) && (jDC == Klkx)) {
                             jNR=1;
                             }*/
                            
                            flagSB=0;
                            sprintf(s0,"mass_unlikeQmixunlike%i_c%i_p%i",jNR,cb,pb);
                            a=(TH1D*) f->Get(s0);
                            if(!a){cerr<<"missing histogram "<<s0<<endl; continue;}
                            for (versionloop=0; versionloop <= (versionmax-1); versionloop ++) {//begin version loop

                                for (typeloop=typestart; typeloop <= typeend; typeloop++) {//begin typeloop for unlike
                                    //if (versionloop == 0) {
                                    if (Use1Type == 1) {
                                        if (typeloop == 8) {
                                            a->Rebin(rb);
                                            //a->GetXaxis()->SetRange(0,50);
                                            a->GetXaxis()->SetRange(a->GetXaxis()->FindBin(1.70),a->GetXaxis()->FindBin(2.10));
                                            Bin1820 = a->GetXaxis()->FindBin(1.82);
                                            Bin1820Width = a->GetXaxis()->GetBinWidth(Bin1820);
                                        }
                                    }
                                    else{
                                        if (typeloop == 1) {
                                            if (rangeloop == 0) {
                                                if (versionloop == 0) {
                                                    //if (jNR == NRstart) {
                                                    a->Rebin(rb);
                                                    a->GetXaxis()->SetRange(a->GetXaxis()->FindBin(1.70),a->GetXaxis()->FindBin(2.10));
                                                    Bin1820 = a->GetXaxis()->FindBin(1.82);
                                                    Bin1820Width = a->GetXaxis()->GetBinWidth(Bin1820);
                                                    //}
                                                }
                                            }
                                        }
                                    }
                                    
                                    a->SetTitle(Form("%s, %1.1f-%1.1f %% , %1.1f<#it{p}_{T}<%1.1f GeV/c",DC1,cbmin[cb],cbmax[cb],ptmin[pb],ptmax[pb]));
                                    a->SetXTitle("Invariant Mass [GeV/#it{c^2}]");
                                    a->GetXaxis()->SetTitleOffset(1.3);
                                    a->SetYTitle("Counts/(0.01 GeV/c^{2})");
                                    a->GetYaxis()->SetTitleOffset(1.3);
                                    a->SetMarkerStyle(1);
                                    a->SetMarkerColor(TColor::GetColor(color));
                                    a->SetLineColor(TColor::GetColor(color));
                                    //}
                                    if (Use1Type != 1) {//not 1 type
                                        AlloftheYield->GetXaxis()->SetBinLabel((typeloop),Form("peak%i",typeloop));
                                        AlloftheMean->GetXaxis()->SetBinLabel((typeloop),Form("peak%i",typeloop));
                                        AlloftheWidth->GetXaxis()->SetBinLabel((typeloop),Form("peak%i",typeloop));
                                        Double_t xarray[versionmax];
                                        Double_t integralarray[2][versionmax];
                                        Double_t meanarray[2][versionmax];
                                        Double_t widtharray[2][versionmax];
                                        Double_t xarrayerror[versionmax];
                                        Double_t integralarrayerror[2][versionmax];
                                        Double_t meanarrayerror[2][versionmax];
                                        Double_t widtharrayerror[2][versionmax];
                                        Double_t systemxarray[2];
                                        Double_t systemxarrayerror[2];
                                        Double_t systemmean[2];
                                        Double_t systemmeanerror[2];
                                        Double_t systemwidth[2];
                                        Double_t systemwidtherror[2];
                                        Double_t yieldarray[3][versionmax];
                                        Double_t yieldplotversion[versionmax];
                                        Double_t yieldplottype[3];
                                        
                                        xarray[0]=0;
                                        xarrayerror[0]=0;
                                        for (int i=0; i <= 1; i++) {
                                            systemxarray[i]=i;
                                            systemxarrayerror[i]=0;
                                        }
                                    }//not 1 type

                                    ALLYield[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                    ALLMean[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                    ALLMeanError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                    ALLWidth[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                    ALLWidthError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                    ALLChi2[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=0;
                                    backnumber=0;
                                    if (fitsubtract == 1) {
                                        htemp = (TH1D*)a->Clone();//new section
                                        htemp2 = (TH1D*)a->Clone();
                                    }
                                    
                                    if (fitflag == 1) {
                                        //double par[8]=0;
                                        if (versionloop == 0) {
                                            polynumber=5;
                                            backnumber=2;
                                        }
                                        else if(versionloop == 1){
                                            polynumber=6;
                                            backnumber=3;
                                        }
                                        else if(versionloop == 2){
                                            polynumber=7;
                                            backnumber=4;
                                        }
                                        else if(versionloop == 3){
                                            polynumber=8;
                                            backnumber=5;
                                        }
                                        
                                        for (int z = 0; z < (polynumber+1); z++) {
                                            par[z]=0;
                                        }

                                        if (versionloop == 0) {
                                            fitfunction = new TF1("fitfunction",fAll,&AllFunctions::function,fitmin,fitmax,polynumber);
                                            //if ((typeloop == 1) || (typeloop == 8)) {
                                            if (UseSigmaRange2 == 1) {
                                                fitline = new TF1("fitline",fAll,&AllFunctions::backline,fitmin,fitmax,backnumber);//new
                                            }
                                            else{
                                                fitline = new TF1("fitline",fAll,&AllFunctions::fline,fitmin,fitmax,backnumber);//new
                                            }
                                            //}
                                        }
                                        else if(versionloop == 1){
                                            fitfunction = new TF1("fitfunction",fAll,&AllFunctions::function2nd,fitmin,fitmax,polynumber);
                                            //if ((typeloop == 1) || (typeloop == 8)) {
                                            if (UseSigmaRange2 == 1) {
                                                fitline = new TF1("fitline",fAll,&AllFunctions::back2nd,fitmin,fitmax,backnumber);//new
                                            }
                                            else{
                                                fitline = new TF1("fitline",fAll,&AllFunctions::f2nd,fitmin,fitmax,backnumber);//new
                                            }
                                            //}
                                        }
                                        else if(versionloop == 2){
                                            fitfunction = new TF1("fitfunction",fAll,&AllFunctions::function3rd,fitmin,fitmax,polynumber);
                                            //if ((typeloop == 1) || (typeloop == 8)) {
                                            if (UseSigmaRange2 == 1) {
                                                fitline = new TF1("fitline",fAll,&AllFunctions::back3rd,fitmin,fitmax,backnumber);//new
                                            }
                                            else{
                                                fitline = new TF1("fitline",fAll,&AllFunctions::f3rd,fitmin,fitmax,backnumber);//new
                                            }
                                            //}
                                        }
                                        else if(versionloop == 3){
                                            fitfunction = new TF1("fitfunction",fAll,&AllFunctions::function4th,fitmin,fitmax,polynumber);
                                            //if ((typeloop == 1) || (typeloop == 8)) {
                                            if (UseSigmaRange2 == 1) {
                                                fitline = new TF1("fitline",fAll,&AllFunctions::back4th,fitmin,fitmax,backnumber);//new
                                            }
                                            else{
                                                fitline = new TF1("fitline",fAll,&AllFunctions::f4th,fitmin,fitmax,backnumber);//new
                                            }
                                            //}
                                        }

                                        fitfunction->SetParLimits(backnumber, 0.0, 500000.0);
                                        fitfunction->ReleaseParameter(backnumber+1);
                                        fitfunction->ReleaseParameter(backnumber+2);
                                        fitfunction->SetParLimits(backnumber+1,0.01,0.06);
                                        fitfunction->SetParLimits(backnumber+2,1.8,1.85);
                                        fitfunction->SetParameter(backnumber+1,0.025);
                                        fitfunction->SetParameter(backnumber+2,1.82);
                                        fitfunction->SetNpx(500);
                                        fitline->SetNpx(500);
                                        //if ((typeloop == 1) || (typeloop == 8)) {
                                        a->Fit("fitline","Q","SAME",fitmin,fitmax);
                                        //}
                                        
                                        if ((fitline->GetChisquare()/fitline->GetNDF()) >= 10) {
                                            cout << "***********************ERROR in FITLine****************************************" << endl;
                                            a->Fit("fitline","Q0","SAME",fitmin,fitmax);
                                            if ((fitline->GetChisquare()/fitline->GetNDF()) >= 10) {
                                                //cout << "refit failed, reject data" << endl;
                                                //continue;
                                            }
                                            else{
                                                cout << "success in refiting" << endl;
                                            }
                                        }
                                        
                                        for (int j=0; j < backnumber; j++) {
                                            if (typeloop <= 7) {
                                                fitfunction->FixParameter(j,fitline->GetParameter(j));
                                            }
                                            else{
                                                fitfunction->SetParameter(j,fitline->GetParameter(j));
                                            }
                                        }
                                        
                                        a->Fit("fitfunction","Q0","SAME",fitmin,fitmax);
                                        if ((fitfunction->GetChisquare()/fitfunction->GetNDF()) >= 10) {
                                            for (int refitloop=0; refitloop <= 10; refitloop++) {
                                                if ((fitfunction->GetChisquare()/fitfunction->GetNDF()) >= 10) {
                                                    cout << "Begin refit" << endl;
                                                    a->Fit("fitfunction","Q0","SAME",fitmin,fitmax);
                                                }
                                            }
                                            if ((fitfunction->GetChisquare()/fitfunction->GetNDF()) >= 10) {
                                                cout << "ERROR With Refit procedure ****************************************************" << endl;
                                                //cout << "reject data" << endl;
                                                //continue;
                                            }
                                        }
                                        
                                        cTotal2->cd();
                                        gPad->SetRightMargin(0.01);
                                        if (typeloop == 1) {
                                            a->Draw();
                                            cTotal2->Update();
                                        }
                                        fitfunction->SetLineColor(typeloop);
                                        fitfunction->SetLineStyle(1);
                                        fitfunction->SetLineWidth(1);
                                        fitfunction->DrawCopy("SAME");
                                        cTotal2->Update();

                                        ALLMean[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fitfunction->GetParameter(backnumber+2);
                                        ALLMeanError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fitfunction->GetParError(backnumber+2);
                                        ALLWidth[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fitfunction->GetParameter(backnumber+1);
                                        ALLWidthError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fitfunction->GetParError(backnumber+1);
                                        ALLChi2[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=(fitfunction->GetChisquare()/fitfunction->GetNDF());
                                        
                                        if (UseNoPoly != 1) {//old subtraction method
                                            if (fitsubtract == 1) {
                                                if (versionloop == 1) {
                                                    backfit = new TF1("backfit",fAll,&AllFunctions::backline,fitmin,fitmax,backnumber);
                                                    fityield = new TF1("fityield",fAll,&AllFunctions::function,1.6,2.4,polynumber);
                                                }
                                                else if(versionloop == 2){
                                                    backfit = new TF1("backfit",fAll,&AllFunctions::back2nd,fitmin,fitmax,backnumber);
                                                    fityield = new TF1("fityield",fAll,&AllFunctions::function2nd,1.6,2.4,polynumber);
                                                    
                                                }
                                                else if(versionloop == 3){
                                                    backfit = new TF1("backfit",fAll,&AllFunctions::back3rd,fitmin,fitmax,backnumber);
                                                    fityield = new TF1("fityield",fAll,&AllFunctions::function3rd,1.6,2.4,polynumber);
                                                }
                                                else if(versionloop == 4){
                                                    backfit = new TF1("backfit",fAll,&AllFunctions::back4th,fitmin,fitmax,backnumber);
                                                    fityield = new TF1("fityield",fAll,&AllFunctions::function4th,1.6,2.4,polynumber);
                                                }
                                                //htemp->Reset();
                                                for (int resetloop=1; resetloop <= htemp->GetNbinsX(); resetloop++) {
                                                    htemp->SetBinContent(resetloop,0);
                                                    htemp->SetBinError(resetloop,0);
                                                }
                                                for (int setloop=0; setloop < backnumber; setloop++) {
                                                    backfit->SetParameter(setloop,fitfunction->GetParameter(setloop));
                                                    backfit->SetParError(setloop,fitfunction->GetParError(setloop));
                                                }
                                                for (int setloop2=0; setloop2 < polynumber; setloop2++) {
                                                    if (setloop2 < backnumber) {
                                                        fityield->SetParameter(setloop2,0);
                                                        fityield->SetParError(setloop2,0);
                                                    }
                                                    else{
                                                        fityield->SetParameter(setloop2,fitfunction->GetParameter(setloop2));
                                                        fityield->SetParError(setloop2,fitfunction->GetParError(setloop2));
                                                    }
                                                }//end of setloop2
                                                if (particleloop == 0) {
                                                    yieldarray[typeloop][0]=((fityield->Integral(1.115683+0.493677,2.4))/(Bin1820Width));//KX
                                                    ALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][(typeloop-1)+(2*cutloop*Nsigmaflag)]=yieldarray[typeloop][0];
                                                }
                                                else if(particleloop == 1){
                                                    yieldarray[typeloop][0]=((fityield->Integral(1.115683+.497611,2.4))/(Bin1820Width));//K0
                                                    ALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][(typeloop-1)+(2*cutloop*Nsigmaflag)]=yieldarray[typeloop][0];
                                                }
                                                AlloftheYield->SetBinContent((typeloop),yieldarray[typeloop][0]);
                                                for (int i=1; i <= ((fitmax-fitmin)/Bin1820Width); i++) {
                                                    double xbin = a->GetBinCenter(i+1+((fitmin-1.6)/Bin1820Width));
                                                    double fval = backfit->Eval(xbin);
                                                    double diff = (a->GetBinContent(i+1+((fitmin-1.6)/Bin1820Width)) - fval);
                                                    htemp->Fill(xbin,diff);
                                                    htemp->SetBinError(i+1+((fitmin-1.6)/Bin1820Width),a->GetBinError(i+1+((fitmin-1.6)/Bin1820Width)));
                                                }
                                            }//end of subtraction
                                        }//old subtraction method
                                        else{//new subtraction method
                                            if (fitsubtract == 1) {
                                                if (UseSigmaRange == 1) {
                                                    if (UseVoigt == 1) {
                                                        if (MCflag == 1) {
                                                            fityield = new TF1("fityield",fAll,&AllFunctions::funVoigtMC,fitmin,fitmax,5);//No polynomial new
                                                        }
                                                        else{
                                                            fityield = new TF1("fityield",fAll,&AllFunctions::funVoigt,fitmin,fitmax,4);//No polynomial new
                                                        }
                                                    }
                                                    else{
                                                        fityield = new TF1("fityield",fAll,&AllFunctions::BWfunction,fitmin,fitmax,3);//No polynomial
                                                    }
                                                }
                                                else{
                                                    fityield = new TF1("fityield",fAll,&AllFunctions::BWfunction,fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),3);//No polynomial
                                                }
                                                
                                                if (versionloop == 0) {
                                                    backfit = new TF1("backfit",fAll,&AllFunctions::backline,fitmin,fitmax,backnumber);
                                                    //backfit = new TF1("backfit",backline,fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),backnumber);
                                                }
                                                else if(versionloop == 1){
                                                    backfit = new TF1("backfit",fAll,&AllFunctions::back2nd,fitmin,fitmax,backnumber);
                                                }
                                                else if(versionloop == 2){
                                                    backfit = new TF1("backfit",fAll,&AllFunctions::back3rd,fitmin,fitmax,backnumber);
                                                }
                                                else if(versionloop == 3){
                                                    backfit = new TF1("backfit",fAll,&AllFunctions::back4th,fitmin,fitmax,backnumber);
                                                }
                                                //htemp->Reset();
                                                for (int resetloop=1; resetloop <= htemp->GetNbinsX(); resetloop++) {
                                                    htemp->SetBinContent(resetloop,0);
                                                    htemp->SetBinError(resetloop,0);
                                                }
                                                for (int setloop=0; setloop < backnumber; setloop++) {
                                                    backfit->SetParameter(setloop,fitfunction->GetParameter(setloop));
                                                    backfit->SetParError(setloop,fitfunction->GetParError(setloop));
                                                }
                                                if (UseVoigt == 1) {
                                                    fityield->SetParameter(0,fitfunction->GetParameter(backnumber));
                                                    fityield->SetParError(0,fitfunction->GetParError(backnumber));
                                                    fityield->SetParameter(1,fitfunction->GetParameter(2+backnumber));
                                                    fityield->SetParError(1,fitfunction->GetParError(2+backnumber));
                                                    if (typeloop == 1) {
                                                        /*if (MCflag == 1) {
                                                            fityield->FixParameter(2,MCSigma[cb][pb-MinPT]);
                                                        }
                                                        else{*/
                                                            fityield->FixParameter(2,0.002);
                                                        //}
                                                    }
                                                    else{
                                                        fityield->SetParameter(2,0.002);
                                                    }
                                                    fityield->SetParameter(3,fitfunction->GetParameter(1+backnumber));
                                                    fityield->SetParError(3,fitfunction->GetParError(1+backnumber));
                                                    //slope
                                                    if (MCflag == 1) {
                                                        fityield->FixParameter(4,AEslope[cb][pb-MinPT]);
                                                    }
                                                    //slope
                                                }
                                                else{
                                                    for (int setloop2=0; setloop2 < (polynumber-backnumber); setloop2++) {
                                                        fityield->SetParameter(setloop2,fitfunction->GetParameter(setloop2+backnumber));
                                                        fityield->SetParError(setloop2,fitfunction->GetParError(setloop2+backnumber));
                                                    }//end of setloop2
                                                }
                                                
                                                if (UseSigmaRange == 1) {//use sigma range
                                                    if (Use1Type == 1) {//use 1 type and No poly
                                                        ALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][0]=(fityield->Integral(fitmin,fitmax)/(Bin1820Width));//KX
                                                        for (int i=1; i <= ((fitmax-fitmin)/(Bin1820Width)); i++) {//definate change
                                                            double xbin = a->GetBinCenter(i+1+((fitmin-1.6)/Bin1820Width));
                                                            double fval = backfit->Eval(xbin);
                                                            double diff = (a->GetBinContent(i+1+(((fitmin)-1.6)/Bin1820Width)) - fval);
                                                            htemp->Fill(xbin,diff);
                                                            htemp->SetBinError(i+1+(((fitmin)-1.6)/Bin1820Width),a->GetBinError(i+1+(((fitmin)-1.6)/Bin1820Width)));
                                                        }
                                                    }//use 1 type and No poly
                                                    else{//else
                                                        yieldarray[typeloop][0]=(fityield->Integral(fitmin,fitmax)/(Bin1820Width));//KX
                                                        ALLYield[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=yieldarray[typeloop][0];
                                                        
                                                        //AlloftheYield->SetBinContent((typeloop),yieldarray[typeloop][0]);
                                                        for (int i=1; i <= ((fitmax-fitmin)/(Bin1820Width)); i++) {//definate change
                                                            double xbin = a->GetBinCenter(i+1+((fitmin-1.6)/Bin1820Width));
                                                            double fval = backfit->Eval(xbin);
                                                            double diff = (a->GetBinContent(i+1+(((fitmin)-1.6)/Bin1820Width)) - fval);
                                                            htemp->Fill(xbin,diff);
                                                            htemp->SetBinError(i+1+(((fitmin)-1.6)/Bin1820Width),a->GetBinError(i+1+(((fitmin)-1.6)/Bin1820Width)));
                                                        }
                                                    }//else
                                                }//use sigma range
                                                else{
                                                    if (Use1Type == 1) {//use 1 type and No poly
                                                        ALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][0]=(fityield->Integral(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))/(Bin1820Width));//KX
                                                        for (int i=1; i <= ((fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))-(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))))/(Bin1820Width)); i++) {//definate change
                                                            double xbin = a->GetBinCenter(i+1+((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))-1.6)/Bin1820Width));
                                                            double fval = backfit->Eval(xbin);
                                                            double diff = (a->GetBinContent(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width)) - fval);
                                                            htemp->Fill(xbin,diff);
                                                            htemp->SetBinError(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width),a->GetBinError(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width)));
                                                        }
                                                    }//use 1 type and No poly
                                                    else{//else
                                                        yieldarray[typeloop][0]=(fityield->Integral(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))/(Bin1820Width));//KX
                                                        ALLYield[cb][pb-MinPT][jNR-NRBegin][0][0][(typeloop-1)+(2*cutloop*Nsigmaflag)]=yieldarray[typeloop][0];
                                                        
                                                        AlloftheYield->SetBinContent((typeloop),yieldarray[typeloop][0]);
                                                        for (int i=1; i <= ((fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))-(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))))/(Bin1820Width)); i++) {//definate change
                                                            double xbin = a->GetBinCenter(i+1+((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548))-1.6)/Bin1820Width));
                                                            double fval = backfit->Eval(xbin);
                                                            double diff = (a->GetBinContent(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width)) - fval);
                                                            htemp->Fill(xbin,diff);
                                                            htemp->SetBinError(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width),a->GetBinError(i+1+(((fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))-1.6)/Bin1820Width)));
                                                        }
                                                    }//else
                                                }
                                            }//end of subtraction
                                        }//new subtraction method
                                        
                                        c->cd();
                                        a->Draw();
                                        fitfunction->SetLineColor(1);
                                        fitfunction->Draw("SAME");
                                        c->Update();
                                        backfit->SetLineColor(3);
                                        backfit->Draw("SAME");
                                        c->Update();
                                        gStyle->SetOptStat("emr");//()
                                        c->Update();
                                        
                                        if (c->GetPrimitive("stats") != NULL) {
                                            ps1 = (TPaveStats*)c->GetPrimitive("stats");
                                            if (!ps1) {
                                                cout << "Error found *111111111111111111" << endl;
                                            }
                                            ps1->SetName("mystats");//mystats
                                            c->Modified();
                                            TList *listofLines2 = ps1->GetListOfLines();
                                            if (!listofLines2) {
                                                cout << "Error here *11" << endl;
                                            }
                                            tconst = ps1->GetLineWith("RMS");
                                            listofLines2->Remove(tconst);
                                            tconst = ps1->GetLineWith("Entries");
                                            listofLines2->Remove(tconst);
                                            tconst = ps1->GetLineWith("Mean");
                                            listofLines2->Remove(tconst);
                                            tconst = ps1->GetLineWith("Std Dev");
                                            listofLines2->Remove(tconst);
                                        }
                                        else{
                                            ps1 = (TPaveStats*)c->GetPrimitive("mystats");
                                            if (!ps1) {
                                                cout << "Error found *111111111111111111" << endl;
                                            }
                                            ps1->SetName("mystats");//mystats
                                            c->Modified();
                                            TList *listofLines2 = ps1->GetListOfLines();
                                            if (!listofLines2) {
                                                cout << "Error here *11" << endl;
                                            }
                                            tconst = ps1->GetLineWith("Integral");
                                            listofLines2->Remove(tconst);
                                            tconst = ps1->GetLineWith("Mean");
                                            listofLines2->Remove(tconst);
                                            tconst = ps1->GetLineWith("Width");
                                            listofLines2->Remove(tconst);
                                            tconst = ps1->GetLineWith("chi2/ndf");
                                            listofLines2->Remove(tconst);
                                        }
                                        a->SetStats(0);
                                        if (UseNoPoly == 1) {//use no poly
                                            //ps1->AddText(Form("Integral = %5.1f",yieldarray[typeloop][0]));
                                            ps1->AddText(Form("Mean = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+2),fitfunction->GetParError(backnumber+2)));
                                            ps1->AddText(Form("Width = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+1),fitfunction->GetParError(backnumber+1)));
                                            ps1->AddText(Form("chi2/ndf = %5.4f",(fitfunction->GetChisquare())/(fitfunction->GetNDF())));
                                            if (BellSwitch == 1) {
                                                ps1->SetX1NDC(0.6);
                                                ps1->SetY1NDC(0.6);
                                                ps1->SetY2NDC(0.9);
                                            }
                                        }//use no poly
                                        else{//use poly
                                            ps1->AddText(Form("Integral = %5.1f",yieldarray[typeloop][0]));
                                            ps1->AddText(Form("Mean = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+2),fitfunction->GetParError(backnumber+2)));
                                            ps1->AddText(Form("Width = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+1),fitfunction->GetParError(backnumber+1)));
                                            ps1->AddText(Form("chi2/ndf = %5.4f",(fitfunction->GetChisquare())/(fitfunction->GetNDF())));
                                            if (BellSwitch == 1) {
                                                ps1->SetX1NDC(0.6);
                                                ps1->SetY1NDC(0.6);
                                                ps1->SetY2NDC(0.9);
                                            }
                                        }//use poly
                                        
                                        c->Modified();//keep here
                                        if (hmtloop == 0) {//Normal
                                            if (cutloop == 0) {//no cut
                                                c->Print(Form("plotsCompleteV3/%s/Type%i/sub_mixunlike_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cb,pb,rangeloop,jNR,versionloop));
                                            }//no cut
                                            else{//cut
                                                c->Print(Form("plotsCompleteV3/%s/Type%i/cut%i_sub_mixunlike_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cutloop,cb,pb,rangeloop,jNR,versionloop));
                                            }//cut
                                        }//Normal
                                        else{//hmt
                                            if (cutloop == 0) {//no cut
                                                c->Print(Form("plotsCompleteV3/%s/Type%i/hmt_sub_mixunlike_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cb,pb,rangeloop,jNR,versionloop));
                                            }//no cut
                                            else{//cut
                                                c->Print(Form("plotsCompleteV3/%s/Type%i/hmtcut%i_sub_mixunlike_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cutloop,cb,pb,rangeloop,jNR,versionloop));
                                            }//cut
                                        }//hmt
                                        c->Clear();
                                        if (typeloop == 1) {
                                            if (rangeloop == 0) {
                                                if (versionloop == 1) {
                                                    if (((jNR == 0) && (particleloop == 0)) || ((jNR == 2) && (particleloop == 1))) {
                                                        if ((cb == 0) && (pb == MinPT)) {
                                                            cSubMix->Clear();
                                                            if (Use1PT == 1) {
                                                                cSubMix->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                                            }
                                                            else{
                                                                cSubMix->Divide(MaxCen,(MaxPT-MinPT));//4,()
                                                            }
                                                            cSubMix->Update();
                                                        }
                                                        cSubMix->cd((cb+1)+(MaxPT-MinPT)*(pb-MinPT));
                                                        a->Draw();
                                                        cSubMix->Update();
                                                        fitfunction->SetLineColor(1);
                                                        fitfunction->DrawCopy("SAME");//draw
                                                        cSubMix->Update();
                                                        backfit->SetLineColor(3);
                                                        backfit->DrawCopy("SAME");//draw
                                                        cSubMix->Update();
                                                        if ((cb == (MaxCen-1)) && (pb == (MaxPT-1))) {
                                                            if (hmtloop == 0) {//Normal
                                                                if (cutloop == 0) {//no cut
                                                                    cSubMix->Print(Form("plotsCompleteV3/%s/Type%i/Sub_Mix_allC_pt.pdf",partname,typeloop));
                                                                }//no cut
                                                                else{//cut
                                                                    cSubMix->Print(Form("plotsCompleteV3/%s/Type%i/cut%i_Sub_Mix_allC_pt.pdf",partname,typeloop,cutloop));
                                                                }//cut
                                                            }//Normal
                                                            else{//hmt
                                                                if (cutloop == 0) {//no cut
                                                                    cSubMix->Print(Form("plotsCompleteV3/%s/Type%i/hmt_Sub_Mix_allC_pt.pdf",partname,typeloop));
                                                                }//no cut
                                                                else{//cut
                                                                    cSubMix->Print(Form("plotsCompleteV3/%s/Type%i/hmtcut%i_Sub_Mix_allC_pt.pdf",partname,typeloop,cutloop));
                                                                }//cut
                                                            }//hmt
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        //if (versionloop == 1) {
                                        c0->Clear();
                                        //c0->Divide(2,2);
                                        //}
                                        c0->cd();
                                        //copy a for htemp2
                                        htemp2->Fit("fitfunction","Q","SAME",fitmin,fitmax);
                                        htemp2->Draw();
                                        c0->Update();
                                        /*if (versionloop == (versionmax-1)) {
                                         c0->Print(Form("plotsCompleteV3/%s/%s/Range%i/Type%i/together_sub_mixunlike%i_c%i_p%i.pdf",normal,partname,rangeloop,typeloop,jNR,cb,pb));
                                         }//complete Together plots*/
                                        c->cd();
                                        if (UseNoPoly == 1) {//use No poly
                                            for (int setloop=0; setloop < backnumber; setloop++) {
                                                fitline->ReleaseParameter(setloop);
                                                backfit->ReleaseParameter(setloop);
                                            }
                                        }//use no poly
                                        else{
                                            fitfunction->ReleaseParameter(backnumber+1);
                                            fitfunction->ReleaseParameter(backnumber+2);
                                            for (int setloop=0; setloop < backnumber; setloop++) {
                                                if (typeloop <= 7) {
                                                    fitfunction->ReleaseParameter(setloop);
                                                    if (typeloop == 7) {
                                                        fitline->ReleaseParameter(setloop);
                                                    }
                                                }
                                                backfit->ReleaseParameter(setloop);
                                            }
                                        }
                                        
                                        if (fitsubtract == 1) {
                                            if (UseNoPoly == 1) {//use No poly
                                                fityield->SetParLimits(0, 0.0, 500000.0);
                                                if (UseVoigt == 1) {
                                                    fityield->SetParLimits(3,0.01,0.06);
                                                    fityield->SetParLimits(1,1.8,1.85);
                                                    fityield->SetParameter(3,0.025);
                                                    fityield->SetParLimits(2,0.001,0.003);
                                                    if (typeloop == 1) {
                                                        /*if (MCflag == 1) {
                                                            fityield->FixParameter(2,MCSigma[cb][pb-MinPT]);
                                                        }
                                                        else{*/
                                                            fityield->FixParameter(2,0.002);
                                                        //}
                                                    }
                                                    else{
                                                        fityield->SetParameter(2,0.002);
                                                    }
                                                    fityield->SetParameter(1,1.82);
                                                    //slope
                                                    if (MCflag == 1) {
                                                        fityield->FixParameter(4,AEslope[cb][pb-MinPT]);
                                                    }
                                                    //slope
                                                }
                                                else{
                                                    fityield->SetParLimits(1,0.01,0.06);
                                                    fityield->SetParLimits(2,1.8,1.85);
                                                    fityield->SetParameter(1,0.025);
                                                    fityield->SetParameter(2,1.82);
                                                }
                                                
                                            }//use no poly
                                            else{
                                                fitfunction->SetParLimits(backnumber, 0.0, 500000.0);
                                                fitfunction->SetParLimits(backnumber+1,0.01,0.06);
                                                fitfunction->SetParLimits(backnumber+2,1.8,1.85);
                                                fitfunction->SetParameter(backnumber+1,0.025);
                                                fitfunction->SetParameter(backnumber+2,1.82);
                                            }
                                        }
                                        
                                        if (UseNoPoly == 1) {
                                            //new linear section
                                            for (int removeloop=0; removeloop < backnumber; removeloop++) {
                                                fitfunction->FixParameter(removeloop,0);
                                            }
                                            //new linear section
                                            if (UseSigmaRange == 1) {//use sigma range
                                                htemp->Fit("fityield","Q","SAME",fitmin2,fitmax2);
                                                for (int setloop=0; setloop < backnumber; setloop++) {
                                                    backfit->FixParameter(setloop,0);
                                                }
                                                backfit->SetLineColor(kBlue);
                                                htemp->Fit("backfit","QR+","SAME",fitmin,fitmax);
                                            }//use sigma range
                                            else{
                                                htemp->Fit("fityield","Q","SAME",fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)));
                                                for (int setloop=0; setloop < backnumber; setloop++) {
                                                    backfit->FixParameter(setloop,0);
                                                }
                                                backfit->SetLineColor(kBlue);
                                                htemp->Fit("backfit","QR+","SAME",fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)));
                                            }
                                        }
                                        else{
                                            htemp->Fit("fitfunction","Q","SAME",fitmin,fitmax);
                                            for (int setloop=0; setloop < backnumber; setloop++) {
                                                backfit->FixParameter(setloop,fitfunction->GetParameter(setloop));
                                            }
                                            backfit->SetLineColor(kBlue);
                                            htemp->Fit("backfit","QR+","SAME",fitmin,fitmax);
                                        }
                                        
                                        c->cd();
                                        gStyle->SetOptStat("emr");//()
                                        htemp->Draw();
                                        c->Update();
                                        
                                        if (c->GetPrimitive("stats") != NULL) {
                                            ps1 = (TPaveStats*)c->GetPrimitive("stats");
                                            if (!ps1) {
                                                cout << "Error found 111111111111111111" << endl;
                                            }
                                            ps1->SetName("mystats");//mystats
                                            c->Modified();
                                            TList *listofLines1 = ps1->GetListOfLines();
                                            if (!listofLines1) {
                                                cout << "Error here 11" << endl;
                                            }
                                            tconst = ps1->GetLineWith("RMS");
                                            listofLines1->Remove(tconst);
                                            tconst = ps1->GetLineWith("Entries");
                                            listofLines1->Remove(tconst);
                                            tconst = ps1->GetLineWith("Mean");
                                            listofLines1->Remove(tconst);
                                            tconst = ps1->GetLineWith("Std Dev");
                                            listofLines1->Remove(tconst);
                                        }
                                        else{
                                            ps1 = (TPaveStats*)c->GetPrimitive("mystats");
                                            if (!ps1) {
                                                cout << "Error found 111111111111111111" << endl;
                                            }
                                            ps1->SetName("mystats");//mystats
                                            c->Modified();
                                            TList *listofLines1 = ps1->GetListOfLines();
                                            if (!listofLines1) {
                                                cout << "Error here 11" << endl;
                                            }
                                            tconst = ps1->GetLineWith("Integral");
                                            listofLines1->Remove(tconst);
                                            tconst = ps1->GetLineWith("Mean");
                                            listofLines1->Remove(tconst);
                                            tconst = ps1->GetLineWith("Width");
                                            listofLines1->Remove(tconst);
                                            if (UseVoigt == 1) {
                                                tconst = ps1->GetLineWith("Sigma");
                                                listofLines1->Remove(tconst);
                                            }
                                            tconst = ps1->GetLineWith("chi2/ndf");
                                            listofLines1->Remove(tconst);
                                            tconst = ps1->GetLineWith("Yield Bin");
                                            listofLines1->Remove(tconst);
                                        }
                                        
                                        axis = htemp->GetXaxis();
                                        if (UseNoPoly == 1) {
                                            bmax = axis->FindBin(1.82);
                                            if (UseSigmaRange == 1) {
                                                integral = fityield->Integral(fitmin2,fitmax2)/axis->GetBinWidth(bmax);
                                                integralbin=0;//
                                                for (i=1; i <= htemp->GetNbinsX(); i++) {
                                                    binvalue=htemp->GetBinContent(i);
                                                    if (binvalue == 0.0) {
                                                        //cout << "i " << i << " = Null bin value " << binvalue << " bin center " << htemp->GetBinCenter(i) << endl;
                                                    }
                                                    else{
                                                        //cout << "i " << i << " not Null bin value " << binvalue << " bin center " << htemp->GetBinCenter(i) << endl;
                                                        if (((i*Bin1820Width+1.60) >= fitmin2) && ((i*Bin1820Width+1.60) <= fitmax2)) {//range2
                                                            integralbin+=binvalue;
                                                        }
                                                    }
                                                }
                                                BinALLYield[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=integralbin;
                                                //insert bin counting method
                                            }
                                            else{
                                                integral = fityield->Integral(fitfunction->GetParameter(backnumber+2)-5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)),fitfunction->GetParameter(backnumber+2)+5.0*(fitfunction->GetParameter(backnumber+1)/(2.3548)))/axis->GetBinWidth(bmax);
                                            }
                                            
                                            htemp->SetStats(0);
                                            ps1->AddText(Form("Integral = %5.1f",integral));
                                            if (UseVoigt == 1) {
                                                ps1->AddText(Form("Mean = %5.4f +- %5.5f",fityield->GetParameter(1),fityield->GetParError(1)));
                                                ps1->AddText(Form("Width = %5.4f +- %5.5f",fityield->GetParameter(3),fityield->GetParError(3)));
                                                ps1->AddText(Form("Sigma = %5.4f +- %5.5f",fityield->GetParameter(2),fityield->GetParError(2)));
                                            }
                                            else{
                                                ps1->AddText(Form("Mean = %5.4f +- %5.5f",fityield->GetParameter(2),fityield->GetParError(2)));
                                                ps1->AddText(Form("Width = %5.4f +- %5.5f",fityield->GetParameter(1),fityield->GetParError(1)));
                                            }
                                            ps1->AddText(Form("chi2/ndf = %5.4f",(fityield->GetChisquare())/(fityield->GetNDF())));
                                            ps1->AddText(Form("Yield Bin = %5.1f",integralbin));
                                            if (BellSwitch == 1) {
                                                ps1->SetX1NDC(0.6);
                                                ps1->SetY1NDC(0.6);
                                                ps1->SetY2NDC(0.9);
                                            }
                                        }
                                        else{
                                            int bmin = axis->FindBin(1.77);
                                            int bmax = axis->FindBin(1.87);
                                            double integral50 = fitfunction->Integral(1.77,1.87)/axis->GetBinWidth(bmax);
                                            if (flag50 == 1) {
                                                double integral = fitfunction->Integral(1.77,1.87)/axis->GetBinWidth(bmax);
                                            }
                                            
                                            double integralback = backfit->Integral(1.77,1.87)/axis->GetBinWidth(bmax);
                                            htemp->SetStats(0);
                                            if (flag50 == 1) {
                                                ps1->AddText(Form("Integral = %5.1f",integral));
                                            }
                                            ps1->AddText(Form("Mean = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+2),fitfunction->GetParError(backnumber+2)));
                                            ps1->AddText(Form("Width = %5.4f +- %5.5f",fitfunction->GetParameter(backnumber+1),fitfunction->GetParError(backnumber+1)));
                                            ps1->AddText(Form("chi2/ndf = %5.4f",(fitfunction->GetChisquare())/(fitfunction->GetNDF())));
                                            if (BellSwitch == 1) {
                                                ps1->SetX1NDC(0.6);
                                                ps1->SetY1NDC(0.6);
                                                ps1->SetY2NDC(0.9);
                                            }
                                        }
                                        
                                        if (UseLinear == 1) {//use linear
                                            if (UseVoigt == 1) {
                                                ALLMean[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fityield->GetParameter(1);
                                                ALLMeanError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fityield->GetParError(1);
                                                ALLWidth[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fityield->GetParameter(3);
                                                ALLWidthError[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=fityield->GetParError(3);
                                            }
                                            else{
                                                ALLMean[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParameter(2);
                                                ALLMeanError[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParError(2);
                                                ALLWidth[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParameter(1);
                                                ALLWidthError[cb][pb-MinPT][jNR-NRBegin][0][0][0]=fityield->GetParError(1);
                                            }
                                            ALLChi2[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=(fityield->GetChisquare()/fityield->GetNDF());
                                            ALLYield[cb][pb-MinPT][jNR-NRBegin][rangeloop][versionloop][(typeloop-1)+(2*cutloop*Nsigmaflag)]=integral;
                                            //}
                                        }//use linear
                                        c->Modified();
                                        if (hmtloop == 0) {//Normal
                                            if (cutloop == 0) {//no cut
                                                c->Print(Form("plotsCompleteV3/%s/Type%i/linear_sub_mixunlike_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cb,pb,rangeloop,jNR,versionloop));
                                            }//no cut
                                            else{//cut
                                                c->Print(Form("plotsCompleteV3/%s/Type%i/cut%i_linear_sub_mixunlike_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cutloop,cb,pb,rangeloop,jNR,versionloop));
                                            }//cut
                                        }//Normal
                                        else{//hmt
                                            if (cutloop == 0) {//no cut
                                                c->Print(Form("plotsCompleteV3/%s/Type%i/hmt_linear_sub_mixunlike_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cb,pb,rangeloop,jNR,versionloop));
                                            }//no cut
                                            else{//cut
                                                c->Print(Form("plotsCompleteV3/%s/Type%i/hmtcut%i_linear_sub_mixunlike_c%i_p%i_R%i_NR%i_v%i.pdf",partname,typeloop,cutloop,cb,pb,rangeloop,jNR,versionloop));
                                            }//cut
                                        }//hmt
                                        if (typeloop == 1) {
                                            if (rangeloop == 0) {
                                                if (versionloop == 1) {
                                                    if (((jNR == 0) && (particleloop == 0)) || ((jNR == 2) && (particleloop == 1))) {
                                                        if ((cb == 0) && (pb == MinPT)) {
                                                            cLinearMix->Clear();
                                                            if (Use1PT == 1) {
                                                                cLinearMix->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                                            }
                                                            else{
                                                                cLinearMix->Divide(MaxCen,(MaxPT-MinPT));//4,()
                                                            }
                                                            cLinearMix->Update();
                                                        }
                                                        cLinearMix->cd((cb+1)+(MaxPT-MinPT)*(pb-MinPT));
                                                        htemp->Draw();
                                                        cLinearMix->Update();
                                                        if ((cb == (MaxCen-1)) && (pb == (MaxPT-1))) {
                                                            if (hmtloop == 0) {//Normal
                                                                if (cutloop == 0) {//no cut
                                                                    cLinearMix->Print(Form("plotsCompleteV3/%s/Type%i/Linear_Mix_allC_pt.pdf",partname,typeloop));
                                                                }//no cut
                                                                else{//cut
                                                                    cLinearMix->Print(Form("plotsCompleteV3/%s/Type%i/cut%i_Linear_Mix_allC_pt.pdf",partname,typeloop,cutloop));
                                                                }//cut
                                                            }//Normal
                                                            else{//hmt
                                                                if (cutloop == 0) {//no cut
                                                                    cLinearMix->Print(Form("plotsCompleteV3/%s/Type%i/hmt_Linear_Mix_allC_pt.pdf",partname,typeloop));
                                                                }//no cut
                                                                else{//cut
                                                                    cLinearMix->Print(Form("plotsCompleteV3/%s/Type%i/hmtcut%i_Linear_Mix_allC_pt.pdf",partname,typeloop,cutloop));
                                                                }//cut
                                                            }//hmt
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        fitfunction->ReleaseParameter(backnumber+1);
                                        fitfunction->ReleaseParameter(backnumber+2);
                                        //new linear section
                                        if (UseNoPoly == 1) {
                                            for (removeloop=0; removeloop <= backnumber; removeloop++) {
                                                fitfunction->ReleaseParameter(backnumber);
                                            }
                                            for (setloop2=0; setloop2 < 3; setloop2++) {
                                                fityield->ReleaseParameter(setloop2);
                                            }
                                            if (UseVoigt == 1) {
                                                fityield->ReleaseParameter(3);
                                                //slope
                                                fityield->ReleaseParameter(4);
                                                //slope
                                            }
                                        }
                                        //new linear section
                                        //if (versionloop == 1) {
                                        c1->Clear();
                                        //c1->Divide(1,3);
                                        //}
                                        c1->cd();
                                        htemp->DrawCopy();
                                        c1->Update();
                                        if (Use1Type != 1) {//not 1 type
                                            if (typeloop == 8) {//loop for types
                                                TMultiGraph *Tmgintegral = new TMultiGraph();//type multigraphs
                                                TMultiGraph *Tmgmean = new TMultiGraph();
                                                TMultiGraph *Tmgwidth = new TMultiGraph();
                                                for (int j=1; j <= 4; j++) {
                                                    if (j==1) {
                                                        systemmean[0]=0;
                                                        systemmeanerror[0]=0;
                                                        systemwidth[0]=0;
                                                        systemwidtherror[0]=0;
                                                    }
                                                    systemmean[0]+=systemmean[j];
                                                    systemwidth[0]+=systemwidth[j];
                                                    
                                                }//loop to assign value and calculate values
                                                systemmean[0]=(systemmean[0]/4.0);
                                                systemwidth[0]=(systemwidth[0]/4.0);
                                                
                                                for (j=1; j <= 4; j++) {
                                                    systemmeanerror[0]+=pow(systemmean[j]-systemmean[0],2);
                                                    systemwidtherror[0]+=pow(systemwidth[j]-systemwidth[0],2);
                                                }
                                                systemmeanerror[0]=sqrt(systemmeanerror[0]/4.0);
                                                systemwidtherror[0]=sqrt(systemwidtherror[0]/4.0);
                                                
                                                //cout << "systematic mean " << systemmean[0] << " error " << systemmeanerror[0] << " width " << systemwidth[0] << " error " << systemwidtherror[0] << endl;
                                                TGraphErrors* gSmean = new TGraphErrors(5,systemxarray,systemmean,systemxarrayerror,systemmeanerror);
                                                TGraphErrors* gSwidth = new TGraphErrors(5,systemxarray,systemwidth,systemxarrayerror,systemwidtherror);
                                                cA->cd();
                                                gSmean->SetTitle(Form("Means %s %s C %i Pt %i",normal,partname,cb,pb));
                                                gSmean->SetLineColor(2);
                                                gSmean->Draw();
                                                cA->Update();
                                                //hmt and cutloop needed
                                                cA->Print(Form("plotsCompleteV3/%s/systemMeans_mixunlike_c%i_p%i_v%i.pdf",partname,cb,pb,versionloop));
                                                cA->Clear();
                                                gSwidth->SetTitle(Form("Widths %s %s C %i Pt %i",normal,partname,cb,pb));
                                                gSwidth->SetLineColor(2);
                                                gSwidth->Draw();
                                                cA->Update();
                                                cA->Print(Form("plotsCompleteV3/%s/systemWidths_mixunlike_c%i_p%i_v%i.pdf",partname,cb,pb,versionloop));
                                                cA->Clear();
                                                
                                                for (int k=0; k <= 8; k++) {
                                                    yieldplottype[k]=yieldarray[k][0];
                                                }
                                                Double_t N[9];
                                                for (int k=0; k <= 8; k++) {
                                                    N[k]=k;
                                                }
                                                TGraph* gyield = new TGraph(9,N,yieldplottype);
                                                cA->cd();
                                                gyield->SetTitle(Form("Yield All Types %s %s C %i Pt %i",normal,partname,cb,pb));
                                                gyield->SetLineColor(2);
                                                gyield->Draw();
                                                cA->Update();
                                                cA->Print(Form("plotsCompleteV3/%s/AllTypeYield_sub_c%i_p%i_v%i.pdf",partname,cb,pb,versionloop));
                                                cA->Clear();
                                                //if ((versionloop == 0) && (rangeloop == 0) && (jNR == NRstart)) {
                                                TMultiGraph *mgyieldtype = new TMultiGraph();
                                                //}//multigrpah for yeilds
                                                gyield->SetLineColor(versionloop+2);
                                                gyield->SetLineStyle(rangeloop);
                                                gyield->SetMarkerStyle(20);
                                                mgyieldtype->Add(gyield);
                                                //if ((versionloop == (versionmax-1)) && (rangeloop == 0) && ((jNR == (nNR-1)) || (Use1NR == 1))) {
                                                cB->cd();
                                                cB->Clear();
                                                cB->Update();
                                                mgyieldtype->SetTitle(Form("All-Yield All-Types %s %s C %i Pt %i",normal,partname,cb,pb));
                                                mgyieldtype->DrawClone("ALP");
                                                cB->Update();
                                                cB->Print(Form("plotsCompleteV3/%s/all_yield_all_type_c%i_p%i.pdf",partname,cb,pb));
                                                cB->Clear();
                                                //}
                                            }//end of if typeloop=8
                                        }//not 1 type

                                        c->cd();
                                        for (int setloop=0; setloop < backnumber; setloop++) {
                                            backfit->ReleaseParameter(setloop);
                                        }
                                    }
                                    gStyle->SetOptStat(0);
                                    fitfunction->ReleaseParameter(backnumber+1);
                                    fitfunction->ReleaseParameter(backnumber+2);
                                }//end of type loop
                            }//end of versionloop
                        }//end of jNR loop and mixunlike
                        
                        
                        if(scheme!=2) continue;
                        
                        
                        for(jNR=0;jNR<nNR;jNR++){//plot sum-mixlike

                            sprintf(s0,"mass_sumQmixlike%i_c%i_p%i",jNR,cb,pb);
                            a=(TH1D*) f->Get(s0);
                            if(!a){cerr<<"missing histogram "<<s0<<endl; continue;}
                            
                            a->Rebin(rb);
                            a->GetXaxis()->SetRange(0,50);
                            a->SetTitle(Form("%s, %1.1f-%1.1f %% , %1.1f<#it{p}_{T}<%1.1f GeV/c",DC1,cbmin[cb],cbmax[cb],ptmin[pb],ptmax[pb]));
                            a->SetXTitle("Invariant Mass [GeV/#it{c}]");
                            a->GetXaxis()->SetTitleOffset(1.3);
                            a->SetYTitle("Counts/(0.01 GeV/c^{2})");
                            a->GetYaxis()->SetTitleOffset(1.3);
                            a->SetMarkerStyle(1);
                            a->SetMarkerColor(TColor::GetColor(color));
                            a->SetLineColor(TColor::GetColor(color));
                            /*if (fitflag == 1) {
                             double par[5]=0;
                             for (int z = 0; z < 5; z++) {
                             par[z]=0;
                             }
                             //if (flagXi1820 == 1) {
                             fitfunction = new TF1("fitfunction",function,fitmin,fitmax,5);
                             fitline = new TF1("fitline",fline,fitmin,fitmax,2);//new
                             fitfunction->SetParLimits(2, 0.0, 50000.0);
                             fitfunction->ReleaseParameter(3);
                             fitfunction->ReleaseParameter(4);
                             fitfunction->SetParLimits(3,0.01,0.06);
                             fitfunction->SetParLimits(4,1.8,1.85);
                             //fitfunction->SetParLimits(3,-2000.,8000.);
                             //fitfunction->SetParLimits(4,-5000.,5000.);
                             //bool reject;
                             //reject = kTRUE;
                             a->Fit("fitline","Q","SAME",fitmin,fitmax);
                             //reject = kFALSE;
                             fitfunction->SetParameter(3,0.025);
                             fitfunction->SetParameter(4,1.82);
                             
                             a->Fit("fitfunction","Q","SAME",fitmin,fitmax);
                             //for (z=0; z < 5; z++) {
                             //    cout << "parmeter " << z << " " << fitfunction->GetParameter(z) << endl;
                             //}
                             //cout << "chi2 Xi1820 " << fitfunction->GetChisquare() << endl;
                             }*/
                            

                            c->cd();
                            a->Draw();
                            if (jDC==Klkx) {
                                if (hmtflag == 1) {
                                    c->Print(Form("plotsCompleteV3/hmt/LambdaKX/sub_mixlike%i_c%i_p%i.gif",jNR,cb,pb));
                                    //c->Print(Form("plotsCompleteV3/hmt/LambdaKX/sub_mixlike%i_c%i_p%i.pdf",jNR,cb,pb));
                                }
                                else
                                {
                                    c->Print(Form("plotsCompleteV3/Normal/LambdaKX/sub_mixlike%i_c%i_p%i.gif",jNR,cb,pb));
                                    //c->Print(Form("plotsCompleteV3/Normal/LambdaKX/sub_mixlike%i_c%i_p%i.eps",jNR,cb,pb));
                                }
                            }
                            else if(jDC==Klk0){
                                if (hmtflag == 1) {
                                    c->Print(Form("plotsCompleteV3/hmt/LambdaK0/sub_mixlike%i_c%i_p%i.gif",jNR,cb,pb));
                                    //c->Print(Form("plotsCompleteV3/hmt/LambdaK0/sub_mixlike%i_c%i_p%i.eps",jNR,cb,pb));
                                }
                                else
                                {
                                    c->Print(Form("plotsCompleteV3/Normal/LambdaK0/sub_mixlike%i_c%i_p%i.gif",jNR,cb,pb));
                                    //c->Print(Form("plotsCompleteV3/Normal/LambdaK0/sub_mixlike%i_c%i_p%i.eps",jNR,cb,pb));
                                }
                            }
                        }
                    }//end of range loop
                }//end cen and pt loop
                
                //return here
                if (cutloop == cutflag) {//cutloop == cutflag
                    cout << "begin major calculations" << endl;
                    if (particleloop <= 1) {
                        //begin ALLYield plots
                        Double_t Xforplot[MaxCen];
                        Double_t Yforplot[MaxCen];
                        Double_t XforplotError[MaxCen];
                        Double_t YforplotError[MaxCen];
                        
                        Double_t Xforplot3[MaxCen];
                        Double_t Yforplot3[MaxCen];
                        Double_t XforplotError3[MaxCen];
                        Double_t YforplotError3[MaxCen];
                        Double_t XforplotError3SYS[MaxCen];
                        Double_t YforplotError3SYS[MaxCen];
                        Double_t YforplotP[2][4][MaxCen];//[gm mixed][pt][cen]
                        Double_t YforplotPError[2][4][MaxCen];//[gm mixed][pt][cen]
                        Double_t YforplotPErrorSYS[2][4][MaxCen];//[gm mixed][pt][cen]
                        Double_t MforplotP[2][4][MaxCen];//[gm mixed][pt][cen]
                        Double_t MforplotPError[2][4][MaxCen];//[gm mixed][pt][cen]
                        Double_t MforplotPErrorSYS[2][4][MaxCen];//[gm mixed][pt][cen]
                        Double_t WforplotP[2][4][MaxCen];//[gm mixed][pt][cen]
                        Double_t WforplotPError[2][4][MaxCen];//[gm mixed][pt][cen]
                        Double_t WforplotPErrorSYS[2][4][MaxCen];//[gm mixed][pt][cen]
                        if (Use1PT == 1) {
                            for (int i=0; i < MaxCen; i++) {
                                Xforplot[i]=i;
                                XforplotError[i]=0;
                            }
                        }
                        
                        if (particleloop == 0) {
                            if (Use1PT != 1) {
                                TMultiGraph *P0MultiAllYield = new TMultiGraph();
                                TMultiGraph *P1MultiAllYield = new TMultiGraph();
                                TMultiGraph *P2MultiAllYield = new TMultiGraph();
                                P0MultiAllYield->SetTitle(Form("All-Yield vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[0],ptmax[0]));
                                P1MultiAllYield->SetTitle(Form("All-Yield vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[1],ptmax[1]));
                                P2MultiAllYield->SetTitle(Form("All-Yield vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[2],ptmax[2]));
                                TMultiGraph *P0MultiAllMean = new TMultiGraph();
                                TMultiGraph *P1MultiAllMean = new TMultiGraph();
                                TMultiGraph *P2MultiAllMean = new TMultiGraph();
                                P0MultiAllMean->SetTitle(Form("All-Mean vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[0],ptmax[0]));
                                P1MultiAllMean->SetTitle(Form("All-Mean vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[1],ptmax[1]));
                                P2MultiAllMean->SetTitle(Form("All-Mean vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[2],ptmax[2]));
                                TMultiGraph *P0MultiAllWidth = new TMultiGraph();
                                TMultiGraph *P1MultiAllWidth = new TMultiGraph();
                                TMultiGraph *P2MultiAllWidth = new TMultiGraph();
                                P0MultiAllWidth->SetTitle(Form("All-Width vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[0],ptmax[0]));
                                P1MultiAllWidth->SetTitle(Form("All-Width vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[1],ptmax[1]));
                                P2MultiAllWidth->SetTitle(Form("All-Width vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[2],ptmax[2]));
                            }
                            P3MultiAllYield = new TMultiGraph();
                            P3MultiAllYield->SetTitle(Form("All-Yield vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[3],ptmax[3]));
                            P3MultiAllMean = new TMultiGraph();
                            P3MultiAllMean->SetTitle(Form("All-Mean vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[3],ptmax[3]));
                            P3MultiAllWidth = new TMultiGraph();
                            P3MultiAllWidth->SetTitle(Form("All-Width vs Cent. for %1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[3],ptmax[3]));
                        }
                        
                        //major calculation section
                        for (int n=0; n < MaxCen; n++) {//cen
                            if ((hmtloop == 1) && (n >= 1)) {
                                continue;
                            }
                            
                            if (n == 0) {
                                P3MultiAllYield = new TMultiGraph();
                                P3MultiAllYield->SetTitle(Form("p_{T} Spectra; p_{T} (GeV/c); d^{2}N/dydp_{T} (Nraw/(N_{ev}*dy*dp_{T}*A*E)) "));
                                
                                P3MultiAllMean = new TMultiGraph();
                                if (BellSwitch == 1) {
                                    P3MultiAllMean->SetTitle(Form("Mean for %s ; p_{T} (GeV/c) ; Mean (GeV/c^{2}) ",partnameComplex));
                                }
                                else{
                                    P3MultiAllMean->SetTitle(Form("Mean for %s ; p_{T} (GeV/c) ; Mean (GeV/c^{2}) ",partname));
                                }
                                
                                P3MultiAllWidth = new TMultiGraph();
                                if (BellSwitch == 1) {
                                    P3MultiAllWidth->SetTitle(Form("Width (FWHM) for %s; p_{T} (GeV/c) ; Width (GeV/c^{2})",partnameComplex));
                                    P3MultiAllWidth->SetMinimum(0.000);
                                }
                                else{
                                    P3MultiAllWidth->SetTitle(Form("Width (FWHM) for X#rightarrow#Lambda(#bar{#Lambda})K; p_{T} (GeV/c) ; Width (GeV/c^{2})"));
                                }
                                
                                //P3MultiAllWidth->GetXaxis()->SetTitle(Form("p_{T} (GeV/c)"));
                                //P3MultiAllWidth->GetYaxis()->SetTitle(Form("Width (GeV/c^{2})"));
                                if (particleloop == 0) {
                                    MultiYieldAverage = new TMultiGraph();
                                    MultiYieldAverage->SetTitle(Form("p_{t} Spectra Average; p_{T} (GeV/c); d^{2}N/dydp_{T} (Nraw/(N_{ev}*dy*dp_{T}*A*E)) "));
                                    //MultiYieldAverage->GetXaxis()->SetTitle(Form("p_{T} (GeV/c)"));
                                    //MultiYieldAverage->GetYaxis()->SetTitle(Form("d^{2}N/dydp_{T} (Nraw/(N_{ev}*dy*dp_{T}*A*E))"));
                                    MultiMeanAverage = new TMultiGraph();
                                    MultiMeanAverage->SetTitle(Form("Mean Average for %s; p_{T} (GeV/c); Mean (GeV/c^{2}) ",partname));
                                    //MultiMeanAverage->GetXaxis()->SetTitle(Form("p_{T} (GeV/c)"));
                                    //MultiMeanAverage->GetYaxis()->SetTitle(Form("Mean (GeV/c^{2})"));
                                    MultiWidthAverage = new TMultiGraph();
                                    MultiWidthAverage->SetTitle(Form("Width Average for %s; p_{T} (GeV/c); Width (GeV/c^{2}) ",partname));
                                    //MultiWidthAverage->GetXaxis()->SetTitle(Form("p_{T} (GeV/c)"));
                                    //MultiWidthAverage->GetYaxis()->SetTitle(Form("Width (GeV/c^{2})"));
                                }
                            }
                            
                            //for (int m=0; m < (NREnd-NRBegin); m++) {//jNR
                            
                            /*//not really used
                             MultiAllYield = new TMultiGraph();
                             for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type l=0; l < 8
                             for (int k=0; k < 1; k++) {//version
                             for (int j=0; j < 1; j++) {//range
                             for (int m=0; m < (NREnd-NRBegin); m++) {//jNR
                             for (int i=0; i < (MaxPT-MinPT); i++) {//pt
                             Xforplot[i]=i;
                             XforplotError[i]=0;
                             Yforplot[i]=0;
                             Yforplot[i]=ALLYield[n][i][m][j][k][l];
                             }
                             cC->cd();
                             TGraph* gAllYield = new TGraph(MaxPT-MinPT,Xforplot,Yforplot);
                             gAllYield->SetLineColor(l+1);//type
                             gAllYield->SetLineStyle(k);//version
                             gAllYield->SetMarkerStyle(20+j);//range
                             MultiAllYield->Add(gAllYield);
                             }
                             }
                             }
                             }
                             cC->cd();
                             cC->Clear();
                             cC->Update();
                             MultiAllYield->SetTitle(Form("All-Yeild All Type vs. pt for Cen %i %s %s",n,partname,normal));
                             MultiAllYield->DrawClone("ALP");
                             cC->Update();
                             if (hmtloop == 0) {//Normal
                             //if (cutloop == 0) {//no cut
                             cC->Print(Form("plotsCompleteV3/%s/AllYield_NR%i_c%i.pdf",partname,m,n));
                             //}//no cut
                             //else{//cut
                             //    cC->Print(Form("plotsCompleteV3/%s/cut1_AllYield_NR%i_c%i.pdf",partname,m,n));
                             //}//cut
                             }//Normal
                             else{//hmt
                             //if (cutloop == 0) {//no cut
                             cC->Print(Form("plotsCompleteV3/%s/hmt_AllYield_NR%i_c%i.pdf",partname,m,n));
                             //}//no cut
                             //else{//cut
                             //    cC->Print(Form("plotsCompleteV3/%s/hmtcut1_AllYield_NR%i_c%i.pdf",partname,m,n));
                             //}//cut
                             }//hmt
                             cC->Clear();
                             *///not really used
                            TH1D* SystematicAllYield = new TH1D(Form("SystematicYield_%s_%s",partname,othernormal),"",4,0,4);//6,0,6
                            TH1D* SystematicAllMean = new TH1D(Form("SystematicMean_%s_%s",partname,othernormal),"",4,0,4);
                            TH1D* SystematicAllWidth = new TH1D(Form("SystematicWidth_%s_%s",partname,othernormal),"",4,0,4);
                            TH1D* StatisticalAllMean = new TH1D(Form("StatisticalMean_%s_%s",partname,othernormal),"",4,0,4);
                            TH1D* StatisticalAllWidth = new TH1D(Form("StatisticalWidth_%s_%s",partname,othernormal),"",4,0,4);
                            calculateSYSmean=0;
                            calculateSYSmeanerror=0;
                            calculateSTAmeanerror=0;
                            calculateSYSyield=0;
                            calculateSYSyielderror=0;
                            calculateSYSyieldBin=0;
                            calculateSYSyielderrorBin=0;
                            calculateSYSwidth=0;
                            calculateSYSwidtherror=0;
                            calculateSTAwidtherror=0;
                            numberofmean=0;
                            numberofwidth=0;
                            numberofyield=0;
                            for (int i=0; i < (MaxPT-MinPT); i++) {//pt
                                numberofyield=0;
                                calculateSYSyield=0;
                                calculateSYSyielderror=0;
                                calculateSYSyieldBin=0;
                                calculateSYSyielderrorBin=0;
                                numberofmean=0;
                                calculateSYSmean=0;
                                calculateSYSmeanerror=0;
                                calculateSTAmeanerror=0;
                                numberofwidth=0;
                                calculateSYSwidth=0;
                                calculateSYSwidtherror=0;
                                calculateSTAwidtherror=0;
                                for (int j=0; j < 3; j++) {//range j=0; j < 3
                                    if ((BellSwitch == 1) && (j != 0)) {
                                        continue;
                                    }
                                    
                                    for (int k=1; k < 4; k++) {//version k=0; k < 3
                                        for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type l=0; l < 8
                                            for (int m=0; m < (NREnd-NRBegin); m++) {//Normalization range (new)

                                                if (BellSwitch == 1) {
                                                    if (particleloop == 0) {
                                                        m=2;
                                                    }
                                                    else if(particleloop == 1){
                                                        m=0;
                                                    }
                                                }
                                                
                                                if (onechange == 1) {
                                                    //default for particle 0 (LAmbda+KX) l=0, k=1, j=0, m=2
                                                    //default for particle 1 (Lambda+K0) l=0, k=1, j=0, m=0
                                                    if ((l%2 == 0) && (l == 0) && (k == 1) && (j == 0)/* && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))*/) {
                                                        //this is changing NR
                                                    }
                                                    else if((l%2 == 0) && (l == 0) && (k == 1) && /*(j == 0) &&*/ (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))) {
                                                        //this is changing range
                                                    }
                                                    else if((l%2 == 0) && (l == 0) && /*(k == 1) &&*/ (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                        //this is changing fit range
                                                    }
                                                    else if((l == 1) && (k == 1) && (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                        //this is changing sigma
                                                    }
                                                    else if(/*(l%2 == 0) &&*/ (l == 0) && (k == 1) && (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                        //this is changing NSigma
                                                    }
                                                    else{
                                                        continue;//this is not needed
                                                    }
                                                }
                                                if ((BellSwitch == 1) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))) {
                                                    if (particleloop == 0) {
                                                        m=0;//change to 0
                                                    }
                                                    else if(particleloop == 1){
                                                        m=0;
                                                    }
                                                }
                                                //cout << "n " << n << " i " << i << " m " << m << " j " << j << " k " << k << " l " << l  << " allyield " << ALLYield[n][i][m][j][k][l] << endl;
                                                calculateSYSyield += ALLYield[n][i][m][j][k][l];
                                                calculateSYSyieldBin += BinALLYield[n][i][m][j][k][l];
                                                numberofyield++;
                                                numberofmean++;
                                                calculateSYSmean += ALLMean[n][i][m][j][k][l];
                                                calculateSTAmeanerror += pow(ALLMeanError[n][i][m][j][k][l],2);
                                                numberofwidth++;
                                                calculateSYSwidth += ALLWidth[n][i][m][j][k][l];
                                                calculateSTAwidtherror += pow(ALLWidthError[n][i][m][j][k][l],2);
                                            }
                                        }
                                    }
                                }
                                if (numberofyield == 0) {
                                    numberofyield=1;
                                    numberofmean=1;
                                    numberofwidth=1;
                                }
                                calculateSYSyieldBin=calculateSYSyieldBin/numberofyield;
                                
                                calculateSYSyield=calculateSYSyield/numberofyield;
                                SystematicAllYield->SetBinContent(i+1,calculateSYSyield);
                                calculateSYSmean=calculateSYSmean/numberofmean;
                                SystematicAllMean->SetBinContent(i+1,calculateSYSmean);
                                calculateSYSwidth=calculateSYSwidth/numberofwidth;
                                SystematicAllWidth->SetBinContent(i+1,calculateSYSwidth);
                                calculateSTAmeanerror=sqrt(calculateSTAmeanerror)/numberofmean;
                                StatisticalAllMean->SetBinContent(i+1,calculateSYSmean);
                                StatisticalAllMean->SetBinError(i+1,calculateSTAmeanerror);
                                calculateSTAwidtherror=sqrt(calculateSTAwidtherror)/numberofwidth;
                                StatisticalAllWidth->SetBinContent(i+1,calculateSYSwidth);
                                StatisticalAllWidth->SetBinError(i+1,calculateSTAwidtherror);
                                
                                //P#Multi
                                if (Use1PT == 1) {//use 1 pt
                                    if (Use1Type == 1) {
                                        //YforplotP[1][3][n]=((calculateSYSyield+calculateSYSyieldBin)/2.0)/((cbmax[n]-cbmin[n]));
                                        YforplotP[1][3][n]=((calculateSYSyield+calculateSYSyieldBin)/2.0);
                                    }
                                    else{
                                        //YforplotP[1][3][n]=calculateSYSyield/(cbmax[n]-cbmin[n]);
                                        YforplotP[1][3][n]=calculateSYSyield;
                                    }
                                    //YforplotP[1][3][n]=((calculateSYSyield+calculateSYSyieldBin)/2.0)/((cbmax[n]-cbmin[n]));//recent
                                    YforplotP[1][3][n]=((calculateSYSyield+calculateSYSyieldBin)/2.0);//recent
                                    //YforplotP3[n]=calculateSYSyield;
                                    MforplotP[1][3][n]=calculateSYSmean;
                                    WforplotP[1][3][n]=calculateSYSwidth;
                                }//use 1 pt
                                else{
                                    YforplotP[1][i][n]=((calculateSYSyield+calculateSYSyieldBin)/2.0);//recent
                                    //YforplotP3[n]=calculateSYSyield;
                                    MforplotP[1][i][n]=calculateSYSmean;
                                    WforplotP[1][i][n]=calculateSYSwidth;
                                    /*
                                     YforplotP[1][i][n]=calculateSYSyield;
                                     MforplotP[1][i][n]=calculateSYSmean;
                                     WforplotP[1][i][n]=calculateSYSwidth;
                                     if(i == 3){
                                     YforplotP[1][i][n]=calculateSYSyield/(cbmax[n]-cbmin[n]);
                                     }*/
                                }
                                //P#Multi
                                for (int j=0; j < 3; j++) {//range j=0; j < 3
                                    if ((BellSwitch == 1) && (j != 0)) {
                                        continue;
                                    }
                                    
                                    for (int k=1; k < 4; k++) {//version k=0; k < 3
                                        for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type l=0; l < 8
                                            for (int m=0; m < (NREnd-NRBegin); m++) {//Normalization range (new)
                                                if (BellSwitch == 1) {
                                                    if (particleloop == 0) {
                                                        m=2;
                                                    }
                                                    else if(particleloop == 1){
                                                        m=0;
                                                    }
                                                }
                                                if (onechange == 1) {
                                                    //default for particle 0 (LAmbda+KX) l=0, k=1, j=0, m=2
                                                    //default for particle 1 (Lambda+K0) l=0, k=1, j=0, m=0
                                                    if ((l%2 == 0) && (l == 0) && (k == 1) && (j == 0)/* && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))*/) {
                                                        //this is changing NR
                                                    }
                                                    else if((l%2 == 0) && (l == 0) && (k == 1) && /*(j == 0) &&*/ (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                        //this is changing range
                                                    }
                                                    else if((l%2 == 0) && (l == 0) && /*(k == 1) &&*/ (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                        //this is changing fit range
                                                    }
                                                    else if((l == 1) && (k == 1) && (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                        //this is changing sigma
                                                    }
                                                    else if(/*(l%2 == 0) &&*/ (l == 0) && (k == 1) && (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                        //this is changing NSigma
                                                    }
                                                    else{
                                                        continue;//this is not needed
                                                    }
                                                }
                                                
                                                if (BellSwitch == 1) {
                                                    if (particleloop == 0) {
                                                        m=0;
                                                    }
                                                    else if(particleloop == 1){
                                                        m=0;
                                                    }
                                                }
                                                calculateSYSyielderror += pow(ALLYield[n][i][m][j][k][l] - calculateSYSyield,2);
                                                calculateSYSmeanerror += pow(ALLMean[n][i][m][j][k][l] - calculateSYSmean,2);
                                                calculateSYSwidtherror += pow(ALLWidth[n][i][m][j][k][l] - calculateSYSwidth,2);
                                            }
                                        }
                                    }
                                }
                                calculateSYSyielderror=sqrt(pow(calculateSYSyield-calculateSYSyieldBin,2)+pow(sqrt(calculateSYSyielderror/numberofyield),2));//sqrt(calculateSYSyielderror/numberofyield)
                                SystematicAllYield->SetBinError(i+1,calculateSYSyielderror);
                                calculateSYSmeanerror=sqrt(calculateSYSmeanerror/numberofmean);
                                SystematicAllMean->SetBinError(i+1,calculateSYSmeanerror);
                                calculateSYSwidtherror=sqrt(calculateSYSwidtherror/numberofwidth);
                                SystematicAllWidth->SetBinError(i+1,calculateSYSwidtherror);
                                if (Use1PT == 1) {//use 1 pt
                                    if (Use1Type == 1) {
                                        //YforplotPError[1][3][n]=((calculateSYSyield-calculateSYSyieldBin)/2.0)/((cbmax[n]-cbmin[n]));
                                        YforplotPError[1][3][n]=TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0));
                                    }
                                    else{
                                        //YforplotPError[1][3][n]=calculateSYSyielderror/(cbmax[n]-cbmin[n]);
                                        YforplotPError[1][3][n]=calculateSYSyielderror;
                                    }
                                    //YforplotPError[1][3][n]=((calculateSYSyield-calculateSYSyieldBin)/2.0)/((cbmax[n]-cbmin[n]));//recent
                                    YforplotPError[1][3][n]=TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0));//previous
                                    //YforplotPError[1][3][n]=sqrt(pow(TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0)),2)+pow(calculateSYSyielderror,2));//current
                                    //YforplotPError3[n]=calculateSYSyielderror;
                                    MforplotPError[1][3][n]=sqrt(pow(calculateSYSmeanerror,2) + pow(calculateSTAmeanerror,2));
                                    WforplotPError[1][3][n]=sqrt(pow(calculateSYSwidtherror,2) + pow(calculateSTAwidtherror,2));
                                    //slight modification
                                    YforplotPErrorSYS[1][3][n]=0;//systematic error is shown above
                                    
                                    MforplotPError[1][3][n]=calculateSTAmeanerror;
                                    MforplotPErrorSYS[1][3][n]=calculateSYSmeanerror;
                                    
                                    WforplotPError[1][3][n]=calculateSTAwidtherror;
                                    WforplotPErrorSYS[1][3][n]=calculateSYSwidtherror;
                                    //slight modification to be removed
                                }//use 1 pt
                                else{
                                    if (Use1Type == 1) {
                                        YforplotPError[1][i][n]=TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0));
                                    }
                                    else{
                                        YforplotPError[1][i][n]=calculateSYSyielderror;
                                    }
                                    YforplotPError[1][i][n]=TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0));//previous
                                    //YforplotPError[1][i][n]=sqrt(pow(TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0)),2)+pow(calculateSYSyielderror,2));//current
                                    YforplotPErrorSYS[1][i][n]=0;//systematic error is shown above
                                    //MforplotPError[1][i][n]=sqrt(pow(calculateSYSmeanerror,2) + pow(calculateSTAmeanerror,2));
                                    MforplotPError[1][i][n]=calculateSTAmeanerror;
                                    MforplotPErrorSYS[1][i][n]=calculateSYSmeanerror;
                                    //WforplotPError[1][i][n]=sqrt(pow(calculateSYSwidtherror,2) + pow(calculateSTAwidtherror,2));
                                    WforplotPError[1][i][n]=calculateSTAwidtherror;
                                    WforplotPErrorSYS[1][i][n]=calculateSYSwidtherror;
                                    /*
                                     YforplotPError[1][i][n]=calculateSYSyielderror;
                                     MforplotPError[1][i][n]=sqrt(pow(calculateSYSmeanerror,2) + pow(calculateSTAmeanerror,2));
                                     WforplotPError[1][i][n]=sqrt(pow(calculateSYSwidtherror,2) + pow(calculateSTAwidtherror,2));
                                     if(i == 3){
                                     //YforplotPError[1][3][n]=calculateSYSyielderror/(cbmax[n]-cbmin[n]);
                                     YforplotPError[1][3][n]=calculateSYSyielderror;
                                     }*/
                                }
                                
                                //systematic mean and width invariant mass plot
                                
                                //if (i == 0) {
                                //    cout << "**ALL systematic caluclations for " << partname << endl;
                                //}
                                //cout << "cen " << n << " pt bin " << i << " mean " << calculateSYSmean << " +- " << calculateSYSmeanerror << " +- Sta " << calculateSTAmeanerror << " width " << calculateSYSwidth << " +- " << calculateSYSwidtherror << " +- Sta " << calculateSTAwidtherror << " yield " << calculateSYSyield << " +- " << calculateSYSyielderror << endl;
                                SystematicAllYield->GetXaxis()->SetBinLabel(i+1,Form("%1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[i+MinPT],ptmax[i+MinPT]));
                                SystematicAllYield->GetYaxis()->SetTitle("Yield");
                                SystematicAllMean->GetXaxis()->SetBinLabel(i+1,Form("%1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[i+MinPT],ptmax[i+MinPT]));
                                SystematicAllMean->GetYaxis()->SetTitle("Mean (GeV/c^{2})");
                                SystematicAllWidth->GetXaxis()->SetBinLabel(i+1,Form("%1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[i+MinPT],ptmax[i+MinPT]));
                                SystematicAllWidth->GetYaxis()->SetTitle("Width (GeV/c^{2})");
                                if (writetxt == 1) {//begin systematic txt file
                                    if ((particleloop == 0) && (n == 0) && (i == 0) /*&& (cutloop == 0)*/ && (hmtloop == 0)) {
                                        //outputFileSTAT.open("plotsCompleteV3/StatisticalDataValues.txt", std::ios::app);
                                        //fopen_s(&Cutoftxt,outputFileSTAT,"wt");
                                        //more work needed
                                        Cutoftxt = fopen("plotsCompleteV3/StatisticalDataValues.txt","w");
                                        if (MCflag == 1) {
                                            if (Cutoftxt!=NULL) {
                                                fprintf(Cutoftxt,"particle | cut/hmt | cb       | pt        | yield/A*E +- sys error | A*E +- error | A*E slope +- error | mean +- sys error +- sta error | width +- sys error +- sta error\n");
                                            }
                                            else{
                                                cout << "error with statistical data values" << endl;
                                            }
                                        }
                                        else{
                                            if (Cutoftxt!=NULL) {
                                                fprintf(Cutoftxt,"particle | cut/hmt | cb       | pt        | yield +- sys error | mean +- sys error +- sta error | width +- sys error +- sta error\n");
                                            }
                                            else{
                                                cout << "error with statistical data values" << endl;
                                            }
                                        }
                                        
                                        Percentage = fopen("plotsCompleteV3/Percentage.txt","w");
                                        if (MCflag == 1) {
                                            if (Percentage!=NULL) {
                                                fprintf(Percentage,"particle | cut/hmt | cb       | pt        | yield/A*E +- sys error | A*E +- error | A*E slope +- error | mean +- sys error +- sta error | width +- sys error +- sta error\n");
                                                fprintf(Percentage,"Nor.Ran.   |range    |version | type    | yield per| mean per| width per|\n");
                                            }
                                            else{
                                                cout << "error with statistical data values" << endl;
                                            }
                                        }
                                        else{
                                            if (Percentage!=NULL) {
                                                fprintf(Percentage,"particle | cut/hmt | cb       | pt        | yield +- sys error | mean +- sys error +- sta error | width +- sys error +- sta error\n");
                                                fprintf(Percentage,"Nor.Ran.   |range    |version | type    | yield per| mean per| width per|\n");
                                                
                                            }
                                            else{
                                                cout << "error with statistical data values" << endl;
                                            }
                                        }
                                        //fprintf(Cutoftxt,"particle | cb | pt | yield +- sys error | mean +- sys error +- sta error | width +- sys error +- sta error\n");
                                    }
                                    
                                    if (MCflag == 1) {
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            fprintf(Cutoftxt,"%s |  Normal |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                            /*}//no cut
                                             else{//cut
                                             fprintf(Cutoftxt,"%s | cut1     |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                             }//cut*/
                                        }//Normal
                                        else{//hmt
                                            //if (cutloop == 0) {//no cut
                                            fprintf(Cutoftxt,"%s | hmt      |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                            /*}//no cut
                                             else{//cut
                                             fprintf(Cutoftxt,"%s | hmtcut1  |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                             }//cut*/
                                        }//hmt
                                    }
                                    else{
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            fprintf(Cutoftxt,"%s |  Normal |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                            /*}//no cut
                                             else{//cut
                                             fprintf(Cutoftxt,"%s | cut1     |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                             }//cut*/
                                        }//Normal
                                        else{//hmt
                                            //if (cutloop == 0) {//no cut
                                            fprintf(Cutoftxt,"%s | hmt      |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                            /*}//no cut
                                             else{//cut
                                             fprintf(Cutoftxt,"%s | hmtcut1  |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                             }//cut*/
                                        }//hmt
                                    }//MC
                                    
                                    //percentage
                                    if (MCflag == 1) {
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            fprintf(Percentage,"%s |  Normal |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                            /*}//no cut
                                             else{//cut
                                             fprintf(Percentage,"%s | cut1     |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                             }//cut*/
                                        }//Normal
                                        else{//hmt
                                            //if (cutloop == 0) {//no cut
                                            fprintf(Percentage,"%s | hmt      |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                            /*}//no cut
                                             else{//cut
                                             fprintf(Percentage,"%s | hmtcut1  |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                             }//cut*/
                                        }//hmt
                                    }
                                    else{
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            fprintf(Percentage,"%s |  Normal |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                            /*}//no cut
                                             else{//cut
                                             fprintf(Percentage,"%s | cut1     |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                             }//cut*/
                                        }//Normal
                                        else{//hmt
                                            //if (cutloop == 0) {//no cut
                                            fprintf(Percentage,"%s | hmt      |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                            /*}//no cut
                                             else{//cut
                                             fprintf(Percentage,"%s | hmtcut1  |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                             }//cut*/
                                        }//hmt
                                    }//MC
                                    //percentage calculation section
                                    tempYield=0;
                                    tempMean=0;
                                    tempWidth=0;
                                    for (int j=0; j < 3; j++) {//range j=0; j < 3
                                        if ((BellSwitch == 1) && (j != 0)) {
                                            continue;
                                        }
                                        for (int k=1; k < 4; k++) {//version k=0; k < 3
                                            for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type l=0; l < 8
                                                for (int m=0; m < (NREnd-NRBegin); m++) {//Normalization range (new)
                                                    if (BellSwitch == 1) {
                                                        if (particleloop == 0) {
                                                            m=2;
                                                        }
                                                        else if(particleloop == 1){
                                                            m=0;
                                                        }
                                                    }
                                                    
                                                    if (onechange == 1) {
                                                        //default for particle 0 (LAmbda+KX) l=0, k=1, j=0, m=2
                                                        //default for particle 1 (Lambda+K0) l=0, k=1, j=0, m=0
                                                        if ((l%2 == 0) && (l == 0) && (k == 1) && (j == 0)/* && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))*/) {
                                                            //this is changing NR
                                                        }
                                                        else if((l%2 == 0) && (l == 0) && (k == 1) && /*(j == 0) &&*/ (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                            //this is changing range
                                                        }
                                                        else if((l%2 == 0) && (l == 0) && /*(k == 1) &&*/ (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                            //this is changing fit range
                                                        }
                                                        else if((l == 1) && (k == 1) && (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                            //this is changing sigma
                                                        }
                                                        else if(/*(l%2 == 0) &&*/ (l == 0) && (k == 1) && (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                            //this is changing NSigma
                                                        }
                                                        else{
                                                            continue;//this is not needed
                                                        }
                                                    }
                                                    
                                                    if (BellSwitch == 1) {
                                                        if (particleloop == 0) {
                                                            m=0;
                                                        }
                                                        else if(particleloop == 1){
                                                            m=0;
                                                        }
                                                    }
                                                    
                                                    tempYield+=TMath::Abs(ALLYield[n][i][m][j][k][l]-calculateSYSyield);
                                                    tempMean+=TMath::Abs(ALLMean[n][i][m][j][k][l]-calculateSYSmean);
                                                    tempWidth+=TMath::Abs(ALLWidth[n][i][m][j][k][l]-calculateSYSwidth);
                                                }
                                            }
                                        }
                                    }
                                    //newer
                                    tempYield+=TMath::Abs(calculateSYSyield-calculateSYSyieldBin);//newer
                                    //newer
                                    for (int j=0; j < 3; j++) {//range j=0; j < 3
                                        if ((BellSwitch == 1) && (j != 0)) {
                                            continue;
                                        }
                                        for (int k=1; k < 4; k++) {//version k=0; k < 3
                                            for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type l=0; l < 8
                                                for (int m=0; m < (NREnd-NRBegin); m++) {//Normalization range (new)
                                                    if (BellSwitch == 1) {
                                                        if (particleloop == 0) {
                                                            m=2;
                                                        }
                                                        else if(particleloop == 1){
                                                            m=0;
                                                        }
                                                    }
                                                    
                                                    if (onechange == 1) {
                                                        //default for particle 0 (LAmbda+KX) l=0, k=1, j=0, m=2
                                                        //default for particle 1 (Lambda+K0) l=0, k=1, j=0, m=0
                                                        if ((l%2 == 0) && (l == 0) && (k == 1) && (j == 0)/* && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))*/) {
                                                            //this is changing NR
                                                        }
                                                        else if((l%2 == 0) && (l == 0) && (k == 1) && /*(j == 0) &&*/ (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                            //this is changing range
                                                        }
                                                        else if((l%2 == 0) && (l == 0) && /*(k == 1) &&*/ (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                            //this is changing fit range
                                                        }
                                                        else if((l == 1) && (k == 1) && (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                            //this is changing sigma
                                                        }
                                                        else if(/*(l%2 == 0) &&*/ (l == 0) && (k == 1) && (j == 0) && (((particleloop == 0) && (m == 2)) || ((particleloop == 1) && (m == 0)))){
                                                            //this is changing NSigma
                                                        }
                                                        else{
                                                            continue;//this is not needed
                                                        }
                                                    }
                                                    
                                                    if (BellSwitch == 1) {
                                                        if (particleloop == 0) {
                                                            m=0;
                                                        }
                                                        else if(particleloop == 1){
                                                            m=0;
                                                        }
                                                    }
                                                    
                                                    //fprintf(Percentage,"range %i |version %i|type %i   | %2.4f | %2.4f | %2.4f | %1.4f | %1.4f |\n",j,k,l,(TMath::Abs((ALLYield[n][i][m][j][k][l]-calculateSYSyield)/sqrt(numberofyield)/calculateSYSyielderror)*100),(TMath::Abs((ALLMean[n][i][m][j][k][l]-calculateSYSmean)/sqrt(numberofmean)/calculateSYSmeanerror)*100),(TMath::Abs((ALLWidth[n][i][m][j][k][l]-calculateSYSwidth)/sqrt(numberofwidth)/calculateSYSwidtherror)*100),ALLMean[n][i][m][j][k][l],ALLWidth[n][i][m][j][k][l]);
                                                    //fprintf(Percentage,"range %i |version %i|type %i   | %2.4f | %2.4f | %2.4f | %1.4f | %1.4f |\n",j,k,l,(TMath::Abs((ALLYield[n][i][m][j][k][l]-calculateSYSyield)/(calculateSYSyielderror*sqrt(numberofyield)))*100),(TMath::Abs((ALLMean[n][i][m][j][k][l]-calculateSYSmean)/(calculateSYSmeanerror*sqrt(numberofmean)))*100),(TMath::Abs((ALLWidth[n][i][m][j][k][l]-calculateSYSwidth)/(calculateSYSwidtherror*sqrt(numberofwidth)))*100),ALLMean[n][i][m][j][k][l],ALLWidth[n][i][m][j][k][l]);
                                                    fprintf(Percentage,"Nor.Ran. %i |range %i |version %i|type %i   | %2.4f | %2.4f | %2.4f | %5.1f | %1.4f | %1.4f |",m,j,k,l,(TMath::Abs((ALLYield[n][i][m][j][k][l]-calculateSYSyield)/(tempYield))*100),(TMath::Abs((ALLMean[n][i][m][j][k][l]-calculateSYSmean)/(tempMean))*100),(TMath::Abs((ALLWidth[n][i][m][j][k][l]-calculateSYSwidth)/(tempWidth))*100),ALLYield[n][i][m][j][k][l]/AEarray[n][i],ALLMean[n][i][m][j][k][l],ALLWidth[n][i][m][j][k][l]);
                                                    if ((TMath::Abs((ALLYield[n][i][m][j][k][l]-calculateSYSyield)/(tempYield))*100) >= 10.0) {//high yield error
                                                        fprintf(Percentage," *Y ");
                                                    }
                                                    if ((TMath::Abs((ALLMean[n][i][m][j][k][l]-calculateSYSmean)/(tempMean))*100) >= 10.0) {//high Mean error
                                                        fprintf(Percentage," *M ");
                                                    }
                                                    if ((TMath::Abs((ALLWidth[n][i][m][j][k][l]-calculateSYSwidth)/(tempWidth))*100) >= 10.0) {//high width error
                                                        fprintf(Percentage," *W ");
                                                    }
                                                    
                                                    if (l == 1) {//only needed once
                                                        fprintf(Percentage," Int-Bin error | %1.4f |", (TMath::Abs(calculateSYSyield-calculateSYSyieldBin)/(tempYield))*100);
                                                    }//only needed once
                                                    fprintf(Percentage,"\n");
                                                }
                                            }
                                        }
                                    }
                                    //percentage calculation section
                                    //percentage
                                    
                                }//end systematic text file
                            }
                            //
                            if (Use1PT == 1) {//use 1 pt
                                if (particleloop == 0) {//temparry remove lambda + k0 section
                                    if (n == (MaxCen-1)) {
                                        
                                        if ((cbmin[0] == 0) && (cbmax[0] == 100)) {
                                            for (int loop3=0; loop3 < (MaxCen-1); loop3++) {
                                                Xforplot3[loop3]=((cbmin[loop3+1]+cbmax[loop3+1])/2.0);
                                                if (BellSwitch == 1) {
                                                    XforplotError3[loop3]=0;
                                                }
                                                else{
                                                    XforplotError3[loop3]=((cbmax[loop3+1]-cbmin[loop3+1])/2.0);
                                                }
                                            }
                                            
                                            for (int y=0; y < (MaxCen-1); y++) {
                                                Yforplot3[y]=YforplotP[1][3][y+1];//yield Vs. centrality
                                                YforplotError3[y]=YforplotPError[1][3][y+1];//yield Vs. centrality
                                            }
                                            gAllYieldP = new TGraphErrors((MaxCen-1),Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                            
                                        }
                                        else{
                                            for (int loop3=0; loop3 < MaxCen; loop3++) {
                                                Xforplot3[loop3]=((cbmin[loop3]+cbmax[loop3])/2.0);
                                                if (BellSwitch == 1) {
                                                    XforplotError3[loop3]=0;
                                                }
                                                else{
                                                    XforplotError3[loop3]=((cbmax[loop3]-cbmin[loop3])/2.0);
                                                }
                                            }
                                            
                                            for (int y=0; y < MaxCen; y++) {
                                                Yforplot3[y]=YforplotP[1][3][y];//yield Vs. centrality
                                                YforplotError3[y]=YforplotPError[1][3][y];//yield Vs. centrality
                                            }
                                            gAllYieldP = new TGraphErrors(MaxCen,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        
                                        if (particleloop == 0) {
                                            gAllYieldP->SetLineColor(1);//type
                                        }
                                        else if (particleloop == 1){
                                            gAllYieldP->SetLineColor(3);//type
                                        }
                                        P3MultiAllYield->Add(gAllYieldP);
                                        
                                        if ((cbmin[0] == 0) && (cbmax[0] == 100)) {
                                            for (int y=0; y < (MaxCen-1); y++) {
                                                Yforplot3[y]=MforplotP[1][3][y+1];//yield Vs. centrality
                                                YforplotError3[y]=MforplotPError[1][3][y+1];//yield Vs. centrality
                                            }
                                            gAllMeanP = new TGraphErrors((MaxCen-1),Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        else{
                                            for (int y=0; y < MaxCen; y++) {
                                                Yforplot3[y]=MforplotP[1][3][y];//yield Vs. centrality
                                                YforplotError3[y]=MforplotPError[1][3][y];//yield Vs. centrality
                                            }
                                            gAllMeanP = new TGraphErrors(MaxCen,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        
                                        cD->cd();
                                        cD->Update();
                                        gAllMeanP->DrawClone("AP");
                                        cD->Modified();
                                        cD->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cD->Print(Form("plotsCompleteV3/TestMean_Pt.pdf"));
                                            /*}//no cut
                                             else{//cut
                                             cD->Print(Form("plotsCompleteV3/cut1_TestMean_Pt.pdf"));
                                             }//cut*/
                                        }//Normal
                                        cD->Clear();
                                        if (particleloop == 0) {
                                            gAllMeanP->SetLineColor(1);//type
                                        }
                                        else if (particleloop == 1){
                                            gAllMeanP->SetLineColor(3);//type
                                        }
                                        P3MultiAllMean->Add(gAllMeanP);
                                        if (particleloop == 0) {
                                            cD->cd();
                                            cD->Update();
                                            P3MultiAllMean->DrawClone("AP");
                                            //P3MultiAllMean->GetXaxis()->SetTitle(Form("p_{T} (GeV/c)"));
                                            //P3MultiAllMean->GetYaxis()->SetTitle(Form("Mean (GeV/c^{2})"));
                                            cD->Modified();
                                            cD->Update();
                                            if (hmtloop == 0) {//Normal
                                                //if (cutloop == 0) {//no cut
                                                cD->Print(Form("plotsCompleteV3/TestMean2_Pt.pdf"));
                                                /*}//no cut
                                                 else{//cut
                                                 cD->Print(Form("plotsCompleteV3/cut1_TestMean2_Pt.pdf"));
                                                 }//cut*/
                                            }//Normal
                                            cD->Clear();
                                        }
                                        
                                        if ((cbmin[0] == 0) && (cbmax[0] == 100)) {
                                            for (int y=0; y < (MaxCen-1); y++) {
                                                Yforplot3[y]=WforplotP[1][3][y+1];//width Vs. centrality
                                                YforplotError3[y]=WforplotPError[1][3][y+1];//width Vs. centrality
                                            }
                                            gAllWidthP = new TGraphErrors((MaxCen-1),Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        else{
                                            for (int y=0; y < MaxCen; y++) {
                                                Yforplot3[y]=WforplotP[1][3][y];//width Vs. centrality
                                                YforplotError3[y]=WforplotPError[1][3][y];//width Vs. centrality
                                            }
                                            gAllWidthP = new TGraphErrors(MaxCen,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        
                                        if (particleloop == 0) {
                                            gAllWidthP->SetLineColor(1);//type
                                        }
                                        else if (particleloop == 1){
                                            gAllWidthP->SetLineColor(3);//type
                                        }
                                        P3MultiAllWidth->Add(gAllWidthP);
                                    }
                                }//tempariary remove kos section
                                
                                if (particleloop == 1) {//ending
                                    if (n == (MaxCen-1)) {
                                        cC->cd();
                                        cC->Update();
                                        //P3MultiAllYield->Fit("pol1");
                                        P3MultiAllYield->Draw("AP");
                                        //P3MultiAllYield->Fit("pol2","F");
                                        //fpol = P3MultiAllYield->GetFunction("pol2");
                                        //fpol->SetLineWidth(1);
                                        //fpol->DrawClone("SAME");
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/AllYield_Pt.pdf"));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_AllYield_Pt.pdf"));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        cC->Update();
                                        
                                        cC->cd();
                                        cC->Update();
                                        P3MultiAllMean->DrawClone("AP");
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/AllMean_Pt.pdf"));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_AllMean_Pt.pdf"));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        cC->Update();
                                        
                                        cC->cd();
                                        cC->Update();
                                        P3MultiAllWidth->DrawClone("AP");
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/AllWidth_Pt.pdf"));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_AllWidth_Pt.pdf"));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        cC->Update();
                                    }
                                    
                                }//ending
                            }//use 1 pt
                            else{//not 1 pt
                                
                                for (int loop3=0; loop3 < 4; loop3++) {
                                    Xforplot3[loop3]=((ptmin[loop3]+ptmax[loop3])/2.0);
                                    if (BellSwitch == 1) {
                                        XforplotError3[loop3]=0;
                                        XforplotError3SYS[loop3]=0.0;//fixed for now
                                    }
                                    else{
                                        XforplotError3[loop3]=((ptmax[loop3]-ptmin[loop3])/2.0);
                                        XforplotError3SYS[loop3]=0.5;//fixed for now
                                    }
                                }
                                
                                for (int y=0; y < 4; y++) {
                                    Yforplot3[y]=(YforplotP[1][y][n]/AEarray[n][y])/(NumberOfEvents*dy*(ptmax[y]-ptmin[y]));//yield Vs. pt
                                    YforplotError3[y]=(Yforplot3[y])*(sqrt(pow((YforplotPError[1][y][n]/YforplotP[1][y][n]),2)+pow((AEarrayerror[n][y]/AEarray[n][y]),2)))/(NumberOfEvents*dy*(ptmax[y]-ptmin[y]));//yield Vs. pt
                                }
                                gAllYieldP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                
                                if (BellSwitch == 1) {
                                    for (int y=0; y < 3; y++) {
                                        Yforplot3[y]=(YforplotP[1][y][n]/AEarray[n][y])/(NumberOfEvents*dy*(ptmax[y]-ptmin[y]));//yield Vs. pt
                                        YforplotError3[y]=(Yforplot3[y])*(sqrt(pow((YforplotPError[1][y][n]/YforplotP[1][y][n]),2)+pow((AEarrayerror[n][y]/AEarray[n][y]),2)))/(NumberOfEvents*dy*(ptmax[y]-ptmin[y]));//yield Vs. pt
                                    }
                                    gAllYieldP = new TGraphErrors(3,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                }
                                
                                if (particleloop == 0) {
                                    gAllYieldP->SetLineStyle(1);//type
                                }
                                else if (particleloop == 1){
                                    gAllYieldP->SetLineStyle(3);//type
                                }
                                gAllYieldP->SetLineColor(1+n);
                                
                                //new to remove
                                if (((n < 2) && (BellSwitch == 1)) || (BellSwitch == 0)) {//remove 10 and 30
                                    if (only1flag == 1) {
                                        if (particleloop == 1) {
                                            P3MultiAllYield->Add(gAllYieldP);
                                        }
                                    }
                                    else{
                                        P3MultiAllYield->Add(gAllYieldP);
                                    }
                                    
                                    if (only1flag == 1) {
                                        if (n == 0) {
                                            legYield= new TLegend(0.7,0.7,0.9,0.9);//(0.55,0.55,0.75,0.75)
                                        }
                                        
                                        if (particleloop == 1) {
                                            if ((n == 0) && (cbmax[0] == 100)) {
                                                sprintf(AEword,"%s (Min.Bias)",partnameComplex);//mixed
                                            }
                                            else{
                                                sprintf(AEword,"%s (%1.1f-%1.1f %%)",partnameComplex,cbmin[n],cbmax[n]);//mixed
                                            }
                                            legYield->AddEntry(gAllYieldP,AEword,"l");
                                        }
                                    }
                                    else{
                                        if (n == 0) {
                                            legYield= new TLegend(0.7,0.7,0.9,0.9);
                                        }
                                        if ((n == 0) && (cbmax[0] == 100)) {
                                            sprintf(AEword,"%s (Min.Bias)",partnameComplex);//mixed
                                        }
                                        else{
                                            sprintf(AEword,"%s (%1.1f-%1.1f %%)",partnameComplex,cbmin[n],cbmax[n]);//mixed
                                        }
                                        legYield->AddEntry(gAllYieldP,AEword,"l");
                                    }
                                }
                                //new to remove
                                
                                
                                for (int y=0; y < 4; y++) {
                                    Yforplot3[y]=MforplotP[1][y][n];//mean Vs. pt
                                    YforplotError3[y]=MforplotPError[1][y][n];//mean Vs. pt
                                    YforplotError3SYS[y]=MforplotPErrorSYS[1][y][n];
                                }
                                gAllMeanP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                
                                if (BellSwitch == 1) {
                                    for (int y=0; y < 3; y++) {
                                        Yforplot3[y]=MforplotP[1][y][n];//width Vs. pt
                                        YforplotError3[y]=MforplotPError[1][y][n];//width Vs. pt
                                        YforplotError3SYS[y]=MforplotPErrorSYS[1][y][n];
                                    }
                                    gAllMeanP = new TGraphErrors(3,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                }
                                
                                if (particleloop == 0) {
                                    gAllMeanP->SetLineStyle(1);//type
                                }
                                else if (particleloop == 1){
                                    gAllMeanP->SetLineStyle(3);//type
                                }
                                gAllMeanP->SetLineColor(1+n);
                                
                                //new to remove
                                if (((n < 2) && (BellSwitch == 1)) || (BellSwitch == 0)) {//remove 10 and 30
                                    if (only1flag == 1) {
                                        if (particleloop == 1) {
                                            P3MultiAllMean->Add(gAllMeanP);
                                        }
                                    }
                                    else{
                                        P3MultiAllMean->Add(gAllMeanP);
                                    }
                                    
                                    //SYS
                                    gAllMeanP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3SYS,YforplotError3SYS);
                                    
                                    if (particleloop == 0) {
                                        gAllMeanP->SetLineStyle(1);//type
                                    }
                                    else if (particleloop == 1){
                                        gAllMeanP->SetLineStyle(3);//type
                                    }
                                    gAllMeanP->SetLineColor(1+n);
                                    //gAllMeanP->SetFillStyle(3001);
                                    ///gAllMeanP->SetFillColorAlpha(1+n,0.10);
                                    gAllMeanP->SetFillStyle(1);
                                    
                                    if (only1flag == 1) {
                                        if (particleloop == 1) {
                                            if (BellSwitch == 1) {
                                                P3MultiAllMean->Add(gAllMeanP,"[]");
                                                gStyle->SetEndErrorSize(4);
                                            }
                                            else{
                                                P3MultiAllMean->Add(gAllMeanP,"5");
                                            }
                                        }
                                    }
                                    else{
                                        if (BellSwitch == 1) {
                                            P3MultiAllMean->Add(gAllMeanP,"[]");
                                            gStyle->SetEndErrorSize(4);
                                        }
                                        else{
                                            P3MultiAllMean->Add(gAllMeanP,"5");
                                        }
                                    }
                                    //SYS
                                    
                                    if (only1flag == 1) {
                                        if (n == 0) {
                                            legMean= new TLegend(0.7,0.7,0.9,0.9);
                                        }
                                        
                                        if (particleloop == 1) {
                                            if ((n == 0) && (cbmax[0] == 100)) {
                                                sprintf(AEword,"%s (Min.Bias)",partnameComplex);//mixed
                                            }
                                            else{
                                                sprintf(AEword,"%s (%1.1f-%1.1f %%)",partnameComplex,cbmin[n],cbmax[n]);//mixed
                                            }
                                            legMean->AddEntry(gAllMeanP,AEword,"l");
                                        }
                                    }
                                    else{
                                        if (n == 0) {
                                            legMean= new TLegend(0.7,0.7,0.9,0.9);
                                        }
                                        
                                        if ((n == 0) && (cbmax[0] == 100)) {
                                            sprintf(AEword,"%s (Min.Bias)",partnameComplex);//mixed
                                        }
                                        else{
                                            sprintf(AEword,"%s (%1.1f-%1.1f %%)",partnameComplex,cbmin[n],cbmax[n]);//mixed
                                        }
                                        legMean->AddEntry(gAllMeanP,AEword,"l");
                                    }
                                }
                                //new to remove
                                
                                for (int y=0; y < 4; y++) {
                                    Yforplot3[y]=WforplotP[1][y][n];//width Vs. pt
                                    YforplotError3[y]=WforplotPError[1][y][n];//width Vs. pt
                                    YforplotError3SYS[y]=WforplotPErrorSYS[1][y][n];
                                }
                                gAllWidthP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                
                                if (BellSwitch == 1) {
                                    for (int y=0; y < 3; y++) {
                                        Yforplot3[y]=WforplotP[1][y][n];//width Vs. pt
                                        YforplotError3[y]=WforplotPError[1][y][n];//width Vs. pt
                                        YforplotError3SYS[y]=WforplotPErrorSYS[1][y][n];
                                    }
                                    gAllWidthP = new TGraphErrors(3,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                }
                                
                                if (particleloop == 0) {
                                    gAllWidthP->SetLineStyle(1);//type
                                }
                                else if (particleloop == 1){
                                    gAllWidthP->SetLineStyle(3);//type
                                }
                                gAllWidthP->SetLineColor(1+n);
                                //new to remove
                                if (((n < 2) && (BellSwitch == 1)) || (BellSwitch == 0)) {//remove 10 and 30
                                    if (only1flag == 1) {
                                        if (particleloop == 1) {
                                            P3MultiAllWidth->Add(gAllWidthP);
                                        }
                                    }
                                    else{
                                        P3MultiAllWidth->Add(gAllWidthP);
                                    }
                                    
                                    if (only1flag == 1) {
                                        if (n == 0) {
                                            if (BellSwitch == 1) {
                                                legWidth= new TLegend(0.7,0.15,0.9,0.35);
                                            }
                                            else{
                                                legWidth= new TLegend(0.7,0.7,0.9,0.9);
                                            }
                                        }
                                        
                                        if (particleloop == 1) {
                                            if ((n == 0) && (cbmax[0] == 100)) {
                                                sprintf(AEword,"%s (Min.Bias)",partnameComplex);//mixed
                                            }
                                            else{
                                                sprintf(AEword,"%s (%1.1f-%1.1f %%)",partnameComplex,cbmin[n],cbmax[n]);//mixed
                                            }
                                            legWidth->AddEntry(gAllWidthP,AEword,"l");
                                        }
                                    }
                                    else{
                                        if (n == 0) {
                                            if (BellSwitch == 1) {
                                                legWidth= new TLegend(0.7,0.15,0.9,0.35);
                                            }
                                            else{
                                                legWidth= new TLegend(0.7,0.7,0.9,0.9);
                                            }
                                        }
                                        
                                        if ((n == 0) && (cbmax[0] == 100)) {
                                            sprintf(AEword,"%s (Min.Bias)",partnameComplex);//mixed
                                        }
                                        else{
                                            sprintf(AEword,"%s (%1.1f-%1.1f %%)",partnameComplex,cbmin[n],cbmax[n]);//mixed
                                        }
                                        legWidth->AddEntry(gAllWidthP,AEword,"l");
                                    }
                                    
                                    //SYS
                                    gAllWidthP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3SYS,YforplotError3SYS);
                                    if (particleloop == 0) {
                                        gAllWidthP->SetLineStyle(1);//type
                                    }
                                    else if (particleloop == 1){
                                        gAllWidthP->SetLineStyle(3);//type
                                    }
                                    gAllWidthP->SetLineColor(1+n);
                                    //gAllWidthP->SetDrawOption("2");
                                    //gAllWidthP->SetFillStyle(3001);
                                    ///gAllWidthP->SetFillColorAlpha(1+n,0.10);
                                    gAllWidthP->SetFillStyle(1);
                                    if (only1flag == 1) {
                                        if (particleloop == 1) {
                                            if (BellSwitch == 1) {
                                                P3MultiAllWidth->Add(gAllWidthP,"[]");
                                                gStyle->SetEndErrorSize(4);
                                            }
                                            else{
                                                P3MultiAllWidth->Add(gAllWidthP,"5");
                                            }
                                        }
                                    }
                                    else{
                                        if (BellSwitch == 1) {
                                            P3MultiAllWidth->Add(gAllWidthP,"[]");
                                            gStyle->SetEndErrorSize(4);
                                        }
                                        else{
                                            P3MultiAllWidth->Add(gAllWidthP,"5");
                                        }
                                    }
                                }
                                //new to remove
                                
                                /*if (only1flag == 1) {
                                 if (particleloop == 1) {
                                 sprintf(AEword,"%s systematic %1.1f-%1.1f %% Cen",partname,cbmin[n],cbmax[n]);
                                 legWidth->AddEntry(gAllWidthP,AEword,"f");
                                 }
                                 }
                                 else{
                                 sprintf(AEword,"%s systematic %1.1f-%1.1f %% Cen",partname,cbmin[n],cbmax[n]);
                                 legWidth->AddEntry(gAllWidthP,AEword,"f");
                                 }*/
                                //SYS
                                
                                if (particleloop == 1) {//ending
                                    if (n == (MaxCen-1)) {
                                        clog->cd();
                                        clog->SetLogy();
                                        clog->Update();
                                        if (BellSwitch == 1) {
                                            P3MultiAllYield->DrawClone("AP");
                                        }
                                        else{
                                            P3MultiAllYield->DrawClone("APL");
                                        }
                                        legYield->Draw("SAME");
                                        clog->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            clog->Print(Form("plotsCompleteV3/AllYield%s_Pt.pdf",partname));
                                            /*}//no cut
                                             else{//cut
                                             clog->Print(Form("plotsCompleteV3/cut%i_AllYield%s_Pt.pdf",cutloop,partname));
                                             }//cut*/
                                        }//Normal
                                        clog->Clear();
                                        clog->Update();
                                        
                                        cC->cd();
                                        cC->Update();
                                        P3MultiAllMean->DrawClone("AP");
                                        legMean->Draw("SAME");
                                        cC->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/AllMean%s_Pt.pdf",partname));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut%i_AllMean%s_Pt.pdf",cutloop,partname));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        cC->Update();
                                        
                                        cC->cd();
                                        cC->Update();
                                        P3MultiAllWidth->DrawClone("AP");
                                        legWidth->Draw("SAME");
                                        cC->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/AllWidth%s_Pt.pdf",partname));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut%i_AllWidth%s_Pt.pdf",cutloop,partname));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        cC->Update();
                                    }
                                    
                                }//ending
                            }//not 1 pt
                            //
                            /*/not used much anymore
                             if (Use1PT != 1) {//not 1 pt
                             cC->cd();
                             cC->Clear();
                             cC->Update();
                             SystematicAllYield->SetTitle(Form("Yields for %s %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                             SystematicAllYield->Draw();
                             cC->Update();
                             if (hmtloop == 0) {//Normal
                             //if (cutloop == 0) {//no cut
                             cC->Print(Form("plotsCompleteV3/%s/SYSYield_NR%i_c%i.pdf",partname,m,n));
                             //}//no cut
                             //else{//cut
                             //    cC->Print(Form("plotsCompleteV3/%s/cut%i_SYSYield_NR%i_c%i.pdf",partname,cutloop,m,n));
                             //}//cut
                             }//Normal
                             cC->Clear();
                             cC->Update();
                             //error starts here
                             SystematicAllMean->SetTitle(Form("Means for %s %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                             SystematicAllMean->DrawClone();
                             cC->Update();
                             if (hmtloop == 0) {//Normal
                             //if (cutloop == 0) {//no cut
                             cC->Print(Form("plotsCompleteV3/%s/SYSMean_NR%i_c%i.pdf",partname,m,n));
                             //}//no cut
                             //else{//cut
                             //    cC->Print(Form("plotsCompleteV3/%s/cut%i_SYSMean_NR%i_c%i.pdf",partname,cutloop,m,n));
                             //}//cut
                             }//Normal
                             cC->Clear();
                             cC->Update();
                             SystematicAllWidth->SetTitle(Form("Widths for %s %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                             SystematicAllWidth->Draw();
                             cC->Update();
                             if (hmtloop == 0) {//Normal
                             //if (cutloop == 0) {//no cut
                             cC->Print(Form("plotsCompleteV3/%s/SYSWidth_NR%i_c%i.pdf",partname,m,n));
                             //}//no cut
                             //else{//cut
                             //    cC->Print(Form("plotsCompleteV3/%s/cut%i_SYSWidth_NR%i_c%i.pdf",partname,cutloop,m,n));
                             //}//cut
                             }//Normal
                             cC->Clear();
                             cC->Update();
                             StatisticalAllMean->SetTitle(Form("Statistical Means for %s %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                             StatisticalAllMean->Draw();
                             cC->Update();
                             if (hmtloop == 0) {//Normal
                             //if (cutloop == 0) {//no cut
                             cC->Print(Form("plotsCompleteV3/%s/STAMean_NR%i_c%i.pdf",partname,m,n));
                             //}//no cut
                             //else{//cut
                             //    cC->Print(Form("plotsCompleteV3/%s/cut%i_STAMean_NR%i_c%i.pdf",partname,cutloop,m,n));
                             //}//cut
                             }//Normal
                             cC->Clear();
                             cC->Update();
                             StatisticalAllWidth->SetTitle(Form("Statistical Widths for %s %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                             StatisticalAllWidth->Draw();
                             cC->Update();
                             if (hmtloop == 0) {//Normal
                             //if (cutloop == 0) {//no cut
                             cC->Print(Form("plotsCompleteV3/%s/STAWidth_NR%i_c%i.pdf",partname,m,n));
                             //}//no cut
                             //else{//cut
                             //    cC->Print(Form("plotsCompleteV3/%s/cut%i_STAWidth_NR%i_c%i.pdf",partname,cutloop,m,n));
                             //}//cut
                             }//Normal
                             cC->Clear();
                             cC->Update();
                             if (n == 0) {
                             cSYSY->Clear();
                             cSYSY->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                             cSYSY->Update();
                             cSYSM->Clear();
                             cSYSM->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                             cSYSM->Update();
                             cSYSW->Clear();
                             cSYSW->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                             cSYSW->Update();
                             cSTAM->Clear();
                             cSTAM->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                             cSTAM->Update();
                             cSTAW->Clear();
                             cSTAW->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                             cSTAW->Update();
                             }
                             cSYSY->cd(n+1);
                             cSYSY->Update();
                             SystematicAllYield->Draw();
                             cSYSY->Update();
                             cSYSM->cd(n+1);
                             cSYSM->Update();
                             SystematicAllMean->Draw();
                             cSYSM->Update();
                             cSYSW->cd(n+1);
                             cSYSW->Update();
                             SystematicAllWidth->Draw();
                             cSYSW->Update();
                             cSTAM->cd(n+1);
                             cSTAM->Update();
                             StatisticalAllMean->Draw();
                             cSTAM->Update();
                             cSTAW->cd(n+1);
                             cSTAW->Update();
                             StatisticalAllWidth->Draw();
                             cSTAW->Update();
                             if (n == (MaxCen-1)) {
                             if (hmtloop == 0) {//Normal
                             //if (cutloop == 0) {//no cut
                             cSYSY->Print(Form("plotsCompleteV3/%s/TogetherSYSYield_NR%i.pdf",partname,m));
                             cSYSM->Print(Form("plotsCompleteV3/%s/TogetherSYSMean_NR%i.pdf",partname,m));
                             cSYSW->Print(Form("plotsCompleteV3/%s/TogetherSYSWidth_NR%i.pdf",partname,m));
                             cSTAM->Print(Form("plotsCompleteV3/%s/TogetherSTAMean_NR%i.pdf",partname,m));
                             cSTAW->Print(Form("plotsCompleteV3/%s/TogetherSTAWidth_NR%i.pdf",partname,m));
                             //}//no cut
                             //else{//cut
                             //    cSYSY->Print(Form("plotsCompleteV3/%s/cut1_TogetherSYSYield_NR%i.pdf",partname,m));
                             //    cSYSM->Print(Form("plotsCompleteV3/%s/cut1_TogetherSYSMean_NR%i.pdf",partname,m));
                             //    cSYSW->Print(Form("plotsCompleteV3/%s/cut1_TogetherSYSWidth_NR%i.pdf",partname,m));
                             //    cSTAM->Print(Form("plotsCompleteV3/%s/cut1_TogetherSTAMean_NR%i.pdf",partname,m));
                             //    cSTAW->Print(Form("plotsCompleteV3/%s/cut1_TogetherSTAWidth_NR%i.pdf",partname,m));
                             //}//cut
                             }//Normal
                             }
                             }//not 1 pt
                             */ //not used much anymore
                            
                            //GM section
                            //return here
                            //not used much anymore
                            
                            if (particleloop == 0) {// use only LambdaKX for GM section
                                
                                /*
                                 MultiAllYield = new TMultiGraph();
                                 for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type
                                 for (int k=0; k < 1; k++) {//version
                                 for (int j=0; j < 1; j++) {//range
                                 //for (int m=0; m < (NREnd-NRBegin); m++) {
                                 for (int i=0; i < (MaxPT-MinPT); i++) {//pt
                                 Xforplot[i]=i;
                                 Yforplot[i]=0;
                                 Yforplot[i]=GMALLYield[n][i][m][j][k][l];
                                 }
                                 //cC->cd();
                                 TGraph* gAllYield = new TGraph((MaxPT-MinPT),Xforplot,Yforplot);
                                 gAllYield->SetLineColor(l+1);//type
                                 gAllYield->SetLineStyle(k);//version
                                 gAllYield->SetMarkerStyle(20+j);//range
                                 MultiAllYield->Add(gAllYield);
                                 }
                                 //}
                                 }
                                 }
                                 cC->cd();
                                 cC->Clear();
                                 cC->Update();
                                 MultiAllYield->SetTitle(Form("All-Yeild All Type GM vs. pt for Cen %i %s %s",n,partname,normal));
                                 MultiAllYield->DrawClone("AP");
                                 cC->Update();
                                 if (hmtloop == 0) {//Normal
                                 //if (cutloop == 0) {//no cut
                                 cC->Print(Form("plotsCompleteV3/%s/AllYield_GM_NR%i_c%i.pdf",partname,m,n));//error here
                                 //}//no cut
                                 //else{//cut
                                 //    cC->Print(Form("plotsCompleteV3/%s/cut1_AllYield_GM_NR%i_c%i.pdf",partname,m,n));//error here
                                 //}//cut
                                 }//Normal
                                 else{//hmt
                                 //if (cutloop == 0) {//no cut
                                 cC->Print(Form("plotsCompleteV3/%s/hmt_AllYield_GM_NR%i_c%i.pdf",partname,m,n));//error here
                                 //}//no cut
                                 //else{//cut
                                 //    cC->Print(Form("plotsCompleteV3/%s/hmtcut1_AllYield_GM_NR%i_c%i.pdf",partname,m,n));//error here
                                 //}//cut
                                 }//hmt
                                 cC->Clear();
                                 */ //not used much anymore
                                
                                //TH1D* SystematicAllYield = new TH1D(Form("SystematicYield_%s_%s",partname,othernormal),"",4,0,4);
                                //TH1D* SystematicAllMean = new TH1D(Form("SystematicMean_%s_%s",partname,othernormal),"",4,0,4);
                                //TH1D* SystematicAllWidth = new TH1D(Form("SystematicWidth_%s_%s",partname,othernormal),"",4,0,4);
                                //TH1D* StatisticalAllMean = new TH1D(Form("StatisticalMean_%s_%s",partname,othernormal),"",4,0,4);
                                //TH1D* StatisticalAllWidth = new TH1D(Form("StatisticalWidth_%s_%s",partname,othernormal),"",4,0,4);
                                calculateSYSmean=0;
                                calculateSYSmeanerror=0;
                                calculateSTAmeanerror=0;
                                calculateSYSyield=0;
                                calculateSYSyielderror=0;
                                calculateSYSwidth=0;
                                calculateSYSwidtherror=0;
                                calculateSTAwidtherror=0;
                                numberofmean=0;
                                numberofwidth=0;
                                numberofyield=0;
                                for (int i=0; i < (MaxPT-MinPT); i++) {//pt
                                    numberofyield=0;
                                    calculateSYSyield=0;
                                    calculateSYSyielderror=0;
                                    calculateSYSyieldBin=0;
                                    calculateSYSyielderrorBin=0;
                                    numberofmean=0;
                                    calculateSYSmean=0;
                                    calculateSYSmeanerror=0;
                                    calculateSTAmeanerror=0;
                                    numberofwidth=0;
                                    calculateSYSwidth=0;
                                    calculateSYSwidtherror=0;
                                    calculateSTAwidtherror=0;
                                    for (int j=0; j < 3; j++) {//range //rangestart-1//3
                                        if ((BellSwitch == 1) && (j != 0)) {
                                            continue;
                                        }
                                        for (int k=1; k < 4; k++) {//version k=0; k < 3
                                            for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type l=0; l < 8
                                                for (int m=0; m < (NREnd-NRBegin); m++) {
                                                    if (BellSwitch == 1) {
                                                        m=2;
                                                    }
                                                    
                                                    if (onechange == 1) {
                                                        //default for particle 0 (LAmbda+KX) l=0, k=1, j=0, m=2
                                                        if ((l%2 == 0) && (l == 0) && (k == 1) && (j == 0) /*&& (m == 2)*/) {
                                                            //this is changing NR
                                                        }
                                                        else if((l%2 == 0) && (l == 0) && (k == 1) && /*(j == 0) &&*/ (m == 2)){
                                                            //this is changing range
                                                        }
                                                        else if((l%2 == 0) && (l == 0) && /*(k == 1) &&*/ (j == 0) && (m == 2)){
                                                            //this is changing fit range
                                                        }
                                                        else if((l == 1) && (k == 1) && (j == 0) && (m == 2)){
                                                            //this is changing sigma
                                                        }
                                                        else if(/*(l%2 == 0) &&*/ (l == 0) && (k == 1) && (j == 0) && (m == 2)){
                                                            //this is changing NSigma
                                                        }
                                                        else{
                                                            continue;//this is not needed
                                                        }
                                                    }
                                                    
                                                    if (BellSwitch == 1) {
                                                        m=0;
                                                    }
                                                    
                                                    calculateSYSyield += GMALLYield[n][i][m][j][k][l];
                                                    calculateSYSyieldBin += BinGMALLYield[n][i][m][j][k][l];
                                                    numberofyield++;
                                                    numberofmean++;
                                                    calculateSYSmean += GMALLMean[n][i][m][j][k][l];
                                                    calculateSTAmeanerror += pow(GMALLMeanError[n][i][m][j][k][l],2);
                                                    numberofwidth++;
                                                    calculateSYSwidth += GMALLWidth[n][i][m][j][k][l];
                                                    calculateSTAwidtherror += pow(GMALLWidthError[n][i][m][j][k][l],2);
                                                }
                                            }
                                        }
                                    }
                                    if (numberofyield == 0) {
                                        numberofyield=1;
                                        numberofmean=1;
                                        numberofwidth=1;
                                    }
                                    calculateSYSyieldBin=calculateSYSyieldBin/numberofyield;
                                    
                                    calculateSYSyield=calculateSYSyield/numberofyield;
                                    SystematicAllYield->SetBinContent(i+1,calculateSYSyield);
                                    calculateSYSmean=calculateSYSmean/numberofmean;
                                    SystematicAllMean->SetBinContent(i+1,calculateSYSmean);
                                    calculateSYSwidth=calculateSYSwidth/numberofwidth;
                                    SystematicAllWidth->SetBinContent(i+1,calculateSYSwidth);
                                    calculateSTAmeanerror=sqrt(calculateSTAmeanerror)/numberofmean;
                                    StatisticalAllMean->SetBinContent(i+1,calculateSYSmean);
                                    StatisticalAllMean->SetBinError(i+1,calculateSTAmeanerror);
                                    calculateSTAwidtherror=sqrt(calculateSTAwidtherror)/numberofwidth;
                                    StatisticalAllWidth->SetBinContent(i+1,calculateSYSwidth);
                                    StatisticalAllWidth->SetBinError(i+1,calculateSTAwidtherror);
                                    //P#Multi
                                    if (Use1PT == 1) {//use 1 pt
                                        //YforplotP[0][3][n]=calculateSYSyield;
                                        MforplotP[0][3][n]=calculateSYSmean;
                                        WforplotP[0][3][n]=calculateSYSwidth;
                                        if (Use1Type == 1) {
                                            //YforplotP[0][3][n]=((calculateSYSyield+calculateSYSyieldBin)/2.0)/((cbmax[n]-cbmin[n])/10.0);
                                            //YforplotP[0][3][n]=((calculateSYSyield+calculateSYSyieldBin)/2.0)/((cbmax[n]-cbmin[n]));
                                            YforplotP[0][3][n]=((calculateSYSyield+calculateSYSyieldBin)/2.0);
                                        }
                                        else{
                                            //YforplotP[0][3][n]=calculateSYSyield/(cbmax[n]-cbmin[n]);
                                            YforplotP[0][3][n]=calculateSYSyield;
                                        }
                                    }//use 1 pt
                                    else{
                                        MforplotP[0][i][n]=calculateSYSmean;
                                        WforplotP[0][i][n]=calculateSYSwidth;
                                        if (Use1Type == 1) {
                                            YforplotP[0][i][n]=((calculateSYSyield+calculateSYSyieldBin)/2.0);
                                        }
                                        else{
                                            YforplotP[0][i][n]=calculateSYSyield;
                                        }
                                        /*
                                         YforplotP[0][i][n]=calculateSYSyield;
                                         MforplotP[0][i][n]=calculateSYSmean;
                                         WforplotP[0][i][n]=calculateSYSwidth;
                                         if(i == 3){
                                         //YforplotP[0][i][n]=calculateSYSyield/((cbmax[n]-cbmin[n]));
                                         YforplotP[0][i][n]=calculateSYSyield;
                                         //YforplotP3GM[n]=calculateSYSyield;
                                         }*/
                                    }
                                    //P#Multi
                                    for (int j=0; j < 3; j++) {//range //rangestart-1 //3
                                        if ((BellSwitch == 1) && (j != 0)) {
                                            continue;
                                        }
                                        for (int k=1; k < 4; k++) {//version k=0; k < 3
                                            for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type l=0; l < 8
                                                for (int m=0; m < (NREnd-NRBegin); m++) {
                                                    if (BellSwitch == 1) {
                                                        m=2;
                                                    }
                                                    
                                                    if (onechange == 1) {
                                                        //default for particle 0 (LAmbda+KX) l=0, k=1, j=0, m=2
                                                        if ((l%2 == 0) && (l == 0) && (k == 1) && (j == 0) /*&& (m == 2)*/){
                                                            //this is changing NR
                                                        }
                                                        else if((l%2 == 0) && (l == 0) && (k == 1) && /*(j == 0) &&*/ (m == 2)){
                                                            //this is changing range
                                                        }
                                                        else if((l%2 == 0) && (l == 0) && /*(k == 1) &&*/ (j == 0) && (m == 2)){
                                                            //this is changing fit range
                                                        }
                                                        else if((l == 1) && (k == 1) && (j == 0) && (m == 2)){
                                                            //this is changing sigma
                                                        }
                                                        else if(/*(l%2 == 0) &&*/ (l == 0) && (k == 1) && (j == 0) && (m == 2)){
                                                            //this is changing NSigma
                                                        }
                                                        else{
                                                            continue;//this is not needed
                                                        }
                                                    }
                                                    
                                                    if (BellSwitch == 1) {
                                                        m=0;
                                                    }
                                                    
                                                    calculateSYSyielderror += pow(GMALLYield[n][i][m][j][k][l] - calculateSYSyield,2);
                                                    calculateSYSmeanerror += pow(GMALLMean[n][i][m][j][k][l] - calculateSYSmean,2);
                                                    calculateSYSwidtherror += pow(GMALLWidth[n][i][m][j][k][l] - calculateSYSwidth,2);
                                                }
                                            }
                                        }
                                    }
                                    calculateSYSyielderror=sqrt(pow(calculateSYSyield-calculateSYSyieldBin,2)+pow(sqrt(calculateSYSyielderror/numberofyield),2));//sqrt(calculateSYSyielderror/numberofyield)
                                    SystematicAllYield->SetBinError(i+1,calculateSYSyielderror);
                                    calculateSYSmeanerror=sqrt(calculateSYSmeanerror/numberofmean);
                                    SystematicAllMean->SetBinError(i+1,calculateSYSmeanerror);
                                    calculateSYSwidtherror=sqrt(calculateSYSwidtherror/numberofwidth);
                                    SystematicAllWidth->SetBinError(i+1,calculateSYSwidtherror);
                                    if (Use1PT == 1) {//use 1 pt
                                        //YforplotPError[0][i][n]=calculateSYSyielderror;
                                        MforplotPError[0][3][n]=sqrt(pow(calculateSYSmeanerror,2) + pow(calculateSTAmeanerror,2));
                                        WforplotPError[0][3][n]=sqrt(pow(calculateSYSwidtherror,2) + pow(calculateSTAwidtherror,2));
                                        if (Use1Type == 1) {
                                            //YforplotPError[0][3][n]=((calculateSYSyield-calculateSYSyieldBin)/2.0)/((cbmax[n]-cbmin[n])/10.0);
                                            //YforplotPError[0][3][n]=((calculateSYSyield-calculateSYSyieldBin)/2.0)/((cbmax[n]-cbmin[n]));
                                            YforplotPError[0][3][n]=TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0));
                                        }
                                        else{
                                            //YforplotPError[0][3][n]=((calculateSYSyield-calculateSYSyieldBin)/2.0)/((cbmax[n]-cbmin[n]));
                                            YforplotPError[0][3][n]=TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0));
                                            //YforplotPError[0][3][n]=calculateSYSyielderror/(cbmax[n]-cbmin[n]);
                                        }
                                        //slight modification
                                        MforplotPError[0][3][n]=calculateSTAmeanerror;
                                        MforplotPErrorSYS[0][3][n]=calculateSYSmeanerror;
                                        
                                        WforplotPError[0][3][n]=calculateSTAwidtherror;
                                        WforplotPErrorSYS[0][3][n]=calculateSYSwidtherror;
                                        //slight modification to be removed
                                    }
                                    else{
                                        //MforplotPError[0][i][n]=sqrt(pow(calculateSYSmeanerror,2) + pow(calculateSTAmeanerror,2));
                                        MforplotPError[0][i][n]=calculateSTAmeanerror;
                                        MforplotPErrorSYS[0][i][n]=calculateSYSmeanerror;
                                        //WforplotPError[0][i][n]=sqrt(pow(calculateSYSwidtherror,2) + pow(calculateSTAwidtherror,2));
                                        WforplotPError[0][i][n]=calculateSTAwidtherror;
                                        WforplotPErrorSYS[0][i][n]=calculateSYSwidtherror;
                                        if (Use1Type == 1) {
                                            YforplotPError[0][i][n]=TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0));
                                        }
                                        else{
                                            YforplotPError[0][i][n]=TMath::Abs(((calculateSYSyield-calculateSYSyieldBin)/2.0));
                                        }
                                        /*
                                         YforplotPError[0][i][n]=calculateSYSyielderror;
                                         MforplotPError[0][i][n]=sqrt(pow(calculateSYSmeanerror,2) + pow(calculateSTAmeanerror,2));
                                         WforplotPError[0][i][n]=sqrt(pow(calculateSYSwidtherror,2) + pow(calculateSTAwidtherror,2));
                                         if(i == 3){
                                         //YforplotPError[0][i][n]=calculateSYSyielderror/(cbmax[n]-cbmin[n]);
                                         YforplotPError[0][i][n]=calculateSYSyielderror;
                                         }*/
                                    }
                                    
                                    SystematicAllYield->GetXaxis()->SetBinLabel(i+1,Form("%1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[i+MinPT],ptmax[i+MinPT]));
                                    SystematicAllYield->GetYaxis()->SetTitle("Yield");
                                    SystematicAllMean->GetXaxis()->SetBinLabel(i+1,Form("%1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[i+MinPT],ptmax[i+MinPT]));
                                    SystematicAllMean->GetYaxis()->SetTitle("Mean (GeV/c^{2})");
                                    SystematicAllWidth->GetXaxis()->SetBinLabel(i+1,Form("%1.1f<#it{p}_{T}<%1.1f GeV/c",ptmin[i+MinPT],ptmax[i+MinPT]));
                                    SystematicAllWidth->GetYaxis()->SetTitle("Width (GeV/c^{2})");
                                    if (writetxt == 1) {//begin systematic txt file
                                        if (MCflag == 1) {
                                            if (hmtloop == 0) {//Normal
                                                //if (cutloop == 0) {//no cut
                                                fprintf(Cutoftxt,"%s GM | Normal   |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                /*}//no cut
                                                 else{//cut
                                                 fprintf(Cutoftxt,"%s GM | cut1     |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                 }//cut*/
                                            }//Normal
                                            else{//hmt
                                                //if (cutloop == 0) {//no cut
                                                fprintf(Cutoftxt,"%s GM | hmt      |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                /*}//no cut
                                                 else{//cut
                                                 fprintf(Cutoftxt,"%s GM | hmtcut1  |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                 }//cut*/
                                            }//hmt
                                        }
                                        else{
                                            if (hmtloop == 0) {//Normal
                                                //if (cutloop == 0) {//no cut
                                                fprintf(Cutoftxt,"%s GM | Normal   |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                /*}//no cut
                                                 else{//cut
                                                 fprintf(Cutoftxt,"%s GM | cut1     |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                 }//cut*/
                                            }//Normal
                                            else{//hmt
                                                //if (cutloop == 0) {//no cut
                                                fprintf(Cutoftxt,"%s GM | hmt      |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                /*}//no cut
                                                 else{//cut
                                                 fprintf(Cutoftxt,"%s GM | hmtcut1  |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                 }//cut*/
                                            }//hmt
                                        }
                                        
                                        //percentage
                                        if (MCflag == 1) {
                                            if (hmtloop == 0) {//Normal
                                                //if (cutloop == 0) {//no cut
                                                fprintf(Percentage,"%s GM | Normal   |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                /*}//no cut
                                                 else{//cut
                                                 fprintf(Percentage,"%s GM | cut1     |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                 }//cut*/
                                            }//Normal
                                            else{//hmt
                                                //if (cutloop == 0) {//no cut
                                                fprintf(Percentage,"%s GM | hmt      |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                /*}//no cut
                                                 else{//cut
                                                 fprintf(Percentage,"%s GM | hmtcut1  |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f |%1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield/AEarray[n][i],calculateSYSyielderror/AEarray[n][i],AEarray[n][i],AEarrayerror[n][i],AEslope[n][i],AEslopeerror[n][i],calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                 }//cut*/
                                            }//hmt
                                        }
                                        else{
                                            if (hmtloop == 0) {//Normal
                                                //if (cutloop == 0) {//no cut
                                                fprintf(Percentage,"%s GM | Normal   |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                /*}//no cut
                                                 else{//cut
                                                 fprintf(Percentage,"%s GM | cut1     |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                 }//cut*/
                                            }//Normal
                                            else{//hmt
                                                //if (cutloop == 0) {//no cut
                                                fprintf(Percentage,"%s GM | hmt      |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                /*}//no cut
                                                 else{//cut
                                                 fprintf(Percentage,"%s GM | hmtcut1  |%1.1f-%1.1f |%1.1f<pt<%1.1f |%5.1f +- %5.1f |%1.4f +- %1.4f +- %1.4f |%1.4f +- %1.4f +- %1.4f| \n",partname,cbmin[n],cbmax[n],ptmin[i+MinPT],ptmax[i+MinPT],calculateSYSyield,calculateSYSyielderror,calculateSYSmean,calculateSYSmeanerror,calculateSTAmeanerror,calculateSYSwidth,calculateSYSwidtherror,calculateSTAwidtherror);
                                                 }//cut*/
                                            }//hmt
                                        }
                                        //percentage calculation section
                                        tempYield=0;
                                        tempMean=0;
                                        tempWidth=0;
                                        for (int j=0; j < 3; j++) {//range j=0; j < 3
                                            if ((BellSwitch == 1) && (j != 0)) {
                                                continue;
                                            }
                                            for (int k=1; k < 4; k++) {//version k=0; k < 3
                                                for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type l=0; l < 8
                                                    for (int m=0; m < (NREnd-NRBegin); m++) {
                                                        if (BellSwitch == 1) {
                                                            m=2;
                                                        }
                                                        
                                                        if (onechange == 1) {
                                                            //default for particle 0 (LAmbda+KX) l=0, k=1, j=0, m=2
                                                            if ((l%2 == 0) && (l == 0) && (k == 1) && (j == 0) /*&& (m == 2)*/){
                                                                //this is changing NR
                                                            }
                                                            else if((l%2 == 0) && (l == 0) && (k == 1) && /*(j == 0) &&*/ (m == 2)){
                                                                //this is changing range
                                                            }
                                                            else if((l%2 == 0) && (l == 0) && /*(k == 1) &&*/ (j == 0) && (m == 2)){
                                                                //this is changing fit range
                                                            }
                                                            else if((l == 1) && (k == 1) && (j == 0) && (m == 2)){
                                                                //this is changing sigma
                                                            }
                                                            else if(/*(l%2 == 0) &&*/ (l == 0) && (k == 1) && (j == 0) && (m == 2)){
                                                                //this is changing NSigma
                                                            }
                                                            else{
                                                                continue;//this is not needed
                                                            }
                                                        }
                                                        
                                                        if (BellSwitch == 1) {
                                                            m=0;
                                                        }
                                                        
                                                        tempYield+=TMath::Abs(GMALLYield[n][i][m][j][k][l]-calculateSYSyield);
                                                        tempMean+=TMath::Abs(GMALLMean[n][i][m][j][k][l]-calculateSYSmean);
                                                        tempWidth+=TMath::Abs(GMALLWidth[n][i][m][j][k][l]-calculateSYSwidth);
                                                    }
                                                }
                                            }
                                        }
                                        //newer
                                        tempYield+=TMath::Abs(calculateSYSyield-calculateSYSyieldBin);//newer
                                        //newer
                                        for (int j=0; j < 3; j++) {//range j=0; j < 3
                                            if ((BellSwitch == 1) && (j != 0)) {
                                                continue;
                                            }
                                            for (int k=1; k < 4; k++) {//version k=0; k < 3
                                                for (int l=0; l < (typeend-typestart+1+(4*Nsigmaflag)); l++) {//type l=0; l < 8
                                                    for (int m=0; m < (NREnd-NRBegin); m++) {
                                                        if (BellSwitch == 1) {
                                                            m=2;
                                                        }
                                                        
                                                        if (onechange == 1) {
                                                            //default for particle 0 (LAmbda+KX) l=0, k=1, j=0, m=2
                                                            if ((l%2 == 0) && (l == 0) && (k == 1) && (j == 0) /*&& (m == 2)*/){
                                                                //this is changing NR
                                                            }
                                                            else if((l%2 == 0) && (l == 0) && (k == 1) && /*(j == 0) &&*/ (m == 2)){
                                                                //this is changing range
                                                            }
                                                            else if((l%2 == 0) && (l == 0) && /*(k == 1) &&*/ (j == 0) && (m == 2)){
                                                                //this is changing fit range
                                                            }
                                                            else if((l == 1) && (k == 1) && (j == 0) && (m == 2)){
                                                                //this is changing sigma
                                                            }
                                                            else if(/*(l%2 == 0) &&*/ (l == 0) && (k == 1) && (j == 0) && (m == 2)){
                                                                //this is changing NSigma
                                                            }
                                                            else{
                                                                continue;//this is not needed
                                                            }
                                                        }
                                                        
                                                        if (BellSwitch == 1) {
                                                            m=0;
                                                        }
                                                        
                                                        //fprintf(Percentage,"range %i |version %i|type %i   | %2.4f | %2.4f | %2.4f | %1.4f | %1.4f |\n",j,k,l,(TMath::Abs((GMALLYield[n][i][m][j][k][l]-calculateSYSyield)/sqrt(numberofyield)/calculateSYSyielderror)*100),(TMath::Abs((GMALLMean[n][i][m][j][k][l]-calculateSYSmean)/sqrt(numberofmean)/calculateSYSmeanerror)*100),(TMath::Abs((GMALLWidth[n][i][m][j][k][l]-calculateSYSwidth)/sqrt(numberofwidth)/calculateSYSwidtherror)*100),GMALLMean[n][i][m][j][k][l],GMALLWidth[n][i][m][j][k][l]);
                                                        //fprintf(Percentage,"range %i |version %i|type %i   | %2.4f | %2.4f | %2.4f | %1.4f | %1.4f |\n",j,k,l,(TMath::Abs((GMALLYield[n][i][m][j][k][l]-calculateSYSyield)/(calculateSYSyielderror*sqrt(numberofyield)))*100),(TMath::Abs((GMALLMean[n][i][m][j][k][l]-calculateSYSmean)/(calculateSYSmeanerror*sqrt(numberofmean)))*100),(TMath::Abs((GMALLWidth[n][i][m][j][k][l]-calculateSYSwidth)/(calculateSYSwidtherror*sqrt(numberofwidth)))*100),GMALLMean[n][i][m][j][k][l],GMALLWidth[n][i][m][j][k][l]);
                                                        fprintf(Percentage,"Nor.Ran. %i |range %i |version %i|type %i   | %2.4f | %2.4f | %2.4f | %5.1f | %1.4f | %1.4f | ",m,j,k,l,(TMath::Abs((GMALLYield[n][i][m][j][k][l]-calculateSYSyield)/(tempYield))*100),(TMath::Abs((GMALLMean[n][i][m][j][k][l]-calculateSYSmean)/(tempMean))*100),(TMath::Abs((GMALLWidth[n][i][m][j][k][l]-calculateSYSwidth)/(tempWidth))*100),GMALLYield[n][i][m][j][k][l]/AEarray[n][i],GMALLMean[n][i][m][j][k][l],GMALLWidth[n][i][m][j][k][l]);
                                                        if ((TMath::Abs((GMALLYield[n][i][m][j][k][l]-calculateSYSyield)/(tempYield))*100) >= 10.0) {//high yield error
                                                            fprintf(Percentage," *Y ");
                                                        }
                                                        if ((TMath::Abs((GMALLMean[n][i][m][j][k][l]-calculateSYSmean)/(tempMean))*100) >= 10.0) {//high Mean error
                                                            fprintf(Percentage," *M ");
                                                        }
                                                        if ((TMath::Abs((GMALLWidth[n][i][m][j][k][l]-calculateSYSwidth)/(tempWidth))*100) >= 10.0) {//high width error
                                                            fprintf(Percentage," *W ");
                                                        }
                                                        
                                                        if (l == 1) {//only needed once
                                                            fprintf(Percentage," Int-Bin error | %1.4f |", (TMath::Abs(calculateSYSyield-calculateSYSyieldBin)/(tempYield))*100);
                                                        }//only needed once
                                                        fprintf(Percentage,"\n");
                                                    }
                                                }
                                            }
                                        }
                                        //percentage calculation section
                                        //percentage
                                    }//end systematic text file
                                }//end of mini pt loop
                                
                                if (Use1PT == 1) {//use 1 pt
                                    if (n == (MaxCen-1)) {
                                        if ((cbmin[0] == 0) && (cbmax[0] == 100)) {
                                            for (int loop3=0; loop3 < (MaxCen-1); loop3++) {
                                                Xforplot3[loop3]=((cbmin[loop3+1]+cbmax[loop3+1])/2.0);
                                                if (BellSwitch == 1) {
                                                    XforplotError3[loop3]=0;
                                                }
                                                else{
                                                    XforplotError3[loop3]=((cbmax[loop3+1]-cbmin[loop3+1])/2.0);
                                                }
                                            }
                                            
                                            for (int y=0; y < (MaxCen-1); y++) {
                                                Yforplot3[y]=YforplotP[0][3][y+1];
                                                YforplotError3[y]=YforplotPError[0][3][y+1];
                                            }
                                            gAllYieldP = new TGraphErrors((MaxCen-1),Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        else{
                                            for (int loop3=0; loop3 < MaxCen; loop3++) {
                                                Xforplot3[loop3]=((cbmin[loop3]+cbmax[loop3])/2.0);
                                                if (BellSwitch == 1) {
                                                    XforplotError3[loop3]=0;
                                                }
                                                else{
                                                    XforplotError3[loop3]=((cbmax[loop3]-cbmin[loop3])/2.0);
                                                }
                                            }
                                            
                                            for (int y=0; y < MaxCen; y++) {
                                                Yforplot3[y]=YforplotP[0][3][y];
                                                YforplotError3[y]=YforplotPError[0][3][y];
                                            }
                                            gAllYieldP = new TGraphErrors(MaxCen,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        
                                        cC->cd();
                                        cC->Update();
                                        gAllYieldP->DrawClone("AP");
                                        cC->Modified();
                                        cC->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/JustYield_Pt.pdf"));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_JustYield_Pt.pdf"));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        
                                        gAllYieldP->SetLineColor(2);//type
                                        P3MultiAllYield->Add(gAllYieldP);
                                        
                                        if ((cbmin[0] == 0) && (cbmax[0] == 100)) {
                                            for (int y=0; y < (MaxCen-1); y++) {
                                                Yforplot3[y]=MforplotP[0][3][y+1];
                                                YforplotError3[y]=MforplotPError[0][3][y+1];
                                            }
                                            gAllMeanP = new TGraphErrors((MaxCen-1),Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        else{
                                            for (int y=0; y < MaxCen; y++) {
                                                Yforplot3[y]=MforplotP[0][3][y];
                                                YforplotError3[y]=MforplotPError[0][3][y];
                                            }
                                            gAllMeanP = new TGraphErrors(MaxCen,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        
                                        gAllMeanP->SetLineColor(2);//type
                                        fit1 = new TF1("fit1","pol1",0,100);
                                        //gAllWidthP->Fit("pol1");
                                        fit1->SetParLimits(0,1.81,1.85);
                                        //gAllWidthP->Fit("fit1","R");
                                        cC->cd();
                                        if (BellSwitch == 0) {
                                            gStyle->SetOptFit(111);
                                        }
                                        cC->Update();
                                        gAllMeanP->Draw("AP");
                                        gAllMeanP->Fit("fit1","B");
                                        //new section, moving stat box
                                        
                                        cC->Update();//update
                                        //fit1->DrawCopy("SAME");
                                        //cC->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/JustMean_Pt.pdf"));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_JustMean_Pt.pdf"));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        P3MultiAllMean->Add(gAllMeanP);
                                        
                                        if ((cbmin[0] == 0) && (cbmax[0] == 100)) {
                                            for (y=0; y < (MaxCen-1); y++) {
                                                Yforplot3[y]=WforplotP[0][3][y+1];
                                                YforplotError3[y]=WforplotPError[0][3][y+1];
                                                YforplotError3SYS[y]=WforplotPErrorSYS[0][3][y+1];
                                            }
                                            gAllWidthP = new TGraphErrors((MaxCen-1),Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        else{
                                            for (y=0; y < MaxCen; y++) {
                                                Yforplot3[y]=WforplotP[0][3][y];
                                                YforplotError3[y]=WforplotPError[0][3][y];
                                                YforplotError3SYS[y]=WforplotPErrorSYS[0][3][y];
                                            }
                                            gAllWidthP = new TGraphErrors(MaxCen,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        
                                        gAllWidthP->SetLineColor(2);//type
                                        gAllWidthP->SetTitle(Form("Width (FWHM) for %s 1.5<p_{t}<20.0 GeV/c",partnameComplex));
                                        gAllWidthP->GetXaxis()->SetTitle("% ");
                                        gAllWidthP->GetYaxis()->SetTitle("Width (GeV/c^{2})");
                                        fit1 = new TF1("fit1","pol1",0,100);
                                        //gAllWidthP->Fit("pol1");
                                        fit1->SetParLimits(0,0.020,0.045);
                                        //gAllWidthP->Fit("fit1","R");
                                        cC->cd();
                                        if (BellSwitch == 0) {
                                            gStyle->SetOptFit(111);
                                        }
                                        cC->Update();
                                        gAllWidthP->Draw("AP");
                                        if (BellSwitch == 0) {
                                            gAllWidthP->Fit("fit1","B");
                                        }
                                        cC->Update();
                                        //fit1->DrawCopy("SAME");
                                        //cC->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/JustWidth_Pt.pdf"));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_JustWidth_Pt.pdf"));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        //new multiplot section
                                        TestAllWidth = new TMultiGraph();
                                        gAllWidthP->SetMarkerStyle(22);
                                        gAllWidthP->SetMarkerColor(2);
                                        gAllWidthP->SetMarkerSize(2.0);
                                        TestAllWidth->Add(gAllWidthP);
                                        
                                        if (BellSwitch == 1) {
                                            legOnePT= new TLegend(0.6,0.70,0.9,0.9);
                                            sprintf(AEword,"#Xi^{**}(1820)");
                                            legOnePT->AddEntry(gAllWidthP,AEword,"l");
                                            legOnePT->SetTextSize(0.03);
                                        }
                                        
                                        //new section
                                        //SYS
                                        if ((cbmin[0] == 0) && (cbmax[0] == 100)) {
                                            gAllWidthP = new TGraphErrors((MaxCen-1),Xforplot3,Yforplot3,XforplotError3,YforplotError3SYS);
                                        }
                                        else{
                                            gAllWidthP = new TGraphErrors(MaxCen,Xforplot3,Yforplot3,XforplotError3,YforplotError3SYS);
                                        }
                                        
                                        //gAllWidthP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3SYS,YforplotError3SYS);
                                        gAllWidthP->SetLineStyle(1);//type
                                        gAllWidthP->SetLineColor(2);
                                        //gAllMeanP->SetFillStyle(3001);
                                        ///gAllMeanP->SetFillColorAlpha(1+n,0.10);
                                        gAllWidthP->SetFillStyle(1);
                                        gAllWidthP->SetMarkerStyle(22);
                                        gAllWidthP->SetMarkerColor(2);
                                        gAllWidthP->SetMarkerSize(2.0);

                                        if (BellSwitch == 1) {
                                            TestAllWidth->Add(gAllWidthP,"[]");
                                            gStyle->SetEndErrorSize(4);
                                            //gStyle->SetMarkerSize(3.0);
                                        }
                                        else{
                                            TestAllWidth->Add(gAllWidthP,"5");
                                        }
                                        //SYS
                                        //new section
                                        
                                        //shading
                                        if ((cbmin[0] == 0) && (cbmax[0] == 100)) {
                                            for (y=0; y < (MaxCen-1); y++) {
                                                Yforplot3[y]=0.024;
                                                YforplotError3[y]=0.006;
                                                //YforplotError3SYS[y]=WforplotPErrorSYS[0][3][y+1];
                                            }
                                            Xforplot3[0]=0;
                                            Xforplot3[1]=50;
                                            Xforplot3[2]=100;
                                            gAllWidthP = new TGraphErrors((MaxCen-1),Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        else{
                                            for (y=0; y < MaxCen; y++) {
                                                Yforplot3[y]=0.24;
                                                YforplotError3[y]=0.006;
                                                //YforplotError3SYS[y]=WforplotPErrorSYS[0][3][y];
                                            }
                                            gAllWidthP = new TGraphErrors(MaxCen,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        }
                                        //gAllWidthP->SetFillColor(2);//type
                                        //gAllWidthP->SetFillStyle(3010);
                                        //gAllMeanP->SetFillStyle(3001);
                                        gAllWidthP->SetFillColorAlpha(2,0.10);
                                        //gAllWidthP->SetFillStyle(1);
                                        
                                        if (BellSwitch == 1) {
                                            TestAllWidth->Add(gAllWidthP,"3");
                                            sprintf(AEword,"#Xi^{**}(1820) PDG Error");
                                            legOnePT->AddEntry(gAllWidthP,AEword,"f");
                                            //gStyle->SetEndErrorSize(4);
                                        }
                                        else{
                                            //TestAllWidth->Add(gAllWidthP,"5");
                                        }
                                        //shading
                                        
                                        Yforplot3[0]=0.00856;
                                        YforplotError3[0]=0.000062;
                                        Yforplot3[1]=0.00980;
                                        YforplotError3[1]=0.000083;
                                        Yforplot3[2]=0.01056;
                                        YforplotError3[2]=0.000169;
                                        Xforplot3[0]=7.5;
                                        Xforplot3[1]=32.5;
                                        Xforplot3[2]=75;
                                        
                                        if (BellSwitch == 0) {
                                            XforplotError3[0]=7.5;
                                            XforplotError3[1]=17.5;
                                            XforplotError3[2]=25;
                                        }
                                        
                                        OtherData = new TGraphErrors(3,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        OtherData->SetLineColor(3);
                                        OtherData->SetMarkerColor(3);
                                        OtherData->SetMarkerStyle(23);
                                        OtherData->SetMarkerSize(2.0);
                                        TestAllWidth->Add(OtherData);
                                        if (BellSwitch == 1) {
                                            sprintf(AEword,"#Xi^{*}(1530) from Bong-Hwi Kim");
                                            legOnePT->AddEntry(OtherData,AEword,"l");
                                            TestAllWidth->SetMinimum(0.000);//0.008
                                            //gStyle->SetMarkerSize(3.0);
                                        }
                                        else{
                                            TestAllWidth->SetMinimum(0.008);//0.008
                                        }
                                        TestAllWidth->SetMaximum(0.035);
                                        //line
                                        Yforplot3[0]=0.0091;
                                        YforplotError3[0]=0.0005;
                                        Yforplot3[1]=0.0091;
                                        YforplotError3[1]=0.0005;
                                        Yforplot3[2]=0.0091;
                                        YforplotError3[2]=0.0005;
                                        Xforplot3[0]=0;
                                        Xforplot3[1]=50;
                                        Xforplot3[2]=100;
                                        OtherData = new TGraphErrors(3,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        //OtherData->SetLineColor(3);
                                        if (BellSwitch == 1) {
                                            OtherData->SetFillColorAlpha(3,0.10);
                                            TestAllWidth->Add(OtherData,"3");
                                            sprintf(AEword,"#Xi^{*}(1530) PDG Error");
                                            legOnePT->AddEntry(OtherData,AEword,"f");
                                        }
                                        
                                        simpleline = new TF1("simpleline","pol0",0,100);
                                        simpleline->SetLineColor(2);
                                        simpleline->SetLineStyle(2);
                                        simpleline->FixParameter(0,0.024);
                                        simpleline2 = new TF1("simpleline2","pol0",0,100);
                                        simpleline2->SetLineColor(3);
                                        simpleline2->SetLineStyle(2);
                                        simpleline2->FixParameter(0,0.0091);
                                        cC->Update();
                                        TestAllWidth->Draw("AP");
                                        TestAllWidth->SetTitle(Form("Width (FWHM) for %s 1.5<p_{t}<20.0 GeV/c",partnameComplex));
                                        TestAllWidth->GetXaxis()->SetTitle("Multiplicity class %");
                                        TestAllWidth->GetYaxis()->SetTitle("Width (GeV/c^{2})");
                                        
                                        simpleline->Draw("SAME");
                                        simpleline2->Draw("SAME");
                                        if (BellSwitch == 1) {
                                            sprintf(AEword,"#Xi^{**}(1820) PDG Value");
                                            legOnePT->AddEntry(simpleline,AEword,"l");
                                            sprintf(AEword,"#Xi^{*}(1530) PDG Value");
                                            legOnePT->AddEntry(simpleline2,AEword,"l");
                                            legOnePT->Draw("SAME");
                                        }
                                        
                                        cC->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/JustWidth2.pdf"));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_JustWidth2.pdf"));
                                             }//cut*/
                                        }//Normal
                                        
                                        cC->Clear();
                                        //simpleline = new TF1("simpleline","pol0",0,100);
                                        //simpleline->SetLineColor(2);
                                        //simpleline->SetLineStyle(2);
                                        //simpleline->FixParameter(0,0.024);
                                        //line
                                        //new multiplot section
                                        P3MultiAllWidth->Add(gAllWidthP);
                                        gStyle->SetOptFit(0);
                                        
                                    }
                                }//use 1 pt
                                else{
                                    for (int loop3=0; loop3 < 4; loop3++) {
                                        Xforplot3[loop3]=((ptmin[loop3]+ptmax[loop3])/2.0);
                                        if (BellSwitch == 1) {
                                            XforplotError3[loop3]=0;
                                            XforplotError3SYS[loop3]=0;
                                        }
                                        else{
                                            XforplotError3[loop3]=((ptmax[loop3]-ptmin[loop3])/2.0);
                                            XforplotError3SYS[loop3]=0.5;
                                        }
                                    }
                                    
                                    if (BellSwitch == 1) {
                                        Xforplot3[0]=1.99;
                                        Xforplot3[1]=3.10;
                                        Xforplot3[2]=5.07;
                                        Xforplot3[3]=2.83;
                                    }
                                    
                                    for (int y=0; y < 4; y++) {
                                        Yforplot3[y]=(YforplotP[0][y][n]/AEarray[n][y])/(NumberOfEvents*dy*(ptmax[y]-ptmin[y]));//yield Vs. pt
                                        YforplotError3[y]=(Yforplot3[y])*(sqrt(pow((YforplotPError[0][y][n]/YforplotP[0][y][n]),2)+pow((AEarrayerror[n][y]/AEarray[n][y]),2)))/(NumberOfEvents*dy*(ptmax[y]-ptmin[y]));//yield Vs. pt
                                    }
                                    gAllYieldP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                    
                                    if (BellSwitch == 1) {
                                        gAllYieldP->SetLineStyle(1);//2//type
                                    }
                                    else{
                                        gAllYieldP->SetLineStyle(2);//type
                                    }
                                    gAllYieldP->SetLineColor(1+n);
                                    //new to remove
                                    if (((n < 2) && (BellSwitch == 1)) || (BellSwitch == 0)) {//remove 10 and 30
                                        P3MultiAllYield->Add(gAllYieldP);
                                        
                                        if (only1flag == 1) {
                                            if (n == 0) {
                                                legYield= new TLegend(0.7,0.7,0.9,0.9);
                                            }
                                        }
                                        
                                        if ((n == 0) && (cbmax[0] == 100)) {
                                            sprintf(AEword,"%s (Min.Bias)",partnameComplex);//same
                                        }
                                        else{
                                            sprintf(AEword,"%s (%1.1f-%1.1f %%)",partnameComplex,cbmin[n],cbmax[n]);//same
                                        }
                                        legYield->AddEntry(gAllYieldP,AEword,"l");
                                    }
                                    //new to remove
                                    
                                    for (int y=0; y < 4; y++) {
                                        Yforplot3[y]=MforplotP[0][y][n];//Mean Vs. pt
                                        YforplotError3[y]=MforplotPError[0][y][n];//mean Vs pt
                                        YforplotError3SYS[y]=MforplotPErrorSYS[0][y][n];
                                    }
                                    gAllMeanP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                    if (BellSwitch == 1) {
                                        gAllMeanP->SetLineStyle(1);//2//type
                                    }
                                    else{
                                        gAllMeanP->SetLineStyle(2);//type
                                    }
                                    gAllMeanP->SetLineColor(1+n);
                                    //new to remove
                                    if (((n < 2) && (BellSwitch == 1)) || (BellSwitch == 0)) {//remove 10 and 30
                                        P3MultiAllMean->Add(gAllMeanP);
                                        
                                        if (only1flag == 1) {
                                            if (n == 0) {
                                                legMean= new TLegend(0.7,0.7,0.9,0.9);
                                            }
                                        }
                                        
                                        if ((n == 0) && (cbmax[0] == 100)) {
                                            sprintf(AEword,"%s (Min.Bias)",partnameComplex);
                                        }
                                        else{
                                            sprintf(AEword,"%s (%1.1f-%1.1f %%)",partnameComplex,cbmin[n],cbmax[n]);
                                        }
                                        legMean->AddEntry(gAllMeanP,AEword,"l");
                                        
                                        //SYS
                                        gAllMeanP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3SYS,YforplotError3SYS);
                                        if (BellSwitch == 1) {
                                            gAllMeanP->SetLineStyle(1);//2//type
                                        }
                                        else{
                                            gAllMeanP->SetLineStyle(2);//type
                                        }
                                        gAllMeanP->SetLineColor(1+n);
                                        //gAllMeanP->SetDrawOption("2");
                                        //gAllMeanP->SetFillStyle(3001);
                                        ///gAllMeanP->SetFillColorAlpha(1+n,0.10);
                                        gAllMeanP->SetFillStyle(1);
                                        if (BellSwitch == 1) {
                                            P3MultiAllMean->Add(gAllMeanP,"[]");
                                            gStyle->SetEndErrorSize(4);
                                        }
                                        else{
                                            P3MultiAllMean->Add(gAllMeanP,"5");
                                        }
                                        //SYS
                                    }
                                    //new to remove
                                    
                                    for (int y=0; y < 4; y++) {
                                        Yforplot3[y]=WforplotP[0][y][n];//width Vs. pt
                                        YforplotError3[y]=WforplotPError[0][y][n];//width Vs. pt
                                        YforplotError3SYS[y]=WforplotPErrorSYS[0][y][n];
                                    }
                                    gAllWidthP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                    if (BellSwitch == 1) {
                                        for (int y=0; y < 3; y++) {
                                            Yforplot3[y]=WforplotP[0][y][n];//width Vs. pt
                                            YforplotError3[y]=WforplotPError[0][y][n];//width Vs. pt
                                            YforplotError3SYS[y]=WforplotPErrorSYS[0][y][n];
                                        }
                                        gAllWidthP = new TGraphErrors(3,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                    }
                                    gAllWidthP->SetLineStyle(1);//2//type
                                    gAllWidthP->SetLineColor(1+n);
                                    //new to remove
                                    if (((n < 2) && (BellSwitch == 1)) || (BellSwitch == 0)) {//remove 10 and 30
                                        P3MultiAllWidth->Add(gAllWidthP);
                                        
                                        if (only1flag == 1) {
                                            if (n == 0) {
                                                if (BellSwitch == 1) {
                                                    legWidth= new TLegend(0.5,0.15,0.9,0.35);
                                                    legWidth->SetTextSize(0.03);
                                                }
                                                else{
                                                    legWidth= new TLegend(0.7,0.7,0.9,0.9);
                                                }
                                            }
                                        }
                                        
                                        if ((n == 0) && (cbmax[0] == 100)) {
                                            sprintf(AEword,"%s (Min.Bias)",partnameComplex);//same
                                        }
                                        else{
                                            sprintf(AEword,"%s (%1.1f-%1.1f %%)",partnameComplex,cbmin[n],cbmax[n]);//same
                                        }
                                        legWidth->AddEntry(gAllWidthP,AEword,"l");
                                        
                                        //SYS
                                        gAllWidthP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3SYS,YforplotError3SYS);
                                        if (BellSwitch == 1) {
                                            gAllWidthP = new TGraphErrors(3,Xforplot3,Yforplot3,XforplotError3SYS,YforplotError3SYS);
                                            gAllWidthP->SetLineStyle(1);//2//type
                                        }
                                        else{
                                            gAllWidthP->SetLineStyle(2);//type
                                        }
                                        gAllWidthP->SetLineColor(1+n);
                                        //gAllWidthP->SetDrawOption("2");
                                        //gAllWidthP->SetFillStyle(3001);
                                        ///gAllWidthP->SetFillColorAlpha(1+n,0.10);
                                        gAllWidthP->SetFillStyle(1);
                                        if (BellSwitch == 1) {
                                            P3MultiAllWidth->Add(gAllWidthP,"[]");
                                            gStyle->SetEndErrorSize(4);
                                        }
                                        else{
                                            P3MultiAllWidth->Add(gAllWidthP,"5");
                                        }
                                        //sprintf(AEword,"%s systematic %1.1f-%1.1f %% Cen",partname,cbmin[n],cbmax[n]);
                                        //legWidth->AddEntry(gAllWidthP,AEword,"f");
                                        
                                        //SYS
                                    }
                                    //new to remove
                                    
                                    if (n == (MaxCen-1)) {
                                        cC->cd();
                                        gAllWidthP->Draw("AP");
                                        cC->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/JustWidth2.pdf"));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_JustWidth2.pdf"));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        cC->Update();
                                    }
                                    
                                    if (n == (MaxCen-1)) {//recent
                                        clog->cd();
                                        clog->SetLogy();
                                        clog->Update();
                                        if (BellSwitch == 1) {
                                            P3MultiAllYield->DrawClone("AP");
                                        }
                                        else{
                                            P3MultiAllYield->DrawClone("APL");
                                        }
                                        legYield->Draw("SAME");
                                        clog->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            clog->Print(Form("plotsCompleteV3/AllYield%s_Pt.pdf",partname));
                                            /*}//no cut
                                             else{//cut
                                             clog->Print(Form("plotsCompleteV3/cut1_AllYield%s_Pt.pdf",partname));
                                             }//cut*/
                                        }//Normal
                                        clog->Clear();
                                        clog->Update();
                                        
                                        cC->cd();
                                        cC->Update();
                                        P3MultiAllMean->DrawClone("AP");
                                        legMean->Draw("SAME");
                                        cC->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/AllMean%s_Pt.pdf",partname));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_AllMean%s_Pt.pdf",partname));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        cC->Update();
                                        
                                        cC->cd();
                                        cC->Update();
                                        P3MultiAllWidth->DrawClone("AP");
                                        legWidth->Draw("SAME");
                                        cC->Update();
                                        if (hmtloop == 0) {//Normal
                                            //if (cutloop == 0) {//no cut
                                            cC->Print(Form("plotsCompleteV3/AllWidth%s_Pt.pdf",partname));
                                            /*}//no cut
                                             else{//cut
                                             cC->Print(Form("plotsCompleteV3/cut1_AllWidth%s_Pt.pdf",partname));
                                             }//cut*/
                                        }//Normal
                                        cC->Clear();
                                        cC->Update();
                                    }//recent
                                    
                                    //average code
                                    if (only1flag != 1) {
                                        for (int loop3=0; loop3 < 4; loop3++) {
                                            Xforplot3[loop3]=((ptmin[loop3]+ptmax[loop3])/2.0);
                                            if (BellSwitch == 1) {
                                                XforplotError3[loop3]=0;
                                                XforplotError3SYS[loop3]=0.0;
                                            }
                                            else{
                                                XforplotError3[loop3]=((ptmax[loop3]-ptmin[loop3])/2.0);
                                                XforplotError3SYS[loop3]=0.5;
                                            }
                                        }
                                        
                                        for (int y=0; y < 4; y++) {
                                            Yforplot3[y]=(((YforplotP[0][y][n]/AEarray[n][y])+(YforplotP[1][y][n]/AEarray[n][y]))/2.0)/(NumberOfEvents*dy*(ptmax[y]-ptmin[y]));//yield Vs. pt
                                            YforplotError3[y]=sqrt(0.25*pow((((YforplotP[0][y][n]/AEarray[n][y]))*(sqrt(pow((YforplotPError[0][y][n]/YforplotP[0][y][n]),2)+pow((AEarrayerror[n][y]/AEarray[n][y]),2)))),2) + 0.25*pow((((YforplotP[1][y][n]/AEarray[n][y]))*(sqrt(pow((YforplotPError[1][y][n]/YforplotP[1][y][n]),2)+pow((AEarrayerror[n][y]/AEarray[n][y]),2)))),2) + pow((YforplotP[0][y][n]/AEarray[n][y])-(YforplotP[1][y][n]/AEarray[n][y]),2))/(NumberOfEvents*dy*(ptmax[y]-ptmin[y]));//yield Vs. pt
                                        }
                                        gAllYieldP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        if (BellSwitch == 1) {
                                            gAllYieldP->SetLineStyle(1);//2//type
                                        }
                                        else{
                                            gAllYieldP->SetLineStyle(2);//type
                                        }
                                        gAllYieldP->SetLineColor(1+n);
                                        MultiYieldAverage->Add(gAllYieldP);
                                        
                                        if (n == 0) {
                                            legYieldAverage= new TLegend(0.7,0.7,0.9,0.9);
                                        }
                                        sprintf(AEword,"%s average %1.1f-%1.1f %% Cen",partname,cbmin[n],cbmax[n]);
                                        legYieldAverage->AddEntry(gAllYieldP,AEword,"l");
                                        
                                        for (int y=0; y < 4; y++) {
                                            Yforplot3[y]=(MforplotP[0][y][n]+MforplotP[1][y][n])/2.0;//Mean Vs. pt
                                            
                                            //YforplotError3[y]=sqrt(0.25*pow(MforplotPError[0][y][n],2) + 0.25*pow(MforplotPError[1][y][n],2) + pow(MforplotP[0][y][n]-MforplotP[1][y][n],2));//mean Vs pt
                                            YforplotError3[y]=sqrt(0.25*pow(MforplotPError[0][y][n],2) + 0.25*pow(MforplotPError[1][y][n],2));//mean Vs pt
                                            YforplotError3SYS[y]=sqrt(0.25*pow(MforplotPErrorSYS[0][y][n],2) + 0.25*pow(MforplotPErrorSYS[1][y][n],2) + pow(MforplotP[0][y][n]-MforplotP[1][y][n],2));
                                        }
                                        gAllMeanP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        gAllMeanP->SetLineStyle(1);//2//type
                                        gAllMeanP->SetLineColor(1+n);
                                        MultiMeanAverage->Add(gAllMeanP);
                                        if (n == 0) {
                                            legMeanAverage= new TLegend(0.7,0.7,0.9,0.9);
                                        }
                                        sprintf(AEword,"%s average %1.1f-%1.1f %% Cen",partname,cbmin[n],cbmax[n]);
                                        legMeanAverage->AddEntry(gAllMeanP,AEword,"l");
                                        
                                        //SYS
                                        gAllMeanP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3SYS,YforplotError3SYS);
                                        if (BellSwitch == 1) {
                                            gAllMeanP->SetLineStyle(1);//2//type
                                        }
                                        else{
                                            gAllMeanP->SetLineStyle(2);//type
                                        }
                                        gAllMeanP->SetLineColor(1+n);
                                        //gAllMeanP->SetDrawOption("2");
                                        //gAllMeanP->SetFillStyle(3001);
                                        ///gAllMeanP->SetFillColorAlpha(1+n,0.10);
                                        gAllMeanP->SetFillStyle(1);
                                        if (BellSwitch == 1) {
                                            MultiMeanAverage->Add(gAllMeanP,"[]");
                                            gStyle->SetEndErrorSize(4);
                                        }
                                        else{
                                            MultiMeanAverage->Add(gAllMeanP,"5");
                                        }
                                        //SYS
                                        
                                        for (int y=0; y < 4; y++) {
                                            Yforplot3[y]=(WforplotP[0][y][n]+WforplotP[1][y][n])/2.0;//Mean Vs. pt
                                            //YforplotError3[y]=sqrt(0.25*pow(WforplotPError[0][y][n],2) + 0.25*pow(WforplotPError[1][y][n],2) + pow(WforplotP[0][y][n]-WforplotP[1][y][n],2));//mean Vs pt
                                            YforplotError3[y]=sqrt(0.25*pow(WforplotPError[0][y][n],2) + 0.25*pow(WforplotPError[1][y][n],2));//mean Vs pt
                                            YforplotError3SYS[y]=sqrt(0.25*pow(WforplotPErrorSYS[0][y][n],2) + 0.25*pow(WforplotPErrorSYS[1][y][n],2) + pow(WforplotP[0][y][n]-WforplotP[1][y][n],2));
                                        }
                                        gAllWidthP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3,YforplotError3);
                                        if (BellSwitch == 1) {
                                            gAllWidthP->SetLineStyle(1);//2//type
                                        }
                                        else{
                                            gAllWidthP->SetLineStyle(2);//type
                                        }
                                        gAllWidthP->SetLineColor(1+n);
                                        MultiWidthAverage->Add(gAllWidthP);
                                        if (n == 0) {
                                            legWidthAverage= new TLegend(0.7,0.7,0.9,0.9);
                                        }
                                        sprintf(AEword,"%s average %1.1f-%1.1f %%",partname,cbmin[n],cbmax[n]);
                                        legWidthAverage->AddEntry(gAllWidthP,AEword,"l");
                                        
                                        //SYS
                                        gAllWidthP = new TGraphErrors(4,Xforplot3,Yforplot3,XforplotError3SYS,YforplotError3SYS);
                                        if (BellSwitch == 1) {
                                            gAllWidthP->SetLineStyle(1);//2//type
                                        }
                                        else{
                                            gAllWidthP->SetLineStyle(2);//type
                                        }
                                        gAllWidthP->SetLineColor(1+n);
                                        //gAllWidthP->SetDrawOption("2");
                                        //gAllWidthP->SetFillStyle(3001);
                                        ///gAllWidthP->SetFillColorAlpha(1+n,0.10);
                                        gAllWidthP->SetFillStyle(1);
                                        if (BellSwitch == 1) {
                                            MultiWidthAverage->Add(gAllWidthP,"[]");
                                            gStyle->SetEndErrorSize(4);
                                        }
                                        else{
                                            MultiWidthAverage->Add(gAllWidthP,"5");
                                        }
                                        //SYS
                                        
                                        if (n == (MaxCen-1)) {//recent
                                            clog->cd();
                                            clog->SetLogy();
                                            clog->Update();
                                            MultiYieldAverage->DrawClone("AP");
                                            legYieldAverage->Draw("SAME");
                                            clog->Update();
                                            if (hmtloop == 0) {//Normal
                                                //if (cutloop == 0) {//no cut
                                                clog->Print(Form("plotsCompleteV3/AllYieldAverage_Pt.pdf"));
                                                /*}//no cut
                                                 else{//cut
                                                 clog->Print(Form("plotsCompleteV3/cut1_AllYieldAverage_Pt.pdf"));
                                                 }//cut*/
                                            }//Normal
                                            clog->Clear();
                                            clog->Update();
                                            
                                            cC->cd();
                                            cC->Update();
                                            MultiMeanAverage->DrawClone("AP");
                                            legMeanAverage->Draw("SAME");
                                            cC->Update();
                                            if (hmtloop == 0) {//Normal
                                                //if (cutloop == 0) {//no cut
                                                cC->Print(Form("plotsCompleteV3/AllMeanAverage_Pt.pdf"));
                                                /*}//no cut
                                                 else{//cut
                                                 cC->Print(Form("plotsCompleteV3/cut1_AllMeanAverage_Pt.pdf"));
                                                 }//cut*/
                                            }//Normal
                                            cC->Clear();
                                            cC->Update();
                                            
                                            cC->cd();
                                            cC->Update();
                                            MultiWidthAverage->DrawClone("AP");
                                            legWidthAverage->Draw("SAME");
                                            cC->Update();
                                            if (hmtloop == 0) {//Normal
                                                //if (cutloop == 0) {//no cut
                                                cC->Print(Form("plotsCompleteV3/AllWidthAverage_Pt.pdf"));
                                                /*}//no cut
                                                 else{//cut
                                                 cC->Print(Form("plotsCompleteV3/cut1_AllWidthAverage_Pt.pdf"));
                                                 }//cut*/
                                            }//Normal
                                            cC->Clear();
                                            cC->Update();
                                        }//recent
                                    }
                                    
                                    //average code
                                    
                                    //
                                }//not 1 pt
                                
                                /*//not used much anymore
                                 if (Use1PT != 1) {//not 1 pt
                                 //
                                 cC->cd();
                                 cC->Clear();
                                 cC->Update();
                                 SystematicAllYield->SetTitle(Form("Yields for %s GM %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                                 SystematicAllYield->Draw();
                                 cC->Update();
                                 if (hmtloop == 0) {//Normal
                                 //if (cutloop == 0) {//no cut
                                 cC->Print(Form("plotsCompleteV3/%s/SYSYield_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//no cut
                                 //else{//cut
                                 //    cC->Print(Form("plotsCompleteV3/%s/cut1_SYSYield_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//cut
                                 }//Normal
                                 cC->Clear();
                                 cC->Update();
                                 SystematicAllMean->SetTitle(Form("Means for %s GM %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                                 SystematicAllMean->Draw();
                                 cC->Update();
                                 if (hmtloop == 0) {//Normal
                                 //if (cutloop == 0) {//no cut
                                 cC->Print(Form("plotsCompleteV3/%s/SYSMean_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//no cut
                                 //else{//cut
                                 //    cC->Print(Form("plotsCompleteV3/%s/cut1_SYSMean_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//cut
                                 }//Normal
                                 cC->Clear();
                                 cC->Update();
                                 SystematicAllWidth->SetTitle(Form("Widths for %s GM %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                                 SystematicAllWidth->Draw();
                                 cC->Update();
                                 if (hmtloop == 0) {//Normal
                                 //if (cutloop == 0) {//no cut
                                 cC->Print(Form("plotsCompleteV3/%s/SYSWidth_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//no cut
                                 //else{//cut
                                 //    cC->Print(Form("plotsCompleteV3/%s/cut1_SYSWidth_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//cut
                                 }//Normal
                                 cC->Clear();
                                 cC->Update();
                                 StatisticalAllMean->SetTitle(Form("Statistical Means for %s GM %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                                 StatisticalAllMean->Draw();
                                 cC->Update();
                                 if (hmtloop == 0) {//Normal
                                 //if (cutloop == 0) {//no cut
                                 cC->Print(Form("plotsCompleteV3/%s/STAMean_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//no cut
                                 //else{//cut
                                 //    cC->Print(Form("plotsCompleteV3/%s/cut1_STAMean_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//cut
                                 }//Normal
                                 cC->Clear();
                                 cC->Update();
                                 StatisticalAllWidth->SetTitle(Form("Statistical Widths for%s GM %1.1f<Cen<%1.1f vs p_{t}",DC1,cbmin[n],cbmax[n]));
                                 StatisticalAllWidth->Draw();
                                 cC->Update();
                                 if (hmtloop == 0) {//Normal
                                 //if (cutloop == 0) {//no cut
                                 cC->Print(Form("plotsCompleteV3/%s/STAWidth_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//no cut
                                 //else{//cut
                                 //    cC->Print(Form("plotsCompleteV3/%s/cut1_STAWidth_NR%i_GM_c%i.pdf",partname,m,n));
                                 //}//cut
                                 }//Normal
                                 cC->Clear();
                                 if (n == 0) {
                                 cSYSYGM->Clear();
                                 cSYSYGM->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                 cSYSYGM->Update();
                                 cSYSMGM->Clear();
                                 cSYSMGM->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                 cSYSMGM->Update();
                                 cSYSWGM->Clear();
                                 cSYSWGM->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                 cSYSWGM->Update();
                                 cSTAMGM->Clear();
                                 cSTAMGM->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                 cSTAMGM->Update();
                                 cSTAWGM->Clear();
                                 cSTAWGM->Divide(MaxCen/2.0,MaxCen/2.0);//2,2
                                 cSTAWGM->Update();
                                 }
                                 cSYSYGM->cd(n+1);
                                 cSYSYGM->Update();
                                 SystematicAllYield->Draw();
                                 cSYSYGM->Update();
                                 cSYSMGM->cd(n+1);
                                 cSYSMGM->Update();
                                 SystematicAllMean->Draw();
                                 cSYSMGM->Update();
                                 cSYSWGM->cd(n+1);
                                 cSYSWGM->Update();
                                 SystematicAllWidth->Draw();
                                 cSYSWGM->Update();
                                 cSTAMGM->cd(n+1);
                                 cSTAMGM->Update();
                                 StatisticalAllMean->Draw();
                                 cSTAMGM->Update();
                                 cSTAWGM->cd(n+1);
                                 cSTAWGM->Update();
                                 StatisticalAllWidth->Draw();
                                 cSTAWGM->Update();
                                 if (n == (MaxCen-1)) {
                                 if (hmtloop == 0) {//Normal
                                 //if (cutloop == 0) {//no cut
                                 cSYSYGM->Print(Form("plotsCompleteV3/%s/TogetherSYSYield_GM_NR%i.pdf",partname,m));
                                 cSYSMGM->Print(Form("plotsCompleteV3/%s/TogetherSYSMean_GM_NR%i.pdf",partname,m));
                                 cSYSWGM->Print(Form("plotsCompleteV3/%s/TogetherSYSWidth_GM_NR%i.pdf",partname,m));
                                 cSTAMGM->Print(Form("plotsCompleteV3/%s/TogetherSTAMean_GM_NR%i.pdf",partname,m));
                                 cSTAWGM->Print(Form("plotsCompleteV3/%s/TogetherSTAWidth_GM_NR%i.pdf",partname,m));
                                 //}//no cut
                                 //else{//cut
                                 //    cSYSYGM->Print(Form("plotsCompleteV3/%s/cut1_TogetherSYSYield_GM_NR%i.pdf",partname,m));
                                 //    cSYSMGM->Print(Form("plotsCompleteV3/%s/cut1_TogetherSYSMean_GM_NR%i.pdf",partname,m));
                                 //    cSYSWGM->Print(Form("plotsCompleteV3/%s/cut1_TogetherSYSWidth_GM_NR%i.pdf",partname,m));
                                 //    cSTAMGM->Print(Form("plotsCompleteV3/%s/cut1_TogetherSTAMean_GM_NR%i.pdf",partname,m));
                                 //    cSTAWGM->Print(Form("plotsCompleteV3/%s/cut1_TogetherSTAWidth_GM_NR%i.pdf",partname,m));
                                 //}//cut
                                 }//Normal
                                 
                                 }
                                 }//not 1 pt
                                 */
                                //not used much anymore
                                
                            }//if particleloop == 0
                            
                            //}//end of jNR
                            
                        }//endof cen loop
                    }//end of particleloop <= 1
                }//cutloop == cutflag
                
                
                /*
                 //find the lowest chi2/ndf
                 cout << "*************************************************************************************************************" << endl;
                 cout << "********************************************Lowest chi2/ndf****************************************************" << endl;
                 for (int g=0; g <= 3; g++) {//Cen
                 for (int i=0; i < (MaxPT-MinPT); i++) {//pt
                 double LOWEST=0;
                 int flag1=0;
                 int TestG=g;
                 int TestI=i;
                 int TestJ=0;
                 int TestK=0;
                 int TestL=0;
                 int TestM=0;
                 for (int j=0; j <= 0; j++) {//NR
                 for (int k=1; k <= 3; k++) {//range k=0; k <= 2
                 for (int l=0; l <= 2 ; l++) {//version
                 for (int m=0; m <= 1; m++) {//type
                 if ((LOWEST == 0) && (flag1 == 0)) {
                 LOWEST=ALLChi2[g][i][j][k][l][m];
                 flag1=1;
                 }
                 else if((TMath::Abs(ALLChi2[g][i][j][k][l][m]-1.0) <= TMath::Abs(LOWEST-1.0))){
                 if((TMath::Abs(ALLChi2[g][i][j][k][l][m]-1.0) == TMath::Abs(LOWEST-1.0)) && (LOWEST != 0)){
                 cout << "Exact same chi2/ndf Cen " << TestG << " NR " << TestJ << " range " << TestK << " version " << TestL << " type " << TestM << " chi2/ndf = " << LOWEST << endl;
                 }
                 //cout << "ALLChi2 " << ALLChi2[i][j][k][l][m] << " LOWEST before " << LOWEST << " ABS chi2-1 " << TMath::Abs(ALLChi2[i][j][k][l][m]-1.0) << " ABS Lowest-1.0 " << TMath::Abs(LOWEST-1.0) << endl;
                 LOWEST=ALLChi2[g][i][j][k][l][m];
                 TestG=g;
                 TestJ=j;
                 TestK=k;
                 TestL=l;
                 TestM=m;
                 }
                 }
                 }
                 }
                 }
                 cout << "pt bin i " << TestI+MinPT << " Cen " << TestG << " NR " << TestJ << " range " << TestK << " version " << TestL << " type " << TestM << " chi2/ndf = " << LOWEST << endl;
                 cout << "Associated Yield " << ALLYield[TestG][TestI][TestJ][TestK][TestL][TestM] << " Mean " << ALLMean[TestG][TestI][TestJ][TestK][TestL][TestM] << " width " << ALLWidth[TestG][TestI][TestJ][TestK][TestL][TestM] << endl;
                 }
                 }//end of chi2 cen loop
                 if (particleloop == 0) {//gm chi2/ndf
                 cout << "*************************************************************************************************************" << endl;
                 cout << "********************************************Lowest chi2/ndf GM****************************************************" << endl;
                 for (int g=0; g <= 3; g++) {//Cen
                 for (int i=0; i < (MaxPT-MinPT); i++) {//pt
                 double LOWEST=0;
                 int flag1=0;
                 int TestG=g;
                 int TestI=i;
                 int TestJ=0;
                 int TestK=0;
                 int TestL=0;
                 int TestM=0;
                 for (int j=0; j <= 0; j++) {//NR
                 for (int k=1; k <= 3; k++) {//range k=0; k <= 2
                 for (int l=0; l <= 2 ; l++) {//version
                 for (int m=0; m <= 1; m++) {//type
                 if ((LOWEST == 0) && (flag1 == 0)) {
                 LOWEST=GMALLChi2[g][i][j][k][l][m];
                 flag1=1;
                 }
                 else if((TMath::Abs(GMALLChi2[g][i][j][k][l][m]-1.0) <= TMath::Abs(LOWEST-1.0))){
                 if((TMath::Abs(GMALLChi2[g][i][j][k][l][m]-1.0) == TMath::Abs(LOWEST-1.0)) && (LOWEST != 0)){
                 cout << "Exact same chi2/ndf GM Cen " << TestG << " NR " << TestJ << " range " << TestK << " version " << TestL << " type " << TestM << " chi2/ndf = " << LOWEST << endl;
                 }
                 //cout << "ALLChi2 " << ALLChi2[i][j][k][l][m] << " LOWEST before " << LOWEST << " ABS chi2-1 " << TMath::Abs(ALLChi2[i][j][k][l][m]-1.0) << " ABS Lowest-1.0 " << TMath::Abs(LOWEST-1.0) << endl;
                 LOWEST=GMALLChi2[g][i][j][k][l][m];
                 TestG=g;
                 TestJ=j;
                 TestK=k;
                 TestL=l;
                 TestM=m;
                 }
                 }
                 }
                 }
                 }
                 cout << "GM pt bin i " << TestI+MinPT << " Cen " << TestG << " NR " << TestJ << " range " << TestK << " version " << TestL << " type " << TestM << " chi2/ndf = " << LOWEST << endl;
                 cout << "GM Associated Yield " << GMALLYield[TestG][TestI][TestJ][TestK][TestL][TestM] << " Mean " << GMALLMean[TestG][TestI][TestJ][TestK][TestL][TestM] << " width " << GMALLWidth[TestG][TestI][TestJ][TestK][TestL][TestM] << endl;
                 }
                 }//end of chi2 cen loop
                 }//gm chi2/ndf
                 */
                
                
                if (writetxt == 1) {//begin output text file
                    if ((particleloop == 0) && (hmtloop == 0) /* && (cutloop == 0)*/) {
                        outputFile.open("plotsCompleteV3/DataValues.txt", std::ios::app);
                    }
                    outputFile << "particle | hmt/cut1 | cb | pt | entries | yield | mean +- error | width +- error | chi2 " << endl;
                    
                    if (particleloop == 0) {
                        outputFile << "Lambda+KX |";
                    }
                    else if(particleloop == 1){
                        outputFile << "Lambda+K0 |";
                    }
                    if (hmtloop == 0) {
                        outputFile << " Normal   |";
                    }
                    else{
                        outputFile << " hmt      |";
                    }
                    
                    for (g=0; g <= 3; g++) {//Cen
                        if ((hmtloop == 1) && (g >= 1)) {
                            continue;
                        }
                        
                        if (g == 0) {
                            outputFile << cbmin[g] << "-" << cbmax[g] << "% |";
                        }
                        else{
                            outputFile << "          |" << cbmin[g] << "-" << cbmax[g] << "% |";
                        }
                        for (i=0; i <= (npb-MinPT-1); i++) {//pt
                            if (i == 0) {
                                outputFile << ptmin[i+MinPT] << "<p_{t}<" << ptmax[i+MinPT] << " |";
                            }
                            else{
                                outputFile << "          |       |" << ptmin[i] << "<p_{t}<" << ptmax[i] << " |";
                            }
                            for (j=0; j <= (nNR-NRBegin-1); j++) {//NR //j <= 0
                                for (k=0; k <= 0; k++) {//range
                                    if ((BellSwitch == 1) && (k != 0)) {
                                        continue;
                                    }
                                    for (l=1; l <= 3 ; l++) {//version l=0; l <= 2
                                        for (m=0; m <= 1+(4.0*Nsigmaflag); m++) {//type
                                            //if (Use1Type == 1) {//use 1 type
                                            outputFile << ALLYield[g][i][j][k][l][m] << " |" << ALLMean[g][i][j][k][l][m] << " +- " << ALLMeanError[g][i][j][k][l][m] << " |" << ALLWidth[g][i][j][k][l][m] << " +- " << ALLWidthError[g][i][j][k][l][m] << endl;
                                            if (particleloop == 0) {//GM
                                                outputFile << "SameEvent |       |              |" << GMALLYield[g][i][j][k][l][m] << " |" << GMALLMean[g][i][j][k][l][m] << " +- " << GMALLMeanError[g][i][j][k][l][m] << " |" << GMALLWidth[g][i][j][k][l][m] << " +- " << GMALLWidthError[g][i][j][k][l][m] << endl;
                                            }//GM
                                            //}//use 1 type
                                        }//type
                                    }//version
                                }//range
                            }//NR
                        }//pt
                    }//Cen
                }//end output text file
                f->Close();
                if (MCflag == 1) {
                    f2->Close();
                }
                fevents->Close();
                //}//end of hmt loop
            }//end of cut loop
        }//end particle loop
    }//end hmt loop
    
    
    /*if (particleloop <= 2) {
     c2->cd();
     c2->Clear();
     c2->Update();
     
     AllMean->SetTitle(Form("All-mean C %i Pt %i",cb,pb));
     AllMean->DrawClone("AP");
     c2->Update();
     c2->Print(Form("plotsCompleteV3/all_mean_c%i_p%i.pdf",cb,pb));
     c2->Clear();
     AllWidth->SetTitle(Form("All-width C %i Pt %i",cb,pb));
     AllWidth->DrawClone("AP");
     c2->Update();
     c2->Print(Form("plotsCompleteV3/all_width_c%i_p%i.pdf",cb,pb));
     c2->Clear();
     AllIntegral->SetTitle(Form("All-integral C %i Pt %i",cb,pb));
     AllIntegral->DrawClone("AP");
     c2->Update();
     c2->Print(Form("plotsCompleteV3/all_integral_c%i_p%i.pdf",cb,pb));
     c2->Clear();
     }*/
    
    if (writetxt == 1) {
        outputFileSTAT.close();
        //Cutoftxt->close();
        outputFile.close();
    }
    
    c0->Close();
    c1->Close();
    c2->Close();
    cA->Close();
    cB->Close();
    cC->Close();
    cD->Close();
    cTotalgm->Close();
    cTotal0->Close();
    cTotal1->Close();
    cTotal2->Close();
    cTotal3->Close();
    cTotalunlike->Close();
    cLinearMix->Close();
    cLinearGM->Close();
    cSubMix->Close();
    cSubGM->Close();
    cSYSY->Close();
    cSYSM->Close();
    cSYSW->Close();
    cSTAM->Close();
    cSTAW->Close();
    cSYSYGM->Close();
    cSYSMGM->Close();
    cSYSWGM->Close();
    cSTAMGM->Close();
    cSTAWGM->Close();
    c->Close();
    clog->Close();
    cMass->Close();
    cMass2->Close();
    cout << "the end" << endl;
    return;
}

