#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TDatabasePDG.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TTree.h>
#include <vector>

#include "AliHFMassFitterVAR.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliVertexingHFUtils.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliEventPlaneResolutionHandler.h"
#include "AliVertexingHFUtils.h"

#endif

//methods for the analysis of AliAnalysisTaskSEHFv2/vn output in case of kEvShape method
//Author: Fabrizio Grosa, INFN Turin grosa@to.infn.it

//*******************************************************//
//                                                       //
//      Main Function: DmesonsFlowEvShapeYieldSyst()      //
//                                                       //
//*******************************************************//

//_____________________________________________________________________________________________
//GLOBAL VARIABLES
//to be set
//input file name
const TString infilename = "$HOME/cernbox/ALICE_WORK/Files/Trains/Run2/LHC15o/AnalysisResults_v2_3050_step2_EvShape_VZERO_CentAxis.root";
const TString suffix = "_Topod0Cut_QoverM_VZERO_EvShape";
const TString partname="Dplus";

const Int_t ncentbins=20;
const Int_t minCent=30;
const Int_t maxCent=50;

const TString reffilename="Cent3050/v2/EvShape/CutCentDep/v2Output_30_50_InOut_Topod0Cut_QoverM_VZERO_EvShape_q2Small60percVsCent_q2Large20percVsCent.root";

const TString outputdir="Cent3050/v2/EvShape/CutCentDep/RawYieldSyst";

//ptbins of the analysis
const Int_t nPtBins=4;
const Int_t nPtLims=nPtBins+1;
//const Double_t PtLims[nPtLims] = {2.,3.,4.,5.,6.,8.,12.,16.,24.};
const Double_t PtLims[nPtLims] = {3.,4.,6.,8.,12.};

//phi bins
const Int_t nPhiBins=4;
const Int_t nPhiLims=nPhiBins+1;
const Double_t PhiLims[nPhiLims]={0,TMath::Pi()/4,TMath::Pi()/2,3*TMath::Pi()/4,TMath::Pi()};

//q2 cut values (absolute cut)
const Double_t q2smalllimit=2.2;
const Double_t q2largelimit=3.2;
//percentage of events with smaller/larger q2
const Double_t q2smallpercevents=0.60;
const Double_t q2largepercevents=0.20;

//EP resolution
const Int_t evPlane=AliEventPlaneResolutionHandler::kVZERO;
Int_t evPlaneRes=AliEventPlaneResolutionHandler::kThreeSub;
const Bool_t useNcollWeight=kFALSE;
const Bool_t useAliHandlerForRes=kFALSE;

// mass fit configuration
const Int_t nReb=5;
const Int_t rebin[nReb]={2,3,4,5,6};
const Int_t nBkgFcn=3;
const Int_t typeb[nBkgFcn]={AliHFMassFitter::kExpo,AliHFMassFitter::kLin,AliHFMassFitter::kPol2};
enum {kGaus=0, kDoubleGaus, kReflTempl};
const Bool_t useTemplD0Refl=kFALSE;
const TString rflFitType="DoubleGaus";
const TString fileNameMCD0refl="./reflections/reflections_fitted_DoubleGaus.root";
const Int_t types=kGaus;//kReflTempl;
const Int_t nMins=5;
const Double_t minMassForFit[nMins]={1.65,1.68,1.70,1.73,1.76};
const Int_t nMaxs=5;
const Double_t maxMassForFit[nMaxs]={2.15,2.10,2.05,2.00,1.95};
const Double_t nSigmaForCounting=3.5;
const Bool_t fixAlsoMass=kFALSE;
const Double_t maxchi=2;

//in-out efficiency
const Double_t effInOverEffOut=1.03;

//not to be set
enum CutMethod{kAbsCut,kPercCut,kPercCutVsCent,kPercentileCut}; //kAbsCut->absolute cut values, kPercCut->cut according to the % of events with smaller/larger q2, kPercCutVsCent->cut according to the % of events with smaller/larger q2 in finer centrality bins, kPercentileCut->cut on percentile (if enabled in the task)
enum SmallOrLarge{kSmall,kLarge,kIntegrated};
enum AnalysisMethod{kEventPlane,kEventPlaneInOut,kScalarProd};

Int_t colors[]={kRed+1,kBlack,kBlue+2,kGreen+2,kMagenta+3,kOrange+7,kCyan+2,kViolet+5,kYellow+2,kBlue-7,kGreen,kMagenta,kAzure,kRed+2};
Int_t markers[]={kFullCircle,kFullSquare,kFullDiamond,kFullTriangleUp,kFullTriangleDown,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenTriangleUp,kOpenTriangleDown};

//_____________________________________________________________________________________________
//FUNCTION PROTOTYPES
Int_t DmesonsFlowEvShapeYieldSyst(Int_t cutmeth=kPercentileCut, Int_t analysismeth=kEventPlaneInOut);
TList* LoadTList();
THnSparseF* LoadSparseFromList(TList* inputlist);
TH2F* GetHistoq2VsCentr(TList* inputlist);
TList* LoadMassHistos(THnSparseF* sparse, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge, Int_t analysismeth);
TList* LoadResolutionHistos(TList *inputlist, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge);
TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Int_t analysismeth, TGraphAsymmErrors *gRelSystEff);
void FillSignalGraph(TList *masslist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2,TGraphAsymmErrors **gSigmaFree,TGraphAsymmErrors **gSigmaFixed, TGraphAsymmErrors **gChiSquareFree, TGraphAsymmErrors **gChiSquareFixed, Int_t analysismeth, Int_t bkgfunc, Double_t minfit, Double_t maxfit, Int_t rebin);
Int_t LoadRefGraphs(TString reffilename, TH1F **hRawYieldRef, TH1F **hRawYieldfsRef, TH1F **hRawYieldBC1Ref, TH1F **hRawYieldBC2Ref, TGraphAsymmErrors *&gv2Ref, TGraphAsymmErrors *&gv2fsRef, TGraphAsymmErrors *&gv2BC1Ref, TGraphAsymmErrors *&gv2BC2Ref, Int_t smallorlarge, Int_t analysismeth);
Double_t GetEventPlaneResolution(Double_t& error, TH1F* hsubev1, TH1F* hsubev2, TH1F* hsubev3);
void ApplyCut(THnSparseF* sparse, Double_t min, Double_t max, UInt_t axnum);
void Applyq2Cut(THnSparseF* sparse, Double_t smallcutvalues, Double_t largecutvalues, Int_t smallorlarge);
Bool_t Defineq2Cuts(TH2F* hq2VsCentr, vector<Double_t> &smallcutvalues, vector<Double_t> &largecutvalues, Int_t cutmeth, Int_t fSparseVers);
Bool_t Getq2CutValuePercEvents(TH1F* hq2, Double_t cutvalues[2]);
void ResetAxes(THnSparseF* sparse, Int_t axnum=-1);
void GetMinMaxHisto(TH1F* histo, Double_t &xmin, Double_t &xmax);
void DivideCanvas(TCanvas* c, const Int_t nPtBins);
void SetStyle();
Bool_t LoadD0toKpiMCHistos(TList *outlist);

//_____________________________________________________________________________________________
//ANALYSIS FUNCTION
Int_t DmesonsFlowEvShapeYieldSyst(Int_t cutmeth, Int_t analysismeth) {
  
  TString workdir=gSystem->pwd();
  if(!gSystem->cd(outputdir.Data())) {gSystem->mkdir(outputdir.Data());}
  else {gSystem->cd(workdir.Data());}
  
  TList* datalist=(TList*)LoadTList();
  if(!datalist){return 1;}
  THnSparseF* hMassPtPhiq2Centr=(THnSparseF*)LoadSparseFromList(datalist);
  if(!hMassPtPhiq2Centr) {return 2;}
  Int_t fSparseVers = 0; //check version of the task output
  if(hMassPtPhiq2Centr->GetAxis(4)) {
    TString centtitle=hMassPtPhiq2Centr->GetAxis(4)->GetTitle();
    if(centtitle.Contains("Centrality") || centtitle.Contains("centrality")) {
      cout << "Version of the sparse with centrality axis!" << endl;
      fSparseVers=1;
    }
    else {
      cout << "Version of the sparse without centrality axis!" << endl;
    }
  }
  TH2F* hq2VsCentr=(TH2F*)GetHistoq2VsCentr(datalist);
  if(!hq2VsCentr) {return 3;}
  
  vector<Double_t> smallcutvalues;
  vector<Double_t> largecutvalues;
  Bool_t defq2cut=Defineq2Cuts(hq2VsCentr,smallcutvalues,largecutvalues,cutmeth,fSparseVers);
  if(!defq2cut) {return 4;}
  
  TString analysismethname="";
  if(analysismeth==kEventPlane) {analysismethname="EP";}
  else if(analysismeth==kEventPlaneInOut) {analysismethname="InOut";}
  else if(analysismeth==kScalarProd) {analysismethname="ScalarProduct";}
  
  Int_t q2region[3] = {kSmall,kLarge,kIntegrated};
  TString q2regionname[3] = {"q2Small","q2Large","q2Int"};
  
  SetStyle();
  
  //Load mass histograms corresponding to the required centrality, pt range and phi binning
  if(analysismeth!=kScalarProd) {
    printf("Load mass histos \n");
    TList *masslist[3];
    masslist[0]=(TList*)LoadMassHistos(hMassPtPhiq2Centr,smallcutvalues,largecutvalues,q2region[0],analysismeth);
    masslist[1]=(TList*)LoadMassHistos(hMassPtPhiq2Centr,smallcutvalues,largecutvalues,q2region[1],analysismeth);
    masslist[2]=(TList*)LoadMassHistos(hMassPtPhiq2Centr,smallcutvalues,largecutvalues,q2region[2],analysismeth);
    
    Int_t nPhi=nPhiBins;
    if(analysismeth==kEventPlaneInOut)nPhi=2;
    
    //Average pt
    Float_t averagePt[3][nPtBins];
    Float_t errorPt[3][nPtBins];
    
    for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
      //average pt for pt bin
      AliVertexingHFUtils *utils=new AliVertexingHFUtils();
      Int_t minCentTimesTen=minCent*10;
      Int_t maxCentTimesTen=maxCent*10;
      TH2F* hmasspt=0x0;
      //apply q2 cut
      for(Int_t iCent=0; iCent<ncentbins; iCent++) {
        if(iq2!=kIntegrated) {Applyq2Cut(hMassPtPhiq2Centr,smallcutvalues[iCent],largecutvalues[iCent],q2region[iq2]);}
        TH2F* htmp=(TH2F*)hMassPtPhiq2Centr->Projection(0,1);
        if(iCent==0) {hmasspt=(TH2F*)htmp->Clone();}
        else {hmasspt->Add(htmp);}
        ResetAxes(hMassPtPhiq2Centr,-1.);
      }
      for(Int_t iPt=0;iPt<nPtBins;iPt++){
        Int_t binMin=hmasspt->FindBin(PtLims[iPt]);
        Int_t binMax=hmasspt->FindBin(PtLims[iPt+1]-0.001);
        if(TMath::Abs(hmasspt->GetXaxis()->GetBinLowEdge(binMin)-PtLims[iPt])>0.001 ||
           TMath::Abs(hmasspt->GetXaxis()->GetBinUpEdge(binMax)-PtLims[iPt+1])>0.001){
          cout << "Error in pt bin limits for projection!\n" << endl;
          return 4;
        }
        TH1F *histtofit = (TH1F*)hmasspt->ProjectionY("_py",binMin,binMax);
        Int_t nMassBins=histtofit->GetNbinsX();
        Double_t hmin=histtofit->GetBinLowEdge(2); // need wide range for <pt>
        Double_t hmax=histtofit->GetBinLowEdge(nMassBins-2); // need wide range for <pt>
        if (partname.Contains("Dstar")) {
          if (hmin < 0.140) hmin=0.140;
          if (hmax > 0.165) hmax=0.165;
        }
        AliHFMassFitterVAR* fitter=new AliHFMassFitterVAR(histtofit,hmin,hmax,1,0,types);
        if(useTemplD0Refl){
          Printf("USE TEMPLATE FOR AVERAGE Pt");
          TH1F *hrflTempl=(TH1F*)(masslist[iq2]->FindObject(Form("histRfl_%d",iPt)))->Clone(Form("histrfl_%d",iPt));
          if(!hrflTempl) {Printf("histRfl_%d not found",iPt); return 15;}
          TH1F *hsigMC=(TH1F*)(masslist[iq2]->FindObject(Form("histSgn_%d",iPt)))->Clone(Form("histsgn_%d",iPt));
          if(!hsigMC) {Printf("histSgn_%d not found",iPt); return 15;}
          fitter->SetTemplateReflections(hrflTempl);
          Float_t sOverRef=(hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999)))/(hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999)));
          Printf("R OVER S = %f",sOverRef);
          fitter->SetFixReflOverS(sOverRef,kTRUE);
        }
        if(partname.Contains("Dstar")) {
          fitter->SetInitialGaussianMean(0.145);
          fitter->SetInitialGaussianSigma(0.0004);
        }
        fitter->SetUseLikelihoodFit();
        fitter->MassFitter(kFALSE);
        Double_t massFromFit=fitter->GetMean();
        Double_t sigmaFromFit=fitter->GetSigma();
        TF1* funcB2=fitter->GetBackgroundRecalcFunc();
        utils->AveragePt(averagePt[iq2][iPt],errorPt[iq2][iPt],PtLims[iPt],PtLims[iPt+1],hmasspt,massFromFit,sigmaFromFit,funcB2,2.5,4.5,0.,3.,1);
        if(averagePt[iq2][iPt]==0) averagePt[iq2][iPt] = (PtLims[iPt+1]+PtLims[iPt])/2;
      }
      cout << Form("Average pt %s region \n",q2regionname[iq2].Data()) << endl;
      for(Int_t iPt=0;iPt<nPtBins;iPt++) cout <<Form("%f +- %f\n",averagePt[iq2][iPt],errorPt[iq2][iPt])<<endl;
    }
    
    printf("Fill TGraphs for signal \n");
    //Fill TGraphs for signal
    TGraphAsymmErrors *gSignal[3][nPtBins];
    TGraphAsymmErrors *gSignalfs[3][nPtBins];
    TGraphAsymmErrors *gSignalBC1[3][nPtBins];
    TGraphAsymmErrors *gSignalBC2[3][nPtBins];
    TGraphAsymmErrors *gSigmaFree[3][nPtBins];
    TGraphAsymmErrors *gSigmaFixed[3][nPtBins];
    TGraphAsymmErrors *gChiSquareFree[3][nPtBins];
    TGraphAsymmErrors *gChiSquareFixed[3][nPtBins];
    for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
      for(Int_t iPt=0;iPt<nPtBins;iPt++){
        gSignal[iq2][iPt]=new TGraphAsymmErrors(nPhi);
        gSignalfs[iq2][iPt]=new TGraphAsymmErrors(nPhi);
        gSignalBC1[iq2][iPt]=new TGraphAsymmErrors(nPhi);
        gSignalBC2[iq2][iPt]=new TGraphAsymmErrors(nPhi);
        
        gSigmaFree[iq2][iPt]=new TGraphAsymmErrors(nPhi);
        gSigmaFixed[iq2][iPt]=new TGraphAsymmErrors(nPhi);
        
        gChiSquareFree[iq2][iPt]=new TGraphAsymmErrors(nPhi);
        gChiSquareFixed[iq2][iPt]=new TGraphAsymmErrors(nPhi);
      }
    }
    
    //EP resolution
    Double_t resol[3] = {-1.,-1.,-1.};
    Double_t errorres[3] = {-1.,-1.,-1.};
    
    if(useAliHandlerForRes) {
      cerr << "Error: AliHandler for resolution not yet implemented. Exit" << endl;
      return 6;
    }
    else {
      for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
        TList* resolist=(TList*)LoadResolutionHistos(datalist,smallcutvalues,largecutvalues,q2region[iq2]);
        TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
        TH1F* hevplresos[3];
        TString namereso[3]={"ResoVsq2","Reso2Vsq2","Reso3Vsq2"};
        Int_t nSubRes=1;
        TH2F* htestversion=(TH2F*)resolist->FindObject(Form("hEvPlane%s%s",namereso[0].Data(),suffixcentr.Data()));
        if(htestversion){
          printf("Old version of the task\n");
        }else{
          printf("New version of the task\n");
          namereso[0]="Reso1Vsq2";
        }
        if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
           evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
        for(Int_t iRes=0;iRes<nSubRes;iRes++){
          hevplresos[iRes]=(TH1F*)resolist->FindObject(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()));
        }
        resol[iq2]=GetEventPlaneResolution(errorres[iq2],hevplresos[0],hevplresos[1],hevplresos[2]);
        cout << "Event plane resolution "<<q2regionname[iq2]<<" region: " << resol[iq2] << " +/- " << errorres[iq2] << endl;
      }
    }
    
    cout << "Compute v2" << endl;
    //compute v2
    
    //load reference graphs
    TGraphAsymmErrors* gv2Ref[3];
    TGraphAsymmErrors* gv2fsRef[3];
    TGraphAsymmErrors* gv2BC1Ref[3];
    TGraphAsymmErrors* gv2BC2Ref[3];
    
    TH1F* hRawYieldRef[3][nPhi];
    TH1F* hRawYieldfsRef[3][nPhi];
    TH1F* hRawYieldBC1Ref[3][nPhi];
    TH1F* hRawYieldBC2Ref[3][nPhi];
    for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
      Int_t loadref=LoadRefGraphs(reffilename, hRawYieldRef[iq2], hRawYieldfsRef[iq2], hRawYieldBC1Ref[iq2], hRawYieldBC2Ref[iq2], gv2Ref[iq2], gv2fsRef[iq2], gv2BC1Ref[iq2], gv2BC2Ref[iq2], q2region[iq2], analysismeth);
      if(loadref>0) {return 7;}
    }
    
    TTree* treeMultiTrial = new TTree("treeMultiTrial","Multi-trial tree");
    //variables for the TTree branches
    const Int_t constnphi = nPhi;
    Double_t v2[3];
    Double_t v2staterr[3];
    Double_t resv2[3];
    Int_t ptbin;
    Double_t extractionmethod;
    Int_t iTrial;
    Double_t rawyield[3][constnphi];
    Double_t sigma[3][constnphi];
    Double_t chi[3][constnphi];
    
    treeMultiTrial->Branch("ptbin",&ptbin);
    treeMultiTrial->Branch("extrmethod",&extractionmethod);
    treeMultiTrial->Branch("TrialNum",&iTrial);
    
    for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
      treeMultiTrial->Branch(Form("v2_%s",q2regionname[iq2].Data()),&v2[iq2]);
      treeMultiTrial->Branch(Form("v2err_%s",q2regionname[iq2].Data()),&v2staterr[iq2]);
      treeMultiTrial->Branch(Form("deltav2_%s",q2regionname[iq2].Data()),&resv2[iq2]);
      for(Int_t iPhi=0; iPhi<constnphi; iPhi++) {
        treeMultiTrial->Branch(Form("rawyield_%d_%s",iPhi,q2regionname[iq2].Data()),&rawyield[iq2][iPhi]);
        treeMultiTrial->Branch(Form("sigma_%d_%s",iPhi,q2regionname[iq2].Data()),&sigma[iq2][iPhi]);
        treeMultiTrial->Branch(Form("chisquare_%d_%s",iPhi,q2regionname[iq2].Data()),&chi[iq2][iPhi]);
      }
    }
    
    TGraphAsymmErrors *gv2[3];
    TGraphAsymmErrors *gv2fs[3];
    TGraphAsymmErrors *gv2BC1[3];
    TGraphAsymmErrors *gv2BC2[3];
    
    iTrial=0;
    for(Int_t iBkgFcn=0; iBkgFcn<nBkgFcn; iBkgFcn++) {
      for(Int_t iReb=0; iReb<nReb; iReb++) {
        for(Int_t iMin=0; iMin<nMins; iMin++) {
          for(Int_t iMax=0; iMax<nMaxs; iMax++) {
            
            Bool_t ok=kTRUE;
            for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
              FillSignalGraph(masslist[iq2],gSignal[iq2],gSignalfs[iq2],gSignalBC1[iq2],gSignalBC2[iq2],gSigmaFree[iq2],gSigmaFixed[iq2],gChiSquareFree[iq2],gChiSquareFixed[iq2],analysismeth,typeb[iBkgFcn],minMassForFit[iMin],maxMassForFit[iMax],rebin[iReb]);
              
              gv2[iq2]=Computev2(gSignal[iq2],resol[iq2],averagePt[iq2],analysismeth,0x0);
              gv2fs[iq2]=Computev2(gSignalfs[iq2],resol[iq2],averagePt[iq2],analysismeth,0x0);
              gv2BC1[iq2]=Computev2(gSignalBC1[iq2],resol[iq2],averagePt[iq2],analysismeth,0x0);
              gv2BC2[iq2]=Computev2(gSignalBC2[iq2],resol[iq2],averagePt[iq2],analysismeth,0x0);
              if(!gv2[iq2] || !gv2fs[iq2] || !gv2BC1[iq2] || !gv2BC2[iq2]) {ok=kFALSE;}
            }
            if(!ok) {continue;}
            
            //fill tree for sigma free entries
            extractionmethod=0;
            for(Int_t iPt=0; iPt<nPtBins; iPt++) {
              ptbin=iPt;
              for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
                Double_t ptcent, v2ref[3];
                gv2[iq2]->GetPoint(iPt,ptcent,v2[iq2]);
                gv2Ref[iq2]->GetPoint(iPt,ptcent,v2ref[iq2]);
                v2staterr[iq2] = gv2[iq2]->GetErrorYhigh(iPt);
                resv2[iq2] = v2[iq2] - v2ref[iq2];
                for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
                  Double_t phi;
                  gSignal[iq2][iPt]->GetPoint(iPhi,phi,rawyield[iq2][iPhi]);
                  gSigmaFree[iq2][iPt]->GetPoint(iPhi,phi,sigma[iq2][iPhi]);
                  gChiSquareFree[iq2][iPt]->GetPoint(iPhi,phi,chi[iq2][iPhi]);
                }
                if(iq2==kIntegrated) {treeMultiTrial->Fill();}
              }
            }
            
            //fill tree for sigma fixed entries
            extractionmethod=1;
            for(Int_t iPt=0; iPt<nPtBins; iPt++) {
              ptbin=iPt;
              for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
                Double_t ptcent, v2ref[3];
                gv2fs[iq2]->GetPoint(iPt,ptcent,v2[iq2]);
                gv2fsRef[iq2]->GetPoint(iPt,ptcent,v2ref[iq2]);
                v2staterr[iq2] = gv2fs[iq2]->GetErrorYhigh(iPt);
                resv2[iq2] = v2[iq2] - v2ref[iq2];
                for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
                  Double_t phi;
                  gSignalfs[iq2][iPt]->GetPoint(iPhi,phi,rawyield[iq2][iPhi]);
                  gSigmaFixed[iq2][iPt]->GetPoint(iPhi,phi,sigma[iq2][iPhi]);
                  gChiSquareFixed[iq2][iPt]->GetPoint(iPhi,phi,chi[iq2][iPhi]);
                }
                if(iq2==kIntegrated) {treeMultiTrial->Fill();}
              }
            }
            
            //fill tree for BC1 entries
            extractionmethod=2;
            for(Int_t iPt=0; iPt<nPtBins; iPt++) {
              ptbin=iPt;
              for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
                Double_t ptcent, v2ref[3];
                gv2BC1[iq2]->GetPoint(iPt,ptcent,v2[iq2]);
                gv2BC1Ref[iq2]->GetPoint(iPt,ptcent,v2ref[iq2]);
                v2staterr[iq2] = gv2BC1[iq2]->GetErrorYhigh(iPt);
                resv2[iq2] = v2[iq2] - v2ref[iq2];
                for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
                  Double_t phi;
                  gSignalBC1[iq2][iPt]->GetPoint(iPhi,phi,rawyield[iq2][iPhi]);
                  gSigmaFree[iq2][iPt]->GetPoint(iPhi,phi,sigma[iq2][iPhi]);
                  gChiSquareFree[iq2][iPt]->GetPoint(iPhi,phi,chi[iq2][iPhi]);
                }
                if(iq2==kIntegrated) {treeMultiTrial->Fill();}
              }
            }
            
            //fill tree for BC2 entries
            extractionmethod=3;
            for(Int_t iPt=0; iPt<nPtBins; iPt++) {
              ptbin=iPt;
              for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
                Double_t ptcent, v2ref[3];
                gv2BC2[iq2]->GetPoint(iPt,ptcent,v2[iq2]);
                gv2BC2Ref[iq2]->GetPoint(iPt,ptcent,v2ref[iq2]);
                v2staterr[iq2] = gv2BC2[iq2]->GetErrorYhigh(iPt);
                resv2[iq2] = v2[iq2] - v2ref[iq2];
                for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
                  Double_t phi;
                  gSignalBC2[iq2][iPt]->GetPoint(iPhi,phi,rawyield[iq2][iPhi]);
                  gSigmaFree[iq2][iPt]->GetPoint(iPhi,phi,sigma[iq2][iPhi]);
                  gChiSquareFree[iq2][iPt]->GetPoint(iPhi,phi,chi[iq2][iPhi]);
                }
                if(iq2==kIntegrated) {treeMultiTrial->Fill();}
              }
            }
            
            for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
              delete gv2[iq2];
              delete gv2fs[iq2];
              delete gv2BC1[iq2];
              delete gv2BC2[iq2];
              gv2[iq2]=0x0;
              gv2fs[iq2]=0x0;
              gv2BC1[iq2]=0x0;
              gv2BC2[iq2]=0x0;
            }
            iTrial++;
          }
        }
      }
    }
    Int_t nTrials=iTrial+1;
    
    //Histos from tree
    TH1F* hResv2[3][nPtBins];
    TH1F* hResv2fs[3][nPtBins];
    TH1F* hResv2BC1[3][nPtBins];
    TH1F* hResv2BC2[3][nPtBins];
    TPaveStats *pv2[3][nPtBins];
    TPaveStats *pv2fs[3][nPtBins];
    TPaveStats *pv2BC1[3][nPtBins];
    TPaveStats *pv2BC2[3][nPtBins];
    
    TH1F* hv2VsTrial[3][nPtBins];
    TH1F* hv2fsVsTrial[3][nPtBins];
    TH1F* hv2BC1VsTrial[3][nPtBins];
    TH1F* hv2BC2VsTrial[3][nPtBins];
    
    TH1F* hv2fsRatioLargeSmallVsTrial[nPtBins];
    TH1F* hv2fsRatioLargeUnbiasedVsTrial[nPtBins];
    TH1F* hv2fsRatioSmallUnbiasedVsTrial[nPtBins];
    
    TH1F* hv2fsRatioLargeSmall[nPtBins];
    TH1F* hv2fsRatioLargeUnbiased[nPtBins];
    TH1F* hv2fsRatioSmallUnbiased[nPtBins];
    
    TH2F* hResv2CorrLargeSmall[nPtBins];
    TH2F* hResv2CorrLargeUnbiased[nPtBins];
    TH2F* hResv2CorrSmallUnbiased[nPtBins];
    
    TString selection;
    TString selectionfs;
    TString selectionBC1;
    TString selectionBC2;
    TString selectioncorr;
    
    Int_t nbins=50;
    
    TCanvas* cResv2CorrLargeSmall = new TCanvas("cResv2CorrLargeSmall","",1920,1080);
    DivideCanvas(cResv2CorrLargeSmall,nPtBins);
    
    TCanvas* cResv2CorrLargeUnbiased = new TCanvas("cResv2CorrLargeUnbiased","",1920,1080);
    DivideCanvas(cResv2CorrLargeUnbiased,nPtBins);
    
    TCanvas* cResv2CorrSmallUnbiased = new TCanvas("cResv2CorrSmallUnbiased","",1920,1080);
    DivideCanvas(cResv2CorrSmallUnbiased,nPtBins);
    
    TCanvas* cv2fsRatioLargeSmallVsTrial = new TCanvas("cv2fsRatioLargeSmallVsTrial","",1920,1080);
    DivideCanvas(cv2fsRatioLargeSmallVsTrial,nPtBins);
    
    TCanvas* cv2fsRatioLargeUnbiasedVsTrial = new TCanvas("cv2fsRatioLargeUnbiasedVsTrial","",1920,1080);
    DivideCanvas(cv2fsRatioLargeUnbiasedVsTrial,nPtBins);
    
    TCanvas* cv2fsRatioSmallUnbiasedVsTrial = new TCanvas("cv2fsRatioSmallUnbiasedVsTrial","",1920,1080);
    DivideCanvas(cv2fsRatioSmallUnbiasedVsTrial,nPtBins);
    
    TCanvas* cv2fsRatioLargeSmall = new TCanvas("cv2fsRatioLargeSmall","",1920,1080);
    DivideCanvas(cv2fsRatioLargeSmall,nPtBins);
    
    TCanvas* cv2fsRatioLargeUnbiased = new TCanvas("cv2fsRatioLargeUnbiased","",1920,1080);
    DivideCanvas(cv2fsRatioLargeUnbiased,nPtBins);
    
    TCanvas* cv2fsRatioSmallUnbiased = new TCanvas("cv2fsRatioSmallUnbiased","",1920,1080);
    DivideCanvas(cv2fsRatioSmallUnbiased,nPtBins);
    
    TCanvas *cv2[3];
    TCanvas *cv2VsTrial[3];
    TCanvas* cv2Syst[3];
    for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
      cv2[iq2] = new TCanvas(Form("cv2_%s",q2regionname[iq2].Data()),q2regionname[iq2].Data(),1920,1080);
      DivideCanvas(cv2[iq2],nPtBins);
      cv2VsTrial[iq2] = new TCanvas(Form("cv2VsTrial_%s",q2regionname[iq2].Data()),q2regionname[iq2].Data(),1920,1080);
      DivideCanvas(cv2VsTrial[iq2],nPtBins);
      cv2Syst[iq2] = new TCanvas(Form("cv2Syst_%s",q2regionname[iq2].Data()),q2regionname[iq2].Data(),1920,1080);
      DivideCanvas(cv2Syst[iq2],nPtBins);
    }
    
    TCanvas* cShiftRMS = new TCanvas("cShiftRMS","",1920,1080);
    
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
        
        Double_t v2refstaterr = gv2Ref[iq2]->GetErrorYhigh(iPt);
        hResv2[iq2][iPt] = new TH1F(Form("hResv2_%s_pt%d",q2regionname[iq2].Data(),iPt),"",nbins,-3*v2refstaterr,3*v2refstaterr);
        hResv2fs[iq2][iPt] = new TH1F(Form("hResv2fs_%s_pt%d",q2regionname[iq2].Data(),iPt),"",nbins,-3*v2refstaterr,3*v2refstaterr);
        hResv2BC1[iq2][iPt] = new TH1F(Form("hResv2BC1_%s_pt%d",q2regionname[iq2].Data(),iPt),"",nbins,-3*v2refstaterr,3*v2refstaterr);
        hResv2BC2[iq2][iPt] = new TH1F(Form("hResv2BC2_%s_pt%d",q2regionname[iq2].Data(),iPt),"",nbins,-3*v2refstaterr,3*v2refstaterr);
        hResv2[iq2][iPt]->GetXaxis()->SetTitle("v_{2}-v_{2}^{ref} {EP}");
        hResv2[iq2][iPt]->GetYaxis()->SetTitle("Entries");
        hResv2[iq2][iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
        hResv2fs[iq2][iPt]->GetXaxis()->SetTitle("v_{2}-v_{2}^{ref} {EP}");
        hResv2fs[iq2][iPt]->GetYaxis()->SetTitle("Entries");
        hResv2fs[iq2][iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
        hResv2BC1[iq2][iPt]->GetXaxis()->SetTitle("v_{2}-v_{2}^{ref} {EP}");
        hResv2BC1[iq2][iPt]->GetYaxis()->SetTitle("Entries");
        hResv2BC1[iq2][iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
        hResv2BC2[iq2][iPt]->GetXaxis()->SetTitle("v_{2}-v_{2}^{ref} {EP}");
        hResv2BC2[iq2][iPt]->GetYaxis()->SetTitle("Entries");
        hResv2BC2[iq2][iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
        hResv2[iq2][iPt]->SetLineColor(colors[0]);
        hResv2fs[iq2][iPt]->SetLineColor(colors[1]);
        hResv2BC1[iq2][iPt]->SetLineColor(colors[2]);
        hResv2BC2[iq2][iPt]->SetLineColor(colors[3]);
        
        hv2VsTrial[iq2][iPt] = new TH1F(Form("hv2VsTrial_%s_pt%d",q2regionname[iq2].Data(),iPt),"",nTrials,0.5,nTrials+0.5);
        hv2fsVsTrial[iq2][iPt] = new TH1F(Form("hv2fsVsTrial_%s_pt%d",q2regionname[iq2].Data(),iPt),"",nTrials,0.5,nTrials+0.5);
        hv2BC1VsTrial[iq2][iPt] = new TH1F(Form("hv2BC1VsTrial_%s_pt%d",q2regionname[iq2].Data(),iPt),"",nTrials,0.5,nTrials+0.5);
        hv2BC2VsTrial[iq2][iPt] = new TH1F(Form("hv2BC2VsTrial_%s_pt%d",q2regionname[iq2].Data(),iPt),"",nTrials,0.5,nTrials+0.5);
        hv2VsTrial[iq2][iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
        hv2VsTrial[iq2][iPt]->GetXaxis()->SetTitle("Trial #");
        hv2VsTrial[iq2][iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
        hv2fsVsTrial[iq2][iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
        hv2fsVsTrial[iq2][iPt]->GetXaxis()->SetTitle("Trial #");
        hv2fsVsTrial[iq2][iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
        hv2BC1VsTrial[iq2][iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
        hv2BC1VsTrial[iq2][iPt]->GetXaxis()->SetTitle("Trial");
        hv2BC1VsTrial[iq2][iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
        hv2BC2VsTrial[iq2][iPt]->GetYaxis()->SetTitle("v_{2} {EP}");
        hv2BC2VsTrial[iq2][iPt]->GetXaxis()->SetTitle("Trial #");
        hv2BC2VsTrial[iq2][iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
        hv2VsTrial[iq2][iPt]->SetLineColor(colors[0]);
        hv2fsVsTrial[iq2][iPt]->SetLineColor(colors[1]);
        hv2BC1VsTrial[iq2][iPt]->SetLineColor(colors[2]);
        hv2BC2VsTrial[iq2][iPt]->SetLineColor(colors[3]);
        hv2VsTrial[iq2][iPt]->SetLineWidth(1);
        hv2fsVsTrial[iq2][iPt]->SetLineWidth(1);
        hv2BC1VsTrial[iq2][iPt]->SetLineWidth(1);
        hv2BC2VsTrial[iq2][iPt]->SetLineWidth(1);
        
        selection = Form("ptbin>%f && ptbin<%f && extrmethod>-0.5 && extrmethod<0.5",iPt-0.01,iPt+0.01);
        for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
          selection+=Form(" && chisquare_%d_%s<%f",iPhi,q2regionname[iq2].Data(),maxchi);
        }
        selectionfs=selection;
        selectionBC1=selection;
        selectionBC2=selection;
        selectionfs.ReplaceAll("extrmethod>-0.5 && extrmethod<0.5","extrmethod>0.5 && extrmethod<1.5");
        selectionBC1.ReplaceAll("extrmethod>-0.5 && extrmethod<0.5","extrmethod>1.5 && extrmethod<2.5");
        selectionBC2.ReplaceAll("extrmethod>-0.5 && extrmethod<0.5","extrmethod>2.5 && extrmethod<3.5");
        treeMultiTrial->Project(Form("hResv2_%s_pt%d",q2regionname[iq2].Data(),iPt),Form("deltav2_%s",q2regionname[iq2].Data()),selection.Data());
        treeMultiTrial->Project(Form("hResv2fs_%s_pt%d",q2regionname[iq2].Data(),iPt),Form("deltav2_%s",q2regionname[iq2].Data()),selectionfs.Data());
        treeMultiTrial->Project(Form("hResv2BC1_%s_pt%d",q2regionname[iq2].Data(),iPt),Form("deltav2_%s",q2regionname[iq2].Data()),selectionBC1.Data());
        treeMultiTrial->Project(Form("hResv2BC2_%s_pt%d",q2regionname[iq2].Data(),iPt),Form("deltav2_%s",q2regionname[iq2].Data()),selectionBC2.Data());
      }
      
      hResv2CorrLargeSmall[iPt] = new TH2F(Form("hResv2CorrLargeSmall_pt%d",iPt),"",nbins,-hResv2[0][iPt]->GetRMS()*5,hResv2[0][iPt]->GetRMS()*5,nbins,-hResv2[0][iPt]->GetRMS()*5,hResv2[0][iPt]->GetRMS()*5);
      hResv2CorrLargeSmall[iPt]->SetStats(kFALSE);
      hResv2CorrLargeSmall[iPt]->GetYaxis()->SetTitle("v_{2} - v_{2}^{ref} (Large-q_{2}) (sigma fix)");
      hResv2CorrLargeSmall[iPt]->GetXaxis()->SetTitle("v_{2} - v_{2}^{ref} (Small-q_{2}) (sigma fix)");
      hResv2CorrLargeSmall[iPt]->GetYaxis()->SetTitleOffset(1.4);
      hResv2CorrLargeSmall[iPt]->GetXaxis()->SetNdivisions(505);
      hResv2CorrLargeSmall[iPt]->GetYaxis()->SetNdivisions(505);
      hResv2CorrLargeSmall[iPt]->GetXaxis()->SetTitleOffset(1.3);
      
      hResv2CorrLargeUnbiased[iPt] = new TH2F(Form("hResv2CorrLargeUnbiased_pt%d",iPt),"",nbins,-hResv2[1][iPt]->GetRMS()*5,hResv2[1][iPt]->GetRMS()*5,nbins,-hResv2[1][iPt]->GetRMS()*5,hResv2[1][iPt]->GetRMS()*5);
      hResv2CorrLargeUnbiased[iPt]->SetStats(kFALSE);
      hResv2CorrLargeUnbiased[iPt]->GetYaxis()->SetTitle("v_{2} - v_{2}^{ref} (Large-q_{2}) (sigma fix)");
      hResv2CorrLargeUnbiased[iPt]->GetXaxis()->SetTitle("v_{2} - v_{2}^{ref} (q_{2}-integrated) (sigma fix)");
      hResv2CorrLargeUnbiased[iPt]->GetYaxis()->SetTitleOffset(1.4);
      hResv2CorrLargeUnbiased[iPt]->GetXaxis()->SetNdivisions(505);
      hResv2CorrLargeUnbiased[iPt]->GetYaxis()->SetNdivisions(505);
      hResv2CorrLargeUnbiased[iPt]->GetXaxis()->SetTitleOffset(1.3);
      
      hResv2CorrSmallUnbiased[iPt] = new TH2F(Form("hResv2CorrSmallUnbiased_pt%d",iPt),"",nbins,-hResv2[0][iPt]->GetRMS()*5,hResv2[0][iPt]->GetRMS()*5,nbins,-hResv2[0][iPt]->GetRMS()*5,hResv2[0][iPt]->GetRMS()*5);
      hResv2CorrSmallUnbiased[iPt]->SetStats(kFALSE);
      hResv2CorrSmallUnbiased[iPt]->GetYaxis()->SetTitle("v_{2} - v_{2}^{ref} (Small-q_{2}) (sigma fix)");
      hResv2CorrSmallUnbiased[iPt]->GetXaxis()->SetTitle("v_{2} - v_{2}^{ref} (q_{2}-integrated) (sigma fix)");
      hResv2CorrSmallUnbiased[iPt]->GetYaxis()->SetTitleOffset(1.4);
      hResv2CorrSmallUnbiased[iPt]->GetXaxis()->SetNdivisions(505);
      hResv2CorrSmallUnbiased[iPt]->GetYaxis()->SetNdivisions(505);
      hResv2CorrSmallUnbiased[iPt]->GetXaxis()->SetTitleOffset(1.3);
      
      selectioncorr = Form("ptbin>%f && ptbin<%f && extrmethod>0.5 && extrmethod<1.5",iPt-0.01,iPt+0.01);
      for(Int_t iq2=kSmall; iq2<=kLarge; iq2++) {
        for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
          selectioncorr+=Form(" && chisquare_%d_%s<%f",iPhi,q2regionname[iq2].Data(),maxchi);
        }
      }
      treeMultiTrial->Project(Form("hResv2CorrLargeSmall_pt%d",iPt),"deltav2_q2Large:deltav2_q2Small",selectioncorr.Data());
      treeMultiTrial->Project(Form("hResv2CorrLargeUnbiased_pt%d",iPt),"deltav2_q2Large:deltav2_q2Int",selectioncorr.Data());
      treeMultiTrial->Project(Form("hResv2CorrSmallUnbiased_pt%d",iPt),"deltav2_q2Small:deltav2_q2Int",selectioncorr.Data());
    }
    
    Int_t itrial;
    Int_t extrmethod;
    
    treeMultiTrial->SetBranchAddress("TrialNum",&itrial);
    treeMultiTrial->SetBranchAddress("ptbin",&ptbin);
    treeMultiTrial->SetBranchAddress("extrmethod",&extractionmethod);
    for(Int_t iq2=kSmall; iq2<=kLarge; iq2++) {
      treeMultiTrial->SetBranchAddress(Form("v2_%s",q2regionname[iq2].Data()),&v2[iq2]);
      treeMultiTrial->SetBranchAddress(Form("v2err_%s",q2regionname[iq2].Data()),&v2staterr[iq2]);
      for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
        treeMultiTrial->SetBranchAddress(Form("chisquare_%d_%s",iPhi,q2regionname[iq2].Data()),&chi[iq2][iPhi]);
      }
    }
    
    for(Int_t iEntry=0; iEntry<treeMultiTrial->GetEntries(); iEntry++) {
      treeMultiTrial->GetEntry(iEntry);
      
      Bool_t isChiOk[3]={kTRUE,kTRUE,kTRUE};
      for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
        for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
          if(chi[iq2][iPhi]>maxchi) {
            isChiOk[iq2] = kFALSE;
            break;
          }
        }
        if(isChiOk[iq2]) {
          if(extractionmethod==0) {
            hv2VsTrial[iq2][ptbin]->SetBinContent(itrial+1,v2[iq2]);
            hv2VsTrial[iq2][ptbin]->SetBinError(itrial+1,v2staterr[iq2]);
          }
          else if(extractionmethod==1) {
            hv2fsVsTrial[iq2][ptbin]->SetBinContent(itrial+1,v2[iq2]);
            hv2fsVsTrial[iq2][ptbin]->SetBinError(itrial+1,v2staterr[iq2]);
          }
          else if(extractionmethod==2) {
            hv2BC1VsTrial[iq2][ptbin]->SetBinContent(itrial+1,v2[iq2]);
            hv2BC1VsTrial[iq2][ptbin]->SetBinError(itrial+1,v2staterr[iq2]);
          }
          else if(extractionmethod==3) {
            hv2BC2VsTrial[iq2][ptbin]->SetBinContent(itrial+1,v2[iq2]);
            hv2BC2VsTrial[iq2][ptbin]->SetBinError(itrial+1,v2staterr[iq2]);
          }
        }
      }
    }
    
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      hv2fsVsTrial[0][iPt]->Sumw2();
      hv2fsVsTrial[1][iPt]->Sumw2();
      hv2fsVsTrial[2][iPt]->Sumw2();
      hv2fsRatioLargeSmallVsTrial[iPt] = (TH1F*)hv2fsVsTrial[1][iPt]->Clone(Form("v2fsRatioLargeSmallVsTrialpt%d",iPt));
      hv2fsRatioLargeSmallVsTrial[iPt]->Divide(hv2fsVsTrial[1][iPt],hv2fsVsTrial[0][iPt],1.,1.);
      hv2fsRatioLargeSmallVsTrial[iPt]->GetYaxis()->SetTitle("v_{2} (Large-q_{2}) / v_{2} (Small-q_{2}) (sigma fix)");
      hv2fsRatioLargeSmallVsTrial[iPt]->SetStats(0);
      hv2fsRatioLargeSmallVsTrial[iPt]->SetLineWidth(1);
      hv2fsRatioLargeSmallVsTrial[iPt]->SetMarkerColor(kRed);
      hv2fsRatioLargeSmallVsTrial[iPt]->SetMarkerSize(0.5);
      hv2fsRatioLargeSmallVsTrial[iPt]->SetMarkerStyle(20);
      hv2fsRatioLargeSmallVsTrial[iPt]->GetYaxis()->SetRangeUser(hv2fsRatioLargeSmallVsTrial[iPt]->GetMinimum()-3*TMath::Abs(hv2fsRatioLargeSmallVsTrial[iPt]->GetBinError(1)),hv2fsRatioLargeSmallVsTrial[iPt]->GetMaximum()+3*TMath::Abs(hv2fsRatioLargeSmallVsTrial[iPt]->GetBinError(1)));
      
      hv2fsRatioLargeUnbiasedVsTrial[iPt] = (TH1F*)hv2fsVsTrial[1][iPt]->Clone(Form("v2fsRatioLargeUnbiasedVsTrial_pt%d",iPt));
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->Divide(hv2fsVsTrial[1][iPt],hv2fsVsTrial[2][iPt],1.,1.);
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetYaxis()->SetTitle("v_{2} (Large-q_{2}) / v_{2} (q_{2}-Integrated) (sigma fix)");
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->SetStats(0);
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->SetLineWidth(1);
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->SetMarkerColor(kRed);
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->SetMarkerSize(0.5);
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->SetMarkerStyle(20);
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetYaxis()->SetRangeUser(hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetMinimum()-3*TMath::Abs(hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetBinError(1)),hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetMaximum()+3*TMath::Abs(hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetBinError(1)));
      
      hv2fsRatioSmallUnbiasedVsTrial[iPt] = (TH1F*)hv2fsVsTrial[0][iPt]->Clone(Form("v2fsRatioSmallUnbiasedVsTrial_pt%d",iPt));
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->Divide(hv2fsVsTrial[0][iPt],hv2fsVsTrial[2][iPt],1.,1.);
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetYaxis()->SetTitle("v_{2} (Small-q_{2}) / v_{2} (q_{2}-Integrated) (sigma fix)");
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->SetStats(0);
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->SetLineWidth(1);
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->SetMarkerColor(kRed);
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->SetMarkerSize(0.5);
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->SetMarkerStyle(20);
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetYaxis()->SetRangeUser(hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetMinimum()-3*TMath::Abs(hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetBinError(1)),hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetMaximum()+3*TMath::Abs(hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetBinError(1)));
    }
    
    TH1F* hv2Syst[3][nPtBins];
    TH1F* hv2Syst2[3][nPtBins];
    TH1F* hv2fsSyst[3][nPtBins];
    TH1F* hv2fsSyst2[3][nPtBins];
    TH1F* hv2BC1Syst[3][nPtBins];
    TH1F* hv2BC1Syst2[3][nPtBins];
    TH1F* hv2BC2Syst[3][nPtBins];
    TH1F* hv2BC2Syst2[3][nPtBins];
    
    TH1F* hShiftRMS[3];
    Int_t q2regioncolor[3] = {kRed+1,kBlue+1,kGray+1};
    for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
      hShiftRMS[iq2] = new TH1F(Form("hShiftRMS_%s",q2regionname[iq2].Data()),"",nPtBins,PtLims);
      hShiftRMS[iq2]->SetLineColor(q2regioncolor[iq2]);
      hShiftRMS[iq2]->GetYaxis()->SetTitle("#sqrt{RMS^{2}+shift^{2}} (sigma fix)");
      hShiftRMS[iq2]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
      hShiftRMS[iq2]->SetStats(kFALSE);
      hShiftRMS[iq2]->GetYaxis()->SetRangeUser(0.,0.05);
    }
    
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
        hv2Syst[iq2][iPt] = new TH1F(Form("hv2Syst_%s_pt%d",q2regionname[iq2].Data(),iPt),"",4,-0.5,3.5);
        hv2Syst[iq2][iPt]->SetStats(kFALSE);
        hv2Syst[iq2][iPt]->GetXaxis()->SetBinLabel(1,"free sigma");
        hv2Syst[iq2][iPt]->GetXaxis()->SetBinLabel(2,"fixed sigma");
        hv2Syst[iq2][iPt]->GetXaxis()->SetBinLabel(3,"BC1");
        hv2Syst[iq2][iPt]->GetXaxis()->SetBinLabel(4,"BC2");
        hv2Syst[iq2][iPt]->GetYaxis()->SetTitle("MEAN v_{2}-v_{2}^{ref}");
        hv2Syst[iq2][iPt]->GetYaxis()->SetTitleOffset(1.5);
        hv2Syst[iq2][iPt]->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
        
        hv2Syst2[iq2][iPt] = new TH1F(Form("hv2Syst2_%s_pt%d",q2regionname[iq2].Data(),iPt),"",4,-0.4,3.6);
        hv2Syst2[iq2][iPt]->SetStats(kFALSE);
        
        hv2fsSyst[iq2][iPt] = (TH1F*)hv2Syst[iq2][iPt]->Clone();
        hv2fsSyst2[iq2][iPt] = (TH1F*)hv2Syst2[iq2][iPt]->Clone();
        hv2BC1Syst[iq2][iPt] = (TH1F*)hv2Syst[iq2][iPt]->Clone();
        hv2BC1Syst2[iq2][iPt] = (TH1F*)hv2Syst2[iq2][iPt]->Clone();
        hv2BC2Syst[iq2][iPt] = (TH1F*)hv2Syst[iq2][iPt]->Clone();
        hv2BC2Syst2[iq2][iPt] = (TH1F*)hv2Syst2[iq2][iPt]->Clone();
        hv2Syst2[iq2][iPt]->SetName(Form("hv2Syst2_%s_pt%d",q2regionname[iq2].Data(),iPt));
        hv2fsSyst[iq2][iPt]->SetName(Form("hv2fsSyst_%s_pt%d",q2regionname[iq2].Data(),iPt));
        hv2fsSyst2[iq2][iPt]->SetName(Form("hv2fsSyst2_%s_pt%d",q2regionname[iq2].Data(),iPt));
        hv2BC1Syst[iq2][iPt]->SetName(Form("hv2BC1Syst_%s_pt%d",q2regionname[iq2].Data(),iPt));
        hv2BC1Syst2[iq2][iPt]->SetName(Form("hv2BC1Syst2_%s_pt%d",q2regionname[iq2].Data(),iPt));
        hv2BC2Syst[iq2][iPt]->SetName(Form("hv2BC2Syst_%s_pt%d",q2regionname[iq2].Data(),iPt));
        hv2BC2Syst2[iq2][iPt]->SetName(Form("hv2BC2Syst2_%s_pt%d",q2regionname[iq2].Data(),iPt));
        
        hv2Syst[iq2][iPt]->SetLineColor(colors[0]);
        hv2fsSyst[iq2][iPt]->SetLineColor(colors[1]);
        hv2BC1Syst[iq2][iPt]->SetLineColor(colors[2]);
        hv2BC2Syst[iq2][iPt]->SetLineColor(colors[3]);
        hv2Syst[iq2][iPt]->SetMarkerColor(colors[0]);
        hv2fsSyst[iq2][iPt]->SetMarkerColor(colors[1]);
        hv2BC1Syst[iq2][iPt]->SetMarkerColor(colors[2]);
        hv2BC2Syst[iq2][iPt]->SetMarkerColor(colors[3]);
        hv2Syst[iq2][iPt]->SetMarkerStyle(markers[0]);
        hv2fsSyst[iq2][iPt]->SetMarkerStyle(markers[1]);
        hv2BC1Syst[iq2][iPt]->SetMarkerStyle(markers[2]);
        hv2BC2Syst[iq2][iPt]->SetMarkerStyle(markers[3]);
        hv2Syst2[iq2][iPt]->SetLineColor(colors[0]);
        hv2fsSyst2[iq2][iPt]->SetLineColor(colors[1]);
        hv2BC1Syst2[iq2][iPt]->SetLineColor(colors[2]);
        hv2BC2Syst2[iq2][iPt]->SetLineColor(colors[3]);
        hv2Syst2[iq2][iPt]->SetMarkerColor(colors[0]);
        hv2fsSyst2[iq2][iPt]->SetMarkerColor(colors[1]);
        hv2BC1Syst2[iq2][iPt]->SetMarkerColor(colors[2]);
        hv2BC2Syst2[iq2][iPt]->SetMarkerColor(colors[3]);
        hv2Syst2[iq2][iPt]->SetMarkerStyle(markers[5]);
        hv2fsSyst2[iq2][iPt]->SetMarkerStyle(markers[6]);
        hv2BC1Syst2[iq2][iPt]->SetMarkerStyle(markers[7]);
        hv2BC2Syst2[iq2][iPt]->SetMarkerStyle(markers[8]);
        
        hv2Syst[iq2][iPt]->SetBinContent(1,hResv2[iq2][iPt]->GetMean());
        hv2Syst[iq2][iPt]->SetBinError(1,hResv2[iq2][iPt]->GetRMS());
        hv2fsSyst[iq2][iPt]->SetBinContent(2,hResv2fs[iq2][iPt]->GetMean());
        hv2fsSyst[iq2][iPt]->SetBinError(2,hResv2fs[iq2][iPt]->GetRMS());
        hv2BC1Syst[iq2][iPt]->SetBinContent(3,hResv2BC1[iq2][iPt]->GetMean());
        hv2BC1Syst[iq2][iPt]->SetBinError(3,hResv2BC1[iq2][iPt]->GetRMS());
        hv2BC2Syst[iq2][iPt]->SetBinContent(4,hResv2BC2[iq2][iPt]->GetMean());
        hv2BC2Syst[iq2][iPt]->SetBinError(4,hResv2BC2[iq2][iPt]->GetRMS());
        hv2Syst2[iq2][iPt]->SetBinContent(1,hResv2[iq2][iPt]->GetMean());
        Double_t xmin, xmax;
        GetMinMaxHisto(hResv2[iq2][iPt],xmin,xmax);
        hv2Syst2[iq2][iPt]->SetBinError(1,(xmax-xmin)/TMath::Sqrt(12));
        hv2fsSyst2[iq2][iPt]->SetBinContent(2,hResv2fs[iq2][iPt]->GetMean());
        GetMinMaxHisto(hResv2fs[iq2][iPt],xmin,xmax);
        hv2fsSyst2[iq2][iPt]->SetBinError(2,(xmax-xmin)/TMath::Sqrt(12));
        hv2BC1Syst2[iq2][iPt]->SetBinContent(3,hResv2BC1[iq2][iPt]->GetMean());
        GetMinMaxHisto(hResv2BC1[iq2][iPt],xmin,xmax);
        hv2BC1Syst2[iq2][iPt]->SetBinError(3,(xmax-xmin)/TMath::Sqrt(12));
        hv2BC2Syst2[iq2][iPt]->SetBinContent(4,hResv2BC2[iq2][iPt]->GetMean());
        GetMinMaxHisto(hResv2BC2[iq2][iPt],xmin,xmax);
        hv2BC2Syst2[iq2][iPt]->SetBinError(4,(xmax-xmin)/TMath::Sqrt(12));
        hShiftRMS[iq2]->SetBinContent(iPt+1,TMath::Sqrt(hResv2fs[iq2][iPt]->GetMean()*hResv2fs[iq2][iPt]->GetMean()+hResv2fs[iq2][iPt]->GetRMS()*hResv2fs[iq2][iPt]->GetRMS()));
      }
    }
    
    TLatex* latex = new TLatex();
    latex->SetTextSize(0.05);
    latex->SetTextFont(42);
    
    TLine* v2refline[3][nPtBins];
    TLine* v2reflinelow[3][nPtBins];
    TLine* v2reflinehigh[3][nPtBins];
    TLegend* leg[3][nPtBins];
    TLegend* methlegend[3];
    
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
        cv2[iq2]->cd(iPt+1);
        hResv2[iq2][iPt]->GetYaxis()->SetRangeUser(0.,hResv2fs[iq2][iPt]->GetMaximum()*2.);
        hResv2[iq2][iPt]->Draw();
        hResv2fs[iq2][iPt]->Draw("sames");
        hResv2BC1[iq2][iPt]->Draw("sames");
        hResv2BC2[iq2][iPt]->Draw("sames");
        
        cv2[iq2]->cd(iPt+1)->Update();
        pv2[iq2][iPt] = (TPaveStats*)hResv2[iq2][iPt]->FindObject("stats");
        pv2[iq2][iPt]->SetTextColor(colors[0]);
        pv2[iq2][iPt]->SetY1NDC(0.74);
        pv2[iq2][iPt]->SetY2NDC(0.89);
        pv2fs[iq2][iPt] = (TPaveStats*)hResv2fs[iq2][iPt]->FindObject("stats");
        pv2fs[iq2][iPt]->SetTextColor(colors[1]);
        pv2fs[iq2][iPt]->SetY1NDC(0.59);
        pv2fs[iq2][iPt]->SetY2NDC(0.74);
        pv2BC1[iq2][iPt] = (TPaveStats*)hResv2BC1[iq2][iPt]->FindObject("stats");
        pv2BC1[iq2][iPt]->SetTextColor(colors[2]);
        pv2BC1[iq2][iPt]->SetY1NDC(0.44);
        pv2BC1[iq2][iPt]->SetY2NDC(0.59);
        pv2BC2[iq2][iPt] = (TPaveStats*)hResv2BC2[iq2][iPt]->FindObject("stats");
        pv2BC2[iq2][iPt]->SetTextColor(colors[3]);
        pv2BC2[iq2][iPt]->SetY1NDC(0.29);
        pv2BC2[iq2][iPt]->SetY2NDC(0.44);
        cv2[iq2]->cd(iPt+1)->Modified();
        
        cv2VsTrial[iq2]->cd(iPt+1);
        hv2VsTrial[iq2][iPt]->GetYaxis()->SetRangeUser(-0.4,0.6);
        hv2VsTrial[iq2][iPt]->SetStats(kFALSE);
        hv2fsVsTrial[iq2][iPt]->SetStats(kFALSE);
        hv2BC1VsTrial[iq2][iPt]->SetStats(kFALSE);
        hv2BC2VsTrial[iq2][iPt]->SetStats(kFALSE);
        if(iPt==0) {
          Double_t y1=0.7;
          Double_t y2=0.89;
          if(iq2==kLarge) {
            y1=0.15;
            y2=0.34;
          }
          methlegend[iq2] = new TLegend(0.5,y1,0.89,y2);
          methlegend[iq2]->SetFillStyle(0);
          methlegend[iq2]->SetFillColor(0);
          methlegend[iq2]->SetTextSize(0.05);
          methlegend[iq2]->AddEntry(hv2VsTrial[iq2][iPt],"free sigma","lpe");
          methlegend[iq2]->AddEntry(hv2fsVsTrial[iq2][iPt],"fixed sigma","lpe");
          methlegend[iq2]->AddEntry(hv2BC1VsTrial[iq2][iPt],"bin counting 1","lpe");
          methlegend[iq2]->AddEntry(hv2BC2VsTrial[iq2][iPt],"bin counting 2","lpe");
        }
        hv2VsTrial[iq2][iPt]->Draw();
        hv2BC1VsTrial[iq2][iPt]->Draw("E1same");
        hv2BC2VsTrial[iq2][iPt]->Draw("E1same");
        hv2fsVsTrial[iq2][iPt]->Draw("E1same");
        methlegend[iq2]->Draw("same");
        
        cv2Syst[iq2]->cd(iPt+1);
        Double_t v2refstaterr = gv2Ref[iq2]->GetErrorYhigh(iPt);
        v2refline[iq2][iPt] = new TLine(-0.5,0,3.5,0);
        v2refline[iq2][iPt]->SetLineWidth(2);
        v2refline[iq2][iPt]->SetLineColor(colors[1]);
        v2reflinelow[iq2][iPt] = new TLine(-0.5,-v2refstaterr,3.5,-v2refstaterr);
        v2reflinelow[iq2][iPt]->SetLineWidth(2);
        v2reflinelow[iq2][iPt]->SetLineColor(colors[1]);
        v2reflinelow[iq2][iPt]->SetLineStyle(7);
        v2reflinehigh[iq2][iPt] = new TLine(-0.5,v2refstaterr,3.5,v2refstaterr);
        v2reflinehigh[iq2][iPt]->SetLineWidth(2);
        v2reflinehigh[iq2][iPt]->SetLineColor(colors[1]);
        v2reflinehigh[iq2][iPt]->SetLineStyle(7);
        leg[iq2][iPt] = new TLegend(0.45,0.7,0.89,0.89);
        leg[iq2][iPt]->SetFillStyle(0);
        leg[iq2][iPt]->AddEntry(v2reflinelow[iq2][iPt],"statistical uncertainty","l");
        leg[iq2][iPt]->AddEntry(hv2Syst[iq2][iPt],"RMS","p");
        leg[iq2][iPt]->AddEntry(hv2Syst2[iq2][iPt],"(x_{max}-x_{min})/#sqrt{12}","p");
        hv2Syst[iq2][iPt]->GetYaxis()->SetRangeUser(-v2refstaterr*2,v2refstaterr*3);
        hv2Syst[iq2][iPt]->Draw("E1X0");
        hv2fsSyst[iq2][iPt]->Draw("E1sameX0");
        hv2BC1Syst[iq2][iPt]->Draw("E1sameX0");
        hv2BC2Syst[iq2][iPt]->Draw("E1sameX0");
        hv2Syst2[iq2][iPt]->Draw("E1sameX0");
        hv2fsSyst2[iq2][iPt]->Draw("E1sameX0");
        hv2BC1Syst2[iq2][iPt]->Draw("E1sameX0");
        hv2BC2Syst2[iq2][iPt]->Draw("E1sameX0");
        v2refline[iq2][iPt]->Draw("same");
        v2reflinelow[iq2][iPt]->Draw("same");
        v2reflinehigh[iq2][iPt]->Draw("same");
        leg[iq2][iPt]->Draw("same");
      }
      
      cResv2CorrLargeSmall->cd(iPt+1);
      hResv2CorrLargeSmall[iPt]->Draw("colz");
      latex->SetTextColor(kBlack);
      latex->DrawLatex(-4*hResv2[0][iPt]->GetRMS(),hResv2[0][iPt]->GetRMS()*4,Form("Correlation coefficient = %0.3f",hResv2CorrLargeSmall[iPt]->GetCorrelationFactor()));
      
      cResv2CorrLargeUnbiased->cd(iPt+1);
      hResv2CorrLargeUnbiased[iPt]->Draw("colz");
      latex->SetTextColor(kBlack);
      latex->DrawLatex(-4*hResv2[1][iPt]->GetRMS(),hResv2[1][iPt]->GetRMS()*4,Form("Correlation coefficient = %0.3f",hResv2CorrLargeUnbiased[iPt]->GetCorrelationFactor()));
      
      cResv2CorrSmallUnbiased->cd(iPt+1);
      hResv2CorrSmallUnbiased[iPt]->Draw("colz");
      latex->SetTextColor(kBlack);
      latex->DrawLatex(-4*hResv2[0][iPt]->GetRMS(),hResv2[0][iPt]->GetRMS()*4,Form("Correlation coefficient = %0.3f",hResv2CorrSmallUnbiased[iPt]->GetCorrelationFactor()));
      
      TF1* fLineLargeSmall = new TF1("fLineLargeSmall","pol1",0,nTrials);
      cv2fsRatioLargeSmallVsTrial->cd(iPt+1);
      hv2fsRatioLargeSmallVsTrial[iPt]->Fit("fLineLargeSmall");
      hv2fsRatioLargeSmallVsTrial[iPt]->Draw("E1");
      latex->SetTextColor(kRed);
      latex->DrawLatex(hv2fsRatioLargeSmallVsTrial[iPt]->GetNbinsX()/10,hv2fsRatioLargeSmallVsTrial[iPt]->GetMaximum()*0.9,Form("slope = %f #pm %f",fLineLargeSmall->GetParameter(1),fLineLargeSmall->GetParError(1)));
      
      TF1* fLineLargeUnbiased = new TF1("fLineLargeUnbiased","pol1",0,nTrials);
      cv2fsRatioLargeUnbiasedVsTrial->cd(iPt+1);
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->Fit("fLineLargeUnbiased");
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->Draw("E1");
      latex->DrawLatex(hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetNbinsX()/10,hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetMaximum()*0.9,Form("slope = %f #pm %f",fLineLargeUnbiased->GetParameter(1),fLineLargeUnbiased->GetParError(1)));
      
      TF1* fLineSmallUnbiased = new TF1("fLineSmallUnbiased","pol1",0,nTrials);
      cv2fsRatioSmallUnbiasedVsTrial->cd(iPt+1);
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->Fit("fLineSmallUnbiased");
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->Draw("E1");
      latex->DrawLatex(hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetNbinsX()/10,hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetMaximum()*0.9,Form("slope = %f #pm %f",fLineSmallUnbiased->GetParameter(1),fLineSmallUnbiased->GetParError(1)));
      
      hv2fsRatioLargeSmall[iPt] = new TH1F(Form("v2fsRatioLargeSmall_pt%d",iPt),";v_{2} (Large-q_{2}) / v_{2} (Small-q_{2}) (sigma fix);Entries",100,-1,-1);
      for(Int_t iTrial=0; iTrial<hv2fsRatioLargeSmallVsTrial[iPt]->GetNbinsX(); iTrial++) {
        if(hv2fsRatioLargeSmallVsTrial[iPt]->GetBinContent(iTrial+1)!=0) hv2fsRatioLargeSmall[iPt]->Fill(hv2fsRatioLargeSmallVsTrial[iPt]->GetBinContent(iTrial+1));
      }
      cv2fsRatioLargeSmall->cd(iPt+1);
      hv2fsRatioLargeSmall[iPt]->Draw();
      
      hv2fsRatioLargeUnbiased[iPt] = new TH1F(Form("v2fsRatioLargeUnbiased_pt%d",iPt),";v_{2} (Large-q_{2}) / v_{2} (q_{2}-Integrated) (sigma fix);Entries",100,-1,-1);
      for(Int_t iTrial=0; iTrial<hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetNbinsX(); iTrial++) {
        if(hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetBinContent(iTrial+1)!=0) hv2fsRatioLargeUnbiased[iPt]->Fill(hv2fsRatioLargeUnbiasedVsTrial[iPt]->GetBinContent(iTrial+1));
      }
      cv2fsRatioLargeUnbiased->cd(iPt+1);
      hv2fsRatioLargeUnbiased[iPt]->Draw();
      
      hv2fsRatioSmallUnbiased[iPt] = new TH1F(Form("v2fsRatioSmallUnbiased_pt%d",iPt),";v_{2} (Small-q_{2}) / v_{2} (q_{2}-Integrated) (sigma fix);Entries",100,-1,-1);
      for(Int_t iTrial=0; iTrial<hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetNbinsX(); iTrial++) {
        if(hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetBinContent(iTrial+1)!=0) hv2fsRatioSmallUnbiased[iPt]->Fill(hv2fsRatioSmallUnbiasedVsTrial[iPt]->GetBinContent(iTrial+1));
      }
      cv2fsRatioSmallUnbiased->cd(iPt+1);
      hv2fsRatioSmallUnbiased[iPt]->Draw();
      
    }
    
    TLegend* q2leg = new TLegend(0.6,0.6,0.89,0.89);
    q2leg->SetFillStyle(0);
    q2leg->SetFillColor(0);
    q2leg->SetTextSize(0.045);
    TString legnames[3] = {"Small-q_{2}","Large-q_{2}","q_{2}-integrated"};
    cShiftRMS->cd();
    for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
      TString drawopt="hist";
      if(iq2>kSmall) {drawopt+=" same";}
      hShiftRMS[iq2]->Draw(drawopt);
      q2leg->AddEntry(hShiftRMS[iq2],legnames[iq2].Data(),"l");
    }
    q2leg->Draw("same");
    
    //Output file
    TString outname=Form("%s/v2RawYieldSyst_%d_%d_%s%s_%s%0.2f_%s%0.2f.root",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2regionname[0].Data(),q2smalllimit,q2regionname[1].Data(),q2largelimit);
    TString percsuffix="perc";
    if(cutmeth==kPercCutVsCent) percsuffix="percVsCent";
    else if(cutmeth==kPercentileCut) percsuffix="percentile";
    if(cutmeth==kPercCut || cutmeth==kPercCutVsCent || cutmeth==kPercentileCut) {outname=Form("%s/v2RawYieldSyst_%d_%d_%s%s_%s%0.f%s_%s%0.f%s.root",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),q2regionname[0].Data(),q2smallpercevents*100,percsuffix.Data(),q2regionname[1].Data(),q2largepercevents*100,percsuffix.Data());}
    TFile *fout=new TFile(outname.Data(),"RECREATE"); //outputfile
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
        hResv2[iq2][iPt]->Write();
        hResv2fs[iq2][iPt]->Write();
        hResv2BC1[iq2][iPt]->Write();
        hResv2BC2[iq2][iPt]->Write();
        hv2VsTrial[iq2][iPt]->Write();
        hv2fsVsTrial[iq2][iPt]->Write();
        hv2BC1VsTrial[iq2][iPt]->Write();
        hv2BC2VsTrial[iq2][iPt]->Write();
      }
      hResv2CorrLargeSmall[iPt]->Write();
      hResv2CorrLargeUnbiased[iPt]->Write();
      hResv2CorrSmallUnbiased[iPt]->Write();
      hv2fsRatioLargeSmallVsTrial[iPt]->Write();
      hv2fsRatioLargeUnbiasedVsTrial[iPt]->Write();
      hv2fsRatioSmallUnbiasedVsTrial[iPt]->Write();
      hv2fsRatioLargeSmall[iPt]->Write();
      hv2fsRatioLargeUnbiased[iPt]->Write();
      hv2fsRatioSmallUnbiased[iPt]->Write();
    }
    for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
      hShiftRMS[iq2]->Write();
      cv2[iq2]->Write();
      cv2Syst[iq2]->Write();
      cv2VsTrial[iq2]->Write();
    }
    cResv2CorrLargeSmall->Write();
    cv2fsRatioLargeSmallVsTrial->Write();
    cShiftRMS->Write();
    cv2fsRatioLargeSmallVsTrial->Write();
    cv2fsRatioLargeUnbiasedVsTrial->Write();
    cv2fsRatioSmallUnbiasedVsTrial->Write();
    cv2fsRatioLargeSmall->Write();
    cv2fsRatioLargeUnbiased->Write();
    cv2fsRatioSmallUnbiased->Write();
    treeMultiTrial->Write();
    fout->Close();
    
    outname.ReplaceAll(".root",".pdf");
    outname.ReplaceAll("/v2RawYieldSyst","/v2RawYieldSyst_RMSshift");
    cShiftRMS->SaveAs(outname.Data());
    outname.ReplaceAll("v2RawYieldSyst_RMSshift","v2Corrq2Regions");
    cResv2CorrLargeSmall->SaveAs(outname.Data());
    outname.ReplaceAll("v2Corrq2Regions","v2CorrLargeUnbiased");
    cResv2CorrLargeUnbiased->SaveAs(outname.Data());
    outname.ReplaceAll("v2CorrLargeUnbiased","v2CorrSmallUnbiased");
    cResv2CorrSmallUnbiased->SaveAs(outname.Data());
    outname.ReplaceAll("v2CorrSmallUnbiased","v2Ratioq2Regions");
    cv2fsRatioLargeSmallVsTrial->SaveAs(outname.Data());
    
    TString plotappendix[3] = {Form("%s%0.2f",q2regionname[0].Data(),q2smalllimit),Form("%s%0.2f",q2regionname[1].Data(),q2largelimit),Form("%s",q2regionname[2].Data())};
    if(cutmeth==kPercCut || cutmeth==kPercCutVsCent || cutmeth==kPercentileCut) {
      plotappendix[0] = Form("%s%0.f%s",q2regionname[0].Data(),q2smallpercevents*100,percsuffix.Data());
      plotappendix[1] = Form("%s%0.f%s",q2regionname[1].Data(),q2largepercevents*100,percsuffix.Data());
    }
    for(Int_t iq2=kSmall; iq2<=kIntegrated; iq2++) {
      cv2[iq2]->SaveAs(Form("%s/v2ResidRawYieldSyst_%d_%d_%s%s_%s.pdf",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),plotappendix[iq2].Data()));
      cv2Syst[iq2]->SaveAs(Form("%s/v2RawYieldSyst_%d_%d_%s%s_%s.pdf",outputdir.Data(),minCent,maxCent,analysismethname.Data(),suffix.Data(),plotappendix[iq2].Data()));
    }
  }
  else {cerr<< "Scalar product not yet implemented. Exit." << endl; return 1;}
  
  return 0;
}

//_____________________________________________________________________________________________
//LOAD DATA LIST
TList* LoadTList() {
  
  cout << "Opening input file " << infilename << "..." <<endl;
  TFile* infile = TFile::Open(infilename.Data(),"READ");
  TDirectoryFile* dir=0x0;
  TList* list=0x0;
  
  TString dirname=Form("PWGHF_D2H_HFv2_%s%s",partname.Data(),suffix.Data());
  TString listname=Form("coutputv2%s%s",partname.Data(),suffix.Data());
  
  if(infile) {dir=(TDirectoryFile*)infile->Get(dirname.Data()); cout << "File opened!" << endl;}
  else {cerr << "Error: File " << infilename << " not found. Exit." << endl; return 0x0;}
  if(dir) list=(TList*)dir->Get(listname.Data());
  else {
    dirname.ReplaceAll("HFv2","HFvn");
    dir=(TDirectoryFile*)infile->Get(dirname.Data()); cout << "File opened!" << endl;
    list=(TList*)dir->Get(listname.Data());
  }
  if(!dir) {cerr << "Error: Wrong TDirectoryFile name " << dirname << ". Exit." << endl; return 0x0;}
  if(!list) {cerr << "Error: Wrong TList name " << listname << ". Exit." << endl; return 0x0;}
  
  infile->Close();
  cout<<"Input file closed."<< endl;
  
  return list;
}

//_____________________________________________________________________________________________
//GET SPARSE
THnSparseF* LoadSparseFromList(TList* inputlist) {
  
  THnSparseF* sparse=(THnSparseF*)inputlist->FindObject("hMassPtPhiq2Centr");
  if(!sparse) {cerr << "Error: No THnSparse found. Check if the name hMassPtPhiq2Centr is right!" << endl; return 0x0;}
  else {cout << "THnSparse got!" << endl;}
  
  return sparse;
}

//_____________________________________________________________________________________________
//GET q2 VS. CENTRALITY HISTO
TH2F* GetHistoq2VsCentr(TList* inputlist) {
  
  Double_t centwidth = (Double_t)(maxCent-minCent)/ncentbins;
  Double_t mincentpermil=minCent*10;
  Double_t maxcentpermil=(minCent+centwidth)*10;
  
  TH2F* hResVsq2=(TH2F*)inputlist->FindObject(Form("hEvPlaneReso1Vsq2centr%0.f_%0.f",mincentpermil,maxcentpermil));
  if(!hResVsq2) {return 0x0;}
  TH1F* hq2=(TH1F*)hResVsq2->ProjectionY();
  
  Int_t nq2bins = hq2->GetNbinsX();
  Int_t q2min = hq2->GetBinLowEdge(1);
  Int_t q2max = hq2->GetBinLowEdge(nq2bins)+hq2->GetBinWidth(nq2bins);
  
  TH2F* hq2VsCentr = new TH2F("hq2VsCentr","q_{2} vs. centrality",ncentbins,minCent,maxCent,nq2bins,q2min,q2max);
  hq2VsCentr->GetXaxis()->SetTitle("Centrality (%)");
  hq2VsCentr->GetYaxis()->SetTitle(hq2->GetXaxis()->GetTitle());
  
  for(Int_t iBin=0; iBin<nq2bins+1; iBin++) {
    hq2VsCentr->SetBinContent(1,iBin+1,hq2->GetBinContent(iBin+1)); //takes also overflow counts
  }
  
  for(Int_t iCent=1; iCent<ncentbins; iCent++) {
    mincentpermil=(minCent+iCent*centwidth)*10;
    maxcentpermil=(minCent+(iCent+1)*centwidth)*10;
    
    hResVsq2=(TH2F*)inputlist->FindObject(Form("hEvPlaneReso1Vsq2centr%0.f_%0.f",mincentpermil,maxcentpermil));
    if(!hResVsq2) return 0x0;
    hq2 = (TH1F*)hResVsq2->ProjectionY();
    for(Int_t iBin=0; iBin<hq2->GetNbinsX()+1; iBin++) {
      hq2VsCentr->SetBinContent(iCent+1,iBin+1,hq2->GetBinContent(iBin+1)); //takes also overflow counts
    }
    hq2VsCentr->GetYaxis()->SetTitle(hq2->GetXaxis()->GetTitle());
  }
  hq2VsCentr->SetStats(0);
  
  return hq2VsCentr;
}

//_____________________________________________________________________________________________
//GET MASS HISTOS IN SELECTED q2 RANGE
TList* LoadMassHistos(THnSparseF* sparse, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge, Int_t analysismeth) {
  
  Int_t ptaxis=1;
  Int_t phiaxis=2;
  Int_t q2axis=3;
  Int_t centaxis=4;
  
  Double_t centwidth = (Double_t)(maxCent-minCent)/ncentbins;
  
  TList *outlist = new TList();
  outlist->SetName("azimuthalhistoslist");
  
  TH1F* hMassTmp[ncentbins];
  TH1F* hMassTmp1[ncentbins];
  TH1F* hMassTmp2[ncentbins];
  
  Int_t fSparseVers = 0; //check version of the task output
  if(sparse->GetAxis(centaxis)) {
    TString centtitle=sparse->GetAxis(centaxis)->GetTitle();
    if(centtitle.Contains("Centrality") || centtitle.Contains("centrality")) {
      cout << "Version of the sparse with centrality axis!" << endl;
      fSparseVers=1;
    }
  }
  
  for(Int_t iPt=0; iPt<nPtBins; iPt++) {
    //get phi- and q2-integrated histos (for fixed sigma)
    ApplyCut(sparse,PtLims[iPt],PtLims[iPt+1],ptaxis); //apply pt cut
    hMassTmp[0]=(TH1F*)sparse->Projection(0);
    hMassTmp[0]->SetName(Form("hMass_pt%d",iPt));
    outlist->Add(hMassTmp[0]->Clone());
    delete hMassTmp[0];
    hMassTmp[0]=0x0;
    
    //get histos in deltaphi bins and with the q2 cut
    if(analysismeth==kEventPlane) {
      for(Int_t iPhi=0; iPhi<nPhiBins; iPhi++) {
        ApplyCut(sparse,PhiLims[iPhi],PhiLims[iPhi+1],phiaxis); //apply deltaphi cut
        if(fSparseVers==0) {
          if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[0],largecutvalues[0],smallorlarge);} //apply q2 cut
          hMassTmp[0]=(TH1F*)sparse->Projection(0);
        }
        else if(fSparseVers==1) {
          for(Int_t iCent=0; iCent<ncentbins; iCent++) {
            ApplyCut(sparse,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,centaxis); //apply centrality cut (to have a q2 cut centrality dependent)
            if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[iCent],largecutvalues[iCent],smallorlarge);} //apply q2 cut centrality dependent
            hMassTmp[iCent]=(TH1F*)sparse->Projection(0);
            if(iCent>0) {hMassTmp[0]->Add(hMassTmp[iCent]);}
            ResetAxes(sparse,centaxis); //release centrality cut
            ResetAxes(sparse,q2axis); //release q2 cut
            if(iCent!=0) {
              delete hMassTmp[iCent];
              hMassTmp[iCent]=0x0;
            }
          }
        }
        hMassTmp[0]->SetName(Form("hMass_pt%d_phi%d",iPt,iPhi));
        outlist->Add(hMassTmp[0]->Clone());
        ResetAxes(sparse,q2axis); //release q2 cut
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        delete hMassTmp[0];
        hMassTmp[0]=0x0;
      }
    }
    else if(analysismeth==kEventPlaneInOut){
      if(fSparseVers==0) {
        if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[0],largecutvalues[0],smallorlarge);} //apply q2 cut
        ApplyCut(sparse,0.,TMath::Pi()/4,phiaxis); //apply deltaphi cut
        hMassTmp[0]=(TH1F*)sparse->Projection(0);
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        ApplyCut(sparse,3./4*TMath::Pi(),TMath::Pi(),phiaxis); //apply deltaphi cut
        hMassTmp1[0]=(TH1F*)sparse->Projection(0);
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        hMassTmp[0]->Add(hMassTmp1[0]);
        hMassTmp[0]->SetName(Form("hMass_pt%d_phi0",iPt));
        outlist->Add(hMassTmp[0]->Clone());
        ApplyCut(sparse,TMath::Pi()/4,3./4*TMath::Pi(),phiaxis); //apply deltaphi cut
        hMassTmp2[0]=(TH1F*)sparse->Projection(0);
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        ResetAxes(sparse,q2axis); //release q2 cut
        hMassTmp2[0]->SetName(Form("hMass_pt%d_phi1",iPt));
        outlist->Add(hMassTmp2[0]->Clone());
      }
      
      else if(fSparseVers==1) {
        ApplyCut(sparse,0.,TMath::Pi()/4,phiaxis); //apply deltaphi cut
        for(Int_t iCent=0; iCent<ncentbins; iCent++) {
          ApplyCut(sparse,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,centaxis); //apply centrality cut (to have a q2 cut centrality dependent)
          if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[iCent],largecutvalues[iCent],smallorlarge);} //apply q2 cut centrality dependent
          hMassTmp[iCent]=(TH1F*)sparse->Projection(0);
          if(iCent>0) {hMassTmp[0]->Add(hMassTmp[iCent]);}
          ResetAxes(sparse,centaxis); //release centrality cut
          ResetAxes(sparse,q2axis); //release q2 cut
          if(iCent!=0) {
            delete hMassTmp[iCent];
            hMassTmp[iCent]=0x0;
          }
        }
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        ApplyCut(sparse,3./4*TMath::Pi(),TMath::Pi(),phiaxis); //apply deltaphi cut
        for(Int_t iCent=0; iCent<ncentbins; iCent++) {
          ApplyCut(sparse,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,centaxis); //apply centrality cut (to have a q2 cut centrality dependent)
          if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[iCent],largecutvalues[iCent],smallorlarge);} //apply q2 cut centrality dependent
          hMassTmp1[iCent]=(TH1F*)sparse->Projection(0);
          if(iCent>0) {hMassTmp1[0]->Add(hMassTmp1[iCent]);}
          ResetAxes(sparse,centaxis); //release centrality cut
          ResetAxes(sparse,q2axis); //release q2 cut
          if(iCent!=0) {
            delete hMassTmp1[iCent];
            hMassTmp1[iCent]=0x0;
          }
        }
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        hMassTmp[0]->Add(hMassTmp1[0]);
        hMassTmp[0]->SetName(Form("hMass_pt%d_phi0",iPt));
        outlist->Add(hMassTmp[0]->Clone());
        ApplyCut(sparse,TMath::Pi()/4,3./4*TMath::Pi(),phiaxis); //apply deltaphi cut
        for(Int_t iCent=0; iCent<ncentbins; iCent++) {
          ApplyCut(sparse,minCent+(iCent*centwidth),minCent+(iCent*centwidth)+1,centaxis); //apply centrality cut (to have a q2 cut centrality dependent)
          if(smallorlarge!=kIntegrated) {Applyq2Cut(sparse,smallcutvalues[iCent],largecutvalues[iCent],smallorlarge);} //apply q2 cut centrality dependent
          hMassTmp2[iCent]=(TH1F*)sparse->Projection(0);
          if(iCent>0) {hMassTmp2[0]->Add(hMassTmp2[iCent]);}
          ResetAxes(sparse,centaxis); //release centrality cut
          ResetAxes(sparse,q2axis); //release q2 cut
          if(iCent!=0) {
            delete hMassTmp2[iCent];
            hMassTmp2[iCent]=0x0;
          }
        }
        ResetAxes(sparse,phiaxis); //release deltaphi cut
        hMassTmp2[0]->SetName(Form("hMass_pt%d_phi1",iPt));
        outlist->Add(hMassTmp2[0]->Clone());
      }
      delete hMassTmp[0];
      delete hMassTmp1[0];
      delete hMassTmp2[0];
      hMassTmp[0]=0x0;
      hMassTmp1[0]=0x0;
      hMassTmp2[0]=0x0;
      ResetAxes(sparse,ptaxis); //release pt cut
    }
  }
  ResetAxes(sparse,-1); //release all cuts
  
  if(useTemplD0Refl){
    Bool_t retCode=LoadD0toKpiMCHistos(outlist);
    if(!retCode)Printf("ERROR: MC histograms loading failed");
    else Printf("******************************************\n MC HISTOGRAMS LOADED\n\n**********************************");
  }
  
  return outlist;
}

//_____________________________________________________________________________________________
//GET RESOLUTION HISTOS IN SELECTED q2 RANGE
TList* LoadResolutionHistos(TList *inputlist, vector<Double_t> smallcutvalues, vector<Double_t> largecutvalues, Int_t smallorlarge) {
  
  TList *outlist = new TList();
  outlist->SetName("eventplanemasslist");
  
  const Int_t nBins=20;
  Double_t ncoll[nBins]={1790.77,1578.44,1394.82,1236.17
    ,1095.08,969.86,859.571,759.959,669.648,589.588,516.039
    ,451.409,392.853,340.493,294.426,252.385,215.484,183.284
    ,155.101,130.963};
  
  Double_t centwidth = (Double_t)(maxCent-minCent)/ncentbins;
  
  Double_t minCentTimesTen=minCent*10;
  Double_t maxCentTimesTen=maxCent*10;
  Double_t centwidthTimesTen=centwidth*10;
  
  TGraphErrors* gResolVsCent=new TGraphErrors(0);
  Int_t iCent=0;
  Int_t nSubRes=1;
  if(evPlaneRes==AliEventPlaneResolutionHandler::kThreeSub||
     evPlaneRes==AliEventPlaneResolutionHandler::kThreeSubTPCGap)nSubRes=3;
  TString namereso[3]={"ResoVsq2","Reso2Vsq2","Reso3Vsq2"};
  TString suffixcentr=Form("centr%d_%d",minCent,maxCent);
  TH2F* htestversion=(TH2F*)inputlist->FindObject(Form("hEvPlaneResocentr%0.f_%0.f",minCentTimesTen,minCentTimesTen+centwidthTimesTen));
  if(htestversion){
    printf("Old version of the task\n");
  }else{
    printf("New version of the task\n");
    namereso[0]="Reso1Vsq2";
  }
  TH1F* hevplresos[3];
  
  for(Int_t iCent=0; iCent<ncentbins; iCent++){
    TH2F* htmpresosvsq2[3];
    TH1F* htmpresos[3];
    
    for(Int_t iRes=0;iRes<nSubRes;iRes++){
      htmpresosvsq2[iRes]=(TH2F*)inputlist->FindObject(Form("hEvPlane%scentr%0.f_%0.f",namereso[iRes].Data(),minCentTimesTen+(iCent*centwidthTimesTen),minCentTimesTen+((iCent+1)*centwidthTimesTen)));
      Int_t binmin=-1;
      Int_t binmax=-1;
      
      //define bins for the q2 cut
      if(smallorlarge==kSmall) {
        binmin=1;
        binmax=htmpresosvsq2[iRes]->GetYaxis()->FindBin(smallcutvalues[iCent]-0.0001);
      }
      else if(smallorlarge==kLarge) {
        binmin=htmpresosvsq2[iRes]->GetYaxis()->FindBin(largecutvalues[iCent]+0.0001);
        binmax=htmpresosvsq2[iRes]->GetYaxis()->FindBin(htmpresosvsq2[iRes]->GetYaxis()->GetNbins()+1); //takes also overflow entries
        Int_t binmaxsmallregion=htmpresosvsq2[iRes]->GetYaxis()->FindBin(smallcutvalues[iCent]-0.0001);
        if(binmaxsmallregion==binmin) {
          binmin += 1; //avoid overlap between the two q2-regions
        }
      }
      else if(smallorlarge==kIntegrated) {
        binmin=1;
        binmax=htmpresosvsq2[iRes]->GetYaxis()->FindBin(htmpresosvsq2[iRes]->GetYaxis()->GetNbins()+1); //takes also overflow entries
      }
      htmpresos[iRes]=(TH1F*)htmpresosvsq2[iRes]->ProjectionX(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()),binmin,binmax);
      if(htmpresos[iRes]){
        TString suffixcentrpart=Form("centr%0.1f_%0.1f",minCent+(iCent*centwidth),minCent+((iCent+1)*centwidth));
        htmpresos[iRes]->SetTitle(Form("Event Plane Resolution %s%s",namereso[iRes].Data(),suffixcentrpart.Data()));
        htmpresos[iRes]->SetName(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentrpart.Data()));
        if(useNcollWeight){
          cout << Form("Centr %0.f Bin %d  Ncoll %f\n",minCentTimesTen,ncentbins,ncoll[ncentbins]) << endl;
          htmpresos[iRes]->Scale(ncoll[ncentbins]);
        }
        outlist->Add(htmpresos[iRes]->Clone());
        if(iCent==0) {
          hevplresos[iRes]=(TH1F*)htmpresos[iRes]->Clone();
        }
        else {
          hevplresos[iRes]->Add(htmpresos[iRes]);
        }
      }
    }
    
    Double_t error;
    Double_t lowestRes=1;
    Double_t highestRes=0;
    Double_t resolBin=GetEventPlaneResolution(error,htmpresos[0],htmpresos[1],htmpresos[2]);
    
    if(resolBin<lowestRes) lowestRes=resolBin;
    if(resolBin>highestRes) highestRes=resolBin;
    
    Double_t binHalfWid=centwidth/2.;
    Double_t binCentr=(Double_t)minCent+(iCent*centwidth)+binHalfWid;
    gResolVsCent->SetPoint(iCent,binCentr,resolBin);
    gResolVsCent->SetPointError(iCent,binHalfWid,error);
  }
  
  for(Int_t iRes=0;iRes<nSubRes;iRes++){
    if(hevplresos[iRes]) {
      hevplresos[iRes]->SetTitle(Form("Event Plane Resolution %s%s",namereso[iRes].Data(),suffixcentr.Data()));
      hevplresos[iRes]->SetName(Form("hEvPlane%s%s",namereso[iRes].Data(),suffixcentr.Data()));
      outlist->Add(hevplresos[iRes]->Clone());
    }
  }
  gResolVsCent->SetName("gResolVsCent");
  gResolVsCent->SetTitle("Resolution vs. Centrality");
  outlist->Add(gResolVsCent->Clone());
  return outlist;
}

TGraphAsymmErrors* Computev2(TGraphAsymmErrors **gSignal, Double_t resol, Float_t *averagePt, Int_t analysismeth, TGraphAsymmErrors *gRelSystEff) {
  
  TGraphAsymmErrors* gv2 = new TGraphAsymmErrors(nPtBins);
  
  if(analysismeth==kEventPlaneInOut) {
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      Double_t *y=gSignal[iPt]->GetY();
      Double_t nIn=y[0];
      Double_t nOut=y[1];
      Double_t enIn=gSignal[iPt]->GetErrorY(0);
      Double_t enOut=gSignal[iPt]->GetErrorY(1);
      Double_t anis=(nIn-nOut)/(nIn+nOut);
      Double_t eAnis=2./((nIn+nOut)*(nIn+nOut))*TMath::Sqrt(nIn*nIn*enOut*enOut+nOut*nOut*enIn*enIn);
      Double_t v2=anis*TMath::Pi()/4./resol;
      Double_t ev2=eAnis*TMath::Pi()/4./resol;
      gv2->SetPoint(iPt,averagePt[iPt],v2);
      gv2->SetPointError(iPt,averagePt[iPt]-PtLims[iPt],PtLims[iPt+1]-averagePt[iPt],ev2,ev2);
      if(gRelSystEff) {
        //systematic uncertainty for in-out efficiency
        Double_t anis1=(nIn-nOut*effInOverEffOut)/(nIn+nOut*effInOverEffOut);
        Double_t anis2=(nIn*effInOverEffOut-nOut)/(nIn*effInOverEffOut+nOut);
        Double_t systEffUp=0.,systEffDown=0.;
        if(anis1>anis && anis1>anis2) systEffUp=anis1/anis;
        if(anis2>anis && anis2>anis1) systEffUp=anis2/anis;
        if(anis1<anis && anis1<anis2) systEffDown=anis1/anis;
        if(anis2<anis && anis2<anis1) systEffDown=anis2/anis;
        cout << Form(" Bin %d <pt>=%.3f  v2=%f+-%f systEff=%f %f\n",iPt,averagePt[iPt],v2,ev2,systEffUp*v2,systEffDown*v2)<<endl;
        gRelSystEff->SetPoint(iPt,averagePt[iPt],v2);
        gRelSystEff->SetPointError(iPt,0.4,0.4,v2*(1-systEffDown),v2*(systEffUp-1));
      }
    }
    return gv2;
  }
  else if(analysismeth==kEventPlane) {
    TF1 *flowFunc = new TF1("flow","[0]*(1.+2.*[1]*TMath::Cos(2.*x))");
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      //v2 from fit to Deltaphi distribution
      gSignal[iPt]->Fit(flowFunc);
      Double_t v2 = flowFunc->GetParameter(1)/resol;
      Double_t ev2=flowFunc->GetParError(1)/resol;
      gv2->SetPoint(iPt,averagePt[iPt],v2);
      gv2->SetPointError(iPt,averagePt[iPt]-PtLims[iPt],PtLims[iPt+1]-averagePt[iPt],ev2,ev2);
    }
    return gv2;
  }
  else {return 0x0;}
}

//_____________________________________________________________________________________________
//FILL SIGNAL GRAPHS
void FillSignalGraph(TList *masslist,TGraphAsymmErrors **gSignal,TGraphAsymmErrors **gSignalfs, TGraphAsymmErrors **gSignalBC1, TGraphAsymmErrors **gSignalBC2,TGraphAsymmErrors **gSigmaFree,TGraphAsymmErrors **gSigmaFixed, TGraphAsymmErrors **gChiSquareFree, TGraphAsymmErrors **gChiSquareFixed, Int_t analysismeth, Int_t bkgfunc, Double_t minfit, Double_t maxfit, Int_t rebin) {
  
  TFile* infile = TFile::Open(infilename.Data(),"READ");
  TDirectoryFile* dir=0x0;
  TList* list=0x0;
  Double_t massD=-1;
  
  TString dirname=Form("PWGHF_D2H_HFv2_%s%s",partname.Data(),suffix.Data());
  
  if(infile) {dir=(TDirectoryFile*)infile->Get(dirname.Data());}
  if(!dir) {
    dirname.ReplaceAll("HFv2","HFvn");
    dir=(TDirectoryFile*)infile->Get(dirname.Data());
  }
  if(dir) {
    if(partname.Contains("Dzero")) {
      massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
    }
    if(partname.Contains("Dplus")){
      massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
    }
    if(partname.Contains("Dstar")) {
      massD=(TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass());
    }
    if(partname.Contains("Ds") && !partname.Contains("Dstar")) {
      massD=(TDatabasePDG::Instance()->GetParticle(431)->Mass());
    }
  }
  else {cerr << "Error: Wrong TDirectoryFile name " << dirname << ". Exit." << endl;}
  
  Int_t nPhi=nPhiBins;
  if(analysismeth==kEventPlaneInOut) nPhi=2;
  
  TH1F *hrflTempl=0x0;
  TH1F *hsigMC=0x0;
  Float_t sOverRef=0.;
  
  Int_t nMassBins;
  Double_t hmin,hmax;
  for(Int_t iPt=0;iPt<nPtBins;iPt++){
    TH1F *histtofitfullsigma=(TH1F*)masslist->FindObject(Form("hMass_pt%d",iPt))->Clone();
    if(useTemplD0Refl){
      hrflTempl=(TH1F*)masslist->FindObject(Form("histRfl_%d",iPt));
      hsigMC=(TH1F*)masslist->FindObject(Form("histSgn_%d",iPt));
    }
    for(Int_t iPhi=0;iPhi<nPhi;iPhi++){
      Double_t signal=0,esignal=0;
      Double_t sigma=0, esigma=0;
      Double_t chisquare=0;
      TH1F *histtofit=(TH1F*)masslist->FindObject(Form("hMass_pt%d_phi%d",iPt,iPhi))->Clone();
      if(!histtofit){
        gSignal[iPt]->SetPoint(iPhi,iPhi,signal);
        gSignal[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
        return;
      }
      histtofit->SetTitle(Form("%.1f < #it{p}_{T} < %.1f, #phi%d",PtLims[iPt],PtLims[iPt+1],iPhi));
      nMassBins=histtofit->GetNbinsX();
      hmin=TMath::Max(minfit,histtofit->GetBinLowEdge(2));
      hmax=TMath::Min(maxfit,histtofit->GetBinLowEdge(nMassBins-2));
      histtofit->Rebin(rebin);
      AliHFMassFitterVAR* fitter=new AliHFMassFitterVAR(histtofit,hmin,hmax,1,bkgfunc,types);
      if(useTemplD0Refl){
        fitter->SetTemplateReflections(hrflTempl);
        sOverRef=(hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999)))/(hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999)));
        fitter->SetFixReflOverS(sOverRef,kTRUE);
      }
      fitter->SetInitialGaussianMean(massD);
      fitter->SetInitialGaussianSigma(0.012);
      if(partname.Contains("Dstar")) {
        fitter->SetInitialGaussianSigma(0.0004);
      }
      fitter->SetUseLikelihoodFit();
      Bool_t ok=fitter->MassFitter(kFALSE);
      Double_t sigmaforcounting=0;
      Double_t meanforcounting=0;
      if(ok){
        signal=fitter->GetRawYield();
        esignal=fitter->GetRawYieldError();
        sigma=fitter->GetSigma();
        esigma=fitter->GetSigmaUncertainty();
        sigmaforcounting=sigma;
        meanforcounting=fitter->GetMean();
        chisquare=fitter->GetReducedChiSquare();
      }
      gSignal[iPt]->SetPoint(iPhi,iPhi,signal);
      gSignal[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
      gSigmaFree[iPt]->SetPoint(iPhi,iPhi,sigma);
      gSigmaFree[iPt]->SetPointError(iPhi,0,0,esigma,esigma);
      gChiSquareFree[iPt]->SetPoint(iPhi,iPhi,chisquare);
      gChiSquareFree[iPt]->SetPointError(iPhi,0,0,0,0);
      TF1* fB1=fitter->GetBackgroundFullRangeFunc();
      TF1* fB2=fitter->GetBackgroundRecalcFunc();
      Double_t minBinSum=histtofit->FindBin(meanforcounting-nSigmaForCounting*sigmaforcounting);
      Double_t maxBinSum=histtofit->FindBin(meanforcounting+nSigmaForCounting*sigmaforcounting);
      Double_t cntSig1=0.;
      Double_t cntSig2=0.;
      Double_t cntErr=0.;
      for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
        Double_t bkg1=fB1 ? fB1->Eval(histtofit->GetBinCenter(iMB)) : 0;
        Double_t bkg2=fB2 ? fB2->Eval(histtofit->GetBinCenter(iMB)) : 0;
        cntSig1+=(histtofit->GetBinContent(iMB)-bkg1);
        cntSig2+=(histtofit->GetBinContent(iMB)-bkg2);
        cntErr+=(histtofit->GetBinContent(iMB));
      }
      cntErr=TMath::Sqrt(cntErr);
      gSignalBC1[iPt]->SetPoint(iPhi,iPhi,cntSig1);
      gSignalBC1[iPt]->SetPointError(iPhi,0,0,cntErr,cntErr);
      gSignalBC2[iPt]->SetPoint(iPhi,iPhi,cntSig2);
      gSignalBC2[iPt]->SetPointError(iPhi,0,0,cntErr,cntErr);
    }
    //fit for fixed sigma
    histtofitfullsigma->SetTitle(Form("%.1f < #it{p}_{T} < %.1f GeV/c",PtLims[iPt],PtLims[iPt+1]));
    nMassBins=histtofitfullsigma->GetNbinsX();
    hmin=TMath::Max(minfit,histtofitfullsigma->GetBinLowEdge(2));
    hmax=TMath::Min(maxfit,histtofitfullsigma->GetBinLowEdge(nMassBins-2));
    histtofitfullsigma->Rebin(rebin);
    AliHFMassFitterVAR* fitter=new AliHFMassFitterVAR(histtofitfullsigma,hmin,hmax,1,bkgfunc,types);
    if(useTemplD0Refl){
      fitter->SetTemplateReflections(hrflTempl);
      sOverRef=hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999))/hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999));
      fitter->SetFixReflOverS(sOverRef,kTRUE);
    }
    fitter->SetInitialGaussianMean(massD);
    fitter->SetUseLikelihoodFit();
    fitter->SetInitialGaussianSigma(0.012);
    if(partname.Contains("Dstar")) {
      fitter->SetInitialGaussianSigma(0.0004);
    }
    Bool_t ok=fitter->MassFitter(kFALSE);
    Double_t sigmatot=fitter->GetSigma();
    Double_t massFromFit=fitter->GetMean();
    for(Int_t iPhi=0;iPhi<nPhi;iPhi++){
      TH1F *histtofit=(TH1F*)masslist->FindObject(Form("hMass_pt%d_phi%d",iPt,iPhi))->Clone();
      histtofit->SetTitle(Form("%.1f < #it{p}_{T} < %.1f, #phi%d",PtLims[iPt],PtLims[iPt+1],iPhi));
      nMassBins=histtofit->GetNbinsX();
      hmin=TMath::Max(minfit,histtofit->GetBinLowEdge(2));
      hmax=TMath::Min(maxfit,histtofit->GetBinLowEdge(nMassBins-2));
      histtofit->Rebin(rebin);
      AliHFMassFitterVAR* fitter2=new AliHFMassFitterVAR(histtofit,hmin,hmax,1,bkgfunc,types);
      if(useTemplD0Refl){
        fitter2->SetTemplateReflections(hrflTempl);
        sOverRef=hrflTempl->Integral(hrflTempl->FindBin(hmin*1.0001),hrflTempl->FindBin(hmax*0.999))/hsigMC->Integral(hsigMC->FindBin(hmin*1.0001),hsigMC->FindBin(hmax*0.999));
        fitter2->SetFixReflOverS(sOverRef,kTRUE);
      }
      fitter2->SetInitialGaussianMean(massD);
      fitter2->SetFixGaussianSigma(sigmatot);
      if(fixAlsoMass) fitter2->SetFixGaussianMean(massFromFit);
      fitter2->SetUseLikelihoodFit();
      Bool_t ok2=fitter2->MassFitter(kFALSE);
      Double_t signal=0,esignal=0;
      Double_t sigma=0, esigma=0;
      Double_t chisquare=0;
      if(ok2){
        signal=fitter2->GetRawYield();
        esignal=fitter2->GetRawYieldError();
        sigma=fitter->GetSigma();
        esigma=fitter->GetSigmaUncertainty();
        chisquare=fitter2->GetReducedChiSquare();
      }
      gSignalfs[iPt]->SetPoint(iPhi,iPhi,signal);
      gSignalfs[iPt]->SetPointError(iPhi,0,0,esignal,esignal);
      gSigmaFixed[iPt]->SetPoint(iPhi,iPhi,sigma);
      gSigmaFixed[iPt]->SetPointError(iPhi,0,0,esigma,esigma);
      gChiSquareFixed[iPt]->SetPoint(iPhi,iPhi,chisquare);
      gChiSquareFixed[iPt]->SetPointError(iPhi,0,0,0,0);
    }
  }//end loop on pt bin
}

//___________________________________________________________
//LOAD REFERENCE GRAPH
Int_t LoadRefGraphs(TString reffilename, TH1F **hRawYieldRef, TH1F **hRawYieldfsRef, TH1F **hRawYieldBC1Ref, TH1F **hRawYieldBC2Ref, TGraphAsymmErrors *&gv2Ref, TGraphAsymmErrors *&gv2fsRef, TGraphAsymmErrors *&gv2BC1Ref, TGraphAsymmErrors *&gv2BC2Ref, Int_t smallorlarge, Int_t analysismeth) {
  
  Int_t nPhi=nPhiBins;
  if(analysismeth==kEventPlaneInOut) nPhi=2;
  
  TString q2regionname="q2Small";
  if(smallorlarge==kLarge) q2regionname="q2Large";
  if(smallorlarge==kIntegrated) q2regionname="q2Int";
  
  TFile* reffile = TFile::Open(reffilename.Data(),"READ");
  if(reffile) {
    gv2Ref=(TGraphAsymmErrors*)reffile->Get(Form("gav2_%s",q2regionname.Data()));
    gv2fsRef=(TGraphAsymmErrors*)reffile->Get(Form("gav2fs_%s",q2regionname.Data()));
    gv2BC1Ref=(TGraphAsymmErrors*)reffile->Get(Form("gav2BC1_%s",q2regionname.Data()));
    gv2BC2Ref=(TGraphAsymmErrors*)reffile->Get(Form("gav2BC2_%s",q2regionname.Data()));
    
    if(!gv2Ref) {return 2;}
    if(!gv2fsRef) {return 3;}
    if(!gv2BC1Ref) {return 4;}
    if(!gv2BC2Ref) {return 5;}
    
    for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
      hRawYieldRef[iPhi] = new TH1F(Form("hRawYieldRef_phi%d",iPhi),"",nPtBins,PtLims);
      hRawYieldfsRef[iPhi]= new TH1F(Form("hRawYieldfsRef_phi%d",iPhi),"",nPtBins,PtLims);
      hRawYieldBC1Ref[iPhi]= new TH1F(Form("hRawYieldBC1Ref_phi%d",iPhi),"",nPtBins,PtLims);
      hRawYieldBC2Ref[iPhi]= new TH1F(Form("hRawYieldBC2Ref_phi%d",iPhi),"",nPtBins,PtLims);
      hRawYieldRef[iPhi]->SetDirectory(0);
      hRawYieldfsRef[iPhi]->SetDirectory(0);
      hRawYieldBC1Ref[iPhi]->SetDirectory(0);
      hRawYieldBC2Ref[iPhi]->SetDirectory(0);
    }
    
    for(Int_t iPt=0; iPt<nPtBins; iPt++) {
      TGraphAsymmErrors* gtmp = (TGraphAsymmErrors*)reffile->Get(Form("gasigpt%d_%s",iPt,q2regionname.Data()));
      TGraphAsymmErrors* gtmpfs = (TGraphAsymmErrors*)reffile->Get(Form("gasigfspt%d_%s",iPt,q2regionname.Data()));
      TGraphAsymmErrors* gtmpBC1 = (TGraphAsymmErrors*)reffile->Get(Form("gasigBC1pt%d_%s",iPt,q2regionname.Data()));
      TGraphAsymmErrors* gtmpBC2 = (TGraphAsymmErrors*)reffile->Get(Form("gasigBC2pt%d_%s",iPt,q2regionname.Data()));
      
      if(gtmp && gtmpfs && gtmpBC1 && gtmpBC2) {
        for(Int_t iPhi=0; iPhi<nPhi; iPhi++) {
          Double_t phi,Y, errY;
          gtmp->GetPoint(iPhi,phi,Y);
          hRawYieldRef[iPhi]->SetBinContent(iPt+1,Y);
          gtmpfs->GetPoint(iPhi,phi,Y);
          hRawYieldfsRef[iPhi]->SetBinContent(iPt+1,Y);
          gtmpBC1->GetPoint(iPhi,phi,Y);
          hRawYieldBC1Ref[iPhi]->SetBinContent(iPt+1,Y);
          gtmpBC2->GetPoint(iPhi,phi,Y);
          hRawYieldBC2Ref[iPhi]->SetBinContent(iPt+1,Y);
        }
      }
    }
    reffile->Close();
  }
  else {return 1;}
  
  return 0;
}

//_____________________________________________________________________________________________
//GET EV PLANE RESOLUTION
Double_t GetEventPlaneResolution(Double_t& error, TH1F* hsubev1, TH1F* hsubev2, TH1F* hsubev3){
  Double_t resolFull=1.;
  if(evPlaneRes==AliEventPlaneResolutionHandler::kTwoRandSub ||
     evPlaneRes==AliEventPlaneResolutionHandler::kTwoChargeSub){
    resolFull=AliVertexingHFUtils::GetFullEvResol(hsubev1);
    error = TMath::Abs(resolFull-AliVertexingHFUtils::GetFullEvResolLowLim(hsubev1));
  }else if(evPlaneRes==AliEventPlaneResolutionHandler::kTwoEtaSub){
    if(evPlane==AliEventPlaneResolutionHandler::kTPCFullEta){
      resolFull=AliVertexingHFUtils::GetFullEvResol(hsubev1);
      error = TMath::Abs(resolFull-AliVertexingHFUtils::GetFullEvResolLowLim(hsubev1));
    }else if(evPlane==AliEventPlaneResolutionHandler::kTPCPosEta){
      resolFull=AliVertexingHFUtils::GetSubEvResol(hsubev1);
      error = TMath::Abs(resolFull-AliVertexingHFUtils::GetSubEvResolLowLim(hsubev1));
    }
  }else{
    Double_t resolSub[3];
    Double_t errors[3];
    TH1F* hevplresos[3] = {hsubev1,hsubev2,hsubev3};
    for(Int_t iRes=0;iRes<3;iRes++){
      resolSub[iRes]=hevplresos[iRes]->GetMean();
      errors[iRes]=hevplresos[iRes]->GetMeanError();
    }
    Double_t lowlim[3];for(Int_t ie=0;ie<3;ie++)lowlim[ie]=TMath::Abs(resolSub[ie]-errors[ie]);
    if(evPlane==AliEventPlaneResolutionHandler::kVZEROC ||
       evPlane==AliEventPlaneResolutionHandler::kVZERO){
      resolFull=TMath::Sqrt(resolSub[1]*resolSub[2]/resolSub[0]);
      error=resolFull-TMath::Sqrt(lowlim[2]*lowlim[1]/lowlim[0]);
    }
    else if(evPlane==AliEventPlaneResolutionHandler::kVZEROA){
      resolFull=TMath::Sqrt(resolSub[0]*resolSub[2]/resolSub[1]);
      error=resolFull-TMath::Sqrt(lowlim[2]*lowlim[0]/lowlim[1]);
    }
    else if(evPlane==AliEventPlaneResolutionHandler::kTPCFullEta ||
            evPlane==AliEventPlaneResolutionHandler::kTPCPosEta){
      resolFull=TMath::Sqrt(resolSub[0]*resolSub[1]/resolSub[2]);
      error=resolFull-TMath::Sqrt(lowlim[0]*lowlim[1]/lowlim[2]);
    }
  }
  return resolFull;
}

//_____________________________________________________________________________________________
//q2 CUT APPLICATION
void Applyq2Cut(THnSparseF* sparse, Double_t smallcutvalues, Double_t largecutvalues, Int_t smallorlarge) {
  
  UInt_t q2axnum=3;
  
  TAxis* q2axis = (TAxis*)sparse->GetAxis(q2axnum);
  Int_t binmin = -1;
  Int_t binmax = -1;
  //define q2 cut bins
  if(smallorlarge==kSmall) {
    binmin=q2axis->FindBin(0.0001);
    binmax=q2axis->FindBin(smallcutvalues-0.0001);
  }
  else if(smallorlarge==kLarge) {
    binmin=q2axis->FindBin(largecutvalues+0.0001);
    binmax=q2axis->FindBin(q2axis->GetNbins()+1); //takes also overflow entries
    Int_t binmaxsmallregion=q2axis->FindBin(smallcutvalues-0.0001);
    if(binmaxsmallregion==binmin) {
      binmin += 1;
    }
  }
  q2axis->SetRange(binmin,binmax);
}

//_____________________________________________________________________________________________
//CUT APPLICATION
void ApplyCut(THnSparseF* sparse, Double_t min, Double_t max, UInt_t axnum) {
  
  TAxis* axis = (TAxis*)sparse->GetAxis(axnum);
  Int_t binmin = axis->FindBin(min+0.0001);
  Int_t binmax = axis->FindBin(max-0.0001);
  axis->SetRange(binmin,binmax);
}

//_____________________________________________________________________________________________
//DEFINE q2 CUTS
Bool_t Defineq2Cuts(TH2F* hq2VsCentr, vector<Double_t> &smallcutvalues, vector<Double_t> &largecutvalues, Int_t cutmeth, Int_t fSparseVers) {
  
  Double_t cutvalues[2] = {-1.,-1.};
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    if(cutmeth==kAbsCut) {smallcutvalues.push_back(q2smalllimit); largecutvalues.push_back(q2largelimit);}
    else if(cutmeth==kPercCut || (cutmeth==kPercCutVsCent && fSparseVers==0)) {
      TH1F* hq2=(TH1F*)hq2VsCentr->ProjectionY();
      Bool_t cuttaken = Getq2CutValuePercEvents(hq2,cutvalues);
      if(cuttaken) {
        smallcutvalues.push_back(cutvalues[0]);
        largecutvalues.push_back(cutvalues[1]);
      }
      else {return kFALSE;}
    }
    else if(cutmeth==kPercCutVsCent && fSparseVers==1) {
      TH1F* hq2centbin=(TH1F*)hq2VsCentr->ProjectionY(Form("hq2centbin%d",iCent+1),iCent+1,iCent+1);
      Bool_t cuttaken = Getq2CutValuePercEvents(hq2centbin,cutvalues);
      if(cuttaken) {
        smallcutvalues.push_back(cutvalues[0]);
        largecutvalues.push_back(cutvalues[1]);
      }
      else {return kFALSE;}
    }
    else if(cutmeth==kPercentileCut) {
      smallcutvalues.push_back(q2smallpercevents*100);
      largecutvalues.push_back(100-q2largepercevents*100);
    }
  }
  
  TString q2name = "q2";
  if(cutmeth==kPercentileCut) {q2name = "q2 percentile";}
  cout << "Cut values for the small-q2 region:"<<endl;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    cout << "Centrality bin "<< iCent << Form(" %s < ",q2name.Data()) << smallcutvalues[iCent] <<endl;
  }
  cout << "\nCut values for the large-q2 region:"<<endl;
  for(Int_t iCent=0; iCent<ncentbins; iCent++) {
    cout << "Centrality bin "<< iCent << Form(" %s > ",q2name.Data()) << largecutvalues[iCent] <<endl;
  }
  cout << endl;
  
  return kTRUE;
}

//_____________________________________________________________________________________________
//q2 CUT VALUE SEARCH (ON A BASIS OF PERCENTAGE OF EVENTS WITH SMALLER/LARGER q2)
Bool_t Getq2CutValuePercEvents(TH1F* hq2, Double_t cutvalues[2]) {
  
  Int_t bincutvalues[2]={-1,-1};
  
  if(hq2) {
    Double_t SmallThresIntegral = q2smallpercevents*hq2->Integral(0.,1000.); //threshold integral value for small-q2 region (including overflow entries)
    Double_t LargeThresIntegral = q2largepercevents*hq2->Integral(0.,1000.); //threshold integral value for large-q2 region (including overflow entries)
    Double_t integral=0;
    Int_t q2bin=0;
    while(integral<SmallThresIntegral) {
      q2bin++;
      integral += hq2->GetBinContent(q2bin);
    }
    if(q2smallpercevents!=1.0) {bincutvalues[0]=q2bin-1;}//rounded down
    else {bincutvalues[0]=q2bin;}//exactly 100%
    if(q2smallpercevents==0.5 && q2largepercevents==0.5) bincutvalues[1]=bincutvalues[0]+1;
    else {
      integral=0;
      q2bin=hq2->GetNbinsX()+2;
      while(integral<LargeThresIntegral) {
        q2bin--;
        integral += hq2->GetBinContent(q2bin); //takes also overflow counts
      }
      if(q2largepercevents!=1.0) {bincutvalues[1]=q2bin+1;}//rounded down
      else {bincutvalues[1]=q2bin;}//exactly 100%
    }
    
    cutvalues[0]=hq2->GetBinLowEdge(bincutvalues[0])+hq2->GetBinWidth(bincutvalues[0]);
    cutvalues[1]=hq2->GetBinLowEdge(bincutvalues[1]);
    
    return kTRUE;
  }
  
  cout << "Warning: hq2 not set, impossible to define the cut values!" << endl;
  return kFALSE;
}

//_____________________________________________________________________________________________
//RESET AXES OF SPARSE
void ResetAxes(THnSparseF* sparse, Int_t axnum) {
  if(axnum<0) {
    for(Int_t iAxis=0; iAxis<sparse->GetNdimensions(); iAxis++) {
      sparse->GetAxis(iAxis)->SetRange(-1,-1);
    }
  }
  else {sparse->GetAxis(axnum)->SetRange(-1,-1);}
}

//_____________________________________________________________________________________________
void GetMinMaxHisto(TH1F* histo, Double_t &xmin, Double_t &xmax) {
  
  Int_t nbins=histo->GetNbinsX();
  xmin=histo->GetBinLowEdge(nbins)+histo->GetBinWidth(nbins);
  Int_t content=0;
  for(Int_t iBin=0; iBin<nbins; iBin++) {
    content=histo->GetBinContent(iBin+1);
    if(content>0) {
      xmin=histo->GetBinLowEdge(iBin+1);
      break;
    }
  }
  content=0;
  xmax=histo->GetBinLowEdge(1);
  for(Int_t iBin=nbins; iBin>0; iBin--) {
    content=histo->GetBinContent(iBin);
    if(content>0) {
      xmax=histo->GetBinLowEdge(iBin)+histo->GetBinWidth(iBin);
      break;
    }
  }
}

//_____________________________________________________________________________________________
void DivideCanvas(TCanvas* c, const Int_t nPtBins) {
  if(nPtBins<2)
    c->Divide(1,1);
  if(nPtBins==2 || nPtBins==3)
    c->Divide(nPtBins,1);
  else if(nPtBins==4 || nPtBins==6 || nPtBins==8)
    c->Divide(nPtBins/2,2);
  else if(nPtBins==5 || nPtBins==7)
    c->Divide((nPtBins+1)/2,2);
  else if(nPtBins==9 || nPtBins==12 || nPtBins==15)
    c->Divide(nPtBins/3,3);
  else if(nPtBins==10 || nPtBins==11)
    c->Divide(4,3);
  else if(nPtBins==13 || nPtBins==14)
    c->Divide(5,3);
  else if(nPtBins>15 && nPtBins<=20 && nPtBins%4==0)
    c->Divide(nPtBins/4,4);
  else if(nPtBins>15 && nPtBins<=20 && nPtBins%4!=0)
    c->Divide(5,4);
  else if(nPtBins==21)
    c->Divide(7,3);
  else if(nPtBins>21 && nPtBins<=25)
    c->Divide(5,5);
  else if(nPtBins>25 && nPtBins%2==0)
    c->Divide(nPtBins/2,2);
  else
    c->Divide((nPtBins+1)/2,2);
}

//_____________________________________________________________________________________________
//DRAWING STYLE
void SetStyle() {
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetTitleFillColor(0);
  gStyle->SetStatColor(0);
  //gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetTextFont(42);
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.3,"y");
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetStatFont(42);
  gStyle->SetStatY(0.89);
  gStyle->SetStatX(0.89);
  gStyle->SetTitleFont(42,"xyzg");
  gStyle->SetHistLineWidth(2);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPalette(53);
}

//_____________________________________________________________________________________________
//LOAD REFLECTION HISTOS FOR D0
Bool_t LoadD0toKpiMCHistos(TList *outlist){
  
  TFile *f=new TFile(fileNameMCD0refl.Data(),"READ");
  if(!f){
    printf("ERROR: file %s does not exist\n",fileNameMCD0refl.Data());
    return kFALSE;
  }
  f->ls();
  TH1F** hsig=new TH1F*[nPtBins];
  TH1F** hrfl=new TH1F*[nPtBins];
  for(Int_t j=0;j<nPtBins;j++){
    hsig[j]=(TH1F*)f->Get(Form("histSgn_%d",j));
    if(!hsig[j]) {Printf("histSgn_%d NOT FOUND",j); return kFALSE;}
    hrfl[j]=(TH1F*)f->Get(Form("histRflFitted%s_ptBin%d",rflFitType.Data(),j));
    if(!hrfl[j]) {Printf("histRflFitted%s_ptBin%d",rflFitType.Data(),j); return kFALSE;}
  }
  for(Int_t k=0;k<nPtBins;k++){
    outlist->Add(hsig[k]->Clone(Form("histSgn_%d",k)));
    outlist->Add(hrfl[k]->Clone(Form("histRfl_%d",k)));
  }
  outlist->ls();
  //for(Int_t p=0;p<nptbinsnew;p++){
  //delete hsig[p];
  //delete hrfl[p];
  //}
  //f->Close();
  return kTRUE;
}
