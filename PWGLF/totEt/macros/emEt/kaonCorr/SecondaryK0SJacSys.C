#include <iostream>
#include <fstream>
#include "TMath.h"

using namespace std;

char myjunk[200];

void SetStyles(TH1 *histo, char *ytitle, char *xtitle, int color, int marker);
void SetCanvasStyle(TCanvas *c);
TLegend *GetLegend(float x1, float y1, float x2, float y2);
Bool_t doLevy = kFALSE;//for debugging
Bool_t doBlast = kTRUE;//for debugging
void Draw(TH1 *data, TF1 *fLevy, TF1 *fBlast,int centbin,char *pid,char *name);
void PrintLatex(int etCutNum, char *det);
void PrintArrays();

// Centrality dependent factors
int nbins = 10;  // number of centrality bins
Double_t etk0[2][10] = {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
Double_t etBlastk0[2][10] = {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
Double_t etBlastk02D[2][10] = {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
//et cuts for kaon energies
//0.00, 0.050, 0.100, 0.150, 0.200
//0.250,0.300, 0.350, 0.400, 0.450
//0.500
//counting is integer and starts are 1

Double_t etCutOffs[11] = {0.00, 0.050, 0.100, 0.150, 0.200,  0.250,0.300, 0.350, 0.400, 0.450,  0.500};

Double_t kaonDeposits[2][10] = {{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};

//note:  Not at all sensitive to which Jacobian is used, which also means the pT used doesn't really matter
void SecondaryK0SJacSys(char *inputfilename = "../workingdir/rootFiles/LHC11a10a_bis/Et.ESD.simPbPb.EMCal.LHC11a10a_bis.Run139465.root", char *det = "EMCal", int whichjacobian = 1, int etCutNum = 7, Bool_t longRun = kFALSE){//Kaon collection: 0=K0S, -1=K-, +1=K+
  Int_t kaonSelection = 1;
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetLineWidth(3);
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gROOT->ProcessLine(".L AliBWFunc.cxx+g");
  gROOT->ProcessLine(".L AliBWTools.cxx+g");
  gSystem->Load("libTree.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWGLFspectra.so");
  gSystem->Load("libPWGTools.so");
  gROOT->LoadMacro("macros/FitParticle.C+");
  gROOT->LoadMacro("macros/GetLevidEtdptTimesPt.C");
  //inputting kaon errors from http://arxiv.org/abs/1303.0737
  //using average of K+ and K-
  Float_t kaonYield[10] = {109,90.5,68,46,30,18.2,10.2,5.1,2.3,0.855};
  Float_t kaonStatErr[10] = {0.3,0.2,0.1,0.1,0.1, 0.06,0.04,0.03,0.02,0.01};
  Float_t kaonSysErr[10] = {9,7,5,4,2, 1.5,0.8,0.4,0.2,0.09};
  Float_t kaonTotErr[10] = {0,0,0,0,0, 0,0,0,0,0};
  Float_t kaonFracErr[10] = {0,0,0,0,0, 0,0,0,0,0};
  for(int i=0;i<10;i++){
    kaonTotErr[i] = TMath::Sqrt(kaonSysErr[i]*kaonSysErr[i]+kaonStatErr[i]*kaonStatErr[i]);
    kaonFracErr[i] = kaonTotErr[i]/kaonYield[i];
    //cout<<"bin "<<i<<" err "<<kaonFracErr[i]<<endl;
  }

  //centrality bin #  (want to loop over 5 bins)
  int centbin = 0;  // can change later
  TString detector = det;
  float etarange = 0.7;
  if(detector.Contains("EMC")){
    etarange = 0.7;
  }
  else{etarange = 0.1;}

  //==================READ IN DATA==========================

  TFile *chargedkaonfile = new TFile("rootFiles/SPECTRA_COMB_20120709.root");
  gROOT->LoadMacro("macros/k0sFinal.C");                // load data spectra
  TClonesArray histosk0 = GetK0S();                     // get K0S histogram
  
  TCanvas *cChK = new TCanvas("cChK","ChargedKaons",700,500);
  //cChK->Divide(4,3);
  TH1D *histoskch[] = {NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL, NULL,NULL};
  if(kaonSelection==1){
    for(int i=0;i<=9;i++){
      TH1D *histoTemp = (TH1D*) chargedkaonfile->Get(Form("cent%i_kaon_plus",i));
      histoskch[i] = (TH1D*)histoTemp->Clone("cent%i_kaon_plus_clone");
      if(i==0){
	((TH1D*)histoskch[i])->Draw();
      }
      else{
	((TH1D*)histoskch[i])->Draw("same");
      }
      //}
    }
  }
  else{
    for(int i=0;i<=9;i++){
      TH1D *histoTemp = (TH1D*) chargedkaonfile->Get(Form("cent%i_kaon_minus",i));
      histoskch[i] = (TH1D*)histoTemp->Clone("cent%i_kaon_minus_clone");
      if(i==0){
	((TH1D*)histoskch[i])->Draw();
      }
      else{
	((TH1D*)histoskch[i])->Draw("same");
      }
      //}
    }
  }
  TCanvas *cEmpty = new TCanvas("cEmpty","Empty",700,500);

  TString particleName = "K0";
  TString particleLatex = "K^{0}_S";
  if(kaonSelection==1){
    particleName = "K+";
    particleLatex = "K^{+}";
  }
  if(kaonSelection==-1){
    particleName = "K-";
    particleLatex = "K^{-}";
  }

  //==================FIT DATA==========================
  TF1 *funck0[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};          
  TF1 *funcBlastk0[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};     // Blastwave fit to K0 pT?
  TF1 *funcetk0[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};         
  TF1 *funcetBlastk0[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};   // Blastwave fit to K0 et
  TF1 *funcetBlastk02D[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};
  TF1 *funcBlastk0Jac[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL}; 
 
  TH1 *histoSysUp[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};  //K0S Clone histogram
  TF1 *funcBlastk0Up[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL}; 
  TH1 *histoSysLow[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};  //K0S Clone histogram
  TF1 *funcBlastk0Low[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL}; 
  
  TH1 *histoSys[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};  //K0S Clone histogram
  
  
  Float_t TInit = 0.1326;    // initialize T parameter of Blastwave fit
  Float_t nInit = 6.6857;    // initialize n parameter of Blastwave fit
  Float_t normInit = 20;     // set normalization parameter of Blastwave fit
  AliBWTools *tools = new AliBWTools();
  Double_t yield;            // yield
  Double_t yieldError;       // yield error
  Double_t partialYields[] = {0,0,0};
  Double_t partialYieldErrors[] = {0,0,0};
  Double_t Ameson = 0;
  Double_t Abaryon = -1;
  Double_t Aantibaryon = 1;
  
  //Systematic erros
  
  TF1 *funcBlastk0Par1[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};
  TF1 *funcBlastk0Par2[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};
  Double_t Par[10] = {0,0,0,0,0, 0,0,0,0,0};
  Double_t ErrPar[10] = {0,0,0,0,0, 0,0,0,0,0};
 
  // Genarate y
   
  TRandom* randy = new TRandom();  
  
  // integral of fit to data
  float fitINTEGRALdata[10] = {NULL,NULL,NULL,NULL,NULL, NULL,NULL,NULL,NULL,NULL};

  if(kaonSelection!=0){//charged kaons
    nbins = 10;
  }

  // loop over centrality bins 0->5 and genarate histos shift up and down to Systematic errors
  for(int i=0;i<nbins;i++){
    //for(int i=centbin;i<centbin+1;i++){
    if(kaonSelection==0){
      histoSys[i] = (TH1F*)histosk0[i]->Clone(Form("Sys%i",i));
      histoSysLow[i] = (TH1F*)histosk0[i]->Clone(Form("SysLow%i",i));
      histoSysUp[i] = (TH1F*)histosk0[i]->Clone(Form("SysUp%i",i));
    }
    else{
      histoSys[i] = (TH1D*)histoskch[i]->Clone(Form("Sys%i",i));
      histoSysLow[i] = (TH1D*)histoskch[i]->Clone(Form("SysLow%i",i));
      histoSysUp[i] = (TH1D*)histoskch[i]->Clone(Form("SysUp%i",i));
    }
    float NbinsSys =  histoSysUp[i]->GetXaxis()->GetNbins();
    
    for(int x = 1; x<=NbinsSys; x++){
		
      float BinContent = histoSysUp[i]->GetBinContent(x);
      float BinError =  histoSysUp[i]->GetBinError(x);
      histoSysUp[i]->SetBinContent(x, BinContent+BinError);
      histoSysLow[i]->SetBinContent(x, BinContent-BinError);
    }

    funcBlastk0Up[i] = FitParticle((TH1*)histoSysUp[i],particleName.Data(),-1,-1,-1,0,0,3,-1,kFitBlastWave);
    funcBlastk0Low[i] = FitParticle((TH1*)histoSysLow[i],particleName.Data(),-1,-1,-1,0,0,3,-1,kFitBlastWave);
    funcBlastk0[i] = FitParticle((TH1*)histoSys[i],particleName.Data(),-1,-1,-1,0,0,3,-1,kFitBlastWave);

    if(kaonSelection==0){
      funck0[i] = FitParticle((TH1F*)histosk0[i],particleName.Data(),TInit,normInit,nInit);
      if(!histosk0[i]) cerr<<"histosk0[i] does not exist!"<<endl;
      if(!funck0[i]) cerr<<"funck0[i] does not exist!"<<endl;
      Draw((TH1F*)histosk0[i],funck0[i],funcBlastk0[i],i,particleLatex.Data(),particleName.Data());
    }
    else{
      if(!histoskch[i]) cerr<<"Warning!  Histogram does not exist!"<<endl;
      funck0[i] = FitParticle((TH1D*)histoskch[i],particleName.Data());
      cerr<<"174"<<endl;
      Draw(histoskch[i],funck0[i],funcBlastk0[i],i,particleLatex.Data(),particleName.Data());
    }
  }

 
  //==================THROWING RANDOM SPECTRA========================== 
  // number of particles to put in histogram (will change later)
  int nParticles = 1e3; // for quick runs
  if(longRun) nParticles = 1e7;   //for statistics

  // set parameters of setstyles
  int marker1 = 21; int marker2 = 24;
  int color1 = 2; int color2 = 4;
  int mixmarker = 22; int mixcolor = 3;
  
  TFile *f = TFile::Open(inputfilename, "READ");
  TList *l = (TList*)f->Get("out1");
  TH3F *fHistK0EDepositsVsPtInAcceptance = l->FindObject("fHistK0EDepositsVsPtInAcceptance");
  TH3F *fHistK0EDepositsVsPtOutOfAcceptance = l->FindObject("fHistK0EDepositsVsPtOutOfAcceptance");
  fHistK0EDepositsVsPtInAcceptance->GetZaxis()->SetRange(etCutNum,etCutNum);
  TH2D *my2DhistoK0S = (TH2D*) fHistK0EDepositsVsPtInAcceptance->Project3D("yx");
  fHistK0EDepositsVsPtOutOfAcceptance->GetZaxis()->SetRange(etCutNum,etCutNum);
  TH2D *my2DhistoK0SOutOfAcceptance = (TH2D*) fHistK0EDepositsVsPtOutOfAcceptance->Project3D("yx");
  TH1F *fHistSimKaonsInAcceptance = l->FindObject("fHistSimKaonsInAcceptance");
  TH1F *fHistSimKaonsOutOfAcceptance = l->FindObject("fHistSimKaonsOutOfAcceptance");
  TH1F *hRatioOutOfAccOverInAcc = fHistSimKaonsOutOfAcceptance->Clone("hRatioOutOfAccOverInAcc");
  hRatioOutOfAccOverInAcc->Divide(fHistSimKaonsInAcceptance);

  //find number of bins of Caio's histogram
  int bincount = 0;
  bincount = my2DhistoK0S->GetXaxis()->GetNbins();
  cout<<"Number of bins in histo: "<<bincount<<endl;
  
  // create array to hold histograms
  TClonesArray histoarrayK0S("TH1D",bincount);  //(nbins);
  TClonesArray histoarrayK0SOutOfAcceptance("TH1D",bincount);  //(nbins);

  //==================PROJECTIONS FOR INPUT ET DEPOSITS==========================
  //loop over bins in x to define projections  (nbins)
  for(int i=0; i<bincount; i++){

    // make projections of 2D histogram my2DhistoK0S
    histoarrayK0S[i] = (TH1D*)my2DhistoK0S->ProjectionY(Form("tmpK0S%i", i+1), i+1, i+1);
    histoarrayK0SOutOfAcceptance[i] = (TH1D*)my2DhistoK0SOutOfAcceptance->ProjectionY(Form("tmpK0SOutOfAcc%i", i+1), i+1, i+1);
  }


  for(int j=0;j<nbins;j++){
    centbin = j;  
    cout<<"Working on centrality bin "<<centbin<<endl;


    // creates histogram of spectra from thrown K0S
    TH1F *histoSpectrumK0S = new TH1F("histoSpectrumK0S","p_{T} spectrum of K^{0}_{S}",100,0.0,5.0);
    SetStyles(histoSpectrumK0S,"dN/dp_{T}","p_{T} (GeV)", color1, marker1);
    TH1F *histoETspectrumK0S = new TH1F("histoETspectrumK0S","transverse energy scaled of K^{0}_{S}",100,0.0,5.0);
    SetStyles(histoETspectrumK0S,"dN/dE_{T} calculated","E_{T} (GeV)", color2, marker2);
    TH1F *histoETCspectrumK0S = new TH1F("histoETCspectrumK0S","transverse energy scaled from data of K^{0}_{S}",100,0.0,5.0);
    SetStyles(histoETCspectrumK0S,"dN/dE_{T} from data","E_{T} (GeV)", mixcolor, mixmarker);

    // initialize totals
    float totET = 0;               // total ET of K0S's calculated
    float totetcK0S = 0;           // total ET of K0S's from data
    float testTOT = 0;
    float testTOTpt = 0;
    float pK0S = 0;  

    // declare edges of first and last bin
    int mybin100 = my2DhistoK0S->GetXaxis()->FindBin(0.1001);
    int mybin0 = my2DhistoK0S->GetXaxis()->FindBin(0.0001);
  
    //Jacobian correction
  
    TH1F *histoSpectrumK0SJac;
    histoSpectrumK0SJac = (TH1F*)histoSys[centbin]->Clone("Jac");
    SetStyles(histoSpectrumK0SJac,"dN/dp_{T}","p_{T} (GeV)", color1, marker1);

    TH1F *histoSpectrumK0SDivide = (TH1F*)histoSys[centbin]->Clone("JacDiv");

  
//     TH1F *histoeta = new TH1F("histoeta","",100,-1,1);
//     SetStyles(histoeta,"N","#eta", color1, marker1);
  
    //Jacobian Correction 
  
    float NbinsJac = histoSys[centbin]->GetXaxis()->GetNbins(); 
    double mK0SJac = 0.497614;            // mass of K0S (GeV)  
    double nK0s = 1000;
    double deltay =0;
    double deltaeta = 0;
  
  
    TF1 *jacobian = new TF1("jacobian","[0]/TMath::Sqrt([1]*[1]*TMath::Power(TMath::CosH(x),2)+[0]*[0])",-etarange,etarange);
    jacobian->SetParameter(1,0.493);
    //here we're going to come up with a simple model that looks at a non-flat eta dependence.
    //if I assume that the pseudorapidity distribution has the shape y = m|x|+b, I have to normalize by the integral of that function after I take the jacobian so I can get the right weighting.  That integral works out to be
    float m = +0.15;
    float b = 1;
    float renormalization = (m*etarange*etarange)/2+b*etarange;
    TF1 *jacobianReweighted = new TF1("jacobianReweighted","([2]*x+[3])*[4]*[0]/TMath::Sqrt([1]*[1]*TMath::Power(TMath::CosH(x),2)+[0]*[0])",-etarange,etarange);
    jacobianReweighted->SetParameter(1,0.493);
    jacobianReweighted->SetParameter(2,m);
    jacobianReweighted->SetParameter(3,b);
    jacobianReweighted->SetParameter(4,renormalization);
    //==================CAIO'S TOY MODEL FOR JACOBIAN==========================
    float overallAverageJac[5] = {0,0,0,0,0};
    int totK0S = 0;
    for(int a = 1; a<NbinsJac; a ++){
		  
      int count = 0;
      double binContent = histoSys[centbin]->GetBinContent(a);
      double binError = histoSys[centbin]->GetBinError(a);
      jacobian->SetParameter(0,histoSys[centbin]->GetXaxis()->GetBinCenter(a));
      float averageJac = jacobian->Integral(-etarange,etarange) / (2*etarange);
      jacobian->SetParameter(0,histoSys[centbin]->GetXaxis()->GetBinLowEdge(a));
      float lowJac = jacobian->Integral(-etarange,etarange) / (2*etarange);
      jacobian->SetParameter(0,histoSys[centbin]->GetXaxis()->GetBinLowEdge(a+1));
      float highJac = jacobian->Integral(-etarange,etarange) / (2*etarange);
      float eta0Jac = jacobian->Eval(0.0);
      jacobianReweighted->SetParameter(0,histoSys[centbin]->GetXaxis()->GetBinCenter(a));
      float altJac = jacobianReweighted->Integral(-etarange,etarange);
      double correc =1;
		
      double mult = (binContent*correc)*deltay*deltaeta;
      float scale = averageJac;
      switch(whichjacobian){
      case 1:
	scale = averageJac;
      case 2:
	scale = lowJac;
      case 3:
	scale = highJac;
      case 4:
	scale = eta0Jac;
      case 5:
	scale = altJac;
      }
      histoSpectrumK0SJac->SetBinContent(a,scale*binContent);
      histoSpectrumK0SJac->SetBinError(a,scale*binError);
      totK0S += binContent;
      overallAverageJac[0] += averageJac*binContent;
      overallAverageJac[1] += lowJac*binContent;
      overallAverageJac[2] += highJac*binContent;
      overallAverageJac[3] += eta0Jac*binContent;
      overallAverageJac[4] += altJac*binContent;
		
    }
    cEmpty->cd();
    histoSys[centbin]->Draw();
    histoSpectrumK0SJac->Draw("same");
    cEmpty->SaveAs(Form("/tmp/KaonCut%i%s.png",etCutNum,det));
    float minJacobian = overallAverageJac[0]/totK0S;
    float maxJacobian = overallAverageJac[0]/totK0S;
    cout<<"Average jacobians: ";
    for(int i=0;i<=4;i++){
      cout<<" "<< overallAverageJac[i]/totK0S;
      if(minJacobian>overallAverageJac[i]/totK0S) minJacobian = overallAverageJac[i]/totK0S;
      if(maxJacobian<overallAverageJac[i]/totK0S) maxJacobian = overallAverageJac[i]/totK0S;
    }
    cout<<endl;
    float errJacobian = (maxJacobian-minJacobian)/2.0;
    float meanJacobian = (maxJacobian+minJacobian)/2.0;
    cout<<"Jacobian with error "<<meanJacobian<<" +/- "<<errJacobian<<endl;
    cEmpty->cd();
    funcBlastk0Jac[centbin] = FitParticle(histoSpectrumK0SJac,"K0",-1,-1,-1,0,0,3,-1,kFitBlastWave);
    funcBlastk0[centbin]->Draw("same");	  
	
    // gets integral in input K0S histograms for different centrality bins   
    fitINTEGRALdata[centbin] = histoSpectrumK0SJac->Integral("width");

    cout<<"integral of fit for centrality bin "<<centbin<<": "<<fitINTEGRALdata[centbin]<<endl;
    cout<<funcBlastk0Jac[centbin]->Integral(0.0,100.0)<<endl;
    
    //==================END CAIO'S TOY MODEL FOR JACOBIAN==========================
     
    //==================THROWING RANDOM K0S==========================
    float totetcK0SOutOfAcc = 0;
    float totetcK0SOutOfAccCorr = 0;
    float avgCorrOutOfAcc = 0;
  
    //THROW random K0S particle 
    for(int i=0; i<nParticles; i++){
      float mK0S = 0.497614;            // mass of K0S (GeV)
      //You can change funBlastk0 to funBlastk0up or funBlastk0low for systematica error analysis 
      float ptK0S = funcBlastk0Jac[centbin]->GetRandom(); // array? centbin =1
      histoSpectrumK0S->Fill(ptK0S);     
      
      // *******************************************
      float etK0S = sqrt(ptK0S*ptK0S + mK0S*mK0S);          // calculate transverse energy of K0S
      histoETspectrumK0S->Fill(etK0S);
      totET += etK0S;

      int mybinK0S = my2DhistoK0S->GetXaxis()->FindBin(ptK0S);  
      int xbins = my2DhistoK0S->GetNbinsX();

      float scaleForOutOfAcc = hRatioOutOfAccOverInAcc->GetBinContent(hRatioOutOfAccOverInAcc->FindBin(ptK0S));
      //cout<<"Scale for out of acc "<<scaleForOutOfAcc<<endl;
      
      float etcK0S = 0;
      float etcK0SOutOfAcc = 0;
      // get ET projection of K0S and sum
      if(mybinK0S>0 && mybinK0S<=bincount){  //nbins  // then this is in our range
	etcK0S = ((TH1D*)histoarrayK0S.At(mybinK0S-1))->GetRandom();
	etcK0SOutOfAcc = ((TH1D*)histoarrayK0SOutOfAcceptance.At(mybinK0S-1))->GetRandom();
	histoETCspectrumK0S->Fill(etcK0S);
	//here we deal with a problem.  The get random function gets confused for zeros.  It gives a small but non-zero value.  But we can never have a deposit less than the minimum energy!
	//this doesn't mess up the kaons in the 
	if(etK0S<etCutOffs[etCutNum-1]){etK0S = 0.0;}
	if(etcK0SOutOfAcc<etCutOffs[etCutNum-1]){etcK0SOutOfAcc = 0.0;}

	//if(etcK0SOutOfAcc>0){cout<<"et dep "<<etcK0SOutOfAcc<<" scale "<<scaleForOutOfAcc<<" pT "<<ptK0S<<endl;}
	//else{cout<<"No energy deposited"<<endl;}
	//else{cout<<"deposit is not small "<<etcK0SOutOfAcc<<endl;}
      }
      else{cerr<<"Did not find pt bin!"<<endl;}    // this should never be reached
      totetcK0S += etcK0S;
      totetcK0SOutOfAccCorr += etcK0SOutOfAcc*scaleForOutOfAcc; 
      //if(etcK0SOutOfAcc>1e-3) cout<<"Deposit is not small "<<etcK0SOutOfAcc<<" corr "<<scaleForOutOfAcc<<" product "<<totetcK0SOutOfAccCorr<<" total "<<endl;
      totetcK0SOutOfAcc +=etcK0SOutOfAcc;
      avgCorrOutOfAcc += scaleForOutOfAcc;
    }         		  

    // print Calculated and Data projection values of scaled Et for K0S
    cout<<"Total ET observed (before renormalization): "<<(totetcK0S)<<endl;//Total ET observed
    cout<<"Total possible  energy (before renormalization):  "<<totET<<endl;
    cout<<"Average ET observed (simple): "<<(totetcK0S)/nParticles<<endl;
    //cout<<"Average ET observed out of acceptance (simple): "<<(totetcK0SOutOfAcc)/nParticles<<endl;
    //cout<<"Average corr for out of acc (simple): "<<(avgCorrOutOfAcc)/nParticles<<endl;
    //cout<<"ET observed out of acceptance corr (simple): "<<(totetcK0SOutOfAccCorr)<<endl;
    cout<<"Average ET observed out of acceptance corr (simple): "<<(totetcK0SOutOfAccCorr)/nParticles<<endl;

    //Estimating the difference between the number of particles in |eta|<0.5 and our eta range
    //assume the distribution of particles is approximately N=m*|eta|+b
    //if N(eta=0) = 1, b=1
    //The percentage drop/increase in eta over one unit is m
    //The integral from 0-0.5 is
    m = 0.15;
    b = 1;
    float etaspectrameas = 0.5;
    float integralSpectraMeas = m*etaspectrameas*etaspectrameas/2+b*etaspectrameas;
    //and the integral over our range is:
    float integralMeasHigh = m*etarange*etarange/2+b*etarange;
    m = -0.15;
    float integralMeasLow = m*etarange*etarange/2+b*etarange;
    float etaRangeWeightMean = (integralMeasHigh/integralSpectraMeas+integralMeasLow/integralSpectraMeas)/2.0;
    float etaRangeWeightErr = (integralMeasHigh/integralSpectraMeas-integralMeasLow/integralSpectraMeas)/2.0;
    //cout<<" range high "<< integralMeasHigh/integralSpectraMeas<< " range low "<< integralMeasLow/integralSpectraMeas<<endl;
  
    //cout<<"Eta range scale: "<<etaRangeWeightMean<<" +/- "<<etaRangeWeightErr<<endl;

    //getting the final ET:
    //multiply by the dN/dy from the spectra paper (with error)
    //scale by Jacobian with the error on the Jacobian (with error)
    //multiply by the 4 kaon species
    //add an error for the fact that we need to extrapolate to our eta range
    float scale = 4*kaonYield[centbin]*meanJacobian*etaRangeWeightMean;
    float meanETperkaon = (totetcK0S+totetcK0SOutOfAccCorr)/nParticles;
    float meanETperkaonInAcc = (totetcK0S)/nParticles;
    float meanETperkaonOutOfAcc = (totetcK0SOutOfAccCorr)/nParticles;
    //now we're going to allow the ratio of kaons in and out of acceptance to vary by another 20%
    float fracerrFromInOutVariance = 0.2*meanETperkaonOutOfAcc/(meanETperkaonInAcc+meanETperkaonOutOfAcc);
    float fracerr = TMath::Sqrt(TMath::Power(kaonFracErr[centbin],2)+TMath::Power(errJacobian/meanJacobian,2)+TMath::Power(etaRangeWeightErr/etaRangeWeightMean,2));
    float fracerrTotal = TMath::Sqrt(TMath::Power(kaonFracErr[centbin],2)+TMath::Power(errJacobian/meanJacobian,2)+TMath::Power(etaRangeWeightErr/etaRangeWeightMean,2)+TMath::Power(fracerrFromInOutVariance,2));
    cout<<"Total ET deposited in acc: "<<meanETperkaonInAcc*scale<<" +/- "<<meanETperkaonInAcc*fracerr*scale<<endl;
    cout<<"Total ET deposited out of acc: "<<meanETperkaonOutOfAcc*scale<<" +/- "<<meanETperkaonOutOfAcc*fracerr*scale<<endl;
    cout<<"Total ET deposited: "<<meanETperkaon*scale<<" +/- "<<meanETperkaon*fracerrTotal*scale<<endl;

    kaonDeposits[0][centbin] = meanETperkaon*scale;
    kaonDeposits[1][centbin] = meanETperkaon*scale*fracerrTotal;

    delete histoSpectrumK0S;
    delete histoETspectrumK0S;
    delete histoETCspectrumK0S;
    delete histoSpectrumK0SJac;
    delete histoSpectrumK0SDivide;
    delete jacobian;
    delete jacobianReweighted;
  }


  PrintLatex(etCutNum,det); 
}


// function to define canvas characteristics globally
void SetCanvasStyle(TCanvas *c){
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  c->SetTopMargin(0.0254237);    // 0.04
  c->SetRightMargin(0.0322581);  // 0.04
  c->SetLeftMargin(0.129032);    // 0.181452
  c->SetBottomMargin(0.134409);  
}

// function to define characteristics globally of plots
void SetStyles(TH1 *histo, char *ytitle, char *xtitle, int color, int marker){
  histo->Sumw2();
  histo->GetYaxis()->SetTitle(ytitle);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->SetMarkerStyle(marker);
}


TLegend *GetLegend(float x1, float y1, float x2, float y2){
  TLegend *leg = new TLegend(x1,y1,x2,y2);
  leg->SetBorderSize(0);
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.0381356);
  return leg;
}

void Draw(TH1 *data, TF1 *fLevy, TF1 *fBlast, int centbin,char *pid,char *name){
  TCanvas *canvas = new TCanvas("canvas","canvas",500,500);
  SetCanvasStyle(canvas);
  //canvas->SetLogy();
  data->SetTitle("p_{T} spectra of K$^{0}_{S}$ centbin=%i");
  data->GetYaxis()->SetTitle("1/2#pi p_{T} d^2N/dp_{T}dy");
  data->GetXaxis()->SetTitle("p_{T} (GeV)");
  data->GetXaxis()->SetTitleSize(0.05);
  data->GetYaxis()->SetTitleSize(0.05);
  data->GetXaxis()->SetRange(1,data->FindBin(3.5));
  //data->SetMarkerStyle(20);
  data->Draw();
  fLevy->Draw("same");
  fBlast->Draw("same");
  fLevy->SetLineColor(2);
  fBlast->SetLineColor(4);
  TLegend *legend = GetLegend(0.118952,0.120763,0.368952,0.271186);
  legend->AddEntry(data,Form("%s centrality bin %i",pid,centbin),"l");
  legend->AddEntry(fLevy,"Levy Fit");
  legend->AddEntry(fBlast,"BlastWave Fit");
  legend->Draw();
  canvas->SaveAs(Form("pics/%s%i.png",name,centbin));
  delete canvas;
}

void PrintLatex(int etCutNum, char *det){
  ofstream myfile;
  myfile.open(Form("datatablesKaonCut%i%s.tex",etCutNum,det));

  myfile<<Form("%3.2f",etCutOffs[etCutNum]);
  myfile<<"& ";
  //myfile<<"K$^{0}_{S}$  ";
  for(int j=0;j<nbins;j++){
    myfile<<Form("%3.1f $\\pm$ %3.1f",kaonDeposits[0][j],kaonDeposits[1][j]);
    if(j!=nbins-1) myfile<<"& ";
  }
  //myfile<<" & ";
  myfile<<"\\\\%data"<<endl;
  myfile<<endl<<endl;

  myfile.close();
}

void PrintArrays(){
  ofstream myfile;
  myfile.open("arraysv0.dat");

  myfile<<"etk0[2][10] = {";
  for(int i=0;i<2;i++){
    myfile<<"{";
    for(int j=0;j<nbins;j++){
      myfile<<etk0[i][j];
      if(j<nbins-1) myfile<<",";
    }
    myfile<<"}";
    if(i<1) myfile<<",";
  }
  myfile<<"};"<<endl;
  myfile<<endl<<endl; 

  myfile<<"etBlastk0[2][10] = {";
  for(int i=0;i<2;i++){
    myfile<<"{";
    for(int j=0;j<nbins;j++){
      myfile<<etBlastk0[i][j];
      if(j<nbins-1) myfile<<",";
    }
    myfile<<"}";
    if(i<1) myfile<<",";
  }
  myfile<<"};"<<endl;
  myfile<<endl<<endl; 

  myfile<<"etBlastk02D[2][10] = {";
  for(int i=0;i<2;i++){
    myfile<<"{";
    for(int j=0;j<nbins;j++){
      myfile<<etBlastk02D[i][j];
      if(j<nbins-1) myfile<<",";
    }
    myfile<<"}";
    if(i<1) myfile<<",";
  }
  myfile<<"};"<<endl;

  myfile.close();
}



 
	  
	   

