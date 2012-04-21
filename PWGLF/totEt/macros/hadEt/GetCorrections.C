//Create by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//This macro takes an input file created by AliAnalysisHadEtMonteCarlo and creates an AliAnalysisHadEtCorrections for the determination of the corrected Et
// #include "TFile.h"
// #include <iostream>
// #include "AliAnalysisHadEtCorrections.h"

// #include <iostream>
// #include <TROOT.h> 
// #include <TSystem.h>
// #include "TStopwatch.h"

//Corrections added in by hand to deal with the inadequacies of PYTHIA and HIJING
TH1D *pp276TPCBkgd();
TH1D *pp276ITSBkgd();
TH1D *PbPb276TPCBkgd();
TH1D *PbPb276ITSBkgd();
   Double_t xAxis1[112] = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 12, 14, 16, 18, 20, 25}; 
   

Float_t CorrNeutral(float ptcut, char *prodname, char *shortprodname, bool ispp = true, bool forSim = true, bool TPC, bool hadronic = false, float etacut = 0.7);
TH1D *GetHistoCorrNeutral(float cut, char *name, bool ispp, bool forSim, int mycase, bool eta, int color, int marker, bool hadronic = false);

Float_t CorrPtCut(float ptcut, char *prodname = "Enter Production Name", char *shortprodname = "EnterProductionName", bool ispp = true, bool forSim = true, int mycase = 0);
TH1D *GetHistoCorrPtCut(float ptcut = 0.15, char *name, bool ispp = true, bool forSim = true, int mycase);

TH1D *GetHistoCorrNotID(float etacut,char *name, bool TPC, bool eta, bool ispp = true, bool forSim = true);
TH1D *CorrNotID(float etacut,char *name, char *prodname, char *shortprodname, bool TPC, bool ispp = true, bool forSim = true);
Float_t CorrNotIDConst(float ptcut, float etacut,char *name, char *prodname, char *shortprodname, bool TPC, bool ispp, bool forSim);

TH1D *GetHistoNoID(float etacut, char *name, bool eta, bool TPC, bool ispp, bool forSim);
TH1D *CorrNoID(float etacut,char *name, char *prodname, char *shortprodname, bool ispp, bool forSim);
Float_t CorrNoIDConst(float etacut, float ptcut,char *name, char *prodname, char *shortprodname, bool ispp, bool forSim);

TH1D* bayneseffdiv(TH1D* numerator, TH1D* denominator,Char_t* name);
TH1D *GetHistoEfficiency(float cut, char *name, int mycase, int color, int marker,bool TPC, bool ITS, int cb = -1, int cblast = -1);
void CorrEfficiencyPlots(bool TPC, char *prodname, char *shortprodname);

TH1D *GetHistoCorrBkgd(float etacut,char *name, bool TPC,bool ispp,bool forSim);
void CorrBkgdPlots(char *prodname, char *shortprodname, bool TPC,bool ispp,bool forSim);

//Some variables that we'll use multiple times.  We'll declare them here since they don't seem to delete right in the functions
char prefix[100];
char histoname[100];
char epsname[100];
char pngname[100];
TFile *file = NULL;//initiated in main function
const char *mynameTPC = "TPC";
const char *mynameITS = "ITS";
const char *mynameTPCITS = "TPCITS";
const char *detectorEMCAL = "EMCAL";
const char *detectorPHOS = "PHOS";
const char *reweightedNo = "";
const char *reweightedYes = "Reweighted";

//===========================================================================================

void GetCorrections(char *prodname = "Enter Production Name", char *shortprodname = "EnterProductionName", bool ispp = true, bool forSim = true, bool TPC = true, char *infilename="Et.ESD.new.sim.merged.root", int dataset = 2009){
    TStopwatch timer;
    timer.Start();
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libXMLIO.so");

    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");

    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");

    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
   gROOT->ProcessLine(".L AliAnalysisEtCuts.cxx+g");
   gROOT->ProcessLine(".L AliAnalysisHadEtCorrections.cxx+g");
   file = new TFile(infilename);

   char outfilename[200];
   char *sim = "ForData";
   if(forSim) sim = "ForSimulations";
   char *system = "PbPb";
   if(ispp) system = "pp";
   sprintf(outfilename,"rootFiles/corrections/corrections.%s.%s.%s.root",shortprodname,system,sim);
   TFile *outfile = new TFile(outfilename,"RECREATE");
   AliAnalysisHadEtCorrections *hadCorrectionEMCAL = new AliAnalysisHadEtCorrections();
   hadCorrectionEMCAL->SetName("hadCorrectionEMCAL");
   float etacut = 0.7;
   hadCorrectionEMCAL->SetEtaCut(etacut);
   hadCorrectionEMCAL->IsData(!forSim);
   hadCorrectionEMCAL->IsEMCal(kTRUE);
   hadCorrectionEMCAL->SetProduction(shortprodname);
   hadCorrectionEMCAL->SetProductionDescription(prodname);
   hadCorrectionEMCAL->SetDataSet(dataset);
   //float etacut = hadCorrectionEMCAL->GetEtaCut();
   //cout<<"eta cut is "<<etacut<<endl;
   hadCorrectionEMCAL->SetAcceptanceCorrectionFull(1.0);
   cout<<"Warning:  Acceptance corrections will have to be updated to include real acceptance maps of the EMCAL"<<endl;
   hadCorrectionEMCAL->SetAcceptanceCorrectionPHOS(360.0/60.0);
   hadCorrectionEMCAL->SetAcceptanceCorrectionEMCAL(360.0/60.0);

   float ptcut = 0.1;
   float neutralCorr = CorrNeutral(ptcut,prodname,shortprodname,ispp,forSim,TPC,false,etacut);
   hadCorrectionEMCAL->SetNeutralCorrection(neutralCorr);
   //Using error from data, see analysis note for details
   if(ispp){
     hadCorrectionEMCAL->SetNeutralCorrectionLowBound(neutralCorr*(1.0-0.013));
     hadCorrectionEMCAL->SetNeutralCorrectionHighBound(neutralCorr*(1.0+0.013));
   }
   else{
     hadCorrectionEMCAL->SetNeutralCorrectionLowBound(neutralCorr*(1.0-0.049));
     hadCorrectionEMCAL->SetNeutralCorrectionHighBound(neutralCorr*(1.0+0.049));
   }

   float hadronicCorr = CorrNeutral(ptcut,prodname,shortprodname,ispp,forSim,TPC,true,etacut);
   hadCorrectionEMCAL->SetNotHadronicCorrection(hadronicCorr);
   if(ispp){
     hadCorrectionEMCAL->SetNotHadronicCorrectionLowBound(hadronicCorr*(1.0-0.008));
     hadCorrectionEMCAL->SetNotHadronicCorrectionHighBound(hadronicCorr*(1.0+0.008));
   }
   else{
     hadCorrectionEMCAL->SetNotHadronicCorrectionLowBound(hadronicCorr*(1.0-0.023));
     hadCorrectionEMCAL->SetNotHadronicCorrectionHighBound(hadronicCorr*(1.0+0.023));
   }

   float ptcutITS = CorrPtCut(0.1,prodname,shortprodname,ispp,forSim);
   hadCorrectionEMCAL->SetpTCutCorrectionITS(ptcutITS);
   float ptcutTPC = CorrPtCut(0.15,prodname,shortprodname,ispp,forSim);
   hadCorrectionEMCAL->SetpTCutCorrectionTPC(ptcutTPC);
   float ptcutITSLow = CorrPtCut(0.1,prodname,shortprodname,ispp,forSim,-1);
   float ptcutTPCLow = CorrPtCut(0.15,prodname,shortprodname,ispp,forSim,-1);
   hadCorrectionEMCAL->SetpTCutCorrectionITSLowBound(ptcutITSLow);
   hadCorrectionEMCAL->SetpTCutCorrectionTPCLowBound(ptcutTPCLow);
   float ptcutITSHigh = CorrPtCut(0.1,prodname,shortprodname,ispp,forSim,1);
   float ptcutTPCHigh = CorrPtCut(0.15,prodname,shortprodname,ispp,forSim,1);
   hadCorrectionEMCAL->SetpTCutCorrectionITSHighBound(ptcutITSHigh);
   hadCorrectionEMCAL->SetpTCutCorrectionTPCHighBound(ptcutTPCHigh);

   TH1D *NotIDTPC = CorrNotID(etacut,"CorrNotIDEMCALTPC",prodname,shortprodname,true,ispp,forSim);
   TH1D *NotIDITS = CorrNotID(etacut,"CorrNotIDEMCALITS",prodname,shortprodname,false,ispp,forSim);
   hadCorrectionEMCAL->SetNotIDCorrectionTPC(NotIDTPC);
   hadCorrectionEMCAL->SetNotIDCorrectionITS(NotIDITS);

   Float_t NotIDConstTPC = CorrNotIDConst(0.15,etacut,"CorrNotIDEMCALTPC2",prodname,shortprodname,true,ispp,forSim);
   Float_t NotIDConstITS = CorrNotIDConst(0.10,etacut,"CorrNotIDEMCALTPC2",prodname,shortprodname,true,ispp,forSim);
   hadCorrectionEMCAL->SetNotIDConstCorrectionTPC(1.0/NotIDConstTPC);
   hadCorrectionEMCAL->SetNotIDConstCorrectionITS(1.0/NotIDConstITS);
   if(ispp){
     hadCorrectionEMCAL->SetNotIDConstCorrectionTPCLowBound(1./NotIDConstTPC*(1.0-0.010));
     hadCorrectionEMCAL->SetNotIDConstCorrectionITSLowBound(1./NotIDConstITS*(1.0-0.010));
     hadCorrectionEMCAL->SetNotIDConstCorrectionTPCHighBound(1./NotIDConstTPC*(1.0+0.010));
     hadCorrectionEMCAL->SetNotIDConstCorrectionITSHighBound(1./NotIDConstITS*(1.0+0.010));
   }
   else{
     hadCorrectionEMCAL->SetNotIDConstCorrectionTPCLowBound(1./NotIDConstTPC*(1.0-0.022));
     hadCorrectionEMCAL->SetNotIDConstCorrectionITSLowBound(1./NotIDConstITS*(1.0-0.022));
     hadCorrectionEMCAL->SetNotIDConstCorrectionTPCHighBound(1./NotIDConstTPC*(1.0+0.022));
     hadCorrectionEMCAL->SetNotIDConstCorrectionITSHighBound(1./NotIDConstITS*(1.0+0.022));
   }

   TH1D *NoID = CorrNoID(etacut,"CorrNoIDEMCAL",prodname,shortprodname,ispp,forSim);
   hadCorrectionEMCAL->SetNotIDCorrectionNoPID(NoID);

   Float_t NoIDTPC = CorrNoIDConst(etacut,0.15,"CorrNoIDEMCAL2",prodname,shortprodname,ispp,forSim);
   Float_t NoIDITS = CorrNoIDConst(etacut,0.1,"CorrNoIDEMCAL2",prodname,shortprodname,ispp,forSim);
   hadCorrectionEMCAL->SetNotIDConstCorrectionTPCNoID(1./NoIDTPC);
   hadCorrectionEMCAL->SetNotIDConstCorrectionITSNoID(1./NoIDITS);
   hadCorrectionEMCAL->SetNotIDConstCorrectionTPCNoIDLowBound(1./NoIDTPC*.98);
   hadCorrectionEMCAL->SetNotIDConstCorrectionITSNoIDLowBound(1./NoIDITS*.98);
   hadCorrectionEMCAL->SetNotIDConstCorrectionTPCNoIDHighBound(1./NoIDTPC*1.02);
   hadCorrectionEMCAL->SetNotIDConstCorrectionITSNoIDHighBound(1./NoIDITS*1.02);
 
   //Here we're going to do a bit of a cheat.  We want the efficiency for ITS standalone tracks + TPC+ITS tracks.  This is returned by the function for the efficiency function if I ask for ITS only efficiency.  Had I known how this worked I probably would have written the code differently...  but...
   //anyhow I left a switch for changing it back.
   bool useITSStandalone = false;
   TH1D *efficiencyPionTPC = GetHistoEfficiency(etacut,"hEfficiencyPionTPC",1,1,20,useITSStandalone,true);
   hadCorrectionEMCAL->SetEfficiencyPionTPC(efficiencyPionTPC);
   TH1D *efficiencyKaonTPC = GetHistoEfficiency(etacut,"hEfficiencyKaonTPC",2,1,20,useITSStandalone,true);
   hadCorrectionEMCAL->SetEfficiencyKaonTPC(efficiencyKaonTPC);
   TH1D *efficiencyProtonTPC = GetHistoEfficiency(etacut,"hEfficiencyProtonTPC",3,1,20,useITSStandalone,true);
   hadCorrectionEMCAL->SetEfficiencyProtonTPC(efficiencyProtonTPC);
   TH1D *efficiencyHadronTPC = GetHistoEfficiency(etacut,"hEfficiencyHadronTPC",0,1,20,useITSStandalone,true);
   hadCorrectionEMCAL->SetEfficiencyHadronTPC(efficiencyHadronTPC);
   TH1D *efficiencyPionITS = GetHistoEfficiency(etacut,"hEfficiencyPionITS",1,1,20,false,true);
   hadCorrectionEMCAL->SetEfficiencyPionITS(efficiencyPionITS);
   TH1D *efficiencyKaonITS = GetHistoEfficiency(etacut,"hEfficiencyKaonITS",2,1,20,false,true);
   hadCorrectionEMCAL->SetEfficiencyKaonITS(efficiencyKaonITS);
   TH1D *efficiencyProtonITS = GetHistoEfficiency(etacut,"hEfficiencyProtonITS",3,1,20,false,true);
   hadCorrectionEMCAL->SetEfficiencyProtonITS(efficiencyProtonITS);
   TH1D *efficiencyHadronITS = GetHistoEfficiency(etacut,"hEfficiencyHadronITS",0,1,20,false,true);
   hadCorrectionEMCAL->SetEfficiencyHadronITS(efficiencyHadronITS);

   if(!ispp){
     TH1D *efficiencyPionTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB0",1,1,20,useITSStandalone,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB5",1,1,20,useITSStandalone,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB10",1,1,20,useITSStandalone,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionEMCAL->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB10->Clone(Form("Test%i",i)),i);

     TH1D *efficiencyKaonTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB0",2,1,20,useITSStandalone,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB5",2,1,20,useITSStandalone,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB10",2,1,20,useITSStandalone,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionEMCAL->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB10->Clone(Form("Test%i",i)),i);//Kaon

     TH1D *efficiencyProtonTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB0",3,1,20,useITSStandalone,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB5",3,1,20,useITSStandalone,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB10",3,1,20,useITSStandalone,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionEMCAL->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB10->Clone(Form("Test%i",i)),i);//Proton

     TH1D *efficiencyHadronTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB0",0,1,20,useITSStandalone,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB5",0,1,20,useITSStandalone,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB10",0,1,20,useITSStandalone,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionEMCAL->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB10->Clone(Form("Test%i",i)),i);//Hadron


     TH1D *efficiencyPionITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB0",1,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB5",1,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB10",1,1,20,false,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionEMCAL->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB10->Clone(Form("Test%i",i)),i);//Pion

     TH1D *efficiencyKaonITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB0",2,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB5",2,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB10",2,1,20,false,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionEMCAL->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB10->Clone(Form("Test%i",i)),i);//Kaon

     TH1D *efficiencyProtonITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB0",3,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB5",3,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB10",3,1,20,false,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionEMCAL->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB10->Clone(Form("Test%i",i)),i);//Proton

     TH1D *efficiencyHadronITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB0",0,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB5",0,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB10",0,1,20,false,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionEMCAL->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB10->Clone(Form("Test%i",i)),i);//Hadron
   }//EMCAL
   hadCorrectionEMCAL->SetEfficiencyErrorLowBound(0.99);
   hadCorrectionEMCAL->SetEfficiencyErrorHighBound(1.01);

   //CorrEfficiencyPlots(true,prodname,shortprodname);
   //CorrEfficiencyPlots(false,prodname,shortprodname,infilename);


   //hadCorrectionEMCAL->GetEfficiencyHadronTPC()->Draw();
   TH1D *backgroundTPC;
   TH1D *backgroundITS;
   if((dataset==20111 || dataset==20100) && !forSim){//2.76 TeV p+p or Pb+Pb
     if(dataset==20111){
       cout<<"Fixing 2.76 TeV p+p background to be average of 900 GeV and 7 TeV scaling"<<endl;
       backgroundTPC = pp276TPCBkgd();
       backgroundTPC->SetName("hBackgroundTPC");
       backgroundITS = pp276ITSBkgd();
       backgroundITS->SetName("hBackgroundITS");
     }
     else{//PbPb
       cout<<"Fixing 2.76 TeV Pb+Pb background to be average of 900 GeV and 7 TeV scaling with baryon enhancement"<<endl;
       backgroundTPC = pp276TPCBkgd();
       backgroundTPC->SetName("hBackgroundTPC");
       //ITS background is currently a placeholder
       backgroundITS = GetHistoCorrBkgd(etacut,"hBackgroundITS",false,ispp,forSim);
     }
   }
   else{
     backgroundTPC = GetHistoCorrBkgd(etacut,"hBackgroundTPC",true,ispp,forSim);
     backgroundITS = GetHistoCorrBkgd(etacut,"hBackgroundITS",false,ispp,forSim);
   }
   float bkgdpcterror = 0.0;
   switch(dataset){
   case 2009:
     bkgdpcterror = 0.37;
     break;
   case 20111:
     bkgdpcterror = 0.38;
     break;
   case 2010:
     bkgdpcterror = 0.13;
     break;
   case 20100:
     bkgdpcterror = 0.76;
     break;
   }
   hadCorrectionEMCAL->SetBackgroundCorrectionTPC(backgroundTPC);
   hadCorrectionEMCAL->SetBackgroundCorrectionITS(backgroundITS);
   hadCorrectionEMCAL->SetBackgroundErrorLowBound(1.0-bkgdpcterror/100.0);
   hadCorrectionEMCAL->SetBackgroundErrorHighBound(1.0+bkgdpcterror/100.0);
   //CorrBkgdPlots(prodname,shortprodname,true,ispp,forSim);
   //CorrBkgdPlots(prodname,shortprodname,false,ispp,forSim);

   hadCorrectionEMCAL->Report();

   outfile->cd();
   hadCorrectionEMCAL->Write();
   outfile->Write();
   delete hadCorrectionEMCAL;

   AliAnalysisHadEtCorrections *hadCorrectionPHOS = new AliAnalysisHadEtCorrections();
   hadCorrectionPHOS->SetName("hadCorrectionPHOS");
   float etacut = 0.12;
   hadCorrectionPHOS->SetEtaCut(etacut);
   hadCorrectionPHOS->IsData(!forSim);
   hadCorrectionPHOS->IsEMCal(kTRUE);
   hadCorrectionPHOS->SetProduction(shortprodname);
   hadCorrectionPHOS->SetProductionDescription(prodname);
   hadCorrectionPHOS->SetDataSet(dataset);
   //float etacut = hadCorrectionPHOS->GetEtaCut();
   //cout<<"eta cut is "<<etacut<<endl;
   hadCorrectionPHOS->SetAcceptanceCorrectionFull(1.0);
   cout<<"Warning:  Acceptance corrections will have to be updated to include real acceptance maps of the PHOS"<<endl;
   hadCorrectionPHOS->SetAcceptanceCorrectionPHOS(360.0/60.0);
   hadCorrectionPHOS->SetAcceptanceCorrectionEMCAL(360.0/60.0);

   float ptcut = 0.1;
   float neutralCorr = CorrNeutral(ptcut,prodname,shortprodname,ispp,forSim,TPC,false,etacut);
   hadCorrectionPHOS->SetNeutralCorrection(neutralCorr);
   //Using error from data, see analysis note for details
   if(ispp){
     hadCorrectionPHOS->SetNeutralCorrectionLowBound(neutralCorr*(1.0-.013));
     hadCorrectionPHOS->SetNeutralCorrectionHighBound(neutralCorr*(1.0+.013));
   }
   else{
     hadCorrectionPHOS->SetNeutralCorrectionLowBound(neutralCorr*(1.0-0.049));
     hadCorrectionPHOS->SetNeutralCorrectionHighBound(neutralCorr*(1.0+0.049));
   }

   float hadronicCorr = CorrNeutral(ptcut,prodname,shortprodname,ispp,forSim,TPC,true,etacut);
   hadCorrectionPHOS->SetNotHadronicCorrection(hadronicCorr);
   if(ispp){
     hadCorrectionPHOS->SetNotHadronicCorrectionLowBound(neutralCorr*(1.0-0.008));
     hadCorrectionPHOS->SetNotHadronicCorrectionHighBound(neutralCorr*(1.0+0.008));
   }
   else{
     hadCorrectionPHOS->SetNotHadronicCorrectionLowBound(neutralCorr*(1.0-0.023));
     hadCorrectionPHOS->SetNotHadronicCorrectionHighBound(neutralCorr*(1.0+0.023));
   }

   float ptcutITS = CorrPtCut(0.1,prodname,shortprodname,ispp,forSim);
   hadCorrectionPHOS->SetpTCutCorrectionITS(ptcutITS);
   float ptcutTPC = CorrPtCut(0.15,prodname,shortprodname,ispp,forSim);
   hadCorrectionPHOS->SetpTCutCorrectionTPC(ptcutTPC);

   float ptcutITSLow = CorrPtCut(0.1,prodname,shortprodname,ispp,forSim,-1);
   float ptcutTPCLow = CorrPtCut(0.15,prodname,shortprodname,ispp,forSim,-1);
   hadCorrectionPHOS->SetpTCutCorrectionITSLowBound(ptcutITSLow);
   hadCorrectionPHOS->SetpTCutCorrectionTPCLowBound(ptcutTPCLow);
   float ptcutITSHigh = CorrPtCut(0.1,prodname,shortprodname,ispp,forSim,1);
   float ptcutTPCHigh = CorrPtCut(0.15,prodname,shortprodname,ispp,forSim,1);
   hadCorrectionPHOS->SetpTCutCorrectionITSHighBound(ptcutITSHigh);
   hadCorrectionPHOS->SetpTCutCorrectionTPCHighBound(ptcutTPCHigh);

   TH1D *NotIDTPC = CorrNotID(etacut,"CorrNotIDPHOSTPC",prodname,shortprodname,true,ispp,forSim);
   TH1D *NotIDITS = CorrNotID(etacut,"CorrNotIDPHOSITS",prodname,shortprodname,false,ispp,forSim);
   hadCorrectionPHOS->SetNotIDCorrectionTPC(NotIDTPC);
   hadCorrectionPHOS->SetNotIDCorrectionITS(NotIDITS);

   Float_t NotIDConstTPC = CorrNotIDConst(0.15,etacut,"CorrNotIDPHOSTPC2",prodname,shortprodname,true,ispp,forSim);
   Float_t NotIDConstITS = CorrNotIDConst(0.10,etacut,"CorrNotIDPHOSTPC2",prodname,shortprodname,true,ispp,forSim);
   hadCorrectionPHOS->SetNotIDConstCorrectionTPC(1./NotIDConstTPC);
   hadCorrectionPHOS->SetNotIDConstCorrectionITS(1./NotIDConstITS);
   if(ispp){
     hadCorrectionPHOS->SetNotIDConstCorrectionTPCLowBound(1./NotIDConstTPC*(1.0-0.010));
     hadCorrectionPHOS->SetNotIDConstCorrectionITSLowBound(1./NotIDConstITS*(1.0-0.010));
     hadCorrectionPHOS->SetNotIDConstCorrectionTPCHighBound(1./NotIDConstTPC*(1.0+0.010));
     hadCorrectionPHOS->SetNotIDConstCorrectionITSHighBound(1./NotIDConstITS*(1.0+0.010));
   }
   else{
     hadCorrectionPHOS->SetNotIDConstCorrectionTPCLowBound(1./NotIDConstTPC*(1.0-0.022));
     hadCorrectionPHOS->SetNotIDConstCorrectionITSLowBound(1./NotIDConstITS*(1.0-0.022));
     hadCorrectionPHOS->SetNotIDConstCorrectionTPCHighBound(1./NotIDConstTPC*(1.0+0.022));
     hadCorrectionPHOS->SetNotIDConstCorrectionITSHighBound(1./NotIDConstITS*(1.0+0.022));
   }


   TH1D *NoID = CorrNoID(etacut,"CorrNoIDPHOS",prodname,shortprodname,ispp,forSim);
   hadCorrectionPHOS->SetNotIDCorrectionNoPID(NoID);


   Float_t NoIDTPC = CorrNoIDConst(etacut,0.15,"CorrNoIDPHOS2",prodname,shortprodname,ispp,forSim);
   Float_t NoIDITS = CorrNoIDConst(etacut,0.1,"CorrNoIDPHOS2",prodname,shortprodname,ispp,forSim);
   hadCorrectionPHOS->SetNotIDConstCorrectionTPCNoID(1./NoIDTPC);
   hadCorrectionPHOS->SetNotIDConstCorrectionITSNoID(1./NoIDITS);
   hadCorrectionPHOS->SetNotIDConstCorrectionTPCNoIDLowBound(1./NoIDTPC*.98);
   hadCorrectionPHOS->SetNotIDConstCorrectionITSNoIDLowBound(1./NoIDITS*.98);
   hadCorrectionPHOS->SetNotIDConstCorrectionTPCNoIDHighBound(1./NoIDTPC*1.02);
   hadCorrectionPHOS->SetNotIDConstCorrectionITSNoIDHighBound(1./NoIDITS*1.02);

   TH1D *efficiencyPionTPC = GetHistoEfficiency(etacut,"hEfficiencyPionTPC",1,1,20,useITSStandalone,true);
   TH1D *efficiencyKaonTPC = GetHistoEfficiency(etacut,"hEfficiencyKaonTPC",2,1,20,useITSStandalone,true);
   TH1D *efficiencyProtonTPC = GetHistoEfficiency(etacut,"hEfficiencyProtonTPC",3,1,20,useITSStandalone,true);
   TH1D *efficiencyHadronTPC = GetHistoEfficiency(etacut,"hEfficiencyHadronTPC",0,1,20,useITSStandalone,true);
   TH1D *efficiencyPionITS = GetHistoEfficiency(etacut,"hEfficiencyPionITS",1,1,20,false,true);
   TH1D *efficiencyKaonITS = GetHistoEfficiency(etacut,"hEfficiencyKaonITS",2,1,20,false,true);
   TH1D *efficiencyProtonITS = GetHistoEfficiency(etacut,"hEfficiencyProtonITS",3,1,20,false,true);
   TH1D *efficiencyHadronITS = GetHistoEfficiency(etacut,"hEfficiencyHadronITS",0,1,20,false,true);
   //CorrEfficiencyPlots(true,prodname,shortprodname);
   //CorrEfficiencyPlots(false,prodname,shortprodname);
   hadCorrectionPHOS->SetEfficiencyPionTPC(efficiencyPionTPC);
   hadCorrectionPHOS->SetEfficiencyKaonTPC(efficiencyKaonTPC);
   hadCorrectionPHOS->SetEfficiencyProtonTPC(efficiencyProtonTPC);
   hadCorrectionPHOS->SetEfficiencyHadronTPC(efficiencyHadronTPC);
   hadCorrectionPHOS->SetEfficiencyPionITS(efficiencyPionITS);
   hadCorrectionPHOS->SetEfficiencyKaonITS(efficiencyKaonITS);
   hadCorrectionPHOS->SetEfficiencyProtonITS(efficiencyProtonITS);
   hadCorrectionPHOS->SetEfficiencyHadronITS(efficiencyHadronITS);


   if(!ispp){
     TH1D *efficiencyPionTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB0",1,1,20,useITSStandalone,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB5",1,1,20,useITSStandalone,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB10",1,1,20,useITSStandalone,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionPHOS->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB10->Clone(Form("Test%i",i)),i);

     TH1D *efficiencyKaonTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB0",2,1,20,useITSStandalone,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB5",2,1,20,useITSStandalone,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB10",2,1,20,useITSStandalone,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionPHOS->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB10->Clone(Form("Test%i",i)),i);//Kaon

     TH1D *efficiencyProtonTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB0",3,1,20,useITSStandalone,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB5",3,1,20,useITSStandalone,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB10",3,1,20,useITSStandalone,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionPHOS->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB10->Clone(Form("Test%i",i)),i);//Proton

     TH1D *efficiencyHadronTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB0",0,1,20,useITSStandalone,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB5",0,1,20,useITSStandalone,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB10",0,1,20,useITSStandalone,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionPHOS->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB10->Clone(Form("Test%i",i)),i);//Hadron


     TH1D *efficiencyPionITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB0",1,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB5",1,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB10",1,1,20,false,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionPHOS->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB10->Clone(Form("Test%i",i)),i);//Pion

     TH1D *efficiencyKaonITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB0",2,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB5",2,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB10",2,1,20,false,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionPHOS->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB10->Clone(Form("Test%i",i)),i);//Kaon

     TH1D *efficiencyProtonITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB0",3,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB5",3,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB10",3,1,20,false,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionPHOS->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB10->Clone(Form("Test%i",i)),i);//Proton

     TH1D *efficiencyHadronITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB0",0,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB5",0,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB10",0,1,20,false,true,10,15);
     for(int i=10;i<=20;i++) hadCorrectionPHOS->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB10->Clone(Form("Test%i",i)),i);//Hadron
   }//EMCAL


   TH1D *backgroundTPC;
   TH1D *backgroundITS;
   if((dataset==20111 || dataset==20100) && !forSim){//2.76 TeV p+p or Pb+Pb
     if(dataset==20111){
       cout<<"Fixing 2.76 TeV p+p background to be average of 900 GeV and 7 TeV scaling"<<endl;
       backgroundTPC = pp276TPCBkgd();
       backgroundTPC->SetName("hBackgroundTPCPHOS");
       backgroundITS = pp276ITSBkgd();
       backgroundITS->SetName("hBackgroundITSPHOS");
     }
     else{//PbPb
       cout<<"Fixing 2.76 TeV Pb+Pb background to be average of 900 GeV and 7 TeV scaling with baryon enhancement"<<endl;
       backgroundTPC = pp276TPCBkgd();
       backgroundTPC->SetName("hBackgroundTPCPHOS");
       //ITS background is currently a placeholder
       backgroundITS = GetHistoCorrBkgd(etacut,"hBackgroundITS",false,ispp,forSim);
     }
   }
   else{
     backgroundTPC = GetHistoCorrBkgd(etacut,"hBackgroundTPC",true,ispp,forSim);
     backgroundITS = GetHistoCorrBkgd(etacut,"hBackgroundITS",false,ispp,forSim);
   }
   hadCorrectionPHOS->SetBackgroundCorrectionTPC(backgroundTPC);
   hadCorrectionPHOS->SetBackgroundCorrectionITS(backgroundITS);
   hadCorrectionPHOS->SetBackgroundErrorLowBound(1.0-0.001);
   hadCorrectionPHOS->SetBackgroundErrorHighBound(1.0+0.001);
   //CorrBkgdPlots(prodname,shortprodname,true,ispp,forSim);
   //CorrBkgdPlots(prodname,shortprodname,false,ispp,forSim);

   hadCorrectionPHOS->Report();
   //Write the output
   outfile->cd();
   hadCorrectionPHOS->Write();
   outfile->Write();
   outfile->Close();


  timer.Stop();
  timer.Print();
}

//==================================CorrNeutral==============================================
Float_t CorrNeutral(float ptcut, char *prodname, char *shortprodname, bool ispp, bool forSim, bool TPC, bool hadronic, float etacut){
  if(!forSim){//for data we have evaluated the neutral correction from ALICE data
    if(hadronic){//for tot et from had et
      if(ispp){
	return 1.0/0.571;
      }
      else{
	return 1.0/0.549;
      }
    }
    else{//for had et only
      if(ispp){
	return 1.0/0.736;
      }
      else{
	return 1.0/0.689;
      }
    }
  }
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",400,400);
  c->SetTopMargin(0.0);
  c->SetRightMargin(0.0);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);

  TPad *ptpad = c->cd(1);
  ptpad->SetTopMargin(0.04);
  ptpad->SetRightMargin(0.04);
  ptpad->SetLeftMargin(0.149288);
  ptpad->SetBorderSize(0);
  ptpad->SetFillColor(0);
  ptpad->SetFillColor(0);
  ptpad->SetBorderMode(0);
  ptpad->SetFrameFillColor(0);
  ptpad->SetFrameBorderMode(0);

  int phosmarker = 20;

  sprintf(prefix,"%s%2.1f",shortprodname,ptcut);

  sprintf(histoname,"%stotal",histoname);
  int colortotal = 1;
  int casetotal = 4;
  if(hadronic) casetotal = 8;
  TH1D *total = GetHistoCorrNeutral(ptcut,histoname,ispp,forSim,casetotal,false,colortotal,phosmarker,hadronic);

  int colorallneutral = 2;
  TH1D *allneutral = GetHistoCorrNeutral(ptcut,"allneutral",ispp,forSim,3,false,colorallneutral,phosmarker,hadronic);

  int colorchargedsecondary = TColor::kViolet-3;
  TH1D *chargedsecondary = GetHistoCorrNeutral(ptcut,"chargedsecondary",ispp,forSim,2,false,colorchargedsecondary,phosmarker,hadronic);

  int colorneutralUndet = 4;
  TH1D *neutralUndet = GetHistoCorrNeutral(ptcut,"neutralUndet",ispp,forSim,1,false,colorneutralUndet,phosmarker,hadronic);

  int colorv0 = TColor::kGreen+2;
  TH1D *v0 = GetHistoCorrNeutral(ptcut,"v0",ispp,forSim,0,false,colorv0,phosmarker,hadronic);

  int colorem = TColor::kCyan;
  TH1D *em = GetHistoCorrNeutral(ptcut,"em",ispp,forSim,9,false,colorem,phosmarker,hadronic);

  TF1 *func = new TF1("func","[0]",-.7,.7);
  func->SetParameter(0,0.2);
  total->Fit(func,"","",-etacut,etacut);

  //total->SetAxisRange(0.0,4);
  total->GetXaxis()->SetLabelSize(0.05);
  total->GetYaxis()->SetLabelSize(0.045);
  total->GetXaxis()->SetTitleSize(0.05);
  total->GetYaxis()->SetTitleSize(0.06);
  if(hadronic){
    total->SetMaximum(0.6);
    total->SetMinimum(0.0);
  }
  else{
    total->SetMaximum(0.3);
    total->SetMinimum(0.0);
  }
  total->Draw();
  allneutral->Draw("same");
  chargedsecondary->Draw("same");
  neutralUndet->Draw("same");
  v0->Draw("same");
  if(hadronic) em->Draw("same");

  TLatex *tex = new TLatex(0.161478,1.0835,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLegend *leg2 = new TLegend(0.518321,0.746873,0.774812,0.955343);
  leg2->AddEntry(total,"Total");
  leg2->AddEntry(allneutral,"#Lambda,#bar{#Lambda},K^{0}_{S},K^{0}_{L},n,#bar{n}");
  leg2->AddEntry(neutralUndet,"K^{0}_{L},n,#bar{n}");
  if(hadronic) leg2->AddEntry(em,"e^{#pm},#gamma,#eta,#pi^{0},#omega");
  leg2->SetFillStyle(0);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.0548607);
  leg2->Draw();
  if(hadronic){
    sprintf(epsname,"pics/%s/fhadronic.eps",shortprodname);
    sprintf(pngname,"pics/%s/fhadronic.png",shortprodname);
  }
  else{
    sprintf(epsname,"pics/%s/fneutral.eps",shortprodname);
    sprintf(pngname,"pics/%s/fneutral.png",shortprodname);
  }
  c->SaveAs(epsname);
  c->SaveAs(pngname);

  delete total;
  delete allneutral;
  delete chargedsecondary;
  delete neutralUndet;
  delete v0;
  delete em;
  delete c;
  float corr = func->GetParameter(0);
  delete func;
  delete tex;
  delete leg2;
  return 1.0/(1.0-corr);

}
TH1D *GetHistoCorrNeutral(float cut, char *name, bool ispp, bool forSim, int mycase, bool eta, int color, int marker, bool hadronic){
  file->cd();
  char *reweightname = "";
  if(!forSim) reweightname = "Reweighted";
  TH2F *numeratorParent; 
  switch(mycase){
  case 0:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedLambda%s",reweightname)))->Clone("v0");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)));
    break;
  case 1:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)))->Clone("Knnbar");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"));
    break;
  case 2:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedOmega"))->Clone("ch2ndary");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"));
    break;
  case 3:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedLambda%s",reweightname)))->Clone("allneutral");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"));
    break;
  case 4:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedLambda%s",reweightname)))->Clone("allneutral");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedOmega"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"));
    break;
  case 5:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedXi"))->Clone("allxi");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"));
    break;
  case 6:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedOmega"))->Clone("allomega");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"));
    break;
  case 7:
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject("EtSimulatedSigma"))->Clone("allsigma");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiSigma"));
    break;
  case 8:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject(Form("EtSimulatedLambda%s",reweightname)))->Clone("allneutral");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedOmega"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedXi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedGamma"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEta"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedPi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedOmega0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEPlus"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEMinus"));
    break;
  case 9:
    numeratorParent= (TH2F*)((TH2F*) out2->FindObject("EtSimulatedGamma"))->Clone("allem");
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEta"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedPi0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedOmega0"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEPlus"));
    numeratorParent->Add((TH2F*) out2->FindObject("EtSimulatedEMinus"));
    break;
  }

  TH2F *allhad;
  //allhad=(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedAllHadron"))->Clone("id");
  allhad=(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedPiPlus"))->Clone("id");
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedPiMinus"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedKMinus"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedKPlus"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedProton"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiProton"));
  allhad->Add((TH2F*) out2->FindObject(Form("EtSimulatedLambda%s",reweightname)));
  allhad->Add((TH2F*) out2->FindObject(Form("EtSimulatedAntiLambda%s",reweightname)));
  allhad->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0S%s",reweightname)));
  allhad->Add((TH2F*) out2->FindObject(Form("EtSimulatedK0L%s",reweightname)));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedNeutron"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiNeutron"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedEPlus"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedEMinus"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedOmega"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiOmega"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedXi"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedSigma"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiSigma"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedXi0"));
  allhad->Add((TH2F*) out2->FindObject("EtSimulatedAntiXi0"));

  if(hadronic){//if we are getting the correction for the hadronic only case...    
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedGamma"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedEta"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedPi0"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedOmega0"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedEPlus"));
    allhad->Add((TH2F*) out2->FindObject("EtSimulatedEMinus"));
  }

  TH1D *denominator;
  TH1D *numerator;
  if(eta){
    int lowbin = numeratorParent->GetYaxis()->FindBin(-cut+.001);//make sure we don't accv0entally get the wrong bin
    int highbin = numeratorParent->GetYaxis()->FindBin(cut-.001);
    //cout<<"Projecting from "<<numeratorParent->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<numeratorParent->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;
    denominator = allhad->ProjectionX("name",lowbin,highbin);
    numerator = numeratorParent->ProjectionX("numerator",lowbin,highbin);
  }
  else{
    int lowbin = allhad->GetXaxis()->FindBin(cut);//make sure we don't accidentally get the wrong bin
    int highbin = allhad->GetXaxis()->GetNbins();
    //cout<<"Projecting from "<<allhad->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<allhad->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
    numerator = numeratorParent->ProjectionY("name",lowbin,highbin);
    denominator = allhad->ProjectionY("denominator",lowbin,highbin);
  }
  numerator->Divide(denominator);
  //numerator->Rebin(2);
  //numerator->Scale(0.5);
  numerator->SetYTitle("E_{T}^{had,sample}/E_{T}^{had,total}");
  numerator->GetYaxis()->SetTitleOffset(1.2);
  numerator->SetMarkerColor(color);
  numerator->SetLineColor(color);
  numerator->SetMarkerStyle(marker);
  delete denominator;
  delete numeratorParent;
  delete allhad;
  //file->Close();
  numerator->SetName(name);
  return numerator;

}

//===============================CorrPtCut=========================================
TH1D *GetHistoCorrPtCut(float ptcut, char *name, bool ispp, bool forSim, int mycase){
  file->cd();
  TH2F *allhad = ((TH2F*) out2->FindObject("EtSimulatedAllHadron"))->Clone("allhad");
  TH2F *ptlow = ((TH2F*) out2->FindObject("EtSimulatedChargedHadronAssumingNoPt"))->Clone("ptlow");
  TH2F *pthigh;
  if(ptcut>0.14){//TPC cut off
    (TH2F*)pthigh =(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedChargedHadronAssumingPtTPCCut"))->Clone("pthigh");
  }
  else{
    (TH2F*)pthigh =(TH2F*) ((TH2F*) out2->FindObject("EtSimulatedChargedHadronAssumingPtITSCut"))->Clone("pthigh");
  }

  int lowbin = allhad->GetXaxis()->FindBin(0.0);//make sure we don't accidentally get the wrong bin
  int highbin = allhad->GetXaxis()->FindBin(ptcut);
  int nbins = allhad->GetXaxis()->GetNbins();
  //cout<<"Projecting from "<<allhad->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<allhad->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
  //cout<<"Projecting from "<<allhad->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<allhad->GetXaxis()->GetBinLowEdge(nbins)<<endl;

  //allhad->Sumw2();
  TH1D *numerator;
  TH1D *denominator;
  switch(mycase){
  case -1:
    numerator = ptlow->ProjectionY("nameLow",lowbin,highbin);
    denominator = allhad->ProjectionY("denominatorLow",highbin,nbins);
    denominator->Add(ptlow);
    break;
  case 1:
    numerator = pthigh->ProjectionY("nameHigh",lowbin,highbin);
    denominator = allhad->ProjectionY("denominatorHigh",highbin,nbins);
    denominator->Add(pthigh);
    break;
  default:
    numerator = allhad->ProjectionY("name",lowbin,highbin);
    denominator = allhad->ProjectionY("denominator",lowbin,nbins);
  }
  numerator->Divide(denominator);
  numerator->SetYTitle("E_{T}^{had, p_{T}<cut-off}/E_{T}^{had, all p_{T}}");
  numerator->GetYaxis()->SetTitleOffset(1.);
  numerator->GetYaxis()->SetTitleSize(0.08);
  numerator->GetYaxis()->SetLabelSize(0.05);
  numerator->GetXaxis()->SetTitleSize(0.08);
  numerator->GetXaxis()->SetLabelSize(0.05);
  numerator->GetXaxis()->SetTitleOffset(.6);
  //numerator->Rebin(2);
  //numerator->Scale(0.5);
  //numerator->Draw("e");
  delete allhad;
  delete denominator;
  delete ptlow;
  delete pthigh;
  numerator->SetName(name);
  return numerator;

}

Float_t CorrPtCut(float ptcut, char *prodname, char *shortprodname, bool ispp, bool forSim, int mycase){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",500,400);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.04);
  c->SetLeftMargin(0.181452);
  c->SetBottomMargin(0.134409);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);



  TH1D *High = GetHistoCorrPtCut(0.15-.001,"High",ispp,forSim,mycase);
  TH1D *Low = GetHistoCorrPtCut(0.1-.001,"Low",ispp,forSim,mycase);
  TH1D *Lowest = GetHistoCorrPtCut(0.05-.001,"Lowest",ispp,forSim,mycase);

  TF1 *func = new TF1("func","[0]",-.7,.7);
  func->SetParameter(0,0.2);
  if(ptcut<.125){//its cuts
    Low->Fit(func);
  }
  else{//tpc cuts
    High->Fit(func);
  }

  High->SetMaximum(0.04);
  High->SetMinimum(0.0);
  High->SetMarkerColor(2);
  Low->SetMarkerColor(4);
  High->SetLineColor(2);
  Low->SetLineColor(4);
  High->SetMinimum(0.0);
  High->SetMarkerStyle(20);
  Low->SetMarkerStyle(21);
  Lowest->SetMarkerStyle(22);
  High->Draw();
  Low->Draw("same");
  Lowest->Draw("same");



  TLatex *tex = new TLatex(-0.723444,0.0373593,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLegend *leg = new TLegend(0.217742,0.66129,0.477823,0.873656);
  leg->AddEntry(High,"p_{T} cut-off = 0.15 GeV/c");
  leg->AddEntry(Low,"p_{T} cut-off = 0.1 GeV/c");
  leg->AddEntry(Lowest,"p_{T} cut-off = 0.05 GeV/c");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.0537634);
  leg->Draw();

  if(ptcut<.125){//its cuts
    c->SaveAs(Form("pics/%s/fptcutITS.eps",shortprodname));
    c->SaveAs(Form("pics/%s/fptcutITS.png",shortprodname));
  }
  else{//tpc cuts
    c->SaveAs(Form("pics/%s/fptcutTPC.eps",shortprodname));
    c->SaveAs(Form("pics/%s/fptcutTPC.png",shortprodname));
  }

  float corr = func->GetParameter(0);
  //cout<<"Pt cut correction: "<<1.0/(1.0-corr)<<endl;
  delete High;
  delete Low;
  delete Lowest;
  delete func;
  delete c;
  delete tex;
  delete leg;
  return 1.0/(1.0-corr);
}



//==================================CorrNotID=================================================
TH1D *GetHistoCorrNotID(float etacut,char *name, bool TPC, bool eta, bool ispp, bool forSim){
  file->cd();
  char *myname = mynameITS;
  if(TPC) myname = mynameTPCITS;
  TH2F *notid = ((TH2F*) out2->FindObject(Form("EtReconstructed%sUnidentifiedAssumingPion",myname)))->Clone("notid");
  TH2F *nNotid = ((TH2F*) out2->FindObject(Form("EtNReconstructed%sUnidentified",myname)))->Clone("nNotid");
  if(!eta){
    //cout<<"Correction determined for all charged hadrons"<<endl;
    notid->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedPiPlus",myname)));
    notid->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedPiMinus",myname)));
    notid->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedKPlus",myname)));
    notid->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedKMinus",myname)));
    notid->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedProton",myname)));
    notid->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedAntiProton",myname)));
    nNotid->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%sPiPlus",myname)));
    nNotid->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%sPiMinus",myname)));
    nNotid->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%sKPlus",myname)));
    nNotid->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%sKMinus",myname)));
    nNotid->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%sProton",myname)));
    nNotid->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%sAntiProton",myname)));
  }

  TH2F *id = ((TH2F*) out2->FindObject(Form("EtReconstructed%sUnidentified",myname)))->Clone("id");
  if(!eta){
    id->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedPiPlus",myname)));
    id->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedPiMinus",myname)));
    id->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedKPlus",myname)));
    id->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedKMinus",myname)));
    id->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedProton",myname)));
    id->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedAntiProton",myname)));
  }

  TH1D *nNotidProj;
  TH1D *denominator;
  TH1D *numerator;
  if(eta){
    int lowbin = notid->GetYaxis()->FindBin(-etacut+.001);//make sure we don't accv0entally get the wrong bin
    int highbin = notid->GetYaxis()->FindBin(etacut-.001);
    //cout<<"Projecting from "<<notid->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<notid->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;
    denominator = id->ProjectionX("name",lowbin,highbin);
    numerator = notid->ProjectionX("numerator",lowbin,highbin);
    nNotidProj = nNotid->ProjectionX("nNotidProj",lowbin,highbin);
  }
  else{
    //cout<<"Getting eta dependence"<<endl;
    int lowbin = id->GetXaxis()->FindBin(etacut);//make sure we don't accidentally get the wrong bin
    int highbin;
    if(etacut<0.15){//then we actually have ITS standalone tracks and we only want this to run from 0.1 to 0.15 because this will be used only for ITS standalone tracks
      highbin = id->GetXaxis()->FindBin(0.15);
    }
    else{
      highbin = id->GetXaxis()->GetNbins();
    }
    //cout<<"Projecting from "<<id->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<id->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
    numerator = notid->ProjectionY("name",lowbin,highbin);
    denominator = id->ProjectionY("denominator",lowbin,highbin);
    nNotidProj = nNotid->ProjectionY("nNotidProj",lowbin,highbin);
  }
  TH1D *result = numerator;
  if(!denominator){
    cerr<<"Uh-oh!  Can't find denominator!!";
    return numerator;
  }
  else{result->Divide(denominator);}
  if(result->GetNbinsX() != nNotidProj->GetNbinsX()){
    cerr<<"Uh-oh!  Can't rescale errors! "<<result->GetNbinsX()<<"!="<<nNotidProj->GetNbinsX()<<endl;
    return result;
  }
  if(!nNotidProj){
    cerr<<"Uh-oh!  Can't find histogram!!";
    return numerator;
  }
  //fixing the errors
  for(int i=1;i<=result->GetNbinsX();i++){
    Float_t value = result->GetBinContent(i);
    Float_t valueerr = result->GetBinError(i);
    Float_t n = nNotidProj->GetBinContent(i);
    Float_t err;
    if(n<=0) err = 0.0;
    else{err= value/TMath::Power(n,0.5);}
    result->SetBinError(i,err);
    //cout<<"Was "<<valueerr<<", setting to "<<err<<endl;
  }
  result->SetYTitle("Ratio of E_{T}^{assuming pion}/E_{T}^{real}");
  result->GetYaxis()->SetTitleOffset(1.2);
  delete denominator;
  delete nNotidProj;
  delete notid;
  delete nNotid;
  delete id;
  result->SetName(name);
  return result;

}

TH1D *CorrNotID(float etacut,char *name, char *prodname, char *shortprodname, bool TPC, bool ispp, bool forSim){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",500,400);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.04);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);

  TH1D *PHOS = GetHistoCorrNotID(etacut,name,TPC,true,ispp,forSim);
  PHOS->SetMarkerColor(2);
  PHOS->SetLineColor(2);
  PHOS->SetAxisRange(0.0,4);
  if(TPC){
    PHOS->SetMaximum(1.1);
    PHOS->SetMinimum(0.85);
  }
  else{
    PHOS->SetAxisRange(0.0,0.5);
  }
  PHOS->SetMarkerStyle(20);
  PHOS->Draw();
  TLatex *tex = new TLatex(0.161478,1.0835,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  char *detector = detectorEMCAL;
  if(etacut<0.2) detector = detectorPHOS;
  if(TPC){
    sprintf(epsname,"pics/%s/fnotidTPC%s.eps",shortprodname,detector);
    sprintf(pngname,"pics/%s/fnotidTPC%s.png",shortprodname,detector);
  }
  else{
    sprintf(epsname,"pics/%s/fnotidITS%s.eps",shortprodname,detector);
    sprintf(pngname,"pics/%s/fnotidITS%s.png",shortprodname,detector);
  }

  c->SaveAs(epsname);
  c->SaveAs(pngname);
  delete c;
  delete tex;
  PHOS->SetName(name);
  return PHOS;
}

Float_t CorrNotIDConst(float ptcut, float etacut,char *name, char *prodname, char *shortprodname, bool TPC, bool ispp, bool forSim){
  if(!forSim){
    if(ispp){
      return 0.996;
    }
    else{
      return 0.976;
    }
  }
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",500,400);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.04);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);

  TH1D *PHOS = GetHistoCorrNotID(ptcut,name,TPC,false,ispp,forSim);
  PHOS->SetMarkerColor(2);
  PHOS->SetLineColor(2);
  PHOS->SetMaximum(1.01);
  PHOS->SetMinimum(0.98);
  TF1 *func = new TF1("func","[0]",-etacut,etacut);
  PHOS->Fit(func,"","",-etacut,etacut);
  PHOS->SetMarkerStyle(20);
  PHOS->Draw();
  TLatex *tex = new TLatex(0.161478,1.0835,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  char *detector = detectorEMCAL;
  if(etacut<0.2) detector = detectorPHOS;
  if(TPC){
    sprintf(epsname,"pics/%s/fnotidConstTPC%s.eps",shortprodname,detector);
    sprintf(pngname,"pics/%s/fnotidConstTPC%s.png",shortprodname,detector);
  }
  else{
    sprintf(epsname,"pics/%s/fnotidConstITS%s.eps",shortprodname,detector);
    sprintf(pngname,"pics/%s/fnotidConstITS%s.png",shortprodname,detector);
  }

  c->SaveAs(epsname);
  c->SaveAs(pngname);
  delete c;
  delete PHOS;
  float value = func->GetParameter(0);
  delete func;
  delete tex;
  return value;
}

//==================================CorrNoID=================================================
TH1D *GetHistoNoID(float etacut, char *name, bool eta, bool TPC, bool ispp, bool forSim){
  file->cd();
  char *myname = mynameITS;
  if(TPC) myname = mynameTPCITS;
  TH2F *notid = ((TH2F*) out2->FindObject(Form("EtReconstructed%sChargedHadronAssumingPion",myname)))->Clone("notid");
  TH2F *nNotid = ((TH2F*) out2->FindObject(Form("EtNReconstructed%sChargedHadron",myname)))->Clone("nNotid");

  TH2F *id = ((TH2F*) out2->FindObject(Form("EtReconstructed%sChargedHadron",myname)))->Clone("id");
  int lowbin = id->GetYaxis()->FindBin(-etacut+.001);//make sure we don't accidentally get the wrong bin
  int highbin = id->GetYaxis()->FindBin(etacut-.001);

  TH1D *nNotidProj;
  TH1D *denominator;
  TH1D *numerator;
  if(eta){
    int lowbin = notid->GetYaxis()->FindBin(-etacut+.001);//make sure we don't accv0entally get the wrong bin
    int highbin = notid->GetYaxis()->FindBin(etacut-.001);
    //cout<<"Projecting from "<<notid->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<notid->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;
    denominator = id->ProjectionX("name",lowbin,highbin);
    numerator = notid->ProjectionX("numerator",lowbin,highbin);
    nNotidProj = nNotid->ProjectionX("nNotidProj",lowbin,highbin);
  }
  else{
    //cout<<"Getting eta dependence"<<endl;
    int lowbin = id->GetXaxis()->FindBin(etacut);//make sure we don't accidentally get the wrong bin
    int highbin;
    if(etacut<0.15){//then we actually have ITS standalone tracks and we only want this to run from 0.1 to 0.15 because this will be used only for ITS standalone tracks
      highbin = id->GetXaxis()->FindBin(0.15);
    }
    else{
      highbin = id->GetXaxis()->GetNbins();
    }
    //cout<<"Projecting from "<<id->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<id->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
    numerator = notid->ProjectionY("name",lowbin,highbin);
    denominator = id->ProjectionY("denominator",lowbin,highbin);
    nNotidProj = nNotid->ProjectionY("nNotidProj",lowbin,highbin);
  }
  if(denominator) numerator->Divide(denominator);

  if(numerator->GetNbinsX() != nNotidProj->GetNbinsX()){
    cerr<<"Uh-oh!  Can't rescale errors! "<<numerator->GetNbinsX()<<"!="<<nNotidProj->GetNbinsX()<<endl;
    return numerator;
  }
  //fixing the errors
  for(int i=1;i<=numerator->GetNbinsX();i++){
    Float_t value = numerator->GetBinContent(i);
    Float_t valueerr = numerator->GetBinError(i);
    Float_t n = nNotidProj->GetBinContent(i);
    Float_t err=0.;
    if(n>0.0){
      err = value/TMath::Power(n,0.5);
    }
    numerator->SetBinError(i,err);;
  }
  numerator->SetYTitle("Ratio of E_{T}^{assuming pion}/E_{T}^{real}");
  numerator->GetYaxis()->SetTitleOffset(1.2);
  delete denominator;
  delete nNotidProj;
  delete notid;
  delete nNotid;
  delete id;
  numerator->SetName(name);
  return numerator;

}

TH1D *CorrNoID(float etacut,char *name, char *prodname, char *shortprodname, bool ispp, bool forSim){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",500,400);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.04);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);

  TH1D *PHOS = GetHistoNoID(etacut,name,true,true,ispp,forSim);
  PHOS->SetMarkerColor(2);
  PHOS->SetLineColor(2);
  PHOS->SetAxisRange(0.0,4);
  PHOS->SetMaximum(1.1);
  PHOS->SetMinimum(0.85);
  PHOS->SetMarkerStyle(20);;
  PHOS->Draw();
  TLatex *tex = new TLatex(0.161478,1.0835,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();


  char *detector = detectorEMCAL;
  if(etacut<0.2) detector = detectorPHOS;
  sprintf(epsname,"pics/%s/fnoid%s.eps",shortprodname,detector);
  sprintf(pngname,"pics/%s/fnoid%s.png",shortprodname,detector);

  c->SaveAs(epsname);
  c->SaveAs(pngname);
  delete c;
  delete tex;
  PHOS->SetName(name);
  return PHOS;

}

Float_t CorrNoIDConst(float etacut, float ptcut,char *name, char *prodname, char *shortprodname, bool ispp, bool forSim){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",500,400);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.04);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);

  bool TPC = true;
  if(ptcut<.15) TPC = false;
  TH1D *PHOS = GetHistoNoID(ptcut,name,false,TPC,ispp,forSim);
  TF1 *func = new TF1("func","[0]",-etacut,etacut);
  PHOS->Fit(func,"","",-etacut,etacut);
  PHOS->SetMarkerColor(2);
  PHOS->SetLineColor(2);
  PHOS->SetMaximum(1.1);
  PHOS->SetMinimum(0.85);
  PHOS->SetMarkerStyle(20);;
  PHOS->Draw();
  TLatex *tex = new TLatex(0.161478,1.0835,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();


  char *detector = detectorEMCAL;
  if(etacut<0.2) detector = detectorPHOS;
  if(TPC){
    sprintf(epsname,"pics/%s/fnoid%sTPC.eps",shortprodname,detector);
    sprintf(pngname,"pics/%s/fnoid%sTPC.png",shortprodname,detector);
  }
  else{
    sprintf(epsname,"pics/%s/fnoid%sITS.eps",shortprodname,detector);
    sprintf(pngname,"pics/%s/fnoid%sITS.png",shortprodname,detector);
  }

  c->SaveAs(epsname);
  c->SaveAs(pngname);
  delete c;
  delete PHOS;
  float value = func->GetParameter(0);
  delete func;
  delete tex;
  return value;

}
//==================================Efficiency=================================================
TH1D* bayneseffdiv(TH1D* numerator, TH1D* denominator,Char_t* name) 
{
    if(!numerator){
      cerr<<"Error:  numerator does not exist!"<<endl;
      return NULL;
    }
    if(!denominator){
      cerr<<"Error:  denominator does not exist!"<<endl;
      return NULL;
    }
    TH1D* result = (TH1D*)numerator->Clone(name);
    Int_t nbins = numerator->GetNbinsX();
    for (Int_t ibin=0; ibin<= nbins+1; ++ibin) {
      Double_t numeratorVal = numerator->GetBinContent(ibin);
      Double_t denominatorVal = denominator->GetBinContent(ibin);
      // Check if the errors are right or the thing is scaled
      Double_t numeratorValErr = numerator->GetBinError(ibin);
      if (!(numeratorValErr==0. || numeratorVal ==0.) ) {
	Double_t rescale = numeratorValErr*numeratorValErr/numeratorVal;
	numeratorVal /= rescale;
      }
      Double_t denominatorValErr = denominator->GetBinError(ibin);
      if (!(denominatorValErr==0. || denominatorVal==0. )) {
	Double_t rescale = denominatorValErr*denominatorValErr/denominatorVal;
	denominatorVal /= rescale;
      }
      Double_t quotient = 0.;
      if (denominatorVal!=0.) {
	quotient = numeratorVal/denominatorVal;
      }
      Double_t quotientErr=0;
      if(numeratorVal>0.0 && denominatorVal>0.0 && (denominatorVal+2.0)>0.0 && (denominatorVal+3.0) >0.0){
	double argument = (numeratorVal+1.0)/(denominatorVal+2.0)*
	  ((numeratorVal+2.0)/(denominatorVal+3.0)-(numeratorVal+1.0)/(denominatorVal+2.0));
	double test = TMath::Power(TMath::Abs(argument),0.5);
	quotientErr = TMath::Sqrt( TMath::Abs(
				  (numeratorVal+1.0)/(denominatorVal+2.0)*
				  ((numeratorVal+2.0)/(denominatorVal+3.0)-(numeratorVal+1.0)/(denominatorVal+2.0))));
      }
      result->SetBinContent(ibin,quotient);
      result->SetBinError(ibin,quotientErr);
      //cout<<"Setting bin "<<ibin<<" to "<<quotient<<" "<<numeratorVal<<"/"<<denominatorVal<<endl;
    }
    return result;
}



TH1D *GetHistoEfficiency(float cut, char *name, int mycase, int color, int marker,bool TPC,bool ITS, int cb, int cblast){
  bool eta = true;
  file->cd();
  char *myname = mynameITS;
  if(TPC && !ITS){ myname = mynameTPCITS;cout<<"I was here I should not be called"<<endl;}
  if(TPC&&ITS) myname = mynameTPCITS;
  char *mynamealt = mynameTPCITS;
  TH2F *numeratorParent; 
  switch(mycase){
  case 0:
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)))->Clone("RecoHadron");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiMinus",i)));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KMinus",i)));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KPlus",i)));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"Proton",i)));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"AntiProton",i)));
	if(!TPC&&ITS){//ITS Standalone tracks - from leftover hits, don't want to double count
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"PiPlus",i)));
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"PiMinus",i)));
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"KMinus",i)));
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"KPlus",i)));
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"Proton",i)));
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"AntiProton",i)));
	}
      }
    }
    else{
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiPlus")))->Clone("RecoHadron");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiMinus")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KMinus")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KPlus")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"Proton")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"AntiProton")));
      if(!TPC&&ITS){//ITS Standalone tracks - from leftover hits, don't want to double count
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"PiPlus")));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"PiMinus")));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"KMinus")));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"KPlus")));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"Proton")));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"AntiProton")));
      }
    }
    break;
  case 1://pion
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)))->Clone("RecoPion");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiMinus",i)));
	if(!TPC&&ITS){//ITS Standalone tracks - from leftover hits, don't want to double count
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"PiMinus",i)));
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"PiPlus",i)));
	}
      }
    }
    else{
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiPlus")))->Clone("RecoPion");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiMinus")));
      if(!TPC&&ITS){//ITS Standalone tracks - from leftover hits, don't want to double count
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"PiPlus")));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"PiMinus")));
      }
    }
    break;
  case 2://kaon
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KPlus",i)))->Clone("RecoKaon");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KPlus",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KMinus",i)));
	if(!TPC&&ITS){//ITS Standalone tracks - from leftover hits, don't want to double count
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"KPlus",i)));
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"KMinus",i)));
	}
      }
    }
    else{
      // cout<<"I am kaoning here !"<<endl;
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KPlus")))->Clone("RecoKaon");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KMinus")));
      if(!TPC&&ITS){//ITS Standalone tracks - from leftover hits, don't want to double count
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"KPlus")));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"KMinus")));
      }
      //cout<<"Numerator "<<numeratorParent->GetEntries()<<endl;
    }
    break;
  case 3://proton
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"Proton",i)))->Clone("RecoProton");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"Proton",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"AntiProton",i)));
	if(!TPC&&ITS){//ITS Standalone tracks - from leftover hits, don't want to double count
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"Proton",i)));
	  numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",mynamealt,"AntiProton",i)));
	}
      }
    }
    else{
      //cout<<"I am protoning here !"<<endl;
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"Proton")))->Clone("RecoProton");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"AntiProton")));
      if(!TPC&&ITS){//ITS Standalone tracks - from leftover hits, don't want to double count
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"Proton")));
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",mynamealt,"AntiProton")));
      }
      //cout<<"Numerator "<<numeratorParent->GetEntries()<<endl;
    }
    break;
  case 4://electron
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s%s",myname,"EPlus",cbname)))->Clone("RecoElectron");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s%s",myname,"EMinus",cbname)));
    if(!TPC&&ITS){//ITS Standalone tracks - from leftover hits, don't want to double count
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s%s",mynamealt,"EPlus",cbname)));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s%s",mynamealt,"EMinus",cbname)));
    }
    break;
  }
  TH2F *denominatorParent; 
  switch(mycase){
  case 0:
    if(cblast != -1){//add more centrality bins
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNSimulatedChargedHadronCB%i",cb)))->Clone("RecoHadron");
      for(int i=cb+1;i<=cblast;i++){
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedChargedHadronCB%i",i)));
      }
    }
    else{
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedChargedHadron"))->Clone("RecoHadron");
    }
    break;
  case 1://pion
    if(cblast != -1){//add more centrality bins
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNSimulatedPiPlusCB%i",cb)))->Clone("RecoPion");
      denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedPiMinusCB%i",cb)));
      for(int i=cb+1;i<=cblast;i++){
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedPiPlusCB%i",i)));
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedPiMinusCB%i",i)));
      }
    }
    else{
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedPiPlus"))->Clone("RecoPion");
      denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedPiMinus"));
    }
    break;
  case 2://kaon
    if(cblast != -1){//add more centrality bins
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNSimulatedKPlusCB%i",cb)))->Clone("RecoKaon");
      denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedKMinusCB%i",cb)));
      for(int i=cb+1;i<=cblast;i++){
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedKPlusCB%i",i)));
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedKMinusCB%i",i)));
      }
    }
    else{
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedKPlus"))->Clone("RecoKaon");
      denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedKMinus"));
    }
    break;
  case 3://proton
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(cb==i)denominatorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNSimulatedProtonCB%i",i)))->Clone("RecoProton");
	else{denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedProtonCB%i",i)));}
	denominatorParent->Add((TH2F*) out2->FindObject(Form("EtNSimulatedAntiProtonCB%i",i)));
      }
    }
    else{
      denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedProton"))->Clone("RecoProton");
      denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedAntiProton"));
    }
    break;
  case 4://electron
    denominatorParent = (TH2F*)((TH2F*) out2->FindObject("EtNSimulatedEPlus"))->Clone("RecoElectron");
    denominatorParent->Add((TH2F*) out2->FindObject("EtNSimulatedEMinus"));
    break;
  }
  TH1D *denominator;
  TH1D *numerator;
  if(eta){
    int lowbin = numeratorParent->GetYaxis()->FindBin(-cut+.001);//make sure we don't accv0entally get the wrong bin
    int highbin = numeratorParent->GetYaxis()->FindBin(cut-.001);
    //cout<<"Projecting from "<<numeratorParent->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<numeratorParent->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;
    denominator = denominatorParent->ProjectionX(Form("garbage%s",name),lowbin,highbin);
    numerator = numeratorParent->ProjectionX(name,lowbin,highbin);
  }
  else{
    int lowbin = denominatorParent->GetXaxis()->FindBin(cut);//make sure we don't accidentally get the wrong bin
    int highbin = denominatorParent->GetXaxis()->GetNbins();
    //cout<<"Projecting from "<<denominatorParent->GetXaxis()->GetBinLowEdge(lowbin)<<" to "<<denominatorParent->GetXaxis()->GetBinLowEdge(highbin+1)<<endl;
    numerator = numeratorParent->ProjectionY(name,lowbin,highbin);
    denominator = denominatorParent->ProjectionY(Form("denominator%s",name),lowbin,highbin);
  }
  delete numeratorParent;
  delete denominatorParent;
  //numerator->Divide(denominator);
  TH1D *result = bayneseffdiv((TH1D*) numerator,(TH1D*)denominator,name);
  //result->Rebin(2);
  //result->Scale(0.5);
  result->SetYTitle("Efficiency");
  result->GetYaxis()->SetTitleOffset(0.8);
  result->GetXaxis()->SetTitleOffset(0.8);
  result->GetYaxis()->SetLabelSize(0.05);
  result->GetXaxis()->SetLabelSize(0.05);
  result->GetYaxis()->SetTitleSize(0.05);
  result->GetXaxis()->SetTitleSize(0.05);
  result->SetMarkerColor(color);
  result->SetLineColor(color);
  result->SetMarkerStyle(marker);
  //result->SetName(name);
  //result->Draw("e");
  delete denominator;
  delete numerator;
  delete numeratorParent;
  delete denominatorParent;
  result->SetName(name);
  return result;

}

void CorrEfficiencyPlots(bool TPC, char *prodname, char *shortprodname){
  bool ITS = true;
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",600,400);
  c->SetTopMargin(0.02);
  c->SetRightMargin(0.02);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);
  //c->SetLogx();

  int colortotal = 1;
  int colorpi = 2;
  int colork = 3;
  int colorp = 4;
  int phosmarker = 20;
  int emcalmarker = 24;
  float ptcut1 = 0.05;
  float ptcut2 = 0.1;
  TH1D *PHOStotal = GetHistoEfficiency(0.12,"PHOStotal",0,colortotal,phosmarker,TPC,ITS);
  TH1D *PHOSpi = GetHistoEfficiency(0.12,"PHOSpi",1,colorpi,phosmarker,TPC,ITS);
  TH1D *PHOSp = GetHistoEfficiency(0.12,"PHOSp",2,colork,phosmarker,TPC,ITS);
  TH1D *PHOSk = GetHistoEfficiency(0.12,"PHOSk",3,colorp,phosmarker,TPC,ITS);
  if(!TPC){PHOStotal->GetXaxis()->SetRange(PHOStotal->GetXaxis()->FindBin(0.05),PHOStotal->GetXaxis()->FindBin(1.0));}
  else{PHOStotal->GetXaxis()->SetRange(PHOStotal->GetXaxis()->FindBin(0.15),PHOStotal->GetXaxis()->FindBin(3.0));}
  PHOStotal->SetMinimum(0.0);
  PHOStotal->SetMaximum(1.0);
  PHOStotal->Draw();
  PHOSpi->Draw("same");
  PHOSp->Draw("same");
  PHOSk->Draw("same");
  TH1D *EMCALtotal = GetHistoEfficiency(0.7,"EMCALtotal",0,colortotal,emcalmarker,TPC,ITS);
  TH1D *EMCALpi = GetHistoEfficiency(0.7,"EMCALpi",1,colorpi,emcalmarker,TPC,ITS);
  TH1D *EMCALp = GetHistoEfficiency(0.7,"EMCALp",2,colork,emcalmarker,TPC,ITS);
  TH1D *EMCALk = GetHistoEfficiency(0.7,"EMCALk",3,colorp,emcalmarker,TPC,ITS);
  EMCALtotal->Draw("same");
  EMCALpi->Draw("same");
  EMCALp->Draw("same");
  EMCALk->Draw("same");


  TLegend *leg = new  TLegend(0.22651,0.247312,0.370805,0.438172);
  leg->AddEntry(PHOStotal,"#pi,K,p");
  leg->AddEntry(PHOSpi,"#pi^{#pm}");
  leg->AddEntry(PHOSk,"K^{#pm}");
  leg->AddEntry(PHOSp,"p,#bar{p}");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
 leg->Draw();

  TLine *line = new TLine(0.2,0.0,0.2,1.0);
  line->Draw();
  line->SetLineWidth(3.0);
  //line->SetLineColor(TColor::kYellow);
  line->SetLineStyle(2);
  TLatex *tex = new TLatex(0.497269,0.0513196,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLatex *tex3 = new TLatex(1.16186,0.28348,"Closed symbols |#eta|<0.12 (PHOS)");
  tex3->SetTextSize(0.0537634);
  tex3->Draw();
  TLatex *tex4 = new TLatex(1.16186,0.213221,"Open symbols |#eta|<0.70 (EMCal)");
  tex4->SetTextSize(0.0537634);
  tex4->Draw();
  TLatex *tex2 = new TLatex(0.241937,0.448436,"Likely TPC cut-off 200 MeV/c");
  tex2->SetTextSize(0.0537634);
  tex2->Draw();
  if(TPC){
    if(ITS){
      sprintf(epsname,"pics/%s/CorrEfficiencyITSTPC.eps",shortprodname);
      sprintf(pngname,"pics/%s/CorrEfficiencyITSTPC.png",shortprodname);
    }
    else{
      sprintf(epsname,"pics/%s/CorrEfficiencyTPC.eps",shortprodname);
      sprintf(pngname,"pics/%s/CorrEfficiencyTPC.png",shortprodname);
    }
  }
  else{
    sprintf(epsname,"pics/%s/CorrEfficiencyITS.eps",shortprodname);
    sprintf(pngname,"pics/%s/CorrEfficiencyITS.png",shortprodname);
  }
  delete PHOStotal;
  delete PHOSpi;
  delete PHOSp;
  delete PHOSk;
  delete EMCALtotal;
  delete EMCALpi;
  delete EMCALp;
  delete EMCALk;
  delete leg;
  c->SaveAs(epsname);
  c->SaveAs(pngname);
  delete c;
  delete line;
  delete tex;
  delete tex3;
}

//==================================CorrBkgd=================================================
TH1D *GetHistoCorrBkgd(float etacut,char *name, bool TPC,bool ispp,bool forSim){
  file->cd();
  char *reweightname = reweightedNo;
  if(!forSim) reweightname = reweightedYes;
  char *myname = mynameITS;
  if(TPC) myname = mynameTPCITS;
  TH2F *signal = ((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedPiPlus",myname)))->Clone("signal");
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedPiMinus",myname)));
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedKMinus",myname)));
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedKPlus",myname)));
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedProton",myname)));
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sIdentifiedAntiProton",myname)));
  signal->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sUnidentifiedAssumingPion",myname)));

  //Et of all unidentified hadrons (plus hadrons identified as pions) calculated assuming their true mass
  TH2F *bkgd = ((TH2F*) out2->FindObject(Form("EtReconstructed%sMisidentifiedElectrons",myname)))->Clone("bkgd");
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sLambdaDaughters%s",myname,reweightname)));
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiLambdaDaughters%s",myname,reweightname)));
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sK0SDaughters%s",myname,reweightname)));
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sXiDaughters",myname)));
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiXiDaughters",myname)));
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sOmegaDaughters",myname)));
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sAntiOmegaDaughters",myname)));
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sConversionElectrons",myname)) );
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sSecondaryMuons",myname)) );
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sSecondaryPions",myname)) );
  bkgd->Add((TH2F*) out2->FindObject(Form("EtReconstructed%sSecondaryProtons",myname)) );
  int lowbin = bkgd->GetYaxis()->FindBin(-etacut+.001);//make sure we don't accidentally get the wrong bin
  int highbin = bkgd->GetYaxis()->FindBin(etacut-.001);
  //cout<<"Projecting from "<<bkgd->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<bkgd->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;


  TH1D *denominator = signal->ProjectionX(name,lowbin,highbin);
  TH1D *numerator = bkgd->ProjectionX("numerator",lowbin,highbin);
  denominator->Add(numerator);
  numerator->Divide(denominator);
  numerator->SetYTitle("Ratio of E_{T}^{background}/E_{T}^{real}");
  numerator->GetYaxis()->SetTitleOffset(1.2);
  delete signal;
  delete bkgd;
  delete denominator;
  numerator->SetName(name);
  return numerator;

}

void CorrBkgdPlots(char *prodname, char *shortprodname, bool TPC,bool ispp,bool forSim){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c = new TCanvas("c","c",500,400);
  c->SetTopMargin(0.04);
  c->SetRightMargin(0.04);
  c->SetBorderSize(0);
  c->SetFillColor(0);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetFrameFillColor(0);
  c->SetFrameBorderMode(0);

  TH1D *PHOS = GetHistoCorrBkgd(0.12,"PHOS2",TPC,ispp,forSim);
  TH1D *EMCAL = GetHistoCorrBkgd(0.7,"EMCAL2",TPC,ispp,forSim);
  PHOS->SetMarkerColor(2);
  EMCAL->SetMarkerColor(4);
  PHOS->SetLineColor(2);
  EMCAL->SetLineColor(4);
  //EMCAL->SetLineWidth(2);
  //PHOS->SetAxisRange(0.0,4);
  PHOS->SetMaximum(0.2);
  PHOS->SetMinimum(0.0);
  PHOS->SetMarkerStyle(20);;
  EMCAL->SetMarkerStyle(21);
  //  TF1 *funcEMCAL = new TF1("funcEMCAL","[0]+0.0*x",0.05,4);
//   funcEMCAL->SetParameter(0,0.95);
//   funcEMCAL->SetParLimits(0,0.9,1.1);
  //EMCAL->Fit(funcEMCAL);//,"","",0.05,3.0);
//   TF1 *funcPHOS = new TF1("funcPHOS","[0]+0.0*x",0.05,4);
//   funcPHOS->SetParameter(0,1.0);
  //PHOS->Fit(funcPHOS);
  if(TPC)  PHOS->GetXaxis()->SetRange(PHOS->GetXaxis()->FindBin(0.0),PHOS->GetXaxis()->FindBin(4.));
  else{  PHOS->GetXaxis()->SetRange(PHOS->GetXaxis()->FindBin(0.0),PHOS->GetXaxis()->FindBin(1.));}
  PHOS->Draw();
  EMCAL->Draw("same");
  TLatex *tex = new TLatex(0.161478,1.0835,prodname);
  tex->SetTextSize(0.0537634);
  tex->Draw();
  TLegend *leg = new TLegend(0.145161,0.604839,0.40121,0.860215);
  leg->AddEntry(PHOS,"|#eta|<0.12");
  leg->AddEntry(EMCAL,"|#eta|<0.70");
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  if(TPC){
    c->SaveAs(Form("pics/%s/bkgdTPC.eps",shortprodname));
    c->SaveAs(Form("pics/%s/bkgdTPC.png",shortprodname));
  }
  else{
    c->SaveAs(Form("pics/%s/bkgdITS.eps",shortprodname));
    c->SaveAs(Form("pics/%s/bkgdITS.png",shortprodname));
  }
  delete c;
  delete PHOS;
  delete EMCAL;
  delete tex;
  delete leg;

}

//TH1D *ppBkgd()
TH1D *pp276TPCBkgd(){
//=========Macro generated from canvas: c4/c4
//=========  (Thu May  5 23:07:14 2011) by ROOT version5.28/00c
   TH1D *Allpt7TeVScaling = new TH1D("Allpt7TeVScaling","Reconstructed E_{T} from misidentified electrons",111, xAxis1);
   Allpt7TeVScaling->SetBinContent(13,0.1673787);
   Allpt7TeVScaling->SetBinContent(14,0.1548408);
   Allpt7TeVScaling->SetBinContent(15,0.1328529);
   Allpt7TeVScaling->SetBinContent(16,0.1220669);
   Allpt7TeVScaling->SetBinContent(17,0.1109226);
   Allpt7TeVScaling->SetBinContent(18,0.0981996);
   Allpt7TeVScaling->SetBinContent(19,0.09152354);
   Allpt7TeVScaling->SetBinContent(20,0.07937801);
   Allpt7TeVScaling->SetBinContent(21,0.07163593);
   Allpt7TeVScaling->SetBinContent(22,0.0562236);
   Allpt7TeVScaling->SetBinContent(23,0.04398756);
   Allpt7TeVScaling->SetBinContent(24,0.039459);
   Allpt7TeVScaling->SetBinContent(25,0.03618901);
   Allpt7TeVScaling->SetBinContent(26,0.0322147);
   Allpt7TeVScaling->SetBinContent(27,0.03294796);
   Allpt7TeVScaling->SetBinContent(28,0.03114461);
   Allpt7TeVScaling->SetBinContent(29,0.03039579);
   Allpt7TeVScaling->SetBinContent(30,0.03261539);
   Allpt7TeVScaling->SetBinContent(31,0.03239832);
   Allpt7TeVScaling->SetBinContent(32,0.03055859);
   Allpt7TeVScaling->SetBinContent(33,0.03181271);
   Allpt7TeVScaling->SetBinContent(34,0.02882665);
   Allpt7TeVScaling->SetBinContent(35,0.03102149);
   Allpt7TeVScaling->SetBinContent(36,0.0272696);
   Allpt7TeVScaling->SetBinContent(37,0.02783727);
   Allpt7TeVScaling->SetBinContent(38,0.02334413);
   Allpt7TeVScaling->SetBinContent(39,0.02357657);
   Allpt7TeVScaling->SetBinContent(40,0.02630686);
   Allpt7TeVScaling->SetBinContent(41,0.02372341);
   Allpt7TeVScaling->SetBinContent(42,0.02279457);
   Allpt7TeVScaling->SetBinContent(43,0.01999143);
   Allpt7TeVScaling->SetBinContent(44,0.02019579);
   Allpt7TeVScaling->SetBinContent(45,0.02029719);
   Allpt7TeVScaling->SetBinContent(46,0.0166231);
   Allpt7TeVScaling->SetBinContent(47,0.01824842);
   Allpt7TeVScaling->SetBinContent(48,0.0188995);
   Allpt7TeVScaling->SetBinContent(49,0.01836067);
   Allpt7TeVScaling->SetBinContent(50,0.02113978);
   Allpt7TeVScaling->SetBinContent(51,0.02236854);
   Allpt7TeVScaling->SetBinContent(52,0.0204874);
   Allpt7TeVScaling->SetBinContent(53,0.02236788);
   Allpt7TeVScaling->SetBinContent(54,0.01862457);
   Allpt7TeVScaling->SetBinContent(55,0.0215448);
   Allpt7TeVScaling->SetBinContent(56,0.01757402);
   Allpt7TeVScaling->SetBinContent(57,0.01899477);
   Allpt7TeVScaling->SetBinContent(58,0.01766989);
   Allpt7TeVScaling->SetBinContent(59,0.01984435);
   Allpt7TeVScaling->SetBinContent(60,0.01842092);
   Allpt7TeVScaling->SetBinContent(61,0.01668053);
   Allpt7TeVScaling->SetBinContent(62,0.01708603);
   Allpt7TeVScaling->SetBinContent(63,0.01557306);
   Allpt7TeVScaling->SetBinContent(64,0.01878942);
   Allpt7TeVScaling->SetBinContent(65,0.01458043);
   Allpt7TeVScaling->SetBinContent(66,0.0151346);
   Allpt7TeVScaling->SetBinContent(67,0.01565407);
   Allpt7TeVScaling->SetBinContent(68,0.01327014);
   Allpt7TeVScaling->SetBinContent(69,0.01284378);
   Allpt7TeVScaling->SetBinContent(70,0.01046102);
   Allpt7TeVScaling->SetBinContent(71,0.009963117);
   Allpt7TeVScaling->SetBinContent(72,0.01359609);
   Allpt7TeVScaling->SetBinContent(73,0.01088263);
   Allpt7TeVScaling->SetBinContent(74,0.01317689);
   Allpt7TeVScaling->SetBinContent(75,0.01433201);
   Allpt7TeVScaling->SetBinContent(76,0.01215939);
   Allpt7TeVScaling->SetBinContent(77,0.01380766);
   Allpt7TeVScaling->SetBinContent(78,0.01499505);
   Allpt7TeVScaling->SetBinContent(79,0.015013);
   Allpt7TeVScaling->SetBinContent(80,0.01140629);
   Allpt7TeVScaling->SetBinContent(81,0.009461979);
   Allpt7TeVScaling->SetBinContent(82,0.007775343);
   Allpt7TeVScaling->SetBinContent(83,0.01092198);
   Allpt7TeVScaling->SetBinContent(84,0.009924767);
   Allpt7TeVScaling->SetBinContent(85,0.006399066);
   Allpt7TeVScaling->SetBinContent(86,0.008370865);
   Allpt7TeVScaling->SetBinContent(87,0.01115192);
   Allpt7TeVScaling->SetBinContent(88,0.003781796);
   Allpt7TeVScaling->SetBinContent(89,0.01596916);
   Allpt7TeVScaling->SetBinContent(90,0.005959723);
   Allpt7TeVScaling->SetBinContent(91,0.006885685);
   Allpt7TeVScaling->SetBinContent(93,0.006011974);
   Allpt7TeVScaling->SetBinContent(94,0.008862645);
   Allpt7TeVScaling->SetBinContent(95,0.009963223);
   Allpt7TeVScaling->SetBinContent(96,0.01457121);
   Allpt7TeVScaling->SetBinContent(97,0.007500135);
   Allpt7TeVScaling->SetBinContent(98,0.02636933);
   Allpt7TeVScaling->SetBinContent(103,0.06353767);
   Allpt7TeVScaling->SetBinError(13,0.02859277);
   Allpt7TeVScaling->SetBinError(14,0.006796474);
   Allpt7TeVScaling->SetBinError(15,0.004234415);
   Allpt7TeVScaling->SetBinError(16,0.0036414);
   Allpt7TeVScaling->SetBinError(17,0.003333835);
   Allpt7TeVScaling->SetBinError(18,0.003009034);
   Allpt7TeVScaling->SetBinError(19,0.002839933);
   Allpt7TeVScaling->SetBinError(20,0.002588262);
   Allpt7TeVScaling->SetBinError(21,0.001703954);
   Allpt7TeVScaling->SetBinError(22,0.001495734);
   Allpt7TeVScaling->SetBinError(23,0.001297475);
   Allpt7TeVScaling->SetBinError(24,0.001213434);
   Allpt7TeVScaling->SetBinError(25,0.001162398);
   Allpt7TeVScaling->SetBinError(26,0.001137158);
   Allpt7TeVScaling->SetBinError(27,0.001208358);
   Allpt7TeVScaling->SetBinError(28,0.001183675);
   Allpt7TeVScaling->SetBinError(29,0.001263691);
   Allpt7TeVScaling->SetBinError(30,0.001320041);
   Allpt7TeVScaling->SetBinError(31,0.001362399);
   Allpt7TeVScaling->SetBinError(32,0.001338259);
   Allpt7TeVScaling->SetBinError(33,0.001369813);
   Allpt7TeVScaling->SetBinError(34,0.001358652);
   Allpt7TeVScaling->SetBinError(35,0.00142864);
   Allpt7TeVScaling->SetBinError(36,0.001378454);
   Allpt7TeVScaling->SetBinError(37,0.001384729);
   Allpt7TeVScaling->SetBinError(38,0.001257903);
   Allpt7TeVScaling->SetBinError(39,0.001295459);
   Allpt7TeVScaling->SetBinError(40,0.001546992);
   Allpt7TeVScaling->SetBinError(41,0.001496938);
   Allpt7TeVScaling->SetBinError(42,0.001486813);
   Allpt7TeVScaling->SetBinError(43,0.00140119);
   Allpt7TeVScaling->SetBinError(44,0.001478331);
   Allpt7TeVScaling->SetBinError(45,0.001566594);
   Allpt7TeVScaling->SetBinError(46,0.001407568);
   Allpt7TeVScaling->SetBinError(47,0.001509742);
   Allpt7TeVScaling->SetBinError(48,0.001583129);
   Allpt7TeVScaling->SetBinError(49,0.001668789);
   Allpt7TeVScaling->SetBinError(50,0.001915902);
   Allpt7TeVScaling->SetBinError(51,0.001968585);
   Allpt7TeVScaling->SetBinError(52,0.00195473);
   Allpt7TeVScaling->SetBinError(53,0.00207896);
   Allpt7TeVScaling->SetBinError(54,0.001937327);
   Allpt7TeVScaling->SetBinError(55,0.002108783);
   Allpt7TeVScaling->SetBinError(56,0.001816997);
   Allpt7TeVScaling->SetBinError(57,0.001892795);
   Allpt7TeVScaling->SetBinError(58,0.001803616);
   Allpt7TeVScaling->SetBinError(59,0.001994761);
   Allpt7TeVScaling->SetBinError(60,0.001964909);
   Allpt7TeVScaling->SetBinError(61,0.001255034);
   Allpt7TeVScaling->SetBinError(62,0.001312232);
   Allpt7TeVScaling->SetBinError(63,0.001358296);
   Allpt7TeVScaling->SetBinError(64,0.001569736);
   Allpt7TeVScaling->SetBinError(65,0.001469072);
   Allpt7TeVScaling->SetBinError(66,0.001588803);
   Allpt7TeVScaling->SetBinError(67,0.001683862);
   Allpt7TeVScaling->SetBinError(68,0.00164699);
   Allpt7TeVScaling->SetBinError(69,0.001760645);
   Allpt7TeVScaling->SetBinError(70,0.00160423);
   Allpt7TeVScaling->SetBinError(71,0.001693146);
   Allpt7TeVScaling->SetBinError(72,0.002042705);
   Allpt7TeVScaling->SetBinError(73,0.001965856);
   Allpt7TeVScaling->SetBinError(74,0.002245841);
   Allpt7TeVScaling->SetBinError(75,0.002520926);
   Allpt7TeVScaling->SetBinError(76,0.002402524);
   Allpt7TeVScaling->SetBinError(77,0.002627656);
   Allpt7TeVScaling->SetBinError(78,0.002963027);
   Allpt7TeVScaling->SetBinError(79,0.003088673);
   Allpt7TeVScaling->SetBinError(80,0.002958055);
   Allpt7TeVScaling->SetBinError(81,0.001434729);
   Allpt7TeVScaling->SetBinError(82,0.001531537);
   Allpt7TeVScaling->SetBinError(83,0.002200666);
   Allpt7TeVScaling->SetBinError(84,0.00242244);
   Allpt7TeVScaling->SetBinError(85,0.002270338);
   Allpt7TeVScaling->SetBinError(86,0.002972849);
   Allpt7TeVScaling->SetBinError(87,0.003966084);
   Allpt7TeVScaling->SetBinError(88,0.002679263);
   Allpt7TeVScaling->SetBinError(89,0.006083804);
   Allpt7TeVScaling->SetBinError(90,0.004227019);
   Allpt7TeVScaling->SetBinError(91,0.004886052);
   Allpt7TeVScaling->SetBinError(93,0.006030514);
   Allpt7TeVScaling->SetBinError(94,0.008892003);
   Allpt7TeVScaling->SetBinError(95,0.01001297);
   Allpt7TeVScaling->SetBinError(96,0.008480241);
   Allpt7TeVScaling->SetBinError(97,0.007527505);
   Allpt7TeVScaling->SetBinError(98,0.01889477);
   Allpt7TeVScaling->SetBinError(103,0.065502);
   Allpt7TeVScaling->SetMinimum(0);
   Allpt7TeVScaling->SetMaximum(0.04);
   Allpt7TeVScaling->SetEntries(1049.145);
   Allpt7TeVScaling->SetStats(0);
   Allpt7TeVScaling->SetMarkerStyle(20);
   Allpt7TeVScaling->GetXaxis()->SetTitle("p_{T}");
   Allpt7TeVScaling->GetXaxis()->SetRange(1,91);
   Allpt7TeVScaling->GetYaxis()->SetTitle("Ratio of E_{T}^{background}/E_{T}^{had, meas.}");
   Allpt7TeVScaling->GetYaxis()->SetTitleOffset(1.2);
   TH1D *Allpt900GeVScaling = new TH1D("Allpt900GeVScaling","Reconstructed E_{T} from misidentified electrons",111, xAxis1);
   Allpt900GeVScaling->SetBinContent(13,0.1950247);
   Allpt900GeVScaling->SetBinContent(14,0.147561);
   Allpt900GeVScaling->SetBinContent(15,0.136378);
   Allpt900GeVScaling->SetBinContent(16,0.1233744);
   Allpt900GeVScaling->SetBinContent(17,0.1126807);
   Allpt900GeVScaling->SetBinContent(18,0.09909555);
   Allpt900GeVScaling->SetBinContent(19,0.09191239);
   Allpt900GeVScaling->SetBinContent(20,0.08438053);
   Allpt900GeVScaling->SetBinContent(21,0.07119319);
   Allpt900GeVScaling->SetBinContent(22,0.05773434);
   Allpt900GeVScaling->SetBinContent(23,0.04292845);
   Allpt900GeVScaling->SetBinContent(24,0.0396127);
   Allpt900GeVScaling->SetBinContent(25,0.03395077);
   Allpt900GeVScaling->SetBinContent(26,0.03305251);
   Allpt900GeVScaling->SetBinContent(27,0.03297093);
   Allpt900GeVScaling->SetBinContent(28,0.031215);
   Allpt900GeVScaling->SetBinContent(29,0.03178372);
   Allpt900GeVScaling->SetBinContent(30,0.03236985);
   Allpt900GeVScaling->SetBinContent(31,0.03271276);
   Allpt900GeVScaling->SetBinContent(32,0.02993124);
   Allpt900GeVScaling->SetBinContent(33,0.03004746);
   Allpt900GeVScaling->SetBinContent(34,0.03077143);
   Allpt900GeVScaling->SetBinContent(35,0.03211214);
   Allpt900GeVScaling->SetBinContent(36,0.02807512);
   Allpt900GeVScaling->SetBinContent(37,0.02651954);
   Allpt900GeVScaling->SetBinContent(38,0.02503885);
   Allpt900GeVScaling->SetBinContent(39,0.02325777);
   Allpt900GeVScaling->SetBinContent(40,0.02773712);
   Allpt900GeVScaling->SetBinContent(41,0.02491155);
   Allpt900GeVScaling->SetBinContent(42,0.02389028);
   Allpt900GeVScaling->SetBinContent(43,0.02085532);
   Allpt900GeVScaling->SetBinContent(44,0.01959721);
   Allpt900GeVScaling->SetBinContent(45,0.01937255);
   Allpt900GeVScaling->SetBinContent(46,0.01745646);
   Allpt900GeVScaling->SetBinContent(47,0.01890853);
   Allpt900GeVScaling->SetBinContent(48,0.02053692);
   Allpt900GeVScaling->SetBinContent(49,0.01813791);
   Allpt900GeVScaling->SetBinContent(50,0.02187495);
   Allpt900GeVScaling->SetBinContent(51,0.02409687);
   Allpt900GeVScaling->SetBinContent(52,0.0206565);
   Allpt900GeVScaling->SetBinContent(53,0.02269165);
   Allpt900GeVScaling->SetBinContent(54,0.01850795);
   Allpt900GeVScaling->SetBinContent(55,0.02115427);
   Allpt900GeVScaling->SetBinContent(56,0.01543302);
   Allpt900GeVScaling->SetBinContent(57,0.01642095);
   Allpt900GeVScaling->SetBinContent(58,0.01542331);
   Allpt900GeVScaling->SetBinContent(59,0.01931254);
   Allpt900GeVScaling->SetBinContent(60,0.01650654);
   Allpt900GeVScaling->SetBinContent(61,0.01827512);
   Allpt900GeVScaling->SetBinContent(62,0.01700125);
   Allpt900GeVScaling->SetBinContent(63,0.01544321);
   Allpt900GeVScaling->SetBinContent(64,0.01830588);
   Allpt900GeVScaling->SetBinContent(65,0.01584369);
   Allpt900GeVScaling->SetBinContent(66,0.0145421);
   Allpt900GeVScaling->SetBinContent(67,0.01382206);
   Allpt900GeVScaling->SetBinContent(68,0.01358825);
   Allpt900GeVScaling->SetBinContent(69,0.01294888);
   Allpt900GeVScaling->SetBinContent(70,0.009346655);
   Allpt900GeVScaling->SetBinContent(71,0.01235781);
   Allpt900GeVScaling->SetBinContent(72,0.0127376);
   Allpt900GeVScaling->SetBinContent(73,0.008968166);
   Allpt900GeVScaling->SetBinContent(74,0.01467765);
   Allpt900GeVScaling->SetBinContent(75,0.01432309);
   Allpt900GeVScaling->SetBinContent(76,0.01204856);
   Allpt900GeVScaling->SetBinContent(77,0.01436071);
   Allpt900GeVScaling->SetBinContent(78,0.01451222);
   Allpt900GeVScaling->SetBinContent(79,0.01378967);
   Allpt900GeVScaling->SetBinContent(80,0.01040148);
   Allpt900GeVScaling->SetBinContent(81,0.01041213);
   Allpt900GeVScaling->SetBinContent(82,0.009971816);
   Allpt900GeVScaling->SetBinContent(83,0.01417668);
   Allpt900GeVScaling->SetBinContent(84,0.008713292);
   Allpt900GeVScaling->SetBinContent(85,0.005425386);
   Allpt900GeVScaling->SetBinContent(86,0.005303683);
   Allpt900GeVScaling->SetBinContent(87,0.006759432);
   Allpt900GeVScaling->SetBinContent(88,0.004445477);
   Allpt900GeVScaling->SetBinContent(89,0.01810702);
   Allpt900GeVScaling->SetBinContent(90,0.0108909);
   Allpt900GeVScaling->SetBinContent(91,0.01186746);
   Allpt900GeVScaling->SetBinContent(92,0.006128581);
   Allpt900GeVScaling->SetBinContent(93,0.007452721);
   Allpt900GeVScaling->SetBinContent(94,0.0259471);
   Allpt900GeVScaling->SetBinContent(95,0.0114321);
   Allpt900GeVScaling->SetBinContent(98,0.03264295);
   Allpt900GeVScaling->SetBinContent(103,0.07829265);
   Allpt900GeVScaling->SetBinError(13,0.03504783);
   Allpt900GeVScaling->SetBinError(14,0.007336667);
   Allpt900GeVScaling->SetBinError(15,0.004765958);
   Allpt900GeVScaling->SetBinError(16,0.004028143);
   Allpt900GeVScaling->SetBinError(17,0.003760166);
   Allpt900GeVScaling->SetBinError(18,0.003357858);
   Allpt900GeVScaling->SetBinError(19,0.003156502);
   Allpt900GeVScaling->SetBinError(20,0.002970822);
   Allpt900GeVScaling->SetBinError(21,0.001886555);
   Allpt900GeVScaling->SetBinError(22,0.001676045);
   Allpt900GeVScaling->SetBinError(23,0.001411794);
   Allpt900GeVScaling->SetBinError(24,0.001346425);
   Allpt900GeVScaling->SetBinError(25,0.001254258);
   Allpt900GeVScaling->SetBinError(26,0.001299384);
   Allpt900GeVScaling->SetBinError(27,0.001328066);
   Allpt900GeVScaling->SetBinError(28,0.001330088);
   Allpt900GeVScaling->SetBinError(29,0.001446767);
   Allpt900GeVScaling->SetBinError(30,0.00145402);
   Allpt900GeVScaling->SetBinError(31,0.001463994);
   Allpt900GeVScaling->SetBinError(32,0.001449209);
   Allpt900GeVScaling->SetBinError(33,0.001414847);
   Allpt900GeVScaling->SetBinError(34,0.001598429);
   Allpt900GeVScaling->SetBinError(35,0.001590854);
   Allpt900GeVScaling->SetBinError(36,0.001526798);
   Allpt900GeVScaling->SetBinError(37,0.001486126);
   Allpt900GeVScaling->SetBinError(38,0.001457142);
   Allpt900GeVScaling->SetBinError(39,0.001474238);
   Allpt900GeVScaling->SetBinError(40,0.001757902);
   Allpt900GeVScaling->SetBinError(41,0.001712601);
   Allpt900GeVScaling->SetBinError(42,0.001671385);
   Allpt900GeVScaling->SetBinError(43,0.001570556);
   Allpt900GeVScaling->SetBinError(44,0.001627725);
   Allpt900GeVScaling->SetBinError(45,0.001637262);
   Allpt900GeVScaling->SetBinError(46,0.001641227);
   Allpt900GeVScaling->SetBinError(47,0.001716345);
   Allpt900GeVScaling->SetBinError(48,0.001907778);
   Allpt900GeVScaling->SetBinError(49,0.001842589);
   Allpt900GeVScaling->SetBinError(50,0.00215445);
   Allpt900GeVScaling->SetBinError(51,0.002310602);
   Allpt900GeVScaling->SetBinError(52,0.002171046);
   Allpt900GeVScaling->SetBinError(53,0.002386521);
   Allpt900GeVScaling->SetBinError(54,0.002149237);
   Allpt900GeVScaling->SetBinError(55,0.002233898);
   Allpt900GeVScaling->SetBinError(56,0.001856697);
   Allpt900GeVScaling->SetBinError(57,0.001891742);
   Allpt900GeVScaling->SetBinError(58,0.001840334);
   Allpt900GeVScaling->SetBinError(59,0.002168741);
   Allpt900GeVScaling->SetBinError(60,0.002044863);
   Allpt900GeVScaling->SetBinError(61,0.001431972);
   Allpt900GeVScaling->SetBinError(62,0.001452706);
   Allpt900GeVScaling->SetBinError(63,0.001477499);
   Allpt900GeVScaling->SetBinError(64,0.001722091);
   Allpt900GeVScaling->SetBinError(65,0.001674803);
   Allpt900GeVScaling->SetBinError(66,0.001737491);
   Allpt900GeVScaling->SetBinError(67,0.001783912);
   Allpt900GeVScaling->SetBinError(68,0.00186443);
   Allpt900GeVScaling->SetBinError(69,0.001965873);
   Allpt900GeVScaling->SetBinError(70,0.001687347);
   Allpt900GeVScaling->SetBinError(71,0.002072538);
   Allpt900GeVScaling->SetBinError(72,0.002233273);
   Allpt900GeVScaling->SetBinError(73,0.001965905);
   Allpt900GeVScaling->SetBinError(74,0.002669571);
   Allpt900GeVScaling->SetBinError(75,0.00278613);
   Allpt900GeVScaling->SetBinError(76,0.002661412);
   Allpt900GeVScaling->SetBinError(77,0.002952605);
   Allpt900GeVScaling->SetBinError(78,0.003268721);
   Allpt900GeVScaling->SetBinError(79,0.003273694);
   Allpt900GeVScaling->SetBinError(80,0.002900898);
   Allpt900GeVScaling->SetBinError(81,0.001677933);
   Allpt900GeVScaling->SetBinError(82,0.001933408);
   Allpt900GeVScaling->SetBinError(83,0.002759328);
   Allpt900GeVScaling->SetBinError(84,0.002526604);
   Allpt900GeVScaling->SetBinError(85,0.002297818);
   Allpt900GeVScaling->SetBinError(86,0.002659213);
   Allpt900GeVScaling->SetBinError(87,0.003392141);
   Allpt900GeVScaling->SetBinError(88,0.003150511);
   Allpt900GeVScaling->SetBinError(89,0.007499872);
   Allpt900GeVScaling->SetBinError(90,0.006322856);
   Allpt900GeVScaling->SetBinError(91,0.00689283);
   Allpt900GeVScaling->SetBinError(92,0.006146806);
   Allpt900GeVScaling->SetBinError(93,0.007481131);
   Allpt900GeVScaling->SetBinError(94,0.01853139);
   Allpt900GeVScaling->SetBinError(95,0.01149767);
   Allpt900GeVScaling->SetBinError(98,0.02345991);
   Allpt900GeVScaling->SetBinError(103,0.08126059);
   Allpt900GeVScaling->SetMinimum(0);
   Allpt900GeVScaling->SetMaximum(0.04);
   Allpt900GeVScaling->SetEntries(724.3671);
   Allpt900GeVScaling->SetStats(0);
   Allpt900GeVScaling->SetMarkerStyle(20);
   Allpt900GeVScaling->GetXaxis()->SetTitle("p_{T}");
   Allpt900GeVScaling->GetXaxis()->SetRange(1,91);
   Allpt900GeVScaling->GetYaxis()->SetTitle("Ratio of E_{T}^{background}/E_{T}^{had, meas.}");
   Allpt900GeVScaling->GetYaxis()->SetTitleOffset(1.2);
   Allpt7TeVScaling->Add(Allpt900GeVScaling);
   Allpt7TeVScaling->Scale(0.5);
   delete Allpt900GeVScaling;
   return Allpt7TeVScaling;
}

TH1D *PbPb276TPCBkgd(){
   TH1D *Allpt7TeVScaling = new TH1D("Allpt7TeVScaling","Reconstructed E_{T} from misidentified electrons",111, xAxis1);
   Allpt7TeVScaling->SetBinContent(13,0.1745755);
   Allpt7TeVScaling->SetBinContent(14,0.1624186);
   Allpt7TeVScaling->SetBinContent(15,0.1431048);
   Allpt7TeVScaling->SetBinContent(16,0.1331614);
   Allpt7TeVScaling->SetBinContent(17,0.121262);
   Allpt7TeVScaling->SetBinContent(18,0.1082923);
   Allpt7TeVScaling->SetBinContent(19,0.1011451);
   Allpt7TeVScaling->SetBinContent(20,0.08875923);
   Allpt7TeVScaling->SetBinContent(21,0.08100933);
   Allpt7TeVScaling->SetBinContent(22,0.06527292);
   Allpt7TeVScaling->SetBinContent(23,0.05220065);
   Allpt7TeVScaling->SetBinContent(24,0.04752904);
   Allpt7TeVScaling->SetBinContent(25,0.04421568);
   Allpt7TeVScaling->SetBinContent(26,0.03879551);
   Allpt7TeVScaling->SetBinContent(27,0.03944679);
   Allpt7TeVScaling->SetBinContent(28,0.03707215);
   Allpt7TeVScaling->SetBinContent(29,0.03663599);
   Allpt7TeVScaling->SetBinContent(30,0.03937437);
   Allpt7TeVScaling->SetBinContent(31,0.03888919);
   Allpt7TeVScaling->SetBinContent(32,0.03646502);
   Allpt7TeVScaling->SetBinContent(33,0.03784992);
   Allpt7TeVScaling->SetBinContent(34,0.03360079);
   Allpt7TeVScaling->SetBinContent(35,0.03761566);
   Allpt7TeVScaling->SetBinContent(36,0.03276653);
   Allpt7TeVScaling->SetBinContent(37,0.03401566);
   Allpt7TeVScaling->SetBinContent(38,0.02803893);
   Allpt7TeVScaling->SetBinContent(39,0.0290018);
   Allpt7TeVScaling->SetBinContent(40,0.03323201);
   Allpt7TeVScaling->SetBinContent(41,0.03006077);
   Allpt7TeVScaling->SetBinContent(42,0.02867674);
   Allpt7TeVScaling->SetBinContent(43,0.02471421);
   Allpt7TeVScaling->SetBinContent(44,0.02610375);
   Allpt7TeVScaling->SetBinContent(45,0.02601621);
   Allpt7TeVScaling->SetBinContent(46,0.02184781);
   Allpt7TeVScaling->SetBinContent(47,0.02479646);
   Allpt7TeVScaling->SetBinContent(48,0.02519956);
   Allpt7TeVScaling->SetBinContent(49,0.02376095);
   Allpt7TeVScaling->SetBinContent(50,0.02836189);
   Allpt7TeVScaling->SetBinContent(51,0.0339595);
   Allpt7TeVScaling->SetBinContent(52,0.03057126);
   Allpt7TeVScaling->SetBinContent(53,0.03238621);
   Allpt7TeVScaling->SetBinContent(54,0.0264272);
   Allpt7TeVScaling->SetBinContent(55,0.03203146);
   Allpt7TeVScaling->SetBinContent(56,0.02554187);
   Allpt7TeVScaling->SetBinContent(57,0.02901947);
   Allpt7TeVScaling->SetBinContent(58,0.02584016);
   Allpt7TeVScaling->SetBinContent(59,0.03055042);
   Allpt7TeVScaling->SetBinContent(60,0.02658193);
   Allpt7TeVScaling->SetBinContent(61,0.02553988);
   Allpt7TeVScaling->SetBinContent(62,0.02787508);
   Allpt7TeVScaling->SetBinContent(63,0.02836785);
   Allpt7TeVScaling->SetBinContent(64,0.03449766);
   Allpt7TeVScaling->SetBinContent(65,0.02710397);
   Allpt7TeVScaling->SetBinContent(66,0.02741966);
   Allpt7TeVScaling->SetBinContent(67,0.03184318);
   Allpt7TeVScaling->SetBinContent(68,0.02993716);
   Allpt7TeVScaling->SetBinContent(69,0.03113033);
   Allpt7TeVScaling->SetBinContent(70,0.01946293);
   Allpt7TeVScaling->SetBinContent(71,0.02890383);
   Allpt7TeVScaling->SetBinContent(72,0.02602078);
   Allpt7TeVScaling->SetBinContent(73,0.02063028);
   Allpt7TeVScaling->SetBinContent(74,0.03186679);
   Allpt7TeVScaling->SetBinContent(75,0.0399062);
   Allpt7TeVScaling->SetBinContent(76,0.02051914);
   Allpt7TeVScaling->SetBinContent(77,0.04454688);
   Allpt7TeVScaling->SetBinContent(78,0.04798449);
   Allpt7TeVScaling->SetBinContent(79,0.0311417);
   Allpt7TeVScaling->SetBinContent(80,0.03774521);
   Allpt7TeVScaling->SetBinContent(81,0.0219128);
   Allpt7TeVScaling->SetBinContent(82,0.02323492);
   Allpt7TeVScaling->SetBinContent(83,0.02412182);
   Allpt7TeVScaling->SetBinContent(84,0.03636555);
   Allpt7TeVScaling->SetBinContent(85,0.007488934);
   Allpt7TeVScaling->SetBinContent(86,0.01676086);
   Allpt7TeVScaling->SetBinContent(87,0.02427021);
   Allpt7TeVScaling->SetBinContent(88,0.02482882);
   Allpt7TeVScaling->SetBinContent(89,0.03680926);
   Allpt7TeVScaling->SetBinContent(90,0.01600685);
   Allpt7TeVScaling->SetBinContent(91,0.006885685);
   Allpt7TeVScaling->SetBinContent(93,0.0195411);
   Allpt7TeVScaling->SetBinContent(94,0.008862645);
   Allpt7TeVScaling->SetBinContent(95,0.009963223);
   Allpt7TeVScaling->SetBinContent(96,0.01399091);
   Allpt7TeVScaling->SetBinContent(97,0.007500135);
   Allpt7TeVScaling->SetBinContent(98,0.01922651);
   Allpt7TeVScaling->SetBinContent(103,0.03511665);
   Allpt7TeVScaling->SetBinError(13,0.03042345);
   Allpt7TeVScaling->SetBinError(14,0.007201474);
   Allpt7TeVScaling->SetBinError(15,0.004623859);
   Allpt7TeVScaling->SetBinError(16,0.004039713);
   Allpt7TeVScaling->SetBinError(17,0.003694238);
   Allpt7TeVScaling->SetBinError(18,0.003377785);
   Allpt7TeVScaling->SetBinError(19,0.00319956);
   Allpt7TeVScaling->SetBinError(20,0.00294949);
   Allpt7TeVScaling->SetBinError(21,0.001965652);
   Allpt7TeVScaling->SetBinError(22,0.001771706);
   Allpt7TeVScaling->SetBinError(23,0.0015726);
   Allpt7TeVScaling->SetBinError(24,0.001495611);
   Allpt7TeVScaling->SetBinError(25,0.00145454);
   Allpt7TeVScaling->SetBinError(26,0.001398941);
   Allpt7TeVScaling->SetBinError(27,0.001495129);
   Allpt7TeVScaling->SetBinError(28,0.001474393);
   Allpt7TeVScaling->SetBinError(29,0.001608602);
   Allpt7TeVScaling->SetBinError(30,0.001684714);
   Allpt7TeVScaling->SetBinError(31,0.001755251);
   Allpt7TeVScaling->SetBinError(32,0.001658064);
   Allpt7TeVScaling->SetBinError(33,0.001734642);
   Allpt7TeVScaling->SetBinError(34,0.001624169);
   Allpt7TeVScaling->SetBinError(35,0.001893747);
   Allpt7TeVScaling->SetBinError(36,0.001748243);
   Allpt7TeVScaling->SetBinError(37,0.001824646);
   Allpt7TeVScaling->SetBinError(38,0.001587514);
   Allpt7TeVScaling->SetBinError(39,0.001697203);
   Allpt7TeVScaling->SetBinError(40,0.002112545);
   Allpt7TeVScaling->SetBinError(41,0.002060819);
   Allpt7TeVScaling->SetBinError(42,0.002026571);
   Allpt7TeVScaling->SetBinError(43,0.001821759);
   Allpt7TeVScaling->SetBinError(44,0.002051793);
   Allpt7TeVScaling->SetBinError(45,0.002193957);
   Allpt7TeVScaling->SetBinError(46,0.001991362);
   Allpt7TeVScaling->SetBinError(47,0.002242771);
   Allpt7TeVScaling->SetBinError(48,0.002305595);
   Allpt7TeVScaling->SetBinError(49,0.002335076);
   Allpt7TeVScaling->SetBinError(50,0.002798601);
   Allpt7TeVScaling->SetBinError(51,0.003536788);
   Allpt7TeVScaling->SetBinError(52,0.003440952);
   Allpt7TeVScaling->SetBinError(53,0.003495971);
   Allpt7TeVScaling->SetBinError(54,0.003067046);
   Allpt7TeVScaling->SetBinError(55,0.003553133);
   Allpt7TeVScaling->SetBinError(56,0.003084681);
   Allpt7TeVScaling->SetBinError(57,0.003207651);
   Allpt7TeVScaling->SetBinError(58,0.002958898);
   Allpt7TeVScaling->SetBinError(59,0.003486212);
   Allpt7TeVScaling->SetBinError(60,0.003160784);
   Allpt7TeVScaling->SetBinError(61,0.00221753);
   Allpt7TeVScaling->SetBinError(62,0.002573063);
   Allpt7TeVScaling->SetBinError(63,0.003014171);
   Allpt7TeVScaling->SetBinError(64,0.003535776);
   Allpt7TeVScaling->SetBinError(65,0.003472028);
   Allpt7TeVScaling->SetBinError(66,0.003649103);
   Allpt7TeVScaling->SetBinError(67,0.004547596);
   Allpt7TeVScaling->SetBinError(68,0.004987629);
   Allpt7TeVScaling->SetBinError(69,0.005726305);
   Allpt7TeVScaling->SetBinError(70,0.003978767);
   Allpt7TeVScaling->SetBinError(71,0.00668463);
   Allpt7TeVScaling->SetBinError(72,0.005573138);
   Allpt7TeVScaling->SetBinError(73,0.005310436);
   Allpt7TeVScaling->SetBinError(74,0.007698246);
   Allpt7TeVScaling->SetBinError(75,0.01002912);
   Allpt7TeVScaling->SetBinError(76,0.006029285);
   Allpt7TeVScaling->SetBinError(77,0.0117225);
   Allpt7TeVScaling->SetBinError(78,0.01333887);
   Allpt7TeVScaling->SetBinError(79,0.009949805);
   Allpt7TeVScaling->SetBinError(80,0.01335105);
   Allpt7TeVScaling->SetBinError(81,0.004966615);
   Allpt7TeVScaling->SetBinError(82,0.006669427);
   Allpt7TeVScaling->SetBinError(83,0.007591192);
   Allpt7TeVScaling->SetBinError(84,0.01217267);
   Allpt7TeVScaling->SetBinError(85,0.002708538);
   Allpt7TeVScaling->SetBinError(86,0.008143581);
   Allpt7TeVScaling->SetBinError(87,0.01149673);
   Allpt7TeVScaling->SetBinError(88,0.01759547);
   Allpt7TeVScaling->SetBinError(89,0.01766708);
   Allpt7TeVScaling->SetBinError(90,0.01339557);
   Allpt7TeVScaling->SetBinError(91,0.004886052);
   Allpt7TeVScaling->SetBinError(93,0.01960311);
   Allpt7TeVScaling->SetBinError(94,0.008892003);
   Allpt7TeVScaling->SetBinError(95,0.01001297);
   Allpt7TeVScaling->SetBinError(96,0.008178781);
   Allpt7TeVScaling->SetBinError(97,0.007527505);
   Allpt7TeVScaling->SetBinError(98,0.01463862);
   Allpt7TeVScaling->SetBinError(103,0.03621621);
   Allpt7TeVScaling->SetMinimum(0);
   Allpt7TeVScaling->SetMaximum(0.04);
   Allpt7TeVScaling->SetEntries(1949.423);
   Allpt7TeVScaling->SetStats(0);
   Allpt7TeVScaling->SetMarkerStyle(20);
   Allpt7TeVScaling->GetXaxis()->SetTitle("p_{T}");
   Allpt7TeVScaling->GetXaxis()->SetRange(1,91);
   Allpt7TeVScaling->GetYaxis()->SetTitle("Ratio of E_{T}^{background}/E_{T}^{had, meas.}");
   Allpt7TeVScaling->GetYaxis()->SetTitleOffset(1.2);

   TH1D *Allpt900GeVScaling = new TH1D("Allpt900GeVScaling","Reconstructed E_{T} from misidentified electrons",111, xAxis1);
   Allpt900GeVScaling->SetBinContent(13,0.2098361);
   Allpt900GeVScaling->SetBinContent(14,0.1557776);
   Allpt900GeVScaling->SetBinContent(15,0.1480364);
   Allpt900GeVScaling->SetBinContent(16,0.1368096);
   Allpt900GeVScaling->SetBinContent(17,0.1244875);
   Allpt900GeVScaling->SetBinContent(18,0.1119436);
   Allpt900GeVScaling->SetBinContent(19,0.104192);
   Allpt900GeVScaling->SetBinContent(20,0.0938359);
   Allpt900GeVScaling->SetBinContent(21,0.08228509);
   Allpt900GeVScaling->SetBinContent(22,0.06783117);
   Allpt900GeVScaling->SetBinContent(23,0.05079209);
   Allpt900GeVScaling->SetBinContent(24,0.04795673);
   Allpt900GeVScaling->SetBinContent(25,0.04181242);
   Allpt900GeVScaling->SetBinContent(26,0.03994469);
   Allpt900GeVScaling->SetBinContent(27,0.03981588);
   Allpt900GeVScaling->SetBinContent(28,0.0377358);
   Allpt900GeVScaling->SetBinContent(29,0.0399667);
   Allpt900GeVScaling->SetBinContent(30,0.04095473);
   Allpt900GeVScaling->SetBinContent(31,0.04058456);
   Allpt900GeVScaling->SetBinContent(32,0.03641218);
   Allpt900GeVScaling->SetBinContent(33,0.03599966);
   Allpt900GeVScaling->SetBinContent(34,0.03853838);
   Allpt900GeVScaling->SetBinContent(35,0.04105409);
   Allpt900GeVScaling->SetBinContent(36,0.03461075);
   Allpt900GeVScaling->SetBinContent(37,0.03478729);
   Allpt900GeVScaling->SetBinContent(38,0.03150168);
   Allpt900GeVScaling->SetBinContent(39,0.03121329);
   Allpt900GeVScaling->SetBinContent(40,0.0370806);
   Allpt900GeVScaling->SetBinContent(41,0.03322445);
   Allpt900GeVScaling->SetBinContent(42,0.03200198);
   Allpt900GeVScaling->SetBinContent(43,0.02753292);
   Allpt900GeVScaling->SetBinContent(44,0.02776155);
   Allpt900GeVScaling->SetBinContent(45,0.02671865);
   Allpt900GeVScaling->SetBinContent(46,0.02612601);
   Allpt900GeVScaling->SetBinContent(47,0.02779079);
   Allpt900GeVScaling->SetBinContent(48,0.02762603);
   Allpt900GeVScaling->SetBinContent(49,0.0288977);
   Allpt900GeVScaling->SetBinContent(50,0.03181282);
   Allpt900GeVScaling->SetBinContent(51,0.04194707);
   Allpt900GeVScaling->SetBinContent(52,0.03287613);
   Allpt900GeVScaling->SetBinContent(53,0.04363778);
   Allpt900GeVScaling->SetBinContent(54,0.03003709);
   Allpt900GeVScaling->SetBinContent(55,0.03286254);
   Allpt900GeVScaling->SetBinContent(56,0.02585759);
   Allpt900GeVScaling->SetBinContent(57,0.03390592);
   Allpt900GeVScaling->SetBinContent(58,0.03212414);
   Allpt900GeVScaling->SetBinContent(59,0.03559123);
   Allpt900GeVScaling->SetBinContent(60,0.02464371);
   Allpt900GeVScaling->SetBinContent(61,0.03680997);
   Allpt900GeVScaling->SetBinContent(62,0.03569024);
   Allpt900GeVScaling->SetBinContent(63,0.03662885);
   Allpt900GeVScaling->SetBinContent(64,0.04690057);
   Allpt900GeVScaling->SetBinContent(65,0.04036927);
   Allpt900GeVScaling->SetBinContent(66,0.03237555);
   Allpt900GeVScaling->SetBinContent(67,0.04285642);
   Allpt900GeVScaling->SetBinContent(68,0.04778185);
   Allpt900GeVScaling->SetBinContent(69,0.0411135);
   Allpt900GeVScaling->SetBinContent(70,0.02465524);
   Allpt900GeVScaling->SetBinContent(71,0.03686565);
   Allpt900GeVScaling->SetBinContent(72,0.03694553);
   Allpt900GeVScaling->SetBinContent(73,0.02206743);
   Allpt900GeVScaling->SetBinContent(74,0.04854334);
   Allpt900GeVScaling->SetBinContent(75,0.04190287);
   Allpt900GeVScaling->SetBinContent(76,0.03255752);
   Allpt900GeVScaling->SetBinContent(77,0.0715847);
   Allpt900GeVScaling->SetBinContent(78,0.07329647);
   Allpt900GeVScaling->SetBinContent(79,0.02727783);
   Allpt900GeVScaling->SetBinContent(80,0.041527);
   Allpt900GeVScaling->SetBinContent(81,0.0400001);
   Allpt900GeVScaling->SetBinContent(82,0.0246785);
   Allpt900GeVScaling->SetBinContent(83,0.03546198);
   Allpt900GeVScaling->SetBinContent(84,0.0244789);
   Allpt900GeVScaling->SetBinContent(85,0.01692187);
   Allpt900GeVScaling->SetBinContent(86,0.01669025);
   Allpt900GeVScaling->SetBinContent(87,0.01506698);
   Allpt900GeVScaling->SetBinContent(88,0.02896953);
   Allpt900GeVScaling->SetBinContent(89,0.04019321);
   Allpt900GeVScaling->SetBinContent(90,0.0108909);
   Allpt900GeVScaling->SetBinContent(91,0.02059646);
   Allpt900GeVScaling->SetBinContent(92,0.01781294);
   Allpt900GeVScaling->SetBinContent(93,0.01930956);
   Allpt900GeVScaling->SetBinContent(94,0.0259471);
   Allpt900GeVScaling->SetBinContent(95,0.0114321);
   Allpt900GeVScaling->SetBinContent(98,0.04961994);
   Allpt900GeVScaling->SetBinContent(103,0.1436233);
   Allpt900GeVScaling->SetBinError(13,0.04135227);
   Allpt900GeVScaling->SetBinError(14,0.007955716);
   Allpt900GeVScaling->SetBinError(15,0.005353158);
   Allpt900GeVScaling->SetBinError(16,0.00465994);
   Allpt900GeVScaling->SetBinError(17,0.004317818);
   Allpt900GeVScaling->SetBinError(18,0.003992759);
   Allpt900GeVScaling->SetBinError(19,0.003794229);
   Allpt900GeVScaling->SetBinError(20,0.003428406);
   Allpt900GeVScaling->SetBinError(21,0.002309662);
   Allpt900GeVScaling->SetBinError(22,0.002082593);
   Allpt900GeVScaling->SetBinError(23,0.001768867);
   Allpt900GeVScaling->SetBinError(24,0.001734781);
   Allpt900GeVScaling->SetBinError(25,0.001668325);
   Allpt900GeVScaling->SetBinError(26,0.001736728);
   Allpt900GeVScaling->SetBinError(27,0.001793642);
   Allpt900GeVScaling->SetBinError(28,0.001824914);
   Allpt900GeVScaling->SetBinError(29,0.002216558);
   Allpt900GeVScaling->SetBinError(30,0.00224095);
   Allpt900GeVScaling->SetBinError(31,0.002244602);
   Allpt900GeVScaling->SetBinError(32,0.002073202);
   Allpt900GeVScaling->SetBinError(33,0.001913981);
   Allpt900GeVScaling->SetBinError(34,0.002456659);
   Allpt900GeVScaling->SetBinError(35,0.002596584);
   Allpt900GeVScaling->SetBinError(36,0.002232065);
   Allpt900GeVScaling->SetBinError(37,0.00251725);
   Allpt900GeVScaling->SetBinError(38,0.002196274);
   Allpt900GeVScaling->SetBinError(39,0.002513529);
   Allpt900GeVScaling->SetBinError(40,0.002973807);
   Allpt900GeVScaling->SetBinError(41,0.002807704);
   Allpt900GeVScaling->SetBinError(42,0.002773487);
   Allpt900GeVScaling->SetBinError(43,0.002527062);
   Allpt900GeVScaling->SetBinError(44,0.002921618);
   Allpt900GeVScaling->SetBinError(45,0.002907984);
   Allpt900GeVScaling->SetBinError(46,0.003115439);
   Allpt900GeVScaling->SetBinError(47,0.003224686);
   Allpt900GeVScaling->SetBinError(48,0.00312276);
   Allpt900GeVScaling->SetBinError(49,0.00393121);
   Allpt900GeVScaling->SetBinError(50,0.003948797);
   Allpt900GeVScaling->SetBinError(51,0.005565144);
   Allpt900GeVScaling->SetBinError(52,0.004615559);
   Allpt900GeVScaling->SetBinError(53,0.006431509);
   Allpt900GeVScaling->SetBinError(54,0.004787345);
   Allpt900GeVScaling->SetBinError(55,0.004506561);
   Allpt900GeVScaling->SetBinError(56,0.004277351);
   Allpt900GeVScaling->SetBinError(57,0.005084143);
   Allpt900GeVScaling->SetBinError(58,0.005085463);
   Allpt900GeVScaling->SetBinError(59,0.005312616);
   Allpt900GeVScaling->SetBinError(60,0.003899113);
   Allpt900GeVScaling->SetBinError(61,0.003977416);
   Allpt900GeVScaling->SetBinError(62,0.004292332);
   Allpt900GeVScaling->SetBinError(63,0.005072534);
   Allpt900GeVScaling->SetBinError(64,0.006271891);
   Allpt900GeVScaling->SetBinError(65,0.006463988);
   Allpt900GeVScaling->SetBinError(66,0.005869634);
   Allpt900GeVScaling->SetBinError(67,0.008405251);
   Allpt900GeVScaling->SetBinError(68,0.009847454);
   Allpt900GeVScaling->SetBinError(69,0.009749354);
   Allpt900GeVScaling->SetBinError(70,0.007179523);
   Allpt900GeVScaling->SetBinError(71,0.01004867);
   Allpt900GeVScaling->SetBinError(72,0.01065873);
   Allpt900GeVScaling->SetBinError(73,0.007849093);
   Allpt900GeVScaling->SetBinError(74,0.01402708);
   Allpt900GeVScaling->SetBinError(75,0.01330162);
   Allpt900GeVScaling->SetBinError(76,0.01245402);
   Allpt900GeVScaling->SetBinError(77,0.02130282);
   Allpt900GeVScaling->SetBinError(78,0.02296044);
   Allpt900GeVScaling->SetBinError(79,0.01175655);
   Allpt900GeVScaling->SetBinError(80,0.0187519);
   Allpt900GeVScaling->SetBinError(81,0.009932965);
   Allpt900GeVScaling->SetBinError(82,0.007798875);
   Allpt900GeVScaling->SetBinError(83,0.01126926);
   Allpt900GeVScaling->SetBinError(84,0.01123177);
   Allpt900GeVScaling->SetBinError(85,0.01013651);
   Allpt900GeVScaling->SetBinError(86,0.01187595);
   Allpt900GeVScaling->SetBinError(87,0.009213867);
   Allpt900GeVScaling->SetBinError(88,0.02054123);
   Allpt900GeVScaling->SetBinError(89,0.01976358);
   Allpt900GeVScaling->SetBinError(90,0.006322856);
   Allpt900GeVScaling->SetBinError(91,0.01400211);
   Allpt900GeVScaling->SetBinError(92,0.01786717);
   Allpt900GeVScaling->SetBinError(93,0.01938449);
   Allpt900GeVScaling->SetBinError(94,0.01853139);
   Allpt900GeVScaling->SetBinError(95,0.01149767);
   Allpt900GeVScaling->SetBinError(98,0.03786105);
   Allpt900GeVScaling->SetBinError(103,0.1493976);
   Allpt900GeVScaling->SetMinimum(0);
   Allpt900GeVScaling->SetMaximum(0.04);
   Allpt900GeVScaling->SetEntries(482.6447);
   Allpt900GeVScaling->SetStats(0);
   Allpt900GeVScaling->SetMarkerStyle(20);
   Allpt900GeVScaling->GetXaxis()->SetTitle("p_{T}");
   Allpt900GeVScaling->GetXaxis()->SetRange(1,91);
   Allpt900GeVScaling->GetYaxis()->SetTitle("Ratio of E_{T}^{background}/E_{T}^{had, meas.}");
   Allpt900GeVScaling->GetYaxis()->SetTitleOffset(1.2);
   Allpt7TeVScaling->Add(Allpt900GeVScaling);
   Allpt7TeVScaling->Scale(0.5);
   delete Allpt900GeVScaling;
   return Allpt7TeVScaling;
}

TH1D *pp276ITSBkgd(){
   TH1D *Allpt7TeVScaling = new TH1D("Allpt7TeVScaling","Reconstructed E_{T} from misidentified electrons",111, xAxis1);
   Allpt7TeVScaling->SetBinContent(11,0.2557976);
   Allpt7TeVScaling->SetBinContent(12,0.1130502);
   Allpt7TeVScaling->SetBinContent(13,0.07504393);
   Allpt7TeVScaling->SetBinContent(14,0.05774832);
   Allpt7TeVScaling->SetBinContent(15,0.06389806);
   Allpt7TeVScaling->SetBinContent(16,0.0468956);
   Allpt7TeVScaling->SetBinContent(17,0.04355167);
   Allpt7TeVScaling->SetBinContent(18,0.04035013);
   Allpt7TeVScaling->SetBinContent(19,0.04359472);
   Allpt7TeVScaling->SetBinContent(20,0.03392261);
   Allpt7TeVScaling->SetBinContent(21,0.0327676);
   Allpt7TeVScaling->SetBinContent(22,0.02916985);
   Allpt7TeVScaling->SetBinContent(23,0.02603138);
   Allpt7TeVScaling->SetBinContent(24,0.02298032);
   Allpt7TeVScaling->SetBinContent(25,0.02366939);
   Allpt7TeVScaling->SetBinContent(26,0.01941947);
   Allpt7TeVScaling->SetBinContent(27,0.02077003);
   Allpt7TeVScaling->SetBinContent(28,0.01985847);
   Allpt7TeVScaling->SetBinContent(29,0.0177459);
   Allpt7TeVScaling->SetBinContent(30,0.01889314);
   Allpt7TeVScaling->SetBinContent(31,0.01783637);
   Allpt7TeVScaling->SetBinContent(32,0.01742894);
   Allpt7TeVScaling->SetBinContent(33,0.01400672);
   Allpt7TeVScaling->SetBinContent(34,0.01509014);
   Allpt7TeVScaling->SetBinContent(35,0.01473409);
   Allpt7TeVScaling->SetBinContent(36,0.01337399);
   Allpt7TeVScaling->SetBinContent(37,0.01458302);
   Allpt7TeVScaling->SetBinContent(38,0.01103552);
   Allpt7TeVScaling->SetBinContent(39,0.01440217);
   Allpt7TeVScaling->SetBinContent(40,0.0147644);
   Allpt7TeVScaling->SetBinContent(41,0.01459072);
   Allpt7TeVScaling->SetBinContent(42,0.01318416);
   Allpt7TeVScaling->SetBinContent(43,0.01364365);
   Allpt7TeVScaling->SetBinContent(44,0.01186023);
   Allpt7TeVScaling->SetBinContent(45,0.01218222);
   Allpt7TeVScaling->SetBinContent(46,0.01015807);
   Allpt7TeVScaling->SetBinContent(47,0.01158454);
   Allpt7TeVScaling->SetBinContent(48,0.009338749);
   Allpt7TeVScaling->SetBinContent(49,0.01320016);
   Allpt7TeVScaling->SetBinContent(50,0.01029357);
   Allpt7TeVScaling->SetBinContent(51,0.01238829);
   Allpt7TeVScaling->SetBinContent(52,0.008714745);
   Allpt7TeVScaling->SetBinContent(53,0.01300408);
   Allpt7TeVScaling->SetBinContent(54,0.009445714);
   Allpt7TeVScaling->SetBinContent(55,0.01301231);
   Allpt7TeVScaling->SetBinContent(56,0.01010809);
   Allpt7TeVScaling->SetBinContent(57,0.01181843);
   Allpt7TeVScaling->SetBinContent(58,0.01037317);
   Allpt7TeVScaling->SetBinContent(59,0.01172911);
   Allpt7TeVScaling->SetBinContent(60,0.007576451);
   Allpt7TeVScaling->SetBinContent(61,0.011185);
   Allpt7TeVScaling->SetBinContent(62,0.01009882);
   Allpt7TeVScaling->SetBinContent(63,0.01042078);
   Allpt7TeVScaling->SetBinContent(64,0.01079128);
   Allpt7TeVScaling->SetBinContent(65,0.008657009);
   Allpt7TeVScaling->SetBinContent(66,0.009641373);
   Allpt7TeVScaling->SetBinContent(67,0.01201938);
   Allpt7TeVScaling->SetBinContent(68,0.01031088);
   Allpt7TeVScaling->SetBinContent(69,0.009905995);
   Allpt7TeVScaling->SetBinContent(70,0.009569216);
   Allpt7TeVScaling->SetBinContent(71,0.009473246);
   Allpt7TeVScaling->SetBinContent(72,0.00825013);
   Allpt7TeVScaling->SetBinContent(73,0.008405926);
   Allpt7TeVScaling->SetBinContent(74,0.006479467);
   Allpt7TeVScaling->SetBinContent(75,0.01135131);
   Allpt7TeVScaling->SetBinContent(76,0.007318432);
   Allpt7TeVScaling->SetBinContent(77,0.01145386);
   Allpt7TeVScaling->SetBinContent(78,0.01425465);
   Allpt7TeVScaling->SetBinContent(79,0.01205772);
   Allpt7TeVScaling->SetBinContent(80,0.009454997);
   Allpt7TeVScaling->SetBinContent(81,0.006863416);
   Allpt7TeVScaling->SetBinContent(82,0.006285407);
   Allpt7TeVScaling->SetBinContent(83,0.004974776);
   Allpt7TeVScaling->SetBinContent(84,0.00748478);
   Allpt7TeVScaling->SetBinContent(85,0.002485323);
   Allpt7TeVScaling->SetBinContent(86,0.009973802);
   Allpt7TeVScaling->SetBinContent(87,0.00543998);
   Allpt7TeVScaling->SetBinContent(89,0.01703395);
   Allpt7TeVScaling->SetBinContent(90,0.004394557);
   Allpt7TeVScaling->SetBinContent(96,0.01833441);
   Allpt7TeVScaling->SetBinContent(98,0.02964621);
   Allpt7TeVScaling->SetBinContent(101,0.0574538);
   Allpt7TeVScaling->SetBinError(11,0.07407797);
   Allpt7TeVScaling->SetBinError(12,0.01275005);
   Allpt7TeVScaling->SetBinError(13,0.00608927);
   Allpt7TeVScaling->SetBinError(14,0.00428747);
   Allpt7TeVScaling->SetBinError(15,0.004071853);
   Allpt7TeVScaling->SetBinError(16,0.003225169);
   Allpt7TeVScaling->SetBinError(17,0.002968269);
   Allpt7TeVScaling->SetBinError(18,0.002761875);
   Allpt7TeVScaling->SetBinError(19,0.002800096);
   Allpt7TeVScaling->SetBinError(20,0.002386005);
   Allpt7TeVScaling->SetBinError(21,0.001589813);
   Allpt7TeVScaling->SetBinError(22,0.001475329);
   Allpt7TeVScaling->SetBinError(23,0.001364666);
   Allpt7TeVScaling->SetBinError(24,0.00126007);
   Allpt7TeVScaling->SetBinError(25,0.001265789);
   Allpt7TeVScaling->SetBinError(26,0.001156208);
   Allpt7TeVScaling->SetBinError(27,0.001182566);
   Allpt7TeVScaling->SetBinError(28,0.001171098);
   Allpt7TeVScaling->SetBinError(29,0.001122203);
   Allpt7TeVScaling->SetBinError(30,0.001173633);
   Allpt7TeVScaling->SetBinError(31,0.001154546);
   Allpt7TeVScaling->SetBinError(32,0.001162445);
   Allpt7TeVScaling->SetBinError(33,0.001067574);
   Allpt7TeVScaling->SetBinError(34,0.001116243);
   Allpt7TeVScaling->SetBinError(35,0.001144851);
   Allpt7TeVScaling->SetBinError(36,0.001110251);
   Allpt7TeVScaling->SetBinError(37,0.001192711);
   Allpt7TeVScaling->SetBinError(38,0.001062774);
   Allpt7TeVScaling->SetBinError(39,0.001255752);
   Allpt7TeVScaling->SetBinError(40,0.001319021);
   Allpt7TeVScaling->SetBinError(41,0.001323284);
   Allpt7TeVScaling->SetBinError(42,0.001303463);
   Allpt7TeVScaling->SetBinError(43,0.001388458);
   Allpt7TeVScaling->SetBinError(44,0.001312344);
   Allpt7TeVScaling->SetBinError(45,0.001383885);
   Allpt7TeVScaling->SetBinError(46,0.001289152);
   Allpt7TeVScaling->SetBinError(47,0.001430287);
   Allpt7TeVScaling->SetBinError(48,0.001302716);
   Allpt7TeVScaling->SetBinError(49,0.001578988);
   Allpt7TeVScaling->SetBinError(50,0.0014581);
   Allpt7TeVScaling->SetBinError(51,0.00165325);
   Allpt7TeVScaling->SetBinError(52,0.001404793);
   Allpt7TeVScaling->SetBinError(53,0.001757303);
   Allpt7TeVScaling->SetBinError(54,0.001562987);
   Allpt7TeVScaling->SetBinError(55,0.001911932);
   Allpt7TeVScaling->SetBinError(56,0.001695077);
   Allpt7TeVScaling->SetBinError(57,0.001932989);
   Allpt7TeVScaling->SetBinError(58,0.001844963);
   Allpt7TeVScaling->SetBinError(59,0.002000117);
   Allpt7TeVScaling->SetBinError(60,0.001689391);
   Allpt7TeVScaling->SetBinError(61,0.001386339);
   Allpt7TeVScaling->SetBinError(62,0.001349763);
   Allpt7TeVScaling->SetBinError(63,0.001454327);
   Allpt7TeVScaling->SetBinError(64,0.001584791);
   Allpt7TeVScaling->SetBinError(65,0.00151686);
   Allpt7TeVScaling->SetBinError(66,0.001668507);
   Allpt7TeVScaling->SetBinError(67,0.001963962);
   Allpt7TeVScaling->SetBinError(68,0.001970941);
   Allpt7TeVScaling->SetBinError(69,0.001992913);
   Allpt7TeVScaling->SetBinError(70,0.002006771);
   Allpt7TeVScaling->SetBinError(71,0.002187254);
   Allpt7TeVScaling->SetBinError(72,0.002140473);
   Allpt7TeVScaling->SetBinError(73,0.002256925);
   Allpt7TeVScaling->SetBinError(74,0.002056219);
   Allpt7TeVScaling->SetBinError(75,0.002854892);
   Allpt7TeVScaling->SetBinError(76,0.002450159);
   Allpt7TeVScaling->SetBinError(77,0.003197387);
   Allpt7TeVScaling->SetBinError(78,0.003840947);
   Allpt7TeVScaling->SetBinError(79,0.003660464);
   Allpt7TeVScaling->SetBinError(80,0.003363956);
   Allpt7TeVScaling->SetBinError(81,0.001624793);
   Allpt7TeVScaling->SetBinError(82,0.00181065);
   Allpt7TeVScaling->SetBinError(83,0.001888966);
   Allpt7TeVScaling->SetBinError(84,0.002671173);
   Allpt7TeVScaling->SetBinError(85,0.001759843);
   Allpt7TeVScaling->SetBinError(86,0.004094911);
   Allpt7TeVScaling->SetBinError(87,0.003871198);
   Allpt7TeVScaling->SetBinError(89,0.008630108);
   Allpt7TeVScaling->SetBinError(90,0.004405367);
   Allpt7TeVScaling->SetBinError(96,0.01308755);
   Allpt7TeVScaling->SetBinError(98,0.02127892);
   Allpt7TeVScaling->SetBinError(101,0.05904744);
   Allpt7TeVScaling->SetMinimum(0);
   Allpt7TeVScaling->SetMaximum(0.04);
   Allpt7TeVScaling->SetEntries(300.6659);
   Allpt7TeVScaling->SetStats(0);
   Allpt7TeVScaling->SetMarkerStyle(20);
   Allpt7TeVScaling->GetXaxis()->SetTitle("p_{T}");
   Allpt7TeVScaling->GetXaxis()->SetRange(1,61);
   Allpt7TeVScaling->GetYaxis()->SetTitle("Ratio of E_{T}^{background}/E_{T}^{had, meas.}");
   Allpt7TeVScaling->GetYaxis()->SetTitleOffset(1.2);
   TH1D *Allpt900GeVScaling = new TH1D("Allpt900GeVScaling","Reconstructed E_{T} from misidentified electrons",111, xAxis1);
   Allpt900GeVScaling->SetBinContent(9,1);
   Allpt900GeVScaling->SetBinContent(11,0.2834258);
   Allpt900GeVScaling->SetBinContent(12,0.1188776);
   Allpt900GeVScaling->SetBinContent(13,0.08085647);
   Allpt900GeVScaling->SetBinContent(14,0.05677858);
   Allpt900GeVScaling->SetBinContent(15,0.0616757);
   Allpt900GeVScaling->SetBinContent(16,0.04844882);
   Allpt900GeVScaling->SetBinContent(17,0.03999142);
   Allpt900GeVScaling->SetBinContent(18,0.03890007);
   Allpt900GeVScaling->SetBinContent(19,0.0406008);
   Allpt900GeVScaling->SetBinContent(20,0.02900627);
   Allpt900GeVScaling->SetBinContent(21,0.03088961);
   Allpt900GeVScaling->SetBinContent(22,0.03033992);
   Allpt900GeVScaling->SetBinContent(23,0.02521969);
   Allpt900GeVScaling->SetBinContent(24,0.02143008);
   Allpt900GeVScaling->SetBinContent(25,0.02266187);
   Allpt900GeVScaling->SetBinContent(26,0.01981764);
   Allpt900GeVScaling->SetBinContent(27,0.02127366);
   Allpt900GeVScaling->SetBinContent(28,0.0192753);
   Allpt900GeVScaling->SetBinContent(29,0.01854744);
   Allpt900GeVScaling->SetBinContent(30,0.01903471);
   Allpt900GeVScaling->SetBinContent(31,0.01819287);
   Allpt900GeVScaling->SetBinContent(32,0.01715784);
   Allpt900GeVScaling->SetBinContent(33,0.01488607);
   Allpt900GeVScaling->SetBinContent(34,0.0156348);
   Allpt900GeVScaling->SetBinContent(35,0.01699235);
   Allpt900GeVScaling->SetBinContent(36,0.01398954);
   Allpt900GeVScaling->SetBinContent(37,0.01435865);
   Allpt900GeVScaling->SetBinContent(38,0.01236315);
   Allpt900GeVScaling->SetBinContent(39,0.01282288);
   Allpt900GeVScaling->SetBinContent(40,0.01668125);
   Allpt900GeVScaling->SetBinContent(41,0.01372736);
   Allpt900GeVScaling->SetBinContent(42,0.01433745);
   Allpt900GeVScaling->SetBinContent(43,0.01477496);
   Allpt900GeVScaling->SetBinContent(44,0.0116776);
   Allpt900GeVScaling->SetBinContent(45,0.01245205);
   Allpt900GeVScaling->SetBinContent(46,0.009475905);
   Allpt900GeVScaling->SetBinContent(47,0.01117031);
   Allpt900GeVScaling->SetBinContent(48,0.009558273);
   Allpt900GeVScaling->SetBinContent(49,0.01404047);
   Allpt900GeVScaling->SetBinContent(50,0.01037464);
   Allpt900GeVScaling->SetBinContent(51,0.01195013);
   Allpt900GeVScaling->SetBinContent(52,0.01155682);
   Allpt900GeVScaling->SetBinContent(53,0.0137127);
   Allpt900GeVScaling->SetBinContent(54,0.009255665);
   Allpt900GeVScaling->SetBinContent(55,0.01133598);
   Allpt900GeVScaling->SetBinContent(56,0.009080946);
   Allpt900GeVScaling->SetBinContent(57,0.01124535);
   Allpt900GeVScaling->SetBinContent(58,0.008398657);
   Allpt900GeVScaling->SetBinContent(59,0.008324869);
   Allpt900GeVScaling->SetBinContent(60,0.007103489);
   Allpt900GeVScaling->SetBinContent(61,0.01142644);
   Allpt900GeVScaling->SetBinContent(62,0.0110906);
   Allpt900GeVScaling->SetBinContent(63,0.01076479);
   Allpt900GeVScaling->SetBinContent(64,0.01168746);
   Allpt900GeVScaling->SetBinContent(65,0.008513672);
   Allpt900GeVScaling->SetBinContent(66,0.00828279);
   Allpt900GeVScaling->SetBinContent(67,0.01140648);
   Allpt900GeVScaling->SetBinContent(68,0.01111323);
   Allpt900GeVScaling->SetBinContent(69,0.007425678);
   Allpt900GeVScaling->SetBinContent(70,0.009847098);
   Allpt900GeVScaling->SetBinContent(71,0.009834419);
   Allpt900GeVScaling->SetBinContent(72,0.008930991);
   Allpt900GeVScaling->SetBinContent(73,0.00669819);
   Allpt900GeVScaling->SetBinContent(74,0.00662618);
   Allpt900GeVScaling->SetBinContent(75,0.01370167);
   Allpt900GeVScaling->SetBinContent(76,0.006451024);
   Allpt900GeVScaling->SetBinContent(77,0.01307605);
   Allpt900GeVScaling->SetBinContent(78,0.01579264);
   Allpt900GeVScaling->SetBinContent(79,0.01255205);
   Allpt900GeVScaling->SetBinContent(80,0.01010045);
   Allpt900GeVScaling->SetBinContent(81,0.01082155);
   Allpt900GeVScaling->SetBinContent(82,0.006495255);
   Allpt900GeVScaling->SetBinContent(83,0.005512525);
   Allpt900GeVScaling->SetBinContent(84,0.003527945);
   Allpt900GeVScaling->SetBinContent(85,0.001806946);
   Allpt900GeVScaling->SetBinContent(86,0.004225982);
   Allpt900GeVScaling->SetBinContent(89,0.01106803);
   Allpt900GeVScaling->SetBinContent(90,0.005044269);
   Allpt900GeVScaling->SetBinContent(98,0.04456838);
   Allpt900GeVScaling->SetBinContent(101,0.07960857);
   Allpt900GeVScaling->SetBinError(9,1.414214);
   Allpt900GeVScaling->SetBinError(11,0.08936516);
   Allpt900GeVScaling->SetBinError(12,0.01474769);
   Allpt900GeVScaling->SetBinError(13,0.006947232);
   Allpt900GeVScaling->SetBinError(14,0.00470276);
   Allpt900GeVScaling->SetBinError(15,0.004394222);
   Allpt900GeVScaling->SetBinError(16,0.003598075);
   Allpt900GeVScaling->SetBinError(17,0.00312581);
   Allpt900GeVScaling->SetBinError(18,0.003002622);
   Allpt900GeVScaling->SetBinError(19,0.002976063);
   Allpt900GeVScaling->SetBinError(20,0.002438669);
   Allpt900GeVScaling->SetBinError(21,0.00172144);
   Allpt900GeVScaling->SetBinError(22,0.001661658);
   Allpt900GeVScaling->SetBinError(23,0.001477901);
   Allpt900GeVScaling->SetBinError(24,0.001329128);
   Allpt900GeVScaling->SetBinError(25,0.00137131);
   Allpt900GeVScaling->SetBinError(26,0.001299862);
   Allpt900GeVScaling->SetBinError(27,0.001332285);
   Allpt900GeVScaling->SetBinError(28,0.001274491);
   Allpt900GeVScaling->SetBinError(29,0.001266741);
   Allpt900GeVScaling->SetBinError(30,0.001305369);
   Allpt900GeVScaling->SetBinError(31,0.001295113);
   Allpt900GeVScaling->SetBinError(32,0.001278961);
   Allpt900GeVScaling->SetBinError(33,0.00122243);
   Allpt900GeVScaling->SetBinError(34,0.001256858);
   Allpt900GeVScaling->SetBinError(35,0.001357681);
   Allpt900GeVScaling->SetBinError(36,0.001267233);
   Allpt900GeVScaling->SetBinError(37,0.00131839);
   Allpt900GeVScaling->SetBinError(38,0.001256835);
   Allpt900GeVScaling->SetBinError(39,0.001306764);
   Allpt900GeVScaling->SetBinError(40,0.001544152);
   Allpt900GeVScaling->SetBinError(41,0.001438378);
   Allpt900GeVScaling->SetBinError(42,0.001507702);
   Allpt900GeVScaling->SetBinError(43,0.001615741);
   Allpt900GeVScaling->SetBinError(44,0.00143817);
   Allpt900GeVScaling->SetBinError(45,0.001536765);
   Allpt900GeVScaling->SetBinError(46,0.001377493);
   Allpt900GeVScaling->SetBinError(47,0.001536932);
   Allpt900GeVScaling->SetBinError(48,0.001466076);
   Allpt900GeVScaling->SetBinError(49,0.001814185);
   Allpt900GeVScaling->SetBinError(50,0.001620601);
   Allpt900GeVScaling->SetBinError(51,0.001814658);
   Allpt900GeVScaling->SetBinError(52,0.001822325);
   Allpt900GeVScaling->SetBinError(53,0.001976947);
   Allpt900GeVScaling->SetBinError(54,0.001729338);
   Allpt900GeVScaling->SetBinError(55,0.001988634);
   Allpt900GeVScaling->SetBinError(56,0.001792157);
   Allpt900GeVScaling->SetBinError(57,0.00206644);
   Allpt900GeVScaling->SetBinError(58,0.001842037);
   Allpt900GeVScaling->SetBinError(59,0.001826512);
   Allpt900GeVScaling->SetBinError(60,0.001732292);
   Allpt900GeVScaling->SetBinError(61,0.001512696);
   Allpt900GeVScaling->SetBinError(62,0.001543542);
   Allpt900GeVScaling->SetBinError(63,0.001633358);
   Allpt900GeVScaling->SetBinError(64,0.001816012);
   Allpt900GeVScaling->SetBinError(65,0.001680175);
   Allpt900GeVScaling->SetBinError(66,0.001692314);
   Allpt900GeVScaling->SetBinError(67,0.00213108);
   Allpt900GeVScaling->SetBinError(68,0.002295308);
   Allpt900GeVScaling->SetBinError(69,0.001926907);
   Allpt900GeVScaling->SetBinError(70,0.002272706);
   Allpt900GeVScaling->SetBinError(71,0.002475109);
   Allpt900GeVScaling->SetBinError(72,0.002399123);
   Allpt900GeVScaling->SetBinError(73,0.002241209);
   Allpt900GeVScaling->SetBinError(74,0.002365641);
   Allpt900GeVScaling->SetBinError(75,0.003450507);
   Allpt900GeVScaling->SetBinError(76,0.002503392);
   Allpt900GeVScaling->SetBinError(77,0.003801781);
   Allpt900GeVScaling->SetBinError(78,0.004417338);
   Allpt900GeVScaling->SetBinError(79,0.00421245);
   Allpt900GeVScaling->SetBinError(80,0.003844397);
   Allpt900GeVScaling->SetBinError(81,0.002274688);
   Allpt900GeVScaling->SetBinError(82,0.002043975);
   Allpt900GeVScaling->SetBinError(83,0.002262563);
   Allpt900GeVScaling->SetBinError(84,0.00204239);
   Allpt900GeVScaling->SetBinError(85,0.00180851);
   Allpt900GeVScaling->SetBinError(86,0.002994568);
   Allpt900GeVScaling->SetBinError(89,0.007939965);
   Allpt900GeVScaling->SetBinError(90,0.005058548);
   Allpt900GeVScaling->SetBinError(98,0.03222111);
   Allpt900GeVScaling->SetBinError(101,0.08264771);
   Allpt900GeVScaling->SetMinimum(0);
   Allpt900GeVScaling->SetMaximum(0.04);
   Allpt900GeVScaling->SetEntries(3.837105);
   Allpt900GeVScaling->SetStats(0);
   Allpt900GeVScaling->SetMarkerStyle(20);
   Allpt900GeVScaling->GetXaxis()->SetTitle("p_{T}");
   Allpt900GeVScaling->GetXaxis()->SetRange(1,61);
   Allpt900GeVScaling->GetYaxis()->SetTitle("Ratio of E_{T}^{background}/E_{T}^{had, meas.}");
   Allpt900GeVScaling->GetYaxis()->SetTitleOffset(1.2);

   Allpt7TeVScaling->Add(Allpt900GeVScaling);
   Allpt7TeVScaling->Scale(0.5);
   delete Allpt900GeVScaling;
   return Allpt7TeVScaling;
}
