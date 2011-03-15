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
   sprintf(outfilename,"corrections.%s.%s.%s.root",shortprodname,system,sim);
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
     hadCorrectionEMCAL->SetNotIDConstCorrectionTPCLowBound(1./NotIDConstTPC*(1.0-0.010)/0.996);
     hadCorrectionEMCAL->SetNotIDConstCorrectionITSLowBound(1./NotIDConstITS*(1.0-0.010)/0.996);
     hadCorrectionEMCAL->SetNotIDConstCorrectionTPCHighBound(1./NotIDConstTPC*(1.0+0.010)/0.996);
     hadCorrectionEMCAL->SetNotIDConstCorrectionITSHighBound(1./NotIDConstITS*(1.0+0.010)/0.996);
   }
   else{
     hadCorrectionEMCAL->SetNotIDConstCorrectionTPCLowBound(1./NotIDConstTPC*(1.0-0.022)/0.976);
     hadCorrectionEMCAL->SetNotIDConstCorrectionITSLowBound(1./NotIDConstITS*(1.0-0.022)/0.976);
     hadCorrectionEMCAL->SetNotIDConstCorrectionTPCHighBound(1./NotIDConstTPC*(1.0+0.022)/0.976);
     hadCorrectionEMCAL->SetNotIDConstCorrectionITSHighBound(1./NotIDConstITS*(1.0+0.022)/0.976);
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

   TH1D *efficiencyPionTPC = GetHistoEfficiency(etacut,"hEfficiencyPionTPC",1,1,20,true,true);
   hadCorrectionEMCAL->SetEfficiencyPionTPC(efficiencyPionTPC);
   TH1D *efficiencyKaonTPC = GetHistoEfficiency(etacut,"hEfficiencyKaonTPC",2,1,20,true,true);
   hadCorrectionEMCAL->SetEfficiencyKaonTPC(efficiencyKaonTPC);
   TH1D *efficiencyProtonTPC = GetHistoEfficiency(etacut,"hEfficiencyProtonTPC",3,1,20,true,true);
   hadCorrectionEMCAL->SetEfficiencyProtonTPC(efficiencyProtonTPC);
   TH1D *efficiencyHadronTPC = GetHistoEfficiency(etacut,"hEfficiencyHadronTPC",0,1,20,true,true);
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
     TH1D *efficiencyPionTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB0",1,1,20,true,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB5",1,1,20,true,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB10",1,1,20,true,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionEMCAL->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB10->Clone(Form("Test%i",i)),i);

     TH1D *efficiencyKaonTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB0",2,1,20,true,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB5",2,1,20,true,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB10",2,1,20,true,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionEMCAL->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB10->Clone(Form("Test%i",i)),i);//Kaon

     TH1D *efficiencyProtonTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB0",3,1,20,true,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB5",3,1,20,true,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB10",3,1,20,true,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionEMCAL->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB10->Clone(Form("Test%i",i)),i);//Proton

     TH1D *efficiencyHadronTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB0",0,1,20,true,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB5",0,1,20,true,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB10",0,1,20,true,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionEMCAL->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB10->Clone(Form("Test%i",i)),i);//Hadron


     TH1D *efficiencyPionITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB0",1,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB5",1,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB10",1,1,20,false,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionEMCAL->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB10->Clone(Form("Test%i",i)),i);//Pion

     TH1D *efficiencyKaonITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB0",2,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB5",2,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB10",2,1,20,false,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionEMCAL->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB10->Clone(Form("Test%i",i)),i);//Kaon

     TH1D *efficiencyProtonITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB0",3,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB5",3,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB10",3,1,20,false,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionEMCAL->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB10->Clone(Form("Test%i",i)),i);//Proton

     TH1D *efficiencyHadronITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB0",0,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionEMCAL->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB5",0,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionEMCAL->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB10",0,1,20,false,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionEMCAL->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB10->Clone(Form("Test%i",i)),i);//Hadron
   }//EMCAL

   //CorrEfficiencyPlots(true,prodname,shortprodname);
   //CorrEfficiencyPlots(false,prodname,shortprodname,infilename);


   //hadCorrectionEMCAL->GetEfficiencyHadronTPC()->Draw();
   TH1D *backgroundTPC = GetHistoCorrBkgd(etacut,"hBackgroundTPC",true,ispp,forSim);
   TH1D *backgroundITS = GetHistoCorrBkgd(etacut,"hBackgroundITS",false,ispp,forSim);
   hadCorrectionEMCAL->SetBackgroundCorrectionTPC(backgroundTPC);
   hadCorrectionEMCAL->SetBackgroundCorrectionITS(backgroundITS);
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
     hadCorrectionPHOS->SetNeutralCorrectionLowBound(neutralCorr*(1.0-.013)/0.736);
     hadCorrectionPHOS->SetNeutralCorrectionHighBound(neutralCorr*(1.0+.013)/0.736);
   }
   else{
     hadCorrectionPHOS->SetNeutralCorrectionLowBound(neutralCorr*(1.0-0.049)/0.689);
     hadCorrectionPHOS->SetNeutralCorrectionHighBound(neutralCorr*(1.0+0.049)/0.689);
   }

   float hadronicCorr = CorrNeutral(ptcut,prodname,shortprodname,ispp,forSim,TPC,true,etacut);
   hadCorrectionPHOS->SetNotHadronicCorrection(hadronicCorr);
   if(ispp){
     hadCorrectionPHOS->SetNotHadronicCorrectionLowBound(neutralCorr*(1.0-0.008)/0.571);
     hadCorrectionPHOS->SetNotHadronicCorrectionHighBound(neutralCorr*(1.0+0.008)/0.571);
   }
   else{
     hadCorrectionPHOS->SetNotHadronicCorrectionLowBound(neutralCorr*(1.0-0.023)/0.549);
     hadCorrectionPHOS->SetNotHadronicCorrectionHighBound(neutralCorr*(1.0+0.023)/0.549);
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
     hadCorrectionPHOS->SetNotIDConstCorrectionTPCLowBound(1./NotIDConstTPC*(1.0-0.010)/0.996);
     hadCorrectionPHOS->SetNotIDConstCorrectionITSLowBound(1./NotIDConstITS*(1.0-0.010)/0.996);
     hadCorrectionPHOS->SetNotIDConstCorrectionTPCHighBound(1./NotIDConstTPC*(1.0+0.010)/0.996);
     hadCorrectionPHOS->SetNotIDConstCorrectionITSHighBound(1./NotIDConstITS*(1.0+0.010)/0.996);
   }
   else{
     hadCorrectionPHOS->SetNotIDConstCorrectionTPCLowBound(1./NotIDConstTPC*(1.0-0.022)/0.976);
     hadCorrectionPHOS->SetNotIDConstCorrectionITSLowBound(1./NotIDConstITS*(1.0-0.022)/0.976);
     hadCorrectionPHOS->SetNotIDConstCorrectionTPCHighBound(1./NotIDConstTPC*(1.0+0.022)/0.976);
     hadCorrectionPHOS->SetNotIDConstCorrectionITSHighBound(1./NotIDConstITS*(1.0+0.022)/0.976);
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

   TH1D *efficiencyPionTPC = GetHistoEfficiency(etacut,"hEfficiencyPionTPC",1,1,20,true,true);
   TH1D *efficiencyKaonTPC = GetHistoEfficiency(etacut,"hEfficiencyKaonTPC",2,1,20,true,true);
   TH1D *efficiencyProtonTPC = GetHistoEfficiency(etacut,"hEfficiencyProtonTPC",3,1,20,true,true);
   TH1D *efficiencyHadronTPC = GetHistoEfficiency(etacut,"hEfficiencyHadronTPC",0,1,20,true,true);
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
     TH1D *efficiencyPionTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB0",1,1,20,true,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB5",1,1,20,true,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyPionTPCCB10",1,1,20,true,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionPHOS->SetEfficiencyPionTPC((TH1D*)efficiencyPionTPCCB10->Clone(Form("Test%i",i)),i);

     TH1D *efficiencyKaonTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB0",2,1,20,true,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB5",2,1,20,true,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyKaonTPCCB10",2,1,20,true,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionPHOS->SetEfficiencyKaonTPC((TH1D*)efficiencyKaonTPCCB10->Clone(Form("Test%i",i)),i);//Kaon

     TH1D *efficiencyProtonTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB0",3,1,20,true,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB5",3,1,20,true,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyProtonTPCCB10",3,1,20,true,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionPHOS->SetEfficiencyProtonTPC((TH1D*)efficiencyProtonTPCCB10->Clone(Form("Test%i",i)),i);//Proton

     TH1D *efficiencyHadronTPCCB0 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB0",0,1,20,true,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronTPCCB5 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB5",0,1,20,true,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronTPCCB10 = GetHistoEfficiency(etacut,"hEfficiencyHadronTPCCB10",0,1,20,true,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionPHOS->SetEfficiencyHadronTPC((TH1D*)efficiencyHadronTPCCB10->Clone(Form("Test%i",i)),i);//Hadron


     TH1D *efficiencyPionITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB0",1,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB5",1,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyPionITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyPionITSCB10",1,1,20,false,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionPHOS->SetEfficiencyPionITS((TH1D*)efficiencyPionITSCB10->Clone(Form("Test%i",i)),i);//Pion

     TH1D *efficiencyKaonITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB0",2,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB5",2,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyKaonITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyKaonITSCB10",2,1,20,false,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionPHOS->SetEfficiencyKaonITS((TH1D*)efficiencyKaonITSCB10->Clone(Form("Test%i",i)),i);//Kaon

     TH1D *efficiencyProtonITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB0",3,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB5",3,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyProtonITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyProtonITSCB10",3,1,20,false,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionPHOS->SetEfficiencyProtonITS((TH1D*)efficiencyProtonITSCB10->Clone(Form("Test%i",i)),i);//Proton

     TH1D *efficiencyHadronITSCB0 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB0",0,1,20,false,true,0,4);
     for(int i=0;i<=4;i++) hadCorrectionPHOS->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB0->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronITSCB5 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB5",0,1,20,false,true,5,9);
     for(int i=5;i<=9;i++) hadCorrectionPHOS->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB5->Clone(Form("Test%i",i)),i);
     TH1D *efficiencyHadronITSCB10 = GetHistoEfficiency(etacut,"hEfficiencyHadronITSCB10",0,1,20,false,true,10,15);
     for(int i=10;i<=19;i++) hadCorrectionPHOS->SetEfficiencyHadronITS((TH1D*)efficiencyHadronITSCB10->Clone(Form("Test%i",i)),i);//Hadron
   }//EMCAL

   TH1D *backgroundTPC = GetHistoCorrBkgd(etacut,"hBackgroundTPC",true,ispp,forSim);
   TH1D *backgroundITS = GetHistoCorrBkgd(etacut,"hBackgroundITS",false,ispp,forSim);
   hadCorrectionPHOS->SetBackgroundCorrectionTPC(backgroundTPC);
   hadCorrectionPHOS->SetBackgroundCorrectionITS(backgroundITS);
   CorrBkgdPlots(prodname,shortprodname,true,ispp,forSim);
   CorrBkgdPlots(prodname,shortprodname,false,ispp,forSim);

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
  if(TPC) myname = mynameTPC;
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
  if(TPC) myname = mynameTPC;
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
  if(TPC&&!ITS) myname = mynameTPC;
  if(TPC&&ITS) myname = mynameTPCITS;
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
      }
    }
    else{
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiPlus")))->Clone("RecoHadron");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiMinus")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KMinus")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KPlus")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"Proton")));
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"AntiProton")));
    }
    break;
  case 1://pion
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)))->Clone("RecoPion");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiPlus",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"PiMinus",i)));
      }
    }
    else{
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiPlus")))->Clone("RecoPion");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"PiMinus")));
    }
    break;
  case 2://kaon
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KPlus",i)))->Clone("RecoKaon");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KPlus",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"KMinus",i)));
      }
    }
    else{
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KPlus")))->Clone("RecoKaon");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"KMinus")));
    }
    break;
  case 3://proton
    if(cblast != -1){//add more centrality bins
      for(int i=cb;i<=cblast;i++){
	if(i==cb) numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"Proton",i)))->Clone("RecoProton");
	else{numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"Proton",i)));}
	numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%sCB%i",myname,"AntiProton",i)));
      }
    }
    else{
      numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"Proton")))->Clone("RecoProton");
      numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s",myname,"AntiProton")));
    }
    break;
  case 4://electron
    numeratorParent = (TH2F*)((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s%s",myname,"EPlus",cbname)))->Clone("RecoElectron");
    numeratorParent->Add((TH2F*) out2->FindObject(Form("EtNReconstructed%s%s%s",myname,"EMinus",cbname)));
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
  if(TPC) myname = mynameTPC;
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
  int lowbin = bkgd->GetYaxis()->FindBin(-etacut+.001);//make sure we don't accidentally get the wrong bin
  int highbin = bkgd->GetYaxis()->FindBin(etacut-.001);
  //cout<<"Projecting from "<<bkgd->GetYaxis()->GetBinLowEdge(lowbin)<<" to "<<bkgd->GetYaxis()->GetBinLowEdge(highbin+1)<<endl;


  TH1D *denominator = signal->ProjectionX(name,lowbin,highbin);
  TH1D *numerator = bkgd->ProjectionX("numerator",lowbin,highbin);
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
