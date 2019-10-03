#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDirectoryFile.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH3F.h"
#include "TH2F.h"


//using namespace std::
TString directory = "";
const Int_t nhistos=2; //pPb !! IMPORTANT to change -->9 for pp
TString *filenames = new TString[nhistos];


Int_t PlottingHFCorrOnFlySim(){
  
  gSystem->Exec("mkdir -p plots/Basic/gif");
  gSystem->Exec("mkdir -p plots/Specific/gif");
  gSystem->Exec("mkdir -p plots/Specific/root");
  gSystem->Exec("mkdir -p plots/Specific/png");
  gSystem->Sleep(100);
  
  //Step 1: Loading Libraries
  LoadLibraries();
  
  //Step 2: Correlations + thr options + name ?
  TString filename="AnalysisResults_new3Dec.root";
  Bool_t Savingfiles= kTRUE; //Want to save your files ?
  
  DoCorreleations(filename, Savingfiles);  
  
}

//______________| Main function for correlation plots
void DoCorreleations(const char *infile="", Bool_t fSave){
  
  TFile* f = new TFile(infile,"READ");
  
  //______________________________________| Fatching Simulatios  task directory
  TDirectoryFile *Simulationsdirectory = (TDirectoryFile*)f->Get("KineSimulations");
  if(!Simulationsdirectory){
    cout<< " No Sim Corr directory "<< Form("KineSimulations") << " found, exiting... "<<endl;
  }
  
  TString objectoutputSpecific ="SpecificPerugia2011";
  TList *SimCorrSpecificlist = (TList*)Simulationsdirectory->Get(objectoutputSpecific);
  if(!SimCorrSpecificlist){
    cout<< " No Sim Specific Plots list  "<< Form(SimCorrSpecificlist) << " found, exiting... "<<endl;
  }
  
    
  TString objectoutputBasic ="BasicPlotsPerugia2011";
  TList *SimCorrBasiclist = (TList*)Simulationsdirectory->Get(objectoutputBasic);
  if(!SimCorrBasiclist){
        cout<< " No Sim Basic Plots list  "<< Form(SimCorrSpecificlist) << " found, exiting... "<<endl;
   }
   
    
  CalculateD0DStarBasicPlpots(SimCorrBasiclist); // PDGofDMeson = #

  for (Double_t pTasso = 0.0; pTasso < 0.5; pTasso+= 0.5){
    
   // CalculateDHadronCorrelations(SimCorrSpecificlist, "DZero", "mid",  pTasso, fSave); // PDGofDMeson = #
   // CalculateDHadronCorrelations(SimCorrSpecificlist, "DStar", "mid",  pTasso, fSave); // PDGofDMeson = #
    
  } 
   
}

void CalculateD0DStarBasicPlpots(TList *CorrelationListBase){

    //Correlations
    TCanvas *MultCanvas = new TCanvas("MultCanvas","Multiplicity",1000,800);
    MultCanvas->SetLogy();
    MultCanvas->cd();
    
    TH1F *MultTot = (TH1F*)CorrelationListBase->FindObject("fMultTot");
    //MultTot->GetYaxis()->SetTitle("entries");
    //MultTot->GetYaxis()->SetTitleOffset(1.3);
    //MultTot->GetXaxis()->SetTitle("Multiplictiy");
    //MultTot->SetName("Various Multiplicties");
    //MultTot->Draw();

    TH1D *fMultPhyPriPart = (TH1D*)CorrelationListBase->FindObject("fMultPhyPriPart");
    fMultPhyPriPart->SetMarkerColor(2);
    fMultPhyPriPart->SetMarkerSize(0.80);
    fMultPhyPriPart->GetYaxis()->SetTitle("entries");
    fMultPhyPriPart->GetYaxis()->SetTitleOffset(1.3);
    fMultPhyPriPart->GetXaxis()->SetTitle("Multiplictiy");
    fMultPhyPriPart->SetTitle("Various Multiplicties");
    fMultPhyPriPart->Draw("ep");

    TH1D *fMultNotPhyPriPart = (TH1D*)CorrelationListBase->FindObject("fMultNotPhyPriPart");
    fMultNotPhyPriPart->SetMarkerColor(4);
    fMultNotPhyPriPart->SetMarkerSize(0.80);
    fMultNotPhyPriPart->Draw("sameep");

    TH1D *fMultPhyPriPartCharge = (TH1D*)CorrelationListBase->FindObject("fMultPhyPriPartCharge");
    fMultPhyPriPartCharge->SetMarkerColor(6);
    fMultPhyPriPartCharge->SetMarkerSize(0.80);
    fMultPhyPriPartCharge->Draw("sameep");
    
    leg = new TLegend(0.65,0.74,0.98,0.94);
    leg->SetHeader("Multiplicity Distributions");
    leg->AddEntry(fMultPhyPriPart,"Physical Primaries","lep");
    leg->AddEntry(fMultNotPhyPriPart,"Not Physical Primaries","lep");
    leg->AddEntry(fMultPhyPriPartCharge,"Charge Physical Primaries","lep");
    leg->Draw();
    MultCanvas->Print("Multiplicities.png");
    MultCanvas->Print("Multiplicities.gif");
    
    

    THnSparseD *ParticlePro = (THnSparseD*)CorrelationListBase->FindObject("fHistNParticle");
    if(!ParticlePro) {cout << "going back" <endl; return;}
    
    //Quarks Level Plots
    // pT
    TCanvas *QuarkpTPlot = new TCanvas("QuarkpTPlot","QuarkpTPlot",1000,800);
    QuarkpTPlot->SetLogy();
    QuarkpTPlot->cd();
    TH1D *hLqpT = GetParticleProp(ParticlePro, "Lquarks", "pT" , "LightQuarkT");
    hLqpT->SetMarkerColor(2);
    hLqpT->Draw("ep");
    hLqpT->GetYaxis()->SetTitle("entries");
    hLqpT->GetYaxis()->SetTitleOffset(1.3);
    hLqpT->GetXaxis()->SetTitle("p_T");
    hLqpT->SetTitle("Quark pT Spectrum");
    
    TH1D *hcpT = GetParticleProp(ParticlePro, "4", "pT" , "CharmpT");
    hcpT->SetMarkerColor(4);
    hcpT->Draw("sameep");
    
    TH1D *hbpT = GetParticleProp(ParticlePro, "5", "pT" , "BeautypT");
    hbpT->SetMarkerColor(6);
    hbpT->Draw("sameep");
    
    legQpT = new TLegend(0.65,0.74,0.98,0.94);
    legQpT->SetHeader("pt Distributions");
    legQpT->AddEntry(hLqpT,"Light Quarks","lep");
    legQpT->AddEntry(hcpT,"Charm Quarks","lep");
    legQpT->AddEntry(hbpT,"Beauty Quarks","lep");
    legQpT->Draw();
    QuarkpTPlot->Print("QuarkpTSpectra.png");
    QuarkpTPlot->Print("QuarkpTSpectra.gif");
    
    
    
    TCanvas *ChargePartpTPlot = new TCanvas("ChargePartpTPlot","ChargePartpTPlot",1000,800);
    ChargePartpTPlot->SetLogy();
    ChargePartpTPlot->cd();
    TH1D *hPionpT = GetParticleProp(ParticlePro, "211pp", "pT" , "Pion");
    hPionpT->SetMarkerColor(2);
    hPionpT->Draw("ep");
    hPionpT->GetYaxis()->SetTitle("entries");
    hPionpT->GetYaxis()->SetTitleOffset(1.3);
    hPionpT->GetXaxis()->SetTitle("p_{T}");
    hPionpT->SetTitle("Charge pT Spectrum");
    
    TH1D *hKaonpT = GetParticleProp(ParticlePro, "321pp", "pT" , "Kaon");
    hKaonpT->SetMarkerColor(4);
    hKaonpT->Draw("sameep");
    
    TH1D *hProtonpT = GetParticleProp(ParticlePro, "2212pp", "pT" , "Proton");
    hProtonpT->SetMarkerColor(6);
    hProtonpT->Draw("sameep");
    
    TH1D *hElecpT = GetParticleProp(ParticlePro, "11pp", "pT" , "electron");
    hElecpT->SetMarkerColor(3);
    hElecpT->Draw("sameep");
    
    TH1D *hMuonpT = GetParticleProp(ParticlePro, "13pp", "pT" , "Muon");
    hMuonpT->SetMarkerColor(9);
    hMuonpT->Draw("sameep");
    
    legchargepT = new TLegend(0.68,0.62,0.94,0.94);
    legchargepT->SetHeader("Charge Particle pt Distributions");
    legchargepT->AddEntry(hPionpT,"Pion","lep");
    legchargepT->AddEntry(hKaonpT,"Kaon","lep");
    legchargepT->AddEntry(hProtonpT,"Protons","lep");
    legchargepT->AddEntry(hElecpT,"Electrons","lep");
    legchargepT->AddEntry(hMuonpT,"Muons","lep");
    legchargepT->Draw();
    ChargePartpTPlot->Print("ChargeParticlepT.png");
    ChargePartpTPlot->Print("ChargeParticlepT.gif");
    

    TCanvas *DMesonBMesonpT = new TCanvas("DMesonBMesonpT","DMesonBMesonpT",1000,800);
    DMesonBMesonpT->SetLogy();
    DMesonBMesonpT->cd();
    TH1D *hDzeropT = GetParticleProp(ParticlePro, "421", "pT" , "DzeropT");
    hDzeropT->SetMarkerColor(2);
    hDzeropT->Draw("ep");
    hDzeropT->GetYaxis()->SetTitle("entries");
    hDzeropT->GetYaxis()->SetTitleOffset(1.3);
    hDzeropT->GetXaxis()->SetTitle("p_T");
    hDzeropT->SetTitle("D and B pT Spectrum");
    
    TH1D *hDpluspT = GetParticleProp(ParticlePro, "411", "pT" , "DpluspT");
    hDpluspT->SetMarkerColor(4);
    hDpluspT->Draw("sameep");
    
    TH1D *hDstarpT = GetParticleProp(ParticlePro, "413", "pT" , "DstarpT");
    hDstarpT->SetMarkerColor(6);
    hDstarpT->Draw("sameep");
    
    TH1D *hBzeropT = GetParticleProp(ParticlePro, "511", "pT" , "BZero");
    hBzeropT->SetMarkerColor(3);
    hBzeropT->Draw("sameep");
    
    TH1D *hBpluspT = GetParticleProp(ParticlePro, "521", "pT" , "BPlus");
    hBpluspT->SetMarkerColor(9);
    hBpluspT->Draw("sameep");
    
    legDBpT = new TLegend(0.78,0.74,0.92,0.92);
    legDBpT->SetHeader("pt Distributions");
    legDBpT->AddEntry(hDzeropT,"DZero","lep");
    legDBpT->AddEntry(hDpluspT,"DPlus","lep");
    legDBpT->AddEntry(hDstarpT,"DStar","lep");
    legDBpT->AddEntry(hBzeropT,"BZero","lep");
    legDBpT->AddEntry(hBpluspT,"BPlus","lep");

    legDBpT->Draw();
    DMesonBMesonpT->Print("DandBMesonpT.png");
    DMesonBMesonpT->Print("DandBMesonpT.gif");
    
    
}

//________________| Projection Function 1D
TH1D *GetParticleProp(THnSparse *sparse, const TString PDGv="error", const TString KineVar="error", const TString hname="some"){
    

    Int_t ParticleType;
    if(PDGv=="Lquarks")     ParticleType = 1;
    else if(PDGv=="4")      ParticleType = 2;
    else if(PDGv=="5")      ParticleType = 3;
    else if(PDGv=="421")    ParticleType = 4;
    else if(PDGv=="411")    ParticleType = 5;
    else if(PDGv=="413")    ParticleType = 6;
    else if(PDGv=="511")    ParticleType = 7;
    else if(PDGv=="521")    ParticleType = 8;
    else if(PDGv=="211pp")  ParticleType = 9;
    else if(PDGv=="321pp")  ParticleType = 10;
    else if(PDGv=="2212pp") ParticleType = 11;
    else if(PDGv=="11pp")   ParticleType = 12;
    else if(PDGv=="13pp")   ParticleType = 13;
    else if(PDGv=="211np")  ParticleType = 14;
    else if(PDGv=="321np")  ParticleType = 15;
    else if(PDGv=="11np")   ParticleType = 16;
    else if(PDGv=="13np")   ParticleType = 17;
    else if(PDGv=="np")     ParticleType = 18;
    else if(PDGv=="rest")   ParticleType = 19;
    else {Printf("Error in Particle Selection.. exiting..");  exit(1);}
    
    Int_t VariableType;
    if(KineVar=="pT")       VariableType = 1;
    else if(KineVar=="eta") VariableType = 2;
    else if(KineVar=="phi") VariableType = 3;
    
    THnSparse *sparseClone  = (THnSparse *)sparse->Clone(hname.Data());
    sparseClone->GetAxis(0)->SetRange(ParticleType, ParticleType);

    TH1D *h = sparseClone->Projection(VariableType);
    h->SetName(hname.Data());
    if(KineVar == "pT")h->SetMarkerStyle(20);
    else if(KineVar == "eta")h->SetMarkerStyle(21);
    else if(KineVar == "phi")h->SetMarkerStyle(22);
    h->SetMarkerSize(0.90);
    return h;
    
}


//________________| Loading Libraries
void LoadLibraries() {
    
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/CDB -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_ROOT/PWGHF/base -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/PWGHF/vertexingHF -I$ALICE_ROOT/PWGHF/correlationHF -I$ALICE_ROOT/PWG/FLOW/Base -I$ALICE_ROOT/PWG/FLOW/Tasks -g");
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libMinuit.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libCORRFW.so");
    gSystem->Load("libPWGflowBase.so");
    gSystem->Load("libPWGflowTasks.so");
    
    //gSystem->Load("libPWGmuon.so");
    return;  
    
}











