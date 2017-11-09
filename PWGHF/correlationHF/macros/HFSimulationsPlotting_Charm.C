
//OFS templates; D <--charm
//Jitendra Kumar
//System: pp or p-Pb

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
Int_t nhistos=12; //pPb2016 = 32, pPb2013 = 6, pp 2013 = 9 - it's N_Dspecies*N_pT(D)
TString * filenames = new TString[nhistos];
TString genName = "Pythia8wBoost";

Int_t HFSimulationsPlotting_Charm(){
    
    gSystem->Exec("mkdir -p plots_charm/Basic/gif");
    gSystem->Exec("mkdir -p plots_charm/Specific/gif");
    gSystem->Exec("mkdir -p plots_charm/Specific/root");
    gSystem->Exec("mkdir -p plots_charm/Specific/png");
    gSystem->Sleep(100);
    
    //Step 1: Loading Libraries
    LoadLibraries();
    
    //Step 2: Correlations + thr options + name ?
    TString filename="AnalysisResults_pPb5TeV_Pythia8_Boost.root";
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
    
    TString objectoutputSpecific ="SpecificPythia8wBoost";
    TList *SimCorrSpecificlist = (TList*)Simulationsdirectory->Get(objectoutputSpecific);
    if(!SimCorrSpecificlist){
        cout<< " No Sim Specific Plots list  "<< Form(SimCorrSpecificlist) << " found, exiting... "<<endl;
    }
    
    
    for (Int_t pTasso = 1; pTasso <= 8; pTasso++){
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dzero", "3To5",  pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dzero", "5To8",  pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dzero", "8To16", pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dzero", "16To24", pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dplus", "3To5",  pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dplus", "5To8",  pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dplus", "8To16", pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dplus", "16To24", pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dstar", "3To5",  pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dstar", "5To8",  pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dstar", "8To16", pTasso, fSave); // PDGofDMeson = #
        CalculateDHadronCorrelations(SimCorrSpecificlist, "Dstar", "16To24", pTasso, fSave); // PDGofDMeson = #
    }
}

//______________________________________| Fatching Simulatios  task directory
void CalculateDHadronCorrelations(TList *CorrelationListtemp, TString DMeson, TString DMesonpTRange, Double_t pTasso, Bool_t SavePlot){
    
    
    TString pTRange = "";
    if(pTasso == 1)     pTRange = "ptAssall0.3to99.0";
    else if(pTasso == 2)pTRange = "ptAssall0.3to1.0";
    else if(pTasso == 3)pTRange = "ptAssall1.0to99.0";
    else if(pTasso == 4)pTRange = "ptAssall2.0to99.0";
    else if(pTasso == 5)pTRange = "ptAssall3.0to99.0";
    else if(pTasso == 6)pTRange = "ptAssall1.0to2.0";
    else if(pTasso == 7)pTRange = "ptAssall1.0to3.0";
    else if(pTasso == 8)pTRange = "ptAssall2.0to3.0";

    Double_t DpTmin = 0, DpTmax = 0;
    if(DMesonpTRange == "3To5")       DpTmin = 3.00, DpTmax = 5.00;
    else if(DMesonpTRange == "5To8")  DpTmin = 5.00, DpTmax = 8.00;
    else if(DMesonpTRange == "8To16") DpTmin = 8.00, DpTmax = 16.00;
    else if(DMesonpTRange == "16To24") DpTmin = 16.00, DpTmax = 24.00;
    else cout << "Select Proper pT of D Meson" << endl;
    
    TString DrwType = "surf1";
    TString FileName = "";
    FileName.Form("ppCorrelationPlots%sPt%sfromC%s_%s_DeltaEta10", genName.Data(), DMeson.Data(), DMesonpTRange.Data(), pTRange.Data());
    cout << FileName <<endl;
    
    
    //Trigger..
    THnSparseD *HFTriggertemp = (THnSparseD*)CorrelationListtemp->FindObject("HFTrgiggerProp");
    THnSparseD *HFTrigger = (THnSparseD*)HFTriggertemp->Clone(FileName.Data());
    if(!HFTrigger) {cout << "going back" << endl; return;}
    HFTrigger->GetAxis(1)->SetRangeUser(DpTmin*1.0001, DpTmax*0.9999);
    TH1D *HistNTrigger = (TH1D *)HFTrigger->Projection(0);
    HistNTrigger->SetName(FileName.Data());
    Int_t nTrg = 0;
    if (DMeson == "Dzero")      nTrg = (TH1D *)HistNTrigger->GetBinContent(9); //9 for charm
    else if (DMeson == "Dplus") nTrg = (TH1D *)HistNTrigger->GetBinContent(10); //10 for charm
    else if (DMeson == "Dstar") nTrg = (TH1D *)HistNTrigger->GetBinContent(11); //11 for charm
    else cout << "Number of Trigger is not found and its.. ->  " << nTrg << endl;
    cout << "IMPORTANT # of Trg are (" <<DMeson.Data() <<", " << DMesonpTRange.Data() << ") -------> " << nTrg << endl;
    
    
    
    //Correlations
    THnSparseD *HFCorrelations = (THnSparseD*)CorrelationListtemp->FindObject("2PCorrBtwn_HF-hadron");
    if(!HFCorrelations)return;
    THnSparseD *HFCorrelationsSoftPi =  (THnSparseD *)HFCorrelations->Clone(FileName.Data());
    THnSparseD *HFCorrelationsSecEl =   (THnSparseD *)HFCorrelations->Clone(FileName.Data());
    THnSparseD *HFCorrelationsCorr   =  (THnSparseD *)HFCorrelations->Clone(FileName.Data());
    
    //Correlations
    TCanvas *Corr2D = new TCanvas("cDeltaPhiEta","cDeltaPhiEta",800,800);
    TH2D *Correlations2D = ComputeCorrelations(HFCorrelationsCorr, DMeson.Data(), DMesonpTRange.Data(), pTasso , FileName);
    if(DMeson == "Dzero"){ //deactivate softpi removal, since we already just project for primary pi,K,p,e,mu!
        //Soft Pions Case of D0
//        TH2D *SoftCorrelations2D = SoftPionRemoval(HFCorrelationsSoftPi, DMeson.Data(), DMesonpTRange.Data(), pTasso , FileName);
//        Correlations2D->Add(SoftCorrelations2D, -1);
    }
    
    //Secondary electrons removal
    TH2D *SecElCorrelations2D = SecondaryElectronRemoval(HFCorrelationsSecEl, DMeson.Data(), DMesonpTRange.Data(), pTasso , FileName);
    Correlations2D->Add(SecElCorrelations2D, -1);
    
    Correlations2D->GetYaxis()->SetTitle("#Delta #eta");
    Correlations2D->GetYaxis()->SetTitleOffset(1.3);
    Correlations2D->GetXaxis()->SetTitle("#Delta #phi");
    Correlations2D->SetName("hCorrDeltaPhi2D");
    Correlations2D->SetTitle(FileName.Data());
    Correlations2D->Draw(DrwType);
    TString namefile_gif2D;
    namefile_gif2D.Form("./plots_charm/Specific/png/2D%s.png", FileName.Data());
    Corr2D->Print(namefile_gif2D);
    Corr2D->Close();
    
    
    Int_t  Eta_binmin =  Correlations2D->GetYaxis()->FindBin(-0.999);// Same for all
    Int_t  Eta_binmax =  Correlations2D->GetYaxis()->FindBin(+0.999); // " - "
    TCanvas *Corr1D = new TCanvas("cDeltaPhi","cDeltaPhi",800,800);
    
    TH1F* Correlations1D = (TH1F*)Correlations2D->ProjectionX(FileName, Eta_binmin, Eta_binmax);
    Correlations1D->GetYaxis()->SetTitle("entries");
    Correlations1D->GetYaxis()->SetTitleOffset(1.3);
    Correlations1D->GetXaxis()->SetTitle("#Delta #varphi");
    Correlations1D->SetName("hCorrDeltaPhi");
    Correlations1D->SetTitle(FileName.Data());
    Double_t BW = Correlations1D->GetBinWidth(5);
    Correlations1D->Scale(1.0/BW); //Bin width corrections..
    Correlations1D->Scale(1.0/nTrg);// Normalization to number of triggers.
    
    Correlations1D->Draw();
    TString namefile_gif1D;
    namefile_gif1D.Form("./plots_charm/Specific/png/%s.png", FileName.Data());
    Corr1D->Print(namefile_gif1D);
    Corr1D->Close();
    
    
    FileName += ".root";
    TFile *fOut=new TFile(FileName.Data(),"RECREATE");
    fOut->cd();
    Correlations2D->Write();
    Corr2D->Write();
    Correlations1D->Write();
    Corr1D->Write();
    Corr2D->Close();
    Corr1D->Close();
    fOut->Close();
    
    HFCorrelationsSoftPi->Delete();
    HFCorrelationsCorr->Delete();
    
}


//________________| Projection Function 1D
TH2D *ComputeCorrelations(THnSparse *sparse, const TString DMesonType, const TString DpTrange, const Double_t pTthrAss, const TString hname = "hnamePro"){
    
    
    Int_t Particle = 0;
    if(DMesonType == "Dzero")       Particle = 9;
    else if(DMesonType == "Dplus")  Particle = 10;
    else if(DMesonType == "Dstar")  Particle = 11;
    else cout << "Select Proper D Meson" << endl;
    
    
    Double_t DpTmin = 0, DpTmax = 0;
    if(DpTrange == "3To5")       DpTmin = 3.00 , DpTmax = 5.00;
    else if(DpTrange == "5To8")  DpTmin = 5.00 , DpTmax = 8.00;
    else if(DpTrange == "8To16") DpTmin = 8.00 , DpTmax = 16.00;
    else if(DpTrange == "16To24") DpTmin = 16.00 , DpTmax = 24.00;
    else cout << "Select Proper pT of D Meson" << endl;
    
    Double_t AssopTmin = 0., AssopTmax = 0.;
    if(pTthrAss == 1)       AssopTmin = 0.00 , AssopTmax = 5.00; //0 means 0.3 since there's an hardcoded lower pt cut! 5 is the upper edge of THnSparse for pTassoc axis!
    else if(pTthrAss == 2)  AssopTmin = 0.00 , AssopTmax = 1.00;
    else if(pTthrAss == 3)  AssopTmin = 1.00 , AssopTmax = 5.00;
    else if(pTthrAss == 4)  AssopTmin = 2.00 , AssopTmax = 5.00;
    else if(pTthrAss == 5)  AssopTmin = 3.00 , AssopTmax = 5.00;
    else if(pTthrAss == 6)  AssopTmin = 1.00 , AssopTmax = 2.00;
    else if(pTthrAss == 7)  AssopTmin = 1.00 , AssopTmax = 3.00;
    else if(pTthrAss == 8)  AssopTmin = 2.00 , AssopTmax = 3.00;
    else {cout << "Select Proper pT of Associated Tracks" << endl;}
    
    sparse->GetAxis(0)->SetRange(Particle, Particle);
    sparse->GetAxis(1)->SetRangeUser(DpTmin*1.001, DpTmax*0.999);
    sparse->GetAxis(3)->SetRangeUser(AssopTmin*1.0001, AssopTmax*0.9999);
    sparse->GetAxis(7)->SetRangeUser(3,10); //take only primary pi,K,p,e,mu, flaggeed as 1-5, so bins 3 to 7!
    
    TH2D *h = sparse->Projection(6,5);
    h->SetName(hname.Data());
    if(DMesonType == "Dzero"){
        
        h->SetMarkerColor(kRed);
        h->SetMarkerStyle(20);
        
    }else if(DMesonType == "Dplus"){
        
        h->SetMarkerColor(kGreen+3);
        h->SetMarkerStyle(21);
        
    }else if(DMesonType == "Dstar"){
        
        h->SetMarkerColor(kAzure-2);
        h->SetMarkerStyle(22);
        
    }
    
    h->SetMarkerSize(1.3);
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

TH2D *SoftPionRemoval(THnSparse *sparseSoftPi, const TString DMesonType, const TString DpTrange, const Double_t pTSoftPi, const TString hname = "hnamePro"){
       
    Double_t DpTmin = 0, DpTmax = 0;
    if(DpTrange == "3To5")       DpTmin = 3.00 , DpTmax = 5.00;
    else if(DpTrange == "5To8")  DpTmin = 5.00 , DpTmax = 8.00;
    else if(DpTrange == "8To16") DpTmin = 8.00 , DpTmax = 16.00;
    else if(DpTrange == "16To24") DpTmin = 16.00 , DpTmax = 24.00;
    else {cout << "Select Proper pT of D Meson" << endl;}
    
    Double_t AssopTmin = 0, AssopTmax = 0;
    if(pTSoftPi == 1)       AssopTmin = 0.00 , AssopTmax = 5.00;
    else if(pTSoftPi == 2)  AssopTmin = 0.00 , AssopTmax = 1.00;
    else if(pTSoftPi == 3)  AssopTmin = 1.00 , AssopTmax = 5.00;
    else if(pTSoftPi == 4)  AssopTmin = 2.00 , AssopTmax = 5.00;
    else if(pTSoftPi == 5)  AssopTmin = 3.00 , AssopTmax = 5.00;
    else if(pTSoftPi == 6)  AssopTmin = 1.00 , AssopTmax = 2.00;
    else if(pTSoftPi == 7)  AssopTmin = 1.00 , AssopTmax = 3.00;
    else if(pTSoftPi == 8)  AssopTmin = 2.00 , AssopTmax = 3.00;
    else {cout << "Select Proper pT of Assocated Tracks " << endl;}
    
    cout << "SoftPion: Assoc --> Min pT = " << AssopTmin<< " and Max pT =  " << AssopTmax<< endl;
    
    sparseSoftPi->GetAxis(0)->SetRange(9, 9); //D0
    sparseSoftPi->GetAxis(1)->SetRangeUser(DpTmin*1.001, DpTmax*0.999);
    sparseSoftPi->GetAxis(3)->SetRangeUser(AssopTmin*1.0001, AssopTmax*0.999);
    sparseSoftPi->GetAxis(7)->SetRange(1, 1); //soft poin
    
    TH2D *hsoftPi = sparseSoftPi->Projection(6,5);
    cout<<"Number of Soft Entries ------> = "<<hsoftPi->GetEntries()<<endl;
    
    hsoftPi->SetName(hname.Data());
    hsoftPi->SetMarkerColor(kRed);
    hsoftPi->SetMarkerStyle(20);
    hsoftPi->SetMarkerSize(1.3);
    return hsoftPi;
}


TH2D *SecondaryElectronRemoval(THnSparse *sparseSecEl, const TString DMesonType, const TString DpTrange, const Double_t pTSecEl, const TString hname = "hnamePro"){
       
    Double_t DpTmin = 0, DpTmax = 0;
    if(DpTrange == "3To5")       DpTmin = 3.00 , DpTmax = 5.00;
    else if(DpTrange == "5To8")  DpTmin = 5.00 , DpTmax = 8.00;
    else if(DpTrange == "8To16") DpTmin = 8.00 , DpTmax = 16.00;
    else if(DpTrange == "16To24") DpTmin = 16.00 , DpTmax = 24.00;
    else {cout << "Select Proper pT of D Meson" << endl;}
    
    Double_t AssopTmin = 0, AssopTmax = 0;
    if(pTSecEl == 1)       AssopTmin = 0.00 , AssopTmax = 5.00;
    else if(pTSecEl == 2)  AssopTmin = 0.00 , AssopTmax = 1.00;
    else if(pTSecEl == 3)  AssopTmin = 1.00 , AssopTmax = 5.00;
    else if(pTSecEl == 4)  AssopTmin = 2.00 , AssopTmax = 5.00;
    else if(pTSecEl == 5)  AssopTmin = 3.00 , AssopTmax = 5.00;
    else if(pTSecEl == 6)  AssopTmin = 1.00 , AssopTmax = 2.00;
    else if(pTSecEl == 7)  AssopTmin = 1.00 , AssopTmax = 3.00;
    else if(pTSecEl == 8)  AssopTmin = 2.00 , AssopTmax = 3.00;
    else {cout << "Select Proper pT of Assocated Tracks " << endl;}
    
    cout << "SecElon: Assoc --> Min pT = " << AssopTmin<< " and Max pT =  " << AssopTmax<< endl;
    
    sparseSecEl->GetAxis(0)->SetRange(9, 11); //D0
    sparseSecEl->GetAxis(1)->SetRangeUser(DpTmin*1.001, DpTmax*0.999);
    sparseSecEl->GetAxis(3)->SetRangeUser(AssopTmin*1.0001, AssopTmax*0.999);
    sparseSecEl->GetAxis(7)->SetRange(9, 9); //gamma conversion
    
    TH2D *hSecEl = sparseSecEl->Projection(6,5);
    cout<<"Number of Soft Entries ------> = "<<hSecEl->GetEntries()<<endl;
    
    hSecEl->SetName(hname.Data());
    hSecEl->SetMarkerColor(kRed);
    hSecEl->SetMarkerStyle(20);
    hSecEl->SetMarkerSize(1.3);
    return hSecEl;
}

// Plot comparison..

////_______________________________________________________________________
void LoadFileNamesDMeson(TString Origin="fromC", TString pTthrAss){
    
    
    filenames[0] = Form("ppCorrelationPlots%sPtDzero%s3To5_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    filenames[1] = Form("ppCorrelationPlots%sPtDplus%s3To5_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    filenames[2] = Form("ppCorrelationPlots%sPtDstar%s3To5_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    
    
    filenames[3] = Form("ppCorrelationPlots%sPtDzero%s5To8_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    filenames[4] = Form("ppCorrelationPlots%sPtDplus%s5To8_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    filenames[5] = Form("ppCorrelationPlots%sPtDstar%s5To8_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    
    
    filenames[6] = Form("ppCorrelationPlots%sPtDzero%s8To16_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    filenames[7] = Form("ppCorrelationPlots%sPtDplus%s8To16_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    filenames[8] = Form("ppCorrelationPlots%sPtDstar%s8To16_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());


    filenames[9] = Form("ppCorrelationPlots%sPtDzero%s16To24_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    filenames[10] = Form("ppCorrelationPlots%sPtDplus%s16To24_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
    filenames[11] = Form("ppCorrelationPlots%sPtDstar%s16To24_ptAssall%s_DeltaEta10.root",genName.Data(),Origin.Data(), pTthrAss.Data());
}


//_______________________________________________________________________
void CompareInSingleCanvas(TString Orgn= "NULL", TString pTth){
    
    
    LoadFileNamesDMeson(Orgn.Data(), pTth.Data());
    
    TString CanName1D = "";
    CanName1D.Form("1DCompare_allpT_CorrFor_AssoParticle_Pt_>_%sGeVc_%s", pTth.Data(),Orgn.Data());
    
    
    TH1D     *ThCorr[12];
    
    const Int_t i=4;
    TCanvas* cnew= new TCanvas(CanName1D.Data(), CanName1D.Data(),1200,400);
    cnew->Divide(4,1);
    
    for(Int_t f=0; f<12; f++){
        TString path = filenames[f];
        cout << "Reading File from path: " << path << endl;
        TFile * file = TFile::Open(path.Data(),"WRITE");
        ThCorr[f] = (TH1D*)file->Get("hCorrDeltaPhi");
        if(!ThCorr[f])cout << ":Something wrong" << endl;
        ThCorr[f]->SetMarkerSize(0.7);
    }
    
    
    cnew->cd(1);
    ThCorr[0]->Draw("p");
    //cout << "GetX bins " <<  ThCorr[1]->GetNbinsX() << endl; exit(1);
    ThCorr[0]->SetMaximum(ThCorr[0]->GetMaximum()*1.15);
    ThCorr[0]->SetMinimum(ThCorr[0]->GetMinimum()*0.85);
    ThCorr[1]->Draw("samep");
    ThCorr[2]->Draw("samep");
    
    cnew->cd(2);
    ThCorr[3]->Draw("p");
    ThCorr[3]->SetMaximum(ThCorr[3]->GetMaximum()*1.15);
    ThCorr[3]->SetMinimum(ThCorr[3]->GetMinimum()*0.85);
    ThCorr[4]->Draw("samep");
    ThCorr[5]->Draw("samep");
    
    cnew->cd(3);
    ThCorr[6]->Draw("p");
    ThCorr[6]->SetMaximum(ThCorr[6]->GetMaximum()*1.15);
    ThCorr[6]->SetMinimum(ThCorr[6]->GetMinimum()*0.85);
    ThCorr[7]->Draw("samep");
    ThCorr[8]->Draw("samep");

    cnew->cd(4);
    ThCorr[9]->Draw("p");
    ThCorr[9]->SetMaximum(ThCorr[9]->GetMaximum()*1.15);
    ThCorr[9]->SetMinimum(ThCorr[9]->GetMinimum()*0.85);
    ThCorr[10]->Draw("samep");
    ThCorr[11]->Draw("samep");    
    
    cnew->SaveAs(Form("plots_Compare/png/1DCompare_allDpT%s_AssoPt_%sGeVc_%s.png",Orgn.Data(), pTth.Data(),genName.Data()));
    
}


void DoCompareInSingleCanvas(){
    
    gSystem->Exec("mkdir -p plots_Compare/png");
    CompareInSingleCanvas("fromC","0.3to99.0");
    CompareInSingleCanvas("fromC","0.3to1.0");
    CompareInSingleCanvas("fromC","1.0to99.0");
    CompareInSingleCanvas("fromC","2.0to99.0");
    CompareInSingleCanvas("fromC","3.0to99.0");
    CompareInSingleCanvas("fromC","1.0to2.0");
    CompareInSingleCanvas("fromC","1.0to3.0");
    CompareInSingleCanvas("fromC","2.0to3.0");

    CompareInSingleCanvas("fromB","0.3to99.0");
    CompareInSingleCanvas("fromB","0.3to1.0");
    CompareInSingleCanvas("fromB","1.0to99.0");
    CompareInSingleCanvas("fromB","2.0to99.0");
    CompareInSingleCanvas("fromB","3.0to99.0");
    CompareInSingleCanvas("fromB","1.0to2.0");
    CompareInSingleCanvas("fromB","1.0to3.0");
    CompareInSingleCanvas("fromB","2.0to3.0");
    
}










