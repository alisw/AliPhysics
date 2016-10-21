/*
 rirusso@cern.ch - last update on 16/02/2014
 Macro to run the ITS QA trending by accessing the std QA output,
 to be mainly used with the automatic scripts to fill the QA repository.
 Launch with
 aliroot -l -b -q "MakeTrendingITSQA.C(\"${fullpath}/QAresults.root\", ${run}, ...)
 The macro produces a file containing the tree of trending variables and the main plots.
 */
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TGrid.h>
#include <TGridResult.h>
#include <TMath.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <Riostream.h>
#include <iostream>
#include <fstream>
//#include <AliITSgeomTGeo.h>
#endif






Double_t LangausFun(Double_t *x, Double_t *par) ;

Int_t MakeTrendingITSQA(TString qafilename ="QAresults.root",       //full path of the QA output; set IsOnGrid to prepend "alien://"
                        Int_t runNumber= 133005,          // run number
                        Bool_t isMC=kFALSE,       //MC flag, to disable meaningless checks
                        Bool_t canvasE = kFALSE,  //enable display plots on canvas and save png
                        Bool_t IsOnGrid = kFALSE, //set to kTRUE to access files on the grid
                        TString ocdbStorage = "raw://") //set the default ocdb storage
{
    // macro to generate tree with ITS QA trending variables
    // access qa PWGPP output files
    if (!qafilename) {
        Printf("Error - Invalid input file");
        return 1;
    }



    char defaultQAoutput[30]="QAresults.root";
    char * treePostFileName="trending.root";

    if (IsOnGrid) TGrid::Connect("alien://");
    TFile * fin = TFile::Open(qafilename,"r");
    if (!fin) {
        Printf("ERROR: QA output not found. Exiting...\n");
        return -1;
    } else {
        Printf("INFO: QA output file %s open. \n",fin->GetName());
    }

    TFile * trendFile = new TFile(treePostFileName,"recreate");


    ///// SDD Variables
    
    Int_t nrun,nEvents, nEventsTriggered;
    Float_t minDrTime,errminDrTime,meanDrTime,errmeanDrTime;
    Float_t fracTrackWithClu1,fracTrackWithClu2,errfracTrackWithClu1,errfracTrackWithClu2;
    Float_t fracTrackWithClu3,fracTrackWithClu4,errfracTrackWithClu3,errfracTrackWithClu4;
    Float_t fracTrackWithClu5,fracTrackWithClu6,errfracTrackWithClu5,errfracTrackWithClu6;
    Float_t fracExtra,errfracExtra;
    Float_t fracT[6]={0.,0.,0.,0.,0.,0.};
    Float_t efracT[6]={0.,0.,0.,0.,0.,0.};
    Int_t nTotEvents;
    Int_t nTrigEvents;
    Double_t averPoints=0.;
    Double_t cntBins=0.;
    Float_t minTime=-999.;
    Float_t errMinTime=0.;
    Float_t MPVdEdxLay3,errMPVdEdxLay3,MPVdEdxLay4,errMPVdEdxLay4;
    Float_t MPVdEdxTB0,errMPVdEdxTB0,MPVdEdxTB5,errMPVdEdxTB5;

    
    //// Vertex Variables
    
    Float_t meanVtxTRKx,meanVtxTRKy,meanVtxTRKz;
    Float_t meanVtxSPDx,meanVtxSPDy,meanVtxSPDz;
    Float_t sigmaVtxTRKx,sigmaVtxTRKy,sigmaVtxTRKz;
    Float_t sigmaVtxSPDx,sigmaVtxSPDy,sigmaVtxSPDz;
    Float_t meanVtxTRKxErr,meanVtxTRKyErr,meanVtxTRKzErr;
    Float_t meanVtxSPDxErr,meanVtxSPDyErr,meanVtxSPDzErr;
    Float_t sigmaVtxTRKxErr,sigmaVtxTRKyErr,sigmaVtxTRKzErr;
    Float_t sigmaVtxSPDxErr,sigmaVtxSPDyErr,sigmaVtxSPDzErr;
    
    
    ///// SSD Variables
    
    
    Float_t MPVL5,MPVErrL5;
    Float_t MPVL6,MPVErrL6;
    Float_t ChargeRatioL5,ChargeRatioErrL5;
    Float_t ChargeRatioL6,ChargeRatioErrL6;
    Float_t EmptyModulesSDD;


    ///// Matching Variables
    


    Float_t FracSPD1;
    Float_t errFracSPD1;
    Float_t FracSPD2;
    Float_t errFracSPD2;
    Float_t Eff6Pt02;
    Float_t errEff6Pt02;
    Float_t Eff6Pt1;
    Float_t errEff6Pt1;
    Float_t Eff6Pt10;
    Float_t errEff6Pt10;
    Float_t Eff5Pt02;
    Float_t errEff5Pt02;
    Float_t Eff5Pt1;
    Float_t errEff5Pt1;
    Float_t Eff5Pt10;
    Float_t errEff5Pt10;
    Float_t Eff4Pt02;
    Float_t errEff4Pt02;
    Float_t Eff4Pt1;
    Float_t errEff4Pt1;
    Float_t Eff4Pt10;
    Float_t errEff4Pt10;
    Float_t Eff3Pt02;
    Float_t errEff3Pt02;
    Float_t Eff3Pt1;
    Float_t errEff3Pt1;
    Float_t Eff3Pt10;
    Float_t errEff3Pt10;
    Float_t Eff2Pt02;
    Float_t errEff2Pt02;
    Float_t Eff2Pt1;
    Float_t errEff2Pt1;
    Float_t Eff2Pt10;
    Float_t errEff2Pt10;
    Float_t EffSPDPt02;
    Float_t errEffSPDPt02;
    Float_t EffSPDPt1;
    Float_t errEffSPDPt1;
    Float_t EffSPDPt10;
    Float_t errEffSPDPt10;
    Float_t EffoneSPDPt02;
    Float_t errEffoneSPDPt02;
    Float_t EffoneSPDPt1;
    Float_t errEffoneSPDPt1;
    Float_t EffoneSPDPt10;
    Float_t errEffoneSPDPt10;
    Float_t EffTOTPt02;
    Float_t errEffTOTPt02;
    Float_t EffTOTPt1;
    Float_t errEffTOTPt1;
    Float_t EffTOTPt10;
    Float_t errEffTOTPt10;
    
    

    
    
    TTree * ttree=new TTree("trending","tree of trending variables");
    ttree->Branch("nrun",&nrun,"nrun/I");
    ttree->Branch("nEvents",&nEvents,"nEvents/I");
    ttree->Branch("nEventsTriggered",&nEventsTriggered,"nEventsTriggered/I");
    ttree->Branch("minDrTime",&minDrTime,"minDrTime/F");
    ttree->Branch("errminDrTime",&errminDrTime,"errminDrTime/F");
    ttree->Branch("meanDrTime",&meanDrTime,"meanDrTime/F");
    ttree->Branch("errmeanDrTime",&errmeanDrTime,"errmeanDrTime/F"); //mean time
    ttree->Branch("fracTrackWithClu1",&fracTrackWithClu1,"fracTrackWithClu1/F"); //fraction of tracks with cluster in layer 1

    ttree->Branch("errfracTrackWithClu1",&errfracTrackWithClu1,"errfracTrackWithClu1/F"); //error fraction of tracks with cluster in layer 1
    ttree->Branch("fracTrackWithClu2",&fracTrackWithClu2,"fracTrackWithClu2/F"); //
    ttree->Branch("errfracTrackWithClu2",&errfracTrackWithClu2,"errfracTrackWithClu2/F"); //
    ttree->Branch("fracTrackWithClu3",&fracTrackWithClu3,"fracTrackWithClu3/F");
    ttree->Branch("errfracTrackWithClu3",&errfracTrackWithClu3,"errfracTrackWithClu3/F");
    ttree->Branch("fracTrackWithClu4",&fracTrackWithClu4,"fracTrackWithClu4/F");
    ttree->Branch("errfracTrackWithClu4",&errfracTrackWithClu4,"errfracTrackWithClu4/F");
    ttree->Branch("fracTrackWithClu5",&fracTrackWithClu5,"fracTrackWithClu5/F");
    ttree->Branch("errfracTrackWithClu5",&errfracTrackWithClu5,"errfracTrackWithClu5/F");
    ttree->Branch("fracTrackWithClu6",&fracTrackWithClu6,"fracTrackWithClu6/F");
    ttree->Branch("errfracTrackWithClu6",&errfracTrackWithClu6,"errfracTrackWithClu6/F");


    
    ttree->Branch("meanVtxTRKx",&meanVtxTRKx,"meanVtxTRKx/F"); // mean of tracks vertex position - x
    ttree->Branch("meanVtxTRKy",&meanVtxTRKy,"meanVtxTRKy/F"); // mean of tracks vertex position - y
    ttree->Branch("meanVtxTRKz",&meanVtxTRKz,"meanVtxTRKz/F"); // mean of tracks vertex position - z
    ttree->Branch("meanVtxTRKxErr",&meanVtxTRKxErr,"meanVtxTRKxErr/F"); // error mean of tracks vertex position - x
    ttree->Branch("meanVtxTRKyErr",&meanVtxTRKyErr,"meanVtxTRKyErr/F"); // error mean of tracks vertex position - y
    ttree->Branch("meanVtxTRKzErr",&meanVtxTRKzErr,"meanVtxTRKzErr/F"); // error mean of tracks vertex position - z

    ttree->Branch("meanVtxSPDx",&meanVtxSPDx,"meanVtxSPDx/F"); // mean of SPD vertex position - x
    ttree->Branch("meanVtxSPDy",&meanVtxSPDy,"meanVtxSPDy/F"); // mean of SPD vertex position - y
    ttree->Branch("meanVtxSPDz",&meanVtxSPDz,"meanVtxSPDz/F"); // mean of SPD vertex position - z
    ttree->Branch("meanVtxSPDxErr",&meanVtxSPDxErr,"meanVtxSPDxErr/F"); // error mean of SPD vertex position - x
    ttree->Branch("meanVtxSPDyErr",&meanVtxSPDyErr,"meanVtxSPDyErr/F"); // error mean of SPD vertex position - y
    ttree->Branch("meanVtxSPDzErr",&meanVtxSPDzErr,"meanVtxSPDzErr/F"); // error mean of SPD vertex position - z

    ttree->Branch("sigmaVtxTRKx",&sigmaVtxTRKx,"sigmaVtxTRKx/F"); // sigma of tracks vertex position - x
    ttree->Branch("sigmaVtxTRKy",&sigmaVtxTRKy,"sigmaVtxTRKy/F"); // sigma of tracks vertex position - y
    ttree->Branch("sigmaVtxTRKz",&sigmaVtxTRKz,"sigmaVtxTRKz/F"); // sigma of tracks vertex position - z
    ttree->Branch("sigmaVtxTRKxErr",&sigmaVtxTRKxErr,"sigmaVtxTRKxErr/F"); // error sigma of tracks vertex position - x
    ttree->Branch("sigmaVtxTRKyErr",&sigmaVtxTRKyErr,"sigmaVtxTRKyErr/F"); // error sigma of tracks vertex position - y
    ttree->Branch("sigmaVtxTRKzErr",&sigmaVtxTRKzErr,"sigmaVtxTRKzErr/F"); // error sigma of tracks vertex position - z

    ttree->Branch("sigmaVtxSPDx",&sigmaVtxSPDx,"sigmaVtxSPDx/F"); // sigma of tracks vertex position - x
    ttree->Branch("sigmaVtxSPDy",&sigmaVtxSPDy,"sigmaVtxSPDy/F"); // sigma of tracks vertex position - y
    ttree->Branch("sigmaVtxSPDz",&sigmaVtxSPDz,"sigmaVtxSPDz/F"); // sigma of tracks vertex position - z
    ttree->Branch("sigmaVtxSPDxErr",&sigmaVtxSPDxErr,"sigmaVtxSPDxErr/F"); // error sigma of tracks vertex position - x
    ttree->Branch("sigmaVtxSPDyErr",&sigmaVtxSPDyErr,"sigmaVtxSPDyErr/F"); // error sigma of tracks vertex position - y
    ttree->Branch("sigmaVtxSPDzErr",&sigmaVtxSPDzErr,"sigmaVtxSPDzErr/F"); // error sigma of tracks vertex position - z
    
    
    
    ttree->Branch("MPVL5",&MPVL5,"MPVL5/F"); // Most Probable Value dEdx Layer 5
    ttree->Branch("MPVErrL5",&MPVErrL5,"MPVErrL5/F"); // Most Probable Value error dEdx Layer 5

    ttree->Branch("MPVL6",&MPVL6,"MPVL6/F"); // Most Probable Value dEdx Layer 6
    ttree->Branch("MPVErrL6",&MPVErrL6,"MPVErrL6/F"); // Most Probable Value error dEdx Layer 6

    ttree->Branch("ChargeRatioL5",&ChargeRatioL5,"ChargeRatioL5/F"); // Charge ratio (2 sides of SSD) Layer 5
    ttree->Branch("ChargeRatioErrL5",&ChargeRatioErrL5,"ChargeRatioErrL5/F"); // Charge ratio error (2 sides of SSD) Layer 5

    ttree->Branch("ChargeRatioL6",&ChargeRatioL6,"ChargeRatioL6/F"); // Charge ratio (2 sides of SSD) Layer 6
    ttree->Branch("ChargeRatioErrL6",&ChargeRatioErrL6,"ChargeRatioErrL6/F"); // Charge ratio error(2 sides of SSD) Layer 6

    ttree->Branch("EmptyModulesSDD",&EmptyModulesSDD,"EmptyModulesSDD/F"); // Number of empty SSD  modules


    
    ttree->Branch("fracExtra",&fracExtra,"fracExtra/F"); // fraction of extra clusters in SDD
    ttree->Branch("errfracExtra",&errfracExtra,"errfracExtra/F"); // fraction of extra clusters in SDD
    
    ttree->Branch("minTime",&minTime,"minTime/F"); // minimum drift time SDD
    ttree->Branch("errMinTime",&errMinTime,"errMinTime/F"); //  error on minimum drift time SDD
    
    ttree->Branch("MPVdEdxLay3",&MPVdEdxLay3,"MPVdEdxLay3/F"); // most probable value of dE/Fx distribution of SDD Layer 3
    ttree->Branch("errMPVdEdxLay3",&errMPVdEdxLay3,"errMPVdEdxLay3/F"); // error  most probable value of dE/Fx distribution of SDD Layer 3

    ttree->Branch("MPVdEdxLay4",&MPVdEdxLay4,"MPVdEdxLay4/F"); // most probable value of dE/Fx distribution of SDD Layer 4
    ttree->Branch("errMPVdEdxLay4",&errMPVdEdxLay4,"errMPVdEdxLay4/F"); // error  most probable value of dE/Fx distribution of SDD Layer 4

    ttree->Branch("MPVdEdxTB0",&MPVdEdxTB0,"MPVdEdxTB0/F"); // most probable value of dE/Fx distribution of SDD - small drift time
    ttree->Branch("errMPVdEdxTB0",&errMPVdEdxTB0,"errMPVdEdxTB0/F"); // most probable value of dE/Fx distribution of SDD - small drift time

    ttree->Branch("MPVdEdxTB5",&MPVdEdxTB5,"MPVdEdxTB5/F"); // most probable value of dE/Fx distribution of SDD - large drift time
    ttree->Branch("errMPVdEdxTB5",&errMPVdEdxTB5,"errMPVdEdxTB5/F"); // most probable value of dE/Fx distribution of SDD - large drift time


    
    
    ttree->Branch("Eff6Pt02",&Eff6Pt02,"Eff6Pt02/F"); // matching efficiency low pt 6 clusters
    ttree->Branch("errEff6Pt02",&errEff6Pt02,"errEff6Pt02/F"); // error matching efficiency low pt 6 clusters
    ttree->Branch("Eff5Pt02",&Eff5Pt02,"Eff5Pt02/F"); // matching efficiency low pt 5 clusters
    ttree->Branch("errEff5Pt02",&errEff5Pt02,"errEff5Pt02/F"); // error matching efficiency low pt 5 clusters
    ttree->Branch("Eff4Pt02",&Eff4Pt02,"Eff4Pt02/F"); // matching efficiency low pt 4 clusters
    ttree->Branch("errEff4Pt02",&errEff4Pt02,"errEff4Pt02/F"); // error matching efficiency low pt 4 clusters
    ttree->Branch("Eff3Pt02",&Eff3Pt02,"Eff3Pt02/F"); // matching efficiency low pt 3 clusters
    ttree->Branch("errEff3Pt02",&errEff3Pt02,"errEff3Pt02/F"); // error matching efficiency low pt 3 clusters
    ttree->Branch("Eff2Pt02",&Eff2Pt02,"Eff2Pt02/F"); // matching efficiency low pt 2 clusters
    ttree->Branch("errEff2Pt02",&errEff2Pt02,"errEff2Pt02/F"); // error matching efficiency low pt 2 clusters
    ttree->Branch("EffSPDPt02",&EffSPDPt02,"EffSPDPt02/F"); // matching efficiency low pt 2 SPD
    ttree->Branch("errEffSPDPt02",&errEffSPDPt02,"errEffSPDPt02/F"); // error matching efficiency low pt 2 SPD
    ttree->Branch("EffoneSPDPt02",&EffoneSPDPt02,"EffoneSPDPt02/F"); // matching efficiency low pt 6 one SPD
    ttree->Branch("errEffoneSPDPt02",&errEffoneSPDPt02,"errEffoneSPDPt02/F"); // error matching efficiency low pt one SPD
    ttree->Branch("EffTOTPt02",&EffTOTPt02,"EffTOTPt02/F"); // matching efficiency low pt
    ttree->Branch("errEffTOTPt02",&errEffTOTPt02,"errEffTOTPt02/F"); // error matching efficiency low pt
    
    

    
    
    ttree->Branch("Eff6Pt1",&Eff6Pt1,"Eff6Pt1/F"); // matching efficiency mid pt 6 clusters
    ttree->Branch("errEff6Pt1",&errEff6Pt1,"errEff6Pt1/F"); // error matching efficiency mid pt 6 clusters
    ttree->Branch("Eff5Pt1",&Eff5Pt1,"Eff5Pt1/F"); // matching efficiency mid pt 5 clusters
    ttree->Branch("errEff5Pt1",&errEff5Pt1,"errEff5Pt1/F"); // error matching efficiency mid pt 5 clusters
    ttree->Branch("Eff4Pt1",&Eff4Pt1,"Eff4Pt1/F"); // matching efficiency mid pt 4 clusters
    ttree->Branch("errEff4Pt1",&errEff4Pt1,"errEff4Pt1/F"); // error matching efficiency mid pt 4 clusters
    ttree->Branch("Eff3Pt1",&Eff3Pt1,"Eff3Pt1/F"); // matching efficiency mid pt 3 clusters
    ttree->Branch("errEff3Pt1",&errEff3Pt1,"errEff3Pt1/F"); // error matching efficiency mid pt 3 clusters
    ttree->Branch("Eff2Pt1",&Eff2Pt1,"Eff2Pt1/F"); // matching efficiency mid pt 2 clusters
    ttree->Branch("errEff2Pt1",&errEff2Pt1,"errEff2Pt1/F"); // error matching efficiency mid pt 2 clusters
    ttree->Branch("EffSPDPt1",&EffSPDPt1,"EffSPDPt1/F"); // matching efficiency mid pt 2 SPD
    ttree->Branch("errEffSPDPt1",&errEffSPDPt1,"errEffSPDPt1/F"); // error matching efficiency mid pt 2 SPD
    ttree->Branch("EffoneSPDPt1",&EffoneSPDPt1,"EffoneSPDPt1/F"); // matching efficiency mid pt 6 one SPD
    ttree->Branch("errEffoneSPDPt1",&errEffoneSPDPt1,"errEffoneSPDPt1/F"); // error matching efficiency mid pt one SPD
    ttree->Branch("EffTOTPt1",&EffTOTPt1,"EffTOTPt1/F"); // matching efficiency mid pt
    ttree->Branch("errEffTOTPt1",&errEffTOTPt1,"errEffTOTPt1/F"); // error matching efficiency mid pt




    ttree->Branch("Eff6Pt10",&Eff6Pt10,"Eff6Pt10/F"); // matching efficiency high pt 6 clusters
    ttree->Branch("errEff6Pt10",&errEff6Pt10,"errEff6Pt10/F"); // error matching efficiency high pt 6 clusters
    ttree->Branch("Eff5Pt10",&Eff5Pt10,"Eff5Pt10/F"); // matching efficiency high pt 5 clusters
    ttree->Branch("errEff5Pt10",&errEff5Pt10,"errEff5Pt10/F"); // error matching efficiency high pt 5 clusters
    ttree->Branch("Eff4Pt10",&Eff4Pt10,"Eff4Pt10/F"); // matching efficiency high pt 4 clusters
    ttree->Branch("errEff4Pt10",&errEff4Pt10,"errEff4Pt10/F"); // error matching efficiency high pt 4 clusters
    ttree->Branch("Eff3Pt10",&Eff3Pt10,"Eff3Pt10/F"); // matching efficiency high pt 3 clusters
    ttree->Branch("errEff3Pt10",&errEff3Pt10,"errEff3Pt10/F"); // error matching efficiency high pt 3 clusters
    ttree->Branch("Eff2Pt10",&Eff2Pt10,"Eff2Pt10/F"); // matching efficiency high pt 2 clusters
    ttree->Branch("errEff2Pt10",&errEff2Pt10,"errEff2Pt10/F"); // error matching efficiency high pt 2 clusters
    ttree->Branch("EffSPDPt10",&EffSPDPt10,"EffSPDPt10/F"); // matching efficiency high pt 2 SPD
    ttree->Branch("errEffSPDPt10",&errEffSPDPt10,"errEffSPDPt10/F"); // error matching efficiency high pt 2 SPD
    ttree->Branch("EffoneSPDPt10",&EffoneSPDPt10,"EffoneSPDPt10/F"); // matching efficiency high pt 6 one SPD
    ttree->Branch("errEffoneSPDPt10",&errEffoneSPDPt10,"errEffoneSPDPt10/F"); // error matching efficiency high pt one SPD
    ttree->Branch("EffTOTPt10",&EffTOTPt10,"EffTOTPt10/F"); // matching efficiency high pt
    ttree->Branch("errEffTOTPt10",&errEffTOTPt10,"errEffTOTPt10/F"); // error matching efficiency high pt

    
    
    ttree->Branch("FracSPD1",&FracSPD1,"FracSPD1/F"); // fraction SPD layers active on 1 layer
    ttree->Branch("errFracSPD1",&errFracSPD1,"errFracSPD1/F");

    ttree->Branch("FracSPD2",&FracSPD2,"FracSPD2/F"); // fraction SPD layers active on 1 layer
    ttree->Branch("errFracSPD2",&errFracSPD2,"errFracSPD2/F");
    
    
    ///////////////   Vertex part
    nrun=runNumber;
    
    
    char VertexDirName[25]="Vertex_Performance";
    
    char VertexListName[25]="cOutputVtxESD";
    
    TDirectoryFile * VertexQAdir=(TDirectoryFile*)fin->Get(VertexDirName);
    if (!VertexQAdir) {
        Printf("ERROR: Vertex QA directory not present in input file.\n");
        return -1;
    }
    TList * VertxList=(TList*)VertexQAdir->Get(VertexListName);
    
    
    if (!VertxList) Printf("WARNING: Vertex QA histograms absent or not accessible\n");
    
    
    
    
    Printf("Vertex - QA");
    
    
    Int_t iRun=runNumber;

    TH1F *xVtxTRK = (TH1F*)VertxList->FindObject("fhTRKVertexX");
    
    
    TH1F *yVtxTRK = (TH1F*)VertxList->FindObject("fhTRKVertexY");
    

    
    TH1F *zVtxTRK = (TH1F*)VertxList->FindObject("fhTRKVertexZ");
    
    
    TH1F *xVtxSPD = (TH1F*)VertxList->FindObject("fhSPDVertexX");
    
    if(xVtxSPD->GetEntries()==0){
        printf("Run %d xVtxSOD EMPTY -> Return\n",iRun);
        

    }
    
    TH1F *yVtxSPD = (TH1F*)VertxList->FindObject("fhSPDVertexY");
    
    if(yVtxSPD->GetEntries()==0){
        printf("Run %d yVtxSPD EMPTY -> Return\n",iRun);
        

    }
    
    TH1F *zVtxSPD = (TH1F*)VertxList->FindObject("fhSPDVertexZ");
    
    if(zVtxSPD->GetEntries()==0){
        printf("Run %d zVtxSPD EMPTY -> Return\n",iRun);
        

    }
    
    TF1 *fxTRK = new TF1("gausx", "gaus", -1, 1);
    xVtxTRK->Fit("gausx", "NQRL");
    
    TF1 *fyTRK = new TF1("gausy", "gaus", -1, 1);
    yVtxTRK->Fit("gausy","NQLR");
    
    TF1 *fzTRK = new TF1("gausz", "gaus", -1, 1);
    zVtxTRK->Fit("gausz","NQRL");
    TF1 *fxSPD = new TF1("gausxSPD", "gaus", -1, 1);
    xVtxSPD->Fit("gausxSPD", "NQRL");
    
    TF1 *fySPD = new TF1("gausySPD", "gaus", -1, 1);
    yVtxSPD->Fit("gausySPD","NQLR");
    
    TF1 *fzSPD = new TF1("gauszSPD", "gaus", -1, 1);
    zVtxSPD->Fit("gauszSPD","NQRL");


    meanVtxTRKx=(Float_t)fxTRK->GetParameter(1);
    meanVtxTRKxErr=(Float_t)fxTRK->GetParError(1);
    sigmaVtxTRKx=(Float_t)fxTRK->GetParameter(2);
    sigmaVtxTRKxErr=(Float_t)fxTRK->GetParError(2);
    meanVtxTRKy=(Float_t)fyTRK->GetParameter(1);
    meanVtxTRKyErr=(Float_t)fyTRK->GetParError(1);
    sigmaVtxTRKy=(Float_t)fyTRK->GetParameter(2);
    sigmaVtxTRKyErr=(Float_t)fyTRK->GetParError(2);
    meanVtxTRKz=(Float_t)fzTRK->GetParameter(1);
    meanVtxTRKzErr=(Float_t)fzTRK->GetParError(1);
    sigmaVtxTRKz=(Float_t)fzTRK->GetParameter(2);
    sigmaVtxTRKzErr=(Float_t)fzTRK->GetParError(2);
    meanVtxSPDx=(Float_t)fxSPD->GetParameter(1);
    meanVtxSPDxErr=(Float_t)fxSPD->GetParError(1);
    sigmaVtxSPDx=(Float_t)fxSPD->GetParameter(2);
    sigmaVtxSPDxErr=(Float_t)fxSPD->GetParError(2);
    meanVtxSPDy=(Float_t)fySPD->GetParameter(1);
    meanVtxSPDyErr=(Float_t)fySPD->GetParError(1);
    sigmaVtxSPDy=(Float_t)fySPD->GetParameter(2);
    sigmaVtxSPDyErr=(Float_t)fySPD->GetParError(2);
    meanVtxSPDz=(Float_t)fzSPD->GetParameter(1);
    sigmaVtxSPDzErr=(Float_t)fzSPD->GetParError(1);
    sigmaVtxSPDz=(Float_t)fzSPD->GetParameter(2);
    sigmaVtxSPDzErr=(Float_t)fzSPD->GetParError(2);
    
    
    /////////// end of vertex part
    
    
    ///////////////////////  SSD Part
    char SSDDirName[25]="PWGPPdEdxSSDQA";
    
    char SSDListName[25]="SSDdEdxQA";
    
    TDirectoryFile * SSDQAdir=(TDirectoryFile*)fin->Get(SSDDirName);
    if (!SSDQAdir) {
        Printf("ERROR: SSD QA directory not present in input file.\n");
        return -1;
    }
    TList * SSDList=(TList*)SSDQAdir->Get(SSDListName);
    
    
    if (!SSDList) Printf("WARNING: SSD QA histograms absent or not accessible\n");
    
    
    
    
    Printf("SSD - QA");
    

    MPVL5=0;
    MPVErrL5=0;
    MPVL6=0;
    MPVErrL6=0;
    ChargeRatioL5=0;
    ChargeRatioErrL5=0;
    ChargeRatioL6=0;
    ChargeRatioErrL6=0;
    EmptyModulesSDD=0;
    
    
    TH2F* QAchargeRatio=(TH2F*)SSDList->FindObject("QAChargeRatio");
    
    if(QAchargeRatio->GetEntries()==0){
        printf("Run %d QAchargeRatio EMPTY -> Return\n",iRun);
    }
    
    TH2F* QAcharge=(TH2F*)SSDList->FindObject("QACharge");
    
    if(QAcharge->GetEntries()==0){
        printf("Run %d QAcharge EMPTY -> Return\n",iRun);
        
    }
    
    if((QAcharge)&&(QAchargeRatio)&&(QAcharge->GetEntries()>10)&&(QAchargeRatio->GetEntries()>10)){

        Int_t biny = QAcharge->GetXaxis()->FindBin(747);
        Int_t maxy = QAcharge->GetXaxis()->GetXmax();

        Int_t  contEmpty=0;
        Int_t  contFull=0;

        TH1D *hChargeL5=QAcharge->ProjectionY("hChargeL5",0,biny);
        TH1D *hChargeL6=QAcharge->ProjectionY("hChargeL6",biny,maxy);



        TH1D *hChargeRatioL5=QAchargeRatio->ProjectionY("hChargeRatioL5",0,biny);
        TH1D *hChargeRatioL6=QAchargeRatio->ProjectionY("hChargeRatioL6",biny,maxy);


        if(QAcharge->GetEntries()< 45000)
            contEmpty=1;

        else{
            for(Int_t i =0;i<1698;i++){

                TString tmpQ("Q");
                tmpQ+=i;

                TH1D* fHist1DQ= QAcharge->ProjectionY(tmpQ,i+1,i+1);
                Double_t mean=fHist1DQ->GetMean();

                if(TMath::Abs(mean)<1.0 ||fHist1DQ->GetEntries()<10)
                    contEmpty++;

                else
                    contFull++;

            }
        }
        
        TF1 *lfunLay5 = new TF1("LangausFunLay5",LangausFun,50.,300.,4);
        lfunLay5->SetParameter(0,5.);
        lfunLay5->SetParameter(1,80.);
        lfunLay5->SetParameter(2,hChargeL5->GetEntries()/10.);
        lfunLay5->SetParameter(3,10.);
        lfunLay5->SetParLimits(3,0.,20);
        hChargeL5->Fit(lfunLay5,"NQLR");
        TF1 *lfunLay6 = new TF1("LangausFunLay6",LangausFun,50.,300.,4);
        lfunLay6->SetParameter(0,5.);
        lfunLay6->SetParameter(1,80.);
        lfunLay6->SetParameter(2,hChargeL6->GetEntries()/10.);
        lfunLay6->SetParameter(3,10.);
        lfunLay6->SetParLimits(3,0.,20);
        hChargeL6->Fit(lfunLay6,"NQLR");
        





        MPVL5=(Float_t)lfunLay5->GetParameter(1);
        MPVErrL5=(Float_t)lfunLay5->GetParError(1);
        MPVL6=(Float_t)lfunLay6->GetParameter(1);
        MPVErrL6=(Float_t)lfunLay6->GetParError(1);
        ChargeRatioL5=(Float_t)hChargeRatioL5->GetMean();
        ChargeRatioErrL5=(Float_t)hChargeRatioL5->GetMeanError();
        ChargeRatioL6=(Float_t)hChargeRatioL6->GetMean();
        ChargeRatioErrL6=(Float_t)hChargeRatioL6->GetMeanError();
        EmptyModulesSDD=(Float_t)contEmpty;

    }

    /////////// end of SSD part
    
    
    
    
    
    
    
    
    
    
    
    
    
    ///////////////////////  SDD Part
    char SDDDirName[25]="SDD_Performance";
    
    char SDDListName[15]="coutputRP";
    
    TDirectoryFile * SDDQAdir=(TDirectoryFile*)fin->Get(SDDDirName);
    if (!SDDQAdir) {
        Printf("ERROR: SDD QA directory not present in input file.\n");
        return -1;
    }
    TList * SDDList=(TList*)SDDQAdir->Get(SDDListName);
    
    
    if (!SDDList) Printf("WARNING: SDD QA histograms absent or not accessible\n");
    
    
    
    
    Printf("SDD - QA");
    
    
    TH1F* hcllay=(TH1F*)SDDList->FindObject("hCluInLay");
    
    


    if(hcllay->GetBinContent(1)>0){
        for(Int_t iLay=0; iLay<6; iLay++){
            fracT[iLay]=hcllay->GetBinContent(iLay+2)/hcllay->GetBinContent(1);
            efracT[iLay]=TMath::Sqrt(fracT[iLay]*(1-fracT[iLay])/hcllay->GetBinContent(1));
        }
    }

    fracTrackWithClu1=fracT[0];
    errfracTrackWithClu1=efracT[0];
    fracTrackWithClu2=fracT[1];
    errfracTrackWithClu2=efracT[1];
    fracTrackWithClu3=fracT[2];
    errfracTrackWithClu3=efracT[2];
    fracTrackWithClu4=fracT[3];
    errfracTrackWithClu4=efracT[3];
    fracTrackWithClu5=fracT[4];
    errfracTrackWithClu5=efracT[4];
    fracTrackWithClu6=fracT[5];
    errfracTrackWithClu6=efracT[5];
    
    cout<<endl<<errfracTrackWithClu6<<endl;

    TH1F* hmodT=(TH1F*)SDDList->FindObject("hTPMod");

    
    TH1F* hgamod=(TH1F*)SDDList->FindObject("hGAMod");
    


    
    TH1F* hev=(TH1F*)SDDList->FindObject("hNEvents");

    nTotEvents=hev->GetBinContent(2);
    nTrigEvents=hev->GetBinContent(3);
    nEvents=nTotEvents;


    
    TH1F* htimT=(TH1F*)SDDList->FindObject("hDrTimTPAll");

    
    TH1F* htimTe=(TH1F*)SDDList->FindObject("hDrTimTPExtra");
    

    if(htimT->GetEntries()>0){
        fracExtra=htimTe->GetEntries()/htimT->GetEntries();
        errfracExtra=TMath::Sqrt(htimTe->GetEntries())/htimT->GetEntries();
    }

    for(Int_t iBin=1; iBin<=htimT->GetNbinsX(); iBin++){
        Float_t tim=htimT->GetBinCenter(iBin);
        if(tim>2000. && tim<4000.){
            averPoints+=htimT->GetBinContent(iBin);
            cntBins+=1;
        }
    }

    if(cntBins>0){
        averPoints/=cntBins;
        for(Int_t iBin=1; iBin<=htimT->GetNbinsX(); iBin++){
            if(htimT->GetBinContent(iBin)>0.5*averPoints){
                minDrTime=htimT->GetBinCenter(iBin);
                errminDrTime=0.5*htimT->GetBinWidth(iBin);
                break;
            }
        }
    }
    meanDrTime=htimT->GetMean();
    errmeanDrTime=htimT->GetMeanError();
    TH2F* hdedxmod=(TH2F*)SDDList->FindObject("hdEdxVsMod");
    
    
    TH1D* hdedxLay3=hdedxmod->ProjectionY("hdedxLay3",1,84);
    TH1D* hdedxLay4=hdedxmod->ProjectionY("hdedxLay4",85,260);
    
    TH1F* hSigTim0=(TH1F*)SDDList->FindObject("hSigTimeInt0");

    
    TH1F* hSigTim5=(TH1F*)SDDList->FindObject("hSigTimeInt5");

    //Fitting the same distributions in order to have the MPV
    TF1 *lfunLay3 = new TF1("LangausFunLay3",LangausFun,50.,300.,4);
    lfunLay3->SetParameter(0,5.);
    lfunLay3->SetParameter(1,80.);
    lfunLay3->SetParameter(2,hdedxLay3->GetEntries()/10.);
    lfunLay3->SetParameter(3,10.);
    lfunLay3->SetParLimits(3,0.,20);
    hdedxLay3->Fit(lfunLay3,"NQLR");
    TF1 *lfunLay4 = new TF1("LangausFunLay4",LangausFun,50.,300.,4);
    lfunLay4->SetParameter(0,5.);
    lfunLay4->SetParameter(1,80.);
    lfunLay4->SetParameter(2,hdedxLay4->GetEntries()/10.);
    lfunLay4->SetParameter(3,10.);
    lfunLay4->SetParLimits(3,0.,20);
    hdedxLay4->Fit(lfunLay4,"NQLR");
    TF1 *lfunTim0 = new TF1("LangausFunTim0",LangausFun,50.,300.,4);
    lfunTim0->SetParameter(0,5.);
    lfunTim0->SetParameter(1,80.);
    lfunTim0->SetParameter(2,hSigTim0->GetEntries()/10.);
    lfunTim0->SetParameter(3,10.);
    lfunTim0->SetParLimits(3,0.,20);
    hSigTim0->Fit(lfunTim0,"NQLR");
    TF1 *lfunTim5 = new TF1("LangausFunTim5",LangausFun,50.,300.,4);
    lfunTim5->SetParameter(0,5.);
    lfunTim5->SetParameter(1,80.);
    lfunTim5->SetParameter(2,hSigTim5->GetEntries()/10.);
    lfunTim5->SetParameter(3,10.);
    lfunTim5->SetParLimits(3,0.,20);
    hSigTim5->Fit(lfunTim5,"NQLR");
    
    




    MPVdEdxLay3=lfunLay3->GetParameter(1);
    errMPVdEdxLay3=lfunLay3->GetParError(1);
    MPVdEdxLay4=lfunLay4->GetParameter(1);
    errMPVdEdxLay4=lfunLay4->GetParError(1);
    MPVdEdxTB0=lfunTim0->GetParameter(1);
    errMPVdEdxTB0=lfunTim0->GetParError(1);
    MPVdEdxTB5=lfunTim5->GetParameter(1);
    errMPVdEdxTB5=lfunTim5->GetParError(1);
    
    

    
    /////////// end of SDD part
    
    
    // Matching Part
    
    cout<<"Tracking"<<endl;
    
    TDirectoryFile *dirMatch=(TDirectoryFile*)fin->GetDirectory("ITS_Performance");
    TList *list=NULL;
    TList *listSPD=NULL;
    
    if(dirMatch) {

        list = (TList*)dirMatch->Get("cOutputITS"); // LHC12e
    }
    dirMatch=(TDirectoryFile*)fin->GetDirectory("SPD_Performance");
    if(dirMatch) listSPD = (TList*)dirMatch->Get("coutput1");
    
    // if(!list) return kFALSE;
    
    Float_t ioValues[30];
    Float_t ioErrors[30];
    for(Int_t jj=0;jj<30;jj++){
        ioValues[jj]=0.;
        ioErrors[jj]=0.;
    }
    
    Float_t ptbin=0;
    
    TH1F *hFiredChip = (TH1F*)listSPD->FindObject("hFiredChip");
    
    if(hFiredChip->GetEntries()==0){
        printf("Run %d hFiredChip EMPTY -> Return\n",iRun);
        
    }
    
    Int_t nHSsInner=0,nHSsOuter=0;
    for(Int_t i=0;i<400;i++) if(hFiredChip->GetBinContent(i)>0) nHSsInner++;
    for(Int_t i=400;i<1200;i++) if(hFiredChip->GetBinContent(i)>0) nHSsOuter++;
    nHSsInner = (Int_t)(nHSsInner/10);
    nHSsOuter = (Int_t)(nHSsOuter/10);
    
    ioValues[0]=(Float_t)nHSsInner/40.;
    ioValues[1]=(Float_t)nHSsOuter/80.;
    
    TH1F *fHistPtTPCInAcc = (TH1F*)list->FindObject("fHistPtTPCInAcc");
    Int_t check1=0;
    
    if(fHistPtTPCInAcc->GetEntries()==0){
        check1=1;
        printf("Run %dfHistPtTPCInAcc  EMPTY -> Return\n",iRun);
    }
    
    TH1F *fHistPtITSMI6InAcc = (TH1F*)list->FindObject("fHistPtITSMI6InAcc");

    
    TH1F *fHistPtITSMI5InAcc = (TH1F*)list->FindObject("fHistPtITSMI5InAcc");
    
    TH1F *fHistPtITSMI4InAcc = (TH1F*)list->FindObject("fHistPtITSMI4InAcc");
    
    TH1F *fHistPtITSMI3InAcc = (TH1F*)list->FindObject("fHistPtITSMI3InAcc");

    TH1F *fHistPtITSMI2InAcc = (TH1F*)list->FindObject("fHistPtITSMI2InAcc");
    
    TH1F *fHistPtITSMISPDInAcc = (TH1F*)list->FindObject("fHistPtITSMISPDInAcc");
    
    TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)list->FindObject("fHistPtITSMIoneSPDInAcc");
    
    
    
    TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
    fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
    fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
    fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
    fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);
    
    
    
    fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
    ptbin=fHistPtITSMI6InAcc->FindBin(0.201);
    ioValues[2]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[2]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(1.001);
    ioValues[3]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[3]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI6InAcc->FindBin(10.001);
    ioValues[4]=fHistPtITSMI6InAcc->GetBinContent(ptbin);
    ioErrors[4]=fHistPtITSMI6InAcc->GetBinError(ptbin);
    
    fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
    
    ptbin=fHistPtITSMI5InAcc->FindBin(0.201);
    ioValues[5]=fHistPtITSMI5InAcc->GetBinContent(ptbin);
    ioErrors[5]=fHistPtITSMI5InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI5InAcc->FindBin(1.001);
    ioValues[6]=fHistPtITSMI5InAcc->GetBinContent(ptbin);
    ioErrors[6]=fHistPtITSMI5InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI5InAcc->FindBin(10.001);
    ioValues[7]=fHistPtITSMI5InAcc->GetBinContent(ptbin);
    ioErrors[7]=fHistPtITSMI5InAcc->GetBinError(ptbin);
    
    fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
    
    ptbin=fHistPtITSMI4InAcc->FindBin(0.201);
    ioValues[9]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
    ioErrors[9]=fHistPtITSMI4InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI4InAcc->FindBin(1.001);
    ioValues[10]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
    ioErrors[10]=fHistPtITSMI4InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI4InAcc->FindBin(10.001);
    ioValues[11]=fHistPtITSMI4InAcc->GetBinContent(ptbin);
    ioErrors[11]=fHistPtITSMI4InAcc->GetBinError(ptbin);
    
    fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");
    
    ptbin=fHistPtITSMI3InAcc->FindBin(0.201);
    ioValues[12]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
    ioErrors[12]=fHistPtITSMI3InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI3InAcc->FindBin(1.001);
    ioValues[13]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
    ioErrors[13]=fHistPtITSMI3InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI3InAcc->FindBin(10.001);
    ioValues[14]=fHistPtITSMI3InAcc->GetBinContent(ptbin);
    ioErrors[14]=fHistPtITSMI3InAcc->GetBinError(ptbin);
    
    fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");
    
    ptbin=fHistPtITSMI2InAcc->FindBin(0.201);
    ioValues[15]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
    ioErrors[15]=fHistPtITSMI2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI2InAcc->FindBin(1.001);
    ioValues[16]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
    ioErrors[16]=fHistPtITSMI2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMI2InAcc->FindBin(10.001);
    ioValues[17]=fHistPtITSMI2InAcc->GetBinContent(ptbin);
    ioErrors[17]=fHistPtITSMI2InAcc->GetBinError(ptbin);
    
    fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
    ptbin=fHistPtITSMISPDInAcc->FindBin(0.201);
    ioValues[18]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[18]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(1.001);
    ioValues[19]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[19]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMISPDInAcc->FindBin(10.001);
    ioValues[20]=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
    ioErrors[20]=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    
    fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");
    
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(0.201);
    ioValues[21]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[21]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(1.001);
    ioValues[22]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[22]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIoneSPDInAcc->FindBin(10.001);
    ioValues[23]=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
    ioErrors[23]=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    
    
    fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
    ptbin=fHistPtITSMIge2InAcc->FindBin(0.201);
    ioValues[24]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[24]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(1.001);
    ioValues[25]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[25]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    ptbin=fHistPtITSMIge2InAcc->FindBin(10.001);
    ioValues[26]=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
    ioErrors[26]=fHistPtITSMIge2InAcc->GetBinError(ptbin);
    
    FracSPD1=ioValues[0];
    errFracSPD1=ioErrors[0];
    FracSPD2=ioValues[1];
    errFracSPD2=ioErrors[1];
    Eff6Pt02=ioValues[2];
    errEff6Pt02=ioErrors[2];
    Eff6Pt1=ioValues[3];
    errEff6Pt1=ioErrors[3];
    Eff6Pt10=ioValues[4];
    errEff6Pt10=ioErrors[4];
    Eff5Pt02=ioValues[5];
    errEff5Pt02=ioErrors[5];
    Eff5Pt1=ioValues[6];
    errEff5Pt1=ioErrors[6];
    Eff5Pt10=ioValues[7];
    errEff5Pt10=ioErrors[7];
    Eff4Pt02=ioValues[8];
    errEff4Pt02=ioErrors[8];
    Eff4Pt1=ioValues[9];
    errEff4Pt1=ioErrors[9];
    Eff4Pt10=ioValues[10];
    errEff4Pt10=ioErrors[10];
    Eff3Pt02=ioValues[11];
    errEff3Pt02=ioErrors[11];
    Eff3Pt1=ioValues[12];
    errEff3Pt1=ioErrors[12];
    Eff3Pt10=ioValues[13];
    errEff3Pt10=ioErrors[13];
    Eff2Pt02=ioValues[14];
    errEff2Pt02=ioErrors[14];
    Eff2Pt1=ioValues[15];
    errEff2Pt1=ioErrors[15];
    Eff2Pt10=ioValues[16];
    errEff2Pt10=ioErrors[16];
    EffSPDPt02=ioValues[17];
    errEffSPDPt02=ioErrors[17];
    EffSPDPt1=ioValues[18];
    errEffSPDPt1=ioErrors[18];
    EffSPDPt10=ioValues[19];
    errEffSPDPt10=ioErrors[19];
    EffoneSPDPt02=ioValues[20];
    errEffoneSPDPt02=ioErrors[20];
    EffoneSPDPt1=ioValues[21];
    errEffoneSPDPt1=ioErrors[21];
    EffoneSPDPt10=ioValues[22];
    errEffoneSPDPt10=ioErrors[22];
    EffTOTPt02=ioValues[23];
    errEffTOTPt02=ioErrors[23];
    EffTOTPt1=ioValues[24];
    errEffTOTPt1=ioErrors[24];
    EffTOTPt10=ioValues[25];
    errEffTOTPt10=ioErrors[25];
    
    
    
    
    
    
    
    


    ttree->Fill();
    printf("==============  Saving trending quantities in tree for run %i ===============\n",runNumber);
    trendFile->cd();
    ttree->Write();
    trendFile->Close();


    return  1;
}




Double_t LangausFun(Double_t *x, Double_t *par) {
    
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation),
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.
    
    // Numeric constants
    Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    Double_t mpshift  = -0.22278298;       // Landau maximum location
    
    // Control constants
    Double_t np = 100.0;      // number of convolution steps
    Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
    
    // Variables
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;
    
    
    // MP shift correction
    mpc = par[1] - mpshift * par[0];
    
    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    
    step = (xupp-xlow) / np;
    
    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
        
        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }
    
    return (par[2] * step * sum * invsq2pi / par[3]);
}


