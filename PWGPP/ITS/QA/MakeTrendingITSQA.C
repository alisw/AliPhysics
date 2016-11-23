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

// global variables for TTree
///// SDD Variables (41 variables)
Int_t nrun,nEvents,nEventsTriggered;
Float_t fracTrackWithClu1,fracTrackWithClu2,errfracTrackWithClu1,errfracTrackWithClu2;
Float_t fracTrackWithClu3,fracTrackWithClu4,errfracTrackWithClu3,errfracTrackWithClu4;
Float_t fracTrackWithClu5,fracTrackWithClu6,errfracTrackWithClu5,errfracTrackWithClu6;
Float_t fracExtra,errfracExtra,fracEvWithSDD,errfracEvWithSDD;
Float_t minDrTime,errminDrTime,meanDrTime,errmeanDrTime;
Float_t MPVdEdxLay3,errMPVdEdxLay3,MPVdEdxLay4,errMPVdEdxLay4;
Float_t MPVdEdxTB0,errMPVdEdxTB0,MPVdEdxTB5,errMPVdEdxTB5;
Float_t fracDead3,errfracDead3,fracDead4,errfracDead4;
Float_t FlagSDD1, FlagSDD2; // flag on fraction of SDD modules ON
Float_t FlagMinTime, FlagMeanTime, FlagdEdx3, FlagdEdx4; // flag on SDD time and charge parameters

///// Vertex Variables (26 variables)
Float_t meanVtxTRKx,meanVtxTRKy,meanVtxTRKz;
Float_t meanVtxSPDx,meanVtxSPDy,meanVtxSPDz;
Float_t sigmaVtxTRKx,sigmaVtxTRKy,sigmaVtxTRKz;
Float_t sigmaVtxSPDx,sigmaVtxSPDy,sigmaVtxSPDz;
Float_t meanVtxTRKxErr,meanVtxTRKyErr,meanVtxTRKzErr;
Float_t meanVtxSPDxErr,meanVtxSPDyErr,meanVtxSPDzErr;
Float_t sigmaVtxTRKxErr,sigmaVtxTRKyErr,sigmaVtxTRKzErr;
Float_t sigmaVtxSPDxErr,sigmaVtxSPDyErr,sigmaVtxSPDzErr;
Float_t pileupSPD,errpileupSPD;

///// SSD Variables (25 variables)
Float_t MPVL5,MPVErrL5;
Float_t MPVL6,MPVErrL6;
Float_t ChargeRatioL5,ChargeRatioErrL5;
Float_t ChargeRatioL6,ChargeRatioErrL6;
Float_t EmptyModulesSSD;
Float_t FracBadn5,errFracBadn5,FracBadp5,errFracBadp5,FracBadn6,errFracBadn6,FracBadp6,errFracBadp6;
Float_t FlagSSD1n,FlagSSD1p,FlagSSD2n,FlagSSD2p; // flag on fraction of SSD bad strips
Float_t FlagChR5,FlagChR6,FlagdEdx5,FlagdEdx6; // flag on SSD CR and MPV

///// Matching Variables

///// TPC-ITS (52+24+2 variables)
Float_t FracSPD1,errFracSPD1,FracSPD2,errFracSPD2,Eff6Pt02,errEff6Pt02,Eff6Pt1,errEff6Pt1,Eff6Pt10,errEff6Pt10;
Float_t Eff5Pt02,errEff5Pt02,Eff5Pt1,errEff5Pt1,Eff5Pt10,errEff5Pt10;
Float_t Eff4Pt02,errEff4Pt02,Eff4Pt1,errEff4Pt1,Eff4Pt10,errEff4Pt10;
Float_t Eff3Pt02,errEff3Pt02,Eff3Pt1,errEff3Pt1,Eff3Pt10,errEff3Pt10;
Float_t Eff2Pt02,errEff2Pt02,Eff2Pt1,errEff2Pt1,Eff2Pt10,errEff2Pt10;
Float_t EffSPDPt02,errEffSPDPt02,EffSPDPt1,errEffSPDPt1,EffSPDPt10,errEffSPDPt10;
Float_t EffoneSPDPt02,errEffoneSPDPt02,EffoneSPDPt1,errEffoneSPDPt1,EffoneSPDPt10,errEffoneSPDPt10;
Float_t EffTOTPt02,errEffTOTPt02,EffTOTPt1,errEffTOTPt1,EffTOTPt10,errEffTOTPt10;
Float_t FracTrackMI1,errFracTrackMI1,FracTrackMI2,errFracTrackMI2,FracTrackMI3,errFracTrackMI3;
Float_t FracTrackMI4,errFracTrackMI4,FracTrackMI5,errFracTrackMI5,FracTrackMI6,errFracTrackMI6;
Float_t FracTrackSA1,errFracTrackSA1,FracTrackSA2,errFracTrackSA2,FracTrackSA3,errFracTrackSA3;
Float_t FracTrackSA4,errFracTrackSA4,FracTrackSA5,errFracTrackSA5,FracTrackSA6,errFracTrackSA6;
Float_t FlagSPD1, FlagSPD2;

///// TOF-ITS (48 variables)
Float_t Eff6Pt02TOF,errEff6Pt02TOF,Eff6Pt1TOF,errEff6Pt1TOF,Eff6Pt10TOF,errEff6Pt10TOF;
Float_t Eff5Pt02TOF,errEff5Pt02TOF,Eff5Pt1TOF,errEff5Pt1TOF,Eff5Pt10TOF,errEff5Pt10TOF;
Float_t Eff4Pt02TOF,errEff4Pt02TOF,Eff4Pt1TOF,errEff4Pt1TOF,Eff4Pt10TOF,errEff4Pt10TOF;
Float_t Eff3Pt02TOF,errEff3Pt02TOF,Eff3Pt1TOF,errEff3Pt1TOF,Eff3Pt10TOF,errEff3Pt10TOF;
Float_t Eff2Pt02TOF,errEff2Pt02TOF,Eff2Pt1TOF,errEff2Pt1TOF,Eff2Pt10TOF,errEff2Pt10TOF;
Float_t EffSPDPt02TOF,errEffSPDPt02TOF,EffSPDPt1TOF,errEffSPDPt1TOF,EffSPDPt10TOF,errEffSPDPt10TOF;
Float_t EffoneSPDPt02TOF,errEffoneSPDPt02TOF,EffoneSPDPt1TOF,errEffoneSPDPt1TOF,EffoneSPDPt10TOF,errEffoneSPDPt10TOF;
Float_t EffTOTPt02TOF,errEffTOTPt02TOF,EffTOTPt1TOF,errEffTOTPt1TOF,EffTOTPt10TOF,errEffTOTPt10TOF;

///// ITSsa (12+18+2+252 variables)
Float_t NITSTPCPtBin0,NITSTPCPtBin1,NITSTPCPtBin2,NITSsaPtBin0,NITSsaPtBin1,NITSsaPtBin2,NITSpureSAPtBin0;
Float_t NITSpureSAPtBin1,NITSpureSAPtBin2,ratioPtBin0,ratioPtBin1,ratioPtBin2,NcluITSpSA,errNcluITSpSA;
Float_t dedx4_3,errdedx4_3,PtpionpSA,errPtpionpSA,NclupSA0,errNclupSA0,NclupSA1,errNclupSA1,NclupSA2;
Float_t errNclupSA2,NclupSA3,errNclupSA3,NclupSA4,errNclupSA4,NclupSA5,errNclupSA5;
Float_t chi2TPCITS,chi2ITSpureSA;
Float_t occ_eta_1[2],occ_eta_2[2],occ_eta_3[2],occ_eta_4[2],occ_eta_5[2],occ_eta_6[2];
Float_t occ_phi_1[40],occ_phi_2[40],occ_phi_3[40],occ_phi_4[40],occ_phi_5[40],occ_phi_6[40];

///// pileup SPD (10 variables)
Float_t npilvtx,errnpilvtx,ntrklpil,errntrklpil,ntrklnopil,errntrklnopil,ncl1pil,errncl1pil,ncl1nopil,errncl1nopil;

///// ITS PID for TPCITS tracks (8 variables)
Float_t nsigmapi02,errnsigmapi02,nsigmapi05,errnsigmapi05,nsigmapi1,errnsigmapi1,nsigmapi3,errnsigmapi3;

///// DCAxy and DCAz for TPCITS tracks (32 variables)
Float_t mdca05,errmdca05,rmsdca05,errrmsdca05,mdca1,errmdca1,rmsdca1,errrmsdca1,mdca5,errmdca5,rmsdca5,errrmsdca5,mdca10,errmdca10,rmsdca10,errrmsdca10;
Float_t mdcaz05,errmdcaz05,rmsdcaz05,errrmsdcaz05,mdcaz1,errmdcaz1,rmsdcaz1,errrmsdcaz1,mdcaz5,errmdcaz5,rmsdcaz5,errrmsdcaz5,mdcaz10,errmdcaz10,rmsdcaz10,errrmsdcaz10;

ofstream myfile("logfile.txt",ios::app);

void FillVertexBranches(TList * VertxList);
void FillSPDBranches(TList * SPDList);
void FillSDDBranches(TList * SDDList);
void FillSSDBranches(TList * SSDList);
void FillMatchingBranches(TList * ITSList);
void FillITSsaBranches(TList * ITSsaList);
void FillPileupBranches(TList * PileUPList);
void FillPIDBranches(TList * PIDList);
void FillDCABranches(TList * DCAList);
Int_t MakeTrendingITSQA(TString qafilename,Int_t runNumber=133505,Bool_t isMC=kFALSE,Bool_t canvasE=kFALSE,Bool_t IsOnGrid=kTRUE,TString ocdbStorage= "raw://");

Double_t LangausFun(Double_t *x, Double_t *par) ;


/////////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t MakeTrendingITSQA(TString qafilename,       // full path of the QA output; set IsOnGrid to prepend "alien://"
                             Int_t runNumber,          // run number
                             Bool_t isMC,       // MC flag, to disable meaningless checks
                             Bool_t canvasE,  // enable display plots on canvas and save on root file
                             Bool_t IsOnGrid, // set to kTRUE to access files on the grid
                             TString ocdbStorage) // set the default ocdb storage
{
    // macro to generate tree with ITS QA trending variables
    // access qa PWGPP output files
    if (!qafilename) {
        printf("Error - Invalid input file");
        myfile << "Error - Invalid input file" << endl;
        return 1;
    }

//    char defaultQAoutput[30]="QAresults.root";
    char * treePostFileName="trending.root";

    if (IsOnGrid) TGrid::Connect("alien://");
    TFile * fin = TFile::Open(qafilename,"r");
    if (!fin) {
        printf("ERROR: QA output not found. Exiting...\n");
        myfile <<"ERROR: QA output not found. Exiting... " << endl;
        return -1;
    } else {
        Printf("INFO: QA output file %s open. \n",fin->GetName());
        myfile <<"INFO: QA output file " << fin->GetName() << endl;
    }

    // initializing TTree variables
    
    nrun=-999;nEvents=-999;nEventsTriggered=-999;
    fracTrackWithClu1=-999.;fracTrackWithClu2=-999.;errfracTrackWithClu1=-999.;errfracTrackWithClu2=-999.;
    fracTrackWithClu3=-999.;fracTrackWithClu4=-999.;errfracTrackWithClu3=-999.;errfracTrackWithClu4=-999.;
    fracTrackWithClu5=-999.;fracTrackWithClu6=-999.;errfracTrackWithClu5=-999.;errfracTrackWithClu6=-999.;
    fracExtra=-999.;errfracExtra=-999.;fracEvWithSDD=-999.;errfracEvWithSDD=-999.;
    minDrTime=-999.;errminDrTime=-999.;meanDrTime=-999.;errmeanDrTime=-999.;
    MPVdEdxLay3=-999.;errMPVdEdxLay3=-999.;MPVdEdxLay4=-999.;errMPVdEdxLay4=-999.;
    MPVdEdxTB0=-999.;errMPVdEdxTB0=-999.;MPVdEdxTB5=-999.;errMPVdEdxTB5=-999.;
    //
    fracDead3=-999.;errfracDead3=-999.;fracDead4=-999.;errfracDead4=-999.;
    FlagSDD1=-1.; FlagSDD2=-1.; // flag on fraction of SDD modules ON
    FlagMinTime=-1.; FlagMeanTime=-1.; FlagdEdx3=-1.; FlagdEdx4=-1.; // flag on SDD time and charge parameters
    
    ///// Vertex Variables (26 variables)
    meanVtxTRKx=-999.;meanVtxTRKy=-999.;meanVtxTRKz=-999.;
    meanVtxSPDx=-999.;meanVtxSPDy=-999.;meanVtxSPDz=-999.;
    sigmaVtxTRKx=-999.;sigmaVtxTRKy=-999.;sigmaVtxTRKz=-999.;
    sigmaVtxSPDx=-999.;sigmaVtxSPDy=-999.;sigmaVtxSPDz=-999.;
    meanVtxTRKxErr=-999.;meanVtxTRKyErr=-999.;meanVtxTRKzErr=-999.;
    meanVtxSPDxErr=-999.;meanVtxSPDyErr=-999.;meanVtxSPDzErr=-999.;
    sigmaVtxTRKxErr=-999.;sigmaVtxTRKyErr=-999.;sigmaVtxTRKzErr=-999.;
    sigmaVtxSPDxErr=-999.;sigmaVtxSPDyErr=-999.;sigmaVtxSPDzErr=-999.;
    //
    pileupSPD=-999.;errpileupSPD=-999.;
    
    
    ///// SSD Variables (25 variables)
    MPVL5=-999.;MPVErrL5=-999.;
    MPVL6=-999.;MPVErrL6=-999.;
    ChargeRatioL5=-999.;ChargeRatioErrL5=-999.;
    ChargeRatioL6=-999.;ChargeRatioErrL6=-999.;
    EmptyModulesSSD=-999.;
    //
    FracBadn5=-999.;errFracBadn5=-999.;FracBadp5=-999.;errFracBadp5=-999.;
    FracBadn6=-999.;errFracBadn6=-999.;FracBadp6=-999.;errFracBadp6=-999.;
    FlagSSD1n=-999.;FlagSSD1p=-999.;FlagSSD2n=-999.;FlagSSD2p=-999.; // flag on fraction of SSD bad strips
    FlagChR5=-999.;FlagChR6=-999.;FlagdEdx5=-999.;FlagdEdx6=-999.; // flag on SSD CR and MPV
    
    
    ///// Matching Variables
    
    ///// TPC-ITS (52+24+2 variables)
    FracSPD1=-999.;errFracSPD1=-999.;FracSPD2=-999.;errFracSPD2=-999.;
    Eff6Pt02=-999.;errEff6Pt02=-999.;Eff6Pt1=-999.;errEff6Pt1=-999.;Eff6Pt10=-999.;errEff6Pt10=-999.;
    Eff5Pt02=-999.;errEff5Pt02=-999.;Eff5Pt1=-999.;errEff5Pt1=-999.;Eff5Pt10=-999.;errEff5Pt10=-999.;
    Eff4Pt02=-999.;errEff4Pt02=-999.;Eff4Pt1=-999.;errEff4Pt1=-999.;Eff4Pt10=-999.;errEff4Pt10=-999.;
    Eff3Pt02=-999.;errEff3Pt02=-999.;Eff3Pt1=-999.;errEff3Pt1=-999.;Eff3Pt10=-999.;errEff3Pt10=-999.;
    Eff2Pt02=-999.;errEff2Pt02=-999.;Eff2Pt1=-999.;errEff2Pt1=-999.;Eff2Pt10=-999.;errEff2Pt10=-999.;
    EffSPDPt02=-999.;errEffSPDPt02=-999.;EffSPDPt1=-999.;errEffSPDPt1=-999.;EffSPDPt10=-999.;errEffSPDPt10=-999.;
    EffoneSPDPt02=-999.;errEffoneSPDPt02=-999.;EffoneSPDPt1=-999.;errEffoneSPDPt1=-999.;EffoneSPDPt10=-999.;errEffoneSPDPt10=-999.;
    EffTOTPt02=-999.;errEffTOTPt02=-999.;EffTOTPt1=-999.;errEffTOTPt1=-999.;EffTOTPt10=-999.;errEffTOTPt10=-999.;
    //
    FracTrackMI1=-999.;errFracTrackMI1=-999.;FracTrackMI2=-999.;errFracTrackMI2=-999.;FracTrackMI3=-999.;errFracTrackMI3=-999.;
    FracTrackMI4=-999.;errFracTrackMI4=-999.;FracTrackMI5=-999.;errFracTrackMI5=-999.;FracTrackMI6=-999.;errFracTrackMI6=-999.;
    FracTrackSA1=-999.;errFracTrackSA1=-999.;FracTrackSA2=-999.;errFracTrackSA2=-999.;FracTrackSA3=-999.;errFracTrackSA3=-999.;
    FracTrackSA4=-999.;errFracTrackSA4=-999.;FracTrackSA5=-999.;errFracTrackSA5=-999.;FracTrackSA6=-999.;errFracTrackSA6=-999.;
    //
    FlagSPD1=-999.;FlagSPD2=-999.;
    
    
    ///// TOF-ITS (48 variables)
    Eff6Pt02TOF=-999.;errEff6Pt02TOF=-999.;Eff6Pt1TOF=-999.;errEff6Pt1TOF=-999.;Eff6Pt10TOF=-999.;errEff6Pt10TOF=-999.;
    Eff5Pt02TOF=-999.;errEff5Pt02TOF=-999.;Eff5Pt1TOF=-999.;errEff5Pt1TOF=-999.;Eff5Pt10TOF=-999.;errEff5Pt10TOF=-999.;
    Eff4Pt02TOF=-999.;errEff4Pt02TOF=-999.;Eff4Pt1TOF=-999.;errEff4Pt1TOF=-999.;Eff4Pt10TOF=-999.;errEff4Pt10TOF=-999.;
    Eff3Pt02TOF=-999.;errEff3Pt02TOF=-999.;Eff3Pt1TOF=-999.;errEff3Pt1TOF=-999.;Eff3Pt10TOF=-999.;errEff3Pt10TOF=-999.;
    Eff2Pt02TOF=-999.;errEff2Pt02TOF=-999.;Eff2Pt1TOF=-999.;errEff2Pt1TOF=-999.;Eff2Pt10TOF=-999.;errEff2Pt10TOF=-999.;
    EffSPDPt02TOF=-999.;errEffSPDPt02TOF=-999.;EffSPDPt1TOF=-999.;errEffSPDPt1TOF=-999.;EffSPDPt10TOF=-999.;errEffSPDPt10TOF=-999.;
    EffoneSPDPt02TOF=-999.;errEffoneSPDPt02TOF=-999.;EffoneSPDPt1TOF=-999.;errEffoneSPDPt1TOF=-999.;EffoneSPDPt10TOF=-999.;errEffoneSPDPt10TOF=-999.;
    EffTOTPt02TOF=-999.;errEffTOTPt02TOF=-999.;EffTOTPt1TOF=-999.;errEffTOTPt1TOF=-999.;EffTOTPt10TOF=-999.;errEffTOTPt10TOF=-999.;
    
    
    ///// ITSsa (12+18+2+252 variables)
    NITSTPCPtBin0=-999.;NITSTPCPtBin1=-999.;NITSTPCPtBin2=-999.;NITSsaPtBin0=-999.;NITSsaPtBin1=-999.;NITSsaPtBin2=-999.;NITSpureSAPtBin0=-999.;
    NITSpureSAPtBin1=-999.;NITSpureSAPtBin2=-999.;ratioPtBin0=-999.;ratioPtBin1=-999.;ratioPtBin2=-999.;NcluITSpSA=-999.;errNcluITSpSA=-999.;
    dedx4_3=-999.;errdedx4_3=-999.;PtpionpSA=-999.;errPtpionpSA=-999.;NclupSA0=-999.;errNclupSA0=-999.;NclupSA1=-999.;errNclupSA1=-999.;NclupSA2=-999.;
    errNclupSA2=-999.;NclupSA3=-999.;errNclupSA3=-999.;NclupSA4=-999.;errNclupSA4=-999.;NclupSA5=-999.;errNclupSA5=-999.;
    chi2TPCITS=-999.;chi2ITSpureSA=-999.;
    // layers occupancy - TPCITS tracks
    for (Int_t ieta=0;ieta<2;ieta++){
        occ_eta_1[ieta]=-999.;occ_eta_2[ieta]=-999.;occ_eta_3[ieta]=-999.;occ_eta_4[ieta]=-999.;occ_eta_5[ieta]=-999.;occ_eta_6[ieta]=-999.;
    }
    for (Int_t iphi=0;iphi<40;iphi++){
        occ_phi_1[iphi]=-999.;occ_phi_2[iphi]=-999.;occ_phi_3[iphi]=-999.;occ_phi_4[iphi]=-999.;occ_phi_5[iphi]=-999.;occ_phi_6[iphi]=-999.;
    }
    
    ///// pileup SPD (10 variables)
    npilvtx=-999.;errnpilvtx=-999.;ntrklpil=-999.;errntrklpil=-999.;ntrklnopil=-999.;errntrklnopil=-999.;
    ncl1pil=-999.;errncl1pil=-999.;ncl1nopil=-999.;errncl1nopil=-999.;
    
    
    ///// ITS PID for TPCITS tracks (8 variables)
    nsigmapi02=-999.;errnsigmapi02=-999.;nsigmapi05=-999.;errnsigmapi05=-999.;nsigmapi1=-999.;errnsigmapi1=-999.;nsigmapi3=-999.;errnsigmapi3=-999.;
    
    
    ///// DCAxy and DCAz for TPCITS tracks (32 variables)
    mdca05=-999.;errmdca05=-999.;rmsdca05=-999.;errrmsdca05=-999.;mdca1=-999.;errmdca1=-999.;rmsdca1=-999.;errrmsdca1=-999.;
    mdca5=-999.;errmdca5=-999.;rmsdca5=-999.;errrmsdca5=-999.;mdca10=-999.;errmdca10=-999.;rmsdca10=-999.;errrmsdca10=-999.;
    mdcaz05=-999.;errmdcaz05=-999.;rmsdcaz05=-999.;errrmsdcaz05=-999.;mdcaz1=-999.;errmdcaz1=-999.;rmsdcaz1=-999.;errrmsdcaz1=-999.;
    mdcaz5=-999.;errmdcaz5=-999.;rmsdcaz5=-999.;errrmsdcaz5=-999.;mdcaz10=-999.;errmdcaz10=-999.;rmsdcaz10=-999.;errrmsdcaz10=-999.;

    
    TFile * trendFile = new TFile(treePostFileName,"recreate");
    
    ///// TTree creation
    TTree * ttree=new TTree("trending","tree of trending variables");
    
    // SDD branches
    ttree->Branch("nrun",&nrun,"nrun/I");
    ttree->Branch("nEvents",&nEvents,"nEvents/I");
    ttree->Branch("nEventsTriggered",&nEventsTriggered,"nEventsTriggered/I");
    ttree->Branch("minDrTime",&minDrTime,"minDrTime/F"); // minimum SDD drift time (50% of plateau)
    ttree->Branch("errminDrTime",&errminDrTime,"errminDrTime/F");
    ttree->Branch("meanDrTime",&meanDrTime,"meanDrTime/F"); // mean SDD drift time
    ttree->Branch("errmeanDrTime",&errmeanDrTime,"errmeanDrTime/F");
    ttree->Branch("MPVdEdxLay3",&MPVdEdxLay3,"MPVdEdxLay3/F"); // most probable value of dE/dx distribution of SDD Layer 3
    ttree->Branch("errMPVdEdxLay3",&errMPVdEdxLay3,"errMPVdEdxLay3/F"); //
    ttree->Branch("MPVdEdxLay4",&MPVdEdxLay4,"MPVdEdxLay4/F"); // most probable value of dE/dx distribution of SDD Layer 4
    ttree->Branch("errMPVdEdxLay4",&errMPVdEdxLay4,"errMPVdEdxLay4/F"); //
    ttree->Branch("MPVdEdxTB0",&MPVdEdxTB0,"MPVdEdxTB0/F"); // most probable value of dE/dx distribution of SDD - small drift time
    ttree->Branch("errMPVdEdxTB0",&errMPVdEdxTB0,"errMPVdEdxTB0/F"); //
    ttree->Branch("MPVdEdxTB5",&MPVdEdxTB5,"MPVdEdxTB5/F"); // most probable value of dE/dx distribution of SDD - large drift time
    ttree->Branch("errMPVdEdxTB5",&errMPVdEdxTB5,"errMPVdEdxTB5/F"); //
    ttree->Branch("fracDead3",&fracDead3,"fracDead3/F"); // fraction of bad SDD modules layer 3
    ttree->Branch("errfracDead3",&errfracDead3,"errfracDead3/F"); // fraction of bad SDD modules layer 3
    ttree->Branch("fracDead4",&fracDead4,"fracDead4/F"); // fraction of bad SDD modules layer 4
    ttree->Branch("errfracDead4",&errfracDead4,"errfracDead4/F"); // fraction of bad SDD modules layer 4
    ttree->Branch("fracTrackWithClu1",&fracTrackWithClu1,"fracTrackWithClu1/F"); //fraction of tracks with cluster in layer 1 (SDD in DAQ)
    ttree->Branch("errfracTrackWithClu1",&errfracTrackWithClu1,"errfracTrackWithClu1/F"); //
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
    ttree->Branch("fracExtra",&fracExtra,"fracExtra/F"); // fraction of extra clusters in SDD
    ttree->Branch("errfracExtra",&errfracExtra,"errfracExtra/F"); // fraction of extra clusters in SDD
    ttree->Branch("FlagSDD1",&FlagSDD1,"FlagSDD1/F"); // Flag for SDD1 ON fraction
    ttree->Branch("FlagSDD2",&FlagSDD2,"FlagSDD2/F"); // Flag for SDD2 ON fraction
    ttree->Branch("FlagMinTime",&FlagMinTime,"FlagMinTime/F"); // Flag for minimum SDD drift time
    ttree->Branch("FlagMeanTime",&FlagMeanTime,"FlagMeanTime/F"); // Flag for mean SDD drift time
    ttree->Branch("FlagdEdx3",&FlagdEdx3,"FlagdEdx3/F"); // Flag for SDD1 dEdx MPV
    ttree->Branch("FlagdEdx4",&FlagdEdx4,"FlagdEdx4/F"); // Flag for SDD2 dEdx MPV
    
    // SSD branches
    ttree->Branch("ChargeRatioL5",&ChargeRatioL5,"ChargeRatioL5/F"); // Charge ratio (2 sides of SSD) Layer 5
    ttree->Branch("ChargeRatioErrL5",&ChargeRatioErrL5,"ChargeRatioErrL5/F"); // Charge ratio error (2 sides of SSD) Layer 5
    ttree->Branch("ChargeRatioL6",&ChargeRatioL6,"ChargeRatioL6/F"); // Charge ratio (2 sides of SSD) Layer 6
    ttree->Branch("ChargeRatioErrL6",&ChargeRatioErrL6,"ChargeRatioErrL6/F"); // Charge ratio error(2 sides of SSD) Layer 6
    ttree->Branch("MPVL5",&MPVL5,"MPVL5/F"); // Most Probable Value dEdx Layer 5
    ttree->Branch("MPVErrL5",&MPVErrL5,"MPVErrL5/F"); // Most Probable Value error dEdx Layer 5
    ttree->Branch("MPVL6",&MPVL6,"MPVL6/F"); // Most Probable Value dEdx Layer 6
    ttree->Branch("MPVErrL6",&MPVErrL6,"MPVErrL6/F"); // Most Probable Value error dEdx Layer 6
    ttree->Branch("FracBadn5",&FracBadn5,"FracBadn5/F"); // fraction of bad n-strips layer 5
    ttree->Branch("errFracBadn5",&errFracBadn5,"errFracBadn5/F"); // fraction of bad n-strips layer 5
    ttree->Branch("FracBadp5",&FracBadp5,"FracBadp5/F"); // fraction of bad p-strips layer 5
    ttree->Branch("errFracBadp5",&errFracBadp5,"errFracBadp5/F"); // fraction of bad p-strips layer 5
    ttree->Branch("FracBadn6",&FracBadn6,"FracBadn6/F"); // fraction of bad n-strips layer 6
    ttree->Branch("errFracBadn6",&errFracBadn6,"errFracBadn6/F"); // fraction of bad n-strips layer 6
    ttree->Branch("FracBadp6",&FracBadp6,"FracBadp6/F"); // fraction of bad p-strips layer 6
    ttree->Branch("errFracBadp6",&errFracBadp6,"errFracBadp6/F"); // fraction of bad p-strips layer 6
    ttree->Branch("EmptyModulesSSD",&EmptyModulesSSD,"EmptyModulesSSD/F"); // Number of empty SSD  modules
    ttree->Branch("FlagSSD1n",&FlagSSD1n,"FlagSSD1n/F"); // flag for run-to-run variation of SSD1 bad n strips
    ttree->Branch("FlagSSD1p",&FlagSSD1p,"FlagSSD1p/F"); // flag for run-to-run variation of SSD1 bad p strips
    ttree->Branch("FlagSSD2n",&FlagSSD2n,"FlagSSD2n/F"); // flag for run-to-run variation of SSD2 bad n strips
    ttree->Branch("FlagSSD2p",&FlagSSD2p,"FlagSSD2p/F"); // flag for run-to-run variation of SSD2 bad p strips
    ttree->Branch("FlagChR5",&FlagChR5,"FlagChR5/F"); // flag for SSD1 charge ratio
    ttree->Branch("FlagChR6",&FlagChR6,"FlagChR6/F"); // flag for SSD2 charge ratio
    ttree->Branch("FlagdEdx5",&FlagdEdx5,"FlagdEdx5/F"); // flag for SSD1 dEdx MPV
    ttree->Branch("FlagdEdx6",&FlagdEdx6,"FlagdEdx6/F"); // flag for SSD2 dEdx MPV

    // vertex branches
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
    ttree->Branch("pileupSPD",&pileupSPD,"pileupSPD/F"); // fraction of events with SPD pileup vertex
    ttree->Branch("errpileupSPD",&errpileupSPD,"errpileupSPD/F"); // fraction of events with SPD pileup vertex

    // TPC-ITS ME branches
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
    ttree->Branch("EffoneSPDPt02",&EffoneSPDPt02,"EffoneSPDPt02/F"); // matching efficiency low pt one SPD
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
    ttree->Branch("FracTrackMI1",&FracTrackMI1,"FracTrackMI1/F"); // fraction of global tracks with hit in ITS layer 1
    ttree->Branch("errFracTrackMI1",&errFracTrackMI1,"errFracTrackMI1/F");
    ttree->Branch("FracTrackMI2",&FracTrackMI2,"FracTrackMI2/F"); // fraction of global tracks with hit in ITS layer 2
    ttree->Branch("errFracTrackMI2",&errFracTrackMI2,"errFracTrackMI2/F");
    ttree->Branch("FracTrackMI3",&FracTrackMI3,"FracTrackMI3/F"); // fraction of global tracks with hit in ITS layer 3
    ttree->Branch("errFracTrackMI3",&errFracTrackMI3,"errFracTrackMI3/F");
    ttree->Branch("FracTrackMI4",&FracTrackMI4,"FracTrackMI4/F"); // fraction of global tracks with hit in ITS layer 4
    ttree->Branch("errFracTrackMI4",&errFracTrackMI4,"errFracTrackMI4/F");
    ttree->Branch("FracTrackMI5",&FracTrackMI5,"FracTrackMI5/F"); // fraction of global tracks with hit in ITS layer 5
    ttree->Branch("errFracTrackMI5",&errFracTrackMI5,"errFracTrackMI5/F");
    ttree->Branch("FracTrackMI6",&FracTrackMI6,"FracTrackMI6/F"); // fraction of global tracks with hit in ITS layer 6
    ttree->Branch("errFracTrackMI6",&errFracTrackMI6,"errFracTrackMI6/F");
    ttree->Branch("FracTrackSA1",&FracTrackSA1,"FracTrackSA1/F"); // fraction of SA tracks with hit in ITS layer 1
    ttree->Branch("errFracTrackSA1",&errFracTrackSA1,"errFracTrackSA1/F");
    ttree->Branch("FracTrackSA2",&FracTrackSA2,"FracTrackSA2/F"); // fraction of SA tracks with hit in ITS layer 2
    ttree->Branch("errFracTrackSA2",&errFracTrackSA2,"errFracTrackSA2/F");
    ttree->Branch("FracTrackSA3",&FracTrackSA3,"FracTrackSA3/F"); // fraction of SA tracks with hit in ITS layer 3
    ttree->Branch("errFracTrackSA3",&errFracTrackSA3,"errFracTrackSA3/F");
    ttree->Branch("FracTrackSA4",&FracTrackSA4,"FracTrackSA4/F"); // fraction of SA tracks with hit in ITS layer 4
    ttree->Branch("errFracTrackSA4",&errFracTrackSA4,"errFracTrackSA4/F");
    ttree->Branch("FracTrackSA5",&FracTrackSA5,"FracTrackSA5/F"); // fraction of SA tracks with hit in ITS layer 5
    ttree->Branch("errFracTrackSA5",&errFracTrackSA5,"errFracTrackSA5/F");
    ttree->Branch("FracTrackSA6",&FracTrackSA6,"FracTrackSA6/F"); // fraction of SA tracks with hit in ITS layer 6
    ttree->Branch("errFracTrackSA6",&errFracTrackSA6,"errFracTrackSA6/F");

    ttree->Branch("FlagSPD1",&FlagSPD1,"FlagSPD1/F"); // flag for run-to-run variation of SPD1 HS ON fraction
    ttree->Branch("FlagSPD2",&FlagSPD2,"FlagSPD2/F"); // flag for run-to-run variation of SPD2 HS ON fraction

    // TOF-ITS ME branches
    ttree->Branch("Eff6Pt02TOF",&Eff6Pt02TOF,"Eff6Pt02TOF/F"); // matching efficiency low pt 6 clusters
    ttree->Branch("errEff6Pt02TOF",&errEff6Pt02TOF,"errEff6Pt02TOF/F"); // error matching efficiency low pt 6 clusters
    ttree->Branch("Eff5Pt02TOF",&Eff5Pt02TOF,"Eff5Pt02TOF/F"); // matching efficiency low pt 5 clusters
    ttree->Branch("errEff5Pt02TOF",&errEff5Pt02TOF,"errEff5Pt02TOF/F"); // error matching efficiency low pt 5 clusters
    ttree->Branch("Eff4Pt02TOF",&Eff4Pt02TOF,"Eff4Pt02TOF/F"); // matching efficiency low pt 4 clusters
    ttree->Branch("errEff4Pt02TOF",&errEff4Pt02TOF,"errEff4Pt02TOF/F"); // error matching efficiency low pt 4 clusters
    ttree->Branch("Eff3Pt02TOF",&Eff3Pt02TOF,"Eff3Pt02TOF/F"); // matching efficiency low pt 3 clusters
    ttree->Branch("errEff3Pt02TOF",&errEff3Pt02TOF,"errEff3Pt02TOF/F"); // error matching efficiency low pt 3 clusters
    ttree->Branch("Eff2Pt02TOF",&Eff2Pt02TOF,"Eff2Pt02TOF/F"); // matching efficiency low pt 2 clusters
    ttree->Branch("errEff2Pt02TOF",&errEff2Pt02TOF,"errEff2Pt02TOF/F"); // error matching efficiency low pt 2 clusters
    ttree->Branch("EffSPDPt02TOF",&EffSPDPt02TOF,"EffSPDPt02TOF/F"); // matching efficiency low pt 2 SPD
    ttree->Branch("errEffSPDPt02TOF",&errEffSPDPt02TOF,"errEffSPDPt02TOF/F"); // error matching efficiency low pt 2 SPD
    ttree->Branch("EffoneSPDPt02TOF",&EffoneSPDPt02TOF,"EffoneSPDPt02TOF/F"); // matching efficiency low pt 6 one SPD
    ttree->Branch("errEffoneSPDPt02TOF",&errEffoneSPDPt02TOF,"errEffoneSPDPt02TOF/F"); // error matching efficiency low pt one SPD
    ttree->Branch("EffTOTPt02TOF",&EffTOTPt02TOF,"EffTOTPt02TOF/F"); // matching efficiency low pt
    ttree->Branch("errEffTOTPt02TOF",&errEffTOTPt02TOF,"errEffTOTPt02TOF/F"); // error matching efficiency low pt
    
    ttree->Branch("Eff6Pt1TOF",&Eff6Pt1TOF,"Eff6Pt1TOF/F"); // matching efficiency mid pt 6 clusters
    ttree->Branch("errEff6Pt1TOF",&errEff6Pt1TOF,"errEff6Pt1TOF/F"); // error matching efficiency mid pt 6 clusters
    ttree->Branch("Eff5Pt1TOF",&Eff5Pt1TOF,"Eff5Pt1TOF/F"); // matching efficiency mid pt 5 clusters
    ttree->Branch("errEff5Pt1TOF",&errEff5Pt1TOF,"errEff5Pt1TOF/F"); // error matching efficiency mid pt 5 clusters
    ttree->Branch("Eff4Pt1TOF",&Eff4Pt1TOF,"Eff4Pt1TOF/F"); // matching efficiency mid pt 4 clusters
    ttree->Branch("errEff4Pt1TOF",&errEff4Pt1TOF,"errEff4Pt1TOF/F"); // error matching efficiency mid pt 4 clusters
    ttree->Branch("Eff3Pt1TOF",&Eff3Pt1TOF,"Eff3Pt1TOF/F"); // matching efficiency mid pt 3 clusters
    ttree->Branch("errEff3Pt1TOF",&errEff3Pt1TOF,"errEff3Pt1TOF/F"); // error matching efficiency mid pt 3 clusters
    ttree->Branch("Eff2Pt1TOF",&Eff2Pt1TOF,"Eff2Pt1TOF/F"); // matching efficiency mid pt 2 clusters
    ttree->Branch("errEff2Pt1TOF",&errEff2Pt1TOF,"errEff2Pt1TOF/F"); // error matching efficiency mid pt 2 clusters
    ttree->Branch("EffSPDPt1TOF",&EffSPDPt1TOF,"EffSPDPt1TOF/F"); // matching efficiency mid pt 2 SPD
    ttree->Branch("errEffSPDPt1TOF",&errEffSPDPt1TOF,"errEffSPDPt1TOF/F"); // error matching efficiency mid pt 2 SPD
    ttree->Branch("EffoneSPDPt1TOF",&EffoneSPDPt1TOF,"EffoneSPDPt1TOF/F"); // matching efficiency mid pt 6 one SPD
    ttree->Branch("errEffoneSPDPt1TOF",&errEffoneSPDPt1TOF,"errEffoneSPDPt1TOF/F"); // error matching efficiency mid pt one SPD
    ttree->Branch("EffTOTPt1TOF",&EffTOTPt1TOF,"EffTOTPt1TOF/F"); // matching efficiency mid pt
    ttree->Branch("errEffTOTPt1TOF",&errEffTOTPt1TOF,"errEffTOTPt1TOF/F"); // error matching efficiency mid pt
    
    ttree->Branch("Eff6Pt10TOF",&Eff6Pt10TOF,"Eff6Pt10TOF/F"); // matching efficiency high pt 6 clusters
    ttree->Branch("errEff6Pt10TOF",&errEff6Pt10TOF,"errEff6Pt10TOF/F"); // error matching efficiency high pt 6 clusters
    ttree->Branch("Eff5Pt10TOF",&Eff5Pt10TOF,"Eff5Pt10TOF/F"); // matching efficiency high pt 5 clusters
    ttree->Branch("errEff5Pt10TOF",&errEff5Pt10TOF,"errEff5Pt10TOF/F"); // error matching efficiency high pt 5 clusters
    ttree->Branch("Eff4Pt10TOF",&Eff4Pt10TOF,"Eff4Pt10TOF/F"); // matching efficiency high pt 4 clusters
    ttree->Branch("errEff4Pt10TOF",&errEff4Pt10TOF,"errEff4Pt10TOF/F"); // error matching efficiency high pt 4 clusters
    ttree->Branch("Eff3Pt10TOF",&Eff3Pt10TOF,"Eff3Pt10TOF/F"); // matching efficiency high pt 3 clusters
    ttree->Branch("errEff3Pt10TOF",&errEff3Pt10TOF,"errEff3Pt10TOF/F"); // error matching efficiency high pt 3 clusters
    ttree->Branch("Eff2Pt10TOF",&Eff2Pt10TOF,"Eff2Pt10TOF/F"); // matching efficiency high pt 2 clusters
    ttree->Branch("errEff2Pt10TOF",&errEff2Pt10TOF,"errEff2Pt10TOF/F"); // error matching efficiency high pt 2 clusters
    ttree->Branch("EffSPDPt10TOF",&EffSPDPt10TOF,"EffSPDPt10TOF/F"); // matching efficiency high pt 2 SPD
    ttree->Branch("errEffSPDPt10TOF",&errEffSPDPt10TOF,"errEffSPDPt10TOF/F"); // error matching efficiency high pt 2 SPD
    ttree->Branch("EffoneSPDPt10TOF",&EffoneSPDPt10TOF,"EffoneSPDPt10TOF/F"); // matching efficiency high pt 6 one SPD
    ttree->Branch("errEffoneSPDPt10TOF",&errEffoneSPDPt10TOF,"errEffoneSPDPt10TOF/F"); // error matching efficiency high pt one SPD
    ttree->Branch("EffTOTPt10TOF",&EffTOTPt10TOF,"EffTOTPt10TOF/F"); // matching efficiency high pt
    ttree->Branch("errEffTOTPt10TOF",&errEffTOTPt10TOF,"errEffTOTPt10TOF/F"); // error matching efficiency high pt
    //
    
    // ITSsa branches
    ttree->Branch("NITSTPCPtBin0",&NITSTPCPtBin0,"NITSTPCPtBin0/F"); // TPCITS tracks/event low pt
    ttree->Branch("NITSTPCPtBin1",&NITSTPCPtBin1,"NITSTPCPtBin1/F"); // TPCITS tracks/event mid pt
    ttree->Branch("NITSTPCPtBin2",&NITSTPCPtBin2,"NITSTPCPtBin2/F"); // TPCITS tracks/event high pt
    ttree->Branch("NITSsaPtBin0",&NITSsaPtBin0,"NITSsaPtBin0/F"); // ITSsa tracks/event low pt
    ttree->Branch("NITSsaPtBin1",&NITSsaPtBin1,"NITSsaPtBin1/F"); // ITSsa tracks/event mid pt
    ttree->Branch("NITSsaPtBin2",&NITSsaPtBin2,"NITSsaPtBin2/F"); // ITSsa tracks/event high pt
    ttree->Branch("NITSpureSAPtBin0",&NITSpureSAPtBin0,"NITSpureSAPtBin0/F"); // ITSpureSA tracks/event low pt
    ttree->Branch("NITSpureSAPtBin1",&NITSpureSAPtBin1,"NITSpureSAPtBin1/F"); // ITSpureSA tracks/event mid pt
    ttree->Branch("NITSpureSAPtBin2",&NITSpureSAPtBin2,"NITSpureSAPtBin2/F"); // ITSpureSA tracks/event high pt
    ttree->Branch("ratioPtBin0",&ratioPtBin0,"ratioPtBin0/F"); // (TPCITS+ITSsa)/ITSpureSA ratio low pt
    ttree->Branch("ratioPtBin1",&ratioPtBin1,"ratioPtBin1/F"); // (TPCITS+ITSsa)/ITSpureSA ratio mid pt
    ttree->Branch("ratioPtBin2",&ratioPtBin2,"ratioPtBin2/F"); // (TPCITS+ITSsa)/ITSpureSA ratio high pt
    ttree->Branch("NcluITSpSA",&NcluITSpSA,"NcluITSpSA/F"); // mean cluster number ITSpureSA, N>=3
    ttree->Branch("errNcluITSpSA",&errNcluITSpSA,"errNcluITSpSA/F"); // mean cluster number ITSpureSA, N>=3
    ttree->Branch("dedx4_3",&dedx4_3,"dedx4_3/F"); // tracks 4cls / tracks 3 cls SDD+SSD
    ttree->Branch("errdedx4_3",&errdedx4_3,"errdedx4_3/F"); // tracks 4cls / tracks 3 cls SDD+SSD
    ttree->Branch("PtpionpSA",&PtpionpSA,"PtpionpSA/F"); // mean pt ITSpureSA tracks
    ttree->Branch("errPtpionpSA",&errPtpionpSA,"errPtpionpSA/F"); // mean pt ITSpureSA tracks
    ttree->Branch("NclupSA0",&NclupSA0,"NclupSA0/F"); // frac. ITSpureSA tracks w/cls SPD1
    ttree->Branch("errNclupSA0",&errNclupSA0,"errNclupSA0/F"); // frac. ITSpureSA tracks w/cls SPD1
    ttree->Branch("NclupSA1",&NclupSA1,"NclupSA1/F"); // frac. ITSpureSA tracks w/cls SPD2
    ttree->Branch("errNclupSA1",&errNclupSA1,"errNclupSA1/F"); // frac. ITSpureSA tracks w/cls SPD2
    ttree->Branch("NclupSA2",&NclupSA2,"NclupSA2/F"); // frac. ITSpureSA tracks w/cls SDD1
    ttree->Branch("errNclupSA2",&errNclupSA2,"errNclupSA2/F"); // frac. ITSpureSA tracks w/cls SDD1
    ttree->Branch("NclupSA3",&NclupSA3,"NclupSA3/F"); // frac. ITSpureSA tracks w/cls SDD2
    ttree->Branch("errNclupSA3",&errNclupSA3,"errNclupSA3/F"); // frac. ITSpureSA tracks w/cls SDD2
    ttree->Branch("NclupSA4",&NclupSA4,"NclupSA4/F"); // frac. ITSpureSA tracks w/cls SSD1
    ttree->Branch("errNclupSA4",&errNclupSA4,"errNclupSA4/F"); // frac. ITSpureSA tracks w/cls SSD1
    ttree->Branch("NclupSA5",&NclupSA5,"NclupSA5/F"); // frac. ITSpureSA tracks w/cls SSD2
    ttree->Branch("errNclupSA5",&errNclupSA5,"errNclupSA5/F"); // frac. ITSpureSA tracks w/cls SSD2
    ttree->Branch("chi1TPCITS",&chi2TPCITS,"chi2TPCITS/F"); // TPCITS tracks mean chi2
    ttree->Branch("chi1ITSpureSA",&chi2ITSpureSA,"chi2ITSpureSA/F"); // TPCITS tracks mean chi2
//
    ttree->Branch("occ_eta_1",occ_eta_1,"occ_eta_1[2]/F"); // eta occupancy layer 1: 2 intervals  to be checked
    ttree->Branch("occ_eta_2",occ_eta_2,"occ_eta_2[2]/F"); // eta occupancy layer 2: 2 intervals  to be checked
    ttree->Branch("occ_eta_3",occ_eta_3,"occ_eta_3[2]/F"); // eta occupancy layer 3: 2 intervals  to be checked
    ttree->Branch("occ_eta_4",occ_eta_4,"occ_eta_4[2]/F"); // eta occupancy layer 4: 2 intervals  to be checked
    ttree->Branch("occ_eta_5",occ_eta_5,"occ_eta_5[2]/F"); // eta occupancy layer 5: 2 intervals  to be checked
    ttree->Branch("occ_eta_6",occ_eta_6,"occ_eta_6[2]/F"); // eta occupancy layer 6: 2 intervals  to be checked
   
    ttree->Branch("occ_phi_1",occ_phi_1,"occ_phi_1[40]/F"); // phi occupancy layer 1: 40 intervals  to be checked
    ttree->Branch("occ_phi_2",occ_phi_2,"occ_phi_2[40]/F"); // phi occupancy layer 2: 40 intervals  to be checked
    ttree->Branch("occ_phi_3",occ_phi_3,"occ_phi_3[40]/F"); // phi occupancy layer 3: 40 intervals  to be checked
    ttree->Branch("occ_phi_4",occ_phi_4,"occ_phi_4[40]/F"); // phi occupancy layer 4: 40 intervals  to be checked
    ttree->Branch("occ_phi_5",occ_phi_5,"occ_phi_5[40]/F"); // phi occupancy layer 5: 40 intervals  to be checked
    ttree->Branch("occ_phi_6",occ_phi_6,"occ_phi_6[40]/F"); // phi occupancy layer 6: 40 intervals  to be checked

    
    // pileup SPD branches
    ttree->Branch("npilvtx",&npilvtx,"npilvtx/F"); // number of pileup vertices/event
    ttree->Branch("errnpilvtx",&errnpilvtx,"errnpilvtx/F"); // number of pileup vertices/event
    ttree->Branch("ntrklpil",&ntrklpil,"ntrklpil/F"); // number of tracklets/pileup vertex
    ttree->Branch("errntrklpil",&errntrklpil,"errntrklpil/F"); // number of tracklets/pileup vertex
    ttree->Branch("ntrklnopil",&ntrklnopil,"ntrklnopil/F"); // number of tracklets/nopileup vertex
    ttree->Branch("errntrklnopil",&errntrklnopil,"errntrklnopil/F"); // number of tracklets/nopileup vertex
    ttree->Branch("ncl1pil",&ncl1pil,"ncl1pil/F"); // number of SPD2 clusters/pileup vertex
    ttree->Branch("errncl1pil",&errncl1pil,"errncl1pil/F"); // number of SPD2 clusters/pileup vertex
    ttree->Branch("ncl1nopil",&ncl1nopil,"ncl1nopil/F"); // number of SPD2 clusters/no pileup vertex
    ttree->Branch("errncl1nopil",&errncl1nopil,"errncl1nopil/F"); // number of SPD2 clusters/no pileup vertex
    
 
    // ITS PID branches (TPCITS tracks)
    ttree->Branch("nsigmapi02",&nsigmapi02,"nsigmapi02/F"); // nsigma for pions - 0.2 GeV/c
    ttree->Branch("errnsigmapi02",&errnsigmapi02,"errnsigmapi02/F"); // nsigma for pions - 0.2 GeV/c
    ttree->Branch("nsigmapi05",&nsigmapi05,"nsigmapi05/F"); // nsigma for pions - 0.5 GeV/c
    ttree->Branch("errnsigmapi05",&errnsigmapi05,"errnsigmapi05/F"); // nsigma for pions - 0.5 GeV/c
    ttree->Branch("nsigmapi1",&nsigmapi1,"nsigmapi1/F"); // nsigma for pions - 1.0 GeV/c
    ttree->Branch("errnsigmapi1",&errnsigmapi1,"errnsigmapi1/F"); // nsigma for pions - 1.0 GeV/c
    ttree->Branch("nsigmapi3",&nsigmapi3,"nsigmapi3/F"); // nsigma for pions - 10 GeV/c
    ttree->Branch("errnsigmapi3",&errnsigmapi3,"errnsigmapi3/F"); // nsigma for pions - 10 GeV/c

 
    // DCAxy and DCAz branches (TPCITS tracks)
    ttree->Branch("mdca05",&mdca05,"mdca05/F"); // mean DCAxy - 0.5 GeV/c
    ttree->Branch("errmdca05",&errmdca05,"errmdca05/F"); // mean DCAxy - 0.5 GeV/c
    ttree->Branch("mdca1",&mdca1,"mdca1/F"); // mean DCAxy - 1.0 GeV/c
    ttree->Branch("errmdca1",&errmdca1,"errmdca1/F"); // mean DCAxy - 1.0 GeV/c
    ttree->Branch("mdca5",&mdca5,"mdca5/F"); // mean DCAxy - 5.0 GeV/c
    ttree->Branch("errmdca5",&errmdca5,"errmdca5/F"); // mean DCAxy - 5.0 GeV/c
    ttree->Branch("mdca10",&mdca10,"mdca10/F"); // mean DCAxy - 10 GeV/c
    ttree->Branch("errmdca10",&errmdca10,"errmdca10/F"); // mean DCAxy - 10 GeV/c
    ttree->Branch("rmsdca05",&rmsdca05,"rmsdca05/F"); // rms DCAxy - 0.5 GeV/c
    ttree->Branch("errrmsdca05",&errrmsdca05,"errrmsdca05/F"); // rms DCAxy - 0.5 GeV/c
    ttree->Branch("rmsdca1",&rmsdca1,"rmsdca1/F"); // rms DCAxy - 1.0 GeV/c
    ttree->Branch("errrmsdca1",&errrmsdca1,"errrmsdca1/F"); // rms DCAxy - 1.0 GeV/c
    ttree->Branch("rmsdca5",&rmsdca5,"rmsdca5/F"); // rms DCAxy - 5.0 GeV/c
    ttree->Branch("errrmsdca5",&errrmsdca5,"errrmsdca5/F"); // rms DCAxy - 5.0 GeV/c
    ttree->Branch("rmsdca10",&rmsdca10,"rmsdca10/F"); // rms DCAxy - 10 GeV/c
    ttree->Branch("errrmsdca10",&errrmsdca10,"errrmsdca10/F"); // rms DCAxy - 10 GeV/c
   
    ttree->Branch("mdcaz05",&mdcaz05,"mdcaz05/F"); // mean DCAz - 0.5 GeV/c
    ttree->Branch("errmdcaz05",&errmdcaz05,"errmdcaz05/F"); // mean DCAz - 0.5 GeV/c
    ttree->Branch("mdcaz1",&mdcaz1,"mdcaz1/F"); // mean DCAz - 1.0 GeV/c
    ttree->Branch("errmdcaz1",&errmdcaz1,"errmdcaz1/F"); // mean DCAz - 1.0 GeV/c
    ttree->Branch("mdcaz5",&mdcaz5,"mdcaz5/F"); // mean DCAz - 5.0 GeV/c
    ttree->Branch("errmdcaz5",&errmdcaz5,"errmdcaz5/F"); // mean DCAz - 5.0 GeV/c
    ttree->Branch("mdcaz10",&mdcaz10,"mdcaz10/F"); // mean DCAz - 10 GeV/c
    ttree->Branch("errmdcaz10",&errmdcaz10,"errmdcaz10/F"); // mean DCAz - 10 GeV/c
    ttree->Branch("rmsdcaz05",&rmsdcaz05,"rmsdcaz05/F"); // rms DCAz - 0.5 GeV/c
    ttree->Branch("errrmsdcaz05",&errrmsdcaz05,"errrmsdcaz05/F"); // rms DCAz - 0.5 GeV/c
    ttree->Branch("rmsdcaz1",&rmsdcaz1,"rmsdcaz1/F"); // rms DCAz - 1.0 GeV/c
    ttree->Branch("errrmsdcaz1",&errrmsdcaz1,"errrmsdcaz1/F"); // rms DCAz - 1.0 GeV/c
    ttree->Branch("rmsdcaz5",&rmsdcaz5,"rmsdcaz5/F"); // rms DCAz - 5.0 GeV/c
    ttree->Branch("errrmsdcaz5",&errrmsdcaz5,"errrmsdcaz5/F"); // rms DCAz - 5.0 GeV/c
    ttree->Branch("rmsdcaz10",&rmsdcaz10,"rmsdcaz10/F"); // rms DCAz - 10 GeV/c
    ttree->Branch("errrmsdcaz10",&errrmsdcaz10,"errrmsdcaz10/F"); // rms DCAz - 10 GeV/c


    ///////////////   TTree filling
    ///////////////   Vertex part
    nrun=runNumber;
    
    TDirectoryFile * VertexQAdir=(TDirectoryFile*)fin->Get("Vertex_Performance");
    if (VertexQAdir) {
        printf("MESSAGE: Vertex QA directory found in input file.\n");
        myfile <<"MESSAGE: Vertex QA directory found in input file" << endl;
  
        TList * VertxList=(TList*)VertexQAdir->Get("cOutputVtxESD");
        if (VertxList) {
            printf("MESSAGE: Vertex QA TList found\n");
            myfile <<"MESSAGE: Vertex QA TList found" << endl;
            FillVertexBranches(VertxList);
            myfile <<"MESSAGE: Vertex branches filled" << endl;
            cout <<"MESSAGE: Vertex branches filled" << endl;
        }
    }
    
   ///////////////////////  SSD Part
    TDirectoryFile * SSDQAdir=(TDirectoryFile*)fin->Get("PWGPPdEdxSSDQA");
    if(!SSDQAdir) SSDQAdir=(TDirectoryFile*)fin->Get("PWG1dEdxSSDQA");
    if (SSDQAdir) {
        Printf("MESSAGE: SSD QA directory found in input file");
        myfile <<"MESSAGE: SSD QA directory found in input file" << endl;

        TList * SSDList=(TList*)SSDQAdir->Get("SSDdEdxQA");
        if (SSDList) {
            printf("MESSAGE: SSD QA TList found\n");
            myfile <<"MESSAGE: SSD QA TList found" << endl;
            FillSSDBranches(SSDList);
            myfile <<"MESSAGE: SSD branches filled" << endl;
            cout <<"MESSAGE: SSD branches filled" << endl;
        }
    
    }

    ///////////////////////  SDD Part
    TDirectoryFile * SDDQAdir=(TDirectoryFile*)fin->Get("SDD_Performance");
    if (SDDQAdir) {
        printf("MESSAGE: SDD QA directory found in input file.\n");
        myfile <<"MESSAGE: SDD QA directory found in input file" << endl;
        TList * SDDList=(TList*)SDDQAdir->Get("coutputRP");
        if(SDDList){
            printf("MESSAGE: SDD QA directory found in input file.\n");
            myfile <<"MESSAGE: SDD QA TList found" << endl;
            cout <<"MESSAGE: SDD QA TList found" << endl;
            FillSDDBranches(SDDList);
            myfile <<"MESSAGE: SDD branches filled" << endl;
            cout <<"MESSAGE: SDD branches filled" << endl;
        }
    }
    
    /////////////////////   SPD Part
    TDirectoryFile *SPDQAdir=(TDirectoryFile*)fin->GetDirectory("SPD_Performance");
    if(SPDQAdir) {
        printf("MESSAGE: SPD QA directory found in input file.\n");
        myfile <<"MESSAGE: SPD QA directory found in input file" << endl;
        TList * SPDList=(TList*)SPDQAdir->Get("coutput1");
        if(SPDList){
            myfile <<"MESSAGE: SPD QA TList found" << endl;
            cout <<"MESSAGE: SPD QA TList found" << endl;
            FillSPDBranches(SPDList);
            myfile <<"MESSAGE: SPD branches filled" << endl;
            cout <<"MESSAGE: SPD branches filled" << endl;
        }
    }

    /////////////////////   Matching Part
    TDirectoryFile *dirMatch=(TDirectoryFile*)fin->GetDirectory("ITS_Performance");
    if(dirMatch) {
        printf("MESSAGE: ITS_Performance QA directory found in input file.\n");
        myfile <<"MESSAGE: ITS_Performance QA directory found in input file" << endl;
        TList * ITSList=(TList*)dirMatch->Get("cOutputITS");
        if(ITSList){
            myfile <<"MESSAGE: ITS_Performance QA TList found" << endl;
            cout <<"MESSAGE: ITS_Performance QA TList found" << endl;
            FillMatchingBranches(ITSList);
            myfile <<"MESSAGE: ITS_Performance branches filled" << endl;
            cout <<"MESSAGE: ITS_Performance branches filled" << endl;
        }
    }

    /////////////////////   ITSsa Part
    TDirectoryFile *ITSsaQAdir=(TDirectoryFile*)fin->GetDirectory("ITSsaTracks");
    if(ITSsaQAdir) {
        printf("MESSAGE: ITSsa QA directory found in input file.\n");
        myfile <<"MESSAGE: ITSsa QA directory found in input file" << endl;
        TList * ITSsaList=(TList*)ITSsaQAdir->Get("clistITSsaTracks");
        if(ITSsaList){
            myfile <<"MESSAGE: ITSsa QA TList found" << endl;
            cout <<"MESSAGE: ITSsa QA TList found" << endl;
            FillITSsaBranches(ITSsaList);
            myfile <<"MESSAGE: ITSsa branches filled" << endl;
            cout <<"MESSAGE: ITSsa branches filled" << endl;
        }
    }
   
    /////////////////////   Pileup  Part
    TDirectoryFile *pileupQAdir=(TDirectoryFile*)fin->GetDirectory("CheckPileupQA");
    if(pileupQAdir) {
        printf("MESSAGE: pileup QA directory found in input file.\n");
        myfile <<"MESSAGE: pileup QA directory found in input file" << endl;
        TList * PuList=(TList*)pileupQAdir->Get("clistPileupSPDQA");
        if(PuList){
            myfile <<"MESSAGE: pileup QA TList found" << endl;
            cout <<"MESSAGE: pileup QA TList found" << endl;
            FillPileupBranches(PuList);
            myfile <<"MESSAGE: pileup branches filled" << endl;
            cout <<"MESSAGE: pileup branches filled" << endl;
        }
    }

    /////////////////////   PID  Part
    TDirectoryFile *PIDQAdir=(TDirectoryFile*)fin->GetDirectory("PIDqa");
    if(PIDQAdir) {
        printf("MESSAGE: PID QA directory found in input file.\n");
        myfile <<"MESSAGE: PID QA directory found in input file" << endl;
        
        TList *PIDList2 = (TList*)PIDQAdir->Get("PIDqa");
        if(PIDList2){
//            cout << "MESSAGE: intermediate QA TList found" << endl;
            TList * PIDList=(TList*)PIDList2->FindObject("ITS");
            if(PIDList){
                myfile <<"MESSAGE: PID QA TList found" << endl;
                cout <<"MESSAGE: PID QA TList found" << endl;
                FillPIDBranches(PIDList);
                myfile <<"MESSAGE: PID branches filled" << endl;
                cout <<"MESSAGE: PID branches filled" << endl;
            }
        }
    }

    /////////////////////   DCA  Part
    TDirectoryFile *DCAQAdir=(TDirectoryFile*)fin->GetDirectory("ImpParRes_Performance");
    if(DCAQAdir) {
        printf("MESSAGE: DCA QA directory found in input file.\n");
        myfile <<"MESSAGE: DCA QA directory found in input file" << endl;
        TList * DCAList=(TList*)DCAQAdir->Get("coutputd0allPointRec_0_1000000");
        if(DCAList){
            myfile <<"MESSAGE: DCA QA TList found" << endl;
            cout <<"MESSAGE: DCA QA TList found" << endl;
            FillDCABranches(DCAList);
            myfile <<"MESSAGE: DCA branches filled" << endl;
            cout <<"MESSAGE: DCA branches filled" << endl;
        }
    }

    ttree->Fill();
    printf("==============  Saving trending quantities in tree for run %i ===============\n",runNumber);
    trendFile->cd();
    ttree->Write();
    trendFile->Close();
  
    
    return  1;
}


/////////////// TTree filling methods
/////////////// Vertex
void FillVertexBranches(TList * VertxList){
    
    Printf("Vertex - QA");
    myfile << "Vertex - QA" << endl;
    
    TH1F *xVtxTRK = (TH1F*)VertxList->FindObject("fhTRKVertexX");
    TH1F *yVtxTRK = (TH1F*)VertxList->FindObject("fhTRKVertexY");
    TH1F *zVtxTRK = (TH1F*)VertxList->FindObject("fhTRKVertexZ");
    TH1F *xVtxSPD = (TH1F*)VertxList->FindObject("fhSPDVertexX");
    TH1F *yVtxSPD = (TH1F*)VertxList->FindObject("fhSPDVertexY");
    TH1F *zVtxSPD = (TH1F*)VertxList->FindObject("fhSPDVertexZ"); // pp runs
    if(zVtxSPD->GetEntries()==0) zVtxSPD = (TH1F*)VertxList->FindObject("fhSPDVertexZonly"); // PbPb runs!!!
    TH1F *zVtxSPDpil = (TH1F*)VertxList->FindObject("fhSPDVertexZPile");
    
    if(xVtxTRK){
     TF1 *fxTRK = new TF1("gausx", "gaus", -1, 1);
        if(xVtxTRK->GetEntries()>0){
            xVtxTRK->Fit("gausx", "NQRL");
            meanVtxTRKx=(Float_t)fxTRK->GetParameter(1);
            meanVtxTRKxErr=(Float_t)fxTRK->GetParError(1);
            sigmaVtxTRKx=(Float_t)fxTRK->GetParameter(2);
            sigmaVtxTRKxErr=(Float_t)fxTRK->GetParError(2);
        }
     }
    else
    {
        myfile << "xVtxTRK not found" << endl;
    }
    
    if(yVtxTRK){
     TF1 *fyTRK = new TF1("gausy", "gaus", -1, 1);
        if(yVtxTRK->GetEntries()>0){
            yVtxTRK->Fit("gausy","NQLR");
            meanVtxTRKy=(Float_t)fyTRK->GetParameter(1);
            meanVtxTRKyErr=(Float_t)fyTRK->GetParError(1);
            sigmaVtxTRKy=(Float_t)fyTRK->GetParameter(2);
            sigmaVtxTRKyErr=(Float_t)fyTRK->GetParError(2);
        }
    }
    else
    {
        myfile << "yVtxTRK not found" << endl;
    }
    
    if(zVtxTRK){
        TF1 *fzTRK = new TF1("gausz", "gaus", -1, 1);
        if(zVtxTRK->GetEntries()>0){
            zVtxTRK->Fit("gausz","NQRL");
            meanVtxTRKz=(Float_t)fzTRK->GetParameter(1);
            meanVtxTRKzErr=(Float_t)fzTRK->GetParError(1);
            sigmaVtxTRKz=(Float_t)fzTRK->GetParameter(2);
            sigmaVtxTRKzErr=(Float_t)fzTRK->GetParError(2);
        }
    }
    else
    {
        myfile << "zVtxTRK not found" << endl;
    }
    
    if(xVtxSPD){
        TF1 *fxSPD = new TF1("gausxSPD", "gaus", -1, 1);
        if(xVtxSPD->GetEntries()>0){
            xVtxSPD->Fit("gausxSPD", "NQRL");
            meanVtxSPDx=(Float_t)fxSPD->GetParameter(1);
            meanVtxSPDxErr=(Float_t)fxSPD->GetParError(1);
            sigmaVtxSPDx=(Float_t)fxSPD->GetParameter(2);
            sigmaVtxSPDxErr=(Float_t)fxSPD->GetParError(2);
        }
    }
    else
    {
        myfile << "xVtxSPD not found" << endl;
    }
    
    if(yVtxSPD){
        TF1 *fySPD = new TF1("gausySPD", "gaus", -1, 1);
        if(yVtxSPD->GetEntries()>0){
            yVtxSPD->Fit("gausySPD","NQLR");
            meanVtxSPDy=(Float_t)fySPD->GetParameter(1);
            meanVtxSPDyErr=(Float_t)fySPD->GetParError(1);
            sigmaVtxSPDy=(Float_t)fySPD->GetParameter(2);
            sigmaVtxSPDyErr=(Float_t)fySPD->GetParError(2);
        }
    }
    else
    {
        myfile << "yVtxSPD not found" << endl;
    }
    
    if(zVtxSPD){
        TF1 *fzSPD = new TF1("gauszSPD", "gaus", -1, 1);
        if(zVtxSPD->GetEntries()>0){
            zVtxSPD->Fit("gauszSPD","NQRL");
            meanVtxSPDz=(Float_t)fzSPD->GetParameter(1);
            sigmaVtxSPDzErr=(Float_t)fzSPD->GetParError(1);
            sigmaVtxSPDz=(Float_t)fzSPD->GetParameter(2);
            sigmaVtxSPDzErr=(Float_t)fzSPD->GetParError(2);
        }
    }
    else
    {
        myfile << "zVtxSPD not found" << endl;
    }

    if(zVtxSPDpil){
        if(zVtxSPDpil->GetEntries()>0){
            pileupSPD=(Float_t)zVtxSPDpil->GetEntries()/(Float_t)zVtxSPD->GetEntries();
            errpileupSPD=TMath::Sqrt((Float_t)zVtxSPDpil->GetEntries())/(Float_t)zVtxSPD->GetEntries();
            }
        }
    
    } /// end void FillVertexBranches(TList * VertxList)


   ///////////////////////  SSD

void FillSSDBranches(TList * SSDList){
 
    Printf("SSD - QA");
    myfile << "SSD - QA" << endl;
    
    TH2F* QAchargeRatio=(TH2F*)SSDList->FindObject("QAChargeRatioSA");
    if(!QAchargeRatio) QAchargeRatio=(TH2F*)SSDList->FindObject("QAChargeRatio");
    if(QAchargeRatio->GetEntries()==0){
        myfile <<"QAchargeRatio EMPTY" << endl;
    }
    
    TH2F* QAcharge=(TH2F*)SSDList->FindObject("QAChargeSA");
    if(!QAcharge) QAcharge=(TH2F*)SSDList->FindObject("QACharge");
    if(QAcharge->GetEntries()==0){
        myfile <<"QAcharge EMPTY" << endl;
        
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
        if(hChargeL5->GetEntries()>0){
            MPVL5=(Float_t)lfunLay5->GetParameter(1);
            MPVErrL5=(Float_t)lfunLay5->GetParError(1);}
        if(hChargeL5->GetEntries()>0){
            MPVL6=(Float_t)lfunLay6->GetParameter(1);
            MPVErrL6=(Float_t)lfunLay6->GetParError(1);}
        ChargeRatioL5=(Float_t)hChargeRatioL5->GetMean();
        ChargeRatioErrL5=(Float_t)hChargeRatioL5->GetMeanError();
        ChargeRatioL6=(Float_t)hChargeRatioL6->GetMean();
        ChargeRatioErrL6=(Float_t)hChargeRatioL6->GetMeanError();
        EmptyModulesSSD=(Float_t)contEmpty;
//
        TH1F* bad_p=(TH1F*)SSDList->FindObject("Bad-p-strips");
        TH1F* bad_n=(TH1F*)SSDList->FindObject("Bad-n-strips");
        // find number of merged subjobs
        Int_t max_p=0;
        Float_t bad_n5=0., bad_p5=0., bad_n6=0., bad_p6=0.;
        
        if(bad_p&&bad_n){
            Float_t n_subjobs=0;
            for(Int_t i=1;i<1699;i++) if(bad_p->GetBinContent(i)>max_p)max_p=bad_p->GetBinContent(i);
            n_subjobs=max_p/768.;
            // find the number of bad n/p strips per layer
            bad_p->Scale(1/n_subjobs);
            bad_n->Scale(1/n_subjobs);
            for(Int_t j=1; j<749; j++) {bad_n5=bad_n5+bad_n->GetBinContent(j); bad_p5=bad_p5+bad_p->GetBinContent(j);}
            for(Int_t j=749; j<1699; j++) {bad_n6=bad_n6+bad_n->GetBinContent(j); bad_p6=bad_p6+bad_p->GetBinContent(j);}
//
            FracBadn5=(Float_t)bad_n5/574464.;
            errFracBadn5=TMath::Sqrt((Float_t)bad_n5)/574464.;
            FracBadp5=(Float_t)bad_p5/574464.;
            errFracBadp5=TMath::Sqrt((Float_t)bad_p5)/574464.;
            FracBadn6=(Float_t)bad_n6/729600.;
            errFracBadn6=TMath::Sqrt((Float_t)bad_n6)/729600.;
            FracBadp6=(Float_t)bad_p6/729600.;
            errFracBadp6=TMath::Sqrt((Float_t)bad_p6)/729600.;
            }
        
        if(TMath::Abs(ChargeRatioL5)<0.01) FlagChR5=1.; else FlagChR5=0.;
        if(TMath::Abs(ChargeRatioL6)<0.01) FlagChR6=1.; else FlagChR6=0.;
        if(MPVL5>80. && MPVL5<85.) FlagdEdx5=1.; else FlagdEdx5=0.;
        if(MPVL6>80. && MPVL6<85.) FlagdEdx6=1.; else FlagdEdx6=0.;
        if(FracBadn5 < 0.2) FlagSSD1n=1.0; else FlagSSD1n=0.;
        if(FracBadn5==-999.)FlagSSD1n=-999.;
        if(FracBadp5 < 0.2) FlagSSD1p=1.0; else FlagSSD1p=0.;
        if(FracBadp5==-999.)FlagSSD1p=-999.;
        if(FracBadn6 < 0.2) FlagSSD2n=1.0; else FlagSSD2n=0.;
        if(FracBadn6==-999.)FlagSSD2n=-999.;
        if(FracBadp6 < 0.2) FlagSSD2p=1.0; else FlagSSD2p=0.;
        if(FracBadp6==-999.)FlagSSD2p=-999.;
        
    } // check on QAcharge and QAchargeRatio

} /// end void FillSSDBranches(TList * SSDList)


///////////////////////  SDD

void FillSDDBranches(TList * SDDList){
    
    Printf("SDD - QA");
    myfile << "SDD - QA" << endl;
    
    //
    Float_t fracT[6]={0.,0.,0.,0.,0.,0.};
    Float_t efracT[6]={0.,0.,0.,0.,0.,0.};
    Int_t nTotEvents;
    Int_t nTrigEvents;
    Double_t averPoints=0.;
    Double_t cntBins=0.;
    Float_t minTime=-999.;
    Float_t errMinTime=0.;

    TH1F* hcllay=(TH1F*)SDDList->FindObject("hCluInLay");
    if(hcllay){
        if(hcllay->GetBinContent(1)>0){
            for(Int_t iLay=0; iLay<6; iLay++){
                fracT[iLay]=hcllay->GetBinContent(iLay+2)/hcllay->GetBinContent(1);
                efracT[iLay]=TMath::Sqrt(fracT[iLay]*(1-fracT[iLay])/hcllay->GetBinContent(1));
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
        }
    }
    
    
    TH1F* hgamod=(TH1F*)SDDList->FindObject("hGAMod");
    if(hgamod){
        if(hgamod->GetEntries()>0){
            Int_t bestMod=0;
            Int_t deadMod3=0, deadMod4=0;
            for(Int_t iMod=0; iMod<260;iMod++){
                Int_t gda=(Int_t)hgamod->GetBinContent(iMod+1);
                if(gda>bestMod) bestMod=gda;
                if(gda == 0) {
                    if(iMod<84) deadMod3+=1;
                    else deadMod4+=1;
                }
            }
    
            Int_t nChunks=1;
            if(bestMod>512){
                nChunks=(Int_t)(bestMod/512.+0.5);
            }
            hgamod->Scale(1./nChunks);

            fracDead3 = 1.-(Float_t)deadMod3/84;
            fracDead4 = 1.-(Float_t)deadMod4/176;
            errfracDead3 = TMath::Sqrt(fracDead3*(1-fracDead3)/84);         //
            errfracDead4 = TMath::Sqrt(fracDead4*(1-fracDead4)/176);
            
            if(fracDead3 > 0.8)FlagSDD1=1.; else FlagSDD1=0.;
            if(fracDead4 > 0.75)FlagSDD2=1.; else FlagSDD2=0.;
        }
    }

    TH1F* hev=(TH1F*)SDDList->FindObject("hNEvents");
    if(hev){
        if(hev->GetEntries()>0){
            nTotEvents=hev->GetBinContent(2);
            nTrigEvents=hev->GetBinContent(3);
            nEvents=nTotEvents;
            nEventsTriggered=nTrigEvents;
        }
    }

    
    TH1F* htimT=(TH1F*)SDDList->FindObject("hDrTimTPAll");
    TH1F* htimTe=(TH1F*)SDDList->FindObject("hDrTimTPExtra");
    if(htimT && htimTe) {
        if(htimT->GetEntries()>0 && htimTe->GetEntries()>0){
            if(htimT->GetEntries()>0){
                fracExtra=htimTe->GetEntries()/htimT->GetEntries();
                errfracExtra=TMath::Sqrt(htimTe->GetEntries())/htimT->GetEntries();
            }
        }
    }
    
    if(htimT){
        if(htimT->GetEntries()>0){
            for(Int_t iBin=1; iBin<=htimT->GetNbinsX(); iBin++){
                Float_t tim=htimT->GetBinCenter(iBin);
                if(tim>2000. && tim<4000.){
                    averPoints+=htimT->GetBinContent(iBin);
                    cntBins+=1;
                }
            }
//            Double_t minTime=-999.;
//            Double_t errMinTime=0.;
            minTime=-999.;
            errMinTime=0.;
            if(cntBins>0){
                averPoints/=cntBins;
                for(Int_t iBin=1; iBin<=htimT->GetNbinsX(); iBin++){
                    if(htimT->GetBinContent(iBin)>0.5*averPoints){
                        minTime=htimT->GetBinCenter(iBin); // minDrTime =-999. sempre ??
                        errMinTime=0.5*htimT->GetBinWidth(iBin);
                        break;
                    }
                }
            }
            minDrTime=minTime;
            errminDrTime=errMinTime;
            if(minDrTime > 485.) FlagMinTime=1.; else FlagMinTime=0.;
            
            meanDrTime=htimT->GetMean();
            errmeanDrTime=htimT->GetMeanError();
            if(meanDrTime > 3150.) FlagMeanTime=1.; else FlagMeanTime=0.;
            
        }
    }
    
    TH2F* hdedxmod=(TH2F*)SDDList->FindObject("hdEdxVsMod");
    if(hdedxmod){
        if(hdedxmod->GetEntries()>0){
            TH1D* hdedxLay3=hdedxmod->ProjectionY("hdedxLay3",1,84);
            TH1D* hdedxLay4=hdedxmod->ProjectionY("hdedxLay4",85,260);
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
        
            MPVdEdxLay3=lfunLay3->GetParameter(1);
            errMPVdEdxLay3=lfunLay3->GetParError(1);
            MPVdEdxLay4=lfunLay4->GetParameter(1);
            errMPVdEdxLay4=lfunLay4->GetParError(1);
            
            if(MPVdEdxLay3 > 80. && MPVdEdxLay3 < 86.) FlagdEdx3=1.; else FlagdEdx3=0.;
            if(MPVdEdxLay4 > 80. && MPVdEdxLay4 < 86.) FlagdEdx4=1.; else FlagdEdx4=0.;
        }
    }
    
    TH1F* hSigTim0=(TH1F*)SDDList->FindObject("hSigTimeInt0");
    if(hSigTim0){
        if(hSigTim0->GetEntries()>0){
            TF1 *lfunTim0 = new TF1("LangausFunTim0",LangausFun,50.,300.,4);
            lfunTim0->SetParameter(0,5.);
            lfunTim0->SetParameter(1,80.);
            lfunTim0->SetParameter(2,hSigTim0->GetEntries()/10.);
            lfunTim0->SetParameter(3,10.);
            lfunTim0->SetParLimits(3,0.,20);
            hSigTim0->Fit(lfunTim0,"NQLR");
            MPVdEdxTB0=lfunTim0->GetParameter(1);
            errMPVdEdxTB0=lfunTim0->GetParError(1);
        }
    }
    
    TH1F* hSigTim5=(TH1F*)SDDList->FindObject("hSigTimeInt5");
    if(hSigTim5){
        if(hSigTim5->GetEntries()>0){
            TF1 *lfunTim5 = new TF1("LangausFunTim5",LangausFun,50.,300.,4);
            lfunTim5->SetParameter(0,5.);
            lfunTim5->SetParameter(1,80.);
            lfunTim5->SetParameter(2,hSigTim5->GetEntries()/10.);
            lfunTim5->SetParameter(3,10.);
            lfunTim5->SetParLimits(3,0.,20);
            hSigTim5->Fit(lfunTim5,"NQLR");
            MPVdEdxTB5=lfunTim5->GetParameter(1);
            errMPVdEdxTB5=lfunTim5->GetParError(1);
        }
    }

} /// end void FillSDDBranches(TList * SDDList)

//
///////////////////////  SPD

void FillSPDBranches(TList * SPDList){
    
    Printf("SPD - QA");
    myfile << "SPD - QA" << endl;
    
    TH1F* hFiredChip=(TH1F*)SPDList->FindObject("hFiredChip");
    if(hFiredChip){
        if(hFiredChip->GetBinContent(1)>0){
        
            Float_t ptbin=0;
            Int_t nHSsInner=0,nHSsOuter=0;
            for(Int_t i=0;i<400;i++) if(hFiredChip->GetBinContent(i)>0) nHSsInner++;
            for(Int_t i=400;i<1200;i++) if(hFiredChip->GetBinContent(i)>0) nHSsOuter++;
            nHSsInner = (Int_t)(nHSsInner/10);
            nHSsOuter = (Int_t)(nHSsOuter/10);
            
            FracSPD1=(Float_t)nHSsInner/40.;
            errFracSPD1=TMath::Sqrt(nHSsInner)/40.;
            FracSPD2=(Float_t)nHSsOuter/80.;
            errFracSPD2=TMath::Sqrt(nHSsOuter)/80.;
            
            if(FracSPD1 > 0.75) FlagSPD1=1.; else FlagSPD1=0.;
            if(FracSPD2 > 0.90) FlagSPD2=1.; else FlagSPD2=0.;
        }
    }

} /// end void FillSPDBranches(TList * SPDList)

/// Matching Part TPCITS && TOFITS

void FillMatchingBranches(TList * ITSList){

    printf("Matching - QA\n");
    myfile << "Matching - QA" << endl;

    //
    TH1F* hclmapMI=(TH1F*)ITSList->FindObject("fHistClusterMapITSMI"); // tracce ricostruite
    TH1F* hclMI=(TH1F*)ITSList->FindObject("fHistNclsITSMI"); //
    
    TH1F* hclmapSA=(TH1F*)ITSList->FindObject("fHistClusterMapITSSA"); // tracce SA
    TH1F* hclSA=(TH1F*)ITSList->FindObject("fHistNclsITSSA"); //
    
    if(hclmapMI->GetEntries()>0 && hclMI->GetEntries()>0){
       Float_t fracTMI[6]={0.,0.,0.,0.,0.,0.};
        Float_t efracTMI[6]={0.,0.,0.,0.,0.,0.};
            for(Int_t iLay=0; iLay<6; iLay++){
                fracTMI[iLay]=hclmapMI->GetBinContent(iLay+1)/(hclMI->GetEntries()-hclMI->GetBinContent(1));
                efracTMI[iLay]=TMath::Sqrt(fracTMI[iLay]*(1-fracTMI[iLay])/(hclMI->GetEntries()-hclMI->GetBinContent(1)));
            }
        FracTrackMI1 = fracTMI[0];
        errFracTrackMI1 = efracTMI[0];
        FracTrackMI2 = fracTMI[1];
        errFracTrackMI2 = efracTMI[1];
        FracTrackMI3 = fracTMI[2];
        errFracTrackMI3 = efracTMI[2];
        FracTrackMI4 = fracTMI[3];
        errFracTrackMI4 = efracTMI[3];
        FracTrackMI5 = fracTMI[4];
        errFracTrackMI5 = efracTMI[4];
        FracTrackMI6 = fracTMI[5];
        errFracTrackMI6 = efracTMI[5];
    }

    if(hclmapSA->GetEntries()>0 && hclSA->GetEntries()>0){
        Float_t fracTSA[6]={0.,0.,0.,0.,0.,0.};
        Float_t efracTSA[6]={0.,0.,0.,0.,0.,0.};
            for(Int_t iLay=0; iLay<6; iLay++){
                fracTSA[iLay]=hclmapSA->GetBinContent(iLay+1)/hclSA->GetEntries();
                efracTSA[iLay]=TMath::Sqrt(fracTSA[iLay]*(1-fracTSA[iLay])/hclSA->GetEntries());
            }
        FracTrackSA1 = fracTSA[0];
        errFracTrackSA1 = efracTSA[0];
        FracTrackSA2 = fracTSA[1];
        errFracTrackSA2 = efracTSA[1];
        FracTrackSA3 = fracTSA[2];
        errFracTrackSA3 = efracTSA[2];
        FracTrackSA4 = fracTSA[3];
        errFracTrackSA4 = efracTSA[3];
        FracTrackSA5 = fracTSA[4];
        errFracTrackSA5 = efracTSA[4];
        FracTrackSA6 = fracTSA[5];
        errFracTrackSA6 = efracTSA[5];
    }
//
    
    Float_t ptbin=0.;
    TH1F *fHistPtTPCInAcc = (TH1F*)ITSList->FindObject("fHistPtTPCInAcc");
    if(fHistPtTPCInAcc){
       if(fHistPtTPCInAcc->GetEntries()>0){
    
           TH1F *fHistPtITSMI6InAcc = (TH1F*)ITSList->FindObject("fHistPtITSMI6InAcc");
           TH1F *fHistPtITSMI5InAcc = (TH1F*)ITSList->FindObject("fHistPtITSMI5InAcc");
           TH1F *fHistPtITSMI4InAcc = (TH1F*)ITSList->FindObject("fHistPtITSMI4InAcc");
           TH1F *fHistPtITSMI3InAcc = (TH1F*)ITSList->FindObject("fHistPtITSMI3InAcc");
           TH1F *fHistPtITSMI2InAcc = (TH1F*)ITSList->FindObject("fHistPtITSMI2InAcc");
           TH1F *fHistPtITSMISPDInAcc = (TH1F*)ITSList->FindObject("fHistPtITSMISPDInAcc");
           TH1F *fHistPtITSMIoneSPDInAcc = (TH1F*)ITSList->FindObject("fHistPtITSMIoneSPDInAcc");
    
           TH1F *fHistPtITSMIge2InAcc = (TH1F*)fHistPtITSMI6InAcc->Clone("fHistPtITSMIge2InAcc");
           fHistPtITSMIge2InAcc->Add(fHistPtITSMI5InAcc);
           fHistPtITSMIge2InAcc->Add(fHistPtITSMI4InAcc);
           fHistPtITSMIge2InAcc->Add(fHistPtITSMI3InAcc);
           fHistPtITSMIge2InAcc->Add(fHistPtITSMI2InAcc);
           
           fHistPtITSMI6InAcc->Divide(fHistPtITSMI6InAcc,fHistPtTPCInAcc,1,1,"B");
           ptbin=fHistPtITSMI6InAcc->FindBin(0.201);
           Eff6Pt02=fHistPtITSMI6InAcc->GetBinContent(ptbin);
           errEff6Pt02=fHistPtITSMI6InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI6InAcc->FindBin(1.001);
           Eff6Pt1=fHistPtITSMI6InAcc->GetBinContent(ptbin);
           errEff6Pt1=fHistPtITSMI6InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI6InAcc->FindBin(10.001);
           Eff6Pt10=fHistPtITSMI6InAcc->GetBinContent(ptbin);
           errEff6Pt10=fHistPtITSMI6InAcc->GetBinError(ptbin);
           
           fHistPtITSMI5InAcc->Divide(fHistPtITSMI5InAcc,fHistPtTPCInAcc,1,1,"B");
           ptbin=fHistPtITSMI5InAcc->FindBin(0.201);
           Eff5Pt02=fHistPtITSMI5InAcc->GetBinContent(ptbin);
           errEff5Pt02=fHistPtITSMI5InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI5InAcc->FindBin(1.001);
           Eff5Pt1=fHistPtITSMI5InAcc->GetBinContent(ptbin);
           errEff5Pt1=fHistPtITSMI5InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI5InAcc->FindBin(10.001);
           Eff5Pt10=fHistPtITSMI5InAcc->GetBinContent(ptbin);
           errEff5Pt10=fHistPtITSMI5InAcc->GetBinError(ptbin);
    
           fHistPtITSMI4InAcc->Divide(fHistPtITSMI4InAcc,fHistPtTPCInAcc,1,1,"B");
           ptbin=fHistPtITSMI4InAcc->FindBin(0.201);
           Eff4Pt02=fHistPtITSMI4InAcc->GetBinContent(ptbin);
           errEff4Pt02=fHistPtITSMI4InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI4InAcc->FindBin(1.001);
           Eff4Pt1=fHistPtITSMI4InAcc->GetBinContent(ptbin);
           errEff4Pt1=fHistPtITSMI4InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI4InAcc->FindBin(10.001);
           Eff4Pt10=fHistPtITSMI4InAcc->GetBinContent(ptbin);
           errEff4Pt10=fHistPtITSMI4InAcc->GetBinError(ptbin);
    
           fHistPtITSMI3InAcc->Divide(fHistPtITSMI3InAcc,fHistPtTPCInAcc,1,1,"B");
           ptbin=fHistPtITSMI3InAcc->FindBin(0.201);
           Eff3Pt02=fHistPtITSMI3InAcc->GetBinContent(ptbin);
           errEff3Pt02=fHistPtITSMI3InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI3InAcc->FindBin(1.001);
           Eff3Pt1=fHistPtITSMI3InAcc->GetBinContent(ptbin);
           errEff3Pt1=fHistPtITSMI3InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI3InAcc->FindBin(10.001);
           Eff3Pt10=fHistPtITSMI3InAcc->GetBinContent(ptbin);
           errEff3Pt10=fHistPtITSMI3InAcc->GetBinError(ptbin);
    
           fHistPtITSMI2InAcc->Divide(fHistPtITSMI2InAcc,fHistPtTPCInAcc,1,1,"B");
           ptbin=fHistPtITSMI2InAcc->FindBin(0.201);
           Eff2Pt02=fHistPtITSMI2InAcc->GetBinContent(ptbin);
           errEff2Pt02=fHistPtITSMI2InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI2InAcc->FindBin(1.001);
           Eff2Pt1=fHistPtITSMI2InAcc->GetBinContent(ptbin);
           errEff2Pt1=fHistPtITSMI2InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMI2InAcc->FindBin(10.001);
           Eff2Pt10=fHistPtITSMI2InAcc->GetBinContent(ptbin);
           errEff2Pt10=fHistPtITSMI2InAcc->GetBinError(ptbin);
    
           fHistPtITSMISPDInAcc->Divide(fHistPtITSMISPDInAcc,fHistPtTPCInAcc,1,1,"B");
           ptbin=fHistPtITSMISPDInAcc->FindBin(0.201);
           EffSPDPt02=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
           errEffSPDPt02=fHistPtITSMISPDInAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMISPDInAcc->FindBin(1.001);
           EffSPDPt1=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
           errEffSPDPt1=fHistPtITSMISPDInAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMISPDInAcc->FindBin(10.001);
           EffSPDPt10=fHistPtITSMISPDInAcc->GetBinContent(ptbin);
           errEffSPDPt10=fHistPtITSMISPDInAcc->GetBinError(ptbin);
    
           fHistPtITSMIoneSPDInAcc->Divide(fHistPtITSMIoneSPDInAcc,fHistPtTPCInAcc,1,1,"B");
           ptbin=fHistPtITSMIoneSPDInAcc->FindBin(0.201);
           EffoneSPDPt02=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
           errEffoneSPDPt02=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMIoneSPDInAcc->FindBin(1.001);
           EffoneSPDPt1=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
           errEffoneSPDPt1=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMIoneSPDInAcc->FindBin(10.001);
           EffoneSPDPt10=fHistPtITSMIoneSPDInAcc->GetBinContent(ptbin);
           errEffoneSPDPt10=fHistPtITSMIoneSPDInAcc->GetBinError(ptbin);
    
           fHistPtITSMIge2InAcc->Divide(fHistPtITSMIge2InAcc,fHistPtTPCInAcc,1,1,"B");
           ptbin=fHistPtITSMIge2InAcc->FindBin(0.201);
           EffTOTPt02=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
           errEffTOTPt02=fHistPtITSMIge2InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMIge2InAcc->FindBin(1.001);
           EffTOTPt1=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
           errEffTOTPt1=fHistPtITSMIge2InAcc->GetBinError(ptbin);
           ptbin=fHistPtITSMIge2InAcc->FindBin(10.001);
           EffTOTPt10=fHistPtITSMIge2InAcc->GetBinContent(ptbin);
           errEffTOTPt10=fHistPtITSMIge2InAcc->GetBinError(ptbin);
        
       }
    }
    /////////// end of TPCITS part
    
    /////////// TOFITS part
    
    TH1F *fHistPtTPCInAccTOF = (TH1F*)ITSList->FindObject("fHistPtTPCInAccTOFbc0");
    if(fHistPtTPCInAccTOF){
        if(fHistPtTPCInAccTOF->GetEntries()>0){
            TH1F *fHistPtITSMI6InAccTOF = (TH1F*)ITSList->FindObject("fHistPtITSMI6InAccTOFbc0");
            TH1F *fHistPtITSMI5InAccTOF = (TH1F*)ITSList->FindObject("fHistPtITSMI5InAccTOFbc0");
            TH1F *fHistPtITSMI4InAccTOF = (TH1F*)ITSList->FindObject("fHistPtITSMI4InAccTOFbc0");
            TH1F *fHistPtITSMI3InAccTOF = (TH1F*)ITSList->FindObject("fHistPtITSMI3InAccTOFbc0");
            TH1F *fHistPtITSMI2InAccTOF = (TH1F*)ITSList->FindObject("fHistPtITSMI2InAccTOFbc0");
            TH1F *fHistPtITSMISPDInAccTOF = (TH1F*)ITSList->FindObject("fHistPtITSMISPDInAccTOFbc0");
            TH1F *fHistPtITSMIoneSPDInAccTOF = (TH1F*)ITSList->FindObject("fHistPtITSMIoneSPDInAccTOFbc0");

            TH1F *fHistPtITSMIge2InAccTOF = (TH1F*)fHistPtITSMI6InAccTOF->Clone("fHistPtITSMIge2InAccTOFbc0");
            fHistPtITSMIge2InAccTOF->Add(fHistPtITSMI5InAccTOF);
            fHistPtITSMIge2InAccTOF->Add(fHistPtITSMI4InAccTOF);
            fHistPtITSMIge2InAccTOF->Add(fHistPtITSMI3InAccTOF);
            fHistPtITSMIge2InAccTOF->Add(fHistPtITSMI2InAccTOF);
    
            fHistPtITSMI6InAccTOF->Divide(fHistPtITSMI6InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
            ptbin=fHistPtITSMI6InAccTOF->FindBin(0.501);
            Eff6Pt02TOF=fHistPtITSMI6InAccTOF->GetBinContent(ptbin);
            errEff6Pt02TOF=fHistPtITSMI6InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI6InAccTOF->FindBin(1.001);
            Eff6Pt1TOF=fHistPtITSMI6InAccTOF->GetBinContent(ptbin);
            errEff6Pt1TOF=fHistPtITSMI6InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI6InAccTOF->FindBin(10.001);
            Eff6Pt10TOF=fHistPtITSMI6InAccTOF->GetBinContent(ptbin);
            errEff6Pt10TOF=fHistPtITSMI6InAccTOF->GetBinError(ptbin);
    
            fHistPtITSMI5InAccTOF->Divide(fHistPtITSMI5InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
            ptbin=fHistPtITSMI5InAccTOF->FindBin(0.501);
            Eff5Pt02TOF=fHistPtITSMI5InAccTOF->GetBinContent(ptbin);
            errEff5Pt02TOF=fHistPtITSMI5InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI5InAccTOF->FindBin(1.001);
            Eff5Pt1TOF=fHistPtITSMI5InAccTOF->GetBinContent(ptbin);
            errEff5Pt1TOF=fHistPtITSMI5InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI5InAccTOF->FindBin(10.001);
            Eff5Pt10TOF=fHistPtITSMI5InAccTOF->GetBinContent(ptbin);
            errEff5Pt10TOF=fHistPtITSMI5InAccTOF->GetBinError(ptbin);
    
            fHistPtITSMI4InAccTOF->Divide(fHistPtITSMI4InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
            ptbin=fHistPtITSMI4InAccTOF->FindBin(0.501);
            Eff4Pt02TOF=fHistPtITSMI4InAccTOF->GetBinContent(ptbin);
            errEff4Pt02TOF=fHistPtITSMI4InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI4InAccTOF->FindBin(1.001);
            Eff4Pt1TOF=fHistPtITSMI4InAccTOF->GetBinContent(ptbin);
            errEff4Pt1TOF=fHistPtITSMI4InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI4InAccTOF->FindBin(10.001);
            Eff4Pt10TOF=fHistPtITSMI4InAccTOF->GetBinContent(ptbin);
            errEff4Pt10TOF=fHistPtITSMI4InAccTOF->GetBinError(ptbin);
    
            fHistPtITSMI3InAccTOF->Divide(fHistPtITSMI3InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
            ptbin=fHistPtITSMI3InAccTOF->FindBin(0.501);
            Eff3Pt02TOF=fHistPtITSMI3InAccTOF->GetBinContent(ptbin);
            errEff3Pt02TOF=fHistPtITSMI3InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI3InAccTOF->FindBin(1.001);
            Eff3Pt1TOF=fHistPtITSMI3InAccTOF->GetBinContent(ptbin);
            errEff3Pt1TOF=fHistPtITSMI3InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI3InAccTOF->FindBin(10.001);
            Eff3Pt10TOF=fHistPtITSMI3InAccTOF->GetBinContent(ptbin);
            errEff3Pt10TOF=fHistPtITSMI3InAccTOF->GetBinError(ptbin);
    
            fHistPtITSMI2InAccTOF->Divide(fHistPtITSMI2InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
            ptbin=fHistPtITSMI2InAccTOF->FindBin(0.501);
            Eff2Pt02TOF=fHistPtITSMI2InAccTOF->GetBinContent(ptbin);
            errEff2Pt02TOF=fHistPtITSMI2InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI2InAccTOF->FindBin(1.001);
            Eff2Pt1TOF=fHistPtITSMI2InAccTOF->GetBinContent(ptbin);
            errEff2Pt1TOF=fHistPtITSMI2InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMI2InAccTOF->FindBin(10.001);
            Eff2Pt10TOF=fHistPtITSMI2InAccTOF->GetBinContent(ptbin);
            errEff2Pt10TOF=fHistPtITSMI2InAccTOF->GetBinError(ptbin);
    
            fHistPtITSMISPDInAccTOF->Divide(fHistPtITSMISPDInAccTOF,fHistPtTPCInAccTOF,1,1,"B");
            ptbin=fHistPtITSMISPDInAccTOF->FindBin(0.501);
            EffSPDPt02TOF=fHistPtITSMISPDInAccTOF->GetBinContent(ptbin);
            errEffSPDPt02TOF=fHistPtITSMISPDInAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMISPDInAccTOF->FindBin(1.001);
            EffSPDPt1TOF=fHistPtITSMISPDInAccTOF->GetBinContent(ptbin);
            errEffSPDPt1TOF=fHistPtITSMISPDInAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMISPDInAccTOF->FindBin(10.001);
            EffSPDPt10TOF=fHistPtITSMISPDInAccTOF->GetBinContent(ptbin);
            errEffSPDPt10TOF=fHistPtITSMISPDInAccTOF->GetBinError(ptbin);
    
            fHistPtITSMIoneSPDInAccTOF->Divide(fHistPtITSMIoneSPDInAccTOF,fHistPtTPCInAccTOF,1,1,"B");
            ptbin=fHistPtITSMIoneSPDInAccTOF->FindBin(0.501);
            EffoneSPDPt02TOF=fHistPtITSMIoneSPDInAccTOF->GetBinContent(ptbin);
            errEffoneSPDPt02TOF=fHistPtITSMIoneSPDInAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMIoneSPDInAccTOF->FindBin(1.001);
            EffoneSPDPt1TOF=fHistPtITSMIoneSPDInAccTOF->GetBinContent(ptbin);
            errEffoneSPDPt1TOF=fHistPtITSMIoneSPDInAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMIoneSPDInAccTOF->FindBin(10.001);
            EffoneSPDPt10TOF=fHistPtITSMIoneSPDInAccTOF->GetBinContent(ptbin);
            errEffoneSPDPt10TOF=fHistPtITSMIoneSPDInAccTOF->GetBinError(ptbin);
    
            fHistPtITSMIge2InAccTOF->Divide(fHistPtITSMIge2InAccTOF,fHistPtTPCInAccTOF,1,1,"B");
            ptbin=fHistPtITSMIge2InAccTOF->FindBin(0.501);
            EffTOTPt02TOF=fHistPtITSMIge2InAccTOF->GetBinContent(ptbin);
            errEffTOTPt02TOF=fHistPtITSMIge2InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMIge2InAccTOF->FindBin(1.001);
            EffTOTPt1TOF=fHistPtITSMIge2InAccTOF->GetBinContent(ptbin);
            errEffTOTPt1TOF=fHistPtITSMIge2InAccTOF->GetBinError(ptbin);
            ptbin=fHistPtITSMIge2InAccTOF->FindBin(10.001);
            EffTOTPt10TOF=fHistPtITSMIge2InAccTOF->GetBinContent(ptbin);
            errEffTOTPt10TOF=fHistPtITSMIge2InAccTOF->GetBinError(ptbin);
        }
    }
    /////////// end of TOFITS part
    
} /// end void FillMatchingBranches(TList * ITSList)


   //////// ITSsa Part
void FillITSsaBranches(TList * ITSsaList){
    
    cout << "ITSsa - QA" << endl;
    myfile << "ITSsa - QA" << endl;
    
    Double_t Lowbin[3]={0.1,0.5,0.9};
    Double_t Upbin[3]={0.2,0.6,1};
    Double_t NTPCITS[3];
    Double_t NITSsa[3];
    Double_t NITSpureSA[3];
    Double_t Ratio[3];
    
    TH1F *hnev =(TH1F*)ITSsaList->FindObject("hNEvents");
    if(hnev){
        Double_t noEvents = hnev->GetBinContent(1);
        if(noEvents<1.)noEvents=1.;   // protection to avoid division by zero

        TH1F *hPtTPCITS=(TH1F*)ITSsaList->FindObject("hPtTPCITS");
        TH1F *hPtITSsa=(TH1F*)ITSsaList->FindObject("hPtITSsa");
        TH1F *hPtITSpureSA=(TH1F*)ITSsaList->FindObject("hPtITSpureSA");
        if(hPtTPCITS && hPtITSsa && hPtITSpureSA){
            for(Int_t ibin=0;ibin<=2;ibin++){
                NTPCITS[ibin]=hPtTPCITS->Integral(hPtTPCITS->FindBin(Lowbin[ibin]),hPtTPCITS->FindBin(Upbin[ibin]))/noEvents;
                NITSsa[ibin]=hPtITSsa->Integral(hPtITSsa->FindBin(Lowbin[ibin]),hPtITSsa->FindBin(Upbin[ibin]))/noEvents;
                NITSpureSA[ibin]=hPtITSpureSA->Integral(hPtITSpureSA->FindBin(Lowbin[ibin]),hPtITSpureSA->FindBin(Upbin[ibin]))/noEvents;
                Double_t totaltrks=NTPCITS[ibin]+NITSsa[ibin];
                if(totaltrks!=0 && NITSpureSA[ibin]!=0 )Ratio[ibin]=totaltrks/NITSpureSA[ibin];
                else Ratio[ibin]=0;
            }
            NITSTPCPtBin0=NTPCITS[0];
            NITSTPCPtBin1=NTPCITS[1];
            NITSTPCPtBin2=NTPCITS[2];
            NITSsaPtBin0=NITSsa[0];
            NITSsaPtBin1=NITSsa[1];
            NITSsaPtBin2=NITSsa[2];
            NITSpureSAPtBin0=NITSpureSA[0];
            NITSpureSAPtBin1=NITSpureSA[1];
            NITSpureSAPtBin2=NITSpureSA[2];
            ratioPtBin0=Ratio[0];
            ratioPtBin1=Ratio[1];
            ratioPtBin2=Ratio[2];
        }
    }
    
    TH1F* hNcluITSpSA =(TH1F*)ITSsaList->FindObject("hNcluITSpureSA"); // numero cluster se N>3
    if(hNcluITSpSA){
        NcluITSpSA=hNcluITSpSA->GetMean();
        errNcluITSpSA=hNcluITSpSA->GetMeanError();
    }
    
    TH1F* hdedx4 =(TH1F*)ITSsaList->FindObject("hdedxvsP4clsITSpureSA");
    TH1F* hdedx3 =(TH1F*)ITSsaList->FindObject("hdedxvsP3clsITSpureSA");
    if(hdedx3 && hdedx4){
        Float_t e4, e3;
        e4 = (Float_t)hdedx4->GetEntries();
        e3 = (Float_t)hdedx3->GetEntries();
        Float_t r43 = 0.;
        Float_t erre43=0.;
        if(hdedx3->GetEntries()>0){
            r43 =e4/e3;
            erre43=TMath::Sqrt(e4/(e3*e3)*(1+r43));
        }
        dedx4_3=r43;
        errdedx4_3=erre43;
    }
    
    TH1F* hPtpSA =(TH1F*)ITSsaList->FindObject("hPtITSpureSA");
    if(hPtpSA){
        PtpionpSA=hPtpSA->GetMean();
        errPtpionpSA=hPtpSA->GetMeanError();
    }
    

//    TH2F* hNcluPion2d = (TH2F*)cOutput->FindObject("hCluInLayITSpureSAPion");
//    TH1D* hNcluPion = hNcluPion2d->ProjectionY();
    TH2F* hNclu2d = (TH2F*)ITSsaList->FindObject("hCluInLayVsPtITSpureSA");
    TH1D* hNclu;
    if(hNclu2d) {
        hNclu = hNclu2d->ProjectionY();
        Float_t fracTpi[6]={0.,0.,0.,0.,0.,0.};
        Float_t efracTpi[6]={0.,0.,0.,0.,0.,0.};
        if(hNclu && hNclu->GetBinContent(1)>0){
            for(Int_t iLay=0; iLay<6; iLay++){
                fracTpi[iLay]=hNclu->GetBinContent(iLay+2)/hNclu->GetBinContent(1);
                efracTpi[iLay]=TMath::Sqrt(fracTpi[iLay]*(1-fracTpi[iLay])/hNclu->GetBinContent(1));
            }
        }
    
        NclupSA0=fracTpi[0];
        errNclupSA0=efracTpi[0];
        NclupSA1=fracTpi[1];
        errNclupSA1=efracTpi[1];
        NclupSA2=fracTpi[2];
        errNclupSA2=efracTpi[2];
        NclupSA3=fracTpi[3];
        errNclupSA3=efracTpi[3];
        NclupSA4=fracTpi[4];
        errNclupSA4=efracTpi[4];
        NclupSA5=fracTpi[5];
        errNclupSA5=efracTpi[5];
    }
    
    ///// chi^2 distribution mean
    TH1F* hChi2TPCITS=(TH1F*)ITSsaList->FindObject("hChi2TPCITS");
    if(hChi2TPCITS){
        chi2TPCITS=hChi2TPCITS->GetMean();
    }
    TH1F* hChi2ITSpureSA=(TH1F*)ITSsaList->FindObject("hChi2ITSpureSA");
    if(hChi2ITSpureSA){
        chi2ITSpureSA=hChi2ITSpureSA->GetMean();
    }
    
    ///// eta-phi TPCITS tracks distributions per layer
    
    TH2F* hEtaPhiTracksLay1TPCITS=(TH2F*)ITSsaList->FindObject("hEtaPhiTracksLay1TPCITS");
    TH2F* hEtaPhiTracksLay2TPCITS=(TH2F*)ITSsaList->FindObject("hEtaPhiTracksLay2TPCITS");
    TH2F* hEtaPhiTracksLay3TPCITS=(TH2F*)ITSsaList->FindObject("hEtaPhiTracksLay3TPCITS");
    TH2F* hEtaPhiTracksLay4TPCITS=(TH2F*)ITSsaList->FindObject("hEtaPhiTracksLay4TPCITS");
    TH2F* hEtaPhiTracksLay5TPCITS=(TH2F*)ITSsaList->FindObject("hEtaPhiTracksLay5TPCITS");
    TH2F* hEtaPhiTracksLay6TPCITS=(TH2F*)ITSsaList->FindObject("hEtaPhiTracksLay6TPCITS");
    TH1D* clus_phi_1, *clus_phi_2, *clus_phi_3, *clus_phi_4, *clus_phi_5, *clus_phi_6;
    TH1D* clus_eta_1, *clus_eta_2, *clus_eta_3, *clus_eta_4, *clus_eta_5, *clus_eta_6;

    if(hEtaPhiTracksLay1TPCITS){
        clus_phi_1 = hEtaPhiTracksLay1TPCITS->ProjectionY();
        clus_phi_1->Rebin(5); // 40 bins in phi, 9 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_1[phibin]=clus_phi_1->GetBinContent(phibin+1)/clus_phi_1->GetEntries();
            //            erocc_phi_1[phibin]= TMath::Sqrt(occ_phi_1[phibin]*(1-occ_phi_1[phibin])/clus_phi_1->GetEntries());
        }
        clus_eta_1 = hEtaPhiTracksLay1TPCITS->ProjectionX();
        clus_eta_1->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_1[etabin]=clus_eta_1->GetBinContent(etabin+1)/clus_eta_1->GetEntries();
            //            erocc_eta_1[etabin]= TMath::Sqrt(occ_eta_1[etabin]*(1-occ_eta_1[etabin])/clus_eta_1->GetEntries());
        }
    }

    if(hEtaPhiTracksLay2TPCITS){
        clus_phi_2 = hEtaPhiTracksLay2TPCITS->ProjectionY();
        clus_phi_2->Rebin(5); // 40 bins in phi, 9 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_2[phibin]=clus_phi_2->GetBinContent(phibin+1)/clus_phi_2->GetEntries();
        }
        clus_eta_2 = hEtaPhiTracksLay2TPCITS->ProjectionX();
        clus_eta_2->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_2[etabin]=clus_eta_2->GetBinContent(etabin+1)/clus_eta_2->GetEntries();
        }
    }
    
    if(hEtaPhiTracksLay3TPCITS){
        clus_phi_3 = hEtaPhiTracksLay3TPCITS->ProjectionY();
        clus_phi_3->Rebin(5); // 40 bins in phi, 9 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_3[phibin]=clus_phi_3->GetBinContent(phibin+1)/clus_phi_3->GetEntries();
        }
        clus_eta_3 = hEtaPhiTracksLay3TPCITS->ProjectionX();
        clus_eta_3->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_3[etabin]=clus_eta_3->GetBinContent(etabin+1)/clus_eta_3->GetEntries();
        }
    }
    
    if(hEtaPhiTracksLay4TPCITS){
        clus_phi_4 = hEtaPhiTracksLay4TPCITS->ProjectionY();
        clus_phi_4->Rebin(5); // 40 bins in phi, 9 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_4[phibin]=clus_phi_4->GetBinContent(phibin+1)/clus_phi_4->GetEntries();
        }
        clus_eta_4 = hEtaPhiTracksLay4TPCITS->ProjectionX();
        clus_eta_4->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_4[etabin]=clus_eta_4->GetBinContent(etabin+1)/clus_eta_4->GetEntries();
        }
    }
    
    if(hEtaPhiTracksLay5TPCITS){
        clus_phi_5 = hEtaPhiTracksLay5TPCITS->ProjectionY();
        clus_phi_5->Rebin(5); // 40 bins in phi, 9 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_5[phibin]=clus_phi_5->GetBinContent(phibin+1)/clus_phi_5->GetEntries();
        }
        clus_eta_5 = hEtaPhiTracksLay5TPCITS->ProjectionX();
        clus_eta_5->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_5[etabin]=clus_eta_5->GetBinContent(etabin+1)/clus_eta_5->GetEntries();
        }
    }
    
    if(hEtaPhiTracksLay6TPCITS){
        clus_phi_6 = hEtaPhiTracksLay6TPCITS->ProjectionY();
        clus_phi_6->Rebin(5); // 40 bins in phi, 9 degrees/bin
        for(Int_t phibin=0;phibin<40;phibin++){
            occ_phi_6[phibin]=clus_phi_6->GetBinContent(phibin+1)/clus_phi_6->GetEntries();
        }
        clus_eta_6 = hEtaPhiTracksLay6TPCITS->ProjectionX();
        clus_eta_6->Rebin(25); // 2 bins in phi, HS, HL
        for(Int_t etabin=0;etabin<2;etabin++){
            occ_eta_6[etabin]=clus_eta_6->GetBinContent(etabin+1)/clus_eta_6->GetEntries();
        }
    }

    
} /////////// end void FillITSsaBranches(TList * ITSsaList)
    
    
void FillPileupBranches(TList *PileUpList){
    
    cout << "Pileup - QA" << endl;
    myfile << "Pileup - QA" << endl;

    TH1F *hnPilVtx = (TH1F*)PileUpList->FindObject("hNOfPileupVertSPD");
    if(hnPilVtx){
    // Pileup vertices from SPD
        if(hnPilVtx->GetEntries()>0){
                Float_t meanPilVtx = 0.0, meanPilVtx_n=0;
                Float_t errmeanPilVtx = 0.0;
                Float_t errmeanPilVtx_n = 0.0;
                Float_t errmeanPilVtx_d = 0.0;
                for(Int_t i=2;i<12;i++){
                    meanPilVtx_n = meanPilVtx_n + (i-1)*hnPilVtx->GetBinContent(i);
                    meanPilVtx = meanPilVtx_n/hnPilVtx->Integral(2,12);
                    errmeanPilVtx_n = meanPilVtx + (i-1)*(i-1)*hnPilVtx->GetBinContent(i);
                }
                errmeanPilVtx_n = TMath::Sqrt(errmeanPilVtx_n);
                errmeanPilVtx_d = TMath::Sqrt(hnPilVtx->Integral(2,12));
                errmeanPilVtx = (errmeanPilVtx_n/meanPilVtx_n)*(errmeanPilVtx_n/meanPilVtx_n)+(errmeanPilVtx_d/hnPilVtx->Integral(2,12))*(errmeanPilVtx_d/hnPilVtx->Integral(2,12));
                errmeanPilVtx = TMath::Sqrt(errmeanPilVtx);

                npilvtx=meanPilVtx;
                errnpilvtx=errmeanPilVtx;
            }
        
        }
    
    TH1F *hnTrklPil = (TH1F*)PileUpList->FindObject("hNtracklPilSPD");
    if(hnTrklPil){
        if(hnTrklPil->GetEntries()>0){
                ntrklpil=hnTrklPil->GetMean();
                errntrklpil=hnTrklPil->GetMeanError();
            }

        }


    TH1F *hnTrklNoPil = (TH1F*)PileUpList->FindObject("hNtracklNoPilSPD");
    if(hnTrklNoPil){
    
    if(hnTrklNoPil->GetEntries()>0){
            ntrklnopil=hnTrklNoPil->GetMean();
            errntrklnopil=hnTrklNoPil->GetMeanError();
        }
    }
    
    TH1F *hnCl1Pil = (TH1F*)PileUpList->FindObject("hNCL1PilSPD");
    if(hnCl1Pil){
        if(hnCl1Pil->GetEntries()>0){
                ncl1pil=hnCl1Pil->GetMean();
                errncl1pil=hnCl1Pil->GetMeanError();
        }
    }
    
    TH1F *hnCl1NoPil = (TH1F*)PileUpList->FindObject("hNCL1NoPilSPD");
    if(hnCl1NoPil){
        if(hnCl1NoPil->GetEntries()>0){
            ncl1nopil=hnCl1NoPil->GetMean();
            errncl1nopil=hnCl1NoPil->GetMeanError();
        }
    }
    
} /////////// end of void FillPileupBranches(TList *PileUpList)
    

void FillPIDBranches(TList *PIDList){

    cout << "PID - QA" << endl;
    myfile << "PID - QA" << endl;

    Int_t bin[4]={25,57,81,105};
    Double_t NsigmaTPCITS[4]={100.,100.,100.,100.};
    Double_t Nsigma_err_TPCITS[4]={0.,0.,0.,0.};

    TH2F *nsigma_pi = (TH2F *)PIDList->FindObject("hNsigmaP_ITS_pion");
    Int_t biny=-1;
    Float_t binycont =-100;
    if(nsigma_pi){
       if(nsigma_pi->GetEntries()>0){
            for(Int_t l=0;l<4;l++){
                TH1D *nsigma_pi_1d=(TH1D *) nsigma_pi->ProjectionY("pro",bin[l],bin[l]);
                biny = nsigma_pi_1d->GetMaximumBin();
                NsigmaTPCITS[l] = nsigma_pi_1d->GetBinLowEdge(biny)+0.5*nsigma_pi_1d->GetBinWidth(biny);
                Nsigma_err_TPCITS[l] = 0.5*nsigma_pi_1d->GetBinWidth(biny);
            }
            nsigmapi02=NsigmaTPCITS[0];
            errnsigmapi02=Nsigma_err_TPCITS[0];
            nsigmapi05=NsigmaTPCITS[1];
            errnsigmapi05=Nsigma_err_TPCITS[1];
            nsigmapi1=NsigmaTPCITS[2];
            errnsigmapi1=Nsigma_err_TPCITS[2];
            nsigmapi3=NsigmaTPCITS[3];
            errnsigmapi3=Nsigma_err_TPCITS[3];
        }
    }
    
} /////////// end of void FillPIDBranches(TList *PIDList)

void FillDCABranches(TList *DCAList){
    
    printf("DCA - QA\n");
    myfile << "DCA - QA" << endl;

    Float_t dcamean[4]={-1000.,-1000.,-1000.,-1000.};
    Float_t err_dcamean[4]={0.01,0.01,0.01,0.01};
    Float_t dcaRMS[4]={0.0,0.0,0.0,0.0};
    Float_t err_dcaRMS[4]={0.0,0.0,0.0,0.0};
    
    Float_t dcazmean[4]={-1000.,-1000.,-1000.,-1000.};
    Float_t err_dcazmean[4]={0.01,0.01,0.01,0.01};
    Float_t dcazRMS[4]={0.0,0.0,0.0,0.0};
    Float_t err_dcazRMS[4]={0.0,0.0,0.0,0.0};
    
    TF1 *h = new TF1("h","gaus",-10000.,10000.);
    Double_t xmin=0.0, xmax=0.0;
    
    TH1F *hDCA_05 = (TH1F *)DCAList->FindObject("d0allpointrphiRec_5");
    if(hDCA_05){
        if(hDCA_05->GetEntries()>0){
            xmin=hDCA_05->GetMean()-3*hDCA_05->GetRMS();
            xmax=hDCA_05->GetMean()+3*hDCA_05->GetRMS();
            hDCA_05->Fit("h","NR,Q","",xmin,xmax);
        //     dcamean[0]=hDCA_05->GetMean();
        //     err_dcamean[0]=hDCA_05->GetMeanError();
        //     dcaRMS[0]=hDCA_05->GetRMS();
        //     err_dcaRMS[0]=hDCA_05->GetRMSError();
            dcamean[0]=h->GetParameter(1);
            err_dcamean[0]=h->GetParError(1);
            dcaRMS[0]=h->GetParameter(2);
            err_dcaRMS[0]=h->GetParError(2);
        }
    }
    
    TH1F *hDCA_1 = (TH1F *)DCAList->FindObject("d0allpointrphiRec_9");
    if(hDCA_1){
        if(hDCA_1->GetEntries()>0){
            xmin=hDCA_1->GetMean()-3*hDCA_1->GetRMS();
            xmax=hDCA_1->GetMean()+3*hDCA_1->GetRMS();
            hDCA_1->Fit("h","NR,Q","",xmin,xmax);
            dcamean[1]=h->GetParameter(1);
            err_dcamean[1]=h->GetParError(1);
            dcaRMS[1]=h->GetParameter(2);
            err_dcaRMS[1]=h->GetParError(2);
        }
    }
    
    TH1F *hDCA_5 = (TH1F *)DCAList->FindObject("d0allpointrphiRec_18");
    if(hDCA_5){
        if(hDCA_5->GetEntries()>0){
            xmin=hDCA_5->GetMean()-3*hDCA_5->GetRMS();
            xmax=hDCA_5->GetMean()+3*hDCA_5->GetRMS();
            hDCA_5->Fit("h","NR,Q","",xmin,xmax);
            dcamean[2]=h->GetParameter(1);
            err_dcamean[2]=h->GetParError(1);
            dcaRMS[2]=h->GetParameter(2);
            err_dcaRMS[2]=h->GetParError(2);
        }
    }
    
    TH1F *hDCA_10 = (TH1F *)DCAList->FindObject("d0allpointrphiRec_20");
    if(hDCA_10){
        if(hDCA_10->GetEntries()>0){
            xmin=hDCA_10->GetMean()-3*hDCA_10->GetRMS();
            xmax=hDCA_10->GetMean()+3*hDCA_10->GetRMS();
            hDCA_10->Fit("h","NR,Q","",xmin,xmax);
            dcamean[3]=h->GetParameter(1);
            err_dcamean[3]=h->GetParError(1);
            dcaRMS[3]=h->GetParameter(2);
            err_dcaRMS[3]=h->GetParError(2);
        }
    }
    
    //
    // DCAz
    //
    TH1F *hDCAz_05 = (TH1F *)DCAList->FindObject("d0allpointzRec_5");
    if(hDCAz_05){
        if(hDCAz_05->GetEntries()>0){
            xmin=hDCAz_05->GetMean()-3*hDCAz_05->GetRMS();
            xmax=hDCAz_05->GetMean()+3*hDCAz_05->GetRMS();
            hDCAz_05->Fit("h","NR,Q","",xmin,xmax);
            dcazmean[0]=h->GetParameter(1);
            err_dcazmean[0]=h->GetParError(1);
            dcazRMS[0]=h->GetParameter(2);
            err_dcazRMS[0]=h->GetParError(2);
        }
    }
    
    TH1F *hDCAz_1 = (TH1F *)DCAList->FindObject("d0allpointzRec_9");
    if(hDCAz_1){
        if(hDCAz_1->GetEntries()>0){
            xmin=hDCAz_1->GetMean()-3*hDCAz_1->GetRMS();
            xmax=hDCAz_1->GetMean()+3*hDCAz_1->GetRMS();
            hDCAz_1->Fit("h","NR,Q","",xmin,xmax);
            dcazmean[1]=h->GetParameter(1);
            err_dcazmean[1]=h->GetParError(1);
            dcazRMS[1]=h->GetParameter(2);
            err_dcazRMS[1]=h->GetParError(2);
        }
    }
    
    TH1F *hDCAz_5 = (TH1F *)DCAList->FindObject("d0allpointzRec_18");
    if(hDCAz_5){
        if(hDCAz_5->GetEntries()>0){
            xmin=hDCAz_5->GetMean()-3*hDCAz_5->GetRMS();
            xmax=hDCAz_5->GetMean()+3*hDCAz_5->GetRMS();
            dcazmean[2]=h->GetParameter(1);
            err_dcazmean[2]=h->GetParError(1);
            dcazRMS[2]=h->GetParameter(2);
            err_dcazRMS[2]=h->GetParError(2);
        }
    }
    
    TH1F *hDCAz_10 = (TH1F *)DCAList->FindObject("d0allpointzRec_20");
    if(hDCAz_10){
        if(hDCAz_10->GetEntries()>0){
            xmin=hDCAz_10->GetMean()-3*hDCAz_10->GetRMS();
            xmax=hDCAz_10->GetMean()+3*hDCAz_10->GetRMS();
            hDCAz_10->Fit("h","NR,Q","",xmin,xmax);
            dcazmean[3]=h->GetParameter(1);
            err_dcazmean[3]=h->GetParError(1);
            dcazRMS[3]=h->GetParameter(2);
            err_dcazRMS[3]=h->GetParError(2);
        }
    }
    
    mdca05=dcamean[0];
    errmdca05=err_dcamean[0];
    rmsdca05=dcaRMS[0];
    errrmsdca05=err_dcaRMS[0];
    mdca1=dcamean[1];
    errmdca1=err_dcamean[1];
    rmsdca1=dcaRMS[1];
    errrmsdca1=err_dcaRMS[1];
    mdca5=dcamean[2];
    errmdca5=err_dcamean[2];
    rmsdca5=dcaRMS[2];
    errrmsdca5=err_dcaRMS[2];
    mdca10=dcamean[3];
    errmdca10=err_dcamean[3];
    rmsdca10=dcaRMS[3];
    errrmsdca10=err_dcaRMS[3];

    mdcaz05=dcazmean[0];
    errmdcaz05=err_dcazmean[0];
    rmsdcaz05=dcazRMS[0];
    errrmsdcaz05=err_dcazRMS[0];
    mdcaz1=dcazmean[1];
    errmdcaz1=err_dcazmean[1];
    rmsdcaz1=dcazRMS[1];
    errrmsdcaz1=err_dcazRMS[1];
    mdcaz5=dcazmean[2];
    errmdcaz5=err_dcazmean[2];
    rmsdcaz5=dcazRMS[2];
    errrmsdcaz5=err_dcazRMS[2];
    mdcaz10=dcazmean[3];
    errmdcaz10=err_dcazmean[3];
    rmsdcaz10=dcazRMS[3];
    errrmsdcaz10=err_dcazRMS[3];

} /////////// end of void FillDCABranches(TList *DCAList)


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


