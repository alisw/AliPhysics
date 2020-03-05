/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */
#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TEventList.h"
#include "TObject.h"
#include "TNamed.h"
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTRDEfficiency.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <typeinfo>
#include "TApplication.h"

#include "AliAODv0.h"
#include "AliAODTrack.h"
/*
////////////////////////////////  don't know if this is necessary
#include "AliConversionPhotonCuts.h"

#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliDataFile.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "AliMCEvent.h"
#include "AliAODConversionMother.h"
#include "TObjString.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "TList.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliOADBContainer.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliTriggerAnalysis.h"
#include "AliV0ReaderV1.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliTRDTriggerAnalysis.h"
#include "AliDalitzAODESDMC.h"
#include "AliDalitzEventMC.h"
//////////////////////////////////////////
*/
#include "AliConversionPhotonCuts.h"
#include "AliConversionSelection.h"
#include "AliConversionTrackCuts.h"
#include "AliTRDTriggerAnalysis.h"
//#include "AliKFParticle.h"
#include "AliV0ReaderV1.h"
//#include "THnSparseD.h"

class AliAnalysisTRDEfficiency;    // your analysis class
//class AliConversionPhotonBase;
class AliConversionPhotonCuts;
//class AliConversionSelection;
class AliV0ReaderV1;

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTRDEfficiency) // classimp: necessary for root

AliAnalysisTRDEfficiency::AliAnalysisTRDEfficiency() : AliAnalysisTaskSE(), 
    fAOD(0), 
    fOutputList(0),
    fHistPt(0),

    file(0),
    fConversionGammas(0),

    fhm1pt2(0), 
    fhNpttr(0), 
    fhNptun(0),
    fhfA(0),
    
    fhR(0),
    fhRpt(0),
    fhMv(0),
    fhvpost(0),
    fhRhqu(0),
    fhRpthqu(0),
    fhMhqu(0),
    fhvposthqu(0),
        
    fhtxv(0),
    fhtyv(0),
    fhtzv(0),
    
    // gamma tracks
    fhgpt(0),
    fhgpttrd(0),
    fhgRpt(0),
    fhgRpttrd(0),
    
    fhgMinvM(0),
    fhgR(0),
    fhgptM(0),
    fhgptMhqu(0),
    
    fhgptQ(0),
    fhgptQhqu(0),
    fhgetaphi(0),
    fhgetaphihqu(0),
    
    fhgxy(0),
    fhgxyhqu(0),

    // v0 daughters
    fhdn(0),
    fhdpt(0),
    
    // n dimensional
    fhna(0),
    fhnp(0),
    fhnhqu(0),
    
    fhgevent2(0),
    fhgevent3(0),
    fhgevent4(0),
    fhgevent5(0),    
    fhgevent6(0),
    fhgevent7(0),
    fhgevent8(0),    
    fhgevent9(0),    
    
    fhgevent(0),
    fhevent(0),
    
    fhtrckvnt(0),   // track event
    fhtrckvnthqu(0), // track hqu event
    fhtrvnt(0),
    fhtrvnthqu(0),
    lsttrckvnt(0),
    lsttrckvnthqu(0),
    
    fhgetaphi1(0),
    fhgR1(0),
    fhgpt1(0),
    fhgetaphi5(0),
    fhgetaphi8(0),
    fhgetaphi9(0)
    
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTRDEfficiency::AliAnalysisTRDEfficiency(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), 
    fOutputList(0), 
    fHistPt(0),

    file(0),
    fConversionGammas(0),

    fhm1pt2(0), 
    fhNpttr(0), 
    fhNptun(0),
    fhfA(0),
    
    fhR(0),
    fhRpt(0),
    fhMv(0),
    fhvpost(0),
    fhRhqu(0),
    fhRpthqu(0),
    fhMhqu(0),
    fhvposthqu(0),
    
    
    fhtxv(0),
    fhtyv(0),
    fhtzv(0),
    
    // gamma tracks
    fhgpt(0),
    fhgpttrd(0),
    fhgRpt(0),
    fhgRpttrd(0),
    
    fhgMinvM(0),
    fhgR(0),
    fhgptM(0),
    fhgptMhqu(0),
    
    fhgptQ(0),
    fhgptQhqu(0),
    fhgetaphi(0),
    fhgetaphihqu(0),
    
    fhgxy(0),
    fhgxyhqu(0),

    // v0 daughters
    fhdn(0),
    fhdpt(0),
    
    // n dimensional
    fhna(0),
    fhnp(0),
    fhnhqu(0),

    fhgevent2(0),
    fhgevent3(0),
    fhgevent4(0),
    fhgevent5(0),
    fhgevent6(0),
    fhgevent7(0),
    fhgevent8(0),
    fhgevent9(0),

    fhgevent(0),
    fhevent(0),
    
    fhtrckvnt(0),   // track event
    fhtrckvnthqu(0), // track event
    fhtrvnt(0),
    fhtrvnthqu(0),
    lsttrckvnt(0),
    lsttrckvnthqu(0),
    
    fhgetaphi1(0),
    fhgR1(0),
    fhgpt1(0),
    fhgetaphi5(0),
    fhgetaphi8(0),
    fhgetaphi9(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTRDEfficiency::~AliAnalysisTRDEfficiency()
{
    // destructor
    if (fConversionGammas){
        cout << "the end of it " << fConversionGammas->GetEntries() << endl;
    }
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTRDEfficiency::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use 
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    cout << "start user create output objects" << endl;
    fOutputList = new TList();          // this is a list which will contain all of your histograms
                                        // at the end of the analysis, the contents of this list are written
                                        // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                        // if requested (dont worry about this now)
    //fConversionGammas = new TClonesArray("AliAODConversionPhoton",100);
    fConversionGammas = new TList();
    //fOutputList->Add(fConversionGammas);
    // example of a histogram
    fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);       // create your histogram
    fhm1pt2 = new TH1F("fhm1pt2","fhm1pt2", 10000, 0, 100);
    fhNpttr = new TH1F("fhNpttr","fhNpttr", 10000, 0, 100);
    fhNptun = new TH1F("fhNptun","fhNptun", 10000, 0, 100);
    fhfA    = new TH1F("fhfA", "fhfA", 200000, -100000, 100000);
    
    fhR     = new TH1F("fhR", "fhR", 5000, 0, 500);
    fhRpt   = new TH2F("fhRpt", "fhRpt", 1000, 0, 100, 5000, 0, 500);
    fhMv    = new TH1F("fhMv", "fhMv", 1000, 0, 1000);
    fhvpost = new TH1F("fhvpost", "fhvpost", 5000, 0, 500);
    fhRhqu  = new TH1F("fhRhqu", "fhRhqu", 5000, 0, 500);
    fhRpthqu= new TH2F("fhRpthqu", "fhRpthqu", 1000, 0, 100, 5000, 0, 500);
    fhMhqu  = new TH1F("fhMhqu", "fhMhqu", 1000, 0, 1000);
    fhvposthqu=new TH1F("fhvposthqu","fhvposthqu",5000,0,500);
    
    fhtxv   = new TH1F("fhtxv", "fhtxv", 4000, -1000, 1000);
    fhtyv   = new TH1F("fhtyv", "fhtyv", 4000, -1000, 1000);
    fhtzv   = new TH1F("fhtzv", "fhtzv", 4000, -1000, 1000);
    
    // gamma
    fhgpt    =new TH1F("fhgpt",     "fhgpt",    1000, 0, 100);
    fhgRpt   =new TH2F("fhgRpt",    "fhgRpt",   1000, 0, 1000, 100, 0, 100);
    fhgpttrd =new TH1F("fhgpttrd",  "fhgpttrd", 1000, 0, 100);
    fhgRpttrd=new TH2F("fhgRpttrd", "fhgRpttrd",1000, 0, 1000, 100, 0, 100);
    
    fhgMinvM = new TH2F("fhgMinvM", "fhgMinvM",     10000, 0, 10000, 10000, 0, 10000);
    fhgR     = new TH1F("fhgR", "fhgR",             1000, 0, 1000);
    fhgptM   = new TH2F("fhgptM", "fhgptM",         1000, 0, 100, 1000, 0, 1000);
    fhgptMhqu= new TH2F("fhgptMhqu", "fhgptMhqu",   1000, 0, 100, 1000, 0, 1000);
    
    fhgptQ      = new TH2F("fhgptQ", "fhgptQ",              1000, 0, 100, 1000, -1, 1);
    fhgptQhqu   = new TH2F("fhgptQhqu", "fhgptQhqu",        1000, 0, 100, 1000, -1, 1);
    fhgetaphi   = new TH2F("fhgetaphi", "fhgetaphi",        1000, -1, 4, 10000, -1, 9);
    fhgetaphihqu= new TH2F("fhgetaphihqu", "fhgetaphihqu",  1000, -1, 4, 10000, -1, 9);

    // v0 daughters
    fhdn        = new TH1F("fhdn", "fhdn",      2000, 0, 2000);
    fhdpt       = new TH1F("fhdpt", "fhdpt",    1000, 0, 200);
    
    // n dimensional histogram ( hqutracks: pt, eta, phi, a, b  
    
    Int_t dim = 6;          // dimensions
    Int_t bins[6]        = {2000, 300, 
                            2000, 1000, 
                            2000, 300};   // # of bins
	Double_t xmin[6]     = {-100, 0,
                            0,    0,
                            -100, 0};           //         min 
	Double_t xmax[6]     = {100, 500,
                            200, 300,
                            100, 300};  // max
    fhgevent2   = new THnSparseD("fhgevent2", "fhgevent2",         dim, bins, xmin, xmax);  // name, title, dim, nbins, xmin, xmax
    fhgevent3   = new THnSparseD("fhgevent3", "fhgevent3",         dim, bins, xmin, xmax);  // name, title, dim, nbins, xmin, xmax
    fhgevent4   = new THnSparseD("fhgevent4", "fhgevent4",         dim, bins, xmin, xmax);  // name, title, dim, nbins, xmin, xmax
    fhgevent5   = new THnSparseD("fhgevent5", "fhgevent5",         dim, bins, xmin, xmax);  // name, title, dim, nbins, xmin, xmax
    fhgevent6   = new THnSparseD("fhgevent6", "fhgevent6",         dim, bins, xmin, xmax);  // name, title, dim, nbins, xmin, xmax
    fhgevent7   = new THnSparseD("fhgevent7", "fhgevent7",         dim, bins, xmin, xmax);  // name, title, dim, nbins, xmin, xmax
    fhgevent8   = new THnSparseD("fhgevent8", "fhgevent8",         dim, bins, xmin, xmax);  // name, title, dim, nbins, xmin, xmax
    fhgevent9   = new THnSparseD("fhgevent9", "fhgevent9",         dim, bins, xmin, xmax);  // name, title, dim, nbins, xmin, xmax
    
    fhgevent2->GetAxis(0)->SetTitle("trd pT (GeV)");        // all 8 THnSparseD graphs have the same axis' names
    fhgevent2->GetAxis(1)->SetTitle("trd pid");
    fhgevent2->GetAxis(2)->SetTitle("gamma pt");
    fhgevent2->GetAxis(3)->SetTitle("gamma R");
    fhgevent2->GetAxis(4)->SetTitle("other trd pT (GeV)");
    fhgevent2->GetAxis(5)->SetTitle("other trd pT (GeV)");

    
    fhgevent    = new TH1F("fhgevent", "fhgevent",  10, 0, 10);
    fhevent     = new TH1F("fhevent", "fhevent",    10, 0, 10);
    
     dim = 4;          // dimensions
//     Int_t bins1[4]     = {10000, 10000};   // # of bins
//	 Int_t xmin1[4]     = {0, 0};           //         min
//	 Int_t xmax1[4]     = {10000, 10000};  // max
    //fhtrckvnt   = new THnSparseD("fhtrckvnt", "fhtrckvnt",          dim, bins1, xmin1, xmax1);
    //fhtrckvnthqu= new THnSparseD("fhtrckvnthqu", "fhtrckvnthqu",    dim, bins1, xmin1, xmax1);
    
    fhtrvnt     = new TH1F("fhtrvnt", "fhtrvnt",        20000, 0, 20000);
    fhtrvnthqu  = new TH1F("fhtrvnthqu", "fhtrvnthqu",  20000, 0, 20000);
    
    lsttrckvnt   = new TList();
    lsttrckvnthqu= new TList();
    
    fhgetaphi1  = new TH2F("fhgetaphi1", "fhgetaphi1",  400, -1, 7, 400, -2, 2);
    fhgR1       = new TH1F("fhgR1", "fhgR1",            2000, 0, 600);
    fhgpt1      = new TH1F("fhgpt1", "fhgpt1",          2000, 0, 200);
    fhgetaphi5  = new TH2F("fhgetaphi5", "fhgetaphi5",  800, -1, 7, 400, -2, 2);
    fhgetaphi8  = new TH2F("fhgetaphi8", "fhgetaphi8",  800, -1, 7, 400, -2, 2);
    fhgetaphi9  = new TH2F("fhgetaphi9", "fhgetaphi9",  800, -1, 7, 400, -2, 2);
    
    //fOutputList->Add(fHistPt);          // don't forget to add it to the list! the list will be written to file, so if you want
                                        // your histogram in the output file, add it to the list!
    //fOutputList->Add(fhm1pt2);
    //fOutputList->Add(fhNpttr);
    //fOutputList->Add(fhNptun);
    //fOutputList->Add(fhfA);
    
    // vertices
    /*fOutputList->Add(fhR);
    fOutputList->Add(fhRpt);
    fOutputList->Add(fhMv);
    fOutputList->Add(fhvpost);
    fOutputList->Add(fhRhqu);
    fOutputList->Add(fhRpthqu);
    fOutputList->Add(fhMhqu);
    fOutputList->Add(fhvposthqu);
    */
    //fOutputList->Add(fhtxv);
    //fOutputList->Add(fhtyv);
    //fOutputList->Add(fhtzv);    
    
    // gamma
    //fOutputList->Add(fhgpt);
    //fOutputList->Add(fhgRpt);
    //fOutputList->Add(fhgpttrd);
    //fOutputList->Add(fhgRpttrd);
    
    //fOutputList->Add(fhgR);
    fOutputList->Add(fhgptM);
    fOutputList->Add(fhgptMhqu);

    //fOutputList->Add(fhgptQ);
    //fOutputList->Add(fhgptQhqu);
    //fOutputList->Add(fhgetaphi);
    //fOutputList->Add(fhgetaphihqu);
    
    // v0 daughters
    //fOutputList->Add(fhdn);
    //fOutputList->Add(fhdpt);
    
    // n dimensional
    //fOutputList->Add(fhna);
    //fOutputList->Add(fhnp);
    //fOutputList->Add(fhnhqu);
    fOutputList->Add(fhgevent2);
    fOutputList->Add(fhgevent3);
    fOutputList->Add(fhgevent4);
    fOutputList->Add(fhgevent5);
    fOutputList->Add(fhgevent6);
    fOutputList->Add(fhgevent7);
    fOutputList->Add(fhgevent8);
    fOutputList->Add(fhgevent9);
    
    fOutputList->Add(fhgevent);
    fOutputList->Add(fhevent);
    
    //fOutputList->Add(fhtrckvnt);
    //fOutputList->Add(fhtrckvnthqu);
    fOutputList->Add(fhtrvnt);
    fOutputList->Add(fhtrvnthqu);
    //fOutputList->Add(lsttrckvnt);
    //fOutputList->Add(lsttrckvnthqu);
    
    fOutputList->Add(fhgetaphi1);
    fOutputList->Add(fhgR1);
    fOutputList->Add(fhgpt1);
    fOutputList->Add(fhgetaphi5);
    fOutputList->Add(fhgetaphi8);
    fOutputList->Add(fhgetaphi9);
    //file = new TFile("AliAODGammaConversion.root", "READ");
    //if (!file) {cout << "well i'm out of ideas" << endl;}
    //TString fDeltaAODFilename("AliAODGammaConversion.root");
    //AddAODBranch("TClonesArray", &fConversionGammas, fDeltaAODFilename.Data());

    //AliAnalysisManager::GetAnalysisManager()->RegisterExtraFile(fDeltaAODFilename.Data());
    
    //cout << "application " << TApplication::InputFiles() << endl;
    
    //cout << "put things into the output list" << endl;
    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
                                        // fOutputList object. the manager will in the end take care of writing your output to file
                                        // so it needs to know what's in the output
    //cout << "finish user stuff" << endl;
}

//void AliAnalysisTRDEfficiency::SetGammaConversionFile(TString *name)
//{
//    return;
//}
/*
 //_______________________________________________________________________
 TFile* AliAnalysisTRDEfficiency::OpenDigitsFile(TString inputfile,
                                             TString digfile,
                                             TString opt)
 {
    // we should check if we are reading ESDs or AODs - for now, only
    //  are supported
  
   if (digfile == "") {
     return NULL;
   }
 
   // construct the name of the digits file from the input file
   inputfile.ReplaceAll("AliESDs.root", digfile);
 
   // open the file
   AliInfo( "opening digits file " + inputfile
            + " with option \"" + opt + "\"");
   cout << inputfile << endl;
   inputfile = "alien://" + inputfile;
   TFile* dfile = new TFile(inputfile, opt);
   if (!dfile) {
     AliWarning("digits file '" + inputfile + "' cannot be opened");
  }
 
   return dfile;
}
*/



//_______________________________________________________________________________

Bool_t AliAnalysisTRDEfficiency::GetAODConversionGammas(AliAODEvent* fAODEvent){
            // yeeted from the https://github.com/alisw/AliPhysics/blob/master/PWGGA/GammaConvBase/AliV0ReaderV1.cxx


    // Get reconstructed conversion photons from satellite AOD file
    Bool_t hasgamma = kFALSE;
    Bool_t hasmatch = kFALSE;
    Bool_t haspt    = kFALSE;
    Bool_t haspid   = kFALSE;
    Bool_t haspitpt = kFALSE;
    Bool_t hashqu   = kFALSE;
    Bool_t hqupt    = kFALSE;
    Bool_t hqupid   = kFALSE;
    Bool_t hquptpid = kFALSE;
    
    Bool_t countedgamma = kFALSE;

    //cout << "start" << endl;
    TString fDeltaAODBranchName("GammaConv_00000003_06000008d00100001100000000_gamma");
 
    if (!fAODEvent->FindListObject(fDeltaAODBranchName.Data())) return kFALSE;
    
    TClonesArray *fInputGammas=dynamic_cast<TClonesArray*>(fAODEvent->FindListObject(fDeltaAODBranchName.Data()));

    if(!fInputGammas){AliError("No Gamma Satellites found");return kFALSE;}

    // Apply Selection Cuts to Gammas and create local working copy
    //if(fInputGammas){
    //cout << "the clones arrays " << fInputGammas->GetEntriesFast() << endl;
    //ofstream ofile;
    //ofile.open("file.txt", ios::app);
    //ofile << "New event" << endl;
    //ofile.close();
    for(Int_t i=0;i<fInputGammas->GetEntriesFast();i++){

        AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fInputGammas->At(i));  // get all of the gammas 
        if (!gamma) continue;
        
        if (!hasgamma){
            hasgamma = kTRUE;
            fhgevent->Fill(1);   
        }
        
        //if(gamma){
            
            //cout << " //////////////////////////////////////////////////////////" << endl;
            //if(fRelabelAODs)RelabelAODPhotonCandidates(gamma);
            //if(fConversionCuts->PhotonIsSelected(gamma,fInputEvent))
            //{new((*fConversionGammas)[fConversionGammas->GetEntriesFast()]) AliAODConversionPhoton(*gamma);}
        //fConversionGammas->Add(gamma);
        //fhgpt->Fill(gamma->Pt());
        //fhgRpt->Fill(gamma->GetConversionRadius(), gamma->Pt());
        //fhgR->Fill(gamma->GetConversionRadius());
                    
        //fhgptM->Fill(gamma->Pt(), gamma->GetPhotonMass());
        //fhgptQ->Fill(gamma->Pt(), gamma->GetPhotonQuality());
        //fhgetaphi->Fill(gamma->GetPhotonEta(), gamma->GetPhotonPhi());
                    
        // get vertex
        AliAODv0 *v0 = fAODEvent->GetV0(gamma->GetV0Index());
                    
            // get the duaghter particles
            //cout << "get number of daughters? " << v0->GetNDaughters() << endl;  // why are there more than two duaghter particles sometimes
        //fhdn->Fill(v0->GetNDaughters());
                    
            // gamma tracks
        Double_t ptg = gamma->Pt();
        Double_t etag= gamma->GetPhotonEta();
        Double_t phig= gamma->GetPhotonPhi();
        Double_t Rg  = gamma->GetConversionRadius();
                    
//        AliConversionPhotonCuts *cuts = new AliConversionPhotonCuts("00000003_06000008d00100001100000000");
                    
        fhgetaphi1->Fill(phig, etag);
        fhgR1->Fill(Rg);
        fhgpt1->Fill(ptg);
        
        for (Int_t j = 0; j < v0->GetNDaughters(); j++) {     // 
                        
                //cout << v0->GetDaughter(j) << endl;
            AliAODTrack *track = (AliAODTrack*)v0->GetDaughter(j);
                        
            // v0 daughters
            if (!track) continue;
                
            //fhdpt->Fill(track->Pt());
                //fhn->FillBin(track->Pt(), 1);

                // aodtracks                        
//            Double_t lst[15] = {0, 0, 0, 0,    // { hqutrack->Pt(), hqutrack->Eta(),        hqutrack->Phi(),        hqutrack->GetA(), hqutrack->GetB(),
//                                0, 0, 0, 0,    //   trdtrack->Pt(), trdtrack->Eta(),        trdtrack->Phi(),        trdtrack->GetA(), trdtrack->GetB(),
//                                0, 0, 0,       //   track->Pt(),    track->Eta(),           track->Phi(),
//                                0, 0, 0, 0};   //   gamma->Pt(),    gamma->GetPhotonEta(),  gamma->GetPhotonPhi(),  gamma->GetConversionRadius() };
                
                
            for (Int_t k = 0; k < fAODEvent->GetNumberOfTrdTracks(); k++){  // get the corresponding trd track
                AliAODTrdTrack *trdtrack = fAODEvent->GetTrdTrack(k);
                            
                if (!trdtrack) continue;
                if (trdtrack->GetTrackMatch() == (AliVTrack*)track){
                    //cout << "+++++++++++++++++++++++++++  we have a match, yippy  +++++++++++++++++++++++++++++++++++++" << endl;
                        //cout << "trdtrack " << trdtrack << endl;
                        //cout << "===================== track PID " << trdtrack->GetPID() << endl;
                    Double_t trdpt = trdtrack->Pt();
                    Double_t trdpid= trdtrack->GetPID();
                    Double_t gpt   = gamma->Pt();
                    Double_t gR    = gamma->GetConversionRadius();
                    Double_t dpt   = 0;
                    Double_t dpid  = 0;
                    
                    countedgamma = haspt + haspid + haspitpt + hashqu + hqupt + hqupid + hquptpid;
                    if (countedgamma){  // prevents double counting; numbers are supposed to be out of range
                        ptg = -1;
                        etag = 100;
                        phig = 100;
                        Rg = -1;
                    }
                    
                    // grab the other daughter track
                    AliAODTrack *other = (AliAODTrack*)v0->GetDaughter(!j); 
                    if (other){
                        for (Int_t l = 0; l < fAODEvent->GetNumberOfTrdTracks(); l++){
                            AliAODTrdTrack *trdother = fAODEvent->GetTrdTrack(l);
                            
                            if (!trdother) continue;
                            if (trdother->GetTrackMatch() == (AliVTrack*)other){
                                dpt = trdother->Pt();
                                dpid= trdother->GetPID();
                                break;  // break out of this loop
                            }
                        }
                    }
                
                    // but all of the variables in the list
                    Double_t lst[6]= {trdpt, trdpid, gpt, gR, dpt, dpid};
                    fhgptM->Fill(gamma->Pt(), gamma->GetPhotonMass());
                    fhgevent2->Fill(lst);
                    if (!hasmatch){    // event has matched trd track
                        hasmatch= kTRUE;
                        //fhgptM->Fill(gamma->Pt(), gamma->GetPhotonMass());
                        fhgevent->Fill(2);
                        //fhgevent2->Fill(lst);
                    }
                    if ((TMath::Abs(trdtrack->Pt()) >= 2.) ){              // e_pt > 2
                        if (!haspt) fhgevent->Fill(3);
                        haspt = kTRUE;
                        //fhgevent->Fill(3);
                        fhgevent3->Fill(lst);
                    }
                    if (trdtrack->GetPID() >= 164 ){     // this number was gotten from  ~PWG/TRD/AliTRDTriggerAnalysis.cxx
                        if (!haspid) fhgevent->Fill(4);
                        haspid = kTRUE;
                        //fhgevent->Fill(4);
                        fhgevent4->Fill(lst);
                            //ofstream ofile;
                            //ofile.open("file.txt", ios::app);
                            //ofile << "RRRR";
                            //ofile << "SSSS";
                            //ofile.close();
                    }
                    if ((TMath::Abs(trdtrack->Pt()) >= 2.) && (trdtrack->GetPID() >= 164) ){ // needs to be modified
                        if (!haspitpt){ // counts event once
                            fhgevent->Fill(5);
                            TString tmp[2] = { ((AliAODHeader*) fAODEvent->GetHeader() )->GetESDFileName(), ((AliAODHeader*) fAODEvent->GetHeader() )->GetEventNumberESDFile() };
                            //fhtrckvnt->Fill( { ((AliAODHeader*) fAODEvent->GetHeader() )->GetESDFileName(), ((AliAODHeader*) fAODEvent->GetHeader() )->GetEventNumberESDFile() } );
                            fhtrvnt->Fill(tmp[0].Append(tmp[1]), 1);
                            //fhgevent5->Fill(lst);
                        }
                        fhgetaphi5->Fill(phig, etag);
                        fhgevent5->Fill(lst);
                        haspitpt = kTRUE;
                        
                    }
                       
                    if ( (fAODEvent->GetFiredTriggerClasses()).Contains("HQU") ){
                                // fill this histogram
                        //cout << index << endl;
                        fhgptMhqu->Fill(gamma->Pt(), gamma->GetPhotonMass());
                        fhgevent6->Fill(lst);

                        if (!hashqu){
                            hashqu = kTRUE;
                            //fhgptMhqu->Fill(gamma->Pt(), gamma->GetPhotonMass());
                            fhgevent->Fill(6);
                            //fhgevent6->Fill(lst);
                        }
                        if ( (TMath::Abs(trdtrack->Pt()) >= 2.)){
                            if (!hqupt) fhgevent->Fill(7);
                            hqupt = kTRUE;
                            //fhgevent->Fill(7);
                            fhgevent7->Fill(lst);
                        }
                        if (trdtrack->GetPID() >= 164){
                            if (!hqupid) fhgevent->Fill(8);
                            hqupid = kTRUE;
                            //fhgevent->Fill(8);
                            fhgevent8->Fill(lst);
                        }
                        if ((TMath::Abs(trdtrack->Pt()) >= 2. && trdtrack->GetPID() >= 164)){
                            if (!hquptpid){ // counts event once
                                TString tmp[2] = { ((AliAODHeader*) fAODEvent->GetHeader() )->GetESDFileName(), ((AliAODHeader*) fAODEvent->GetHeader() )->GetEventNumberESDFile() };
                                //fhtrckvnthqu->Fill( tmp[0]+tmp[1], 1 );
                                fhtrvnthqu->Fill(tmp[0].Append(tmp[1]), 1);
                                //lsttrckvnthqu->Add(tmp[0], tmp[0].Data());
                                fhgevent->Fill(9);
                            }
                            fhgetaphi9->Fill(phig, etag);
                            hquptpid = kTRUE;
                            //cout << "a unicorn" << endl;    // print special track more than once
                            fhgevent9->Fill(lst);
                            
                        }
                        
                        //cout << "also has hqu" << endl;
                                
    
                        //fhnhqu->Fill(lst);
                    } // hqu
                    //fhnp->Fill(lst);
                } // match tracks
            } // trd tracks
                        
            //fhna->Fill(lst);
                        
        } // aod tracks
                    
    }
    //}
    //}
    //cout << "??????????????????????????????? " << fConversionGammas->GetEntries() << endl;
    //if(fConversionGammas->GetEntries()){return kTRUE;}

    return kFALSE;
}


//Bool_t AliAnalysisTRDEfficiency::GetSimilarVertex();
//{
//    return kTRUE;
//}


//_____________________________________________________________________________
void AliAnalysisTRDEfficiency::UserExec(Option_t *)
{
    //cout << "so the debugging begins" << endl;
    
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you 
    // have access to the current event. 
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
                                                        // there's another event format (ESD) which works in a similar wya
                                                        // but is more cpu/memory unfriendly. for now, we'll stick with aod's
    if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
        // example part: i'll show how to loop over the tracks in an event 
        // and extract some information from them which we'll store in a histogram
    //cout << "got events or whatever" << endl;
    //AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    //if (!fESD) return;
    //if (fESD) cout << "oh lawd this actually worked" << endl;
    //cout << "esd file name    " << ((AliAODHeader*) fAOD->GetHeader() )->GetESDFileName() << endl;
    //cout << "esd event number " << ((AliAODHeader*) fAOD->GetHeader() )->GetEventNumberESDFile() << endl;
    //cout << "run number       " << fAOD->GetRunNumber() << endl;
    //cout << "peroid number    " << fAOD->GetPeriodNumber() << endl;
    //cout << "orbit number     " << fAOD->GetOrbitNumber() << endl;
    //cout << "event type       " << fAOD->GetEventType() << endl;
    TString clss = fAOD->GetFiredTriggerClasses();
    if (!clss.Contains("CINT7-T")) return;          // Only grab these min. bias events
    //if (clss.Contains("HQU")){ cout << "hqu " << endl; fhgevent->Fill(6); }
    fhgevent->Fill(0);       // all events are inputted
    fhevent->Fill(0);
    
    GetAODConversionGammas(fAOD);
    
    //AliKFConversionPhoton *tmp = new AliKFConversionPhoton();

    //AliConversionSelection* acs = new AliConversionSelection("00010113", "00200009227300008250400000", "0152101500000000");  // aliconversion, will this work?
    
    //cout << "Does this thing exist " << acs << " and event number " << acs->GetEventNumber(fAOD) << endl;
    //cout << "get cut number" << acs->GetCutString() << endl;
    //cout << "number of good photons " << acs->GetGoodGammas() << endl; //->GetNumberOfPhotons() << endl;
   
//    TClonesArray *photons = new TClonesArray();
    //photons->Add();
    //Bool_t bl = acs->ProcessEvent( photons, fAOD, NULL);
    //cout << "ending boolean " << bl << endl;
    
    Int_t tracks    = fAOD->GetNumberOfTracks();
    Int_t trdtracks = fAOD->GetNumberOfTrdTracks();
//    Int_t v0tracks  = fAOD->GetNumberOfV0s();
   
    //cout << "v0tracks  " << v0tracks << endl;
    //cout << "vzero     " << fAOD->GetVZEROData() << endl;
    //cout << "trdtracks " << trdtracks << endl;
   
    //GetNumberOfPhotons();
   
    //Photons(fAOD);   // vertices
   

//    Bool_t hasgamma = kFALSE;
    Bool_t hasmatch = kFALSE;
    Bool_t haspt    = kFALSE;
    Bool_t haspid   = kFALSE;
    Bool_t haspitpt = kFALSE;
    Bool_t hashqu   = kFALSE;
    Bool_t hqupt    = kFALSE;
    Bool_t hqupid   = kFALSE;
    Bool_t hquptpid = kFALSE;
    
    /////////////////////  tracks  ////////////////////////
    for (Int_t i = 0; i < tracks; i++){   // all events
       
       AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
       if (!track) continue;
       
       //cout << "++++++++++++++++++== track pid " << track->PID() << "  " << track->GetID() << endl;
       //cout << "mass " << track->GetMassForTracking() << endl;
       //cout << "poisition " << track->Xv() << ", " << track->Yv() << ", " << track->Zv() << endl;
       
       fhtxv->Fill(track->Xv());
       fhtyv->Fill(track->Yv());
       fhtzv->Fill(track->Zv());
       
        for (Int_t j = 0; j < trdtracks; j++){      // trd track matching
            //cout << "sup " << j << endl;
            AliAODTrdTrack* trdtrack = static_cast<AliAODTrdTrack*>(fAOD->GetTrdTrack(j));
            if(!trdtrack) continue;
        
            // match the tracks
            AliVTrack* match = trdtrack->GetTrackMatch();  // work???
            if (!match) continue;
            
            if (!hasmatch){
                hasmatch = kTRUE;
                fhevent->Fill(2);
            }
            if ((TMath::Abs(trdtrack->Pt()) >= 2.) && !haspt){              // e_pt > 2
                haspt = kTRUE;
                fhevent->Fill(3);
            }
            if (trdtrack->GetPID() >= 164 && !haspid){     // this number was gotten from  ~PWG/TRD/AliTRDTriggerAnalysis.cxx
                haspid = kTRUE;
                fhevent->Fill(4);
            }
            if ((TMath::Abs(trdtrack->Pt()) >= 2.) && (trdtrack->GetPID() >= 164) && !haspitpt){
                haspitpt = kTRUE;
                fhevent->Fill(5);
            }
            
            if (clss.Contains("HQU")) {
                if (!hashqu){
                    hashqu = kTRUE;
                    fhevent->Fill(6);
                }
                if (!hqupt && (TMath::Abs(trdtrack->Pt()) >= 2.)){
                    hqupt = kTRUE;
                    fhevent->Fill(7);
                }
                if (!hqupid && trdtrack->GetPID() >= 164){
                    hqupid = kTRUE;
                    fhevent->Fill(8);
                }
                if (!hquptpid && (TMath::Abs(trdtrack->Pt()) >= 2. && trdtrack->GetPID() >= 164)){
                    hquptpid = kTRUE;
                    fhevent->Fill(9);
                }
            }
            
            fhfA->Fill(trdtrack->GetA());
            
        }
   }
   
   //TString clss = fAOD->GetFiredTriggerClasses();
    if (clss.Contains("HQU")) { //|| clss.Contains("HJE") || clss.Contains("HSE")) { 
        //cout << "have trigger hqu or the hje or hse triggers " << clss << endl; 
        Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
            for(Int_t i(0); i < iTracks; i++) {                 // loop ove rall these tracks
            AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
            if(!track || !track->TestFilterBit(1)) continue;                            // if we failed, skip this track
            
            if (track->Pt() < 2 && track->Pt() > 1)
                fhm1pt2->Fill(track->M());
            fhNpttr->Fill(track->Pt());
            
            
            //cout << track->GetA() << "" << endl;
        }                                 
    }
    else{
        Int_t iTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
        for(Int_t i(0); i < iTracks; i++) {                 // loop ove rall these tracks
            AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
            if(!track || !track->TestFilterBit(1)) continue;                            // if we failed, skip this track
            fhNptun->Fill(track->Pt());
        }
    }                                               // continue until all the tracks are processed

    PostData(1, fOutputList);                           // stream the results the analysis of this event to
                                                        // the output manager which will take care of writing
                                                        // it to a file
}




//________________________________________________________________________________
void AliAnalysisTRDEfficiency::Photons(AliAODEvent* fAOD)
{
    Int_t count = 0;
    Int_t ecount= 0;
    Int_t v0tracks  = fAOD->GetNumberOfV0s();
    
    //AliV0ReaderV1::ProcessESDV0s();

    
    for (Int_t i = 0; i < v0tracks; i++){  // all photons?
       
       AliAODv0* v0 = fAOD->GetV0(i);// = fAOD->GetVZERO(i);
       if (!v0) return;
       
       AliAODTrack* d0 = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(0));  // positive daughter
       AliAODTrack* d1 = (AliAODTrack *) (v0->GetSecondaryVtx()->GetDaughter(1));  // negative daughter
       
       //cout << d0->PID() << endl;
       //cout << v0->GetNegID() << "  " << v0->GetPosID() << endl;
       if (d0->PID() == 0 && d1->PID() == 0) ecount++;
       else count++;
       
       Double_t Vertptpos = (v0->DecayVertexV0X())*(v0->DecayVertexV0X()) + (v0->DecayVertexV0Y())*(v0->DecayVertexV0Y());
       Vertptpos = sqrt(Vertptpos);
       //cout << "legit v0 " << Vertptpos << endl;
       
       
       fhR->Fill(v0->RadiusV0());
       fhRpt->Fill(v0->Pt2V0(), v0->RadiusV0());
       //fhMv->Fill(v0->GetEffMass());
       fhvpost->Fill(Vertptpos);
       
       // check if has hqu trigger (would be more efficient if i did this outside the loop)
       // make a turn on curve
   }
   //cout << "electron counts " << ecount << " counts " << count << endl;
   
   TString clss = fAOD->GetFiredTriggerClasses();
   //cout << "trd trigger " << fAOD->GetTrdMask() << endl;
   if (clss.Contains("HQU")) { //|| clss.Contains("HJE") || clss.Contains("HSE")) { 
      for (Int_t i = 0; i < v0tracks; i++){
          
         AliAODv0* v0 = fAOD->GetV0(i);// = fAOD->GetVZERO(i);
         if (!v0) return;
         
         Double_t Vertptpos = (v0->DecayVertexV0X())*(v0->DecayVertexV0X()) + (v0->DecayVertexV0Y())*(v0->DecayVertexV0Y());
         Vertptpos = sqrt(Vertptpos);
         //cout << "legit v0 " << Vertptpos << endl;
         
         fhRhqu->Fill(v0->RadiusV0());
         fhRpthqu->Fill(v0->Pt2V0(), v0->RadiusV0());
         //fhMhqu->Fill(v0->M());
         fhvposthqu->Fill(Vertptpos);
         
      }
   }
   
   
}


//_____________________________________________________________________________
void AliAnalysisTRDEfficiency::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
