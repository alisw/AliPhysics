/////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
/////////////////////////////////////////////////////////////
//
// Description: AliRoot macro to make AliFlowEvents from KineTree (new way) 
//
/////////////////////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH2.h"
#include "TFile.h"
#include "TObjArray"
#include "TStopwatch.h"
//#include "AliESD.h"

using namespace std; //required for resolving the 'cout' symbol

//////////////////////////////////////////////////////////////////////////////////////////////////////

const char* aDataDir = "./" ; 
Int_t aRuns = -1 ; Int_t offset = -1 ;

AliRunLoader*  rl = 0 ;
AliRun*    gAlice = 0 ;
TTree*      treeK = 0 ;
AliStack*  kStack = 0 ;
TChain*  esdChain = 0 ;

TTree* kOpen(int evtN=0) ;
Int_t  kInit(const char* dir="./") ;

Double_t Norm(Double_t nu[3]) ; 
Double_t Phi(Double_t nu[3]) ;
Double_t Pt(Double_t nu[3]) ;
Double_t Eta(Double_t nu[3]) ;

//////////////////////////////////////////////////////////////////////////////////////////////////////

int efficiency2s()
{
 cout << " . Here the efficiency2 things ... " << endl ;
 cout << endl ;

 TStopwatch timer ;
 timer.Start() ;

// global settings, counters and pointers ...

 int binPt = 100 ;                // bins in pT 
 int binEta = 42 ;                // bins in eta 
 float maxPt  = 5. ;              // max pT
 float maxEta = 2.1 ;             // max eta
 float pTbinWidth = maxPt/binPt ; // bin width (to quantify the smearing of pT)

 const Int_t kCode[6] = {11,13,211,321,2212,10010020} ;  // stable particles
  
 int countsim = 0 ;               // number of simulations (esd files or folders)
 int countev = 0 ;                // total number of looped events
 int tot = 0 ;                    // total number of looped esd tracks
 int totPerf = 0 ;                // total number of perfectly reconstructed tracks (in cut + mc primary + same bin)

 AliESD*       esd = 0 ;          // Event Summary Data
 AliESDtrack*  esdTrack = 0 ;     // reconstructed track
 TParticle*    particle = 0 ;     // TParticle (kineTree)
 TParticlePDG* particlePDG = 0 ;  // TParticlePDG (type, charge, etc.)
 //TString       esdFileName("AliESDs.root") ;
 //TString       gAliceName("galice.root") ;
  
// define the structs ...

 struct structSameBin 
 { 
  TH2F* hYield2Desd ;
  TH2F* hYield2Dmc ;
  
  TH2F* hPtPt ;
  TH2F* hEtaEta ;
  TH2F* hPhiPhi ;
 } ;
 struct structMcPrim { struct structSameBin histSameBin[2] ; } ;
 struct structEsdIn  { struct structMcPrim histMcPrim[2] ; } ;
 struct structEsdIn  histEsdIn[2] ;
 
// define the histograms ...

 TString histTitle ;
 for(int i = 0; i < 2; i++)
 {
  for(int j = 0; j < 2; j++) 
  {
   for(int k = 0; k < 2; k++)
   {
    TString histDetails ;
    if(i) { histDetails = "_inCut"    ; } else { histDetails = "_outCut"   ; }
    if(j) { histDetails += "_mcPrim"  ; } else { histDetails += "_mcSec"   ; }
    if(k) { histDetails += "_sameBin" ; } else { histDetails += "_diffBin" ; }
    
   // Yield pt - eta
    histTitle = "hYield2Desd" ; histTitle += histDetails ;
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Desd = new TH2F(histTitle.Data(), histTitle.Data(), binEta, -maxEta, maxEta, binPt, 0., maxPt);
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Desd->Sumw2();
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Desd->SetXTitle("#eta");
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Desd->SetYTitle("p_{T} (GeV/c)");

   // Yield pt - eta
    histTitle = "hYield2Dmc" ; histTitle += histDetails ;
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Dmc = new TH2F(histTitle.Data(), histTitle.Data(), binEta, -maxEta, maxEta, binPt, 0., maxPt);
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Dmc->Sumw2();
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Dmc->SetXTitle("#eta");
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Dmc->SetYTitle("p_{T} (GeV/c)");

   // comparison (esd vs kine)
    histTitle = "hPtPt" ; histTitle += histDetails ;
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPtPt = new TH2F(histTitle.Data(), histTitle.Data(), binPt, 0., maxPt, binPt, 0., maxPt);
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPtPt->Sumw2();
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPtPt->SetXTitle("p_{T} (GeV/c) kine");
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPtPt->SetYTitle("p_{T} (GeV/c) esd");

    histTitle = "hEtaEta" ; histTitle += histDetails ;
    histEsdIn[i].histMcPrim[j].histSameBin[k].hEtaEta = new TH2F(histTitle.Data(), histTitle.Data(), binEta, -maxEta, maxEta, binEta, -maxEta, maxEta);
    histEsdIn[i].histMcPrim[j].histSameBin[k].hEtaEta->Sumw2();
    histEsdIn[i].histMcPrim[j].histSameBin[k].hEtaEta->SetXTitle("#eta kine");
    histEsdIn[i].histMcPrim[j].histSameBin[k].hEtaEta->SetYTitle("#eta esd");

    histTitle = "hPhiPhi" ; histTitle += histDetails ;
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPhiPhi = new TH2F(histTitle.Data(), histTitle.Data(), 180, 0., 360., 180, 0., 360.);
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPhiPhi->Sumw2();
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPhiPhi->SetXTitle("#phi kine");
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPhiPhi->SetYTitle("#phi esd");

    // delete histTitle ; delete histDetails ;
   }
  }
 }

//  TH2F* pTpT = new TH2F("", "", 20, -2., 2., 30, 0., 3.);
//  ->Sumw2();
//  ->SetXTitle("Pt (GeV/c)");
//  ->SetYTitle("Pt (GeV/c)");

// loop (more simulations) ...

 TString execDir(gSystem->pwd());
 TSystemDirectory* baseDir = new TSystemDirectory(".", aDataDir);
 TList* dirList 	   = baseDir->GetListOfFiles();
 Int_t nDirs		   = dirList->GetEntries();
 gSystem->cd(execDir);
 
 for (Int_t iDir=0; iDir<nDirs; ++iDir)
 {
  TSystemFile* presentDir = (TSystemFile*)dirList->At(iDir);
  if(!presentDir || !presentDir->IsDirectory() || strcmp(presentDir->GetName(), ".") == 0 || strcmp(presentDir->GetName(), "..") == 0) 
  { continue ; }
  if(offset > 0)  { --offset ; continue ; }
  if(countsim == aRuns) { break ; }
 
  TString presentDirName(aDataDir); 
  presentDirName += presentDir->GetName();
  presentDirName += "/"; 

// initialise esdChain and kineTree ... 
  
  Int_t fNumberOfESDs      = eInit(presentDirName.Data()) ;
  Int_t fNumberOfKineTrees = kInit(presentDirName.Data()) ;

  Int_t fNumberOfEvents = -1 ;
  if(fNumberOfKineTrees == fNumberOfESDs)  { fNumberOfEvents = fNumberOfKineTrees ; }
  if(fNumberOfEvents < 0) 		   { return fNumberOfEvents ; }

  esdChain->SetBranchAddress("ESD",&esd) ;

  cout << " Directory :  " << presentDirName << " ;  n. of events : " << fNumberOfEvents << endl ;
  
// loop (events) ...

  for(int evtN=0;evtN<fNumberOfEvents;evtN++)   // use (1) for testing the folder loop
  {
   int countrk = 0 ;
   int nConstrainable = 0 ;
   int nPrimaries = 0 ;
   int nPerfect = 0 ;
   
   treeK = kOpen(evtN) ;
   treeK->SetBranchAddress("Particles",&particle) ;    // set the place where to put the TParticle from the tree (filled every time you call GetEntry() )
  
   esdChain->GetEntry(evtN) ;

   Int_t fNumberOfParticles = treeK->GetEntries() ;
   Int_t fNumberOfPrimaries = kStack->GetNprimary() ;   
   Int_t fNumberOfTracks = esd->GetNumberOfTracks() ;

   cout << " Event " << evtN << " has " << fNumberOfTracks << " tracks in the ESD, and " << fNumberOfParticles << " particles in the kinetree (" << fNumberOfPrimaries << " primaries) ." << endl ;
  
// loop (tracks) ...

   float ptTr, ptP ;
   float phiTr, phiP ;
   float etaTr, etaP ;
   Double_t cD[3] ; // for(Int_t iii=0;iii<3;iii++) { cD[iii] = 0 ; }				
   for(Int_t nTracks = 0; nTracks < fNumberOfTracks; nTracks++)  
   {
    esdTrack = esd->GetTrack(nTracks) ;
    Int_t label = esdTrack->GetLabel() ;

    Int_t idStack = TMath::Abs(label) ; 	      // (idStack) is just the label of the particle (refers to the stack)
    Int_t idTree = kStack->TreeKEntry(idStack) ;      // the stack knows the position (idTree) of the particle (idStack) in the tree
    treeK->GetEntry(idTree) ;			      // get the particle  (idTree) out of the tree
    
    // cout << "  Track " << nTracks << " has label = " << label << " and idTree = " << idTree << endl ;

    Int_t incut   = 0 ; if(esdTrack->GetPxPyPz(cD))	      { incut = 1 ; } 
    Int_t primary = 0 ; if(particle->IsPrimary())	      { primary = 1 ; }
    Int_t perfect = 0 ; if(TMath::Abs(ptTr-ptP) < pTbinWidth) { perfect = 1 ; }

    if(incut) { nConstrainable++ ; if(primary) { nPrimaries++ ; if(perfect) { nPerfect++ ; } } }

   // Constrained pt, eta, phi from the ESDtrack ...
    phiTr = (Float_t)Phi(cD) ; if(phiTr<0) { phiTr += 2*TMath::Pi() ; }
    ptTr  = (Float_t)Pt(cD) ;  if(ptTr<=0) { cout << " !!! pt = " << ptTr << endl ; continue ; }
    etaTr = (Float_t)Eta(cD) ;     
   // Particle kinematic ...
    phiP = particle->Phi() ;   if(phiP<0) { phiP += 2*TMath::Pi() ; }
    ptP = particle->Pt();
    etaP = particle->Eta();
   // Pid ...
    // particlePDG = particle->GetPDG() ;
    // ...

    histEsdIn[incut].histMcPrim[primary].histSameBin[perfect].hYield2Desd->Fill(etaTr,ptTr) ; 
    histEsdIn[incut].histMcPrim[primary].histSameBin[perfect].hYield2Dmc->Fill(etaP,ptP) ; 
    histEsdIn[incut].histMcPrim[primary].histSameBin[perfect].hPtPt->Fill(ptP,ptTr) ;
    histEsdIn[incut].histMcPrim[primary].histSameBin[perfect].hEtaEta->Fill(etaP,etaTr) ;
    histEsdIn[incut].histMcPrim[primary].histSameBin[perfect].hPhiPhi->Fill((180./TMath::Pi())*phiP,(180./TMath::Pi())*phiTr) ;
    
    // cout << " particle n. " << i << " mothers: " << particle->GetMother(0) << " , " << particle->GetMother(1) <<  endl ;	
    //cout << "  Track " << nTracks << " has label = " << label << " and idTree = " << idTree << " is primary && constrainable " << endl ;
    countrk++ ;
   }

   cout << " Event looped " << evtN << " - " << countrk << " tracks found: " << nConstrainable << " constrainable, of which " << nPrimaries << " primaries, and " << nPerfect << " perfectly reconstructed . " << endl ; cout << endl ; 
   tot += countrk ;
   totPerf += nPerfect ;	     
   countev++ ;
  }	
  esdChain->Reset() ; // delete esdChain ;
  
  //rl->UnloadgAlice() ; 
  //rl->UnloadHeader() ; 
  //rl->UnloadKinematics() ; 
  rl->UnloadAll() ;
  rl->RemoveEventFolder() ;
  rl->Clear() ; // delete rl ;
  
  delete presentDir ;
  countsim++ ;
 }
 
 cout <<  endl ;
 cout << " Finished ... " << endl ;
 cout << "  nEvts     :  " << countev << " in " << countsim << " simulations " << endl ;  	     
 cout << "  nParticles:  " << tot << " (" << totPerf << " perfect) " << endl ;   
 cout << " . " << endl ; 

//  // plots ...
// 
//  TCanvas* vs = new TCanvas("vs","vs",900,400) ;
//  vs->Divide(3,1) ;
//  vs->cd(1) ; ...->Draw("COLZ") ;
//  vs->cd(2) ; ...->Draw("COLZ") ;
//  vs->cd(3) ; ...->Draw("COLZ") ;

// break ;

// output file ...

 TString output("efficiency2s.root") ;
 TFile * file = new TFile(output.Data(),"RECREATE") ;
 file->cd() ; 
 for(int i = 0; i < 2; i++)
 {
  for(int j = 0; j < 2; j++) 
  {
   for(int k = 0; k < 2; k++)
   {
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Desd->Write() ;
    histEsdIn[i].histMcPrim[j].histSameBin[k].hYield2Dmc->Write() ; 
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPtPt->Write() ;
    histEsdIn[i].histMcPrim[j].histSameBin[k].hEtaEta->Write() ;
    histEsdIn[i].histMcPrim[j].histSameBin[k].hPhiPhi->Write() ;
   }
  }
 }
 file->Close() ;
 cout << " file " << output.Data() << " saved ... " << endl ;

 timer.Stop() ;
 cout << endl ;
 timer.Print() ;
 cout << " . here it was (efficiency2s) ... " << endl ;  //juice!
 cout << endl ;

 cout << endl ; cout << " Memory Check (from Paul)" << endl ; 
 gObjectTable->Print();
 cout << endl ; cout << endl ;

 return 2 ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t eInit(const char* dir)
{
 TString fileName(dir) ; 
 fileName += "AliESDs.root" ;
 
 Long_t *id, *size, *flags, *modtime ;
 if(gSystem->GetPathInfo(fileName.Data(),id,size,flags,modtime)) 
 { 
  cout << " File : " << fileName << " does NOT exist ! - Skipping ... " << endl ; 
  return 0 ; 
 }
 esdChain = new TChain("esdTree") ;
 esdChain->Add(fileName.Data()) ;
 Int_t nEvents = esdChain->GetEntries() ;

 return nEvents ;
} 

//////////////////////////////////////////////////////////////////////////////////////////////////////

Int_t kInit(const char* dir)
{
 TString fileName(dir) ; 
 fileName += "galice.root" ;
 
 Long_t *id, *size, *flags, *modtime ;
 if(gSystem->GetPathInfo(fileName.Data(),id,size,flags,modtime)) 
 { 
  cout << " File : " << fileName << " does NOT exist ! - Skipping ... " << endl ; 
  return 0 ; 
 }
 
 rl = AliRunLoader::Open(fileName.Data(),"MyEvent","read");  // AliRunLoader* 
 rl->LoadgAlice();
 gAlice = rl->GetAliRun();  // AliRun*
 rl->LoadHeader();
 rl->LoadKinematics();
 Int_t nEvents = rl->GetNumberOfEvents() ;
 
 return nEvents ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

TTree* kOpen(int evtN)
{
 Int_t exitStatus = rl->GetEvent(evtN) ; if(exitStatus!=0) { return 0 ; }

 TTree* kTree = (TTree*)rl->TreeK();  // Particles TTree (KineTree)
 kStack = gAlice->Stack();	      // Particles Stack (use "Label()" to get the number in the stack)

 Int_t fNumberOfParticles = kTree->GetEntries() ;
 Int_t nPart = kStack->GetNtrack() ;
 Int_t nPrim = kStack->GetNprimary() ;
 if(fNumberOfParticles !=  nPart) { cout << " Something is wrong !!! " << fNumberOfParticles << " !=  " << nPart << " . " << endl ; }
 // else { cout << " Event n. " << evtN << "  contains :  " << fNumberOfParticles << "  particles in the TTree  ( =  " << nPart << "  in the stack, " << nPrim << " ) . " << endl ; }

 return kTree ; 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------
//*** USEFULL METHODS for a 3-array of double (~ TVector3) ***
//-----------------------------------------------------------------------
Double_t Norm(Double_t nu[3])
{ 
 // returns the norm of a double[3] 

 Double_t norm2 = nu[0]*nu[0] + nu[1]*nu[1] + nu[2]*nu[2] ;
 return TMath::Sqrt(norm2) ; 
}
//-----------------------------------------------------------------------
Double_t Phi(Double_t nu[3])
{
 // returns the azimuthal angle of a double[3] 

 if(nu[0]==0 && nu[1]==0) { return 0. ; }
 else 			  { return TMath::ATan2(nu[1],nu[0]) ; }
}
//-----------------------------------------------------------------------
Double_t Pt(Double_t nu[3])
{
 // returns the transvers momentum of a double[3] 

 Double_t trans = nu[0]*nu[0] + nu[1]*nu[1] ;
 return TMath::Sqrt(trans) ; 
}
//-----------------------------------------------------------------------
Double_t Eta(Double_t nu[3])
{
 // returns the PseudoRapidity of a double[3] 
 // if transvers momentum = 0 --> returns +/- 1.000

 Double_t m = Norm(nu) ;
 if(nu[0]!=0 || nu[1]!=0) { return 0.5*TMath::Log((m+nu[2])/(m-nu[2])) ; }
 else     	 	  { return TMath::Sign((Double_t)1000.,nu[2]) ; }
}
//-----------------------------------------------------------------------
