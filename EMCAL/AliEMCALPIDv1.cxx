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

/* $Id$ */

//_________________________________________________________________________
// Implementation version v1 of the EMCAL particle identifier 

//*-- Author: Yves Schutz (SUBATECH) 
//
// --- ROOT system ---
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TSystem.h"
#include "TBenchmark.h"

#include "TSystem.h"

// --- Standard library ---

#include <Riostream.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliGenerator.h"
#include "AliEMCAL.h"
#include "AliEMCALPIDv1.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALTrackSegment.h"
#include "AliEMCALTrackSegmentMakerv1.h"
#include "AliEMCALRecParticle.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGetter.h"

ClassImp( AliEMCALPIDv1) 

//____________________________________________________________________________
AliEMCALPIDv1::AliEMCALPIDv1():AliEMCALPID()
{ 
  // default ctor
 
  InitParameters() ; 
  fDefaultInit = kTRUE ; 

}

//____________________________________________________________________________
AliEMCALPIDv1::AliEMCALPIDv1(const char * headerFile,const char * name, const Bool_t toSplit)
:AliEMCALPID(headerFile, name,toSplit)
{ 
  //ctor with the indication on where to look for the track segments
 
  InitParameters() ; 

  Init() ;
  fDefaultInit = kFALSE ; 

}

//____________________________________________________________________________
AliEMCALPIDv1::~AliEMCALPIDv1()
{ 
  // dtor

  if (!fDefaultInit) {  
    fSplitFile = 0 ; 
  }
}

//____________________________________________________________________________
const TString AliEMCALPIDv1::BranchName() const 
{  
  TString branchName(GetName() ) ;
  branchName.Remove(branchName.Index(Version())-1) ;
  return branchName ;
}
 
//____________________________________________________________________________
void AliEMCALPIDv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  // Add the PID task to the list of EMCAL tasks

  if ( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;

  TString branchname(GetName()) ;
  branchname.Remove(branchname.Index(Version())-1) ;    
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance(GetTitle(),branchname.Data(),fToSplit ) ; 

  //  gime->SetRecParticlesTitle(BranchName()) ;
  if ( gime == 0 ) {
    Error("Init", "Could not obtain the Getter object !" ) ;  
    return ;
  } 

  fSplitFile = 0 ;
  if(fToSplit){
    //First - extract full path if necessary
    TString fileName(GetTitle()) ;
    Ssiz_t islash = fileName.Last('/') ;
    if(islash<fileName.Length())
      fileName.Remove(islash+1,fileName.Length()) ;
    else
      fileName="" ;
    fileName+="EMCAL.RecData." ;
    if((strcmp(branchname.Data(),"Default")!=0)&&(strcmp(branchname.Data(),"")!=0)){
      fileName+=branchname.Data() ;
      fileName+="." ;
    }
    fileName+="root" ;
    fSplitFile = static_cast<TFile*>(gROOT->GetFile(fileName.Data()));   
    if(!fSplitFile)
      fSplitFile =  TFile::Open(fileName.Data(),"update") ;
  }
  
  gime->PostPID(this) ;
  gime->PostRecParticles(branchname) ; 
  
}

//____________________________________________________________________________
void AliEMCALPIDv1::InitParameters()
{

  fRecParticlesInRun = 0 ; 
  fNEvent            = 0 ;            
  fRecParticlesInRun = 0 ;
  TString pidName( GetName()) ;
  if (pidName.IsNull() ) 
    pidName = "Default" ; 
  pidName.Append(":") ; 
  pidName.Append(Version()) ; 
  SetName(pidName) ;
  fPi0Analysis = kFALSE ;
}

//____________________________________________________________________________

void  AliEMCALPIDv1::Exec(Option_t * option) 
{
  //Steering method
  
  if( strcmp(GetName(), "")== 0 ) 
    Init() ;
  
  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALPID");
  
  if(strstr(option,"print")) {
    Print("") ; 
    return ; 
  }
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  if(gime->BranchExists("RecParticles") )
    return ;
  Int_t nevents = gime->MaxEvent() ;       //(Int_t) gAlice->TreeE()->GetEntries() ;
  Int_t ievent ;


  for(ievent = 0; ievent < nevents; ievent++){
    gime->Event(ievent,"R") ;
 
    MakeRecParticles() ;
    
    WriteRecParticles(ievent);
    
    if(strstr(option,"deb"))
      PrintRecParticles(option) ;

    //increment the total number of rec particles per run 
    fRecParticlesInRun += gime->RecParticles(BranchName())->GetEntriesFast() ; 

  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALPID");
    Info("Exec", "took %f seconds for PID %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALPID"),  
	 gBenchmark->GetCpuTime("EMCALPID")/nevents) ;
  } 
}

//____________________________________________________________________________
void  AliEMCALPIDv1::MakeRecParticles(){

  // Makes a RecParticle out of a TrackSegment
  
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  TObjArray * aECRecPoints = gime->ECALRecPoints() ; 
  TObjArray * aPRRecPoints = gime->PRERecPoints() ; 
  TObjArray * aHCRecPoints = gime->HCALRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments() ; 
  if ( !aECRecPoints || !aPRRecPoints || !aHCRecPoints || !trackSegments ) {
    Fatal("MakeRecParticles", "RecPoints or TrackSegments not found !") ;  
  }
  TClonesArray * recParticles  = gime->RecParticles() ; 
  recParticles->Clear();

  TIter next(trackSegments) ; 
  AliEMCALTrackSegment * ts ; 
  Int_t index = 0 ; 
  AliEMCALRecParticle * rp ; 
  while ( (ts = (AliEMCALTrackSegment *)next()) ) {
    
    new( (*recParticles)[index] ) AliEMCALRecParticle() ;
    rp = (AliEMCALRecParticle *)recParticles->At(index) ; 
    rp->SetTrackSegment(index) ;
    rp->SetIndexInList(index) ;
    	
    AliEMCALTowerRecPoint * ecal = 0 ;
    if(ts->GetECIndex()>=0)
      ecal = dynamic_cast<AliEMCALTowerRecPoint *>(aECRecPoints->At(ts->GetECIndex())) ;
    
    AliEMCALTowerRecPoint    * pre = 0 ;
    if(ts->GetPREIndex()>=0)
      pre = dynamic_cast<AliEMCALTowerRecPoint *>(aPRRecPoints->At(ts->GetPREIndex())) ;
    
    AliEMCALTowerRecPoint    * hcal = 0 ;
    if(ts->GetHCIndex()>=0)
      hcal = dynamic_cast<AliEMCALTowerRecPoint *>(aHCRecPoints->At(ts->GetHCIndex())) ;

    // Now set type (reconstructed) of the particle

    // Choose the cluster energy range
    
    if (!ecal) {
      Fatal("MakeRecParticles", "-> emcal(%d) = %d", ts->GetECIndex(), ecal ) ;
    }

    Float_t    e = ecal->GetEnergy() ;   
    
    Float_t  lambda[2] ;
    ecal->GetElipsAxis(lambda) ;
    
    if((lambda[0]>0.01) && (lambda[1]>0.01)){
      // Looking PCA. Define and calculate the data (X),
      // introduce in the function X2P that gives the components (P).  

      Float_t  Spher = 0. ;
      Float_t  Emaxdtotal = 0. ; 
      
      if((lambda[0]+lambda[1])!=0) 
	Spher=fabs(lambda[0]-lambda[1])/(lambda[0]+lambda[1]); 
      
      Emaxdtotal=ecal->GetMaximalEnergy()/ecal->GetEnergy(); 
    }
    
    //    Float_t time = ecal->GetTime() ;
      
    //Set momentum, energy and other parameters 
    Float_t  encal = e;
    TVector3 dir(0., 0., 0.) ; 
    dir.SetMag(encal) ;
    rp->SetMomentum(dir.X(),dir.Y(),dir.Z(),encal) ;
    rp->SetCalcMass(0);
    rp->Name(); //If photon sets the particle pdg name to gamma
    rp->SetProductionVertex(0,0,0,0);
    rp->SetFirstMother(-1);
    rp->SetLastMother(-1);
    rp->SetFirstDaughter(-1);
    rp->SetLastDaughter(-1);
    rp->SetPolarisation(0,0,0);
    index++ ; 
  }
  
}

//____________________________________________________________________________
void  AliEMCALPIDv1:: Print()
{
  // Print the parameters used for the particle type identification

  TString message ; 
    message  = "\n=============== AliEMCALPID1 ================\n" ;
    message += "Making PID\n";
    message += "    Pricipal analysis file from 0.5 to 100 %s\n" ; 
    message += "    Name of parameters file     %s\n" ;
    message += "    Matrix of Parameters: 9x4\n" ;
    message += "        RCPV 2x3 rows x and z, columns function cut parameters\n" ;
    message += "        TOF  1x3 [High Eff-Low Pur,Medium Eff-Pur, Low Eff-High Pur]\n" ;
    message += "        PCA  5x4 [5 ellipse parametres and 4 parametres to calculate them: A/Sqrt(E) + B* E + C * E^2 + D]\n" ;
    message += "        Energy Calibration  1x3 [3 parametres to calibrate energy: A + B* E + C * E^2]\n" ;
    Info("Print", message.Data(), fFileName.Data(), fFileNamePar.Data() ) ; 
}

//____________________________________________________________________________
void  AliEMCALPIDv1::WriteRecParticles(Int_t event)
{
 
  AliEMCALGetter *gime = AliEMCALGetter::GetInstance() ; 

  TClonesArray * recParticles = gime->RecParticles() ; 
  recParticles->Expand(recParticles->GetEntriesFast() ) ;
  TTree * treeR ;

  if(fToSplit){
    if(!fSplitFile)
      return ;
    fSplitFile->cd() ;
    char name[10] ;
    sprintf(name,"%s%d", "TreeR",event) ;
    treeR = dynamic_cast<TTree*>(fSplitFile->Get(name)); 
  }
  else{
    treeR = gAlice->TreeR();
  }
  
  if(!treeR){
    gAlice->MakeTree("R", fSplitFile);
    treeR = gAlice->TreeR() ;
  }
  
  //First rp
  Int_t bufferSize = 32000 ;    
  TBranch * rpBranch = treeR->Branch("EMCALRP",&recParticles,bufferSize);
  rpBranch->SetTitle(BranchName());

  
  //second, pid
  Int_t splitlevel = 0 ; 
  AliEMCALPIDv1 * pid = this ;
  TBranch * pidBranch = treeR->Branch("AliEMCALPID","AliEMCALPIDv1",&pid,bufferSize,splitlevel);
  pidBranch->SetTitle(BranchName());
  
  rpBranch->Fill() ;
  pidBranch->Fill() ; 
  
  treeR->AutoSave() ; //Write(0,kOverwrite) ;  
  if(gAlice->TreeR()!=treeR){
    treeR->Delete();
  }
}

//____________________________________________________________________________
void AliEMCALPIDv1::PrintRecParticles(Option_t * option)
{
  // Print table of reconstructed particles

  AliEMCALGetter *gime = AliEMCALGetter::GetInstance() ; 

  TClonesArray * recParticles = gime->RecParticles(BranchName()) ; 

  TString message ; 
  message  = "\nevent " ;
  message += gAlice->GetEvNumber() ; 
  message += "       found " ; 
  message += recParticles->GetEntriesFast(); 
  message += " RecParticles\n" ; 

  if(strstr(option,"all")) {  // printing found TS
    message += "\n  PARTICLE         Index    \n" ; 
    
    Int_t index ;
    for (index = 0 ; index < recParticles->GetEntries() ; index++) {
      AliEMCALRecParticle * rp = (AliEMCALRecParticle * ) recParticles->At(index) ;       
      message += "\n" ;
      message += rp->Name().Data() ;  
      message += " " ;
      message += rp->GetIndexInList() ;  
      message += " " ;
      message += rp->GetType()  ;
    }
  }
  Info("Print", message.Data() ) ; 
}



