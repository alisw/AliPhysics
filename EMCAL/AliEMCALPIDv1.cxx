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
//#include "TROOT.h"
#include "TTree.h"
#include "TBenchmark.h"
#include "TSystem.h"
  
  // --- Standard library ---
  
  // --- AliRoot header files ---
#include "AliGenerator.h"
#include "AliEMCALPIDv1.h"
#include "AliEMCALRecParticle.h"
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
AliEMCALPIDv1::AliEMCALPIDv1(const AliEMCALPIDv1 & pid ):AliEMCALPID(pid)
{ 
  // ctor
  InitParameters() ; 
  Init() ;
  
}

//____________________________________________________________________________
AliEMCALPIDv1::AliEMCALPIDv1(const TString alirunFileName, const TString eventFolderName):AliEMCALPID(alirunFileName, eventFolderName)
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
}

//____________________________________________________________________________
const TString AliEMCALPIDv1::BranchName() const 
{  

  return GetName() ;
}
 
//____________________________________________________________________________
Float_t  AliEMCALPIDv1::GetCalibratedEnergy(Float_t e) const
{
//      It calibrates Energy depending on the recpoint energy.
//      The energy of the reconstructed cluster is corrected with 
//      the formula A + B* E  + C* E^2, whose parameters where obtained 
//      through the study of the reconstructed energy distribution of 
//      monoenergetic photons.
 
  //Float_t p[]={0.,0.,0.};
  //for (Int_t i=0; i<3; i++) p[i] = GetParameterCalibration(i);
  Float_t enerec = e ; // p[0] +  p[1]*e + p[2]*e*e;
  return enerec ;

}
//____________________________________________________________________________
TVector3 AliEMCALPIDv1::GetMomentumDirection(AliEMCALRecPoint * emc)const 
{ 
  // Calculates the momentum direction:
  // direction is given by IP and this RecPoint
  

  TVector3 dir(0,0,0) ; 
  TVector3 emcglobalpos ;
  // TMatrix  dummy ;
  
  emc->GetGlobalPosition(emcglobalpos) ;
  

  dir = emcglobalpos ;  
  // dir.SetMag(1.) ; Removed to avoid warings !!!!!!!!!!!!!! TO BE REVISED

  //account correction to the position of IP
  Float_t xo,yo,zo ; //Coordinates of the origin
  gAlice->Generator()->GetOrigin(xo,yo,zo) ;
  TVector3 origin(xo,yo,zo);
  dir = dir - origin ;

  return dir ;  
}

//____________________________________________________________________________
void AliEMCALPIDv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  // Add the PID task to the list of EMCAL tasks


  AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName.Data()) ; 

  if ( !gime->PID() ) 
    gime->PostPID(this) ;
}

//____________________________________________________________________________
void AliEMCALPIDv1::InitParameters()
{
  // Initialize the parameters
  fRecParticlesInRun = 0 ; 
  fNEvent            = 0 ;            
  fRecParticlesInRun = 0 ;
}

//____________________________________________________________________________

void  AliEMCALPIDv1::Exec(Option_t * option) 
{
  //Steering method

  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALPID");
  
  if(strstr(option,"print")) {
    Print("") ; 
    return ; 
  }
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  
  if (fLastEvent == -1) 
    fLastEvent = gime->MaxEvent() - 1 ;
  else 
    fLastEvent = TMath::Min(fLastEvent,gime->MaxEvent());
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;

  Int_t ievent ;

  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    gime->Event(ievent,"R") ;
    MakeRecParticles() ;
    WriteRecParticles();
    if(strstr(option,"deb"))
      PrintRecParticles(option) ;
    //increment the total number of rec particles per run 
    fRecParticlesInRun += gime->RecParticles()->GetEntriesFast() ;  
  }
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALPID");
    printf("Exec: took %f seconds for PID %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALPID"),  
	 gBenchmark->GetCpuTime("EMCALPID")/nEvents) ;
  } 

  Unload();
}

//____________________________________________________________________________
void  AliEMCALPIDv1::MakeRecParticles(){

  // Makes a RecParticle out of a TrackSegment
  
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  TObjArray * aECARecPoints = gime->ECARecPoints() ; 
  if ( !aECARecPoints ) {
    Fatal("MakeRecParticles", "RecPoints or TrackSegments not found !") ;  
  }
  TClonesArray * recParticles  = gime->RecParticles() ; 
  recParticles->Clear();

  TIter next(aECARecPoints) ; 
  AliEMCALRecPoint * eca ; 
  Int_t index = 0 ; 
  AliEMCALRecParticle * rp ; 
  while ( (eca = (AliEMCALRecPoint *)next()) ) {
    
    new( (*recParticles)[index] ) AliEMCALRecParticle() ;
    rp = (AliEMCALRecParticle *)recParticles->At(index) ; 
    rp->SetRecPoint(index) ;
    rp->SetIndexInList(index) ;
    	
    // Now set type (reconstructed) of the particle

    // Choose the cluster energy range
    
    Float_t  lambda[2] ;
    eca->GetElipsAxis(lambda) ;
    
    if((lambda[0]>0.01) && (lambda[1]>0.01)){
      // Looking PCA. Define and calculate the data (X),
      // introduce in the function X2P that gives the components (P).  

      Float_t  spher = 0. ;
      Float_t  emaxdtotal = 0. ; 
      
      if((lambda[0]+lambda[1])!=0) 
	spher=fabs(lambda[0]-lambda[1])/(lambda[0]+lambda[1]); 
      
      emaxdtotal=eca->GetMaximalEnergy()/eca->GetEnergy(); 
    }
    
    //    Float_t time = eca->GetTime() ;
      
    //Set momentum, energy and other parameters 
    Float_t  encal = GetCalibratedEnergy(eca->GetEnergy());
    TVector3 dir   = GetMomentumDirection(eca) ; 
    // dir.SetMag(encal) ;Removed to avoid warings !!!!!!!!!!!!!! TO BE REVISED
    rp->SetMomentum(dir.X(),dir.Y(),dir.Z(),encal) ;
    rp->SetCalcMass(0);
    rp->Name(); //If photon sets the particle pdg name to gamma
    rp->SetProductionVertex(0,0,0,0);
    rp->SetFirstMother(-1);
    rp->SetLastMother(-1);
    rp->SetFirstDaughter(-1);
    rp->SetLastDaughter(-1);
    rp->SetPolarisation(0,0,0);
    //Set the position in global coordinate system from the RecPoint
    //AliEMCALGeometry * geom = gime->EMCALGeometry() ; 
    //AliEMCALTowerRecPoint  * erp = gime->ECARecPoint(rp->GetEMCALRPIndex()) ; 
    TVector3 pos ; 
    //geom->GetGlobal(erp, pos) ; !!!!!!!!!! to check 
    rp->SetPos(pos);


    index++ ; 
  }
  
}

//____________________________________________________________________________
void  AliEMCALPIDv1:: Print(Option_t * /*option*/) const
{
  // Print the parameters used for the particle type identification

    printf("Print: =============== AliEMCALPID1 ================") ;
    printf("Making PID\n");
    printf("    Pricipal analysis file from 0.5 to 100 %s\n", fFileName.Data() ) ; 
    printf("    Name of parameters file     %s\n", fFileNamePar.Data() )  ;
}

//____________________________________________________________________________
void  AliEMCALPIDv1::Print() const
{
  // Print the parameters used for the particle type identification

    Info("Print", "=============== AliEMCALPIDv1 ================") ;
}

//____________________________________________________________________________
void AliEMCALPIDv1::PrintRecParticles(Option_t * option)
{
  // Print table of reconstructed particles

  AliEMCALGetter *gime = AliEMCALGetter::Instance() ; 

  TClonesArray * recParticles = gime->RecParticles() ; 

  printf("\nevent %i", gAlice->GetEvNumber()); 
  printf("       found %i", recParticles->GetEntriesFast()); 
  printf(" RecParticles\n"); 

  if(strstr(option,"all")) {  // printing found TS
    printf("\n  PARTICLE         Index    \n"); 
    
    Int_t index ;
    for (index = 0 ; index < recParticles->GetEntries() ; index++) {
      AliEMCALRecParticle * rp = (AliEMCALRecParticle * ) recParticles->At(index) ;       
      printf("\n");
      printf(rp->Name().Data());  
      printf(" %i", rp->GetIndexInList());  
      printf(" %i", rp->GetType());
    }
  }
}

//____________________________________________________________________________
void AliEMCALPIDv1::Unload() 
{
  // Unloads RecPoints, TrackSegments and RecParticles from the folder 
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ;  
  gime->EmcalLoader()->UnloadRecPoints() ;
  gime->EmcalLoader()->UnloadTracks() ;
  gime->EmcalLoader()->UnloadRecParticles() ;
}

//____________________________________________________________________________
void  AliEMCALPIDv1::WriteRecParticles()
{
  // Write RecParticles array to file
  AliEMCALGetter *gime = AliEMCALGetter::Instance() ; 

  TClonesArray * recParticles = gime->RecParticles() ; 
  recParticles->Expand(recParticles->GetEntriesFast() ) ;
  TTree * treeP =  gime->TreeP() ;

  
  
  //First rp
  Int_t bufferSize = 32000 ;    
  TBranch * rpBranch = treeP->Branch("EMCALRP",&recParticles,bufferSize);
  rpBranch->SetTitle(BranchName());

  rpBranch->Fill() ;
 
  gime->WriteRecParticles("OVERWRITE");
  gime->WritePID("OVERWRITE");

}

