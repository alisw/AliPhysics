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
// Implementation version v0 of the PHOS particle identifier 
// Particle identification based on the 
//     - CPV information, 
//     - Preshower information (in MIXT or GPS2 geometries)
//     - shower width.
//
// CPV or Preshower clusters should be closer in PHOS plane than fCpvEmcDistance (in cm).
// This parameter can be set by method SetCpvtoEmcDistanceCut(Float_t cut)  
//
// One can set desirable ID method by the function SetIdentificationMethod(option).
// Presently the following options can be used together or separately :
//     - "disp": use dispersion cut on shower width 
//               (width can be set by method SetDispersionCut(Float_t cut)
//     - "ell" : use cut on the axis of the ellipse, drawn around shower 
//       (this cut can be changed by SetShowerProfileCut(char* formula), 
//        where formula - any function of two variables f(lambda[0],lambda[1]).
//        Shower is considered as EM if f() > 0 )
// One can visualize current cuts calling method PlotDispersionCuts().    
//
// use case:
//  root [0] AliPHOSPIDv0 * p1 = new AliPHOSPIDv0("galice.root")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//  root [1] p1->SetIdentificationMethod("disp ellipse")
//  root [2] p1->ExecuteTask()
//  root [3] AliPHOSPIDv0 * p2 = new AliPHOSPIDv0("galice1.root","ts1")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                // reading headers from file galice1.root and TrackSegments 
//                // with title "ts1"
//  root [4] p2->SetRecParticlesBranch("rp1")
//                // set file name for the branch RecParticles
//  root [5] p2->ExecuteTask("deb all time")
//                // available options
//                // "deb" - prints # of reconstructed particles
//                // "deb all" -  prints # and list of RecParticles
//                // "time" - prints benchmarking results
//                  
//*-- Author: Yves Schutz (SUBATECH)  & Gines Martinez (SUBATECH) & 
//            Dmitri Peressounko (SUBATECH & Kurchatov Institute)
//            Completely redesined by Dmitri Peressounko, March 2001

// --- ROOT system ---
#include "TTree.h"
#include "TF2.h"
#include "TFormula.h"
#include "TCanvas.h"

#include "TBenchmark.h"
// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliRun.h"
#include "AliGenerator.h"
#include "AliPHOSPIDv0.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSLoader.h"

ClassImp( AliPHOSPIDv0) 

//____________________________________________________________________________
AliPHOSPIDv0::AliPHOSPIDv0():AliPHOSPID()
{ 
  // default ctor
  fFormula           = 0 ;
  fDispersion        = 0. ; 
  fCpvEmcDistance    = 0 ; 
  fTimeGate          = 2.e-9 ;
  fEventFolderName    = "" ; 
  fTrackSegmentsTitle= "" ; 
  fRecPointsTitle    = "" ; 
  fRecParticlesTitle = "" ; 
  fIDOptions         = "dis time" ; 
  fRecParticlesInRun = 0 ;
  fClusterizer = 0;
  fTSMaker = 0;
}

//____________________________________________________________________________
AliPHOSPIDv0::AliPHOSPIDv0(const char * evFolderName,const char * name) : AliPHOSPID(evFolderName, name)
{ 
  //ctor with the indication on where to look for the track segments

  fFormula        = new TFormula("LambdaCuts","(x>1)*(x<2.5)*(y>0)*(y<x)") ;   
  fDispersion     = 2.0 ; 
  fCpvEmcDistance = 3.0 ;
  fTimeGate          = 2.e-9 ;
 
  fEventFolderName     = GetTitle() ; 
  fTrackSegmentsTitle = GetName();
  fRecPointsTitle     = GetName();
  fRecParticlesTitle  = GetName();
  fIDOptions          = "dis time" ;
    
  fRecParticlesInRun = 0 ; 

  Init() ;

}

//____________________________________________________________________________
AliPHOSPIDv0::~AliPHOSPIDv0()
{ 
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv0::GetDistance(AliPHOSEmcRecPoint * emc,AliPHOSRecPoint * cpv, Option_t *  Axis)const
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
 
  const AliPHOSGeometry * geom = AliPHOSLoader::GetPHOSGeometry() ; 
  TVector3 vecEmc ;
  TVector3 vecCpv ;
  
  emc->GetLocalPosition(vecEmc) ;
  cpv->GetLocalPosition(vecCpv) ; 
  if(emc->GetPHOSMod() == cpv->GetPHOSMod()){ 
    
    // Correct to difference in CPV and EMC position due to different distance to center.
    // we assume, that particle moves from center
    Float_t dCPV = geom->GetIPtoOuterCoverDistance();
    Float_t dEMC = geom->GetIPtoCrystalSurface() ;
    dEMC         = dEMC / dCPV ;
    vecCpv = dEMC * vecCpv  - vecEmc ; 
    if (Axis == "X") return vecCpv.X();
    if (Axis == "Y") return vecCpv.Y();
    if (Axis == "Z") return vecCpv.Z();
    if (Axis == "R") return vecCpv.Mag();
  } 
 
  return 100000000 ;
}

//____________________________________________________________________________
void  AliPHOSPIDv0::Exec(Option_t * option) 
{
  //Steering method
  
  if( strcmp(GetName(), "")== 0 ) 
    Init() ;
  
  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSPID");
  
  if(strstr(option,"print")) {
    Print() ; 
    return ; 
  }

  AliRunLoader* runget = AliRunLoader::GetRunLoader(GetTitle());
  if(runget == 0x0) 
   {
     AliError(Form("Can not find run getter in event folder \"%s\"",
		   GetTitle()));
     return;
   }
  
  AliPHOSLoader* gime = dynamic_cast<AliPHOSLoader*>(runget->GetLoader("PHOSLoader"));
  if ( gime == 0 ) 
   {
     AliError("Could not obtain the Loader object !"); 
     return ;
   } 

  if(gime->BranchExists("RecParticles") )
    return ;

  
  Int_t nevents = runget->GetNumberOfEvents() ;       //(Int_t) gAlice->TreeE()->GetEntries() ;

  Int_t ievent ;
  
  for(ievent = 0; ievent < nevents; ievent++){
    runget->GetEvent(ievent);
    AliInfo(Form("event %d %d %d", 
		 ievent, gime->EmcRecPoints(), 
		 gime->TrackSegments())) ;
    MakeRecParticles() ;
    
    WriteRecParticles();
    
    if(strstr(option,"deb"))
      PrintRecParticles(option) ;

    //increment the total number of rec particles per run 
    fRecParticlesInRun += gime->RecParticles()->GetEntriesFast() ; 

  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSPID");
    AliInfo(Form("took %f seconds for PID %f seconds per event", 
		 gBenchmark->GetCpuTime("PHOSPID"), 
		 gBenchmark->GetCpuTime("PHOSPID")/nevents)) ; 
  }
  
}
//____________________________________________________________________________
void AliPHOSPIDv0::Init()
{
  // Make all memory allocations that are not possible in default constructor
  // Add the PID task to the list of PHOS tasks
  
  AliRunLoader* runget = AliRunLoader::GetRunLoader(GetTitle());
  if(runget == 0x0) 
   {
     AliError(Form("Can not find run getter in event folder \"%s\"",
		   GetTitle()));
     return;
   }
  
  AliPHOSLoader* gime = dynamic_cast<AliPHOSLoader*>(runget->GetLoader("PHOSLoader"));
  if ( gime == 0 ) 
   {
     AliError("Could not obtain the Loader object !"); 
     return ;
   } 
   
  gime->PostPID(this);
  gime->LoadRecParticles("UPDATE");
  
}

//____________________________________________________________________________
void  AliPHOSPIDv0::MakeRecParticles()
{
  // Reconstructs the particles from the tracksegments

  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;

  AliRunLoader* runget = AliRunLoader::GetRunLoader(GetTitle());
  if(runget == 0x0) 
   {
     AliError(Form("Can not find run getter in event folder \"%s\"",
		   GetTitle()));
     return;
   }
  
  AliPHOSLoader* gime = dynamic_cast<AliPHOSLoader*>(runget->GetLoader("PHOSLoader"));
  if ( gime == 0 ) 
   {
     AliError("Could not obtain the Loader object !"); 
     return ;
   } 

  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments() ; 
  TClonesArray * recParticles  = gime->RecParticles() ; 

  recParticles->Clear();
  
  TIter next(trackSegments) ; 
  AliPHOSTrackSegment * ts ; 
  Int_t index = 0 ; 
  AliPHOSRecParticle * rp ; 
  
  Bool_t ellips = fIDOptions.Contains("ell",TString::kIgnoreCase ) ;
  Bool_t disp   = fIDOptions.Contains("dis",TString::kIgnoreCase ) ;
  Bool_t time   = fIDOptions.Contains("tim",TString::kIgnoreCase ) ;
  
  while ( (ts = (AliPHOSTrackSegment *)next()) ) {
    
    new( (*recParticles)[index] ) AliPHOSRecParticle() ;
    rp = (AliPHOSRecParticle *)recParticles->At(index) ; 
    rp->SetTrackSegment(index) ;
    rp->SetIndexInList(index) ;
    
    AliPHOSEmcRecPoint * emc = 0 ;
    if(ts->GetEmcIndex()>=0)
      emc = (AliPHOSEmcRecPoint *) emcRecPoints->At(ts->GetEmcIndex()) ;
    
    AliPHOSRecPoint    * cpv = 0 ;
    if(ts->GetCpvIndex()>=0)
      cpv = (AliPHOSRecPoint *)   cpvRecPoints->At(ts->GetCpvIndex()) ;
    
    //set momentum and energy first
    Float_t    e = emc->GetEnergy() ;
    TVector3 dir = GetMomentumDirection(emc,cpv) ; 
    dir.SetMag(e) ;

    rp->SetMomentum(dir.X(),dir.Y(),dir.Z(),e) ;
    rp->SetCalcMass(0);

    //now set type (reconstructed) of the particle    
    Int_t showerprofile = 0;  // 0 narrow and 1 wide
    
    if(ellips){
      Float_t lambda[2] ;
      emc->GetElipsAxis(lambda) ;
      if(fFormula->Eval(lambda[0],lambda[1]) <= 0 )
	showerprofile = 1 ;  // not narrow
    }
    
    if(disp)
      if(emc->GetDispersion() > fDispersion )
	showerprofile = 1 ;  // not narrow
    
    Int_t slow = 0 ;
    if(time)
      if(emc->GetTime() > fTimeGate )
	slow = 0 ; 
        
    // Looking at the CPV detector
    Int_t cpvdetector= 0 ;  //1 hit and 0 no hit     
    if(cpv)
      if(GetDistance(emc, cpv,  "R") < fCpvEmcDistance) 
	cpvdetector = 1 ;  
    
    Int_t type = showerprofile + 2 * slow  + 4 * cpvdetector ;
    rp->SetType(type) ; 
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
void  AliPHOSPIDv0:: Print() const
{
  // Print the parameters used for the particle type identification
  TString message ; 
  message  = "=============== AliPHOSPIDv0 ================\n" ;
  message += "Making PID\n" ;
  message += "    Headers file:               %s\n" ; 
  message += "    RecPoints branch title:     %s\n" ;
  message += "    TrackSegments Branch title: %s\n" ; 
  message += "    RecParticles Branch title   %s\n" ;  
  message += "with parameters:\n"  ;
  message += "    Maximal EMC - CPV  distance (cm) %f\n" ;
  AliInfo(Form( message.Data(),  
       GetTitle(), 
       fRecPointsTitle.Data(), 
       fTrackSegmentsTitle.Data(), 
       fRecParticlesTitle.Data(), 
       fCpvEmcDistance ));

  if(fIDOptions.Contains("dis",TString::kIgnoreCase ))
    AliInfo(Form("                    dispersion cut %f",  fDispersion )) ;
  if(fIDOptions.Contains("ell",TString::kIgnoreCase ))
    AliInfo(Form("             Eliptic cuts function: %s",  
		 fFormula->GetTitle() )) ;
  if(fIDOptions.Contains("tim",TString::kIgnoreCase ))
    AliInfo(Form("             Time Gate used: %f",  fTimeGate)) ;
}

//____________________________________________________________________________
void  AliPHOSPIDv0::SetShowerProfileCut(const char * formula)
{
  //set shape of the cut on the axis of ellipce, drown around shouer
  //shower considered "narrow" if Formula(lambda[0],lambda[1]) > 0.
  if(fFormula) 
    delete fFormula; 
  fFormula = new TFormula("Lambda Cut",formula) ;
}
//____________________________________________________________________________
void  AliPHOSPIDv0::WriteRecParticles()
{
  // Saves the reconstructed particles too a file
 
  AliRunLoader* runget = AliRunLoader::GetRunLoader(GetTitle());
  if(runget == 0x0) 
   {
     AliError(Form("Can not find run getter in event folder \"%s\"",
		   GetTitle()));
     return;
   }
  
  AliPHOSLoader* gime = dynamic_cast<AliPHOSLoader*>(runget->GetLoader("PHOSLoader"));
  if ( gime == 0 ) 
   {
     AliError("Could not obtain the Loader object !"); 
     return ;
   } 

  TClonesArray * recParticles = gime->RecParticles() ; 
  recParticles->Expand(recParticles->GetEntriesFast() ) ;

  TTree * treeR = gime->TreeR();
  
  if(!treeR){
    gime->MakeTree("R");
    treeR = gime->TreeR() ;
  }
  
  //First rp
  Int_t bufferSize = 32000 ;    
  TBranch * rpBranch = treeR->Branch("PHOSRP",&recParticles,bufferSize);
  rpBranch->SetTitle(fRecParticlesTitle);
  
  //second, pid
  Int_t splitlevel = 0 ; 
  AliPHOSPIDv0 * pid = this ;
  TBranch * pidBranch = treeR->Branch("AliPHOSPID","AliPHOSPIDv0",&pid,bufferSize,splitlevel);
  pidBranch->SetTitle(fRecParticlesTitle.Data());
  
  rpBranch->Fill() ;
  pidBranch->Fill() ;

  gime->WriteRecParticles("OVERWRITE");
  gime->WritePID("OVERWRITE");

}

//____________________________________________________________________________
void  AliPHOSPIDv0::PlotDispersionCuts()const
{
  // produces a plot of the dispersion cut
  TCanvas*  lambdas = new TCanvas("lambdas","Cuts on the ellipse axis",200,10,700,500);
 
  if(fIDOptions.Contains("ell",TString::kIgnoreCase ) ){
    TF2 * ell = new TF2("Elliptic Cuts",fFormula->GetName(),0,3,0,3) ;
    ell->SetMinimum(0.0000001) ;
    ell->SetMaximum(0.001) ;
    ell->SetLineStyle(1) ;
    ell->SetLineWidth(2) ;
    ell->Draw() ;
  }
  
  if( fIDOptions.Contains("dis",TString::kIgnoreCase ) ){
    TF2 * dsp = new TF2("dispersion","(y<x)*(x*x+y*y < [0]*[0])",0,3,0,3) ;
    dsp->SetParameter(0,fDispersion) ;
    dsp->SetMinimum(0.0000001) ;
    dsp->SetMaximum(0.001) ;
    dsp->SetLineStyle(1) ;
    dsp->SetLineColor(2) ;
    dsp->SetLineWidth(2) ;
    dsp->SetNpx(200) ;
    dsp->SetNpy(200) ;
    if(fIDOptions.Contains("ell",TString::kIgnoreCase ) )
      dsp->Draw("same") ;
    else
      dsp->Draw() ;
  }
  lambdas->Update();
}

//____________________________________________________________________________
TVector3 AliPHOSPIDv0::GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * )const 
{ 
  // Calculates the momentum direction:
  //   1. if only a EMC RecPoint, direction is given by IP and this RecPoint
  //   2. if a EMC RecPoint and CPV RecPoint, direction is given by the line through the 2 recpoints 
  //  However because of the poor position resolution of PPSD the direction is always taken as if we were 
  //  in case 1.

  TVector3 dir(0,0,0) ; 
  
  TVector3 emcglobalpos ;
  TMatrix  dummy ;
  
  emc->GetGlobalPosition(emcglobalpos, dummy) ;
  
 
  // The following commented code becomes valid once the PPSD provides 
  // a reasonable position resolution, at least as good as EMC ! 
  //   TVector3 ppsdlglobalpos ;
  //   TVector3 ppsduglobalpos ;
  //   if( fPpsdLowRecPoint ){ // certainly a photon that has concerted
  //     fPpsdLowRecPoint->GetGlobalPosition(ppsdlglobalpos, mdummy) ; 
  //     dir = emcglobalpos -  ppsdlglobalpos ; 
  //     if( fPpsdUpRecPoint ){ // not looks like a charged       
  //        fPpsdUpRecPoint->GetGlobalPosition(ppsduglobalpos, mdummy) ; 
  //        dir = ( dir +  emcglobalpos -  ppsduglobalpos ) * 0.5 ; 
  //      }
  //   }
  //   else { // looks like a neutral
  //    dir = emcglobalpos ;  
  //  }
  
  dir = emcglobalpos ;  
  dir.SetZ( -dir.Z() ) ;   // why ?  
  dir.SetMag(1.) ;

  //account correction to the position of IP
  Float_t xo,yo,zo ; //Coordinates of the origin
  gAlice->Generator()->GetOrigin(xo,yo,zo) ;
  TVector3 origin(xo,yo,zo);
  dir = dir - origin ;

  return dir ;  
}
//____________________________________________________________________________
void AliPHOSPIDv0::PrintRecParticles(Option_t * option)
{
  // Print table of reconstructed particles

  AliRunLoader* runget = AliRunLoader::GetRunLoader(GetTitle());
  if(runget == 0x0) 
   {
     AliError(Form("Can not find run getter in event folder \"%s\"",
		   GetTitle()));
     return;
   }
  
  AliPHOSLoader* gime = dynamic_cast<AliPHOSLoader*>(runget->GetLoader("PHOSLoader"));
  if ( gime == 0 ) 
   {
     AliError("Could not obtain the Loader object !"); 
     return ;
   } 

  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;
  TClonesArray * recParticles = gime->RecParticles() ; 
  
  TString message ; 
  message  = "event %d\n" ; 
  message += "       found %d RecParticles\n" ; 
  AliInfo(Form(message.Data(), 
	       gAlice->GetEvNumber(), recParticles->GetEntriesFast() )) ;   

  if(strstr(option,"all")) {  // printing found TS
    AliInfo("  PARTICLE   Index"  ) ; 
   
    Int_t index ;
    for (index = 0 ; index < recParticles->GetEntries() ; index++) {
      AliPHOSRecParticle * rp = (AliPHOSRecParticle * ) recParticles->At(index) ;       
      
      Text_t particle[11];
      switch(rp->GetType()) {
      case  AliPHOSFastRecParticle::kNEUTRALEMFAST:
	strcpy( particle, "NEUTRAL EM FAST");
	break;
      case  AliPHOSFastRecParticle::kNEUTRALHAFAST:
	strcpy(particle, "NEUTRAL HA FAST");
	break;
      case  AliPHOSFastRecParticle::kNEUTRALEMSLOW:
	strcpy(particle, "NEUTRAL EM SLOW");
	break ;
      case  AliPHOSFastRecParticle::kNEUTRALHASLOW: 
	strcpy(particle, "NEUTRAL HA SLOW");
	break ;
      case  AliPHOSFastRecParticle::kCHARGEDEMFAST:
	strcpy(particle, "CHARGED EM FAST") ;
	break ;
      case  AliPHOSFastRecParticle::kCHARGEDHAFAST:
	strcpy(particle, "CHARGED HA FAST") ;
	break ;	
      case  AliPHOSFastRecParticle::kCHARGEDEMSLOW:
	strcpy(particle, "CHARGEDEMSLOW") ;
	break ;
      case  AliPHOSFastRecParticle::kCHARGEDHASLOW:
	strcpy(particle, "CHARGED HA SLOW") ;
	break ; 
      }
      
      //    Int_t * primaries; 
      //    Int_t nprimaries;
      //    primaries = rp->GetPrimaries(nprimaries);
      
      AliInfo(Form("          %s     %d",  
		   particle, rp->GetIndexInList())) ;
    }
  }  
}



