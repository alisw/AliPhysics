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
// Implementation version v1 of the PHOS particle identifier 
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
//  root [0] AliPHOSPIDv1 * p1 = new AliPHOSPIDv1("galice.root")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//  root [1] p1->SetIdentificationMethod("disp ellipse")
//  root [2] p1->ExecuteTask()
//  root [3] AliPHOSPIDv1 * p2 = new AliPHOSPIDv1("galice1.root","ts1")
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
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TF2.h"
#include "TFormula.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TSystem.h"
#include "TBenchmark.h"
// --- Standard library ---

#include <iostream.h>
#include <iomanip.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliGenerator.h"
#include "AliPHOS.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSRecParticle.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSGetter.h"

ClassImp( AliPHOSPIDv1) 

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1():AliPHOSPID()
{ 
  // default ctor
  fFormula           = 0 ;
  fDispersion        = 0. ; 
  fCpvEmcDistance    = 0 ; 
  fTimeGate          = 2.e-9 ;
  fHeaderFileName    = "" ; 
  fTrackSegmentsTitle= "" ; 
  fRecPointsTitle    = "" ; 
  fRecParticlesTitle = "" ; 
  fIDOptions         = "dis time" ; 
  fRecParticlesInRun = 0 ;
  fClusterizer = 0;
  fTSMaker = 0;
}

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1(const char * headerFile,const char * name) : AliPHOSPID(headerFile, name)
{ 
  //ctor with the indication on where to look for the track segments

  fFormula        = new TFormula("LambdaCuts","(x>1)*(x<2.5)*(y>0)*(y<x)") ;   
  fDispersion     = 2.0 ; 
  fCpvEmcDistance = 3.0 ;
  fTimeGate          = 2.e-9 ;
 
  fHeaderFileName     = GetTitle() ; 
  fTrackSegmentsTitle = GetName() ; 
  fRecPointsTitle     = GetName() ; 
  fRecParticlesTitle  = GetName() ; 
  fIDOptions          = "dis time" ;
    
  TString tempo(GetName()) ; 
  tempo.Append(":") ;
  tempo.Append(Version()) ; 
  SetName(tempo) ; 
  fRecParticlesInRun = 0 ; 

  Init() ;

}

//____________________________________________________________________________
AliPHOSPIDv1::~AliPHOSPIDv1()
{ 
}


//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetDistance(AliPHOSEmcRecPoint * emc,AliPHOSRecPoint * cpv, Option_t *  Axis)const
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
 
  const AliPHOSGeometry * geom = AliPHOSGetter::GetInstance()->PHOSGeometry() ; 
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
void  AliPHOSPIDv1::Exec(Option_t * option) 
{
  //Steering method
  
  if( strcmp(GetName(), "")== 0 ) 
    Init() ;
  
  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSPID");
  
  if(strstr(option,"print")) {
    Print("") ; 
    return ; 
  }

  gAlice->GetEvent(0) ;
  //check, if the branch with name of this" already exits?
  TObjArray * lob = (TObjArray*)gAlice->TreeR()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ;  
  Bool_t phospidfound = kFALSE, pidfound = kFALSE ; 
  
  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;

  while ( (branch = (TBranch*)next()) && (!phospidfound || !pidfound) ) {
    if ( (strcmp(branch->GetName(), "PHOSPID")==0) && (strcmp(branch->GetTitle(), taskName.Data())==0) ) 
      phospidfound = kTRUE ;
    
    else if ( (strcmp(branch->GetName(), "AliPHOSPID")==0) && (strcmp(branch->GetTitle(), taskName.Data())==0) ) 
      pidfound = kTRUE ; 
  }

  if ( phospidfound || pidfound ) {
    cerr << "WARNING: AliPHOSPIDv1::Exec -> RecParticles and/or PIDtMaker branch with name " 
	 << taskName.Data() << " already exits" << endl ;
    return ; 
  }       
  
  Int_t nevents = (Int_t) gAlice->TreeE()->GetEntries() ;
  Int_t ievent ;
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  
  for(ievent = 0; ievent < nevents; ievent++){
    gime->Event(ievent,"R") ;
    
    MakeRecParticles() ;
    
    WriteRecParticles(ievent);
    
    if(strstr(option,"deb"))
      PrintRecParticles(option) ;

    //increment the total number of rec particles per run 
    fRecParticlesInRun += gime->RecParticles()->GetEntriesFast() ; 

  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSPID");
    cout << "AliPHOSPID:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("PHOSPID") << " seconds for PID " 
	 <<  gBenchmark->GetCpuTime("PHOSPID")/nevents << " seconds per event " << endl ;
    cout << endl ;
  }
  
}
//____________________________________________________________________________
void AliPHOSPIDv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  // Add the PID task to the list of PHOS tasks
  
  if ( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;
  
  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;
  
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(GetTitle(), taskName.Data()) ; 
  if ( gime == 0 ) {
    cerr << "ERROR: AliPHOSPIDv1::Init -> Could not obtain the Getter object !" << endl ; 
    return ;
  } 
   
  gime->PostPID(this) ;
  // create a folder on the white board //YSAlice/WhiteBoard/RecParticles/PHOS/recparticlesName
  gime->PostRecParticles(taskName.Data() ) ; 
  
}

//____________________________________________________________________________
void  AliPHOSPIDv1::MakeRecParticles(){

  // Makes a RecParticle out of a TrackSegment
  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  TObjArray * emcRecPoints = gime->EmcRecPoints(taskName) ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints(taskName) ; 
  TClonesArray * trackSegments = gime->TrackSegments(taskName) ; 
  TClonesArray * recParticles  = gime->RecParticles(taskName) ; 
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
    rp->SetTraskSegment(index) ;
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
void  AliPHOSPIDv1:: Print(Option_t * option) const
{
  // Print the parameters used for the particle type identification
    cout <<  "=============== AliPHOSPID1 ================" << endl ;
    cout <<  "Making PID "<< endl ;
    cout <<  "    Headers file:               " << fHeaderFileName.Data() << endl ;
    cout <<  "    RecPoints branch title:     " << fRecPointsTitle.Data() << endl ;
    cout <<  "    TrackSegments Branch title: " << fTrackSegmentsTitle.Data() << endl ;
    cout <<  "    RecParticles Branch title   " << fRecParticlesTitle.Data() << endl;
    cout <<  "with parameters: " << endl ;
    cout <<  "    Maximal EMC - CPV  distance (cm) " << fCpvEmcDistance << endl ;
    if(fIDOptions.Contains("dis",TString::kIgnoreCase ))
      cout <<  "                    dispersion cut " << fDispersion << endl ;
    if(fIDOptions.Contains("ell",TString::kIgnoreCase )){
      cout << "             Eliptic cuts function: " << endl ;
      cout << fFormula->GetTitle() << endl ;
    }
    if(fIDOptions.Contains("tim",TString::kIgnoreCase ))
      cout << "             Time Gate uzed: " << fTimeGate <<  endl ;
    cout <<  "============================================" << endl ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetShowerProfileCut(char * formula)
{
  //set shape of the cut on the axis of ellipce, drown around shouer
  //shower considered "narrow" if Formula(lambda[0],lambda[1]) > 0.
  if(fFormula) 
    delete fFormula; 
  fFormula = new TFormula("Lambda Cut",formula) ;
}
//____________________________________________________________________________
void  AliPHOSPIDv1::WriteRecParticles(Int_t event)
{
 
  AliPHOSGetter *gime = AliPHOSGetter::GetInstance() ; 
  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;
  TClonesArray * recParticles = gime->RecParticles(taskName) ; 
  recParticles->Expand(recParticles->GetEntriesFast() ) ;

  //Make branch in TreeR for RecParticles 
  char * filename = 0;
  if(gSystem->Getenv("CONFIG_SPLIT_FILE")!=0){   //generating file name
    filename = new char[strlen(gAlice->GetBaseFile())+20] ;
    sprintf(filename,"%s/PHOS.Reco.root",gAlice->GetBaseFile()) ; 
  }
  
  TDirectory *cwd = gDirectory;
  
  //First rp
  Int_t bufferSize = 32000 ;    
  TBranch * rpBranch = gAlice->TreeR()->Branch("PHOSRP",&recParticles,bufferSize);
  rpBranch->SetTitle(fRecParticlesTitle);
  if (filename) {
    rpBranch->SetFile(filename);
    TIter next( rpBranch->GetListOfBranches());
    TBranch * sb ;
    while ((sb=(TBranch*)next())) {
      sb->SetFile(filename);
    }   
    cwd->cd();
  }
  
  //second, pid
  Int_t splitlevel = 0 ; 
  AliPHOSPIDv1 * pid = this ;
  TBranch * pidBranch = gAlice->TreeR()->Branch("AliPHOSPID","AliPHOSPIDv1",&pid,bufferSize,splitlevel);
  pidBranch->SetTitle(fRecParticlesTitle.Data());
  if (filename) {
    pidBranch->SetFile(filename);
    TIter next( pidBranch->GetListOfBranches());
    TBranch * sb ;
    while ((sb=(TBranch*)next())) {
      sb->SetFile(filename);
    }   
    cwd->cd();
  }    
  
  rpBranch->Fill() ;
  pidBranch->Fill() ;
  
  gAlice->TreeR()->Write(0,kOverwrite) ;  
  
}

//____________________________________________________________________________
void  AliPHOSPIDv1::PlotDispersionCuts()const
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
TVector3 AliPHOSPIDv1::GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv)const 
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
void AliPHOSPIDv1::PrintRecParticles(Option_t * option)
{
  // Print table of reconstructed particles

  AliPHOSGetter *gime = AliPHOSGetter::GetInstance() ; 

  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;
  TClonesArray * recParticles = gime->RecParticles(taskName) ; 
  
  cout << "AliPHOSPIDv1: event "<<gAlice->GetEvNumber()  << endl ;
  cout << "       found " << recParticles->GetEntriesFast() << " RecParticles " << endl ;
  
  if(strstr(option,"all")) {  // printing found TS
    
    cout << "  PARTICLE "   
	 << "  Index    "  << endl ;
      //	 << "  X        "     
      //	 << "  Y        " 
      //	 << "  Z        "    
      //	 << " # of primaries "          
      //	 << " Primaries list "    <<  endl;      
    
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
      
      cout << setw(10) << particle << "  "
	   << setw(5) <<  rp->GetIndexInList() << " "  ;
	//	   << setw(4) <<  nprimaries << "  ";
	//      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++)
	//	cout << setw(4)  <<  primaries[iprimary] << " ";
      cout << endl;  	 
    }
    cout << "-------------------------------------------" << endl ;
  }
  
}



