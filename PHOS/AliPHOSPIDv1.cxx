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
//     - Preshower information (in MIX or GPS2 geometries)
//     - shower width.

// CPV or Preshower cluster should be clother in PHOS plane than fCpvEmcDistance (in cm).
// This variable can be set by method SetCpvtoEmcDistanceCut(Float_t cut)  
//
// One can set desirable ID method by the function SetIdentificationMethod(option).
// Now the following options can be used together or separately :
//     - "disp": use dispersion cut on shower width 
//               (width can be set by method SetDispersionCut(Float_t cut)
//     - "ell" : use cut on the axis of the ellipse, drown around shower 
//       (this cut can be changed by SetShowerProfileCut(char* formula), 
//        where formula - any function of two variables f(lambda[0],lambda[1]).
//        Shower is considered as EM if f() > 0 )
// One can visualize current cuts calling method PlotDispersionCuts().    
//
// Below we present usercase:
// root [0] AliPHOSPIDv1 * p1 = new AliPHOSPIDv1("galice.root")
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
// root [1] p1->SetIdentificationMethod("disp ellipse")
// root [2] p1->ExecuteTask()
// root [3] AliPHOSPIDv1 * p2 = new AliPHOSPIDv1("galice1.root","ts1")
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                // reading headers from file galice1.root and TrackSegments 
//                // with title "ts1"
// root [4] p2->SetRecParticlesBranch("rp1")
//                // set file name for the branch RecParticles
// root [5] p2->ExecuteTask("deb all time")
//                // available options
//                // "deb" - prints # of reconstructed particles
//                // "deb all" -  prints # and list of RecParticles
//                // "time" - prints benchmarking results
//                  
//*-- Author: Yves Schutz (SUBATECH)  & Gines Martinez (SUBATECH) & 
//            Dmitri Peressounko (SUBATECH & Kurchatov Institute)
//            Complitely redesined by Dmitri Peressounko, March 2001

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

ClassImp( AliPHOSPIDv1) 

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1():AliPHOSPID()
{ 
  fIsInitialized = kFALSE ;
}

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1(const char * headeFile,const char * tsBranchTitle):AliPHOSPID()
{ 
  
  fHeaderFileName = headeFile ;

  fTSTitle = tsBranchTitle ;

  SetName("AliPHOSPID") ;
  SetTitle("version1") ;

  TFile * file = (TFile*) gROOT->GetFile(fHeaderFileName.Data() ) ;

  if(file == 0){
    file = new TFile(fHeaderFileName.Data(),"update") ;
    gAlice = (AliRun *) file->Get("gAlice") ;
  }
  
  AliPHOS * phos = (AliPHOS *) gAlice->GetDetector("PHOS") ;    
  fGeom  = AliPHOSGeometry::GetInstance(phos->GetGeometry()->GetName(),phos->GetGeometry()->GetTitle() );
  
  fTrackSegments = new TClonesArray("AliPHOSTrackSegment",1) ;
  fTSMaker       = 0 ;
  fEmcRecPoints  = new TObjArray(1) ;
  fCpvRecPoints  = new TObjArray(1) ;
  fClusterizer   = 0 ;
  fRecParticles  = new TClonesArray("AliPHOSRecParticle",100) ;

  fFormula = new TFormula("LambdaCuts","(x>1)*(x<3)*(y>0)*(y<x)") ;
  
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 

  fDispersion = 2.0; 
  fCpvEmcDistance = 3.0 ;
  fIsInitialized = kTRUE ;

}
//____________________________________________________________________________
AliPHOSPIDv1::~AliPHOSPIDv1()
{ 

}
//____________________________________________________________________________
void AliPHOSPIDv1::Init()
{
  if(!fIsInitialized){
    if(fHeaderFileName.IsNull())
      fHeaderFileName = "galice.root" ;
    
    TFile * file = (TFile*) gROOT->GetFile(fHeaderFileName.Data() ) ;

    if(file == 0){
      file = new TFile(fHeaderFileName.Data(),"update") ;
      gAlice = (AliRun *) file->Get("gAlice") ;
    }

    AliPHOS * phos = (AliPHOS *) gAlice->GetDetector("PHOS") ;    
    fGeom  = AliPHOSGeometry::GetInstance(phos->GetGeometry()->GetName(),phos->GetGeometry()->GetTitle() );

    fTrackSegments = new TClonesArray("AliPHOSTrackSegment",1) ;
    fTSMaker       = new AliPHOSTrackSegmentMakerv1() ;
    fEmcRecPoints  = new TObjArray(1) ;
    fCpvRecPoints  = new TObjArray(1) ;
    fClusterizer   = new AliPHOSClusterizerv1() ;
    fRecParticles  = new TClonesArray("AliPHOSRecParticle",100) ;
    
    fFormula = new TFormula("LambdaCuts","(x>1)*(x<2.5)*(y>0)*(y<x)") ;
    
    // add Task to //root/Tasks folder
    TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
    roottasks->Add(this) ; 

    fDispersion = 2.0; 
    fCpvEmcDistance = 3.0 ;

    fIsInitialized = kTRUE ;
  }


}
//____________________________________________________________________________
Bool_t AliPHOSPIDv1::ReadTrackSegments()
{
  //Fist read Track Segment Branch and extract RecPointsBranch from fTSMaker

  fTrackSegments->Clear() ; 
  fEmcRecPoints->Clear() ;
  fCpvRecPoints->Clear() ;
  fRecParticles->Clear() ;

  gAlice->GetEvent(fNEvent) ;

  TTree * treeR = gAlice->TreeR()  ; 

  if(treeR==0){
    char treeName[20]; 
    sprintf(treeName,"TreeR%d",fNEvent);
    cout << "Error in AliPHOSClusterizerv1 : no "<<treeName << endl  ;
    cout << "   Do nothing " << endl ;
    return kFALSE ;
  }

  //first read TSMaker branch and extract information about RecPoints Branches
  TBranch * tsMakerBranch = 0;
  TBranch * tsBranch = 0;

  TObjArray * branches = treeR->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t tsMakerNotFound = kTRUE ;
  Bool_t tsNotFound = kTRUE ;
  
  for(ibranch = 0;(ibranch <branches->GetEntries())&&(tsMakerNotFound||tsNotFound);ibranch++){
    if(tsMakerNotFound){
      tsMakerBranch=(TBranch *) branches->At(ibranch) ;
      if( fTSTitle.CompareTo(tsMakerBranch->GetTitle())==0 )
	if( strcmp(tsMakerBranch->GetName(),"AliPHOSTrackSegmentMaker") == 0) 
	  tsMakerNotFound = kFALSE ;
    }
    if(tsNotFound){
      tsBranch=(TBranch *) branches->At(ibranch) ;
      if( fTSTitle.CompareTo(tsBranch->GetTitle())==0 )
	if( strcmp(tsBranch->GetName(),"PHOSTS") == 0) 
	  tsNotFound = kFALSE ;
    }
  }
  
  if(tsMakerNotFound ||tsNotFound ){
    cout << "Can't find Branch with TrackSegmentMaker and TrackSegments " ;
    cout << "Do nothing" <<endl  ;
    return kFALSE ;
  }

  tsMakerBranch->SetAddress(&fTSMaker) ;
  tsBranch->SetAddress(&fTrackSegments) ;

  treeR->GetEvent(0) ;

  fRecPointsTitle = fTSMaker->GetRecPointsBranch() ;

  //reading now recponts branches
  TBranch * emcBranch = 0;
  TBranch * cpvBranch = 0;
  TBranch * cluBranch = 0;

  Bool_t emcNotFound = kTRUE ;
  Bool_t cpvNotFound = kTRUE ;
  Bool_t cluNotFound = kTRUE ;
 
  for(ibranch = 0;(ibranch <branches->GetEntries())&&(emcNotFound||cpvNotFound||cluNotFound);ibranch++){
    if(emcNotFound){
      emcBranch=(TBranch *) branches->At(ibranch) ;
      if( fRecPointsTitle.CompareTo(emcBranch->GetTitle())==0 )
	if( strcmp(emcBranch->GetName(),"PHOSEmcRP") == 0) 
	  emcNotFound = kFALSE ;
    }
    if(cpvNotFound){
      cpvBranch=(TBranch *) branches->At(ibranch) ;
      if( fRecPointsTitle.CompareTo(cpvBranch->GetTitle())==0 )
	if( strcmp(cpvBranch->GetName(),"PHOSCpvRP") == 0) 
	  cpvNotFound = kFALSE ;
    }
    if(cluNotFound){
      cluBranch=(TBranch *) branches->At(ibranch) ;
      if( fRecPointsTitle.CompareTo(cluBranch->GetTitle())==0 )
	if( strcmp(cluBranch->GetName(),"AliPHOSClusterizer") == 0) 
	  cluNotFound = kFALSE ;
    }
  }
  
  if(emcNotFound ||cpvNotFound ||cluNotFound ){
    cout << "Can't find Branch with RecPoints or AliPHOSClusterizer " ;
    cout << "Do nothing" <<endl  ;
    return kFALSE ;
  }
  
  emcBranch->SetAddress(&fEmcRecPoints) ;
  cpvBranch->SetAddress(&fCpvRecPoints) ;
  cluBranch->SetAddress(&fClusterizer) ;

  treeR->GetEvent(0) ;
  return kTRUE ;



}
//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetDistance(AliPHOSEmcRecPoint * emc,AliPHOSRecPoint * cpv, Option_t *  Axis)const
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
 
  TVector3 vecEmc ;
  TVector3 vecCpv ;
  
  emc->GetLocalPosition(vecEmc) ;
  cpv->GetLocalPosition(vecCpv) ; 
  if(emc->GetPHOSMod() == cpv->GetPHOSMod()){ 
    
    // Correct to difference in CPV and EMC position due to different distance to center.
    // we assume, that particle moves from center
    Float_t dCPV = fGeom->GetIPtoOuterCoverDistance();
    Float_t dEMC = fGeom->GetIPtoCrystalSurface() ;
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
  if(!fIsInitialized) 
    Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSPID");


  Int_t nEvents = (Int_t) gAlice->TreeE()->GetEntries() ;

  for(fNEvent = 0 ;fNEvent <nEvents; fNEvent++){
    if(!ReadTrackSegments())
      return ;
    MakeRecParticles() ;
    WriteRecParticles();
    if(strstr(option,"deb"))
      PrintRecParticles(option) ;
  }

  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSPID");
    cout << "AliPHOSPID:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("PHOSPID") << " seconds for PID " 
	 <<  gBenchmark->GetCpuTime("PHOSPID")/nEvents << " seconds per event " << endl ;
    cout << endl ;
  }

}
//____________________________________________________________________________
void  AliPHOSPIDv1::MakeRecParticles(){

  // Makes a RecParticle out of a TrackSegment

  TIter next(fTrackSegments) ; 
  AliPHOSTrackSegment * ts ; 
  Int_t index = 0 ; 
  AliPHOSRecParticle * rp ; 
  
  Bool_t ellips = fIDOptions.Contains("ell",TString::kIgnoreCase ) ;
  Bool_t disp   = fIDOptions.Contains("dis",TString::kIgnoreCase ) ;
  
  while ( (ts = (AliPHOSTrackSegment *)next()) ) {
    
    new( (*fRecParticles)[index] ) AliPHOSRecParticle() ;
    rp = (AliPHOSRecParticle *)fRecParticles->At(index) ; 
    rp->SetTraskSegment(index) ;

    AliPHOSEmcRecPoint * emc = 0 ;
    if(ts->GetEmcIndex()>=0)
      emc = (AliPHOSEmcRecPoint *) fEmcRecPoints->At(ts->GetEmcIndex()) ;
    
    AliPHOSRecPoint    * cpv = 0 ;
    if(ts->GetCpvIndex()>=0)
      cpv = (AliPHOSRecPoint *)   fCpvRecPoints->At(ts->GetCpvIndex()) ;
    
    AliPHOSRecPoint    * ppsd = 0 ;
    if(ts->GetPpsdIndex()>=0)
      ppsd= (AliPHOSRecPoint *)   fCpvRecPoints->At(ts->GetPpsdIndex()) ;

    //set momentum and energy first
    Float_t    e = emc->GetEnergy() ;
    TVector3 dir = GetMomentumDirection(emc,cpv,ppsd) ; 
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
    
    
    // Looking at the photon conversion detector
    Int_t pcdetector= 0 ;  //1 hit and 0 no hit
    if(ppsd)
      if(GetDistance(emc, ppsd, "R") < fCpvEmcDistance) 
	pcdetector = 1 ;  
    
    // Looking at the CPV detector
    Int_t cpvdetector= 0 ;  //1 hit and 0 no hit     
    if(cpv)
      if(GetDistance(emc, cpv,  "R") < fCpvEmcDistance) 
	cpvdetector = 1 ;  
    
    Int_t type = showerprofile + 2 * pcdetector + 4 * cpvdetector ;
    rp->SetType(type) ; 
    index++ ; 
  }
  
}

//____________________________________________________________________________
void  AliPHOSPIDv1:: Print(Option_t * option) const
{
  // Print the parameters used for the particle type identification
  
  cout << "AliPHOSPIDv1 : cuts for the particle idendification based on the shower profile " << endl ;

  cout << "Eliptic cuts function " << endl ;
  cout << "    " << fFormula->GetTitle() << endl ;

}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetShowerProfileCut(char * formula){
  //set shape of the cut on the axis of ellipce, drown around shouer
  //shower considered "narrow" if Formula(lambda[0],lambda[1]) > 0.
  if(fFormula) 
    delete fFormula; 
  fFormula = new TF2("Lambda Cut",formula,0,3,0,3) ;
}
//____________________________________________________________________________
void  AliPHOSPIDv1::WriteRecParticles()
{

  //check, if these branches already exist  
  TBranch * pidBranch = 0;
  TBranch * rpBranch = 0;

  TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t pidNotFound = kTRUE ;
  Bool_t rpNotFound = kTRUE ;
  
  for(ibranch = 0;(ibranch <branches->GetEntries())&& pidNotFound && rpNotFound;ibranch++){
    if(pidNotFound){
      pidBranch=(TBranch *) branches->At(ibranch) ;
      if( (strcmp(pidBranch->GetName(),"PHOSPID") == 0) &&
	  (fRecparticlesTitle.CompareTo(pidBranch->GetTitle()) == 0) )
	pidNotFound = kFALSE ;
    }
    if(rpNotFound){
      rpBranch=(TBranch *) branches->At(ibranch) ;
      if( (strcmp(rpBranch->GetName(),"PHOSRP") == 0) &&
	  (fRecparticlesTitle.CompareTo(rpBranch->GetTitle())==0 ))
	rpNotFound = kFALSE ;
    }
  }
  
  if(!pidNotFound || !rpNotFound) {
    cout << "AliPHOSPIDv1 error: " << endl ;
    cout << "       Branch PHOSRP and PHOSPID with title '"<<fRecparticlesTitle.Data()<<"' already exist "<< endl ;
    cout << "       can not overwrite " << endl ;
    return ;
  }

  //Make branch in TreeR for TrackSegments 
  char * filename = 0;
  if(gSystem->Getenv("CONFIG_SPLIT_FILE")!=0){   //generating file name
    filename = new char[strlen(gAlice->GetBaseFile())+20] ;
    sprintf(filename,"%s/PHOS.Reco.root",gAlice->GetBaseFile()) ; 
  }

  TDirectory *cwd = gDirectory;
  
  //First rp
  Int_t bufferSize = 32000 ;    
  rpBranch = gAlice->TreeR()->Branch("PHOSRP",&fRecParticles,bufferSize);
  rpBranch->SetTitle(fRecparticlesTitle.Data());
  if (filename) {
    rpBranch->SetFile(filename);
    TIter next( rpBranch->GetListOfBranches());
    while ((rpBranch=(TBranch*)next())) {
      rpBranch->SetFile(filename);
    }   
    cwd->cd();
  }

  //second, pid
  Int_t splitlevel = 0 ; 
  AliPHOSPIDv1 * pid = this ;
  pidBranch = gAlice->TreeR()->Branch("AliPHOSPID","AliPHOSPIDv1",&pid,bufferSize,splitlevel);
  pidBranch->SetTitle(fRecparticlesTitle.Data());
  if (filename) {
    pidBranch->SetFile(filename);
    TIter next( pidBranch->GetListOfBranches());
    while ((pidBranch=(TBranch*)next())) {
      pidBranch->SetFile(filename);
    }   
    cwd->cd();
  }    
  
  gAlice->TreeR()->Fill() ;    
  gAlice->TreeR()->Write(0,kOverwrite) ;  
  
}
//____________________________________________________________________________
void  AliPHOSPIDv1::PlotDispersionCuts()const
{
  TCanvas*  lambdas = new TCanvas("lambdas","Cuts on the elipse axise",200,10,700,500);
  
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
TVector3 AliPHOSPIDv1::GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv,AliPHOSRecPoint * ppsd)const 
{ 
  // Calculates the momentum direction:
  //   1. if only a EMC RecPoint, direction is given by IP and this RecPoint
  //   2. if a EMC RecPoint and one PPSD RecPoint, direction is given by the line through the 2 recpoints 
  //   3. if a EMC RecPoint and two PPSD RecPoints, dirrection is given by the average line through 
  //      the 2 pairs of recpoints  
  // However because of the poor position resolution of PPSD the direction is always taken as if we were 
  //  in case 1.

  TVector3 dir(0,0,0) ; 
  
  TVector3 emcglobalpos ;
  TMatrix  dummy ;
  
  emc->GetGlobalPosition(emcglobalpos, dummy) ;
  
 
  // The following commeneted code becomes valid once the PPSD provides 
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
void AliPHOSPIDv1::PrintRecParticles(Option_t * option){

  cout << "AliPHOSPIDv1: " << endl ;
  cout << "       found " << fRecParticles->GetEntriesFast() << " RecParticles " << endl ;

  if(strstr(option,"all")) {  // printing found TS
    
    cout << "  PARTICLE "   
	 << "  Index    "  << endl ;
      //	 << "  X        "     
      //	 << "  Y        " 
      //	 << "  Z        "    
      //	 << " # of primaries "          
      //	 << " Primaries list "    <<  endl;      
    
    Int_t index ;
    for (index = 0 ; index < fRecParticles->GetEntries() ; index++) {
      AliPHOSRecParticle * rp = (AliPHOSRecParticle * ) fRecParticles->At(index) ;       
      
      Text_t particle[11];
      switch(rp->GetType()) {
      case  AliPHOSFastRecParticle::kNEUTRALEM:
	strcpy( particle, "NEUTRAL_EM");
	break;
      case  AliPHOSFastRecParticle::kNEUTRALHA:
	strcpy(particle, "NEUTRAL_HA");
	break;
      case  AliPHOSFastRecParticle::kGAMMA:
	strcpy(particle, "GAMMA");
	break ;
      case  AliPHOSFastRecParticle::kGAMMAHA: 
	strcpy(particle, "GAMMA_H");
	break ;
      case  AliPHOSFastRecParticle::kABSURDEM:
	strcpy(particle, "ABSURD_EM") ;
	break ;
      case  AliPHOSFastRecParticle::kABSURDHA:
	strcpy(particle, "ABSURD_HA") ;
	break ;	
      case  AliPHOSFastRecParticle::kELECTRON:
	strcpy(particle, "ELECTRON") ;
	break ;
      case  AliPHOSFastRecParticle::kCHARGEDHA:
	strcpy(particle, "CHARGED_HA") ;
	break ; 
      }
      
      //    Int_t * primaries; 
      //    Int_t nprimaries;
      //    primaries = rp->GetPrimaries(nprimaries);
      
      cout << setw(15) << particle << "  "
	   << setw(3) <<  rp->GetIndexInList() << " "  ;
	//	   << setw(4) <<  nprimaries << "  ";
	//      for (Int_t iprimary=0; iprimary<nprimaries; iprimary++)
	//	cout << setw(4)  <<  primaries[iprimary] << " ";
      cout << endl;  	 
    }
    cout << "-------------------------------------------" << endl ;
  }
  
}



