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
// Implementation version v2 of the PHOS particle identifier 
// Particle identification based on the 
//     - RCPV: distance from CPV recpoint to EMCA recpoint.
//     - TOF 
//     - PCA: Principal Components Analisis..
// The identified particle has an identification number corresponding 
// to a 3 bits number:
//     -Bit 0: bit set if RCPV > fCpvEmcDistance 
//     -Bit 1: bit set if TOF  < fTimeGate
//     -Bit 2: bit set if Principal Components are 
//      inside an ellipse defined by fX_center, fY_center, fA, fB, fAngle
//
//
// use case:
//  root [0] AliPHOSPIDv2 * p1 = new AliPHOSPIDv2("galice.root")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//  root [1] p1->SetIdentificationMethod("disp ellipse")
//  root [2] p1->ExecuteTask()
//  root [3] AliPHOSPIDv2 * p2 = new AliPHOSPIDv2("galice1.root","ts1")
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
//            Gustavo Conesa April 2002

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
#include "TEllipse.h"
#include "TPrincipal.h"
// --- Standard library ---

#include <iostream.h>
#include <iomanip.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliGenerator.h"
#include "AliPHOS.h"
#include "AliPHOSPIDv2.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSRecParticle.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSGetter.h"

ClassImp( AliPHOSPIDv2) 

//____________________________________________________________________________
AliPHOSPIDv2::AliPHOSPIDv2():AliPHOSPID()
{ 
  // default ctor
 
  fHeaderFileName    = "" ; 
  fTrackSegmentsTitle= "" ; 
  fRecPointsTitle    = "" ; 
  fRecParticlesTitle = "" ; 
  fRecParticlesInRun = 0 ;
  fClusterizer       = 0 ; 
  fTSMaker           = 0 ;
  
}

//____________________________________________________________________________
AliPHOSPIDv2::AliPHOSPIDv2(const char * headerFile,const char * name) : AliPHOSPID(headerFile, name)
{ 
  //ctor with the indication on where to look for the track segments
 
  fHeaderFileName     = GetTitle() ; 
  fTrackSegmentsTitle = GetName() ; 
  fRecPointsTitle     = GetName() ; 
  fRecParticlesTitle  = GetName() ;
  TString tempo(GetName()) ; 
  tempo.Append(":") ;
  tempo.Append(Version()) ; 
  SetName(tempo) ; 
  fRecParticlesInRun = 0 ; 

  Init() ;

}

//____________________________________________________________________________
AliPHOSPIDv2::~AliPHOSPIDv2()
{ 
  delete [] fX;
  delete [] fP; 
}


//____________________________________________________________________________
Float_t  AliPHOSPIDv2::GetDistance(AliPHOSEmcRecPoint * emc,AliPHOSRecPoint * cpv, Option_t *  Axis)const
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
 void  AliPHOSPIDv2::SetEllipseParameters(Float_t x, Float_t y,Float_t a, Float_t b,Float_t angle)
{
  fX_center = x ;
  fY_center = y ;
  fA = a  ;
  fB = b ;
  fAngle = angle ;
}

//____________________________________________________________________________
Int_t  AliPHOSPIDv2::GetPrincipalSign(Double_t* P )const
{
  //This method gives if the PCA of the particle are inside a defined ellipse
 
  Int_t      prinsign;
  Double_t   fDx        = 0. ; 
  Double_t   fDelta     = 0. ; 
  Double_t   fY         = 0. ; 
  Double_t   fY_1       = 0. ; 
  Double_t   fY_2       = 0. ;
  Double_t   fPi        = TMath::Pi() ;
  Double_t   fCos_Theta = TMath::Cos(fPi*fAngle/180.) ;
  Double_t   fSin_Theta = TMath::Sin(fPi*fAngle/180.) ;   

  fDx = P[0] - fX_center ; 
  fDelta = 4.*fA*fA*fA*fB* (fA*fA*fCos_Theta*fCos_Theta + fB*fB*fSin_Theta*fSin_Theta - fDx*fDx) ; 
  if (fDelta < 0.) 
    {prinsign=0;} 
  
  else if (fDelta == 0.) 
    { 
      fY = fCos_Theta*fSin_Theta*(fA*fA - fB*fB)*fDx / (fA*fA*fCos_Theta*fCos_Theta +
							fB*fB*fSin_Theta*fSin_Theta) ; 
      fY += fY_center ; 
      if(P[1]==fY ) 
	{prinsign=1;} 
      else 
	{prinsign=0;} 
    } 
  else 
    { 
      fY_1 = (fCos_Theta*fSin_Theta*(fA*fA - fB*fB) *fDx +
	     TMath::Sqrt(fDelta)/2.)/(fA*fA*fCos_Theta*fCos_Theta + fB*fB*fSin_Theta*fSin_Theta) ; 
      fY_2 = (fCos_Theta*fSin_Theta*(fA*fA - fB*fB) *fDx -
	     TMath::Sqrt(fDelta)/2.)/(fA*fA*fCos_Theta*fCos_Theta + fB*fB*fSin_Theta*fSin_Theta) ; 
      fY_1 += fY_center ; 
      fY_2 += fY_center ; 
      if ((P[1]<=fY_1) && (P[1]>=fY_2)) 
	{prinsign=1;} 
      else 
	{prinsign=0;}  
    } 
  return prinsign;
}
//____________________________________________________________________________
void  AliPHOSPIDv2::Exec(Option_t * option) 
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
    cerr << "WARNING: AliPHOSPIDv2::Exec -> RecParticles and/or PIDtMaker branch with name " 
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
void AliPHOSPIDv2::Init()
{
  // Make all memory allocations that are not possible in default constructor
  // Add the PID task to the list of PHOS tasks
  
  if ( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;
  
  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;

  fCpvEmcDistance = 4.0 ; 
  fTimeGate       = 0.162e-7 ;

  //PCA 
  fX              = new double[7]; // Data for the PCA 
  fP              = new double[7]; // Eigenvalues of the PCA  
  fFileName       = "$ALICE_ROOT/PHOS/PCA8pa15_0.5-100.root" ; 
  TFile f( fFileName.Data(), "read" ) ;
  fPrincipal      = dynamic_cast<TPrincipal*> (f.Get("principal")) ; 
 
  // Ellipse parameters 
  fX_center          = 2.0 ; 
  fY_center          = -0.35 ; 

  fA                 = 1.9 ; 
  fB                 = 1.0 ; 
  fAngle             = -60. ;
 
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(GetTitle(), taskName.Data()) ; 
  if ( gime == 0 ) {
    cerr << "ERROR: AliPHOSPIDv2::Init -> Could not obtain the Getter object !" << endl ; 
    return ;
  } 
   
  gime->PostPID(this) ;
  // create a folder on the white board //YSAlice/WhiteBoard/RecParticles/PHOS/recparticlesName
  gime->PostRecParticles(taskName.Data() ) ; 
  
}

//____________________________________________________________________________
void  AliPHOSPIDv2::MakeRecParticles(){

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
    
    // Looking at the CPV detector. If RCPV grater than fCpvEmcDistance, first bit 1.
    if(cpv)
      if(GetDistance(emc, cpv,  "R") > fCpvEmcDistance )  
	rp->SetPIDBit(0);
    
    //Looking the TOF. If TOF smaller than gate, second bit 1.
    if(emc->GetTime()< fTimeGate) 
      rp->SetPIDBit(1);				    
    
    
    //Loking PCA. Define and calculate the data (X) introduce in the function 
    //X2P that gives the components (P).  
    Float_t    fSpher = 0. ;
    Float_t    fEmaxdtotal = 0. ; 
    Float_t lambda[2] ;
    
    emc->GetElipsAxis(lambda) ;
    if((lambda[0]+lambda[1])!=0) fSpher=fabs(lambda[0]-lambda[1])/(lambda[0]+lambda[1]); 
    
    fEmaxdtotal=emc->GetMaximalEnergy()/emc->GetEnergy(); 
    
    fX[0] = lambda[0]; //cout <<" lambda0 "<<fX[0]; 
    fX[1] = lambda[1]; //cout <<" lambda1 "<<fX[1]; 
    fX[2] = emc->GetDispersion(); // cout <<" disper "<<fX[2]; 
    fX[3] = fSpher; // cout <<" spher "<<fX[3]; 
    fX[4] = emc->GetMultiplicity(); // cout <<" mult "<<fX[4]; 
    fX[5] = fEmaxdtotal; // cout <<" emax "<<fX[5]; 
    fX[6] = emc->GetCoreEnergy(); //  cout <<" core "<<fX[5]<< endl ; 
    
    // cout<<"Principal "<<fPrincipal<<endl;
    fPrincipal->X2P(fX,fP);
    
    //If we are inside the ellipse, third bit 1
    if(GetPrincipalSign(fP)== 1) 
      rp->SetPIDBit(2) ;
    
    
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
void  AliPHOSPIDv2:: Print(Option_t * option) const
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
    cout <<  "    Time Gate used:                  " << fTimeGate <<  endl ;
    cout <<  "    Principal Ellipse Parameters     " << endl ;
    cout <<  "       Ellipse center   (x,y)           (" << fX_center<<","<<fY_center<<")"<< endl;
    cout <<  "       Ellipse focus    (a,b)           (" << fA<<","<<fB<<")"<< endl;
    cout <<  "       Ellipse angle                     " << fAngle<< endl;        
    cout <<  "============================================" << endl ;
}

//____________________________________________________________________________
void  AliPHOSPIDv2::WriteRecParticles(Int_t event)
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
  AliPHOSPIDv2 * pid = this ;
  TBranch * pidBranch = gAlice->TreeR()->Branch("AliPHOSPID","AliPHOSPIDv2",&pid,bufferSize,splitlevel);
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
  //pidBranch->Write(0,kOverwrite) ;  

  delete [] filename ; 
}

//____________________________________________________________________________
TVector3 AliPHOSPIDv2::GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSRecPoint * cpv)const 
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
void AliPHOSPIDv2::PrintRecParticles(Option_t * option)
{
  // Print table of reconstructed particles

  AliPHOSGetter *gime = AliPHOSGetter::GetInstance() ; 

  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;
  TClonesArray * recParticles = gime->RecParticles(taskName) ; 
  
  cout << "AliPHOSPIDv2: event "<<gAlice->GetEvNumber()  << endl ;
  cout << "       found " << recParticles->GetEntriesFast() << " RecParticles " << endl ;
  
  if(strstr(option,"all")) {  // printing found TS
    
    cout << "  PARTICLE "   
	 << "  Index    "  << endl ;
    
    Int_t index ;
    for (index = 0 ; index < recParticles->GetEntries() ; index++) {
       AliPHOSRecParticle * rp = (AliPHOSRecParticle * ) recParticles->At(index) ;       
      
       Text_t particle[11];
       switch(rp->GetType()) {
       case  AliPHOSFastRecParticle::kCHARGEDHASLOW:
	 strcpy(particle, "CHARGED HA SLOW") ;
	 break ; 
       case  AliPHOSFastRecParticle::kNEUTRALHASLOW: 
	 strcpy(particle, "NEUTRAL HA SLOW");
	 break ;    
       case  AliPHOSFastRecParticle::kCHARGEDHAFAST:
	 strcpy(particle, "CHARGED HA FAST") ;
	 break ;	
       case  AliPHOSFastRecParticle::kNEUTRALHAFAST:
	 strcpy(particle, "NEUTRAL HA FAST");
	 break;
       case  AliPHOSFastRecParticle::kCHARGEDEMSLOW:
	 strcpy(particle, "CHARGED EM SLOW") ;
	 break ;
       case  AliPHOSFastRecParticle::kNEUTRALEMSLOW:
	 strcpy(particle, "NEUTRAL EM SLOW");
	 break ;
       case  AliPHOSFastRecParticle::kCHARGEDEMFAST:
	 strcpy(particle, "CHARGED EM FAST") ;
	 break ;
       case  AliPHOSFastRecParticle::kNEUTRALEMFAST:
	 strcpy( particle, "NEUTRAL EM FAST");
	 break;
      }

      cout << setw(10) << particle << "  "
	   << setw(5) <<  rp->GetIndexInList() << " " <<endl;
   
    }
    cout << "-------------------------------------------" << endl ;
  }
  
}



