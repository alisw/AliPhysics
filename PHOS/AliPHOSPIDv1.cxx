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
//     - RCPV: distance from CPV recpoint to EMCA recpoint.
//     - TOF 
//     - PCA: Principal Components Analysis..
// The identified particle has an identification number corresponding 
// to a 9 bits number:
//     -Bit 0 to 2: bit set if RCPV > CpvEmcDistance (each bit corresponds
//      to a different efficiency-purity point of the photon identification) 
//     -Bit 3 to 5: bit set if TOF  < TimeGate (each bit corresponds
//      to a different efficiency-purity point of the photon identification) 
//     -Bit 6 to 9: bit set if Principal Components are 
//      inside an ellipse defined by the parameters a, b, c, x0 and y0.
//      (each bit corresponds to a different efficiency-purity point of the 
//      photon identification)
//      The PCA (Principal components analysis) needs a file that contains
//      a previous analysis of the correlations between the particles. This 
//      file is $ALICE_ROOT/PHOS/PCA8pa15_0.5-100.root. Analysis done for 
//      energies between 0.5 and 100 GeV.
//      A calibrated energy is calculated. The energy of the reconstructed
//      cluster is corrected with the formula A + B * E  + C * E^2, whose 
//      parameters where obtained through the study of the reconstructed 
//      energy distribution of monoenergetic photons. 
//
//      All the parameters (RCPV(2 rows-3 columns),TOF(1r-3c),PCA(5r-4c) 
//      and calibration(1r-3c))are stored in a file called 
//      $ALICE_ROOT/PHOS/Parameters.dat. Each time that AliPHOSPIDv1 is 
//      initialized, this parameters are copied to a Matrix (9,4), a 
//      TMatrixD object.  
//
// use case:
//  root [0] AliPHOSPIDv1 * p = new AliPHOSPIDv1("galice1.root")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//          // reading headers from file galice1.root and create  RecParticles 
            // TrackSegments and RecPoints are used 
//          // set file name for the branch RecParticles
//  root [1] p->ExecuteTask("deb all time")
//          // available options
//          // "deb" - prints # of reconstructed particles
//          // "deb all" -  prints # and list of RecParticles
//          // "time" - prints benchmarking results
//                  
//  root [2] AliPHOSPIDv1 * p2 = new AliPHOSPIDv1("galice1.root","v1",kTRUE)
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                //Split mode.  
//  root [3] p2->ExecuteTask()
//


//*-- Author: Yves Schutz (SUBATECH)  & Gines Martinez (SUBATECH) & 
//            Gustavo Conesa April 2002
//            PCA redesigned by Gustavo Conesa October 2002:
//            The way of using the PCA has changed. Instead of 2
//            files with the PCA, each one with different energy ranges 
//            of application, we use the wide one (0.5-100 GeV), and instead
//            of fixing 3 ellipses for different ranges of energy, it has been
//            studied the dependency of the ellipses parameters with the 
//            energy, and they are implemented in the code as a funtion 
//            of the energy. 
//
//
//
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
#include "TMatrixD.h"
#include "TPrincipal.h"
#include "TSystem.h"

// --- Standard library ---


// --- AliRoot header files ---

#include "AliGenerator.h"
#include "AliPHOS.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSEmcRecPoint.h"
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
 
  InitParameters() ; 
  fDefaultInit = kTRUE ; 
}

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1(const AliPHOSPIDv1 & pid ):AliPHOSPID(pid)
{ 
  // ctor
  InitParameters() ; 
  Init() ;

}

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1(const TString alirunFileName, const TString eventFolderName):AliPHOSPID(alirunFileName, eventFolderName)
{ 
  //ctor with the indication on where to look for the track segments
 
  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 
}

//____________________________________________________________________________
AliPHOSPIDv1::~AliPHOSPIDv1()
{ 
  // dtor

  delete [] fX ;       // Principal input 
  delete [] fPPhoton ; // Photon Principal components
  delete [] fPPi0 ;    // Pi0 Principal components
}
//____________________________________________________________________________
const TString AliPHOSPIDv1::BranchName() const 
{  

  return GetName() ;
}
 
//____________________________________________________________________________
void AliPHOSPIDv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  // Add the PID task to the list of PHOS tasks

  AliPHOSGetter * gime = AliPHOSGetter::Instance(GetTitle(), fEventFolderName.Data()) ; 

  if ( !gime->PID() ) 
    gime->PostPID(this) ;
}

//____________________________________________________________________________
void AliPHOSPIDv1::InitParameters()
{
  // Initialize PID parameters
  fRecParticlesInRun = 0 ; 
  fNEvent            = 0 ;            
  fRecParticlesInRun = 0 ;
  SetParameters() ; // fill the parameters matrix from parameters file
  SetEventRange(0,-1) ;
}

//________________________________________________________________________
void  AliPHOSPIDv1::Exec(Option_t *option)
{
  // Steering method to perform particle reconstruction and identification
  // for the event range from fFirstEvent to fLastEvent.
  // This range is optionally set by SetEventRange().
  // if fLastEvent=-1 (by default), then process events until the end.
  
  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSPID");
  
  if(strstr(option,"print")) {
    Print() ; 
    return ; 
  }


  AliPHOSGetter * gime = AliPHOSGetter::Instance(GetTitle()) ; 
 
  if (fLastEvent == -1) 
    fLastEvent = gime->MaxEvent() - 1 ;
  else 
    fLastEvent = TMath::Min(fLastEvent,gime->MaxEvent());
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;

  Int_t ievent ; 
  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    gime->Event(ievent,"TR") ;
    if(gime->TrackSegments() && //Skip events, where no track segments made
       gime->TrackSegments()->GetEntriesFast()) {
      MakeRecParticles() ;
      WriteRecParticles();
      if(strstr(option,"deb"))
	PrintRecParticles(option) ;
      //increment the total number of rec particles per run 
      fRecParticlesInRun += gime->RecParticles()->GetEntriesFast() ; 
    }
  }
  if(strstr(option,"deb"))
      PrintRecParticles(option);
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSPID");
    Info("Exec", "took %f seconds for PID %f seconds per event", 
	 gBenchmark->GetCpuTime("PHOSPID"),  
	 gBenchmark->GetCpuTime("PHOSPID")/nEvents) ;
  }
  Unload();
}

//____________________________________________________________________________
const TString AliPHOSPIDv1::GetFileNamePrincipal(TString particle) const
{
  //Get file name that contains the PCA for a particle ("photon or pi0")
  particle.ToLower();
  TString name;
  if      (particle=="photon") name = fFileNamePrincipalPhoton ;
  else if (particle=="pi0"   ) name = fFileNamePrincipalPi0    ;
  else    Error("GetFileNamePrincipal","Wrong particle name: %s (choose from pi0/photon)\n",particle.Data());
  return name;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterCalibration(Int_t i) const 
{
  // Get the i-th parameter "Calibration"
  Float_t param = 0.;
  if (i>2 || i<0)
    Error("GetParameterCalibration","Invalid parameter number: %d",i);
  else
    param = (*fParameters)(0,i);
  return param;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetCalibratedEnergy(Float_t e) const
{
//      It calibrates Energy depending on the recpoint energy.
//      The energy of the reconstructed cluster is corrected with 
//      the formula A + B* E  + C* E^2, whose parameters where obtained 
//      through the study of the reconstructed energy distribution of 
//      monoenergetic photons.
 
  Float_t p[]={0.,0.,0.};
  for (Int_t i=0; i<3; i++) p[i] = GetParameterCalibration(i);
  Float_t enerec = p[0] +  p[1]*e + p[2]*e*e;
  return enerec ;

}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterCpv2Emc(Int_t i, TString axis) const 
{
  // Get the i-th parameter "CPV-EMC distance" for the specified axis
  Float_t param = 0.;
  if(i>2 || i<0)
    Error("GetParameterCpv2Emc","Invalid parameter number: %d",i);
  else {
    axis.ToLower();
    if      (axis == "x") param = (*fParameters)(1,i);
    else if (axis == "z") param = (*fParameters)(2,i);
    else Error("GetParameterCpv2Emc","Invalid axis name: %s",axis.Data());
  }
  return  param;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetCpv2EmcDistanceCut(TString axis, Float_t e) const
{
  // Get CpvtoEmcDistance Cut depending on the cluster energy, axis and 
  // Purity-Efficiency point 

  axis.ToLower();
  Float_t p[]={0.,0.,0.};
  for (Int_t i=0; i<3; i++) p[i] = GetParameterCpv2Emc(i,axis);
  Float_t sig = p[0] + TMath::Exp(p[1] - p[2]*e);
  return sig;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetEllipseParameter(TString particle, TString param, Float_t e) const 
{
  // Calculates the parameter param of the ellipse

  particle.ToLower();
  param.   ToLower();
  Float_t p[4]={0.,0.,0.,0.};
  Float_t value = 0.0;
  for (Int_t i=0; i<4; i++) p[i] = GetParameterToCalculateEllipse(particle,param,i);
  if (particle == "photon") {
    if      (param.Contains("a"))  e = TMath::Min((Double_t)e,70.);
    else if (param.Contains("b"))  e = TMath::Min((Double_t)e,70.);
    else if (param.Contains("x0")) e = TMath::Max((Double_t)e,1.1);
  }

  value = p[0]/TMath::Sqrt(e) + p[1]*e + p[2]*e*e + p[3];
  return value;
}

//_____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterPhotonBoundary (Int_t i) const
{ 
  // Get the parameter "i" to calculate the boundary on the moment M2x
  // for photons at high p_T
  Float_t param = 0;
  if (i>3 || i<0)
    Error("GetParameterPhotonBoundary","Wrong parameter number: %d\n",i);
  else
    param = (*fParameters)(14,i) ;
  return param;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterPi0Boundary (Int_t i) const
{ 
  // Get the parameter "i" to calculate the boundary on the moment M2x
  // for pi0 at high p_T
  Float_t param = 0;
  if (i>2 || i<0)
    Error("GetParameterPi0Boundary","Wrong parameter number: %d\n",i);
  else
    param = (*fParameters)(15,i) ;
  return param;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterTimeGate(Int_t i) const
{
  // Get TimeGate parameter depending on Purity-Efficiency i:
  // i=0 - Low purity, i=1 - Medium purity, i=2 - High purity
  Float_t param = 0.;
  if(i>2 || i<0)
    Error("GetParameterTimeGate","Invalid Efficiency-Purity choice %d",i);
  else
    param = (*fParameters)(3,i) ; 
  return param;
}

//_____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetParameterToCalculateEllipse(TString particle, TString param, Int_t i) const
{ 
  // Get the parameter "i" that is needed to calculate the ellipse 
  // parameter "param" for the particle "particle" ("photon" or "pi0")

  particle.ToLower();
  param.   ToLower();
  Int_t offset = -1;
  if      (particle == "photon") offset=0;
  else if (particle == "pi0")    offset=5;
  else
    Error("GetParameterToCalculateEllipse","Wrong particle name: %s (choose from pi0/photon)\n",particle.Data());

  Int_t p= -1;
  Float_t par = 0;

  if     (param.Contains("a")) p=4+offset; 
  else if(param.Contains("b")) p=5+offset; 
  else if(param.Contains("c")) p=6+offset; 
  else if(param.Contains("x0"))p=7+offset; 
  else if(param.Contains("y0"))p=8+offset;

  if      (i>4 || i<0)
    Error("GetParameterToCalculateEllipse", "No parameter with index", i) ; 
  else if (p==-1)
    Error("GetParameterToCalculateEllipse", "No parameter with name %s", param.Data() ) ; 
  else
    par = (*fParameters)(p,i) ;
  
  return par;
}


//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetDistance(AliPHOSEmcRecPoint * emc,AliPHOSCpvRecPoint * cpv, Option_t *  axis)const
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
  
  const AliPHOSGeometry * geom = AliPHOSGetter::Instance()->PHOSGeometry() ; 
  TVector3 vecEmc ;
  TVector3 vecCpv ;
  if(cpv){
    emc->GetLocalPosition(vecEmc) ;
    cpv->GetLocalPosition(vecCpv) ; 
    
    if(emc->GetPHOSMod() == cpv->GetPHOSMod()){      
      // Correct to difference in CPV and EMC position due to different distance to center.
      // we assume, that particle moves from center
      Float_t dCPV = geom->GetIPtoOuterCoverDistance();
      Float_t dEMC = geom->GetIPtoCrystalSurface() ;
      dEMC         = dEMC / dCPV ;
      vecCpv = dEMC * vecCpv  - vecEmc ; 
      if (axis == "X") return vecCpv.X();
      if (axis == "Y") return vecCpv.Y();
      if (axis == "Z") return vecCpv.Z();
      if (axis == "R") return vecCpv.Mag();
    }
    return 100000000 ;
  }
  return 100000000 ;
}
//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetCPVBit(AliPHOSEmcRecPoint * emc,AliPHOSCpvRecPoint * cpv, Int_t effPur, Float_t e) const
{
  if(effPur>2 || effPur<0)
    Error("GetCPVBit","Invalid Efficiency-Purity choice %d",effPur);
  
  Float_t sigX = GetCpv2EmcDistanceCut("X",e);
  Float_t sigZ = GetCpv2EmcDistanceCut("Z",e);
  
  Float_t deltaX = TMath::Abs(GetDistance(emc, cpv,  "X"));
  Float_t deltaZ = TMath::Abs(GetDistance(emc, cpv,  "Z"));

  if((deltaX>sigX*(effPur+1))|(deltaZ>sigZ*(effPur+1)))
    return 1;//Neutral
  else
    return 0;//Charged
}

//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetPrincipalBit(TString particle, const Double_t* p, Int_t effPur, Float_t e)const
{
  //Is the particle inside de PCA ellipse?
  
  particle.ToLower();
  Int_t    prinbit  = 0 ;
  Float_t a        = GetEllipseParameter(particle,"a" , e); 
  Float_t b        = GetEllipseParameter(particle,"b" , e);
  Float_t c        = GetEllipseParameter(particle,"c" , e);
  Float_t x0 = GetEllipseParameter(particle,"x0", e); 
  Float_t y0 = GetEllipseParameter(particle,"y0", e);
  
  Float_t r = TMath::Power((p[0] - x0)/a,2) + 
              TMath::Power((p[1] - y0)/b,2) +
            c*(p[0] -  x0)*(p[1] - y0)/(a*b) ;
  //3 different ellipses defined
  if((effPur==2) && (r<1./2.)) prinbit= 1;
  if((effPur==1) && (r<2.   )) prinbit= 1;
  if((effPur==0) && (r<9./2.)) prinbit= 1;

  if(r<0)
    Error("GetPrincipalBit", "Negative square?") ;

  return prinbit;

}
//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetHardPhotonBit(AliPHOSEmcRecPoint * emc) const
{
  // Set bit for identified hard photons (E > 30 GeV)
  // if the second moment M2x is below the boundary

  Float_t e   = emc->GetEnergy();
  if (e < 30.0) return 0;
  Float_t m2x = emc->GetM2x();
  Float_t m2xBoundary = GetParameterPhotonBoundary(0) *
    TMath::Exp(-TMath::Power(e-GetParameterPhotonBoundary(1),2)/2.0/
	        TMath::Power(GetParameterPhotonBoundary(2),2)) +
    GetParameterPhotonBoundary(3);
  //Info("GetHardPhotonBit","E=%f, m2x=%f, boundary=%f",e,m2x,m2xBoundary);
  if (m2x < m2xBoundary)
    return 1;// A hard photon
  else
    return 0;// Not a hard photon
}

//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetHardPi0Bit(AliPHOSEmcRecPoint * emc) const
{
  // Set bit for identified hard pi0  (E > 30 GeV)
  // if the second moment M2x is above the boundary

  Float_t e   = emc->GetEnergy();
  if (e < 30.0) return 0;
  Float_t m2x = emc->GetM2x();
  Float_t m2xBoundary = GetParameterPi0Boundary(0) +
                    e * GetParameterPi0Boundary(1);
  //Info("GetHardPi0Bit","E=%f, m2x=%f, boundary=%f",e,m2x,m2xBoundary);
  if (m2x > m2xBoundary)
    return 1;// A hard pi0
  else
    return 0;// Not a hard pi0
}

//____________________________________________________________________________
TVector3 AliPHOSPIDv1::GetMomentumDirection(AliPHOSEmcRecPoint * emc, AliPHOSCpvRecPoint * )const 
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
  dir.SetMag(1.) ;

  //account correction to the position of IP
  Float_t xo,yo,zo ; //Coordinates of the origin
  gAlice->Generator()->GetOrigin(xo,yo,zo) ;
  TVector3 origin(xo,yo,zo);
  dir = dir - origin ;

  return dir ;  
}

//____________________________________________________________________________
void  AliPHOSPIDv1::MakeRecParticles()
{
  // Makes a RecParticle out of a TrackSegment
  
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments() ; 
  if ( !emcRecPoints || !cpvRecPoints || !trackSegments ) {
    Fatal("MakeRecParticles", "RecPoints or TrackSegments not found !") ;  
  }
  TClonesArray * recParticles  = gime->RecParticles() ; 
  recParticles->Clear();

  TIter next(trackSegments) ; 
  AliPHOSTrackSegment * ts ; 
  Int_t index = 0 ; 
  AliPHOSRecParticle * rp ; 
  while ( (ts = (AliPHOSTrackSegment *)next()) ) {
    
    new( (*recParticles)[index] ) AliPHOSRecParticle() ;
    rp = (AliPHOSRecParticle *)recParticles->At(index) ; 
    rp->SetTrackSegment(index) ;
    rp->SetIndexInList(index) ;
    	
    AliPHOSEmcRecPoint * emc = 0 ;
    if(ts->GetEmcIndex()>=0)
      emc = (AliPHOSEmcRecPoint *) emcRecPoints->At(ts->GetEmcIndex()) ;
    
    AliPHOSCpvRecPoint * cpv = 0 ;
    if(ts->GetCpvIndex()>=0)
      cpv = (AliPHOSCpvRecPoint *) cpvRecPoints->At(ts->GetCpvIndex()) ;
    
    Int_t track = 0 ; 
    track = ts->GetTrackIndex() ; 
      
    // Now set type (reconstructed) of the particle

    // Choose the cluster energy range
    
    if (!emc) {
      Fatal("MakeRecParticles", "-> emc(%d) = %d", ts->GetEmcIndex(), emc ) ;
    }

    Float_t e = emc->GetEnergy() ;   
    
    Float_t  lambda[2] ;
    emc->GetElipsAxis(lambda) ;
    
    if((lambda[0]>0.01) && (lambda[1]>0.01)){
      // Looking PCA. Define and calculate the data (X),
      // introduce in the function X2P that gives the components (P).  

      Float_t  Spher = 0. ;
      Float_t  Emaxdtotal = 0. ; 
      
      if((lambda[0]+lambda[1])!=0) 
	Spher=fabs(lambda[0]-lambda[1])/(lambda[0]+lambda[1]); 
      
      Emaxdtotal=emc->GetMaximalEnergy()/emc->GetEnergy(); 
      
      fX[0] = lambda[0] ;  
      fX[1] = lambda[1] ; 
      fX[2] = emc->GetDispersion() ; 
      fX[3] = Spher ; 
      fX[4] = emc->GetMultiplicity() ;  
      fX[5] = Emaxdtotal ;  
      fX[6] = emc->GetCoreEnergy() ;  
      
      fPrincipalPhoton->X2P(fX,fPPhoton);
      fPrincipalPi0   ->X2P(fX,fPPi0);

    }
    else{
      fPPhoton[0]=-100.0;  //We do not accept clusters with 
      fPPhoton[1]=-100.0;  //one cell as a photon-like
      fPPi0[0]   =-100.0;
      fPPi0[1]   =-100.0;
    }
    
    Float_t time =emc->GetTime() ;
    
    // Loop of Efficiency-Purity (the 3 points of purity or efficiency 
    // are taken into account to set the particle identification)
    for(Int_t effPur = 0; effPur < 3 ; effPur++){
      
      // Looking at the CPV detector. If RCPV greater than CpvEmcDistance, 
      // 1st,2nd or 3rd bit (depending on the efficiency-purity point )
      // is set to 1
      if(GetCPVBit(emc, cpv, effPur,e) == 1 )  
	rp->SetPIDBit(effPur) ;
      
      // Looking the TOF. If TOF smaller than gate,  4th, 5th or 6th 
      // bit (depending on the efficiency-purity point )is set to 1             
      if(time< (*fParameters)(2,effPur)) 
	rp->SetPIDBit(effPur+3) ;		    
      
      //Photon PCA
      //If we are inside the ellipse, 7th, 8th or 9th 
      // bit (depending on the efficiency-purity point )is set to 1 
      if(GetPrincipalBit("photon",fPPhoton,effPur,e) == 1) 
	rp->SetPIDBit(effPur+6) ;

      //Pi0 PCA
      //If we are inside the ellipse, 10th, 11th or 12th 
      // bit (depending on the efficiency-purity point )is set to 1 
      if(GetPrincipalBit("pi0"   ,fPPi0   ,effPur,e) == 1) 
	rp->SetPIDBit(effPur+9) ;
    }
    if(GetHardPhotonBit(emc))
      rp->SetPIDBit(12) ;
    if(GetHardPi0Bit   (emc))
      rp->SetPIDBit(13) ;
    
    if(track >= 0) 
      rp->SetPIDBit(14) ; 

    //Set momentum, energy and other parameters 
    Float_t  encal = GetCalibratedEnergy(e);
    TVector3 dir   = GetMomentumDirection(emc,cpv) ; 
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
    //Set the position in global coordinate system from the RecPoint
    AliPHOSGeometry * geom = gime->PHOSGeometry() ; 
    AliPHOSTrackSegment * ts  = gime->TrackSegment(rp->GetPHOSTSIndex()) ; 
    AliPHOSEmcRecPoint  * erp = gime->EmcRecPoint(ts->GetEmcIndex()) ; 
    TVector3 pos ; 
    geom->GetGlobal(erp, pos) ; 
    rp->SetPos(pos);

    index++ ; 
  }
}
  
//____________________________________________________________________________
void  AliPHOSPIDv1::Print() const
{
  // Print the parameters used for the particle type identification

    Info("Print", "=============== AliPHOSPIDv1 ================") ;
    printf("Making PID\n") ;
    printf("    Pricipal analysis file from 0.5 to 100 %s\n", fFileNamePrincipalPhoton.Data() )   ; 
    printf("    Name of parameters file     %s\n", fFileNameParameters.Data() )  ;
    printf("    Matrix of Parameters: 14x4\n") ;
    printf("        Energy Calibration  1x3 [3 parametres to calibrate energy: A + B* E + C * E^2]\n") ;
    printf("        RCPV 2x3 rows x and z, columns function cut parameters\n") ;
    printf("        TOF  1x3 [High Eff-Low Pur,Medium Eff-Pur, Low Eff-High Pur]\n") ;
    printf("        PCA  5x4 [5 ellipse parametres and 4 parametres to calculate them: A/Sqrt(E) + B* E + C * E^2 + D]\n") ;
    Printf("    Pi0 PCA  5x3 [5 ellipse parametres and 3 parametres to calculate them: A + B* E + C * E^2]\n") ;
    fParameters->Print() ;
}



//____________________________________________________________________________
void AliPHOSPIDv1::PrintRecParticles(Option_t * option)
{
  // Print table of reconstructed particles

  AliPHOSGetter *gime = AliPHOSGetter::Instance() ; 

  TClonesArray * recParticles = gime->RecParticles() ; 

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
      AliPHOSRecParticle * rp = (AliPHOSRecParticle * ) recParticles->At(index) ;       
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

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameters() 
{
  // PCA : To do the Principal Components Analysis it is necessary 
  // the Principal file, which is opened here
  fX       = new double[7]; // Data for the PCA 
  fPPhoton = new double[7]; // Eigenvalues of the PCA
  fPPi0    = new double[7]; // Eigenvalues of the Pi0 PCA

  // Read photon principals from the photon file
  
  fFileNamePrincipalPhoton = "$ALICE_ROOT/PHOS/PCA8pa15_0.5-100.root" ; 
  TFile f( fFileNamePrincipalPhoton.Data(), "read" ) ;
  fPrincipalPhoton = dynamic_cast<TPrincipal*> (f.Get("principal")) ; 
  f.Close() ; 

  // Read pi0 principals from the pi0 file

  fFileNamePrincipalPi0    = "$ALICE_ROOT/PHOS/PCA_pi0_40-120.root" ;
  TFile fPi0( fFileNamePrincipalPi0.Data(), "read" ) ;
  fPrincipalPi0    = dynamic_cast<TPrincipal*> (fPi0.Get("principal")) ; 
  fPi0.Close() ;

  // Open parameters file and initialization of the Parameters matrix. 
  // In the File Parameters.dat are all the parameters. These are introduced 
  // in a matrix of 16x4  
  // 
  // All the parameters defined in this file are, in order of row: 
  // line   0   : calibration 
  // lines  1,2 : CPV rectangular cat for X and Z
  // line   3   : TOF cut
  // lines  4-8 : parameters to calculate photon PCA ellipse
  // lines  9-13: parameters to calculate pi0 PCA ellipse
  // lines 14-15: parameters to calculate border for high-pt photons and pi0

  fFileNameParameters = gSystem->ExpandPathName("$ALICE_ROOT/PHOS/Parameters.dat");
  fParameters = new TMatrix(16,4) ;
  const Int_t maxLeng=255;
  char string[maxLeng];

  // Open a text file with PID parameters
  FILE *fd = fopen(fFileNameParameters.Data(),"r");
  if (!fd)
    Fatal("SetParameter","File %s with a PID parameters cannot be opened\n",
	  fFileNameParameters.Data());

  Int_t i=0;
  // Read parameter file line-by-line and skip empty line and comments
  while (fgets(string,maxLeng,fd) != NULL) {
    if (string[0] == '\n' ) continue;
    if (string[0] == '!'  ) continue;
    sscanf(string, "%f %f %f %f",
	   &(*fParameters)(i,0), &(*fParameters)(i,1), 
	   &(*fParameters)(i,2), &(*fParameters)(i,3));
    i++;
    //printf("line %d: %s",i,string);
  }
  fclose(fd);
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterCalibration(Int_t i,Float_t param) 
{
  // Set parameter "Calibration" i to a value param
  if(i>2 || i<0)
    Error("SetParameterCalibration","Invalid parameter number: %d",i);
  else
    (*fParameters)(0,i) = param ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterCpv2Emc(Int_t i, TString axis, Float_t cut) 
{
  // Set the parameters to calculate Cpv-to-Emc Distance Cut depending on 
  // Purity-Efficiency point i

  if(i>2 || i<0)
    Error("SetParameterCpv2Emc","Invalid parameter number: %d",i);
  else {
    axis.ToLower();
    if      (axis == "x") (*fParameters)(1,i) = cut;
    else if (axis == "z") (*fParameters)(2,i) = cut;
    else Error("SetParameterCpv2Emc","Invalid axis name: %s",axis.Data());
  }
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterPhotonBoundary(Int_t i,Float_t param) 
{
  // Set parameter "Hard photon boundary" i to a value param
  if(i>4 || i<0)
    Error("SetParameterPhotonBoundary","Invalid parameter number: %d",i);
  else
    (*fParameters)(14,i) = param ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterPi0Boundary(Int_t i,Float_t param) 
{
  // Set parameter "Hard pi0 boundary" i to a value param
  if(i>1 || i<0)
    Error("SetParameterPi0Boundary","Invalid parameter number: %d",i);
  else
    (*fParameters)(15,i) = param ;
}

//_____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterTimeGate(Int_t i, Float_t gate) 
{
  // Set the parameter TimeGate depending on Purity-Efficiency point i 
  if (i>2 || i<0)
    Error("SetParameterTimeGate","Invalid Efficiency-Purity choice %d",i);
  else
    (*fParameters)(3,i)= gate ; 
} 

//_____________________________________________________________________________
void  AliPHOSPIDv1::SetParameterToCalculateEllipse(TString particle, TString param, Int_t i, Float_t par) 
{  
  // Set the parameter "i" that is needed to calculate the ellipse 
  // parameter "param" for a particle "particle"
  
  particle.ToLower();
  param.   ToLower();
  Int_t p= -1;
  Int_t offset=0;

  if      (particle == "photon") offset=0;
  else if (particle == "pi0")    offset=5;
  else
    Error("SetParameterToCalculateEllipse","Wrong particle name: %s (choose from pi0/photon)\n",particle.Data());

  if     (param.Contains("a")) p=4+offset; 
  else if(param.Contains("b")) p=5+offset; 
  else if(param.Contains("c")) p=6+offset; 
  else if(param.Contains("x0"))p=7+offset; 
  else if(param.Contains("y0"))p=8+offset;
  if((i>4)||(i<0))
    Error("SetEllipseParameter", "No parameter with index %d", i) ; 
  else if(p==-1)
    Error("SetEllipseParameter", "No parameter with name %s", param.Data() ) ; 
  else
    (*fParameters)(p,i) = par ;
} 

//____________________________________________________________________________
void AliPHOSPIDv1::Unload() 
{
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ;  
  gime->PhosLoader()->UnloadRecPoints() ;
  gime->PhosLoader()->UnloadTracks() ;
  gime->PhosLoader()->UnloadRecParticles() ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::WriteRecParticles()
{
 
  AliPHOSGetter *gime = AliPHOSGetter::Instance() ; 

  TClonesArray * recParticles = gime->RecParticles() ; 
  recParticles->Expand(recParticles->GetEntriesFast() ) ;
  TTree * treeP =  gime->TreeP();
  
  //First rp
  Int_t bufferSize = 32000 ;
  TBranch * rpBranch = treeP->Branch("PHOSRP",&recParticles,bufferSize);
  rpBranch->SetTitle(BranchName());
  
  rpBranch->Fill() ;
 
  gime->WriteRecParticles("OVERWRITE");
  gime->WritePID("OVERWRITE");
}

