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

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  if(!gime)
    gime = AliPHOSGetter::Instance(GetTitle(), fEventFolderName.Data()) ; 

  if ( !gime->PID() ) 
    gime->PostPID(this) ;
}

//____________________________________________________________________________
void AliPHOSPIDv1::InitParameters()
{
  // Initialize PID parameters
  fWrite                   = kTRUE ;
  fRecParticlesInRun = 0 ; 
  fNEvent            = 0 ;            
  fRecParticlesInRun = 0 ;
  fBayesian          = kTRUE ;
  SetParameters() ; // fill the parameters matrix from parameters file
  SetEventRange(0,-1) ;

  // initialisation of response function parameters
  // Tof
  // Photons
  fTphoton[0] = 0.218    ;
  //fTphoton[0] = 1.    ;
  fTphoton[1] = 1.55E-8  ; 
  fTphoton[2] = 5.05E-10 ;
  fTFphoton = new TFormula("ToF response to photons" , "gaus") ; 
  fTFphoton->SetParameters( fTphoton[0], fTphoton[1], fTphoton[2]) ; 
//   // Electrons
//   fTelectron[0] = 0.2      ;
//   fTelectron[1] = 1.55E-8  ; 
//   fTelectron[2] = 5.35E-10 ;
//   fTFelectron = new TFormula("ToF response to electrons" , "gaus") ; 
//   fTFelectron->SetParameters( fTelectron[0], fTelectron[1], fTelectron[2]) ; 
//   // Muons
//   fTmuon[0] = 0.2     ;
//   fTmuon[1] = 1.55E-8 ; 
//   fTmuon[2] = 5.1E-10 ;
//   fTFmuon = new TFormula("ToF response to muons" , "gaus") ; 
//   fTFmuon->SetParameters( fTmuon[0], fTmuon[1], fTmuon[2]) ; 

  // Pions
  //Gaus (0 to max probability)
  fTpiong[0] = 0.0971    ; 
  //fTpiong[0] = 1.       ;
  fTpiong[1] = 1.58E-8  ; 
  fTpiong[2] = 5.69E-10 ;
  fTFpiong = new TFormula("ToF response to pions" , "gaus") ; 
  fTFpiong->SetParameters( fTpiong[0], fTpiong[1], fTpiong[2]) ; 
  // Landau (max probability to inf) 
//   fTpionl[0] = 0.05     ; 
//   //fTpionl[0] = 5.53     ; 
//   fTpionl[1] = 1.68E-8  ; 
//   fTpionl[2] = 5.38E-10 ;
//   fTFpionl = new TFormula("ToF response to pions" , "landau") ; 
//   fTFpionl->SetParameters( fTpionl[0], fTpionl[1], fTpionl[2]) ; 


  // Kaons
  //Gaus (0 to max probability)
  fTkaong[0] = 0.0542  ; 
  //fTkaong[0] = 1.      ;
  fTkaong[1] = 1.64E-8 ; 
  fTkaong[2] = 6.07-10 ;
  fTFkaong = new TFormula("ToF response to kaon" , "gaus") ; 
  fTFkaong->SetParameters( fTkaong[0], fTkaong[1], fTkaong[2]) ; 
  //Landau (max probability to inf) 
  fTkaonl[0] = 0.264   ;
  //fTkaonl[0] = 5.53     ;
  fTkaonl[1] = 1.68E-8  ; 
  fTkaonl[2] = 4.10E-10 ;
  fTFkaonl = new TFormula("ToF response to kaon" , "landau") ; 
  fTFkaonl->SetParameters( fTkaonl[0], fTkaonl[1], fTkaonl[2]) ; 

  //Heavy Hadrons
  //Gaus (0 to max probability)
  fThhadrong[0] = 0.0302   ;  
  //fThhadrong[0] = 1.       ; 
  fThhadrong[1] = 1.73E-8  ; 
  fThhadrong[2] = 9.52E-10 ;
  fTFhhadrong = new TFormula("ToF response to heavy hadrons" , "gaus") ; 
  fTFhhadrong->SetParameters( fThhadrong[0], fThhadrong[1], fThhadrong[2]) ; 
  //Landau (max probability to inf) 
  fThhadronl[0] = 0.139    ;  
  //fThhadronl[0] = 5.53      ;   
  fThhadronl[1] = 1.745E-8  ; 
  fThhadronl[2] = 1.00E-9  ;
  fTFhhadronl = new TFormula("ToF response to heavy hadrons" , "landau") ; 
  fTFhhadronl->SetParameters( fThhadronl[0], fThhadronl[1], fThhadronl[2]) ; 

/// /gaussian parametrization for pions
//   fTpion[0] = 3.93E-2 ;  fTpion[1] = 0.130   ; fTpion[2] =-6.37E-2 ;//constant
//   fTpion[3] = 1.65E-8 ;  fTpion[4] =-1.40E-9 ; fTpion[5] = 5.96E-10;//mean
//   fTpion[6] = 8.09E-10;  fTpion[7] =-4.65E-10; fTpion[8] = 1.50E-10;//sigma

// //landau parametrization for kaons
//   fTkaon[0] = 0.107   ;  fTkaon[1] = 0.166   ; fTkaon[2] = 0.243   ;//constant
//   fTkaon[3] = 1.80E-8 ;  fTkaon[4] =-2.96E-9 ; fTkaon[5] = 9.60E-10;//mean
//   fTkaon[6] = 1.37E-9 ;  fTkaon[7] =-1.80E-9 ; fTkaon[8] = 6.74E-10;//sigma

// //landau parametrization for nucleons
//   fThhadron[0] = 6.33E-2 ;  fThhadron[1] = 2.52E-2 ; fThhadron[2] = 2.16E-2 ;//constant
//   fThhadron[3] = 1.94E-8 ;  fThhadron[4] =-7.06E-10; fThhadron[5] =-4.69E-10;//mean
//   fThhadron[6] = 2.55E-9 ;  fThhadron[7] =-1.90E-9 ; fThhadron[8] = 5.41E-10;//sigma


  // Shower shape: dispersion gaussian parameters
  // Photons

  fDphoton[0] = 0.1    ;  fDphoton[1] = 0.      ; fDphoton[2] = 0.     ;//constant
  //fDphoton[0] = 1.0    ;  fDphoton[1] = 0.      ; fDphoton[2] = 0.     ;//constant
  fDphoton[3] = 1.55   ;  fDphoton[4] =-0.0863  ; fDphoton[5] = 0.287  ;//mean
  fDphoton[6] = 0.0451 ;  fDphoton[7] =-0.0803  ; fDphoton[8] = 0.314  ;//sigma

   fDpi0[0] = 0.0586  ;  fDpi0[1] = 1.06E-3 ; fDpi0[2] = 0.      ;//constant
   //fDpi0[0] = 1.0     ;  fDpi0[1] = 0.0     ; fDpi0[2] = 0.      ;//constant
  fDpi0[3] = 2.67    ;  fDpi0[4] =-2.00E-2 ; fDpi0[5] = 9.37E-5 ;//mean
  fDpi0[6] = 0.153   ;  fDpi0[7] = 9.34E-4 ; fDpi0[8] =-1.49E-5 ;//sigma
  //landau
//   fDhadron[0] = 0.007  ;  fDhadron[1] = 0.      ; fDhadron[2] = 0.    ;//constant
//   //fDhadron[0] = 5.53   ;  fDhadron[1] = 0.      ; fDhadron[2] = 0.    ;//constant
//   fDhadron[3] = 3.38   ;  fDhadron[4] = 0.0833  ; fDhadron[5] =-0.845 ;//mean
//   fDhadron[6] = 0.627  ;  fDhadron[7] = 0.012   ; fDhadron[8] =-0.170 ;//sigma

  fDhadron[0] =-5.10E-3 ;  fDhadron[1] =-5.35E-3 ; fDhadron[2] = 3.77E-2 ;//constant
  fDhadron[3] = 4.03    ;  fDhadron[4] = 0.292   ; fDhadron[5] =-1.50    ;//mean
  fDhadron[6] = 0.958   ;  fDhadron[7] = 0.117   ; fDhadron[8] =-0.598   ;//sigma
  // Muons
  fDmuon[0] = 0.0631     ;
  fDmuon[1] = 1.4 ; 
  fDmuon[2] = 0.0557 ;
  fDFmuon = new TFormula("Shower shape response to muons" , "landau") ; 
  fDFmuon->SetParameters( fDmuon[0], fDmuon[1], fDmuon[2]) ; 

  // CPV-EMC distance gaussian parameters

  fCPVelectron[0] = 0.0     ;  fCPVelectron[1] = 0.0160 ; fCPVelectron[2] = 0.    ;//constant
  //fCPVelectron[0] = 1.0     ;  fCPVelectron[1] = 0.     ; fCPVelectron[2] = 0.    ;//constant
  fCPVelectron[3] = 0.0682  ;  fCPVelectron[4] =-1.32   ; fCPVelectron[5] = 6.67  ;//mean
  fCPVelectron[6] = 0.276   ;  fCPVelectron[7] = 0.234  ; fCPVelectron[8] = 0.356 ;//sigma

  //all charged landau
 //  fCPVcharged[0] = 0.0     ;  fCPVcharged[1] = 0.0464  ; fCPVcharged[2] = 0.      ;//constant
//   //fCPVcharged[0] = 5.53    ;  fCPVcharged[1] = 0.      ; fCPVcharged[2] = 0.      ;//constant
//   fCPVcharged[3] = 0.297   ;  fCPVcharged[4] =-0.652   ; fCPVcharged[5] = 1.91    ;//mean
//   fCPVcharged[6] = 0.0786  ;  fCPVcharged[7] =-0.237   ; fCPVcharged[8] = 0.752   ;//sigma

// //charged hadrons landau
//   fCPVchargedl[0] = 0.103   ;  fCPVchargedl[1] = 8.84E-3 ; fCPVchargedl[2] =-2.40E-2 ;//constant
//   fCPVchargedl[3] = 2.86    ;  fCPVchargedl[4] =-0.214   ; fCPVchargedl[5] = 0.817   ;//mean
//   fCPVchargedl[6] = 0.182   ;  fCPVchargedl[7] =-0.0693  ; fCPVchargedl[8] = 0.319   ;//sigma
//   //charged hadrons gaus
//   fCPVchargedg[0] = 0.0135  ;  fCPVchargedg[1] = 2.43E-5 ; fCPVchargedg[2] = 3.01E-3 ;//constant
//   fCPVchargedg[3] = 2.37    ;  fCPVchargedg[4] =-0.181   ; fCPVchargedg[5] = 0.726   ;//mean
//   fCPVchargedg[6] = 0.908   ;  fCPVchargedg[7] =-0.0558  ; fCPVchargedg[8] = 0.219   ;//sigma


  //charged hadrons landau
  fCPVcharged[0] = 6.06E-2 ;  fCPVcharged[1] = 3.80E-3 ; fCPVcharged[2] =-1.40E-2 ;//constant
  fCPVcharged[3] = 1.15    ;  fCPVcharged[4] =-0.563   ; fCPVcharged[5] = 2.63    ;//mean
  fCPVcharged[6] = 0.915   ;  fCPVcharged[7] =-0.0790  ; fCPVcharged[8] = 0.307   ;//sigma

  for (Int_t i =0; i<  AliESDtrack::kSPECIESN ; i++)
    fInitPID[i] = 1.;
  
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


  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
 
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

      if(fBayesian)
	MakePID() ; 
      
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
  if(fWrite)
    Unload();
}

//________________________________________________________________________
Double_t  AliPHOSPIDv1::GausF(Double_t  x, Double_t  y, Double_t * par)
{
  Double_t cnt    = par[2] * (x*x) + par[1] * x + par[0] ;
  Double_t mean   = par[4] / (x*x) + par[5] / x + par[3] ;
  Double_t sigma  = par[7] / (x*x) + par[8] / x + par[6] ;
  //cout<<"c "<< cnt << " mean "<<mean<<" sigma "<<sigma<<endl;
  //  Double_t arg    = - (y-mean) * (y-mean) / (2*sigma*sigma) ;
  //  return cnt * TMath::Exp(arg) ;
  if(mean != 0. && sigma/mean > 1.e-4 ){
    TF1 * f = new TF1("gaus","gaus",0,100);
    f->SetParameters(cnt,mean,sigma);
    //cout<<"gaus value "<<f->Eval(y)<<endl ;
    Double_t arg  = f->Eval(y) ;
    return arg;
  }
  else
    return 0.;

}
//________________________________________________________________________
Double_t  AliPHOSPIDv1::GausPol2(Double_t  x, Double_t y, Double_t * par)
{
  Double_t cnt    = par[0] + par[1] * x + par[2] * x * x ;
  Double_t mean   = par[3] + par[4] * x + par[5] * x * x ;
  Double_t sigma  = par[6] + par[7] * x + par[8] * x * x ;

  if(mean != 0. && sigma/mean > 1.e-4 ){
    TF1 * f = new TF1("gaus","gaus",0,100);
    f->SetParameters(cnt,mean,sigma);
    Double_t arg  = f->Eval(y) ;
    return arg;
  }
  else
    return 0.;
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

 if (particle == "photon")
    value = p[0]/TMath::Sqrt(e) + p[1]*e + p[2]*e*e + p[3];
  else if (particle == "pi0")
    value = p[0] + p[1]*e + p[2]*e*e;

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
  //Info("GetCPVBit"," xdist %f, sigx %f, zdist %f, sigz %f",deltaX, sigX, deltaZ,sigZ ) ;
  if((deltaX>sigX*(effPur+1))&&(deltaZ>sigZ*(effPur+1)))
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
  if(gAlice && gAlice->GetMCApp() && gAlice->Generator())
    gAlice->Generator()->GetOrigin(xo,yo,zo) ;
  else{
    xo=yo=zo=0.;
  }
  TVector3 origin(xo,yo,zo);
  dir = dir - origin ;

  return dir ;  
}

//________________________________________________________________________
Double_t  AliPHOSPIDv1::LandauF(Double_t  x, Double_t y, Double_t * par)
{
  Double_t cnt    = par[2] * (x*x) + par[1] * x + par[0] ;
  Double_t mean   = par[4] / (x*x) + par[5] / x + par[3] ;
  Double_t sigma  = par[7] / (x*x) + par[8] / x + par[6] ;
 //  Double_t mean   = par[3] + par[4] * x + par[5] * x * x ;
//   Double_t sigma  = par[6] + par[7] * x + par[8] * x * x ;
  
  //Double_t arg    = -(y-mean)*(y-mean)/(2*sigma*sigma) ;
  //return cnt * TMath::Exp(arg) ;
  if(mean != 0. && sigma/mean > 1.e-4 ){
    TF1 * f = new TF1("landau","landau",0.,100.);
    f->SetParameters(cnt,mean,sigma);
    Double_t arg  = f->Eval(y) ;
    return arg;
  }
  else
    return 0.;

}
//________________________________________________________________________
Double_t  AliPHOSPIDv1::LandauPol2(Double_t  x, Double_t y, Double_t * par)
{
  Double_t cnt    = par[2] * (x*x) + par[1] * x + par[0] ;
  Double_t mean   = par[4] * (x*x) + par[5] * x + par[3] ;
  Double_t sigma  = par[7] * (x*x) + par[8] * x + par[6] ;

  if(mean != 0. && sigma/mean > 1.e-4 ){
    TF1 * f = new TF1("landau","landau",0.,100.);
    f->SetParameters(cnt,mean,sigma);
    Double_t arg  = f->Eval(y) ;
    return arg;
  }
  else
    return 0.;
}
// //________________________________________________________________________
// Double_t  AliPHOSPIDv1::ChargedHadronDistProb(Double_t  x, Double_t y, Double_t * parg, Double_t * parl)
// {
//   Double_t cnt   = 0.0 ;
//   Double_t mean  = 0.0 ;
//   Double_t sigma = 0.0 ;
//   Double_t arg   = 0.0 ;
//   if (y < parl[4] / (x*x) + parl[5] / x + parl[3]){
//     cnt    = parg[1] / (x*x) + parg[2] / x + parg[0] ;
//     mean   = parg[4] / (x*x) + parg[5] / x + parg[3] ;
//     sigma  = parg[7] / (x*x) + parg[8] / x + parg[6] ;
//     TF1 * f = new TF1("gaus","gaus",0.,100.);
//     f->SetParameters(cnt,mean,sigma);
//     arg  = f->Eval(y) ;
//   }
//   else{
//     cnt    = parl[1] / (x*x) + parl[2] / x + parl[0] ;
//     mean   = parl[4] / (x*x) + parl[5] / x + parl[3] ;
//     sigma  = parl[7] / (x*x) + parl[8] / x + parl[6] ;
//     TF1 * f = new TF1("landau","landau",0.,100.);
//     f->SetParameters(cnt,mean,sigma);
//     arg  = f->Eval(y) ;
//   }
//   //  Double_t mean   = par[3] + par[4] * x + par[5] * x * x ;
//   //   Double_t sigma  = par[6] + par[7] * x + par[8] * x * x ;
  
//   //Double_t arg    = -(y-mean)*(y-mean)/(2*sigma*sigma) ;
//   //return cnt * TMath::Exp(arg) ;
  
//   return arg;
  
// }
//____________________________________________________________________________
void  AliPHOSPIDv1::MakePID()
{
  // construct the PID weight from a Bayesian Method

  Int_t index ;
  const Int_t kSPECIES = AliESDtrack::kSPECIESN ;
  Int_t nparticles = AliPHOSGetter::Instance()->RecParticles()->GetEntriesFast() ;

//   const Int_t kMAXPARTICLES = 2000 ; 
//   if (nparticles >= kMAXPARTICLES) 
//     Error("MakePID", "Change size of MAXPARTICLES") ; 
//   Double_t stof[kSPECIES][kMAXPARTICLES] ;


  Double_t * stof[kSPECIES] ;
  Double_t * sdp [kSPECIES]  ;
  Double_t * scpv[kSPECIES] ;
  
  //Info("MakePID","Begin MakePID"); 
  
  for (Int_t i =0; i< kSPECIES; i++){
    stof[i] = new Double_t[nparticles] ;
    sdp [i] = new Double_t[nparticles] ;
    scpv[i] = new Double_t[nparticles] ;
  }
  
  // make the normalized distribution of pid for this event 
  // w(pid) in the Bayesian formulation
  for(index = 0 ; index < nparticles ; index ++) {
    
    AliPHOSRecParticle * recpar  = AliPHOSGetter::Instance()->RecParticle(index) ;  
    AliPHOSEmcRecPoint * emc     = AliPHOSGetter::Instance()->EmcRecPoint(index) ;
    AliPHOSCpvRecPoint * cpv     = AliPHOSGetter::Instance()->CpvRecPoint(index) ;
    
    Float_t en   = emc->GetEnergy(); 
    
    // Tof
    Double_t time = recpar->ToF() ;
    //cout<<">>>>>>>Energy "<<en<<"Time "<<time<<endl;
    //Electrons initial population to be removed
    fInitPID[AliESDtrack::kEleCon]   = 0. ;
    
    // now get the signals probability
    // s(pid) in the Bayesian formulation
    //     stof[AliESDtrack::kPhoton][index]   = fTFphoton     ->Eval(time) ; 
    //     stof[AliESDtrack::kElectron][index] =  stof[AliESDtrack::kPhoton][index] ;
    //     if(time < fTpionl[1])
    //       stof[AliESDtrack::kPion][index]     = fTFpiong      ->Eval(time) ; //gaus distribution
    //     else
    //       stof[AliESDtrack::kPion][index]     = fTFpionl      ->Eval(time) ; //landau distribution
    //     if(time < fTkaonl[1])
    //       stof[AliESDtrack::kKaon][index]     = fTFkaong      ->Eval(time) ; //gaus distribution
    //     else 
    //       stof[AliESDtrack::kKaon][index]     = fTFkaonl      ->Eval(time) ; //landau distribution
    //     if(time < fThhadronl[1])
    //       stof[AliESDtrack::kProton][index]   = fTFhhadrong   ->Eval(time) ; //gaus distribution
    //     else
    //       stof[AliESDtrack::kProton][index]   = fTFhhadronl   ->Eval(time) ; //landau distribution
    
    //     stof[AliESDtrack::kNeutron][index]  = stof[AliESDtrack::kProton][index] ;
    //     stof[AliESDtrack::kEleCon][index]   = stof[AliESDtrack::kPhoton][index] ;
    //     // a conversion electron has the photon ToF
    //     stof[AliESDtrack::kKaon0][index]    = stof[AliESDtrack::kKaon][index] ;
    //     stof[AliESDtrack::kMuon][index]     = stof[AliESDtrack::kPhoton][index] ;
    if(en < 2.) {
      //       cout<<"TOF: pi "<< GausPol2(en, time, fTpion)<<endl;
      //       cout<<"TOF: k  "<< LandauPol2(en, time, fTkaon)<<endl;
      //       cout<<"TOF: N  "<< LandauPol2(en, time, fThhadron)<<endl;
      stof[AliESDtrack::kPhoton][index]   = fTFphoton     ->Eval(time) ; 
      stof[AliESDtrack::kElectron][index] =  stof[AliESDtrack::kPhoton][index] ;
//       stof[AliESDtrack::kPion][index]     = GausPol2(en, time, fTpion) ; //gaus distribution
//       stof[AliESDtrack::kKaon][index]     = LandauPol2(en, time, fTkaon) ; //gaus distribution
//       stof[AliESDtrack::kProton][index]   = LandauPol2(en, time, fThhadron); //gaus distribution
      stof[AliESDtrack::kPion][index]     = fTFpiong      ->Eval(time) ; //landau distribution
      if(time < fTkaonl[1])
	stof[AliESDtrack::kKaon][index]     = fTFkaong      ->Eval(time) ; //gaus distribution
      else 
	stof[AliESDtrack::kKaon][index]     = fTFkaonl      ->Eval(time) ; //landau distribution
      if(time < fThhadronl[1])
	stof[AliESDtrack::kProton][index]   = fTFhhadrong   ->Eval(time) ; //gaus distribution
      else
	stof[AliESDtrack::kProton][index]   = fTFhhadronl   ->Eval(time) ; //landau distribution

      stof[AliESDtrack::kNeutron][index]  = stof[AliESDtrack::kProton][index] ;
      stof[AliESDtrack::kEleCon][index]   = stof[AliESDtrack::kPhoton][index] ;
      // a conversion electron has the photon ToF
      stof[AliESDtrack::kKaon0][index]    = stof[AliESDtrack::kKaon][index] ;
      stof[AliESDtrack::kMuon][index]     = stof[AliESDtrack::kPhoton][index] ;
    } 
    else {
      stof[AliESDtrack::kPhoton][index]   = 1.; 
      stof[AliESDtrack::kElectron][index] = 1.;
      stof[AliESDtrack::kPion][index]     = 1.; 
      stof[AliESDtrack::kKaon][index]     = 1.; 
      stof[AliESDtrack::kProton][index]   = 1.;
      stof[AliESDtrack::kNeutron][index]  = 1.;
      stof[AliESDtrack::kEleCon][index]   = 1.;
      stof[AliESDtrack::kKaon0][index]    = 1.;
      stof[AliESDtrack::kMuon][index]     = 1.; 
    }
    //    Info("MakePID", "TOF passed");
    
    // Shower shape: Dispersion
    Float_t dispersion = emc->GetDispersion();
    //dispersion is not well defined if the cluster is only in few crystals
    
    //     Info("MakePID","multiplicity %d, dispersion %f", emc->GetMultiplicity(), 
    // 	 dispersion);
    //     Info("MakePID","ss: photon %f, hadron %f ", GausF   (en , dispersion, fDphoton), 
    //       LandauF(en , dispersion, fDhadron ) );
    if(emc->GetMultiplicity() > 4){ 
      sdp[AliESDtrack::kPhoton][index]   = GausF   (en , dispersion, fDphoton) ;
      sdp[AliESDtrack::kElectron][index] = sdp[AliESDtrack::kPhoton][index] ;
      sdp[AliESDtrack::kPion][index]     = LandauF(en , dispersion, fDhadron ) ; 
      sdp[AliESDtrack::kKaon][index]     = sdp[AliESDtrack::kPion][index]  ; 
      sdp[AliESDtrack::kProton][index]   = sdp[AliESDtrack::kPion][index]  ;
      sdp[AliESDtrack::kNeutron][index]  = sdp[AliESDtrack::kPion][index]  ;
      sdp[AliESDtrack::kEleCon][index]   = sdp[AliESDtrack::kPhoton][index]; 
      sdp[AliESDtrack::kKaon0][index]    = sdp[AliESDtrack::kPion][index]  ; 
      sdp[AliESDtrack::kMuon][index]     = fDFmuon ->Eval(dispersion) ; //landau distribution
    }
    else{
      sdp[AliESDtrack::kPhoton][index]   = 1. ;
      sdp[AliESDtrack::kElectron][index] = 1. ;
      sdp[AliESDtrack::kPion][index]     = 1. ; 
      sdp[AliESDtrack::kKaon][index]     = 1. ; 
      sdp[AliESDtrack::kProton][index]   = 1. ;
      sdp[AliESDtrack::kNeutron][index]  = 1. ; 
      sdp[AliESDtrack::kEleCon][index]   = 1. ; 
      sdp[AliESDtrack::kKaon0][index]    = 1. ; 
      sdp[AliESDtrack::kMuon][index]     = 1. ; 
    }
    
    
    // CPV-EMC  Distance
    Float_t distance = GetDistance(emc, cpv,  "R") ;
    //    Info("MakePID", "Distance %f", distance);
    Float_t pcpv = 0 ;
    Float_t pcpvneutral  = 0. ;
    Float_t pcpvelectron = GausF  (en , distance, fCPVelectron) ;
    Float_t pcpvcharged  = LandauF(en , distance, fCPVcharged) ;
    //Float_t pcpvcharged  = ChargedHadronDistProb(en , distance, fCPVchargedg, fCPVchargedl) ;
    //    Info("MakePID", "CPV: electron %f, hadron %f", pcpvelectron, pcpvcharged);
    if(pcpvelectron >= pcpvcharged)  
      pcpv = pcpvelectron ;
    else
      pcpv = pcpvcharged ;
    
    if(pcpv < 1e-4)
      {
	pcpvneutral  = 1. ;
	pcpvcharged  = 0. ;
	pcpvelectron = 0. ;
      }
    
    scpv[AliESDtrack::kPion][index]     =  pcpvcharged  ; 
    scpv[AliESDtrack::kKaon][index]     =  pcpvcharged  ; 
    scpv[AliESDtrack::kProton][index]   =  pcpvcharged  ;
    scpv[AliESDtrack::kPhoton][index]   =  pcpvneutral  ;
    scpv[AliESDtrack::kElectron][index] =  pcpvelectron ;
    scpv[AliESDtrack::kNeutron][index]  =  pcpvneutral  ; 
    scpv[AliESDtrack::kEleCon][index]   =  pcpvelectron ; 
    scpv[AliESDtrack::kKaon0][index]    =  pcpvneutral  ; 
    scpv[AliESDtrack::kMuon][index]     =  pcpvelectron ; 
    
    //   Info("MakePID", "CPV passed");
    
    if(en > 30.){
      // pi0 are detected via decay photon
      stof[AliESDtrack::kPi0][index]  = fTFphoton  ->Eval(time) ;
      scpv[AliESDtrack::kPi0][index]  = pcpvneutral  ;
      if(emc->GetMultiplicity() > 4)
	sdp [AliESDtrack::kPi0][index]  = GausPol2(en , dispersion, fDpi0) ;
      else
	sdp [AliESDtrack::kPi0][index]  = 1. ;
    }
    else{
      stof[AliESDtrack::kPi0][index]      = 0. ;  
      scpv[AliESDtrack::kPi0][index]      = 0. ;
      sdp [AliESDtrack::kPi0][index]      = 0. ;
      fInitPID[AliESDtrack::kPi0]         = 0. ;
    }
    
    if(en > 0.5){
      //Muons deposit few energy
      scpv[AliESDtrack::kMuon][index]     =  0;
      stof[AliESDtrack::kMuon][index]     =  0;
      sdp [AliESDtrack::kMuon][index]     =  0;
    }
//     cout<<"MakePID: energy "<<en<<", tof "<<time<<", distance "<<distance<<", dispersion "<<dispersion<<endl ;
//     cout<<"Photon   , pid "<< fInitPID[AliESDtrack::kPhoton]<<" tof "<<stof[AliESDtrack::kPhoton][index]
//      	<<", cpv "<<scpv[AliESDtrack::kPhoton][index]<<", ss "<<sdp[AliESDtrack::kPhoton][index]<<endl;
//     cout<<"EleCon   , pid "<< fInitPID[AliESDtrack::kEleCon]<<", tof "<<stof[AliESDtrack::kEleCon][index]
//      	<<", cpv "<<scpv[AliESDtrack::kEleCon][index]<<" ss "<<sdp[AliESDtrack::kEleCon][index]<<endl;
//     cout<<"Electron , pid "<< fInitPID[AliESDtrack::kElectron]<<", tof "<<stof[AliESDtrack::kElectron][index]
//      	<<", cpv "<<scpv[AliESDtrack::kElectron][index]<<" ss "<<sdp[AliESDtrack::kElectron][index]<<endl;
//     cout<<"Muon     , pid "<< fInitPID[AliESDtrack::kMuon]<<", tof "<<stof[AliESDtrack::kMuon][index]
//      	<<", cpv "<<scpv[AliESDtrack::kMuon][index]<<" ss "<<sdp[AliESDtrack::kMuon][index]<<endl;
//     cout<<"Pi0      , pid "<< fInitPID[AliESDtrack::kPi0]<<", tof "<<stof[AliESDtrack::kPi0][index]
//      	<<", cpv "<<scpv[AliESDtrack::kPi0][index]<<" ss "<<sdp[AliESDtrack::kPi0][index]<<endl;
//     cout<<"Pion     , pid "<< fInitPID[AliESDtrack::kPion]<<", tof "<<stof[AliESDtrack::kPion][index]
//      	<<", cpv "<<scpv[AliESDtrack::kPion][index]<<" ss "<<sdp[AliESDtrack::kPion][index]<<endl;
//     cout<<"Kaon0    , pid "<< fInitPID[AliESDtrack::kKaon0]<<", tof "<<stof[AliESDtrack::kKaon0][index]
//      	<<", cpv "<<scpv[AliESDtrack::kKaon0][index]<<" ss "<<sdp[AliESDtrack::kKaon0][index]<<endl;
//     cout<<"Kaon     , pid "<< fInitPID[AliESDtrack::kKaon]<<", tof "<<stof[AliESDtrack::kKaon][index]
//      	<<", cpv "<<scpv[AliESDtrack::kKaon][index]<<" ss "<<sdp[AliESDtrack::kKaon][index]<<endl;
//     cout<<"Neutron  , pid "<< fInitPID[AliESDtrack::kNeutron]<<", tof "<<stof[AliESDtrack::kNeutron][index]
//      	<<", cpv "<<scpv[AliESDtrack::kNeutron][index]<<" ss "<<sdp[AliESDtrack::kNeutron][index]<<endl;
//     cout<<"Proton   , pid "<< fInitPID[AliESDtrack::kProton]<<", tof "<<stof[AliESDtrack::kProton][index]
//      	<<", cpv "<<scpv[AliESDtrack::kProton][index]<<" ss "<<sdp[AliESDtrack::kProton][index]<<endl;
  }
  
  //for (index = 0 ; index < kSPECIES ; index++) 
  // pid[index] /= nparticles ; 
  
  //  Info("MakePID", "Total Probability calculation");
  
  for(index = 0 ; index < nparticles ; index ++) {
    // calculates the Bayesian weight
    Int_t jndex ;
    Double_t wn = 0.0 ; 
    for (jndex = 0 ; jndex < kSPECIES ; jndex++) 
      //wn += stof[jndex][index] * pid[jndex] ;
      wn += stof[jndex][index] * sdp[jndex][index]  * scpv[jndex][index] * fInitPID[jndex] ;
    //cout<<"*************wn "<<wn<<endl;
    AliPHOSRecParticle * recpar = AliPHOSGetter::Instance()->RecParticle(index) ;  
    if (TMath::Abs(wn)>0)
      for (jndex = 0 ; jndex < kSPECIES ; jndex++) {
	//cout<<"jndex "<<jndex<<" wn "<<wn<<" SetPID * wn"
	//<<stof[jndex][index] * sdp[jndex][index] * pid[jndex]  << endl;
	//cout<<" tof "<<stof[jndex][index] << " disp " <<sdp[jndex][index] << " pid "<< fInitPID[jndex] << endl;
// 	cout<<"Particle "<<jndex<<"  final prob * wn   "
// 	    <<stof[jndex][index] * sdp[jndex][index] * scpv[jndex][index] * fInitPID[jndex] <<"  wn  "<< wn<<endl;
	recpar->SetPID(jndex, stof[jndex][index] * sdp[jndex][index] * 
		       scpv[jndex][index] * fInitPID[jndex] / wn) ; 
// 	cout<<"final prob "<<stof[jndex][index] * sdp[jndex][index] * scpv[jndex][index] * fInitPID[jndex] / wn<<endl;
	//recpar->SetPID(jndex, stof[jndex][index] * fInitPID[jndex] / wn) ; 
	//cout<<"After SetPID"<<endl;
	//recpar->Print();
      }
  }
  //  Info("MakePID", "Delete");
  
  //   for (Int_t i =0; i< kSPECIES; i++){
  //     delete [] stof[i];
  //     cout<<i<<endl;
  //     delete [] sdp[i];
  //     cout<<i<<endl;
  //     delete [] scpv[i];
  //     cout<<i<<endl;
  //   }
  
  //  Info("MakePID","End MakePID"); 
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
    
    Float_t time = emc->GetTime() ;
    rp->SetTof(time) ; 
    
    // Loop of Efficiency-Purity (the 3 points of purity or efficiency 
    // are taken into account to set the particle identification)
    for(Int_t effPur = 0; effPur < 3 ; effPur++){
      
      // Looking at the CPV detector. If RCPV greater than CpvEmcDistance, 
      // 1st,2nd or 3rd bit (depending on the efficiency-purity point )
      // is set to 1
      if(GetCPVBit(emc, cpv, effPur,e) == 1 ){  
	rp->SetPIDBit(effPur) ;
	//cout<<"CPV bit "<<effPur<<endl;
      }
      // Looking the TOF. If TOF smaller than gate,  4th, 5th or 6th 
      // bit (depending on the efficiency-purity point )is set to 1             
      if(time< (*fParameters)(3,effPur)) 
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
    //Info("SetParameters", "line %d: %s",i,string);
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
  if(fWrite){
    TTree * treeP =  gime->TreeP();
    
    //First rp
    Int_t bufferSize = 32000 ;
    TBranch * rpBranch = treeP->Branch("PHOSRP",&recParticles,bufferSize);
    rpBranch->SetTitle(BranchName());
    
    rpBranch->Fill() ;
    
    gime->WriteRecParticles("OVERWRITE");
    gime->WritePID("OVERWRITE");
  }
}


//_______________________________________________________________________
void AliPHOSPIDv1::SetInitPID(const Double_t *p) {
  // Sets values for the initial population of each particle type 
  for (Int_t i=0; i<AliESDtrack::kSPECIESN; i++) fInitPID[i] = p[i];
}
//_______________________________________________________________________
void AliPHOSPIDv1::GetInitPID(Double_t *p) const {
  // Gets values for the initial population of each particle type 
  for (Int_t i=0; i<AliESDtrack::kSPECIESN; i++) p[i] = fInitPID[i];
}
