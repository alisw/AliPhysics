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
//     -Bit 0 to 2: bit set if RCPV > fCpvEmcDistance (each bit corresponds
//      to a different efficiency-purity point of the photon identification) 
//     -Bit 3 to 5: bit set if TOF  < fTimeGate (each bit corresponds
//      to a different efficiency-purity point of the photon identification) 
//     -Bit 6 to 9: bit set if Principal Components are 
//      inside an ellipse defined by fX_center, fY_center, fA, fB, fAngle
//      (each bit corresponds to a different efficiency-purity point of the 
//      photon identification) 
//
//
// use case:
//  root [0] AliPHOSPIDv1 * p1 = new AliPHOSPIDv1("galice.root")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//  root [1] p1->SetIdentificationMethod("disp ellipse")
//  root [2] p1->ExecuteTask()
//  root [3] AliPHOSPIDv1 * p2 = new AliPHOSPIDv1("galice1.root","v1")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                // reading headers from file galice1.root and create  RecParticles with title v1
                  // TrackSegments and RecPoints with title "v1" are used 
//                // set file name for the branch RecParticles
//  root [4] p2->ExecuteTask("deb all time")
//                // available options
//                // "deb" - prints # of reconstructed particles
//                // "deb all" -  prints # and list of RecParticles
//                // "time" - prints benchmarking results
//                  
//  root [5] AliPHOSPIDv1 * p3 = new AliPHOSPIDv1("galice1.root","v1","v0")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                // reading headers from file galice1.root and create  RecParticles with title v1
                  // RecPoints and TrackSegments with title "v0" are used 
//  root [6] p3->ExecuteTask()
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
#include "TMatrixD.h"
#include "TPrincipal.h"
#include "TSystem.h"

// --- Standard library ---

#include <iostream.h>
#include <fstream.h>
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
 
  fFileName           = "" ; 
  fFileNamePar        = "" ;      
  fFrom               = "" ;
  fHeaderFileName     = "" ; 
  fOptFileName        = "Default" ;
  fTrackSegmentsTitle = "" ; 
  fRecPointsTitle     = "" ;   
  fRecParticlesTitle  = "" ; 
  
  fNEvent            = 0 ;            
  fClusterizer       = 0 ;      
  fTSMaker           = 0 ;        
  fRecParticlesInRun = 0 ;
  fX                 = 0 ;
  fP                 = 0 ; 
  fParameters        = 0 ;     

}

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1(const char * headerFile,const char * name, const char * from) : AliPHOSPID(headerFile, name)

			  
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
  if ( from == 0 ) 
    fFrom = name ; 
  else
    fFrom = from ; 
  fOptFileName        = "Default" ;
  Init() ;

}

//____________________________________________________________________________
AliPHOSPIDv1::~AliPHOSPIDv1()
{ 

  delete [] fX ; // Principal input 
  delete [] fP ; // Principal components
  delete fParameters ; // Matrix of Parameters 
}

//____________________________________________________________________________
const TString AliPHOSPIDv1::BranchName() const 
{  
  TString branchName(GetName() ) ;
  branchName.Remove(branchName.Index(Version())-1) ;
  return branchName ;
}
 
//____________________________________________________________________________
void AliPHOSPIDv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  // Add the PID task to the list of PHOS tasks

  if ( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;
  
  SetParameters() ; // fill the parameters matrix from parameters file
  
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(GetTitle(), fFrom.Data()) ; 

  gime->SetRecParticlesTitle(BranchName()) ;
  if ( gime == 0 ) {
    cerr << "ERROR: AliPHOSPIDv1::Init -> Could not obtain the Getter object !" << endl ; 
    return ;
  } 
  
  gime->PostPID(this) ;
  // create a folder on the white board //YSAlice/WhiteBoard/RecParticles/PHOS/recparticlesName
  gime->PostRecParticles(BranchName()) ; 
  
}
//____________________________________________________________________________
Double_t  AliPHOSPIDv1::GetCpvtoEmcDistanceCut(const Float_t Cluster_En, const TString Eff_Pur)const 
{
  // Get CpvtoEmcDistanceCut parameter depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
 
  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;

  // Check the cluster energy range  
  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster = 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;
  
  if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
    cout<<"Invalid Efficiency-Purity option"<<endl;
  }
  else{
    return (*fParameters)(cluster,eff_pur) ;
  }
}
//____________________________________________________________________________

Double_t  AliPHOSPIDv1::GetTimeGate(const Float_t Cluster_En, const TString Eff_Pur) const 
{
  // Get TimeGate parameter depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;
  
  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster = 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;
 
   if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
    cout<<"Invalid Efficiency-Purity option"<<endl; 
  }
  else{
    return (*fParameters)(cluster+3,eff_pur) ; 
  }
}
//_____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetDistance(AliPHOSEmcRecPoint * emc,AliPHOSRecPoint * cpv, Option_t *  Axis)const
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
  
  const AliPHOSGeometry * geom = AliPHOSGetter::GetInstance()->PHOSGeometry() ; 
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
      if (Axis == "X") return vecCpv.X();
      if (Axis == "Y") return vecCpv.Y();
      if (Axis == "Z") return vecCpv.Z();
      if (Axis == "R") return vecCpv.Mag();
  } 
    
    return 100000000 ;
  }
  return 100000000 ;
}

//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetPrincipalSign(Double_t* P, Int_t cluster, Int_t eff_pur )const
{
  //This method gives if the PCA of the particle are inside a defined ellipse
  
  // Get the parameters that define the ellipse stored in the 
  // fParameters matrix.
  Double_t X_center = (*fParameters)(cluster+6,eff_pur) ; 
  Double_t Y_center = (*fParameters)(cluster+9,eff_pur) ; 
  Double_t A        = (*fParameters)(cluster+12,eff_pur) ; 
  Double_t B        = (*fParameters)(cluster+15,eff_pur) ; 
  Double_t Angle    = (*fParameters)(cluster+18,eff_pur) ;

  Int_t      prinsign;
  Double_t   Dx        = 0. ; 
  Double_t   Delta     = 0. ; 
  Double_t   Y         = 0. ; 
  Double_t   Y_1       = 0. ; 
  Double_t   Y_2       = 0. ;
  Double_t   Pi        = TMath::Pi() ;
  Double_t   Cos_Theta = TMath::Cos(Pi*Angle/180.) ;
  Double_t   Sin_Theta = TMath::Sin(Pi*Angle/180.) ;   

  Dx = P[0] - X_center ; 
  Delta = 4.*A*A*A*B* (A*A*Cos_Theta*Cos_Theta 
			+ B*B*Sin_Theta*Sin_Theta - Dx*Dx) ; 
  if (Delta < 0.) 
    {prinsign=0;} 
  
  else if (Delta == 0.) 
    { 
      Y = Cos_Theta*Sin_Theta*(A*A - B*B)*Dx / 
	(A*A*Cos_Theta*Cos_Theta + B*B*Sin_Theta*Sin_Theta) ; 
      Y += Y_center ; 
      if(P[1]==Y ) 
	{prinsign=1;} 
      else 
	{prinsign=0;} 
    } 
  else 
    { 
      Y_1 = (Cos_Theta*Sin_Theta*(A*A - B*B) *Dx +
	     TMath::Sqrt(Delta)/2.)/(A*A*Cos_Theta*Cos_Theta + 
				     B*B*Sin_Theta*Sin_Theta) ; 
      Y_2 = (Cos_Theta*Sin_Theta*(A*A - B*B) *Dx -
	     TMath::Sqrt(Delta)/2.)/(A*A*Cos_Theta*Cos_Theta 
				      + B*B*Sin_Theta*Sin_Theta) ; 
      Y_1 += Y_center ; 
      Y_2 += Y_center ; 
      if ((P[1]<=Y_1) && (P[1]>=Y_2)) 
	{prinsign=1;} 
      else 
	{prinsign=0;}  
    } 
  return prinsign;
}

//____________________________________________________________________________
void   AliPHOSPIDv1::SetPrincipalFileOptions(TString OptFileName) {

  if(OptFileName.Contains("Small energy range")||OptFileName.Contains("Default")){
    fFileName  = "$ALICE_ROOT/PHOS/PCA8pa15_0.5-5.root" ;
    fFileNamePar = gSystem->ExpandPathName("$ALICE_ROOT/PHOS/Parameters_0.5_5.dat"); 
  }
  
  if(OptFileName.Contains("Wide energy range")){
    fFileName  = "$ALICE_ROOT/PHOS/PCA8pa15_0.5-100.root" ;
    fFileNamePar = gSystem->ExpandPathName("$ALICE_ROOT/PHOS/Parameters_0.5_100.dat"); 
  }
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetEllipseParameters(Float_t Cluster_En, TString Eff_Pur, Float_t x, Float_t y,Float_t a, Float_t b,Float_t angle)
{

  // Set all ellipse parameters depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)

  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;
  
  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster= 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;
  
  if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
    cout<<"Invalid Efficiency-Purity option"<<endl;
  }
  else{    
    (*fParameters)(cluster+6,eff_pur) = x ;
    (*fParameters)(cluster+9,eff_pur) = y ;
    (*fParameters)(cluster+12,eff_pur) = a ;
    (*fParameters)(cluster+15,eff_pur) = b ;
    (*fParameters)(cluster+18,eff_pur) = angle ;
  }

}
//__________________________________________________________________________ 
void  AliPHOSPIDv1::SetEllipseXCenter(Float_t Cluster_En, TString Eff_Pur, Float_t x) 
{
  // Set the ellipse parameter x_center depending on the custer energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)

  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;
  
  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster = 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;
   
  if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
    cout<<"Invalid Efficiency-Purity option"<<endl;
  }
  else{
    (*fParameters)(cluster+6,eff_pur) = x ; 
  }
}
//_________________________________________________________________________    
void  AliPHOSPIDv1::SetEllipseYCenter(Float_t Cluster_En, TString Eff_Pur, Float_t y) 
{

  // Set the ellipse parameter y_center depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
  
  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;
  
  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster = 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;

  if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
    cout<<"Invalid Efficiency-Purity option"<<endl ; 
  }
  else{
    (*fParameters)(cluster+9,eff_pur) = y ;
  }
}
//_________________________________________________________________________
void  AliPHOSPIDv1::SetEllipseAParameter(Float_t Cluster_En, TString Eff_Pur, Float_t a) 
{
  // Set the ellipse parameter a depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)

  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;

  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster = 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;
  
  if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
    cout<<"Invalid Efficiency-Purity option"<<endl; 
  }
  else{ 
    (*fParameters)(cluster+12,eff_pur) = a ;    
  }
}
//________________________________________________________________________
void  AliPHOSPIDv1::SetEllipseBParameter(Float_t Cluster_En, TString Eff_Pur, Float_t b) 
{
  // Set the ellipse parameter b depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)

  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;
  
  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster = 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;
  
  if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
   cout<<"Invalid Efficiency-Purity option"<<endl; 
  }
  else{
    (*fParameters)(cluster+15,eff_pur) = b ;
  }
}
//________________________________________________________________________
void  AliPHOSPIDv1::SetEllipseAngle(Float_t Cluster_En, TString Eff_Pur, Float_t angle) 
{

  // Set the ellipse parameter angle depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)

  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;

  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster = 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;

  if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
    cout<<"Invalid Efficiency-Purity option"<<endl; 
  }
  else{
    (*fParameters)(cluster+18,eff_pur) = angle ;
  }
} 
//_____________________________________________________________________________
void  AliPHOSPIDv1::SetCpvtoEmcDistanceCut(Float_t Cluster_En, TString Eff_Pur, Float_t cut) 
{

  // Set the parameter Cpvto EmcDistanceCut depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)

  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;

  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster = 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;

  if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
    cout<<"Invalid Efficiency-Purity option"<<endl; 
  }
  else{
    (*fParameters)(cluster,eff_pur) = cut ;
  }
}  

//_____________________________________________________________________________
void  AliPHOSPIDv1::SetParameters() 
{
  // PCA : To do the Principal Components Analysis it is necessary 
  // the Principal file, which is opened here
  fX         = new double[7]; // Data for the PCA 
  fP         = new double[7]; // Eigenvalues of the PCA  
  
  SetPrincipalFileOptions(fOptFileName);
  TFile f( fFileName.Data(), "read" ) ;
  fPrincipal = dynamic_cast<TPrincipal*> (f.Get("principal")) ; 
  f.Close() ; 

  // Initialization of the Parameters matrix. In the File Parameters.dat
  // are all the parameters. These are introduced in a matrix of 21x3 elements.
  // All the parameters defined in this file are, in order of row (there are
  // 3 rows per parameter): CpvtoEmcDistanceCut, TimeGate (and the ellipse 
  // parameters), X_center, Y_center, a, b, angle. Each row of a given parameter
  // depends on the cluster energy range (0.3-1,1-2, >2 )
  // Each column designes the parameters for a point in the Efficiency-Purity
  // of the photon identification P1(0.959,0.625), P2(0.919,0.835), P3(0.833,0.901).
 
  fParameters = new TMatrixD(21,3) ; 

  ifstream paramFile(fFileNamePar, ios::in) ; 
  
  Int_t i,j ;
  
  for(i = 0; i< 21; i++){
    for(j = 0; j< 3; j++){
      paramFile >> (*fParameters)(i,j) ;
    }
  }
  paramFile.close();  
}

//_____________________________________________________________________________
void  AliPHOSPIDv1::SetTimeGate(Float_t Cluster_En, TString Eff_Pur, Float_t gate) 
{

  // Set the parameter TimeGate depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)

  Int_t cluster = -1 ;
  Int_t eff_pur = -1 ;

  if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)) cluster = 0 ;
  if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)) cluster = 1 ;
  if( Cluster_En > 2.0) cluster = 2 ;
  
  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") ) eff_pur = 0 ;
  if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) eff_pur = 1 ;
  if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) eff_pur = 2 ;

  if(cluster ==-1){
    cout<<"Invalid Cluster Energy option"<<endl;
  }
  else if(eff_pur ==-1){
    cout<<"Invalid Efficiency-Purity option"<<endl; 
  }
  else{
    (*fParameters)(cluster+3,eff_pur) = gate ;
  }
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

  //cout << gDirectory->GetName() << endl ; 

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
    fRecParticlesInRun += gime->RecParticles(BranchName())->GetEntriesFast() ; 

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
void  AliPHOSPIDv1::MakeRecParticles(){

  // Makes a RecParticle out of a TrackSegment
  
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  TObjArray * emcRecPoints = gime->EmcRecPoints(fFrom) ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints(fFrom) ; 
  TClonesArray * trackSegments = gime->TrackSegments(fFrom) ; 
  if ( !emcRecPoints || !cpvRecPoints || !trackSegments ) {
    cerr << "ERROR:  AliPHOSPIDv1::MakeRecParticles -> RecPoints or TrackSegments with name " 
	 << fFrom << " not found ! " << endl ; 
    abort() ; 
  }
  TClonesArray * recParticles  = gime->RecParticles(BranchName()) ; 
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
    
    // Now set type (reconstructed) of the particle

    // Choose the cluster energy range
    Int_t cluster = 0 ; // Ellipse and rcpv cut in function of the cluster energy
    if((e > 0.3)&&(e <= 1.0)) cluster = 0 ;
    if((e > 1.0)&&(e <= 2.0)) cluster = 1 ;
    if( e > 2.0) cluster = 2 ;

    // Loop of Efficiency-Purity (the 3 points of purity or efficiency are taken 
    // into account to set the particle identification)
    for(Int_t eff_pur = 0; eff_pur < 3 ; eff_pur++){

    // Looking at the CPV detector. If RCPV greater than CpvEmcDistance, 1st, 
    // 2nd or 3rd bit (depending on the efficiency-purity point )is set to 1 .  
      if(GetDistance(emc, cpv,  "R") > (*fParameters)(cluster,eff_pur) )  
	rp->SetPIDBit(eff_pur) ;
    
    // Looking the TOF. If TOF smaller than gate,  4th, 5th or 6th 
    // bit (depending on the efficiency-purity point )is set to 1             
      if(emc->GetTime()< (*fParameters)(cluster+3,eff_pur))  
	rp->SetPIDBit(eff_pur+3) ;		    
     
    // Looking PCA. Define and calculate the data (X), introduce in the function 
    // X2P that gives the components (P).  
    Float_t    fSpher = 0. ;
    Float_t    fEmaxdtotal = 0. ; 
    Float_t lambda[2] ;
    
    emc->GetElipsAxis(lambda) ;
    if((lambda[0]+lambda[1])!=0) fSpher=fabs(lambda[0]-lambda[1])/(lambda[0]+lambda[1]); 
    
    fEmaxdtotal=emc->GetMaximalEnergy()/emc->GetEnergy(); 
    
    fX[0] = lambda[0] ;  
    fX[1] = lambda[1] ; 
    fX[2] = emc->GetDispersion() ; 
    fX[3] = fSpher ; 
    fX[4] = emc->GetMultiplicity() ;  
    fX[5] = fEmaxdtotal ;  
    fX[6] = emc->GetCoreEnergy() ;  
    
    fPrincipal->X2P(fX,fP);
    
    //If we are inside the ellipse, 7th, 8th or 9th 
    // bit (depending on the efficiency-purity point )is set to 1 
    if(GetPrincipalSign(fP,cluster,eff_pur) == 1) 
      rp->SetPIDBit(eff_pur+6) ;
 
  }
  
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
void  AliPHOSPIDv1:: Print()
{
  // Print the parameters used for the particle type identification
    cout <<  "=============== AliPHOSPID1 ================" << endl ;
    cout <<  "Making PID "<< endl ;
    cout <<  "    Headers file:               " << fHeaderFileName.Data() << endl ;
    cout <<  "    RecPoints branch title:     " << fRecPointsTitle.Data() << endl ;
    cout <<  "    TrackSegments Branch title: " << fTrackSegmentsTitle.Data() << endl ;
    cout <<  "    RecParticles Branch title   " << fRecParticlesTitle.Data() << endl;
    cout <<  "    Matrix of Parameters: "<<endl;
    cout <<  "           3  Columns [High Eff-Low Pur,Medium Eff-Pur, Low Eff-High Pur]"<<endl;
    cout <<  "           21 Rows, each 3 [ RCPV, TOF, X_Center, Y_Center, A, B, Angle ]"<<endl;
    SetParameters() ;
    fParameters->Print() ; 
    cout <<  "============================================" << endl ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::WriteRecParticles(Int_t event)
{
 
  AliPHOSGetter *gime = AliPHOSGetter::GetInstance() ; 

  TClonesArray * recParticles = gime->RecParticles(BranchName()) ; 
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

  delete [] filename ; 
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

  TClonesArray * recParticles = gime->RecParticles(BranchName()) ; 
  
  cout << "AliPHOSPIDv1: event "<<gAlice->GetEvNumber()  << endl ;
  cout << "       found " << recParticles->GetEntriesFast() << " RecParticles " << endl ;
  
  if(strstr(option,"all")) {  // printing found TS
    
    cout << "  PARTICLE       "   
	 << "  Index    "  << endl ;
    
    Int_t index ;
    for (index = 0 ; index < recParticles->GetEntries() ; index++) {
       AliPHOSRecParticle * rp = (AliPHOSRecParticle * ) recParticles->At(index) ;       

       cout << setw(10) << rp->Name() << "  "
	    << setw(5) <<  rp->GetIndexInList() << " " <<endl;
       cout << "Type "<<  rp->GetType() << endl;
    }
    cout << "-------------------------------------------" << endl ;
  }
  
}



