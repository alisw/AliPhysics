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
//      A calibrated energy is calculated. The energy of the reconstructed
//      cluster is corrected with the formula A + B * E  + C * E^2, whose parameters
//      where obtained thourgh the study of the reconstructed energy 
//      distribution of monoenergetic photons. 
//
//
//
// use case:
//  root [0] AliPHOSPIDv1 * p = new AliPHOSPIDv1("galice1.root","v1")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                // reading headers from file galice1.root and create  RecParticles with title v1
                  // TrackSegments and RecPoints with title "v1" are used 
//                // set file name for the branch RecParticles
//  root [1] p->ExecuteTask("deb all time")
//                // available options
//                // "deb" - prints # of reconstructed particles
//                // "deb all" -  prints # and list of RecParticles
//                // "time" - prints benchmarking results
//                  
//  root [2] AliPHOSPIDv1 * p2 = new AliPHOSPIDv1("galice1.root","v1","v0")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                // reading headers from file galice1.root and create  RecParticles with title v1
                  // RecPoints and TrackSegments with title "v0" are used 
//  root [3] p2->ExecuteTask()
//
//  There are two possible principal files available to do the analysis. 
//  One for energy ranges from 0.5 to 5 GeV, and another 
//  one from 5 to 100 GeV. This files are automatically called in function
//  of the cluster energy.

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

#include <iostream>
#include <fstream>
#include <iomanip>

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
 
  InitParameters() ; 
  fDefaultInit = kTRUE ; 

}

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1(const char * headerFile,const char * name, const Bool_t toSplit)
:AliPHOSPID(headerFile, name,toSplit)
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
  // fDefaultInit = kTRUE if PID created by default ctor (to get just the parameters)

  delete [] fX ; // Principal input 
  delete [] fP ; // Principal components
//  delete fParameters ; // Matrix of Parameters 
//  delete fParameters5 ; // Matrix of Parameters 
//  delete fParameters100 ; // Matrix of Parameters 
 

  if (!fDefaultInit) {  
//    AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
    // remove the task from the folder list
//     gime->RemoveTask("P",GetName()) ;
//     TString name(GetName()) ; 
//     name.ReplaceAll("pid", "clu") ; 
//     gime->RemoveTask("C",name) ;
    
//     // remove the data from the folder list
//     name = GetName() ; 
//     name.Remove(name.Index(":")) ; 
//     gime->RemoveObjects("RE", name) ; // EMCARecPoints
//     gime->RemoveObjects("RC", name) ; // CPVRecPoints
//     gime->RemoveObjects("T", name) ;  // TrackSegments
//     gime->RemoveObjects("P", name) ;  // RecParticles
    
//     // Delete gAlice
//     gime->CloseFile() ; 
    
    fSplitFile = 0 ; 
  }
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

  TString branchname(GetName()) ;
  branchname.Remove(branchname.Index(Version())-1) ;    
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(GetTitle(),branchname.Data(),fToSplit ) ; 

  //  gime->SetRecParticlesTitle(BranchName()) ;
  if ( gime == 0 ) {
    cerr << "ERROR: AliPHOSPIDv1::Init -> Could not obtain the Getter object !" << endl ; 
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
    fileName+="PHOS.RecData." ;
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
  // create a folder on the white board //YSAlice/WhiteBoard/RecParticles/PHOS/recparticlesName
  gime->PostRecParticles(branchname) ; 
  
}

//____________________________________________________________________________
void AliPHOSPIDv1::InitParameters()
{
//   fFrom               = "" ;
//   fHeaderFileName     = GetTitle() ; 
//   TString name(GetName()) ; 
//   if (name.IsNull()) 
//     name = "Default" ;
//   fTrackSegmentsTitle = name ; 
//   fRecPointsTitle     = name ; 
//   fRecParticlesTitle  = name ;
//   name.Append(":") ;
//   name.Append(Version()) ; 
//   SetName(name) ; 
  fRecParticlesInRun = 0 ; 
  fNEvent            = 0 ;            
  //  fClusterizer       = 0 ;      
  //  fTSMaker           = 0 ;        
  fRecParticlesInRun = 0 ;
  TString pidName( GetName()) ;
  if (pidName.IsNull() ) 
    pidName = "Default" ; 
  pidName.Append(":") ; 
  pidName.Append(Version()) ; 
  SetName(pidName) ;
  SetParameters() ; // fill the parameters matrix from parameters file
}

//____________________________________________________________________________
Double_t  AliPHOSPIDv1::GetCpvtoEmcDistanceCut(const Float_t Cluster_En, const TString Eff_Pur)
{
  // Get CpvtoEmcDistanceCut parameter depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
  
  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  
  GetAnalysisParameters(Cluster_En) ;
  if((fClusterrcpv!= -1)&&(eff_pur != -1))
    return (*fParameters)(fClusterrcpv,eff_pur) ;
  else
    return 0.0;
}
//____________________________________________________________________________

Double_t  AliPHOSPIDv1::GetTimeGate(const Float_t Cluster_En, const TString Eff_Pur)  
{
  // Get TimeGate parameter depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
 
  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  GetAnalysisParameters(Cluster_En) ;

  if((fCluster!= -1)&&(eff_pur != -1))
    return (*fParameters)(fCluster+3+fMatrixExtraRow,eff_pur) ; 
  else
    return 0.0;

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
Double_t  AliPHOSPIDv1::CalibratedEnergy(Float_t e){
  //It calibrates Energy depending on the recpoint energy.
//      The energy of the reconstructed
//      cluster is corrected with the formula A + B* E  + C* E^2, whose parameters
//      where obtained through the study of the reconstructed energy 
//      distribution of monoenergetic photons. 
  Double_t enerec; 
    enerec = fACalParameter + fBCalParameter * e+ fCCalParameter * e * e;
  return enerec ;

}
//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetPrincipalSign(Double_t* P, Int_t cluster, Int_t eff_pur)const
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
  Delta = 4.*A*A*B*B* (A*A*Cos_Theta*Cos_Theta 
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
void  AliPHOSPIDv1::SetEllipseParameters(Float_t Cluster_En, TString Eff_Pur, Float_t x, Float_t y,Float_t a, Float_t b,Float_t angle)
{

  // Set all ellipse parameters depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
  
  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  GetAnalysisParameters(Cluster_En) ;
  if((fCluster!= -1)&&(eff_pur != -1)){   
    (*fParameters)(fCluster+6 +fMatrixExtraRow,eff_pur) = x ;
    (*fParameters)(fCluster+9 +fMatrixExtraRow,eff_pur) = y ;
    (*fParameters)(fCluster+12+fMatrixExtraRow,eff_pur) = a ;
    (*fParameters)(fCluster+15+fMatrixExtraRow,eff_pur) = b ;
    (*fParameters)(fCluster+18+fMatrixExtraRow,eff_pur) = angle ;
  }
  
}
//__________________________________________________________________________ 
void  AliPHOSPIDv1::SetEllipseXCenter(Float_t Cluster_En, TString Eff_Pur, Float_t x) 
{
  // Set the ellipse parameter x_center depending on the custer energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  GetAnalysisParameters(Cluster_En) ;
  if((fCluster!= -1)&&(eff_pur != -1))
    (*fParameters)(fCluster+6+fMatrixExtraRow,eff_pur) = x ; 
}
//_________________________________________________________________________    
void  AliPHOSPIDv1::SetEllipseYCenter(Float_t Cluster_En, TString Eff_Pur, Float_t y) 
{

  // Set the ellipse parameter y_center depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
 
  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  GetAnalysisParameters(Cluster_En) ;
  if((fCluster!= -1)&&(eff_pur != -1))
    (*fParameters)(fCluster+9+fMatrixExtraRow,eff_pur) = y ;
}
//_________________________________________________________________________
void  AliPHOSPIDv1::SetEllipseAParameter(Float_t Cluster_En, TString Eff_Pur, Float_t a) 
{
  // Set the ellipse parameter a depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
  
  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  GetAnalysisParameters(Cluster_En) ;
  if((fCluster!= -1)&&(eff_pur != -1)) 
    (*fParameters)(fCluster+12+fMatrixExtraRow,eff_pur) = a ;    
}
//________________________________________________________________________
void  AliPHOSPIDv1::SetEllipseBParameter(Float_t Cluster_En, TString Eff_Pur, Float_t b) 
{
  // Set the ellipse parameter b depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
  
  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  GetAnalysisParameters(Cluster_En) ;
  if((fCluster!= -1)&&(eff_pur != -1))
    (*fParameters)(fCluster+15+fMatrixExtraRow,eff_pur) = b ;
}
//________________________________________________________________________
void  AliPHOSPIDv1::SetEllipseAngle(Float_t Cluster_En, TString Eff_Pur, Float_t angle) 
{

  // Set the ellipse parameter angle depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
 
  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  GetAnalysisParameters(Cluster_En) ;
  if((fCluster!= -1)&&(eff_pur != -1))
    (*fParameters)(fCluster+18+fMatrixExtraRow,eff_pur) = angle ;
} 
//_____________________________________________________________________________
void  AliPHOSPIDv1::SetCpvtoEmcDistanceCut(Float_t Cluster_En, TString Eff_Pur, Float_t cut) 
{

  // Set the parameter Cpvto EmcDistanceCut depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)


  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  GetAnalysisParameters(Cluster_En) ;
  if((fClusterrcpv!= -1)&&(eff_pur != -1))
    (*fParameters)(fClusterrcpv,eff_pur) = cut ;
}
//_____________________________________________________________________________
void  AliPHOSPIDv1::SetTimeGate(Float_t Cluster_En, TString Eff_Pur, Float_t gate) 
{

  // Set the parameter TimeGate depending on the cluster energy and 
  // Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)
    
  Int_t eff_pur = GetEffPurOption(Eff_Pur);
  GetAnalysisParameters(Cluster_En) ;
  if((fCluster!= -1)&&(eff_pur != -1))
    (*fParameters)(fCluster+3+fMatrixExtraRow,eff_pur) = gate ;
} 
//_____________________________________________________________________________
void  AliPHOSPIDv1::SetParameters()
				  //TString OptFileName) 
{
  // PCA : To do the Principal Components Analysis it is necessary 
  // the Principal file, which is opened here
  fX         = new double[7]; // Data for the PCA 
  fP         = new double[7]; // Eigenvalues of the PCA
  

  // Set the principal and parameters files to be used
  fFileName5  = "$ALICE_ROOT/PHOS/PCA8pa15_0.5-5.root" ;
  fFileNamePar5 = gSystem->ExpandPathName("$ALICE_ROOT/PHOS/Parameters_0.5_5.dat"); 
  fFileName100  = "$ALICE_ROOT/PHOS/PCA8pa15_0.5-100.root" ;
  fFileNamePar100 = gSystem->ExpandPathName("$ALICE_ROOT/PHOS/Parameters_0.5_100.dat"); 

  //SetPrincipalFileOptions();
  //fOptFileName);
  TFile f5( fFileName5.Data(), "read" ) ;
  fPrincipal5 = dynamic_cast<TPrincipal*> (f5.Get("principal")) ; 
  f5.Close() ; 
  TFile f100( fFileName100.Data(), "read" ) ;
  fPrincipal100 = dynamic_cast<TPrincipal*> (f100.Get("principal")) ; 
  f100.Close() ; 
  TFile f( fFileName100.Data(), "read" ) ;
  fPrincipal = dynamic_cast<TPrincipal*> (f.Get("principal")) ; 
  f.Close() ; 
  // Initialization of the Parameters matrix. In the File ParametersXX.dat
  // are all the parameters. These are introduced in a matrix of 21x3 or 22x3 
  // elements (depending on the principal file 21 rows for 0.5-5 GeV and 22 
  // rows for 5-100).
  // All the parameters defined in this file are, in order of row (there are
  // 3 rows per parameter): CpvtoEmcDistanceCut(if the principal file is 5-100 
  // GeV then 4 rows), TimeGate and the ellipse parameters, X_center, Y_center,
  // a, b, angle. Each row of a given parameter depends on the cluster energy range 
  // (wich depends on the chosen principal file)
  // Each column designs the parameters for a point in the Efficiency-Purity
  // of the photon identification P1(96%,63%), P2(87%,0.88%) and P3(68%,94%) 
  // for the principal file from 0.5-5 GeV and for the other one P1(95%,79%),
  // P2(89%,90%) and P3(72%,96%)

  fEnergyAnalysisCut = 5.; // Energy cut to change PCA

  fParameters5 = new TMatrixD(21,3) ; 
  fParameters100 = new TMatrixD(22,3) ; 
  fParameters = new TMatrixD(22,3) ;
 
  ifstream paramFile5(fFileNamePar5) ; 

  Int_t i,j ;
  
  for(i = 0; i< 21; i++){
    for(j = 0; j< 3; j++){
      paramFile5 >> (*fParameters5)(i,j) ;
    }
  }
  paramFile5.close();
 
  ifstream paramFile100(fFileNamePar100) ; 
  
  Int_t l,k ;
  
  for(l = 0; l< 22; l++){
    for(k = 0; k< 3; k++){
      paramFile100 >> (*fParameters100)(l,k) ;
    }
  }
  paramFile100.close();
 
  ifstream paramFile(fFileNamePar100) ; 
  Int_t h,n;
  for(h = 0; h< 22; h++){
    for(n = 0; n< 3; n++){
      paramFile >> (*fParameters)(h,n) ;
    }
  }
  paramFile.close();

  fCluster = -1;
  fClusterrcpv = -1;
  fMatrixExtraRow = 0;

  //Calibration parameters Encal = C * E^2 + B * E + A  (E is the energy from cluster)
  fACalParameter = 0.0241  ;
  fBCalParameter = 1.0504  ;
  fCCalParameter = 0.000249 ;
 
  // fParameters->Print();
}
//_____________________________________________________________________________
void  AliPHOSPIDv1::GetAnalysisParameters(Float_t Cluster_En) 
{
  if(Cluster_En <=  fEnergyAnalysisCut){
    fPrincipal  = fPrincipal5;
    fParameters = fParameters5;
    fMatrixExtraRow = 0;
    GetClusterOption(Cluster_En,kFALSE) ;
  }
  else{
    fPrincipal  = fPrincipal100;
    fParameters = fParameters100;
    fMatrixExtraRow = 1;
    GetClusterOption(Cluster_En,kTRUE) ;
  }
}

//_____________________________________________________________________________
void  AliPHOSPIDv1::GetClusterOption(const Float_t Cluster_En, const Bool_t range) 
{

  // Gives the cluster energy range.
  // range = kFALSE Default analysis range from 0.5 to 5 GeV
  // range = kTRUE  analysis range from 0.5 to 100 GeV

  
  //Int_t cluster = -1 ;
  
  if((range == kFALSE)){
    if((Cluster_En > 0.3)&&(Cluster_En <= 1.0)){
      fCluster = 0 ;
      fClusterrcpv = 0 ;
    }
    if((Cluster_En > 1.0)&&(Cluster_En <= 2.0)){
      fCluster = 1 ;
      fClusterrcpv = 1 ;
    }
    if( Cluster_En > 2.0){
      fCluster = 2 ;
      fClusterrcpv = 2 ;
    }
  }
  else if(range == kTRUE){
    if((Cluster_En > 0.5 )&&(Cluster_En <= 20.0)) fCluster = 0 ;
    if((Cluster_En > 20.0)&&(Cluster_En <= 50.0)) fCluster = 1 ;
    if( Cluster_En > 50.0)                        fCluster = 2 ;
    if((Cluster_En > 5.0 )&&(Cluster_En <= 10.0)) fClusterrcpv = 0 ;
    if((Cluster_En > 10.0)&&(Cluster_En <= 20.0)) fClusterrcpv = 1 ;
    if((Cluster_En > 20.0)&&(Cluster_En <= 30.0)) fClusterrcpv = 2 ;
    if( Cluster_En > 30.0)                        fClusterrcpv = 3 ;
  }
  else {
    fCluster = -1 ;
    fClusterrcpv = -1;
    cout<<"Invalid Energy option"<<endl;
  }
  
  //return cluster;
}
//____________________________________________________________________________
Int_t  AliPHOSPIDv1::GetEffPurOption(const TString Eff_Pur) const
{

  // Looks for the Purity-Efficiency point (possible options "HIGH EFFICIENCY" 
  // "MEDIUM EFFICIENCY" "LOW EFFICIENCY" and 3 more options changing 
  // EFFICIENCY by PURITY)

  Int_t eff_pur = -1 ;

  if(Eff_Pur.Contains("HIGH EFFICIENCY") ||Eff_Pur.Contains("LOW PURITY") )
    eff_pur = 0 ;
  else if(Eff_Pur.Contains("MEDIUM EFFICIENCY") ||Eff_Pur.Contains("MEDIUM PURITY") ) 
    eff_pur = 1 ;
  else if(Eff_Pur.Contains("LOW EFFICIENCY")||Eff_Pur.Contains("HIGH PURITY") ) 
    eff_pur = 2 ;
  else{
    eff_pur = -1;
    cout<<"Invalid Efficiency-Purity option"<<endl;
    cout<<"Possible options: HIGH EFFICIENCY =    LOW PURITY"<<endl;
    cout<<"                MEDIUM EFFICIENCY = MEDIUM PURITY"<<endl;
    cout<<"                   LOW EFFICIENCY =   HIGH PURITY"<<endl;
  }

  return eff_pur;
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


//   gAlice->GetEvent(0) ;

//   //check, if the branch with name of this" already exits?
//   if (gAlice->TreeR()) {
//     TObjArray * lob = (TObjArray*)gAlice->TreeR()->GetListOfBranches() ;
//     TIter next(lob) ; 
//     TBranch * branch = 0 ;  
//     Bool_t phospidfound = kFALSE, pidfound = kFALSE ; 
    
//     TString taskName(GetName()) ; 
//     taskName.Remove(taskName.Index(Version())-1) ;
    
//     while ( (branch = (TBranch*)next()) && (!phospidfound || !pidfound) ) {
//       if ( (strcmp(branch->GetName(), "PHOSPID")==0) && (strcmp(branch->GetTitle(), taskName.Data())==0) ) 
// 	phospidfound = kTRUE ;
      
//       else if ( (strcmp(branch->GetName(), "AliPHOSPID")==0) && (strcmp(branch->GetTitle(), taskName.Data())==0) ) 
// 	pidfound = kTRUE ; 
//     }
    
//     if ( phospidfound || pidfound ) {
//       cerr << "WARNING: AliPHOSPIDv1::Exec -> RecParticles and/or PIDtMaker branch with name " 
// 	   << taskName.Data() << " already exits" << endl ;
//       return ; 
//     }       
//   }

//   Int_t nevents = (Int_t) gAlice->TreeE()->GetEntries() ;
//   Int_t ievent ;
//   AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;  

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
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
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments() ; 
  if ( !emcRecPoints || !cpvRecPoints || !trackSegments ) {
    cerr << "ERROR:  AliPHOSPIDv1::MakeRecParticles -> RecPoints or TrackSegments not found ! " << endl ; 
    abort() ; 
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
    
    AliPHOSRecPoint    * cpv = 0 ;
    if(ts->GetCpvIndex()>=0)
      cpv = (AliPHOSRecPoint *)   cpvRecPoints->At(ts->GetCpvIndex()) ;
    
    // Now set type (reconstructed) of the particle

    // Choose the cluster energy range
    
    // YK: check if (emc != 0) !!!
    if (!emc) {
      cerr << "ERROR:  AliPHOSPIDv1::MakeRecParticles -> emc("
	   <<ts->GetEmcIndex()<<") = "         <<emc<< endl;
      abort();
    }
    Float_t    e = emc->GetEnergy() ;   
    
    GetAnalysisParameters(e);// Gives value to fCluster, fClusterrcpv, fMatrixExtraRow, and to fPrincipal and fParameters depending on the energy.
    
    if((fCluster== -1)||(fClusterrcpv == -1)) continue ;
    
    Float_t  lambda[2] ;
    emc->GetElipsAxis(lambda) ;
    Float_t time =emc->GetTime() ;
    
    if((lambda[0]>0.01) && (lambda[1]>0.01) && time > 0.){
      
      // Loop of Efficiency-Purity (the 3 points of purity or efficiency are taken 
      // into account to set the particle identification)
      for(Int_t eff_pur = 0; eff_pur < 3 ; eff_pur++){
	
	// Looking at the CPV detector. If RCPV greater than CpvEmcDistance, 1st, 
	// 2nd or 3rd bit (depending on the efficiency-purity point )is set to 1 . 
	
	if(GetDistance(emc, cpv,  "R") > (*fParameters)(fClusterrcpv,eff_pur) )  
	  rp->SetPIDBit(eff_pur) ;
	
	// Looking the TOF. If TOF smaller than gate,  4th, 5th or 6th 
	// bit (depending on the efficiency-purity point )is set to 1             
	if(time< (*fParameters)(fCluster+3+fMatrixExtraRow,eff_pur))  
	  rp->SetPIDBit(eff_pur+3) ;		    
	
	// Looking PCA. Define and calculate the data (X), introduce in the function 
	// X2P that gives the components (P).  
	Float_t  Spher = 0. ;
	Float_t  Emaxdtotal = 0. ; 
	
	if((lambda[0]+lambda[1])!=0) Spher=fabs(lambda[0]-lambda[1])/(lambda[0]+lambda[1]); 
	
	Emaxdtotal=emc->GetMaximalEnergy()/emc->GetEnergy(); 
	
	fX[0] = lambda[0] ;  
	fX[1] = lambda[1] ; 
	fX[2] = emc->GetDispersion() ; 
	fX[3] = Spher ; 
	fX[4] = emc->GetMultiplicity() ;  
	fX[5] = Emaxdtotal ;  
	fX[6] = emc->GetCoreEnergy() ;  
	
	fPrincipal->X2P(fX,fP);
	
	//If we are inside the ellipse, 7th, 8th or 9th 
	// bit (depending on the efficiency-purity point )is set to 1 
	if(GetPrincipalSign(fP,fCluster+fMatrixExtraRow,eff_pur) == 1) 
	  rp->SetPIDBit(eff_pur+6) ;
	
      }
    }
    
    //Set momentum, energy and other parameters 
    Float_t  encal = CalibratedEnergy(e);
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
    index++ ; 
  }
  
}

//____________________________________________________________________________
void  AliPHOSPIDv1:: Print()
{
  // Print the parameters used for the particle type identification
    cout <<  "=============== AliPHOSPID1 ================" << endl ;
    cout <<  "Making PID "<< endl ;
//     cout <<  "    Headers file:               " << fHeaderFileName.Data() << endl ;
//     cout <<  "    RecPoints branch title:     " << fRecPointsTitle.Data() << endl ;
//     cout <<  "    TrackSegments Branch title: " << fTrackSegmentsTitle.Data() << endl ;
//     cout <<  "    RecParticles Branch title   " << fRecParticlesTitle.Data() << endl;

    cout <<  "    Pricipal analysis file from 0.5 to 5 " << fFileName5.Data() << endl;
    cout <<  "    Name of parameters file     "<<fFileNamePar5.Data() << endl ;
    cout <<  "    Matrix of Parameters: "<<endl;
    cout <<  "           3  Columns [High Eff-Low Pur,Medium Eff-Pur, Low Eff-High Pur]"<<endl;
    cout <<  "           21 Rows, each 3 [ RCPV, TOF, X_Center, Y_Center, A, B, Angle ]"<<endl;
    fParameters5->Print() ;

    cout <<  "    Pricipal analysis file from 5 to 100 " << fFileName100.Data() << endl;
    cout <<  "    Name of parameters file     "<<fFileNamePar100.Data() << endl ;
    cout <<  "    Matrix of Parameters: "<<endl;
    cout <<  "           3  Columns [High Eff-Low Pur,Medium Eff-Pur, Low Eff-High Pur]"<<endl;
    cout <<  "           22 Rows, [ 4 RCPV, 3 TOF, 3 X_Center, 3 Y_Center, 3 A, 3 B, 3 Angle ]"<<endl;
    fParameters100->Print() ;

    cout <<  "    Energy Calibration Parameters  A + B* E + C * E^2"<<endl;
    cout <<  "    E is the energy from the cluster "<<endl;
    cout <<  "           A = "<< fACalParameter  << endl;
    cout <<  "           B = "<< fBCalParameter  << endl;   
    cout <<  "           C = "<< fCCalParameter  << endl; 
    cout <<  "============================================" << endl ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::WriteRecParticles(Int_t event)
{
 
  AliPHOSGetter *gime = AliPHOSGetter::GetInstance() ; 

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
  TBranch * rpBranch = treeR->Branch("PHOSRP",&recParticles,bufferSize);
  rpBranch->SetTitle(BranchName());

  
  //second, pid
  Int_t splitlevel = 0 ; 
  AliPHOSPIDv1 * pid = this ;
  TBranch * pidBranch = treeR->Branch("AliPHOSPID","AliPHOSPIDv1",&pid,bufferSize,splitlevel);
  pidBranch->SetTitle(BranchName());
  
  rpBranch->Fill() ;
  pidBranch->Fill() ; 
  
  treeR->AutoSave() ; //Write(0,kOverwrite) ;  
  if(gAlice->TreeR()!=treeR){
    treeR->Delete();
  }
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



