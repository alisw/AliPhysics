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
#include "TSystem.h"
#include "TBenchmark.h"
#include "TMatrixD.h"
#include "TPrincipal.h"

// --- Standard library ---

//#include <Riostream.h>

// --- AliRoot header files ---

#include "AliGenerator.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSTrackSegment.h"
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
AliPHOSPIDv1::AliPHOSPIDv1(AliPHOSPIDv1 & pid ):AliPHOSPID(pid)
{ 
  InitParameters() ; 

  Init() ;
  fDefaultInit = kFALSE ; 

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

  delete [] fX ;    // Principal input 
  delete [] fP ;    // Principal components
  delete [] fPPi0 ; // Pi0 Principal components

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
  // gives the name of the current branch
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
  fPi0Analysis = kFALSE ;
  SetParameters() ; // fill the parameters matrix from parameters file
}

//____________________________________________________________________________
const Float_t  AliPHOSPIDv1::GetCpvtoEmcDistanceCut(const Float_t e, const TString Axis) const
{
  // Get CpvtoEmcDistance Cut depending on the cluster energy, axis and 
  // Purity-Efficiency point 

  Int_t i = -1;
  if      (Axis.Contains("X")) i = 1;
  else if (Axis.Contains("Z")) i = 2;
  else
    Error("GetCpvtoEmcDistanceCut"," Invalid axis option ");
   
  Float_t a = (*fParameters)(i,0) ;
  Float_t b = (*fParameters)(i,1) ;
  Float_t c = (*fParameters)(i,2) ;

  Float_t sig = a + TMath::Exp(b-c*e);
  return sig;
}
//____________________________________________________________________________

const Double_t  AliPHOSPIDv1::GetTimeGate(const Int_t effpur) const
{
  // Get TimeGate parameter depending on Purity-Efficiency point 
 
   if(effpur>2 || effpur<0)
    Error("GetTimeGate","Invalid Efficiency-Purity choice %d",effpur);
    return (*fParameters)(3,effpur) ; 

}
//_____________________________________________________________________________
const Float_t  AliPHOSPIDv1::GetDistance(AliPHOSEmcRecPoint * emc,AliPHOSRecPoint * cpv, Option_t *  Axis)const
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
const Int_t  AliPHOSPIDv1::GetCPVBit(AliPHOSEmcRecPoint * emc,AliPHOSRecPoint * cpv,const Int_t EffPur, const Float_t e) const
{
  if(EffPur>2 || EffPur<0)
    Error("GetCPVBit","Invalid Efficiency-Purity choice %d",EffPur);
  
  Float_t sigX = GetCpvtoEmcDistanceCut(e,"X");
  Float_t sigZ = GetCpvtoEmcDistanceCut(e,"Z");
  
  Float_t deltaX = TMath::Abs(GetDistance(emc, cpv,  "X"));
  Float_t deltaZ = TMath::Abs(GetDistance(emc, cpv,  "Z"));
       
  if((deltaX>sigX*(EffPur+1)) || (deltaZ>sigZ*(EffPur+1)))
    return 1;//Neutral
  else
    return 0;//Charged
  
}

//____________________________________________________________________________
const Double_t  AliPHOSPIDv1::GetCalibratedEnergy(const Float_t e) const
{
//      It calibrates Energy depending on the recpoint energy.
//      The energy of the reconstructed cluster is corrected with 
//      the formula A + B* E  + C* E^2, whose parameters where obtained 
//      through the study of the reconstructed energy distribution of 
//      monoenergetic photons.
 
  Double_t p[]={0.,0.,0.};
  Int_t i;
  for(i=0;i<3;i++) p[i]= (*fParameters)(0,i);
  Double_t  enerec = p[0] +  p[1]* e+ p[2] * e * e;
  return enerec ;

}
//____________________________________________________________________________
const Int_t  AliPHOSPIDv1::GetPrincipalBit(const Double_t* p ,const Int_t effpur, const Float_t e)const
{
  //Is the particle inside de PCA ellipse?

  Int_t    prinbit  = 0 ;
  Double_t a        = GetEllipseParameter("a", e); 
  Double_t b        = GetEllipseParameter("b", e);
  Double_t c        = GetEllipseParameter("c", e);
  Double_t xCenter = GetEllipseParameter("x0", e); 
  Double_t yCenter = GetEllipseParameter("y0", e);
  
  Double_t r = TMath::Power((p[0] - xCenter)/a,2) + 
      TMath::Power((p[1] - yCenter)/b,2) +
     c*(p[0] -  xCenter)*(p[1] - yCenter)/(a*b) ;
  //3 different ellipses defined
  if((effpur==2)&&(r <1./2.)) prinbit= 1;
  if((effpur==1)&&(r <2.   )) prinbit= 1;
  if((effpur==0)&&(r <9./2.)) prinbit= 1;

  if(r<0)
    Error("GetPrincipalBit", "Negative square? R=%f \n",r) ;

  return prinbit;

}
//____________________________________________________________________________
const Int_t  AliPHOSPIDv1::GetPrincipalPi0Bit(const Double_t* p, const Int_t effpur, const Float_t e)const
{
  //Is the particle inside de Pi0 PCA ellipse?

  Int_t    prinbit  = 0 ;
  Double_t a        = GetEllipseParameterPi0("a", e); 
  Double_t b        = GetEllipseParameterPi0("b", e);
  Double_t c        = GetEllipseParameterPi0("c", e);
  Double_t xCenter = GetEllipseParameterPi0("x0", e); 
  Double_t yCenter = GetEllipseParameterPi0("y0", e);
  
  Double_t r = TMath::Power((p[0] - xCenter)/a,2) + 
      TMath::Power((p[1] - yCenter)/b,2) +
      c*(p[0] -  xCenter)*(p[1] - yCenter)/(a*b) ;
  //3 different ellipses defined
  if((effpur==2)&&(r <1./2.)) prinbit= 1;
  if((effpur==1)&&(r <2.   )) prinbit= 1;
  if((effpur==0)&&(r <9./2.)) prinbit= 1;

  if(r<0)
    Error("GetPrincipalPi0Bit", "Negative square?") ;

  return prinbit;

}
//_____________________________________________________________________________
void  AliPHOSPIDv1::SetCpvtoEmcDistanceCutParameters(Float_t e, Int_t effpur, TString Axis,Float_t cut) 
{
  // Set the parameters to calculate Cpvto EmcDistanceCut depending on the cluster energy and 
  // Purity-Efficiency point 

  if(effpur>2 || effpur<0)
     Error("SetCpvtoEmcDistanceCutParameters","Invalid Efficiency-Purity choice %d",effpur);

  Int_t i = -1;
  if     (Axis.Contains("X")) i = 1;
  else if(Axis.Contains("Z")) i = 2;
  else
    Error("SetCpvtoEmcDistanceCutParameters"," Invalid axis option");
  
  (*fParameters)(i,effpur) = cut ;
}
//_____________________________________________________________________________
void  AliPHOSPIDv1::SetTimeGate(Int_t effpur, Float_t gate) 
{
  // Set the parameter TimeGate depending on the cluster energy and 
  // Purity-Efficiency point 
  if(effpur>2 || effpur<0)
    Error("SetTimeGate","Invalid Efficiency-Purity choice %d",effpur);
  
  (*fParameters)(3,effpur)= gate ; 
} 
//_____________________________________________________________________________
void  AliPHOSPIDv1::SetParameters() 
{
  // PCA : To do the Principal Components Analysis it is necessary 
  // the Principal file, which is opened here
  fX         = new double[7]; // Data for the PCA 
  fP         = new double[7]; // Eigenvalues of the PCA
  fPPi0      = new double[7]; // Eigenvalues of the Pi0 PCA

  // Read photon principals from the photon file
  
  fFileName     = "$ALICE_ROOT/PHOS/PCA8pa15_0.5-100.root" ; 
  TFile f( fFileName.Data(), "read" ) ;
  fPrincipal = dynamic_cast<TPrincipal*> (f.Get("principal")) ; 
  f.Close() ; 

  // Read pi0 principals from the pi0 file

  fFileNamePi0  = "$ALICE_ROOT/PHOS/PCA_pi0_40-120.root" ;
  TFile fPi0( fFileNamePi0.Data(), "read" ) ;
  fPrincipalPi0 = dynamic_cast<TPrincipal*> (fPi0.Get("principal")) ; 
  fPi0.Close() ;

  // Open parameters file and initialization of the Parameters matrix. 
  // In the File Parameters.dat are all the parameters. These are introduced 
  // in a matrix of 9x4  
  // 
  // All the parameters defined in this file are, in order of row: 
  // -CpvtoEmcDistanceCut (2 row (x and z) and 3 columns, each one depending 
  // on the parameter of the funtion that sets the cut in x or z.   
  // -TimeGate, 1 row and 3 columns (3 efficiency-purty cuts) 
  // -PCA, parameters of the functions that 
  // calculate the ellipse parameters, x0,y0,a, b, c. These 5 parameters 
  // (5 rows) depend on 4 parameters (columns). 
  // -Finally there is a row with the energy calibration parameters, 
  // 3 parameters. 

  fFileNamePar = gSystem->ExpandPathName("$ALICE_ROOT/PHOS/Parameters.dat");
  fParameters = new TMatrixD(14,4) ;
  const Int_t kmaxLeng=255;
  char string[kmaxLeng];

  // Open a text file with PID parameters
  FILE *fd = fopen(fFileNamePar.Data(),"r");
  if (!fd)
    Fatal("SetParameter","File %s with a PID parameters cannot be opened\n",
	  fFileNamePar.Data());

  Int_t i=0;
  // Read parameter file line-by-line and skip empty line and comments
  while (fgets(string,kmaxLeng,fd) != NULL) {
    if (string[0] == '\n' ) continue;
    if (string[0] == '!'  ) continue;
    sscanf(string, "%lf %lf %lf %lf",
	   &(*fParameters)(i,0), &(*fParameters)(i,1), 
	   &(*fParameters)(i,2), &(*fParameters)(i,3));
    i++;
  }
  fclose(fd);
}


//________________________________________________________________________
void  AliPHOSPIDv1::SetEllipseParameter(TString Param, Int_t i, Double_t par) 
{  
  // Set the parameter "i" that is needed to calculate the ellipse 
  // parameter "Param".
  
  Int_t p= -1;
  if     (Param.Contains("a")) p=4; 
  else if(Param.Contains("b")) p=5; 
  else if(Param.Contains("c")) p=6; 
  else if(Param.Contains("x0"))p=7; 
  else if(Param.Contains("y0"))p=8;
  if((i>4)||(i<0))
    Error("SetEllipseParameter", "No parameter with index %d", i) ; 
  else if(p==-1)
     Error("SetEllipseParameter", "No parameter with name %s", Param.Data() ) ; 
  else
    (*fParameters)(p,i) = par ;
} 
//________________________________________________________________________
void  AliPHOSPIDv1::SetEllipseParameterPi0(TString Param, Int_t i, Double_t par) 
{  
  // Set the parameter "i" that is needed to calculate the ellipse 
  // parameter "Param".
  if(!fPi0Analysis) Error("SetPi0EllipseParameter", "Pi 0 Analysis is off") ; 
  Int_t p= -1;
  if     (Param.Contains("a")) p=9; 
  else if(Param.Contains("b")) p=10; 
  else if(Param.Contains("c")) p=11; 
  else if(Param.Contains("x0"))p=12; 
  else if(Param.Contains("y0"))p=13;
  if((i>4)||(i<0))
    Error("SetPi0EllipseParameter", "No parameter with index %d", i) ; 
  else if(p==-1)
     Error("SetPi0EllipseParameter", "No parameter with name %s", Param.Data() ) ; 
  else
    (*fParameters)(p,i) = par ;
} 
//________________________________________________________________________
const Double_t  AliPHOSPIDv1::GetParameterToCalculateEllipse(const TString Param, const Int_t i) const
{ 
  // Get the parameter "i" that is needed to calculate the ellipse 
  // parameter "Param".

  Int_t p= -1;
  Double_t par = -1;

  if     (Param.Contains("a")) p=4; 
  else if(Param.Contains("b")) p=5; 
  else if(Param.Contains("c")) p=6; 
  else if(Param.Contains("x0"))p=7; 
  else if(Param.Contains("y0"))p=8;

  if((i>4)||(i<0))
    Error("GetParameterToCalculateEllipse", "No parameter with index", i) ; 
  else if(p==-1)
    Error("GetParameterToCalculateEllipse", "No parameter with name %s", Param.Data() ) ; 
  else
    par = (*fParameters)(p,i) ;
  
  return par;

} 
//____________________________________________________________________________
const Double_t  AliPHOSPIDv1::GetParameterToCalculatePi0Ellipse(const TString Param, const Int_t i) const
{ 
  // Get the parameter "i" that is needed to calculate the ellipse 
  // parameter "Param".

  if(!fPi0Analysis) Error("GetParameterToCalculatePi0Ellipse", "Pi 0 Analysis is off") ;

  Int_t p= -1;
  Double_t par = -1;

  if(Param.Contains("a")) p=9; 
  if(Param.Contains("b")) p=10; 
  if(Param.Contains("c")) p=11; 
  if(Param.Contains("x0"))p=12; 
  if(Param.Contains("y0"))p=13;

  if((i>4)||(i<0))
    Error("GetParameterToCalculatePi0Ellipse", "No parameter with index", i) ; 
  else if(p==-1)
    Error("GetParameterToCalculatePi0Ellipse", "No parameter with name %s", Param.Data() ) ; 
  else
    par = (*fParameters)(p,i) ;
  
  return par;

} 
//____________________________________________________________________________
void  AliPHOSPIDv1::SetCalibrationParameter(Int_t i,Double_t param) const
{
  (*fParameters)(0,i) = param ;
}
//____________________________________________________________________________
const Double_t  AliPHOSPIDv1::GetCalibrationParameter(const Int_t i) const 
{
  Float_t param = (*fParameters)(0,i);
  return param;
}
//____________________________________________________________________________
const Double_t  AliPHOSPIDv1::GetEllipseParameter(const TString Param,Float_t E) const 
{
  // Calculates the parameter Param of the ellipse
  
  Double_t p[4]={0.,0.,0.,0.};
  Double_t value = 0.0;
  Int_t i;

  if(Param.Contains("a")){
    for(i=0;i<4;i++)p[i]=(*fParameters)(4,i);
    if(E>70.)E=70.;
  }
  
  else if(Param.Contains("b")){
    for(i=0;i<4;i++)p[i]=(*fParameters)(5,i);
    if(E>70.)E=70.;
  }
  
  else if(Param.Contains("c"))
    for(i=0;i<4;i++)p[i]=(*fParameters)(6,i);
  
  else if(Param.Contains("x0")){
    for(i=0;i<4;i++)p[i]=(*fParameters)(7,i);
    if(E<1.)E=1.1;
  }
  else if(Param.Contains("y0"))
    for(i=0;i<4;i++)p[i]=(*fParameters)(8,i);
  
  value = p[0]/TMath::Sqrt(E)+p[1]*E+p[2]*E*E+p[3];
  return value;
}

//____________________________________________________________________________
// const Double_t  AliPHOSPIDv1::GetEllipseParameter(const TString Param,Float_t E) const 
// {
//   // Calculates the parameter Param of the pi0 ellipse
  
//   Double_t p[3]  = {0.,0.,0.};
//   Double_t value = 0.0;
//   Int_t    i;

//   if(Param.Contains("a"))
//     for(i=0;i<3;i++)p[i]=(*fParameters)(4,i);
//   else if(Param.Contains("b"))
//     for(i=0;i<3;i++)p[i]=(*fParameters)(5,i);
//   else if(Param.Contains("c"))
//     for(i=0;i<3;i++)p[i]=(*fParameters)(6,i);
//   else if(Param.Contains("x0"))
//     for(i=0;i<3;i++)p[i]=(*fParameters)(7,i);
//   else if(Param.Contains("y0"))
//     for(i=0;i<3;i++)p[i]=(*fParameters)(8,i);
  
//   value = p[0] + p[1]*E + p[2]*E*E;
//   return value;
// }
//____________________________________________________________________________
const Double_t  AliPHOSPIDv1::GetEllipseParameterPi0(const TString Param,Float_t E) const 
{
  // Calculates the parameter Param of the pi0 ellipse
  
  Double_t p[3]  = {0.,0.,0.};
  Double_t value = 0.0;
  Int_t    i;

  if(Param.Contains("a"))
    for(i=0;i<3;i++)p[i]=(*fParameters)(9,i);
  else if(Param.Contains("b"))
    for(i=0;i<3;i++)p[i]=(*fParameters)(10,i);
  else if(Param.Contains("c"))
    for(i=0;i<3;i++)p[i]=(*fParameters)(11,i);
  else if(Param.Contains("x0"))
    for(i=0;i<3;i++)p[i]=(*fParameters)(12,i);
  else if(Param.Contains("y0"))
    for(i=0;i<3;i++)p[i]=(*fParameters)(13,i);
  
  value = p[0] + p[1]*E + p[2]*E*E;
  return value;
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


  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  if(gime->BranchExists("RecParticles") )
    return ;
  Int_t nevents = gime->MaxEvent() ;  
  Int_t ievent ;

  for(ievent = 0; ievent < nevents; ievent++){
    gime->Event(ievent,"R") ;
 
    if(gime->TrackSegments() && //Skip events, where no track segments made
       gime->TrackSegments()->GetEntriesFast()) {
      MakeRecParticles() ;
      WriteRecParticles(ievent);
      if(strstr(option,"deb"))
        PrintRecParticles(option) ;
      //increment the total number of rec particles per run 
      fRecParticlesInRun+=gime->RecParticles(BranchName())->GetEntriesFast() ; 
    }
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSPID");
    Info("Exec", "took %f seconds for PID %f seconds per event", 
	 gBenchmark->GetCpuTime("PHOSPID"),  
	 gBenchmark->GetCpuTime("PHOSPID")/nevents) ;
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
    
    AliPHOSRecPoint    * cpv = 0 ;
    if(ts->GetCpvIndex()>=0)
      cpv = (AliPHOSRecPoint *)   cpvRecPoints->At(ts->GetCpvIndex()) ;
    
    // Now set type (reconstructed) of the particle

    // Choose the cluster energy range
    
    if (!emc) {
      Fatal("MakeRecParticles", "-> emc(%d) = %d", ts->GetEmcIndex(), emc ) ;
    }

    Float_t    e = emc->GetEnergy() ;   
    
    Float_t  lambda[2] ;
    emc->GetElipsAxis(lambda) ;
    
    if((lambda[0]>0.01) && (lambda[1]>0.01)){
      // Looking PCA. Define and calculate the data (X),
      // introduce in the function X2P that gives the components (P).  

      Float_t  spher = 0. ;
      Float_t  emaxdTotal = 0. ; 
      
      if((lambda[0]+lambda[1])!=0) 
	spher=fabs(lambda[0]-lambda[1])/(lambda[0]+lambda[1]); 
      
      emaxdTotal=emc->GetMaximalEnergy()/emc->GetEnergy(); 
      
      fX[0] = lambda[0] ;  
      fX[1] = lambda[1] ; 
      fX[2] = emc->GetDispersion() ; 
      fX[3] = spher ; 
      fX[4] = emc->GetMultiplicity() ;  
      fX[5] = emaxdTotal ;  
      fX[6] = emc->GetCoreEnergy() ;  
      
      fPrincipal->X2P(fX,fP);
      if(fPi0Analysis)
	fPrincipalPi0->X2P(fX,fPPi0);

    }
    else{
      fP[0]=-100.0;  //We do not accept clusters with 
      fP[1]=-100.0;  //one cell as a photon-like
      if(fPi0Analysis){
	fPPi0[0]=-100.0;
	fPPi0[1]=-100.0;
      }
    }
    
    Float_t time =emc->GetTime() ;
    
    // Loop of Efficiency-Purity (the 3 points of purity or efficiency 
    // are taken into account to set the particle identification)
    for(Int_t eff_pur = 0; eff_pur < 3 ; eff_pur++){
      
      // Looking at the CPV detector. If RCPV greater than CpvEmcDistance, 
      // 1st,2nd or 3rd bit (depending on the efficiency-purity point )
      // is set to 1
      if(GetCPVBit(emc, cpv, eff_pur,e) == 1 )  
	rp->SetPIDBit(eff_pur) ;
      
      // Looking the TOF. If TOF smaller than gate,  4th, 5th or 6th 
      // bit (depending on the efficiency-purity point )is set to 1             
      if(time< (*fParameters)(2,eff_pur)) 
	rp->SetPIDBit(eff_pur+3) ;		    
      
      //If we are inside the ellipse, 7th, 8th or 9th 
      // bit (depending on the efficiency-purity point )is set to 1 
      if(GetPrincipalBit(fP,eff_pur,e) == 1) 
	rp->SetPIDBit(eff_pur+6) ;

      //Pi0 analysis
      //If we are inside the ellipse, 10th, 11th or 12th 
      // bit (depending on the efficiency-purity point )is set to 1 
      if(fPi0Analysis){
	if(GetPrincipalPi0Bit(fPPi0,eff_pur,e) == 1) 
	  rp->SetPIDBit(eff_pur+9) ;
      }
    }
    
    
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
    index++ ; 
  }
  
}

//____________________________________________________________________________
void  AliPHOSPIDv1::Print()
{
  // Print the parameters used for the particle type identification

  TString message ; 
    message  = "\n=============== AliPHOSPID1 ================\n" ;
    message += "Making PID\n";
    message += "    Pricipal analysis file from 0.5 to 100 %s\n" ; 
    message += "    Name of parameters file     %s\n" ;
    message += "    Matrix of Parameters: 14x4\n" ;
    message += "        Energy Calibration  1x3 [3 parametres to calibrate energy: A + B* E + C * E^2]\n" ;
    message += "        RCPV 2x3 rows x and z, columns function cut parameters\n" ;
    message += "        TOF  1x3 [High Eff-Low Pur,Medium Eff-Pur, Low Eff-High Pur]\n" ;
    message += "        PCA  5x4 [5 ellipse parametres and 4 parametres to calculate them: A/Sqrt(E) + B* E + C * E^2 + D]\n" ;
    message += "    Pi0 PCA  5x3 [5 ellipse parametres and 3 parametres to calculate them: A + B* E + C * E^2]\n" ;
    Info("Print", message.Data(), fFileName.Data(), fFileNamePar.Data() ) ; 
    fParameters->Print() ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::WriteRecParticles(Int_t event)
{
  // writes the reconstructed particles to file
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

  //account correction to the position of IP
  Float_t xo,yo,zo ; //Coordinates of the origin
  gAlice->Generator()->GetOrigin(xo,yo,zo) ;
  TVector3 origin(xo,yo,zo);
  dir = dir - origin ;
  dir.SetMag(1.) ;
  return dir ;  
}
//____________________________________________________________________________
void AliPHOSPIDv1::PrintRecParticles(Option_t * option)
{
  // Print table of reconstructed particles

  AliPHOSGetter *gime = AliPHOSGetter::GetInstance() ; 

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



