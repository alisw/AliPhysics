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


///////////////////////////////////////////////////////////////////////
//                                                                   //
//  		AliZDCv4 --- new ZDC geometry		     	     //
//  	    with both ZDC arms geometry implemented 		     //
//                                                                   //  
///////////////////////////////////////////////////////////////////////

// --- Standard libraries
#include "stdio.h"

// --- ROOT system
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoShape.h>
#include <TGeoScaledShape.h>
#include <TGeoCompositeShape.h>
#include <TParticle.h>

// --- AliRoot classes
#include "AliLog.h"
#include "AliConst.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliZDCv4.h"
#include "AliMC.h"
#include "AliMCParticle.h"
#include "AliTrackReference.h"
 
class  AliZDCHit;
class  AliPDG;
class  AliDetector;
 
 
ClassImp(AliZDCv4)

//_____________________________________________________________________________
AliZDCv4::AliZDCv4() : 
  AliZDC(),
  fMedSensF1(0),
  fMedSensF2(0),
  fMedSensZP(0),
  fMedSensZN(0),
  fMedSensZEM(0),
  fMedSensGR(0),
  fMedSensPI(0),
  fMedSensTDI(0),
  fMedSensVColl(0),
  fMedSensLumi(0),
  fNalfan(0),
  fNalfap(0),
  fNben(0),  
  fNbep(0),
  fZEMLength(0),
  fpLostITC(0), 
  fpLostD1C(0), 
  fpcVCollC(0),
  fpDetectedC(0),
  fnDetectedC(0),
  fpLostITA(0), 
  fpLostD1A(0), 
  fpLostTDI(0), 
  fpcVCollA(0),
  fpDetectedA(0),
  fnDetectedA(0),
  fVCollSideCAperture(7./2.),
  fVCollSideCApertureNeg(7./2.),
  fVCollSideCCentreY(0.),
  fTCDDAperturePos(2.0),
  fTCDDApertureNeg(2.0),
  fTDIAperturePos(5.5),
  fTDIApertureNeg(5.5),
  fLumiLength(15.),
  fSwitchOnTrackRef(kFALSE)
{
  //
  // Default constructor for Zero Degree Calorimeter
  //
  for(Int_t i=0; i<3; i++){
     fDimZN[i] = fDimZP[i] = 0.;
     fPosZNC[i] = fPosZNA[i] = fPosZPC[i]= fPosZPA[i] = fPosZEM[i] = 0.;
     fFibZN[i] = fFibZP[i] = 0.;
  }
}
 
//_____________________________________________________________________________
AliZDCv4::AliZDCv4(const char *name, const char *title) : 
  AliZDC(name,title),
  fMedSensF1(0),
  fMedSensF2(0),
  fMedSensZP(0),
  fMedSensZN(0),
  fMedSensZEM(0),
  fMedSensGR(0),
  fMedSensPI(0),
  fMedSensTDI(0),
  fMedSensVColl(0),
  fMedSensLumi(0),
  fNalfan(90),
  fNalfap(90),
  fNben(18),  
  fNbep(28), 
  fZEMLength(0),
  fpLostITC(0), 
  fpLostD1C(0), 
  fpcVCollC(0),
  fpDetectedC(0),
  fnDetectedC(0),
  fpLostITA(0), 
  fpLostD1A(0), 
  fpLostTDI(0), 
  fpcVCollA(0),
  fpDetectedA(0),
  fnDetectedA(0),
  fVCollSideCAperture(7./2.),
  fVCollSideCApertureNeg(7./2.),
  fVCollSideCCentreY(0.),
  fTCDDAperturePos(2.0),
  fTCDDApertureNeg(2.0),
  fTDIAperturePos(5.5),
  fTDIApertureNeg(5.5),
  fLumiLength(15.),
  fSwitchOnTrackRef(kFALSE)  
{
  //
  // Standard constructor for Zero Degree Calorimeter 
  //
  //
  // Check that DIPO, ABSO, DIPO and SHIL is there (otherwise tracking is wrong!!!)
  
  AliModule* pipe=gAlice->GetModule("PIPE");
  AliModule* abso=gAlice->GetModule("ABSO");
  AliModule* dipo=gAlice->GetModule("DIPO");
  AliModule* shil=gAlice->GetModule("SHIL");
  if((!pipe) || (!abso) || (!dipo) || (!shil)) {
    Error("Constructor","ZDC needs PIPE, ABSO, DIPO and SHIL!!!\n");
    exit(1);
  } 
  //
  Int_t ip,jp,kp;
  for(ip=0; ip<4; ip++){
     for(kp=0; kp<fNalfap; kp++){
        for(jp=0; jp<fNbep; jp++){
           fTablep[ip][kp][jp] = 0;
        } 
     }
  }
  Int_t in,jn,kn;
  for(in=0; in<4; in++){
     for(kn=0; kn<fNalfan; kn++){
        for(jn=0; jn<fNben; jn++){
           fTablen[in][kn][jn] = 0;
        } 
     }
  }
  //
  // Parameters for hadronic calorimeters geometry
  // Positions updated after post-installation measurements
  fDimZN[0] = 3.52;
  fDimZN[1] = 3.52;
  fDimZN[2] = 50.;  
  fDimZP[0] = 11.2;
  fDimZP[1] = 6.;
  fDimZP[2] = 75.;    
  fPosZNC[0] = 0.;
  fPosZNC[1] = 0.;
  fPosZNC[2] = -11397.3+136; 
  fPosZPC[0] = 24.35;
  fPosZPC[1] = 0.;
  fPosZPC[2] = -11389.3+136; 
  fPosZNA[0] = 0.;
  fPosZNA[1] = 0.;
  fPosZNA[2] = 11395.8-136;  
  fPosZPA[0] = 24.35;
  fPosZPA[1] = 0.;
  fPosZPA[2] = 11387.8-136; 
  fFibZN[0] = 0.;
  fFibZN[1] = 0.01825;
  fFibZN[2] = 50.;
  fFibZP[0] = 0.;
  fFibZP[1] = 0.0275;
  fFibZP[2] = 75.;
  // Parameters for EM calorimeter geometry
  fPosZEM[0] = 8.5;
  fPosZEM[1] = 0.;
  fPosZEM[2] = 735.;
  Float_t kDimZEMPb  = 0.15*(TMath::Sqrt(2.));  // z-dimension of the Pb slice
  Float_t kDimZEMAir = 0.001; 			// scotch
  Float_t kFibRadZEM = 0.0315; 			// External fiber radius (including cladding)
  Int_t   kDivZEM[3] = {92, 0, 20}; 		// Divisions for EM detector
  Float_t kDimZEM0 = 2*kDivZEM[2]*(kDimZEMPb+kDimZEMAir+kFibRadZEM*(TMath::Sqrt(2.)));
  fZEMLength = kDimZEM0;
  
}
 
//_____________________________________________________________________________
void AliZDCv4::CreateGeometry()
{
  //
  // Create the geometry for the Zero Degree Calorimeter version 2
  //* Initialize COMMON block ZDC_CGEOM
  //*

  CreateBeamLine();
  CreateZDC();
}
  
//_____________________________________________________________________________
void AliZDCv4::CreateBeamLine()
{
  //
  // Create the beam line elements
  //
  if(fOnlyZEM) printf("\n  Only ZEM configuration requested: no side-C beam pipe, no side-A hadronic ZDCs\n\n");
  
  Double_t zd1=0., zd2=0., zCorrDip=0., zInnTrip=0., zD1=0.;
  Double_t tubpar[3]={0.,0.,0}, boxpar[3]={0.,0.,0};
  Double_t tubspar[5]={0.,0.,0.,0.,0.};
  Double_t conpar[15];
  for(int i=0; i<15; i++) conpar[i]=0.;

  //-- rotation matrices for the legs
  Int_t irotpipe1, irotpipe2;
  TVirtualMC::GetMC()->Matrix(irotpipe1,90.-1.0027,0.,90.,90.,1.0027,180.);      
  TVirtualMC::GetMC()->Matrix(irotpipe2,90.+1.0027,0.,90.,90.,1.0027,0.);

  Int_t *idtmed = fIdtmed->GetArray();
      
  ////////////////////////////////////////////////////////////////
  //								//
  //                SIDE C - RB26 (dimuon side)			//
  //								//
  ////////////////////////////////////////////////////////////////
  
//if(!fOnlyZEM){  
  // -- Mother of the ZDCs (Vacuum PCON)
  //zd1 = 1921.6;
  //const Double_t kZComDip = 1972.5;
  // New -> to accomodate AD detector (by A. Morsch)
  zd1 = 1947.2;
  const Double_t kZComDip = 1974.0;
  //
  conpar[ 0] = 0.;
  conpar[ 1] = 360.;
  conpar[ 2] = 4.;      // Num radius specifications: 4
  conpar[ 3] = -13500.; // (1) end of mother vol
  conpar[ 4] = 0.;
  conpar[ 5] = 55.;
  conpar[ 6] = -kZComDip; // (2) Beginning of Compensator Dipole
  conpar[ 7] = 0.;
  conpar[ 8] = 55.;
  conpar[ 9] = -kZComDip; // (3) Beginning of Compensator Dipole
  conpar[10] = 0.;
  conpar[11] = 6.7/2.;
  conpar[12] = -zd1;    // (4) Beginning of ZDCC mother volume
  conpar[13] = 0.;
  conpar[14] = 6.7/2.; 
  TVirtualMC::GetMC()->Gsvolu("ZDCC", "PCON", idtmed[10], conpar, 15);
  TVirtualMC::GetMC()->Gspos("ZDCC", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  

  // -- BEAM PIPE from compensator dipole to the beginning of D1) 
  tubpar[0] = 6.3/2.;
  tubpar[1] = 6.7/2.;
  // From beginning of ZDC volumes to beginning of D1
  tubpar[2] = (5838.3-zd1)/2.;
  TVirtualMC::GetMC()->Gsvolu("QT01", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT01", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("  QT01 TUBE pipe from z = %1.2f to z = %1.2f (D1 begin)\n",-zd1,-2*tubpar[2]-zd1);
  
  //-- BEAM PIPE from the end of D1 to the beginning of D2) 
  
  //-- FROM MAGNETIC BEGINNING OF D1 TO MAGNETIC END OF D1
  //-- 	Cylindrical pipe (r = 3.47) + conical flare  
  // -> Beginning of D1
  zd1 += 2.*tubpar[2];
  
  tubpar[0] = 6.94/2.;
  tubpar[1] = 7.34/2.;
  tubpar[2] = (6909.8-zd1)/2.;
  TVirtualMC::GetMC()->Gsvolu("QT02", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT02", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT02 TUBE pipe from z = %1.2f to z = %1.2f (D1 magnetic end)\n",-zd1,-2*tubpar[2]-zd1);

  zd1 += 2.*tubpar[2];
  
  tubpar[0] = 8./2.;
  tubpar[1] = 8.6/2.;
  tubpar[2] = (6958.3-zd1)/2.;
  TVirtualMC::GetMC()->Gsvolu("QT0B", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT0B", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT0B TUBE pipe from z = %1.2f to z = %1.2f \n",-zd1,-2*tubpar[2]-zd1);
 
  zd1 += 2.*tubpar[2];
  
  tubpar[0] = 9./2.;
  tubpar[1] = 9.6/2.;
  tubpar[2] = (7022.8-zd1)/2.;
  TVirtualMC::GetMC()->Gsvolu("QT03", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT03", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT03 TUBE pipe from z = %1.2f to z = %1.2f (D1 end)\n",-zd1,-2*tubpar[2]-zd1);

  zd1 += 2.*tubpar[2];
  
  conpar[0] = 39.2/2.;
  conpar[1] = 18./2.;
  conpar[2] = 18.6/2.;
  conpar[3] = 9./2.;
  conpar[4] = 9.6/2.;
  TVirtualMC::GetMC()->Gsvolu("QC01", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QC01", 1, "ZDCC", 0., 0., -conpar[0]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QC01 CONE pipe from z = %1.2f to z= %1.2f (VCTCQ-I)\n",-zd1,-2*conpar[0]-zd1);
  
  zd1 += conpar[0] * 2.;
  
  // ******************************************************
  // N.B.-> according to last vacuum layout 
  // private communication by D. Macina, mail 27/1/2009
  // updated to new ZDC installation (Janiary 2012) 
  // ****************************************************** 
  // 2nd section of    VCTCQ+VAMTF+TCLIA+VAMTF+1st part of VCTCP
  Float_t totLength1 = 160.8 + 78. + 148. + 78. + 9.3;
  //
  tubpar[0] = 18.6/2.;
  tubpar[1] = 7.6/2.;
  tubpar[2] = totLength1/2.;
//  TVirtualMC::GetMC()->Gsvolu("QE01", "ELTU", idtmed[7], tubpar, 3);  
  // temporary replace with a scaled tube (AG)
  TGeoTube *tubeQE01 = new TGeoTube(0.,tubpar[0],tubpar[2]);
  TGeoScale *scaleQE01 = new TGeoScale(1., tubpar[1]/tubpar[0], 1.);
  TGeoScaledShape *sshapeQE01 = new TGeoScaledShape(tubeQE01, scaleQE01);
  new TGeoVolume("QE01", sshapeQE01, gGeoManager->GetMedium(idtmed[7]));

  tubpar[0] = 18.0/2.;
  tubpar[1] = 7.0/2.;
  tubpar[2] = totLength1/2.;
//  TVirtualMC::GetMC()->Gsvolu("QE02", "ELTU", idtmed[10], tubpar, 3);  
  // temporary replace with a scaled tube (AG)
  TGeoTube *tubeQE02 = new TGeoTube(0.,tubpar[0],tubpar[2]);
  TGeoScale *scaleQE02 = new TGeoScale(1., tubpar[1]/tubpar[0], 1.);
  TGeoScaledShape *sshapeQE02 = new TGeoScaledShape(tubeQE02, scaleQE02);
  new TGeoVolume("QE02", sshapeQE02, gGeoManager->GetMedium(idtmed[10]));

  TVirtualMC::GetMC()->Gspos("QE01", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("QE02", 1, "QE01", 0., 0., 0., 0, "ONLY");  
  // Ch.debug
  //printf("	QE01 ELTU from z = %1.2f to z = %1.2f (VCTCQ-II+VAMTF+TCLIA+VAMTF+VCTCP-I)\n",-zd1,-2*tubpar[2]-zd1);
  
  // TCLIA collimator jaws (defined ONLY if fVCollAperture<3.5!)
  if(fVCollSideCAperture<3.5){
    boxpar[0] = 5.4/2.;
    boxpar[1] = (3.5-fVCollSideCAperture-fVCollSideCCentreY-0.7)/2.;
    if(boxpar[1]<0.) boxpar[1]=0.;
    boxpar[2] = 124.4/2.;
    printf("  AliZDCv4 -> C side injection collimator (TCLIA) jaws: apertures +%1.2f/-%1.2f center %1.2f [cm]\n", 
    	fVCollSideCAperture, fVCollSideCApertureNeg,fVCollSideCCentreY);
    TVirtualMC::GetMC()->Gsvolu("QCVC" , "BOX ", idtmed[14], boxpar, 3); 
    TVirtualMC::GetMC()->Gspos("QCVC", 1, "QE02", -boxpar[0],  fVCollSideCAperture+fVCollSideCCentreY+boxpar[1], -totLength1/2.+160.8+78.+148./2., 0, "ONLY");  
    TVirtualMC::GetMC()->Gspos("QCVC", 2, "QE02", -boxpar[0], -fVCollSideCApertureNeg+fVCollSideCCentreY-boxpar[1], -totLength1/2.+160.8+78.+148./2., 0, "ONLY");  
  }
  
  zd1 += tubpar[2] * 2.;
  
  // 2nd part of VCTCP
  conpar[0] = 31.5/2.;
  conpar[1] = 21.27/2.;
  conpar[2] = 21.87/2.;
  conpar[3] = 18.0/2.;
  conpar[4] = 18.6/2.;
  TVirtualMC::GetMC()->Gsvolu("QC02", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QC02", 1, "ZDCC", 0., 0., -conpar[0]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QC02 CONE pipe from z = %1.2f to z= %1.2f (VCTCP-II)\n",-zd1,-2*conpar[0]-zd1);
  
  zd1 += conpar[0] * 2.;

  // 3rd section of VCTCP+VCDWC+VMLGB	
  //Float_t totLenght2 = 9.2 + 530.5+40.;
  Float_t totLenght2 = (8373.3-zd1);
  tubpar[0] = 21.2/2.;
  tubpar[1] = 21.9/2.;
  tubpar[2] = totLenght2/2.;
  TVirtualMC::GetMC()->Gsvolu("QT04", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT04", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT04 TUBE pipe from z = %1.2f to z= %1.2f (VCTCP-III)\n",-zd1,-2*tubpar[2]-zd1);
  
  zd1 += tubpar[2] * 2.;
  
  // First part of VCTCD
  // skewed transition cone from ID=212.7 mm to ID=797 mm
  conpar[0] = 121./2.;
  conpar[1] = 79.7/2.;
  conpar[2] = 81.3/2.;
  conpar[3] = 21.27/2.;
  conpar[4] = 21.87/2.;
  TVirtualMC::GetMC()->Gsvolu("QC03", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QC03", 1, "ZDCC", 0., 0., -conpar[0]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QC03 CONE pipe from z = %1.2f to z = %1.2f (VCTCD-I)\n",-zd1,-2*conpar[0]-zd1);
  
  zd1 += 2.*conpar[0];
  
  // VCDGB + 1st part of VCTCH
  // Modified according to 2012 ZDC installation
  tubpar[0] = 79.7/2.;
  tubpar[1] = 81.3/2.;
  tubpar[2] = (5*475.2+97.-136)/2.;
  TVirtualMC::GetMC()->Gsvolu("QT05", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT05", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT05 TUBE pipe from z = %1.2f to z = %1.2f (VCDGB+VCTCH-I)\n",-zd1,-2*tubpar[2]-zd1);
  
  zd1 += 2.*tubpar[2];
     
  // 2nd part of VCTCH
  // Transition from ID=797 mm to ID=196 mm:
  // in order to simulate the thin window opened in the transition cone
  // we divide the transition cone in three cones:
  // (1) 8 mm thick (2) 3 mm thick (3) the third 8 mm thick
  
  // (1) 8 mm thick
  conpar[0] = 9.09/2.; // 15 degree
  conpar[1] = 74.82868/2.;
  conpar[2] = 76.42868/2.; // thickness 8 mm 
  conpar[3] = 79.7/2.;
  conpar[4] = 81.3/2.; // thickness 8 mm  
  TVirtualMC::GetMC()->Gsvolu("QC04", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QC04", 1, "ZDCC", 0., 0., -conpar[0]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QC04 CONE pipe from z = %1.2f to z = %1.2f (VCTCH-II)\n",-zd1,-2*conpar[0]-zd1);

  zd1 += 2.*conpar[0];  

  // (2) 3 mm thick
  conpar[0] = 96.2/2.; // 15 degree
  conpar[1] = 23.19588/2.;
  conpar[2] = 23.79588/2.; // thickness 3 mm 
  conpar[3] = 74.82868/2.;
  conpar[4] = 75.42868/2.; // thickness 3 mm  
  TVirtualMC::GetMC()->Gsvolu("QC05", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QC05", 1, "ZDCC", 0., 0., -conpar[0]-zd1, 0, "ONLY");  
  // Ch.debug
  //printf("	QC05 CONE pipe from z = %1.2f to z = %1.2f (VCTCH-III)\n",-zd1,-2*conpar[0]-zd1);

  zd1 += 2.*conpar[0];
  
  // (3) 8 mm thick
  conpar[0] = 6.71/2.; // 15 degree
  conpar[1] = 19.6/2.;
  conpar[2] = 21.2/2.;// thickness 8 mm 
  conpar[3] = 23.19588/2.;
  conpar[4] = 24.79588/2.;// thickness 8 mm 
  TVirtualMC::GetMC()->Gsvolu("QC06", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QC06", 1, "ZDCC", 0., 0., -conpar[0]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QC06 CONE pipe from z = %1.2f to z = %1.2f (VCTCH-III)\n",-zd1,-2*conpar[0]-zd1);

  zd1 += 2.*conpar[0];
  
  // VMZAR (5 volumes)  
  tubpar[0] = 20.2/2.;
  tubpar[1] = 20.6/2.;
  tubpar[2] = 2.15/2.;
  TVirtualMC::GetMC()->Gsvolu("QT06", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT06", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT06 TUBE pipe from z = %1.2f to z = %1.2f (VMZAR-I)\n",-zd1,-2*tubpar[2]-zd1);

  zd1 += 2.*tubpar[2];
  
  conpar[0] = 6.9/2.;
  conpar[1] = 23.9/2.;
  conpar[2] = 24.3/2.;
  conpar[3] = 20.2/2.;
  conpar[4] = 20.6/2.;
  TVirtualMC::GetMC()->Gsvolu("QC07", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QC07", 1, "ZDCC", 0., 0., -conpar[0]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QC07 CONE pipe from z = %1.2f to z = %1.2f (VMZAR-II)\n",-zd1,-2*conpar[0]-zd1);

  zd1 += 2.*conpar[0];

  tubpar[0] = 23.9/2.;
  tubpar[1] = 25.5/2.;
  tubpar[2] = 17.0/2.;
  TVirtualMC::GetMC()->Gsvolu("QT07", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT07", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT07 TUBE pipe from z = %1.2f to z = %1.2f (VMZAR-III)\n",-zd1,-2*tubpar[2]-zd1);
 
  zd1 += 2.*tubpar[2];
  
  conpar[0] = 6.9/2.;
  conpar[1] = 20.2/2.;
  conpar[2] = 20.6/2.;
  conpar[3] = 23.9/2.;
  conpar[4] = 24.3/2.;
  TVirtualMC::GetMC()->Gsvolu("QC08", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QC08", 1, "ZDCC", 0., 0., -conpar[0]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QC08 CONE pipe from z = %1.2f to z = %1.2f (VMZAR-IV)\n",-zd1,-2*conpar[0]-zd1);

  zd1 += 2.*conpar[0];
  
  tubpar[0] = 20.2/2.;
  tubpar[1] = 20.6/2.;
  tubpar[2] = 2.15/2.;
  TVirtualMC::GetMC()->Gsvolu("QT08", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT08", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT08 TUBE pipe from z = %1.2f to z = %1.2f (VMZAR-V)\n",-zd1,-2*tubpar[2]-zd1);

  zd1 += 2.*tubpar[2];
  
  // Flange (ID=196 mm)(last part of VMZAR and first part of VCTYB)
  tubpar[0] = 19.6/2.;
  tubpar[1] = 25.3/2.;
  tubpar[2] = 4.9/2.;
  TVirtualMC::GetMC()->Gsvolu("QT09", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT09", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT09 TUBE pipe from z = %1.2f to z = %1.2f (VMZAR-VI+VCTYB-I)\n",-zd1,-2*tubpar[2]-zd1);
 
  zd1 += 2.*tubpar[2];
  // Ch.debug
  //printf("	Beginning of VCTYB volume @ z = %1.2f \n",-zd1);
  
  // simulation of the trousers (VCTYB)     
  tubpar[0] = 19.6/2.;
  tubpar[1] = 20.0/2.;
  tubpar[2] = 3.9/2.;
  TVirtualMC::GetMC()->Gsvolu("QT10", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT10", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT10 TUBE pipe from z = %1.2f to z = %1.2f (VCTYB-II)\n",-zd1,-2*tubpar[2]-zd1);

  zd1 += 2.*tubpar[2];

  // transition cone from ID=196. to ID=216.6
  conpar[0] = 32.55/2.;
  conpar[1] = 21.66/2.;
  conpar[2] = 22.06/2.;
  conpar[3] = 19.6/2.;
  conpar[4] = 20.0/2.;
  TVirtualMC::GetMC()->Gsvolu("QC09", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QC09", 1, "ZDCC", 0., 0., -conpar[0]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QC09 CONE pipe from z = %1.2f to z= %1.2f\n",-zd1,-2*conpar[0]-zd1);

  zd1 += 2.*conpar[0]; 
  
  // tube  
  tubpar[0] = 21.66/2.;
  tubpar[1] = 22.06/2.;
  tubpar[2] = 28.6/2.;
  TVirtualMC::GetMC()->Gsvolu("QT11", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT11", 1, "ZDCC", 0., 0., -tubpar[2]-zd1, 0, "ONLY");
  // Ch.debug
  //printf("	QT11 TUBE pipe from z = %1.2f to z= %1.2f\n",-zd1,-2*tubpar[2]-zd1);

  zd1 += 2.*tubpar[2];
  // Ch.debug
  //printf("	Beginning of C side recombination chamber @ z = %f \n",-zd1);

  // --------------------------------------------------------
  // RECOMBINATION CHAMBER IMPLEMENTED USING TGeo CLASSES!!!!
  // author: Chiara (August 2008)
  // --------------------------------------------------------
  // TRANSFORMATION MATRICES
  // Combi transformation: 
  Double_t dx = -3.970000;
  Double_t dy = 0.000000;
  Double_t dz = 0.0;
  // Rotation: 
  Double_t thx = 84.989100;   Double_t phx = 180.000000;
  Double_t thy = 90.000000;   Double_t phy = 90.000000;
  Double_t thz = 185.010900;  Double_t phz = 0.000000;
  TGeoRotation *rotMatrix1c = new TGeoRotation("c",thx,phx,thy,phy,thz,phz);
  // Combi transformation: 
  dx = -3.970000;
  dy = 0.000000;
  dz = 0.0;
  TGeoCombiTrans *rotMatrix2c = new TGeoCombiTrans("ZDCC_c1", dx,dy,dz,rotMatrix1c);
  rotMatrix2c->RegisterYourself();
  // Combi transformation: 
  dx = 3.970000;
  dy = 0.000000;
  dz = 0.0;
  // Rotation: 
  thx = 95.010900;   phx = 180.000000;
  thy = 90.000000;   phy = 90.000000;
  thz = 180.-5.010900;    phz = 0.000000;
  TGeoRotation *rotMatrix3c = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans *rotMatrix4c = new TGeoCombiTrans("ZDCC_c2", dx,dy,dz,rotMatrix3c);
  rotMatrix4c->RegisterYourself();

  // VOLUMES DEFINITION
  // Volume: ZDCC
  TGeoVolume *pZDCC = gGeoManager->GetVolume("ZDCC");
  
  conpar[0] = (90.1-0.95-0.26-0.0085)/2.;
  conpar[1] = 0.0/2.;
  conpar[2] = 21.6/2.;
  conpar[3] = 0.0/2.;
  conpar[4] = 5.8/2.;
  new TGeoCone("QCLext", conpar[0],conpar[1],conpar[2],conpar[3],conpar[4]);
  
  conpar[0] = (90.1-0.95-0.26-0.0085)/2.;
  conpar[1] = 0.0/2.;
  conpar[2] = 21.2/2.;
  conpar[3] = 0.0/2.;
  conpar[4] = 5.4/2.;
  new TGeoCone("QCLint", conpar[0],conpar[1],conpar[2],conpar[3],conpar[4]);

  // Outer trousers
  TGeoCompositeShape *pOutTrousersC = new TGeoCompositeShape("outTrousersC", "QCLext:ZDCC_c1+QCLext:ZDCC_c2");
  
  // Volume: QCLext
  TGeoMedium *medZDCFe = gGeoManager->GetMedium("ZDC_ZIRON");
  TGeoVolume *pQCLext = new TGeoVolume("QCLext",pOutTrousersC, medZDCFe);
  pQCLext->SetLineColor(kGreen);
  pQCLext->SetVisLeaves(kTRUE);
  //
  TGeoTranslation *tr1c = new TGeoTranslation(0., 0., (Double_t) -conpar[0]-0.95-zd1);
  //printf("	C side recombination chamber from z = %1.2f to z= %1.2f\n",-zd1,-2*conpar[0]-0.95-zd1);
  //
  pZDCC->AddNode(pQCLext, 1, tr1c);
  // Inner trousers
  TGeoCompositeShape *pIntTrousersC = new TGeoCompositeShape("intTrousersC", "QCLint:ZDCC_c1+QCLint:ZDCC_c2");
  // Volume: QCLint
  TGeoMedium *medZDCvoid = gGeoManager->GetMedium("ZDC_ZVOID");
  TGeoVolume *pQCLint = new TGeoVolume("QCLint",pIntTrousersC, medZDCvoid);
  pQCLint->SetLineColor(kTeal);
  pQCLint->SetVisLeaves(kTRUE);
  pQCLext->AddNode(pQCLint, 1);
    
  zd1 += 90.1;
  Double_t offset = 0.5;
  zd1 = zd1+offset;
  
  //  second section : 2 tubes (ID = 54. OD = 58.)  
  tubpar[0] = 5.4/2.;
  tubpar[1] = 5.8/2.;
  tubpar[2] = 40.0/2.;
  TVirtualMC::GetMC()->Gsvolu("QT12", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT12", 1, "ZDCC", -15.8/2., 0., -tubpar[2]-zd1, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("QT12", 2, "ZDCC",  15.8/2., 0., -tubpar[2]-zd1, 0, "ONLY");  
  // Ch.debug
  //printf("	QT12 TUBE from z = %1.2f to z = %1.2f (separate beam pipes)\n",-zd1,-2*tubpar[2]-zd1);
  
  zd1 += 2.*tubpar[2];
  
  // transition x2zdc to recombination chamber : skewed cone  
  conpar[0] = (10.-0.2-offset)/2.;
  conpar[1] = 6.3/2.;
  conpar[2] = 7.0/2.;
  conpar[3] = 5.4/2.;
  conpar[4] = 5.8/2.;
  TVirtualMC::GetMC()->Gsvolu("QC10", "CONE", idtmed[7], conpar, 5); 
  TVirtualMC::GetMC()->Gspos("QC10", 1, "ZDCC", -7.9-0.175, 0., -conpar[0]-0.1-zd1, irotpipe1, "ONLY");
  TVirtualMC::GetMC()->Gspos("QC10", 2, "ZDCC", 7.9+0.175, 0., -conpar[0]-0.1-zd1, irotpipe2, "ONLY");
  //printf("	QC10 CONE from z = %1.2f to z = %1.2f (transition X2ZDC)\n",-zd1,-2*conpar[0]-0.2-zd1);

  zd1 += 2.*conpar[0]+0.2;
  
  // 2 tubes (ID = 63 mm OD=70 mm)      
  tubpar[0] = 6.3/2.;
  tubpar[1] = 7.0/2.;
  tubpar[2] = 639.8/2.;
  TVirtualMC::GetMC()->Gsvolu("QT13", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QT13", 1, "ZDCC", -16.5/2., 0., -tubpar[2]-zd1, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("QT13", 2, "ZDCC",  16.5/2., 0., -tubpar[2]-zd1, 0, "ONLY");
  //printf("	QT13 TUBE from z = %1.2f to z = %1.2f (separate beam pipes)\n",-zd1,-2*tubpar[2]-zd1);  

  zd1 += 2.*tubpar[2];
  //printf("	END OF C SIDE BEAM PIPE DEFINITION @ z = %f m from IP2\n\n",-zd1/100.);

	   
  // -- Luminometer (Cu box) in front of ZN - side C
  if(fLumiLength>0.){
    boxpar[0] = 8.0/2.;
    boxpar[1] = 8.0/2.;
    boxpar[2] = fLumiLength/2.;
    TVirtualMC::GetMC()->Gsvolu("QLUC", "BOX ", idtmed[9], boxpar, 3);
    TVirtualMC::GetMC()->Gspos("QLUC", 1, "ZDCC", 0., 0.,  fPosZNC[2]+66.+boxpar[2], 0, "ONLY");
    printf("	C SIDE LUMINOMETER %1.2f < z < %1.2f\n",  fPosZNC[2]+66., fPosZNC[2]+66.+2*boxpar[2]);
  }
//}  	         
  // --  END OF BEAM PIPE VOLUME DEFINITION FOR SIDE C (RB26 SIDE) 
  // ----------------------------------------------------------------

  ////////////////////////////////////////////////////////////////
  //								//
  //                SIDE A - RB24 				//
  //								//
  ///////////////////////////////////////////////////////////////

  // Rotation Matrices definition
  Int_t irotpipe3, irotpipe4, irotpipe5;
  //-- rotation matrices for the tilted cone after the TDI to recenter vacuum chamber      
  TVirtualMC::GetMC()->Matrix(irotpipe3,90.-1.8934,0.,90.,90.,1.8934,180.);    
  //-- rotation matrices for the tilted tube before and after the TDI 
  TVirtualMC::GetMC()->Matrix(irotpipe4,90.-3.8,0.,90.,90.,3.8,180.);       
  //-- rotation matrix for the tilted cone after the TDI
  TVirtualMC::GetMC()->Matrix(irotpipe5,90.+9.8,0.,90.,90.,9.8,0.);     

  // -- Mother of the ZDCs (Vacuum PCON)		
  zd2 = 1910.22;// zd2 initial value
  
  conpar[0] = 0.;
  conpar[1] = 360.;
  conpar[2] = 2.;
  conpar[3] = zd2;
  conpar[4] = 0.;
  conpar[5] = 55.;
  conpar[6] = 13500.;
  conpar[7] = 0.;
  conpar[8] = 55.;
  TVirtualMC::GetMC()->Gsvolu("ZDCA", "PCON", idtmed[10], conpar, 9);
  TVirtualMC::GetMC()->Gspos("ZDCA", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  
  // To avoid overlaps 1 micron are left between certain volumes!
  Double_t dxNoOverlap = 0.0;
  //zd2 += dxNoOverlap;  
  
  // BEAM PIPE from 19.10 m to inner triplet beginning (22.965 m)  
  tubpar[0] = 6.0/2.;
  tubpar[1] = 6.4/2.;
  tubpar[2] = 386.28/2. - dxNoOverlap; 
  TVirtualMC::GetMC()->Gsvolu("QA01", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA01", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	QA01 TUBE centred in %f from z = %1.2f to z = %1.2f (IT begin)\n",tubpar[2]+zd2,zd2,2*tubpar[2]+zd2);
  
  zd2 += 2.*tubpar[2];  

  // -- FIRST SECTION OF THE BEAM PIPE (from beginning of inner triplet to
  //    beginning of D1)  
  tubpar[0] = 6.3/2.;
  tubpar[1] = 6.7/2.;
  tubpar[2] = 3541.8/2. - dxNoOverlap;
  TVirtualMC::GetMC()->Gsvolu("QA02", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA02", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	QA02 TUBE from z = %1.2f to z= %1.2f (D1 begin)\n",zd2,2*tubpar[2]+zd2);
  
  zd2 += 2.*tubpar[2]; 
  
    
  // -- SECOND SECTION OF THE BEAM PIPE (from the beginning of D1 to the beginning of D2)
  //
  //  FROM (MAGNETIC) BEGINNING OF D1 TO THE (MAGNETIC) END OF D1 + 126.5 cm
  //  CYLINDRICAL PIPE of diameter increasing from 6.75 cm up to 8.0 cm
  //  from magnetic end :
  //  1) 80.1 cm still with ID = 6.75 radial beam screen
  //  2) 2.5 cm conical section from ID = 6.75 to ID = 8.0 cm
  //  3) 43.9 cm straight section (tube) with ID = 8.0 cm

  tubpar[0] = 6.75/2.;
  tubpar[1] = 7.15/2.;
  tubpar[2] = (945.0+80.1)/2.;
  TVirtualMC::GetMC()->Gsvolu("QA03", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA03", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	QA03 TUBE from z = %1.2f to z = %1.2f (D1 end)\n",zd2,2*tubpar[2]+zd2);
  
  zd2 += 2.*tubpar[2];

  // Transition Cone from ID=67.5 mm  to ID=80 mm
  conpar[0] = 2.5/2.;
  conpar[1] = 6.75/2.;
  conpar[2] = 7.15/2.;
  conpar[3] = 8.0/2.;
  conpar[4] = 8.4/2.;
  TVirtualMC::GetMC()->Gsvolu("QA04", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA04", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  //printf("	QA04 CONE from z = %1.2f to z = %1.2f (transition cone)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];
  
  tubpar[0] = 8.0/2.;
  tubpar[1] = 8.4/2.;
  tubpar[2] = (43.9+20.+28.5+28.5)/2.;
  TVirtualMC::GetMC()->Gsvolu("QA05", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA05", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	QA05 TUBE from z = %1.2f to z = %1.2f\n",zd2,2*tubpar[2]+zd2);
  
  zd2 += 2.*tubpar[2];

  // Second section of VAEHI (transition cone from ID=80mm to ID=98mm)
  conpar[0] = 4.0/2.;
  conpar[1] = 8.0/2.;
  conpar[2] = 8.4/2.;
  conpar[3] = 9.8/2.;
  conpar[4] = 10.2/2.;
  TVirtualMC::GetMC()->Gsvolu("QAV1", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QAV1", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  //printf("	QAV1 CONE from z = %1.2f to z = %1.2f (VAEHI-I)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];
  
  //Third section of VAEHI (transition cone from ID=98mm to ID=90mm)
  conpar[0] = 1.0/2.;
  conpar[1] = 9.8/2.;
  conpar[2] = 10.2/2.;
  conpar[3] = 9.0/2.;
  conpar[4] = 9.4/2.;
  TVirtualMC::GetMC()->Gsvolu("QAV2", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QAV2", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  //printf("	QAV2 CONE from z = %1.2f to z = %1.2f (VAEHI-II)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];
 
  // Fourth section of VAEHI (tube ID=90mm)    
  tubpar[0] = 9.0/2.;
  tubpar[1] = 9.4/2.;
  tubpar[2] = 31.0/2.;
  TVirtualMC::GetMC()->Gsvolu("QAV3", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QAV3", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	QAV3 TUBE from z = %1.2f to z = %1.2f (VAEHI-III)\n",zd2,2*tubpar[2]+zd2);
  
  zd2 += 2.*tubpar[2]; 

  //---------------------------- TCDD beginning ----------------------------------    
  // space for the insertion of the collimator TCDD (2 m)
  // TCDD ZONE - 1st volume
  conpar[0] = 1.3/2.;
  conpar[1] = 9.0/2.;
  conpar[2] = 13.0/2.;
  conpar[3] = 9.6/2.;
  conpar[4] = 13.0/2.;
  TVirtualMC::GetMC()->Gsvolu("Q01T", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("Q01T", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  //printf("	Q01T CONE from z = %1.2f to z = %1.2f (TCDD-I)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];  

  // TCDD ZONE - 2nd volume    
  tubpar[0] = 9.6/2.;
  tubpar[1] = 10.0/2.;
  tubpar[2] = 1.0/2.;
  TVirtualMC::GetMC()->Gsvolu("Q02T", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("Q02T", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	Q02T TUBE from z = %1.2f to z= %1.2f (TCDD-II)\n",zd2,2*tubpar[2]+zd2);
  
  zd2 += 2.*tubpar[2]; 

  // TCDD ZONE - third volume
  conpar[0] = 9.04/2.;
  conpar[1] = 9.6/2.;
  conpar[2] = 10.0/2.;
  conpar[3] = 13.8/2.;
  conpar[4] = 14.2/2.;
  TVirtualMC::GetMC()->Gsvolu("Q03T", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("Q03T", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  //printf("	Q03T CONE from z = %1.2f to z= %1.2f (TCDD-III)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];  

  // TCDD ZONE - 4th volume    
  tubpar[0] = 13.8/2.;
  tubpar[1] = 14.2/2.;
  tubpar[2] = 38.6/2.;
  TVirtualMC::GetMC()->Gsvolu("Q04T", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("Q04T", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	Q04T TUBE from z = %1.2f to z= %1.2f (TCDD-IV)\n",zd2,2*tubpar[2]+zd2);
  
  zd2 += 2.*tubpar[2]; 

  // TCDD ZONE - 5th volume    
  tubpar[0] = 21.0/2.;
  tubpar[1] = 21.4/2.;
  tubpar[2] = 100.12/2.;
  TVirtualMC::GetMC()->Gsvolu("Q05T", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("Q05T", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	Q05T TUBE from z = %1.2f to z= %1.2f (TCDD-V)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2]; 
 
  // TCDD ZONE - 6th volume    
  tubpar[0] = 13.8/2.;
  tubpar[1] = 14.2/2.;
  tubpar[2] = 38.6/2.;
  TVirtualMC::GetMC()->Gsvolu("Q06T", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("Q06T", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	Q06T TUBE from z = %1.2f to z= %1.2f (TCDD-VI)\n",zd2,2*tubpar[2]+zd2);
  
  zd2 += 2.*tubpar[2];

  // TCDD ZONE - 7th volume
  conpar[0] = 11.34/2.;
  conpar[1] = 13.8/2.;
  conpar[2] = 14.2/2.;
  conpar[3] = 18.0/2.;
  conpar[4] = 18.4/2.;
  TVirtualMC::GetMC()->Gsvolu("Q07T", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("Q07T", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  //printf("	Q07T CONE from z = %1.2f to z= %1.2f (TCDD-VII)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];

  // Upper section : one single phi segment of a tube 
  //  5 parameters for tubs: inner radius = 0.,
  //	outer radius = 7. cm, half length = 50 cm
  //	phi1 = 0., phi2 = 180. 
  tubspar[0] = 0.0/2.;
  tubspar[1] = 14.0/2.;
  tubspar[2] = 100.0/2.;
  tubspar[3] = 0.;
  tubspar[4] = 180.;  
  TVirtualMC::GetMC()->Gsvolu("Q08T", "TUBS", idtmed[7], tubspar, 5);
  
  // rectangular beam pipe inside TCDD upper section (Vacuum)  
  boxpar[0] = 7.0/2.;
  boxpar[1] = 2.2/2.;
  boxpar[2] = 100./2.;
  TVirtualMC::GetMC()->Gsvolu("Q09T", "BOX ", idtmed[10], boxpar, 3);
  // positioning vacuum box in the upper section of TCDD
  TVirtualMC::GetMC()->Gspos("Q09T", 1, "Q08T", 0., 1.1,  0., 0, "ONLY");
  
  // lower section : one single phi segment of a tube       
  tubspar[0] = 0.0/2.;
  tubspar[1] = 14.0/2.;
  tubspar[2] = 100.0/2.;
  tubspar[3] = 180.;
  tubspar[4] = 360.;  
  TVirtualMC::GetMC()->Gsvolu("Q10T", "TUBS", idtmed[7], tubspar, 5);
  // rectangular beam pipe inside TCDD lower section (Vacuum)  
  boxpar[0] = 7.0/2.;
  boxpar[1] = 2.2/2.;
  boxpar[2] = 100./2.;
  TVirtualMC::GetMC()->Gsvolu("Q11T", "BOX ", idtmed[10], boxpar, 3);
  // positioning vacuum box in the lower section of TCDD
  TVirtualMC::GetMC()->Gspos("Q11T", 1, "Q10T", 0., -1.1,  0., 0, "ONLY");  
  
  // positioning  TCDD elements in ZDCA, (inside TCDD volume)
  TVirtualMC::GetMC()->Gspos("Q08T", 1, "ZDCA", 0., fTCDDAperturePos, -100.+zd2, 0, "ONLY");  
  TVirtualMC::GetMC()->Gspos("Q10T", 1, "ZDCA", 0., -fTCDDApertureNeg, -100.+zd2, 0, "ONLY");  
  printf("  AliZDCv4 -> TCDD apertures +%1.2f/-%1.2f cm\n", 
    	fTCDDAperturePos, fTCDDApertureNeg);
    
  // RF screen 
  boxpar[0] = 0.2/2.;
  boxpar[1] = 4.0/2.;
  boxpar[2] = 100./2.;
  TVirtualMC::GetMC()->Gsvolu("Q12T", "BOX ", idtmed[7], boxpar, 3);  
  // positioning RF screen at both sides of TCDD
  TVirtualMC::GetMC()->Gspos("Q12T", 1, "ZDCA", tubspar[1]+boxpar[0], 0., -100.+zd2, 0, "ONLY");  
  TVirtualMC::GetMC()->Gspos("Q12T", 2, "ZDCA", -tubspar[1]-boxpar[0], 0., -100.+zd2, 0, "ONLY");      
  //---------------------------- TCDD end ---------------------------------------    

  // The following elliptical tube 180 mm x 70 mm
  // (obtained positioning the void QA06 in QA07)
  // represents VAMTF + first part of VCTCP (93 mm)
  // updated according to 2012 new ZDC installation

  tubpar[0] = 18.4/2.;
  tubpar[1] = 7.4/2.;
  tubpar[2] = (78+9.3)/2.;
//  TVirtualMC::GetMC()->Gsvolu("QA06", "ELTU", idtmed[7], tubpar, 3);  
  // temporary replace with a scaled tube (AG)
  TGeoTube *tubeQA06 = new TGeoTube(0.,tubpar[0],tubpar[2]);
  TGeoScale *scaleQA06 = new TGeoScale(1., tubpar[1]/tubpar[0], 1.);
  TGeoScaledShape *sshapeQA06 = new TGeoScaledShape(tubeQA06, scaleQA06);
  new TGeoVolume("QA06", sshapeQA06, gGeoManager->GetMedium(idtmed[7]));
  //printf("	QA06 TUBE from z = %1.2f to z = %1.2f (VAMTF+VCTCP-I)\n",zd2,2*tubpar[2]+zd2);

  tubpar[0] = 18.0/2.;
  tubpar[1] = 7.0/2.;
  tubpar[2] = (78+9.3)/2.;
//  TVirtualMC::GetMC()->Gsvolu("QA07", "ELTU", idtmed[10], tubpar, 3);  
  // temporary replace with a scaled tube (AG)
  TGeoTube *tubeQA07 = new TGeoTube(0.,tubpar[0],tubpar[2]);
  TGeoScale *scaleQA07 = new TGeoScale(1., tubpar[1]/tubpar[0], 1.);
  TGeoScaledShape *sshapeQA07 = new TGeoScaledShape(tubeQA07, scaleQA07);
  new TGeoVolume("QA07", sshapeQA07, gGeoManager->GetMedium(idtmed[10]));
  //printf("	QA07 TUBE from z = %1.2f to z= %1.2f\n",zd2,2*tubpar[2]+zd2);
  TVirtualMC::GetMC()->Gspos("QA06", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("QA07", 1, "QA06", 0., 0., 0., 0, "ONLY");  
    
  zd2 += 2.*tubpar[2];
      
  // VCTCP second part: transition cone from ID=180 to ID=212.7 
  conpar[0] = 31.5/2.;
  conpar[1] = 18.0/2.;
  conpar[2] = 18.6/2.;
  conpar[3] = 21.27/2.;
  conpar[4] = 21.87/2.;
  TVirtualMC::GetMC()->Gsvolu("QA08", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA08", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QA08 CONE from z = %f to z = %f (VCTCP-II)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];
  
  // Tube ID 212.7 mm
  // Represents VCTCP third part (92 mm) + VCDWB (765 mm) + VMBGA (400 mm) +
  //            VCDWE (300 mm) + VMBGA (400 mm)
  // + TCTVB space + VAMTF space (new installation Jan 2012)
  tubpar[0] = 21.27/2.;
  tubpar[1] = 21.87/2.;
  tubpar[2] = (195.7+148.+78.)/2.;
  TVirtualMC::GetMC()->Gsvolu("QA09", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA09", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  //printf("	QA09 TUBE from z = %1.2f to z= %1.2f (VCTCP-III+VCDWB+VMBGA+VCDWE+VMBGA)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];

  // skewed transition piece (ID=212.7 mm to 332 mm) (before TDI)   
  conpar[0] = (50.0-0.73-1.13)/2.;
  conpar[1] = 21.27/2.;
  conpar[2] = 21.87/2.;
  conpar[3] = 33.2/2.;
  conpar[4] = 33.8/2.;
  TVirtualMC::GetMC()->Gsvolu("QA10", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA10", 1, "ZDCA", -1.66, 0., conpar[0]+0.73+zd2, irotpipe4, "ONLY");
  // Ch.debug  
  //printf("	QA10 skewed CONE from z = %1.2f to z= %1.2f\n",zd2,2*conpar[0]+0.73+1.13+zd2);

  zd2 += 2.*conpar[0]+0.73+1.13;
      
  // Vacuum chamber containing TDI  
  tubpar[0] = 0.;
  tubpar[1] = 54.6/2.;
  tubpar[2] = 540.0/2.;
  TVirtualMC::GetMC()->Gsvolu("Q13TM", "TUBE", idtmed[10], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("Q13TM", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  tubpar[0] = 54.0/2.;
  tubpar[1] = 54.6/2.;
  tubpar[2] = 540.0/2.;
  TVirtualMC::GetMC()->Gsvolu("Q13T", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("Q13T", 1, "Q13TM", 0., 0., 0., 0, "ONLY");
  // Ch.debug
  //printf("	Q13T TUBE from z = %1.2f to z= %1.2f (TDI vacuum chamber)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];
  
  //---------------- INSERT TDI INSIDE Q13T -----------------------------------    
  boxpar[0] = 11.0/2.;
  boxpar[1] = 9.0/2.;
  boxpar[2] = 418.5/2.;
  TVirtualMC::GetMC()->Gsvolu("QTD1", "BOX ", idtmed[7], boxpar, 3);
  TVirtualMC::GetMC()->Gspos("QTD1", 1, "Q13TM", -3.8, boxpar[1]+fTDIAperturePos,  0., 0, "ONLY");
  boxpar[0] = 11.0/2.;
  boxpar[1] = 9.0/2.;
  boxpar[2] = 418.5/2.;
  TVirtualMC::GetMC()->Gsvolu("QTD2", "BOX ", idtmed[7], boxpar, 3);
  TVirtualMC::GetMC()->Gspos("QTD2", 1, "Q13TM", -3.8, -boxpar[1]-fTDIApertureNeg,  0., 0, "ONLY");  
  boxpar[0] = 5.1/2.;
  boxpar[1] = 0.2/2.;
  boxpar[2] = 418.5/2.;
  TVirtualMC::GetMC()->Gsvolu("QTD3", "BOX ", idtmed[7], boxpar, 3);
  TVirtualMC::GetMC()->Gspos("QTD3", 1, "Q13TM", -3.8+5.5+boxpar[0], fTDIAperturePos,  0., 0, "ONLY");  
  TVirtualMC::GetMC()->Gspos("QTD3", 2, "Q13TM", -3.8+5.5+boxpar[0], -fTDIApertureNeg,  0., 0, "ONLY"); 
  TVirtualMC::GetMC()->Gspos("QTD3", 3, "Q13TM", -3.8-5.5-boxpar[0], fTDIAperturePos,  0., 0, "ONLY");  
  TVirtualMC::GetMC()->Gspos("QTD3", 4, "Q13TM", -3.8-5.5-boxpar[0], -fTDIApertureNeg,  0., 0, "ONLY");  
  printf("  AliZDCv4 -> TDI apertures +%1.2f/-%1.2f cm\n", fTDIAperturePos, fTDIApertureNeg);
  //
  tubspar[0] = 12.0/2.;
  tubspar[1] = 12.4/2.;
  tubspar[2] = 418.5/2.;
  tubspar[3] = 90.;
  tubspar[4] = 270.;  
  TVirtualMC::GetMC()->Gsvolu("QTD4", "TUBS", idtmed[6], tubspar, 5);
  TVirtualMC::GetMC()->Gspos("QTD4", 1, "Q13TM", -3.8-10.6, 0.,  0., 0, "ONLY");
  tubspar[0] = 12.0/2.;
  tubspar[1] = 12.4/2.;
  tubspar[2] = 418.5/2.;
  tubspar[3] = -90.;
  tubspar[4] = 90.;  
  TVirtualMC::GetMC()->Gsvolu("QTD5", "TUBS", idtmed[6], tubspar, 5);
  TVirtualMC::GetMC()->Gspos("QTD5", 1, "Q13TM", -3.8+10.6, 0.,  0., 0, "ONLY"); 
  //---------------- END DEFINING TDI INSIDE Q13T -------------------------------
  
  // VCTCG skewed transition piece (ID=332 mm to 212.7 mm) (after TDI)
  conpar[0] = (50.0-2.92-1.89)/2.;
  conpar[1] = 33.2/2.;
  conpar[2] = 33.8/2.;
  conpar[3] = 21.27/2.;
  conpar[4] = 21.87/2.;
  TVirtualMC::GetMC()->Gsvolu("QA11", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA11", 1, "ZDCA", 4.32-3.8, 0., conpar[0]+2.92+zd2, irotpipe5, "ONLY");
  // Ch.debug  
  //printf("	QA11 skewed CONE from z = %f to z =%f (VCTCG)\n",zd2,2*conpar[0]+2.92+1.89+zd2);

  zd2 += 2.*conpar[0]+2.92+1.89;
  
  // The following tube ID 212.7 mm  
  // represents VMBGA (400 mm) + VCDWE (300 mm) + VMBGA (400 mm) +
  //            BTVTS (600 mm) + VMLGB (400 mm)  
  tubpar[0] = 21.27/2.;
  tubpar[1] = 21.87/2.;
  tubpar[2] = 210.0/2.;
  TVirtualMC::GetMC()->Gsvolu("QA12", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA12", 1, "ZDCA", 4., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	QA12 TUBE from z = %1.2f to z= %1.2f (VMBGA+VCDWE+VMBGA+BTVTS+VMLGB)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];  
  
  // First part of VCTCC
  // skewed transition cone from ID=212.7 mm to ID=797 mm
  conpar[0] = (121.0-0.37-1.35)/2.;
  conpar[1] = 21.27/2.;
  conpar[2] = 21.87/2.;
  conpar[3] = 79.7/2.;
  conpar[4] = 81.3/2.;
  TVirtualMC::GetMC()->Gsvolu("QA13", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA13", 1, "ZDCA", 4.-2., 0., conpar[0]+0.37+zd2, irotpipe3, "ONLY");
  // Ch.debug  
  //printf("	QA13 CONE from z = %1.2f to z = %1.2f (VCTCC-I)\n",zd2,2*conpar[0]+0.37+1.35+zd2);

  zd2 += 2.*conpar[0]+0.37+1.35;
  
  // The following tube ID 797 mm  
  // represents the second part of VCTCC (4272 mm) + 
  //            4 x VCDGA (4 x 4272 mm) + 
  //            the first part of VCTCR (850 mm)
  // updated according to 2012 ZDC installation
  tubpar[0] = 79.7/2.;
  tubpar[1] = 81.3/2.;
  tubpar[2] = (2221.-136.)/2.;
  TVirtualMC::GetMC()->Gsvolu("QA14", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA14", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QA14 TUBE from z = %1.2f to z = %1.2f (VCTCC-II)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];
        
  // Second part of VCTCR
  // Transition from ID=797 mm to ID=196 mm:
  // in order to simulate the thin window opened in the transition cone
  // we divide the transition cone in three cones:
  // (1) 8 mm thick (2) 3 mm thick (3) the third 8 mm thick
  
  // (1) 8 mm thick
  conpar[0] = 9.09/2.; // 15 degree
  conpar[1] = 79.7/2.;
  conpar[2] = 81.3/2.; // thickness 8 mm  
  conpar[3] = 74.82868/2.;
  conpar[4] = 76.42868/2.; // thickness 8 mm 
  TVirtualMC::GetMC()->Gsvolu("QA15", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA15", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  //printf("	QA15 CONE from z = %1.2f to z= %1.2f (VCTCR-I)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];  

  // (2) 3 mm thick
  conpar[0] = 96.2/2.; // 15 degree
  conpar[1] = 74.82868/2.;
  conpar[2] = 75.42868/2.; // thickness 3 mm  
  conpar[3] = 23.19588/2.;
  conpar[4] = 23.79588/2.; // thickness 3 mm 
  TVirtualMC::GetMC()->Gsvolu("QA16", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA16", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");  
  //printf("	QA16 CONE from z = %1.2f to z= %1.2f\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];
  
  // (3) 8 mm thick
  conpar[0] = 6.71/2.; // 15 degree
  conpar[1] = 23.19588/2.;
  conpar[2] = 24.79588/2.;// thickness 8 mm 
  conpar[3] = 19.6/2.;
  conpar[4] = 21.2/2.;// thickness 8 mm 
  TVirtualMC::GetMC()->Gsvolu("QA17", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA17", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  //printf("	QA17 CONE from z = %1.2f to z= %1.2f (VCTCR-II)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];
 
  // Third part of VCTCR: tube (ID=196 mm)  
  tubpar[0] = 19.6/2.;
  tubpar[1] = 21.2/2.;
  tubpar[2] = 9.55/2.;
  TVirtualMC::GetMC()->Gsvolu("QA18", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA18", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QA18 TUBE from z = %1.2f to z= %1.2f (VCTCR-III)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];  
  
  // Flange (ID=196 mm) (last part of VCTCR and first part of VMZAR)
  tubpar[0] = 19.6/2.;
  tubpar[1] = 25.3/2.;
  tubpar[2] = 4.9/2.;
  TVirtualMC::GetMC()->Gsvolu("QF01", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QF01", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QF01  TUBE from z = %1.2f to z= %1.2f (VMZAR-I)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];
  
  // VMZAR (5 volumes)  
  tubpar[0] = 20.2/2.;
  tubpar[1] = 20.6/2.;
  tubpar[2] = 2.15/2.;
  TVirtualMC::GetMC()->Gsvolu("QA19", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA19", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QA19  TUBE from z = %1.2f to z = %1.2f (VMZAR-II)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];
  
  conpar[0] = 6.9/2.;
  conpar[1] = 20.2/2.;
  conpar[2] = 20.6/2.;
  conpar[3] = 23.9/2.;
  conpar[4] = 24.3/2.;
  TVirtualMC::GetMC()->Gsvolu("QA20", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA20", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QA20 CONE from z = %1.2f to z = %1.2f (VMZAR-III)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];

  tubpar[0] = 23.9/2.;
  tubpar[1] = 25.5/2.;
  tubpar[2] = 17.0/2.;
  TVirtualMC::GetMC()->Gsvolu("QA21", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA21", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QA21  TUBE from z = %1.2f to z = %1.2f  (VMZAR-IV)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];
  
  conpar[0] = 6.9/2.;
  conpar[1] = 23.9/2.;
  conpar[2] = 24.3/2.;
  conpar[3] = 20.2/2.;
  conpar[4] = 20.6/2.;
  TVirtualMC::GetMC()->Gsvolu("QA22", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA22", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QA22 CONE from z = %1.2f to z = %1.2f (VMZAR-V)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0];
  
  tubpar[0] = 20.2/2.;
  tubpar[1] = 20.6/2.;
  tubpar[2] = 2.15/2.;
  TVirtualMC::GetMC()->Gsvolu("QA23", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA23", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QA23  TUBE from z = %1.2f to z= %1.2f (VMZAR-VI)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];
  
  // Flange (ID=196 mm)(last part of VMZAR and first part of VCTYD)
  tubpar[0] = 19.6/2.;
  tubpar[1] = 25.3/2.;
  tubpar[2] = 4.9/2.;
  TVirtualMC::GetMC()->Gsvolu("QF02", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QF02", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QF02 TUBE from z = %1.2f to z= %1.2f (VMZAR-VII)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];
  
  // simulation of the trousers (VCTYB)     
  tubpar[0] = 19.6/2.;
  tubpar[1] = 20.0/2.;
  tubpar[2] = 3.9/2.;
  TVirtualMC::GetMC()->Gsvolu("QA24", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA24", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug
  //printf("	QA24  TUBE from z = %1.2f to z= %1.2f (VCTYB)\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];

  // transition cone from ID=196. to ID=216.6
  conpar[0] = 32.55/2.;
  conpar[1] = 19.6/2.;
  conpar[2] = 20.0/2.;
  conpar[3] = 21.66/2.;
  conpar[4] = 22.06/2.;
  TVirtualMC::GetMC()->Gsvolu("QA25", "CONE", idtmed[7], conpar, 5);
  TVirtualMC::GetMC()->Gspos("QA25", 1, "ZDCA", 0., 0., conpar[0]+zd2, 0, "ONLY");
  // Ch.debug  
  //printf("	QA25 CONE from z = %1.2f to z= %1.2f (transition cone)\n",zd2,2*conpar[0]+zd2);

  zd2 += 2.*conpar[0]; 
  
  // tube  
  tubpar[0] = 21.66/2.;
  tubpar[1] = 22.06/2.;
  tubpar[2] = 28.6/2.;
  TVirtualMC::GetMC()->Gsvolu("QA26", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA26", 1, "ZDCA", 0., 0., tubpar[2]+zd2, 0, "ONLY");
  // Ch.debug 
  //printf("	QA26  TUBE from z = %1.2f to z= %1.2f\n",zd2,2*tubpar[2]+zd2);

  zd2 += 2.*tubpar[2];
  // Ch.debug
  //printf("	Begin of recombination chamber z = %1.2f\n",zd2);

  // --------------------------------------------------------
  // RECOMBINATION CHAMBER IMPLEMENTED USING TGeo CLASSES!!!!
  // author: Chiara (June 2008)
  // --------------------------------------------------------
  // TRANSFORMATION MATRICES
  // Combi transformation: 
  dx = -3.970000;
  dy = 0.000000;
  dz = 0.0;
  // Rotation: 
  thx = 84.989100;   phx = 0.000000;
  thy = 90.000000;   phy = 90.000000;
  thz = 5.010900;    phz = 180.000000;
  TGeoRotation *rotMatrix1 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  // Combi transformation: 
  dx = -3.970000;
  dy = 0.000000;
  dz = 0.0;
  TGeoCombiTrans *rotMatrix2 = new TGeoCombiTrans("ZDC_c1", dx,dy,dz,rotMatrix1);
  rotMatrix2->RegisterYourself();
  // Combi transformation: 
  dx = 3.970000;
  dy = 0.000000;
  dz = 0.0;
  // Rotation: 
  thx = 95.010900;   phx = 0.000000;
  thy = 90.000000;   phy = 90.000000;
  thz = 5.010900;    phz = 0.000000;
  TGeoRotation *rotMatrix3 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans *rotMatrix4 = new TGeoCombiTrans("ZDC_c2", dx,dy,dz,rotMatrix3);
  rotMatrix4->RegisterYourself();
  
  
  // VOLUMES DEFINITION
  // Volume: ZDCA
  TGeoVolume *pZDCA = gGeoManager->GetVolume("ZDCA");
  
  conpar[0] = (90.1-0.95-0.26)/2.;
  conpar[1] = 0.0/2.;
  conpar[2] = 21.6/2.;
  conpar[3] = 0.0/2.;
  conpar[4] = 5.8/2.;
  new TGeoCone("QALext", conpar[0],conpar[1],conpar[2],conpar[3],conpar[4]);
  
  conpar[0] = (90.1-0.95-0.26)/2.;
  conpar[1] = 0.0/2.;
  conpar[2] = 21.2/2.;
  conpar[3] = 0.0/2.;
  conpar[4] = 5.4/2.;
  new TGeoCone("QALint", conpar[0],conpar[1],conpar[2],conpar[3],conpar[4]);

  // Outer trousers
  TGeoCompositeShape *pOutTrousers = new TGeoCompositeShape("outTrousers", "QALext:ZDC_c1+QALext:ZDC_c2");
  
  // Volume: QALext
  TGeoVolume *pQALext = new TGeoVolume("QALext",pOutTrousers, medZDCFe);
  pQALext->SetLineColor(kBlue);
  pQALext->SetVisLeaves(kTRUE);
  //
  TGeoTranslation *tr1 = new TGeoTranslation(0., 0., (Double_t) conpar[0]+0.95+zd2);
  pZDCA->AddNode(pQALext, 1, tr1);
  // Inner trousers
  TGeoCompositeShape *pIntTrousers = new TGeoCompositeShape("intTrousers", "QALint:ZDC_c1+QALint:ZDC_c2");
  // Volume: QALint
  //TGeoMedium *medZDCvoid = gGeoManager->GetMedium("ZDC_ZVOID");
  TGeoVolume *pQALint = new TGeoVolume("QALint",pIntTrousers, medZDCvoid);
  pQALint->SetLineColor(kAzure);
  pQALint->SetVisLeaves(kTRUE);
  pQALext->AddNode(pQALint, 1);
    
  zd2 += 90.1;
  // Ch.debug
  //printf("	End of recombination chamber z = %1.2f\n",zd2);
  
  
  //  second section : 2 tubes (ID = 54. OD = 58.)  
  tubpar[0] = 5.4/2.;
  tubpar[1] = 5.8/2.;
  tubpar[2] = 40.0/2.;
  TVirtualMC::GetMC()->Gsvolu("QA27", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA27", 1, "ZDCA", -15.8/2., 0., tubpar[2]+zd2, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("QA27", 2, "ZDCA",  15.8/2., 0., tubpar[2]+zd2, 0, "ONLY");  
  // Ch.debug
  //printf("	QA27 TUBE from z = %1.2f to z= %1.2f (separate pipes)\n",zd2,2*tubpar[2]+zd2);
  
  zd2 += 2.*tubpar[2];
 
  // transition x2zdc to recombination chamber : skewed cone  
  conpar[0] = (10.-1.)/2.;
  conpar[1] = 5.4/2.;
  conpar[2] = 5.8/2.;
  conpar[3] = 6.3/2.;
  conpar[4] = 7.0/2.;
  TVirtualMC::GetMC()->Gsvolu("QA28", "CONE", idtmed[7], conpar, 5); 
  TVirtualMC::GetMC()->Gspos("QA28", 1, "ZDCA", -7.9-0.175, 0., conpar[0]+0.5+zd2, irotpipe1, "ONLY");
  TVirtualMC::GetMC()->Gspos("QA28", 2, "ZDCA", 7.9+0.175, 0., conpar[0]+0.5+zd2, irotpipe2, "ONLY");
  //printf("	QA28 CONE from z = %1.2f to z= %1.2f (transition X2ZDC)\n",zd2,2*conpar[0]+0.2+zd2);

  zd2 += 2.*conpar[0]+1.;
  
  // 2 tubes (ID = 63 mm OD=70 mm)      
  tubpar[0] = 6.3/2.;
  tubpar[1] = 7.0/2.;
  tubpar[2] = (342.5+498.3)/2.;
  TVirtualMC::GetMC()->Gsvolu("QA29", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("QA29", 1, "ZDCA", -16.5/2., 0., tubpar[2]+zd2, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("QA29", 2, "ZDCA",  16.5/2., 0., tubpar[2]+zd2, 0, "ONLY");
  //printf("	QA29 TUBE from z = %1.2f to z= %1.2f (separate pipes)\n",zd2,2*tubpar[2]+zd2);  

  zd2 += 2.*tubpar[2];
	   
  // -- Luminometer (Cu box) in front of ZN - side A
  if(fLumiLength>0.){
    boxpar[0] = 8.0/2.;
    boxpar[1] = 8.0/2.;
    boxpar[2] = fLumiLength/2.;
    TVirtualMC::GetMC()->Gsvolu("QLUA", "BOX ", idtmed[9], boxpar, 3);
    TVirtualMC::GetMC()->Gspos("QLUA", 1, "ZDCA", 0., 0.,  fPosZNA[2]-66.-boxpar[2], 0, "ONLY");
    printf("	A SIDE LUMINOMETER %1.2f < z < %1.2f\n\n",  fPosZNA[2]-66., fPosZNA[2]-66.-2*boxpar[2]);
  }
  //printf("	END OF A SIDE BEAM PIPE VOLUME DEFINITION AT z = %f m from IP2\n",zd2/100.);
  

  // ----------------------------------------------------------------
  // --  MAGNET DEFINITION  -> LHC OPTICS 6.5  
  // ----------------------------------------------------------------      
  // ***************************************************************  
  //		SIDE C - RB26  (dimuon side) 
  // ***************************************************************   
  // --  COMPENSATOR DIPOLE (MBXW)
  zCorrDip = 1972.5;   
  
  // --  GAP (VACUUM WITH MAGNETIC FIELD)
  tubpar[0] = 0.;
  tubpar[1] = 3.14;
  // New -> Added to accomodate AD (A. Morsch)
  tubpar[2] = 150./2.;
  //tubpar[2] = 153./2.;
  TVirtualMC::GetMC()->Gsvolu("MBXW", "TUBE", idtmed[11], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("MBXW", 1, "ZDCC", 0., 0., -tubpar[2]-zCorrDip, 0, "ONLY");
  // Ch.debug
  //printf("	MBXW volume: %1.2f < z < %1.2f\n\n",  -zCorrDip, -zCorrDip-2*tubpar[2]);
  // --  YOKE 
  tubpar[0] = 4.5;
  tubpar[1] = 55.;
  tubpar[2] = 153./2.;
  TVirtualMC::GetMC()->Gsvolu("YMBX", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("YMBX", 1, "ZDCC", 0., 0., -tubpar[2]-zCorrDip, 0, "ONLY");
  // Ch.debug
  //printf("	MBXW yoke: %1.2f < z < %1.2f\n\n",  -zCorrDip, -zCorrDip-2*tubpar[2]);
  
  
  // -- INNER TRIPLET 
  zInnTrip = 2296.5; 

  // -- DEFINE MQXL AND MQX QUADRUPOLE ELEMENT 
  // --  MQXL 
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  tubpar[0] = 0.;
  tubpar[1] = 3.14;
  tubpar[2] = 637./2.;
  TVirtualMC::GetMC()->Gsvolu("MQXL", "TUBE", idtmed[11], tubpar, 3);
    
  // --  YOKE 
  tubpar[0] = 3.5;
  tubpar[1] = 22.;
  tubpar[2] = 637./2.;
  TVirtualMC::GetMC()->Gsvolu("YMQL", "TUBE", idtmed[7], tubpar, 3);
  
  TVirtualMC::GetMC()->Gspos("MQXL", 1, "ZDCC", 0., 0., -tubpar[2]-zInnTrip, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("YMQL", 1, "ZDCC", 0., 0., -tubpar[2]-zInnTrip, 0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("MQXL", 2, "ZDCC", 0., 0., -tubpar[2]-zInnTrip-2400., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("YMQL", 2, "ZDCC", 0., 0., -tubpar[2]-zInnTrip-2400., 0, "ONLY");
  
  // --  MQX 
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  tubpar[0] = 0.;
  tubpar[1] = 3.14;
  tubpar[2] = 550./2.;
  TVirtualMC::GetMC()->Gsvolu("MQX ", "TUBE", idtmed[11], tubpar, 3);
  
  // --  YOKE 
  tubpar[0] = 3.5;
  tubpar[1] = 22.;
  tubpar[2] = 550./2.;
  TVirtualMC::GetMC()->Gsvolu("YMQ ", "TUBE", idtmed[7], tubpar, 3);
  
  TVirtualMC::GetMC()->Gspos("MQX ", 1, "ZDCC", 0., 0., -tubpar[2]-zInnTrip-908.5,  0, "ONLY");
  TVirtualMC::GetMC()->Gspos("YMQ ", 1, "ZDCC", 0., 0., -tubpar[2]-zInnTrip-908.5,  0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("MQX ", 2, "ZDCC", 0., 0., -tubpar[2]-zInnTrip-1558.5, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("YMQ ", 2, "ZDCC", 0., 0., -tubpar[2]-zInnTrip-1558.5, 0, "ONLY");
  
  // -- SEPARATOR DIPOLE D1 
  zD1 = 5838.3001;
  
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  tubpar[0] = 0.;
  tubpar[1] = 3.46;
  tubpar[2] = 945./2.;
  TVirtualMC::GetMC()->Gsvolu("MD1 ", "TUBE", idtmed[11], tubpar, 3);
  
  // --  Insert horizontal Cu plates inside D1 
  // --   (to simulate the vacuum chamber)
  boxpar[0] = TMath::Sqrt(tubpar[1]*tubpar[1]-(2.98+0.2)*(2.98+0.2)) - 0.05;
  boxpar[1] = 0.2/2.;
  boxpar[2] = 945./2.;
  TVirtualMC::GetMC()->Gsvolu("MD1V", "BOX ", idtmed[6], boxpar, 3);
  TVirtualMC::GetMC()->Gspos("MD1V", 1, "MD1 ", 0., 2.98+boxpar[1], 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("MD1V", 2, "MD1 ", 0., -2.98-boxpar[1], 0., 0, "ONLY");
    
  // --  YOKE 
  tubpar[0] = 3.68;
  tubpar[1] = 110./2.;
  tubpar[2] = 945./2.;
  TVirtualMC::GetMC()->Gsvolu("YD1 ", "TUBE", idtmed[7], tubpar, 3);
  
  TVirtualMC::GetMC()->Gspos("YD1 ", 1, "ZDCC", 0., 0., -tubpar[2]-zD1, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("MD1 ", 1, "ZDCC", 0., 0., -tubpar[2]-zD1, 0, "ONLY");
  // Ch debug
  //printf("	MD1 from z = %1.2f to z= %1.2f cm\n",-zD1, -zD1-2*tubpar[2]); 
  
  // -- DIPOLE D2 
  Float_t zD2 = 12167.8;
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  tubpar[0] = 0.;
  tubpar[1] = 7.5/2.;
  tubpar[2] = 945./2.;
  TVirtualMC::GetMC()->Gsvolu("MD2 ", "TUBE", idtmed[11], tubpar, 3);
  
  // --  YOKE 
  tubpar[0] = 0.;
  tubpar[1] = 55.;
  tubpar[2] = 945./2.;
  TVirtualMC::GetMC()->Gsvolu("YD2 ", "TUBE", idtmed[7], tubpar, 3);
  
  TVirtualMC::GetMC()->Gspos("YD2 ", 1, "ZDCC", 0., 0., -tubpar[2]-zD2, 0, "ONLY");
  // Ch debug
  printf("	YD2 from z = %1.2f to z= %1.2f cm\n",-zD2, -zD2-2*tubpar[2]); 
  
  TVirtualMC::GetMC()->Gspos("MD2 ", 1, "YD2 ", -9.4, 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("MD2 ", 2, "YD2 ",  9.4, 0., 0., 0, "ONLY");
  
  // ***************************************************************  
  //		SIDE A - RB24 
  // ***************************************************************
  
  // COMPENSATOR DIPOLE (MCBWA) (2nd compensator)
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  tubpar[0] = 0.;
  tubpar[1] = 3.;  
  tubpar[2] = 153./2.;
  TVirtualMC::GetMC()->Gsvolu("MCBW", "TUBE", idtmed[11], tubpar, 3);  
  TVirtualMC::GetMC()->Gspos("MCBW", 1, "ZDCA", 0., 0., tubpar[2]+zCorrDip, 0, "ONLY");
  // Ch.debug
  //printf("	MCBWA volume: %1.2f < z < %1.2f\n\n",  zCorrDip, zCorrDip+2*tubpar[2]);
    
   // --  YOKE 
  tubpar[0] = 4.5;
  tubpar[1] = 55.;
  tubpar[2] = 153./2.;
  TVirtualMC::GetMC()->Gsvolu("YMCB", "TUBE", idtmed[7], tubpar, 3);
  TVirtualMC::GetMC()->Gspos("YMCB", 1, "ZDCA", 0., 0., tubpar[2]+zCorrDip, 0, "ONLY");  
  // Ch.debug
  //printf("	MCBWA volume: %1.2f < z < %1.2f\n\n",  zCorrDip, zCorrDip+2*tubpar[2]);
  
   // -- INNER TRIPLET 
  // -- DEFINE MQX1 AND MQX2 QUADRUPOLE ELEMENT 
  // --  MQX1 
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  tubpar[0] = 0.;
  tubpar[1] = 3.14;
  tubpar[2] = 637./2.;
  TVirtualMC::GetMC()->Gsvolu("MQX1", "TUBE", idtmed[11], tubpar, 3);
  TVirtualMC::GetMC()->Gsvolu("MQX4", "TUBE", idtmed[11], tubpar, 3);
    
  // --  YOKE 
  tubpar[0] = 3.5;
  tubpar[1] = 22.;
  tubpar[2] = 637./2.;
  TVirtualMC::GetMC()->Gsvolu("YMQ1", "TUBE", idtmed[7], tubpar, 3);

  // -- Q1
  TVirtualMC::GetMC()->Gspos("MQX1", 1, "ZDCA", 0., 0., tubpar[2]+zInnTrip, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("YMQ1", 1, "ZDCA", 0., 0., tubpar[2]+zInnTrip, 0, "ONLY");

   // -- BEAM SCREEN FOR Q1
   tubpar[0] = 4.78/2.;
   tubpar[1] = 5.18/2.;
   tubpar[2] = 637./2.;
   TVirtualMC::GetMC()->Gsvolu("QBS1", "TUBE", idtmed[6], tubpar, 3);
   TVirtualMC::GetMC()->Gspos("QBS1", 1, "MQX1", 0., 0., 0., 0, "ONLY");
   // INSERT VERTICAL PLATE INSIDE Q1
   boxpar[0] = 0.2/2.0;
   boxpar[1] = TMath::Sqrt(tubpar[0]*tubpar[0]-(1.9+0.2)*(1.9+0.2));
   boxpar[2] = 637./2.;
   TVirtualMC::GetMC()->Gsvolu("QBS2", "BOX ", idtmed[6], boxpar, 3);
   TVirtualMC::GetMC()->Gspos("QBS2", 1, "MQX1", 1.9+boxpar[0], 0., 0., 0, "ONLY");
   TVirtualMC::GetMC()->Gspos("QBS2", 2, "MQX1", -1.9-boxpar[0], 0., 0., 0, "ONLY");

   // -- Q3   
   TVirtualMC::GetMC()->Gspos("MQX4", 1, "ZDCA", 0., 0., tubpar[2]+zInnTrip+2400., 0, "ONLY");
   TVirtualMC::GetMC()->Gspos("YMQ1", 2, "ZDCA", 0., 0., tubpar[2]+zInnTrip+2400., 0, "ONLY");

   // -- BEAM SCREEN FOR Q3
   tubpar[0] = 5.79/2.;
   tubpar[1] = 6.14/2.;
   tubpar[2] = 637./2.;
   TVirtualMC::GetMC()->Gsvolu("QBS3", "TUBE", idtmed[6], tubpar, 3);
   TVirtualMC::GetMC()->Gspos("QBS3", 1, "MQX4", 0., 0., 0., 0, "ONLY");
   // INSERT VERTICAL PLATE INSIDE Q3
   boxpar[0] = 0.2/2.0;
   boxpar[1] = TMath::Sqrt(tubpar[0]*tubpar[0]-(2.405+0.2)*(2.405+0.2));
   boxpar[2] =637./2.;
   TVirtualMC::GetMC()->Gsvolu("QBS4", "BOX ", idtmed[6], boxpar, 3);
   TVirtualMC::GetMC()->Gspos("QBS4", 1, "MQX4", 2.405+boxpar[0], 0., 0., 0, "ONLY");
   TVirtualMC::GetMC()->Gspos("QBS4", 2, "MQX4", -2.405-boxpar[0], 0., 0., 0, "ONLY");
    
  
  
  // --  MQX2
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  tubpar[0] = 0.;
  tubpar[1] = 3.14;
  tubpar[2] = 550./2.;
  TVirtualMC::GetMC()->Gsvolu("MQX2", "TUBE", idtmed[11], tubpar, 3);
  TVirtualMC::GetMC()->Gsvolu("MQX3", "TUBE", idtmed[11], tubpar, 3);
  
  // --  YOKE 
  tubpar[0] = 3.5;
  tubpar[1] = 22.;
  tubpar[2] = 550./2.;
  TVirtualMC::GetMC()->Gsvolu("YMQ2", "TUBE", idtmed[7], tubpar, 3);

   // -- BEAM SCREEN FOR Q2
   tubpar[0] = 5.79/2.;
   tubpar[1] = 6.14/2.;
   tubpar[2] = 550./2.;
   TVirtualMC::GetMC()->Gsvolu("QBS5", "TUBE", idtmed[6], tubpar, 3);
   //    VERTICAL PLATE INSIDE Q2
   boxpar[0] = 0.2/2.0;
   boxpar[1] = TMath::Sqrt(tubpar[0]*tubpar[0]-(2.405+0.2)*(2.405+0.2));
   boxpar[2] =550./2.;
   TVirtualMC::GetMC()->Gsvolu("QBS6", "BOX ", idtmed[6], boxpar, 3);

  // -- Q2A
  TVirtualMC::GetMC()->Gspos("MQX2", 1, "ZDCA", 0., 0., tubpar[2]+zInnTrip+908.5,  0, "ONLY");
  TVirtualMC::GetMC()->Gspos("QBS5", 1, "MQX2", 0., 0., 0., 0, "ONLY");  
  TVirtualMC::GetMC()->Gspos("QBS6", 1, "MQX2", 2.405+boxpar[0], 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("QBS6", 2, "MQX2", -2.405-boxpar[0], 0., 0., 0, "ONLY");  
  TVirtualMC::GetMC()->Gspos("YMQ2", 1, "ZDCA", 0., 0., tubpar[2]+zInnTrip+908.5,  0, "ONLY");

  
  // -- Q2B
  TVirtualMC::GetMC()->Gspos("MQX3", 1, "ZDCA", 0., 0., tubpar[2]+zInnTrip+1558.5, 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("QBS5", 2, "MQX3", 0., 0., 0., 0, "ONLY");  
  TVirtualMC::GetMC()->Gspos("QBS6", 3, "MQX3", 2.405+boxpar[0], 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("QBS6", 4, "MQX3", -2.405-boxpar[0], 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("YMQ2", 2, "ZDCA", 0., 0., tubpar[2]+zInnTrip+1558.5, 0, "ONLY");

  // -- SEPARATOR DIPOLE D1 
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  tubpar[0] = 0.;
  tubpar[1] = 6.75/2.;//3.375
  tubpar[2] = 945./2.;
  TVirtualMC::GetMC()->Gsvolu("MD1L", "TUBE", idtmed[11], tubpar, 3);

  // --  The beam screen tube is provided by the beam pipe in D1 (QA03 volume)
  // --  Insert the beam screen horizontal Cu plates inside D1  
  // --   (to simulate the vacuum chamber)
  boxpar[0] = TMath::Sqrt(tubpar[1]*tubpar[1]-(2.885+0.2)*(2.885+0.2));
  boxpar[1] = 0.2/2.;
  boxpar[2] =945./2.;  
  TVirtualMC::GetMC()->Gsvolu("QBS7", "BOX ", idtmed[6], boxpar, 3);
  TVirtualMC::GetMC()->Gspos("QBS7", 1, "MD1L", 0., 2.885+boxpar[1],0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("QBS7", 2, "MD1L", 0., -2.885-boxpar[1],0., 0, "ONLY");  
    
  // --  YOKE 
  tubpar[0] = 3.68;
  tubpar[1] = 110./2;
  tubpar[2] = 945./2.;
  TVirtualMC::GetMC()->Gsvolu("YD1L", "TUBE", idtmed[7], tubpar, 3);
  
  TVirtualMC::GetMC()->Gspos("YD1L", 1, "ZDCA", 0., 0., tubpar[2]+zD1, 0, "ONLY");  
  TVirtualMC::GetMC()->Gspos("MD1L", 1, "ZDCA", 0., 0., tubpar[2]+zD1, 0, "ONLY");  
  
  // -- DIPOLE D2 
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  tubpar[0] = 0.;
  tubpar[1] = 7.5/2.; // this has to be checked
  tubpar[2] = 945./2.;
  TVirtualMC::GetMC()->Gsvolu("MD2L", "TUBE", idtmed[11], tubpar, 3);
  
  // --  YOKE 
  tubpar[0] = 0.;
  tubpar[1] = 55.;
  tubpar[2] = 945./2.;
  TVirtualMC::GetMC()->Gsvolu("YD2L", "TUBE", idtmed[7], tubpar, 3);
  
  TVirtualMC::GetMC()->Gspos("YD2L", 1, "ZDCA", 0., 0., tubpar[2]+zD2, 0, "ONLY");
  
  TVirtualMC::GetMC()->Gspos("MD2L", 1, "YD2L", -9.4, 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("MD2L", 2, "YD2L",  9.4, 0., 0., 0, "ONLY");
  
  // -- END OF MAGNET DEFINITION     
}
  
//_____________________________________________________________________________
void AliZDCv4::CreateZDC()
{
 //
 // Create the various ZDCs (ZN + ZP)
 //
  
  Float_t dimPb[6], dimVoid[6];
  
  Int_t *idtmed = fIdtmed->GetArray();

  // Parameters for EM calorimeter geometry
  // NB -> parameters used ONLY in CreateZDC()
  Float_t kDimZEMPb  = 0.15*(TMath::Sqrt(2.));  // z-dimension of the Pb slice
  Float_t kFibRadZEM = 0.0315; 			// External fiber radius (including cladding)
  Int_t   fDivZEM[3] = {92, 0, 20}; 		// Divisions for EM detector
  Float_t fDimZEM[6] = {fZEMLength, 3.5, 3.5, 45., 0., 0.}; // Dimensions of EM detector
  Float_t fFibZEM2 = fDimZEM[2]/TMath::Sin(fDimZEM[3]*kDegrad)-kFibRadZEM;
  Float_t fFibZEM[3] = {0., 0.0275, fFibZEM2};  // Fibers for EM calorimeter

//if(!fOnlyZEM){
  // Parameters for hadronic calorimeters geometry
  // NB -> parameters used ONLY in CreateZDC()
  Float_t fGrvZN[3] = {0.03, 0.03, 50.};  // Grooves for neutron detector
  Float_t fGrvZP[3] = {0.04, 0.04, 75.};  // Grooves for proton detector
  Int_t   fDivZN[3] = {11, 11, 0};  	  // Division for neutron detector
  Int_t   fDivZP[3] = {7, 15, 0};  	  // Division for proton detector
  Int_t   fTowZN[2] = {2, 2};  		  // Tower for neutron detector
  Int_t   fTowZP[2] = {4, 1};  		  // Tower for proton detector


  
  //-- Create calorimeters geometry
  
  // -------------------------------------------------------------------------------
  //--> Neutron calorimeter (ZN) 
  
  TVirtualMC::GetMC()->Gsvolu("ZNEU", "BOX ", idtmed[1], fDimZN, 3); // Passive material  
  TVirtualMC::GetMC()->Gsvolu("ZNF1", "TUBE", idtmed[3], fFibZN, 3); // Active material
  TVirtualMC::GetMC()->Gsvolu("ZNF2", "TUBE", idtmed[4], fFibZN, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZNF3", "TUBE", idtmed[4], fFibZN, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZNF4", "TUBE", idtmed[3], fFibZN, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZNG1", "BOX ", idtmed[12], fGrvZN, 3); // Empty grooves 
  TVirtualMC::GetMC()->Gsvolu("ZNG2", "BOX ", idtmed[12], fGrvZN, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZNG3", "BOX ", idtmed[12], fGrvZN, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZNG4", "BOX ", idtmed[12], fGrvZN, 3); 
  
  // Divide ZNEU in towers (for hits purposes) 
  
  TVirtualMC::GetMC()->Gsdvn("ZNTX", "ZNEU", fTowZN[0], 1); // x-tower 
  TVirtualMC::GetMC()->Gsdvn("ZN1 ", "ZNTX", fTowZN[1], 2); // y-tower
  
  //-- Divide ZN1 in minitowers 
  //  fDivZN[0]= NUMBER OF FIBERS PER TOWER ALONG X-AXIS, 
  //  fDivZN[1]= NUMBER OF FIBERS PER TOWER ALONG Y-AXIS
  //  (4 fibres per minitower) 
  
  TVirtualMC::GetMC()->Gsdvn("ZNSL", "ZN1 ", fDivZN[1], 2); // Slices 
  TVirtualMC::GetMC()->Gsdvn("ZNST", "ZNSL", fDivZN[0], 1); // Sticks
  
  // --- Position the empty grooves in the sticks (4 grooves per stick)
  Float_t dx = fDimZN[0] / fDivZN[0] / 4.;
  Float_t dy = fDimZN[1] / fDivZN[1] / 4.;
  
  TVirtualMC::GetMC()->Gspos("ZNG1", 1, "ZNST", 0.-dx, 0.+dy, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZNG2", 1, "ZNST", 0.+dx, 0.+dy, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZNG3", 1, "ZNST", 0.-dx, 0.-dy, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZNG4", 1, "ZNST", 0.+dx, 0.-dy, 0., 0, "ONLY");
  
  // --- Position the fibers in the grooves 
  TVirtualMC::GetMC()->Gspos("ZNF1", 1, "ZNG1", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZNF2", 1, "ZNG2", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZNF3", 1, "ZNG3", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZNF4", 1, "ZNG4", 0., 0., 0., 0, "ONLY");
  
  // --- Position the neutron calorimeter in ZDC 
  // -- Rotation of ZDCs
  Int_t irotzdc;
  TVirtualMC::GetMC()->Matrix(irotzdc, 90., 180., 90., 90., 180., 0.);
  //
  TVirtualMC::GetMC()->Gspos("ZNEU", 1, "ZDCC", fPosZNC[0], fPosZNC[1], fPosZNC[2]-fDimZN[2], irotzdc, "ONLY");
  //Ch debug
  if(TMath::Abs(fPosZNC[1])>0.) printf("\n ZNC placed at  y %f\n",fPosZNC[1]);
  printf("\n ZNC -> %f < z < %f cm\n",fPosZNC[2],fPosZNC[2]-2*fDimZN[2]);

  // --- Position the neutron calorimeter in ZDC2 (left line) 
  // -- No Rotation of ZDCs
  TVirtualMC::GetMC()->Gspos("ZNEU", 2, "ZDCA", fPosZNA[0], fPosZNA[1], fPosZNA[2]+fDimZN[2], 0, "ONLY");
  //Ch debug
  if(TMath::Abs(fPosZNA[1])>0.) printf("\n ZNA placed at  y %f\n",fPosZNA[1]);
  printf("\n ZNA -> %f < z < %f cm\n",fPosZNA[2],fPosZNA[2]+2*fDimZN[2]);


  // -------------------------------------------------------------------------------
  //--> Proton calorimeter (ZP)  
  
  TVirtualMC::GetMC()->Gsvolu("ZPRO", "BOX ", idtmed[2], fDimZP, 3); // Passive material
  TVirtualMC::GetMC()->Gsvolu("ZPF1", "TUBE", idtmed[3], fFibZP, 3); // Active material
  TVirtualMC::GetMC()->Gsvolu("ZPF2", "TUBE", idtmed[4], fFibZP, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZPF3", "TUBE", idtmed[4], fFibZP, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZPF4", "TUBE", idtmed[3], fFibZP, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZPG1", "BOX ", idtmed[12], fGrvZP, 3); // Empty grooves 
  TVirtualMC::GetMC()->Gsvolu("ZPG2", "BOX ", idtmed[12], fGrvZP, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZPG3", "BOX ", idtmed[12], fGrvZP, 3); 
  TVirtualMC::GetMC()->Gsvolu("ZPG4", "BOX ", idtmed[12], fGrvZP, 3); 
    
  //-- Divide ZPRO in towers(for hits purposes) 
  
  TVirtualMC::GetMC()->Gsdvn("ZPTX", "ZPRO", fTowZP[0], 1); // x-tower 
  TVirtualMC::GetMC()->Gsdvn("ZP1 ", "ZPTX", fTowZP[1], 2); // y-tower
  
  
  //-- Divide ZP1 in minitowers 
  //  fDivZP[0]= NUMBER OF FIBERS ALONG X-AXIS PER MINITOWER, 
  //  fDivZP[1]= NUMBER OF FIBERS ALONG Y-AXIS PER MINITOWER
  //  (4 fiber per minitower) 
  
  TVirtualMC::GetMC()->Gsdvn("ZPSL", "ZP1 ", fDivZP[1], 2); // Slices 
  TVirtualMC::GetMC()->Gsdvn("ZPST", "ZPSL", fDivZP[0], 1); // Sticks
  
  // --- Position the empty grooves in the sticks (4 grooves per stick)
  dx = fDimZP[0] / fTowZP[0] / fDivZP[0] / 2.;
  dy = fDimZP[1] / fTowZP[1] / fDivZP[1] / 2.;
  
  TVirtualMC::GetMC()->Gspos("ZPG1", 1, "ZPST", 0.-dx, 0.+dy, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZPG2", 1, "ZPST", 0.+dx, 0.+dy, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZPG3", 1, "ZPST", 0.-dx, 0.-dy, 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZPG4", 1, "ZPST", 0.+dx, 0.-dy, 0., 0, "ONLY");
  
  // --- Position the fibers in the grooves 
  TVirtualMC::GetMC()->Gspos("ZPF1", 1, "ZPG1", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZPF2", 1, "ZPG2", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZPF3", 1, "ZPG3", 0., 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZPF4", 1, "ZPG4", 0., 0., 0., 0, "ONLY");
  

  // --- Position the proton calorimeter in ZDCC
  TVirtualMC::GetMC()->Gspos("ZPRO", 1, "ZDCC", fPosZPC[0], fPosZPC[1], fPosZPC[2]-fDimZP[2], irotzdc, "ONLY");
  //Ch debug
  printf("\n ZPC -> %f < z < %f cm\n",fPosZPC[2],fPosZPC[2]-2*fDimZP[2]);
  if(TMath::Abs(fPosZPC[1])>0.) printf("\n ZNA placed at  y %f\n",fPosZPC[1]);
  
  // --- Position the proton calorimeter in ZDCA
  // --- No rotation 
  TVirtualMC::GetMC()->Gspos("ZPRO", 2, "ZDCA", fPosZPA[0], fPosZPA[1], fPosZPA[2]+fDimZP[2], 0, "ONLY");
  //Ch debug
  printf("\n ZPA -> %f < z < %f cm\n",fPosZPA[2],fPosZPA[2]+2*fDimZP[2]);  
  if(TMath::Abs(fPosZPA[1])>0.) printf("\n ZNA placed at  y %f\n",fPosZPA[1]);
//}    
  
  // -------------------------------------------------------------------------------
  // -> EM calorimeter (ZEM)  
  
  TVirtualMC::GetMC()->Gsvolu("ZEM ", "PARA", idtmed[10], fDimZEM, 6);

  Int_t irot1, irot2;
  TVirtualMC::GetMC()->Matrix(irot1,0.,0.,90.,90.,-90.,0.); 		       // Rotation matrix 1  
  TVirtualMC::GetMC()->Matrix(irot2,180.,0.,90.,fDimZEM[3]+90.,90.,fDimZEM[3]);// Rotation matrix 2
  //printf("irot1 = %d, irot2 = %d \n", irot1, irot2);
  
  TVirtualMC::GetMC()->Gsvolu("ZEMF", "TUBE", idtmed[3], fFibZEM, 3); 	// Active material

  TVirtualMC::GetMC()->Gsdvn("ZETR", "ZEM ", fDivZEM[2], 1); 	     	// Tranches 
  
  dimPb[0] = kDimZEMPb;					// Lead slices 
  dimPb[1] = fDimZEM[2];
  dimPb[2] = fDimZEM[1];
  //dimPb[3] = fDimZEM[3]; //controllare
  dimPb[3] = 90.-fDimZEM[3]; //originale
  dimPb[4] = 0.;
  dimPb[5] = 0.;
  TVirtualMC::GetMC()->Gsvolu("ZEL0", "PARA", idtmed[5], dimPb, 6);
  TVirtualMC::GetMC()->Gsvolu("ZEL1", "PARA", idtmed[5], dimPb, 6);
  TVirtualMC::GetMC()->Gsvolu("ZEL2", "PARA", idtmed[5], dimPb, 6);
  
  // --- Position the lead slices in the tranche 
  Float_t zTran = fDimZEM[0]/fDivZEM[2]; 
  Float_t zTrPb = -zTran+kDimZEMPb;
  TVirtualMC::GetMC()->Gspos("ZEL0", 1, "ZETR", zTrPb, 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZEL1", 1, "ZETR", kDimZEMPb, 0., 0., 0, "ONLY");
  
  // --- Vacuum zone (to be filled with fibres)
  dimVoid[0] = (zTran-2*kDimZEMPb)/2.;
  dimVoid[1] = fDimZEM[2];
  dimVoid[2] = fDimZEM[1];
  dimVoid[3] = 90.-fDimZEM[3];
  dimVoid[4] = 0.;
  dimVoid[5] = 0.;
  TVirtualMC::GetMC()->Gsvolu("ZEV0", "PARA", idtmed[10], dimVoid,6);
  TVirtualMC::GetMC()->Gsvolu("ZEV1", "PARA", idtmed[10], dimVoid,6);
  
  // --- Divide the vacuum slice into sticks along x axis
  TVirtualMC::GetMC()->Gsdvn("ZES0", "ZEV0", fDivZEM[0], 3); 
  TVirtualMC::GetMC()->Gsdvn("ZES1", "ZEV1", fDivZEM[0], 3); 
  
  // --- Positioning the fibers into the sticks
  TVirtualMC::GetMC()->Gspos("ZEMF", 1,"ZES0", 0., 0., 0., irot2, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZEMF", 1,"ZES1", 0., 0., 0., irot2, "ONLY");
  
  // --- Positioning the vacuum slice into the tranche
  //Float_t displFib = fDimZEM[1]/fDivZEM[0];
  TVirtualMC::GetMC()->Gspos("ZEV0", 1,"ZETR", -dimVoid[0], 0., 0., 0, "ONLY");
  TVirtualMC::GetMC()->Gspos("ZEV1", 1,"ZETR", -dimVoid[0]+zTran, 0., 0., 0, "ONLY");

  // --- Positioning the ZEM into the ZDC - rotation for 90 degrees  
  // NB -> ZEM is positioned in ALIC (instead of in ZDC) volume
  TVirtualMC::GetMC()->Gspos("ZEM ", 1,"ALIC", -fPosZEM[0], fPosZEM[1], fPosZEM[2]+fDimZEM[0], irot1, "ONLY");
  
  // Second EM ZDC (same side w.r.t. IP, just on the other side w.r.t. beam pipe)
  TVirtualMC::GetMC()->Gspos("ZEM ", 2,"ALIC", fPosZEM[0], fPosZEM[1], fPosZEM[2]+fDimZEM[0], irot1, "ONLY");
  
  // --- Adding last slice at the end of the EM calorimeter 
  Float_t zLastSlice = fPosZEM[2]+kDimZEMPb+2*fDimZEM[0];
  TVirtualMC::GetMC()->Gspos("ZEL2", 1,"ALIC", fPosZEM[0], fPosZEM[1], zLastSlice, irot1, "ONLY");
  //Ch debug
  //printf("\n ZEM lenght = %f cm\n",2*fZEMLength);
  printf("\n ZEM -> %f < z < %f cm\n\n",fPosZEM[2],fPosZEM[2]+2*fZEMLength+kDimZEMPb);
  
}
 
//_____________________________________________________________________________
void AliZDCv4::CreateMaterials()
{
  //
  // Create Materials for the Zero Degree Calorimeter
  //
  Float_t ubuf[1]={0.};
  Float_t wmat[3]={0.,0,0}, a[3]={0.,0,0}, z[3]={0.,0,0};

  // --- W alloy -> ZN passive material
  a[0] = 183.85;
  a[1] = 55.85;
  a[2] = 58.71;
  z[0] = 74.;
  z[1] = 26.;
  z[2] = 28.;
  wmat[0] = .93;
  wmat[1] = .03;
  wmat[2] = .04;
  AliMixture(1, "WALL", a, z, 17.6, 3, wmat);

  // --- Brass (CuZn)  -> ZP passive material
  a[0] = 63.546;
  a[1] = 65.39;
  z[0] = 29.;
  z[1] = 30.;
  wmat[0] = .63;
  wmat[1] = .37;
  AliMixture(2, "BRASS", a, z, 8.48, 2, wmat);
  
  // --- SiO2 
  a[0] = 28.086;
  a[1] = 15.9994;
  z[0] = 14.;
  z[1] = 8.;
  wmat[0] = 1.;
  wmat[1] = 2.;
  AliMixture(3, "SIO2", a, z, 2.64, -2, wmat);  
  
  // --- Lead 
  ubuf[0] = 1.12;
  AliMaterial(5, "LEAD", 207.19, 82., 11.35, .56, 0., ubuf, 1);

  // --- Copper (energy loss taken into account)
  ubuf[0] = 1.10;
  AliMaterial(6, "COPP0", 63.54, 29., 8.96, 1.43, 0., ubuf, 1);

  // --- Copper 
  AliMaterial(9, "COPP1", 63.54, 29., 8.96, 1.43, 0., ubuf, 1);
  
  // --- Iron (energy loss taken into account)
  AliMaterial(7, "IRON0", 55.85, 26., 7.87, 1.76, 0., ubuf, 1);
  
  // --- Iron (no energy loss)
  ubuf[0] = 1.1;
  AliMaterial(8, "IRON1", 55.85, 26., 7.87, 1.76, 0., ubuf, 1);
  
  // --- Tatalum 
  ubuf[0] = 1.1;
  AliMaterial(13, "TANT", 183.84, 74., 19.3, 0.35, 0., ubuf, 1);
    
  // ---------------------------------------------------------  
  Float_t aResGas[3]={1.008,12.0107,15.9994};
  Float_t zResGas[3]={1.,6.,8.};
  Float_t wResGas[3]={0.28,0.28,0.44};
  Float_t dResGas = 3.2E-14;

  // --- Vacuum (no magnetic field) 
  AliMixture(10, "VOID", aResGas, zResGas, dResGas, 3, wResGas);
  
  // --- Vacuum (with magnetic field) 
  AliMixture(11, "VOIM", aResGas, zResGas, dResGas, 3, wResGas);
  
  // --- Air (no magnetic field)
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  //
  AliMixture(12, "Air    $", aAir, zAir, dAir, 4, wAir);
  
  // --- Aluminum 
  AliMaterial(14, "ALUM", 26.98, 13., 2.7, 8.9, 0., ubuf, 1);
  
  // --- Carbon 
  AliMaterial(15, "GRAPH", 12.011, 6., 2.265, 18.8, 49.9);
  
  // ---  Definition of tracking media: 
  
  // --- Tantalum = 1 ; 
  // --- Brass = 2 ; 
  // --- Fibers (SiO2) = 3 ; 
  // --- Fibers (SiO2) = 4 ; 
  // --- Lead = 5 ; 
  // --- Copper (with high thr.)= 6 ;
  // --- Copper (with low thr.)=  9;
  // --- Iron (with energy loss) = 7 ; 
  // --- Iron (without energy loss) = 8 ; 
  // --- Vacuum (no field) = 10 
  // --- Vacuum (with field) = 11 
  // --- Air (no field) = 12 
  
  // **************************************************** 
  //     Tracking media parameters
  //
  Float_t epsil  = 0.01;   // Tracking precision, 
  Float_t stmin  = 0.01;   // Min. value 4 max. step (cm)
  Float_t stemax = 1.;     // Max. step permitted (cm) 
  Float_t tmaxfd = 0.;     // Maximum angle due to field (degrees) 
  Float_t tmaxfdv = 0.1;   // Maximum angle due to field (degrees) 
  Float_t deemax = -1.;    // Maximum fractional energy loss
  Float_t nofieldm = 0.;   // Max. field value (no field)
  Float_t fieldm = 45.;    // Max. field value (with field)
  Int_t isvol = 0;         // ISVOL =0 -> not sensitive volume
  Int_t isvolActive = 1;   // ISVOL =1 -> sensitive volume
  Int_t inofld = 0;        // IFIELD=0 -> no magnetic field
  Int_t ifield = 2;        // IFIELD=2 -> magnetic field defined in AliMagFC.h
  // *****************************************************
  
  AliMedium(1, "ZWALL", 1, isvolActive, inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(2, "ZBRASS",2, isvolActive, inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(3, "ZSIO2", 3, isvolActive, inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(4, "ZQUAR", 3, isvolActive, inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(5, "ZLEAD", 5, isvolActive, inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(6, "ZCOPP", 6, isvol, 	inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(7, "ZIRON", 7, isvol, 	inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(8, "ZIRONN",8, isvol, 	inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(9, "ZCOPL", 6, isvol, 	inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(10,"ZVOID",10, isvol, 	inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(11,"ZVOIM",11, isvol, 	ifield, fieldm,   tmaxfdv,stemax, deemax, epsil, stmin);
  AliMedium(12,"ZAIR", 12, isvolActive, inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(13,"ZALUM",13, isvol, inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(14,"ZGRAPH",14,isvolActive, inofld, nofieldm, tmaxfd, stemax, deemax, epsil, stmin);

} 

//_____________________________________________________________________________
void AliZDCv4::AddAlignableVolumes() const
{
 //
 // Create entries for alignable volumes associating the symbolic volume
 // name with the corresponding volume path. Needs to be syncronized with
 // eventual changes in the geometry.
 //
 //if(fOnlyZEM) return;
 
 TString volpath1 = "ALIC_1/ZDCC_1/ZNEU_1";
 TString volpath2 = "ALIC_1/ZDCC_1/ZPRO_1";
 TString volpath3 = "ALIC_1/ZDCA_1/ZNEU_2";
 TString volpath4 = "ALIC_1/ZDCA_1/ZPRO_2";

 TString symname1="ZDC/NeutronZDC_C";
 TString symname2="ZDC/ProtonZDC_C";
 TString symname3="ZDC/NeutronZDC_A";
 TString symname4="ZDC/ProtonZDC_A";

 if(!gGeoManager->SetAlignableEntry(symname1.Data(),volpath1.Data()))
     AliFatal(Form("Alignable entry %s not created. Volume path %s not valid",   symname1.Data(),volpath1.Data()));

 if(!gGeoManager->SetAlignableEntry(symname2.Data(),volpath2.Data()))
     AliFatal(Form("Alignable entry %s not created. Volume path %s not valid",   symname2.Data(),volpath2.Data()));

 if(!gGeoManager->SetAlignableEntry(symname3.Data(),volpath3.Data()))
     AliFatal(Form("Alignable entry %s not created. Volume path %s not valid",   symname1.Data(),volpath1.Data()));

 if(!gGeoManager->SetAlignableEntry(symname4.Data(),volpath4.Data()))
     AliFatal(Form("Alignable entry %s not created. Volume path %s not valid",   symname2.Data(),volpath2.Data()));

}


//_____________________________________________________________________________
void AliZDCv4::Init()
{
 InitTables();
  Int_t *idtmed = fIdtmed->GetArray();  
  //
  fMedSensZN     = idtmed[1];  // Sensitive volume: ZN passive material
  fMedSensZP     = idtmed[2];  // Sensitive volume: ZP passive material
  fMedSensF1     = idtmed[3];  // Sensitive volume: fibres type 1
  fMedSensF2     = idtmed[4];  // Sensitive volume: fibres type 2
  fMedSensZEM    = idtmed[5];  // Sensitive volume: ZEM passive material
  fMedSensTDI    = idtmed[6];  // Sensitive volume: TDI Cu shield
  fMedSensPI     = idtmed[7];  // Sensitive volume: beam pipes
  fMedSensLumi   = idtmed[9];  // Sensitive volume: luminometer
  fMedSensGR     = idtmed[12]; // Sensitive volume: air into the grooves
  fMedSensVColl  = idtmed[14]; // Sensitive volume: collimator vertical jaws
}

//_____________________________________________________________________________
void AliZDCv4::InitTables()
{
 //
 // Read light tables for Cerenkov light production parameterization 
 //

  Int_t k, j;
  int read=1;

  //  --- Reading light tables for ZN 
  char *lightfName1 = gSystem->ExpandPathName("$ALICE_ROOT/ZDC/light22620362207s");
  FILE *fp1 = fopen(lightfName1,"r");
  if(fp1 == NULL){
     printf("Cannot open light table from file %s \n",lightfName1);
     return;
  }
  else{
    for(k=0; k<fNalfan; k++){
      for(j=0; j<fNben; j++){
       read = fscanf(fp1,"%f",&fTablen[0][k][j]);
       if(read==0) AliDebug(3, " Error in reading light table 1");
      }
    }
    fclose(fp1);
  }
  char *lightfName2 = gSystem->ExpandPathName("$ALICE_ROOT/ZDC/light22620362208s");
  FILE *fp2 = fopen(lightfName2,"r");
  if(fp2 == NULL){
     printf("Cannot open light table from file %s \n",lightfName2);
     return;
  }  
  else{
    for(k=0; k<fNalfan; k++){
      for(j=0; j<fNben; j++){
       read = fscanf(fp2,"%f",&fTablen[1][k][j]);
       if(read==0) AliDebug(3, " Error in reading light table 2");
      }
    }
    fclose(fp2);
  }
  char *lightfName3 = gSystem->ExpandPathName("$ALICE_ROOT/ZDC/light22620362209s");
  FILE *fp3 = fopen(lightfName3,"r");
  if(fp3 == NULL){
     printf("Cannot open light table from file %s \n",lightfName3);
     return;
  }
  else{
    for(k=0; k<fNalfan; k++){
      for(j=0; j<fNben; j++){
       read = fscanf(fp3,"%f",&fTablen[2][k][j]);
       if(read==0) AliDebug(3, " Error in reading light table 3");
      }
    }
    fclose(fp3);
  }
  char *lightfName4 = gSystem->ExpandPathName("$ALICE_ROOT/ZDC/light22620362210s");
  FILE *fp4 = fopen(lightfName4,"r");
  if(fp4 == NULL){
     printf("Cannot open light table from file %s \n",lightfName4);
     return;
  }
  else{
    for(k=0; k<fNalfan; k++){
      for(j=0; j<fNben; j++){
       read = fscanf(fp4,"%f",&fTablen[3][k][j]);
       if(read==0) AliDebug(3, " Error in reading light table 4");
      }
    }
    fclose(fp4);
  }
    
  //  --- Reading light tables for ZP and ZEM
  char *lightfName5 = gSystem->ExpandPathName("$ALICE_ROOT/ZDC/light22620552207s");
  FILE *fp5 = fopen(lightfName5,"r");
  if(fp5 == NULL){
     printf("Cannot open light table from file %s \n",lightfName5);
     return;
  }
  else{
    for(k=0; k<fNalfap; k++){
      for(j=0; j<fNbep; j++){
       read = fscanf(fp5,"%f",&fTablep[0][k][j]);
       if(read==0) AliDebug(3, " Error in reading light table 5");
      }
    }
    fclose(fp5);
  }
  char *lightfName6 = gSystem->ExpandPathName("$ALICE_ROOT/ZDC/light22620552208s");
  FILE *fp6 = fopen(lightfName6,"r");
  if(fp6 == NULL){
     printf("Cannot open light table from file %s \n",lightfName6);
     return;
  }
  else{
    for(k=0; k<fNalfap; k++){
      for(j=0; j<fNbep; j++){
       read = fscanf(fp6,"%f",&fTablep[1][k][j]);
       if(read==0) AliDebug(3, " Error in reading light table 6");
      }
    }
    fclose(fp6);
  }
  char *lightfName7 = gSystem->ExpandPathName("$ALICE_ROOT/ZDC/light22620552209s");
  FILE *fp7 = fopen(lightfName7,"r");
  if(fp7 == NULL){
     printf("Cannot open light table from file %s \n",lightfName7);
     return;
  }
  else{
    for(k=0; k<fNalfap; k++){
      for(j=0; j<fNbep; j++){
       read = fscanf(fp7,"%f",&fTablep[2][k][j]);
       if(read==0) AliDebug(3, " Error in reading light table 7");
      }
    }
   fclose(fp7);
  }
  char *lightfName8 = gSystem->ExpandPathName("$ALICE_ROOT/ZDC/light22620552210s");
  FILE *fp8 = fopen(lightfName8,"r");
  if(fp8 == NULL){
     printf("Cannot open light table from file %s \n",lightfName8);
     return;
  }
  else{
    for(k=0; k<fNalfap; k++){
      for(j=0; j<fNbep; j++){
       read = fscanf(fp8,"%f",&fTablep[3][k][j]);
       if(read==0) AliDebug(3, " Error in reading light table 8");
      }
    }
   fclose(fp8);
  }

}
//_____________________________________________________________________________
void AliZDCv4::StepManager()
{
  //
  // Routine called at every step in the Zero Degree Calorimeters
  //
  Int_t   vol[2]={0,0}, ibeta=0, ialfa=0, ibe=0, nphe=0;
  Float_t x[3]={0.,0.,0.}, xdet[3]={999.,999.,999.}, um[3]={0.,0.,0.}, ud[3]={0.,0.,0.};
  Double_t s[3]={0.,0.,0.}, p[4]={0.,0.,0.,0.};
  Float_t destep=0., be=0., out=0.;
  //
  Float_t hits[14];
  for(int j=0; j<14; j++) hits[j]=-999.;
  const char *knamed = (fMC)->CurrentVolName();
  Int_t  mid = fMC->CurrentMedium();
  fMC->TrackPosition(s[0],s[1],s[2]);
  //printf("\tZDC::StepManager\t volume %s medium %d (x,y,z) = (%f, %f, %f)\n", knamed, mid, s[0], s[1], s[2]);
  
  // Study spectator protons distributions at TDI z
  /*fMC->TrackPosition(s[0],s[1],s[2]);
  if(s[2]>=7813.30 && s[2]<=8353.30){
     printf(" \t**** particle in vol. %s\n ",knamed);  
     fMC->TrackMomentum(p[0], p[1], p[2], p[3]);
     Int_t ctrack = gAlice->GetMCApp()->GetCurrentTrackNumber();
     TParticle *cpart = gAlice->GetMCApp()->Particle(ctrack);
     printf("\t TDIpc  %d %f %f %f %f \n", cpart->GetPdgCode(), s[0],s[1],s[2],p[3]);
  }
  else if(s[2]>=8353.30 && s[2]<=8403.30){
     fMC->TrackMomentum(p[0], p[1], p[2], p[3]);
     Int_t ctrack = gAlice->GetMCApp()->GetCurrentTrackNumber();
     TParticle *cpart = gAlice->GetMCApp()->Particle(ctrack);
     printf("\t TDIpc  %d %f %f %f %f \n", cpart->GetPdgCode(), s[0],s[1],s[2],p[3]);
  }
  else if(s[2]>8403.30){ 
     fMC->StopTrack();
     return;
  }*/
  //
  // --- This part is for no shower developement in beam pipe, TDI, VColl
  // If particle interacts with beam pipe, TDI, VColl -> return
  if(fNoShower==1 && ((mid == fMedSensPI) || (mid == fMedSensTDI) ||  
  	(mid == fMedSensVColl) || (mid == fMedSensLumi))){ 
   
   // Avoid to stop track in skewed cones between recombination chambers or separate beam pipes and ZDC (Jan 2015)
   if((strncmp(knamed,"QA27",4)) && (strncmp(knamed,"QA28",4)) &&
   	(strncmp(knamed,"QA29",4))){ // true if it is NOT in QA27 || QA28 || QA29
    
    // If option NoShower is set -> StopTrack
    //printf(" \t**** particle in vol. %s\n ",knamed);  
    
    Int_t ipr = 0; 
      fMC->TrackPosition(s[0],s[1],s[2]);
      //printf("\t\t(x,y,z) = (%f, %f, %f)\n", s[0], s[1], s[2]);
      fMC->TrackMomentum(p[0], p[1], p[2], p[3]);
      
      if(mid == fMedSensPI){
        if(!strncmp(knamed,"YMQ",3)){
	  if(s[2]<0) fpLostITC += 1;
	  else fpLostITA += 1;
	  ipr=1;
        }
	else if(!strncmp(knamed,"QA02",4)){
	  if((s[2]>26.15 && s[2]<32.52) || (s[2]>34.80 && s[2]<40.30) || 
	     (s[2]>41.30 && s[2]<46.80) || (s[2]>50.15 && s[2]<56.52)) fpLostITA += 1;
	}
	else if(!strncmp(knamed,"YD1",3)){
	  if(s[2]<0) fpLostD1C += 1;
	  else fpLostD1A += 1;
	  ipr=1;
	}
	else if(!strncmp(knamed,"QA03",4)) fpLostD1A += 1;
	else if(!strncmp(knamed,"QT02",4)) fpLostD1C += 1;
	else if(!strncmp(knamed,"QTD",3) || strncmp(knamed,"Q13T",4)) fpLostTDI += 1;
      }
      else if(mid == fMedSensTDI){  // fMedSensTDI also involves beam screen inside IT and D1
        if(!strncmp(knamed,"QBS1",4) || !strncmp(knamed,"QBS2",4) || // beam screens inside Q1
	   !strncmp(knamed,"QBS3",4) || !strncmp(knamed,"QBS4",4) || // beam screens inside Q3
	   !strncmp(knamed,"QBS5",4) || !strncmp(knamed,"QBS6",4)    // beam screens inside Q2A/Q2B
	){
	  if(s[2]<0) fpLostITC += 1;
	  else fpLostITA += 1;
        }
        else if(!strncmp(knamed,"MD1",3)){
          if(s[2]<0) fpLostD1C += 1;
          else  fpLostD1A += 1;
        }
	else if(!strncmp(knamed,"QTD",3)) fpLostTDI += 1;
        ipr=1;
      }
      else if(mid == fMedSensVColl){ 
        if(!strncmp(knamed,"QCVC",4)) fpcVCollC++;
 	else if(!strncmp(knamed,"QCVA",4))  fpcVCollA++;
	ipr=1;
      }
      //
      //printf("\t Particle: mass = %1.3f, E = %1.3f GeV, pz = %1.2f GeV -> stopped in volume %s\n", fMC->TrackMass(), p[3], p[2], knamed);
      //
      if(ipr<0){
        printf("\n\t **********************************\n");
        printf("\t ********** Side C **********\n");
        printf("\t # of particles in IT = %d\n",fpLostITC);
        printf("\t # of particles in D1 = %d\n",fpLostD1C);
        printf("\t # of particles in VColl = %d\n",fpcVCollC);
        printf("\t ********** Side A **********\n");
        printf("\t # of particles in IT = %d\n",fpLostITA);
        printf("\t # of particles in D1 = %d\n",fpLostD1A);
        printf("\t # of particles in TDI = %d\n",fpLostTDI);
        printf("\t # of particles in VColl = %d\n",fpcVCollA);
        printf("\t **********************************\n");
      }
      fMC->StopTrack();
      return;
     }
  }
  
  if((mid == fMedSensZN) || (mid == fMedSensZP) ||
     (mid == fMedSensGR) || (mid == fMedSensF1) ||
     (mid == fMedSensF2) || (mid == fMedSensZEM)){

    //Ch. debug
    //printf(" ** pc. track %d in vol. %s \n",gAlice->GetMCApp()->GetCurrentTrackNumber(), knamed);
    
  //Particle coordinates 
    fMC->TrackPosition(s[0],s[1],s[2]);
    for(int j=0; j<=2; j++) x[j] = s[j];
    hits[0] = x[0];
    hits[1] = x[1];
    hits[2] = x[2];

  // Determine in which ZDC the particle is
    if(!strncmp(knamed,"ZN",2)){
          if(x[2]<0.) vol[0]=1; // ZNC (dimuon side)
	  else if(x[2]>0.) vol[0]=4; //ZNA
    }
    else if(!strncmp(knamed,"ZP",2)){ 
          if(x[2]<0.) vol[0]=2; //ZPC (dimuon side)
	  else if(x[2]>0.) vol[0]=5; //ZPA  
    }
    else if(!strncmp(knamed,"ZE",2)) vol[0]=3; //ZEM
    
    // February 2015: Adding TrackReference
    if(fSwitchOnTrackRef==kTRUE && (fMC->IsTrackEntering() || fMC->IsTrackExiting())) {
       AliTrackReference* trackRef = AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kZDC);
       if(vol[0]>0){
         trackRef->SetUserId(vol[0]);
         //printf("Adding track reference for track %d in vol. %d\n", gAlice->GetMCApp()->GetCurrentTrackNumber(), vol[0]);
       }
    }
  
  // Determine in which quadrant the particle is
    if(vol[0]==1){	//Quadrant in ZNC
      // Calculating particle coordinates inside ZNC
      xdet[0] = x[0]-fPosZNC[0];
      xdet[1] = x[1]-fPosZNC[1];
      // Calculating quadrant in ZN
      if(xdet[0]<=0.){
        if(xdet[1]<=0.) vol[1]=1;
	else vol[1]=3;
      }
      else if(xdet[0]>0.){
        if(xdet[1]<=0.) vol[1]=2;
        else vol[1]=4;
      }
    }
    else if(vol[0]==2){	//Quadrant in ZPC
      // Calculating particle coordinates inside ZPC
      xdet[0] = x[0]-fPosZPC[0];
      xdet[1] = x[1]-fPosZPC[1];
      if(xdet[0]>=fDimZP[0])  xdet[0]=fDimZP[0]-0.01;
      if(xdet[0]<=-fDimZP[0]) xdet[0]=-fDimZP[0]+0.01;
      // Calculating tower in ZP
      Float_t xqZP = xdet[0]/(fDimZP[0]/2.);
      for(int i=1; i<=4; i++){
         if(xqZP>=(i-3) && xqZP<(i-2)){
 	   vol[1] = i;
 	   break;
 	 }
      }
    }
    //
    // Quadrant in ZEM: vol[1] = 1 -> particle in 1st ZEM (placed at x = 8.5 cm)
    // 		 	vol[1] = 2 -> particle in 2nd ZEM (placed at x = -8.5 cm)
    else if(vol[0] == 3){	
      if(x[0]>0.){
        vol[1] = 1;
        // Particle x-coordinate inside ZEM1
        xdet[0] = x[0]-fPosZEM[0];
      }
      else{
   	vol[1] = 2;
        // Particle x-coordinate inside ZEM2
        xdet[0] = x[0]+fPosZEM[0];
      }
      xdet[1] = x[1]-fPosZEM[1];
    }
    //
    else if(vol[0]==4){	//Quadrant in ZNA
      // Calculating particle coordinates inside ZNA
      xdet[0] = x[0]-fPosZNA[0];
      xdet[1] = x[1]-fPosZNA[1];
      // Calculating quadrant in ZNA
      if(xdet[0]>=0.){
        if(xdet[1]<=0.) vol[1]=1;
	else vol[1]=3;
      }
      else if(xdet[0]<0.){
        if(xdet[1]<=0.) vol[1]=2;
        else vol[1]=4;
      }
    }    
    //
    else if(vol[0]==5){	//Quadrant in ZPA
      // Calculating particle coordinates inside ZPA
      xdet[0] = x[0]-fPosZPA[0];
      xdet[1] = x[1]-fPosZPA[1];
      if(xdet[0]>=fDimZP[0])  xdet[0]=fDimZP[0]-0.01;
      if(xdet[0]<=-fDimZP[0]) xdet[0]=-fDimZP[0]+0.01;
      // Calculating tower in ZP
      Float_t xqZP = -xdet[0]/(fDimZP[0]/2.);
      for(int i=1; i<=4; i++){
         if(xqZP>=(i-3) && xqZP<(i-2)){
 	   vol[1] = i;
 	   break;
 	 }
      }
    }    
    if((vol[1]!=1) && (vol[1]!=2) && (vol[1]!=3) && (vol[1]!=4))
      AliError(Form(" WRONG tower for det %d: tow %d with xdet=(%f, %f)\n",
		vol[0], vol[1], xdet[0], xdet[1]));
    // Ch. debug
    //printf("\t *** det %d vol %d xdet(%f, %f)\n",vol[0], vol[1], xdet[0], xdet[1]);
    
    // Store impact point and kinetic energy of the ENTERING particle
    
    if(fMC->IsTrackEntering()){
      //Particle energy
      fMC->TrackMomentum(p[0],p[1],p[2],p[3]);
      hits[3] = p[3];
      
      // Impact point on ZDC
      // X takes into account the LHC x-axis sign
      // which is opposite to positive x on detector front face
      // for side A detectors (ZNA and ZPA)  
      if(vol[0]==4 || vol[0]==5){
        hits[4] = -xdet[0];
      }
      else{
        hits[4] = xdet[0];
      }
      hits[5] = xdet[1];
      hits[6] = 0;
      hits[7] = 0;
      hits[8] = 0;
      hits[9] = 0;
      //
      Int_t curTrackN = gAlice->GetMCApp()->GetCurrentTrackNumber();
      TParticle *part = gAlice->GetMCApp()->Particle(curTrackN);
      hits[10] = part->GetPdgCode();
      hits[11] = 0;
      hits[12] = 1.0e09*fMC->TrackTime(); // in ns!
      hits[13] = part->Eta();
      //
      if(fFindMother){
         Int_t imo = part->GetFirstMother();
         //printf(" tracks: pc %d -> mother %d \n", curTrackN,imo); 
         
         int trmo = imo;
         TParticle *pmot = 0x0;
         Bool_t isChild = kFALSE;
         if(imo>-1){
           pmot = gAlice->GetMCApp()->Particle(imo);
      	   trmo = pmot->GetFirstMother();
      	   isChild = kTRUE;
           while(trmo!=-1){
      	      pmot = gAlice->GetMCApp()->Particle(trmo);
              //printf("  **** pc %d -> mother %d \n", trch,trmo); 
      	      trmo = pmot->GetFirstMother();
      	   }
         }
      
         if(isChild && pmot){
             hits[6]  = 1;
             hits[11] = pmot->GetPdgCode();
      	     hits[13] = pmot->Eta();
         }
      }
      

      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);

      if(fNoShower==1){
        if(vol[0]==1){
          fnDetectedC += 1;
          printf("	### Particle in ZNC\n\n");
        }
        else if(vol[0]==2){
          fpDetectedC += 1;
          printf("	### Particle in ZPC\n\n");
        }
        else if(vol[0]==3) printf("	### Particle in ZEM\n\n");	  
        else if(vol[0]==4){
          fnDetectedA += 1;
          printf("	### Particle in ZNA\n\n");	  
        }
        else if(vol[0]==5){
          fpDetectedA += 1;
          printf("	### Particle in ZPA\n\n"); 	 
        }
    	//
        //printf("\t Track %d: x %1.2f y %1.2f z %1.2f  E %1.2f GeV pz = %1.2f GeV in volume %s -> det %d\n", 
           //gAlice->GetMCApp()->GetCurrentTrackNumber(),x[0],x[1],x[2],p[3],p[2],knamed, vol[0]);
        //printf("\t Track %d: pc %d  E %1.2f GeV pz = %1.2f GeV in volume %s -> det %d\n", 
           //gAlice->GetMCApp()->GetCurrentTrackNumber(),part->GetPdgCode(),p[3],p[2],knamed, vol[0]);
        //
        fMC->StopTrack();
        return;
      }
    }
    	   
    // Particle energy loss
    if(fMC->Edep() != 0){
      hits[9] = fMC->Edep();
      hits[7] = 0.;
      hits[8] = 0.;
      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
    }
 

  // *** Light production in fibres 
  if((mid == fMedSensF1) || (mid == fMedSensF2)){

     //Select charged particles
     if((destep=fMC->Edep())){

       // Particle velocity
       Float_t beta = 0.;
       fMC->TrackMomentum(p[0],p[1],p[2],p[3]);
       Float_t ptot=TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
       if(p[3] > 0.00001) beta =  ptot/p[3];
       else return;
       if(beta<0.67)return;
       else if((beta>=0.67) && (beta<=0.75)) ibeta = 0;
       else if((beta>0.75)  && (beta<=0.85)) ibeta = 1;
       else if((beta>0.85)  && (beta<=0.95)) ibeta = 2;
       else if(beta>0.95) ibeta = 3;
 
       // Angle between particle trajectory and fibre axis
       // 1 -> Momentum directions
       um[0] = p[0]/ptot;
       um[1] = p[1]/ptot;
       um[2] = p[2]/ptot;
       fMC->Gmtod(um,ud,2);
       // 2 -> Angle < limit angle
       Double_t alfar = TMath::ACos(ud[2]);
       Double_t alfa = alfar*kRaddeg;
       if(alfa>=110.) return;
       //
       ialfa = Int_t(1.+alfa/2.);
 
       // Distance between particle trajectory and fibre axis
       fMC->TrackPosition(s[0],s[1],s[2]);
       for(int j=0; j<=2; j++){
   	  x[j] = s[j];
       }
       fMC->Gmtod(x,xdet,1);
       if(TMath::Abs(ud[0])>0.00001){
         Float_t dcoeff = ud[1]/ud[0];
         be = TMath::Abs((xdet[1]-dcoeff*xdet[0])/TMath::Sqrt(dcoeff*dcoeff+1.));
       }
       else{
         be = TMath::Abs(ud[0]);
       }
 
       ibe = Int_t(be*1000.+1);
  
       //Looking into the light tables 
       Float_t charge = 0.;
       Int_t curTrackN = gAlice->GetMCApp()->GetCurrentTrackNumber();
       TParticle *part = gAlice->GetMCApp()->Particle(curTrackN);
       Int_t pdgCode = part->GetPdgCode();
       if(pdgCode<10000) charge = fMC->TrackCharge();
       else{
          float z = (pdgCode/10000-100000);
          charge = TMath::Abs(z);
          //printf(" PDG %d   charge %f\n",pdgCode,charge);
       } 
       
       if(vol[0]==1 || vol[0]==4) {	// (1)  ZN fibres
         if(ibe>fNben) ibe=fNben;
         out =  charge*charge*fTablen[ibeta][ialfa][ibe];
	 nphe = gRandom->Poisson(out);
	 // Ch. debug
         //if(ibeta==3) printf("\t %f \t %f \t %f\n",alfa, be, out);
	 //printf("\t ibeta = %d, ialfa = %d, ibe = %d -> nphe = %d\n\n",ibeta,ialfa,ibe,nphe);
	 if(mid == fMedSensF1){
	   hits[7] = nphe;  	//fLightPMQ
	   hits[8] = 0;
	   hits[9] = 0;
	   AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
	 }
	 else{
	   hits[7] = 0;
	   hits[8] = nphe;	//fLightPMC
	   hits[9] = 0;
	   AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
	 }
       } 
       else if(vol[0]==2 || vol[0]==5) {// (2) ZP fibres
         if(ibe>fNbep) ibe=fNbep;
         out =  charge*charge*fTablep[ibeta][ialfa][ibe];
	 nphe = gRandom->Poisson(out);
	 if(mid == fMedSensF1){
	   hits[7] = nphe;  	//fLightPMQ
	   hits[8] = 0;
	   hits[9] = 0;
	   AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
	 }
	 else{
	   hits[7] = 0;
	   hits[8] = nphe;	//fLightPMC
	   hits[9] = 0;
	   AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
	 }
       } 
       else if(vol[0]==3) {	// (3) ZEM fibres
         if(ibe>fNbep) ibe=fNbep;
         out =  charge*charge*fTablep[ibeta][ialfa][ibe];
     fMC->TrackPosition(s[0],s[1],s[2]);
	 Float_t xalic[3];
         for(int j=0; j<3; j++){
            xalic[j] = s[j];
         }
	 // z-coordinate from ZEM front face 
	 // NB-> fPosZEM[2]+fZEMLength = -1000.+2*10.3 = 979.69 cm
	 Float_t z = -xalic[2]+fPosZEM[2]+2*fZEMLength-xalic[1];
	 //z = xalic[2]-fPosZEM[2]-fZEMLength-xalic[1]*(TMath::Tan(45.*kDegrad));
         //printf("	fPosZEM[2]+2*fZEMLength = %f", fPosZEM[2]+2*fZEMLength);
         //
	 // Parametrization for light guide uniformity
         // NEW!!! Light guide tilted @ 51 degrees
         Float_t guiPar[4]={0.31,-0.0006305,0.01337,0.8895};
	 Float_t guiEff = guiPar[0]*(guiPar[1]*z*z+guiPar[2]*z+guiPar[3]);
	 out = out*guiEff;
	 nphe = gRandom->Poisson(out);
         //printf("	out*guiEff = %f	nphe = %d", out, nphe);
	 if(vol[1] == 1){
	   hits[7] = 0;  	
	   hits[8] = nphe;	//fLightPMC (ZEM1)
	   hits[9] = 0;
	   AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
	 }
	 else{
	   hits[7] = nphe;  	//fLightPMQ (ZEM2)
	   hits[8] = 0;		
	   hits[9] = 0;
	   AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
	 }
       }
     }
   }
  }
}
