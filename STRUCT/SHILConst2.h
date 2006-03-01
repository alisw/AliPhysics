#ifndef SHILCONST2_H
#define SHILCONST2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */

/* $Id$ */

// angle of 2nd cone
  const Float_t kThetaOpen2   = 0.83*kDegrad;
// angle of 3rd cone
  const Float_t kThetaOpen3   = 0.83*kDegrad;
// angle of beam tube in second cone
  const Float_t kThetaOpenB   = 0.83*kDegrad;
//  const Float_t thetaOpenB = 0.84*kDegrad;
// inner lead Cone opening angle
//  const Float_t kThetaOpenPb  = 0.6*kDegrad;
// Outer Pb Cone opening angle
  const Float_t kThetaOpenPbO = 1.6*kDegrad;
// Start of lead cone
  const Float_t kZPb = 720.;   
// y-position of trigger wall
  const Float_t kZFilterIn       = 1471.;
  const Float_t kZFilterOut      = kZFilterIn+120.;
// end of 2-degree outer cone
//  const Float_t zConeE    = 30./TMath::Tan(accMin);
  const Float_t kZConeE    = 859.;
//
//
  Float_t dTubeS = 0.1;
  const Float_t kDInsuS=1.5;
  const Float_t kDEnveS=0.1;
  const Float_t kDProtS=0.1;
  const Float_t kDFreeS=0.00;

  Float_t dVacuS=dTubeS+kDInsuS+kDEnveS+kDProtS+kDFreeS;

//
// Radii and z-positions imposed by vacuum chamber layout
//
// FIRST SECTION
// delta_R for bellows
//  const Float_t dr11=2.65;
  const Float_t kDr11=4.1;
// delta_R for flange
  const Float_t kDr12=0.01;
// delta_R to catch up with cone
  const Float_t kDr13=2.0;

// flange length
  const Float_t kDF1=2.*3.9;
// bellow length
  const Float_t kDB1=18.482-kDF1/2-kDr11/6.;
// Flange position
  const Float_t kZvac2=577.;
  const Float_t kZvac1=kZvac2-kDF1/2-kDB1-kDr11/10.-kDr12;
  const Float_t kZvac3=kZvac2+kDF1/2+kDB1+kDr12+kDr13;
//  const Float_t kZvac4=648.;
  const Float_t kZvac4 = 618.;
  const Float_t kZvac41= 558.;
// Outer shield dimensions
  const Float_t kR11=15.4;
// Steel Envelope
  const Float_t kDRSteel2=4.;
//  According to design
//  const Float_t R21=20.3;
//  to avoid overlap with 2deg line
  const Float_t kR21= 19.5;
//
// 2nd Section
// 
//  const Float_t kZvac5=kZvac4+4.;
//
// 3rd Section
// 
  const Float_t kZvac6=711.;
  const Float_t kZvac8=1274.;
//const Float_t kDr21=2.263;
  const Float_t kDr21=.5;
  const Float_t kDr22=1.3;
  const Float_t kDr23=2.263;
  const Float_t kDB2=24.118;
  const Float_t kDF2=10.6;
  const Float_t kZvac7=kZvac8-kDF2/2-kDB2-kDr22-kDr21;
  const Float_t kZvac9=kZvac8+kDF2/2+kDB2+kDr22+kDr23;
//
// 4th Section
// 
  const Float_t kZvac10=1466.;
  const Float_t kZvac11=1800.;
  const Float_t kZvac12=1900.;

  const Float_t kR41=35.;
  const Float_t kR42=48.0;
  const Float_t kR43=110.;

//
// Vacuum System
//

// Bellow1
//
  const Float_t kRB1=5.97+0.4;
  const Float_t kHB1=2.25;
  const Float_t kLB1=1.2;
  const Float_t kEB1=0.04;
//
// Flange1
//
  const Float_t kRF1=8.87+0.4;
  const Float_t kDFlange=0.1;

//
// Bellow2
//
  const Float_t kRB2=16.35;
  const Float_t kHB2=2.25;
  const Float_t kLB2=2.32;
  const Float_t kEB2=0.05;
//
// Flange2
//
  const Float_t kRF2=18.5;

//
// Chamber positions
//

const Float_t kZch11 =  517.70;
const Float_t kZch12 =  553.70;
const Float_t kZch21 =  661.30;
const Float_t kZch22 =  709.90;
const Float_t kZch31 =  951.25;
const Float_t kZch32 = 1014.75;
const Float_t kZch41 = 1259.90;
const Float_t kZch42 = 1324.10;
const Float_t kZch51 = 1390.00;
const Float_t kZch52 = 1454.20;
#endif

