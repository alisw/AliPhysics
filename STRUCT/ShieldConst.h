#ifndef SHIELDCONST_H
#define SHIELDCONST_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// z-positions defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// start of the absorber 
    const Float_t abs_d   = 90.;    
// end of the W-nose
    const Float_t z_nose  = 102.;
// end of the 5deg line below the TPC field cage
    const Float_t z_cone  = 285.;
//
    const Float_t abs_cc  = 315.; 
// start of inner opening cone
    const Float_t abs_c   = 362.;
//    const Float_t abs_c   = 300.;
// rear end of the absorber
    const Float_t abs_l   = 503.;
// thickness of rear shield
    const Float_t d_rear  =  35.;
//
// angles defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// angle of nose
    const Float_t theta1  = 24.*kDegrad;
// angle of second outer cone below field cage
    const Float_t theta2  = 5. *kDegrad;
// outer angler of W rear shield  
    const Float_t theta_r = 3. *kDegrad;
// max acceptance angle
    const Float_t acc_max = 10.*kDegrad;
// min acceptance angle
    const Float_t acc_min = 2. *kDegrad;     
// opening angle of inner shielding cone
    const Float_t theta_open  = 0.7*kDegrad;
    const Float_t theta_open1 = 1.1*kDegrad;
//    const Float_t theta_open1 = 0.70*kDegrad;
    const Float_t theta_open2 = 0.8*kDegrad;
    const Float_t theta_open3 = 0.9*kDegrad;
//
// thicknesses defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// steel envelope
    const Float_t d_steel = 1.;  
// poly-ethylene layer
    const Float_t d_poly  = 5.0;
//
// radii defining the absorber
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^
// absorber inner radius
    const Float_t r_abs    = 4.95;

    const Float_t epsilon = .01;
    Float_t zr;
// start of 2 deg cone defined by absorber inner radius
    Float_t z_2deg = r_abs/tan(acc_min);
    

// Index of heavy shield material 
    Int_t idHeavy=1660;

// y-position of trigger wall
  const Float_t zfil_in    = 1471.;
  const Float_t zfil_out   = zfil_in+120.;
// end of 2-degree outer cone
  const Float_t zcone_e    = 30./TMath::Tan(acc_min);
// end of opening cone
  const Float_t zcone_c    = 1800.;
  const Float_t z_out      = 1900.;
//
// Chamber positions
  const Float_t cz1= 511;
  const Float_t cz2= 686;
//
// Radii and z-positions imposed by vacuum chamber layout
//
// FIRST SECTION
// delta_R for bellows
  const Float_t dr11=1.8;
// delta_R for flange
  const Float_t dr12=0.875;
// delta_R to catch up with cone
  const Float_t dr13=0.525;
// flange length
  const Float_t dF1=2.*3.9;
// bellow length
  const Float_t dB1=18.482-dF1/2-dr11;
// Flange position
  const Float_t zvac2=518.;
  const Float_t zvac1=zvac2-dF1/2-dB1-dr11-dr12;
  const Float_t zvac3=zvac2+dF1/2+dB1+dr11+dr13;
  const Float_t zvac4=558.;
// Outer shield dimensions
  const Float_t R11=15.45;
// Steel Envelope
  const Float_t dRSteel1=2.;
  const Float_t dRSteel2=4.;
//  According to design
//  const Float_t R21=20.3;
//  to avoid overlap with 2deg line
  const Float_t R21=19.4;
//
// 2nd Section
// 

  const Float_t zvac5=zvac4+4.;

//
// 3rd Section
// 
  const Float_t zvac6=711.;
  const Float_t zvac8=1274.;
  const Float_t dr21=2.263;
  const Float_t dr22=1.3;
  const Float_t dr23=0.1;
  const Float_t dB2=24.118;
  const Float_t dF2=10.6;
  const Float_t zvac7=zvac8-dF2/2-dB2-dr22-dr21;
  const Float_t zvac9=zvac8+dF2/2+dB2+dr22+dr23;

//
// 4th Section
// 
  const Float_t zvac10=1466.;
  const Float_t zvac11=1800.;
  const Float_t zvac12=1900.;

  const Float_t R41=35.;
  const Float_t R42=50.;
  const Float_t R43=110.;

//
// Vacuum System
//
// beam pipe outer radius
// absorber inner radius
 
  const Float_t d_tube=0.1;
  const Float_t d_insu=0.9;
  const Float_t d_enve=0.1;
  const Float_t d_prot=0.2;
  const Float_t d_free=0.5;

  const Float_t d_vacu=d_tube+d_insu+d_enve+d_prot+d_free;
  const Float_t r_vacu=r_abs-d_vacu;


//
// Bellow1
//
  const Float_t rB1=5.5;
  const Float_t hB1=2.25;
  const Float_t lB1=0.77;
  const Float_t eB1=0.04;
//
// Flange1
//
  const Float_t rF1=8.5;
  const Float_t d_flange=0.1;

//
// Bellow2
//
  const Float_t rB2=15.35;
  const Float_t hB2=2.25;
  const Float_t lB2=2.32;
  const Float_t eB2=0.05;
//
// Flange2
//
  const Float_t rF2=18.8;

//
// Chamber positions
//
const Float_t dzch=10.;
const Float_t zch1=532.5;
const Float_t zch2=690.;

#endif
