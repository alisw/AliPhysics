#ifndef _GenTypeDefs_H
#define _GenTypeDefs_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

typedef enum {Pion, Kaon, Eta, Omega, Etaprime, Phi, Baryon, pion_p, kaon_p, phi_p, jpsi_p, upsilon_p, charm_p, beauty_p}
Param_t;

typedef enum
{ charm, beauty, charm_unforced, beauty_unforced, jpsi, jpsi_chi, mb}
Process_t;

typedef enum
{ semielectronic, dielectron, semimuonic, dimuon,
  b_jpsi_dimuon, b_jpsi_dielectron, 
  b_psip_dimuon, b_psip_dielectron, pitomu, katomu, nodecay, hadronicD, all}
Decay_t;

typedef enum
{
    DO_Set_1=1006,
    GRV_LO=5005,
    GRV_HO=5006,
    MRS_D_minus=3031,
    MRS_D0=3030,
    MRS_G=3041,
    CTEQ_2pM=4024,
    CTEQ_4M=4034
}
StrucFunc_t;

typedef enum
{
    analog,
    non_analog
}
Weighting_t;
#endif


