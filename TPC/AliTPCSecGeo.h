#ifndef TPCSecGeo_H


//Some things from the old AliTPCSecGeo

const Float_t z_end = 250.; 
const Float_t alpha_low=0.523598775; // 30 degrees
const Float_t alpha_up=0.261799387; //  15 degrees




const Float_t q_el = 1.602e-19; // elementary charge
const Float_t adc_sat = 1023; // dynamic range (10 bits)
const Float_t dyn_range = 2000.; // output dynamic range (mV)

/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE TPC Cluster Parameters
//--------------------------------------------------------------------
//
//
// Sigma rphi
const Float_t a_rphi=0.41818e-2;
const Float_t b_rphi=0.17460e-4;
const Float_t c_rphi=0.30993e-2;
const Float_t d_rphi=0.41061e-3;
// Sigma z
const Float_t a_z=0.39614e-2;
const Float_t b_z=0.22443e-4;
const Float_t c_z=0.51504e-1;
// Cluster width in rphi
const Float_t ac_rphi=0.18322;
const Float_t bc_rphi=0.59551e-3;
const Float_t cc_rphi=0.60952e-1;
// Cluster width in z
const Float_t ac_z=0.19081;
const Float_t bc_z=0.55938e-3;
const Float_t cc_z=0.30428;
//

#define TPCSecGeo_H

#endif
