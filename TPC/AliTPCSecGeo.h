#ifndef TPCSecGeo_H
#define TPCSecGeo_H
/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE TPC sector geometry
//--------------------------------------------------------------------
//
const Float_t z_end = 250.; // position of the readout chamber
//
//   Lower sectors (numbers 1-12), opening angle in radians
//
const Float_t alpha_low=0.523598775; // 30 degrees
//
//   Upper sectors (numbers 13-36), opening angle in radians
//
const Float_t alpha_up=0.261799387; //  15 degrees
//
// Pad size 2.05 x 0.35 cm
//
const Float_t pad_pitch_l=2.05;
const Float_t pad_pitch_w=0.35;
//
//  number of pad rows per sector
//
const Int_t nrow_low = 23;
const Int_t nrow_up = 52;
//
// Lower sector, pad row radii
//
const Float_t pad_row_low[23]={
   89.445,   91.495,   93.545,   95.595,  97.645,  99.695,  101.745,  103.795,  
  105.845,  107.895,  109.945,  111.995, 114.045, 116.095,  118.145,  120.195,  
  122.245,  124.295,  126.345,  128.395, 130.445, 132.495,  134.545}; 
//
// Upper sector, pad row radii
//
const Float_t pad_row_up[52]={
143.72,  145.77,  147.82,  149.87,  151.92,  153.97,  156.02, 158.07,  
160.12,  162.17,  164.22,  166.27,  168.32,  170.37,  172.42, 174.47,  
176.52,  178.57,  180.62,  182.67,  184.72,  186.77,  188.82, 190.87,  
192.92,  194.97,  197.02,  199.07,  201.12,  203.17,  205.22, 207.27,  
209.32,  211.37,  213.42,  215.47,  217.52,  219.57,  221.62, 223.67,  
225.72,  227.77,  229.82,  231.87,  233.92,  235.97,  238.02, 240.07,  
242.12,  244.17,  246.22,  248.27};
//
// Lower sector, number of pads per row
//
const Int_t npads_low[23]={
 129,  133,  135,  139,  143,  145,  149,  151,  155,  157,  161,  165,
 167,  171,  173,  177,  179,  183,  187,  189,  193,  195,  199};
//
//  Upper sector, number of pads per row
//
const Int_t npads_up[52]={
  101,  103,  103,  105,  107,  109,  111,  111,  113,  115,  117,  117,
  119,  121,  123,  123,  125,  127,  129,  131,  131,  133,  135,  137,
  137,  139,  141,  143,  145,  145,  147,  149,  151,  151,  153,  155,
  157,  157,  159,  161,  163,  165,  165,  167,  169,  171,  171,  173,
  175,  177,  177,  179};

//
//  Number of wires per pad and wire-wire pitch
//
const Int_t nwires = 5;
const Float_t ww_pitch = 0.41;

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



/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE TPC Gas Parameters
//--------------------------------------------------------------------
//
//
//   Diffusion constants
//
const Float_t diff_t=2.2e-2; 
const Float_t diff_l=2.2e-2;  
//
//  Lorentz angle (omega_tau)
//
const Float_t omega_tau = 0.125;
//
//  Electron drift velocity
//
const Float_t v_drift=2.85e6;
/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE TPC  Electronics Parameters
//--------------------------------------------------------------------
//
//

const Float_t t_sample = 2.e-7; // sampling time
const Float_t fwhm = 2.5e-7; // width of the Preamp/Shaper function

//

const Float_t gas_gain = 1.e4; // gas gain
const Float_t q_el = 1.602e-19; // elementary charge

//

const Float_t adc_sat = 1023; // dynamic range (10 bits)
const Float_t zero_sup = 5.; // zero suppression
const Float_t sigma_noise = 500.; // electronics noise (no. of electrons)
const Float_t chip_gain = 24.; // conversion gain (mV/fC)
const Float_t dyn_range = 2000.; // output dynamic range (mV)

#endif
