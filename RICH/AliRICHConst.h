#ifndef RICHConst_H
#define RICHConst_H
/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE RICH chambers geometry
//--------------------------------------------------------------------
//
const Float_t zend = 511.+0.15-2*0.001; // z-out position of first chamber


/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE RICH  Electronics Parameters
//--------------------------------------------------------------------
//
//

const Float_t adc_satm  = 4096; // dynamic range (10 bits)
const Float_t zero_supm = 6.; // zero suppression
const Float_t sig_noise = 500.; // electronics noise (no. of electrons)

const Int_t kMaxNeighbours = 24; // max number of neighbours
#endif
