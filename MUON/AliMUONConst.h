#ifndef MUONConst_H
#define MUONConst_H
/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE MUON chambers geometry
//--------------------------------------------------------------------
//
const Float_t zend = 511.+0.15-2*0.001; // z-out position of first chamber


/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE MUON  Electronics Parameters
//--------------------------------------------------------------------
//
//

const Float_t adc_satm  = 1024; // dynamic range (10 bits)
const Int_t zero_supm = 6; // zero suppression
const Float_t sig_noise = 500.; // electronics noise (no. of electrons)
const Float_t kScale=1.;


/////////////////////////////////////////////////////////////////////////////
//
//---------------------------------------------------------------------
//   ALICE MUON  segmentation Parameters
//--------------------------------------------------------------------
//
//
const Int_t kMaxNeighbours = 24; // max number of neighbours

#endif

