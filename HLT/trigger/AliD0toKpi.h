// @(#) $Id$

#ifndef AliD0toKpi_H
#define AliD0toKpi_H
  /*******************************************************************
   *                                                                 *
   *  Reconstructed D^0 -> K^- pi^+ candidate clss                   * 
   *                                                                 *
   *  origin: A. Dainese    andrea.dainese@pd.infn.it                *  
   *******************************************************************/

#include <TMath.h>
#include <TVector3.h>

// particle masses
const Double_t kMD0 = 1.8645;  // D0  mass
const Double_t kMK  = 0.49368; // K+  mass
const Double_t kMPi = 0.13957; // Pi+ mass

//  --- TOF tagging probabilities --- 
//  central HIJING
//  B = 0.4 T
//  tracking errors in TPC included
//  With TRD
//
//  *** TOF default cuts for PbPb *** 
//
/*
//  PIONS
const Int_t kPiBins = 20;
const Double_t kPiBinWidth = 0.125;
const Double_t kPiTagPi[kPiBins] = {0.000000,0.591070,0.908567,0.848575,0.770351,0.722956,0.706888,0.695148,0.687808,0.686050,0.689962,0.690553,0.687770,0.687936,0.691099,0.696442,0.704000,0.708352,0.709448,0.712480};
const Double_t kPiTagNid[kPiBins] = {1.000000,0.408930,0.091433,0.151425,0.229649,0.277044,0.293112,0.304852,0.312192,0.313950,0.310038,0.309447,0.312230,0.312064,0.308901,0.303558,0.296000,0.291648,0.290552,0.287520};
// KAONS
const Int_t kKBins = 20;
const Double_t kKBinWidth = 0.125;
const Double_t kKTagK[kKBins] = {0.000000,0.000000,0.039577,0.164515,0.308500,0.366434,0.412601,0.450548,0.472678,0.497862,0.517206,0.529481,0.550700,0.522383,0.491087,0.472946,0.380244,0.360244,0.314122,0.221509};
const Double_t kKTagPi[kKBins] = {0.000000,0.204099,0.576011,0.617346,0.390463,0.282998,0.221205,0.178188,0.161715,0.145718,0.139147,0.140506,0.133797,0.165070,0.187582,0.216733,0.340557,0.377610,0.426489,0.521731};
const Double_t kKTagNid[kKBins] = {1.000000,0.795901,0.384412,0.218139,0.301037,0.350568,0.366194,0.371264,0.365607,0.356420,0.343647,0.330012,0.315503,0.312547,0.321331,0.310321,0.279199,0.262147,0.259390,0.256760};
// PROTONS
const Int_t kPBins = 36;
const Double_t kPBinWidth = 0.125;
const Double_t kPTagP[kPBins] = {0.000000,0.002889,0.038156,0.085494,0.192829,0.297104,0.377435,0.429866,0.486685,0.501615,0.511019,0.524895,0.541178,0.559191,0.577053,0.585619,0.574745,0.572519,0.580219,0.586295,0.591467,0.601517,0.609334,0.611916,0.606260,0.596428,0.576680,0.556772,0.551359,0.547600,0.532903,0.506974,0.468212,0.309770,0.150187,0.000000};
const Double_t kPTagPi[kPBins] = {0.000000,0.200825,0.576272,0.691605,0.494897,0.341871,0.254065,0.204454,0.160916,0.152080,0.141468,0.127294,0.111641,0.098054,0.088475,0.077363,0.074684,0.074395,0.075396,0.076185,0.076013,0.068483,0.060666,0.064324,0.072480,0.084812,0.107060,0.129468,0.137381,0.143640,0.160837,0.189266,0.230528,0.391470,0.553553,0.706240};
const Double_t kPTagNid[kPBins] = {1.000000,0.796286,0.385572,0.222901,0.312273,0.361025,0.368500,0.365680,0.352399,0.346305,0.347514,0.347811,0.347181,0.342755,0.334471,0.337018,0.350571,0.353086,0.344385,0.337520,0.332520,0.330000,0.330000,0.323760,0.321260,0.318760,0.316260,0.313760,0.311260,0.308760,0.306260,0.303760,0.301260,0.298760,0.296260,0.293760};
//
// ***** Cuts_1 ******
//
//  PIONS
const Int_t kPiBins = 10;
const Double_t kPiBinWidth = 0.250;
const Double_t kPiTagPi[kPiBins] = {0.040920,0.402873,0.561528,0.549381,0.481544,0.445760,0.393588,0.369643,0.345340,0.323165};
const Double_t kPiTagNid[kPiBins] = {0.959080,0.597127,0.438472,0.450619,0.518456,0.554240,0.606412,0.630357,0.654660,0.676835};
//  KAONS
const Int_t kKBins = 10;
const Double_t kKBinWidth = 0.250;
const Double_t kKTagK[kKBins] = {0.000000,0.101255,0.397662,0.467586,0.517008,0.555023,0.601888,0.589805,0.609983,0.412181};
const Double_t kKTagPi[kKBins] = {0.102049,0.066665,0.031082,0.018655,0.016786,0.014284,0.017703,0.023592,0.039782,0.123654};
const Double_t kKTagNid[kKBins] = {0.897951,0.832080,0.571256,0.513759,0.466206,0.430694,0.380410,0.386603,0.350236,0.464165};
//  PROTONS
const Int_t kPBins = 9;
const Double_t kPBinWidth = 0.500;
const Double_t kPTagP[kPBins] = {0.017940,0.350681,0.535286,0.583264,0.562935,0.560524,0.545992,0.598060,0.351245};
const Double_t kPTagPi[kPBins] = {0.045999,0.034981,0.016826,0.013887,0.003778,0.008493,-0.000000,0.000000,0.000000};
const Double_t kPTagNid[kPBins] = {0.936062,0.614337,0.447888,0.402849,0.433287,0.430983,0.454008,0.401940,0.648755};
//
//  **** Cuts_2 ****
//
//  PIONS
const Int_t kPiBins = 10;
const Double_t kPiBinWidth = 0.250;
const Double_t kPiTagPi[kPiBins] = {0.041204,0.405645,0.563019,0.550815,0.485707,0.454333,0.406930,0.387574,0.376332,0.349019};
const Double_t kPiTagNid[kPiBins] = {0.958796,0.594355,0.436981,0.449185,0.514293,0.545667,0.593070,0.612426,0.623668,0.650981};
//  KAONS
const Int_t kKBins = 10;
const Double_t kKBinWidth = 0.250;
const Double_t kKTagK[kKBins] = {0.000000,0.101255,0.397662,0.467586,0.517008,0.555023,0.588611,0.581941,0.596722,0.370963};
const Double_t kKTagPi[kKBins] = {0.102049,0.071067,0.032453,0.021062,0.017905,0.014284,0.017703,0.031456,0.053042,0.206090};
const Double_t kKTagNid[kKBins] = {0.897951,0.827677,0.569885,0.511352,0.465087,0.430694,0.393687,0.386603,0.350236,0.422947};
//  PROTONS
const Int_t kPBins = 9;
const Double_t kPBinWidth = 0.500;
const Double_t kPTagP[kPBins] = {0.017940,0.350681,0.535286,0.583264,0.562935,0.560524,0.545992,0.598060,0.351245};
const Double_t kPTagPi[kPBins] = {0.046459,0.036068,0.017878,0.013887,0.003778,0.008493,0.015166,0.000000,0.000000};
const Double_t kPTagNid[kPBins] = {0.935602,0.613251,0.446836,0.402849,0.433287,0.430983,0.438842,0.401940,0.648755};
//
//  **** Cuts_3 ****
//
//  PIONS
const Int_t kPiBins = 10;
const Double_t kPiBinWidth = 0.250;
const Double_t kPiTagPi[kPiBins] = {0.185562,0.587715,0.605042,0.606765,0.600542,0.597826,0.552951,0.531017,0.482590,0.423347};
const Double_t kPiTagNid[kPiBins] = {0.814438,0.412285,0.394958,0.393235,0.399458,0.402174,0.447049,0.468983,0.517410,0.576653};
//  KAONS
const Int_t kKBins = 10;
const Double_t kKBinWidth = 0.250;
const Double_t kKTagK[kKBins] = {0.000000,0.101255,0.397662,0.467586,0.517008,0.555023,0.584185,0.542621,0.517159,0.267917};
const Double_t kKTagPi[kKBins] = {0.102049,0.220749,0.078618,0.046939,0.032453,0.022446,0.035405,0.055048,0.066302,0.247308};
const Double_t kKTagNid[kKBins] = {0.897951,0.677996,0.523719,0.485475,0.450539,0.422532,0.380410,0.402331,0.416538,0.484774};
//  PROTONS
const Int_t kPBins = 9;
const Double_t kPBinWidth = 0.500;
const Double_t kPTagP[kPBins] = {0.017940,0.350681,0.535286,0.583264,0.562935,0.560524,0.545992,0.598060,0.351245};
const Double_t kPTagPi[kPBins] = {0.136617,0.076481,0.032075,0.020831,0.003778,0.016986,0.030333,0.000000,0.000000};
const Double_t kPTagNid[kPBins] = {0.845444,0.572838,0.432639,0.395905,0.433287,0.422491,0.423675,0.401940,0.648755};
//
// *** Cuts_4 ***
//
//  PIONS
const Int_t kPiBins = 10;
const Double_t kPiBinWidth = 0.250;
const Double_t kPiTagPi[kPiBins] = {0.211421,0.652184,0.624421,0.614727,0.610777,0.628015,0.631520,0.630324,0.637551,0.575235};
const Double_t kPiTagNid[kPiBins] = {0.788579,0.347816,0.375579,0.385273,0.389223,0.371985,0.368480,0.369676,0.362449,0.424765};
//  KAONS
const Int_t kKBins = 10;
const Double_t kKBinWidth = 0.250;
const Double_t kKTagK[kKBins] = {0.000000,0.101255,0.397662,0.467586,0.517008,0.555023,0.584185,0.519029,0.464117,0.247308};
const Double_t kKTagPi[kKBins] = {0.102049,0.289930,0.101930,0.057771,0.040286,0.028567,0.053108,0.094369,0.066302,0.247308};
const Double_t kKTagNid[kKBins] = {0.897951,0.608815,0.500408,0.474643,0.442705,0.416410,0.362707,0.386603,0.469580,0.505383};
//  PROTONS
const Int_t kPBins = 9;
const Double_t kPBinWidth = 0.500;
const Double_t kPTagP[kPBins] = {0.017940,0.350681,0.535286,0.583264,0.562935,0.560524,0.545992,0.598060,0.351245};
const Double_t kPTagPi[kPBins] = {0.195955,0.094949,0.039962,0.026039,0.007556,0.016986,0.030333,0.000000,0.000000};
const Double_t kPTagNid[kPBins] = {0.786105,0.554370,0.424751,0.390697,0.429508,0.422491,0.423675,0.401940,0.648755};
//
// pp PYTHIA 
//
// *** Same Cuts as PbPb_4
//
//  PIONS
const Int_t kPiBins = 10;
const Double_t kPiBinWidth = 0.250;
const Double_t kPiTagPi[kPiBins] = {0.194528,0.447097,0.603364,0.646413,0.647125,0.669157,0.688139,0.682564,0.689910,0.665710};
const Double_t kPiTagNid[kPiBins] = {0.805472,0.552903,0.396636,0.353587,0.352875,0.330843,0.311861,0.317436,0.310090,0.334290};
//  KAONS
const Int_t kKBins = 10;
const Double_t kKBinWidth = 0.250;
const Double_t kKTagK[kKBins] = {0.000000,0.115097,0.419704,0.516093,0.585443,0.599174,0.586021,0.646101,0.434783,0.273885};
const Double_t kKTagPi[kKBins] = {0.000000,0.001495,0.000000,0.000000,-0.000000,0.000000,0.032258,0.060572,0.101449,0.242038};
const Double_t kKTagNid[kKBins] = {1.000000,0.883408,0.580296,0.483907,0.414557,0.400826,0.381720,0.293327,0.463768,0.484076};
//  PROTONS
const Int_t kPBins = 9;
const Double_t kPBinWidth = 0.500;
const Double_t kPTagP[kPBins] = {0.006918,0.386968,0.613710,0.665152,0.634961,0.657711,0.703704,0.685714,0.235294};
const Double_t kPTagPi[kPBins] = {-0.000000,0.000000,0.000000,-0.000000,0.000000,0.000000,-0.000000,0.014286,-0.000000};
const Double_t kPTagNid[kPBins] = {0.993082,0.613032,0.386290,0.334848,0.365039,0.342289,0.296296,0.300000,0.764706};
*/
//
// *** cuts for pp ***
//
//  PIONS
const Int_t kPiBins = 10;
const Double_t kPiBinWidth = 0.250;
const Double_t kPiTagPi[kPiBins] = {0.194528,0.447097,0.603364,0.646413,0.647125,0.669157,0.688139,0.682564,0.689910,0.665710};
const Double_t kPiTagNid[kPiBins] = {0.805472,0.552903,0.396636,0.353587,0.352875,0.330843,0.311861,0.317436,0.310090,0.334290};
//  KAONS
const Int_t kKBins = 10;
const Double_t kKBinWidth = 0.250;
const Double_t kKTagK[kKBins] = {0.000000,0.173393,0.439690,0.519423,0.587025,0.605372,0.586021,0.650139,0.444444,0.299363};
const Double_t kKTagPi[kKBins] = {0.000000,0.001495,0.000000,-0.000000,-0.000000,0.000000,0.032258,0.060572,0.101449,0.242038};
const Double_t kKTagNid[kKBins] = {1.000000,0.825112,0.560310,0.480577,0.412975,0.394628,0.381720,0.289289,0.454106,0.458599};
//  PROTONS
const Int_t kPBins = 9;
const Double_t kPBinWidth = 0.500;
const Double_t kPTagP[kPBins] = {0.029404,0.438640,0.613710,0.665152,0.634961,0.657711,0.703704,0.685714,0.235294};
const Double_t kPTagPi[kPBins] = {0.000000,0.000000,0.000000,-0.000000,0.000000,0.000000,-0.000000,0.014286,-0.000000};
const Double_t kPTagNid[kPBins] = {0.970596,0.561360,0.386290,0.334848,0.365039,0.342289,0.296296,0.300000,0.764706};



//-----------------------------------------------------------------------------
class AliD0toKpi : public TObject {
 public:
  //
  AliD0toKpi();
  AliD0toKpi(Bool_t sgn,Int_t ev,  
	     Double_t v1[3],Double_t v2[3],Double_t dca,
	     Double_t mom[6],Double_t d0[2],
	     Int_t pdg[2],Int_t mum[2]);


  Double_t PtChild(Int_t child)  const {return TMath::Sqrt(fPx[child]*fPx[child]+fPy[child]*fPy[child]);} 
  Double_t PChild(Int_t child)   const {return TMath::Sqrt(fPx[child]*fPx[child]+fPy[child]*fPy[child]+fPz[child]*fPz[child]);} 
  Double_t EtaChild(Int_t child)    const;
  Int_t    GetPdgChild(Int_t child) const {return fPdg[child];}
  Int_t    GetMumId(Int_t child)    const {return fMum[child];}
  Double_t Getd0Child(Int_t child)  const {return fd0[child];}
  Double_t P()  const {return TMath::Sqrt(Pt()*Pt()+Pz()*Pz());} 
  Double_t Pt() const {return TMath::Sqrt(Px()*Px()+Py()*Py());}
  Double_t Px() const {return (fPx[0]+fPx[1]);}
  Double_t Py() const {return (fPy[0]+fPy[1]);}
  Double_t Pz() const {return (fPz[0]+fPz[1]);}
  Double_t Energy() const {return TMath::Sqrt(P()*P()+kMD0*kMD0);}
  Double_t Rapidity() const {return 0.5*TMath::Log((Energy()+Pz())/(Energy()-Pz()+1.e-13));}
  Double_t Eta() const;
  Double_t ChildrenRelAngle()       const; 
  Double_t d0d0() const { return fd0[0]*fd0[1]; } 
  Double_t GetDCA()  const {return 10000.*fDCA;}
  Double_t CPta() const;
  Double_t CPtaXY() const;
  void     DrawPIDinTOF() const;
  void     InvMass(Double_t&,Double_t&) const;
  Double_t qT() const;
  Double_t qL(Int_t child) const;
  Double_t Alpha() const {return (qL(0)-qL(1))/(qL(0)+qL(1));}
  void     CosThetaStar(Double_t&,Double_t&) const;

  Bool_t   Select(const Double_t* cuts,Int_t&,Int_t&) const;
  void     ComputeWgts();
  void     CorrectWgt4BR(Double_t factor);
  void     SetPtWgts();
  void     GetWgts(Double_t&,Double_t&,Option_t* sample) const;
  //
 private:
  //
  Bool_t   fSignal; // TRUE if signal, FALSE if background
  Int_t    fEvent;  // number of the event this D0 comes from
  Double_t fV1x; //
  Double_t fV1y; // position of the primary vertex of the event
  Double_t fV1z; //
  Double_t fV2x; //
  Double_t fV2y; // position of the reconstructed secondary vertex
  Double_t fV2z; //
  Double_t fDCA; // DCA of the two tracks

  Double_t fPx[2];  // 
  Double_t fPy[2];  // momenta of the two tracks
  Double_t fPz[2];  // at the reconstructed vertex  

  Double_t fd0[2];  //  impact parameters in the bending plane

  Int_t fPdg[2];  // PDG codes 
  Int_t fMum[2];  // PDG codes of the mothers


  Double_t fWgtAD0,fWgtAD0bar; //
  Double_t fWgtBD0,fWgtBD0bar; // weights for the 3 samples 
  Double_t fWgtCD0,fWgtCD0bar; // A: (K,Pi)+(K,?) B: (?,Pi) C: (?,?)



  Double_t LinearInterpolation(Double_t p,Int_t nBins,Double_t Bin,const Double_t *values) const;

  void     SetZeroWgts(Option_t* opt);

  ClassDef(AliD0toKpi,1)  // reconstructed D0 candidate class
};


#endif








