/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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

//===============================================================
// To be modified for hadron correction using particle ID and track merging
// Author : magali.estienne@ires.in2p3.fr
//===============================================================


// --- AliRoot header files ---
#include "AliJetHadronCorrectionv0.h"

const Int_t maxVariant = 8;  // size eta grid
const Int_t nVec = 10;       // size momentum grid
const Int_t nPol = 4;        // number coefficients of polinom
static Double_t etaGrid[maxVariant]={ 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.67};
static Double_t pGrid[nVec]={ 0.2, 0.5, 1.0, 2.0, 3.0, 5.0,10.0,15.0,25.0,40.0};
// c[][] - first index for eta, second for momentum 
static Double_t c[maxVariant][nPol] ={
    {1.305705e-01, 3.725653e-01, -1.219962e-02, 1.806235e-04},
    {1.296153e-01, 3.827408e-01, -1.238640e-02, 1.822804e-04},
    {1.337690e-01, 3.797454e-01, -1.245227e-02, 1.848243e-04},
    {1.395796e-01, 3.623994e-01, -9.196803e-03, 1.243278e-04},
    {1.457184e-01, 3.753655e-01, -1.035324e-02, 1.473447e-04},
    {1.329164e-01, 4.219044e-01, -1.310515e-02, 1.948883e-04},
    {8.136581e-02, 4.646087e-01, -1.531917e-02, 2.274749e-04},
    {1.119836e-01, 4.262497e-01, -1.160125e-02, 1.628738e-04} };

ClassImp(AliJetHadronCorrectionv0)

AliJetHadronCorrectionv0* AliJetHadronCorrectionv0::fHadrCorr = 0;

AliJetHadronCorrectionv0::AliJetHadronCorrectionv0(const char *name,const char *title) 
                           :AliJetHadronCorrection(name, title)
{
  fHadrCorr = this;
}

AliJetHadronCorrectionv0*
AliJetHadronCorrectionv0::Instance()
{
  fHadrCorr = new AliJetHadronCorrectionv0();
  return fHadrCorr;
}

Double_t 
AliJetHadronCorrectionv0::GetEnergy(Double_t pmom, Double_t eta, Int_t /*gid*/)
{
  Int_t iEta=0; // index 
  Double_t etaw = TMath::Abs(eta);
  if(etaw > etaGrid[maxVariant-1]) etaw = etaGrid[maxVariant-1];
  for(Int_t i=0; i<maxVariant; i++) if(eta>=etaGrid[i]) {iEta = i; break;}

  Double_t e[2], y, pw = pmom;
  if(pmom > pGrid[nVec-1]) pw = pGrid[nVec-1];
  for(Int_t i=0; i<2; i++){ // e for two eta value
    e[i] = c[iEta][0];
    y = 1.;
    for(Int_t j=1; j<nPol; j++){
      y *= pw;
      e[i] += c[iEta][j]*y;
    }
    if(i==0) iEta ++;
  }

  Double_t deta = etaGrid[iEta] - etaGrid[iEta-1];
  Double_t a = (e[1] - e[0])/deta; // slope 
  Double_t energy = e[0] + a*(eta-etaGrid[iEta-1]);
  return energy;
}
