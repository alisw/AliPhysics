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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//   The TRD particle identification class                                   //
//                                                                           //
//   Its main purposes are:                                                  //
//      - Creation and bookkeeping of the propability distributions          //
//      - Assignment of a e/pi propability to a given track based on         //
//        the LQ method                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <math.h>

#include <TROOT.h>
#include <TH1.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>

#include "AliRun.h"
#include "AliTRDpidLQ.h"
#include "AliTRDcluster.h"
#include "AliTRDtrack.h"
#include "AliTRDtracker.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDpidLQ)

//_____________________________________________________________________________
AliTRDpidLQ::AliTRDpidLQ():AliTRDpid()
{
  //
  // AliTRDpidLQ default constructor
  // 

  fNMom   = 0;
  fMinMom = 0;
  fMaxMom = 0;
  fWidMom = 0;

  fNLh    = 0;
  fMinLh  = 0;
  fMaxLh  = 0;

  fHist   = NULL;

  fChargeMin   = 0;
  fNClusterMin = 0;

}

//_____________________________________________________________________________
AliTRDpidLQ::AliTRDpidLQ(const char* name, const char* title)
            :AliTRDpid(name,title)
{
  //
  // AliTRDpidLQ constructor
  //

  fNMom   = 0;
  fMinMom = 0;
  fMaxMom = 0;
  fWidMom = 0;

  fNLh    = 0;
  fMinLh  = 0;
  fMaxLh  = 0;

  fHist   = NULL;

  Init();

}

//_____________________________________________________________________________
AliTRDpidLQ::AliTRDpidLQ(const AliTRDpidLQ &p):AliTRDpid(p)
{
  //
  // AliTRDpidLQ copy constructor
  //

  ((AliTRDpidLQ &) p).Copy(*this);

}

//_____________________________________________________________________________
AliTRDpidLQ::~AliTRDpidLQ()
{
  //
  // AliTRDpidLQ destructor
  //

  if (fHist) {
    fHist->Delete();
    delete fHist;
  }

}

//_____________________________________________________________________________
AliTRDpidLQ &AliTRDpidLQ::operator=(const AliTRDpidLQ &p)
{
  //
  // Assignment operator
  //

  if (this != &p) ((AliTRDpidLQ &) p).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDpidLQ::Copy(TObject &p) const
{
  //
  // Copy function
  //

  fHist->Copy(*((AliTRDpidLQ &) p).fHist);

  ((AliTRDpidLQ &) p).fNMom        = fNMom;
  ((AliTRDpidLQ &) p).fMinMom      = fMinMom;
  ((AliTRDpidLQ &) p).fMaxMom      = fMaxMom;
  ((AliTRDpidLQ &) p).fWidMom      = fWidMom;
  ((AliTRDpidLQ &) p).fChargeMin   = fChargeMin;
  ((AliTRDpidLQ &) p).fNClusterMin = fNClusterMin;
  ((AliTRDpidLQ &) p).fNLh         = fNLh;
  ((AliTRDpidLQ &) p).fMinLh       = fMinLh;
  ((AliTRDpidLQ &) p).fMaxLh       = fMaxLh;

  AliTRDpid::Copy(p);

}

//_____________________________________________________________________________
Bool_t AliTRDpidLQ::Init()
{
  //
  // Initializes the PID object 
  //

  fChargeMin   = 0.0;
  fNClusterMin = 10;

  fNLh         = 100;
  fMinLh       =   0.0;
  fMaxLh       = 500.0;

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDpidLQ::AssignLikelihood(AliTRDtrack *t)
{
  //
  // Assigns the e / pi likelihood to a given track
  //

  const Int_t kNpla = AliTRDgeometry::Nplan();
  Float_t * charge = new Float_t[kNpla];
  Int_t   * nCluster = new Int_t[kNpla];

  Float_t  lhPi = 0;
  Float_t  lhEl = 0;

  Double_t pPi  = 1;
  Double_t pEl  = 1;
  Double_t pSum = 0;

  Int_t   indexEl;
  Int_t   indexPi;
  TH1F   *hTmpEl;
  TH1F   *hTmpPi;

  t->SetLikelihoodElectron(-1.);
  if (TMath::IsNaN(t->GetP())) return kFALSE;
  Float_t mom = t->GetP();

  // Calculate the total charge in each plane
  if (!SumCharge(t,charge,nCluster)) return kFALSE;

  indexEl = GetIndex(mom,kElectron);
  indexPi = GetIndex(mom,kPion);
  if ((indexEl > -1) && (indexPi > -1)) {
    hTmpEl = (TH1F *) fHist->UncheckedAt(indexEl);
    hTmpPi = (TH1F *) fHist->UncheckedAt(indexPi);
    for (Int_t ipla = 0; ipla < kNpla; ipla++) {
      Float_t chargePlane   = charge[ipla];
      Int_t   nClusterPlane = nCluster[ipla];
      if ((chargePlane   >   fChargeMin) &&
          (nClusterPlane > fNClusterMin)){
        Float_t chargeNorm = chargePlane / ((Float_t) nClusterPlane);
        if (chargeNorm < fMaxLh) {
          Int_t   ibinEl     = hTmpEl->FindBin(chargeNorm);
          Float_t pElPlane   = hTmpEl->GetBinContent(ibinEl);
          if (pElPlane > 0) { 
            pEl  = pEl * pElPlane;
          } 
          Int_t   ibinPi     = hTmpPi->FindBin(chargeNorm);
          Float_t pPiPlane   = hTmpPi->GetBinContent(ibinPi);
          if (pPiPlane > 0) { 
            pPi  = pPi * pPiPlane;
          } 
//           printf("          ipla = %d, nCluster = %d, charge = %f\n"
//                  ,ipla,nClusterPlane,chargeNorm);
//           printf("electron: pElPlane = %f, ibinEl = %d, pEl = %f\n"
// 	         ,pElPlane,ibinEl,pEl);
//           printf("    pion: pPiPlane = %f, ibinPi = %d, pPi = %f\n"
// 	         ,pPiPlane,ibinPi,pPi);
	}
      }
    }
  }
  else {
    delete [] charge;
    delete [] nCluster;
    return kTRUE;
  }

  pSum = pEl + pPi;
  if (pSum <= 0) {
    delete [] charge;
    delete [] nCluster;
    return kFALSE;
  }
  lhEl = pEl / pSum;
  lhPi = pPi / pSum;

//   printf("---- mom = %f, lhEl = %f, lhPi = %f\n",mom,lhEl,lhPi);

  // Assign the likelihoods 
  t->SetLikelihoodElectron(lhEl);

  delete [] charge;
  delete [] nCluster;
  return kTRUE;  

}

//_____________________________________________________________________________
Bool_t AliTRDpidLQ::CreateHistograms(Int_t   nmom, Float_t minmom, Float_t maxmom)
{
  //
  // Creates the histograms
  //

  Int_t imom;
  Int_t ipid;

  fNMom   = nmom;
  fMinMom = minmom;
  fMaxMom = maxmom;
  fWidMom = (maxmom - minmom) / ((Float_t) nmom);

  fHist = new TObjArray(kNpid * nmom);
  for (imom = 0; imom < nmom;  imom++) {
    for (ipid = 0; ipid < kNpid; ipid++) {

      Int_t index = GetIndex(imom,ipid);
      Char_t name[10];
      Char_t title[80];
      sprintf(name ,"hQ%03d",index);
      if (ipid == kElectron) {
        sprintf(title,"Q electron p-bin %03d",imom);
      }
      else {
        sprintf(title,"Q pion p-bin %03d",imom);
      }
      TH1F *hTmp = new TH1F(name,title,fNLh,fMinLh,fMaxLh);
      fHist->AddAt(hTmp,index);

    }
  }

  return kTRUE;

}

// //_____________________________________________________________________________
// Bool_t AliTRDpidLQ::FillSpectra(const AliTRDtrack *t)
// {
//   //
//   // Fills the energy loss distributions with track <t>.
//   //

//   Bool_t status = kTRUE;

//   if (TMath::IsNaN(t->GetP())) return kFALSE;

//   Float_t        mom     = t->GetP();
//   Int_t          ipid    = MCpid(t);
//   TH1F          *hTmp    = NULL;
//   AliTRDcluster *cluster = NULL;

//   Int_t index = GetIndex(mom,ipid);
//   if (index > -1) {
//     hTmp = (TH1F *) fHist->UncheckedAt(index);
//     // Loop through all clusters associated to this track
//     Int_t nClus = t->GetNclusters();
//     for (Int_t iClus = 0; iClus < nClus; iClus++) {
//       // Get a cluster
//       Int_t idxClus = t->GetClusterIndex(iClus);
//       if (!(cluster = (AliTRDcluster *) fClusterArray->UncheckedAt(idxClus))) {
//         status = kFALSE;
//         break;
//       } 
//       hTmp->Fill(cluster->GetQ());
//     }
//   }  

//   return status;

// }

//_____________________________________________________________________________
Bool_t AliTRDpidLQ::FillSpectra(const AliTRDtrack *t)
{
  //
  // Fills the energy loss distributions with track <t>.
  //

  const Int_t kNpla = AliTRDgeometry::Nplan();

  if (TMath::IsNaN(t->GetP())) return kFALSE;

  Float_t * charge = new Float_t[kNpla];
  Int_t   * nCluster = new Int_t[kNpla];
  Float_t mom  = t->GetP();
  Int_t   ipid = MCpid(t);
  TH1F   *hTmp = NULL;

  if (!SumCharge(t,charge,nCluster)) {
    delete [] charge;
    delete [] nCluster;
    return kFALSE;
  }

  Int_t index = GetIndex(mom,ipid);
  if (index > -1) {
    hTmp = (TH1F *) fHist->UncheckedAt(index);
    for (Int_t ipla = 0; ipla < kNpla; ipla++) {
      if ((charge[ipla]   >   fChargeMin) &&
          (nCluster[ipla] > fNClusterMin)){
        hTmp->Fill(charge[ipla] / ((Float_t) nCluster[ipla]));
      }
    }
  }  

  delete [] charge;
  delete [] nCluster;
  return kTRUE;

}

//_____________________________________________________________________________
Int_t AliTRDpidLQ::GetIndex(const AliTRDtrack *t)
{
  //
  // Returns the histogram index
  //

  if (TMath::IsNaN(t->GetP())) return -1;
  Float_t mom  = t->GetP();
  Int_t   ipid = MCpid(t);

  return GetIndex(mom,ipid);

}

//_____________________________________________________________________________
Int_t AliTRDpidLQ::GetIndex(Float_t mom, Int_t ipid)
{
  //
  // Returns the histogram index
  //

  if ((mom < fMinMom) || (mom > fMaxMom))  return -1;
  Int_t imom = ((Int_t) ((mom - fMinMom) / fWidMom));
  return GetIndex(imom,ipid);

}

//_____________________________________________________________________________
Int_t AliTRDpidLQ::GetIndex(Int_t imom, Int_t ipid)
{
  //
  // Returns the histogram index
  //

  if ((ipid < 0) || (ipid >= kNpid)) return -1;
  return imom * kNpid + ipid; 

}
