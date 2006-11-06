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

//-------------------------------------------------------
//          Implementation of the TPC clusterer
//
//   Origin: Jouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include "AliTPCclusterer.h"
#include "AliTPCcluster.h"
#include <TObjArray.h>
#include <TFile.h>
#include "AliDigits.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliTPCClustersRow.h"
#include <TTree.h>

ClassImp(AliTPCclusterer)


//____________________________________________________
AliTPCclusterer::AliTPCclusterer(const AliTPCclusterer &param)
                :TObject(),
                 fPar(0)
{
  //
  // dummy
  //
  fPar = param.fPar;
}
//-----------------------------------------------------
AliTPCclusterer & AliTPCclusterer::operator =(const AliTPCclusterer & param)
{
  //
  // assignment operator - dummy
  //
    fPar = param.fPar;
  return (*this);
}
//____________________________________________________
void AliTPCclusterer::FindPeaks(Int_t k,Int_t max,
AliBin *b,Int_t *idx,UInt_t *msk,Int_t& n) {
  //find local maxima
  if (n<31)
  if (IsMaximum(k,max,b)) {
    idx[n]=k; msk[n]=(2<<n);
    n++;
  }
  b[k].SetMask(0);
  if (b[k-max].GetMask()&1) FindPeaks(k-max,max,b,idx,msk,n);
  if (b[k-1  ].GetMask()&1) FindPeaks(k-1  ,max,b,idx,msk,n);
  if (b[k+max].GetMask()&1) FindPeaks(k+max,max,b,idx,msk,n);
  if (b[k+1  ].GetMask()&1) FindPeaks(k+1  ,max,b,idx,msk,n);
}

void AliTPCclusterer::MarkPeak(Int_t k, Int_t max, AliBin *bins, UInt_t m) {
  //mark this peak
  UShort_t q=bins[k].GetQ();

  bins[k].SetMask(bins[k].GetMask()|m); 

  if (bins[k-max].GetQ() <= q)
     if ((bins[k-max].GetMask()&m) == 0) MarkPeak(k-max,max,bins,m);
  if (bins[k-1  ].GetQ() <= q)
     if ((bins[k-1  ].GetMask()&m) == 0) MarkPeak(k-1  ,max,bins,m);
  if (bins[k+max].GetQ() <= q)
     if ((bins[k+max].GetMask()&m) == 0) MarkPeak(k+max,max,bins,m);
  if (bins[k+1  ].GetQ() <= q)
     if ((bins[k+1  ].GetMask()&m) == 0) MarkPeak(k+1  ,max,bins,m);
}

void AliTPCclusterer::MakeCluster(Int_t k,Int_t max,AliBin *bins,UInt_t m,
AliTPCcluster &c) {
  //make cluster using digits of this peak
  Float_t q=(Float_t)bins[k].GetQ();
  Int_t i=k/max, j=k-i*max;

  c.SetQ(c.GetQ()+q);
  c.SetY(c.GetY()+i*q); 
  c.SetZ(c.GetZ()+j*q); 
  c.SetSigmaY2(c.GetSigmaY2()+i*i*q);
  c.SetSigmaZ2(c.GetSigmaZ2()+j*j*q);

  bins[k].SetMask(0xFFFFFFFE);
  
  if (bins[k-max].GetMask() == m) MakeCluster(k-max,max,bins,m,c);
  if (bins[k-1  ].GetMask() == m) MakeCluster(k-1  ,max,bins,m,c);
  if (bins[k+max].GetMask() == m) MakeCluster(k+max,max,bins,m,c);
  if (bins[k+1  ].GetMask() == m) MakeCluster(k+1  ,max,bins,m,c);
}

//_____________________________________________________________________________
Int_t AliTPCclusterer::Digits2Clusters(TTree *dTree, TTree *cTree) {
  //-----------------------------------------------------------------
  // This is a simple cluster finder.
  //-----------------------------------------------------------------
  TBranch *branch=dTree->GetBranch("Segment");
  if (!branch) {
     Error("Digits2Cluster","Can't get the branch !");
     return 1;
  }
  AliSimDigits digarr, *dummy=&digarr;
  branch->SetAddress(&dummy);
  
  AliTPCClustersRow ddd,*clrow=&ddd;
  clrow->SetClass("AliTPCcluster"); clrow->SetArray(1);
  cTree->Branch("Segment","AliTPCClustersRow",&clrow,32000,200);    

  const Int_t kMAXZ=fPar->GetMaxTBin()+2;

  Int_t nclusters=0;

  Int_t nentries = (Int_t)dTree->GetEntries();
  for (Int_t n=0; n<nentries; n++) {
   
    Int_t sec, row;
    dTree->GetEvent(n);

    if (!fPar->AdjustSectorRow(digarr.GetID(),sec,row)) {
       Error("Digits2Clusters","!nvalid segment ID ! %d",digarr.GetID());
       continue;
    }

    clrow=new AliTPCClustersRow();

    clrow->SetClass("AliTPCcluster"); clrow->SetArray(1);
    clrow->SetID(digarr.GetID());

    cTree->GetBranch("Segment")->SetAddress(&clrow);

    Float_t rx=fPar->GetPadRowRadii(sec,row);

    Int_t npads, sign;
    {
       const Int_t kNIS=fPar->GetNInnerSector(), kNOS=fPar->GetNOuterSector();
       if (sec < kNIS) {
          npads = fPar->GetNPadsLow(row);
          sign = (sec < kNIS/2) ? 1 : -1;
       } else {
          npads = fPar->GetNPadsUp(row);
          sign = ((sec-kNIS) < kNOS/2) ? 1 : -1;
       }
    }

    const Int_t kMAXBIN=kMAXZ*(npads+2);
    AliBin *bins=new AliBin[kMAXBIN];
    for (Int_t ii=0;ii<kMAXBIN;ii++) {
       bins[ii].SetQ(0); bins[ii].SetMask(0xFFFFFFFE);
    }

    digarr.First();
    do {
       Short_t dig=digarr.CurrentDigit();
       if (dig<=fPar->GetZeroSup()) continue;
       Int_t j=digarr.CurrentRow()+1, i=digarr.CurrentColumn()+1;
       bins[i*kMAXZ+j].SetQ(dig);
       bins[i*kMAXZ+j].SetMask(1);
    } while (digarr.Next());

    Int_t ncl=0;
    for (Int_t i=0; i<kMAXBIN; i++) {
      if ((bins[i].GetMask()&1) == 0) continue;
      Int_t idx[32]; UInt_t msk[32]; Int_t npeaks=0;
      FindPeaks(i, kMAXZ, bins, idx, msk, npeaks);

      if (npeaks>30) continue;

      Int_t k,l;
      for (k=0; k<npeaks-1; k++){//mark adjacent peaks
        if (idx[k] < 0) continue; //this peak is already removed
        for (l=k+1; l<npeaks; l++) {
           if (idx[l] < 0) continue; //this peak is already removed
           Int_t ki=idx[k]/kMAXZ, kj=idx[k] - ki*kMAXZ;
           Int_t li=idx[l]/kMAXZ, lj=idx[l] - li*kMAXZ;
           Int_t di=TMath::Abs(ki - li);
           Int_t dj=TMath::Abs(kj - lj);
           if (di>1 || dj>1) continue;
           if (bins[idx[k]].GetQ() > bins[idx[l]].GetQ()) {
              msk[l]=msk[k];
              idx[l]*=-1;
           } else {
              msk[k]=msk[l];
              idx[k]*=-1;
              break;
           } 
        }
      }

      for (k=0; k<npeaks; k++) {
        MarkPeak(TMath::Abs(idx[k]), kMAXZ, bins, msk[k]);
      }
        
      for (k=0; k<npeaks; k++) {
         if (idx[k] < 0) continue; //removed peak
         AliTPCcluster c;
         MakeCluster(idx[k], kMAXZ, bins, msk[k], c);
         if (c.GetQ() < 5) continue; //noise cluster
         c.SetY(c.GetY()/c.GetQ());
         c.SetZ(c.GetZ()/c.GetQ());

         Float_t s2 = c.GetSigmaY2()/c.GetQ() - c.GetY()*c.GetY();
         Float_t w=fPar->GetPadPitchWidth(sec);
         c.SetSigmaY2((s2 + 1./12.)*w*w);
         if (s2 != 0.) {
	   c.SetSigmaY2(c.GetSigmaY2()*0.108);
	   if (sec<fPar->GetNInnerSector()) c.SetSigmaY2(c.GetSigmaY2()*2.07);
         }

         s2 = c.GetSigmaZ2()/c.GetQ() - c.GetZ()*c.GetZ();
         w=fPar->GetZWidth();
         c.SetSigmaZ2((s2 + 1./12.)*w*w);
         if (s2 != 0.) {
	   c.SetSigmaZ2(c.GetSigmaZ2()*0.169);
	   if (sec<fPar->GetNInnerSector()) c.SetSigmaZ2(c.GetSigmaZ2()*1.77);
         }

         c.SetY((c.GetY() - 0.5 - 0.5*npads)*fPar->GetPadPitchWidth(sec));
         c.SetZ(fPar->GetZWidth()*(c.GetZ()-1)); 
         c.SetZ(c.GetZ() - 3.*fPar->GetZSigma()); // PASA delay 
         c.SetZ(sign*(fPar->GetZLength() - c.GetZ()));

         if (rx<230./250.*TMath::Abs(c.GetZ())) continue;

         Int_t ki=idx[k]/kMAXZ, kj=idx[k] - ki*kMAXZ;
         c.SetLabel(digarr.GetTrackID(kj-1,ki-1,0),0);
         c.SetLabel(digarr.GetTrackID(kj-1,ki-1,1),1);
         c.SetLabel(digarr.GetTrackID(kj-1,ki-1,2),2);

         c.SetQ(bins[idx[k]].GetQ());

         if (ki==1 || ki==npads || kj==1 || kj==kMAXZ-2) {
           c.SetSigmaY2(c.GetSigmaY2()*25.);
           c.SetSigmaZ2(c.GetSigmaZ2()*4.);
         }
         clrow->InsertCluster(&c); ncl++;
      }
    }
    cTree->Fill();

    delete clrow;

    nclusters+=ncl;

    delete[] bins;  
  }

  Info("Digits2Cluster","Number of found clusters : %d",nclusters);

  return 0;
}

