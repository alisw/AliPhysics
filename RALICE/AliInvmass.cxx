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

/*
$Log$
Revision 1.2  1999/09/29 09:24:28  fca
Introduction of the Copyright and cvs Log

*/

////////////////////////////////////////////////////////////////////////////////
// Class AliInvmass
// Construction of invariant mass and combinatorial background.
//
// Example :
// ---------
//
// TObjArray* photons=new TObjArray(); // Array with photon tracks for pi0 rec.
//
// // Code to create some photon tracks from pi0 decays
// Int_t ntracks=200;
// for (Int_t i=0; i<ntracks; i++)
// {
//  photons->Add(new Alitrack);
//  ...
//  ...
//  ...
// }  
//
// // Perform the invariant mass and comb. bkg. reconstruction
//
// TObjArray* allm=q.Invmass(photons,photons); // All reconstructed invariant masses
//
// TH1F* hall=new TH1F("hall","hall",200,0,2); // Histo with M_inv of all combinations
//
// Int_t nall=0;
// if (allm) nall=allm->GetEntries();
//
// AliTrack* t;
// Float_t minv;
// for (Int_t j=0; j<nall; j++)
// {
//  t=(AliTrack*)allm->At(j);
//  if (t)
//  {
//   minv=t->GetMass();
//   hall->Fill(minv);
//  }
// }
//
// TObjArray* bkgm=q.CombBkg(photons,photons); // Reconstructed comb. background
//
// TH1F* hbkg=new TH1F("hbkg","hbkg",200,0,2); // Histo with M_inv. of comb. background
//
// Int_t nbkg=0;
// if (bkgm) nbkg=bkgm->GetEntries();
//
// for (Int_t j=0; j<nbkg; j++)
// {
//  t=(AliTrack*)bkgm->At(j);
//  if (t)
//  {
//   minv=t->GetMass();
//   hbkg->Fill(minv);
//  }
// }
//
// TH1F* hsig=new TH1F("sig","sig",200,0,2);   // Histo with the bkg. subtracted signal
// hsig->Sumw2();
// hsig->Add(hall,hbkg,1,-1);
//
//
// Note : By default the storage of the reconstructed information is performed
//        in separate TObjArrays for the signal and comb. background resp.
//        In order to limit the memory usage, AliInvmass::SetStorageMode(1) may be
//        used to activate only a single TObjArray to store the reconstructed information.
//        Consequently, the following statements 
//
//         TObjArray* allm=q.Invmass(photons,photons);
//         TObjArray* bkgm=q.CombBkg(photons,photons);
//
//        will result in the fact that after he invokation of CombBkg
//        the information of "allm" is lost due to the fact that the storage is
//        is re-used for "bkgm" in case the "single storage" option has been selected.
//        Usage of the, in that case invalid, pointer "allm" may cause your
//        program to crash.
//
//        * Thus : In case of single storage usage, all invokations of the returned
//                 array pointer have to be completed before invoking any memberfunction
//                 of the same AliInvmass object again.
//        
//        
//
//--- Author: Nick van Eijndhoven 12-apr-1999 UU-SAP Utrecht
////////////////////////////////////////////////////////////////////////////////

#include "AliInvmass.h"
 
ClassImp(AliInvmass) // Class implementation to enable ROOT I/O
 
AliInvmass::AliInvmass()
{
// Creation of an AliInvmass object and initialisation of parameters
 fPi=acos(-1.);
 fMode=2;
 fBkg=0;
 fNewtheta=1;
 fNewphi=1;
 fMinv=0;
 fMbkg=0;
}
////////////////////////////////////////////////////////////////////////////////
AliInvmass::~AliInvmass()
{
// Destructor to delete dynamically allocated memory
 if (fMinv)
 {
  fMinv->Delete();
  delete fMinv;
  fMinv=0;
 }

 if (fMbkg)
 {
  fMbkg->Delete();
  delete fMbkg;
  fMbkg=0;
 }
}
////////////////////////////////////////////////////////////////////////////////
void AliInvmass::SetStorageMode(Int_t m)
{
// Set storage mode for the result arrays for inv. mass and comb. background
 fMode=2;
 if (m==1) fMode=1;
}
////////////////////////////////////////////////////////////////////////////////
void AliInvmass::SetThetaSwitch(Int_t i)
{
// Enable/Disable (1/0) switching of theta angle in comb. bkg. reconstruction.
// Default : Switching of theta is enabled.
 fNewtheta=1;
 if (i==0) fNewtheta=0;
}
////////////////////////////////////////////////////////////////////////////////
void AliInvmass::SetPhiSwitch(Int_t i)
{
// Enable/Disable (1/0) switching of phi angle in comb. bkg. reconstruction.
// Default : Switching of phi is enabled.
 fNewphi=1;
 if (i==0) fNewphi=0;
}
////////////////////////////////////////////////////////////////////////////////
Int_t AliInvmass::GetStorageMode()
{
// Provide mode of storage for the result arrays for inv. mass and comb. background
 return fMode;
}
////////////////////////////////////////////////////////////////////////////////
Int_t AliInvmass::GetThetaSwitch()
{
// Provide the theta switching flag
 return fNewtheta;
}
////////////////////////////////////////////////////////////////////////////////
Int_t AliInvmass::GetPhiSwitch()
{
// Provide the phi switching flag
 return fNewphi;
}
////////////////////////////////////////////////////////////////////////////////
TObjArray* AliInvmass::Invmass(TObjArray* a1,TObjArray* a2)
{
// Perform two-particle invariant mass reconstruction
 fBkg=0;
 Combine(a1,a2);
 return fMinv;
}
////////////////////////////////////////////////////////////////////////////////
TObjArray* AliInvmass::CombBkg(TObjArray* a1,TObjArray* a2)
{
// Perform two-particle combinatorial background reconstruction
 fBkg=1;
 Combine(a1,a2);
 if (fMode != 1)
 {
  return fMbkg;
 }
 else
 {
  return fMinv;
 }
}
////////////////////////////////////////////////////////////////////////////////
void AliInvmass::Combine(TObjArray* a1,TObjArray* a2)
{
// Perform two-particle invariant mass reconstruction
 
 if ((!fBkg || fMode==1) && fMinv)
 {
  fMinv->Delete();
  delete fMinv;
  fMinv=0;
 }
 
 if (fBkg && (fMode !=1) && fMbkg)
 {
  fMbkg->Delete();
  delete fMbkg;
  fMbkg=0;
 }

 Int_t isame; // Indicates whether both lists are identical
 isame=0;
 if (a1==a2) isame=1;

 // Index i must loop over the shortest of a1 and a2
 TObjArray* listi=a1;
 TObjArray* listj=a2;
 Int_t ni=a1->GetEntries();
 Int_t nj=a2->GetEntries();
 if (nj < ni)
 {
  ni=a2->GetEntries();
  nj=a1->GetEntries();
  listi=a2;
  listj=a1;
 }

 AliTrack* p1=0;  
 AliTrack* p2=0;
 AliTrack* px=0;
 Ali4Vector ptot;
 AliTrack* t=0;
 Double_t v2[4],vx[4];
 Float_t q1,q2;
 
 Int_t jmin; // Start index for list j
 Int_t jx;   // Index for randomly picked particle for comb. bkg. reconstruction 
 
 for (Int_t i=0; i<ni; i++) // Select first a particle from list i
 {
  p1=(AliTrack*)listi->At(i);
  p2=0;

  if (!p1) continue;

  jmin=0;
  if (isame) jmin=i+1;
  for (Int_t j=jmin; j<nj; j++) // Select also a particle from list j
  {
   p2=(AliTrack*)listj->At(j);
   if (p1==p2) p2=0; // Don't combine particle with itself

   if (!p2) continue;

   p2->GetVector(v2,"sph");

   // Take theta and phi from randomly chosen other list j particle for bkg. reconstr.
   if (fBkg)
   {
    px=0; 
    if ((!isame && nj>1) || (isame && nj>2))
    {
     jx=int(fRndm.Uniform(0,float(nj)));
     px=(AliTrack*)listj->At(jx);
     
     while (!px || px==p2 || px==p1)
     {
      jx++;
      if (jx >= nj) jx=0;
      px=(AliTrack*)listj->At(jx);
     }

     px->GetVector(vx,"sph");
     if (fNewtheta) v2[2]=vx[2]; // Replace the theta angle in the v2 vector
     if (fNewphi) v2[3]=vx[3]; // Replace the phi angle in the v2 vector
    }
   }
 
   if ((!fBkg && p2) || (fBkg && px))
   {
    // Store the data of this two-particle combination
    ptot.SetVector(v2,"sph");
    ptot=(Ali4Vector)(ptot+(*p1));
    q1=p1->GetCharge();
    q2=p2->GetCharge();
    t=new AliTrack;
    t->Set4Momentum(ptot);
    t->SetCharge(q1+q2);
    if (!fBkg || fMode==1) 
    {
     if (!fMinv) fMinv=new TObjArray();
     fMinv->Add(t);
    }
    else
    {
     if (!fMbkg) fMbkg=new TObjArray();
     fMbkg->Add(t);
    }
   }

  } // End of second particle loop

 } // End of first particle loop
}
////////////////////////////////////////////////////////////////////////////////
