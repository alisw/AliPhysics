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
//____________________________________________________________________
//
// Utility class for analysing ESD data. 
// This class does sharing and background correction 
// It can form a base class for other things too.
// This line is here to make the silly code checker happy!
// This line is here to make the silly code checker happy!
//
#include "AliFMDAnaRing.h"
#include <TH2.h>
#include <TMath.h>
#include <TBrowser.h>
// #include <AliLog.h>
#define ne    80
#define emin  -3.7
#define emax   5.3
#define nm    100
#define mmin  -.5
#define mmax  9.5
// namespace {
//   Int_t   ne   = 80;
//   Float_t emin = -3.7;
//   Float_t emax =  5.3;
//   Int_t   nm   = 100;
//   Int_t   mmin = -.5;
//   Int_t   mmax = 9.5;
// }

//____________________________________________________________________
AliFMDAnaRing::AliFMDAnaRing()
  : fDet(0), 
    fRing('\0'), 
    fBg(0), 
    fCut0(0), 
    fCut1(0), 
    fUseBgCor(kFALSE), 
    fNSeq(0),
    fBareMult(),
    fMergedMult(),
    fRemovedMult(),
    fMult(),
    fStep1(),
    fStep2(),
    fStep3(),
    fNEvents(0)
{
  // Default constructor.  
  // Do not use.  
  // Used by ROOT I/O
}

//____________________________________________________________________
AliFMDAnaRing::AliFMDAnaRing(UShort_t det, Char_t ring, 
			     TH2* bg, Float_t c0, Float_t c1) 
  : fDet(det), 
    fRing(ring), 
    fBg(bg), 
    fCut0(c0), 
    fCut1(c1), 
    fUseBgCor(bg ? kTRUE : kFALSE), 
    fNSeq(ring == 'I' || ring == 'i' ? 20 : 40), 
    fBareMult(Form("bareMult%d%c", fDet, fRing), "Bare multiplicity",
	      ne, emin, emax, fNSeq, 0, 2*TMath::Pi()),
    fMergedMult(Form("mergedMult%d%c", fDet, fRing), "Merged multiplicity",
		ne, emin, emax, fNSeq, 0, 2*TMath::Pi()),
    fRemovedMult(Form("removedMult%d%c", fDet, fRing), "Removed multiplicity",
		 ne, emin, emax, fNSeq, 0, 2*TMath::Pi()),
    fMult(Form("mult%d%c", fDet, fRing), "Multiplicity", 
	  ne, emin, emax, fNSeq, 0, 2*TMath::Pi()), 
    fStep1(Form("step1_%d%c",fDet,fRing), "Step 1", nm,mmin,mmax,nm,mmin,mmax),
    fStep2(Form("step2_%d%c",fDet,fRing), "Step 2", nm,mmin,mmax,nm,mmin,mmax),
    fStep3(Form("step3_%d%c",fDet,fRing), "Step 3", nm,mmin,mmax,nm,mmin,mmax),
    fNEvents(0)
{
  // Constructor 
  // 
  // Parameters
  //     det   Detector 
  //     ring  Ring 
  //     bg    Background 
  //     c0    Lower cut 
  //     c1    higher cut */
  fName[0] = 'F'; fName[1] = 'M'; fName[2] = 'D'; 
  fName[3] = Char_t(det+48);      fName[4] = ring; 
  fName[5] = '\0';
  fBareMult.SetDirectory(0);
  fBareMult.SetXTitle("#eta");
  fBareMult.SetYTitle("#varphi");
  fBareMult.SetZTitle("d^2M_{ch}'/d#etad#varphi");
  fMergedMult.SetDirectory(0);
  fMergedMult.SetXTitle("#eta");
  fMergedMult.SetYTitle("#varphi");
  fMergedMult.SetZTitle("d^2M_{ch}''/d#etad#varphi");
  fRemovedMult.SetDirectory(0);
  fRemovedMult.SetXTitle("#eta");
  fRemovedMult.SetYTitle("#varphi");
  fRemovedMult.SetZTitle("d^2M_{ch}'''/d#etad#varphi");
  fMult.SetDirectory(0);
  fMult.SetXTitle("#eta");
  fMult.SetYTitle("#varphi");
  fMult.SetZTitle("d^2M_{ch}/d#etad#varphi");
  fStep1.SetDirectory(0);
  fStep1.SetXTitle("Bare multiplicity");
  fStep1.SetYTitle("Merged multiplicity");
  fStep2.SetDirectory(0);
  fStep2.SetXTitle("Merged multiplicity");
  fStep2.SetYTitle("Cut multiplicity");
  fStep3.SetDirectory(0);
  fStep3.SetXTitle("Cut multiplicity");
  fStep3.SetYTitle("Multiplicity");
}

//____________________________________________________________________
AliFMDAnaRing::AliFMDAnaRing(const AliFMDAnaRing& o)
  : TObject(o),
    fDet(o.fDet), 
    fRing(o.fRing), 
    fBg(o.fBg), 
    fCut0(o.fCut0), 
    fCut1(o.fCut1), 
    fUseBgCor(o.fUseBgCor), 
    fNSeq(o.fNSeq), 
    fBareMult(o.fBareMult),
    fMergedMult(o.fMergedMult),
    fRemovedMult(o.fRemovedMult),
    fMult(o.fMult), 
    fStep1(o.fStep1),
    fStep2(o.fStep2),
    fStep3(o.fStep3),
    fNEvents(o.fNEvents)
{
  // Copy constructor.  
  // Do not use.  
  // Used by ROOT I/O
  fBareMult.SetDirectory(0);
  fBareMult.SetXTitle("#eta");
  fBareMult.SetYTitle("#varphi");
  fBareMult.SetZTitle("dM_{ch}'/d#eta");
  fMergedMult.SetDirectory(0);
  fMergedMult.SetXTitle("#eta");
  fMergedMult.SetYTitle("#varphi");
  fMergedMult.SetZTitle("dM_{ch}''/d#eta");
  fRemovedMult.SetDirectory(0);
  fRemovedMult.SetXTitle("#eta");
  fRemovedMult.SetYTitle("#varphi");
  fRemovedMult.SetZTitle("dM_{ch}'''/d#eta");
  fMult.SetDirectory(0);
  fMult.SetXTitle("#eta");
  fMult.SetYTitle("#varphi");
  fMult.SetZTitle("dM_{ch}/d#eta");
}

//____________________________________________________________________
AliFMDAnaRing&
AliFMDAnaRing::operator=(const AliFMDAnaRing& o)
{
  // Assignment operator 
  // 
  // Parameters: 
  //   o   Object to assign from 
  // Returns: 
  //   Reference to this object. 
  this->fDet      = o.fDet;
  this->fRing     = o.fRing;
  this->fBg       = o.fBg;
  this->fCut0     = o.fCut0; 
  this->fCut1     = o.fCut1;
  this->fUseBgCor = o.fUseBgCor;
  this->fNSeq     = o.fNSeq;
  fBareMult.Reset();	fBareMult.Add(&o.fBareMult);
  fMergedMult.Reset();	fMergedMult.Add(&o.fMergedMult);
  fRemovedMult.Reset();	fRemovedMult.Add(&o.fRemovedMult);
  fMult.Reset();	fMult.Add(&o.fMult);
  return *this;
}


//____________________________________________________________________
Bool_t 
AliFMDAnaRing::ProcessESD(Float_t phi, Float_t eta, Float_t& m1, Float_t m2) 
{ 
  // Process ESD 
  // Parameters
  //     phi Azimuthal angle @f$ \varphi @f$ of the hit 
  //     eta Psuedo-rapidity @f$ \eta@f$ of hit
  //     m1  Multiplicity of this strip. C ontains the corrected 
  //         value on return 
  //     m2  Multiplicity of neighbor strip 
  // Return 
  //     true if hits are merged */

  // Merge shared hits 
  Bool_t  merged = false;
  Float_t sig    = m1;
  fBareMult.Fill(eta, phi, m1);
  if ((m1 < fCut0 || m2 < fCut0) && m2 > 0.001) { 
    sig = m1 + m2;
    merged = true;
    fRemovedMult.Fill(eta, phi, TMath::Min(m1,m2));
  }
  fMergedMult.Fill(eta, phi, sig);
  fStep1.Fill(m1, sig);

  // Real multiplicity 
  Double_t cmult = 1;
  if (sig < fCut0 / 2)  cmult = 0; // return merged; // cmult = 0;
  if (sig > fCut1)      cmult = 2;
  fStep2.Fill(sig, cmult);

  // Background correction
  Float_t bmult = cmult;
  if (fUseBgCor && fBg) { 
    // Float_t bgPhi = phi + (phi >  TMath::Pi() ? -2*TMath::Pi() : 0);
    Float_t bgPhi = phi - TMath::Pi();
    Int_t   bgBin = fBg->FindBin(eta, bgPhi);
    Float_t bgCor = (bgBin > 0 ? fBg->GetBinContent(bgBin) : 0);
    bmult         = (bgCor < 0.01 ? 0 : cmult / bgCor);
  }
  fStep3.Fill(cmult, bmult);
  fMult.Fill(eta, phi, bmult);

  // Call user code 
  Fill(phi, eta, bmult);
  m1 = bmult;

  return merged;
}
//____________________________________________________________________
void
AliFMDAnaRing::Finish()
{
  // Finish this task 
  // 
  // Parameters: 
  //    none
  // Returns: 
  //    nothing
  Double_t de    = (fMult.GetXaxis()->GetXmax() - 
		    fMult.GetXaxis()->GetXmin())/fMult.GetNbinsX(); 
  Double_t dp    = 2*TMath::Pi() / fNSeq;
  Double_t scale = ((fNEvents == 0 ? 1. : 1. / fNEvents) * 1 / de * 1 / dp);
  fBareMult.Scale(scale);
  fMergedMult.Scale(scale);
  fRemovedMult.Scale(scale);
  fMult.Scale(scale);
}
//____________________________________________________________________
Int_t  
AliFMDAnaRing::Color() const 
{ 
  // Get the ring specific color 
  // 
  // Returns:
  //   The color of the current ring
  switch (fDet) { 
  case 1: return kGreen + 2; 
  case 2: return kRed  + (fRing == 'I' || fRing == 'i' ? 2 : -7);
  case 3: return kBlue + (fRing == 'I' || fRing == 'i' ? 2 : -7);
  }
  return 0;
}
//____________________________________________________________________
void
AliFMDAnaRing::Browse(TBrowser* b) 
{
  // Browse this object
  // 
  // Parameters: 
  //    b  Browser to use
  // Returns: 
  //    nothing
  b->Add(&fBareMult);
  b->Add(&fMergedMult);
  b->Add(&fRemovedMult);
  b->Add(&fMult);
  b->Add(&fStep1);
  b->Add(&fStep2);
  b->Add(&fStep3);
  if (fBg) b->Add(fBg);
}

//____________________________________________________________________
//
// EOF
//

