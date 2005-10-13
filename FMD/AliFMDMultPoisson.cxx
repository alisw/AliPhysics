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

//____________________________________________________________________
// 
// Reconstruct charged particle multiplicity in the FMD 
// 
// [See also the AliFMDReconstructor class]
// 
// This class reconstructs the muliplicity in regions based on the
// ratio of empty to full strips. 
//
#include "AliFMD.h"			// ALIFMD_H
#include "AliFMDGeometry.h"		// ALIFMDGEOMETRY_H
#include "AliFMDDetector.h"		// ALIFMDDETECTOR_H
#include "AliFMDRing.h"			// ALIFMDRING_H
#include "AliFMDMultPoisson.h"		// ALIFMDMULTPOISSON_H
#include "AliFMDMultRegion.h"		// ALIFMDMULTREGION_H
#include "AliFMDDigit.h"		// ALIFMDDIGIT_H
#include "AliLog.h"			// ALILOG_H
#include <TClonesArray.h>               // ROOT_TClonesArray
#include <TTree.h>               	// ROOT_TTree

//____________________________________________________________________
ClassImp(AliFMDMultPoisson)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDMultPoisson::AliFMDMultPoisson()
  : AliFMDMultAlgorithm("Poisson", "Poisson"),
    fDeltaEta(0), 
    fDeltaPhi(0), 
    fThreshold(0)
{
  SetDeltaEta();
  SetDeltaPhi();
  SetThreshold();
  fMult = new TClonesArray("AliFMDMultRegion", 1000);
}

//____________________________________________________________________
void
AliFMDMultPoisson::PreEvent(TTree* tree, Float_t ipZ) 
{
  // Reset internal data
  AliFMDMultAlgorithm::PreEvent(tree, ipZ);
  fCurrentVertexZ = ipZ;
  fEmpty.Reset(kFALSE);

  // Make a branch in the reconstruction tree. 
  const Int_t kBufferSize = 16000;
  fTreeR->Branch("FMDPoisson", &fMult, kBufferSize);  
  
}

//____________________________________________________________________
void
AliFMDMultPoisson::ProcessDigit(AliFMDDigit*  digit, 
				Float_t       /* eta */, 
				Float_t       /* phi */, 
				UShort_t      count)
{
  // Process one digit. 
  // 
  // Parameters: 
  //    
  //   digit		Digit to process 
  //   ipZ		Z--coordinate of the primary interaction
  //                    vertex of this event 
  //
  if (!digit) return;
  if (count < fThreshold) fEmpty(digit->Detector() - 1, 
				 digit->Ring(), 
				 digit->Sector(), 
				 digit->Strip()) = kTRUE;
}

//____________________________________________________________________
void
AliFMDMultPoisson::PostEvent() 
{
  // Fill the branch 
  // Based on the information in the cache, do the reconstruction. 

  // Loop over the detectors 
  for (Int_t i = 1; i <= 3; i++) {
    AliFMDGeometry* fmd = AliFMDGeometry::Instance();
    AliFMDDetector* sub = fmd->GetDetector(i);
    if (!sub) continue;
	
    // Loop over the rings in the detector
    for (Int_t j = 0; j < 2; j++) {
      AliFMDRing* r  = sub->GetRing((j == 0 ? 'I' : 'O'));
      Float_t     rZ = sub->GetRingZ((j == 0 ? 'I' : 'O'));
      if (!r) continue;
      
      // Calculate low/high theta and eta 
      // FIXME: Is this right? 
      Float_t realZ    = fCurrentVertexZ + rZ;
      Float_t thetaOut = TMath::ATan2(r->GetHighR(), realZ);
      Float_t thetaIn  = TMath::ATan2(r->GetLowR(), realZ);
      Float_t etaOut   = - TMath::Log(TMath::Tan(thetaOut / 2));
      Float_t etaIn    = - TMath::Log(TMath::Tan(thetaIn / 2));
      if (TMath::Abs(etaOut) > TMath::Abs(etaIn)) {
	Float_t tmp = etaIn;
	etaIn       = etaOut;
	etaOut      = tmp;
      }

      //-------------------------------------------------------------
      //
      // Here starts poisson method 
      //
      // Calculate eta step per strip, number of eta steps, number of
      // phi steps, and check the sign of the eta increment 
      Float_t stripEta = (Float_t(r->GetNStrips()) / (etaIn - etaOut));
      Int_t   nEta     = Int_t(TMath::Abs(etaIn - etaOut) / fDeltaEta); 
      Int_t   nPhi     = Int_t(360. / fDeltaPhi);
      Float_t sign     = TMath::Sign(Float_t(1.), etaIn);

      AliDebug(10, Form("FMD%d%c Eta range: %f, %f %d Phi steps",
			sub->GetId(), r->GetId(), etaOut, etaIn, nPhi));

      // Loop over relevant phi values 
      for (Int_t p = 0; p < nPhi; p++) {
	Float_t  minPhi    = p * fDeltaPhi;
	Float_t  maxPhi    = minPhi + fDeltaPhi;
	UShort_t minSector = UShort_t(minPhi / 360) * r->GetNSectors();
	UShort_t maxSector = UShort_t(maxPhi / 360) * r->GetNSectors();
	
	//	AliDebug(10, Form(" Now in phi range %f, %f (sectors %d,%d)",
	//		  minPhi, maxPhi, minSector, maxSector));
	// Loop over relevant eta values 
	for (Int_t e = nEta; e >= 0; --e) {
	  Float_t  maxEta   = etaIn  - sign * e * fDeltaEta;
	  Float_t  minEta   = maxEta - sign * fDeltaEta;
	  if (sign > 0)  minEta = TMath::Max(minEta, etaOut);
	  else           minEta = TMath::Min(minEta, etaOut);
	  Float_t  theta1   = 2 * TMath::ATan(TMath::Exp(-minEta));
	  Float_t  theta2   = 2 * TMath::ATan(TMath::Exp(-maxEta));
	  Float_t  minR     = TMath::Abs(realZ * TMath::Tan(theta2));
	  Float_t  maxR     = TMath::Abs(realZ * TMath::Tan(theta1));
	  // Calculate the weighted mean eta of the region
	  Float_t  minW2    = TMath::Power(minR * 2 * TMath::Pi() * 
					   ((maxPhi - minPhi)/360),2);
	  Float_t  maxW2    = TMath::Power(minR * 2 * TMath::Pi() * 
					   ((maxPhi - minPhi)/360), 2);
	  Float_t  meanEta  = ((minEta / minW2 + maxEta / maxW2) / 
			       (1 / (minW2 + maxW2)));
	  //UShort_t minStrip = UShort_t((etaIn - maxEta) * stripEta + 0.5);
	  // UShort_t maxStrip = UShort_t((etaIn - minEta) * stripEta + 0.5);

	  UShort_t minStrip = UShort_t(r->GetNStrips() - 
				       (etaIn - minEta) * stripEta + 0.5);
	  UShort_t maxStrip = UShort_t(r->GetNStrips() - 
				       (etaIn - maxEta) * stripEta + 0.5);

	  AliDebug(10, Form("  Now in eta range %f, %f (strips %d, %d)\n"
	  		    "    [radii %f, %f, thetas %f, %f, sign %d]", 
	    		    minEta, maxEta, minStrip, maxStrip,
			    minR, maxR, theta1, theta2, sign));
	  
	  // Count number of empty strips
	  Int_t   emptyStrips = 0;
	  for (Int_t sector = minSector; sector < maxSector; sector++) 
	    for (Int_t strip = minStrip; strip < maxStrip; strip++) 
	      if (fEmpty(sub->GetId() - 1, r->GetId(), sector, strip)) 
		emptyStrips++;
	  
	  // The total number of strips 
	  Float_t nTotal = (maxSector - minSector) * (maxStrip - minStrip);
	  // Log ratio of empty to total number of strips 
	  
	  Double_t lambda = (emptyStrips > 0 ? 
			     - TMath::Log(Double_t(emptyStrips) / nTotal) :
			     1);

	  // The reconstructed number of particles is then given by 
	  Int_t reconstructed = Int_t(lambda * nTotal + 0.5);
	  AliDebug(10, Form("Lambda= %d / %f = %f particles %d", 
			    emptyStrips, nTotal, 
			    Float_t(emptyStrips) / nTotal , reconstructed));
	    
	  // Add a AliFMDMultRegion to the reconstruction tree. 
	  AliFMDMultRegion* m = new((*fMult)[fNMult])   
	    AliFMDMultRegion(sub->GetId(), r->GetId(),
			     minSector, maxSector, minStrip, maxStrip,
			     minEta, maxEta, meanEta, minPhi, maxPhi,
			     reconstructed, AliFMDMultRegion::kPoission);
	  (void)m;
	  fNMult++;
	} // phi 
      } // eta
    } // ring 
  } // detector 
}


//____________________________________________________________________
// 
// EOF
//
