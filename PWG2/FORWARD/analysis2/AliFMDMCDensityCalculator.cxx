#include "AliFMDMCDensityCalculator.h"
#include <TMath.h>
#include "AliForwardCorrectionManager.h"
#include "AliFMDStripIndex.h"
#include "AliMCEvent.h"
// #include "AliFMDAnaParameters.h"
#include "AliLog.h"
#include <TH2D.h>

ClassImp(AliFMDMCDensityCalculator)
#if 0
; // For Emacs
#endif 


//____________________________________________________________________
AliFMDMCDensityCalculator&
AliFMDMCDensityCalculator::operator=(const AliFMDMCDensityCalculator& o)
{
  AliFMDDensityCalculator::operator=(o);
  return *this;
}

    
//____________________________________________________________________
Bool_t
AliFMDMCDensityCalculator::CalculateMC(const AliMCEvent&       event,
				       AliForwardUtil::Histos& hists,
				       Double_t                vz,
				       UShort_t                vtxbin)
{
  Int_t nTracks = event.GetNumberOfTracks();
  for (Int_t iTr = 0; iTr < nTracks; iTr++) { 
    AliMCParticle* particle = 
      static_cast<AliMCParticle*>(event.GetTrack(iTr));
    
    // Check the returned particle 
    if (!particle) continue;
    
    // Check if this charged and a primary 
    Bool_t isCharged = particle->Charge() != 0;
    if (!isCharged) continue;

    Int_t nTrRef = particle->GetNumberOfTrackReferences();
    for (Int_t iTrRef = 0; iTrRef < nTrRef; iTrRef++) { 
      AliTrackReference* ref = particle->GetTrackReference(iTrRef);
      
      // Check existence 
      if (!ref) continue;

      // Check that we hit an FMD element 
      if (ref->DetectorId() != AliTrackReference::kFMD) 
	continue;

      // Get the detector coordinates 
      UShort_t d, s, t;
      Char_t r;
      AliFMDStripIndex::Unpack(ref->UserId(), d, r, s, t);

      Double_t x    = ref->X();
      Double_t y    = ref->Y();
      Double_t z    = ref->Z()-vz;
      Double_t rr   = TMath::Sqrt(x*x+y*y);
      Double_t phi  = TMath::ATan2(y,x);
      Double_t theta= TMath::ATan2(rr,z);
      Double_t eta  = -TMath::Log(TMath::Tan(theta/2));

      Float_t  c    = Correction(d,r,s,t,vtxbin,eta,false);
      fCorrections->Fill(c);
      
      TH2D*    h    = hists.Get(d,r);
      h->Fill(eta,phi, 1 * c);
    }
  }
  return kTRUE;
}

//____________________________________________________________________
Bool_t
AliFMDMCDensityCalculator::Calculate(const AliESDFMD&,
				     AliForwardUtil::Histos&,
				     UShort_t,
				     Bool_t)
{
  AliWarning("Method Calculate disabled for this class. If you need this, "
	     "make an AliFMDDensityCalculator object instead");
  return kFALSE;
}

//____________________________________________________________________
//
// EOF
//
	  


