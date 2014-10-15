#include <AliFMDESDRevertexer.h>
#include <AliFMDGeometry.h>
#include <AliESDFMD.h>
#include <TMath.h>
#include <AliLog.h>

ClassImp(AliFMDESDRevertexer)
#if 0 // for emacs 
;
#endif

//____________________________________________________________________
AliFMDESDRevertexer::AliFMDESDRevertexer()
{
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  geom->Init();
  geom->InitTransformations();
}

//____________________________________________________________________
Bool_t
AliFMDESDRevertexer::Revertex(AliESDFMD* fmdEsd, Double_t vz) const
{
  // Recalculate the various quantities based on updated 
  // primary vertex position. 
  // 
  // Parameters: 
  //    fmdEsd    FMD ESD object 
  //    vz        New vertex location (along the z-axis)
  //
  // Return:
  //    true on success, false if there was an error during the 
  //    recalculations.   Please inspect log output for details. 
  // 
  if (!fmdEsd) return kFALSE;
  
  Bool_t ret = kTRUE;
  for (UShort_t det = 1; det <= 3; det++) { 
    UShort_t nrng = (det == 1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nrng; ir++) {
      Char_t   rng  = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ?  20 :  40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for (UShort_t str = 0; str < nstr; str++) { 
	Double_t phi, r, theta;
	Double_t eta      = AliESDFMD::kInvalidEta;
	Double_t oldEta   = fmdEsd->Eta(det, rng, 0, str);
	// if (oldEta == AliESDFMD::kInvalidEta) continue;

	Double_t oldTheta = Eta2Theta(oldEta);
	Bool_t   ret1     = PhysicalCoordinates(det, rng, 0, str, vz, 
						eta, phi, r, theta);
	fmdEsd->SetEta(det, rng, 0, str, eta);

	if (!ret1) {
	  // If the was an error, then there's no reason to go on with
	  // this strip-ring.  Note, that the eta is correctly set to
	  // AliESDFMD::kInvalidMult. 
	  AliWarning(Form("Failed to calculate eta, phi for "
			  "FMD%d%c[%02d,%03d] with v_z=%9.4f",
			  det, rng, 0, str, vz));
	  ret = kFALSE;
	  continue;
	}

	Double_t corr = TMath::Abs(TMath::Cos(theta));
	if (fmdEsd->IsAngleCorrected()) {
	  if (oldEta != AliESDFMD::kInvalidMult)
	    corr /= TMath::Abs(TMath::Cos(oldTheta));
	  for (UShort_t sec = 0; sec < nsec; sec++) { 
	    Double_t mult = fmdEsd->Multiplicity(det, rng, sec, str);
	    if (mult == AliESDFMD::kInvalidMult) continue;
	    fmdEsd->SetMultiplicity(det, rng, sec, str, corr * mult);
	  }
	}
      }
    }
  }

  return ret;
}

//____________________________________________________________________
Double_t
AliFMDESDRevertexer::Eta2Theta(Double_t eta) const
{
  if (eta == AliESDFMD::kInvalidEta) return 0;
  return 2 * TMath::ATan(TMath::Exp(-eta));
}


//____________________________________________________________________
Bool_t
AliFMDESDRevertexer::PhysicalCoordinates(UShort_t  det, 
					 Char_t    rng, 
					 UShort_t  sec, 
					 UShort_t  str,
					 Double_t  vz,
					 Double_t& eta, 
					 Double_t& phi,
					 Double_t& r,
					 Double_t& theta) const
{
  // Get the eta and phi of a digit 
  // 
  // Get geometry. 
  AliFMDGeometry* geom = AliFMDGeometry::Instance();
  Double_t x=0, y=0, z=0;
  geom->Detector2XYZ(det, rng, sec, str, x, y, z);

  return AliFMDGeometry::XYZ2REtaPhiTheta(x, y, z-vz, r, eta, phi, theta);
}


