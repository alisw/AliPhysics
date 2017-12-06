#ifndef AliTPCLookUpTable3DInterpolatorD_H
#define AliTPCLookUpTable3DInterpolatorD_H


/// \class AliTPCCorrectionLookUpTableFunctor
/// \brief Wrap up look-up table for correction/distortion integral or derivative (electric field)
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date Mar 4, 2015

#include <TMath.h>
#include "TMatrixD.h"
#include "AliTPC3DCylindricalInterpolator.h"

class AliTPCLookUpTable3DInterpolatorD {
public:

	void SetNR(Int_t nr) {fNR = nr;}
	void SetNPhi(Int_t nphi) {fNPhi = nphi;}
	void SetNZ(Int_t nz) {fNZ = nz;}
	
	Int_t GetNR() {return fNR;}
	Int_t GetNPhi() {return fNPhi;}
	Int_t GetNZ() {return fNZ;}
	
	
	void SetRList(Double_t *rlist) {fRlist = rlist; }
	void SetPhiList(Double_t *philist) {fPhilist = philist; }
	void SetZList(Double_t *zlist) {fZlist = zlist; }
	
	void SetLookUpR(TMatrixD **lookupr) {fLookUpR = lookupr;}
	void SetLookUpPhi(TMatrixD **lookupphi) {fLookUpPhi = lookupphi;}
	void SetLookUpZ(TMatrixD **lookupz) {fLookUpZ = lookupz;}
	
	
  AliTPCLookUpTable3DInterpolatorD();
  AliTPCLookUpTable3DInterpolatorD(Int_t nr, Double_t rmin, Double_t rmax, Int_t nphi, Double_t phimin, Double_t phimax, Int_t nz , Double_t zmin, Double_t zmax );
  AliTPCLookUpTable3DInterpolatorD
  (
		Int_t nr, 
		TMatrixD**lookupr, 
		Double_t *r, 
		Int_t nphi, 
		TMatrixD**lookupphi, 
		Double_t *philist,
		Int_t nz, 
		TMatrixD**lookupz, 
		Double_t *zlist,
		Int_t order
	);
  virtual ~AliTPCLookUpTable3DInterpolatorD();
  
  void GetValue
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Double_t &vr,
		Double_t &vphi,
		Double_t &vz
	);
	
	
  void GetValue
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Float_t &vr,
		Float_t &vphi,
		Float_t &vz
	);
	
	void SetOrder(Int_t order) {
		fOrder = order;
	}
  
  Double_t SplineInt3
	(	
		const Double_t xArray[], 
		const Double_t yArray[],
		Double_t y2Array[],	
		const Int_t n,
		Double_t x
	);
	
	void Spline3
	(	
		const Double_t xArray[], 
		const Double_t yArray[],
		const Int_t n,
		Double_t y2Array[]
	);
	
	 Double_t SplineInt3
	(	
		const Double_t xArray[], 
		const Float_t yArray[],
		Double_t y2Array[],	
		const Int_t n,
		Double_t x
	);
	
	void Spline3
	(	
		const Double_t xArray[], 
		const Float_t yArray[],
		const Int_t n,
		Double_t y2Array[]
	);
	
	void CopyVals();
	
private:


	
	Int_t fOrder;

	Int_t fNR;
	Int_t fNPhi;
	Int_t fNZ;
	
	
	TMatrixD **fLookUpR;   ///< Array to store distortion followind the drift
  TMatrixD **fLookUpPhi;   ///< Array to store distortion followind the drift
  TMatrixD **fLookUpZ;   ///< Array to store distortion followind the drift
  
  AliTPC3DCylindricalInterpolator *fInterpolatorR;
  AliTPC3DCylindricalInterpolator *fInterpolatorPhi;
  AliTPC3DCylindricalInterpolator *fInterpolatorZ;

	
	Double_t *fRlist;
	Double_t *fPhilist;
	Double_t *fZlist;
	
	Bool_t fIsAllocatingLookUp;
	
	
	Double_t InterpolatePhi
	( 
		const Double_t xArray[], 
		const Int_t ilow,
		const Int_t nx,
		const Float_t yArray[],
		Int_t order, 
		Double_t x 
	);
	
	
	Double_t Interpolate3DTableCyl
	( 
		Int_t order, 
		Double_t r,   
		Double_t z,
		Double_t phi,   	
		Int_t  nr,    
		Int_t  nz,
		Int_t  nphi,    	
		const Double_t rlist[], 
		const Double_t zlist[],
		const Double_t philist[], 
	
		TMatrixD **arrayofArrays 
	);
	
	void Search
	( 
		Int_t n, 
		const Double_t xArray[], 
		Double_t x, 
		Int_t &low 
	);
	
	
	Float_t Interpolate
	( 
		const Double_t xArray[], 
		const Float_t yArray[],
		Int_t order, 
		Double_t x 
	);
	
	Double_t Interpolate
	( 
		const Double_t xArray[], 
		const Double_t yArray[],
		Int_t order, 
		Double_t x 
	);
	
	
	
	
/// \cond CLASSIMP
	ClassDef(AliTPCLookUpTable3DInterpolatorD,1);
/// \endcond
};


#endif
