#ifndef AliTPCLookUpTable3DInterpolatorDFull_H
#define AliTPCLookUpTable3DInterpolatorDFull_H


/// \class AliTPCCorrectionLookUpTableFunctor
/// \brief Wrap up look-up table for correction/distortion integral or derivative (electric field)
///
/// \author Rifki Sadikin <rifki.sadikin@cern.ch>, Indonesian Institute of Sciences
/// \date Mar 4, 2015

#include <TMath.h>
#include "TMatrixD.h"
#include "AliTPC3DCylindricalInterpolatorFull.h"

class AliTPCLookUpTable3DInterpolatorDFull {
public:

	void SetNR(Int_t nr) {fNR = nr;}
	void SetNPhi(Int_t nphi) {fNPhi = nphi;}
	void SetNZ(Int_t nz) {fNZ = nz;}
	
	Int_t GetNR() {return fNR;}
	Int_t GetNPhi() {return fNPhi;}
	Int_t GetNZ() {return fNZ;}
	
	
	void SetRList(TMatrixD **rlist) {fRlist = rlist; }
	void SetPhiList(TMatrixD **philist) {fPhilist = philist; }
	void SetZList(TMatrixD **zlist) {fZlist = zlist; }
	
	void SetLookUpR(TMatrixD **lookupr) {fLookUpR = lookupr;}
	void SetLookUpPhi(TMatrixD **lookupphi) {fLookUpPhi = lookupphi;}
	void SetLookUpZ(TMatrixD **lookupz) {fLookUpZ = lookupz;}
	
	
  AliTPCLookUpTable3DInterpolatorDFull();
  AliTPCLookUpTable3DInterpolatorDFull
    (
     Int_t nr,
     Double_t rmin,Double_t rmax, Int_t nphi, Double_t phimin, Double_t phimax, Int_t nz , Double_t zmin, Double_t zmax
    );
  
  AliTPCLookUpTable3DInterpolatorDFull
  (
		Int_t nr, 
		TMatrixD**lookupr, 
		TMatrixD**r, 
		Double_t *rn,
		Int_t nphi, 
		TMatrixD**lookupphi, 
		TMatrixD**philist,
		Double_t *phin,
		Int_t nz, 
		TMatrixD**lookupz, 
		TMatrixD**zlist,
		Double_t * zlistn,
		Int_t order,
		Int_t stepR,
		Int_t stepZ,
		Int_t stepPhi,
		Int_t type
  );
  
  virtual ~AliTPCLookUpTable3DInterpolatorDFull();
  
  void GetValue
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Double_t &vr,
		Double_t &vphi,
		Double_t &vz,
		Int_t rindex,
		Int_t phiindex,
		Int_t zindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ
	);
	
  void GetValueCartesian
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Double_t &vr,
		Double_t &vphi,
		Double_t &vz,
		Int_t rindex,
		Int_t phiindex,
		Int_t zindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ
	);


  void GetValue
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Double_t &vr,
		Double_t &vphi,
		Double_t &vz,
		Int_t rindex,
		Int_t phiindex,
		Int_t zindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Int_t minzindex
	);
	
	
  void GetValueCartesian
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Double_t &vr,
		Double_t &vphi,
		Double_t &vz,
		Int_t rindex,
		Int_t phiindex,
		Int_t zindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ,
		Int_t minzindex
	);
	
  void GetValue
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Float_t &vr,
		Float_t &vphi,
		Float_t &vz,
		Int_t rindex,
		Int_t phiindex,
		Int_t zindex,
		Int_t stepR,
		Int_t stepPhi,
		Int_t stepZ
	);



        void GetValue
	(
		Double_t r, 
		Double_t phi, 
		Double_t z,
		Double_t &vr,
		Double_t &vphi,
		Double_t &vz
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
	void CopyValsCartesian();
	void CopyVals(Int_t j);


	Int_t GetIrregularGridSize() {return fInterpolatorR->GetIrregularGridSize();}
	void SetIrregularGridSize(Int_t size) 
	{	
		fInterpolatorR->SetIrregularGridSize(size);
		fInterpolatorPhi->SetIrregularGridSize(size);
		fInterpolatorZ->SetIrregularGridSize(size);
	}
	
	void SetKernelType(Int_t kernelType) 
	{
  		fInterpolatorR->SetKernelType(kernelType);
  		fInterpolatorPhi->SetKernelType(kernelType);
  		fInterpolatorZ->SetKernelType(kernelType);

	}
	
	Int_t GetKernelType() {return fInterpolatorR->GetKernelType();}	
private:


	
	Int_t fOrder;
	Int_t fIrregularGridSize;

	Int_t fNR;
	Int_t fNPhi;
	Int_t fNZ;
	
	
	TMatrixD **fLookUpR;   ///< Array to store distortion followind the drift
  TMatrixD **fLookUpPhi;   ///< Array to store distortion followind the drift
  TMatrixD **fLookUpZ;   ///< Array to store distortion followind the drift
  
  AliTPC3DCylindricalInterpolatorFull *fInterpolatorR;
  AliTPC3DCylindricalInterpolatorFull *fInterpolatorPhi;
  AliTPC3DCylindricalInterpolatorFull *fInterpolatorZ;

	
	TMatrixD **fRlist;
	TMatrixD **fPhilist;
	TMatrixD **fZlist;
	
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
	ClassDef(AliTPCLookUpTable3DInterpolatorDFull,1);
/// \endcond
};


#endif
