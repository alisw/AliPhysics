#ifndef ALIITSGEOMMATRIX_H
#define ALIITSGEOMMATRIX_H
/* Copyright(c) 2000, ALICE Experiment at CERN, All rights reserved. *
 * see cxx source for full Copyright notice.                         */
/* $Id: */
////////////////////////////////////////////////////////////////////////
// ITS geometry manipulation routines on the module level. This class is
// to replace the structure ITS_geom in the class AliITSgeom.
// Created May 30 2000.
// version 0.0.0
// By Bjorn S. Nilsen
////////////////////////////////////////////////////////////////////////

class AliITSgeomMatrix{
 public:
	AliITSgeomMatrix(); // Default constructor
	AliITSgeomMatrix(const Int_t idt,const Int_t id[3],
			 const Double_t rot[3],const Double_t tran[3]);
        AliITSgeomMatrix(const Int_t idt,const Int_t id[3],
		         const Double_t matrix[3][3],const Double_t tran[3]);
        AliITSgeomMatrix(const Double_t rotd[6]/*degrees Geant angles*/,
                         const Int_t idt,const Int_t id[3],
                         const Double_t tran[3]);
	AliITSgeomMatrix(const AliITSgeomMatrix &source);
	void operator=(const AliITSgeomMatrix &sourse); // copy
	virtual ~AliITSgeomMatrix(){};
	void print(ostream *os);
	void PrintTitles(ostream *os);
	void read(istream *is);

	void SetAngles(const Double_t rot[3]){// [radians]
              for(Int_t i=0;i<3;i++)frot[i] = rot[i];this->MatrixFromAngle();}
	void SetTranslation(const Double_t tran[3]){
	                    for(Int_t i=0;i<3;i++) ftran[i] = tran[i];}
	void SetMatrix(const Double_t matrix[3][3]){ for(Int_t i=0;i<3;i++)
	 for(Int_t j=0;j<3;j++) fm[i][j]=matrix[i][j];this->AngleFromMatrix();}
	void SetDetectorIndex(const Int_t idt) {fDetectorIndex = idt;}
	void SetIndex(const Int_t id[3]){
	                   for(Int_t i=0;i<3;i++) fid[i] = id[i];}
	void GetAngles(Double_t rot[3]){// [radians]
	                   for(Int_t i=0;i<3;i++)  rot[i] = frot[i];}
	void GetTranslation(Double_t tran[3]){
	                    for(Int_t i=0;i<3;i++) tran[i] = ftran[i];}
	void GetMatrix(Double_t matrix[3][3]){for(Int_t i=0;i<3;i++)
		         for(Int_t j=0;j<3;j++) matrix[i][j] = fm[i][j];}
	Int_t GetDetectorIndex() {return fDetectorIndex;}
	void  GetIndex(Int_t id[3]){for(Int_t i=0;i<3;i++) id[i] = fid[i];}
	void  MatrixFromSixAngles(const Double_t *ang);
	void  SixAnglesFromMatrix(Double_t *ang);

	void GtoLPosition(const Double_t g[3],Double_t l[3]);
	void LtoGPosition(const Double_t l[3],Double_t g[3]);
	void GtoLMomentum(const Double_t g[3],Double_t l[3]);
	void LtoGMomentum(const Double_t l[3],Double_t g[3]);
	void GtoLPositionError(const Double_t g[3][3],Double_t l[3][3]);
	void LtoGPositionError(const Double_t l[3][3],Double_t g[3][3]);
	// Tracking Related Routines
	void GtoLPositionTracking(const Double_t g[3],Double_t l[3]);
	void LtoGPositionTracking(const Double_t l[3],Double_t g[3]);
	void GtoLMomentumTracking(const Double_t g[3],Double_t l[3]);
	void LtoGMomentumTracking(const Double_t l[3],Double_t g[3]);
	void GtoLPositionErrorTracking(const Double_t g[3][3],
				       Double_t l[3][3]);
	void LtoGPositionErrorTracking(const Double_t l[3][3],
				       Double_t g[3][3]);
	Double_t Distance2(const Double_t t[3]){Double_t d=0.0,q;
                 for(Int_t i=0;i<3;i++){q = t[i]-ftran[i]; d += q*q;}
                 return d;}
 private: // private functions
	void MatrixFromAngle();
	void AngleFromMatrix();
 private: // Data members.
	Int_t    fDetectorIndex; // Detector type index (like fShapeIndex was)
	Int_t    fid[3];         // layer, ladder, detector numbers.
	Double_t frot[3];        // vector of rotations about x,y,z [radians].
	Double_t ftran[3];       // Translation vector of module x,y,z.
	Double_t fm[3][3];       // Rotation matrix based on frot.

	ClassDef(AliITSgeomMatrix,1) // Matrix class used by AliITSgeom.
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomMatrix &source);
istream &operator>>(istream &os,AliITSgeomMatrix &source);

#endif
