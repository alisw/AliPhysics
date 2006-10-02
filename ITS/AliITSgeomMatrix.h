#ifndef ALIITSGEOMMATRIX_H
#define ALIITSGEOMMATRIX_H
/* Copyright(c) 2000, ALICE Experiment at CERN, All rights reserved. *
 * see cxx source for full Copyright notice.                         */

//////////////////////////////////////////////////////////////////////////
//  ITS geometry manipulation routines on the module level. This class is
//  to replace the structure ITS_geom in the class AliITSgeom.
//  Created May 30 2000.
//   version 0.0.0
// By Bjorn S. Nilsen
////////////////////////////////////////////////////////////////////////////
#include <TObject.h>
#include <TString.h>
class TPolyLine3D;
class TNode;
class TShape;

class AliITSgeomMatrix : public TObject {
 public:
	AliITSgeomMatrix(); // Default constructor
	// Standard constructor #1
	AliITSgeomMatrix(Int_t idt,const Int_t id[3],
			 const Double_t rot[3],const Double_t tran[3]);
	// Standard constructor #2
        AliITSgeomMatrix(Int_t idt,const Int_t id[3],
		         Double_t matrix[3][3],const Double_t tran[3]);
	// Standard constructor #3
        AliITSgeomMatrix(const Double_t rotd[6]/*degrees GEANT angles*/,
                         Int_t idt,const Int_t id[3],
                         const Double_t tran[3]);
	// Copy constructor
	AliITSgeomMatrix(const AliITSgeomMatrix &source);
	// Assignment operator
	AliITSgeomMatrix& operator=(const AliITSgeomMatrix &source); 
	virtual ~AliITSgeomMatrix(){}; // default constructor.
	// Prints a line describing the output format of the function Print.
	void PrintComment(ostream *os) const;
	// Prints out the content of this class in ASCII format.
	void Print(ostream *os)const;
	// Prints out the content of this class in ASCII format but includes
	// formating and strings that make it more humanly readable.
	void PrintTitles(ostream *os) const;
	// Reads in the content of this class in the format of Print
	void Read(istream *is);
     virtual void Print(Option_t *option="") const {
                       TObject::Print(option);}
     virtual Int_t Read(const char *name) {return TObject::Read(name);}

     // Returns the geometry path corresponding to this transformation
     TString& GetPath(){return fPath;}
     // Sets the geometry path
     void SetPath(const Char_t *p){fPath = p;}
     void SetPath(const TString &p){fPath = p;}
	// Given the rotation angles [radians] it fills frot and computes
	// the rotation matrix fm.
	void SetAngles(const Double_t rot[3]){// [radians]
              for(Int_t i=0;i<3;i++)frot[i] = rot[i];this->MatrixFromAngle();}
	// Sets the translation vector and computes fCylR and fCylPhi.
	void SetTranslation(const Double_t tran[3]);
	// sets the rotation matrix and computes the rotation angles [radians]
	void SetMatrix(const Double_t matrix[3][3]){ for(Int_t i=0;i<3;i++)
	 for(Int_t j=0;j<3;j++) fm[i][j]=matrix[i][j];this->AngleFromMatrix();}
	// Sets the detector index value
	void SetDetectorIndex(Int_t idt) {fDetectorIndex = idt;}
	// Sets the detector layer, ladder, detector (id) values.
	void SetIndex(const Int_t id[3]){
	                   for(Int_t i=0;i<3;i++) fid[i] = id[i];}
	// Returns the rotation angles [radians]
	void GetAngles(Double_t rot[3]) const {// [radians]
	                   for(Int_t i=0;i<3;i++)  rot[i] = frot[i];}
	// Returns the translation vector [cm]
	void GetTranslation(Double_t tran[3]) const {
	                    for(Int_t i=0;i<3;i++) tran[i] = ftran[i];}
	// Returns the translation vector in cylindrical
	// coordinates [cm,radians]
	void GetTranslationCylinderical (Double_t tran[3]) const {
	                    tran[0] = fCylR;
			    tran[1] = fCylPhi;
			    tran[2] = ftran[2];}
	// Returns the values of the rotation matrix
	void GetMatrix(Double_t matrix[3][3]) const {for(Int_t i=0;i<3;i++)
		         for(Int_t j=0;j<3;j++) matrix[i][j] = fm[i][j];}
	// Returns the detector index value.
	Int_t GetDetectorIndex() const {return fDetectorIndex;}
	// returns the modules index layer, ladder, detector
	void  GetIndex(Int_t id[3]) const {for(Int_t i=0;i<3;i++) id[i] = fid[i];}
	// Sets the rotation matrix based on the 6 GEANT rotation
	// angles [radian]
	void  MatrixFromSixAngles(const Double_t *ang);
	// Returns the 6 GEANT rotation angles [radians] from the
	// existing rotation matrix.
	void  SixAnglesFromMatrix(Double_t *ang)const;

	// Given a position in Cartesian ALICE global coordinates [cm]
	// returns the position in Cartesian detector/module local
	//coordinates [cm]
	void GtoLPosition(const Double_t g[3],Double_t l[3]) const;
	// Given a position in Cartesian detector/module local coordinates [cm]
	// returns the position in Cartesian ALICE global
	//coordinates [cm]
	void LtoGPosition(const Double_t l[3],Double_t g[3]) const;
	// Given a momentum in Cartesian ALICE global coordinates
	// returns the momentum in Cartesian detector/module local
	//coordinates
	void GtoLMomentum(const Double_t g[3],Double_t l[3]) const;
	// Given a momentum in Cartesian detector/module local coordinates 
	// returns the momentum in Cartesian ALICE global coordinates
	void LtoGMomentum(const Double_t l[3],Double_t g[3]) const;
	// given a position error matrix in ALICE Cartesian global
	// coordinates [cm] returns a position error matrix in detector/
	// module local Cartesian local coordinates [cm]
	void GtoLPositionError(const Double_t g[3][3],Double_t l[3][3]) const;
	// given a position error matrix in detector/module Cartesian local
	// coordinates [cm] returns a position error matrix in ALICE
	// Cartesian global coordinates [cm]
	void LtoGPositionError(const Double_t l[3][3],Double_t g[3][3]) const;
	// Tracking Related Routines
	void GtoLPositionTracking(const Double_t g[3],Double_t l[3]) const;
	// Given a position in Cartesian Tracking global coordinates [cm]
	// returns the position in Cartesian detector/module local
	// coordinates [cm]
	void LtoGPositionTracking(const Double_t l[3],Double_t g[3]) const;
	// Given a position in Cartesian detector/module local coordinates [cm]
	// returns the position in Cartesian Tracking global
	//coordinates [cm]
	void GtoLMomentumTracking(const Double_t g[3],Double_t l[3]) const;
	// Given a momentum in Cartesian detector/module local coordinates 
	// returns the momentum in Cartesian Tracking global coordinates
	void LtoGMomentumTracking(const Double_t l[3],Double_t g[3]) const;
	// given a position error matrix in Tracking Cartesian global
	// coordinates [cm] returns a position error matrix in detector/
	// module local Cartesian local coordinates [cm]
	void GtoLPositionErrorTracking(const Double_t g[3][3],
				       Double_t l[3][3]) const;
	// given a position error matrix in detector/module Cartesian local
	// coordinates [cm] returns a position error matrix in Tracking
	// Cartesian global coordinates [cm]
	void LtoGPositionErrorTracking(const Double_t l[3][3],
				       Double_t g[3][3]) const;
	// Computes the distance squared [cm^2] between a point t[3] and
	// this module/detector
	Double_t Distance2(const Double_t t[3]) const {Double_t d=0.0,q;
                 for(Int_t i=0;i<3;i++){q = t[i]-ftran[i]; d += q*q;}
                 return d;}
     //
     // Documentation related Class
     TPolyLine3D* CreateLocalAxis() const;
     TPolyLine3D* CreateLocalAxisTracking() const;
     TNode* CreateNode(const Char_t *nodeName,const Char_t *nodeTitle,
                       TNode *mother,TShape *shape,Bool_t axis=kTRUE) const;
     void MakeFigures() const;
     //
 private: // private functions
	// Given the rotation matrix fm it fills the rotation angles frot
	void MatrixFromAngle();
	// Given the rotation angles frot it fills the rotation matrix fm
	void AngleFromMatrix();
 private: // Data members.
	Int_t    fDetectorIndex; // Detector type index (like fShapeIndex was)
	Int_t    fid[3];         // layer, ladder, detector numbers.
	Double_t frot[3];        //! vector of rotations about x,y,z [radians].
	Double_t ftran[3];       // Translation vector of module x,y,z.
	Double_t fCylR,fCylPhi;  //! Translation vector in Cylindrical coord.
	Double_t fm[3][3];       // Rotation matrix based on frot.
	TString  fPath;          // Path within Geometry to this volume

	// Note, fCylR and fCylPhi are added as data members because it costs
	// about a factor of 10 to compute them over looking them up. Since
	// they are used in some tracking algorithms this can be a large cost
	// in computing time. They are not written out but computed.

	ClassDef(AliITSgeomMatrix,2) // Matrix class used by AliITSgeom.
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,AliITSgeomMatrix &source);
istream &operator>>(istream &os,AliITSgeomMatrix &source);

#endif
