#ifndef ALIITSNEURALPOINT_H
#define ALIITSNEURALPOINT_H
///////////////////////////////////////////////////////////////
// AliITSneuralPoint                                         //
//                                                           //
// A class which resumes the information of ITS clusters     //
// in the global reference frame.                            //
// Author: A. Pulvirenti                                     //
///////////////////////////////////////////////////////////////
class AliITSgeom;
class AliITSgeomMatrix;
class AliITSRecPoint;
class AliITSRecPoint;

class AliITSNeuralPoint : public TObject {

public:

	AliITSNeuralPoint();
	AliITSNeuralPoint(AliITSNeuralPoint *p);
	AliITSNeuralPoint(AliITSRecPoint *rp, AliITSgeomMatrix *gm);
	AliITSNeuralPoint(AliITSRecPoint *rp, AliITSgeom *geom, Short_t module, Short_t index);

	virtual ~AliITSNeuralPoint() { }

	Double_t& X()    {return fX;}   // reference to X coord
	Double_t& Y()    {return fY;}   // reference to Y coord
	Double_t& Z()    {return fZ;}   // reference to Z coord
	Double_t& ErrX() {return fEX;}  // reference to X error
	Double_t& ErrY() {return fEY;}  // reference to Y error
	Double_t& ErrZ() {return fEZ;}  // reference to Z error
		
	Double_t  GetR2()    const  {return TMath::Sqrt(GetR2sq());} // xy radius
	Double_t  GetR3()    const  {return TMath::Sqrt(GetR3sq());} // 3D radius
	Double_t  GetR2sq()  const  {return fX*fX+fY*fY;}            // xy rad. square
	Double_t  GetR3sq()  const  {return GetR2sq()+fZ*fZ;}        // 3D rad. square
	Double_t  GetPhi()   const;
	Double_t  GetTheta() const  {return TMath::ATan2(GetR2(),fZ);} // polar angle
	Double_t  GetConfX() const  {return fConfX;}
	Double_t  GetConfY() const  {return fConfY;}
	Double_t  GetError(Option_t *opt);
	void      ConfMap(Double_t vx, Double_t vy);

	Double_t  GetCharge()       const {return fCharge;}        // ADC signal
	Short_t   GetIndex()        const {return fIndex;}         // Reference in TreeR
	Long_t    GetLabel(Int_t i) const {return fLabel[Chk(i)];} // GEANT owner particle
	Short_t   GetLayer()        const {return fLayer;}         // ITS layer
	Short_t   GetModule()       const {return fModule;}        // ITS module 
	Short_t   GetUser()         const {return fUser;}          // Found track owner
				
	void      SetCharge(Double_t val)       {fCharge = val;}
	void      SetIndex(Short_t val)         {fIndex = val;}
	void      SetLabel(Int_t i, Long_t val) {fLabel[Chk(i)] = val;}
	void      SetLayer(Short_t val)         {fLayer = val;}
	void      SetModule(Short_t val)        {fModule = val;}
	void      SetUser(Short_t val)          {fUser = val;}
	
	Bool_t    HasID (Int_t ID) const;
	Int_t*    SharedID(AliITSNeuralPoint *p) const;

protected:
	
	Int_t     Chk(Int_t i) const {if(i<0)i=0;if(i>=3)i=3;return i;}

	Double_t  fX;   // position 
	Double_t  fY;   // position
	Double_t  fZ;   // position
	
	Double_t  fConfX; // conformal mapping X
	Double_t  fConfY; // conformal mapping Y
		
	Double_t  fEX;  // position error
	Double_t  fEY;  // position error
	Double_t  fEZ;  // position error

	Double_t  fCharge;   // total charge signal in cluster

	Short_t   fModule;   // ITS module containing the point (0 - 2197)
	Short_t   fIndex;    // index as TClonesArray entry in TreeR (usually not > 600)
	Short_t   fLayer;    // ITS layer containing the point
	Short_t   fUser;     // owner recognized track or flag to evidence using
	Short_t   fZSort;    // order as a function of local Z

	Int_t     fLabel[3]; // GEANT labels of the owner tracks

	ClassDef(AliITSNeuralPoint, 1) // AliITSNeuralPoints class
};

#endif
