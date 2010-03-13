#ifndef ALI_TPC_DISTORTIONS_H
#define ALI_TPC_DISTORTIONS_H

#include <TObject.h>

class AliMagF;

class AliTPCDistortions : public TObject {
public:
  AliTPCDistortions();
  virtual ~AliTPCDistortions() {};
  
  void UndoTwistDistortion(        const Double_t x[],Double_t xprime[],Int_t roc);
  void UndoIFCShiftDistortion(     const Double_t x[],Double_t xprime[],Int_t roc);
  void UndoGGVoltErrorDistortion(  const Double_t x[],Double_t xprime[],Int_t roc);
  void UndoExBShapeDistortion(     const Double_t x[],Double_t xprime[],Int_t roc);

  // common setters and getters for ExB
  void SetOmegaTau(Double_t omegaTau) {fOmegaTau=omegaTau;RecalculateCs();}
  void SetT1(Double_t t1) {fT1=t1;RecalculateCs();}
  void SetT2(Double_t t2) {fT2=t2;RecalculateCs();}
  void SetC0(Double_t c0) {fC0=c0;}
  void SetC1(Double_t c1) {fC1=c1;}
  void SetC2(Double_t c2) {fC2=c2;}
  Double_t GetOmegaTau() const {return fOmegaTau;}
  Double_t GetT1() const {return fT1;}
  Double_t GetT2() const {return fT2;}
  Double_t GetC0() const {return fC0;}
  Double_t GetC1() const {return fC1;}
  Double_t GetC2() const {return fC2;}

  // setters and getters for twist
  void SetXTwist(Double_t xTwist) {fXTwist=xTwist;}
  void SetYTwist(Double_t yTwist) {fYTwist=yTwist;}
  Double_t GetXTwist() const {return fXTwist;}
  Double_t GetYTwist() const {return fYTwist;}

  // setter and getter for IFC
  void SetIFCShift(Double_t ifcShift) {fIFCShift=ifcShift;}
  Double_t GetIFCshift() const {return fIFCShift;}

  // setters and getters for GG
  void SetDeltaVGGA(Double_t deltaVGGA) {fDeltaVGGA=deltaVGGA;}
  void SetDeltaVGGC(Double_t deltaVGGC) {fDeltaVGGC=deltaVGGC;}
  Double_t GetDeltaVGGA() const {return fDeltaVGGA;}
  Double_t GetDeltaVGGC() const {return fDeltaVGGC;}

  // setter and getter for B field shape
  void SetBField(AliMagF *bField) {fBField=bField;}
  AliMagF* GetBField() const {return fBField;}
  static AliTPCDistortions*   Instance();
  static AliTPCDistortions*   fgInstance; //! Instance of this class (singleton implementation)
private:
  enum {kNR=   92};              // Number of R points in table
  enum {kNZ=  270};              // Number of Z points in the tables for interpolating distortion data
  enum {kNPhi= 19};              // Number of Phi points in table ( plus one extra for 360 == 0 )
  void RecalculateCs();
  void InitIFCShiftDistortion();
  void InitGGVoltErrorDistortion();

  void Interpolate2DEdistortion( const Int_t order, const Double_t r, const Double_t z, 
				 const Double_t er[kNZ][kNR], Double_t &er_value );
  Double_t Interpolate( const Double_t xArray[], const Double_t yArray[], 
		       const Int_t order, const Double_t x );
  void Search( const Int_t n, const Double_t xArray[], const Double_t x, Int_t &low );

  static const Double_t fgkIFCRadius;   // Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
  static const Double_t fgkOFCRadius;   // Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
  static const Double_t fgkTPC_Z0;      // Z location of TPC Gated Grid (cm)
  static const Double_t fgkZOffSet;     // Offset from CE: calculate all distortions closer to CE as if at this point
  static const Double_t fgkCathodeV;    // Cathode Voltage (volts)
  static const Double_t fgkGG;          // Gating Grid voltage (volts)
  static const Double_t fgkAliceDriftV; // Drift Velocity (cm/microSec) Magnitude
  static const Double_t fgkRList[kNR];
  static const Double_t fgkPhiList[kNPhi];
  static const Double_t fgkZList[kNZ];

  Double_t fOmegaTau;             // omega tau factor               (compare Jim Thomas's notes for definitions)
  Double_t fT1;                   // 1st tensor term for omega tau  (compare Jim Thomas's notes for definitions)
  Double_t fT2;                   // 2nd tensor term for omega tau  (compare Jim Thomas's notes for definitions)
  Double_t fC0;                   // coefficient C0                 (compare Jim Thomas's notes for definitions)
  Double_t fC1;                   // coefficient C1                 (compare Jim Thomas's notes for definitions)
  Double_t fC2;                   // coefficient C2                 (compare Jim Thomas's notes for definitions)
  AliMagF *fBField;              // Magnetic field map for ExB shape distortions due to B field
  Double_t fXTwist;               // Twist of E to B filed in X-Z [rad]
  Double_t fYTwist;               // Twist of E to B filed in Y-Z [rad]
  Double_t fIFCShift;             // Shift of inner field cage [cm]
  Double_t fDeltaVGGA;            // Missmatch of gating grid voltage on A-side [V]
  Double_t fDeltaVGGC;            // Missmatch of gating grid voltage on C-side [V]
  Double_t  fShiftER[kNZ][kNR];       // Array to store electric field for IFCShift calcualtion
  Double_t  fSpaceER[kNZ][kNR];       // Array to store electric field for SpaceCharge calculation
  Double_t  fShortER[kNZ][kNR];       // Array to store electric field for ShortedRing calculation
  Double_t  fGGVoltErrorER[kNZ][kNR]; // Array to store electric field for GGVoltError calculation
  Int_t fJLow;
  Int_t fKLow;

  ClassDef(AliTPCDistortions,1);
};

#endif
