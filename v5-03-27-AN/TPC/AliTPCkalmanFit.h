#ifndef ALITPCKALMANFIT_H
#define ALITPCKALMANFIT_H



#include "TNamed.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjArray.h"
class TTreeSRedirector;
class AliTrackPointArray;
class AliTrackPoint;
class TFormula;
class TBits;
class THnSparse;
//

class AliTPCkalmanFit: public TNamed{
public:
  AliTPCkalmanFit();
  void Init();
  void InitTransformation();
  void Add(const AliTPCkalmanFit * kalman);

  void  AddCalibration(AliTPCTransformation * calib);
  AliTPCTransformation * GetTransformation(Int_t i){return (fCalibration)? (AliTPCTransformation *)fCalibration->At(i):0;}  
  Int_t GetTransformationIndex(const char * trName);
  //
  void SetStatus(const char * name, Bool_t setOn, Bool_t isOr=kTRUE);
  //
  void FitTrackLinear(AliTrackPointArray& points,  TTreeSRedirector *debug=0, Float_t scalingRMSY=1., Float_t scalingRMSZ=1.);
  void DumpTrackLinear(AliTrackPointArray& points, TTreeSRedirector *debug);
  void UpdateLinear(AliTrackPoint &point, TTreeSRedirector *debug=0);

  void Propagate(TTreeSRedirector *debug=0);
  void PropagateTime(Int_t time);
  void AddCovariance(const char * varName, Double_t sigma);
  void Update(const AliTPCkalmanFit * kalman);

  static AliTrackPointArray * SortPoints(AliTrackPointArray &points);
  static AliTrackPointArray * MakePointArrayLinear(Double_t alpha, Double_t y0, Double_t z0, Double_t ky, Double_t kz, Double_t err=0.02); 
  void  ApplyCalibration(AliTrackPointArray *array, Double_t csign);
  Bool_t  CheckCovariance(TMatrixD &covar, Float_t maxEl);
  
  Bool_t DumpCorelation(Double_t threshold, const char *mask0=0, const char *mask1=0);
  Bool_t DumpCalib(const char *mask=0, Float_t correlationCut=-1);
  //
  //Double_t GetTPCDeltaXYZ(Int_t coord, Int_t volID, Double_t x, Double_t y, Double_t z);
  //static Double_t SGetTPCDeltaXYZ(Int_t coord, Int_t volID, Double_t x, Double_t y, Double_t z);

  Double_t GetTPCDeltaXYZ(Int_t coord, Int_t volID, Int_t icoordsys, Double_t x, Double_t y, Double_t z);
  static Double_t SGetTPCDeltaXYZ(Int_t coord, Int_t volID, Int_t icoordsys, Double_t x, Double_t y, Double_t z);
  Double_t GetTPCtransXYZ(Int_t coord, Int_t volID, Int_t calibID, Int_t icoordsys, Double_t x, Double_t y, Double_t z);
  static Double_t SGetTPCtransXYZ(Int_t coord, Int_t volID, Int_t calibID, Int_t icoordsys, Double_t x, Double_t y, Double_t z);
  void MakeTreeTrans(TTreeSRedirector *debug, const char *treeName="all");

  AliTPCkalmanFit *Test(Int_t ntracks);
 public:
  //
  // Calibration parameters
  //
  TObjArray *fCalibration;  // array of calibrations
  TMatrixD  *fCalibParam;   // calibration parameters 
  TMatrixD  *fCalibCovar;   // calibration parameters 
  //
  // Linear track
  //
  TMatrixD  *fLinearParam;      // linear parameters
  TMatrixD  *fLinearCovar;      // linear covariance
  THnSparse *fLinearTrackDelta[12];   // linear tracks matching residuals - delta 
  THnSparse *fLinearTrackPull[12];    // linear tracks matching residuals  - pull
  //
  //
  //
  Int_t      fLastTimeStamp; // last time stamp - used for propagation of parameters
  //static AliTPCkalmanFit* Instance();
  void SetInstance(AliTPCkalmanFit*param){fgInstance = param;}
  static AliTPCkalmanFit*   fgInstance; //! Instance of this class (singleton implementation)
 private:  
  Double_t   fCurrentAlpha; //! current rotation frame
  Double_t   fCA;           //! cosine of current angle
  Double_t   fSA;           //! sinus of current angle  
  AliTPCkalmanFit&  operator=(const AliTPCkalmanFit&);// not implemented
  AliTPCkalmanFit(const AliTPCkalmanFit&):TNamed(),
                                  fCalibration(0),
                                  fCalibParam(0),
                                  fCalibCovar(0),
                                  fLinearParam(0),
                                  fLinearCovar(0),
                                  fLastTimeStamp(-1),
                                  fCurrentAlpha(0),  //! current rotation frame
                                  fCA(0),     //! cosine of current angle
                                  fSA(0)     //! sinus of current angle  
                                          
{;} //not implemented
  ClassDef(AliTPCkalmanFit,3);
};



#endif

