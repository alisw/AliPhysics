#ifndef ALIL3_HOUGHTRANSFORMERVHDL
#define ALIL3_HOUGHTRANSFORMERVDHL

#include <stream.h>

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

#ifdef USEFFLOAT
#include "AliL3FFloat.h"
#else
typedef Float_t AliL3FFloat;
//#define AliL3FFloat float
#endif

class AliL3Histogram;

class AliL3HoughTransformerVhdl : public AliL3HoughBaseTransformer 
{

 private:

  AliL3Histogram **fParamSpace; //!
#ifdef VHDLVERSION
  Int_t fMinRow;
  Int_t fMaxRow;
  Int_t fNRows;
  Int_t fNEtas;
  Int_t fNPhi0;
  Int_t fSector;
  Int_t fSectorRow;
  Int_t fZSign;
  AliL3FFloat fZLengthPlusOff;
  AliL3FFloat fTimeWidth;
  AliL3FFloat fPadPitch;
  AliL3FFloat fEtaSlice;

  //VESTBO: Here is the problem: it seems that
  //values into the arrays somehow overwrite 
  //other data member. And that even in the
  //static version, where I dont use heap!!!
#ifdef VHDLSTATIC
  float fLUTX[32]; 
  float fLUTY[32]; 
  float fLUTEta[256]; 
  float fLUTphi0[256]; 
  float fLUT2sinphi0[256];   
  float fLUT2cosphi0[256];
#else
  AliL3FFloat *fLUTX; //!
  AliL3FFloat *fLUTY; //!
  AliL3FFloat *fLUTEta; //!
  AliL3FFloat *fLUTphi0; //!
  AliL3FFloat *fLUT2sinphi0; //!   
  AliL3FFloat *fLUT2cosphi0; //!
#endif
  
  Float_t CalcRoverZ2(Float_t eta);
  Float_t CalcEta(Float_t roverz2);
  Float_t CalcX(Int_t row);
  Float_t CalcY(Int_t pad, Int_t row);
  Float_t CalcZ(Int_t time);  

  Int_t FindIndex(Double_t rz2);
#endif
  void DeleteHistograms();


 public:

  AliL3HoughTransformerVhdl(); 
  AliL3HoughTransformerVhdl(Int_t slice,Int_t patch,Int_t n_eta_segments);
  virtual ~AliL3HoughTransformerVhdl();

  void CreateHistograms(Int_t nxbin,Double_t ptmin,Int_t nybin,Double_t phimin,Double_t phimax);
  void CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,
			Int_t nybin,Double_t ymin,Double_t ymax);
  void Reset();
  void TransformCircle();
  void TransformCircleC(Int_t row_range) {cerr<<"TransformCircleC is not defined!"<<endl;}
  void TransformLine() {cerr<<"TransformLine is not defined!"<<endl;}

  Int_t GetEtaIndex(Double_t eta);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  Double_t GetEta(Int_t eta_index,Int_t slice);

#ifdef VHDLVERSION
  void Print();
  void Init(Int_t slice=0,Int_t patch=0,Int_t n_eta_segments=100);
#endif

  ClassDef(AliL3HoughTransformerVhdl,1) //Normal Hough transformation class

};

#endif






