// @(#) $Id$

#ifndef ALIL3_HOUGHTRANSFORMERLUT
#define ALIL3_HOUGHTRANSFORMERLUT

#include "AliHLTRootTypes.h"
#include "AliHLTHoughBaseTransformer.h"
#include "AliHLTHistogram.h"

/* use if you also want the eta index 
   to be looked up and linear searched */
//#define FULLLUT

class AliHLTHoughTransformerLUT : public AliHLTHoughBaseTransformer {
  
 public:

  AliHLTHoughTransformerLUT(); 
  AliHLTHoughTransformerLUT(Int_t slice,Int_t patch,Int_t nEtaSegments);
  virtual ~AliHLTHoughTransformerLUT();
  
  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t pres,Int_t nybin,Float_t psi);
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();

  void TransformCircle();
  void TransformCircle(Int_t *row_range,Int_t every) {
    AliHLTHoughBaseTransformer::TransformCircle(row_range,every);
  }

  Int_t GetEtaIndex(Double_t eta) const;
  AliHLTHistogram *GetHistogram(Int_t etaIndex);
  Double_t GetEta(Int_t etaIndex,Int_t slice) const;
  
  void Print();
  void Init(Int_t slice=0,Int_t patch=0,Int_t nEtaSegments=100,Int_t nSeqs=-1);

 protected:
  AliHLTHistogram **fParamSpace; //!

  void DeleteHistograms();
  
  Int_t fMinRow; // Min row (= first row)
  Int_t fMaxRow; // Max row (= last row)
  Int_t fNRows;  // Number of rows
  Int_t fNEtas;  // Number of eta slices
  Int_t fNPhi0;  // Number of phi bins
  Int_t fSlice;  // Current slice
  Int_t fSector; // Current sector
  Int_t fSectorRow; // Sector row (?)
  Int_t fZSign; // Z sign
  Float_t fZLengthPlusOff; // Z lenght plus offset
  Float_t fTimeWidth; // Time width
  Float_t fPadPitch; // Pad pitch
  Float_t fEtaSlice; // Eta slice

  Float_t *fLUTX; //! LUT for X
  Float_t *fLUTY; //! LUT for Y
  Float_t *fLUTEta; //! LUT for eta
  Float_t *fLUTEtaReal; //! LUT for real eta (?)
  Float_t *fLUTphi0; //! LUT for phi0
  Float_t *fLUT2sinphi0; //! LUT for sin(phi0)
  Float_t *fLUT2cosphi0; //! LUT for cos(phi0)
  //not used but need for VHDL version
  Float_t *fLUTKappa; //! LUT for kappa
  
  Int_t fLastPad; // Last pad
  Int_t fLastIndex; // Last index
  Int_t fAccCharge; // Accepted charge
  Float_t fX,fY; //trafo values per pad

  Float_t CalcRoverZ2(Float_t eta) const;
  Float_t CalcEta(Float_t roverz2) const;
  Float_t CalcX(Int_t row) const;
  Float_t CalcY(Int_t pad, Int_t row) const;
  Float_t CalcZ(Int_t time) const;  
  Int_t FindIndex(Float_t rz2, Int_t start=-100) const;
  void AddCurveToHistogram(Int_t newEtaIndex=-1);

  ClassDef(AliHLTHoughTransformerLUT,1) //LUT Hough transformation class

};

typedef AliHLTHoughTransformerLUT AliL3HoughTransformerLUT; // for backward comaptibility

inline Float_t AliHLTHoughTransformerLUT::CalcRoverZ2(Float_t eta) const
{
  Float_t e=exp(2*eta);
  Float_t ret=(e+1)/(e-1);
  ret*=ret;
  return ret;
}

inline Float_t AliHLTHoughTransformerLUT::CalcEta(Float_t roverz2) const
{
  Float_t rz=sqrt(roverz2);
  if(fZSign<0) rz=-rz;
  Float_t ret=(1+rz)/(rz-1);
  ret=0.5*log(ret);
  return ret;
}

inline Float_t AliHLTHoughTransformerLUT::CalcX(Int_t row) const
{
  return fLUTX[row];
}

inline Float_t AliHLTHoughTransformerLUT::CalcY(Int_t pad,Int_t row) const
{
  return pad*fPadPitch-fLUTY[row];
}

inline Float_t AliHLTHoughTransformerLUT::CalcZ(Int_t time) const
{
  Float_t ret=time*fTimeWidth;
  if(fZSign>0) ret=fZLengthPlusOff-ret;
  else ret=ret-fZLengthPlusOff;
  return ret;
}

inline void AliHLTHoughTransformerLUT::AddCurveToHistogram(Int_t newEtaIndex)
{
  //get correct histogrampointer
  AliHLTHistogram *hist = fParamSpace[fLastIndex];

  //Fill the histogram along the phirange
  Float_t r2=fX*fX+fY*fY;
  for(Int_t b=0; b<fNPhi0; b++){
    Float_t kappa=(fY*fLUT2cosphi0[b]-fX*fLUT2sinphi0[b])/r2;
    hist->Fill(kappa,b+1,fAccCharge);
    //cout << kappa << " " << fLUTphi0[b] << " " << fAccCharge << endl;
  }

  fAccCharge=0;
  fLastIndex=newEtaIndex;
}

inline Int_t AliHLTHoughTransformerLUT::FindIndex(Float_t rz2, Int_t start) const
{
  //could improve search through devide and conquere strategy
  
  Int_t index=start; 
  if(index==-100){
    index=0;
    while((index<fNEtas)&&(rz2<=fLUTEta[index])){
      index++;
    }
  } else {
    while((index>=0)&&(rz2>fLUTEta[index])){
      index--;
    }
    index++;
  }
  //cout << start << " - " << index << " " << ": " << rz2 << " " << fLUTEta[index] << endl;

  return index;
}




#endif
