// @(#) $Id$

#ifndef ALIL3_HOUGHTRANSFORMERLUT
#define ALIL3_HOUGHTRANSFORMERLUT

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"
#include "AliL3Histogram.h"

/* use if you also want the eta index 
   to be looked up and linear searched */
//#define FULLLUT

class AliL3HoughTransformerLUT : public AliL3HoughBaseTransformer {
  
 protected:
  AliL3Histogram **fParamSpace; //!

  void DeleteHistograms();
  
  Int_t fMinRow;
  Int_t fMaxRow;
  Int_t fNRows;
  Int_t fNEtas;
  Int_t fNPhi0;
  Int_t fSlice;
  Int_t fSector;
  Int_t fSectorRow;
  Int_t fZSign;
  Float_t fZLengthPlusOff;
  Float_t fTimeWidth;
  Float_t fPadPitch;
  Float_t fEtaSlice;

  Float_t *fLUTX; //!
  Float_t *fLUTY; //!
  Float_t *fLUTEta; //!
  Float_t *fLUTEtaReal; //!
  Float_t *fLUTphi0; //!
  Float_t *fLUT2sinphi0; //!   
  Float_t *fLUT2cosphi0; //!
  //not used but need for VHDL version
  Float_t *fLUTKappa; //!
  
  Int_t fLastPad;
  Int_t fLastIndex;
  Int_t fAccCharge;
  Float_t fX,fY; //trafo values per pad

  Float_t CalcRoverZ2(Float_t eta);
  Float_t CalcEta(Float_t roverz2);
  Float_t CalcX(Int_t row);
  Float_t CalcY(Int_t pad, Int_t row);
  Float_t CalcZ(Int_t time);  
  Int_t FindIndex(Float_t rz2, Int_t start=-100);
  void AddCurveToHistogram(Int_t new_eta_index=-1);

 public:

  AliL3HoughTransformerLUT(); 
  AliL3HoughTransformerLUT(Int_t slice,Int_t patch,Int_t n_eta_segments);
  virtual ~AliL3HoughTransformerLUT();
  
  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t pres,Int_t nybin,Float_t psi);
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,Int_t nybin,Float_t ymin,Float_t ymax);
  void Reset();

  void TransformCircle();

  Int_t GetEtaIndex(Double_t eta);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  Double_t GetEta(Int_t eta_index,Int_t slice);
  
  void Print();
  void Init(Int_t slice=0,Int_t patch=0,Int_t n_eta_segments=100,Int_t n_seqs=-1);

  ClassDef(AliL3HoughTransformerLUT,1) //LUT Hough transformation class

};

inline Float_t AliL3HoughTransformerLUT::CalcRoverZ2(Float_t eta)
{
  Float_t e=exp(2*eta);
  Float_t ret=(e+1)/(e-1);
  ret*=ret;
  return ret;
}

inline Float_t AliL3HoughTransformerLUT::CalcEta(Float_t roverz2)
{
  Float_t rz=sqrt(roverz2);
  if(fZSign<0) rz=-rz;
  Float_t ret=(1+rz)/(rz-1);
  ret=0.5*log(ret);
  return ret;
}

inline Float_t AliL3HoughTransformerLUT::CalcX(Int_t row)
{
  return fLUTX[row];
}

inline Float_t AliL3HoughTransformerLUT::CalcY(Int_t pad,Int_t row)
{
  return pad*fPadPitch-fLUTY[row];
}

inline Float_t AliL3HoughTransformerLUT::CalcZ(Int_t time)
{
  Float_t ret=time*fTimeWidth;
  if(fZSign>0) ret=fZLengthPlusOff-ret;
  else ret=ret-fZLengthPlusOff;
  return ret;
}

inline void AliL3HoughTransformerLUT::AddCurveToHistogram(Int_t new_eta_index)
{
  //get correct histogrampointer
  AliL3Histogram *hist = fParamSpace[fLastIndex];

  //Fill the histogram along the phirange
  Float_t r2=fX*fX+fY*fY;
  for(Int_t b=0; b<fNPhi0; b++){
    Float_t kappa=(fY*fLUT2cosphi0[b]-fX*fLUT2sinphi0[b])/r2;
    hist->Fill(kappa,fLUTphi0[b],fAccCharge);
    //cout << kappa << " " << fLUTphi0[b] << " " << fAccCharge << endl;
  }

  fAccCharge=0;
  fLastIndex=new_eta_index;
}

inline Int_t AliL3HoughTransformerLUT::FindIndex(Float_t rz2, Int_t start)
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
