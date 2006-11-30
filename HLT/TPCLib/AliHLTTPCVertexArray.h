// @(#) $Id$
// Original: AliHLTVertexArray.h,v 1.4 2004/06/11 16:06:33 loizides 

#ifndef AliHLTTPCVERTEXARRAY_H
#define AliHLTTPCVERTEXARRAY_H

#include <math.h>
#include "AliHLTTPCRootTypes.h"
 
class AliHLTTPCVertexArray {
  private:

  Char_t fArray[8320][8][8]; //array
  Double_t fZSector;    //sector
  Double_t fZSectorErr; //sector error
  Int_t fMaxSeed;  //max seed
  Int_t fNSeed;    //number of seeds
  Float_t fZSeed[400]; //seed in Z
  Float_t fRSeed[400]; //seed in XY
  Int_t fSecSeed[400]; //seed for sectors

  void FindMean(Float_t *vertex,Int_t *array, Int_t len);
  void AnalyzeSector(Float_t *vertex, Int_t *array, Int_t len);

  public:
  AliHLTTPCVertexArray(){fNSeed=0;fMaxSeed=400;}
  AliHLTTPCVertexArray(AliHLTTPCVertexArray&){fNSeed=0;fMaxSeed=400;}
  AliHLTTPCVertexArray(Int_t maxseed){fNSeed=0;fMaxSeed=maxseed;}
  virtual ~AliHLTTPCVertexArray(){;}
  
  Int_t GetContent(Float_t z,Float_t r,Int_t sec);
  Int_t Trace(Float_t z,Float_t r,Int_t sec,Float_t vertex);

  Double_t GetZSectorErr(){return fZSectorErr;}
  Double_t GetZSector(){return fZSector;}
  Int_t GetMaxSeed(){return fMaxSeed;}
  Int_t GetNSeed(){return fNSeed;}
  Float_t GetRSeed(Int_t i){if(i<400) return fRSeed[i]; else return -999;}
  Float_t GetZSeed(Int_t i){if(i<400) return fZSeed[i]; else return -999;} 
  Float_t GetSecSeed(Int_t i){if(i<400) return fSecSeed[i]; else return -999;}

  void FillSector2D(Float_t z,Float_t r,Int_t sec);
  void FillSectorSeed2D(Float_t z,Float_t r,Int_t sec);
  void FillSector3D(Float_t x,Float_t y, Float_t z);
  void FillSectorSeed3D(Float_t x,Float_t y, Float_t z);
  void FindSectorVertex(Double_t pos = 0,Double_t range = 60,Int_t nbin = 60);
  void ResetSector();

  ClassDef(AliHLTTPCVertexArray,1)  //The HLTTPC Fast Vertex Finder Base Class

};

inline void AliHLTTPCVertexArray::FillSector3D(Float_t x, Float_t y, Float_t z)
{
  // Filling routine in coordinates
  Int_t sec = Int_t( (y+.168*x)/(.336*x)*8); // 8 subsec!!
  Float_t r = sqrt(pow(y,2)+pow(x,2));
  FillSector2D(z,r,sec); 
}

inline void AliHLTTPCVertexArray:: FillSectorSeed3D(Float_t x,Float_t y, Float_t z)
{
  // Filling routine for seeds in coordinates
  Int_t sec = Int_t( (y+.168*x)/(.336*x)*8); // 8 subsec!!
  Float_t r = sqrt(pow(y,2)+pow(x,2));
  FillSectorSeed2D(z,r,sec);    
}

inline void AliHLTTPCVertexArray::FillSectorSeed2D(Float_t z,Float_t r,Int_t sec)
{
  // Filling routine in r,z coordinates 
  if(fNSeed>=400) return;
  fZSeed[fNSeed] = z; fRSeed[fNSeed] = r; fSecSeed[fNSeed] = sec;
  fNSeed++; 
}

inline void AliHLTTPCVertexArray::FillSector2D(Float_t z,Float_t r,Int_t sec)
{
  // Filling routine for seeds in r,z coordinates
  if(z>r||z<=0||r<220||r>=252) return;
  fArray[Int_t(z/r*32*260)][(Int_t(r-220))/4][sec] += 1;
}

inline Int_t AliHLTTPCVertexArray::GetContent(Float_t z,Float_t r,Int_t sec)
{
  // Return content of array in r,z coordinates
  if(z>r||z<=0||r<220||r>=252) return 0;
  return  fArray[Int_t(z/r*32*260)][(Int_t(r-220))/4][sec];
}

inline void AliHLTTPCVertexArray::ResetSector()
{
  // do it!
  fZSector=0;
  fZSectorErr=0;
  fNSeed=0;
  for(Int_t i =0;i<400;i++)
    fZSeed[i] =  fRSeed[i] =  fSecSeed[i] = 0;
  for(Int_t z=0;z<8320;z++)
    for(Int_t r=0;r<8;r++)
      for(Int_t sec =0;sec<8;sec++)
        fArray[z][r][sec] = 0;
}

inline Int_t AliHLTTPCVertexArray::Trace(Float_t z,Float_t r,Int_t sec,Float_t vertex)
{
// count the number of entries along starting from z,r to vertex,0
  Int_t cont=0;
  for(Int_t i = 0;i<8;i++){
    Float_t ry = 222 +(i*4);
    Float_t zx = (ry/r)*(z-vertex)+vertex;
    cont += GetContent(zx,ry,sec);
  }
  if(cont < 5) return 0;
  return cont;
}

#endif
