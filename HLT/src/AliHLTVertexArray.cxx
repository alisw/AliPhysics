// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTVertexArray.h"


/** \class AliHLTVertexArray
<pre>
//_____________________________________________________________
// AliHLTVertexArray
//
// The L3 Fast Vertex Finder Base Class
//
// usage:
//
//for(Int_t sec=0;sec<NSEC;sec++){
//  ResetSector();
//  FillSectorSeed3D(x,y,z);
//  FillSector3D(x,y,z);
//  FindSectorVertex();
//  Double_t z = GetZSector();
//  Double_t zerr = GetZSectorErr();
//  // do somethink with z, zerr
//}
</pre>
*/

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTVertexArray)

void AliHLTVertexArray::AnalyzeSector(Float_t *vertex, Int_t *array, Int_t len)
{
  //loop over all seeds and all vertex position
  LOG(AliHLTLog::kInformational,"AliHLTVertexArray::AnalyzeSector","Analyze")
  <<AliHLTLog::kDec<<"Number of Seeds: "<<fNSeed<<ENDLOG;
  for(Int_t i =0;i<fNSeed;i++)
    for(Int_t bin = 0;bin<len;bin++)
      array[bin] += Trace(fZSeed[i],fRSeed[i],fSecSeed[i],vertex[bin]);
}

void AliHLTVertexArray::FindSectorVertex(Double_t pos, Double_t range, Int_t nbin)
{
  //define position and range for search and
  //loop over all seeds
  //and find position of vertex and error 
  const Int_t klen = nbin;
  const Double_t kwidth = range; 
  const Double_t kxmin = pos - kwidth/2;
  const Double_t kstep = kwidth/klen;
  const Double_t kstart = kxmin + kstep/2.;
  Int_t * array = new Int_t[klen];
  Float_t * ver   = new Float_t[klen];
  for(Int_t i=0;i<klen;i++){
    ver[i] =  kstart + kstep * i;
    array[i] = 0;
  }
  AnalyzeSector(ver,array,klen);
  FindMean(ver,array,klen);
  delete[] array;
  delete[] ver;
}

void AliHLTVertexArray::FindMean(Float_t *vertex,Int_t *array, Int_t len)
{
  //find mean and error of array and store it in
  //fZSector and fZSectorErr
  const Int_t knbin = len;
  Int_t xbin =0;
  Int_t max=0;
  for(Int_t i = 0;i<knbin;i++){
    if(array[i]>max){
      max = array[i];
      xbin =i;
    }
  }
  Int_t hmax = max/2;
  Int_t xmin,xmax;
  Int_t ops = 0;
  xmin = xbin;
  while(xmin--){
    if(xmin<0) {ops++;break;}
    if(array[xmin]<hmax) {
      break;
    }
  }
  xmax = xbin;
  while(xmax++){
    if(xmax>=knbin) {ops++;break;}
    if(array[xmax]<hmax){
      break;
    }
  }
  if(ops){
    if(xbin >= knbin/2){xmin = 2 * xbin - knbin +1;xmax = knbin-1;}
    else{xmin = 0;xmax = 2 * xbin;}
  }
  Double_t sumw=0;
  Double_t sumw2=0;
  Double_t sumwx=0;
  Double_t sumwx2=0;
  for(Int_t bin = xmin;bin<=xmax;bin++){
    sumw   += array[bin];
    sumw2  += array[bin] * array[bin]; 
    sumwx  += array[bin] * vertex[bin];
    sumwx2 += array[bin] * vertex[bin] * vertex[bin];
  }
  if(sumw){
    Double_t mean = sumwx/sumw;
    Double_t rms2 = fabs(sumwx2/sumw - mean*mean);
    fZSectorErr = sqrt(rms2/sumw);
    fZSector = mean;
  }
  else{fZSectorErr = fZSector = 0;}
  sumw=sumw2=sumwx=sumwx2=0;
  for(Int_t bin = 0;bin<knbin;bin++){
    sumw   += array[bin];
    sumw2  += array[bin] * array[bin];
    sumwx  += array[bin] * vertex[bin];
    sumwx2 += array[bin] * vertex[bin] * vertex[bin];
  }
  if(sumw){
    Double_t mean = fZSector;
    Double_t rms2 = fabs(sumwx2/sumw - mean*mean);
    fZSectorErr = sqrt(rms2/sumw);
  }
}

