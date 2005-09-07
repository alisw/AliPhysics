// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCVertexArray.h"


/** \class AliHLTTPCVertexArray
<pre>
//_____________________________________________________________
// AliHLTTPCVertexArray
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

ClassImp(AliHLTTPCVertexArray)
void AliHLTTPCVertexArray::AnalyzeSector(Float_t *vertex, Int_t *array, Int_t len){
//loop over all seeds and all vertex position
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCVertexArray::AnalyzeSector","Analyze")
  <<AliHLTTPCLog::kDec<<"Number of Seeds: "<<fNSeed<<ENDLOG;
  for(Int_t i =0;i<fNSeed;i++)
    for(Int_t bin = 0;bin<len;bin++)
      array[bin] += Trace(fZSeed[i],fRSeed[i],fSecSeed[i],vertex[bin]);
}

void AliHLTTPCVertexArray::FindSectorVertex(Double_t pos, Double_t range, Int_t nbin){
//define position and range for search and
//loop over all seeds
//and find position of vertex and error 
  const Int_t len = nbin;
  const Double_t width = range; 
  const Double_t  xmin = pos - width/2;
  const Double_t step = width/len;
  const Double_t start = xmin + step/2.;
  Int_t * array = new Int_t[len];
  Float_t * ver   = new Float_t[len];
  for(Int_t i=0;i<len;i++){
    ver[i] =  start + step * i;
    array[i] = 0;
  }
  AnalyzeSector(ver,array,len);
  FindMean(ver,array,len);
  delete[] array;
  delete[] ver;
}

void AliHLTTPCVertexArray::FindMean(Float_t *vertex,Int_t *array, Int_t len){
//find mean and error of array and store it in
//fZSector and fZSectorErr
  const Int_t nbin = len;
  Int_t xbin =0;
  Int_t max=0;
  for(Int_t i = 0;i<nbin;i++){
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
    if(xmax>=nbin) {ops++;break;}
    if(array[xmax]<hmax){
      break;
    }
  }
  if(ops){
    if(xbin >= nbin/2){xmin = 2 * xbin - nbin +1;xmax = nbin-1;}
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
  for(Int_t bin = 0;bin<nbin;bin++){
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

