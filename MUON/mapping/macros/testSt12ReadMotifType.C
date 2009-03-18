// $Id$
// $MpId: testReadMotifType.C,v 1.12 2005/09/26 16:05:25 ivana Exp $
//
// Test macro for reading motif type data.

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpMotifReader.h"
#include "AliMpVMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifMap.h"
#include "AliMpConstants.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TH2.h>

#endif


void testReadMotifType(AliMq::Station12Type station, AliMp::PlaneType plane)
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpMotifReader r(dataStreams, AliMp::kStation12, station, plane);
  //r.SetVerboseLevel(2);

  TString names;
  TString names2;
  Int_t nv =0;
  if ( station == AliMq::kStation1 ) {
    if ( plane == AliMp::kBendingPlane ) 
      names ="ABCDEFGHI";
    else
      names = "ABCDEFGHIJKLMN";
  }    
  else if ( station == AliMq::kStation2 ) {
    if ( plane == AliMp::kBendingPlane ) {
      names ="ABCDEFGHIJKLMNOPQRSTUVWXY";
      names2 ="abcdefghimnptuv";
      nv = 5;
    }  
    else {
      names = "ABCEFGHIJKLMN";
      names2 ="abcdefgijklmnopqrstuwv";
      nv = 5;
    }  
  }  
    
  for (Int_t i=0;i<names.Length();++i){
     r.BuildMotifType(names[i])->Print("G");
  }
  
  // motifs a1, b1, ..., u1, v1
  for (Int_t i2=0;i2<names2.Length();++i2){
    TString mtName = names2[i2];
    mtName += "1";
    r.BuildMotifType(mtName)->Print("G");
  }
  
  // motifs v2, ..., v5, v6
  TString names4="v";
  for (Int_t i3=2;i3<nv+1;++i3) { 
      TString mtName = "v";
      mtName += i3;
      r.BuildMotifType(mtName)->Print("G");
  }
}

void testSt12ReadMotifType()
{
  AliMq::Station12Type  station[2] = { AliMq::kStation1, AliMq::kStation2 }; 
  AliMp::PlaneType      plane[2]   = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };
  
  for ( Int_t is = 0; is < 2; is++ ) {
    for ( Int_t ip = 0; ip < 2; ip++ ) {
    
      cout << "Running testReadMotifType for " 
           << AliMq::Station12TypeName(station[is]) << "  "
           << AliMp::PlaneTypeName(plane[ip])  << " ... " << endl;
       
      testReadMotifType(station[is], plane[ip]);
    
      cout << "... end running " << endl << endl;
    }  
  }   
}  
  
