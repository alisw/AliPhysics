/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-------------------------------------------------------
//          Implementation of the TPC cluster debug information
//
//   Origin: Marian Ivanov   Marian.Ivanov@cern.ch
//   
//   Additional cluster information to monitor clustering performance
//   and to extract a features of detector response
//   Information attached to the default cluster
//   ONLY in DEBUG MODE
//    
//-------------------------------------------------------

/* $Id$ */

#include "AliTPCclusterInfo.h"
#include "AliLog.h"

ClassImp(AliTPCclusterInfo);


AliTPCclusterInfo::AliTPCclusterInfo():
  fNPads(0),
  fNTimeBins(0),
  fNBins(0),
  fGraph(0)
{
  //
  // default constructor
  //  
  for (Int_t i=0; i<25;i++){
    fMatrix[i] = i; 
  }
}

AliTPCclusterInfo::AliTPCclusterInfo(const  AliTPCclusterInfo & info):
  TObject(info),
  fNPads(info.fNPads),
  fNTimeBins(info.fNTimeBins),
  fNBins(info.fNBins),
  fGraph(0)    
{
  //
  // copy constructor
  //
  //  AliInfo("Copy constructor\n");
  for (Int_t i=0; i<25;i++){
    fMatrix[i] = info.fMatrix[i]; 
  }
  if (info.fGraph) fGraph = new Float_t[fNBins];
  for (Int_t i=0;i<fNBins; i++){
    fGraph[i] = info.fGraph[i];
  }
  
}


AliTPCclusterInfo::AliTPCclusterInfo(Bool_t extend):
  fNPads(0),
  fNTimeBins(0),
  fNBins(0),
  fGraph(0)
{
  //
  // allocate dummy graph - neccessary for IO part to use automatic branching
  //  
  for (Int_t i=0; i<25;i++){
    fMatrix[i] = i; 
  }
  if (extend){
    fNBins = 1;
    fGraph  = new Float_t[1];
    fGraph[0]=-1;
  }
}

AliTPCclusterInfo::AliTPCclusterInfo(Float_t *matrix, Int_t nbins, Float_t* graph){
  //
  // constructor of the info
  //
  for (Int_t i=0;i<25;i++){
    fMatrix[i]=matrix[i];
  }
  fNPads=0;
  fNTimeBins=0;
  Int_t center = 5+5+2;
  for (Int_t i=-2; i<=2;i++) if (matrix[center+i]>0) fNTimeBins++;
  for (Int_t i=-2; i<=2;i++) if (matrix[center+i*5]>0) fNPads++;
  fNBins = nbins;
  fGraph = 0;
  if (fNBins>0) {
    fGraph = new Float_t[fNBins];
    for (Int_t i=0;i<fNBins; i++){
      fGraph[i] = graph[i];
    }
  }
}

UChar_t AliTPCclusterInfo::GetNPads(Float_t threshold) const { 
  //
  //
  Int_t nPads=0;
  Int_t center = 5+5+2;
  for (Int_t i=-2; i<=2;i++) if (fMatrix[center+i*5]>threshold) nPads++;
  return nPads;
}

UChar_t AliTPCclusterInfo::GetNTimeBins(Float_t threshold) const { 
  //
  //
  //
  Int_t nTimeBins=0;
  Int_t center = 5+5+2;
  for (Int_t i=-2; i<=2;i++) if (fMatrix[center+i]>0) nTimeBins++;
  return nTimeBins;
}




AliTPCclusterInfo::~AliTPCclusterInfo(){
  //
  // destructor
  //
  if (fGraph)  delete [] fGraph;
}

