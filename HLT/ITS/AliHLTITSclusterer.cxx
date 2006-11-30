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

//-------------------------------------------------------------------------
//               Implementation of the HLT ITS clusterer class
//    The class derives from AliITSclustererV2.
//    There is one new method added which allows to read ITS raw data
//    and store the clusters in a tree without using runloaders.
//    In this case, the labels filling is skipped.
//          Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch
//-------------------------------------------------------------------------

#include "AliHLTITSclusterer.h"
#include "AliRawReader.h"
#include "AliITSgeom.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSRawStreamSSD.h"
#include <TTree.h>
#include <TClonesArray.h>

ClassImp(AliHLTITSclusterer)

AliHLTITSclusterer::AliHLTITSclusterer(const AliITSgeom *geom):AliITSclustererV2(geom)
{
fNModule = geom->GetIndexMax();
}

void AliHLTITSclusterer::Digits2Clusters(AliRawReader* rawReader,TTree *cTree)
{

  TClonesArray *array=new TClonesArray("AliITSclusterV2",1000);
  cTree->Branch("Clusters",&array);
  delete array;

  TClonesArray** clusters = new TClonesArray*[fNModule]; 
  for (Int_t iModule = 0; iModule < fNModule; iModule++) {
    clusters[iModule] = NULL;
  }

  rawReader->Reset();
  AliITSRawStreamSPD inputSPD(rawReader);
  FindClustersSPD(&inputSPD, clusters);

  rawReader->Reset();
  AliITSRawStreamSDD inputSDD(rawReader);
  FindClustersSDD(&inputSDD, clusters);

  rawReader->Reset();
  AliITSRawStreamSSD inputSSD(rawReader);
  FindClustersSSD(&inputSSD, clusters);

  // write all clusters to the tree
  Int_t nClusters = 0;
  for (Int_t iModule = 0; iModule < fNModule; iModule++) {
    array = clusters[iModule];
    if (!array) {
      Error("Digits2Clusters", "data for module %d missing!", iModule);
      array = new TClonesArray("AliITSclusterV2");
    }
    cTree->SetBranchAddress("Clusters", &array);
    cTree->Fill();
    nClusters += array->GetEntriesFast();
    delete array;
  }

  delete[] clusters;

  Info("Digits2Clusters", "total number of found clusters in ITS: %d\n", 
       nClusters);
}
