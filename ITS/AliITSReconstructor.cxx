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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ITS reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliITSReconstructor.h"
#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliITSclustererV2.h"
#include "AliITStrackerMI.h"
#include "AliITStrackerSA.h"
#include "AliITSVertexerIons.h"
#include "AliITSVertexerFast.h"
#include "AliITSVertexerPPZ.h"
#include "AliITSVertexerZ.h"
#include "AliESD.h"
#include "AliITSpidESD.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliRun.h"
#include "AliITS.h"


ClassImp(AliITSReconstructor)


//_____________________________________________________________________________
void AliITSReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
// reconstruct clusters

  AliLoader* loader = runLoader->GetLoader("ITSLoader");
  if (!loader) {
    Error("Reconstruct", "ITS loader not found");
    return;
  }
  loader->LoadRecPoints("recreate");
  loader->LoadDigits("read");
  runLoader->LoadKinematics();

  AliITSgeom* geom = GetITSgeom(runLoader);
  if (!geom) return;
  AliITSclustererV2 clusterer(geom);
  Int_t nEvents = runLoader->GetNumberOfEvents();

  for (Int_t iEvent = 0; iEvent < nEvents; iEvent++) {
    runLoader->GetEvent(iEvent);

    TTree* treeClusters = loader->TreeR();
    if (!treeClusters) {
      loader->MakeTree("R");
      treeClusters = loader->TreeR();
    }
    TTree* treeDigits = loader->TreeD();
    if (!treeDigits) {
      Error("Reconstruct", "Can't get digits tree !");
      return;
    }

    clusterer.Digits2Clusters(treeDigits, treeClusters);
         
    loader->WriteRecPoints("OVERWRITE");
  }

  loader->UnloadRecPoints();
  loader->UnloadDigits();
  runLoader->UnloadKinematics();
}

//_____________________________________________________________________________
void AliITSReconstructor::Reconstruct(AliRunLoader* runLoader,
				      AliRawReader* rawReader) const
{
// reconstruct clusters from raw data

  AliLoader* loader = runLoader->GetLoader("ITSLoader");
  if (!loader) {
    Error("Reconstruct", "ITS loader not found");
    return;
  }
  loader->LoadRecPoints("recreate");

  AliITSgeom* geom = GetITSgeom(runLoader);
  if (!geom) return;
  AliITSclustererV2 clusterer(geom);

  Int_t iEvent = 0;
  while (rawReader->NextEvent()) {
    runLoader->GetEvent(iEvent++);

    TTree* treeClusters = loader->TreeR();
    if (!treeClusters) {
      loader->MakeTree("R");
      treeClusters = loader->TreeR();
    }

    clusterer.Digits2Clusters(rawReader);
         
    loader->WriteRecPoints("OVERWRITE");
  }

  loader->UnloadRecPoints();
}

//_____________________________________________________________________________
AliTracker* AliITSReconstructor::CreateTracker(AliRunLoader* runLoader) const
{
// create a ITS tracker

  AliITSgeom* geom = GetITSgeom(runLoader);
  if (!geom) return NULL; 
  TString selectedTracker = GetOption();
  if (selectedTracker.Contains("MI")) return new AliITStrackerMI(geom);
  return new AliITStrackerSA(geom);
}

//_____________________________________________________________________________
AliVertexer* AliITSReconstructor::CreateVertexer(AliRunLoader* /*runLoader*/) const
{
// create a ITS vertexer

  TString selectedVertexer = GetOption();
  if(selectedVertexer.Contains("ions") || selectedVertexer.Contains("IONS")){
    Info("CreateVertexer","a AliITSVertexerIons object has been selected\n");
    return new AliITSVertexerIons("null");
  }
  if(selectedVertexer.Contains("smear") || selectedVertexer.Contains("SMEAR")){
    Double_t smear[3]={0.005,0.005,0.01};
    Info("CreateVertexer","a AliITSVertexerFast object has been selected\n"); 
    return new AliITSVertexerFast(smear);
  }
  if(selectedVertexer.Contains("ppz") || selectedVertexer.Contains("PPZ")){
    Info("CreateVertexer","a AliITSVertexerPPZ object has been selected\n");
    return new AliITSVertexerPPZ("null");
  }
  // by default an AliITSVertexerZ object is instatiated
  Info("CreateVertexer","a AliITSVertexerZ object has been selected\n");
  return new AliITSVertexerZ("null");
}

//_____________________________________________________________________________
void AliITSReconstructor::FillESD(AliRunLoader* /*runLoader*/, 
				  AliESD* esd) const
{
// make PID, find V0s and cascades

  Double_t parITS[] = {34., 0.15, 10.};
  AliITSpidESD itsPID(parITS);
  itsPID.MakePID(esd);

  // V0 finding
  Double_t cuts[]={33,  // max. allowed chi2
		   0.16,// min. allowed negative daughter's impact parameter 
		   0.05,// min. allowed positive daughter's impact parameter 
		   0.08,// max. allowed DCA between the daughter tracks
		   0.99,// max. allowed cosine of V0's pointing angle
		   0.9,  // min. radius of the fiducial volume
		   2.9   // max. radius of the fiducial volume
  };
  AliV0vertexer vtxer(cuts);
  Double_t vtx[3], cvtx[6];
  esd->GetVertex()->GetXYZ(vtx);
  esd->GetVertex()->GetSigmaXYZ(cvtx);
  vtxer.SetVertex(vtx);
  vtxer.Tracks2V0vertices(esd);

  // cascade finding
  Double_t cts[]={33.,    // max. allowed chi2
		  0.05,   // min. allowed V0 impact parameter 
		  0.008,  // window around the Lambda mass 
		  0.035,  // min. allowed bachelor's impact parameter 
		  0.10,   // max. allowed DCA between a V0 and a track
		  0.9985, //max. allowed cosine of the cascade pointing angle
		  0.9,    // min. radius of the fiducial volume
		  2.9     // max. radius of the fiducial volume
  };
  AliCascadeVertexer cvtxer=AliCascadeVertexer(cts);
  cvtxer.SetVertex(vtx);
  cvtxer.V0sTracks2CascadeVertices(esd);
}


//_____________________________________________________________________________
AliITSgeom* AliITSReconstructor::GetITSgeom(AliRunLoader* runLoader) const
{
// get the ITS geometry

  if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
  if (!runLoader->GetAliRun()) {
    Error("GetITSgeom", "couldn't get AliRun object");
    return NULL;
  }
  AliITS* its = (AliITS*) runLoader->GetAliRun()->GetDetector("ITS");
  if (!its) {
    Error("GetITSgeom", "couldn't get ITS detector");
    return NULL;
  }
  if (!its->GetITSgeom()) {
    Error("GetITSgeom", "no ITS geometry available");
    return NULL;
  }
  return its->GetITSgeom();
}
