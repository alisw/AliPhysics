/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliLog.h"
#include "AliPHOSLogbackCluster.h"
#include "AliPHOSGeometry.h"
#include "AliVCluster.h"
#include "AliOADBContainer.h"
#include "AliVEvent.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMultiInputEventHandler.h"

ClassImp(AliPHOSLogbackCluster)

//______________________________________________________________________________
AliPHOSLogbackCluster::AliPHOSLogbackCluster(AliVCluster* cluster)
: fE( cluster->E() ), fCoreE(0)
{
  //TODO: fCoreE

  cluster->GetPosition( fPosition );
}

AliPHOSLogbackCluster::~AliPHOSLogbackCluster()
{
}


TLorentzVector AliPHOSLogbackCluster::GetMomentum(Double_t * vertex)
{
   // Returns TLorentzVector with momentum of the cluster. Only valid for clusters
  // identified as photons or pi0 (overlapped gamma) produced on the vertex
  // Vertex can be recovered with esd pointer doing:
  //" Double_t vertex[3] ; esd->GetVertex()->GetXYZ(vertex) ; "

  Double_t pos[3]={ fPosition[0] - vertex[0],
		    fPosition[1] - vertex[1],
		    fPosition[2] - vertex[2]  };

  Double_t r = TMath::Sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]   ) ;

  return TLorentzVector( fE*pos[0]/r,  fE*pos[1]/r,  fE*pos[2]/r,  fE) ;
}


AliPHOSGeometry* AliPHOSLogbackCluster::GetGeometry() const
{
  // Gets Geometry, from OADB

  static AliPHOSGeometry* sGeometry = 0x0;

  AliVEvent* vevent = AliPHOSLogbackCluster::GetCurrentEvent();
  Int_t runNumber = vevent->GetRunNumber();

  //Init geometry
  if(!sGeometry){
    AliOADBContainer geomContainer("phosGeo_AliPHOSLogbackCluster");
    geomContainer.InitFromFile("$ALICE_ROOT/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(runNumber,"PHOSRotationMatrixes");
    sGeometry =  AliPHOSGeometry::GetInstance("IHEP") ;
    for(Int_t mod=0; mod<5; mod++) {
      if(!matrixes->At(mod)) {
	// Please note that not all modules are present as of 2014
	// 	if( fDebug )
	// 	  AliInfo(Form("No PHOS Matrix for mod:%d, geo=%p\n", mod, sGeometry));
	AliDebug(2, Form("No PHOS Matrix for mod:%d, geo=%p\n", mod, sGeometry) );
	continue;
      }
      else {
	sGeometry->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
// 	if( fDebug >1 )
// 	  AliInfo(Form("Adding PHOS Matrix for mod:%d, geo=%p\n", mod, sGeometry));
      }
    }
  }

  return sGeometry;
}


AliVEvent* AliPHOSLogbackCluster::GetCurrentEvent() const
{
  // Hackish way of getting the current event.
  // Its probably not appropriate to call this function outside of
  // AliAnalysisTaskSE::UserExec

  AliAnalysisManager* analysisManager = dynamic_cast<AliAnalysisManager*>(AliAnalysisManager::GetAnalysisManager());
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(analysisManager->GetInputEventHandler());
  AliMultiInputEventHandler *multiInputHandler = dynamic_cast<AliMultiInputEventHandler *>(inputHandler);
  if (multiInputHandler)
    inputHandler = dynamic_cast<AliInputEventHandler *>(multiInputHandler->GetFirstInputEventHandler());

  AliVEvent* inputEvent = dynamic_cast<AliVEvent*>(inputHandler->GetEvent());
  if( ! inputEvent )
    AliError("Was not able to retrieve event!");

  return inputEvent;

}
