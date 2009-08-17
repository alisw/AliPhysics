// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Ingrid Kielen                                         *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTITSSSDQARecPointsComponent.cxx
    @author Ingrid Kielen
    @brief  Component for the SSD clusters QA
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTITSSSDQARecPointsComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliITSRecPoint.h"
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include "TObjString.h"
#include "TObjArray.h"
#include "AliITSgeomTGeo.h"

//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSSSDQARecPointsComponent)

AliHLTITSSSDQARecPointsComponent::AliHLTITSSSDQARecPointsComponent()
  :AliHLTProcessor(),
fHistArray(NULL),
fHistSSDModuleIdLayer5(NULL),
fHistSSDLocalXLayer5(NULL),
fHistSSDLocalZLayer5(NULL),
fHistSSDGlobalXLayer5(NULL),
fHistSSDGlobalYLayer5(NULL),
fHistSSDGlobalZLayer5(NULL),
fHistSSDPhiLayer5(NULL),
fHistSSDThetaLayer5(NULL),
fHistSSDRadiusLayer5(NULL),
fHistSSDClusterTypeLayer5(NULL),
fHistSSDChargeRatioLayer5(NULL),
fHistSSDChargekeVLayer5(NULL),
fHistSSDChargePSideLayer5(NULL),
fHistSSDChargeNSideLayer5(NULL),
fHistSSDChargeRatio2Layer5(NULL),
fHistSSDChargePNSideLayer5(NULL),
fHistSSDChargeMapLayer5(NULL),
fHistSSDClusterMapLayer5(NULL),
fHistSSDModuleIdLayer6(NULL),
fHistSSDLocalXLayer6(NULL),
fHistSSDLocalZLayer6(NULL),
fHistSSDGlobalXLayer6(NULL),
fHistSSDGlobalYLayer6(NULL),
fHistSSDGlobalZLayer6(NULL),
fHistSSDPhiLayer6(NULL),
fHistSSDThetaLayer6(NULL),
fHistSSDRadiusLayer6(NULL),
fHistSSDClusterTypeLayer6(NULL),
fHistSSDChargeRatioLayer6(NULL),
fHistSSDChargekeVLayer6(NULL),
fHistSSDChargePSideLayer6(NULL),
fHistSSDChargeNSideLayer6(NULL),
fHistSSDChargeRatio2Layer6(NULL),
fHistSSDChargePNSideLayer6(NULL),
fHistSSDChargeMapLayer6(NULL),
fHistSSDClusterMapLayer6(NULL) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTITSSSDQARecPointsComponent::~AliHLTITSSSDQARecPointsComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTITSSSDQARecPointsComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "ITSSSDQARecPoints";
}

void AliHLTITSSSDQARecPointsComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list) {
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD );
}

AliHLTComponentDataType AliHLTITSSSDQARecPointsComponent::GetOutputDataType() {
  // see header file for class documentation
  return kAliHLTDataTypeTObjArray;
}

void AliHLTITSSSDQARecPointsComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 10;
}

AliHLTComponent* AliHLTITSSSDQARecPointsComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTITSSSDQARecPointsComponent;
}

int AliHLTITSSSDQARecPointsComponent::DoInit(int argc, const char** argv) {
  //Inititalization of histograms for the SSD QA component

  Int_t nModuleOffset = 500;
  Int_t nITSTotalModules = AliITSgeomTGeo::GetNModules();  
  fHistSSDModuleIdLayer5 = new TH1F("SSD/Statistics/fHistSSDModuleIdLayer5",
				    "Module Id - Layer 5;Module Id;Entries",
				    fgkSSDMODULESLAYER5,
				    nModuleOffset - 0.5,
				    nITSTotalModules-fgkSSDMODULESLAYER6+0.5);
  
  fHistSSDModuleIdLayer6 = new TH1F("SSD/Statistics/fHistSSDModuleIdLayer6",
				    "Module Id - Layer 6;Module Id;Entries",
				    fgkSSDMODULESLAYER6,
				    nModuleOffset+fgkSSDMODULESLAYER5 - 0.5,
				    nITSTotalModules + 0.5);
  
  fHistSSDClusterPerEventLayer5 = new TH1F("SSD/Statistics/fHistSSDClusterPerEventLayer5",
					   "N_{clusters} - Layer 5;N_{clusters};Entries;",
					   100,0.1,5000);
  
  fHistSSDClusterPerEventLayer6 = new TH1F("SSD/Statistics/fHistSSDClusterPerEventLayer6",
					   "N_{clusters} - Layer 6;N_{clusters};Entries;",
					   100,0.1,5000);
  
  fHistSSDLocalXLayer5 = new TH1F("SSD/Coordinates/fHistSSDLocalXLayer5",
				  "Local x coord.- Layer 5;x [cm];Entries;",
				  100,-4.,4.);
  
  fHistSSDLocalXLayer6 = new TH1F("SSD/Coordinates/fHistSSDLocalXLayer6",
				  "Local x coord.- Layer 6;x [cm];Entries;",
				  100,-4.,4.);
  
  fHistSSDLocalZLayer5 = new TH1F("SSD/Coordinates/fHistSSDLocalZLayer5",
				  "Local z coord.- Layer 5;z [cm];Entries;",
				  100,-4.,4.);
  
  fHistSSDLocalZLayer6 = new TH1F("SSD/Coordinates/fHistSSDLocalZLayer6",
				  "Local z coord.- Layer 6;z [cm];Entries;",
				  100,-4.,4.);
  
  fHistSSDGlobalXLayer5 = new TH1F("SSD/Coordinates/fHistSSDGlobalXLayer5",
				   "Global x - Layer 5;x [cm];Entries;",
				   100,-40.,40.);
  
  fHistSSDGlobalXLayer6 = new TH1F("SSD/Coordinates/fHistSSDGlobalXLayer6",
				   "Global x - Layer 6;x [cm];Entries;",
				   100,-45.,45.);
  
  fHistSSDGlobalYLayer5 = new TH1F("SSD/Coordinates/fHistSSDGlobalYLayer5",
				   "Global y - Layer 5;y [cm];Entries;",
				   100,-40.,40);
  
  fHistSSDGlobalYLayer6 = new TH1F("SSD/Coordinates/fHistSSDGlobalYLayer6",
				   "Global y - Layer 6;y [cm];Entries;",
				   100,-45.,45.);
  
  fHistSSDGlobalZLayer5 = new TH1F("SSD/Coordinates/fHistSSDGlobalZLayer5",
				   "Global z - Layer 5;z [cm];Entries;",
				   100,-45.,45);
  
  fHistSSDGlobalZLayer6 = new TH1F("SSD/Coordinates/fHistSSDGlobalZLayer6",
				   "Global z - Layer 6;z [cm];Entries;",
				   100,-55.,55.);
  
  fHistSSDPhiLayer5 = new TH1F("SSD/Coordinates/fHistSSDPhiLayer5",
			       "#phi - Layer 5;#phi [rad];Entries;",
			       100,-TMath::Pi(),TMath::Pi());
  
  fHistSSDPhiLayer6 = new TH1F("SSD/Coordinates/fHistSSDPhiLayer6",
			       "#phi - Layer 6;#phi [rad];Entries;",
			       100,-TMath::Pi(),TMath::Pi());
  
  fHistSSDThetaLayer5 = new TH1F("SSD/Coordinates/fHistSSDThetaLayer5",
				 "#theta - Layer 5;#theta [rad];Entries;",
				 100,-TMath::Pi(),TMath::Pi());
  
  fHistSSDThetaLayer6 = new TH1F("SSD/Coordinates/fHistSSDThetaLayer6",
				 "#theta - Layer 6;#theta [rad];Entries;",
				 100,-TMath::Pi(),TMath::Pi());
  
  fHistSSDRadiusLayer5 = new TH1F("SSD/Coordinates/fHistSSDRadiusLayer5",
				  "r - Layer 5;r [cm];Entries;",
				  100,35.,50.);
  
  fHistSSDRadiusLayer6 = new TH1F("SSD/Coordinates/fHistSSDRadiusLayer6",
				  "r - Layer 6;r [cm];Entries;",
				  100,35.,50.);
  
  fHistSSDClusterTypeLayer5 = new TH1F("SSD/ClusterCharge/Layer5/fHistSSDClusterTypeLayer5",
				       "CL type - Layer 5;Cluster type;Entries;",
				       150,0,150);
  
  fHistSSDClusterTypeLayer6 = new TH1F("SSD/ClusterCharge/Layer6/fHistSSDClusterTypeLayer6",
				       "CL type - Layer 6;Cluster type;Entries;",
				       150,0,150);
  
  fHistSSDChargeRatioLayer5 = new TH1F("SSD/ClusterCharge/Layer5/fHistSSDChargeRatioLayer5",
				       "Charge ratio - Layer 5;q_{ratio};Entries;",
				       100,-2.0,2.0);
  
  fHistSSDChargeRatioLayer6 = new TH1F("SSD/ClusterCharge/Layer6/fHistSSDChargeRatioLayer6",
				       "Charge ratio - Layer 6;q_{ratio};Entries;",
				       100,-2.0,2.0);
  
  fHistSSDChargeRatio2Layer5 = new TH1F("fHistSSDChargeRatio2Layer5",
					"Charge Ratio qN/qP - Layer 5;q_{N}/q_{P};Entries;",
					100,0,2);
  fHistSSDChargeRatio2Layer6 = new TH1F("fHistSSDChargeRatio2Layer6",
					"Charge Ratio qN/qP - Layer 6;q_{N}/q_{P};Entries;",
					100,0,2);
  fHistSSDChargekeVLayer5 = new TH1F("SSD/ClusterCharge/Layer5/fHistSSDChargekeVLayer5",
				     "Charge - Layer 5;q [keV];Entries;",
				     100,0.,300.);
  
  fHistSSDChargekeVLayer6 = new TH1F("SSD/ClusterCharge/Layer6/fHistSSDChargekeVLayer6",
				     "Charge - Layer 6;q [keV];Entries;",
				     100,0.,300.);
  
  fHistSSDChargePSideLayer5 = new TH1F("SSD/ClusterCharge/Layer5/fHistSSDChargePSideLayer5",
				       "Charge P- Layer 5;q_{P} [keV];Entries;",
				       100,0.,300.);
  
  fHistSSDChargePSideLayer6 = new TH1F("SSD/ClusterCharge/Layer6/fHistSSDChargePSideLayer6",
				       "Charge P- Layer 6;q_{P} [keV];Entries;",
				       100,0.,300.);
  
  fHistSSDChargeNSideLayer5 = new TH1F("SSD/ClusterCharge/Layer5/fHistSSDChargeNSideLayer5",
				       "Charge N- Layer 5;q_{N} [keV];Entries;",
				       100,0.,300.);
  
  fHistSSDChargeNSideLayer6 = new TH1F("SSD/ClusterCharge/Layer6/fHistSSDChargeNSideLayer6",
				       "Charge N- Layer 6;q_{N} [keV];Entries;",
				       100,0.,300.);
  
  fHistSSDChargeRatio2Layer5 = new TH1F("SSD/ClusterCharge/Layer5/fHistSSDChargeRatio2Layer5",
					"Charge Ratio qN/qP - Layer 5;q_{N}/q_{P};Entries;",
					100,0,2);
  
  fHistSSDChargeRatio2Layer6 = new TH1F("SSD/ClusterCharge/Layer6/fHistSSDChargeRatio2Layer6",
					"Charge Ratio qN/qP - Layer 6;q_{N}/q_{P};Entries;",
					100,0,2);
  
  fHistSSDChargePNSideLayer5 = new TH2F("SSD/ClusterCharge/Layer5/fHistSSDChargePNSideLayer5",
					"Charge correlation - Layer 5;q_{P} [keV];q_{N} [keV]",
					100,0.,300.,
					100,0.,300.);
  
  fHistSSDChargePNSideLayer6 = new TH2F("SSD/ClusterCharge/Layer6/fHistSSDChargePNSideLayer6",
					"Charge correlation - Layer 6;q_{P} [keV];q_{N} [keV]",
					100,0.,300.,
					100,0.,300.);
  
  fHistSSDChargeMapLayer5 = new TH2F("SSD/ClusterCharge/Layer5/fHistSSDChargeMapLayer5",
				     "Charge map;N_{modules};N_{Ladders}",
				     fgkSSDMODULESPERLADDERLAYER5,
				     -0.5,fgkSSDMODULESPERLADDERLAYER5+0.5,
				     3*fgkSSDLADDERSLAYER5,
				     -0.5,fgkSSDLADDERSLAYER5+0.5);
  
  fHistSSDChargeMapLayer6 = new TH2F("SSD/ClusterCharge/Layer6/fHistSSDChargeMapLayer6",
				     "Charge map;N_{modules};N_{Ladders}",
				     fgkSSDMODULESPERLADDERLAYER6,
				     -0.5,fgkSSDMODULESPERLADDERLAYER6+0.5,
				     3*fgkSSDLADDERSLAYER6,
				     -0.5,fgkSSDLADDERSLAYER6+0.5);
  
  fHistSSDClusterMapLayer5 = new TH2F("SSD/Statistics/Layer5/fHistSSDClusterMapLayer5",
				      "Layer 5;N_{module};N_{ladder}",
				      22,1,23,
				      34,500,534);
  fHistSSDClusterMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistSSDClusterMapLayer5->SetStats(kFALSE);
  fHistSSDClusterMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistSSDClusterMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistSSDClusterMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistSSDClusterMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistSSDClusterMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistSSDClusterMapLayer5->GetZaxis()->SetTitleOffset(1.4);
  fHistSSDClusterMapLayer5->GetZaxis()->SetTitle("N_{clusters}");
  
  fHistSSDClusterMapLayer6 = new TH2F("SSD/Statistics/Layer6/fHistSSDClusterMapLayer6",
				      "Layer 6;N_{module};N_{ladder}",
				      25,1,26,
				      38,600,638);
  fHistSSDClusterMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistSSDClusterMapLayer6->SetStats(kFALSE);
  fHistSSDClusterMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistSSDClusterMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistSSDClusterMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistSSDClusterMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistSSDClusterMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistSSDClusterMapLayer6->GetZaxis()->SetTitleOffset(1.4);
  fHistSSDClusterMapLayer6->GetZaxis()->SetTitle("N_{clusters}");

  HLTInfo("Finished initialization of SSD histograms");  

  return 0; 
 
}
  
int AliHLTITSSSDQARecPointsComponent::DoDeinit() {
  // see header file for class documentation
  if(fHistSSDModuleIdLayer5) delete fHistSSDModuleIdLayer5;
  if(fHistSSDModuleIdLayer6) delete fHistSSDModuleIdLayer6;    
  if(fHistSSDClusterPerEventLayer5) delete fHistSSDClusterPerEventLayer5;
  if(fHistSSDClusterPerEventLayer6) delete fHistSSDClusterPerEventLayer6;   
  if(fHistSSDLocalXLayer5) delete fHistSSDLocalXLayer5;
  if(fHistSSDLocalXLayer6) delete fHistSSDLocalXLayer6;   
  if(fHistSSDLocalZLayer5) delete fHistSSDLocalZLayer5;
  if(fHistSSDLocalZLayer6) delete fHistSSDLocalZLayer6;   
  if(fHistSSDGlobalXLayer5) delete fHistSSDGlobalXLayer5;
  if(fHistSSDGlobalXLayer6) delete fHistSSDGlobalXLayer6;   
  if(fHistSSDGlobalYLayer5) delete fHistSSDGlobalYLayer5;
  if(fHistSSDGlobalYLayer6) delete fHistSSDGlobalYLayer6;   
  if(fHistSSDGlobalZLayer5) delete fHistSSDGlobalZLayer5;
  if(fHistSSDGlobalZLayer6) delete fHistSSDGlobalZLayer6;   
  if(fHistSSDPhiLayer5) delete fHistSSDPhiLayer5;
  if(fHistSSDPhiLayer6) delete fHistSSDPhiLayer6;      
  if(fHistSSDThetaLayer5) delete fHistSSDThetaLayer5;
  if(fHistSSDThetaLayer6) delete fHistSSDThetaLayer6;      
  if(fHistSSDRadiusLayer5) delete fHistSSDRadiusLayer5;
  if(fHistSSDRadiusLayer6) delete fHistSSDRadiusLayer6;         
  if(fHistSSDClusterTypeLayer5) delete fHistSSDClusterTypeLayer5;
  if(fHistSSDClusterTypeLayer6) delete fHistSSDClusterTypeLayer6;         
  if(fHistSSDChargeRatioLayer5) delete fHistSSDChargeRatioLayer5;
  if(fHistSSDChargeRatioLayer6) delete fHistSSDChargeRatioLayer6;         
  if(fHistSSDChargekeVLayer5) delete fHistSSDChargekeVLayer5;
  if(fHistSSDChargekeVLayer6) delete fHistSSDChargekeVLayer6;         
  if(fHistSSDChargePSideLayer5) delete fHistSSDChargePSideLayer5;
  if(fHistSSDChargePSideLayer6) delete fHistSSDChargePSideLayer6;         
  if(fHistSSDChargeNSideLayer5) delete fHistSSDChargeNSideLayer5;
  if(fHistSSDChargeNSideLayer6) delete fHistSSDChargeNSideLayer6;         
  if(fHistSSDChargeRatio2Layer5) delete fHistSSDChargeRatio2Layer5;
  if(fHistSSDChargeRatio2Layer6) delete fHistSSDChargeRatio2Layer6;               
  if(fHistSSDChargePNSideLayer5) delete fHistSSDChargePNSideLayer5;
  if(fHistSSDChargePNSideLayer6) delete fHistSSDChargePNSideLayer6;         
  if(fHistSSDChargeMapLayer5) delete fHistSSDChargeMapLayer5;
  if(fHistSSDChargeMapLayer6) delete fHistSSDChargeMapLayer6;         
  if(fHistSSDClusterMapLayer5) delete fHistSSDClusterMapLayer5;
  if(fHistSSDClusterMapLayer6) delete fHistSSDClusterMapLayer6; 
  return 0;
}

int AliHLTITSSSDQARecPointsComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  
  int TotalSpacePoint = 0;
    
  const AliHLTComponentBlockData* iter = NULL;
  
  if(!IsDataEvent())
    return 0;
   
 for ( iter = GetFirstInputBlock(kAliHLTDataTypeClusters); iter != NULL; iter = GetNextInputBlock() ) {
    
    if(iter->fDataType!=(kAliHLTAnyDataType|kAliHLTDataOriginITSSSD))
       continue;

    const AliHLTITSClusterData* clusterData = (const AliHLTITSClusterData*) iter->fPtr;
    Int_t nSpacepoint = (Int_t) clusterData->fSpacePointCnt;
    TotalSpacePoint += nSpacepoint;
    AliHLTITSSpacePointData *clusters = (AliHLTITSSpacePointData*) clusterData->fSpacePoints;

    for(int i=0;i<nSpacepoint;i++){ //cluster loop
      Int_t lab[4]={0,0,0,0};
      Float_t hit[6]={0,0,0,0,0,0};
      Int_t info[3]={0,0,0};
      
      lab[0]=clusters[i].fTracks[0];
      lab[1]=clusters[i].fTracks[1];
      lab[2]=clusters[i].fTracks[2];
      lab[3]=clusters[i].fIndex;
      hit[0]=clusters[i].fY;
      hit[1]=clusters[i].fZ;
      hit[2]=clusters[i].fSigmaY2;
      hit[3]=clusters[i].fSigmaZ2;
      hit[4]=clusters[i].fQ;
      hit[5]=clusters[i].fSigmaYZ;
      info[0]=clusters[i].fNy;
      info[1]=clusters[i].fNz;
      info[2]=clusters[i].fLayer;

      AliITSRecPoint recp(lab,hit,info);
     
     Int_t module = 0;
     Int_t gLayer = 0, gLadder = 0, gModule = 0;
     Int_t lLadderLocationY = 0;
      
     Int_t nClustersLayer5 = 0, nClustersLayer6 = 0;   
     Float_t cluglo[3]={0.,0.,0.}; 
     
     Int_t layer = recp.GetLayer();
     if (layer == 4) module = recp.GetDetectorIndex() + 500;
     if (layer == 5) module = recp.GetDetectorIndex() + 1248;

     AliITSgeomTGeo::GetModuleId(module,gLayer,gLadder,gModule);
     lLadderLocationY = 3*gLadder; 

     recp.GetGlobalXYZ(cluglo);
     Float_t radius = TMath::Sqrt(cluglo[0]*cluglo[0]+cluglo[1]*cluglo[1]); 
     Float_t phi = TMath::ATan2(cluglo[1],cluglo[0]);
     Float_t theta = TMath::ATan2(radius,cluglo[2]);
     Double_t chargeRatio = recp.GetChargeRatio();
     Double_t clusterCharge = recp.GetQ();
     Double_t chargePSide = clusterCharge*(1. + chargeRatio);
     Double_t chargeNSide = clusterCharge*(1. - chargeRatio);
     
      if(layer == 4) {
	fHistSSDModuleIdLayer5->Fill(module);
	fHistSSDLocalXLayer5->Fill(recp.GetDetLocalX());
	fHistSSDLocalZLayer5->Fill(recp.GetDetLocalZ());
	fHistSSDGlobalXLayer5->Fill(cluglo[0]);
	fHistSSDGlobalYLayer5->Fill(cluglo[1]);
	fHistSSDGlobalZLayer5->Fill(cluglo[2]);
	fHistSSDPhiLayer5->Fill(phi);
	fHistSSDThetaLayer5->Fill(theta);
	fHistSSDRadiusLayer5->Fill(radius);
	fHistSSDClusterTypeLayer5->Fill(recp.GetType());
	fHistSSDChargeRatioLayer5->Fill(recp.GetChargeRatio());
	fHistSSDChargekeVLayer5->Fill(recp.GetQ());
	fHistSSDChargePSideLayer5->Fill(chargePSide);
	fHistSSDChargeNSideLayer5->Fill(chargeNSide);
	if(chargePSide != 0.) fHistSSDChargeRatio2Layer5->Fill(chargeNSide/chargePSide);
	fHistSSDChargePNSideLayer5->Fill(chargePSide,chargeNSide);
	fHistSSDChargeMapLayer5->SetBinContent(gModule,lLadderLocationY,recp.GetQ());
	((TH2F *)fHistSSDClusterMapLayer5)->Fill(gModule,499+gLadder,1);
	nClustersLayer5 += 1;
      }//layer 5 histograms
      if(layer == 5) {
	fHistSSDModuleIdLayer6->Fill(module);
	fHistSSDLocalXLayer6->Fill(recp.GetDetLocalX());
	fHistSSDLocalZLayer6->Fill(recp.GetDetLocalZ());
	fHistSSDGlobalXLayer6->Fill(cluglo[0]);
	fHistSSDGlobalYLayer6->Fill(cluglo[1]);
	fHistSSDGlobalZLayer6->Fill(cluglo[2]);
	fHistSSDPhiLayer6->Fill(phi);
	fHistSSDThetaLayer6->Fill(theta);
	fHistSSDRadiusLayer6->Fill(radius);
	fHistSSDClusterTypeLayer6->Fill(recp.GetType());
	fHistSSDChargeRatioLayer6->Fill(recp.GetChargeRatio());
	fHistSSDChargekeVLayer6->Fill(recp.GetQ());
	fHistSSDChargePSideLayer6->Fill(chargePSide);
	fHistSSDChargeNSideLayer6->Fill(chargeNSide);
	if(chargePSide != 0.) fHistSSDChargeRatio2Layer6->Fill(chargeNSide/chargePSide);
	fHistSSDChargePNSideLayer6->Fill(chargePSide,chargeNSide);
	fHistSSDChargeMapLayer6->SetBinContent(gModule,lLadderLocationY,recp.GetQ());
	((TH2F *)fHistSSDClusterMapLayer6)->Fill(gModule,599+gLadder,1);
	nClustersLayer6 += 1;
      }//layer 6 histograms
    }//cluster loop
  }//input block
      
    //Write all histograms in TObjArray
    fHistArray = new TObjArray; 
    fHistArray->AddLast(fHistSSDModuleIdLayer5); 
    fHistArray->AddLast(fHistSSDLocalXLayer5); 
    fHistArray->AddLast(fHistSSDLocalZLayer5); 
    fHistArray->AddLast(fHistSSDGlobalXLayer5); 
    fHistArray->AddLast(fHistSSDGlobalYLayer5); 
    fHistArray->AddLast(fHistSSDGlobalZLayer5); 
    fHistArray->AddLast(fHistSSDPhiLayer5); 
    fHistArray->AddLast(fHistSSDThetaLayer5); 
    fHistArray->AddLast(fHistSSDRadiusLayer5); 
    fHistArray->AddLast(fHistSSDClusterTypeLayer5); 
    fHistArray->AddLast(fHistSSDChargeRatioLayer5); 
    fHistArray->AddLast(fHistSSDChargekeVLayer5); 
    fHistArray->AddLast(fHistSSDChargePSideLayer5); 
    fHistArray->AddLast(fHistSSDChargeNSideLayer5); 
    fHistArray->AddLast(fHistSSDChargeRatio2Layer5); 
    fHistArray->AddLast(fHistSSDChargePNSideLayer5); 
    fHistArray->AddLast(fHistSSDChargeMapLayer5); 
    fHistArray->AddLast(fHistSSDClusterMapLayer5);  
    fHistArray->AddLast(fHistSSDModuleIdLayer6); 
    fHistArray->AddLast(fHistSSDLocalXLayer6); 
    fHistArray->AddLast(fHistSSDLocalZLayer6); 
    fHistArray->AddLast(fHistSSDGlobalXLayer6); 
    fHistArray->AddLast(fHistSSDGlobalYLayer6); 
    fHistArray->AddLast(fHistSSDGlobalZLayer6); 
    fHistArray->AddLast(fHistSSDPhiLayer6); 
    fHistArray->AddLast(fHistSSDThetaLayer6); 
    fHistArray->AddLast(fHistSSDRadiusLayer6); 
    fHistArray->AddLast(fHistSSDClusterTypeLayer6); 
    fHistArray->AddLast(fHistSSDChargeRatioLayer6); 
    fHistArray->AddLast(fHistSSDChargekeVLayer6); 
    fHistArray->AddLast(fHistSSDChargePSideLayer6); 
    fHistArray->AddLast(fHistSSDChargeNSideLayer6); 
    fHistArray->AddLast(fHistSSDChargeRatio2Layer6); 
    fHistArray->AddLast(fHistSSDChargePNSideLayer6); 
    fHistArray->AddLast(fHistSSDChargeMapLayer6); 
    fHistArray->AddLast(fHistSSDClusterMapLayer6); 

    //Publish array
    AliHLTUInt32_t fSpecification = 0x0;
    PushBack( (TObjArray*) fHistArray,kAliHLTDataTypeTObjArray,fSpecification);
  return 0;
}
