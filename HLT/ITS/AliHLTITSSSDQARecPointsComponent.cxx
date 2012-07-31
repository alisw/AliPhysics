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
    @author Ingrid Kielen - Panos Christakoglou (Panos.Christakoglou@cern.ch)
    @brief  Component for the SSD clusters QA
*/

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
#include "TH1.h"
#include "TH2.h"
#include "AliITSgeomTGeo.h"
#include "AliGeomManager.h"

//#include <stdlib.h>
//#include <cerrno>

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTITSSSDQARecPointsComponent)

AliHLTITSSSDQARecPointsComponent::AliHLTITSSSDQARecPointsComponent()
  :AliHLTProcessor(),
   fHistSSDArray(NULL) {
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

int AliHLTITSSSDQARecPointsComponent::DoInit(int /*argc*/, const char** /*argv*/) {
  //Inititalization of histograms for the SSD QA component
  
  
  if(AliGeomManager::GetGeometry()==NULL){
     AliGeomManager::LoadGeometry();
  }
  
  
  fHistSSDArray = new TObjArray();  
  fHistSSDArray->SetName("ssdArray");

  Int_t nModuleOffset = 500;
  Int_t nITSTotalModules = AliITSgeomTGeo::GetNModules();  
  TH1F *gHistSSDModuleIdLayer5 = new TH1F("fHistSSDModuleIdLayer5",
					  "Module Id - Layer 5;Module Id;Entries",
					  fgkSSDMODULESLAYER5,
					  nModuleOffset - 0.5,
					  nITSTotalModules-fgkSSDMODULESLAYER6+0.5);
  fHistSSDArray->Add(gHistSSDModuleIdLayer5);//entry 0

  TH1F *gHistSSDModuleIdLayer6 = new TH1F("fHistSSDModuleIdLayer6",
					  "Module Id - Layer 6;Module Id;Entries",
					  fgkSSDMODULESLAYER6,
					  nModuleOffset+fgkSSDMODULESLAYER5 - 0.5,
					  nITSTotalModules + 0.5);
  fHistSSDArray->Add(gHistSSDModuleIdLayer6);//entry 1
  
  TH1F *gHistSSDClusterPerEventLayer5 = new TH1F("fHistSSDClusterPerEventLayer5",
						 "N_{clusters} - Layer 5;N_{clusters};Entries;",
						 100,0.1,5000);
  fHistSSDArray->Add(gHistSSDClusterPerEventLayer5);//entry 2
  
  TH1F *gHistSSDClusterPerEventLayer6 = new TH1F("fHistSSDClusterPerEventLayer6",
						 "N_{clusters} - Layer 6;N_{clusters};Entries;",
						 100,0.1,5000);
  fHistSSDArray->Add(gHistSSDClusterPerEventLayer6);//entry 3
  
  TH1F *gHistSSDLocalXLayer5 = new TH1F("fHistSSDLocalXLayer5",
					"Local x coord.- Layer 5;x [cm];Entries;",
					100,-4.,4.);
  fHistSSDArray->Add(gHistSSDLocalXLayer5);//entry 4
  
  TH1F *gHistSSDLocalXLayer6 = new TH1F("fHistSSDLocalXLayer6",
					"Local x coord.- Layer 6;x [cm];Entries;",
					100,-4.,4.);
  fHistSSDArray->Add(gHistSSDLocalXLayer6);//entry 5
  
  TH1F *gHistSSDLocalZLayer5 = new TH1F("fHistSSDLocalZLayer5",
					"Local z coord.- Layer 5;z [cm];Entries;",
					100,-4.,4.);
  fHistSSDArray->Add(gHistSSDLocalZLayer5);//entry 6
  
  TH1F *gHistSSDLocalZLayer6 = new TH1F("fHistSSDLocalZLayer6",
					"Local z coord.- Layer 6;z [cm];Entries;",
					100,-4.,4.);
  fHistSSDArray->Add(gHistSSDLocalZLayer6);//entry 7
  
  TH1F *gHistSSDGlobalXLayer5 = new TH1F("fHistSSDGlobalXLayer5",
					 "Global x - Layer 5;x [cm];Entries;",
					 100,-40.,40.);
  fHistSSDArray->Add(gHistSSDGlobalXLayer5);//entry 8
  
  TH1F *gHistSSDGlobalXLayer6 = new TH1F("fHistSSDGlobalXLayer6",
					 "Global x - Layer 6;x [cm];Entries;",
					 100,-45.,45.);
  fHistSSDArray->Add(gHistSSDGlobalXLayer6);//entry 9
  
  TH1F *gHistSSDGlobalYLayer5 = new TH1F("fHistSSDGlobalYLayer5",
					 "Global y - Layer 5;y [cm];Entries;",
					 100,-40.,40);
  fHistSSDArray->Add(gHistSSDGlobalYLayer5);//entry 10
  
  TH1F *gHistSSDGlobalYLayer6 = new TH1F("fHistSSDGlobalYLayer6",
					 "Global y - Layer 6;y [cm];Entries;",
					 100,-45.,45.);
  fHistSSDArray->Add(gHistSSDGlobalYLayer6);//entry 11
  
  TH1F *gHistSSDGlobalZLayer5 = new TH1F("fHistSSDGlobalZLayer5",
					 "Global z - Layer 5;z [cm];Entries;",
					 100,-45.,45);
  fHistSSDArray->Add(gHistSSDGlobalZLayer5);//entry 12
  
  TH1F *gHistSSDGlobalZLayer6 = new TH1F("fHistSSDGlobalZLayer6",
					 "Global z - Layer 6;z [cm];Entries;",
					 100,-55.,55.);
  fHistSSDArray->Add(gHistSSDGlobalZLayer6);//entry 13
  
  TH1F *gHistSSDPhiLayer5 = new TH1F("fHistSSDPhiLayer5",
				     "#phi - Layer 5;#phi [rad];Entries;",
				     100,-TMath::Pi(),TMath::Pi());
  fHistSSDArray->Add(gHistSSDPhiLayer5);//entry 14
  
  TH1F *gHistSSDPhiLayer6 = new TH1F("fHistSSDPhiLayer6",
				     "#phi - Layer 6;#phi [rad];Entries;",
				     100,-TMath::Pi(),TMath::Pi());
  fHistSSDArray->Add(gHistSSDPhiLayer6);//entry 15
  
  TH1F *gHistSSDThetaLayer5 = new TH1F("fHistSSDThetaLayer5",
				       "#theta - Layer 5;#theta [rad];Entries;",
				       100,-TMath::Pi(),TMath::Pi());
  fHistSSDArray->Add(gHistSSDThetaLayer5);//entry 16
  
  TH1F *gHistSSDThetaLayer6 = new TH1F("fHistSSDThetaLayer6",
				       "#theta - Layer 6;#theta [rad];Entries;",
				       100,-TMath::Pi(),TMath::Pi());
  fHistSSDArray->Add(gHistSSDThetaLayer6);//entry 17
  
  TH1F *gHistSSDRadiusLayer5 = new TH1F("fHistSSDRadiusLayer5",
					"r - Layer 5;r [cm];Entries;",
					100,35.,50.);
  fHistSSDArray->Add(gHistSSDRadiusLayer5);//entry 18
  
  TH1F *gHistSSDRadiusLayer6 = new TH1F("fHistSSDRadiusLayer6",
					"r - Layer 6;r [cm];Entries;",
					100,35.,50.);
  fHistSSDArray->Add(gHistSSDRadiusLayer6);//entry 19
  
  TH1F *gHistSSDClusterTypeLayer5 = new TH1F("Layer5/fHistSSDClusterTypeLayer5",
					     "CL type - Layer 5;Cluster type;Entries;",
					     150,0,150);
  fHistSSDArray->Add(gHistSSDClusterTypeLayer5);//entry 20
  
  TH1F *gHistSSDClusterTypeLayer6 = new TH1F("fHistSSDClusterTypeLayer6",
					     "CL type - Layer 6;Cluster type;Entries;",
					     150,0,150);
  fHistSSDArray->Add(gHistSSDClusterTypeLayer6);//entry 21
  
  TH1F *gHistSSDChargeRatioLayer5 = new TH1F("fHistSSDChargeRatioLayer5",
					     "Charge ratio - Layer 5;q_{ratio};Entries;",
					     100,-2.0,2.0);
  fHistSSDArray->Add(gHistSSDChargeRatioLayer5);//entry 22
  
  TH1F *gHistSSDChargeRatioLayer6 = new TH1F("fHistSSDChargeRatioLayer6",
					     "Charge ratio - Layer 6;q_{ratio};Entries;",
					     100,-2.0,2.0);
  fHistSSDArray->Add(gHistSSDChargeRatioLayer6);//entry 23
  
  TH1F *gHistSSDChargekeVLayer5 = new TH1F("fHistSSDChargekeVLayer5",
					   "Charge - Layer 5;q [keV];Entries;",
					   100,0.,300.);
  fHistSSDArray->Add(gHistSSDChargekeVLayer5);//entry 24
  
  TH1F *gHistSSDChargekeVLayer6 = new TH1F("fHistSSDChargekeVLayer6",
					   "Charge - Layer 6;q [keV];Entries;",
					   100,0.,300.);
  fHistSSDArray->Add(gHistSSDChargekeVLayer6);//entry 25
  
  TH1F *gHistSSDChargePSideLayer5 = new TH1F("fHistSSDChargePSideLayer5",
					     "Charge P- Layer 5;q_{P} [keV];Entries;",
					     100,0.,300.);
  fHistSSDArray->Add(gHistSSDChargePSideLayer5);//entry 26
  
  TH1F *gHistSSDChargePSideLayer6 = new TH1F("fHistSSDChargePSideLayer6",
					     "Charge P- Layer 6;q_{P} [keV];Entries;",
					     100,0.,300.);
  fHistSSDArray->Add(gHistSSDChargePSideLayer6);//entry 27
  
  TH1F *gHistSSDChargeNSideLayer5 = new TH1F("fHistSSDChargeNSideLayer5",
					     "Charge N- Layer 5;q_{N} [keV];Entries;",
					     100,0.,300.);
  fHistSSDArray->Add(gHistSSDChargeNSideLayer5);//entry 28
  
  TH1F *gHistSSDChargeNSideLayer6 = new TH1F("fHistSSDChargeNSideLayer6",
					     "Charge N- Layer 6;q_{N} [keV];Entries;",
					     100,0.,300.);
  fHistSSDArray->Add(gHistSSDChargeNSideLayer6);//entry 29
  
  TH1F *gHistSSDChargeRatio2Layer5 = new TH1F("fHistSSDChargeRatio2Layer5",
					      "Charge Ratio qN/qP - Layer 5;q_{N}/q_{P};Entries;",
					      100,0,2);
  fHistSSDArray->Add(gHistSSDChargeRatio2Layer5);//entry 30
  
  TH1F *gHistSSDChargeRatio2Layer6 = new TH1F("fHistSSDChargeRatio2Layer6",
					      "Charge Ratio qN/qP - Layer 6;q_{N}/q_{P};Entries;",
					      100,0,2);
  fHistSSDArray->Add(gHistSSDChargeRatio2Layer6);//entry 31
  
  TH2F *gHistSSDChargePNSideLayer5 = new TH2F("fHistSSDChargePNSideLayer5",
					      "Charge correlation - Layer 5;q_{P} [keV];q_{N} [keV]",
					      100,0.,300.,
					      100,0.,300.);
  fHistSSDArray->Add(gHistSSDChargePNSideLayer5);//entry 32
  
  TH2F *gHistSSDChargePNSideLayer6 = new TH2F("fHistSSDChargePNSideLayer6",
					      "Charge correlation - Layer 6;q_{P} [keV];q_{N} [keV]",
					      100,0.,300.,
					      100,0.,300.);
  fHistSSDArray->Add(gHistSSDChargePNSideLayer6);//entry 33
  
  TH2F *gHistSSDChargeMapLayer5 = new TH2F("fHistSSDChargeMapLayer5",
					   "Charge map;N_{modules};N_{Ladders}",
					   fgkSSDMODULESPERLADDERLAYER5,
					   -0.5,fgkSSDMODULESPERLADDERLAYER5+0.5,
					   3*fgkSSDLADDERSLAYER5,
					   -0.5,fgkSSDLADDERSLAYER5+0.5);
  fHistSSDArray->Add(gHistSSDChargeMapLayer5);//entry 34
  
  TH2F *gHistSSDChargeMapLayer6 = new TH2F("fHistSSDChargeMapLayer6",
					   "Charge map;N_{modules};N_{Ladders}",
					   fgkSSDMODULESPERLADDERLAYER6,
					   -0.5,fgkSSDMODULESPERLADDERLAYER6+0.5,
					   3*fgkSSDLADDERSLAYER6,
					   -0.5,fgkSSDLADDERSLAYER6+0.5);
  fHistSSDArray->Add(gHistSSDChargeMapLayer6);//entry 35
  
  TH2F *gHistSSDClusterMapLayer5 = new TH2F("fHistSSDClusterMapLayer5",
					    "Layer 5;N_{module};N_{ladder}",
					    22,1,23,
					    34,500,534);
  gHistSSDClusterMapLayer5->GetXaxis()->SetTitleColor(1);
  gHistSSDClusterMapLayer5->SetStats(kFALSE);
  gHistSSDClusterMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  gHistSSDClusterMapLayer5->GetXaxis()->SetNdivisions(22);
  gHistSSDClusterMapLayer5->GetYaxis()->SetNdivisions(34);
  gHistSSDClusterMapLayer5->GetXaxis()->SetLabelSize(0.03);
  gHistSSDClusterMapLayer5->GetYaxis()->SetLabelSize(0.03);
  gHistSSDClusterMapLayer5->GetZaxis()->SetTitleOffset(1.4);
  gHistSSDClusterMapLayer5->GetZaxis()->SetTitle("N_{clusters}");
  fHistSSDArray->Add(gHistSSDClusterMapLayer5); //entry 36
  
  TH2F *gHistSSDClusterMapLayer6 = new TH2F("fHistSSDClusterMapLayer6",
					    "Layer 6;N_{module};N_{ladder}",
					    25,1,26,
					    38,600,638);
  gHistSSDClusterMapLayer6->GetXaxis()->SetTitleColor(1);
  gHistSSDClusterMapLayer6->SetStats(kFALSE);
  gHistSSDClusterMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  gHistSSDClusterMapLayer6->GetXaxis()->SetNdivisions(25);
  gHistSSDClusterMapLayer6->GetYaxis()->SetNdivisions(38);
  gHistSSDClusterMapLayer6->GetXaxis()->SetLabelSize(0.03);
  gHistSSDClusterMapLayer6->GetYaxis()->SetLabelSize(0.03);
  gHistSSDClusterMapLayer6->GetZaxis()->SetTitleOffset(1.4);
  gHistSSDClusterMapLayer6->GetZaxis()->SetTitle("N_{clusters}");
  fHistSSDArray->Add(gHistSSDClusterMapLayer6);//entry 37

  HLTInfo("Finished initialization of SSD histograms");  

  return 0;  
}
  
int AliHLTITSSSDQARecPointsComponent::DoDeinit() {
  // see header file for class documentation
  if(fHistSSDArray) delete fHistSSDArray;

  return 0;
}

int AliHLTITSSSDQARecPointsComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  
  Int_t nClustersLayer5 = 0, nClustersLayer6 = 0;   
  
  const AliHLTComponentBlockData* iter = NULL;
  
  if(!IsDataEvent())
    return 0;
  
  for ( iter = GetFirstInputBlock(kAliHLTDataTypeClusters); iter != NULL; iter = GetNextInputBlock() ) {
    
    if(iter->fDataType!=(kAliHLTAnyDataType|kAliHLTDataOriginITSSSD))
      continue;
    
    const AliHLTITSClusterData* clusterData = (const AliHLTITSClusterData*) iter->fPtr;
    Int_t nSpacepoint = (Int_t) clusterData->fSpacePointCnt;
    AliHLTITSSpacePointData *clusters = (AliHLTITSSpacePointData*) clusterData->fSpacePoints;
    
    for(int i = 0; i < nSpacepoint; i++) { //cluster loop
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
	dynamic_cast<TH1F *>(fHistSSDArray->At(0))->Fill(module);
	dynamic_cast<TH1F *>(fHistSSDArray->At(4))->Fill(recp.GetDetLocalX());
	dynamic_cast<TH1F *>(fHistSSDArray->At(6))->Fill(recp.GetDetLocalZ());
	dynamic_cast<TH1F *>(fHistSSDArray->At(8))->Fill(cluglo[0]);
	dynamic_cast<TH1F *>(fHistSSDArray->At(10))->Fill(cluglo[1]);
	dynamic_cast<TH1F *>(fHistSSDArray->At(12))->Fill(cluglo[2]);
	dynamic_cast<TH1F *>(fHistSSDArray->At(14))->Fill(phi);
	dynamic_cast<TH1F *>(fHistSSDArray->At(16))->Fill(theta);
	dynamic_cast<TH1F *>(fHistSSDArray->At(18))->Fill(radius);
	dynamic_cast<TH1F *>(fHistSSDArray->At(20))->Fill(recp.GetType());
	dynamic_cast<TH1F *>(fHistSSDArray->At(22))->Fill(recp.GetChargeRatio());
	dynamic_cast<TH1F *>(fHistSSDArray->At(24))->Fill(recp.GetQ());
	dynamic_cast<TH1F *>(fHistSSDArray->At(26))->Fill(chargePSide);
	dynamic_cast<TH1F *>(fHistSSDArray->At(28))->Fill(chargeNSide);
	if(chargePSide != 0.)
	  dynamic_cast<TH1F *>(fHistSSDArray->At(30))->Fill(chargeNSide/chargePSide);
	dynamic_cast<TH2F *>(fHistSSDArray->At(32))->Fill(chargePSide,
							  chargeNSide);
	dynamic_cast<TH2F *>(fHistSSDArray->At(34))->SetBinContent(gModule,lLadderLocationY,recp.GetQ());
	dynamic_cast<TH2F *>(fHistSSDArray->At(36))->Fill(gModule,499+gLadder,1);
	nClustersLayer5 += 1;
      }//layer 5 histograms
      if(layer == 5) {
	dynamic_cast<TH1F *>(fHistSSDArray->At(1))->Fill(module);
	dynamic_cast<TH1F *>(fHistSSDArray->At(5))->Fill(recp.GetDetLocalX());
	dynamic_cast<TH1F *>(fHistSSDArray->At(6))->Fill(recp.GetDetLocalZ());
	dynamic_cast<TH1F *>(fHistSSDArray->At(9))->Fill(cluglo[0]);
	dynamic_cast<TH1F *>(fHistSSDArray->At(11))->Fill(cluglo[1]);
	dynamic_cast<TH1F *>(fHistSSDArray->At(13))->Fill(cluglo[2]);
	dynamic_cast<TH1F *>(fHistSSDArray->At(15))->Fill(phi);
	dynamic_cast<TH1F *>(fHistSSDArray->At(17))->Fill(theta);
	dynamic_cast<TH1F *>(fHistSSDArray->At(19))->Fill(radius);
	dynamic_cast<TH1F *>(fHistSSDArray->At(21))->Fill(recp.GetType());
	dynamic_cast<TH1F *>(fHistSSDArray->At(23))->Fill(recp.GetChargeRatio());
	dynamic_cast<TH1F *>(fHistSSDArray->At(25))->Fill(recp.GetQ());
	dynamic_cast<TH1F *>(fHistSSDArray->At(27))->Fill(chargePSide);
	dynamic_cast<TH1F *>(fHistSSDArray->At(29))->Fill(chargeNSide);
	if(chargePSide != 0.)
	  dynamic_cast<TH1F *>(fHistSSDArray->At(31))->Fill(chargeNSide/chargePSide);
	dynamic_cast<TH2F *>(fHistSSDArray->At(33))->Fill(chargePSide,
							  chargeNSide);
	dynamic_cast<TH2F *>(fHistSSDArray->At(35))->SetBinContent(gModule,lLadderLocationY,recp.GetQ());
	dynamic_cast<TH2F *>(fHistSSDArray->At(37))->Fill(gModule,599+gLadder,1);
	nClustersLayer6 += 1;
      }//layer 6 histograms
    }//cluster loop
    dynamic_cast<TH1F *>(fHistSSDArray->At(2))->Fill(nClustersLayer5);
    dynamic_cast<TH1F *>(fHistSSDArray->At(3))->Fill(nClustersLayer6);
 }//input block
 
 //Publish array
 AliHLTUInt32_t fSpecification = 0x0;
 PushBack( (TObjArray*) fHistSSDArray,kAliHLTDataTypeTObjArray|kAliHLTDataOriginITSSSD,fSpecification);
 
 HLTInfo("ITSSSDQARecPoints found %d SSD clusters", nClustersLayer5+nClustersLayer6);

 return 0;
}
