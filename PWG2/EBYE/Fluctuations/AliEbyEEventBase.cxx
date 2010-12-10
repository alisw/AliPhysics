/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Satyajit Jena.                                                 *
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

//-----------------------------------------------------------------
//                 AliEbyEEventBase class
//   This is the class to deal with the EbyE Fluctuation  analysis
//               Origin: Satyajit Jena, sjena@cern.ch
//-----------------------------------------------------------------
#include <Riostream.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TH2F.h>
#include <TList.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <AliExternalTrackParam.h>
#include <AliESDEvent.h>
#include <AliPID.h>
#include <AliVertexerTracks.h>
#include <AliESDpid.h>
#include <AliTPCPIDResponse.h>
#include <AliAODTrack.h>
#include <AliMCParticle.h>
#include <AliESDCentrality.h>
#include <AliCentralitySelectionTask.h>
#include <AliESDVZERO.h>
class AliLog;
class AliESDVertex;

#include "AliEbyEEventBase.h"

ClassImp(AliEbyEEventBase)

//____________________________________________________________________//
AliEbyEEventBase::AliEbyEEventBase() : TObject(),  
  fAnalysisLevel("ESD"), 
  fAnalysisMode(kTPC),
  fDebugMode(kFALSE),
  fPhySel(0),
  fListQA(0),
  fNBinsX(0), 
  fMinX(0), 
  fMaxX(0),
  fNBinsY(0), 
  fMinY(0), 
  fMaxY(0),
  fVxMax(20.), 
  fVyMax(20.), 
  fVzMax(20.),
  fCentralityType(kFlat),
  fCentralityBin(50),  
  fCentralityEstimator("V0M"),
  fCentralityPointer(0),
  fFile1(""),
  fFile2("") {
  //Default constructor

  fListQA = new TList();
  fListQA->SetName("fListQA"); 
  TH1F *gHistVx = new TH1F("gHistVx",
			   "Vx distribution;V_{x} [cm];Entries",
			   100,-5.,5.);
  gHistVx->SetFillColor(kRed-2);
  fListQA->Add(gHistVx);//0
  TH1F *gHistVxAccepted = new TH1F("gHistVxaccepted",
				   "Vx distribution;V_{x} [cm];Entries",
				   100,-5.,5.);
  fListQA->Add(gHistVxAccepted); //1
  TH1F *gHistVy = new TH1F("gHistVy",
			   "Vy distribution;V_{y} [cm];Entries",
			   100,-5.,5.);
  gHistVy->SetFillColor(kRed-2);
  fListQA->Add(gHistVy);//2
  TH1F *gHistVyAccepted = new TH1F("gHistVyaccepted",
				   "Vy distribution;V_{y} [cm];Entries",
				   100,-5.,5.);
  fListQA->Add(gHistVyAccepted);//3
  TH1F *gHistVz = new TH1F("gHistVz",
			   "Vz distribution;V_{z} [cm];Entries",
			   100,-25.,25.);
  gHistVz->SetFillColor(kRed-2);
  fListQA->Add(gHistVz);//4
  TH1F *gHistVzAccepted = new TH1F("gHistVzaccepted",
				   "Vz distribution;V_{z} [cm];Entries",
				   100,-25.,25.);
  fListQA->Add(gHistVzAccepted);//5
 

  TH2F *hCentralityQA = new TH2F("HCentralityQA","Tracks In Centrality",fCentralityBin+1,0,(Double_t)fCentralityBin+1,1000,0,20000);
  fListQA->Add(hCentralityQA);//6
}

//____________________________________________________________________//
AliEbyEEventBase::~AliEbyEEventBase() {
  if(fListQA) delete fListQA;

}

//____________________________________________________________________//
const AliESDVertex* AliEbyEEventBase::GetVertex(AliESDEvent* esd, 
						     AnalysisMode mode,
						     Double_t gVxMax,
						     Double_t gVyMax,
						     Double_t gVzMax) {
  // Vertex selection
  const AliESDVertex* vertex = 0;
  if((mode == kITS)||(mode == kTPCnITS))
    vertex = esd->GetPrimaryVertexSPD();
  else if(mode == kTPC){
    Double_t kBz = esd->GetMagneticField();
    AliVertexerTracks vertexer(kBz);
    vertexer.SetTPCMode();
    AliESDVertex *vTPC = vertexer.FindPrimaryVertex(esd);
    esd->SetPrimaryVertexTPC(vTPC);
    for (Int_t i=0; i<esd->GetNumberOfTracks(); i++) {
      AliESDtrack *t = esd->GetTrack(i);
      t->RelateToVertexTPC(vTPC, kBz, kVeryBig);
    }
    delete vTPC;
    vertex = esd->GetPrimaryVertexTPC();
  }
  else if(mode == kGlobal)
    vertex = esd->GetPrimaryVertex();
  else
    Printf("GetVertex: ERROR: Invalid second argument %d", mode);
  
  if(!vertex) {
    //  if(fDebugMode)
      Printf("GetVertex: Event rejected because there is no valid vertex object");
    return 0;
  }
  
    
  // check resolution
  Double_t zRes = vertex->GetZRes();
  if(zRes == 0) {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because the value of the vertex resolution in z is 0");
    return 0;
  }
  ((TH1F *)(fListQA->At(0)))->Fill(vertex->GetXv());
  ((TH1F *)(fListQA->At(2)))->Fill(vertex->GetYv());
  ((TH1F *)(fListQA->At(4)))->Fill(vertex->GetZv());

  //check position
  if(TMath::Abs(vertex->GetXv()) > gVxMax) {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because it has a Vx value of %lf cm (accepted interval: -%lf - %lf)",TMath::Abs(vertex->GetXv()),gVxMax,gVxMax);
    return 0;
  }
  if(TMath::Abs(vertex->GetYv()) > gVyMax)  {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because it has a Vy value of %lf cm (accepted interval: -%lf - %lf)",TMath::Abs(vertex->GetYv()),gVyMax,gVyMax);
    return 0;
  }
  if(TMath::Abs(vertex->GetZv()) > gVzMax)  {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because it has a Vz value of %lf cm (accepted interval: -%lf - %lf)",TMath::Abs(vertex->GetZv()),gVzMax,gVzMax);
    return 0;
  }
  ((TH1F *)(fListQA->At(1)))->Fill(vertex->GetXv());
  ((TH1F *)(fListQA->At(3)))->Fill(vertex->GetYv());
  ((TH1F *)(fListQA->At(5)))->Fill(vertex->GetZv());
   
  return vertex;
}

//____________________________________________________________________//
Int_t AliEbyEEventBase::GetCentrality(AliESDEvent *esd) const {
  //Get the centrality
  AliESDCentrality *centrality = esd->GetCentrality();
  Int_t cent = -1;

  if(fCentralityType == kHardFlat) {
    AliESDVZERO* esdV0 = esd->GetVZEROData();
    Double_t mult = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();

    ((TH2F *)(fListQA->At(6)))->Fill(50,mult);
      
    if(fAnalysisLevel == "ESD") cent = FindCentralityESD(mult);
    if(fAnalysisLevel == "MC") cent = FindCentralityMC(mult);
    if(cent > -1)  ((TH2F *)(fListQA->At(6)))->Fill(cent,mult);
    return cent;
  }

 if(fCentralityType == kHardAcc) {
   AliESDVZERO* esdV0 = esd->GetVZEROData();
   Double_t mult = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();
   ((TH2F *)(fListQA->At(6)))->Fill(50,mult);
   if(fAnalysisLevel == "ESD") cent = FindCentralityESD(mult);
   if(fAnalysisLevel == "MC") cent = FindCentralityMC(mult);
   if(cent > -1)  ((TH2F *)(fListQA->At(6)))->Fill(cent,mult);
   return cent;
 }

  if(fCentralityType == kAll) {
    if(centrality->IsEventInCentralityClass(0,100,fCentralityEstimator)) return 0;
  }
  if(fCentralityType == kFlat) {
    for(int i = 0; i < 50; i++) {
	if(centrality->IsEventInCentralityClass(i*2,i*2+2,fCentralityEstimator)) {
	    cent = i;
	    break;
	}
    }
  }
  
  if(fCentralityType == kAcc) {
    for(int i = 0; i < 50; i++) {
      if(centrality->IsEventInCentralityClass(0,i*2+2,fCentralityEstimator)) {
	cent = i;
	break;
      }
    }
  }
  
  return cent;
}

//____________________________________________________________________//
Int_t AliEbyEEventBase::FindCentralityESD(Double_t mult) const {
  //hardcoded centrality bins (to be removed)
  Double_t data[101] = { 250000.0,
			 17760.0,  16870.0,  16070.0,  15360.0,  14660.0, 
			 13970.0,  13320.0,  12700.0,  12140.0,  11570.0, 
			 11020.0,  10530.0,  10030.0,   9560.0,   9110.0, 
			 8680.0,   8270.0,   7870.0,   7490.0,   7120.0, 
			 6780.0,   6440.0,   6120.0,   5810.0,   5510.0, 
			 5230.0,   4950.0,   4690.0,   4440.0,   4190.0, 
			 3960.0,   3740.0,   3530.0,   3330.0,   3130.0, 
			 2940.0,   2770.0,   2600.0,   2440.0,   2290.0, 
			 2150.0,   2010.0,   1880.0,   1760.0,   1640.0, 
			 1540.0,   1440.0,   1350.0,   1260.0,   1180.0, 
			 1110.0,   1030.0,    970.0,    900.0,    840.0, 
			 790.0,    740.0,    700.0,    650.0,    610.0, 
			 580.0,    540.0,    510.0,    480.0,    450.0, 
			 430.0,    400.0,    380.0,    360.0,    340.0, 
			 320.0,    300.0,    290.0,    270.0,    250.0, 
			 240.0,    230.0,    210.0,    200.0,    190.0, 
			 180.0,    170.0,    160.0,    150.0,    150.0, 
			 140.0,    130.0,    120.0,    120.0,    110.0, 
			 100.0,    100.0,     90.0,     90.0,     80.0, 
			 80.0,     80.0,     70.0,     70.0,     70.0 };
  

 if(fCentralityType == kHardFlat) { 
    for(Int_t i = 0; i < 100; i+=2) {
      
      if( (mult < data[i]) && (mult >= data[i+2]) ) 
	return i/2;
      }
  }
  if(fCentralityType == kHardAcc) {
    for(int i = 0; i < 100; i+=2) {
      if((mult < data[0]) && (mult >data[i+2]) ) {
	return i/2; 
      }
    }
  } 
  
return -2;
  
  
}

//____________________________________________________________________//
Int_t AliEbyEEventBase::FindCentralityMC(Double_t mult) const {
  //hardcoded centrality bins (MC) ==> to be removed
  Double_t data[101] = { 250000.0, 17760.0,  16870.0,  16070.0,  15360.0,  
			 14660.0, 13970.0,  13320.0,  12700.0,  12140.0,  
			 11570.0, 
      11020.0,  10530.0,  10030.0,   9560.0,   9110.0, 
      8680.0,   8270.0,   7870.0,   7490.0,   7120.0, 
      6780.0,   6440.0,   6120.0,   5810.0,   5510.0, 
      5230.0,   4950.0,   4690.0,   4440.0,   4190.0, 
      3960.0,   3740.0,   3530.0,   3330.0,   3130.0, 
      2940.0,   2770.0,   2600.0,   2440.0,   2290.0, 
      2150.0,   2010.0,   1880.0,   1760.0,   1640.0, 
      1540.0,   1440.0,   1350.0,   1260.0,   1180.0, 
      1110.0,   1030.0,    970.0,    900.0,    840.0, 
      790.0,    740.0,    700.0,    650.0,    610.0, 
      580.0,    540.0,    510.0,    480.0,    450.0, 
      430.0,    400.0,    380.0,    360.0,    340.0, 
      320.0,    300.0,    290.0,    270.0,    250.0, 
      240.0,    230.0,    210.0,    200.0,    190.0, 
      180.0,    170.0,    160.0,    150.0,    150.0, 
      140.0,    130.0,    120.0,    120.0,    110.0, 
      100.0,    100.0,     90.0,     90.0,     80.0, 
      80.0,     80.0,     70.0,     70.0,     70.0 
    };

 if(fCentralityType == kHardFlat) { 
    for(Int_t i = 0; i < 100; i+=2) {
      if( (mult < data[i]) && (mult >= data[i+2]) ) 
	return i/2;
      }
  }
  if(fCentralityType == kHardAcc) {
    for(int i = 0; i < 100; i+=2) {
      if((mult < data[0]) && (mult >= data[i+2]) ) {
	return i/2; 
      }
    }
  } 
  
  return -3;
}

