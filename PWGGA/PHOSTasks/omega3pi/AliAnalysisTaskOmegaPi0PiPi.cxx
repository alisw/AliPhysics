/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Boris Polishchuk (Boris.Polishchuk@cern.ch)                    *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//------------------------------------------------------------------------// 
// Class used to do omega(782) -> pi0 pi+ pi- analysis.                   //
//------------------------------------------------------------------------//
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskOmegaPi0PiPi.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "AliESDEvent.h"
#include "TRefArray.h"
#include "TLorentzVector.h"
#include "AliESDCaloCluster.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "TString.h"
#include "AliPHOSGeoUtils.h"

ClassImp(AliAnalysisTaskOmegaPi0PiPi)

AliAnalysisTaskOmegaPi0PiPi::AliAnalysisTaskOmegaPi0PiPi() :
  AliAnalysisTaskSE(),fOutputContainer(0x0),fhM2piSel(0x0),fhDxy(0x0),fhMggSel(0x0),
  fhImpXY(0x0),fhM3pi0to2(0x0),fhM3pi2to4(0x0),fhM3pi4to6(0x0),fhM3pi6to8(0x0),
  fhM3pi8to10(0x0),fhM3pi10to12(0x0),fhM3pi12(0x0),fhM3pi0to8(0x0),fModules("123")
{
  //Default constructor.
}


AliAnalysisTaskOmegaPi0PiPi::AliAnalysisTaskOmegaPi0PiPi(const char* name) :
  AliAnalysisTaskSE(name),fOutputContainer(0x0),fhM2piSel(0x0),fhDxy(0x0),fhMggSel(0x0),
  fhImpXY(0x0),fhM3pi0to2(0x0),fhM3pi2to4(0x0),fhM3pi4to6(0x0),fhM3pi6to8(0x0),
  fhM3pi8to10(0x0),fhM3pi10to12(0x0),fhM3pi12(0x0),fhM3pi0to8(0x0),fModules("123")
{
  //Constructor used to create a named task.

  DefineOutput(1,TList::Class());
}


AliAnalysisTaskOmegaPi0PiPi::~AliAnalysisTaskOmegaPi0PiPi()
{
  //Destructor
  
  if(fOutputContainer){
    fOutputContainer->Delete() ; 
    delete fOutputContainer ;
  }
}


void AliAnalysisTaskOmegaPi0PiPi::UserCreateOutputObjects()
{
  //Allocates histograms and puts it to the output container.
  
  fOutputContainer = new TList();
  
  fhM2piSel = new TH1F("hM2pi_sel","V0 inv. mass, Dxy<2cm",100,0.,1.);
  fOutputContainer->Add(fhM2piSel);

  fhDxy = new TH1F("hDxy","Distance from V0 to primary vertex in XY-plane", 500, 0.,500.);
  fOutputContainer->Add(fhDxy);

  fhMggSel = new TH1F("hMgg_sel","Inv. mass of two clusters, pT>1GeV",100,0.,0.3);
  fOutputContainer->Add(fhMggSel);

  fhImpXY = new TH1F("hImpXY","Impact parameters in XY-plane",1000,0.,1.);
  fOutputContainer->Add(fhImpXY);


  //M(3pi) vs pT(3pi), ptrack < 2GeV
  fhM3pi0to2 = new TH2F("hM3pi_0_2","pTrack: 0-2GeV",100,0.,10., 200,0.4,1.);
  fOutputContainer->Add(fhM3pi0to2);

  //M(3pi) vs pT(3pi), 2GeV < ptrack < 4GeV
  fhM3pi2to4 = new TH2F("hM3pi_2_4","pTrack: 2-4GeV",100,0.,10., 200,0.4,1.);
  fOutputContainer->Add(fhM3pi2to4);

  //M(3pi) vs pT(3pi), 4GeV < ptrack < 6GeV
  fhM3pi4to6 = new TH2F("hM3pi_4_6","pTrack: 4-6GeV",100,0.,10., 200,0.4,1.);
  fOutputContainer->Add(fhM3pi4to6);

  //M(3pi) vs pT(3pi), 6GeV < ptrack < 8GeV
  fhM3pi6to8 = new TH2F("hM3pi_6_8","pTrack: 6-8GeV",100,0.,10., 200,0.4,1.);
  fOutputContainer->Add(fhM3pi6to8);

  //M(3pi) vs pT(3pi), 8GeV < ptrack < 10GeV
  fhM3pi8to10 = new TH2F("hM3pi_8_10","pTrack: 8-10GeV",100,0.,10., 200,0.4,1.);
  fOutputContainer->Add(fhM3pi8to10);

  //M(3pi) vs pT(3pi), 10GeV < ptrack < 12GeV
  fhM3pi10to12 = new TH2F("hM3pi_10_12","pTrack: 10-12GeV",100,0.,10., 200,0.4,1.);
  fOutputContainer->Add(fhM3pi10to12);

  //M(3pi) vs pT(3pi), ptrack > 12GeV
  fhM3pi12 = new TH2F("hM3pi_12","pTrack>12GeV",100,0.,10., 200,0.4,1.);
  fOutputContainer->Add(fhM3pi12);

  //M(3pi) vs pT(3pi), ptrack < 8GeV
  fhM3pi0to8 = new TH2F("hM3pi_0_8","pTrack: 0-8GeV",100,0.,10., 200,0.4,1.);
  fOutputContainer->Add(fhM3pi0to8);

}


void AliAnalysisTaskOmegaPi0PiPi::UserExec(Option_t* /* option */)
{
  //Main function that does all the job.
  
  AliVEvent* event = InputEvent();
  AliESDEvent* esd = (AliESDEvent*)event ;

  Double_t v[3] ; //vertex ;
  esd->GetVertex()->GetXYZ(v) ;
  TVector3 vtx(v);

  AliDebug(2,Form("Vertex: (%.3f,%.3f,%.3f)\n",vtx.X(),vtx.Y(),vtx.Z()));
  
  TRefArray * caloClustersArr  = new TRefArray();
  esd->GetPHOSClusters(caloClustersArr);

  const Int_t kNumberOfTracks = esd->GetNumberOfTracks();
  const Int_t kNumberOfV0s = esd->GetNumberOfV0s();

  AliDebug(2,Form("Tracks: %d. V0s: %d\n",kNumberOfTracks,kNumberOfV0s));
  
  Float_t xyImp,zImp,imp1,imp2;
  
  const Int_t kNumberOfPhosClusters   = caloClustersArr->GetEntries() ;
  AliDebug(2,Form("CaloClusters: %d\n", kNumberOfPhosClusters));
  if(kNumberOfPhosClusters<2){
    delete caloClustersArr;
    return;
  }  
  TLorentzVector pc1; //4-momentum of PHOS cluster 1
  TLorentzVector pc2; //4-momentum of PHOS cluster 2
  TLorentzVector pc12;

  TLorentzVector ppos; //4-momentum of positive track
  TLorentzVector pneg; //4-momentum of negative track
  TLorentzVector ptr12;

  TLorentzVector pomega;
  Double_t etrack;
  Double_t p1,p2; // 3-momentum of tracks
  
  Int_t relid[4] ;
  UShort_t *absIds;

  Int_t minCells=3; //Ncells>2
  const Double_t kMpi = 0.1396;
  AliPHOSGeoUtils* geom = new AliPHOSGeoUtils("IHEP","");
  
  for(Int_t iClu=0; iClu<kNumberOfPhosClusters; iClu++) {
    AliESDCaloCluster *c1 = (AliESDCaloCluster *) caloClustersArr->At(iClu);
    if(c1->GetNCells()<minCells) continue;
 
    absIds = c1->GetCellsAbsId();
    geom->AbsToRelNumbering(absIds[0], relid) ;
    TString phosMod1;
    phosMod1 += relid[0];

    if(!fModules.Contains(phosMod1.Data())) continue;
    c1->GetMomentum(pc1,v);

    for (Int_t jClu=iClu; jClu<kNumberOfPhosClusters; jClu++) {
      AliESDCaloCluster *c2 = (AliESDCaloCluster *) caloClustersArr->At(jClu);
      if(c2->IsEqual(c1)) continue;
      if(c2->GetNCells()<minCells) continue;

      absIds = c2->GetCellsAbsId();
      geom->AbsToRelNumbering(absIds[0], relid) ; 
      TString phosMod2;
      phosMod2 += relid[0];

      if(!fModules.Contains(phosMod2.Data())) continue;
      c2->GetMomentum(pc2,v);

      pc12 = pc1+pc2;
      AliDebug(2,Form("pc12.M(): %.3f\n",pc12.M()));

      if(pc12.M()<0.100 || pc12.M()>0.160) continue; //not a pi0 candidate!
      if(pc12.Pt()<1.5) continue; //pT(pi0) > 1.5GeV
      
      fhMggSel->Fill(pc12.M());
      
      for(Int_t iTrack=0; iTrack<kNumberOfTracks; iTrack++) {
	AliESDtrack* tr1 = esd->GetTrack(iTrack);
	p1 = tr1->P();
	tr1->GetImpactParameters(xyImp,zImp);
	imp1 = TMath::Abs(xyImp);
	fhImpXY->Fill(imp1);
	Short_t charge = tr1->Charge();
	if(imp1>0.004) continue; // not from the primary vertex!

	etrack = TMath::Sqrt(p1*p1 + kMpi*kMpi);
	ppos.SetPxPyPzE(tr1->Px(),tr1->Py(),tr1->Pz(),etrack);
	
	for(Int_t jTrack=iTrack; jTrack<kNumberOfTracks; jTrack++) {
	  AliESDtrack* tr2 = esd->GetTrack(jTrack);
	  p2 = tr2->P();
	  if(tr2->Charge()==charge) continue; //same sign track (WRONG!!)
	  tr2->GetImpactParameters(xyImp,zImp);
	  imp2=TMath::Abs(xyImp);
	  if(imp2>0.004) continue; // not from the primary vertex!
	  
	  etrack = TMath::Sqrt(p2*p2 + kMpi*kMpi);
	  pneg.SetPxPyPzE(tr2->Px(),tr2->Py(),tr2->Pz(),etrack);

	  ptr12 = ppos + pneg;

// 	  printf("ptr12.M()=%.3f, xyImp1=%f, xyImp2=%f, ch1=%d, ch2=%d, 
//                  p1=%.3f, p2=%.3f, m1=%.3f, m2=%.3f\n",
// 		 ptr12.M(),imp1,imp2,tr1->Charge(),tr2->Charge(),
// 		 tr1->P(),tr2->P(),tr1->M(),tr2->M());
	  
	  pomega = pc12 + ptr12;
	  AliDebug(2,Form("pomega.M(): %f\n",pomega.M()));

	  if( p1<2. && p2<2. )
	    fhM3pi0to2->Fill(pomega.Pt(),pomega.M());
	  
	  if( (p1>=2. && p1<4.) && (p2>=2. && p2<4.) ) 
	    fhM3pi2to4->Fill(pomega.Pt(),pomega.M());	  
	  
	  if( (p1>=4. && p1<6.) && (p2>=4. && p2<6.) ) 
	    fhM3pi4to6->Fill(pomega.Pt(),pomega.M());
	  
	  if( (p1>=6. && p1<8.) && (p2>=6. && p2<8.) ) 
	    fhM3pi6to8->Fill(pomega.Pt(),pomega.M());
	  
	  if( (p1>=8. && p1<10.) && (p2>=8. && p2<10.) ) 
	    fhM3pi8to10->Fill(pomega.Pt(),pomega.M());
	  
	  if( (p1>=10. && p1<12.) && (p2>=10. && p2<12.) ) 
	    fhM3pi10to12->Fill(pomega.Pt(),pomega.M());
	  
	  if( (p1>=12.) && (p2>=12.) ) 
	    fhM3pi12->Fill(pomega.Pt(),pomega.M());
	  
	  if( p1<8. && p2<8. )
	    fhM3pi0to8->Fill(pomega.Pt(),pomega.M());

	}
      }
      
      for(Int_t iVert=0; iVert<kNumberOfV0s; iVert++) {
	AliESDv0* v0 = esd->GetV0(iVert);
	AliESDtrack* ptrack = esd->GetTrack(v0->GetPindex()); //positive track
	AliESDtrack* ntrack = esd->GetTrack(v0->GetNindex()); //negative track

	etrack = TMath::Sqrt(ptrack->P()*ptrack->P() + kMpi*kMpi);
	ppos.SetPxPyPzE(ptrack->Px(),ptrack->Py(),ptrack->Pz(),etrack);

	etrack = TMath::Sqrt(ntrack->P()*ntrack->P() + kMpi*kMpi);
	pneg.SetPxPyPzE(ntrack->Px(),ntrack->Py(),ntrack->Pz(),etrack);

	ptr12 = ppos + pneg;

	Double_t dx = vtx.X() - v0->Xv();
	Double_t dy = vtx.Y() - v0->Yv();
	Double_t dxy = TMath::Sqrt(dx*dx + dy*dy);
	fhDxy->Fill(dxy);

	AliDebug(2,Form("V0: dxy=%.3f\n",dxy));

	if(dxy<2.) fhM2piSel->Fill(ptr12.M());
      } 
    
    }
  }
  
  delete caloClustersArr;
  PostData(1,fOutputContainer);
}
