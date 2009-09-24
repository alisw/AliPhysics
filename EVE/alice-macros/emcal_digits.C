// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef __CINT__

#include <TEveManager.h>
#include <TEveQuadSet.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TStyle.h>
#include <TEveTrans.h>
#include <TClonesArray.h>

#include <EveBase/AliEveEventManager.h>

#include <AliRunLoader.h>
#include <AliCluster.h>
#include <EMCAL/AliEMCALGeometry.h>
#include <EMCAL/AliEMCALDigit.h>

#include <iostream>
#endif

void emcal_digits()
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();

  rl->LoadDigits("EMCAL");
  TTree* dt = rl->GetTreeD("EMCAL", kFALSE);
  if (!dt) return;

  AliEveEventManager::AssertGeometry();

  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");

  TGeoBBox* bbbox = (TGeoBBox*) node->GetDaughter(0) ->GetVolume()->GetShape();
  TGeoBBox* sbbox = (TGeoBBox*) node->GetDaughter(10)->GetVolume()->GetShape();

  TEveElementList* l = new TEveElementList("EMCAL");
  l->SetTitle("Tooltip");
  gEve->AddElement(l);

  TEveFrameBox* frame_big = new TEveFrameBox();
  frame_big->SetFrameColorRGBA(200,200,0,50);
  frame_big->SetAABoxCenterHalfSize(0, 0, 0, bbbox->GetDX(), bbbox->GetDY(), bbbox->GetDZ());

  TEveFrameBox* frame_sml = new TEveFrameBox();
  frame_sml->SetFrameColorRGBA(200,200,0,50);
  frame_sml->SetAABoxCenterHalfSize(0, 0, 0, sbbox->GetDX(), sbbox->GetDY(), sbbox->GetDZ());

  gStyle->SetPalette(1, 0);
  TEveRGBAPalette* pal = new TEveRGBAPalette(0, 512);
  pal->SetLimits(0, 1024);

  TEveQuadSet* smodules[12];


  AliEMCALGeometry * geom  = AliEMCALGeometry::GetInstance();  
  if (!geom) geom = AliEMCALGeometry::GetInstance("","");

  for (Int_t sm=0; sm<12; ++sm)
  {
    TEveQuadSet* q = new TEveQuadSet(Form("SM %d", sm+1));
    q->SetOwnIds(kTRUE);
    q->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    q->SetDefWidth (geom->GetPhiTileSize());
    q->SetDefHeight(geom->GetEtaTileSize());

    q->RefMainTrans().SetFrom(*node->GetDaughter(sm)->GetMatrix());

    q->SetFrame(sm < 10 ? frame_big : frame_sml);
    q->SetPalette(pal);

    gEve->AddElement(q, l);
    smodules[sm] = q;
  }

  TClonesArray *digits = 0;
  dt->SetBranchAddress("EMCAL", &digits);
  dt->GetEntry(0);
  Int_t nEnt = digits->GetEntriesFast();
  AliEMCALDigit * dig;

  Float_t amp   = -1 ;
  Float_t time  = -1 ;
  Int_t id      = -1 ;
  Int_t iSupMod =  0 ;
  Int_t iTower  =  0 ;
  Int_t iIphi   =  0 ;
  Int_t iIeta   =  0 ;
  Int_t iphi    =  0 ;
  Int_t ieta    =  0 ;
  Double_t x, y, z;

  for (Int_t idig = 0; idig < nEnt; ++idig)
  {
    dig = static_cast<AliEMCALDigit *>(digits->At(idig));

    if(dig != 0) {
      id   = dig->GetId() ; //cell (digit) label
      amp  = dig->GetAmp(); //amplitude in cell (digit)
      time = dig->GetTime();//time of creation of digit after collision

      std::cout<<"Cell ID "<<id<<" Amp "<<amp<<std::endl;//" time "<<time<<endl;

      //Geometry methods
      geom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
      //Gives SuperModule and Tower numbers
      geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
					iIphi, iIeta,iphi,ieta);
      //Gives label of cell in eta-phi position per each supermodule

      std::cout<< "SModule "<<iSupMod<<"; Tower "<<iTower
	  <<"; Eta "<<iIeta<<"; Phi "<<iIphi
	  <<"; Cell Eta "<<ieta<<"; Cell Phi "<<iphi<<std::endl;

      geom->RelPosCellInSModule(id, x, y, z);
      std::cout << x <<" "<< y <<" "<< z <<std::endl;

      TEveQuadSet* q = smodules[iSupMod];
      q->AddQuad(y, z);
      q->QuadValue(TMath::Nint(amp));
      q->QuadId(new AliEMCALDigit(*dig));
    } else {
      std::cout<<"Digit pointer 0x0"<<std::endl;
    }
  }

  rl->UnloadDigits("EMCAL");


  for (Int_t sm = 0; sm < 12; ++sm)
  {
    smodules[iSupMod]->RefitPlex();
  }

  gEve->Redraw3D();
}
