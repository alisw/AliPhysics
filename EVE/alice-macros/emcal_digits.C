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
#include <TTree.h>

#include <EveBase/AliEveEventManager.h>

#include <AliRunLoader.h>
#include <AliCluster.h>
#include <EMCAL/AliEMCALGeometry.h>
#include <EMCAL/AliEMCALDigit.h>
#include <AliLog.h>

// #include <Riostream.h>
#endif

void emcal_digits()
{
  AliEveEventManager::AssertGeometry();

  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");
  if (!node) return;

  Int_t nModules = node->GetNdaughters();

  TEveElementList* l = new TEveElementList("EMCAL");
  l->SetTitle("Tooltip");
  gEve->AddElement(l);

  TGeoBBox* bbbox = (TGeoBBox*) node->GetDaughter(0) ->GetVolume()->GetShape();
  TEveFrameBox* frame_big = new TEveFrameBox();
  frame_big->SetFrameColorRGBA(200,200,0,50);
  frame_big->SetAABoxCenterHalfSize(0, 0, 0, bbbox->GetDX(), bbbox->GetDY(), bbbox->GetDZ());

  TEveFrameBox* frame_sml = 0x0;

  if (nModules==12) {
    TGeoBBox* sbbox = (TGeoBBox*) node->GetDaughter(10)->GetVolume()->GetShape();
    frame_sml = new TEveFrameBox();
    frame_sml->SetFrameColorRGBA(200,200,0,50);
    frame_sml->SetAABoxCenterHalfSize(0, 0, 0, sbbox->GetDX(), sbbox->GetDY(), sbbox->GetDZ());
  }

  gStyle->SetPalette(1, 0);
  TEveRGBAPalette* pal = new TEveRGBAPalette(0, 512);
  pal->SetLimits(0, 1024);

  TEveQuadSet* smodules[12];
  memset(smodules,0,12*sizeof(TEveQuadSet*));


  AliEMCALGeometry * geom  = AliEMCALGeometry::GetInstance();  
  if (!geom) geom = AliEMCALGeometry::GetInstance("","");

  for (Int_t sm=0; sm<nModules; ++sm)
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

  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();

  rl->LoadDigits("EMCAL");
  TTree* dt = rl->GetTreeD("EMCAL", kFALSE);
  if (!dt) return;

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

      AliDebugGeneral("emcal_digits", 5, Form("Cell ID %3d, Amplitude: %f", id, amp));
      // cout<<"Cell ID "<<id<<" Amp "<<amp<<endl;//" time "<<time<<endl;

      //Geometry methods
      geom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
      //Gives SuperModule and Tower numbers
      geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
					iIphi, iIeta,iphi,ieta);
      //Gives label of cell in eta-phi position per each supermodule

      AliDebugGeneral("emcal_digits", 5, Form("SModule %3d; Tover %3d; Eta %3d; Phi %3d; Cell Eta %3d; Cell Phi %3d",
		       iSupMod, iTower, iIeta, iIphi, ieta, iphi));
      // cout<< "SModule "<<iSupMod<<"; Tower "<<iTower
      //     <<"; Eta "<<iIeta<<"; Phi "<<iIphi
      //     <<"; Cell Eta "<<ieta<<"; Cell Phi "<<iphi<<endl;

      geom->RelPosCellInSModule(id, x, y, z);
      // cout << x <<" "<< y <<" "<< z <<endl;
      AliDebugGeneral("emcal_digits", 5, Form("(x,y,z)=(%8.3f,%8.3f,%8.3f)", x, y, z));

      TEveQuadSet* q = smodules[iSupMod];
      if (q) {
	q->AddQuad(y, z);
	q->QuadValue(TMath::Nint(amp));
	q->QuadId(new AliEMCALDigit(*dig));
      }
    } else {
      AliDebugGeneral("emcal_digits", 1, Form("Digit pointer 0x0"));
      // cout<<"Digit pointer 0x0"<<endl;
    }
  }

  rl->UnloadDigits("EMCAL");


  for (Int_t sm = 0; sm < nModules; ++sm)
  {
    smodules[iSupMod]->RefitPlex();
  }

  gEve->Redraw3D();
}
