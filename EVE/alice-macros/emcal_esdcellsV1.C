// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGLViewer.h>
#include <TGLWidget.h>
#include <TH2.h>
#include <TMath.h>
#include <TTree.h>
#include <TEveBrowser.h>
#include <TEveCalo.h>
#include <TEveCaloData.h>
#include <TEveCaloLegoOverlay.h>
#include <TEveLegoEventHandler.h>
#include <TEveManager.h>
#include <TEveScene.h>
#include <TEveTrans.h>
#include <TEveViewer.h>
#include <TEveWindow.h>

#include <AliESDEvent.h>
#include <AliEveEventManager.h>
#include <AliEveMultiView.h>
#endif

double pi = TMath::Pi();
TEveViewer *g_histo2d_v   = 0;
TEveScene  *g_histo2d_s   = 0;
TEveScene  *g_histo2d_s2  = 0;
TEveCaloLegoOverlay* g_histo2d_lego_overlay = 0;

Double_t GetPhi(Double_t phi);
TEveCaloLego* CreateHistoLego(TEveCaloData* data);
TEveCalo3D* Create3DView(TEveCaloData* data);
AliEveMultiView* CreateProjections(TEveCaloData* data, TEveCalo3D *calo3d);

Bool_t fMatrixEMSet = kFALSE;
Bool_t fMatrixPHSet = kFALSE;

TH2F* fHistoEM  = 0;
TH2F* fHistoPH  = 0;
TH2F* fHistoneg = 0;
TH2F* fHistopos = 0;

TEveCaloDataHist* emcal_esdcellsV1()
{ 
  
  // Access to esdTree
  AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
  //TTree* t = AliEveEventManager::GetMaster()->GetESDTree();
  
  TEveCaloDataHist* data = new TEveCaloDataHist();
  AliEveEventManager::GetMaster->RegisterTransient(data);
  
  // Get the EMCAL geometry
  //
  AliEMCALGeometry * geomEM  = AliEMCALGeometry::GetInstance();  
  if (!geomEM) 
  {
    printf("xxx Set default geo as Run2 xxx\n");
    geomEM  = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  }
  
  // Set the geometry alignment matrices from the ones stored in the data.
  if(!fMatrixEMSet)
  {
    printf("LOAD EMCAL MATRICES\n");
    
    for(Int_t mod = 0; mod < (geomEM->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
    { 
      printf("Load EMCAL ESD matrix %d, %p\n",mod,esd->GetEMCALMatrix(mod));
      
      if( esd->GetEMCALMatrix(mod) ) 
        geomEM->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
      else // set default identity matrix
        geomEM->SetMisalMatrix((new TGeoHMatrix),mod) ;
    }// loop over super modules	     
    
    fMatrixEMSet = kTRUE;
  }
  
  // Get the PHOS geometry
  //
  AliPHOSGeometry * geomPH  = AliPHOSGeometry::GetInstance();  
  if (!geomPH) 
  {
    printf("xxx Set PHOS default geo as Run2 xxx\n");
    geomPH  = AliPHOSGeometry::GetInstance("Run2");
  }
  
  if(!fMatrixPHSet)
  {
    for(Int_t mod = 0; mod < 5; mod++)
    { 
      printf("Load PHOS ESD matrix %d, %p\n",mod,esd->GetPHOSMatrix(mod));
      
      if(esd->GetPHOSMatrix(mod)) 
        geomPH->SetMisalMatrix(esd->GetPHOSMatrix(mod),mod) ;
      else // set default identity matrix
        geomPH->SetMisalMatrix((new TGeoHMatrix),mod) ;
    }// loop over modules	
    
    fMatrixPHSet = kTRUE;
  }
  
  // Creating 2D histograms
  
  if(!fHistoEM)
  {
    fHistoEM  = new TH2F("histoEMcell","EMCal Cell #eta vs #phi vs E",
                         100,-1.5,1.5,80,-pi,pi);
    fHistoPH  = new TH2F("histoPHcell","PHOS Cell #eta vs #phi vs E",
                         100,-1.5,1.5,80,-pi,pi);
    fHistopos = new TH2F("histopos_t","Histo 2d positive",
                         100,-1.5,1.5,80,-pi,pi);
    fHistoneg = new TH2F("histoneg_t","Histo 2d negative",
                         100,-1.5,1.5,80,-pi,pi);
  }
  
  // Getting current tracks for each event, filling histograms 
  //for ( int event = 0; event < t->GetEntries(); event++ ) 
  //{	
    // t->GetEntry(event);
    
    
  //
  // Tracks
  //
  
  printf("Number of Tracks %d\n",esd->GetNumberOfTracks());
  
  for ( int n = 0; n < esd->GetNumberOfTracks(); ++n ) 
  {    
    if ( esd->GetTrack(n)->GetSign() > 0 )
    {
      fHistopos->Fill(esd->GetTrack(n)->Eta(),
                      GetPhi(esd->GetTrack(n)->Phi()),
                      fabs(esd->GetTrack(n)->Pt()));
    } 
    else 
    {
      fHistoneg->Fill(esd->GetTrack(n)->Eta(),
                     GetPhi(esd->GetTrack(n)->Phi()),
                     fabs(esd->GetTrack(n)->Pt()));
    }
  }
  
    
  //
  // Calorimeter
  //
  
  AliESDCaloCells &cellsEM= *(esd->GetEMCALCells());
  AliESDCaloCells &cellsPH= *(esd->GetPHOSCells());
  
  Int_t ncellEM = cellsEM.GetNumberOfCells() ;  
  Int_t ncellPH = cellsPH.GetNumberOfCells() ;  
  
  printf("Number of ESD CaloCells: EM %d PH %d\n",ncellEM,ncellPH);
  
  Float_t amp   = -1 ;
  Float_t time  = -1 ;
  Int_t id      = -1 ;
  
  Float_t phi   =  0 ;
  Float_t eta   =  0 ;
  
  // Extract EMCAL cell information from the ESDs
  for (Int_t icell=  0; icell <  ncellEM; icell++) 
  {
    id  = cellsEM.GetCellNumber(icell);
    amp = cellsEM.GetAmplitude (icell); // GeV
    
    //if(amp < 0.1) continue ;
    
    geomEM->EtaPhiFromIndex(id,eta,phi);
    
    //printf("CaloCell %d, ID %d, energy %2.3f,eta %f, phi %f\n",
    //       icell,id,amp,eta,GetPhi(phi));
    
    fHistoEM->Fill(eta,GetPhi(phi),amp);
  }
  
  // Extract PHOS cells information from the ESDs
//  for (Int_t icell =  0; icell <  ncellPH; icell++) 
//  {
//    printf("cell %d, ncell %d\n",icell,ncellPH);
//    id  = cellsPH.GetCellNumber(icell);
//    amp = cellsPH.GetAmplitude (icell); // GeV
//    
//    //if(amp < 0.1) continue ;
//    
//    TVector3 xyz;
//    Int_t relId[4], module;
//    Float_t xCell, zCell;
//    
//    geomPH->AbsToRelNumbering(id,relId);
//    printf("PHOS, mod %d\n",module);
//    module = relId[0];
//    geomPH->RelPosInModule(relId,xCell,zCell);
//    geomPH->Local2Global(module,xCell,zCell,xyz);
//    
//    printf("PHOS CaloCell %d, ID %d, energy %2.3f,eta %f, phi %f\n",
//           icell,id,amp,xyz.Eta(),xyz.Phi()*TMath::RadToDeg());        
//    
//    fHistoPH->Fill(xyz.Eta(),GetPhi(xyz.Phi()),amp);
//  }
  
  
  //} event loop
  
  data->AddHistogram(fHistoneg);
  data->RefSliceInfo(0).Setup("NegCg:", 0, kBlue);
  data->AddHistogram(fHistopos);
  data->RefSliceInfo(1).Setup("PosCg:", 0, kRed);
  data->AddHistogram(fHistoEM);
  data->RefSliceInfo(2).Setup("EMCell:", 0, kViolet);
  data->AddHistogram(fHistoPH);
  data->RefSliceInfo(3).Setup("PHCell:", 0, kOrange);
  
  data->GetEtaBins()->SetTitleFont(120);
  data->GetEtaBins()->SetTitle("h");
  data->GetPhiBins()->SetTitleFont(120);
  data->GetPhiBins()->SetTitle("f");
  data->IncDenyDestroy();
  
  // Plotting the lego histogram in a new tab
  CreateHistoLego(data);
  
  // Plotting the 3D histogram
  TEveCalo3D *calo3d = Create3DView(data);
  
  // Plotting projections RPhi and RhoZ
  CreateProjections(data, calo3d);
  
  gEve->Redraw3D(kTRUE);
  
  return data;
}

//______________________________________________________________________________
Double_t GetPhi(Double_t phi)
{
  if (phi > pi) {
    phi -= 2*pi;
  }
  return phi;
}


//______________________________________________________________________________
TEveCaloLego* CreateHistoLego(TEveCaloData* data){

   TGLViewer* glv;

   // Viewer initialization, tab creation
   if ( g_histo2d_v == 0 ) {
      TEveWindowSlot *slot    = 0;
      TEveBrowser    *browser = gEve->GetBrowser();

      slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
      slot->MakeCurrent();
      g_histo2d_v = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
      g_histo2d_s = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
      g_histo2d_v->AddScene(g_histo2d_s);
      g_histo2d_v->SetElementName("2D Lego Viewer");
      g_histo2d_s->SetElementName("2D Lego Scene");

      glv = g_histo2d_v->GetGLViewer();
      g_histo2d_lego_overlay = new TEveCaloLegoOverlay();
      glv->AddOverlayElement(g_histo2d_lego_overlay);
      glv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
   } else {
      glv = g_histo2d_v->GetGLViewer(); 
   }
   
   //plotting histo
   TEveCaloLego* lego = new TEveCaloLego(data);
   g_histo2d_s->AddElement(lego);
   AliEveEventManager::GetMaster()->RegisterTransient(lego);

   // move to real world coordinates
   lego->InitMainTrans();
   Float_t sc = TMath::Min(lego->GetEtaRng(), lego->GetPhiRng());
   lego->RefMainTrans().SetScale(sc, sc, sc);

   // set event handler to move from perspective to orthographic view.
   glv->SetEventHandler(new TEveLegoEventHandler(glv->GetGLWidget(), glv, lego));

   g_histo2d_lego_overlay->SetCaloLego(lego);

   return lego;
}

//______________________________________________________________________________
TEveCalo3D* Create3DView(TEveCaloData* data){
  
   //initialization
   if ( g_histo2d_s2 == 0 ) {
      g_histo2d_s2 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
      gEve->GetDefaultViewer()->AddScene(g_histo2d_s2);
      g_histo2d_s2->SetElementName("3D Histogram Scene");
   }
 
   TEveCalo3D* calo3d = new TEveCalo3D(data);
   AliEveEventManager::GetMaster()->RegisterTransient(calo3d);
   
   calo3d->SetBarrelRadius(550);
   calo3d->SetEndCapPos(550);
   g_histo2d_s2->AddElement(calo3d);
 
   return calo3d;
}

//______________________________________________________________________________
AliEveMultiView* CreateProjections(TEveCaloData* data, TEveCalo3D *calo3d){

   AliEveMultiView *al = AliEveMultiView::Instance();
   al->ImportEventRPhi(calo3d);
   al->ImportEventRhoZ(calo3d);

   return al;
}
