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
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEveScene.h>
#include <TEveTrans.h>
#include <TEveViewer.h>
#include <TEveWindow.h>

#include <AliESDEvent.h>
#include <AliESDCaloCluster.h>
#include <AliESDCaloCells.h>
#include <AliEveEventManager.h>
#endif

double pi = TMath::Pi();
TEveViewer *g_histo2d_all_events_v0 = 0;
TEveViewer *g_histo2d_all_events_v1 = 0;
TEveViewer *g_histo2d_all_events_v2 = 0;
TEveViewer *g_histo2d_all_events_v3 = 0;
TEveScene  *g_histo2d_all_events_s0 = 0;
TEveScene  *g_histo2d_all_events_s1 = 0;
TEveScene  *g_histo2d_all_events_s2 = 0;
TEveScene  *g_histo2d_all_events_s3 = 0;
TEveCaloLegoOverlay* g_histo2d_all_events_lego_overlay = 0;
TEveWindowSlot* g_histo2d_all_events_slot = 0;

Double_t GetPhi(Double_t phi);
TEveCaloLego* CreateHistoLego(TEveCaloData* data, TEveWindowSlot* slot);
TEveCalo3D* Create3DView(TEveCaloData* data, TEveWindowSlot* slot);
void CreateProjections(TEveCaloData* data, TEveCalo3D *calo3d, TEveWindowSlot* slot1, TEveWindowSlot* slot2);

Float_t GetECross(Bool_t isEMCAL, Int_t iSM, Int_t icol, Int_t irow/*, Float_t dtcut*/) ;
void    GetModuleNumberColAndRow(Int_t absId, Bool_t isEMCAL, Int_t & iSM, Int_t & icol, Int_t & irow)  ;
void    GetMaxEnergyCellAbsId(AliESDCaloCluster* clu, Int_t & absId, Float_t & eMax)   ;

// Cluster cuts {PHOS,EMCAL}
Int_t   ncellsCut[] = {3   , 2   };
Float_t energyCut[] = {0.30, 0.30};
Float_t m02lowCut[] = {0.01, 0.10};
Float_t m20lowCut[] = {0.01, 0.01};
Float_t exoCut = 0.95;

TLorentzVector momentum;

Double_t vertex[] = {0.,0.,0.};

AliEMCALGeometry * fGeomEM  = 0;
AliPHOSGeometry  * fGeomPH  = 0;

AliESDCaloCells  * fCellsEM = 0;
AliESDCaloCells  * fCellsPH = 0;

Bool_t fMatrixEMSet = kFALSE;
Bool_t fMatrixPHSet = kFALSE;

TEveCaloDataHist* emcal_esdclustercellsV2()
{

  TEveCaloDataHist* data_t;
  
  if ( g_histo2d_all_events_slot == 0 )
  {
    Info("histo2d_all_events", "Filling histogram...");
    
    // Access to esdTree
    AliESDEvent* esd = AliEveEventManager::GetMaster()->AssertESD();
    
    TTree* t = AliEveEventManager::GetMaster()->GetESDTree();
    
    // Get the EMCAL geometry
    //
    if(!fGeomEM) 
      fGeomEM  = AliEMCALGeometry::GetInstance();  
    if (!fGeomEM) 
    {
      printf("xxx Set EMCal default geo as Run2 xxx\n");
      fGeomEM  = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
    }
    
    // Set the geometry alignment matrices from the ones stored in the data.
    if(!fMatrixEMSet)
    {
      printf("LOAD EMCAL MATRICES\n");
      
      for(Int_t mod = 0; mod < (fGeomEM->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
      { 
        printf("Load EMCAL ESD matrix %d, %p\n",mod,esd->GetEMCALMatrix(mod));
        
        if( esd->GetEMCALMatrix(mod) ) 
          fGeomEM->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
        else // set default identity matrix
          fGeomEM->SetMisalMatrix((new TGeoHMatrix),mod) ;
      }// loop over super modules	     
      
      fMatrixEMSet = kTRUE;
    }
    
    // Get the PHOS geometry
    //
    if(!fGeomPH)
      fGeomPH  = AliPHOSGeometry::GetInstance();  
    if (!fGeomPH) 
    {
      printf("xxx Set PHOS default geo as Run2 xxx\n");
      fGeomPH  = AliPHOSGeometry::GetInstance("Run2");
    }
    
    if(!fMatrixPHSet)
    {
      for(Int_t mod = 0; mod < 5; mod++)
      { 
        printf("Load PHOS ESD matrix %d, %p\n",mod,esd->GetPHOSMatrix(mod));
        if(esd->GetPHOSMatrix(mod)) 
          fGeomPH->SetMisalMatrix(esd->GetPHOSMatrix(mod),mod) ;
        else // set default identity matrix
          fGeomPH->SetMisalMatrix((new TGeoHMatrix),mod) ;
      }// loop over modules	
      
      fMatrixPHSet = kTRUE;
    }
    
    // Creating 2D histograms
    
    
    TH2F *histoEM_t  = new TH2F("histoEMcell_t","EMCal Cell #eta vs #phi vs E",
                                //214,-1.5,1.5,448,-pi,pi);
                                100,-1.5,1.5,80,-pi,pi);

    TH2F *histoPH_t  = new TH2F("histoPHcell_t","PHOS Cell #eta vs #phi vs E",
                                //584,-1.5,1.5,1224,-pi,pi);
                                100,-1.5,1.5,80,-pi,pi);

    TH2F *histopos_t = new TH2F("histopos_t","Histo 2d positive",
                                100,-1.5,1.5,80,-pi,pi);
    TH2F *histoneg_t = new TH2F("histoneg_t","Histo 2d negative",
                                100,-1.5,1.5,80,-pi,pi);
    
    // Getting current tracks for each event, filling histograms 
    for ( int event = 0; event < t->GetEntries(); event++ ) 
    {	
      t->GetEntry(event);
      
      //
      // Tracks
      //
      
      //printf("Number of Tracks %d\n",esd->GetNumberOfTracks());

      for ( int n = 0; n < esd->GetNumberOfTracks(); ++n ) 
      {    
        if ( esd->GetTrack(n)->GetSign() > 0 )
        {
          histopos_t->Fill(esd->GetTrack(n)->Eta(),
                           GetPhi(esd->GetTrack(n)->Phi()),
                           fabs(esd->GetTrack(n)->Pt()));
        } 
        else 
        {
          histoneg_t->Fill(esd->GetTrack(n)->Eta(),
                           GetPhi(esd->GetTrack(n)->Phi()),
                           fabs(esd->GetTrack(n)->Pt()));
        }
      }
      
      //
      // Calorimeter
      //
      
      fCellsEM = esd->GetEMCALCells();
      fCellsPH = esd->GetPHOSCells();
      
      Int_t ncellEM = fCellsEM->GetNumberOfCells() ;  
      Int_t ncellPH = fCellsPH->GetNumberOfCells() ;  
      
      //printf("Number of ESD CaloCells: EM %d PH %d\n",ncellEM,ncellPH);
      
      Float_t amp   = -1 ;
      Float_t time  = -1 ;
      Int_t   id    = -1 ;
      
      Float_t phi   =  0 ;
      Float_t eta   =  0 ;
      
      // Extract EMCAL/PHOS cell-cluster information from the ESDs
      printf("Number of ESD CaloClusters: %d\n",esd->GetNumberOfCaloClusters());
      for (Int_t iclus =  0; iclus < esd->GetNumberOfCaloClusters(); iclus++) 
      {
        AliESDCaloCluster * cluster = esd->GetCaloCluster(iclus);
        
        Bool_t isEMCAL = cluster->IsEMCAL();
        
        if(cluster->GetNCells() < ncellsCut[isEMCAL]) continue ;
        
        if(cluster->E()         < energyCut[isEMCAL]) continue ;

        if(cluster->GetM02()    < m02lowCut[isEMCAL]) continue ;
        
        if(cluster->GetM20()    < m20lowCut[isEMCAL]) continue ;
        
        cluster->GetMomentum(momentum,vertex);
        
        //if(isEMCAL) printf("Cluster %d, E %2.2f, eta %2.2f, phi %2.2f\n",
        //                    iclus, momentum.E(),momentum.Eta(),momentum.Phi()*TMath::RadToDeg());
         
        Int_t   absId  = -1;
        Float_t eMax   = 0.;
        GetMaxEnergyCellAbsId(cluster,absId,eMax) ;

        if ( eMax < 0 || absId < 0 ) continue;

        Int_t ism  =  -1, icol = -1, irow = -1;
        GetModuleNumberColAndRow(absId, isEMCAL, ism,icol,irow)  ;
        
        Float_t eCross = GetECross(isEMCAL, ism, icol, irow);
        
        //printf("\t SM %d, Id %d, EMax %2.2f, Ecross %2.2f, fraction %2.2f\n",ism,absId,eMax,eCross,1-eCross/eMax);
        
        if(1-eCross/eMax > exoCut) 
        {
          //printf("Remove SM %d, Id %d, EMax %2.2f, Ecross %2.2f, fraction %2.2f\n",ism,absId,eMax,eCross,1-eCross/eMax);

          continue;
        }
        //printf("\t Fill\n");
        
        for(Int_t icell = 0; icell < cluster->GetNCells(); icell++)
        {
          id  = cluster->GetCellAbsId(icell);
          
          if(isEMCAL)
          {
            amp = fCellsEM->GetCellAmplitude(id); // GeV

            //if(amp < 0.1) continue ;
            
            fGeomEM->EtaPhiFromIndex(id,eta,phi);
            
//            printf("CaloCell %d, ID %d, energy %2.2f,eta %2.2f, phi %2.2f\n",
//                   icell,id,amp,eta,GetPhi(phi)*TMath::RadToDeg());
            
            histoEM_t->Fill(eta,GetPhi(phi),amp);
          }
//          else // PHOS
//          {
//            amp = fCellsPH->GetCellAmplitude(id); // GeV
//                                                  //if(amp < 0.1) continue ;
//            
//            TVector3 xyz;
//            Int_t relId[4], module;
//            Float_t xCell, zCell;
//            
//            fGeomPH->AbsToRelNumbering(id,relId);
//            printf("PHOS, mod %d\n",module);
//            module = relId[0];
//            fGeomPH->RelPosInModule(relId,xCell,zCell);
//            fGeomPH->Local2Global(module,xCell,zCell,xyz);
//            
//            //printf("PHOS CaloCell %d, ID %d, energy %2.3f,eta %f, phi %f\n",
//            //       icell,id,amp,xyz.Eta(),xyz.Phi()*TMath::RadToDeg());        
//            
//            histoPH_t->Fill(xyz.Eta(),GetPhi(xyz.Phi()),amp);
//          }
          
        }
      }
    } // event loop
    
    data_t = new TEveCaloDataHist();
    
    data_t->AddHistogram(histoneg_t);
    data_t->RefSliceInfo(0).Setup("NegCg:", 0, kBlue);
    data_t->AddHistogram(histopos_t);
    data_t->RefSliceInfo(1).Setup("PosCg:", 0, kRed);
    data_t->AddHistogram(histoEM_t);
    data_t->RefSliceInfo(2).Setup("EMCell:", 0, kViolet);
    data_t->AddHistogram(histoPH_t);
    data_t->RefSliceInfo(3).Setup("PHCell:", 0, kOrange);
    
    data_t->GetEtaBins()->SetTitleFont(120);
    data_t->GetEtaBins()->SetTitle("h");
    data_t->GetPhiBins()->SetTitleFont(120);
    data_t->GetPhiBins()->SetTitle("f");
    data_t->IncDenyDestroy();
    
    // Creating frames
    g_histo2d_all_events_slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
    TEveWindowPack* packH = g_histo2d_all_events_slot->MakePack();
    packH->SetElementName("Projections");
    packH->SetHorizontal();
    packH->SetShowTitleBar(kFALSE);
    
    g_histo2d_all_events_slot = packH->NewSlot();
    TEveWindowPack* pack0 = g_histo2d_all_events_slot->MakePack();
    pack0->SetShowTitleBar(kFALSE);
    TEveWindowSlot*  slotLeftTop   = pack0->NewSlot();
    TEveWindowSlot* slotLeftBottom = pack0->NewSlot();
    
    g_histo2d_all_events_slot = packH->NewSlot();
    TEveWindowPack* pack1 = g_histo2d_all_events_slot->MakePack();
    pack1->SetShowTitleBar(kFALSE);
    TEveWindowSlot* slotRightTop    = pack1->NewSlot();
    TEveWindowSlot* slotRightBottom = pack1->NewSlot();
    
    // Creating viewers and scenes   
    TEveCalo3D* calo3d = Create3DView(data_t, slotLeftTop);
    CreateHistoLego(data_t, slotLeftBottom);
    CreateProjections(data_t, calo3d, slotRightTop, slotRightBottom);
    
    gEve->Redraw3D(kTRUE);
    
    Info("emcal_esdcellsV2", "...Finished");
   }

   return data_t;
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
TEveCaloLego* CreateHistoLego(TEveCaloData* data, TEveWindowSlot* slot){

   TEveCaloLego* lego;

   // Viewer initialization, tab creation
   if ( g_histo2d_all_events_v0 == 0 ) {

      TEveBrowser *browser = gEve->GetBrowser();
      slot->MakeCurrent();
      g_histo2d_all_events_v0 = gEve->SpawnNewViewer("2D Lego Histogram", "2D Lego Histogram");
      g_histo2d_all_events_s0 = gEve->SpawnNewScene("2D Lego Histogram", "2D Lego Histogram");
      g_histo2d_all_events_v0->AddScene(g_histo2d_all_events_s0);
      g_histo2d_all_events_v0->SetElementName("2D Lego Viewer");
      g_histo2d_all_events_s0->SetElementName("2D Lego Scene");

      TGLViewer* glv = g_histo2d_all_events_v0->GetGLViewer();
      g_histo2d_all_events_lego_overlay = new TEveCaloLegoOverlay();
      glv->AddOverlayElement(g_histo2d_all_events_lego_overlay);
      glv->SetCurrentCamera(TGLViewer::kCameraPerspXOY);

      // Plotting histogram lego
      lego = new TEveCaloLego(data);
      g_histo2d_all_events_s0->AddElement(lego);

      // Move to real world coordinates
      lego->InitMainTrans();
      Float_t sc = TMath::Min(lego->GetEtaRng(), lego->GetPhiRng());
      lego->RefMainTrans().SetScale(sc, sc, sc);

      // Set event handler to move from perspective to orthographic view.
      glv->SetEventHandler(new TEveLegoEventHandler(glv->GetGLWidget(), glv, lego));

      g_histo2d_all_events_lego_overlay->SetCaloLego(lego);
   }

   return lego;
}

//______________________________________________________________________________
TEveCalo3D* Create3DView(TEveCaloData* data, TEveWindowSlot* slot){

   TEveCalo3D* calo3d;

   if ( g_histo2d_all_events_v1 == 0 ) {
      
      TEveBrowser *browser = gEve->GetBrowser();
      slot->MakeCurrent();
      g_histo2d_all_events_v1 = gEve->SpawnNewViewer("3D Histogram", "3D Histogram");
      g_histo2d_all_events_s1 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
      g_histo2d_all_events_v1->AddScene(g_histo2d_all_events_s1);
      g_histo2d_all_events_v1->SetElementName("3D Histogram Viewer");
      g_histo2d_all_events_s1->SetElementName("3D Histogram Scene");

      calo3d = new TEveCalo3D(data);
   
      calo3d->SetBarrelRadius(550);
      calo3d->SetEndCapPos(550);
      g_histo2d_all_events_s1->AddElement(calo3d);
   } 

   return calo3d;
}

//______________________________________________________________________________
void CreateProjections(TEveCaloData* data, TEveCalo3D *calo3d, TEveWindowSlot* slot1, TEveWindowSlot* slot2){

   if ( g_histo2d_all_events_v2 == 0 ) {
      
      TEveBrowser *browser = gEve->GetBrowser();
      slot1->MakeCurrent();
      g_histo2d_all_events_v2 = gEve->SpawnNewViewer("RPhi projection", "RPhi projection");
      g_histo2d_all_events_s2 = gEve->SpawnNewScene("RPhi projection", "RPhi projection");
      g_histo2d_all_events_v2->AddScene(g_histo2d_all_events_s2);
      g_histo2d_all_events_v2->SetElementName("RPhi Projection Viewer");
      g_histo2d_all_events_s2->SetElementName("RPhi Projection Scene");

      TEveProjectionManager* mng1 = new TEveProjectionManager();
      mng1->SetProjection(TEveProjection::kPT_RPhi);

      TEveProjectionAxes* axeg_histo2d_all_events_s1 = new TEveProjectionAxes(mng1);
      g_histo2d_all_events_s2->AddElement(axeg_histo2d_all_events_s1);
      TEveCalo2D* calo2d1 = (TEveCalo2D*) mng1->ImportElements(calo3d);
      g_histo2d_all_events_s2->AddElement(calo2d1);

      g_histo2d_all_events_v2->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   } 

   if ( g_histo2d_all_events_v3 == 0 ) {
      
      TEveBrowser *browser = gEve->GetBrowser();
      slot2->MakeCurrent();
      g_histo2d_all_events_v3 = gEve->SpawnNewViewer("RhoZ projection", "RhoZ projection");
      g_histo2d_all_events_s3 = gEve->SpawnNewScene("RhoZ projection", "RhoZ projection");
      g_histo2d_all_events_v3->AddScene(g_histo2d_all_events_s3);
      g_histo2d_all_events_v3->SetElementName("RhoZ Projection Viewer");
      g_histo2d_all_events_s3->SetElementName("RhoZ Projection Viewer");

      TEveProjectionManager* mng2 = new TEveProjectionManager();
      mng2->SetProjection(TEveProjection::kPT_RhoZ);

      TEveProjectionAxes* axeg_histo2d_all_events_s2 = new TEveProjectionAxes(mng2);
      g_histo2d_all_events_s3->AddElement(axeg_histo2d_all_events_s2);
      TEveCalo2D* calo2d2 = (TEveCalo2D*) mng2->ImportElements(calo3d);
      g_histo2d_all_events_s3->AddElement(calo2d2);

      g_histo2d_all_events_v3->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   } 

   return;
}

//______________________________________________________________________________
Float_t GetECross(Bool_t isEMCAL, Int_t imod, Int_t icol, Int_t irow)
{  
  Float_t  ecell1 =  0, ecell2  = 0, ecell3  = 0, ecell4  = 0;
  
  Int_t    absId1 = -1, absId2 = -1, absId3 = -1, absId4 = -1;
  
  AliESDCaloCells * cells = 0;
  
  if ( isEMCAL )
  {
    // Get close cells index, energy and time, not in corners
    
    cells = fCellsEM;
    
    Int_t rowMax = AliEMCALGeoParams::fgkEMCALRows;
    Int_t colMax = AliEMCALGeoParams::fgkEMCALCols;
    
    if(imod == 11 || imod == 10 || imod == 18 || imod == 19) 
      rowMax = AliEMCALGeoParams::fgkEMCALRows/3;
    
    if(imod > 12 && imod < 18)
      colMax = AliEMCALGeoParams::fgkEMCALCols*2/3;

    if( irow < rowMax) absId1 = fGeomEM->GetAbsCellIdFromCellIndexes(imod, irow+1, icol);
    if( irow > 0     ) absId2 = fGeomEM->GetAbsCellIdFromCellIndexes(imod, irow-1, icol);
    
    if( icol < colMax-1 )
      absId3 = fGeomEM->GetAbsCellIdFromCellIndexes(imod, irow, icol+1);
    if( icol > 0 )    
      absId4 = fGeomEM->GetAbsCellIdFromCellIndexes(imod, irow, icol-1);

    // In case of cell in eta = 0 border, depending on SM shift the cross cell index
    if(imod > 11 && imod < 18)
    {
      if     ( icol == colMax - 1 && !(imod%2) )
      {
        absId3 = fGeomEM->GetAbsCellIdFromCellIndexes(imod+1, irow, 0);
        absId4 = fGeomEM->GetAbsCellIdFromCellIndexes(imod  , irow, icol-1); 
      }
      else if( icol == 0 && imod%2 )
      {
        absId3 = fGeomEM->GetAbsCellIdFromCellIndexes(imod  , irow, icol+1);
        absId4 = fGeomEM->GetAbsCellIdFromCellIndexes(imod-1, irow, colMax-1); 
      }
    }
  }
  else // PHOS
  {       
    cells = fCellsPH;
    
    Int_t relId1[] = { imod+1, 0, irow+1, icol   };
    Int_t relId2[] = { imod+1, 0, irow-1, icol   };
    Int_t relId3[] = { imod+1, 0, irow  , icol+1 };
    Int_t relId4[] = { imod+1, 0, irow  , icol-1 };
    
    fGeomPH->RelToAbsNumbering(relId1, absId1);
    fGeomPH->RelToAbsNumbering(relId2, absId2);
    fGeomPH->RelToAbsNumbering(relId3, absId3);
    fGeomPH->RelToAbsNumbering(relId4, absId4);
  }
  
  if(absId1 > 0 ) ecell1 = cells->GetCellAmplitude(absId1);
  if(absId2 > 0 ) ecell2 = cells->GetCellAmplitude(absId2);
  if(absId3 > 0 ) ecell3 = cells->GetCellAmplitude(absId3);
  if(absId4 > 0 ) ecell4 = cells->GetCellAmplitude(absId4);
  
  //    if(TMath::Abs(tcell-tcell1)*1.e9 > dtcut) ecell1 = 0 ;
  //    if(TMath::Abs(tcell-tcell2)*1.e9 > dtcut) ecell2 = 0 ;
  //    if(TMath::Abs(tcell-tcell3)*1.e9 > dtcut) ecell3 = 0 ;
  //    if(TMath::Abs(tcell-tcell4)*1.e9 > dtcut) ecell4 = 0 ;
  
  return ecell1+ecell2+ecell3+ecell4;

}

//___________________________________________________________________________________________________
void GetModuleNumberColAndRow(Int_t absId, Bool_t isEMCAL,
                              Int_t & iSM, Int_t & icol, Int_t & irow) 
{
  Int_t imod = -1;
  
  if ( absId < 0) return -1 ;
  
  if ( isEMCAL )
  {
    Int_t iTower = -1, iIphi = -1, iIeta = -1;
    fGeomEM->GetCellIndex(absId,iSM,iTower,iIphi,iIeta);
    fGeomEM->GetCellPhiEtaIndexInSModule(iSM,iTower,iIphi, iIeta,irow,icol);
  } // EMCAL
  else //PHOS
  {
    Int_t    relId[4];
    fGeomPH->AbsToRelNumbering(absId,relId);
    irow = relId[2];
    icol = relId[3];
    iSM  = relId[0]-1;
   
  }//PHOS
 
  if(iSM < 0 )
    AliFatal(Form("Negative value for super module: %d",imod));
}

//________________________________________________________________________________________
void GetMaxEnergyCellAbsId(AliESDCaloCluster* clu, Int_t & absId, Float_t & eMax) 
{  
  Double_t eCell       =-1.;
  Int_t    cellAbsId   =-1 ; 
  Int_t    iSupMod     =-1 ;
  
  AliESDCaloCells*  cells = 0;
  
  if(clu->IsEMCAL()) cells = fCellsEM;
  else               cells = fCellsPH;
  
  for (Int_t iDig=0; iDig < clu->GetNCells(); iDig++) 
  {
    cellAbsId = clu->GetCellAbsId(iDig);
    
    eCell  = cells->GetCellAmplitude(cellAbsId);
    
    if(eCell > eMax)  
    { 
      eMax  = eCell; 
      absId = cellAbsId;
    }
  }// cell loop
}

