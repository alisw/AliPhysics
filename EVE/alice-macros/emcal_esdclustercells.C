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

#include <TGeoNode.h>
#include <TGeoManager.h>
#include <TEveFrameBox.h>
#include <TGeoBBox.h>
#include <TStyle.h>
#include <TEveRGBAPalette.h>
#include <TEveQuadSet.h>

#include <AliESDEvent.h>
#include <AliEveEventManager.h>
#include <AliEveMultiView.h>

#include "AliEMCALGeometry.h"
#include "AliPHOSGeometry.h"
#endif

double pi = TMath::Pi();
TEveViewer *g_histo2d_v   = 0;
TEveScene  *g_histo2d_s   = 0;
TEveScene  *g_histo2d_s2  = 0;
TEveCaloLegoOverlay* g_histo2d_lego_overlay = 0;

Double_t GetPhi(Double_t phi);
TEveCaloLego* CreateHistoLego(TEveCaloData* data);
TEveCalo3D* Create3DView(TEveCaloData* data);
AliEveMultiView* CreateProjections(TEveCalo3D *calo3d);

/// Calorimeter methods
void GetMaxEnergyCellAbsId(AliESDCaloCluster* clu, Int_t & absId, Float_t & eMax);
void GetModuleNumberColAndRow(Int_t absId, Bool_t isEMCAL,
                              Int_t & iSM, Int_t & icol, Int_t & irow);
Float_t GetECross(Bool_t isEMCAL, Int_t imod, Int_t icol, Int_t irow);

/// Cluster cuts {PHOS,EMCAL}
Int_t   ncellsCut[] = { 3   , 2    }; /// Number of cells in cluster must be larger than this value.
Float_t energyCut[] = { 0.30, 0.30 }; /// Cluster energy must be larger than this value.
Float_t m02lowCut[] = { 0.01, 0.10 }; /// Cluster shower shape major axis must be larger than this value.
Float_t m20lowCut[] = { 0.01, 0.01 }; /// Cluster shower shape lower axis must be larger than this value.
Float_t exoCut = 0.95;                /// Reject clusters with this exoticity value.
Bool_t  plotcells[] = {kTRUE , kFALSE}; /// Display cells in clusters or cluster center.

// Declare once different parameters/histograms
TLorentzVector momentum;

Double_t vertex[] = {0.,0.,0.};

AliEMCALGeometry * fGeomEM  = 0;
AliPHOSGeometry  * fGeomPH  = 0;

AliESDCaloCells  * fCellsEM = 0;
AliESDCaloCells  * fCellsPH = 0;

Bool_t fMatrixEMSet = kFALSE;
Bool_t fMatrixPHSet = kFALSE;

TH2F* fHistoEM  = 0;
TH2F* fHistoPH  = 0;
//TH2F* fHistoneg = 0;
//TH2F* fHistopos = 0;

//TEveCaloDataHist* 
void emcal_esdclustercells()
{ 
  printf("Execute emcal_esdclustercells\n");
  
  AliEveEventManager::AssertGeometry();

  // Access to esdTree
  AliESDEvent* esd = AliEveEventManager::AssertESD();
  //TTree* t = AliEveEventManager::GetMaster()->GetESDTree();
  
  TEveCaloDataHist* data = new TEveCaloDataHist();
  AliEveEventManager::RegisterTransient(data);
  
  //-------------------------
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
  if( !fMatrixEMSet )
  {
    //printf("LOAD EMCAL MATRICES\n");
    
    for(Int_t mod = 0; mod < (fGeomEM->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
    { 
      //printf("Load EMCAL ESD matrix %d, %p\n",mod,esd->GetEMCALMatrix(mod));
      
      if( esd->GetEMCALMatrix(mod) ) 
        fGeomEM->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
      else // set default identity matrix
        fGeomEM->SetMisalMatrix((new TGeoHMatrix),mod) ;
    }// loop over super modules	     
    
    fMatrixEMSet = kTRUE;
  }
  
  //-------------------------
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
      //printf("Load PHOS ESD matrix %d, %p\n",mod,esd->GetPHOSMatrix(mod));
      if(esd->GetPHOSMatrix(mod)) 
        fGeomPH->SetMisalMatrix(esd->GetPHOSMatrix(mod),mod) ;
      else // set default identity matrix
        fGeomPH->SetMisalMatrix((new TGeoHMatrix),mod) ;
    }// loop over modules	
    
    fMatrixPHSet = kTRUE;
  }
  
  //-------------------------
  // EMCAL volumes 
  
  TGeoNode* node = gGeoManager->GetTopVolume()->FindNode("XEN1_1");
  if (!node) return ;//0;
  
  Int_t nModules = node->GetNdaughters();
  
  // Check that the EMCAL geo and the nodes from EMCAL have the same number of entries
  //
  if(nModules != fGeomEM->GetNumberOfSuperModules())
    printf("*** === EMCAL DIGITS - N Daughter Nodes %d - N super mod %d === ***\n", 
           node->GetNdaughters(), fGeomEM->GetNumberOfSuperModules());
  
  // Get the EMCAL bounding boxes for the super modules.
  // 4 kind of SM: 10 Full EMCal, 2 1/3 EMCal, 6 DCal (2/3 EMCal) and 2 1/3 EMCal in DCal region.
  //
  TGeoBBox* bbbox = (TGeoBBox*) node->GetDaughter(0) ->GetVolume()->GetShape();
  TEveFrameBox* frame_big = new TEveFrameBox();
  frame_big->SetFrameColorRGBA(200,200,0,50);
  frame_big->SetAABoxCenterHalfSize(0, 0, 0, bbbox->GetDX(), bbbox->GetDY(), bbbox->GetDZ());
  
  TEveFrameBox* frame_sml  = 0x0;
  TEveFrameBox* frame_dcl  = 0x0;
  TEveFrameBox* frame_smld = 0x0;
  
  if (nModules > 10) 
  {
    TGeoBBox* sbbox = (TGeoBBox*) node->GetDaughter(10)->GetVolume()->GetShape();
    frame_sml = new TEveFrameBox();
    frame_sml->SetFrameColorRGBA(200,200,0,50);
    frame_sml->SetAABoxCenterHalfSize(0, 0, 0, sbbox->GetDX(), sbbox->GetDY(), sbbox->GetDZ());
  }
  
  if (nModules > 12) 
  {
    TGeoBBox* dbbox = (TGeoBBox*) node->GetDaughter(12)->GetVolume()->GetShape();
    frame_dcl = new TEveFrameBox();
    frame_dcl->SetFrameColorRGBA(200,200,0,50);
    frame_dcl->SetAABoxCenterHalfSize(0, 0, 0, dbbox->GetDX(), dbbox->GetDY(), dbbox->GetDZ());
    
    TGeoBBox* sdbbox = (TGeoBBox*) node->GetDaughter(18)->GetVolume()->GetShape();
    frame_smld = new TEveFrameBox();
    frame_smld->SetFrameColorRGBA(200,200,0,50);
    frame_smld->SetAABoxCenterHalfSize(0, 0, 0, sdbbox->GetDX(), sdbbox->GetDY(), sdbbox->GetDZ());
  }

  //--------------------------------
  // Define EVE stuff for EMCAL 3d view
  //
  TEveElementList* l = new TEveElementList("EMCAL");
  l->SetTitle("Tooltip");
  gEve->AddElement(l);
  
  gStyle->SetPalette(1, 0);
  TEveRGBAPalette* pal = new TEveRGBAPalette(0, 10);
//  pal->SetLimits(0, 1024);
  
  // Here we will store the EMCAL data that will be treated by EVE
  // per each super-module.
  // Pass the SM bounding boxes (frames).
  //
  const Int_t nSM = nModules;
  TEveQuadSet* smodules[nSM];
  memset(smodules,0,nModules*sizeof(TEveQuadSet*));
  
  for (Int_t sm = 0; sm < nModules; ++sm)
  {
    TEveQuadSet* q = new TEveQuadSet(Form("SM %d", sm+1));
    q->SetOwnIds(kTRUE);
    
    // Type of object to be displayed, rectangle with cell size
    q->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    q->SetDefWidth (fGeomEM->GetPhiTileSize());
    q->SetDefHeight(fGeomEM->GetEtaTileSize());
    
    //printf(">>>> Tile size eta %f  phi %f <<<<\n", fGeomEM->GetEtaTileSize(),fGeomEM->GetPhiTileSize());
    
    q->RefMainTrans().SetFrom(*node->GetDaughter(sm)->GetMatrix());
    
//    if     (sm < 10) q->SetFrame(frame_big );
//    else if(sm < 12) q->SetFrame(frame_sml );
//    else if(sm < 18) q->SetFrame(frame_dcl );
//    else if(sm < 20) q->SetFrame(frame_smld);
    
    q->SetPalette(pal);
    
    gEve->AddElement(q, l);
    smodules[sm] = q;
  }

  //-------------------------
  // PHOS volumes 
  
  // Define EVE stuff for PHOS 3d view
  //
  TEveElementList* lPH = new TEveElementList("PHOS");
  lPH->SetTitle("Tooltip");
  gEve->AddElement(lPH);
  
  TEveRGBAPalette* palPH = new TEveRGBAPalette(0, 100);
  palPH->SetLimits(0, 200);
  
  // Here we will store the PHOS data that will be treated by EVE
  // per each module.
  // Pass the module bounding boxes (frames).
  //
  TEveQuadSet* smodulesPH[4];
  memset(smodulesPH,0,4*sizeof(TEveQuadSet*));
  
  for (Int_t mod = 0; mod < 4; ++mod)
  {
    TGeoNode* nodePH = gGeoManager->GetTopVolume()->FindNode(Form("PHOS_%d",mod+1));
        
    if(!nodePH) 
    {
      printf("PHOS node %d not found\n",mod+1);
      continue;
    }
    
    // Get the PHOS bounding boxes for the modules.
    //
    TGeoBBox* bbboxPH = (TGeoBBox*) nodePH->GetDaughter(0)->GetVolume()->GetShape();
    TEveFrameBox* frame_PH = new TEveFrameBox();
    frame_PH->SetFrameColorRGBA(200,200,0,50);
    frame_PH->SetAABoxCenterHalfSize(0, 0, 0, bbboxPH->GetDX(), bbboxPH->GetDY(), bbboxPH->GetDZ());
    
    TEveQuadSet* qPH = new TEveQuadSet(Form("Mod %d", mod+1));
    qPH->SetOwnIds(kTRUE);
    
    // Type of object to be displayed, rectangle with cell size
    qPH->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    qPH->SetDefWidth (2.2);
    qPH->SetDefHeight(2.2);
    
    qPH->RefMainTrans().SetFrom(*nodePH->GetDaughter(0)->GetMatrix());
    
    qPH->SetFrame(frame_PH);
 
    qPH->SetPalette(pal);
    
    gEve->AddElement(qPH, lPH);
    smodulesPH[mod] = qPH;
  }

  //-------------------------
  // Creating 2D histograms for projection plots
  
  if(!fHistoEM)
  {
    fHistoEM  = new TH2F("histoEMcell","EMCal Cell #eta vs #phi vs E",
                              100,-1.5,1.5,80,-pi,pi);
    fHistoPH  = new TH2F("histoPHcell","PHOS Cell #eta vs #phi vs E",
                              100,-1.5,1.5,80,-pi,pi);
//    fHistopos = new TH2F("histopos_t","Histo 2d positive",
//                              100,-1.5,1.5,80,-pi,pi);
//    fHistoneg = new TH2F("histoneg_t","Histo 2d negative",
//                              100,-1.5,1.5,80,-pi,pi);
  }
  else
  {
    fHistoEM ->Reset();
    fHistoPH ->Reset();   
//    fHistopos->Reset();
//    fHistoneg->Reset();
  }
  
  //-------------------------
  // Data analysis
    
  // Getting current tracks for each event, filling histograms 
  //for ( int event = 0; event < t->GetEntries(); event++ ) 
  //{	
    // t->GetEntry(event);
  
  //
  // Tracks
  //
  
//  printf("Number of Tracks %d\n",esd->GetNumberOfTracks());
//  
//  for ( int n = 0; n < esd->GetNumberOfTracks(); ++n ) 
//  {    
//    if ( esd->GetTrack(n)->GetSign() > 0 )
//    {
//      fHistopos->Fill(esd->GetTrack(n)->Eta(),
//                      GetPhi(esd->GetTrack(n)->Phi()),
//                      fabs(esd->GetTrack(n)->Pt()));
//    } 
//    else 
//    {
//      fHistoneg->Fill(esd->GetTrack(n)->Eta(),
//                      GetPhi(esd->GetTrack(n)->Phi()),
//                      fabs(esd->GetTrack(n)->Pt()));
//    }
//  }
  
  
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
  //printf("Number of ESD CaloClusters: %d\n",esd->GetNumberOfCaloClusters());
  for (Int_t iclus =  0; iclus < esd->GetNumberOfCaloClusters(); iclus++) 
  {
    AliESDCaloCluster * cluster = esd->GetCaloCluster(iclus);
    
    Bool_t isEMCAL = cluster->IsEMCAL();
    
    //if(!isEMCAL) printf("PHOS cluster  \n");
    //else         printf("EMCAL cluster \n");

    if(cluster->GetNCells() < ncellsCut[isEMCAL]) continue ;
    
    if(cluster->E()         < energyCut[isEMCAL]) continue ;
    
    if(cluster->GetM02()    < m02lowCut[isEMCAL]) continue ;
    
    if(cluster->GetM20()    < m20lowCut[isEMCAL]) continue ;
    
    cluster->GetMomentum(momentum,vertex);
    
    //printf("Cluster %d, E %2.2f, eta %2.2f, phi %2.2f\n",
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
    
    // Plot EMCal
    //
    if(isEMCAL)
    {
      // Plot clusters or cluster cells?
      Double_t x=0., y=0., z=0.;
      Int_t iSupMod =  0 ;
      Float_t amp   = -1 ;
      Float_t time  = -1 ;
      Int_t id      = -1 ;
      Int_t iTower  =  0 ;
      Int_t iIphi   =  0 ;
      Int_t iIeta   =  0 ;
      Int_t iphi    =  0 ;
      Int_t ieta    =  0 ;

      if( plotcells[1] ) // cells
      {
        //printf("*** EMCAL cluster cells ***\n");
        
        for(Int_t icell = 0; icell < cluster->GetNCells(); icell++)
        {
          id  = cluster->GetCellAbsId(icell);
          
          amp = fCellsEM->GetCellAmplitude(id); // GeV
          
          //if(amp < 0.1) continue ;
          
          // 2d projection view
          
          fGeomEM->EtaPhiFromIndex(id,eta,phi);
          
          //            printf("CaloCell %d, ID %d, energy %2.2f,eta %2.2f, phi %2.2f\n",
          //                   icell,id,amp,eta,GetPhi(phi)*TMath::RadToDeg());
          
          fHistoEM->Fill(eta,GetPhi(phi),amp);
          
          // 3d view
          
          fGeomEM->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
          //Gives SuperModule and Tower numbers
          fGeomEM->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
          //Gives label of cell in eta-phi position per each supermodule
          
          fGeomEM->RelPosCellInSModule(id, x, y, z);
          
          // It should not happen, but in case the OCDB file is not the
          // correct one.
          if(iSupMod >= nModules) continue;

          // Push the data to the 3dvisualization tools
          TEveQuadSet* q = smodules[iSupMod];
          if (q) 
          {
            q->AddQuad(y, z);
            q->QuadValue(amp);
          }
        }
      } // cells in cluster
      else // clusters
      {
        //printf("*** EMCAL cluster position ***\n");

        // 2d projection
        fHistoEM->Fill(momentum.Eta(),GetPhi(momentum.Phi()),momentum.E());
        
        // Push the data to the 3d visualization tools
        
        fGeomEM->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta);
        fGeomEM->RelPosCellInSModule(absId, x, y, z);
        
        TEveQuadSet* q = smodules[iSupMod];
        if (q) 
        {
          q->AddQuad(y, z);
          q->QuadValue(momentum.E());
        }
      }
    } // EMCAL
    else // PHOS
    {
      TVector3 xyz;
      Int_t relId[4], module;
      Float_t xCell, zCell;

      // Plot clusters or cluster cells?
      if( plotcells[0] ) // cells
      {       
        //printf("*** PHOS cluster cells ***\n");

        for(Int_t icell = 0; icell < cluster->GetNCells(); icell++)
        {
          id  = cluster->GetCellAbsId(icell);
          
          amp = fCellsPH->GetCellAmplitude(id); // GeV
          
          //if(amp < 0.1) continue ;
          
          amp = fCellsPH->GetCellAmplitude(id); // GeV
                                                //if(amp < 0.1) continue ;
          
          // 2d view
          TVector3 xyz;
          Int_t relId[4], module;
          Float_t xCell, zCell;
          
          fGeomPH->AbsToRelNumbering(id,relId);
          //        printf("PHOS, mod %d\n",module);
          module = relId[0];
          fGeomPH->RelPosInModule(relId,xCell,zCell);
          
          // One should use the following lines but it crashes in the geometry method
          //fGeomPH->Local2Global(module,xCell,zCell,xyz); <<-- It crashes here?
          //fHistoPH->Fill(xyz.Eta(),GetPhi(xyz.Phi()),amp);
          
          //        printf("PHOS CaloCell %d, ID %d, energy %2.3f,eta %f, phi %f\n",
          //               icell,id,amp,xyz.Eta(),xyz.Phi()*TMath::RadToDeg());     
          
          // Temporary fix, one should use the previous lines
          if(id==absId) fHistoPH->Fill(momentum.Eta(),GetPhi(momentum.Phi()),momentum.Energy());
          
          // 3d view
             
          // Push the data to the 3dvisualization tools
          TEveQuadSet* q = smodulesPH[module-1];
          if (q) 
          {
            q->AddQuad(xCell, zCell);
            q->QuadValue(amp);
          }

        } // loop
      } // cells in cluster
      else // cluster
      {
        // printf("*** PHOS cluster position ***\n");

        fHistoPH->Fill(momentum.Eta(),GetPhi(momentum.Phi()),momentum.E());
                
        // 3d view
        
        fGeomPH->AbsToRelNumbering(absId,relId);
        module = relId[0];
        
        printf("\t mod %d\n",module);
        
        fGeomPH->RelPosInModule(relId,xCell,zCell);
            
        // Push the data to the 3dvisualization tools
        TEveQuadSet* q = smodulesPH[module-1];
        if (q) 
        {
          printf("Add quad\n");
          q->AddQuad(xCell, zCell);
          q->QuadValue(amp);
        }
      }
    }
  } // cluster loop
  //} event loop
  
  //==================================================
  // 3d display
  // EMCAL
  for (Int_t sm = 0; sm < nModules; ++sm)
  {
    smodules[sm]->RefitPlex();
  }

//  // PHOS
  for (Int_t sm = 0; sm < 4; ++sm)
  {
    smodulesPH[sm]->RefitPlex();
  }
  
  //==================================================
  // 2d display
  data->AddHistogram(fHistoEM);
  data->RefSliceInfo(0).Setup("EMCell:", 0, kViolet);
  data->AddHistogram(fHistoPH);
  data->RefSliceInfo(1).Setup("PHCell:", 0, kOrange);
  
//  data->AddHistogram(fHistoneg);
//  data->RefSliceInfo(3).Setup("NegCg:", 0, kBlue);
//  data->AddHistogram(fHistopos);
//  data->RefSliceInfo(4).Setup("PosCg:", 0, kRed);
   
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
  CreateProjections(calo3d);
  
  
  gEve->Redraw3D();
  
  //return data;
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
   AliEveEventManager::RegisterTransient(lego);

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
   AliEveEventManager::RegisterTransient(calo3d);
   
   calo3d->SetBarrelRadius(550);
   calo3d->SetEndCapPos(550);
   g_histo2d_s2->AddElement(calo3d);
 
   return calo3d;
}

//______________________________________________________________________________
AliEveMultiView* CreateProjections(TEveCalo3D *calo3d){

   AliEveMultiView *al = AliEveMultiView::Instance();
   al->ImportEventRPhi(calo3d);
   al->ImportEventRhoZ(calo3d);

   return al;
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
  
  if ( absId < 0) return ;
  
  if ( isEMCAL )
  {
    Int_t iTower = -1, iIphi = -1, iIeta = -1;
    fGeomEM->GetCellIndex(absId,iSM,iTower,iIphi,iIeta);
    fGeomEM->GetCellPhiEtaIndexInSModule(iSM,iTower,iIphi, iIeta,irow,icol);
  } // EMCAL
  else // PHOS
  {
    Int_t    relId[4];
    fGeomPH->AbsToRelNumbering(absId,relId);
    irow = relId[2];
    icol = relId[3];
    iSM  = relId[0]-1;
    
  }// PHOS  
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


