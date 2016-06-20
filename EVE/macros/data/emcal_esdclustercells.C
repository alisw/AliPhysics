
/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

///
/// \file emcal_esdclustercells.C
/// \brief event display for EMCal and PHOS
///
/// Two kind of displays:
///  * 3D: PHOS and EMCAL cells/clusters fired in the event are set in boxes TEveQuadSet,
///    the color of the box will depend in the amount of signal
///  * 2D: PHOS and EMCAL cells/clusters eta and phi location fills a histogram, its projection is shown.
///
/// In order to avoid the display of too many bad channels (not completelly possible to remove at this stage)
/// a set of cuts are implemented in this macro:
///   * Cluster energy, see energyCut[].
///   * Number of cells in cluster, see ncellsCut[].
///   * Shower shape, see m02lowCut[] and m20lowCut[]
///   * Exoticity, see exoCut
/// The cuts are implemented as global variables, with 2 entries one set for PHOS and another for EMCal.
///
/// What can be displayed are either the energy deposit in the cells of the clusters or just the total cluster
/// energy deposit and location, assuming that it has 3x3 (EMCAL) or 5x5 (PHOS) cells,
/// which only happens at high energy but for display it might be enough. See plotcells[]
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS
///

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

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliEMCALGeometry.h"
#include "AliPHOSGeometry.h"
#endif

#include <iostream>
using namespace std;

double pi = TMath::Pi();

Double_t GetPhi(Double_t phi);

//
// 2D view stuff
//
TEveViewer *g_histo2d_v   = 0;
TEveScene  *g_histo2d_s   = 0;
TEveScene  *g_histo2d_s2  = 0;
TEveCaloLegoOverlay* g_histo2d_lego_overlay = 0;

TEveCaloLego* CreateHistoLego(TEveCaloData* data);
TEveCalo3D* Create3DView(TEveCaloData* data);
//

//
// Calorimeter methods
//
void FillEMCALClusters(Int_t absIdEMaxCell);
void FillPHOSClusters (Int_t absIdEMaxCell);

void GetMaxEnergyCellAbsId(Int_t & absId, Float_t & eMax);
void GetModuleNumberColAndRow(Int_t absId, Bool_t isEMCAL,
                              Int_t & iSM, Int_t & icol, Int_t & irow);
Float_t GetECross(Bool_t isEMCAL, Int_t imod, Int_t icol, Int_t irow);

Bool_t IsBadCluster(Int_t absId, Float_t eMax);

void SetUpEMCALGeometry(AliESDEvent * esd);
void SetUpPHOSGeometry (AliESDEvent * esd);

void SetEMCALMatrices(AliESDEvent * esd, Bool_t & ok);

void SetUpEMCALQuads();
void SetUpPHOSQuads();

//void AnalyzeTracks(AliESDEvent * esd);
void AnalyzeClusters(AliESDEvent * esd);
//

//
// Declare once different commonly used parameters/histograms
//
TLorentzVector fClusterMomentum;  /// Current cluster kinematics in TLorentzVector
AliESDCaloCluster * fCaloCluster; /// Pointer with current cluster info

Double_t vertex[] = {0.,0.,0.}; /// Vertex container

AliEMCALGeometry * fGeomEM  = 0; /// EMCal geometry
AliPHOSGeometry  * fGeomPH  = 0; /// PHOS geometry

AliESDCaloCells  * fCellsEM = 0; /// List with EMCAL cells
AliESDCaloCells  * fCellsPH = 0; /// List with PHOS cells

Bool_t fGeoSet = kFALSE; /// Check if EMCAL/PHOS geometry, matrices, volumes etc already set once

TH2F* fHistoEM  = 0; /// Histogram with EMCal signals and location
TH2F* fHistoPH  = 0; /// Histogram with PHOS signals and location
//TH2F* fHistoneg = 0;
//TH2F* fHistopos = 0;

Int_t debug = 10;

TGeoNode* fNodeEM    ; /// EMCAL volumes node
TGeoNode* fNodePH[4] ; /// PHOS volumes nodes

TEveQuadSet* fQuadsEMCAL[20]; /// EMCAL list of quads
TEveQuadSet* fQuadsPHOS [4];  /// PHOS list of quads
//

//------------------------------------------------------------------------------
//
// Cluster cuts          { PHOS, EMCAL}
//
Int_t   nMinCellsCut[] = { 3   , 2    };    /// Number of cells in cluster must be larger than this value.
Int_t   nMaxCellsCut[] = { 60  , 30   };    /// Number of cells in cluster must be smaller than this value.
Float_t    energyCut[] = { 0.30, 0.30 };    /// Cluster energy must be larger than this value.
Float_t    m02lowCut[] = { 0.70, 0.10 };    /// Cluster shower shape major axis must be larger than this value.
Float_t    m02higCut[] = { 7.00, 7.00 };    /// Cluster shower shape major axis must be larger than this value.
Float_t    m20lowCut[] = { 0.50, -1.0 };    /// Cluster shower shape lower axis must be larger than this value.
Float_t       exoCut   = 0.95;              /// Reject clusters with this exoticity value.
Bool_t     plotcells[] = {kFALSE , kTRUE};  /// Display cells in clusters or cluster center.
//
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
///
/// Main execution method
///
void emcal_esdclustercells()
{
    if ( debug > 0 ) printf("Execute emcal_esdclustercells\n");
    if ( debug > 9 ) AliLog::SetModuleDebugLevel("EMCAL",100);
  
    //--------------------
    // Access to esd event
    //--------------------
  
    AliESDEvent* esd = AliEveEventManager::AssertESD();
          
    if(esd->GetNumberOfCaloClusters() <= 0 )
    {
      if ( debug > 0 ) printf("emcal_esdclustercells(): No calorimeter info, nclusters 0, skip!\n");
      return;
    }
  
    //-----------------------------------
    // Set geometry, volumes, alignment,
    // 2d histograms once for first event
    //-----------------------------------
    
    if( !fGeoSet )
    {
        AliEveEventManager::AssertGeometry();
        
        SetUpEMCALGeometry(esd);
        
        SetUpPHOSGeometry(esd);
        
        fHistoEM  = new TH2F("histoEMcell","EMCal Cell #eta vs #phi vs E",
                             100,-1.5,1.5,80,-pi,pi);
        fHistoPH  = new TH2F("histoPHcell","PHOS Cell #eta vs #phi vs E",
                             100,-1.5,1.5,80,-pi,pi);
        //    fHistopos = new TH2F("histopos_t","Histo 2d positive",
        //                              100,-1.5,1.5,80,-pi,pi);
        //    fHistoneg = new TH2F("histoneg_t","Histo 2d negative",
        //                              100,-1.5,1.5,80,-pi,pi);
        
        fGeoSet = kTRUE ;
    }
    else
    {
        fHistoEM ->Reset();
        fHistoPH ->Reset();
        //    fHistopos->Reset();
        //    fHistoneg->Reset();
    }
  
    // Do here the setting of the EMCal matrices, 
    // in case the first event did not have them in the ESD
    // Do it once.
    // Get first EMCal/DCal SM matrix in geometry if non null skip.  
    if ( !fGeomEM->GetMatrixForSuperModuleFromArray(0) || 
         !fGeomEM->GetMatrixForSuperModuleFromArray(12) ) 
    {
      Bool_t ok = kFALSE;

      SetEMCALMatrices(esd,ok);
      
      if(!ok) 
      {
        printf("emcal_esdclustercells: Alignment Matrices not available, skip!\n");
        return;
      }
    }
  
    // Matrix debugging
    if ( debug > 9 )
    {
      printf("***>>> Begin EMCal ESD alignment matrices ***<<<\n");
    
      for(Int_t mod = 0; mod < fGeomEM->GetNumberOfSuperModules(); mod++)
      {
       if ( esd->GetEMCALMatrix(mod) ) 
         esd->GetEMCALMatrix(mod)->Print();
       else 
         printf("No esd matrix for SM %d\n",mod);
      }
      
      printf("***>>> End EMCal ESD alignment matrices ***<<<\n\n");
    
      printf("***>>> Begin EMCal Geometry alignment matrices ***<<<\n");
    
      for(Int_t mod = 0; mod < fGeomEM->GetNumberOfSuperModules(); mod++)
      {
        if ( fGeomEM->GetMatrixForSuperModule(mod) )
          fGeomEM->GetMatrixForSuperModule(mod)->Print();
        else 
          printf("No geo matrix for SM %d\n",mod);
      }
      
      printf("***>>> End EMCal Geometry alignment matrices ***<<<\n\n");
    
      printf("***>>> Begin EMCal OCDB alignment matrices ***<<<\n");
    
      AliCDBManager* man = AliCDBManager::Instance();
      AliCDBEntry *cdb = (AliCDBEntry*)  man->Get("EMCAL/Align/Data");
    
      cdb->Print();
    
      printf("***>>> End EMCal OCDB alignment matrices ***<<<\n\n");      
    }
  
    //---------------------------
    // Set up 3d quads
    //---------------------------

    SetUpEMCALQuads();
    
    SetUpPHOSQuads();
    
    //-------------------------
    // Data analysis
    //-------------------------
    
    // Tracks
    //
    //AnalyzeTracks(esd);
    
    // Calorimeter
    //
    AnalyzeClusters(esd);

    //---------------------------
    // 2d display
    //---------------------------
    
    TEveCaloDataHist* data = new TEveCaloDataHist();
    AliEveEventManager *manager = AliEveEventManager::Instance();
    cout<<"\n\n adding 2d data to event manager"<<endl;
    manager->RegisterTransient(data);
    
    data->AddHistogram(fHistoEM);
    data->RefSliceInfo(0).Setup("EMCell:", 0, kOrange+7);
    data->AddHistogram(fHistoPH);
    data->RefSliceInfo(1).Setup("PHCell:", 0, kYellow);
    
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
//    CreateHistoLego(data); // this function breaks Event Display closing

    // Plotting the 3D histogram and projections RPhi and RhoZ
    TEveCalo3D *calo3d = Create3DView(data);
    
    //-----------------
    // Send to EVE
    //-----------------
    
    gEve->Redraw3D();
    
}

////______________________________________________________________________________
/////
///// Fill 2D histograms for tracks
/////
/////
//void AnalyzeTracks(AliESDEvent* esd)
//{
//  //printf("Number of Tracks %d\n",esd->GetNumberOfTracks());
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
//}

//______________________________________________________________________________
///
/// Access calorimeter clusters and cells, loop over the clusters
/// and fill the 2d histograms or the Quads for the 3d view.
/// Clusters will be selected depending on different quality criteria in method
/// IsBadCluster.
///
/// \param esd ESD event pointer
///
void AnalyzeClusters(AliESDEvent * esd)
{
    // Extract EMCAL/PHOS cell-cluster information from the ESDs
    //
    fCellsEM = esd->GetEMCALCells();
    fCellsPH = esd->GetPHOSCells();
    
    Int_t ncellEM = 0; 
    if(fCellsEM) fCellsEM->GetNumberOfCells() ;
    Int_t ncellPH = 0; 
    if(fCellsPH) fCellsPH->GetNumberOfCells() ;
    
    if ( debug > 0 ) 
    {
      printf("Number of ESD CaloCells: EM %d PH %d\n",ncellEM,ncellPH);
      printf("Number of ESD CaloClusters: %d\n",esd->GetNumberOfCaloClusters());
    }
  
    Int_t   absIdEMaxCell = -1;
    Float_t eMaxCell      = 0.;
    
    for (Int_t iclus =  0; iclus < esd->GetNumberOfCaloClusters(); iclus++)
    {
        fCaloCluster = esd->GetCaloCluster(iclus);
        
        fCaloCluster->GetMomentum(fClusterMomentum,vertex);
        
        absIdEMaxCell = -1;
        eMaxCell      = 0.;
        GetMaxEnergyCellAbsId(absIdEMaxCell, eMaxCell) ;
        
        if( IsBadCluster(absIdEMaxCell, eMaxCell) ) continue ;
        
        if ( debug > 0 ) 
        {
          printf("Cluster %d, is EMCAL %d, E %2.2f, eta %2.2f, phi %2.2f\n",
                 iclus, fCaloCluster->IsEMCAL(),fClusterMomentum.E(),
                 fClusterMomentum.Eta(),fClusterMomentum.Phi()*TMath::RadToDeg());
        }
      
        // Plot clusters or cells in clusters
        //
        if(fCaloCluster->IsEMCAL()) FillEMCALClusters(absIdEMaxCell);
        else                        FillPHOSClusters (absIdEMaxCell);
    } // cluster loop
}

//______________________________________________________________________________
///
/// Plot EMCal clusters or cluster cells
/// Fill the histograms for the 2D or the Quads for the 3D
///
/// \param ID number of cell with highest energy in the cluster
///
void FillEMCALClusters(Int_t absIdEMaxCell)
{  
    if(!fCellsEM) return;
  
    // If any of the matrices is missing, skip analysis 
    for(Int_t mod = 0; mod < fGeomEM->GetNumberOfSuperModules(); mod++)
    {
      if ( !fGeomEM->GetMatrixForSuperModule(mod) )  
        printf("emcal_esdclustercells.C::FillEMCALClusters() No geo matrix for SM %d, skip this event for EMCal!!!\n",mod);
    }
  
    Double_t x=0., y=0., z=0.;
    Int_t iSupMod =  -1 ;
    Int_t iTower  =  -1 ;
    Int_t iIphi   =  -1 ;
    Int_t iIeta   =  -1 ;
  
    if( plotcells[1] ) // cells in cluster
    {
        if ( debug > 1 ) printf("\t *** EMCAL cluster cells ***\n");
        
        Float_t amp   = -1 ;
        Int_t   id    = -1 ;
        Float_t phi   =  0 ;
        Float_t eta   =  0 ;
        Int_t iphi    = -1 ;
        Int_t ieta    = -1 ;
        
        for(Int_t icell = 0; icell < fCaloCluster->GetNCells(); icell++)
        {
            id  = fCaloCluster->GetCellAbsId(icell);
            
            amp = fCellsEM->GetCellAmplitude(id); // GeV
            
            //if(amp < 0.1) continue ;
            
            //
            fGeomEM->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
            //Gives SuperModule and Tower numbers
            fGeomEM->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
            //Gives label of cell in eta-phi position per each supermodule

            // 2d projection view
            
            fGeomEM->EtaPhiFromIndex(id,eta,phi);
            
            if ( debug > 1 ) 
              printf("\t CaloCell %d, ID %d, energy %2.2f,eta %2.2f, phi %2.2f\n",
                     icell,id,amp,eta,GetPhi(phi)*TMath::RadToDeg());
            
            if(TMath::Abs(eta) < 0.7)
            {
              fHistoEM->Fill(eta,GetPhi(phi),amp);
            }
            else // Workaround
            {
              printf("Wrong eta value for EMCal clusters, active workaround!!!!");
              Float_t etaCluster = fClusterMomentum.Eta();  
              Float_t phiCluster = fClusterMomentum.Phi();  
              
              Int_t iSupModMax =  -1 ;
              Int_t iTowerMax  =  -1 ;
              Int_t iIphiMax   =  -1 ;
              Int_t iIetaMax   =  -1 ;
              Int_t iphiMax    =  -1 ;
              Int_t ietaMax    =  -1 ;
              
              fGeomEM->GetCellIndex(absIdEMaxCell,iSupModMax,iTowerMax,iIphiMax,iIetaMax);
              //Gives SuperModule and Tower numbers
              fGeomEM->GetCellPhiEtaIndexInSModule(iSupModMax,iTowerMax,iIphiMax, iIetaMax,iphiMax,ietaMax);
              //Gives label of cell in eta-phi position per each supermodule

              Int_t iPhiDiff = iphiMax - iphi;
              Int_t iEtaDiff = ietaMax - ieta;
              
              Float_t etaWA = etaCluster + iEtaDiff*0.014;
              Float_t phiWA = phiCluster + iPhiDiff*0.014;
              printf("\t \t Workaround Cell eta = (std %2.2f, wa %2.2f); phi = (std %2.2f, wa %2.2f)\n",eta,etaWA,phi,phiWA);
              
              fHistoEM->Fill(etaWA,GetPhi(phiWA),amp);
            }
          
            // 3d view
            
            fGeomEM->RelPosCellInSModule(id, x, y, z);
            
            if ( debug > 1 ) 
              printf("\t , SM %d (Nodes %d),  iEta %d,  iPhi %d, x %3.3f, y %3.3f, z %3.3f \n",
                     iSupMod,fNodeEM->GetNdaughters(),ieta,iphi,x,y,z);

          
            // It should not happen, but in case the OCDB file is not the
            // correct one.
            if(iSupMod >= fNodeEM->GetNdaughters()) continue;
            
            // Push the data to the 3dvisualization tools
            //
            if (fQuadsEMCAL[iSupMod])
            {
                fQuadsEMCAL[iSupMod]->AddQuad(y, z);
                fQuadsEMCAL[iSupMod]->QuadValue(amp*1000);
            }
        }
    } // cells in cluster
    else // clusters
    {
        // 2d projection
        //
        fHistoEM->Fill(fClusterMomentum.Eta(),GetPhi(fClusterMomentum.Phi()),fClusterMomentum.E());
        
        // 3d view
        //
        fGeomEM->GetCellIndex(absIdEMaxCell,iSupMod,iTower,iIphi,iIeta); // needed to get iSupMod
        
        fGeomEM->RelPosCellInSModule(absIdEMaxCell, x, y, z);
        
        Float_t phi = fClusterMomentum.Phi()*TMath::RadToDeg();
        if ( phi < 0 ) phi+=360;
      
        if(debug > 1)
          printf("\t *** EMCAL cluster position, eta = %f, phi = %f, E = %f, SM = %d, location (x=%f, y=%f, z=%f) ***\n",
                 fClusterMomentum.Eta(),phi,fClusterMomentum.E(),iSupMod, x, y ,x);
        
        // Push the data to the 3d visualization tools
        //
        if (fQuadsEMCAL[iSupMod])
        {
            fQuadsEMCAL[iSupMod]->AddQuad(y, z);
            fQuadsEMCAL[iSupMod]->QuadValue(fClusterMomentum.E()*1000);
        }
    }
}

//______________________________________________________________________________
///
/// Plot EMCal clusters or cluster cells
/// Fill the histograms for the 2D or the Quads for the 3D
///
/// \param ID number of cell with highest energy in the cluster
///
void FillPHOSClusters(Int_t absIdEMaxCell)
{
    if(!fCellsPH) return;

    TVector3 xyz;
    Int_t relId[4], module;
    Float_t xCell, zCell;
    
    // Plot clusters or cluster cells?
    if( plotcells[0] ) // cells
    {
      if ( debug > 0 ) printf("\t *** PHOS cluster cells ***\n");
        
        Float_t amp   = -1 ;
        Int_t   id    = -1 ;
        
        for(Int_t icell = 0; icell < fCaloCluster->GetNCells(); icell++)
        {
            id  = fCaloCluster->GetCellAbsId(icell);
            
            amp = fCellsPH->GetCellAmplitude(id); // GeV
            
            //if(amp < 0.1) continue ;
            
            amp = fCellsPH->GetCellAmplitude(id); // GeV
            //if(amp < 0.1) continue ;
            
            
            TVector3 xyz;
            Int_t relId[4], module;
            Float_t xCell, zCell;
            
            fGeomPH->AbsToRelNumbering(id,relId);
            //        printf("PHOS, mod %d\n",module);
            
            module = relId[0];
            fGeomPH->RelPosInModule(relId,xCell,zCell);
            
            // 2d view
            //
            
            // -->>One should use the following lines but it crashes in the geometry method<<--
            //
            //fGeomPH->Local2Global(module,xCell,zCell,xyz); <<-- It crashes here?
            //fHistoPH->Fill(xyz.Eta(),GetPhi(xyz.Phi()),amp);
            
            if ( debug > 1 ) 
            printf("\t PHOS CaloCell %d, ID %d, energy %2.3f,eta %f, phi %f\n",
                    icell,id,amp,xyz.Eta(),xyz.Phi()*TMath::RadToDeg());
            
            // -->>Temporary fix, one should use the previous lines<<--
            if(id == absIdEMaxCell) fHistoPH->Fill(fClusterMomentum.Eta(),GetPhi(fClusterMomentum.Phi()),fClusterMomentum.Energy());
            
            // Push the data to the 3d visualization tools
            //
            if (fQuadsPHOS[module-1])
            {
                fQuadsPHOS[module-1]->AddQuad(xCell, zCell);
                fQuadsPHOS[module-1]->QuadValue(amp*1000);
            }
            
        } // loop
    } // cells in cluster
    else // cluster
    {
        // 2d view
        //
        fHistoPH->Fill(fClusterMomentum.Eta(),GetPhi(fClusterMomentum.Phi()),fClusterMomentum.E());
        
        // 3d view
        //
        fGeomPH->AbsToRelNumbering(absIdEMaxCell,relId);
        module = relId[0];
        
        if(debug > 1)
        {
          Float_t phi = fClusterMomentum.Phi()*TMath::RadToDeg();
          if ( phi < 0 ) phi+=360;
        
          printf("*** PHOS cluster position, eta = %f, phi = %f, E = %f, mod %d, location (x=%f, z=%f) ***\n",
               fClusterMomentum.Eta(),phi,fClusterMomentum.E(),module,xCell,zCell);
        }
      
        fGeomPH->RelPosInModule(relId,xCell,zCell);
        
        // Push the data to the 3d visualization tools
        //
        if (fQuadsPHOS[module-1])
        {
            fQuadsPHOS[module-1]->AddQuad(xCell, zCell);
            fQuadsPHOS[module-1]->QuadValue(fClusterMomentum.E()*1000);
        }
    }
}

//______________________________________________________________________________
///
/// \return sum of energy in neighbor cells (cross) to the reference cell
///
/// \param isEMCAL bool true for EMCAL, false for PHOS
/// \param imod reference cell module
/// \param icol reference cell column
/// \param irow reference cell row
///
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

//_____________________________________________________________________
///
/// Get for a given cell absolute ID the (super)module number, column and row
///
/// \param absId absolute cell identity number
/// \param isEMCAL bool true for EMCAL, false for PHOS
/// \param iSM reference cell module
/// \param icol reference cell column
/// \param irow reference cell row
///
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

//______________________________________________________________________________
///
/// Get the cell with highest energy in the cluster and its energy
///
/// \param absId absolute cell identity number with highest energy in cluster
/// \param eMax energy of most significant cell in cluster
///
void GetMaxEnergyCellAbsId(Int_t & absId, Float_t & eMax)
{
    Double_t eCell       =-1.;
    Int_t    cellAbsId   =-1 ;
    Int_t    iSupMod     =-1 ;
    
    AliESDCaloCells*  cells = 0;
    
    if(fCaloCluster->IsEMCAL()) cells = fCellsEM;
    else                        cells = fCellsPH;
    
    for (Int_t iDig=0; iDig < fCaloCluster->GetNCells(); iDig++)
    {
        cellAbsId = fCaloCluster->GetCellAbsId(iDig);
        
        eCell  = cells->GetCellAmplitude(cellAbsId);
        //printf("\t \t \t ecell %f\n",eCell);
        if(eCell > eMax)
        {
            eMax  = eCell;
            absId = cellAbsId;
        }
    }// cell loop
}

//______________________________________________________________________________
///
/// Check goodness of clusters
///
/// \return true for bad clusters
///
/// \param absId   ID number of cell with highest energy in cluster
/// \param eMax    Energy of cell with highest energy in cluster
///
Bool_t IsBadCluster(Int_t absId, Float_t eMax)
{
    Bool_t isEMCAL = fCaloCluster->IsEMCAL();
    
    //  if(!isEMCAL) printf("PHOS cluster  \n");
    //  else         printf("EMCAL cluster \n");
    //
    //  printf("n cells %d, E %f, m02 %f, m20 %f\n",
    //         fCaloCluster->GetNCells(),fCaloCluster->E(),fCaloCluster->GetM02(),fCaloCluster->GetM20());
    
    if(fCaloCluster->GetNCells() < nMinCellsCut[isEMCAL]) return kTRUE ;
    if(fCaloCluster->GetNCells() > nMaxCellsCut[isEMCAL]) return kTRUE ;
    
    //printf("n Cells OK\n");
    
    if(fCaloCluster->E()         < energyCut[isEMCAL]) return kTRUE ;
    
    //printf("E OK\n");
    
    if(fCaloCluster->GetM02()    < m02lowCut[isEMCAL]) return kTRUE ;
    if(fCaloCluster->GetM02()    > m02higCut[isEMCAL]) return kTRUE ;
    
    //printf("M02 OK\n");
    
    if(fCaloCluster->GetM20()    < m20lowCut[isEMCAL]) return kTRUE ;
    
    //printf("M20 OK\n");
    
    if ( eMax < energyCut[isEMCAL]/2. || absId < 0 )
    {
        //printf("\t Remove Emax %f, absId %d\n",eMax,absId);
        return kTRUE;
    }
    
    Int_t ism  =  -1, icol = -1, irow = -1;
    GetModuleNumberColAndRow(absId, isEMCAL, ism,icol,irow)  ;
    
    Float_t eCross = GetECross(isEMCAL, ism, icol, irow);
    
    //printf("\t SM %d, Id %d, Emax %2.2f, Ecross %2.2f, fraction %2.2f\n",ism,absId,eMax,eCross,1-eCross/eMax);
    
    if(1-eCross/eMax > exoCut)
    {
        //printf("\t \t Remove SM %d, Id %d, EMax %2.2f, Ecross %2.2f, fraction %2.2f\n",ism,absId,eMax,eCross,1-eCross/eMax);
        return kTRUE;
    }
    
    //printf("Cluster accepted!\n");
    
    return kFALSE;
}

//______________________________________________________________________________
///
/// Get the EMCal geometry, alignment matrices and volumes.
/// Initialize data members fGeomEM and fNodeEM.
/// Set it at least for the first event.
///
/// \param esd ESD event pointer needed for matrix retrieval
///
void SetUpEMCALGeometry(AliESDEvent * esd)
{
    //
    // Set the geometry
    //
    if(!fGeomEM)
        fGeomEM  = AliEMCALGeometry::GetInstance();
    
    if (!fGeomEM)
    {
        printf("xxx Set EMCal default geo as Run2 xxx\n");
        fGeomEM  = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
    }
    
    Bool_t ok = kFALSE;
    SetEMCALMatrices(esd,ok); // Do it outside also if we could not set them here.
  
    if(!ok) 
    {
      printf("emcal_esdclustercells::SetUpEMCALGeometry() - Alignment Matrices not available, skip!");
      return;
    }    //
    // EMCAL volumes
    //
    fNodeEM = gGeoManager->GetTopVolume()->FindNode("XEN1_1");
    
    //Int_t nModules = fNodeEM->GetNdaughters();
    
    // Check that the EMCAL geo and the nodes from EMCAL have the same number of entries
    //
    if(fNodeEM->GetNdaughters() != fGeomEM->GetNumberOfSuperModules())
        printf("*** === EMCAL DIGITS - N Daughter Nodes %d - N super mod %d === ***\n",
               fNodeEM->GetNdaughters(), fGeomEM->GetNumberOfSuperModules());
    
    //  // Get the EMCAL bounding boxes for the super modules.
    //  // 4 kind of SM: 10 Full EMCal, 2 1/3 EMCal, 6 DCal (2/3 EMCal) and 2 1/3 EMCal in DCal region.
    //  //
    //  TGeoBBox* bbbox = (TGeoBBox*) node->GetDaughter(0) ->GetVolume()->GetShape();
    //  TEveFrameBox* frame_big = new TEveFrameBox();
    //  frame_big->SetFrameColorRGBA(200,200,0,50);
    //  frame_big->SetAABoxCenterHalfSize(0, 0, 0, bbbox->GetDX(), bbbox->GetDY(), bbbox->GetDZ());
    //
    //  TEveFrameBox* frame_sml  = 0x0;
    //  TEveFrameBox* frame_dcl  = 0x0;
    //  TEveFrameBox* frame_smld = 0x0;
    //
    //  if (nModules > 10)
    //  {
    //    TGeoBBox* sbbox = (TGeoBBox*) node->GetDaughter(10)->GetVolume()->GetShape();
    //    frame_sml = new TEveFrameBox();
    //    frame_sml->SetFrameColorRGBA(200,200,0,50);
    //    frame_sml->SetAABoxCenterHalfSize(0, 0, 0, sbbox->GetDX(), sbbox->GetDY(), sbbox->GetDZ());
    //  }
    //
    //  if (nModules > 12)
    //  {
    //    TGeoBBox* dbbox = (TGeoBBox*) node->GetDaughter(12)->GetVolume()->GetShape();
    //    frame_dcl = new TEveFrameBox();
    //    frame_dcl->SetFrameColorRGBA(200,200,0,50);
    //    frame_dcl->SetAABoxCenterHalfSize(0, 0, 0, dbbox->GetDX(), dbbox->GetDY(), dbbox->GetDZ());
    //
    //    TGeoBBox* sdbbox = (TGeoBBox*) node->GetDaughter(18)->GetVolume()->GetShape();
    //    frame_smld = new TEveFrameBox();
    //    frame_smld->SetFrameColorRGBA(200,200,0,50);
    //    frame_smld->SetAABoxCenterHalfSize(0, 0, 0, sdbbox->GetDX(), sdbbox->GetDY(), sdbbox->GetDZ());
    //  }
}



///
/// Set the geometry alignment matrices from the ones stored in the data.
/// Ideally it should be done just in first event, but it has been observed
/// that the first event does not always contain them (???)
///
void SetEMCALMatrices(AliESDEvent * esd, Bool_t & ok)
{  
  ok = kTRUE;
  // Set all the matrices
  for(Int_t mod = 0; mod < fGeomEM->GetNumberOfSuperModules(); mod++)
  {
    if( debug > 1 ) printf("Load EMCAL ESD matrix %d, %p\n",mod,esd->GetEMCALMatrix(mod));
    
    if( esd->GetEMCALMatrix(mod) )
    {
      fGeomEM->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
    }
    else // set default identity matrix
    {
      printf("Could not set EMCal geo matrix for SM %d",mod);
      ok = kFALSE;
      //fGeomEM->SetMisalMatrix((new TGeoHMatrix),mod) ;
    }
  } // loop over super modules
  
  if(ok) return;
  
  // Re-set the bool
  ok = kTRUE;
  // Try to get now the matrices from OCDB
  printf("Get alignment matrices form OCDB\n");
  
  AliCDBManager* man = AliCDBManager::Instance();
  //man->SetDefaultStorage("raw://");
  man->SetRun(esd->GetRunNumber());
  
  AliCDBEntry *cdb = (AliCDBEntry*)  man->Get("EMCAL/Align/Data");
  TClonesArray * matrixArr = (TClonesArray*) cdb->GetObject();
  
  for(Int_t mod = 0; mod < fGeomEM->GetNumberOfSuperModules(); mod++)
  {
    if( debug > 1 ) printf("Load EMCAL OCDB matrix %d, %p\n",mod,matrixArr->At(mod));
    
    if( matrixArr->At(mod) )
    {
      fGeomEM->SetMisalMatrix(((TGeoHMatrix *) matrixArr->At(mod)),mod) ;
    }
    else // set default identity matrix
    {
      printf("Could not set EMCal geo matrix for SM %d from OCDB",mod);
      ok = kFALSE;
      //fGeomEM->SetMisalMatrix((new TGeoHMatrix),mod) ;
    }
  } // loop over super modules
  printf("ok? %d\n",ok);
  
  for(Int_t mod = 0; mod < fGeomEM->GetNumberOfSuperModules(); mod++) 
    printf("Matrix in geometry: imod %d, %p\n",mod,fGeomEM->GetMatrixForSuperModuleFromArray(mod));

  
}

//______________________________________________________________________________
///
/// Get the PHOS geometry, alignment matrices and volumes.
/// Initialize data members fGeomPH and fNodePH.
/// Set it at least for the first event.
///
/// \param esd ESD event pointer needed for matrix retrieval
///
void SetUpPHOSGeometry(AliESDEvent * esd)
{
    //
    // Set the geometry
    //
    if(!fGeomPH)
        fGeomPH  = AliPHOSGeometry::GetInstance();
    
    if (!fGeomPH)
    {
        printf("xxx Set PHOS default geo as Run2 xxx\n");
        fGeomPH  = AliPHOSGeometry::GetInstance("Run2");
    }
    
    //
    // Set the geometry alignment matrices from the ones stored in the data.
    //
    for(Int_t mod = 0; mod < 5; mod++)
    {
        //printf("Load PHOS ESD matrix %d, %p\n",mod,esd->GetPHOSMatrix(mod));
        if(esd->GetPHOSMatrix(mod))
            fGeomPH->SetMisalMatrix(esd->GetPHOSMatrix(mod),mod) ;
        else // set default identity matrix
            fGeomPH->SetMisalMatrix((new TGeoHMatrix),mod) ;
    }// loop over modules
    
    //
    // PHOS volumes
    //
    for (Int_t mod = 1; mod < 5; ++mod)
    {
        fNodePH[mod-1] = gGeoManager->GetTopVolume()->FindNode(Form("PHOS_%d",mod));
    }
    
}


//______________________________________________________________________________
///
/// Set here the elements for the EMCAL 3D view
/// The basic unit are squares with the size of an EMCAL tower or a combination or towers
/// depending on the plotting of towers or clusters
/// Initialize data members fQuadsEMCAL.
/// Set it at least for the first event.
///
///
void SetUpEMCALQuads()
{
    // Define EVE stuff for EMCAL 3d view
    //
    TEveElementList* l = new TEveElementList("EMCAL");
    l->SetTitle("Tooltip");
    gEve->AddElement(l);
    
    //
    // Define energy range for the color palette
    //
    Int_t maxEMCalE = 10000; // MeV
    if(plotcells[1]) maxEMCalE = 2000; // MeV
    
    Int_t minEMCalE = Int_t(energyCut[1]*1000); // MeV
    if(plotcells[1]) minEMCalE = 100; // MeV
    
    //printf("EMCAL Energy limit: max %d MeV, min %d MeV\n",minEMCalE, maxEMCalE);
    
    gStyle->SetPalette(1, 0);
    TEveRGBAPalette* pal = new TEveRGBAPalette(minEMCalE, maxEMCalE);
    //pal->SetLimits(0, 1024);
    
    //
    // Here we will store the EMCAL data that will be treated by EVE
    // per each super-module.
    //
    
    memset(fQuadsEMCAL,0,fGeomEM->GetNumberOfSuperModules()*sizeof(TEveQuadSet*));
    
    // Quad size
    Float_t quadSizeEMCal = 6  ; // cm, tower side size
    if(!plotcells[1]) quadSizeEMCal *=3; // in case of plotting just clusters, increase to 3x3
    
    for (Int_t sm = 0; sm < fNodeEM->GetNdaughters(); ++sm)
    {
        fQuadsEMCAL[sm] = new TEveQuadSet(Form("SM %d", sm+1));
        fQuadsEMCAL[sm]->SetOwnIds(kTRUE);
        
        // Type of object to be displayed, rectangle with cell size
        fQuadsEMCAL[sm]->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
        
        //    fQuadsEMCAL[sm]->SetDefWidth (fGeomEM->GetPhiTileSize());
        //    fQuadsEMCAL[sm]->SetDefHeight(fGeomEM->GetEtaTileSize());
        fQuadsEMCAL[sm]->SetDefWidth (quadSizeEMCal);
        fQuadsEMCAL[sm]->SetDefHeight(quadSizeEMCal);
        
        //printf(">>>> Tile size eta %f  phi %f, set %f <<<<\n", fGeomEM->GetEtaTileSize(),fGeomEM->GetPhiTileSize(), quadSizeEMCal);
        
        fQuadsEMCAL[sm]->RefMainTrans().SetFrom(*fNodeEM->GetDaughter(sm)->GetMatrix());
        
        //    if     (sm < 10) q->SetFrame(frame_big );
        //    else if(sm < 12) q->SetFrame(frame_sml );
        //    else if(sm < 18) q->SetFrame(frame_dcl );
        //    else if(sm < 20) q->SetFrame(frame_smld);
        
        fQuadsEMCAL[sm]->SetPalette(pal);
        
        fQuadsEMCAL[sm]->RefitPlex();
        
        gEve->AddElement(fQuadsEMCAL[sm], l);
    }
}

//______________________________________________________________________________
///
/// Set here the elements for PHOS the 3D view
/// The basic unit are squares with the size of an PHOS crystal or a combination or crystals
/// depending on the plotting of towers or clusters
/// Initialize data members fQuadsPHOS.
/// Set it at least for the first event.
///
///
void SetUpPHOSQuads()
{
    // Define EVE stuff for PHOS 3d view
    //
    TEveElementList* lPH = new TEveElementList("PHOS");
    lPH->SetTitle("Tooltip");
    gEve->AddElement(lPH);
    
    
    // Define energy range for the color palette
    Int_t maxPHOSE = 10000; // MeV
    if(plotcells[0]) maxPHOSE = 2000; // MeV
    
    Int_t minPHOSE = Int_t(energyCut[0]*1000); // MeV
    if(plotcells[1]) minPHOSE = 100; // MeV
    
    //printf("PHOS Energy limit: max %d MeV, min %d MeV\n",minPHOSE, maxPHOSE);
    
    gStyle->SetPalette(1, 0);
    TEveRGBAPalette* palPH = new TEveRGBAPalette(minPHOSE, maxPHOSE);
    //palPH->SetLimits(0, 1024);
    
    // Here we will store the PHOS data that will be treated by EVE
    // per each module.
    //
    TEveQuadSet* smodulesPH[4];
    memset(smodulesPH,0,4*sizeof(TEveQuadSet*));
    
    // Quad size
    Float_t quadSizePHOS = 2.2; // cm, crystal side size
    if(!plotcells[0]) quadSizePHOS*=5; // in case of plotting just clusters, increase to 5x5
    
    for (Int_t mod = 0; mod < 4; ++mod)
    {
        if(!fNodePH[mod])
        {
            printf("PHOS node %d not found\n",mod);
            continue;
        }
        
        fQuadsPHOS[mod] = new TEveQuadSet(Form("Mod %d", mod));
        fQuadsPHOS[mod]->SetOwnIds(kTRUE);
        
        // Type of object to be displayed, rectangle with cell size
        fQuadsPHOS[mod]->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
        fQuadsPHOS[mod]->SetDefWidth (quadSizePHOS);
        fQuadsPHOS[mod]->SetDefHeight(quadSizePHOS);
        //printf("quad size phos %f\n",quadSizePHOS);
        
        
        fQuadsPHOS[mod] ->RefMainTrans().SetFrom(*fNodePH[mod]->GetMatrix());
        //printf("%p\n",*nodePH->GetMatrix());
        
        // Get the PHOS bounding boxes for the modules.
        //
        //    TGeoBBox* bbboxPH = (TGeoBBox*) nodePH->GetDaughter(0)->GetVolume()->GetShape();
        //    TEveFrameBox* frame_PH = new TEveFrameBox();
        //    frame_PH->SetFrameColorRGBA(200,200,0,50);
        //    frame_PH->SetAABoxCenterHalfSize(0, 0, 0, bbboxPH->GetDX(), bbboxPH->GetDY(), bbboxPH->GetDZ());
        
        //    fQuadsPHOS[mod] ->SetFrame(frame_PH);
        
        fQuadsPHOS[mod]->SetPalette(palPH);
        
        fQuadsPHOS[mod]->RefitPlex();
        
        gEve->AddElement(fQuadsPHOS[mod] , lPH);
    }
}

//______________________________________________________________________________
///
/// Shift phi angle for values larger than pi to -pi
///
Double_t GetPhi(Double_t phi)
{
    if (phi > pi)
        phi -= 2*pi;
    
    return phi;
}

//______________________________________________________________________________
TEveCaloLego* CreateHistoLego(TEveCaloData* data)
{
    TGLViewer* glv;
    
    // Viewer initialization, tab creation
    if ( g_histo2d_v == 0 )
    {
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
    } else
    {
        glv = g_histo2d_v->GetGLViewer(); 
    }
    
    // plotting histo
    TEveCaloLego* lego = new TEveCaloLego(data);
    
    g_histo2d_s->AddElement(lego);
    AliEveEventManager *manager = AliEveEventManager::Instance();
    cout<<"\n\nadding lego to event manager"<<endl;
    manager->RegisterTransient(lego);
    
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
TEveCalo3D* Create3DView(TEveCaloData* data)
{  
    // initialization
    if ( g_histo2d_s2 == 0 ) 
    {
        g_histo2d_s2 = gEve->SpawnNewScene("3D Histogram", "3D Histogram");
        gEve->GetDefaultViewer()->AddScene(g_histo2d_s2);
        AliEveMultiView::Instance()->Get3DView()->AddScene(g_histo2d_s2);
        g_histo2d_s2->SetElementName("3D Histogram Scene");
    }
    
    TEveCalo3D* calo3d = new TEveCalo3D(data);
    AliEveEventManager *manager = AliEveEventManager::Instance();
    cout<<"Adding 3d view to event manager"<<endl;
    manager->RegisterTransient(calo3d);
    
    calo3d->SetBarrelRadius(600);
    calo3d->SetEndCapPos(550);
    calo3d->SetMaxTowerH(300);
    calo3d->SetFrameTransparency(100);
    g_histo2d_s2->AddElement(calo3d);
    
    AliEveMultiView::Instance()->ImportEvent(calo3d);
    
    return calo3d;
}

