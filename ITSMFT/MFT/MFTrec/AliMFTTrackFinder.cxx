/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

//============================================================================================
//
//      Class for finding Tracks from Cluster of the ALICE Muon Forward Tracker
//
//      Contact author: raphael.tieulent@cern.ch
//
//============================================================================================

#include <fstream>
#include <iostream>

#include "TMath.h"
#include "TF1.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMFTCAHit.h"
#include "AliMFTCARoad.h"
#include "AliMFTCACell.h"
#include "AliMFTCALayer.h"
#include "AliMFTCATrack.h"
#include "AliMFTCluster.h"
#include "AliMFTTrackFinder.h"


ClassImp(AliMFTTrackFinder)

//___________________________________________________________________________
AliMFTTrackFinder::AliMFTTrackFinder() :
TObject(),
fXCut(0.0025),
fYCut(0.0025),
fMaxSegAngle(180.-15.),
fCellGID(0),
fMaxCellStatus(0),
fNlayers(fNDetMax),
fNtracks(0),
fCalcVertex(kFALSE),
fZVertCalc(0.),
fZVertDet(0.),
fZVertRange(),
fMinTrackLength(fNDetMax),
fPlaneDetEff(),
fAddNoise(kTRUE),
fPixelNoise(1.E-5),
fAddQED(kFALSE),
fMBRate(0.),
fCMOSIntTime(0.),
fThick(0.004),
fReadGeom(kTRUE),
fUseTF(kFALSE),
fDebug(0),
fCPUTime(0.),
fRealTime(0.),
fNDifTracks(0),
fPlanesZ(),
fZGap(1.0),
fNRoads(0),
fErrX(0.),
fErrY(0.)
{
  
#ifdef OLDGEOM
  fZGap = 0.2;
#else
  fZGap = 1.0;
#endif
  
  fZVertRange[0] = fZVertRange[1] = 0.;
  
  fHistList = new TList();
  
  fTracks = new TClonesArray("AliMFTCATrack", 1000);
  
  fRoads = new TClonesArray("AliMFTCARoad", 1000);
  
  for (Int_t iL = 0; iL < fNlayers; iL++) {
    
    fACutV[iL] = 0.20;
    fACutN[iL] = 0.10;
    
    fLayers[iL] = new AliMFTCALayer();
    fLayers[iL]->SetID(iL);
    fPlaneDetEff[iL] = 1.00;
    fPlanesZ[iL] = 0.;
    
    hDA[iL] = new TH1F(Form("hDA%d",iL),Form("#Delta #theta between two segments in Layer %d; #Delta #theta (deg)",iL),200,0.,0.5);
    hDAv[iL] = new TH1F(Form("hDAv%d",iL),Form("#Delta #theta with respect to the vertex in Layer %d; #Delta #theta (deg)",iL),200,0.,3.0);
    fHistList->Add(hDA[iL]);
    fHistList->Add(hDAv[iL]);
    hDXY[iL] = new TH2F(Form("hDXY%d",iL),Form("Distance between 2 hits on layer %d ; X (microns) ; Y (microns)",iL),100,0.,10.,100,0.,+10.);
    //hDXY[iL] = new TH2F(Form("hDXY%d",iL),Form("Distance between 2 hits on layer %d ; X (cm) ; Y (cm)",iL),100,0.,+1.0,100,0.,+1.0);
    fHistList->Add(hDXY[iL]);
    
  }
  hNGoodCell = new TH1F("hNGoodCell","Number of good cell in the track",5,-0.5,4.5);
  fHistList->Add(hNGoodCell);
  hTrackType = new TH1F("hTrackType","Track type",4,-0.5,3.5);
  fHistList->Add(hTrackType);
  hAngleCells = new TH1F("hAngleCell1gap2gaps","#Delta #theta between two segments one 1 gap and one 2 layers gap; #Delta #theta (deg)",200,0.,45.);
  fHistList->Add(hAngleCells);
  hThetaCells = new TH1F("hThetaCells","#theta of the cells; #theta (deg)",200,0.,15.);
  fHistList->Add(hThetaCells);
  
  // limited range for possible z vertex
  hVertZ = new TH1F("hVertZ","hVertZ",400,-10.,+10.);
  hVertZa = new TH1F("hVertZa","hVertZa",400,-10.,+10.);
  
  fHistList->Add(hVertZ);
  fHistList->Add(hVertZa);
  
  // temporary
  hHitDifXY = new TH2F("hHitDifXY","hHitDifXY",100,-0.04,+0.04,100,-0.04,+0.04);
  fHistList->Add(hHitDifXY);
  hAngDifAll = new TH1F("hAngDifAll","hAngDifAll",1000,0.,+0.5);
  fHistList->Add(hAngDifAll);
  hAngDifDup = new TH1F("hAngDifDup","hAngDifDup",1000,0.,+0.5);
  fHistList->Add(hAngDifDup);
  hIntDifXYAll = new TH2F("hIntDifXYAll","hIntDifXYAll",100,-0.01,+0.01,100,-0.01,+0.01);
  fHistList->Add(hIntDifXYAll);
#ifdef OLDGEOM
  hIntDifXYDup = new TH2F("hIntDifXYDup","hIntDifXYDup",100,-0.01,+0.01,100,-0.01,+0.01);
#else
  hIntDifXYDup = new TH2F("hIntDifXYDup","hIntDifXYDup",100,-0.1,+0.1,100,-0.1,+0.1);
#endif
  fHistList->Add(hIntDifXYDup);
  
}

//___________________________________________________________________________
AliMFTTrackFinder::~AliMFTTrackFinder()
{
  
}
//___________________________________________________________________________
void AliMFTTrackFinder::Init(Char_t *parfile)
{

  AliInfo(Form("Parameter file  set to  %s ",parfile));
  
  ReadParam(parfile);
  
  ReadGeom();
  
  SetDebug(1);
  

  
}
//___________________________________________________________________________
void AliMFTTrackFinder::LoadClusters( TClonesArray *clusterArrayFront[AliMFTConstants::fNMaxPlanes], TClonesArray *clusterArrayBack[AliMFTConstants::fNMaxPlanes] ){
  AliCodeTimerAuto("",0);
 AliMFTCALayer *caLayer = NULL;
  AliMFTCAHit   *caHit = NULL;

  for (int iPlane = 0; iPlane<AliMFTConstants::kNDisks; iPlane++) {
    Int_t nClusterFront = clusterArrayFront[iPlane]->GetEntriesFast();
    Int_t nClusterBack  = clusterArrayBack[iPlane]->GetEntriesFast();
    AliDebug(1,Form("Loading %d + %d = %d Clusters for Plane %d",nClusterFront , nClusterBack,nClusterFront + nClusterBack, iPlane));

    // Treating FRONT face of the plane
    caLayer = GetLayer(iPlane*2);
    for (Int_t iCluster=0; iCluster<nClusterFront; iCluster++) {
      caHit = caLayer->AddHit();
      AliMFTCluster * cluster = (AliMFTCluster *)clusterArrayFront[iPlane]->At(iCluster);
      caHit->SetPos(cluster->GetX(), cluster->GetY(), cluster->GetZ());
      caHit->SetTrackGID(cluster->GetMCLabel(0),iPlane*2,caLayer->GetNhits()-1,0);
      caHit->SetMFTClsId(iCluster);
    }
    // Treating BACK face of the plane
    caLayer = GetLayer(iPlane*2+1);
    for (Int_t iCluster=0; iCluster<nClusterBack; iCluster++) {
      caHit = caLayer->AddHit();
      AliMFTCluster * cluster = (AliMFTCluster *)clusterArrayBack[iPlane]->At(iCluster);
      caHit->SetPos(cluster->GetX(), cluster->GetY(), cluster->GetZ());
      caHit->SetTrackGID(cluster->GetMCLabel(0),iPlane*2+1,caLayer->GetNhits()-1,0);
      caHit->SetMFTClsId(iCluster);
    }

  }
}

//___________________________________________________________________________
void AliMFTTrackFinder::ReadParam(Char_t *parfile)
{
  
  AliInfo(Form("Reading Parameter File %s ",parfile));
  
  std::ifstream in;
  in.open(parfile,std::ios::in);
  
  std::string line;
  
  getline(in,line);
  in >> fUseTF;
  printf("Use the TrackFinder: %d \n",fUseTF);
  
  getline(in,line); getline(in,line);
  in >> fNlayers;
  printf("Number of detecting planes: %d \n",fNlayers);
  
  getline(in,line); getline(in,line);
  in >> fThick;
  printf("Layer thickness in X0: %5.3f \n",fThick);
  
  getline(in,line); getline(in,line);
  in >> fPixelNoise;
  printf("Pixel noise: %4.2e \n",fPixelNoise);
  in >> fAddNoise;
  printf("Add noise: %d \n",fAddNoise);
  
  getline(in,line); getline(in,line);
  in >> fMinTrackLength;
  printf("fMinTrackLength: %d \n",fMinTrackLength);
  
  getline(in,line); getline(in,line);
  in >> fXCut;
  printf("fXCut: %6.4f cm\n",fXCut);
  in >> fYCut;
  printf("fYCut: %6.4f cm\n",fYCut);
  
  getline(in,line); getline(in,line);
  in >> fMaxSegAngle;
  printf("fMaxSegAngle: %4.1f deg\n",fMaxSegAngle);
  fMaxSegAngle = 180. - fMaxSegAngle;
  
  getline(in,line); getline(in,line);
  for (Int_t i = 0; i < fNlayers; i++) {
    in >> fACutV[i];
    printf("fACutV[%d]: %4.2f \n",i,fACutV[i]);
  }
  
  getline(in,line); getline(in,line);
  for (Int_t i = 0; i < fNlayers; i++) {
    in >> fACutN[i];
    printf("fACutN[%d]: %4.2f \n",i,fACutN[i]);
  }
  
  getline(in,line); getline(in,line);
  for (Int_t i = 0; i < fNlayers; i++) {
    in >> fPlaneDetEff[i];
    printf("PlaneDetEff[%d]: %4.2f \n",i,fPlaneDetEff[i]);
  }
  
  getline(in,line); getline(in,line);
  in >> fCalcVertex;
  printf("Calculate vertex: %d \n",fCalcVertex);
  
  getline(in,line); getline(in,line);
  in >> fAddQED;
  printf("Add QED hits: %d \n",fAddQED);
  in >> fMBRate;
  printf("Hadronic MB Rate: %5.1f kHz\n",fMBRate);
  in >> fCMOSIntTime;
  printf("CMOS integration time: %5.1f microsec \n",fCMOSIntTime);
  
  getline(in,line); getline(in,line);
  in >> fReadGeom;
  printf("Read geometry: %d \n",fReadGeom);
  in >> fGeomName;
  if (fReadGeom) printf("... from file: %s \n",fGeomName.Data());
  
  in.close();
  printf("==================================\n");
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::Clear(Option_t *)
{
  
  if (fTracks) fTracks->Clear("C");
  
  for (Int_t iL = 0; iL < fNlayers; iL++) {
    fLayers[iL]->Clear("");
  }
  
  fCellGID       = 0;
  fMaxCellStatus = 0;
  fNtracks       = 0;
  
  if (fRoads) fRoads->Clear("C");
  
  fNRoads = 0;
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::ClearCells()
{
  
  for (Int_t iL = 0; iL < fNlayers; iL++) {
    fLayers[iL]->ClearCells();
  }
  
  fCellGID       = 0;
  fMaxCellStatus = 0;
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::ResetCells()
{
  
  AliMFTCALayer *layer;
  AliMFTCACell *cell;
  
  for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
    layer = GetLayer(iL);
    for (Int_t iC = 0; iC < layer->GetNcells(); iC++) {
      cell = layer->GetCell(iC);
      cell->Reset();
    }
  }
  
}

//___________________________________________________________________________
Double_t AliMFTTrackFinder::GetCellAngleDif(AliMFTCACell *cell1, AliMFTCACell *cell2) {
  
  Double_t *c1h1 = cell1->GetHit1();
  Double_t *c1h2 = cell1->GetHit2();
  Double_t *c2h1 = cell2->GetHit1();
  Double_t *c2h2 = cell2->GetHit2();
  
  TVector3 seg1 = TVector3(c1h2[0]-c1h1[0],c1h2[1]-c1h1[1],c1h2[2]-c1h1[2]);
  TVector3 seg2 = TVector3(c2h2[0]-c2h1[0],c2h2[1]-c2h1[1],c2h2[2]-c2h1[2]);
  
  return (seg1.Angle(seg2))*TMath::RadToDeg();
  
}

//___________________________________________________________________________
Double_t AliMFTTrackFinder::GetCellInterceptX(AliMFTCACell *cell, Double_t z) {
  
  Double_t *ch1 = cell->GetHit1();
  Double_t *ch2 = cell->GetHit2();
  
  return ch1[0]+(ch2[0]-ch1[0])/(ch2[2]-ch1[2])*(z-ch1[2]);
  
}

//___________________________________________________________________________
Double_t AliMFTTrackFinder::GetCellInterceptY(AliMFTCACell *cell, Double_t z) {
  
  Double_t *ch1 = cell->GetHit1();
  Double_t *ch2 = cell->GetHit2();
  
  return ch1[1]+(ch2[1]-ch1[1])/(ch2[2]-ch1[2])*(z-ch1[2]);
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::AnalyzeCells() {
  
  printf("Total number of cells = %d \n",GetNcells());
  
  // analyze cell for aliroot input
  
  Double_t cellAngDifCut; // [deg]
  Double_t cellIntDifCut; // [cm]
  
  cellAngDifCut = 0.03;
  cellIntDifCut = 0.01;
  //cellIntDifCut = fXCut; // old
  
  AliMFTCALayer *layer;
  AliMFTCACell *cell1, *cell2, *cellM, *cell;
  Int_t trkid1h1, trkid1h2, trkid2h1, trkid2h2;
  Double_t xc1h1, xc1h2, xc2h1, xc2h2;
  Double_t yc1h1, yc1h2, yc2h1, yc2h2;
  Double_t zc1h1, zc1h2, zc2h1, zc2h2, zMin, zMax;
  Double_t cellAngDif, cellIntDifX, cellIntDifY, zint;
  Double_t cellDifX1, cellDifX2, cellDifY1, cellDifY2;
  Double_t hit1[3], hit2[3];
  Bool_t multCell1;
  TList *multCell = new TList();
  
  const Int_t nMaxh = 100;
  Double_t xTr[nMaxh], yTr[nMaxh], zTr[nMaxh];
  Double_t ax, axe, bx, bxe, ay, aye, by, bye;
  Double_t xTrErrDet = 0.0025/TMath::Sqrt(12.);
  Double_t yTrErrDet = 0.0025/TMath::Sqrt(12.);
  Double_t xTrErrMS = 0.00055; // estimated at p = 5.5 GeV/c
  Double_t yTrErrMS = 0.00055; // estimated at p = 5.5 GeV/c
  Double_t xTrErr[nMaxh], yTrErr[nMaxh];
  for (Int_t i = 0; i < nMaxh; i++) {
    xTrErr[i] = TMath::Sqrt(xTrErrDet*xTrErrDet+xTrErrMS*xTrErrMS);
    yTrErr[i] = TMath::Sqrt(yTrErrDet*yTrErrDet+yTrErrMS*yTrErrMS);
  }
  Int_t nTr, nCells;
  
  // print info on multiple cells
  //cell = GetCellByGID(3091);
  //cell->PrintCell("MC");
  if (kFALSE) {
    Bool_t recTrack;
    Int_t nCells, nCells1, nDiffTracks = 0;
    Int_t nTrackID[10000], TrackID[10000];
    for (Int_t i = 0; i < 10000; i++) {
      TrackID[i]  = -1;
      nTrackID[i] =  0;
    }
    for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
      layer = GetLayer(iL);
      nCells = layer->GetNcells();
      
      nDiffTracks = 0;
      for (Int_t i = 0; i < 10000; i++) {
        TrackID[i] = -1;
        nTrackID[i] = 0;
      }
      nCells1 = 0;
      
      for (Int_t iC = 0; iC < nCells; iC++) {
        cell = layer->GetCell(iC);
        if (cell->GetTrackGID(0) == cell->GetTrackGID(1)) {
          nCells1++;
          recTrack = kTRUE;
          for (Int_t idt = 0; idt < nDiffTracks; idt++) {
            if (TrackID[idt] == cell->GetTrackGID(0)) {
              nTrackID[idt]++;
              recTrack = kFALSE;
              break;
            }
          }
          if (recTrack) {
            TrackID[nDiffTracks] = cell->GetTrackGID(0);
            nTrackID[nDiffTracks]++;
            nDiffTracks++;
          }
        }
      }
      
      if (nDiffTracks < nCells1) {
        printf("AnalyzeCells: L %d cells %d , %d diff tracks %d \n",iL,nCells,nCells1,nDiffTracks);
      }
      for (Int_t idt = 0; idt < nDiffTracks; idt++) {
        if (nTrackID[idt] > 1) {
          for (Int_t iC = 0; iC < nCells; iC++) {
            cell = layer->GetCell(iC);
            if (cell->GetTrackGID(0) == TrackID[idt] &&
                cell->GetTrackGID(1) == TrackID[idt]) {
              printf("Cell:  \n"); cell->PrintCell("MC");
            }
          }
        }
      }
      
    } // end loop layers
      //return;
  } // end print info on multiple cells
  
  for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
    
    layer = GetLayer(iL);
    
    nCells = layer->GetNcells();
    
    for (Int_t iC1 = 0; iC1 < nCells; iC1++) {
      
      cell1 = layer->GetCell(iC1);
      if (cell1->GetStatus() < 1 || cell1->IsMerged()) continue;
      
      zc1h1 = cell1->GetHitp1()[2];
      zc1h2 = cell1->GetHitp2()[2];
      trkid1h1 = cell1->GetTrackGID(0);
      trkid1h2 = cell1->GetTrackGID(1);
      
      multCell1 = kFALSE;
      for (Int_t iC2 = (iC1+1); iC2 < nCells; iC2++) {
        
        cell2 = layer->GetCell(iC2);
        if (cell2->GetStatus() < 1 || cell2->IsMerged()) continue;
        
        if (cell2->GetLength() != cell1->GetLength()) continue;
        
        zc2h1 = cell2->GetHitp1()[2];
        zc2h2 = cell2->GetHitp2()[2];
        trkid2h1 = cell2->GetTrackGID(0);
        trkid2h2 = cell2->GetTrackGID(1);
        
        if ((TMath::Abs(zc1h1-zc2h1) > 0.5*fZGap) ||
            (TMath::Abs(zc1h2-zc2h2) > 0.5*fZGap)) {
          
          xc1h1 = cell1->GetHit1()[0];
          xc1h2 = cell1->GetHit2()[0];
          xc2h1 = cell2->GetHit1()[0];
          xc2h2 = cell2->GetHit2()[0];
          yc1h1 = cell1->GetHit1()[1];
          yc1h2 = cell1->GetHit2()[1];
          yc2h1 = cell2->GetHit1()[1];
          yc2h2 = cell2->GetHit2()[1];
          /*
           // angle between the two cells
           cellAngDif = GetCellAngleDif(cell1,cell2);
           // intercept point at the mid distance between layers
           zint = 0.5*(fPlanesZ[cell1->GetLayers()[0]]+fPlanesZ[cell1->GetLayers()[1]]);
           cellIntDifX = GetCellInterceptX(cell1,zint)-GetCellInterceptX(cell2,zint);
           cellIntDifY = GetCellInterceptY(cell1,zint)-GetCellInterceptY(cell2,zint);
           */
          cellDifX1 = xc1h1-xc2h1;
          cellDifX2 = xc1h2-xc2h2;
          cellDifY1 = yc1h1-yc2h1;
          cellDifY2 = yc1h2-yc2h2;
          
          if (trkid1h1 == trkid2h1 && trkid1h2 == trkid2h2) {
            cellAngDif = GetCellAngleDif(cell1,cell2);
            hAngDifDup->Fill(cellAngDif);
            hIntDifXYDup->Fill(cellDifX1,cellDifY1);
            hIntDifXYDup->Fill(cellDifX2,cellDifY2);
            //printf("CellIntDif %f %f %f %f \n",cellDifX1,cellDifX2,cellDifY1,cellDifY2);
          }
          
          //if ((cellAngDif < cellAngDifCut) &&
          //    (TMath::Abs(cellIntDifX) < fXCut) &&
          //    (TMath::Abs(cellIntDifY) < fYCut)) {
          if (TMath::Abs(cellDifX1) < cellIntDifCut &&
              TMath::Abs(cellDifX2) < cellIntDifCut &&
              TMath::Abs(cellDifY1) < cellIntDifCut &&
              TMath::Abs(cellDifY2) < cellIntDifCut) {
            if (!multCell1) {
              //printf("Add mult cell1 %d \n",cell1->GetGID());
              //cell1->PrintCell("MC");
              multCell1 = kTRUE;
              multCell->Add(cell1);
            }
            //printf("Add mult cell2 %d \n",cell2->GetGID());
            //cell2->PrintCell("MC");
            multCell->Add(cell2);
          }
        }
        /*
         if (trkid1h1 != trkid2h1 && trkid1h2 != trkid2h2) {
         // angle between the two cells
         cellAngDif = GetCellAngleDif(cell1,cell2);
         // intercept point at the mid distance between layers
         zint = 0.5*(fPlanesZ[cell1->GetLayers()[0]]+fPlanesZ[cell1->GetLayers()[1]]);
         cellIntDifX = GetCellInterceptX(cell1,zint)-GetCellInterceptX(cell2,zint);
         cellIntDifY = GetCellInterceptY(cell1,zint)-GetCellInterceptY(cell2,zint);
         cellIntDif = TMath::Sqrt(cellIntDifX*cellIntDifX+cellIntDifY*cellIntDifY);
         hAngDifAll->Fill(cellAngDif);
         hIntDifXYAll->Fill(cellIntDifX,cellIntDifY);
         }
         */
      } // end cell2 loop
      
      // merge multi cells
      if (multCell->GetSize() > 0) {
        
        //printf("multCell size %d \n",multCell->GetSize());
        nTr = 0;
        zMin = +9999.;
        zMax = -9999.;
        //printf("Found %d multi cells.\n",multCell->GetSize());
        for (Int_t imc = 0; imc < multCell->GetSize(); imc++) {
          cellM = (AliMFTCACell*)multCell->At(imc);
          xTr[nTr] = cellM->GetHit1()[0];
          yTr[nTr] = cellM->GetHit1()[1];
          zTr[nTr] = cellM->GetHit1()[2];
          zMin = TMath::Min(zMin,zTr[nTr]);
          zMax = TMath::Max(zMax,zTr[nTr]);
          nTr++;
          xTr[nTr] = cellM->GetHit2()[0];
          yTr[nTr] = cellM->GetHit2()[1];
          zTr[nTr] = cellM->GetHit2()[2];
          zMin = TMath::Min(zMin,zTr[nTr]);
          zMax = TMath::Max(zMax,zTr[nTr]);
          nTr++;
          cellM->SetStatus(0);
          //printf("Rm cell %d \n",cellM->GetGID());
          //cellM->PrintCell("MC");
        }
        /*
         printf("Number of points: %d \n",nTr);
         for (Int_t iTr = 0; iTr < nTr; iTr++) {
         printf("%d %f %f %f \n",iTr,xTr[iTr],yTr[iTr],zTr[iTr]);
         }
         */
        if (LinFit(nTr,zTr,xTr,xTrErr,ax,axe,bx,bxe) &&
            LinFit(nTr,zTr,yTr,yTrErr,ay,aye,by,bye)) {
          hit1[0] = ax*zMin+bx;
          hit1[1] = ay*zMin+by;
          hit1[2] = zMin;
          hit2[0] = ax*zMax+bx;
          hit2[1] = ay*zMax+by;
          hit2[2] = zMax;
          // add a new cell
          cell = layer->AddCell();
          cell->SetHits(hit1,hit2,fPlanesZ[cell1->GetLayers()[0]],fPlanesZ[cell1->GetLayers()[1]]);
          cell->SetLayers(cell1->GetLayers()[0],cell1->GetLayers()[1]);
          cell->SetStatus(1);
          cell->SetGID(fCellGID++,trkid1h1,trkid1h2);
          cell->SetIsMerged();
          //printf("Add merged cell %d \n",cell->GetGID());
          //printf("H1: %f %f %f \n",hit1[0],hit1[1],hit1[2]);
          //printf("H2: %f %f %f \n",hit2[0],hit2[1],hit2[2]);
          /*
           for (Int_t imc = 0; imc < multCell->GetSize(); imc++) {
           cellM = (AliMFTCACell*)multCell->At(imc);
           cellAngDif = GetCellAngleDif(cell,cellM);
           hAngDifDup->Fill(cellAngDif);
           }
           */
        } else {
          printf("No line fit possible!\n");
        }
        multCell->Clear();
      } // end merge multi cells
      
    } // end cell1 loop
    
  } // end layer loop
  
  delete multCell;
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::CreateGapCells() {
  
  // create long cells (with one missing layer in the middle) starting from the
  // short cells which did not find neighbours during the first RunForward()
  
  Bool_t prn = kFALSE;
  
  AliMFTCALayer *layerL;
  AliMFTCACell *cellL;
  Double_t h1[3], h2[3];
  Int_t iL2, nH2, trackID1, trackID2;
  
  Int_t nGapCells = 0;
  
  for (Int_t iL1 = 0; iL1 < (fNlayers-1); iL1++) {
    
    layerL = GetLayer(iL1);
    
    for (Int_t iCL = 0; iCL < layerL->GetNcells(); iCL++) { // Loop over cell in layer iL
      
      // short cell
      cellL = layerL->GetCell(iCL);
      
      // attach at right a long cell
      if ((iL1 < (fNlayers-3)) && (cellL->GetNNbR() == 0)) {
        
        h1[0] = cellL->GetHit2()[0];
        h1[1] = cellL->GetHit2()[1];
        h1[2] = cellL->GetHit2()[2];
        
        iL2 = iL1 + 3;
        nH2 = GetLayer(iL2)->GetNhits();
        
        for (Int_t iH2 = 0; iH2 < nH2; iH2++) {
          
          h2[0] = GetLayer(iL2)->GetHit(iH2)->GetPos()[0];
          h2[1] = GetLayer(iL2)->GetHit(iH2)->GetPos()[1];
          h2[2] = GetLayer(iL2)->GetHit(iH2)->GetPos()[2];
          
          TVector3 vec(h2[0]-h1[0],h2[1]-h1[1],h2[2]-h1[2]);
          if (vec.Theta() < fMaxSegAngle*TMath::DegToRad()) continue;
          
          if (!RuleSelect2LayersGap(iL1+1,iL1+2,h1,h2)) continue;
          
          if (!RuleSelectCell(h1,h2,iL1+1)) continue;
          
          nGapCells++;
          
          AliMFTCACell *cell = GetLayer(iL1+1)->AddCell();
          cell->SetHits(h1,h2,fPlanesZ[iL1+1],fPlanesZ[iL1+3]);
          cell->SetLayers(iL1+1,iL1+3);
          cell->SetStatus(1);
          trackID1 = cellL->GetTrackGID(1);
          trackID2 = GetLayer(iL1+3)->GetHit(iH2)->GetTrackGID();
          cell->SetGID(fCellGID++,trackID1,trackID2);
          
          if (prn) printf("Create gap (L) cell %d in layer %d. \n",cell->GetGID(),iL1+1);
          
        } // end loop hits
        
      } // end attach at right a long cell
      
      // attach at left a long cell; 1 <> 2
      if ((iL1 > 1) && (cellL->GetNNbL() == 0)) {
        
        h1[0] = cellL->GetHit1()[0];
        h1[1] = cellL->GetHit1()[1];
        h1[2] = cellL->GetHit1()[2];
        
        iL2 = iL1 - 2;
        nH2 = GetLayer(iL2)->GetNhits();
        
        for (Int_t iH2 = 0; iH2 < nH2; iH2++) {
          
          h2[0] = GetLayer(iL2)->GetHit(iH2)->GetPos()[0];
          h2[1] = GetLayer(iL2)->GetHit(iH2)->GetPos()[1];
          h2[2] = GetLayer(iL2)->GetHit(iH2)->GetPos()[2];
          
          TVector3 vec(h1[0]-h2[0],h1[1]-h2[1],h1[2]-h2[2]);
          if (vec.Theta() < fMaxSegAngle*TMath::DegToRad()) continue;
          
          if (!RuleSelect2LayersGap(iL1-2,iL1-1,h2,h1)) continue;
          
          if (!RuleSelectCell(h2,h1,iL1-2)) continue;
          
          nGapCells++;
          
          AliMFTCACell *cell = GetLayer(iL1-2)->AddCell();
          cell->SetHits(h2,h1,fPlanesZ[iL1-2],fPlanesZ[iL1]);
          cell->SetLayers(iL1-2,iL1);
          cell->SetStatus(1);
          trackID2 = cellL->GetTrackGID(0);
          trackID1 = GetLayer(iL1-2)->GetHit(iH2)->GetTrackGID();
          cell->SetGID(fCellGID++,trackID1,trackID2);
          
          if (prn) printf("Create gap (R) cell %d in layer %d. \n",cell->GetGID(),iL1-2);
          
        } // end loop hits
        
      } // end attach at left a long cell
      
    } // end loop cells
    
  } // end loop layers
  
  printf("Found %d gap cells.\n",nGapCells);
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::CreateCellsR(AliMFTCARoad *road) {
  
  // create cells from hits selected in roads
  
  Bool_t prn = kFALSE;
  
  AliMFTCACell *cell;
  AliMFTCAHit *hit1, *hit2;
  Int_t iL1, iL2, nH1, nH2;
  Int_t iL1min, iL1max;
  Int_t iL2min, iL2max;
  Int_t trackID1, trackID2, detElemID1, detElemID2;
  Bool_t noCell;
  Double_t h1[3], h2[3], h[3], hx, hy, dR;
  Int_t mftClsId1, mftClsId2;
  Int_t nCombi = 0;
  
  iL1min = road->GetLayer1();
  iL1max = road->GetLayer2()-1;
  //printf("R%d iL1 %d %d \n",ir,iL1min,iL1max);
  
  for (iL1 = iL1min; iL1 <= iL1max; iL1++) {
    //printf("iL1 %d \n",iL1);
    iL2min = iL1+1;
    nH1 = road->GetNhitsInLayer(iL1);
    for (Int_t iH1 = 0; iH1 < nH1; iH1++) {
      //printf("iH1 %d \n",iH1);
      hit1 = road->GetHitInLayer(iL1,iH1);
      iL2max = TMath::Min((iL1+(5-hit1->IsFace())),fNlayers-1);
      h1[0] = hit1->GetPos()[0];
      h1[1] = hit1->GetPos()[1];
      h1[2] = hit1->GetPos()[2];
      mftClsId1 = hit1->GetMFTClsId();
      noCell = kTRUE;
      iL2 = iL2min;
      while (noCell && (iL2 <= iL2max)) {
        //printf("iL2 %d \n",iL2);
        nH2 = road->GetNhitsInLayer(iL2);
        for (Int_t iH2 = 0; iH2 < nH2; iH2++) {
          hit2 = road->GetHitInLayer(iL2,iH2);
          h2[0] = hit2->GetPos()[0];
          h2[1] = hit2->GetPos()[1];
          h2[2] = hit2->GetPos()[2];
	  mftClsId2 = hit2->GetMFTClsId();
          nCombi++;
          if (RuleSelectCell(h1,h2,iL1)) {
            noCell = kFALSE;
            cell = GetLayer(iL1)->AddCell();
            cell->SetHits(h1,h2,fPlanesZ[iL1],fPlanesZ[iL2]);
	    cell->SetMFTClsId(mftClsId1,mftClsId2);
            cell->SetLayers(iL1,iL2);
            cell->SetStatus(1);
            trackID1 = hit1->GetTrackGID();
            trackID2 = hit2->GetTrackGID();
            detElemID1 = hit1->GetDetElemID();
            detElemID2 = hit2->GetDetElemID();
            cell->SetGID(fCellGID++,trackID1,trackID2);
            cell->SetDetElemID(detElemID1,detElemID2);
            road->AddCell(cell);
            //printf("Cell %5d: L %d-%d H %03d-%03d MC %04d-%04d \n",cell->GetGID(),iL1,iL2,iH1,iH2,trackID1,trackID2);
          } // end create cell
        } // end loop iH2
        iL2++;
      } // end loop iL2
    } // end loop iH1
  } // end loop iL1
  /*
   for (Int_t iL = 0; iL < fNlayers; iL++) {
   printf("%5d %2d %5d \n",ir,iL,road->GetNhitsInLayer(iL));
   }
   */
  
  if (kFALSE) {
    for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
      for (Int_t iC = 0; iC < road->GetNcellsInLayer(iL); iC++) {
        cell = road->GetCellInLayer(iL,iC);
        printf("L%d,%d-%d CellGID %d MC %d %d \n",iL,cell->GetLayers()[0],cell->GetLayers()[1],cell->GetGID(),cell->GetTrackGID(0),cell->GetTrackGID(1));
      }
    }
  }
  
  if (kFALSE) {
    printf("From %d combinations: \n",nCombi);
    Long_t nTotCell =0;
    for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
      printf("Layer %d nr of cells: %d \n",iL,GetLayer(iL)->GetNcells());
      nTotCell += GetLayer(iL)->GetNcells();
    }
    printf("Tot cells: %ld \n",nTotCell);
    
    for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
      Int_t nc = road->GetNcellsInLayer(iL);
      printf("L%1d C%2d \n",iL,nc);
    }
    
  }
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::CreateCells(Bool_t cv /*Calculate Vertex*/) {
  
  Bool_t prn = kFALSE;
  
  AliMFTCACell *cell;
  AliMFTCAHit *hit1, *hit2;
  Int_t iL1, iL2, nH1, nH2;
  Int_t iL1min, iL1max;
  Int_t iL2min, iL2max;
  Int_t trackID1, trackID2, detElemID1, detElemID2;
  Bool_t noCell;
  Double_t h1[3], h2[3], h[3], hx, hy, dR;
  Int_t nCombi = 0;
  
  iL1min = 0;
  iL1max = fNlayers-2;
  
  for (iL1 = iL1min; iL1 <= iL1max; iL1++) {
    //printf("iL1 %d \n",iL1);
    iL2min = iL1+1;
    nH1 = GetLayer(iL1)->GetNhits();
    for (Int_t iH1 = 0; iH1 < nH1; iH1++) {
      //printf("iH1 %d \n",iH1);
      hit1 = GetLayer(iL1)->GetHit(iH1);
      iL2max = TMath::Min((iL1+(5-hit1->IsFace())),fNlayers-1);
      h1[0] = hit1->GetPos()[0];
      h1[1] = hit1->GetPos()[1];
      h1[2] = hit1->GetPos()[2];
      noCell = kTRUE;
      iL2 = iL2min;
      while (noCell && (iL2 <= iL2max)) {
        //printf("iL2 %d \n",iL2);
        nH2 = GetLayer(iL2)->GetNhits();
        for (Int_t iH2 = 0; iH2 < nH2; iH2++) {
          hit2 = GetLayer(iL2)->GetHit(iH2);
          h2[0] = hit2->GetPos()[0];
          h2[1] = hit2->GetPos()[1];
          h2[2] = hit2->GetPos()[2];
          nCombi++;
          if (RuleSelectCell(h1,h2,iL1)) {
            noCell = kFALSE;
            cell = GetLayer(iL1)->AddCell();
            cell->SetHits(h1,h2,fPlanesZ[iL1],fPlanesZ[iL2]);
            cell->SetLayers(iL1,iL2);
            cell->SetStatus(1);
            trackID1 = hit1->GetTrackGID();
            trackID2 = hit2->GetTrackGID();
            detElemID1 = hit1->GetDetElemID();
            detElemID2 = hit2->GetDetElemID();
            cell->SetGID(fCellGID++,trackID1,trackID2);
            cell->SetDetElemID(detElemID1,detElemID2);
            //printf("Cell %5d: L %d-%d H %03d-%03d MC %04d-%04d \n",cell->GetGID(),iL1,iL2,iH1,iH2,trackID1,trackID2);
          } // end create cell
        } // end loop iH2
        iL2++;
      } // end loop iL2
    } // end loop iH1
  } // end loop iL1
  /*
   for (Int_t iL = 0; iL < fNlayers; iL++) {
   printf("%5d %2d %5d \n",ir,iL,road->GetNhitsInLayer(iL));
   }
   */
  
  if (kFALSE) {
    printf("From %d combinations: \n",nCombi);
    Long_t nTotCell =0;
    for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
      printf("Layer %d nr of cells: %d \n",iL,GetLayer(iL)->GetNcells());
      nTotCell += GetLayer(iL)->GetNcells();
    }
    printf("Tot cells: %ld \n",nTotCell);
    
  }
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::CreateCellsOld(Bool_t cv /*Calculate Vertex*/) {
  
  // create only short cells (between two subsequent layers)
  
  // print info on multiple cells
  if (kFALSE) {
    Bool_t recTrack;
    Int_t nHits, nDiffTracks = 0;
    Int_t nTrackID[10000], TrackID[10000];
    for (Int_t i = 0; i < 10000; i++) {
      TrackID[i] = -1;
      nTrackID[i] = 0;
    }
    for (Int_t iL = 0; iL < fNlayers; iL++) {
      nHits = GetLayer(iL)->GetNhits();
      nDiffTracks = 0;
      for (Int_t i = 0; i < 10000; i++) {
        TrackID[i] = -1;
        nTrackID[i] = 0;
      }
      for (Int_t iH = 0; iH < nHits; iH++) {
        recTrack = kTRUE;
        for (Int_t idt = 0; idt < nDiffTracks; idt++) {
          if (TrackID[idt] == GetLayer(iL)->GetHit(iH)->GetTrackGID()) {
            nTrackID[idt]++;
            recTrack = kFALSE;
            break;
          }
        }
        if (recTrack) {
          TrackID[nDiffTracks] = GetLayer(iL)->GetHit(iH)->GetTrackGID();
          nTrackID[nDiffTracks]++;
          nDiffTracks++;
        }
      }
      printf("CreateCells: L %d hits %d diff tracks %d \n",iL,nHits,nDiffTracks);
    }
  }
  
  Bool_t prn = kFALSE;
  
  // loop over the layers
  Int_t nH1, nH2, trackID1, trackID2, detElemID1, detElemID2;
  Double_t h1[3], h2[3];
  Double_t nCombi = 0.;
  
  Int_t iL1, iL2;
  
  iL1 = 0;
  while (iL1 < (fNlayers-1)) {
    
    iL2 = iL1 + 1;
    
    nH1 = GetLayer(iL1)->GetNhits();
    if (prn) printf("---> 1st Layer L %d with %d hits.\n",iL1,nH1);
    
    nH2 = GetLayer(iL2)->GetNhits();
    if (prn) printf("---> 2nd Layer R %d with %d hits.\n",iL2,nH2);
    
    for (Int_t iH1 = 0; iH1 < nH1; iH1++) {
      
      h1[0] = GetLayer(iL1)->GetHit(iH1)->GetPos()[0];
      h1[1] = GetLayer(iL1)->GetHit(iH1)->GetPos()[1];
      h1[2] = GetLayer(iL1)->GetHit(iH1)->GetPos()[2];
      
      for (Int_t iH2 = 0; iH2 < nH2; iH2++) {
        
        h2[0] = GetLayer(iL2)->GetHit(iH2)->GetPos()[0];
        h2[1] = GetLayer(iL2)->GetHit(iH2)->GetPos()[1];
        h2[2] = GetLayer(iL2)->GetHit(iH2)->GetPos()[2];
        
        TVector3 vec(h2[0]-h1[0],h2[1]-h1[1],h2[2]-h1[2]);
        if (vec.Theta() < fMaxSegAngle*TMath::DegToRad()) continue;
        
        nCombi += 1.;
        
        if (!cv && !RuleSelectCell(h1,h2,iL1)) continue;
        
        AliMFTCACell *cell = GetLayer(iL1)->AddCell();
        cell->SetHits(h1,h2,fPlanesZ[iL1],fPlanesZ[iL2]);
        cell->SetLayers(iL1,iL2);
        cell->SetStatus(1);
        trackID1 = GetLayer(iL1)->GetHit(iH1)->GetTrackGID();
        trackID2 = GetLayer(iL2)->GetHit(iH2)->GetTrackGID();
        detElemID1 = GetLayer(iL1)->GetHit(iH1)->GetDetElemID();
        detElemID2 = GetLayer(iL2)->GetHit(iH2)->GetDetElemID();
        cell->SetGID(fCellGID++,trackID1,trackID2);
        cell->SetDetElemID(detElemID1,detElemID2);
        //cell->PrintCell("FULL");
        //printf("Cell nr: %d \n",GetLayer(iL1)->GetNcells());
        //cell = GetLayer(iL1)->GetCell(GetLayer(iL1)->GetNcells()-1);
        //cell->PrintCell("FULL");
        
      } // end loop hit in layer 2
      
    } // end loop hit in layer 1
    
    if (prn) printf("Create cell %d in layer %d-%d. \n",GetLayer(iL1)->GetNcells(),iL1,GetLayer(iL1)->GetID());
    
    if (cv) break;
    
    iL1++;
    
  } // end loop layer 1
  
  if (cv) CalculateVertex();
  
  if (kTRUE || prn) {
    printf("From %.0f combinations: \n",nCombi);
    Long_t nTotCell =0;
    for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
      printf("Layer %d nr of cells: %d \n",iL,GetLayer(iL)->GetNcells());
      nTotCell += GetLayer(iL)->GetNcells();
    }
    printf("Tot cells: %ld \n",nTotCell);
  }
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::CalculateVertex() {
  
  hVertZ->Reset();
  hVertZa->Reset();
  
  AliMFTCALayer *layer = GetLayer(0);
  AliMFTCACell *cell;
  const Double_t *h1, *h2;
  Double_t a, b, c, x0, y0, z0;
  Double_t n1, n2, n3, n4;
  Double_t zmin;
  
  for (Int_t iC = 0; iC < layer->GetNcells(); iC++) {
    
    cell = layer->GetCell(iC);
    
    h1 = cell->GetHit1();
    h2 = cell->GetHit2();
    
    x0 = h1[0];
    y0 = h1[1];
    z0 = h1[2];
    a = (h2[0]-h1[0])/(h2[2]-h1[2]);
    b = (h2[1]-h1[1])/(h2[2]-h1[2]);
    c = 1.;
    
    n1 = (c*x0-a*z0)/c;
    n2 = a/c;
    n3 = (c*y0-b*z0)/c;
    n4 = b/c;
    
    zmin = -(n1*n2+n3*n4)/(n2*n2+n4*n4);
    
    hVertZ->Fill(zmin);
    hVertZa->Fill(zmin);
    
  }
  
  Float_t zvert = 0., sum = 0., maxsum = 0.;
  Int_t bin1, bin2, binW = 2;
  Int_t binMin = 1;
  Int_t binMax = hVertZ->GetNbinsX();
  
  hVertZ->Fit("pol1","QW0");
  TF1 *f = hVertZ->GetFunction("pol1");
  
  for (Int_t i = binMin; i <= binMax; i++) {
    hVertZ->SetBinContent(i,TMath::Max(0.,hVertZ->GetBinContent(i)-f->Eval(hVertZ->GetBinCenter(i))));
  }
  
  for (Int_t i = binMin; i <= binMax; i++) {
    bin1 = TMath::Max(binMin,i-binW);
    bin2 = TMath::Min(binMax,i+binW);
    sum = hVertZ->Integral(bin1,bin2);
    if (sum > maxsum) {
      maxsum = sum;
      zvert = hVertZ->GetBinCenter(i);
    }
    //printf("%d %d %d %f %f %f \n",bin1,i,bin2,sum,maxsum,zvert);
  }
  fZVertCalc = zvert;
  printf("Fit vertex z = %f cm \n",fZVertCalc);
  
  // range = zvert +/- 3 cm
  fZVertRange[0] = fZVertCalc - 3.;
  fZVertRange[1] = fZVertCalc + 3.;
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::RunForwardR(AliMFTCARoad *road, Int_t& trackGID) {
  
  Bool_t prn = kFALSE;
  
  if (prn) AliInfo("Run forward (roads) ==================================== \n");
  
  AliMFTCALayer *layerL;
  AliMFTCALayer *layerR;
  AliMFTCACell *cellL;
  AliMFTCACell *cellR;
  Bool_t stch;
  Int_t iL, iR, iter;
  Double_t nCombiTot = 0;
  Double_t nCombiMatch = 0;
  Int_t cellLayers[2];
  
  fMaxCellStatus = 0;
  
  iter = 0;
  
  stch = kTRUE;
  while (stch) {
    
    stch = kFALSE;
    iter++;
    
    for (iL = 0; iL < (fNlayers-2); iL++) { // loop over layers
      
      for (Int_t iCL = 0; iCL < road->GetNcellsInLayer(iL); iCL++) {
        
        cellL = road->GetCellInLayer(iL,iCL);
        if (cellL->GetStatus() == 0) continue;
        
        for(Int_t i = 0; i < 2; i++)  cellLayers[i] = cellL->GetLayers()[i];
        
        iR = iL+(cellLayers[1]-cellLayers[0]);
        if (iR >= (fNlayers-1))
        continue;
        
        for (Int_t iCR = 0; iCR < road->GetNcellsInLayer(iR); iCR++) {
          
          cellR = road->GetCellInLayer(iR,iCR);
          if (cellR->GetStatus() == 0) continue;
          
          nCombiTot += 1;
          
          if ((cellL->GetStatus() == cellR->GetStatus()) &&
              RuleSelect(cellL,cellR)) {
            
            if (prn){
              AliInfo(Form("Matching cells: L(%d) cellGID(%d) %d  R(%d) cellGID(%d) %d  \n",iL,iCL,cellL->GetGID(),iR,iCR,cellR->GetGID()));
            }
            
            nCombiMatch += 1;
            
            if (iter == 1) {
              cellL->AddRightNeighbour(cellR->GetGID());
              cellR->AddLeftNeighbour(cellL->GetGID());
            }
            
            cellR->IncrStatus();
            
            stch = kTRUE;
            
          } // END : matching cells
          
        } // END : loop over cells in layer iR
        
      } // END : loop over cells in layer iL
      
    } // END : loop over layer iL
    
    UpdateCellStatusR();
    
    if (kFALSE || prn) {
      AliInfo(Form("Iteration: %5d ----------------- \n",iter));
      for (iL = 0; iL < (fNlayers-1); iL++) {
        for (Int_t iCL = 0; iCL < road->GetNcellsInLayer(iL); iCL++) {
          cellL = road->GetCellInLayer(iL,iCL);
          if (cellL->HasNbL() || cellL->HasNbR()) {
            printf("L%1d C%03d S%1d GID%03d NNb %d %d \n",iL,iCL,cellL->GetStatus(),cellL->GetGID(),cellL->GetNNbL(),cellL->GetNNbR());
          }
        }
      }
    }
    
  } // end status change
  
  if (kFALSE || prn) {
    printf("End iteration: ----------------- \n");
    for (iL = 0; iL < (fNlayers-1); iL++) {
      for (Int_t iCL = 0; iCL < road->GetNcellsInLayer(iL); iCL++) {
        cellL = road->GetCellInLayer(iL,iCL);
        if (cellL->HasNbL() || cellL->HasNbR()) {
          printf("L%1d C%03d S%1d GID%03d NNb %d %d \n",iL,iCL,cellL->GetStatus(),cellL->GetGID(),cellL->GetNNbL(),cellL->GetNNbR());
        }
      }
    }
  }
  
  if (kFALSE || prn) {
    
    printf("RunForward after %d iterations, nr of Total combinations %.0f , nr of matched combinations  %.0f\n",iter,nCombiTot, nCombiMatch);
    
    printf("After RunForward max cell status = %d \n",fMaxCellStatus);
    
  }
  
  RunBackwardR(road,trackGID);
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::RunForward() {
  
  Bool_t prn = kFALSE;
  
  if (prn) printf("Run forward ==================================== \n");
  
  AliMFTCALayer *layerL;
  AliMFTCALayer *layerR;
  AliMFTCACell *cellL;
  AliMFTCACell *cellR;
  Bool_t stch = kTRUE;
  Int_t iL, iR, iter = 0;
  Double_t nCombiTot = 0;
  Double_t nCombiMatch = 0;
  
  Double_t nCombiIter, nCombiIterMatch;
  
  while (stch) {
    
    nCombiIter = 0;
    nCombiIterMatch = 0;
    
    stch = kFALSE;
    iter++;
    
    for (iL = 0; iL < (fNlayers-2); iL++) { // loop over layers
      
      layerL = GetLayer(iL);
      
      if (prn) printf("L %d cells %d  \n",iL,layerL->GetNcells());
      
      for (Int_t iCL = 0; iCL < layerL->GetNcells(); iCL++) { // Loop over cell in layer iL
        
        cellL = layerL->GetCell(iCL);
        if (cellL->GetStatus() == 0) continue;
        
        Int_t cellLayers[2];
        for(Int_t i = 0; i < 2; i++)  cellLayers[i] = cellL->GetLayers()[i];
        
        iR = iL+(cellLayers[1] - cellLayers[0]);
        if (iR < (fNlayers-1))
        layerR = GetLayer(iR);
        else
        continue;
        
        // vertex selection here ?
        //if (!RuleSelectCell(cellL)) continue;
        
        for (Int_t iCR = 0; iCR < layerR->GetNcells(); iCR++) { // Loop over cell in layer iL
          
          cellR = layerR->GetCell(iCR);
          if (cellR->GetStatus() == 0) continue;
          
          // vertex selection here ?
          //if (!RuleSelectCell(cellR)) continue;
          
          nCombiTot += 1;
          nCombiIter += 1;
          
          if ((cellL->GetStatus() == cellR->GetStatus()) && RuleSelect(cellL,cellR)) { // If Cells are matching in angles and have equal status values
            
            if (prn){
              printf("Matching cells: L cellGID %d  R cellGID %d  \n",cellL->GetGID(),cellR->GetGID());
              printf("Layer L %d cell %d - Layer R %d cell %d \n",iL,iCL,iR,iCR);
            }
            nCombiMatch += 1;
            nCombiIterMatch += 1;
            
            if (iter == 1) {
              cellL->AddRightNeighbour(cellR->GetGID());
              cellR->AddLeftNeighbour(cellL->GetGID());
            }
            
            cellR->IncrStatus();
            
            stch = kTRUE;
            
          } // END : matching cells
          
        } // END : loop over cell in layer iL-1
        
      } // END : loop over cell layer iL
    } // END : loop over layers
    
    if (prn) {
      printf("Iteration: %5d nr of combinations %.0f match %.0f \n",iter,nCombiIter,nCombiIterMatch);
    }
    
    UpdateCellStatus();
    
    if (prn) {
      printf("Iteration: %5d ----------------- \n",iter);
      for (iL = 0; iL < (fNlayers-1); iL++) {
        layerL = GetLayer(iL);
        for (Int_t iCL = 0; iCL < layerL->GetNcells(); iCL++) {
          cellL = layerL->GetCell(iCL);
          if (cellL->HasNbL() || cellL->HasNbR()) {
            printf("L%1d C%03d S%1d GID%03d NNb %d %d \n",iL,iCL,cellL->GetStatus(),cellL->GetGID(),cellL->GetNNbL(),cellL->GetNNbR());
          }
        }
      }
    }
    
  } // end status change
  
  if (prn) {
    printf("End iteration: ----------------- \n");
    for (iL = 0; iL < (fNlayers-1); iL++) {
      layerL = GetLayer(iL);
      for (Int_t iCL = 0; iCL < layerL->GetNcells(); iCL++) {
        cellL = layerL->GetCell(iCL);
        if (cellL->HasNbL() || cellL->HasNbR()) {
          printf("L%1d C%03d S%1d GID%03d NNb %d %d \n",iL,iCL,cellL->GetStatus(),cellL->GetGID(),cellL->GetNNbL(),cellL->GetNNbR());
        }
      }
    }
  }
  
  if (kTRUE || prn) {
    
    printf("RunForward after %d iterations, nr of Total combinations %.0f , nr of matched combinations  %.0f\n",iter,nCombiTot, nCombiMatch);
    
    printf("After RunForward max cell status = %d \n",fMaxCellStatus);
    
  }
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::RunBackwardR(AliMFTCARoad *road, Int_t& trackGID) {
  
  Bool_t prn = kFALSE;
  
  if (prn) printf("Run backward road %d  ==================================== \n",road->GetID());
  
  if (fMaxCellStatus == 1) return; // only isolated cells
  
  Double_t chisqNbLprev, cellAngDif, cellAngDifPrev;
  Int_t addCellIdToTrack, iNbLchiSqMin, iCellAngDifMin, nHitSta;
  
  AliMFTCALayer *layerR;
  AliMFTCACell *cellR;
  AliMFTCACell *cellL;
  AliMFTCACell *cell;
  AliMFTCATrack *track;
  
  Bool_t addCellToTrack, hitSta[5];
  
  Int_t minStartLayer = 6;
  Int_t maxStartLayer = 8;
  
  for (Int_t startLayer = maxStartLayer; startLayer >= minStartLayer; startLayer--) {
    
    if (prn) printf("Start layer %d \n",startLayer);
    
    for (Int_t iCR = 0; iCR < road->GetNcellsInLayer(startLayer); iCR++) {
      
      cellR = road->GetCellInLayer(startLayer,iCR);
      
      if (cellR->GetStatus() == 0) continue;
      if (cellR->IsUsed()) continue;
      if (cellR->GetStatus() < (fMinTrackLength-1)) continue;
      
      if (prn) printf("Create new track %d \n",trackGID);
      
      track = AddTrack(trackGID++);
      track->SetStartLayer(startLayer);
      track->AddCell(cellR);
      cellR->SetUsed(kTRUE);
      
      // add cells to new track
      addCellToTrack = kTRUE;
      while (addCellToTrack) {
        
        cellR = road->GetCellByGID(track->GetLastCellGID()); // !!!
        
        addCellToTrack = kFALSE;
        
        // find the left neighbor giving the smalles chisquare
        iNbLchiSqMin = 0;
        chisqNbLprev = 0.;
        
        // find the left neighbor giving the smallest deviation
        cellAngDifPrev = -1.;
        iCellAngDifMin = 0;
        
        // do a first loop to check all possible associations
        for (Int_t iNNbL = 0; iNNbL < cellR->GetNNbL(); iNNbL++) {
          
          cellL = road->GetCellByGID(cellR->GetNbLgid(iNNbL));
          
          if (kFALSE && prn) {
            printf("To track %d attach cell GID {L}: %d  - TrackID of this cell : %d - %d, Status %d \n",track->GetGID(),cellL->GetGID(),cellL->GetTrackGID(0),cellL->GetTrackGID(1),cellL->GetStatus());
          }
          
          if (cellL->GetStatus() == 0) continue;
          if (cellL->IsUsed()) continue;
          if (cellL->GetStatus() != (cellR->GetStatus()-1)) continue;
          
          // ... smallest deviation
          cellAngDif = GetCellAngleDif(cellL,cellR);
          if (cellAngDifPrev < 0.) {
            cellAngDifPrev = cellAngDif;
          } else {
            if (cellAngDif < cellAngDifPrev) {
              iCellAngDifMin = iNNbL;
            }
          }
          
          // ... smallest chisquare
          if (track->AddCellToChiSq(cellL) < chisqNbLprev) {
            iNbLchiSqMin = iNNbL;
          }
          
          chisqNbLprev = track->AddCellToChiSq(cellL);
          
        } // END : left neighbour loop
        
        //if (cellR->GetNNbL() > 1) {
        //  printf("%d %d \n",iNbLchiSqMin,iCellAngDifMin);
        //}
        
        // use the angular deviation instead of the chisquare
        //iNbLchiSqMin = iCellAngDifMin;
        
        // do a second loop and take the good association of cells
        if (cellR->GetNNbL() > 0) {
          cellL = road->GetCellByGID(cellR->GetNbLgid(iNbLchiSqMin));
          addCellToTrack = kTRUE;
          addCellIdToTrack = cellR->GetNbLgid(iNbLchiSqMin);
          cellL = road->GetCellByGID(addCellIdToTrack);
          track->AddCell(cellL);
          cellL->SetUsed(kTRUE);
          if (prn) {
            printf("To track %d attach cell GID {L}: %d  - TrackID of this cell : %d - %d, Status %d \n",track->GetGID(),cellL->GetGID(),cellL->GetTrackGID(0),cellL->GetTrackGID(1),cellL->GetStatus());
          }
        }
        
      } // END : addCellToTrack
      
      // check again track length
      for (Int_t j = 0; j < 5; j++) hitSta[j] = kFALSE;
      for (Int_t iCell = 0; iCell < track->GetNcells(); iCell++) {
        cell = GetCellByGID(track->GetCellGID(iCell));
        hitSta[cell->GetLayers()[0]/2] = kTRUE;
        hitSta[cell->GetLayers()[1]/2] = kTRUE;
      }
      nHitSta = 0;
      for (Int_t i = 0; i < 5; i++) {
        if (hitSta[i]) nHitSta++;
      }
      if (nHitSta < fMinTrackLength) {
        for (Int_t iCell = 0; iCell < track->GetNcells(); iCell++) {
          cell = track->GetCell(iCell);
          cell->SetUsed(kFALSE);
        }
        RemoveLastTrack();
        trackGID--;
      } else {
        road->SetGood();
      }
      
    } // END : startLayer cells loop
    
  } // END : startLayer loop
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::RunBackward() {
  
  Bool_t prn = kFALSE;
  
  if (prn) printf("Run backward ==================================== \n");
  
  if (fMaxCellStatus == 1) return; // only isolated cells
  
  Double_t chisqNbLprev, cellAngDif, cellAngDifPrev;
  Int_t addCellIdToTrack, iNbLchiSqMin, iCellAngDifMin, trackGID = 0;
  
  AliMFTCALayer *layerR;
  AliMFTCACell *cellR;
  AliMFTCACell *cellL;
  AliMFTCACell *cell;
  AliMFTCATrack *track;
  
  Bool_t addCellToTrack;
  
  Int_t minStartLayer = fMinTrackLength-2;
  Int_t maxStartLayer = (fNlayers-1)-1;
  
  for (Int_t startLayer = maxStartLayer; startLayer >= minStartLayer; startLayer--) {
    
    if (prn) printf("Start layer %d \n",startLayer);
    
    layerR = GetLayer(startLayer);
    
    for (Int_t iCR = 0; iCR < layerR->GetNcells(); iCR++) {
      cellR = layerR->GetCell(iCR);
      if (cellR->GetStatus() == 0) continue;
      if (cellR->IsUsed()) continue;
      if (cellR->GetStatus() >= (fMinTrackLength-1)) {
        if (prn) printf("Create new track %d \n",trackGID);
        track = AddTrack(trackGID++);
        track->SetStartLayer(startLayer);
        //track->AddCellGID(cellR->GetGID());
        track->AddCell(cellR);
        cellR->SetUsed(kTRUE);
        if (prn) {
          printf("To track %d attach cell GID: %d  - TrackID of this cell : %d - %d\n",track->GetGID(),cellR->GetGID(),cellR->GetTrackGID(0),cellR->GetTrackGID(1));
        }
        // add cells to new track
        addCellToTrack = kTRUE;
        while (addCellToTrack) {
          cellR = GetCellByGID(track->GetLastCellGID());
          addCellToTrack = kFALSE;
          iNbLchiSqMin = 0;
          chisqNbLprev = 0.;
          cellAngDifPrev = -1.;
          iCellAngDifMin = 0;
          for (Int_t iNNbL = 0; iNNbL < cellR->GetNNbL(); iNNbL++) {
            cellL = GetCellByGID(cellR->GetNbLgid(iNNbL));
            if (cellL->GetStatus() == 0) continue;
            if (cellL->IsUsed()) continue;
            if (cellL->GetStatus() == (cellR->GetStatus()-1)) {
              cellAngDif = GetCellAngleDif(cellL,cellR);
              if (cellAngDifPrev < 0.) {
                cellAngDifPrev = cellAngDif;
              } else {
                if (cellAngDif < cellAngDifPrev) {
                  iCellAngDifMin = iNNbL;
                }
              }
              if (track->AddCellToChiSq(cellL) < chisqNbLprev) {
                iNbLchiSqMin = iNNbL;
              }
            }
            chisqNbLprev = track->AddCellToChiSq(cellL);
          } // END : left neighbour loop
            //if (cellR->GetNNbL() > 1) {
            //  printf("%d %d \n",iNbLchiSqMin,iCellAngDifMin);
            //}
            // use the angular deviation instead of the chisquare
            //iNbLchiSqMin = iCellAngDifMin;
          if (cellR->GetNNbL() > 0) {
            cellL = GetCellByGID(cellR->GetNbLgid(iNbLchiSqMin));
            addCellToTrack = kTRUE;
            addCellIdToTrack = cellR->GetNbLgid(iNbLchiSqMin);
            //track->AddCellGID(addCellIdToTrack);
            cellL = GetCellByGID(addCellIdToTrack);
            track->AddCell(cellL);
            cellL->SetUsed(kTRUE);
            if (prn) {
              printf("To track %d attach cell GID: %d  - TrackID of this cell : %d - %d\n",track->GetGID(),cellL->GetGID(),cellL->GetTrackGID(0),cellL->GetTrackGID(1));
            }
          }
        } // END : addCellToTrack
        
        // check again track length
        if (track->GetNcells() < (fMinTrackLength-1)) {
          for (Int_t iCell = 0; iCell < track->GetNcells(); iCell++) {
            cell = track->GetCell(iCell);
            cell->SetUsed(kFALSE);
          }
          RemoveLastTrack();
          trackGID--;
        }
        
      } // END : create new track
      
    } // END : startLayer cells loop
    
  } // END : startLayer loop
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::FilterTracks() {
  AliCodeTimerAuto("",0);

  Bool_t prn = kFALSE;
  
  AliInfo(Form("Filtering %d tracks",GetNtracks()));
  
  Int_t nTrkC = 0, nTrkG = 0, nTrkF = 0, nTrkN = 0;
  
  AliMFTCATrack *track;
  AliMFTCACell *cell, *celln;
  Int_t ndof, nptr, cellGID, cellGIDn, cellGIDprev = -1, nTrkSplitEnd = 0;
  const Int_t nMaxh = 100;
  Double_t xTr[nMaxh], yTr[nMaxh], zTr[nMaxh];
  Double_t a, ae, b, be, x0, xS, y0, yS, zmin, chisqx, chisqy;
  Double_t trkPhi, trkThe;
  Bool_t splitTrack, recTrack, cleanTrack, goodTrack, fakeTrack, noiseTrack;
  const Int_t nTrackMax = 10000;
  Int_t nrHitTrackID[nTrackMax], idHitTrackID[nTrackMax], nDiffTracks;
  Int_t idMaxHitsTrack, nMaxHits, nNoisyPix;
  for (Int_t i = 0; i < nTrackMax; i++) {
    nrHitTrackID[i] =  0;
    idHitTrackID[i] = -2;
  }
  Int_t nrAliMFTCATrackID[nTrackMax], idAliMFTCATrackID[nTrackMax], nDiffAliMFTCATracks = 0;
  for (Int_t i = 0; i < nTrackMax; i++) {
    nrAliMFTCATrackID[i] =  0;
    idAliMFTCATrackID[i] = -2;
  }
  Double_t xTrErrDet = 0.0028/TMath::Sqrt(12.);
  Double_t yTrErrDet = 0.0028/TMath::Sqrt(12.);
  Double_t xTrErrMS = 0.00055; // estimated at p = 5.5 GeV/c
  Double_t yTrErrMS = 0.00055; // estimated at p = 5.5 GeV/c
  Double_t xTrErr[nMaxh], yTrErr[nMaxh];
  for (Int_t i = 0; i < nMaxh; i++) {
    xTrErr[i] = TMath::Sqrt(xTrErrDet*xTrErrDet+xTrErrMS*xTrErrMS);
    yTrErr[i] = TMath::Sqrt(yTrErrDet*yTrErrDet+yTrErrMS*yTrErrMS);
  }
  fErrX = xTrErr[0];
  fErrY = yTrErr[0];
  
  Int_t nTotalHits = 0, nCleanTotalHits = 0;
  
  for (Int_t iTrack = 0; iTrack < GetNtracks(); iTrack++) {
    track = GetTrack(iTrack);
    nDiffTracks = 0;
    nptr = 0;
    for (Int_t iCell = 0; iCell < track->GetNcells(); iCell++) {
      cellGID = track->GetCellGID(iCell);
      cell = GetCellByGID(cellGID);
      //printf("Track %3d Cell %5d \n",iTrack,cellGID);
      // tracks split from the first (highest status) cell
      if (iCell == 0) {
        if (cellGIDprev >= 0) {
          if (cellGID == cellGIDprev) {
            // this is a split track
            nTrkSplitEnd++;
            splitTrack = kTRUE;
          } else {
            splitTrack = kFALSE;
          }
        }
        cellGIDprev = cellGID;
      }
      if (prn) {
        printf("Cell %4d from track %d ",cellGID,track->GetGID());
        printf("(%4d   %4d   %d) \n",cell->GetTrackGID(0),cell->GetTrackGID(1),cell->GetGID());
        printf("2: %f %f %f \n",cell->GetHit2()[0],cell->GetHit2()[1],cell->GetHit2()[2]);
        printf("1: %f %f %f \n",cell->GetHit1()[0],cell->GetHit1()[1],cell->GetHit1()[2]);
      }
      // extract hit x,y,z
      if (nptr == 0) {
        xTr[nptr] = cell->GetHit2()[0];
        yTr[nptr] = cell->GetHit2()[1];
        zTr[nptr] = cell->GetHit2()[2];
        nptr++;
        xTr[nptr] = cell->GetHit1()[0];
        yTr[nptr] = cell->GetHit1()[1];
        zTr[nptr] = cell->GetHit1()[2];
        nptr++;
      } else {
        xTr[nptr] = cell->GetHit1()[0];
        yTr[nptr] = cell->GetHit1()[1];
        zTr[nptr] = cell->GetHit1()[2];
        nptr++;
      }
      // count the generated tracks which contribute to this reconstructed track
      for (Int_t ihc = 0; ihc < 2; ihc++) {
        recTrack = kTRUE;
        for (Int_t idt = 0; idt < nDiffTracks; idt++) {
          if (idHitTrackID[idt] == cell->GetTrackGID(ihc)) {
            nrHitTrackID[idt]++;
            recTrack = kFALSE;
            break;
          }
        }
        if (recTrack) {
          idHitTrackID[nDiffTracks] = cell->GetTrackGID(ihc);
          nrHitTrackID[nDiffTracks] = 1;
          nDiffTracks++;
        }
      }	// cell hits
      /*
       // debug
       if (kFALSE || pTot[cell->GetTrackGID(0)] > 4.0) {
       RuleSelectCell(cell);
       if (iCell < (track->GetNcells()-1)) {
       cellGIDn = track->GetCellGID(iCell+1);
       celln = GetCellByGID(cellGIDn);
       SetDebug(1);
       RuleSelect(celln,cell);
       SetDebug(0);
       }
       }
       //
       */
    } // end cell loop
      // assert quality of the track
    cleanTrack = goodTrack = fakeTrack = noiseTrack = kFALSE;
    if (nDiffTracks == 1 && idHitTrackID[0] >= 0) {
      cleanTrack = kTRUE;
      idMaxHitsTrack = idHitTrackID[0];
    } else {
      nNoisyPix = 0;
      nMaxHits  = 0;
      for (Int_t idt = 0; idt < nDiffTracks; idt++) {
        if (idHitTrackID[idt] == -1) {
          nNoisyPix = nrHitTrackID[idt];
        } else if (nMaxHits < nrHitTrackID[idt]) {
          nMaxHits = nrHitTrackID[idt];
          idMaxHitsTrack = idHitTrackID[idt];
        }
      }
      if (GetNDet() == 5) {
        // allow one fake hit
        if (nMaxHits >= (2*track->GetNcells())-1) {
          goodTrack = kTRUE;
        } else {
          if (nNoisyPix > 0) noiseTrack = kTRUE;
          else fakeTrack = kTRUE;
        }
      }
      if (GetNDet() == 5*2) {
        // allow two fake hits
        if (nMaxHits >= (2*track->GetNcells())-2) {
          goodTrack = kTRUE;
        } else {
          if (nNoisyPix > 0) noiseTrack = kTRUE;
          else fakeTrack = kTRUE;
        }
      }
      if (GetNDet() == 6) {
        // allow one fake hit
        if (nMaxHits >= (2*track->GetNcells())-1) {
          goodTrack = kTRUE;
        } else {
          if (nNoisyPix > 0) noiseTrack = kTRUE;
          else fakeTrack = kTRUE;
        }
      }
    } // end assert track quality
    nTotalHits += nptr;
    if (cleanTrack) {
      // count the duplicated clean tracks
      recTrack = kTRUE;
      for (Int_t idt = 0; idt < nDiffAliMFTCATracks; idt++) {
        if (idAliMFTCATrackID[idt] == idMaxHitsTrack) {
          nrAliMFTCATrackID[idt]++;
          recTrack = kFALSE;
          break;
        }
      }
      if (recTrack) {
        idAliMFTCATrackID[nDiffAliMFTCATracks] = idMaxHitsTrack;
        nrAliMFTCATrackID[nDiffAliMFTCATracks] = 1;
        nDiffAliMFTCATracks++;
      }
      track->SetMCflag(1);
      track->SetMCindex(idMaxHitsTrack);
      nTrkC++;
      nCleanTotalHits += nptr;
    }
    if (goodTrack) {
      track->SetMCflag(2);
      nTrkG++;
    }
    if (fakeTrack) {
      track->SetMCflag(3);
      nTrkF++;
    }
    if (noiseTrack) {
      track->SetMCflag(4);
      nTrkN++;
    }
    // linear regression
    //printf("Fit line with %d points.\n",nptr);
    if (LinFit(nptr,zTr,xTr,xTrErr,a,ae,b,be)) {
      x0 = b; xS = a;
      if (LinFit(nptr,zTr,yTr,yTrErr,a,ae,b,be)) {
        y0 = b; yS = a;
        chisqx = 0.;
        chisqy = 0.;
        for (Int_t iptr = 0; iptr < nptr; iptr++) {
          //printf("%d  %f  %f  %f  \n",iptr,xTr[iptr],yTr[iptr],zTr[iptr]);
          chisqx += (xTr[iptr]-(xS*zTr[iptr]+x0))*(xTr[iptr]-(xS*zTr[iptr]+x0))/(xTrErr[iptr]*xTrErr[iptr]);
          chisqy += (yTr[iptr]-(yS*zTr[iptr]+y0))*(yTr[iptr]-(yS*zTr[iptr]+y0))/(yTrErr[iptr]*yTrErr[iptr]);
        }
        // track phi and theta
        trkPhi = trkThe = 0.;
        if (TMath::Abs(xS) > 0.) {
          trkPhi = TMath::ATan(yS/xS);
          // put the correct signs
          if (xS < 0. && yS > 0.) {
            trkPhi = TMath::Pi()+trkPhi;
          }
          if (xS < 0. && yS < 0.) {
            trkPhi = TMath::Pi()+trkPhi;
          }
          if (xS > 0. && yS < 0.) {
            trkPhi = TMath::TwoPi()+trkPhi;
          }
          if (TMath::Abs(TMath::Sin(trkPhi)) > 0.) {
            trkThe = TMath::ATan(yS/TMath::Sin(trkPhi));
            // put the correct signs
            trkThe = TMath::Abs(trkThe); // is always smaller than 90deg
            trkPhi *= TMath::RadToDeg();
            trkThe *= TMath::RadToDeg();
          }
        }
        // calculate DCA with the beam axis
        zmin = -(x0*xS+y0*yS)/(xS*xS+yS*yS);
        track->SetTheta(trkThe);
        track->SetPhi(trkPhi);
        if (fCalcVertex) {
          track->SetVertX(x0+xS*fZVertCalc);
          track->SetVertY(y0+yS*fZVertCalc);
	  track->SetVertZ(fZVertCalc);
        } else {
          track->SetVertX(x0+xS*fZVertDet);
          track->SetVertY(y0+yS*fZVertDet);
	  track->SetVertZ(fZVertDet);
        }
        ndof = nptr-2;
        track->SetChiSqX(chisqx/(Double_t)ndof);
        track->SetChiSqY(chisqy/(Double_t)ndof);
      } // yz fit
    } // xz fit
      //printf("End fit.\n");
  } // end track loop
  
  fNDifTracks = nDiffAliMFTCATracks;
  
  AliInfo(Form("Track found -> C: %5d G: %5d F: %5d N: %5d Dif: %5d   TotalHits: %d Clean: %d",nTrkC,nTrkG,nTrkF,nTrkN,nDiffAliMFTCATracks,nTotalHits,nCleanTotalHits));
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::DrawTracks(Double_t *pTot, Double_t *Theta) {
  
  Bool_t prn = kFALSE;
  
  AliMFTCATrack *track;
  AliMFTCACell *cell;
  Int_t cellGID, trackGID, nTrackS = 0;
  Bool_t single[10000], hitFromNoisyPix, hitFromDiffTrack;
  Int_t nGoodCell, firstHitTrackID, idTrack[10000], nGoodTracks = 0;
  for (Int_t i = 0; i < 10000; i++) idTrack[i] = -1;
  
  printf("Draw %d tracks. \n",GetNtracks());
  /*
   Double_t x[nDet], y[nDet], z[nDet], a, ae, b, be;
   Double_t errx[nDet], erry[nDet];
   for (Int_t i = 0; i < nDet; i++) {
   errx[i] = det[i]->GetSigmaX();
   erry[i] = det[i]->GetSigmaY();
   }
   */
  //Double_t chisq;
  for (Int_t iT = 0; iT < GetNtracks(); iT++) {
    // cell content in a track
    single[iT] = kTRUE;
    track = GetTrack(iT);
    cellGID = track->GetCellGID(0);
    cell = GetCellByGID(cellGID);
    trackGID = cell->GetTrackGID(0);
    single[iT] &= (cell->GetTrackGID(1) == trackGID);
    nGoodCell = 0;
    for (Int_t iC = 0; iC < track->GetNcells(); iC++) {
      cellGID = track->GetCellGID(iC);
      cell = GetCellByGID(cellGID);
      single[iT] &= (cell->GetTrackGID(0) == trackGID);
      single[iT] &= (cell->GetTrackGID(1) == trackGID);
      if (single[iT]) nGoodCell++;
      //x[iC] = cell->GetHit1()[0];
      //y[iC] = cell->GetHit1()[1];
      //z[iC] = cell->GetHit1()[2];
    }
    //____________________________________________________________________
    /*
     chisq = 0.;
     if (LinFit(track->GetNcells(),z,x,errx,a,ae,b,be)) {
     for (Int_t iC = 0; iC < track->GetNcells(); iC++) {
     chisq += (x[iC]-(a*z[iC]+b))*(x[iC]-(a*z[iC]+b));
     }
     #ifdef HOUGH
     hRTheXZ->Fill(90.-TMath::ATan(a)*TMath::RadToDeg(),b*TMath::Sin(TMath::ATan(a)));
     #endif
     }
     if (LinFit(track->GetNcells(),z,y,erry,a,ae,b,be)) {
     for (Int_t iC = 0; iC < track->GetNcells(); iC++) {
     chisq += (y[iC]-(a*z[iC]+b))*(y[iC]-(a*z[iC]+b));
     }
     #ifdef HOUGH
     hRTheYZ->Fill(90.-TMath::ATan(a)*TMath::RadToDeg(),b*TMath::Sin(TMath::ATan(a)));
     #endif
     }
     if (prn) printf("ChiSq = %e \n",chisq);
     */
    //____________________________________________________________________
    if (fDebug > 0) {
      hNGoodCell->Fill(nGoodCell);
    }
    if (single[iT]) {
      nTrackS++;
    } else {
      if (prn || fDebug > 0) {
        for (Int_t iC = 0; iC < track->GetNcells(); iC++) {
          cellGID = track->GetCellGID(iC);
          cell = GetCellByGID(cellGID);
          if (prn) printf("Track %4d Cell %2d TrackGID %5d %5d \n",iT,iC,cell->GetTrackGID(0),cell->GetTrackGID(1));
        }
      }
    }
    if (prn) printf("Track %d with %d cells: \n",iT,track->GetNcells());
    hitFromNoisyPix = kFALSE;
    hitFromDiffTrack = kFALSE;
    firstHitTrackID = -1;
    for (Int_t iC = 0; iC < track->GetNcells(); iC++) {
      cellGID = track->GetCellGID(iC);
      cell = GetCellByGID(cellGID);
      if (iC == 0) firstHitTrackID = cell->GetTrackGID(0);
      if (prn) printf("\t%5d from track: ",cellGID);
      if (prn) printf("%5d   %5d \n",cell->GetTrackGID(0),cell->GetTrackGID(1));
      if (cell->GetTrackGID(0) == -1 || cell->GetTrackGID(1) == -1)
      hitFromNoisyPix = kTRUE;
      else if (cell->GetTrackGID(0) != firstHitTrackID ||
               cell->GetTrackGID(1) != firstHitTrackID)
      hitFromDiffTrack = kTRUE;
      
      //cell->PrintCell("FULL");
    }
    if (hitFromDiffTrack) {
      hTrackType->Fill(1);
      //PrintTrack(iT);
    } else if (hitFromNoisyPix) {
      hTrackType->Fill(2);
      //PrintTrack(iT);
    } else {
      if (idTrack[firstHitTrackID] >= 0) {
        hTrackType->Fill(3);
        printf("Double track %6d from %6d and %6d \n",firstHitTrackID,iT,idTrack[firstHitTrackID]);
        PrintTrack(iT);
        PrintTrack(idTrack[firstHitTrackID]);
      } else {
        idTrack[firstHitTrackID] = iT;
        hTrackType->Fill(0);
        if (track->GetNcells() >= (fMinTrackLength-1)) nGoodTracks++;
        //PrintTrack(iT);
        //printf("%d \n",track->GetNcells());
      }
    }
  } // end loop tracks
  
  printf("... %d single. \n",nTrackS);
  printf("... %d (%d) good \n",(Int_t)hTrackType->GetBinContent(1),nGoodTracks);
  
}

//___________________________________________________________________________
AliMFTCACell *AliMFTTrackFinder::GetCellByGID(Int_t gid) {
  
  AliMFTCALayer *layer;
  AliMFTCACell *cell;
  
  for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
    layer = GetLayer(iL);
    for (Int_t iC = 0; iC < layer->GetNcells(); iC++) {
      cell = layer->GetCell(iC);
      if (gid == cell->GetGID()) return cell;
    }
  }
  
  return 0;
  
}

//___________________________________________________________________________
Bool_t AliMFTTrackFinder::RuleSelect(AliMFTCACell *cellL, AliMFTCACell *cellR) { // Are cells matching in angle ?
  /*
   if (0) {
   // ideal
   if (cellL->GetTrackGID(0) != cellL->GetTrackGID(1)) return kFALSE;
   if (cellR->GetTrackGID(0) != cellR->GetTrackGID(1)) return kFALSE;
   if (cellL->GetTrackGID(1) != cellR->GetTrackGID(0)) return kFALSE;
   }
   */
  //if (cellL->GetLayers()[1] != cellR->GetLayers()[0]) {
  //  printf("Neighbor cells do not touch the same common layer!\n");
  //}
  Int_t layer = cellL->GetLayers()[1];
  
  TVector3 *segL_ = cellL->GetSeg();
  TVector3 *segR_ = cellR->GetSeg();
  
  TVector3 segL = TVector3(*segL_);
  TVector3 segR = TVector3(*segR_);
  
  const Double_t *hitL[2];
  const Double_t *hitR[2];
  
  hitL[0] = cellL->GetHit1();
  hitL[1] = cellL->GetHit2();
  hitR[0] = cellR->GetHit1();
  hitR[1] = cellR->GetHit2();
  
  //printf("%f %f %f - %f %f %f \n",hitL[1][0],hitL[1][1],hitL[1][2],hitR[0][0],hitR[0][1],hitR[0][2]);
  
  Double_t dx, dy, a;
  dx = (hitL[1][0]-hitR[0][0])*(hitL[1][0]-hitR[0][0]); //  Distance in X direction between the 2 hits of the two different cells in the same layer
  dy = (hitL[1][1]-hitR[0][1])*(hitL[1][1]-hitR[0][1]); //  Distance in Y direction between the 2 hits of the two different cells in the same layer
  dx = (dx > 0.) ? (TMath::Sqrt(dx)) : 0.;
  dy = (dy > 0.) ? (TMath::Sqrt(dy)) : 0.;
  a = (segL.Angle(segR))*TMath::RadToDeg(); // Angle between the segments of each cell
  
  if (fDebug > 0) {
    hDA[cellL->GetLayers()[1]]->Fill(a);
    hDXY[cellL->GetLayers()[1]]->Fill(dx*1.E4,dy*1.E4); // in microns
                                                        //hDXY[cellL->GetLayers()[1]]->Fill(dx,dy);
    /*
     printf("--------------------\n");
     segL.Print();
     segR.Print();
     TVector3 segL1  = TVector3(hitL[1][0]-hitL[0][0],
     hitL[1][1]-hitL[0][1],
     hitL[1][2]-hitL[0][2]);
     TVector3 segR1  = TVector3(hitR[1][0]-hitR[0][0],
     hitR[1][1]-hitR[0][1],
     hitR[1][2]-hitR[0][2]);
     segL1.Print();
     segR1.Print();
     printf("%f %f %f \n",hitL[0][0],hitL[0][1],hitL[0][2]);
     printf("%f %f %f - %f %f %f \n",hitL[1][0],hitL[1][1],hitL[1][2],hitR[0][0],hitR[0][1],hitR[0][2]);
     printf("%f %f %f \n",hitR[1][0],hitR[1][1],hitR[1][2]);
     */
    /*
     if (cellL->GetTrackGID(0) == cellL->GetTrackGID(1) &&
     cellL->GetTrackGID(1) == cellR->GetTrackGID(0) &&
     cellR->GetTrackGID(0) == cellR->GetTrackGID(1)) {
     hDA[cellL->GetLayers()[1]]->Fill(a);
     hDXY[cellL->GetLayers()[1]]->Fill(dx,dy);
     }
     */
  }
  /*
   if (a < fACutN[layer]) {
   printf("dx, dy, a : %f %f %f \n",dx,dy,a);
   }
   */
  //printf("RS: %f %f %f %f %f %f \n",dx,fXCut,dy,fYCut,a,fACutN[layer]);
  if ((dx > fXCut) || (dy > fYCut) || (a > fACutN[layer])) return kFALSE;
  
  return kTRUE;
  
}

//___________________________________________________________________________
Bool_t AliMFTTrackFinder::RuleSelect2LayersGap(Int_t iL1, Int_t iL2, Double_t *hit1, Double_t *hit2) { // Are Cell formed by hit1 and hit2 compatible with any cell of layers 1 and 2 ?
  
  Bool_t prn = kFALSE;
  
  Bool_t findCompatibleCells = kFALSE;
  
  AliMFTCALayer *layer1 = 0;
  if (iL1 >= 0) layer1 = GetLayer(iL1);
  else return kFALSE;
  
  AliMFTCALayer *layer2 = 0;
  if (iL2 >= 0) layer2 = GetLayer(iL2);
  
  for (Int_t iCell = 0; iCell < layer1->GetNcells(); iCell++) {
    
    AliMFTCACell * cell = layer1->GetCell(iCell);
    const Double_t *hit = cell->GetHit1();
    
    // Distance in X direction between the 2 hits of the two different cells in the same layer
    Double_t dx = (hit[0]-hit1[0])*(hit[0]-hit1[0]);
    
    // Distance in Y direction between the 2 hits of the two different cells in the same layer
    Double_t dy = (hit[1]-hit1[1])*(hit[1]-hit1[1]);
    
    //Double_t radius = (dx+dy > 0.) ? (TMath::Sqrt(dx+dy)) : 0.;
    
    if ((TMath::Abs(dx) > fXCut) || (TMath::Abs(dy) > fYCut)) continue;
    
    TVector3 *seg1_ = cell->GetSeg();
    TVector3 seg1 = TVector3(*seg1_);
    
    TVector3 seg2(hit2[0]-hit1[0],hit2[1]-hit1[1],hit2[2]-hit1[2]);
    Double_t a = (seg1.Angle(seg2))*TMath::RadToDeg(); // Angle between the segments of each cell
    
    //printf(Form("Angle = %f\n",a));
    if (fDebug > 0) hAngleCells->Fill(a);
    
    if (prn) {
      printf("Compare with cell %d in layer %d \n",iCell,layer1->GetID());
    }
    
    if (a < 0.1) findCompatibleCells = kTRUE;
    
  }
  
  if (iL2 < 0) return (!findCompatibleCells);
  
  for (Int_t iCell = 0; iCell < layer2->GetNcells(); iCell++) {
    
    AliMFTCACell * cell = layer2->GetCell(iCell);
    const Double_t *hit = cell->GetHit2();
    
    // Distance in X direction between the 2 hits of the two different cells in the same layer
    Double_t dx = (hit[0]-hit2[0])*(hit[0]-hit2[0]);
    
    // Distance in Y direction between the 2 hits of the two different cells in the same layer
    Double_t dy = (hit[1]-hit2[1])*(hit[1]-hit2[1]);
    
    //Double_t radius = (dx+dy > 0.) ? (TMath::Sqrt(dx+dy)) : 0.;
    
    if ((TMath::Abs(dx) > fXCut) || (TMath::Abs(dy) > fYCut)) continue;
    
    TVector3 *seg1_ = cell->GetSeg();
    TVector3 seg1 = TVector3(*seg1_);
    
    TVector3 seg2(hit2[0]-hit1[0],hit2[1]-hit1[1],hit2[2]-hit1[2]);
    Double_t a = (seg1.Angle(seg2))*TMath::RadToDeg(); // Angle between the segments of each cell
    
    //printf(Form("Angle = %f\n",a));
    if (fDebug > 0) hAngleCells->Fill(a);
    
    if (prn) {
      printf("Compare with cell %d in layer %d \n",iCell,layer2->GetID());
    }
    
    if (a < 0.1) findCompatibleCells = kTRUE;
    
  }
  
  return (!findCompatibleCells);
  
}



//___________________________________________________________________________
Bool_t AliMFTTrackFinder::RuleSelectCell(AliMFTCACell *cell) { // Look if segment pointing to the vertex
  
  Int_t layer = cell->GetLayers()[0];
  
  TVector3 *seg_ = cell->GetSeg();
  
  TVector3 seg = TVector3(*seg_);
  
  const Double_t *hit;
  
  hit = cell->GetHit1();
  
  Double_t vert[3] = { 0., 0., 0. };
  if (fCalcVertex) vert[2] = fZVertCalc;
  else vert[2] = fZVertDet;
  
  Double_t a;
  TVector3 segV;
  
  segV = TVector3(hit[0]-vert[0],hit[1]-vert[1],hit[2]-vert[2]);
  a = (seg.Angle(segV))*TMath::RadToDeg();
  
  hThetaCells->Fill(seg.Theta()*TMath::RadToDeg());
  hDAv[layer]->Fill(a);
  
  if ( a > fACutV[layer]) return kFALSE;
  
  return kTRUE;
  
}

//___________________________________________________________________________
Bool_t AliMFTTrackFinder::RuleSelectCell(Double_t *h1, Double_t *h2, Int_t iL1, TF1 *f, Bool_t acalc) {

  // Look if segment pointing to the vertex (using directly the hits)
  
  Double_t a, av[2], acut;
  TVector3 segV, segVdet;
  
  TVector3 seg  = TVector3(h2[0]-h1[0],  h2[1]-h1[1],  h2[2]-h1[2]);
  
  Double_t vert[3] = { 0., 0., 0. };
  if (fCalcVertex) vert[2] = fZVertCalc;
  else vert[2] = fZVertDet;
  
  segVdet = TVector3(h1[0]-vert[0],h1[1]-vert[1],h1[2]-vert[2]);
  a = (seg.Angle(segVdet))*TMath::RadToDeg();
  
  //hDAv[iL1]->Fill(a);
  
  if (acalc) {
    acut = f->Eval(180.-segVdet.Theta()*TMath::RadToDeg());
  } else {
    acut = fACutV[iL1];
  }
  //printf("%f %f %d \n",180.-segVdet.Theta()*TMath::RadToDeg(),acut,acalc);
  //printf("RSC: %f %f \n",a,acut);
  
  if ( a > acut) {
    AliDebug(2,"Cell NOT Selected");
    return kFALSE;
  }
  AliDebug(2,"Cell Selected");

  return kTRUE;
  
}

//___________________________________________________________________________
AliMFTCATrack* AliMFTTrackFinder::AddTrack(Int_t gid) {
  
  new ((*fTracks)[fNtracks++]) AliMFTCATrack();
  AliMFTCATrack *track = (AliMFTCATrack*)fTracks->At(fTracks->GetLast());
  track->SetGID(gid);
  
  return track;
  
}

//___________________________________________________________________________
AliMFTCATrack* AliMFTTrackFinder::AddTrack(Int_t gid, const AliMFTCATrack& trk) {
  
  // create new track and copy
  
  new ((*fTracks)[fNtracks++]) AliMFTCATrack();
  AliMFTCATrack *track = (AliMFTCATrack*)fTracks->At(fTracks->GetLast());
  track->SetGID(gid);
  
  track->SetStartLayer(trk.GetStartLayer());
  
  for (Int_t iC = 0; iC < trk.GetNcells(); iC++) {
    //track->AddCellGID(trk.GetCellGID(iC));
    track->AddCell(trk.GetCell(iC));
  }
  
  return track;
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::UpdateCellStatusR() {
  
  AliMFTCARoad *road;
  AliMFTCACell *cell;
  
  for (Int_t ir = 0; ir < fNRoads; ir++) {
    road = GetRoad(ir);
    for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
      for (Int_t iC = 0; iC < road->GetNcellsInLayer(iL); iC++) {
        cell = road->GetCellInLayer(iL,iC);
        cell->UpdateStatus();
        fMaxCellStatus = TMath::Max(fMaxCellStatus,cell->GetStatus());
      }
    }
  }
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::UpdateCellStatus() {
  
  AliMFTCALayer *layer;
  AliMFTCACell *cell;
  
  for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
    layer = GetLayer(iL);
    for (Int_t iC = 0; iC < layer->GetNcells(); iC++) {
      cell = layer->GetCell(iC);
      cell->UpdateStatus();
      fMaxCellStatus = TMath::Max(fMaxCellStatus,cell->GetStatus());
    }
  }
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::PrintTrack(Int_t id) {
  
  AliMFTCATrack *track = GetTrack(id);
  Int_t cellGID;
  AliMFTCACell *cell;
  printf("Track:\t%6d \n",id);
  for (Int_t iC = 0; iC < track->GetNcells(); iC++) {
    cellGID = track->GetCellGID(iC);
    cell = GetCellByGID(cellGID);
    printf("cell\t%6d gid\t%6d \tfrom track: ",iC,cellGID);
    printf("\t%6d\t%6d \tin layers: ",cell->GetTrackGID(0),cell->GetTrackGID(1));
    printf("\t%1d\t%1d\n",cell->GetLayers()[0],cell->GetLayers()[1]);
  }
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::PrintAll() {
  
  AliMFTCALayer *layer;
  AliMFTCACell *cell;
  
  for (Int_t iL = 0; iL < (fNlayers-1); iL++) {
    layer = GetLayer(iL);
    printf("LayerID %d \n",layer->GetID());
    for (Int_t iC = 0; iC < layer->GetNcells(); iC++) {
      cell = layer->GetCell(iC);
      cell->PrintCell("FULL");
    }
  }
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::DrawHisto() {
  
//  TCanvas *ca1 = new TCanvas("CA1","",50,50,800,800);
//  TCanvas *ca2 = new TCanvas("CA2","",50,50,800,800);
//  TCanvas *ca3 = new TCanvas("CA3","",50,50,800,800);
//  ca1->Divide(3,2);
//  ca2->Divide(3,2);
//  for (Int_t i = 0; i<fNlayers; i++) {
//    ca1->cd(i+1);
//    hDA[i]->Draw();
//    ca2->cd(i+1);
//    hDXY[i]->Draw("colz");
//  }
//  ca3->Divide(1,2);
//  ca3->cd(1);
//  hNGoodCell->Draw();
//  ca3->cd(2);
//  hTrackType->Draw();
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::PrintParam() {
  
  printf("==== AliMFTTrackFinder::PrintParam ====\n");
  
  printf("Number of layers: %d \n",fNlayers);
  printf("Use the TrackFinder: %d \n",fUseTF);
  printf("Layer thickness in X0: %5.3f \n",fThick);
  printf("Pixel noise: %4.2e \n",fPixelNoise);
  printf("Add noise: %d \n",fAddNoise);
  printf("fMinTrackLength: %d \n",fMinTrackLength);
  printf("fXCut: %6.4f cm\n",fXCut);
  printf("fYCut: %6.4f cm\n",fYCut);
  printf("fMaxSegAngle: %4.1f deg\n",fMaxSegAngle);
  for (Int_t i = 0; i < fNlayers; i++) {
    printf("fACutV[%d]: %4.2f \n",i,fACutV[i]);
  }
  for (Int_t i = 0; i < fNlayers; i++) {
    printf("fACutN[%d]: %4.2f \n",i,fACutN[i]);
  }
  for (Int_t i = 0; i < fNlayers; i++) {
    printf("PlaneDetEff[%d]: %4.2f \n",i,fPlaneDetEff[i]);
  }
  printf("Calculate vertex: %d \n",fCalcVertex);
  
  printf("CPU time: %f seconds\n",fCPUTime);
  printf("Real time: %f seconds\n",fRealTime);
  
  printf("============================\n");
  
}

//___________________________________________________________________________
AliMFTCARoad* AliMFTTrackFinder::AddRoad() {
  
  new ((*fRoads)[fNRoads++]) AliMFTCARoad();
  AliMFTCARoad *road = (AliMFTCARoad*)fRoads->At(fRoads->GetLast());
  
  return road;
  
}

//___________________________________________________________________________
void AliMFTTrackFinder::BuildRoads() {
  AliCodeTimerAuto("",0);

  Bool_t prn = kFALSE;
  
  // planes: 0, 1, 2, ..., 9 (layer ... becomes plane)
  // rules for combining first/last plane in a road:
  // 0 with 6, 7, 8, 9
  // 1 with 6, 7, 8, 9
  // 2 with 8, 9
  // 3 with 8, 9
  
  Int_t iPla1Min = 0, iPla1Max = 3;
  Int_t iPla2Min[4] = { 6, 6, 6, 6 };
  Int_t iPla2Max[4] = { 9, 9, 7, 7 };
  
  Int_t nH1, nH2, nH;
  Int_t roadLen, nRoads = 0, trackGID = 0;
  Double_t h1[3], h2[3], h[3], hx, hy, dR;
  
  Double_t dRmin = 0.0400;
  
  AliMFTCAHit *hit1, *hit2, *hit, *htmp;
  AliMFTCARoad *road, *road1;
  Bool_t hitSta[5] = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
  Bool_t isUsed, newRoad;
  
  for (Int_t iL1 = iPla1Min; iL1 <= iPla1Max; iL1++) {
    
    for (Int_t iL2 = iPla2Max[iL1]; iL2 >= iPla2Min[iL1]; iL2--) {
      
      // see AliMFTCARoad::SetLength
      roadLen = iL2/2-iL1/2;
      
      nH1 = GetLayer(iL1)->GetNhits();
      if(prn) AliInfo(Form("nH1 = %d ",nH1));
      for (Int_t iH1 = 0; iH1 < nH1; iH1++) {
        //AliInfo(Form("iH1 = %d ",iH1));

        hit1 = GetLayer(iL1)->GetHit(iH1);
        
        if (hit1->IsUsed()) continue;
        /*
         // check if it belongs to a good longer road previously found
         isUsed = kFALSE;
         for (Int_t i = 0; i < hit1->GetNRoads(); i++) {
         road1 = GetRoad(hit1->GetInRoad(i));
         if (road1->IsGood()) {
         if (road1->GetLength() >= roadLen) {
         isUsed = kTRUE;
         break;
         }
         }
         }
         //if (isUsed) continue;
         */
        if (prn)
        printf("Hit1: %d %d %d %d \n",hit1->GetLayer(),hit1->GetID(),hit1->GetDetElemID(),hit1->GetTrackGID());
        
        h1[0] = hit1->GetPos()[0];
        h1[1] = hit1->GetPos()[1];
        h1[2] = hit1->GetPos()[2];
        
        nH2 = GetLayer(iL2)->GetNhits();
        if(prn) AliInfo(Form("nH2 = %d ",nH2));

        for (Int_t iH2 = 0; iH2 < nH2; iH2++) {
          
          hit2 = GetLayer(iL2)->GetHit(iH2);
          
          if (hit2->IsUsed()) continue;
          /*
           // check if it belongs to a good longer road previously found
           isUsed = kFALSE;
           for (Int_t i = 0; i < hit2->GetNRoads(); i++) {
           road1 = GetRoad(hit2->GetInRoad(i));
           if (road1->IsGood()) {
           if (road1->GetLength() >= roadLen) {
           isUsed = kTRUE;
           break;
           }
           }
           }
           //if (isUsed) continue;
           */
          if (prn)
          printf("Hit2: %d %d %d %d \n",hit2->GetLayer(),hit2->GetID(),hit2->GetDetElemID(),hit2->GetTrackGID());
          
          h2[0] = hit2->GetPos()[0];
          h2[1] = hit2->GetPos()[1];
          h2[2] = hit2->GetPos()[2];
          
          TVector3 vec(h2[0]-h1[0],h2[1]-h1[1],h2[2]-h1[2]);
          if (vec.Theta() < fMaxSegAngle*TMath::DegToRad()) continue;
          
          if (!RuleSelectCell(h1,h2,iL1)) continue;
          
          newRoad = kTRUE;
          for (Int_t i = 0; i < 5; i++) hitSta[i] = kFALSE;
          
          for (Int_t iL = (iL1+1); iL <= (iL2-1); iL++) {
            
            nH = GetLayer(iL)->GetNhits();
            
            for (Int_t iH = 0; iH < nH; iH++) {
              
              hit = GetLayer(iL)->GetHit(iH);
              
              if (hit->IsUsed()) continue;
              
              h[0] = hit->GetPos()[0];
              h[1] = hit->GetPos()[1];
              h[2] = hit->GetPos()[2];
              
              hx = h1[0] + (h2[0]-h1[0])*(h[2]-h1[2])/(h2[2]-h1[2]);
              hy = h1[1] + (h2[1]-h1[1])*(h[2]-h1[2])/(h2[2]-h1[2]);
              
              dR = TMath::Sqrt((hx-h[0])*(hx-h[0])+(hy-h[1])*(hy-h[1]));
              if (dR >= dRmin) continue;
              /*
               // check if it belongs to a good longer road previously found
               isUsed = kFALSE;
               for (Int_t i = 0; i < hit->GetNRoads(); i++) {
               AliMFTCARoad *road1 = GetRoad(hit->GetInRoad(i));
               if (road1->IsGood()) {
               if (road1->GetLength() > roadLen) {
               isUsed = kTRUE;
               break;
               }
               }
               }
               //if (isUsed) continue;
               */
              if (prn)
              printf("Hit: %d %d %d %d \n",hit->GetLayer(),hit->GetID(),hit->GetDetElemID(),hit->GetTrackGID());
              
              if (newRoad) {
                //AliInfo("Adding new road");
                road = AddRoad();
                road->SetID(nRoads++);
                road->AddHit(hit1);
                road->AddHit(hit2);
                
                hit1->SetInRoad(road->GetID());
                hit2->SetInRoad(road->GetID());
                
                hitSta[iL1/2] = kTRUE;
                hitSta[iL2/2] = kTRUE;
                
                road->SetLength(iL1,iL2);
                
                newRoad = kFALSE;
                
              }
              
              road->AddHit(hit);
              hit->SetInRoad(road->GetID());
              hitSta[iL/2] = kTRUE;
              
            } // end loop hits middle plane
            
          } // end loop middle plane
          
          // count the number of hit stations (disks)
          if (!newRoad) {
            Int_t nHitSta = 0;
            for (Int_t i = 0; i < 5; i++) 
            if (hitSta[i]) nHitSta++;
            road->SetNHitSta(nHitSta);
            if (nHitSta >= fMinTrackLength) {
              CreateCellsR(road);
              RunForwardR(road,trackGID);
              if (road->IsGood()) {
                for (Int_t j = 0; j < fNlayers; j++) {
                  for (Int_t k = 0; k < road->GetNhitsInLayer(j); k++) {
                    hit = road->GetHitInLayer(j,k);
                    //printf("%d %d %d \n",hit->GetLayer(),hit->GetID(),GetLayer(hit->GetLayer())->GetNhits());
                    htmp = GetLayer(hit->GetLayer())->GetHit(hit->GetID());
                    htmp->SetUsed();
                    if (prn)
                    printf("Hit used: %d %d %d %d \n",hit->GetLayer(),hit->GetID(),hit->GetDetElemID(),hit->GetTrackGID());
                  }
                }
              }
            }
          }
          
        } // end loop hits last plane
        
      } // end loop hits first plane
      
    } // end loop last plane
    
  } // end loop first plane
  /*
   Int_t nRoadsGood = 0, nTotalHits = 0;
   for (Int_t i = 0; i < nRoads; i++) {
   road = GetRoad(i);
   if (road->IsGood()) {
   //printf("Road %d with %d hits.\n",i,road->GetNhits());
   nRoadsGood++;
   nTotalHits += road->GetNhits();
   }
   }
   printf("Found %d roads %d good with %d hits. \n",nRoads,nRoadsGood,nTotalHits);
   */
}

//___________________________________________________________________________
void AliMFTTrackFinder::FindTracks() {
  AliCodeTimerAuto("",0);

  Bool_t prn = kFALSE;
  
  Double_t h1[3], h2[3], h[3], hx, hy;
  Double_t htr1[3], htr2[3];
  const Int_t nMaxh = 100;
  Double_t trX[nMaxh], trY[nMaxh], trZ[nMaxh];
  Int_t lay[nMaxh], trkid[nMaxh], hit[nMaxh], detid[nMaxh];
  Int_t nH1, nH2, nH, iL1, iL2, lay1, lay2, trkid1, trkid2, nptr;
  Int_t hit1, hit2, detid1, detid2, nHitSta;
  Int_t mftClsId1, mftClsId2;
  Int_t trackGID = 0;
  AliMFTCACell *cell;
  AliMFTCATrack *track;
  
  Double_t dR, dRmin, dRcut = 0.0100;
  
  TF1 *fACutF = new TF1("fACutF","[0]+x*[1]",0.,1.);
  Float_t cut170 = 0.6; // cut [deg] at theta = 170 deg
  Float_t cut177 = 0.3; // cut [deg] at theta = 177 deg
  Float_t par[2];
  par[1] = (cut177-cut170)/(177.-170.);
  par[0] = cut170-par[1]*170.;
  fACutF->SetParameter(0,par[0]);
  fACutF->SetParameter(1,par[1]);
  
  Bool_t addHit, selCAgapCell;
  Bool_t hitSta[5];
  Bool_t seed = kTRUE;
  
  Int_t step = 0;
  
  iL1 = 0;
  
  while (seed) {
    
    if (step == 0) {
      iL2 = fNlayers-1;
    } else {
      iL2--;
    }
    
    step++;
    
    if (iL2 < iL1 + (fMinTrackLength-1)) {
      iL1++;
      if (iL1 > fNlayers - (fMinTrackLength-1)) break;
      step = 0;
      continue;
    }
    
    if (prn) printf("iL1,iL2 %d %d \n",iL1,iL2);
    //continue;
    
    nH1 = GetLayer(iL1)->GetNhits();
    nH2 = GetLayer(iL2)->GetNhits();
    if (prn) printf("nH1 %d nH2 %d\n",nH1,nH2);
    for (Int_t iH1 = 0; iH1 < nH1; iH1++) {
      
      if (GetLayer(iL1)->GetHit(iH1)->IsUsed()) continue;
      
      
      h1[0] = GetLayer(iL1)->GetHit(iH1)->GetPos()[0];
      h1[1] = GetLayer(iL1)->GetHit(iH1)->GetPos()[1];
      h1[2] = GetLayer(iL1)->GetHit(iH1)->GetPos()[2];
      
      if (prn) printf("H1 %d (%.1f,%.1f,%.1f)\n",iH1,h1[0],h1[1],h1[2]);

      for (Int_t iH2 = 0; iH2 < nH2; iH2++) {
        
        if (GetLayer(iL2)->GetHit(iH2)->IsUsed()) continue;
        
        
        h2[0] = GetLayer(iL2)->GetHit(iH2)->GetPos()[0];
        h2[1] = GetLayer(iL2)->GetHit(iH2)->GetPos()[1];
        h2[2] = GetLayer(iL2)->GetHit(iH2)->GetPos()[2];
       
        TVector3 vec(h2[0]-h1[0],h2[1]-h1[1],h2[2]-h1[2]);

        if (vec.Theta() < fMaxSegAngle*TMath::DegToRad()) continue;
        //if (prn) printf("Theta =%f  (max = %f)\n",vec.Theta(),fMaxSegAngle*TMath::DegToRad());

        //if (!RuleSelectCell(h1,h2,iL1,fACutF,kTRUE)) continue;
        if (!RuleSelectCell(h1,h2,iL1)) continue;
        if (prn) printf("H2 %d (%.1f,%.1f,%.1f)\n",iH2,h2[0],h2[1],h2[2]);

        // this is a seed connecting the first and the last layers
        
        for (Int_t i = 0; i < 5; i++) hitSta[i] = kFALSE;
        nHitSta = 0;
        
        hitSta[iL1/2] = kTRUE;
        hitSta[iL2/2] = kTRUE;
        
        nptr = 0;
        
        trX[nptr] = h1[0];
        trY[nptr] = h1[1];
        trZ[nptr] = h1[2];
        lay[nptr] = iL1;
        hit[nptr] = iH1;
        trkid[nptr] = GetLayer(iL1)->GetHit(iH1)->GetTrackGID();
        detid[nptr] = GetLayer(iL1)->GetHit(iH1)->GetDetElemID();
        nptr++;
        
        for (Int_t iL = (iL1+1); iL <= (iL2-1); iL++) {
          
          nH = GetLayer(iL)->GetNhits();
          
          if (prn) printf("L %d nH %d \n",iL,nH);
          
          dRmin = dRcut;
          addHit = kFALSE;
          
          for (Int_t iH = 0; iH < nH; iH++) {
            
            //if (prn) printf("H %d \n",iH);
            
            if (GetLayer(iL)->GetHit(iH)->IsUsed()) continue;
            
            h[0] = GetLayer(iL)->GetHit(iH)->GetPos()[0];
            h[1] = GetLayer(iL)->GetHit(iH)->GetPos()[1];
            h[2] = GetLayer(iL)->GetHit(iH)->GetPos()[2];
            
            hx = h1[0] + (h2[0]-h1[0])*(h[2]-h1[2])/(h2[2]-h1[2]);
            hy = h1[1] + (h2[1]-h1[1])*(h[2]-h1[2])/(h2[2]-h1[2]);
            
            dR = TMath::Sqrt((hx-h[0])*(hx-h[0])+(hy-h[1])*(hy-h[1]));
            if (dR >= dRmin) continue;
            AliDebug(1,Form("Hit %d added",iH));
            hitSta[iL/2] = kTRUE;
            
            dRmin = dR;
            
            trX[nptr] = h[0];
            trY[nptr] = h[1];
            trZ[nptr] = h[2];
            lay[nptr] = iL;
            hit[nptr] = iH;
            trkid[nptr] = GetLayer(iL)->GetHit(iH)->GetTrackGID();
            detid[nptr] = GetLayer(iL)->GetHit(iH)->GetDetElemID();
            
            addHit = kTRUE;
            
          } // loop hits intermediate layer
          
          if (addHit) nptr++;
          
        } // loop intermediate layers
        
        trX[nptr] = h2[0];
        trY[nptr] = h2[1];
        trZ[nptr] = h2[2];
        lay[nptr] = iL2;
        hit[nptr] = iH2;
        trkid[nptr] = GetLayer(iL2)->GetHit(iH2)->GetTrackGID();
        detid[nptr] = GetLayer(iL2)->GetHit(iH2)->GetDetElemID();
        nptr++;
        
        // create cells and tracks
        //
        // calculate how many stations have hit(s)
        AliDebug(2,Form("nptr %d fMinTrackLength %d",nptr,fMinTrackLength));
        if (nptr < fMinTrackLength) continue;
        nHitSta = 0;
        for (Int_t i = 0; i < 5; i++) {
          if (hitSta[i]) nHitSta++;
        }
        AliDebug(2,Form("nHitSta %d fMinTrackLength %d",nHitSta,fMinTrackLength));
       if (nHitSta < fMinTrackLength) continue;
        /*      
         // as in CA impose gap cells with no more than one layer missing
         selCAgapCell = kTRUE;
         for (Int_t i = 0; i < (nptr-1); i++) {
         if (lay[i+1]-lay[i] > 2) {
         selCAgapCell = kFALSE;
         break;
         }
         }
         if (!selCAgapCell) continue;
         */
        // create track
        AliDebug(1,Form("Adding Track %d ",trackGID));
        track = AddTrack(trackGID++);
        track->SetStartLayer(iL1);
        
        // loop over hits in reverse order, like in RunBackward
        for (Int_t iptr = (nptr-1); iptr >= 1; iptr--) {
          
          trkid1 = trkid[iptr-1];
          lay1   = lay[iptr-1];
          hit1   = hit[iptr-1];
          detid1 = detid[iptr-1];
          
          htr1[0] = GetLayer(lay1)->GetHit(hit1)->GetPos()[0];
          htr1[1] = GetLayer(lay1)->GetHit(hit1)->GetPos()[1];
          htr1[2] = GetLayer(lay1)->GetHit(hit1)->GetPos()[2];
	  mftClsId1 = GetLayer(lay1)->GetHit(hit1)->GetMFTClsId();
          
          GetLayer(lay1)->GetHit(hit1)->SetUsed();
          
          trkid2 = trkid[iptr];
          lay2   = lay[iptr];
          hit2   = hit[iptr];
          detid2 = detid[iptr];
          
          htr2[0] = GetLayer(lay2)->GetHit(hit2)->GetPos()[0];
          htr2[1] = GetLayer(lay2)->GetHit(hit2)->GetPos()[1];
          htr2[2] = GetLayer(lay2)->GetHit(hit2)->GetPos()[2];
	  mftClsId2 = GetLayer(lay2)->GetHit(hit2)->GetMFTClsId();
          
          GetLayer(lay2)->GetHit(hit2)->SetUsed();

          // create cells
          cell = GetLayer(lay[iptr-1])->AddCell();
          cell->SetHits(htr1,htr2,fPlanesZ[lay1],fPlanesZ[lay2]);
	  cell->SetMFTClsId(mftClsId1,mftClsId2);
          cell->SetLayers(lay1,lay2);
          cell->SetStatus(1);
          cell->SetGID(fCellGID++,trkid1,trkid2);
          cell->SetDetElemID(detid1,detid2);
          
          // add cell to track
          track->AddCell(cell);
          cell->SetUsed(kTRUE);
          
        }
        
      } // loop hits last layer
    } // loop hits first layer
    
  } // end seed
  
}




//___________________________________________________________________________
Bool_t AliMFTTrackFinder::LinFit(Int_t nDet, Double_t *xcl,
              Double_t *ycl, Double_t *yerr,
              Double_t &a, Double_t &ae, Double_t &b, Double_t &be,
              Int_t skip) {
  
  
  // y=a*x+b
  
  const Int_t nMaxh = 100;
  Double_t xCl[nMaxh], yCl[nMaxh], yErr[nMaxh];
  Int_t idet = 0;
  for (Int_t i = 0; i < nDet; i++) {
    if (i == skip) continue;
    xCl[idet] = xcl[i];
    yCl[idet] = ycl[i];
    yErr[idet] = yerr[i];
    idet++;
  }
  
  Double_t S1, SXY, SX, SY, SXX, SsXY, SsXX, SsYY, Xm, Ym, s, delta, difx;
  
  S1 = SXY = SX = SY = SXX = 0.0;
  SsXX = SsYY = SsXY = Xm = Ym = 0.;
  difx = 0.;
  for (Int_t i = 0; i < idet; i++) {
    S1  += 1.0/(yErr[i]*yErr[i]);
    SXY += xCl[i]*yCl[i]/(yErr[i]*yErr[i]);
    SX  += xCl[i]/(yErr[i]*yErr[i]);
    SY  += yCl[i]/(yErr[i]*yErr[i]);
    SXX += xCl[i]*xCl[i]/(yErr[i]*yErr[i]);
    if (i > 0) difx += TMath::Abs(xCl[i]-xCl[i-1]);
    Xm  += xCl[i];
    Ym  += yCl[i];
    SsXX += xCl[i]*xCl[i];
    SsYY += yCl[i]*yCl[i];
    SsXY += xCl[i]*yCl[i];
  }
  delta = SXX*S1 - SX*SX;
  if (delta == 0.) {
    return kFALSE;
  }
  a = (SXY*S1 - SX*SY)/delta;
  b = (SY*SXX - SX*SXY)/delta;
  
  Ym /= (Double_t)idet;
  Xm /= (Double_t)idet;
  SsYY -= (Double_t)idet*(Ym*Ym);
  SsXX -= (Double_t)idet*(Xm*Xm);
  SsXY -= (Double_t)idet*(Ym*Xm);
  Double_t eps = 1.E-24;
  if ((idet > 2) && (TMath::Abs(difx) > eps) && ((SsYY-(SsXY*SsXY)/SsXX) > 0.)) {
    s = TMath::Sqrt((SsYY-(SsXY*SsXY)/SsXX)/(idet-2));
    be = s*TMath::Sqrt(1./(Double_t)idet+(Xm*Xm)/SsXX);
    ae = s/TMath::Sqrt(SsXX);
  } else {
    be = 0.;
    ae = 0.;
  }
  
  return kTRUE;
  
}

