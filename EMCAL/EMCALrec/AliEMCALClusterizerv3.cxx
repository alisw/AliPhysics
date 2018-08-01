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

// --- ROOT system ---
#include <TFile.h> 
#include <TMath.h> 
#include <TMinuit.h>
#include <TTree.h> 
#include <TBenchmark.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TList.h>
#include <TClonesArray.h>

// --- Standard library ---
#include <cstring>
#include <cassert>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliEMCALClusterizerv3.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALCalibTime.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALUnfolding.h"

/// \cond CLASSIMP
ClassImp(AliEMCALClusterFinder);
/// \endcond


///
/// Constructor
//____________________________________________________________________________
AliEMCALClusterFinder::AliEMCALClusterFinder(TObjArray* outputArray, AliEMCALGeometry* geometry, Double_t timeCut, Double_t timeMin, Double_t timeMax, Double_t gradientCut, Bool_t doEnergyGradientCut, Double_t thresholdSeedE, Double_t thresholdCellE) :
    fEMCALGeometry(geometry), fSeedList(), fDigitMap(), fCellMask(),
    fFoundClusters(outputArray), fNumFoundClusters(0), fTimeCut(timeCut), fTimeMin(timeMin), fTimeMax(timeMax), fGradientCut(gradientCut), fDoEnergyGradientCut(doEnergyGradientCut), fThresholdSeedEnergy(thresholdSeedE),fThresholdCellEnergy(thresholdCellE)
{
}

///
/// Destructor
//____________________________________________________________________________
AliEMCALClusterFinder::~AliEMCALClusterFinder()
{
  if(fFoundClusters) fFoundClusters->Delete();
}

///
/// Recursively search for neighbours (EMCAL)
//____________________________________________________________________________
AliEMCALRecPoint* AliEMCALClusterFinder::GetClusterFromNeighbours(AliEMCALRecPoint* recPoint, Int_t row, Int_t column)
{
  // Cluster/recpoint is already initialized with seed digit
  Bool_t shared = kFALSE;

  // Recursion 0
  if(!recPoint)
  {
    recPoint = new AliEMCALRecPoint("");
    recPoint->SetClusterType(AliVCluster::kEMCALClusterv1);
    recPoint->AddDigit(*fDigitMap[row][column], fDigitMap[row][column]->GetCalibAmp(), kFALSE);
  }
  Int_t currentSM = fEMCALGeometry->GetSuperModuleNumber(fDigitMap[row][column]->GetId());

  // Mark the current cell as clustered
  fCellMask[row][column] = kTRUE;

  // Now go recursively to the next 4 neighbours and add them to the cluster if they fulfill the conditions
  for(Int_t dir=0; dir<4; dir++)
  {
    Int_t rowDiff = 0;
    Int_t colDiff = 0;
    if(dir==0)      {rowDiff = -1; colDiff =  0;}
    else if(dir==1) {rowDiff =  0; colDiff = -1;}
    else if(dir==2) {rowDiff =  0; colDiff = +1;}
    else if(dir==3) {rowDiff = +1; colDiff =  0;}

    if( (row+rowDiff < 0) || (row+rowDiff >= EMCALClusterFinder::kNrows) ) continue;
    if( (column+colDiff < 0) || (column+colDiff >= EMCALClusterFinder::kNcolumns) ) continue;

    if(fDigitMap[row+rowDiff][column+colDiff])
      if(!fCellMask[row+rowDiff][column+colDiff])
        if (fDoEnergyGradientCut && not (fDigitMap[row+rowDiff][column+colDiff]->GetCalibAmp()>fDigitMap[row][column]->GetCalibAmp()+fGradientCut))
          if (not (TMath::Abs(fDigitMap[row+rowDiff][column+colDiff]->GetTime() - fDigitMap[row][column]->GetTime()) > fTimeCut))
          {
            // Check if cluster extends over two SMs (= shared cluster)
            shared = (currentSM != fEMCALGeometry->GetSuperModuleNumber(fDigitMap[row+rowDiff][column+colDiff]->GetId()));
            (void) GetClusterFromNeighbours(recPoint, row+rowDiff, column+colDiff);
            // Add the digit to the current cluster -- if we end up here, the selected cluster fulfills the condition
            recPoint->AddDigit(*fDigitMap[row+rowDiff][column+colDiff], fDigitMap[row+rowDiff][column+colDiff]->GetCalibAmp(), shared);
            // Mask the cell as already used for clustering
          }
  }
  return recPoint;
}


///
/// Get row (phi) and column (eta) of a digit, values corresponding to topology
///
//____________________________________________________________________________
void AliEMCALClusterFinder::GetTopologicalRowColumn(AliEMCALDigit* digit, Int_t& row, Int_t& column)
{
  Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  row=0;
  column=0;
  // Get SM number and relative row/column for SM
  fEMCALGeometry->GetCellIndex(digit->GetId(), nSupMod,nModule,nIphi,nIeta);
  fEMCALGeometry->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, row,column);

  // Add shifts wrt. supermodule and type of calorimeter
  // NOTE:
  // * Rows (phi) are arranged that one space is left empty between supermodules in phi
  //   This is due to the physical gap that forbids clustering
  // * For DCAL, there is an additional empty column between two supermodules in eta
  //   Again, this is to account for the gap in DCAL

  row    += nSupMod/2 * (24+1);
  // In DCAL, leave a gap between two SMs with same phi
  if(!fEMCALGeometry->IsDCALSM(nSupMod)) // EMCAL
    column += nSupMod%2 * 48;
  else
    column += nSupMod%2 * (48+1);
}


///
/// Return number of found clusters. Start clustering from highest energy cell.
//____________________________________________________________________________
Int_t AliEMCALClusterFinder::FindClusters(TClonesArray* digitArray)
{  
  fFoundClusters->Delete();

  // Algorithm
  // - Fill digits in 2D topological map
  // - Fill struct arrays (energy,x,y)  (to get mapping energy -> (x,y))
  // - Create 2D bitmap (digit is already clustered or not)
  // - Sort struct arrays with descending energy
  //
  // - Loop over arrays:
  // --> Check 2D bitmap (don't use digit which are already clustered)
  // --> Take valid digit with highest energy as seed (they are already sorted)
  // --> Recursive to neighboughs and create cluster
  // --> Seed cell and all neighbours belonging to cluster will be put in 2D bitmap

  // Reset digit maps and cell masks
  std::memset(fCellMask, 0, sizeof(Bool_t) * EMCALClusterFinder::kNrows*EMCALClusterFinder::kNcolumns);
  std::memset(fDigitMap, 0, sizeof(AliEMCALDigit*) * EMCALClusterFinder::kNrows*EMCALClusterFinder::kNcolumns);

  // Calibrate digits and fill the maps/arrays
  Int_t nCells = 0;
  Double_t ehs = 0.0;
  AliEMCALDigit *digit = 0;
  TIter nextdigit(digitArray);
  while ( (digit = static_cast<AliEMCALDigit*>(nextdigit())) ) 
  {
    Float_t dEnergyCalibrated = digit->GetAmplitude();
    Float_t time              = digit->GetTime();

    if (dEnergyCalibrated < fThresholdCellEnergy || time > fTimeMax || time < fTimeMin) 
      continue;
    if (!fEMCALGeometry->CheckAbsCellId(digit->GetId()))
      continue;

    ehs += dEnergyCalibrated;

    // Put digit to 2D map
    Int_t row = 0, column = 0;
    GetTopologicalRowColumn(digit, row, column);
    fDigitMap[row][column] = digit;
    fSeedList[nCells].energy = dEnergyCalibrated;
    fSeedList[nCells].row = row;
    fSeedList[nCells].column = column;
    nCells++;
  }

  // Sort struct arrays with ascending energy
  std::sort(fSeedList, fSeedList+nCells);

  // Take next valid digit in calorimeter as seed (in descending energy order)
  fNumFoundClusters = 0;
  for(Int_t i=nCells-1; i>=0; i--)
  {
    Int_t row = fSeedList[i].row, column = fSeedList[i].column;
    // Continue if the cell is already masked (i.e. was already clustered)
    if(fCellMask[row][column])
      continue;
    // Continue if energy constraints are not fulfilled
    if (fSeedList[i].energy<=fThresholdSeedEnergy)
      continue;

    // Seed is found, form cluster recursively
    if (fNumFoundClusters>=fFoundClusters->GetSize()) 
      fFoundClusters->Expand(2*fNumFoundClusters+1);

    AliEMCALRecPoint* recPoint = GetClusterFromNeighbours(0, row, column);
    fFoundClusters->AddAt(recPoint, fNumFoundClusters++);
  }

  AliDebug(1,Form("%d clusters found from %d digits (total=%d)-> ehs %f (minE %f)\n",fNumFoundClusters,nCells,digitArray->GetEntriesFast(), ehs,fThresholdCellEnergy));
  return fNumFoundClusters;
}



/// \cond CLASSIMP
ClassImp(AliEMCALClusterizerv3) ;
/// \endcond


///
/// Default constructor
//____________________________________________________________________________
AliEMCALClusterizerv3::AliEMCALClusterizerv3() 
  : AliEMCALClusterizerv1(), fDoEnGradCut(1), fClusterFinder(0)
{ }

///
/// Constructor 
/// 
/// \param geometry: EMCal geometry pointer
//____________________________________________________________________________
AliEMCALClusterizerv3::AliEMCALClusterizerv3(AliEMCALGeometry* geometry)
  : AliEMCALClusterizerv1(geometry), fDoEnGradCut(1), fClusterFinder(0)
{ }

///
/// Constructor, geometry and calibration are initialized elsewhere.
///
/// \param geometry: EMCal geometry pointer
/// \param calib: EMCal energy calibration container
/// \param calibt: EMCal time calibration container
/// \param caloped: EMCal bad map container
//____________________________________________________________________________
AliEMCALClusterizerv3::AliEMCALClusterizerv3(AliEMCALGeometry* geometry, 
                                             AliEMCALCalibData* calib, 
                                             AliEMCALCalibTime* calibt, 
                                             AliCaloCalibPedestal* caloped)
  : AliEMCALClusterizerv1(geometry, calib, calibt, caloped), fDoEnGradCut(1), fClusterFinder(0)
{ }

///
/// Destructtor
//____________________________________________________________________________
AliEMCALClusterizerv3::~AliEMCALClusterizerv3()
{
  if(fClusterFinder) delete fClusterFinder;
}

///
/// Make list of clusters. Use AliEMCALClusterFinder class
//____________________________________________________________________________
void AliEMCALClusterizerv3::MakeClusters()
{  
  if(fGeom==0)
    AliFatal("AliEMCALGeometry object not properly loaded.");

  // Create cluster finder only once, at execution time (cannot be done in constructor, because settings could change afterwards)
  // fRecPoints is the TObjArray we will write to
  if(!fClusterFinder)
    fClusterFinder = new AliEMCALClusterFinder(fRecPoints, fGeom, fTimeCut, fTimeMin, fTimeMax, fECALocMaxCut, fDoEnGradCut, fECAClusteringThreshold, fMinECut);

  // Do calibration of cells and sort out some cells already
  AliEMCALDigit *digit = 0;
  TIter nextdigit(fDigitsArr);
  while ( (digit = static_cast<AliEMCALDigit*>(nextdigit())) ) 
  {
    Float_t dEnergyCalibrated = digit->GetAmplitude();
    Float_t time              = digit->GetTime();
    Calibrate(dEnergyCalibrated, time, digit->GetId());
    
    digit->SetCalibAmp(dEnergyCalibrated);
    digit->SetTime(time);
  }

  // Perform cluster finding. Output is written to fRecPoints
  fNumberOfECAClusters = fClusterFinder->FindClusters(fDigitsArr);
}

