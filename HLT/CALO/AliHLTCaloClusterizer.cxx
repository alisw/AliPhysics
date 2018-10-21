// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/**
 * @file   AliHLTCaloClusterizer.cxx
 * @author Rudiger Haake (Yale)
 * @date
 * @brief  Clusterizer for EMCAL/DCAL (HLT implementation)
 */

// see header file for class documentation

#include "AliHLTCaloClusterizer.h"
#include "AliHLTLogging.h"
#include "TMath.h"
#include "AliHLTCaloRecPointDataStruct.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloConstantsHandler.h"
#include "AliEMCALGeometry.h"
#include "AliHLTEMCALGeometry.h"
#include "AliEMCALCalibBadChannels.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "TH2.h"
#include "TMath.h"


/// \cond CLASSIMP
ClassImp(AliHLTClusterFinder);
/// \endcond


///
/// Constructor
//____________________________________________________________________________
AliHLTClusterFinder::AliHLTClusterFinder(AliHLTCaloRecPointDataStruct** outputArray, Int_t maxNumClusters, AliHLTEMCALGeometry* geometry, Double_t timeCut, Double_t timeMin, Double_t timeMax, Double_t gradientCut, Bool_t doEnergyGradientCut, Double_t thresholdSeedE, Double_t thresholdCellE) :
    fSeedList(), fDigitMap(), fCellMask(), fBadCellMask(), fCurrentClusterDigits(), fFoundClusters(outputArray), fNumFoundClusters(0), fMaxNumClusters(maxNumClusters), fGeometry(geometry), fTimeCut(timeCut), fTimeMin(timeMin), fTimeMax(timeMax), fGradientCut(gradientCut), fDoEnergyGradientCut(doEnergyGradientCut), fThresholdSeedEnergy(thresholdSeedE), fThresholdCellEnergy(thresholdCellE)
{
  // 
}

///
/// Destructor
//____________________________________________________________________________
AliHLTClusterFinder::~AliHLTClusterFinder()
{
  // 
}

///
/// Recursively search for neighbours
//____________________________________________________________________________
AliHLTCaloRecPointDataStruct* AliHLTClusterFinder::GetClusterFromNeighbours(AliHLTCaloRecPointDataStruct* recPoint, Int_t row, Int_t column)
{
  // Cluster/recpoint is already initialized with seed digit
  Bool_t shared = kFALSE;

  // Recursion 0
  if(!recPoint)
  {
    recPoint = new AliHLTCaloRecPointDataStruct();
    recPoint->fMultiplicity = 0;
    fCurrentClusterDigits[recPoint->fMultiplicity] = fDigitMap[row][column];
    recPoint->fMultiplicity++;
    recPoint->fAmp = fDigitMap[row][column]->fEnergy;
    recPoint->fLeadingCellID  = fDigitMap[row][column]->fID;
    recPoint->fModule = fDigitMap[row][column]->fModule;
    recPoint->fTime= fDigitMap[row][column]->fTime;
  }
  Int_t currentSM = fGeometry->GetGeometryPtr()->GetSuperModuleNumber(fDigitMap[row][column]->fID);
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

    if( (row+rowDiff < 0) || (row+rowDiff >= HLTClusterFinder::kNrows) ) continue;
    if( (column+colDiff < 0) || (column+colDiff >= HLTClusterFinder::kNcolumns) ) continue;

    if(fDigitMap[row+rowDiff][column+colDiff])
      if(!fCellMask[row+rowDiff][column+colDiff])
        if (fDoEnergyGradientCut && not (fDigitMap[row+rowDiff][column+colDiff]->fEnergy>fDigitMap[row][column]->fEnergy+fGradientCut))
          if (not (TMath::Abs(fDigitMap[row+rowDiff][column+colDiff]->fTime - fDigitMap[row][column]->fTime) > fTimeCut))
          {
            // Check if cluster extends over two SMs (= shared cluster)
            shared = (currentSM != fGeometry->GetGeometryPtr()->GetSuperModuleNumber(fDigitMap[row+rowDiff][column+colDiff]->fID));
            (void) GetClusterFromNeighbours(recPoint, row+rowDiff, column+colDiff);
            // Add the digit to the current cluster -- if we end up here, the selected cluster fulfills the condition
            if(recPoint->fMultiplicity < 500) fCurrentClusterDigits[recPoint->fMultiplicity] = fDigitMap[row+rowDiff][column+colDiff];
            recPoint->fMultiplicity++;
            recPoint->fAmp += fDigitMap[row+rowDiff][column+colDiff]->fEnergy;
            // Mask the cell as already used for clustering
          }
  }
  return recPoint;
}

///
/// Add bad cells for supermodule to bad cell mask
///
//____________________________________________________________________________
void AliHLTClusterFinder::AddBadCellsForSM(Int_t sm, TH2* hBadCells)
{
  // Loop over bad cell histogram for supermodule sm and add cells to bad cell mask
  for(Int_t i=1; i<=hBadCells->GetNbinsX(); i++)
    for(Int_t j=1; j<=hBadCells->GetNbinsY(); j++)
    {
      if (hBadCells->GetBinContent(i,j) > 0)
      {
        Int_t row = 0;
        Int_t col = 0;
        Int_t absID = fGeometry->GetGeometryPtr()->GetAbsCellIdFromCellIndexes(sm, j-1, i-1);
        GetTopologicalRowColumn(absID, row, col);
        fBadCellMask[row][col] = kTRUE;
      }
    }

}

///
/// Get row (phi) and column (eta) of a digit, values corresponding to topology
///
//____________________________________________________________________________
void AliHLTClusterFinder::GetTopologicalRowColumn(AliHLTCaloDigitDataStruct* digit, Int_t& row, Int_t& column)
{
  GetTopologicalRowColumn(digit->fID, row, column);
}

///
/// Get row (phi) and column (eta) for a digit ID, values corresponding to topology
///
//____________________________________________________________________________
void AliHLTClusterFinder::GetTopologicalRowColumn(UInt_t digitID, Int_t& row, Int_t& column)
{
  Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  row=0;
  column=0;
  // Get SM number and relative row/column for SM
  fGeometry->GetGeometryPtr()->GetCellIndex(digitID, nSupMod,nModule,nIphi,nIeta);
  fGeometry->GetGeometryPtr()->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, row,column);

  // Add shifts wrt. supermodule and type of calorimeter
  // NOTE:
  // * Rows (phi) are arranged that one space is left empty between supermodules in phi
  //   This is due to the physical gap that forbids clustering
  // * For DCAL, there is an additional empty column between two supermodules in eta
  //   Again, this is to account for the gap in DCAL

  row    += nSupMod/2 * (24+1);
  // In DCAL, leave a gap between two SMs with same phi
  if(!fGeometry->GetGeometryPtr()->IsDCALSM(nSupMod)) // EMCAL
    column += nSupMod%2 * 48;
  else
    column += nSupMod%2 * (48+1);
}

///
/// Return number of found clusters. Start clustering from highest energy cell.
//____________________________________________________________________________
Int_t AliHLTClusterFinder::FindClusters(AliHLTCaloDigitDataStruct** digitArray, Int_t numDigits)
{  
  // 
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
  std::memset(fCellMask, 0, sizeof(Bool_t) * HLTClusterFinder::kNrows*HLTClusterFinder::kNcolumns);
  std::memset(fDigitMap, 0, sizeof(AliHLTCaloDigitDataStruct*) * HLTClusterFinder::kNrows*HLTClusterFinder::kNcolumns);

  // Calibrate digits and fill the maps/arrays
  Int_t nCells = 0;
  Double_t ehs = 0.0;
  for (Int_t iDigit = 0; iDigit<numDigits; iDigit++)
  {
    AliHLTCaloDigitDataStruct* digit = digitArray[iDigit];

    Float_t dEnergyCalibrated = digit->fEnergy;
    Float_t time              = digit->fTime;

    if (dEnergyCalibrated < fThresholdCellEnergy || time > fTimeMax || time < fTimeMin) 
      continue;

    ehs += dEnergyCalibrated;

    // Put digit to 2D map
    Int_t row = 0, column = 0;
    GetTopologicalRowColumn(digit, row, column);
    if(row >= HLTClusterFinder::kNrows) {
      HLTError("Row exceed range: %d\n", row);
      continue;
    }
    if(column >= HLTClusterFinder::kNcolumns) {
      HLTError("Column exceed range: %d\n", column);
      continue;
    }
    if(nCells >= HLTClusterFinder::kNcolumns * HLTClusterFinder::kNrows) {
      HLTError("Too may seeds\n");
      break;
    }

    // Do not use digit if it is masked as bad cell
    if(fBadCellMask[row][column])
      continue;

    fDigitMap[row][column] = digit;
    // NOTE: Here we misuse fAssociatedCluster to hold the original digit index in digitArray
    //       This index is needed for the rec point member fDigits
    fDigitMap[row][column]->fAssociatedCluster = iDigit;
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
    if(fNumFoundClusters>=fMaxNumClusters)
    {
      HLTWarning("Buffer exceeded after %d clusters", fNumFoundClusters);
      break;
    }
    Int_t row = fSeedList[i].row, column = fSeedList[i].column;
    // Continue if the cell is already masked (i.e. was already clustered)
    if(fCellMask[row][column])
      continue;
    // Continue if energy constraints are not fulfilled
    if (fSeedList[i].energy<=fThresholdSeedEnergy)
      continue;
    // Seed is found, form cluster recursively
    AliHLTCaloRecPointDataStruct* recPoint = GetClusterFromNeighbours(0, row, column);
    CalculateCenterOfGravity(recPoint);
    fFoundClusters[fNumFoundClusters++] = recPoint;
  }

  HLTDebug(Form("%d clusters found from %d digits (total=%d)-> ehs %f (minE %f)\n",fNumFoundClusters,nCells,numDigits, ehs,fThresholdCellEnergy));
  return fNumFoundClusters;
}

///
/// Calculate cluster center-of-gravity (code taken from AliHLTCaloClusterAnalyser)
//____________________________________________________________________________
void AliHLTClusterFinder::CalculateCenterOfGravity(AliHLTCaloRecPointDataStruct* recPoint)
{
  Float_t wtot = 0.;
  Float_t x = 0.;
  Float_t z = 0.;
  Float_t xi = 0.;
  Float_t zi = 0.;

  Float_t delta_1 = 0.;
  Float_t delta_2 = 0.;
  Float_t delta_3 = 0.;
  Float_t delta_4 = 0.;
  Float_t delta_5 = 0.;

  wtot = 0;
  UInt_t iDigit = 0;
  Double_t logWeight = 4.5; // TODO: Verify this value

  for(UInt_t iDigit = 0; iDigit < TMath::Min(int(recPoint->fMultiplicity), 500); iDigit++)
  {
    AliHLTCaloDigitDataStruct* digit = fCurrentClusterDigits[iDigit];
    xi = digit->fX;
    zi = digit->fZ;
    if (recPoint->fAmp > 0 && digit->fEnergy > 0) 
    {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( digit->fEnergy / recPoint->fAmp));
      x    += xi * w;
      z    += zi * w;
      wtot += w;
      delta_1 += w*xi*xi;
      delta_2 += w*xi;
      delta_3 += w*zi*zi;
      delta_4 += w*zi;
      delta_5 += w*xi*zi;
    }
  }
  if (wtot>0) 
  {
    recPoint->fX = x/wtot;
    recPoint->fZ = z/wtot;

    Double_t sEta = delta_1/wtot - 1/wtot/wtot * delta_2*delta_2;
    Double_t sPhi = delta_3/wtot - 1/wtot/wtot * delta_4*delta_4;
    Double_t sEtaPhi = delta_5/wtot - 1/wtot/wtot * delta_2*delta_4;

    Double_t l0 = (0.5 * (sEta + sPhi) + TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
    Double_t l1 = (0.5 * (sEta + sPhi) - TMath::Sqrt( 0.25 * (sEta - sPhi) * (sEta - sPhi) + sEtaPhi * sEtaPhi ));
    recPoint->fM2x = l0;
    recPoint->fM2z = l1;
  }
  else
  {
    recPoint->fX = -9999;
    recPoint->fZ =-9999;
    recPoint->fM2x = -9999;
    recPoint->fM2z = -9999;
  }
}

ClassImp(AliHLTCaloClusterizer);


///
/// Clusterize event using the AliHLTClusterFinder
//____________________________________________________________________________
Int_t AliHLTCaloClusterizer::ClusterizeEvent(Int_t nDigits)
{
  Int_t nRecPoints = 0;
  
  // Create cluster finder only once, at execution time (cannot be done in constructor, because settings could change afterwards)
  if(!fClusterFinder)
  {
    fClusterFinder = new AliHLTClusterFinder(fRecPointArray, fArraySize, fGeometry, fEmcTimeGate, fCellTimeMin, fCellTimeMax, fGradientCut, fUseGradientCut, fEmcClusteringThreshold, fEmcMinEnergyThreshold);

    // Add bad channel map
    AliCDBManager* cdb = AliCDBManager::Instance();
    AliCDBEntry*   cdbEntry = cdb->Get("EMCAL/Calib/BadChannels");
    if(cdbEntry)
    {
      AliEMCALCalibBadChannels* badChannels = static_cast<AliEMCALCalibBadChannels*>(cdbEntry->GetObject());
      for(Int_t iSM=0; iSM<20; iSM++)
      {
        fClusterFinder->AddBadCellsForSM(iSM, badChannels->GetBadChannelHistForSupermodule(iSM));
      }
    }
  }
  // Perform cluster finding. Output is written to fRecPointArray
  nRecPoints = fClusterFinder->FindClusters(fDigitsPointerArray, nDigits);
  return  nRecPoints;
}


///
/// Constructor
//____________________________________________________________________________
AliHLTCaloClusterizer::AliHLTCaloClusterizer(TString det):
        AliHLTCaloConstantsHandler(det),
        fRecPointArray(0),
        fArraySize(0),
        fEmcClusteringThreshold(0),
        fEmcMinEnergyThreshold(0),
        fEmcTimeGate(0),
        fCellTimeMin(-1.),
        fCellTimeMax(1.),
        fUseGradientCut(kTRUE),
        fGradientCut(0),
        fDigitsPointerArray(0),
        fGeometry(0),
        fClusterFinder(0)
{
    fEmcClusteringThreshold = 0.1;
    fEmcMinEnergyThreshold = 0.01;
    fEmcTimeGate = 1.e-6 ;

    // Reserve memory for buffer 
    fArraySize = 10000; // max num clusters
    fRecPointArray = new AliHLTCaloRecPointDataStruct*[fArraySize];
}

///
/// Destructor
//____________________________________________________________________________
AliHLTCaloClusterizer::~AliHLTCaloClusterizer()
{
  // TODO: Is this already deleted before by analyzer?
  //  for(Int_t i=0; i<fArraySize; i++)
  //    if(fRecPointArray[i]) delete fRecPointArray[i];
  delete [] fRecPointArray;
}
