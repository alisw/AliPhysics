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
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TGeoManager.h>

// --- Standard library --

// --- AliRoot header files ---
#include <AliEMCALAfterBurnerUF.h>
#include <AliEMCALGeometry.h>
#include <AliEMCALUnfolding.h>
#include <AliAODCaloCluster.h>
#include <AliVCaloCells.h>
#include <AliEMCALRecPoint.h>
#include <AliEMCALDigit.h>

ClassImp(AliEMCALAfterBurnerUF)

//------------------------------------------------------------------------
AliEMCALAfterBurnerUF::AliEMCALAfterBurnerUF():
  fGeom(NULL),
  fLogWeight(4.5),      // correct?
  fECALocMaxCut(0.03),  // value suggested by Adam
  fMinECut(0.01),
  fRecPoints(NULL),
  fDigitsArr(NULL),
  fClusterUnfolding(NULL)
{
  // Use this constructor, if unsure

  Init();
}

//------------------------------------------------------------------------
AliEMCALAfterBurnerUF::AliEMCALAfterBurnerUF(Float_t logWeight, Float_t ecaLocMaxCut, Float_t minECut):
  fGeom(NULL),
  fLogWeight(logWeight),
  fECALocMaxCut(ecaLocMaxCut),
  fMinECut(minECut),
  fRecPoints(NULL),
  fDigitsArr(NULL),
  fClusterUnfolding(NULL)
{
  // This constructor allows to set parameters
  // Recommended values:
  //   Float_t logWeight = 4.5, ECALocMaxCut = 0.03

  Init();
}

//------------------------------------------------------------------------
void AliEMCALAfterBurnerUF::Init()
{
  // After-burner initialization
  // Imports geometry.root (if required), creates unfolding class instance
  //
  // TODO: geometry.root does not allow to use the method AliEMCALRecPoint::EvalAll();
  //       for this to work, the OCDB geometry must be imported instead

  if (!gGeoManager)
    Warning("AliEMCALAfterBurnerUF::Init","GeoManager not found, please import the geometry.root file or pass to the geometry the misalignment matrices");
//    TGeoManager::Import("geometry.root");

  // required for global cluster position recalculation
  if (!gGeoManager)
    Info("AliEMCALAfterBurnerUF::Init", "gGeoManager was not set, be careful");

  // initialize geometry, if not yet initialized
  if (!AliEMCALGeometry::GetInstance()) {
    Warning("AliEMCALAfterBurnerUF::Init", "AliEMCALGeometry is not yet initialized. Initializing with EMCAL_COMPLETE12SMV1_DCAL_8SM");
    AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  }

  // AliEMCALRecPoint is using exactly this call
  fGeom = AliEMCALGeometry::GetInstance();
  if (!fGeom)
    Fatal("AliEMCALAfterBurnerUF::AliEMCALAfterBurnerUF", "could not get EMCAL geometry");

  fClusterUnfolding = new AliEMCALUnfolding(fGeom);
  fClusterUnfolding->SetECALocalMaxCut(fECALocMaxCut);
  fClusterUnfolding->SetThreshold(fMinECut); 
  
  // clusters --> recPoints, cells --> digits and back ;)
  fRecPoints = new TObjArray(100);
  fDigitsArr = new TClonesArray("AliEMCALDigit",1152);
}

//------------------------------------------------------------------------
AliEMCALAfterBurnerUF::~AliEMCALAfterBurnerUF()
{
  //
  // destructor
  //

  if (fClusterUnfolding) delete fClusterUnfolding;

  if (fRecPoints) {
    fRecPoints->Delete();
    delete fRecPoints;
  }
  if (fDigitsArr) {
    fDigitsArr->Clear("C");
    delete fDigitsArr;
  }
}

//------------------------------------------------------------------------
void AliEMCALAfterBurnerUF::Clear()
{
  //Clean the arrays
  
  if (fRecPoints) fRecPoints->Delete();  // do not Clear(), it leaks, why?
  if (fDigitsArr) fDigitsArr->Clear("C");
  
}
//------------------------------------------------------------------------
void AliEMCALAfterBurnerUF::RecPoints2Clusters(TObjArray *clusArray)
{
  // Restore clusters from recPoints
  // Cluster energy, global position, cells and their amplitude fractions are restored
  //
  // TODO: restore time and other parameters

  for(Int_t i = 0; i < fRecPoints->GetEntriesFast(); i++)
  {
    AliEMCALRecPoint *recPoint = (AliEMCALRecPoint*) fRecPoints->At(i);

    const Int_t ncells = recPoint->GetMultiplicity();
    Int_t ncellsTrue = 0;

    // cells and their amplitude fractions
    UShort_t absIds[ncells];  // NOTE: unfolding must not give recPoints with no cells!
    Double32_t ratios[ncells];

    for (Int_t c = 0; c < ncells; c++) {
      AliEMCALDigit *digit = (AliEMCALDigit*) fDigitsArr->At(recPoint->GetDigitsList()[c]);

      absIds[ncellsTrue] = digit->GetId();
      ratios[ncellsTrue] = recPoint->GetEnergiesList()[c]/digit->GetAmplitude();

      if (ratios[ncellsTrue] > 0.001) ncellsTrue++;
    }

    if (ncellsTrue < 1) {
      Warning("AliEMCALAfterBurnerUF::RecPoints2Clusters", "skipping cluster with no cells");
      continue;
    }

    TVector3 gpos;
    Float_t g[3];

    // calculate new cluster position
    recPoint->EvalGlobalPosition(fLogWeight, fDigitsArr);
    recPoint->GetGlobalPosition(gpos);
    gpos.GetXYZ(g);

    // create a new cluster
    AliAODCaloCluster *clus = new AliAODCaloCluster();
    clus->SetType(AliVCluster::kEMCALClusterv1);
    clus->SetE(recPoint->GetEnergy());
    clus->SetPosition(g);
    clus->SetNCells(ncellsTrue);
    clus->SetCellsAbsId(absIds);
    clus->SetCellsAmplitudeFraction(ratios);
    // TODO: time not stored
    // TODO: some other properties not stored

    clusArray->Add(clus);
  } // recPoints loop

}

//------------------------------------------------------------------------
void AliEMCALAfterBurnerUF::UnfoldClusters(TObjArray *clusArray, AliVCaloCells *cellsEMCAL)
{
  // Unfolds clusters.
  //
  // Input: TObjArray of clusters, EMCAL cells.
  // Output is added to the same array, original clusters are _deleted_ or moved to another position.

  Int_t ndigits = 0;

  Int_t nclus = clusArray->GetEntriesFast();

  /* Fill recPoints with digits
  */
  for (Int_t i = 0; i < nclus; i++)
  {
    AliVCluster *clus = (AliVCluster*) clusArray->At(i);
    if (!clus->IsEMCAL()) continue;

    // new recPoint
    AliEMCALRecPoint *recPoint = new AliEMCALRecPoint("");
    recPoint->SetClusterType(AliVCluster::kEMCALClusterv1);
    fRecPoints->Add(recPoint);

    // fill digits
    for (Int_t c = 0; c < clus->GetNCells(); c++) {
      Int_t absId = clus->GetCellAbsId(c);
      Double_t amp = cellsEMCAL->GetCellAmplitude(absId);
      Double_t time = cellsEMCAL->GetCellTime(absId);

      // NOTE: it is easy to include cells recalibration here:
      // amp *= factor;

      AliEMCALDigit *digit = (AliEMCALDigit*) fDigitsArr->New(ndigits);

      digit->SetId(absId);
      digit->SetAmplitude(amp);
      digit->SetTime(time);
      digit->SetTimeR(time);
      digit->SetIndexInList(ndigits);

      recPoint->AddDigit(*digit, amp, kFALSE);

      ndigits++;
    }

    // this cluster will be substituted with the result of unfolding
    clusArray->RemoveAt(i);
    delete clus;
  } // cluster loop


  /* Peform unfolding
  */
  fClusterUnfolding->SetInput(fRecPoints->GetEntriesFast(), fRecPoints, fDigitsArr);
  fClusterUnfolding->MakeUnfolding();

  /* Restore clusters from recPoints.
  */
  RecPoints2Clusters(clusArray);

  // clean up
  fRecPoints->Delete(); // do not Clear(), it leaks, why?
  fDigitsArr->Clear("C");

  clusArray->Compress();

}
