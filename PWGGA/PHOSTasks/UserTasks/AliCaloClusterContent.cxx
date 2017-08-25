#include "AliCaloClusterContent.h"
#include <vector>
#include <Rtypes.h>
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliPHOSGeometry.h"

class TObject;

ClassImp(AliCaloClusterContent)

//________________________________________________________________________
AliCaloClusterContent::AliCaloClusterContent():
  fType(0),
  fLabel(0),
  fNCells(0),
  fNTracksMatched(0),
  fIsFilled(0),
  fIsExotic(0),
  fIsEMCAL(0),
  fIsPHOS(0),
  fCoreEnergy(0),
  fDispersion(0),
  fDistanceToBadChannel(0),
  fEmcCpvDistance(0),
  fEnergy(0),
  fM02(0),
  fM20(0),
  fTOF(0),
  fTrackDx(0),
  fTrackDz(0),
  fCellAbsID(0),
  fCellDetector(0),
  fCellMod(0),
  fCellRelIDX(0),
  fCellRelIDZ(0),
  fCellEnergy(0),
  fCellTime(0)
{
// allocate memory
  fCellAbsID.reserve(30);
  fCellMod.reserve(30);
  fCellDetector.reserve(30);
  fCellRelIDX.reserve(30);
  fCellRelIDZ.reserve(30);
  fCellEnergy.reserve(30);
  fCellTime.reserve(30);
}

//________________________________________________________________________
AliCaloClusterContent::AliCaloClusterContent(const AliVCluster* clust, AliVCaloCells* cells, const AliPHOSGeometry* fgeom):
  fType(0),
  fLabel(0),
  fNCells(0),
  fNTracksMatched(0),
  fIsFilled(0),
  fIsExotic(0),
  fIsEMCAL(0),
  fIsPHOS(0),
  fCoreEnergy(0),
  fDispersion(0),
  fDistanceToBadChannel(0),
  fEmcCpvDistance(0),
  fEnergy(0),
  fM02(0),
  fM20(0),
  fTOF(0),
  fTrackDx(0),
  fTrackDz(0),
  fCellAbsID(0),
  fCellDetector(0),
  fCellMod(0),
  fCellRelIDX(0),
  fCellRelIDZ(0),
  fCellEnergy(0),
  fCellTime(0)
{

  Int_t NCells = clust ->GetNCells();
// allocate memory
  fCellAbsID.reserve(2*NCells);
  fCellMod.reserve(2*NCells);
  fCellDetector.reserve(2*NCells);
  fCellRelIDX.reserve(2*NCells);
  fCellRelIDZ.reserve(2*NCells);
  fCellEnergy.reserve(2*NCells);
  fCellTime.reserve(2*NCells);

  this->SetClusterAndCells(clust, cells, fgeom);

/*
// read out cluster information
  fType                 = clust ->GetType();
  fLabel                = clust ->GetLabel();
  fNCells               = clust ->GetNCells();
  fNTracksMatched       = clust ->GetNTracksMatched();
  fIsExotic             = clust ->GetIsExotic();
  fIsEMCAL              = clust ->IsEMCAL();
  fIsPHOS               = clust ->IsPHOS();
  fCoreEnergy           = clust ->GetCoreEnergy();
  fDispersion           = clust ->GetDispersion();
  fDistanceToBadChannel = clust ->GetDistanceToBadChannel();
  fEmcCpvDistance       = clust ->GetEmcCpvDistance();
  fEnergy               = clust ->E();
  fM02                  = clust ->GetM02();
  fM20                  = clust ->GetM20();
  fTOF                  = clust ->GetTOF();
  fTrackDx              = clust ->GetTrackDx();
  fTrackDz              = clust ->GetTrackDz();

// allocate memory
  fCellAbsID.reserve(2*fNCells);
  fCellMod.reserve(2*fNCells);
  fCellDetector.reserve(2*fNCells);
  fCellRelIDX.reserve(2*fNCells);
  fCellRelIDZ.reserve(2*fNCells);
  fCellEnergy.reserve(2*fNCells);
  fCellTime.reserve(2*fNCells);




// Loop over all cells of current cluster
  for(int icell=0; icell<fNCells; icell++){
    Int_t iCellAbsId = clust->GetCellAbsId(icell);
    Int_t irelIDCell[4] = {0};
    fgeom ->AbsToRelNumbering(iCellAbsId, irelIDCell);

    Int_t iMod, iDetector, icellX, icellZ;
    iMod      = irelIDCell[0];
    iDetector = irelIDCell[1];
    icellX    = irelIDCell[2];
    icellZ    = irelIDCell[3];

    Double_t dCellEnergy = cells ->GetCellAmplitude(iCellAbsId);
    Double_t dCellTime   = cells ->GetCellTime(iCellAbsId);

    fCellAbsID.push_back(iCellAbsId);
    fCellMod.push_back(iMod);
    fCellDetector.push_back(iDetector);
    fCellRelIDX.push_back(icellX);
    fCellRelIDZ.push_back(icellZ);
    fCellEnergy.push_back(dCellEnergy);
    fCellTime.push_back(dCellTime);
    }

  fIsFilled = kTRUE;
  */
}

//________________________________________________________________________
AliCaloClusterContent::~AliCaloClusterContent()
{

}



//________________________________________________________________________
void AliCaloClusterContent::SetClusterAndCells(const AliVCluster* clust, AliVCaloCells* cells, const AliPHOSGeometry* fgeom){
// read out cluster information
  fType                 = clust ->GetType();
  fLabel                = clust ->GetLabel();
  fNCells               = clust ->GetNCells();
  fNTracksMatched       = clust ->GetNTracksMatched();
  fIsExotic             = clust ->GetIsExotic();
  fIsEMCAL              = clust ->IsEMCAL();
  fIsPHOS               = clust ->IsPHOS();
  fCoreEnergy           = clust ->GetCoreEnergy();
  fDispersion           = clust ->GetDispersion();
  fDistanceToBadChannel = clust ->GetDistanceToBadChannel();
  fEmcCpvDistance       = clust ->GetEmcCpvDistance();
  fEnergy               = clust ->E();
  fM02                  = clust ->GetM02();
  fM20                  = clust ->GetM20();
  fTOF                  = clust ->GetTOF();
  fTrackDx              = clust ->GetTrackDx();
  fTrackDz              = clust ->GetTrackDz();

// Loop over all cells of current cluster
  for(int icell=0; icell<fNCells; icell++){
    Int_t iCellAbsId = clust->GetCellAbsId(icell);
    Int_t irelIDCell[4] = {0};
    fgeom ->AbsToRelNumbering(iCellAbsId, irelIDCell);

    Int_t iMod, iDetector, icellX, icellZ;
    iMod      = irelIDCell[0];
    iDetector = irelIDCell[1];
    icellX    = irelIDCell[2];
    icellZ    = irelIDCell[3];


    Double_t dCellEnergy = cells ->GetCellAmplitude(iCellAbsId);
    Double_t dCellTime   = cells ->GetCellTime(iCellAbsId);

    fCellAbsID.push_back(iCellAbsId);
    fCellMod.push_back(iMod);
    fCellDetector.push_back(iDetector);
    fCellRelIDX.push_back(icellX);
    fCellRelIDZ.push_back(icellZ);
    fCellEnergy.push_back(dCellEnergy);
    fCellTime.push_back(dCellTime);
    }

  fIsFilled = kTRUE;
}






//________________________________________________________________________
void AliCaloClusterContent::Reset(){

  fType                 = 0;
  fLabel                = 0;
  fNCells               = 0;
  fNTracksMatched       = 0;
  fIsExotic             = 0;
  fIsEMCAL              = 0;
  fIsPHOS               = 0;
  fCoreEnergy           = 0;
  fDispersion           = 0;
  fDistanceToBadChannel = 0;
  fEmcCpvDistance       = 0;
  fEnergy               = 0;
  fM02                  = 0;
  fM20                  = 0;
  fTOF                  = 0;
  fTrackDx              = 0;
  fTrackDz              = 0;

  fCellAbsID.clear();
  fCellMod.clear();
  fCellDetector.clear();
  fCellRelIDX.clear();
  fCellRelIDZ.clear();
  fCellEnergy.clear();
  fCellTime.clear();


  fIsFilled = kFALSE;
}





