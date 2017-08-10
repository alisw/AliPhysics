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
  fCellsAbsId(0),
//  frelID(),
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

}

//________________________________________________________________________
AliCaloClusterContent::AliCaloClusterContent(const AliVCluster* clust, AliVCaloCells* cells, const AliPHOSGeometry* fgeom):
  fType(0),
  fCellsAbsId(0),
//  relID(),
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
  tempVecSize(0),
  fCellTime(0)
{
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
 /// fCellType.clear();

  tempVecSize = 0;



  fIsFilled = kFALSE;
}





