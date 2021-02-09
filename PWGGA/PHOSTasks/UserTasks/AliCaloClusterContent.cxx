#include "AliCaloClusterContent.h"
#include <vector>
#include <Rtypes.h>
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliPHOSGeometry.h"


/*
2017-10-14: added cluster position
*/








class TObject;

ClassImp(AliCaloClusterContent)

//________________________________________________________________________
AliCaloClusterContent::AliCaloClusterContent():
  fType(-1),
  fLabel(-1),
  fNCells(-1),
  fNTracksMatched(-1),
  fIsFilled(-1),
  fIsExotic(-1),
  fIsEMCAL(-1),
  fIsPHOS(-1),
  fPosition(),  //TODO move down and initialize all entries with -1
  fCoreEnergy(-1),
  fDispersion(-1),
  fDistanceToBadChannel(-1),
  fEmcCpvDistance(-1),
  fEnergy(-1),
  fM02(-1),
  fM20(-1),
  fTOF(-1),
  fTrackDx(-1),
  fTrackDz(-1),
  fCellAbsID(0),
  fCellDetector(0),
  fCellMod(0),
  fCellRelIDX(0),
  fCellRelIDZ(0),
  fCellEnergy(0),
  fCellTime(0),
  fCellAmpFrac(0)
{
// allocate memory
  fCellAbsID.reserve(30);
  fCellMod.reserve(30);
  fCellDetector.reserve(30);
  fCellRelIDX.reserve(30);
  fCellRelIDZ.reserve(30);
  fCellEnergy.reserve(30);
  fCellTime.reserve(30);
  fCellAmpFrac.reserve(30);

}

//________________________________________________________________________
AliCaloClusterContent::AliCaloClusterContent(const AliVCluster* clust, AliVCaloCells* cells, const AliPHOSGeometry* fgeom):
  fType(-1),
  fLabel(-1),
  fNCells(-1),
  fNTracksMatched(-1),
  fIsFilled(-1),
  fIsExotic(-1),
  fIsEMCAL(-1),
  fIsPHOS(-1),
  fPosition(),  //TODO move down and initialize all entries with -1
  fCoreEnergy(-1),
  fDispersion(-1),
  fDistanceToBadChannel(-1),
  fEmcCpvDistance(-1),
  fEnergy(-1),
  fM02(-1),
  fM20(-1),
  fTOF(-1),
  fTrackDx(-1),
  fTrackDz(-1),
  fCellAbsID(0),
  fCellDetector(0),
  fCellMod(0),
  fCellRelIDX(0),
  fCellRelIDZ(0),
  fCellEnergy(0),
  fCellTime(0),
  fCellAmpFrac(0)
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
  fCellAmpFrac.reserve(2*NCells);

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
AliCaloClusterContent::AliCaloClusterContent(const AliCaloClusterContent &src ):
  fType(src.fType),
  fLabel(src.fLabel),
  fNCells(src.fNCells),
  fNTracksMatched(src.fNTracksMatched),
  fIsFilled(src.fIsFilled),
  fIsExotic(src.fIsExotic),
  fIsEMCAL(src.fIsEMCAL),
  fIsPHOS(src.fIsPHOS),
  fCoreEnergy(src.fCoreEnergy),
  fDispersion(src.fDispersion),
  fDistanceToBadChannel(src.fDistanceToBadChannel),
  fEmcCpvDistance(src.fEmcCpvDistance),
  fEnergy(src.fEnergy),
  fM02(src.fM02),
  fM20(src.fM20),
  fTOF(src.fTOF),
  fTrackDx(src.fTrackDx),
  fTrackDz(src.fTrackDz),
  fCellAbsID(src.fCellAbsID),
  fCellDetector(src.fCellDetector),
  fCellMod(src.fCellMod),
  fCellRelIDX(src.fCellRelIDX),
  fCellRelIDZ(src.fCellRelIDZ),
  fCellEnergy(src.fCellEnergy),
  fCellTime(src.fCellTime),
  fCellAmpFrac(src.fCellAmpFrac)

{
  for(int i=0; i<3; i++){
    fPosition[i] = src.fPosition[i];
  }
}
//________________________________________________________________________
AliCaloClusterContent& AliCaloClusterContent::operator=(const AliCaloClusterContent &src)
{
  if(&src == this)  return *this;

  fType                 = src.fType;
  fLabel                = src.fLabel;
  fNCells               = src.fNCells;
  fNTracksMatched       = src.fNTracksMatched;
  fIsFilled             = src.fIsFilled;
  fIsExotic             = src.fIsExotic;
  fIsEMCAL              = src.fIsEMCAL;
  fIsPHOS               = src.fIsPHOS;
  fCoreEnergy           = src.fCoreEnergy;
  fDispersion           = src.fDispersion;
  fDistanceToBadChannel = src.fDistanceToBadChannel;
  fEmcCpvDistance       = src.fEmcCpvDistance;
  fEnergy               = src.fEnergy;
  fM02                  = src.fM02;
  fM20                  = src.fM20;
  fTOF                  = src.fTOF;
  fTrackDx              = src.fTrackDx;
  fTrackDz              = src.fTrackDz;
  fCellAbsID            = src.fCellAbsID;
  fCellDetector         = src.fCellDetector;
  fCellMod              = src.fCellMod;
  fCellRelIDX           = src.fCellRelIDX;
  fCellRelIDZ           = src.fCellRelIDZ;
  fCellEnergy           = src.fCellEnergy;
  fCellTime             = src.fCellTime;
  fCellAmpFrac          = src.fCellAmpFrac;


  for(int i=0; i<3; i++){
    fPosition[i] = src.fPosition[i];
  }


  return *this;

}






//________________________________________________________________________
AliCaloClusterContent::~AliCaloClusterContent()
{

}



//________________________________________________________________________
void AliCaloClusterContent::SetClusterAndCells(const AliVCluster* clust, AliVCaloCells* cells, const AliPHOSGeometry* fgeom){
// read out cluster information
  clust ->GetPosition(fPosition);
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


    Double_t dCellEnergy  = cells ->GetCellAmplitude(iCellAbsId);
    Double_t dCellTime    = cells ->GetCellTime(iCellAbsId);
    Double_t dCellAmpFrac = clust ->GetCellAmplitudeFraction(iCellAbsId);

    fCellAbsID.push_back(iCellAbsId);
    fCellMod.push_back(iMod);
    fCellDetector.push_back(iDetector);
    fCellRelIDX.push_back(icellX);
    fCellRelIDZ.push_back(icellZ);
    fCellEnergy.push_back(dCellEnergy);
    fCellTime.push_back(dCellTime);
    fCellAmpFrac.push_back(dCellAmpFrac);

    }

  fIsFilled = kTRUE;
}

//________________________________________________________________________
//void AliCaloClusterContent::SetMatchedTrack(Double_t Distance, Double_t DistanceX, Double DistanceZ){
//
//
//  fDistanceToTrack  = Distance;
//  fDistanceXToTrack = DistanceX;
//  fDistanceZToTrack = DistanceZ;
//
//
//
//  //fIsTrackMatched = kTRUE;
//}

//________________________________________________________________________
void AliCaloClusterContent::Reset(){

  fType                 = -1;
  fLabel                = -1;
  fNCells               = -1;
  fNTracksMatched       = -1;
  fIsExotic             = -1;
  fIsEMCAL              = -1;
  fIsPHOS               = -1;
  fCoreEnergy           = -1;
  fDispersion           = -1;
  fDistanceToBadChannel = -1;
  fEmcCpvDistance       = -1;
  fEnergy               = -1;
  fM02                  = -1;
  fM20                  = -1;
  fTOF                  = -1;
  fTrackDx              = -1;
  fTrackDz              = -1;

  fCellAbsID.clear();
  fCellMod.clear();
  fCellDetector.clear();
  fCellRelIDX.clear();
  fCellRelIDZ.clear();
  fCellEnergy.clear();
  fCellTime.clear();
  fCellAmpFrac.clear();


  fIsFilled = kFALSE;
}





