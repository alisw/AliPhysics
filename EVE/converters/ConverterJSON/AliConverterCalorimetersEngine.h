/// \class AliConverterCalorimetersEngine
/// AliConverterCalorimetersEngine is responsible for
/// extracting data which refers to Calorimeters
/// Usage example:
///     AliMinimalisticEvent &event;
///     AliConverterCalorimetersEngine *caloEngine = new AliConverterCalorimetersEngine(fESDEvent); /// fESDEvent is a member of ExternalFormatConverter
///     caloEngine->PopulateEventWithCaloClusters(event);
//
/// \author Maciej Grochowicz <maciej.aleksander.grochowicz@cern.ch>, Warsaw University of Technology


#ifndef ALIROOT_ALICONVERTERCALORIMETERSENGINE_H
#define ALIROOT_ALICONVERTERCALORIMETERSENGINE_H

#include <AliMinimalisticEvent.h>

#include <AliESDEvent.h>
#include <TEveVector.h>

#include <AliEMCALGeometry.h>
#include <AliPHOSGeometry.h>

#include <AliESDCaloCells.h>
#include <AliESDCaloCluster.h>

class AliConverterCalorimetersEngine
{
public:
    AliConverterCalorimetersEngine(AliESDEvent *event);
    ~AliConverterCalorimetersEngine();

    void PopulateEventWithCaloClusters(AliMinimalisticEvent &event);
    
private:
    void AssertGeometry();
    void AssertMagField();

    float GetECross(bool isEMCAL, int imod, int icol, int irow);
    void GetModuleNumberColAndRow(int absId, bool isEMCAL,int &iSM, int &icol, int &irow);
    void GetMaxEnergyCellAbsId(int &absId, float &eMax);
    bool IsBadCluster(int absId, float eMax);
    void SetUpEMCALGeometry();
    void SetUpPHOSGeometry();
    double GetPhi(double phi);

    void AddEMCALClustersToEvent(AliMinimalisticEvent &event, Int_t id, Float_t amp, Float_t phi, Float_t eta);
    void AddPHOSCalClusterToEvent(AliMinimalisticEvent &event, Int_t &id, Float_t &amp);

    AliESDEvent *fESDEvent;

    TLorentzVector fClusterMomentum;  /// Current cluster kinematics in TLorentzVector
    AliESDCaloCluster *fCaloCluster; /// Pointer with current cluster info

    AliEMCALGeometry *fGeomEM; /// EMCal geometry
    AliPHOSGeometry  *fGeomPH; /// PHOS geometry

    AliESDCaloCells  *fCellsEM; /// List with EMCAL cells
    AliESDCaloCells  *fCellsPH; /// List with PHOS cells

    /// Copy ctor -- it is not recommended to copy CalorimetersEngine therefore it is private
    AliConverterCalorimetersEngine(const AliConverterCalorimetersEngine&) {};
    /// Assignment operator -- it is not recommended to assign CalorimetersEngine therefore it is private
    AliConverterCalorimetersEngine& operator=(const AliConverterCalorimetersEngine&) {};


};

#endif
