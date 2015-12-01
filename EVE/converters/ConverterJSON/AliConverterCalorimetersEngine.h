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

    AliESDEvent *fESDEvent;
    
    TLorentzVector fClusterMomentum;  /// Current cluster kinematics in TLorentzVector
    AliESDCaloCluster *fCaloCluster; /// Pointer with current cluster info
    
    AliEMCALGeometry *fGeomEM; /// EMCal geometry
    AliPHOSGeometry  *fGeomPH; /// PHOS geometry
    
    AliESDCaloCells  *fCellsEM; /// List with EMCAL cells
    AliESDCaloCells  *fCellsPH; /// List with PHOS cells
    
    AliConverterCalorimetersEngine(const AliConverterCalorimetersEngine&) {};
    AliConverterCalorimetersEngine& operator=(const AliConverterCalorimetersEngine&) {};
};

#endif
