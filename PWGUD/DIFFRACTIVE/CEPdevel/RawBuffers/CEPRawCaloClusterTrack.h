#ifndef CEPRawCaloClusterTrack_H
#define CEPRawCaloClusterTrack_H

#include "AliESDCaloCluster.h"

class CEPRawCaloClusterTrack : public TObject {

  private:
    /// Cluster variables
    /// E measured by calo [GeV]
    Double_t        fE;
    /// Cluster shape dispersion
    Double_t        fShapeDisp;
    /// Chi2 of cluster fit
    Double_t        fChi2;
    /// Distance from PHOS/EMC rec. point to the closest CPV rec. point
    Double_t        fCaloCpvDist;
    /// check wether the cluster is found in the EMCal or PHOS 
    Bool_t          fIsEMCAL;
    Bool_t          fIsPHOS;


  public:
                    CEPRawCaloClusterTrack();
    virtual         ~CEPRawCaloClusterTrack() {};
    /// Modifiers
    void            Reset();
 
    /// Setter functions
    void            SetCaloClusterE          (Double_t e)              { fE           = e;      }
    void            SetCaloClusterShapeDisp  (Double_t shDisp)         { fShapeDisp   = shDisp; }
    void            SetCaloClusterChi2       (Double_t chi2)           { fChi2        = chi2;   }
    void            SetCaloClusterCPVDist    (Double_t cpvD)           { fCaloCpvDist = cpvD;   }
    void            SetCaloClusterIsEMCAL    (Bool_t   isEMC)          { fIsEMCAL     = isEMC;  }
    void            SetCaloClusterIsPHOS     (Bool_t   isPHOS)         { fIsPHOS      = isPHOS; }


    /// Global Setter
    void            SetCaloClusterVariables(AliESDCaloCluster* ClusterObj);

    /// Accessors
    Float_t         GetCaloClusterE()                 const { return fE;           }
    Float_t         GetCaloClusterShapeDispersion()   const { return fShapeDisp;   }
    Float_t         GetCaloClusterChi2()              const { return fChi2;        }
    Float_t         GetCaloClusterCPVDist()           const { return fCaloCpvDist; }
    Float_t         GetCaloClusterIsEMCAL()           const { return fIsEMCAL;     }
    Float_t         GetCaloClusterIsPHOS()            const { return fIsPHOS;      }


    ClassDef(CEPRawCaloClusterTrack,1);
};

#endif
