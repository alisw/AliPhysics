#ifndef CEPRawCaloClusterTrack_H
#define CEPRawCaloClusterTrack_H

#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"

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
    /// Global position of the cluster 
    Double_t        fGlobalPos[3];   /// position in global coordinate system (cm)
    /// Moments along eigenaxis
    Double_t        fM02;            /// 2nd moment along the main eigen axis
    Double_t        fM20;            /// 2nd moment along the second eigen axis
    /// Time in the higest amplitude cell
    Double_t        fTime;

  public:
                    CEPRawCaloClusterTrack();
    virtual         ~CEPRawCaloClusterTrack() {};
    /// Modifiers
    void            Reset();
 
    /// Setter functions
    void            SetCaloClusterE          (Double_t e)       { fE           = e;      }
    void            SetCaloClusterShapeDisp  (Double_t shDisp)  { fShapeDisp   = shDisp; }
    void            SetCaloClusterChi2       (Double_t chi2)    { fChi2        = chi2;   }
    void            SetCaloClusterCPVDist    (Double_t cpvD)    { fCaloCpvDist = cpvD;   }
    void            SetCaloClusterIsEMCAL    (Bool_t   isEMC)   { fIsEMCAL     = isEMC;  }
    void            SetCaloClusterIsPHOS     (Bool_t   isPHOS)  { fIsPHOS      = isPHOS; }

    void            SetCaloClusterM20        (Double_t m20)     { fM20         = m20;    }
    void            SetCaloClusterM02        (Double_t m02)     { fM02         = m02;    }
    void            SetCaloClusterTime       (Double_t tme)     { fTime        = tme;    }

    void            SetCaloClusterGlobalPosition(Float_t *x);


    /// Global Setter
    void            SetCaloClusterVariables(AliESDCaloCluster* ClusterObj, AliESDCaloCells* CaloCells);

    /// Accessors
    Float_t         GetCaloClusterE()                 const { return fE;           }
    Float_t         GetCaloClusterShapeDispersion()   const { return fShapeDisp;   }
    Float_t         GetCaloClusterChi2()              const { return fChi2;        }
    Float_t         GetCaloClusterCPVDist()           const { return fCaloCpvDist; }
    Float_t         GetCaloClusterIsEMCAL()           const { return fIsEMCAL;     }
    Float_t         GetCaloClusterIsPHOS()            const { return fIsPHOS;      }

    Float_t         GetCaloClusterM20()   const  { return fM20; }
    Float_t         GetCaloClusterM02()   const  { return fM02; }
    Float_t         GetCaloClusterTime()  const  { return fTime; }
    
    void            GetCaloClusterGlobalPosition(Float_t &x, Float_t &y, Float_t &z);


    ClassDef(CEPRawCaloClusterTrack,1)
};

#endif
