// Class for reduced calorimeter cluster information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 

#ifndef ALIREDUCEDCALOCLUSTERINFO_H
#define ALIREDUCEDCALOCLUSTERINFO_H

#include <TObject.h>

//_________________________________________________________________________
class AliReducedCaloClusterInfo : public TObject {
  
  friend class AliAnalysisTaskReducedTreeMaker;         // friend analysis task which fills the object
  
 public:
  enum ClusterType {
    kUndefined=0, kEMCAL, kPHOS  
  };
   
  AliReducedCaloClusterInfo();
  virtual ~AliReducedCaloClusterInfo();
  
  Bool_t  IsEMCAL()    const {return (fType==kEMCAL ? kTRUE : kFALSE);}
  Bool_t  IsPHOS()     const {return (fType==kPHOS ? kTRUE : kFALSE);}
  Float_t Energy()     const {return fEnergy;}
  Float_t Dx()         const {return fTrackDx;}
  Float_t Dz()         const {return fTrackDz;}
  Float_t M20()        const {return fM20;}
  Float_t M02()        const {return fM02;}
  Float_t Dispersion() const {return fDispersion;}
  Float_t X()          const {return fPosition[0];}
  Float_t Y()          const {return fPosition[1];}
  Float_t Z()          const {return fPosition[2];}
  Float_t TOF()        const {return fTOF;}
  Short_t NCells()     const {return fNCells;}
  
 protected:
  Char_t  fType;         // cluster type (EMCAL/PHOS)
  Float_t fEnergy;       // cluster energy
  Float_t fTrackDx;      // distance to closest track in phi
  Float_t fTrackDz;      // distance to closest track in z
  Float_t fM20;          // short axis
  Float_t fM02;          // long axis
  Float_t fDispersion;   // dispersion
  Float_t fPosition[3];  // cluster position
  Float_t fTOF;          // time of flight
  Short_t fNCells;       // number of cells
  //---------------------------------------------------
  
  AliReducedCaloClusterInfo(const AliReducedCaloClusterInfo &c);
  AliReducedCaloClusterInfo& operator= (const AliReducedCaloClusterInfo &c);

  ClassDef(AliReducedCaloClusterInfo, 1);
};

#endif