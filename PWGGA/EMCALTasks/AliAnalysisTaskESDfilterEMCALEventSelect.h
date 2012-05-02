#ifndef ALIANALYSISTASKESDFILTEREMCALEVENTSELECT_H
#define ALIANALYSISTASKESDFILTEREMCALEVENTSELECT_H

// $Id$

#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliAnalysisTaskESDfilter.h"

class AliAnalysisTaskESDfilterEMCALEventSelect : public AliAnalysisTaskESDfilter 
{
public:
  AliAnalysisTaskESDfilterEMCALEventSelect();                   // default constructor
  AliAnalysisTaskESDfilterEMCALEventSelect(const char *name);   // named constructor
  virtual ~AliAnalysisTaskESDfilterEMCALEventSelect() { ; }     // destructor
  
  void    UserExec(Option_t *option);               
  Bool_t  AcceptEventEMCAL() ;
  void    AccessBadMap();
  void    SetGeometryName(TString name)  { fGeoName = name   ; } 
  TString GetGeometryName()        const { return fGeoName   ; } 
  void    SetEnergyCut(Float_t cut)      { fEnergyCut = cut  ; }
  Float_t GetEnergyCut()           const { return fEnergyCut ; }
  void    SetNcellsCut(Int_t cut)        { fNcellsCut = cut  ; }
  Int_t   GetNcellsCut()           const { return fNcellsCut ; }
  
  AliEMCALRecoUtils* GetRecoUtils()      { return fRecoUtils ; }
  
private:
  Float_t             fEnergyCut;       //  At least a cluster with this energy in the event
  Int_t               fNcellsCut;       //  At least a cluster with fNCellsCut cells over fEnergyCut
  AliEMCALRecoUtils * fRecoUtils;       //  RecoUtils
  AliEMCALGeometry  * fGeometry;        //  Access to EMCAL geometry utils
  TString             fGeoName;         //  Name of geometry used
    
  AliAnalysisTaskESDfilterEMCALEventSelect(           const AliAnalysisTaskESDfilterEMCALEventSelect&); // not implemented
  AliAnalysisTaskESDfilterEMCALEventSelect& operator=(const AliAnalysisTaskESDfilterEMCALEventSelect&); // not implemented
  
  ClassDef(AliAnalysisTaskESDfilterEMCALEventSelect, 1);  
};
#endif 
