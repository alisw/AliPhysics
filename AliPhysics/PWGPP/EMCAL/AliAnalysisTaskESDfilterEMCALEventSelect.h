#ifndef ALIANALYSISTASKESDFILTEREMCALEVENTSELECT_H
#define ALIANALYSISTASKESDFILTEREMCALEVENTSELECT_H

//______________________________________________________
/// \class AliAnalysisTaskESDfilterEMCALEventSelect
/// \ingroup EMCALPerformance 
/// \brief Filter ESDs events into AODs with some significant calorimeter signal
///
/// Filter the ESD Events to AODs, only those events with
/// some signal in EMCAL, righ now at least a
/// cluster of high energy.
/// Class derived from AliAnalysisTaskESDfilter.
///
/// Something similar is done in AliAnalysisTaskCaloFilter but more complete (?).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//______________________________________________________

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
    
  Float_t             fEnergyCut;       ///<  At least a cluster with this energy in the event.
    
  Int_t               fNcellsCut;       ///<  At least a cluster with fNCellsCut cells over fEnergyCut.
    
  AliEMCALRecoUtils * fRecoUtils;       ///<  AliEMCALRecoUtils pointer.
    
  AliEMCALGeometry  * fGeometry;        ///<  Access to EMCAL geometry utils.
    
  TString             fGeoName;         ///<  Name of geometry used.
    
  /// Copy constructor not implemented.
  AliAnalysisTaskESDfilterEMCALEventSelect(           const AliAnalysisTaskESDfilterEMCALEventSelect&) ;
    
  /// Assignment operator not implemented.
  AliAnalysisTaskESDfilterEMCALEventSelect& operator=(const AliAnalysisTaskESDfilterEMCALEventSelect&) ;
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskESDfilterEMCALEventSelect, 1) ;
  /// \endcond

};
#endif 
