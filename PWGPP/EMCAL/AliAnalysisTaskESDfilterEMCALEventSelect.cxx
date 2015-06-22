/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliESDCaloCluster.h"

#include "AliAnalysisTaskESDfilterEMCALEventSelect.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskESDfilterEMCALEventSelect) ;
/// \endcond

///
/// Default constructor.
///
//____________________________________________________________________________________
AliAnalysisTaskESDfilterEMCALEventSelect::AliAnalysisTaskESDfilterEMCALEventSelect() :
AliAnalysisTaskESDfilter("ESD Filte : EMCAL selected events"), 
fEnergyCut(10),            fNcellsCut (2),
fRecoUtils(0x0),           
fGeometry(0),              fGeoName("EMCAL_COMPLETE12SMV1_DCAL_8SM")        
{
  fRecoUtils = new AliEMCALRecoUtils;
}		      

///
/// Constructor.
///
//____________________________________________________________________________________________________
AliAnalysisTaskESDfilterEMCALEventSelect::AliAnalysisTaskESDfilterEMCALEventSelect(const char *name) :
AliAnalysisTaskESDfilter(name), 
fEnergyCut(10),            fNcellsCut (2),
fRecoUtils(0x0),           
fGeometry(0),              fGeoName("EMCAL_COMPLETE12SMV1_DCAL_8SM")
{
  fRecoUtils = new AliEMCALRecoUtils;
}

///
/// \return True if there is signal in EMCal
/// Accept event given there is a EMCAL cluster with
/// enough energy and cells.
///
//_________________________________________________________________
Bool_t AliAnalysisTaskESDfilterEMCALEventSelect::AcceptEventEMCAL()
{
  if(!fGeometry)  fGeometry  = AliEMCALGeometry::GetInstance(fGeoName);

  Int_t           nCluster = InputEvent() -> GetNumberOfCaloClusters();
  AliVCaloCells * caloCell = InputEvent() -> GetEMCALCells();
  Int_t           bc       = InputEvent() -> GetBunchCrossNumber();
  
  for(Int_t icalo = 0; icalo < nCluster; icalo++)
  {
    AliESDCaloCluster *clus = (AliESDCaloCluster*) (InputEvent()->GetCaloCluster(icalo));
    
    if( ( clus->IsEMCAL() ) && ( clus->GetNCells() > fNcellsCut ) && ( clus->E() > fEnergyCut ) &&
        fRecoUtils->IsGoodCluster(clus,fGeometry,caloCell,bc))
    {
      
     // printf("Accept event %d, E %2.2f > %2.2f, nCells %d > %d \n", 
     //        (Int_t) Entry(),clus->E(), fEnergyCut, clus->GetNCells(), fNcellsCut);
      
      return kTRUE;
    }
    
  }// loop
  
  return kFALSE;
}

///
/// Main method, execute per event:
/// * Check if the events contains what we want in EMCAL, if not, do not copy the ESD into AOD
/// * Continue the processing in the same way as in the ESD filter.
///
//_________________________________________________________________
void AliAnalysisTaskESDfilterEMCALEventSelect::UserExec(Option_t *)
{
  if(!AcceptEventEMCAL()) return ;
    
  AliAnalysisTaskESDfilter::UserExec("");
}
