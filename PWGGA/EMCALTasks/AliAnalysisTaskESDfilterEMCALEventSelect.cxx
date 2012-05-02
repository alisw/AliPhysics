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

//////////////////////////////////////////////////////////
// Calls derived from AliAnalysisTaskESDfilter
// Filter the ESD Events to AODs, only those events with
// some signal in EMCAL, righ now at least a 
// cluster of high energy
//
// Author: Gustavo Conesa Balbastre (LPSC - Grenoble)
//
// $Id$
//
//////////////////////////////////////////////////////////

#include "AliESDCaloCluster.h"

#include "AliAnalysisTaskESDfilterEMCALEventSelect.h"

ClassImp(AliAnalysisTaskESDfilterEMCALEventSelect)

//____________________________________________________________________________________
AliAnalysisTaskESDfilterEMCALEventSelect::AliAnalysisTaskESDfilterEMCALEventSelect() : 
AliAnalysisTaskESDfilter("ESD Filte : EMCAL selected events"), 
fEnergyCut(10),            fNcellsCut (2),
fRecoUtils(0x0),           
fGeometry(0),              fGeoName("EMCAL_COMPLETE12SMV1")        
{
  // Default constructor
  
  fRecoUtils = new AliEMCALRecoUtils;

}		      

//____________________________________________________________________________________________________
AliAnalysisTaskESDfilterEMCALEventSelect::AliAnalysisTaskESDfilterEMCALEventSelect(const char *name) : 
AliAnalysisTaskESDfilter(name), 
fEnergyCut(10),            fNcellsCut (2),
fRecoUtils(0x0),           
fGeometry(0),              fGeoName("EMCAL_COMPLETE12SMV1")
{
  // Constructor
  
  fRecoUtils = new AliEMCALRecoUtils;

}

//_________________________________________________________________
Bool_t AliAnalysisTaskESDfilterEMCALEventSelect::AcceptEventEMCAL()
{
  // Accept event given there is a cluster with enough energy
    
  if(!fGeometry)  fGeometry  = AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1");

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

//_________________________________________________________________
void AliAnalysisTaskESDfilterEMCALEventSelect::UserExec(Option_t *) 
{
  // Main loop

  // Check if the events contains what we want in EMCAL, if not, 
  // do not copy the ESD into AOD
  
  if(!AcceptEventEMCAL()) return ;
  
  // Continue the processing in the same way as in the ESD filter
  
  AliAnalysisTaskESDfilter::UserExec("");
  
}
