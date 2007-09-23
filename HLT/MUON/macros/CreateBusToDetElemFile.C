/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Indranil Das <indra.das@saha.ac.in>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 
/*********************************************
Purpose:  A macro to generate buspatch to ddl mapping file 

Created:  7/10/2005
Modified: 22/12/2005
Modified: 09/02/2006
Modified: 09/04/2007
Modified: 24/08/2007 (To adopt to AliRoot v4-06-Release)

Author:   Indranil Das, HEP, SINP, Kolkata
Email:    indra.das@saha.ac.in
***********************************************/


#include <TArrayI.h>
#include <AliMpDDLStore.h>
#include <AliMpDetElement.h>
#include <AliMpSegmentation.h>

int CreateBusToDetElemFile()
{
  FILE *fp = fopen("BusToDetElem.dat","w");

  AliMpDetElement* fDetElement;
  AliMpSegmentation::ReadData(); 
  AliMpDDLStore::ReadData();

  fprintf(fp,"#DE\tBusPatch\tDDL\n");
  for(int ch=7;ch<=10;ch++){
    fprintf(fp,"# Chamber %d\n",ch);
    for(int i=0; i<26 ; i++){
      fDetElement = AliMpDDLStore::Instance()->GetDetElement(ch*100 + i);
       fprintf(fp,"%d\t%d - %d\t%d\n",ch*100 + i,fDetElement->GetBusPatchId(0),
	       fDetElement->GetBusPatchId(fDetElement->GetNofBusPatches()-1),fDetElement->GetDdlId());
    }
  }
  //delete fDDLStore;
  delete fDetElement;

  fclose(fp);
  return 0;
}
