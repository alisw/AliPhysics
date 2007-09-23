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

/**********************************************************************
 Created on : 08/09/2007
 Purpose    : TriggerDDL Regional to LocalBoard Mapping.
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com
**********************************************************************/

#include <iostream>

#include "AliMpDDLStore.h"
#include "AliMpSegmentation.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"

using namespace std;

bool CreateRegionalToLocalCardMapping()
{
  AliMpSegmentation::ReadData(); 
  AliMpDDLStore::ReadData();

  FILE *fout = fopen("RegionalToLocalCard.dat","w");

  for(int iTrigDDL = 0; iTrigDDL<2 ; iTrigDDL++){

    for(int iReg = 0; iReg < 8 ; iReg++){
    
      AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->GetTriggerCrate(iTrigDDL,iReg);
      cout<<"Nof Local Board : "<<crate->GetNofLocalBoards()<<endl;
      AliMpLocalBoard* locBoard ;
      
      for(int i = 0 ; i < crate->GetNofLocalBoards() ; i++){
      locBoard = AliMpDDLStore::Instance()->GetLocalBoard(crate->GetLocalBoardId(i));

      cout<<"iTrigDDL : "<<iTrigDDL<<", iReg : "<<iReg<<", iLoc : "<<i
	  <<", locID : "<<locBoard->GetId()<<", Switch : (" ;
      fprintf(fout,"%1d\t%1d\t%-2d\t%-3d\t",iTrigDDL,iReg,i,locBoard->GetId());

      int switch_packed ;
      switch_packed &= 0x0 ;
      for(int j = 0 ; j < locBoard->GetNofSwitches() ; j++){
	cout<<locBoard->GetSwitch(j);
	
	if(j<(locBoard->GetNofSwitches()-1))
	  switch_packed = (switch_packed | locBoard->GetSwitch(j))<<1 ;
	else
	  switch_packed |= locBoard->GetSwitch(j) ;
      }

      cout<<"), packed switched : "<<switch_packed<<", DDL : ";
      fprintf(fout,"%-4d\t",switch_packed);

      for(int j = 0 ; j < 4 ; j++){ 
	// intensionally made 4 instead of locBoard->GetNofDEs() to give identical structure for the fake local IDs
	if(locBoard->GetNofDEs()<4){
	  cout<<locBoard->GetDEId(0)<<", ";
	  fprintf(fout,"%-4d\t",locBoard->GetDEId(0));
	}else{
	  cout<<locBoard->GetDEId(j)<<", ";
	  fprintf(fout,"%-4d\t",locBoard->GetDEId(j));
	}//if detelem condn
      }// detelem loop

      cout<<endl;
      fprintf(fout,"\n");

      }// loc board loop
      delete crate;
      //delete locBoard;

    }//iReg loop
  }//iTrigDDL loop

  fclose(fout);
  return true;
}
