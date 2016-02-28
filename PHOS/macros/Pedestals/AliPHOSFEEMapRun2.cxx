#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <TObject.h>
#include <TString.h>
#include <TMath.h>
#include "AliPHOSFEEMapRun2.h"

//This class is developed for calibration of PHOS in Run-2 (2015-2017).
//This class is based on hard ware numbering.
//So, cellx should be [0,63] and cellz should be [0,55].
//Module ID should be  [0-4], hard ware numbering. M1:half module,M4:closet module to HMPID.
//Author : Daiki Sekihata (Hiroshima University), 15.August.2015
//daiki.sekihata@cern.ch
//last update : 15.Sep.2015, Daiki Sekihata


//usage 1: if one wants to convert a cell ID to all electronics numbering at once.
//root [0] AliPHOSFEEMapRun2 *map = new AliPHOSFEEMapRun2(2,23,42)
//root [1] map->Print()
//root [2] Int_t sru, fee, altro, csp, hvch, altlg, althg
//root [3] map->GetElectronicsMap(sru,fee,altro,csp,hvch,altlg,althg)

//usage 2: if one wants to convert to a cell ID to an explicit device ID,
//for examle in a for-loop over cells.
//AliPHOSFEEMapRun2 *map = new AliPHOSFEEMapRun2()
//for(Int_t im=0;im<5;im++){
//  for(Int_t ix=0;ix<64;ix++){
//    for(Int_t iz=0;iz<64;iz++){
//
//      map->CellToCSPID(im,ix,iz);
//
//    }
//  }
//} 


using namespace std;

ClassImp(AliPHOSFEEMapRun2)
//________________________________________________________________________
AliPHOSFEEMapRun2::AliPHOSFEEMapRun2()
:fSRUID(-1),
 fFEEID(-1),
 fALTRO(-1),
 fCSPID(-1),
 fHVID(-1),
 fALTCH_LG(-1),
 fALTCH_HG(-1),
 fTRUID(-1),
 fTRUCH(-1)
{
	//default constructor

}
//________________________________________________________________________
AliPHOSFEEMapRun2::AliPHOSFEEMapRun2(Int_t module,Int_t cellx,Int_t cellz)
:fSRUID(-1),
 fFEEID(-1),
 fALTRO(-1),
 fCSPID(-1),
 fHVID(-1),
 fALTCH_LG(-1),
 fALTCH_HG(-1),
 fTRUID(-1),
 fTRUCH(-1)
{
	fSRUID = CellToSRUID(cellx);
	fFEEID = CellToFEEID(cellz);
	fALTRO = CellToALTRO(cellx,cellz);
	fCSPID = CellToCSPID(module,cellx,cellz);
  fHVID  = CSPToHVID(fCSPID);
  fALTCH_LG = CSPToALTROChannel(fCSPID,"Low");
  fALTCH_HG = CSPToALTROChannel(fCSPID,"High");
  fTRUID = CellToTRUID(module,cellx,cellz);
  fTRUCH = CellToTRUChannel(cellx,cellz);

}
//________________________________________________________________________
void AliPHOSFEEMapRun2::Print(Option_t *) const
{

  printf("***** %s Print *****\n SRU ID : %d.\n FEE ID : %d.\n ALTRO ID :%d.\n CSP ID : %d.\n HV ID : %d. (This HV ID should be converted to HEX)\n ALTRO Channel LG : %d.\n ALTRO Channel HG : %d.\n",
  GetName(),fSRUID,fFEEID,fALTRO,fCSPID,fHVID,fALTCH_LG,fALTCH_HG);
  printf("***** %s Print end *****\n",GetName());

}
//________________________________________________________________________
AliPHOSFEEMapRun2::~AliPHOSFEEMapRun2(){

  //destructor


}
//________________________________________________________________________
Int_t AliPHOSFEEMapRun2::CellToSRUID(Int_t cellx)
{
	//cellx should be [0,63].

  Int_t sru = -1;
	if(cellx < 0 || 63 < cellx) return -1;
	else{
    sru = cellx/16;
    return sru;
	}

}
//________________________________________________________________________
Int_t AliPHOSFEEMapRun2::CellToFEEID(Int_t cellz)
{
	//cellz should be [0,55].

	Int_t fee = -1;

	if(cellz < 0 || 55 < cellz) return -1;
	else{
		if(cellz > 27) fee = 1 + (cellz - 27 - 1)/2; //for branch 0
		else           fee = 34 - (cellz)/2; //for branch 1
		return fee;
	}

}
//________________________________________________________________________
Int_t AliPHOSFEEMapRun2::CellToCSPID(Int_t module, Int_t cellx, Int_t cellz)
{
	//cellx should be [0,63].
	//cellz should be [0,55].

	//CSP mapping on 1 "new" FEE card for module 1,2,3.
	//
	//   x=     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
	//     |---------------------------------------------------
	//z=1  |CSP16|17|18|19|20|21|22|23| 8| 9|10|11|12|13|14|15|
	//top  |--------------------------------------------------| bottom
	//z=0  |CSP 0| 1| 2| 3| 4| 5| 6| 7|24|25|26|27|28|29|30|31|
	//     |---------------------------------------------------
	//     |   ALTRO 2    |  ALTRO 3  |  ALTRO 0  |  ALTRO 4  |
	//

	//CSP mapping on 1 "old" FEE card for module 4.
	//
	//   x=     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
	//     |---------------------------------------------------
	//z=1  |CSP16|17|18|19|20|21|22|23|24|25|26|27|28|29|30|31|
	//top  |--------------------------------------------------| bottom
	//z=0  |CSP 0| 1| 2| 3| 4| 5| 6| 7| 8| 9|10|11|12|13|14|15|
	//     |---------------------------------------------------
	//     |   ALTRO 2    |  ALTRO 3  |  ALTRO 0  |  ALTRO 4  |
	//

	if(module == 1 && cellx < 32) return -1;//module 1 is half module.
	if(cellx < 0 || 63 < cellx || cellz < 0 || 55 < cellz) return -1;

	Int_t x = -1;
	Int_t z = -1;
	Int_t csp = -1;

	if(module==4){//module 4
		x = cellx % 16;//x = [0,15]
		z = cellz % 2; //z = [0,1]

		if(z==0) csp = x;
		else if(z==1) csp = x + 16;
		else{
			printf("calculation of z is wrong. return -1.");
			return -1;
		}

		return csp;
	}
	else{//for module 1,2,3
		x = cellx % 16;//x = [0,15]
		z = cellz % 2; //z = [0,1]

		if(z==0){
			if(x < 8) csp = x;
			else csp = x + 16;
		}
		else if(z==1){
			if(x > 7) csp = x;
			else csp = x + 16;
		}
		else{
			printf("calculation of z is wrong. return -1.");
			return -1;
		}

		return csp;
	}

}
//________________________________________________________________________
Int_t AliPHOSFEEMapRun2::CellToALTRO(Int_t cellx, Int_t cellz)
{

	if(cellx < 0 || 63 < cellx || cellz < 0 || 55 < cellz) return -1;

	Int_t x =  cellx % 16;//x = [0,15]
	//Int_t z =  cellz % 2; //z = [0,1]
	
	Int_t altro = -1;
	if(x<4)       altro = 2;
	else if(x<8)  altro = 3;
	else if(x<12) altro = 0;
	else if(x<16) altro = 4;
	else{
		printf("calculation of x is wrong. return -1.");
		return -1;
	}

	return altro;
}
//________________________________________________________________________
Int_t AliPHOSFEEMapRun2::CSPToHVID(Int_t csp)
{

	//csp ID has to be 0-31.
	if(csp < 0 || 31 < csp) return -1;
	else{
		Int_t hvid = -1;

		if(csp < 16)      hvid = 104 + csp;
		else if(csp < 24) hvid = 104 - csp + 15;
		else if(csp < 32) hvid = 127 - csp + 24;

		return hvid;
	}

}
//________________________________________________________________________
Int_t AliPHOSFEEMapRun2::CSPToALTROChannel(Int_t csp, TString gain)
{

  Int_t channel = -1;

  if(gain=="High" || gain=="Low"){
    //csp ID has to be 0-31.
    if(csp < 0 || 31 < csp) return -1;
    else{
      if(gain=="Low"){
        switch(csp){
          case 0   : return 11;
          case 1   : return 15;
          case 2   : return 4;
          case 3   : return 0;
          case 4   : return 0;
          case 5   : return 4;
          case 6   : return 15;
          case 7   : return 11;

          case 8   : return 11;
          case 9   : return 15;
          case 10  : return 4;
          case 11  : return 0;
          case 12  : return 0;
          case 13  : return 4;
          case 14  : return 15;
          case 15  : return 11;

          case 16  : return 9;
          case 17  : return 13;
          case 18  : return 6;
          case 19  : return 2;
          case 20  : return 2;
          case 21  : return 6;
          case 22  : return 13;
          case 23  : return 9;

          case 24  : return 9;
          case 25  : return 13;
          case 26  : return 6;
          case 27  : return 2;
          case 28  : return 2;
          case 29  : return 6;
          case 30  : return 13;
          case 31  : return 9;
          default : return -1;
        }
      }
      else if(gain=="High"){
        switch(csp){
          case 0   : return 10;
          case 1   : return 14;
          case 2   : return 5;
          case 3   : return 1;
          case 4   : return 1;
          case 5   : return 5;
          case 6   : return 14;
          case 7   : return 10;

          case 8   : return 10;
          case 9   : return 14;
          case 10  : return 5;
          case 11  : return 1;
          case 12  : return 1;
          case 13  : return 5;
          case 14  : return 14;
          case 15  : return 10;

          case 16  : return 8;
          case 17  : return 12;
          case 18  : return 7;
          case 19  : return 3;
          case 20  : return 3;
          case 21  : return 7;
          case 22  : return 12;
          case 23  : return 8;

          case 24  : return 8;
          case 25  : return 12;
          case 26  : return 7;
          case 27  : return 3;
          case 28  : return 3;
          case 29  : return 7;
          case 30  : return 12;
          case 31  : return 8;
          default : return -1;
        }

      }

      return channel;
    }
  }
  else{
    printf("gain must be High or Low! return -1.");
    return -1;
  }

}
//________________________________________________________________________
Int_t AliPHOSFEEMapRun2::CellToTRUID(Int_t module, Int_t cellx, Int_t cellz)
{
  if(module < 1 || 5 < module) return -1;// this is based on hardware numbering.
	if(cellx < 0 || 63 < cellx || cellz < 0 || 55 < cellz) return -1;

  if(module==1 && cellx<32) return -1;

  Int_t tru = -1;

  Int_t x = cellx / 16;//x = [0,1,2,3]
  Int_t z = 1 - cellz / 28;//z = [0,1]
  Int_t offset = (module-1)*8-4;

  tru = 2*x + z + offset;

	return tru;
}
//________________________________________________________________________
Int_t AliPHOSFEEMapRun2::CellToTRUChannel(Int_t cellx, Int_t cellz)
{
	if(cellx < 0 || 63 < cellx || cellz < 0 || 55 < cellz) return -1;
  Int_t truch = -1;

  //if(cellz>27) cellz -= 27;

  Int_t z = (27 - (cellz % 28)) / 2;//z = [0-13]
  Int_t x = (cellx%16) / 2;
  truch = z*8 + x;

	return truch;
}
//________________________________________________________________________
void AliPHOSFEEMapRun2::TRUHWToCellID(Int_t ddl, Int_t hwaddress, Int_t &cellx, Int_t &cellz)
{
  if(hwaddress < 0 || (111 < hwaddress && hwaddress < 2048 )|| 2160 < hwaddress ){
    cellx = -1;
    cellz = -1;
    return;
  }
  else{
    Int_t offset_x = 16*((ddl-4)%4);

    if(hwaddress>=2048){// this is for cellz<28
      cellx = 2*((hwaddress-2048) % 8) + offset_x;
      cellz = 28 - 2*((hwaddress-2048) / 8 + 1);
      //cout << "ddl = " << ddl << " , HW = " << hwaddress << " , cellx = " << cellx << " , cellz = " << cellz << endl;
    }
    else{// this is cell for cellz>=28
      cellx = 2 * (hwaddress % 8) + offset_x;
      cellz = 28 - 2 * (hwaddress / 8 + 1) + 28;
    }

    return;
  }

}
//________________________________________________________________________

