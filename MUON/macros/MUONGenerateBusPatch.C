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

// $Id$

#if !defined(__CINT__) || defined(__MAKECINT__)
// ROOT includes
#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>

// MUON includes
#include "AliMpBusPatch.h"

#endif

/// \ingroup macros
/// \file MUONGenerateBusPatch.C
/// \brief Generates buspatch id and DDL id for given detection element id
///
/// - station 1 & 2 assuming 24 buspatches per quadrant         
/// - station345, reading from DetElemIdToSlatType.dat file and calculates
///   the number of bus patches per slat (and also number of Translator and Bridge Boards).
/// Generates an output file DetElemIdToBusPatch.dat.out, preserve from overwriting
///
/// List of changes:
/// - July 05, first version
/// - Nov. 05, added DDL
/// - June 06, correction for St123
/// - Feb. 07, add 1st manu list for St12 (starting on NB !) and new ddl sharing for station 3)
/// - June. 07, new numbering of station 345, and buspatch starts at 1)
///
/// \author Ch. Finck, July 05


void MUONGenerateBusPatch()
{

  TString dirPath2 = gSystem->Getenv("ALICE_ROOT");
  dirPath2 += "/MUON/mapping/data/"; 

  TString dirPath1 = dirPath2 + "station345/";

  TString infile =  dirPath1 + "DetElemIdToSlatType.dat";
  TString outfile = dirPath2 + "DetElemIdToBusPatch.dat.out";

  ifstream in(infile, ios::in);
  ofstream out(outfile, ios::out);
  char line[80];
  char nameSlat[82];
  Int_t idDE;
  Int_t i;
  Int_t cursor = 0;
  Int_t begin[900];
  Int_t end[900];
  Int_t nbBusPatch;
  Int_t nbTB = 0;
  Int_t nbBB = 0;
  Int_t idSt12[] = {100, 101, 102, 103, 
		    200, 201, 202, 203, 
		    300, 301, 302, 303,
		    400, 401, 402, 403};

  Char_t manuListSt1[] 
      = " 1,27,53,79,105,131,157,183,201,214,224,232,1025,1051,1077,1103,1129,1155,1181,1207,1225,1238,1249,1257";

  Char_t manuListSt2[] 
      = " 1,27,53,79,105,131,157,183,201,214,226,246,1025,1051,1077,1103,1129,1155,1181,1207,1225,1238,1251,1269";



  Int_t iDDL = 0;
  // station 1 & 2
  nbBusPatch = 24;
  Int_t nbHalfBusPatch =  nbBusPatch/2;
  cout << "#DE BusPatch DDL SlatName" << endl;
  out << "#DE BusPatch DDL  1st manu in buspatch" << endl;

  for (Int_t j = 0; j < 16; j++) {

    idDE = idSt12[j];
    if (idDE % 100 == 0) {
      iDDL++;
      begin[cursor] = AliMpBusPatch::GetGlobalBusID(1, iDDL-1);
      cout << "# Chamber " << idDE/100 << endl;
      out  << "# Chamber " << idDE/100 << endl;
    }
    if (idDE % 100 == 1) { 
      iDDL++;
      begin[cursor] = AliMpBusPatch::GetGlobalBusID(1, iDDL-1);
    }
   
    if (idDE % 100 == 3) {
      iDDL--;
      begin[cursor] = AliMpBusPatch::GetGlobalBusID(1, iDDL-1) + nbBusPatch;
    }

    end[cursor]     = begin[cursor] + nbBusPatch - 1;
    ++cursor; // it seems that vec[++i] does not work in Cint
    begin[cursor] = end[cursor-1] + 1;

    cout << idDE << " " << begin[cursor-1] + nbHalfBusPatch << "-" <<end[cursor-1] << ";" <<
	begin[cursor-1]	<< "-" << end[cursor-1] - nbHalfBusPatch << " " << iDDL-1 << endl;

    if (idDE < 300 )
	out << idDE  << " " << begin[cursor-1] + nbHalfBusPatch << "-" <<end[cursor-1] << ";" <<
	    begin[cursor-1]	<< "-" << end[cursor-1] - nbHalfBusPatch << " " << iDDL-1 <<
	    manuListSt1 << endl;
    else
	out << idDE  << " " << begin[cursor-1] + nbHalfBusPatch << "-" <<end[cursor-1] << ";" <<
	    begin[cursor-1]	<< "-" << end[cursor-1] - nbHalfBusPatch << " " <<  iDDL-1 <<
	    manuListSt2  << endl;
    if (idDE % 100 == 3) iDDL++;

  }
  
  // station 345
  nbBusPatch = 0;
  Int_t idCh5swp1 = 501; // 1/4 chamber for DDL on horizontal
  Int_t idCh5swp2 = 505; // 1/4 chamber for DDL on horizontal
  Int_t idCh5swp3 = 510; // 1/4 chamber for DDL on horizontal
  Int_t idCh5swp4 = 514; // 1/4 chamber for DDL on horizontal

  Int_t idCh6swp1 = 600; // 1/4 chamber for DDL on horizontal
  Int_t idCh6swp2 = 605; // 1/4 chamber for DDL on horizontal
  Int_t idCh6swp3 = 609; // 1/4 chamber for DDL on horizontal
  Int_t idCh6swp4 = 614; // 1/4 chamber for DDL on horizontal

  Int_t idSt45swp1 = 7;  // half chamber for DDL in vertical cutting twice the official numbering
  Int_t idSt45swp2 = 20;

  Int_t offsetBusCh5swp1 = 13; // number of buspatches between DE 501-504
  Int_t offsetBusCh5swp2 = 17; // number of buspatches before DE 510

  Int_t offsetBusCh6swp1 = 17; // number of buspatches before DE 600-604
  Int_t offsetBusCh6swp2 = 13; // number of buspatches between DE 614-617

  Int_t offsetBusSt4 = 25;     // number of buspatches before DE 700-800
  Int_t offsetBusSt4swp = 46;  // number of buspatches between 801-806 (701-706)

  Int_t offsetBusSt5 = 27;     // number of buspatches before DE 900-1000
  Int_t offsetBusSt5swp = 50;  // number of buspatches between 901-906 (1001-1006)

  // reads from file
  while ( in.getline(line,80) ) {

    if ( line[0] == '#' || line[0] == '\0' ) continue;
    i = 0;
    sscanf(line, "%d %s", &idDE, nameSlat);

    // number of PCB's
    while(nameSlat[i] < '4' && nameSlat[i] !='0') i++;
    Int_t len = strlen(nameSlat);

    // number of buspatch 
    if (i == 2 || i == 3)
      nbBusPatch = 2;
    else
      nbBusPatch = 4;

    // calculate the number of TB & BB
    // R1: 2 bridges more per PCB
    // R2 & R3: 1 translator more per PCB
    if (i == 2) nbBB += 2; 
    if (i == 3) nbBB += 4; 
    if (i == 4) nbBB += 4; 
    if (i == 5) nbBB += 6;
    if (i == 6) nbBB += 8;

    nbTB += nbBusPatch;
    if (nameSlat[len-1] == '1') {
      nbTB += 2;
      nbBusPatch+=2;
    }
    if (nameSlat[len-1] == '2' || nameSlat[len-1] == '3') {
      nbTB += 1;
      nbBusPatch++;
    }

    // station 3
    // for buspatch length reasons, one ddl connects one 1/4 of two chambers
    // really messy isn't it ? 
    // much more with the new DDL sharing for station 3

    if (idDE < 700 ) {
  
      if (idDE % 100 == 0) {
	cout << "# Chamber " << idDE/100 << endl;
	out  << "# Chamber " << idDE/100 << endl;
      }

      // chamber 5
      if (idDE  == 500 ) {
	iDDL = (idDE/100)*2+2;
	begin[cursor] = AliMpBusPatch::GetGlobalBusID(1, iDDL-1);
	
	end[cursor]   = begin[cursor] + nbBusPatch - 1;
	++cursor;
	begin[cursor] = end[cursor-1] + 1;

	cout << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] <<" " << iDDL-1 <<endl;
      }

      if (idDE  >= idCh5swp1 && idDE < idCh5swp2) {
	
	if ( idDE ==  idCh5swp1) {
	  iDDL = (idDE/100)*2 ;
	  begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusCh5swp1, iDDL-1);
	}
	
	end[cursor]   = begin[cursor] - nbBusPatch + 1;
	++cursor;
	begin[cursor] = end[cursor-1] - 1;

	cout << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] <<" " << iDDL-1 <<endl;
      }

      if (idDE >=  idCh5swp2  && idDE  < idCh5swp3) {

	if (idDE ==  idCh5swp2) {
	  iDDL = (idDE/100)*2-1;
	      begin[cursor] = AliMpBusPatch::GetGlobalBusID(1, iDDL-1);	  
	}

	end[cursor]   = begin[cursor] + nbBusPatch - 1;
	++cursor;
	begin[cursor] = end[cursor-1] + 1;

	cout << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] <<" " << iDDL-1 <<endl;
      }

      if (idDE >=  idCh5swp3  && idDE  < idCh5swp4) {

	if (idDE ==  idCh5swp3) {
	  iDDL = (idDE/100)*2+1;
	      begin[cursor] = AliMpBusPatch::GetGlobalBusID(1, iDDL-1);	  
	}

	end[cursor]   = begin[cursor] + nbBusPatch - 1;
	++cursor;
	begin[cursor] = end[cursor-1] + 1;

	cout << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] <<" " << iDDL-1 <<endl;
      }

     if (idDE  >= idCh5swp4 && idDE < 600) {
	
	if ( idDE ==  idCh5swp4) {
	  iDDL = (idDE/100)*2+2 ;
	  begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusCh5swp2, iDDL-1);
	}
	
	end[cursor]   = begin[cursor] - nbBusPatch + 1;
	++cursor;
	begin[cursor] = end[cursor-1] - 1;

	cout << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] <<" " << iDDL-1 <<endl;
      }

     // chamber 6
      if (idDE  >= idCh6swp1 && idDE < idCh6swp2) {
	
	if ( idDE ==  idCh6swp1) {
	  iDDL = (idDE/100)*2-2 ;
	  begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusCh5swp1+offsetBusCh6swp1, iDDL-1);
	}
	
	end[cursor]   = begin[cursor] - nbBusPatch + 1;
	++cursor;
	begin[cursor] = end[cursor-1] - 1;

	cout << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] <<" " << iDDL-1 <<endl;
      }

      if (idDE >=  idCh6swp2  && idDE  < idCh6swp3) {

	if (idDE ==  idCh6swp2) {
	  iDDL = (idDE/100)*2-3;
	      begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusCh5swp2+1, iDDL-1);	  
	}

	end[cursor]   = begin[cursor] + nbBusPatch - 1;
	++cursor;
	begin[cursor] = end[cursor-1] + 1;

	cout << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] <<" " << iDDL-1 <<endl;
      }

      if (idDE >=  idCh6swp3  && idDE  < idCh6swp4) {

	if (idDE ==  idCh6swp3) {
	  iDDL = (idDE/100)*2-1;
	      begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusCh5swp1+1, iDDL-1);	  
	}

	end[cursor]   = begin[cursor] + nbBusPatch - 1;
	++cursor;
	begin[cursor] = end[cursor-1] + 1;

	cout << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] <<" " << iDDL-1 <<endl;
      }

     if (idDE  >= idCh6swp4 && idDE < 700) {
	
	if ( idDE ==  idCh6swp4) {
	  iDDL = (idDE/100)*2 ;
	  begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusCh5swp2+offsetBusCh6swp2, iDDL-1);
	}
	
	end[cursor]   = begin[cursor] - nbBusPatch + 1;
	++cursor;
	begin[cursor] = end[cursor-1] - 1;

	cout << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] <<" " << iDDL-1 <<endl;
      }

    }
 
    // station 4 & 5
    // back to normal one ddl connects one 1/2 chamber
    // finally !
 
    if (idDE >= 700) {

      if (idDE % 100 == 0) {
	cout << "# Chamber " << idDE/100 << endl;
	out  << "# Chamber " << idDE/100 << endl;
      }

      if ((idDE % 100) >= idSt45swp1 && (idDE % 100) < idSt45swp2) {
	if ((idDE % 100) == idSt45swp1) {
	  iDDL = (idDE/100)*2-1;
	  begin[cursor] = AliMpBusPatch::GetGlobalBusID(1, iDDL-1);
	}
	end[cursor]   = begin[cursor] + nbBusPatch - 1;
	++cursor;
	begin[cursor] = end[cursor-1] + 1;

	cout << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] <<" " << iDDL-1 <<endl;
      }

      if (idDE % 100 >= idSt45swp2) {

	if ((idDE % 100) == idSt45swp2 ) {
	  iDDL = (idDE/100)*2 ;
	  if (idDE/100 == 7 || idDE/100 == 8)
	      begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusSt4swp, iDDL-1);
	  else
	      begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusSt5swp, iDDL-1);

	}
	end[cursor]   = begin[cursor] - nbBusPatch + 1;
	++cursor;
	begin[cursor] = end[cursor-1] - 1;

	cout << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] <<" " << iDDL-1 <<endl;
      }

      if ((idDE % 100) >= 0 && (idDE % 100) < idSt45swp1) {

	if ((idDE % 100) == 0 ) {
	  iDDL = (idDE/100)*2;
	  if (idDE/100 == 7 || idDE/100 == 8)
	      begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusSt4, iDDL-1);
	  else 
	      begin[cursor] = AliMpBusPatch::GetGlobalBusID(offsetBusSt5, iDDL-1);
	  
	}
	end[cursor]   = begin[cursor] - nbBusPatch + 1;
	++cursor;
	begin[cursor] = end[cursor-1] - 1;

	cout << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] << " " << iDDL-1 << " " <<nameSlat <<endl;
	out << idDE << " " << end[cursor-1]<<"-"<<begin[cursor-1] <<" " << iDDL-1 <<endl;
      }
 
    }

  }
  printf(" Slat: number of TB %d and BB %d\n", nbTB, nbBB);
  in.close();
  out.close();  
}
