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


void MUONGenerateBusPatch()
{
  // Generates buspatch id and DDL id for given detection element id
  // station 1 & 2 assuming 24 buspatches per quadrant
  // station345, reading from DetElemIdToSlatType.dat file and calculates
  // the number of bus patches per slat (and also number of Translator and Bridge Boards).
  // Generates an output file DetElemIdToBusPatch.dat.out, preserve from overwriting
  // (Ch. Finck, July 05)
  // (Nov. 05: added DDL)

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
  Int_t begin[800];
  Int_t end[800];
  Int_t nbBusPatch;
  Int_t listDE[800];
  Int_t nbTB = 0;
  Int_t nbBB = 0;
  Int_t idSt12[] = {100, 101, 102, 103, 
		    200, 201, 202, 203, 
		    300, 301, 302, 303,
		    400, 401, 402, 403};

  Int_t idSt3swp = 9; // half chamber for DDL on horizontal
  Int_t idSt45swp1 = 7; // half chamber for DDL in vertical cutting twice the official numbering
  Int_t idSt45swp2 = 20;

  Int_t iDDL = 0;
  // station 1 & 2
  nbBusPatch = 24;
  cout << "#DE BusPatch DDL SlatName" << endl;
  out << "#DE BusPatch DDL " << endl;

  for (Int_t j = 0; j < 16; j++) {

    idDE = idSt12[j];
    if (idDE % 100 == 0) {
      begin[cursor] = idDE;
      cout << "# Chamber " << idDE/100 << endl;
      out  << "# Chamber " << idDE/100 << endl;
    }
    if (idDE % 2 == 0) iDDL++;
    end[cursor]     = begin[cursor] + nbBusPatch - 1;
    begin[++cursor] = end[cursor] + 1;

    cout << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1]  <<" " << iDDL << endl;
    out << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] <<" " << iDDL  <<endl;
  }
  
  // station 345
  nbBusPatch = 0;
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
    if (nameSlat[len-1] == '1')
      nbBB += 4;
    if (nameSlat[len-1] == '2' || nameSlat[len-1] == '3') {
      nbTB += 1;
      nbBusPatch++;
    }
    listDE[cursor] = idDE;
    if (idDE % 100 == 0) {
      iDDL++;
      if (idDE >=800)
	iDDL++;
      begin[cursor] = idDE;
      cout << "# Chamber " << idDE/100 << endl;
      out  << "# Chamber " << idDE/100 << endl;
    }
    if (idDE == 500+idSt3swp || idDE == 600+idSt3swp )
	iDDL++;

    if (idDE == 700+idSt45swp1 || idDE == 800+idSt45swp1 || idDE == 900+idSt45swp1 || idDE == 1000+idSt45swp1 )
	iDDL++;

    if (idDE == 700+idSt45swp2 || idDE == 800+idSt45swp2 || idDE == 900+idSt45swp2 || idDE == 1000+idSt45swp2 )
	iDDL--;

    end[cursor]     = begin[cursor] + nbBusPatch - 1;
    begin[++cursor] = end[cursor] + 1;

    cout << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] << " " << iDDL << " " <<nameSlat <<endl;
    out << idDE << " " << begin[cursor-1]<<"-"<<end[cursor-1] <<" " << iDDL <<endl;
  }
  printf(" Slat: number of TB %d and BB %d\n", nbTB, nbBB);
  in.close();
  out.close();  
}
