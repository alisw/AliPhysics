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

//-------------------------------------------------------------------------
//           AliTPCmapper
//  Origin: Jens Wiechula    (j.wiechula@gsi.de)
//
//  class to map (row,pad,sector) to (RCU, Branch, FEC, Altro, Altro channel)
//  create Altro address 


#include <stdio.h>
#include <TMath.h>
#include <TSystem.h>

#include "AliTPCmapper.h"


ClassImp(AliTPCmapper)


//______________________________________________________________
AliTPCmapper::AliTPCmapper()
{
    Init();
    ReadMapping();
}

//______________________________________________________________
void AliTPCmapper::Init()
{

    //Initialize arrays
    for (Int_t i=0; i<kNrcu; i++){
	for (Int_t j=0; j<kNbranch; j++){
	    for (Int_t k=0; k<kNfecMax; k++){
		for (Int_t l=0; l<kNaltro; l++){
		    for (Int_t m=0; m<kNchannel; m++){
			fAddressToRow[i][j][k][l][m] = -1;
			fAddressToPad[i][j][k][l][m] = -1;
		    }
		}
	    }
	}
    }

    for (Int_t i=0; i<kNpadrow; i++){
	for (Int_t j=0; j<kNpadMax; j++){
	    fRowPadToRCU[i][j] = -1;
	    fRowPadToBranch[i][j] = -1;
	    fRowPadToFEC[i][j] = -1;
	    fRowPadToAltro[i][j] = -1;
	    fRowPadToChannel[i][j] = -1;
            fRowPadToPadsec[i][j] = -1;
	}
    }

    for (Int_t i=0; i<kNpadSector; i++){
	fPadsecToRow[i]=-1;
	fPadsecToPad[i]=-1;
    }

//    for (Int_t i=0; i<kNaddrSize; i++)
//        fAddressArray[i] = 0;

    strcpy(fMapfileName,gSystem->ExpandPathName("$ALICE_ROOT/TPC/final_mapping.txt"));
}

//______________________________________________________________
Int_t AliTPCmapper::ReadMapping()
{
    FILE *fin;
    char line[255];

    int pad = -1, row = -1;
    int rcu = -1, bra = -1, fec = -1, alt = -1, chn = -1;
    int padsec = 0;


    fin = fopen(fMapfileName,"r");
    if (!fin){
	fprintf(stderr, "cannot open file '%s'!\n",fMapfileName);
	return 1;
    }

    fgets(line,256,fin);
    while (sscanf(line,"%d %d %d %d %d %d %d %d",
		  &padsec, &row, &pad,
		  &rcu, &bra, &fec, &alt, &chn
		 )!=8){
        fgets(line, 256, fin);
	fprintf(stderr,"%s",line);
    }

    while (!feof(fin)){
        sscanf(line,"%d %d %d %d %d %d %d %d",
	       &padsec, &row, &pad,
	       &rcu, &bra, &fec, &alt, &chn
	      );
	fAddressToRow[rcu][bra][fec][alt][chn] = row;
	fAddressToPad[rcu][bra][fec][alt][chn] = pad;

	fRowPadToRCU[row][pad]     = rcu;
	fRowPadToBranch[row][pad]  = bra;
	fRowPadToFEC[row][pad]     = fec;
	fRowPadToAltro[row][pad]   = alt;
	fRowPadToChannel[row][pad] = chn;
	fRowPadToPadsec[row][pad]  = padsec;

	fPadsecToRow[padsec]        = row;
        fPadsecToPad[padsec]        = pad;



	fgets(line, 256, fin);
    }

    fclose(fin);

    return 0;
}

//______________________________________________________________
void AliTPCmapper::PrintRBFACinfo(Int_t row, Int_t pad)
{
    fprintf(stderr,"RCU: %d, Branch: %d, FEC: %2d, Altro: %d, Channel: %2d\n",
	    GetRCUs(row,pad), GetBranchS(row,pad),
	    GetFECs(row,pad), GetAltroS(row,pad),
	    GetChannelS(row,pad));

}

//______________________________________________________________
Int_t AliTPCmapper::GetPadsInRowS(Int_t row) const{
  //
  //GetPadsInRowS
  //
    if ( row == 0 )
	return 68;

    if ( row < 63 )
	return 2*(Int_t)(row/3.+33.67);

    row -= 63;
    Double_t k1=10./6.*tan(TMath::Pi()/18.);

    if ( row < 64 )
	return 2*(Int_t)(k1 * row + 37.75);

    Double_t k2=15./6.*tan(TMath::Pi()/18.);
    return 2*(Int_t)(k2 * (row-64) + 56.66);
}

//______________________________________________________________
Double_t AliTPCmapper::GetPadXlocalS(Int_t row, Int_t pad) const {
  //
  //GetPadXlocalS
  //
    if ( row < 63 ) //IROC
	return (852.25 + 7.5 * (Double_t)row)*1.e-1; //divide by 10 to get cm

    row -= 63;

    if ( row < 64 ) //OROC inner part
	return (10.* row + 1351.)*1.e-1;  //divide by 10 to get cm


    //OROC outer part
    return (15.*(row - 64) + 1993.5)*1.e-1;  //divide by 10 to get cm
}

//______________________________________________________________
Double_t AliTPCmapper::GetPadXlocalS(Int_t padsector) const{
  //
  //GetPadXlocalS
  //
    Int_t row=GetRowFromPadSector(padsector);
    Int_t pad=GetPadFromPadSector(padsector);
    return GetPadXlocalS(row,pad);
}

//______________________________________________________________
Double_t AliTPCmapper::GetPadYlocalS(Int_t row, Int_t pad) const{
  //
  //:GetPadYlocalS
  //
    Int_t padsInRow = GetPadsInRowS(row);

    if ( row < 63 ) //IROC
	return (2.* padsInRow - 4.*pad - 2.)*1.e-1;  //divide by 10 to get cm

    //OROC
        return (3.* padsInRow -6.*pad - 3.)*1.e-1;  //divide by 10 to get cm
}

//______________________________________________________________
Double_t AliTPCmapper::GetPadYlocalS(Int_t padsector) const{
  //
  //:GetPadYlocalS
  //
    Int_t row = GetRowFromPadSector(padsector);
    Int_t pad = GetPadFromPadSector(padsector);
    return GetPadYlocalS(row,pad);
}

//______________________________________________________________
Double_t AliTPCmapper::GetPadXglobalS(Int_t row, Int_t pad,Int_t sector) const{
  //
  // GetPadXglobalS
  //
    Double_t angle = (Double_t)(( sector * 20. ) +10. ) * TMath::DegToRad();
    return GetPadXlocalS(row,pad)*TMath::Cos(angle) -
	GetPadYlocalS(row,pad)*TMath::Sin(angle);
}

//______________________________________________________________
Double_t AliTPCmapper::GetPadYglobalS(Int_t row, Int_t pad,Int_t sector) const{
  //
  // GetPadYglobalS
  //
    Double_t angle = (Double_t)(( sector * 20. ) + 10. ) * TMath::DegToRad();
    return GetPadXlocalS(row,pad)*TMath::Sin(angle) +
	GetPadYlocalS(row,pad)*TMath::Cos(angle);
}

//______________________________________________________________
Double_t AliTPCmapper::GetPadWidthS(Int_t row) const
{
  //
  // :GetPadWidthS
  //
    if (row < 63 ) return .4;
    return .6;
}

//______________________________________________________________
Double_t AliTPCmapper::GetPadLengthS(Int_t row) const
{
  //
  // GetPadLengthS
  //

    if ( row < 63  ) return  .75;
    if ( row < 127 ) return 1.;
    return 1.5;
}

//______________________________________________________________
Int_t AliTPCmapper::GetAltroAddrwPatch(const Int_t row,const Int_t pad) const
{
  //
  // :GetAltroAddrwPatch
  //
    return GetChannelS(row,pad)+
	(GetAltroS (row,pad) <<4)+
	(GetFECs   (row,pad) <<7)+
	(GetBranchS(row,pad) <<11)+
        (GetRCUs   (row,pad) <<12);
}

//______________________________________________________________
Int_t AliTPCmapper::GetAltroAddrwPatch(const Int_t padsector) const
{
  //
  // GetAltroAddrwPatch
  //
    return GetAltroAddrwPatch(GetRowFromPadSector(padsector),GetPadFromPadSector(padsector));
}

//______________________________________________________________
Int_t AliTPCmapper::GetRow(Int_t altroaddr) const
{
  //
  // GetRow
  //
    Int_t rcu = (altroaddr>>12)&0x07;
    Int_t bra = (altroaddr>>11)&0x01;
    Int_t fec = (altroaddr>>7)&0x0F;
    Int_t alt = (altroaddr>>4)&0x07;
    Int_t chn = (altroaddr)&0x0F;
    return GetRow(rcu,bra,fec,alt,chn);
}

//______________________________________________________________
Int_t AliTPCmapper::GetPad(Int_t altroaddr) const
{
  //
  // GetPad
  //
    Int_t rcu = (altroaddr>>12)&0x07;
    Int_t bra = (altroaddr>>11)&0x01;
    Int_t fec = (altroaddr>>7)&0x0F;
    Int_t alt = (altroaddr>>4)&0x07;
    Int_t chn = (altroaddr)&0x0F;
    return GetPad(rcu,bra,fec,alt,chn);
}

//______________________________________________________________
void AliTPCmapper::PrintAddressArray(Int_t row, Int_t pad)
{
  //
  // PrintAddressArray
  //
    Bool_t a[kNaddrSize];

    Int_t addr = GetAltroAddrwPatch(row,pad);
    for (Int_t i=0; i<kNaddrSize; i++)
	a[i] = addr&(Int_t)TMath::Power(2,i);

    fprintf(stderr,"Par|Bro|BC/|Bra|   FEC HW   |  Altro  |Altro Chann.|  Instruction  \n");
    fprintf(stderr,"ity|adc|AL |nch|  Address   |Chip Addr|   Address  |      Code     \n");
    fprintf(stderr,"-------------------------------------------------------------------\n");
    fprintf(stderr," %d | %d | %d | %d | %d  %d  %d  %d | %d  %d  %d | %d  %d  %d  %d | %d  %d  %d  %d  %d\n",
	    a[19],a[18],a[17],a[16],a[15],a[14],a[13],a[12],a[11],
	    a[10],a[9],a[8],a[7],a[6],a[5],a[4],a[3],a[2],a[1],a[0]);

    fprintf(stderr,"-------------------------------------------------------------------\n");
    fprintf(stderr,"19 |18 |17 |16 |15 14 13 12 |11 10  9 | 8  7  6  5 | 4  3  2  1  0\n\n");

}

//______________________________________________________________
AliTPCmapper::~AliTPCmapper()
{
//    delete fAddressArray;

}
