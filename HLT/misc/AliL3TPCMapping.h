// @(#) $Id$

#ifndef ALIL3TPCMAPPING_H
#define ALIL3TPCMAPPING_H

#include "AliL3RootTypes.h"

#define fNIROC_def 5504 //see ulis file
#define fNOROC_def 9984 //see ulis file


class AliL3TPCMapping {

 private:
  static const Int_t fNIROC; 
  static const Int_t fIRows[fNIROC_def];
  static const Int_t fIPad[fNIROC_def];
  static const Int_t fICon[fNIROC_def];
  static const Int_t fIPin[fNIROC_def];
  static const Int_t fIFec[fNIROC_def];
  static const Int_t fIFecChannel[fNIROC_def];
  static const Int_t fIFecCon[fNIROC_def];

  static const Int_t fNOROC; 
  static const Int_t fORows[fNOROC_def];
  static const Int_t fOPad[fNOROC_def];
  static const Int_t fOCon[fNOROC_def];
  static const Int_t fOPin[fNOROC_def];
  static const Int_t fOFec[fNOROC_def];
  static const Int_t fOFecChannel[fNOROC_def];
  static const Int_t fOFecCon[fNOROC_def];

 public:
  //taken from GSI TPC numbering document
  static Int_t GetRealNPads(Int_t slicerow);           //Number of pads per row
  static Double_t GetRealX(Int_t slicerow);            //Local X in cm for modules 0,36
  static Double_t GetRealY(Int_t slicerow, Int_t pad); //Local Y in cm for modules 0,36

  ClassDef(AliL3TPCMapping,1)
};
#endif


#if __old__

#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
/*
g++ -O3 -c -o TPCMapping.o TPCMapping.C
g++ -O3 TPCMapping.o -o TPCMappingmain.app TPCMappingmain.C
*/

class TPCMapping
{
	public:
	TPCMapping(char *file);
	~TPCMapping();
	void open();
	void isOpen();
	void read();
	void read(int* listofRCUs, int numofRCU);
	void print();
	void print(int index);
	void print(int start, int end);
	void myprint();
	void myprint1();
	void myprint2();
	//COLUMN 0 -> INDEX (0 - 9983)
	//COLUMN 1 -> PADROW (0 - 95)
	//COLUMN 2 -> PAD (0 - (Np-1))
	//COLUMN 3 -> Connector (1 - 468)
	//COLUMN 4 -> Pin (0 - 22)
	//COLUMN 5 -> FEC (0 - 77)
	//COLUMN 6 -> FEC Channel (0 - 127)
	//COLUMN 7 -> FEC Connector (0 - 5)
	int getIndex(int index);
	int	getPadrow(int index);
	int getPad(int index);
	int getConnector(int index);
	int getPin(int index);
	int getFEC(int index);
	int getFECchannel(int index);
	int getFECconnector(int index);
	int getAltroChannel(int index);
	int getAltro(int index);

	int getPadsperRow(int row);
	
	private:
	int kreadfile;
//	int kpartialread;
	ifstream *fin;
	char* ffile;
	short fIRORC [5504][8];
	int fsizeoffIRORC;
};


#endif
