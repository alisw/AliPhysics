// @(#) $Id$

// Author: Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"
#include "AliHLTLogging.h"
#include "AliHLTTransform.h"
#include "AliHLTTPCMapping.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//generated from Ulis txt file with read-roc.sh
#include "AliHLTTPCMapping-iroc.generated"
#include "AliHLTTPCMapping-oroc.generated"


Int_t AliHLTTPCMapping::GetRealNPads(Int_t slicerow)
{ 
  //see tpc numbering doc

  if(slicerow<0 || slicerow >= AliHLTTransform::GetNRows()){
    LOG(AliHLTLog::kError,"AliHLTTPCMapping::GetRealNPads","Slicerow")
      <<"Wrong slicerow "<<slicerow<<ENDLOG;
    return -1; 
  }

  const Float_t k1=0.293878;// == 10/6*tan(10)
  const Float_t k2=0.440817;// == 15/6*tan(10)
  const Int_t knRowLow=AliHLTTransform::GetNRowLow();
  const Int_t knRowUp1=AliHLTTransform::GetNRowUp1();

  if(slicerow==0) return 68;
  else if(slicerow<knRowLow){
    Double_t dummy=slicerow/3+33.67;
    return (2*Int_t(dummy));
  }

  Int_t rowup=slicerow-knRowLow;
  if(rowup<knRowUp1){
    Double_t dummy=k1*rowup+37.75;
    return (2*Int_t(dummy));
  }
  else
  {
    Double_t dummy=k2*(rowup-knRowUp1+56.66);
    return (2*Int_t(dummy));
  }
}

Double_t AliHLTTPCMapping::GetRealX(Int_t slicerow)
{ 
  //see tpc numbering doc
  if(slicerow<0 || slicerow >= AliHLTTransform::GetNRows()){
    LOG(AliHLTLog::kError,"AliHLTTPCMapping::GetRealX","Slicerow")
      <<"Wrong slicerow "<<slicerow<<ENDLOG;
    return -1; 
  }

  const Int_t knRowLow=AliHLTTransform::GetNRowLow();
  const Int_t knRowUp1=AliHLTTransform::GetNRowUp1();

  if(slicerow<knRowLow){
    return (85.225+0.75*slicerow);
  }

  Int_t rowup=slicerow-knRowLow;
  if(rowup<knRowUp1){
    return (rowup+135.1);
  }
  else
  {
    return (1.5*(rowup-knRowUp1)+199.35);
  }
}

Double_t AliHLTTPCMapping::GetRealY(Int_t slicerow, Int_t pad)
{ 
  //see tpc numbering doc
  if(slicerow<0 || slicerow >= AliHLTTransform::GetNRows()){
    LOG(AliHLTLog::kError,"AliHLTTransform::GetRealY","Slicerow")
      <<"Wrong slicerow "<<slicerow<<ENDLOG;
    return -1; 
  }

  Int_t npads=GetRealNPads(slicerow);
  if(pad<0 || pad >= npads){
    LOG(AliHLTLog::kError,"AliHLTTransform::GetRealY","pad")
      <<"Wrong pad "<<pad<<" npads " <<npads<<ENDLOG;
    return 0.; 
  }

  const Int_t knRowLow=AliHLTTransform::GetNRowLow();

  if(slicerow<knRowLow){
    return (0.4*pad+0.2-0.2*npads);
    //== (pad-0.5*(npads-1))*fPadPitchWidthLow;
  }
  else
  {
    return (0.6*pad+0.3-0.3*npads);
    // == (pad-0.5*(npads-1))*fPadPitchWidthUp;
  }
}


#if __old__
#include "TPCMapping.h"

TPCMapping::TPCMapping(char* file){
  // constructor
	fin = new ifstream();
	ffile = file;
	kreadfile = 0;
}

TPCMapping::~TPCMapping(){
  // destructor
	fin->close();
	delete fin;
}
void TPCMapping::open(){
  // opens file
	fin->open(ffile);
}
void TPCMapping::isOpen(){
  // dummy
}
void TPCMapping::read(){
  // reads the RORC
	for(int i = 0; i < 5504 ; i++){
		for(int j = 0 ; j < 8 ; j++){
			*fin >> fIRORC[i][j];
		}
	}
	fsizeoffIRORC = 5504;
	kreadfile = 1;
}
void TPCMapping::read(int* listofRCUs, int numofRCU){
  // reads (?)
	int pos = 0;
	for(int i = 0; i < 5504 ; i++){
		for(int j = 0 ; j < 8 ; j++){
			*fin >> fIRORC[pos][j];
		}
		for(int j = 0 ; j < numofRCU ; j++){
			if( fIRORC[pos][5] == listofRCUs[j]){
				pos++;
			}
		}
	}
	fsizeoffIRORC = pos;
	kreadfile = 1;
}

void TPCMapping::print(){
  // debug printout
	cout << " Index  |  Row  |  Pad  |Connect|  Pin  |  FEC  |channel| FECcon|AltroCH| Altro |" << endl;
	for(int i = 0; i < fsizeoffIRORC ; i++){
		print(i);
	}	
}
void TPCMapping::print(int start, int end){
  // debug printout
	cout << " Index  |  Row  |  Pad  |Connect|  Pin  |  FEC  |channel| FECcon|AltroCH| Altro |" << endl;
	for(int i = start; i <= end ; i++){
		print(i);
	}	
}
void TPCMapping::print(int index){
  // debug printout
	cout <<
	getIndex(index) << "\t|" <<
	getPadrow(index) << "\t|" <<
	getPad(index) << "\t|" <<
	getConnector(index) << "\t|" <<
	getPin(index) << "\t|" <<
	getFEC(index) << "\t|" <<
	getFECchannel(index) << "\t|" <<
	getFECconnector(index) << "\t|" <<
	getAltroChannel(index) << "\t|" <<
	getAltro(index) << "\t|" << endl;
}
int TPCMapping::getIndex(int index){
	//COLUMN 0 -> INDEX (0 - 9983)
	int retval;
	if(kreadfile == 1)
		retval = fIRORC[index][0];
	else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}
int	TPCMapping::getPadrow(int index){
	//COLUMN 1 -> PADROW (0 - 95)
	int retval;
	if(kreadfile == 1)
		retval = fIRORC[index][1];
	else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}
int TPCMapping::getPad(int index){
	//COLUMN 2 -> PAD (0 - (Np-1))
	int retval;
	if(kreadfile == 1)
		retval = fIRORC[index][2];
	else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}
int TPCMapping::getConnector(int index){
	//COLUMN 3 -> Connector (1 - 468)
	int retval;
	if(kreadfile == 1)
		retval = fIRORC[index][3];
	else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}
int TPCMapping::getPin(int index){
	//COLUMN 4 -> Pin (0 - 22)
	int retval;
	if(kreadfile == 1)
		retval = fIRORC[index][4];
	else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}
int TPCMapping::getFEC(int index){
	//COLUMN 5 -> FEC (0 - 77)
	int retval;
	if(kreadfile == 1)
		retval = fIRORC[index][5];
	else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}
int TPCMapping::getFECchannel(int index){
	//COLUMN 6 -> FEC Channel (0 - 127)
	int retval;
	if(kreadfile == 1)
		retval = fIRORC[index][6];
	else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}
int TPCMapping::getFECconnector(int index){
	//COLUMN 7 -> FEC Connector (0 - 5)
	int retval;
	if(kreadfile == 1)
		retval = fIRORC[index][7];
	else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}

int TPCMapping::getAltroChannel(int index){
// every the channel to Altro ordering is 0,1,2,3,4,5,6,7,16,15,14,13,12,11,10,9,8
	int retval;
	int channel = getFECchannel(index);
	if(kreadfile == 1){
		if( (channel/8)%2 == 0 ){ //every even block hast to be changed
			retval = channel;
		}else{
			retval = (channel/8)*8 + 7-channel%8; //Blocknumber + reversed internel Number
		}
		if( (channel/16)%2 != 0){
//			retval = 999;
			if( (channel/8)%2 == 0 )
				retval += 8;
			else
				retval -= 8;
		}
	}else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}
int TPCMapping::getAltro(int index){
  // gets Altro 
	int retval;
	int channel = getFECchannel(index);
	if(kreadfile == 1)
		retval = channel/16;
	else{
		retval = 0;
		cout << "ERROR: Array Empty!" << endl;
	}
	return retval;
}

int TPCMapping::getPadsperRow(int row){
	//IRORC Formula !!!!!!!!!!!!!!!!!!!!!!!!!
	int retval;
	if(row == 0){
		retval = 68;
	}else{
		retval = 2*((int)(row/3.0+33.67));
	}
	return retval;
}
void TPCMapping::myprint(){
  // debug printout
	int FEC = 0;
	int channel = 0;
	cout << "channel Index  |  Row  |  Pad  |Connect|  Pin  |  FEC  |channel| FECcon|AltroCH| Altro |" << endl;

	for(int j = 0; j < 512 ; j++){
//		cout << "j: " << j << " fec " << (j/128+29) << " channel: " << j%128 << endl; 
		FEC = (32-j/128);
		channel = (j%128);
		for(int i = 0; i < fsizeoffIRORC ; i++){
			if((getFEC(i) == FEC) && (getAltroChannel(i) == channel)) {
				cout << j << "\t" << getIndex(i) << "\t" << getPadrow(i) << "\t" <<
				getPad(i) << "\t" << getConnector(i) << "\t" << getPin(i) << "\t" <<
				getFEC(i) << "\t" << getFECchannel(i) << "\t" << getFECconnector(i) << "\t" <<
				getAltroChannel(i) << "\t" << getAltro(i) << "\t" << endl;
			}
		}
	}
}
void TPCMapping::myprint1(){
  // debug printout
	int FEC = 0;
	int channel = 0;
	cout << "channel Index  |  Row  |  Pad  |Connect|  Pin  |  FEC  |channel| FECcon|AltroCH| Altro |" << endl;

	for(int j = 0; j < 95 ; j++){
		for(int k = 0; k < getPadsperRow(j) ; k++){
			for(int i = 0; i < fsizeoffIRORC ; i++){
				if(((getFEC(i) >= 29) && (getFEC(i) <= 32)) && (getPadrow(i) == j) && (getPad(i) == k)) {
					cout << j << "\t" << k << "\t" << getIndex(i) << "\t" << getPadrow(i) << "\t" <<
					getPad(i) << "\t" << getConnector(i) << "\t" << getPin(i) << "\t" <<
					getFEC(i) << "\t" << getFECchannel(i) << "\t" << getFECconnector(i) << "\t" <<
					getAltroChannel(i) << "\t" << getAltro(i) << "\t" << endl;
				}
			}
		}
	}
}
void TPCMapping::myprint2(){
  // debug printout
	int FEC = 0;
	int channel = 0;
	int rownum = 0;

	for(int j = 0; j < 95 ; j++){
		for(int i = 0; i < fsizeoffIRORC ; i++){
			if(((getFEC(i) >= 29) && (getFEC(i) <= 32)) && (getPadrow(i) == j) ){
				rownum++;
			}
		}		
		if(rownum != 0){ 
			cout << rownum << "\t" << j <<  "\t";
			for(int k = 0; k < getPadsperRow(j) ; k++){
				for(int i = 0; i < fsizeoffIRORC ; i++){
					if(((getFEC(i) >= 29) && (getFEC(i) <= 32)) && (getPadrow(i) == j) && (getPad(i) == k)) {
						cout << getAltroChannel(i)+(32-getFEC(i))*128 << "\t";
//						cout << getAltroChannel(i) << "\t";
					}
				}
			}
			cout << endl;
		}
		rownum = 0;
	}
}

#endif
