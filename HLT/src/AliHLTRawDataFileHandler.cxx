// @(#) $Id$

// Author: C. Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTRootTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTTransform.h"
#include "AliHLTMemHandler.h"
#include "AliHLTDigitData.h"

#include "AliHLTRawDataFileHandler.h"

#if __GNUC__ >= 3
using namespace std;
#endif

/** \class AliHLTRawDataFileHandler 
<pre>
//_____________________________________________________________
// AliHLTRawDataFileHandler
//
</pre>
*/

ClassImp(AliHLTRawDataFileHandler)

AliHLTRawDataFileHandler::AliHLTRawDataFileHandler()
{
  //constructor
  fConvert=kTRUE;
  fInRaw = 0;
  fInRawPed = 0;
  fMapping = 0;
  fPedestals=0;
  fCharges=0;
  fOutRaw = 0;
  fRow=0;
  fPad=0;
  fRowPad=0;

  FreeAll();

  if((sizeof(Int_t) != 4) || (sizeof(Short_t) != 2)){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::AliHLTRawDataFileHandler","Constructor")
      <<"Check architecture to run the conversion on! Int_t should be 32 and Short_t should be 16 bit."<<ENDLOG;
  }
}

AliHLTRawDataFileHandler::~AliHLTRawDataFileHandler()
{
  //destructor
  FreeAll();
}

void AliHLTRawDataFileHandler::FreeAll()
{
  //free all heap
  if(fInRaw) CloseRawInput();
  if(fInRawPed) CloseRawPedestalsInput();
  if(fMapping) CloseMappingFile();
  if(fNChannels){
    delete[] fRow;
    delete[] fPad;
    delete[] fRowPad;
  }
  if(fPedestals) delete[] fPedestals;
  if(fCharges) delete[] fCharges;
  if(fOutRaw) CloseRawOutput();
  fConvert=kTRUE;
  fInRaw = 0;
  fInRawPed = 0;
  fMapping = 0;
  fPedestals=0;
  fCharges=0;
  fOutRaw = 0;
  fRow=0;
  fPad=0;
  fRowPad=0;
  fNChannels=0;
  fRowMinUsed=AliHLTTransform::GetNRows();
  fRowMaxUsed=0;
  fPadMinUsed=255;
  fPadMaxUsed=0;
  fNTimeBins=0;
  for(Int_t i=0;i<AliHLTTransform::GetNRows();i++) fNPads[i]=0;
  fPedVal=0;
}

Bool_t AliHLTRawDataFileHandler::SetRawInput(Char_t *name)
{
  //set raw input
  if(fInRaw){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawInput","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  //Open the raw data file with name.
  fInRaw = new ifstream();
#ifndef __DECCXX
  fInRaw->open(name,fstream::binary);
#else
  fInRaw->open(name);
#endif

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fInRaw->rdbuf()->is_open()){
#else
  if(!fInRaw->is_open()){
#endif
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawInput","File Open")
      <<"Pointer to ifstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t AliHLTRawDataFileHandler::SetRawInput(ifstream *file)
{
  //set raw input
  if(fInRaw){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawInput","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  //Open the raw data file with given file.
  fInRaw = file;
#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fInRaw->rdbuf()->is_open()){
#else
  if(!fInRaw->is_open()){
#endif
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawInput","File Open")
      <<"Pointer to ifstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

void AliHLTRawDataFileHandler::CloseRawInput()
{
  //close raw input
  if(!fInRaw){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::CloseRawInput","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
#if defined(__HP_aCC) || defined(__DECCXX)
  if(fInRaw->rdbuf()->is_open()) fInRaw->close();
#else
  if(fInRaw->is_open()) fInRaw->close();
#endif
  delete fInRaw;
  fInRaw = 0;
}

Bool_t AliHLTRawDataFileHandler::SetRawOutput(Char_t *name)
{
  //set raw output
  if(fOutRaw){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawOutput","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  fOutRaw = new ofstream();
#ifndef __DECCXX
  fOutRaw->open(name,fstream::binary);
#else
  fOutRaw->open(name);
#endif

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fOutRaw->rdbuf()->is_open()){
#else
  if(!fOutRaw->is_open()){
#endif
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawOutput","File Open")
      <<"Pointer to ofstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t AliHLTRawDataFileHandler::SetRawOutput(ofstream *file)
{
  //set raw output
  if(fOutRaw){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawOutput","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  //Open the raw data file with given file.
  fOutRaw = file;

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fOutRaw->rdbuf()->is_open()){
#else
  if(!fOutRaw->is_open()){
#endif
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawOutput","File Open")
      <<"Pointer to ofstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

void AliHLTRawDataFileHandler::CloseRawOutput()
{
  //close raw output
  if(!fOutRaw){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::CloseRawOutput","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
#if defined(__HP_aCC) || defined(__DECCXX)
  if(fOutRaw->rdbuf()->is_open()) fOutRaw->close();
#else
  if(fOutRaw->is_open()) fOutRaw->close();
#endif
  delete fOutRaw;
  fOutRaw = 0;
}


Bool_t AliHLTRawDataFileHandler::SetRawPedestalsInput(Char_t *name)
{
  //set raw pedestals
  if(fInRawPed){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawPedestalsInput","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  //Open the raw data file with name.
  fInRawPed = new ifstream();
#ifndef __DECCXX
  fInRawPed->open(name,fstream::binary);
#else
  fInRawPed->open(name);
#endif

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fInRawPed->rdbuf()->is_open()){
#else
  if(!fInRawPed->is_open()){
#endif
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawPedestalsInput","File Open")
      <<"Pointer to ifstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t AliHLTRawDataFileHandler::SetRawPedestalsInput(ifstream *file)
{
  //set raw pedestals input
  if(fInRawPed){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawPedestalsInput","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  //Open the raw data file with given file.
  fInRawPed = file;
#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fInRawPed->rdbuf()->is_open()){
#else
  if(!fInRawPed->is_open()){
#endif
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetRawPedestalsInput","File Open")
      <<"Pointer to ifstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

void AliHLTRawDataFileHandler::CloseRawPedestalsInput()
{
  //close raw pedestals input
  if(!fInRawPed){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::CloseRawPedestalsInput","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
#if defined(__HP_aCC) || defined(__DECCXX)
  if(fInRawPed->rdbuf()->is_open()) fInRawPed->close();
#else
  if(fInRawPed->is_open()) fInRawPed->close();
#endif
  delete fInRawPed;
  fInRaw = 0;
}

Bool_t AliHLTRawDataFileHandler::SetMappingFile(Char_t *name)
{
  //set mapping file
  if(fMapping){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetMapping","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  fMapping = fopen(name,"r");
  if(!fMapping){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::SetMappingFile","File Open")
      <<"Pointer to file = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t AliHLTRawDataFileHandler::SetMappingFile(FILE *file)
{
  //set mapping file
  if(fMapping){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::SetMapping","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  fMapping = file;
  if(!fMapping){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::SetMappingFile","File Open")
      <<"Pointer to file = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

void AliHLTRawDataFileHandler::CloseMappingFile()
{
  //close mapping file
  if(!fMapping){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::CloseMappingFile","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fMapping);
  fMapping = 0;
}

Int_t AliHLTRawDataFileHandler::ReadMappingFile()
{
  //read mapping file
  if(!fMapping){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadMappingFile","File Open")
      <<"Pointer to file = 0x0"<<ENDLOG;
    return -1;
  }

  Char_t dummy[100];
  fgets(dummy,80,fMapping);

  Int_t nboard,nadc;
  fscanf(fMapping,"%s %d %s %d",dummy,&nboard,dummy,&nadc); 

  fNChannels=nboard*nadc;
  if(fNChannels<=0){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::ReadMappingFile","Data Inconsistency")
      <<"fNChannels must be greater than 0"<<ENDLOG;
    return -1;
  }

  fRow=new Byte_t[fNChannels];
  fPad=new Byte_t[fNChannels];
  fRowPad=new Short_t*[AliHLTTransform::GetNRows()];
  for(Int_t i=0; i < AliHLTTransform::GetNRows(); i++){
    fRowPad[i]=new Short_t[AliHLTTransform::GetNPads(i)];
    for(Int_t j=0; j < AliHLTTransform::GetNPads(i); j++) fRowPad[i][j]=-1;
  }

  for(UInt_t i=0;i<fNChannels;i++){
    Int_t board,adc,row,pad;
    if(fscanf(fMapping,"%d %d %d %d",&board,&adc,&row,&pad)!=4) break; 
    //store the mapping
    fRow[i]=(Byte_t)row;
    fPad[i]=(Byte_t)pad;
    fRowPad[row][pad]=i;
    if(row>fRowMaxUsed) fRowMaxUsed=row;
    if(row<fRowMinUsed) fRowMinUsed=row;
    if(pad>fPadMaxUsed) fPadMaxUsed=pad;
    if(pad<fPadMinUsed) fPadMinUsed=pad;

    fNPads[row]++;
    //cout << i << " " << row << " " << pad << endl;
  }

  CloseMappingFile();
  return fNChannels;
}

inline Int_t AliHLTRawDataFileHandler::Convert4(Int_t i) const
{ 
  //BigEndian i0i1i2i3 -> LittleEndian i3i2i1i0
  if(!fConvert) return i;
  Char_t *p=(Char_t*)&i;
  Char_t temp[4];
  temp[0]=p[3];
  temp[1]=p[2];
  temp[2]=p[1];
  temp[3]=p[0];
  return (*(Int_t*)temp);
}

inline Short_t AliHLTRawDataFileHandler::Convert2(Short_t s) const
{ 
  //BigEndian i0i1 -> LittleEndian i1i0
  if(!fConvert) return s;
  Char_t *p=(Char_t*)&s;
  Char_t temp[2];
  temp[0]=p[1];
  temp[1]=p[0];
  return (*(Short_t*)temp);
}

Int_t AliHLTRawDataFileHandler::ReadRawInput()
{
  //Read data from cosmics file into memory
  if(!fInRaw){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawInput","File Open")
      <<"No Input avalible: no object ifstream"<<ENDLOG;
    return 0; 
  }

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fInRaw->rdbuf()->is_open()){
#else
  if(!fInRaw->is_open()){
#endif
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawInput","File Open")
      <<"No Input avalible: ifstream not opened"<<ENDLOG;
    return 0;
  }

  Int_t dummy4;
  Short_t dummy2;
  fInRaw->read((Char_t*)&dummy4,sizeof(dummy4));
  if(dummy4==(Int_t)fNChannels) fConvert=kFALSE;
  else {
    Int_t knumofChannels = Convert4(dummy4);    
    if(knumofChannels!=(Int_t)fNChannels){
      LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawInput","Data Inconsistency")
	<<"Number of Channels should be equal to fNChannels "<<knumofChannels<<" "<<fNChannels<<ENDLOG;
      return 0;
    }
  }
  
  //read used altrochannels (for the moment
  //this information is not really needed as
  //all channels per FEC are used
  for(UInt_t i = 0 ; i < fNChannels ; i++){
    fInRaw->read((Char_t*)&dummy2,sizeof(dummy2));
    UShort_t channel = Convert2(dummy2);
    if(channel>fNChannels){
      LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::ReadRawInput","Data Inconsistency")
	<<AliHLTLog::kDec<<"Channel number must be smaller then fNChannels "<<channel<<" "<<fNChannels<<ENDLOG;
      return 0;
    }
  }
  
   fInRaw->read((Char_t*)&dummy4,sizeof(dummy4));
   Int_t numofChannelsTest = Convert4(dummy4);

  if (numofChannelsTest != (Int_t)fNChannels){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::ReadRawInput","Data Inconsistency")
      <<AliHLTLog::kDec<<"Number of test channels should be equal to fNChannels "<<numofChannelsTest<<" "<<fNChannels<<ENDLOG;
    return 0;
  }

  //Timebins
  fInRaw->read((Char_t*)&dummy4,sizeof(dummy4));
  fNTimeBins=Convert4(dummy4);

  if(fNTimeBins!=AliHLTTransform::GetNTimeBins()){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawInput","Data Inconsistency")
      <<AliHLTLog::kDec<<"fNTimeBins does not match AliHLTTransformer, check AliHLTTransform::Init() "<<fNTimeBins<<" "<<AliHLTTransform::GetNTimeBins()<<ENDLOG;
  }

  //assign array
  if(fCharges) delete[] fCharges;
  fCharges=new Short_t*[fNChannels];
  for(UInt_t c=0;c<fNChannels;c++) fCharges[c]=new Short_t[fNTimeBins];

  //read data
  for(UInt_t channel = 0; channel < fNChannels; channel++){
    for(Int_t timebin = 0 ; timebin < fNTimeBins ; timebin++){
      Short_t dummy2;
      fInRaw->read((Char_t*)&dummy2,sizeof(dummy2));//1024012));
      Short_t charge = Convert2(dummy2);

      //Pedestal substraction
      if(fPedestals) charge-=fPedestals[channel][timebin];
      else charge-=fPedVal;
      if(charge<0) charge=0;

      fCharges[channel][timebin]=charge;
    }
  }
  return fNChannels;
}

Int_t AliHLTRawDataFileHandler::ReadRawInputPointer(const Char_t *ptr)
{
  //Read data from cosmics pointer into memory
  if(!ptr){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawInputPointer","Pointer")
      <<"Pointer equals 0x0!"<<ENDLOG;
    return 0; 
  }
  Int_t dummy4;
  Short_t dummy2;
  dummy4=*(Int_t*)ptr; ptr+=sizeof(dummy4);
  if(dummy4==(Int_t)fNChannels) fConvert=kFALSE;
  else {
    Int_t knumofChannels = Convert4(dummy4);    
    if(knumofChannels!=(Int_t)fNChannels){
      LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawInputPointer","Data Inconsistency")
	<<"Number of Channels should be equal to fNChannels "<<knumofChannels<<" "<<fNChannels<<ENDLOG;
      return 0;
    }
  }
  //read used altrochannels (for the moment
  //this information is not really needed as
  //all channels per FEC are used
  for(UInt_t i = 0 ; i < fNChannels ; i++){
    dummy2=*(Short_t*)ptr; ptr+=sizeof(dummy2);
    UShort_t channel = Convert2(dummy2);
    if(channel>fNChannels){
      LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::ReadRawInputPointer","Data Inconsistency")
	<<AliHLTLog::kDec<<"Channel number must be smaller then fNChannels "<<channel<<" "<<fNChannels<<ENDLOG;
      return 0;
    }
  }
  dummy4=*(Int_t*)ptr; ptr+=sizeof(dummy4);
  Int_t numofChannelsTest = Convert4(dummy4);
  if (numofChannelsTest != (Int_t)fNChannels){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::ReadRawInputPointer","Data Inconsistency")
      <<AliHLTLog::kDec<<"Number of test channels should be equal to fNChannels "<<numofChannelsTest<<" "<<fNChannels<<ENDLOG;
    return 0;
  }
  //Timebins
  dummy4=*(Int_t*)ptr; ptr+=sizeof(Int_t);
  fNTimeBins=Convert4(dummy4);
  if(fNTimeBins!=AliHLTTransform::GetNTimeBins()){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawInputPointer","Data Inconsistency")
      <<AliHLTLog::kDec<<"fNTimeBins does not match AliHLTTransformer, check AliHLTTransform::Init() "<<fNTimeBins<<" "<<AliHLTTransform::GetNTimeBins()<<ENDLOG;
  }
  //assign array
  if(fCharges) delete[] fCharges;
  fCharges=new Short_t*[fNChannels];
  for(UInt_t c=0;c<fNChannels;c++) fCharges[c]=new Short_t[fNTimeBins];
  //read data
  for(UInt_t channel = 0; channel < fNChannels; channel++){
    for(Int_t timebin = 0 ; timebin < fNTimeBins ; timebin++){
      Short_t dummy2=*(Short_t*)ptr;
      Short_t charge = Convert2(dummy2);

      //Pedestal substraction
      if(fPedestals) charge-=fPedestals[channel][timebin];
      else charge-=fPedVal;
      if(charge<0) charge=0;

      fCharges[channel][timebin]=charge;
      ptr+=sizeof(dummy2);
    }
  }
  
  return fNChannels;
}


Short_t** AliHLTRawDataFileHandler::GetRawData(Int_t &channels, Int_t &timebins)
{
  //get raw data
  Short_t **charges=0;
  channels=0;
  timebins=0;

  if(fNTimeBins==0){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::GetRawData","Data Inconsistency")
      <<"Call AliHLTRawDataFileHandler::RawReadInput() first"<<ENDLOG;
    if(!ReadRawInput()){
      LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::GetRawData","Data Inconsistency")
	<<"Something went wrong reading data header"<<ENDLOG;
      return 0;
    }
  }

  charges=new Short_t*[fNChannels];
  for(UInt_t c=0;c<fNChannels;c++) charges[c]=new Short_t[fNTimeBins];

  for(UInt_t channel = 0; channel < fNChannels; channel++){
    for(Int_t timebin = 0 ; timebin < fNTimeBins ; timebin++){
      Short_t dummy2;
      fInRaw->read((Char_t*)&dummy2,sizeof(dummy2));
      Short_t charge = Convert2(dummy2);
      charges[channel][timebin]=charge;
    }
  }

  channels=fNChannels;
  timebins=fNTimeBins;

  return charges;
}

Int_t AliHLTRawDataFileHandler::StoreRawData(Short_t **charges)
{
  //store charges in the raw data format
  if(!fOutRaw){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::StoreRawData","File Open")
      <<"No Output avalible: no object ofstream"<<ENDLOG;
    return 0; 
  }

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fOutRaw->rdbuf()->is_open()){
#else
  if(!fOutRaw->is_open()){
#endif
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::StoreRawData","File Open")
      <<"No Output avalible: ofstream not opened"<<ENDLOG;
    return 0;
  }

  Int_t dummy4;
  Short_t dummy2;

  dummy4=Convert4(fNChannels);

  fOutRaw->write((Char_t*)&dummy4,sizeof(dummy4));
  for(UInt_t i = 0 ; i < fNChannels ; i++){
    dummy2 = Convert2(Short_t(i));
    fOutRaw->write((Char_t*)&dummy2,sizeof(dummy2));
  }

  dummy4=Convert4(fNChannels);
  fOutRaw->write((Char_t*)&dummy4,sizeof(dummy4));

  //Timebins
  dummy4=Convert4(fNTimeBins);
  fOutRaw->write((Char_t*)&dummy4,sizeof(dummy4));

  for(UInt_t channel = 0; channel < fNChannels; channel++){
    for(Int_t timebin = 0 ; timebin < fNTimeBins ; timebin++){
      Short_t charge=charges[channel][timebin];
      dummy2 = Convert2(charge);
      fOutRaw->write((Char_t*)&dummy2,sizeof(dummy2));
    }
  }

  return fNChannels;
}

Int_t AliHLTRawDataFileHandler::ReadRawPedestalsInput()
{
  //read raw pedestals input
  if(!fInRawPed){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawPedestalsInput","File Open")
      <<"No Input avalible: no object ifstream"<<ENDLOG;
    return 0; 
  }

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fInRawPed->rdbuf()->is_open()){
#else
  if(!fInRawPed->is_open()){
#endif
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawPedestalsInput","File Open")
      <<"No Input avalible: ifstream not opened"<<ENDLOG;
    return 0;
  }

  Int_t dummy4;
  Short_t dummy2;
  fInRawPed->read((Char_t*)&dummy4,sizeof(dummy4));
  if(dummy4==(Int_t)fNChannels) fConvert=kFALSE;
  else {
    Int_t knumofChannels = Convert4(dummy4);    
    if(knumofChannels!=(Int_t)fNChannels){
      LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawPedestalsInput","Data Inconsistency")
	<<AliHLTLog::kDec<<"Number of Channels should be equal to fNChannels "<<knumofChannels<<" "<<fNChannels<<ENDLOG;
      return 0;
    }
  }
  
  //read used altrochannels (for the moment
  for(UInt_t i = 0 ; i < fNChannels ; i++){
    fInRawPed->read((Char_t*)&dummy2,sizeof(dummy2));
    //UShort_t channel = Convert2(dummy2);
  }

  fInRawPed->read((Char_t*)&dummy4,sizeof(dummy4));
  //Int_t numofChannelsTest = Convert4(dummy4);

  //Timebins
  fInRawPed->read((Char_t*)&dummy4,sizeof(dummy4));
  fNTimeBins=Convert4(dummy4);

  if(fNTimeBins!=AliHLTTransform::GetNTimeBins()){
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::ReadRawPedestalsInput","Data Inconsistency")
      <<AliHLTLog::kDec<<"fNTimeBins does not match AliHLTTransformer, check AliHLTTransform::Init() "<<fNTimeBins<<" "<<AliHLTTransform::GetNTimeBins()<<ENDLOG;
  }

  //Read the data
  fPedestals=new Short_t*[fNChannels];
  for(UInt_t c=0;c<fNChannels;c++) fPedestals[c]=new Short_t[fNTimeBins];

  for(UInt_t channel = 0; channel < fNChannels; channel++){
    for(Int_t timebin = 0 ; timebin < fNTimeBins ; timebin++){
      Short_t dummy2;
      fInRawPed->read((Char_t*)&dummy2,sizeof(dummy2));
      Short_t charge = Convert2(dummy2);
      fPedestals[channel][timebin]=charge;
    }
  }
  CloseRawPedestalsInput();
  return fNChannels;
}

AliHLTDigitRowData * AliHLTRawDataFileHandler::RawData2Memory(UInt_t &nrow,Int_t /*event*/)
{
  //convert raw data to memory
  AliHLTDigitRowData *data = 0;
  nrow=0;

  if(fNTimeBins==0){
    LOG(AliHLTLog::kWarning,"AliHLTRawDataFileHandler::RawData2Memory","Data Inconsistency")
      <<"Call AliHLTRawDataFileHandler::RawReadInput() first"<<ENDLOG;
    if(!ReadRawInput()){
      LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::RawData2Memory","Data Inconsistency")
	<<"Something went wrong reading data header"<<ENDLOG;
      return 0;
    }
  }

  //get data size
  Int_t nrows=0;
  Int_t ndigitcount=0;
  Int_t *ndigits=new Int_t[AliHLTTransform::GetNRows()];
  for(Int_t i=0;i<AliHLTTransform::GetNRows();i++) ndigits[i]=0;

  //no need to search for slice/sector given by init
  //but check for row/patch boundaries
  //assume slice 0
  for(Int_t slrow=0;slrow<AliHLTTransform::GetNRows();slrow++){
    
    if(slrow<fRowMin) continue;
    if(slrow>fRowMax) break;

    for(Int_t pad=0;pad<AliHLTTransform::GetNPads(slrow);pad++){
      Short_t channel=fRowPad[slrow][pad];
      if(channel==-1) continue; //no data on that channel;

      for(Int_t timebin = 0 ; timebin < fNTimeBins ; timebin++){
	Int_t dig=fCharges[channel][timebin];
	
	if(dig <= AliHLTTransform::GetZeroSup()) continue;
	if(dig >= AliHLTTransform::GetADCSat())
	  dig = AliHLTTransform::GetADCSat();

	ndigits[slrow]++; //for this row only
	ndigitcount++;  //total number of digits to be published
      }
    }
    
    //count number of rows
    nrows++;
  }

  //test data consistency
  Int_t ndigitcounttest=0;
  for(Int_t slrow=0;slrow<AliHLTTransform::GetNRows();slrow++)
    ndigitcounttest+=ndigits[slrow];
  if(ndigitcount!=ndigitcounttest)
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::RawData2Memory","Digits")
      <<AliHLTLog::kDec<<"Found Inconsistency "<<ndigitcount<<" != "<<ndigitcounttest<<ENDLOG;
    
  Int_t size = sizeof(AliHLTDigitData)*ndigitcount
    + nrows*sizeof(AliHLTDigitRowData);
  LOG(AliHLTLog::kDebug,"AliHLTRawDataFileHandler::RawData2Memory","Digits")
    <<AliHLTLog::kDec<<"Found "<<ndigitcount<<" Digits on "<<nrows<<" rows"<<ENDLOG;

  //now copy data
  data=(AliHLTDigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  //memset(data,1,size); //for debugging

  Int_t ndigitcounttest2=0;
  AliHLTDigitRowData *tempPt = data;
  for(Int_t slrow=0;slrow<AliHLTTransform::GetNRows();slrow++){
    
    if(slrow<fRowMin) continue;
    if(slrow>fRowMax) break;
    
    tempPt->fRow = slrow;
    tempPt->fNDigit = ndigits[slrow];

    Int_t localcount=0;
    for(Int_t pad=0;pad<AliHLTTransform::GetNPads(slrow);pad++){
      Short_t channel=fRowPad[slrow][pad];
      if(channel==-1) continue; //no data on that channel;

      for(Int_t timebin = 0 ; timebin < fNTimeBins ; timebin++){
	Int_t dig=fCharges[channel][timebin];
	
	if(dig <= AliHLTTransform::GetZeroSup()) continue;
	if(dig >= AliHLTTransform::GetADCSat())
	  dig = AliHLTTransform::GetADCSat();

	//Exclude data outside cone:
	//AliHLTTransform::Raw2Local(xyz,sector,row,pad,time);
	//if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2])) continue;

	tempPt->fDigitData[localcount].fCharge=(UShort_t)dig;
	tempPt->fDigitData[localcount].fPad=(UChar_t)pad;
	tempPt->fDigitData[localcount].fTime=(UShort_t)timebin;
#ifdef do_mc
	tempPt->fDigitData[localcount].fTrackID[0] = 0;
	tempPt->fDigitData[localcount].fTrackID[1] = 0;
	tempPt->fDigitData[localcount].fTrackID[2] = 0;
#endif
	localcount++;
	ndigitcounttest2++;
      } //time
    } //pad 

    if(localcount != ndigits[slrow])
      LOG(AliHLTLog::kFatal,"AliHLTRawDataFileHandler::RawData2Memory","Memory")
	<<AliHLTLog::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
	<<ndigits[slrow]<<ENDLOG;

    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t size = sizeof(AliHLTDigitRowData)
      + ndigits[slrow]*sizeof(AliHLTDigitData);
    tmp += size;
    tempPt = (AliHLTDigitRowData*)tmp;
  }//row

  if(ndigitcount!=ndigitcounttest2)
    LOG(AliHLTLog::kError,"AliHLTRawDataFileHandler::RawData2Memory","Digits")
      <<AliHLTLog::kDec<<"Found Inconsistency "<<ndigitcount<<" != "<<ndigitcounttest2<<ENDLOG;

  delete [] ndigits;
  return data;
}

Bool_t AliHLTRawDataFileHandler::RawData2CompBinary(Int_t event)
{
  //raw data to binary
  Bool_t out = kTRUE;
  UInt_t ndigits=0;
  AliHLTDigitRowData *digits=0;
  digits = RawData2Memory(ndigits,event);
  out = Memory2CompBinary(ndigits,digits);
  Free();
  return out;
}
