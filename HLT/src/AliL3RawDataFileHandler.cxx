// @(#) $Id$

// Author: C. Loizides <loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3RootTypes.h"
#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3MemHandler.h"
#include "AliL3DigitData.h"

#include "AliL3RawDataFileHandler.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** \class AliL3RawDataFileHandler 
<pre>
//_____________________________________________________________
// AliL3RawDataFileHandler
//
</pre>
*/

ClassImp(AliL3RawDataFileHandler)

AliL3RawDataFileHandler::AliL3RawDataFileHandler()
{
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
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::AliL3RawDataFileHandler","Constructor")
      <<"Check architecture to run the conversion on! Int_t should be 32 and Short_t should be 16 bit."<<ENDLOG;
  }
}

AliL3RawDataFileHandler::~AliL3RawDataFileHandler()
{
  FreeAll();
}

void AliL3RawDataFileHandler::FreeAll()
{
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
  fRowMinUsed=AliL3Transform::GetNRows();
  fRowMaxUsed=0;
  fPadMinUsed=255;
  fPadMaxUsed=0;
  fNTimeBins=0;
  for(Int_t i=0;i<AliL3Transform::GetNRows();i++) fNPads[i]=0;
  fPedVal=0;
}

Bool_t AliL3RawDataFileHandler::SetRawInput(Char_t *name)
{
  if(fInRaw){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawInput","File Open")
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
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawInput","File Open")
      <<"Pointer to ifstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t AliL3RawDataFileHandler::SetRawInput(ifstream *file)
{
  if(fInRaw){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawInput","File Open")
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
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawInput","File Open")
      <<"Pointer to ifstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

void AliL3RawDataFileHandler::CloseRawInput()
{
  if(!fInRaw){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::CloseRawInput","File Close")
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

Bool_t AliL3RawDataFileHandler::SetRawOutput(Char_t *name)
{
  if(fOutRaw){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawOutput","File Open")
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
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawOutput","File Open")
      <<"Pointer to ofstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t AliL3RawDataFileHandler::SetRawOutput(ofstream *file)
{
  if(fOutRaw){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawOutput","File Open")
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
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawOutput","File Open")
      <<"Pointer to ofstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

void AliL3RawDataFileHandler::CloseRawOutput()
{
  if(!fOutRaw){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::CloseRawOutput","File Close")
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


Bool_t AliL3RawDataFileHandler::SetRawPedestalsInput(Char_t *name)
{
  if(fInRawPed){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawPedestalsInput","File Open")
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
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawPedestalsInput","File Open")
      <<"Pointer to ifstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t AliL3RawDataFileHandler::SetRawPedestalsInput(ifstream *file)
{
  if(fInRawPed){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawPedestalsInput","File Open")
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
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetRawPedestalsInput","File Open")
      <<"Pointer to ifstream = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

void AliL3RawDataFileHandler::CloseRawPedestalsInput()
{
  if(!fInRawPed){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::CloseRawPedestalsInput","File Close")
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

Bool_t AliL3RawDataFileHandler::SetMappingFile(Char_t *name)
{
  if(fMapping){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetMapping","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  fMapping = fopen(name,"r");
  if(!fMapping){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::SetMappingFile","File Open")
      <<"Pointer to file = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t AliL3RawDataFileHandler::SetMappingFile(FILE *file)
{
  if(fMapping){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::SetMapping","File Open")
      <<"File ptr is already in use, close file first"<<ENDLOG;
    return kFALSE;
  }

  fMapping = file;
  if(!fMapping){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::SetMappingFile","File Open")
      <<"Pointer to file = 0x0"<<ENDLOG;
    return kFALSE;
  }
  
  return kTRUE;
}

void AliL3RawDataFileHandler::CloseMappingFile()
{
  if(!fMapping){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::CloseMappingFile","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fMapping);
  fMapping = 0;
}

Int_t AliL3RawDataFileHandler::ReadMappingFile()
{
  if(!fMapping){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadMappingFile","File Open")
      <<"Pointer to file = 0x0"<<ENDLOG;
    return -1;
  }

  Char_t dummy[100];
  fgets(dummy,80,fMapping);

  Int_t nboard,nadc;
  fscanf(fMapping,"%s %d %s %d",dummy,&nboard,dummy,&nadc); 

  fNChannels=nboard*nadc;
  if(fNChannels<=0){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::ReadMappingFile","Data Inconsistency")
      <<"fNChannels must be greater than 0"<<ENDLOG;
    return -1;
  }

  fRow=new Byte_t[fNChannels];
  fPad=new Byte_t[fNChannels];
  fRowPad=new Short_t*[AliL3Transform::GetNRows()];
  for(Int_t i=0; i < AliL3Transform::GetNRows(); i++){
    fRowPad[i]=new Short_t[AliL3Transform::GetNPads(i)];
    for(Int_t j=0; j < AliL3Transform::GetNPads(i); j++) fRowPad[i][j]=-1;
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

inline Int_t AliL3RawDataFileHandler::Convert4(Int_t i)
{ //BigEndian i0i1i2i3 -> LittleEndian i3i2i1i0
  if(!fConvert) return i;
  Char_t *p=(Char_t*)&i;
  Char_t temp[4];
  temp[0]=p[3];
  temp[1]=p[2];
  temp[2]=p[1];
  temp[3]=p[0];
  return (*(Int_t*)temp);
}

inline Short_t AliL3RawDataFileHandler::Convert2(Short_t s)
{ //BigEndian i0i1 -> LittleEndian i1i0
  if(!fConvert) return s;
  Char_t *p=(Char_t*)&s;
  Char_t temp[2];
  temp[0]=p[1];
  temp[1]=p[0];
  return (*(Short_t*)temp);
}

Int_t AliL3RawDataFileHandler::ReadRawInput()
{
  //Read data from cosmics file into memory
  if(!fInRaw){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawInput","File Open")
      <<"No Input avalible: no object ifstream"<<ENDLOG;
    return 0; 
  }

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fInRaw->rdbuf()->is_open()){
#else
  if(!fInRaw->is_open()){
#endif
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawInput","File Open")
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
      LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawInput","Data Inconsistency")
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
      LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::ReadRawInput","Data Inconsistency")
	<<AliL3Log::kDec<<"Channel number must be smaller then fNChannels "<<channel<<" "<<fNChannels<<ENDLOG;
      return 0;
    }
  }
  
   fInRaw->read((Char_t*)&dummy4,sizeof(dummy4));
   Int_t numofChannelsTest = Convert4(dummy4);

  if (numofChannelsTest != (Int_t)fNChannels){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::ReadRawInput","Data Inconsistency")
      <<AliL3Log::kDec<<"Number of test channels should be equal to fNChannels "<<numofChannelsTest<<" "<<fNChannels<<ENDLOG;
    return 0;
  }

  //Timebins
  fInRaw->read((Char_t*)&dummy4,sizeof(dummy4));
  fNTimeBins=Convert4(dummy4);

  if(fNTimeBins!=AliL3Transform::GetNTimeBins()){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawInput","Data Inconsistency")
      <<AliL3Log::kDec<<"fNTimeBins does not match AliL3Transformer, check AliL3Transform::Init() "<<fNTimeBins<<" "<<AliL3Transform::GetNTimeBins()<<ENDLOG;
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

Int_t AliL3RawDataFileHandler::ReadRawInputPointer(const Char_t *ptr)
{
  //Read data from cosmics pointer into memory
  if(!ptr){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawInputPointer","Pointer")
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
      LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawInputPointer","Data Inconsistency")
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
      LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::ReadRawInputPointer","Data Inconsistency")
	<<AliL3Log::kDec<<"Channel number must be smaller then fNChannels "<<channel<<" "<<fNChannels<<ENDLOG;
      return 0;
    }
  }
  dummy4=*(Int_t*)ptr; ptr+=sizeof(dummy4);
  Int_t numofChannelsTest = Convert4(dummy4);
  if (numofChannelsTest != (Int_t)fNChannels){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::ReadRawInputPointer","Data Inconsistency")
      <<AliL3Log::kDec<<"Number of test channels should be equal to fNChannels "<<numofChannelsTest<<" "<<fNChannels<<ENDLOG;
    return 0;
  }
  //Timebins
  dummy4=*(Int_t*)ptr; ptr+=sizeof(Int_t);
  fNTimeBins=Convert4(dummy4);
  if(fNTimeBins!=AliL3Transform::GetNTimeBins()){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawInputPointer","Data Inconsistency")
      <<AliL3Log::kDec<<"fNTimeBins does not match AliL3Transformer, check AliL3Transform::Init() "<<fNTimeBins<<" "<<AliL3Transform::GetNTimeBins()<<ENDLOG;
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


Short_t** AliL3RawDataFileHandler::GetRawData(Int_t &channels, Int_t &timebins)
{
  Short_t **charges=0;
  channels=0;
  timebins=0;

  if(fNTimeBins==0){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::GetRawData","Data Inconsistency")
      <<"Call AliL3RawDataFileHandler::RawReadInput() first"<<ENDLOG;
    if(!ReadRawInput()){
      LOG(AliL3Log::kError,"AliL3RawDataFileHandler::GetRawData","Data Inconsistency")
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

Int_t AliL3RawDataFileHandler::StoreRawData(Short_t **charges)
{
  //store charges in the raw data format

  if(!fOutRaw){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::StoreRawData","File Open")
      <<"No Output avalible: no object ofstream"<<ENDLOG;
    return 0; 
  }

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fOutRaw->rdbuf()->is_open()){
#else
  if(!fOutRaw->is_open()){
#endif
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::StoreRawData","File Open")
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

Int_t AliL3RawDataFileHandler::ReadRawPedestalsInput()
{
  if(!fInRawPed){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawPedestalsInput","File Open")
      <<"No Input avalible: no object ifstream"<<ENDLOG;
    return 0; 
  }

#if defined(__HP_aCC) || defined(__DECCXX)
  if(!fInRawPed->rdbuf()->is_open()){
#else
  if(!fInRawPed->is_open()){
#endif
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawPedestalsInput","File Open")
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
      LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawPedestalsInput","Data Inconsistency")
	<<AliL3Log::kDec<<"Number of Channels should be equal to fNChannels "<<knumofChannels<<" "<<fNChannels<<ENDLOG;
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

  if(fNTimeBins!=AliL3Transform::GetNTimeBins()){
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::ReadRawPedestalsInput","Data Inconsistency")
      <<AliL3Log::kDec<<"fNTimeBins does not match AliL3Transformer, check AliL3Transform::Init() "<<fNTimeBins<<" "<<AliL3Transform::GetNTimeBins()<<ENDLOG;
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

AliL3DigitRowData * AliL3RawDataFileHandler::RawData2Memory(UInt_t &nrow,Int_t /*event*/)
{
  AliL3DigitRowData *data = 0;
  nrow=0;

  if(fNTimeBins==0){
    LOG(AliL3Log::kWarning,"AliL3RawDataFileHandler::RawData2Memory","Data Inconsistency")
      <<"Call AliL3RawDataFileHandler::RawReadInput() first"<<ENDLOG;
    if(!ReadRawInput()){
      LOG(AliL3Log::kError,"AliL3RawDataFileHandler::RawData2Memory","Data Inconsistency")
	<<"Something went wrong reading data header"<<ENDLOG;
      return 0;
    }
  }
  

  //get data size
  Int_t nrows=0;
  Int_t ndigitcount=0;
  Int_t *ndigits=new Int_t[AliL3Transform::GetNRows()];
  for(Int_t i=0;i<AliL3Transform::GetNRows();i++) ndigits[i]=0;

  //no need to search for slice/sector given by init
  //but check for row/patch boundaries
  //assume slice 0
  for(Int_t slrow=0;slrow<AliL3Transform::GetNRows();slrow++){
    
    if(slrow<fRowMin) continue;
    if(slrow>fRowMax) break;

    for(Int_t pad=0;pad<AliL3Transform::GetNPads(slrow);pad++){
      Short_t channel=fRowPad[slrow][pad];
      if(channel==-1) continue; //no data on that channel;

      for(Int_t timebin = 0 ; timebin < fNTimeBins ; timebin++){
	Int_t dig=fCharges[channel][timebin];
	
	if(dig <= AliL3Transform::GetZeroSup()) continue;
	if(dig >= AliL3Transform::GetADCSat())
	  dig = AliL3Transform::GetADCSat();

	ndigits[slrow]++; //for this row only
	ndigitcount++;  //total number of digits to be published
      }
    }
    
    //count number of rows
    nrows++;
  }

  //test data consistency
  Int_t ndigitcounttest=0;
  for(Int_t slrow=0;slrow<AliL3Transform::GetNRows();slrow++)
    ndigitcounttest+=ndigits[slrow];
  if(ndigitcount!=ndigitcounttest)
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::RawData2Memory","Digits")
      <<AliL3Log::kDec<<"Found Inconsistency "<<ndigitcount<<" != "<<ndigitcounttest<<ENDLOG;
    
  Int_t size = sizeof(AliL3DigitData)*ndigitcount
    + nrows*sizeof(AliL3DigitRowData);
  LOG(AliL3Log::kDebug,"AliL3RawDataFileHandler::RawData2Memory","Digits")
    <<AliL3Log::kDec<<"Found "<<ndigitcount<<" Digits on "<<nrows<<" rows"<<ENDLOG;

  //now copy data
  data=(AliL3DigitRowData*) Allocate(size);
  nrow = (UInt_t)nrows;
  //memset(data,1,size); //for debugging

  Int_t ndigitcounttest2=0;
  AliL3DigitRowData *tempPt = data;
  for(Int_t slrow=0;slrow<AliL3Transform::GetNRows();slrow++){
    
    if(slrow<fRowMin) continue;
    if(slrow>fRowMax) break;
    
    tempPt->fRow = slrow;
    tempPt->fNDigit = ndigits[slrow];

    Int_t localcount=0;
    for(Int_t pad=0;pad<AliL3Transform::GetNPads(slrow);pad++){
      Short_t channel=fRowPad[slrow][pad];
      if(channel==-1) continue; //no data on that channel;

      for(Int_t timebin = 0 ; timebin < fNTimeBins ; timebin++){
	Int_t dig=fCharges[channel][timebin];
	
	if(dig <= AliL3Transform::GetZeroSup()) continue;
	if(dig >= AliL3Transform::GetADCSat())
	  dig = AliL3Transform::GetADCSat();

	//Exclude data outside cone:
	//AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
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
      LOG(AliL3Log::kFatal,"AliL3RawDataFileHandler::RawData2Memory","Memory")
	<<AliL3Log::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
	<<ndigits[slrow]<<ENDLOG;

    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t size = sizeof(AliL3DigitRowData)
      + ndigits[slrow]*sizeof(AliL3DigitData);
    tmp += size;
    tempPt = (AliL3DigitRowData*)tmp;
  }//row

  if(ndigitcount!=ndigitcounttest2)
    LOG(AliL3Log::kError,"AliL3RawDataFileHandler::RawData2Memory","Digits")
      <<AliL3Log::kDec<<"Found Inconsistency "<<ndigitcount<<" != "<<ndigitcounttest2<<ENDLOG;

  delete [] ndigits;
  return data;
}

Bool_t AliL3RawDataFileHandler::RawData2CompBinary(Int_t event)
{
  Bool_t out = kTRUE;
  UInt_t ndigits=0;
  AliL3DigitRowData *digits=0;
  digits = RawData2Memory(ndigits,event);
  out = Memory2CompBinary(ndigits,digits);
  Free();
  return out;
}
