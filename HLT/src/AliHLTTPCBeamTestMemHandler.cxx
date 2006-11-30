// @(#) $Id$

/** \class AliHLTTPCBeamTestMemHandler 
<pre>
//_____________________________________________________________
// AliHLTTPCBeamTestMemHandler
//
// Class for converting the test beam data of May 2004 
// to the HLT file format using R. Bramms tables.
//
// Author: C. Loizides <loizides@ikf.uni-frankfurt.de>
// -- Copyright &copy ALICE HLT Group
</pre>
*/

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTTransform.h"
#include "AliHLTMemHandler.h"
#include "AliHLTDigitData.h"
#include "AliHLTTPCBeamTestMemHandler.h"

#if __GNUC__ >= 3
using namespace std;
#endif


ClassImp(AliHLTTPCBeamTestMemHandler)

AliHLTTPCBeamTestMemHandler::AliHLTTPCBeamTestMemHandler(Char_t *fPathToMappingFile) : AliHLTMemHandler()
{ 
  //constructor
  fMinTimeBin=1;
  fNumOfChannels=7807+1; //must be big enough to contain all channels (per patch)

  Char_t readcarry[255];
  Int_t actPos=0;
  Int_t oldPos=0;
  ifstream *in = new ifstream();
  in->open(fPathToMappingFile); 
#if defined(__HP_aCC) || defined(__DECCXX)
  if(!in->rdbuf()->is_open()){
#else
  if(!in->is_open()){
#endif
    LOG(AliHLTLog::kFatal,"AliHLTTPCBeamTestMemHandler","Mapping File")
	<<"Can't open file " << fPathToMappingFile << " !!!" <<ENDLOG;
  }
  fMapping = new short*[fNumOfChannels];
  for ( Int_t i = 0; i < fNumOfChannels; i++ )
      fMapping[i] = 0;
  fMappingEmptyRow = new short[11]; //11 colums per row in mapping file
  for(Int_t i = 0; i < 11 ; i++){
      fMappingEmptyRow[i] = 0;
  }
  Short_t *mappingRow;
  for(Int_t i = 0; i < 5504 ; i++) { //5504 is size of irorc mapping at the moment only for irorc
      mappingRow = new Short_t[11];
      for(Int_t j = 0 ; j < 11 ; j++) {
	  *in >> readcarry;
	  mappingRow[j] = atoi(readcarry);
     }
      actPos = mappingRow[0];
      fMapping[actPos] = mappingRow;
      if( (actPos - oldPos) > 1){
	  for(Int_t j = (oldPos+1); j < actPos; j++){
	      fMapping[j] = fMappingEmptyRow;
	  }
      }
      oldPos = actPos;
  }
  in->close();
  delete in;
}

AliHLTTPCBeamTestMemHandler::~AliHLTTPCBeamTestMemHandler()
{ 
  //destructor
  for(Int_t i = 0; i < 5504 ; i++) { 
	if(fMapping[i] != fMappingEmptyRow && fMapping[i]) delete[] fMapping[i];
  }
  delete [] fMappingEmptyRow;
 delete[] fMapping;
}

AliHLTDigitRowData* AliHLTTPCBeamTestMemHandler::RawData2Memory(UInt_t &nrow,Int_t /*event*/)
{ 
  //convert the raw data
  AliHLTDigitRowData *data = 0;
  nrow=0;

  Int_t nrowsdummy=AliHLTTransform::GetNRows(fPatch);
  fRows = new AliRowStructure[nrowsdummy];
  for(Int_t i=0;i<nrowsdummy;i++){
   fRows[i].fRow=-1;
   fRows[i].fNDigits=0;
   fRows[i].fPadPos= new Int_t[AliHLTTransform::GetNPads(i+fRowMin)];
   for(Int_t p=0;p<AliHLTTransform::GetNPads(i+fRowMin);p++)
     fRows[i].fPadPos[p]=-1;
  }

  Int_t ntimebins=AliHLTTransform::GetNTimeBins();
  Int_t npads=fInputSize/(ntimebins+1);
  Int_t ndigitcount=0; //total number of digits to be published
  for(Int_t i=0;i<npads;i++){
    Int_t pos=i*(ntimebins+1);
    Short_t hw=fInputPtr[pos];
    if(hw>=fNumOfChannels) continue;
    Int_t pad=MappingGetPad(hw);
    Int_t lrow=MappingGetPadRow(hw);
    if((lrow<0) ||(pad<0)) continue;
    if(lrow+fRowMin>fRowMax) continue;
    fRows[lrow].fRow=lrow+fRowMin;
    if(fRows[lrow].fPadPos[pad]!=-1){
      continue;
    }
    Bool_t isThereDataOnThisPad=kFALSE;
    Int_t digmean=0;
#if 1
    for(Int_t timebin = fMinTimeBin ; timebin <= ntimebins ; timebin++){
      Int_t dig=fInputPtr[pos+timebin];
      digmean+=dig;
    }
    digmean/=(ntimebins-fMinTimeBin+1);
#else
    digmean = 40;
#endif
    for(Int_t timebin = fMinTimeBin ; timebin <= ntimebins ; timebin++){
      Int_t dig=fInputPtr[pos+timebin]-digmean;
    
      if(dig <= AliHLTTransform::GetZeroSup()) continue;
      if(dig >= AliHLTTransform::GetADCSat())
        dig = AliHLTTransform::GetADCSat();

	fRows[lrow].fNDigits++; //for this row only
	ndigitcount++;  
	isThereDataOnThisPad=kTRUE;
      }

    if(isThereDataOnThisPad) {
      fRows[lrow].fPadPos[pad]=pos;
    }
  }

  Int_t nrows=0;
  for(Int_t i=0;i<AliHLTTransform::GetNRows(fPatch);i++){
      if(fRows[i].fRow!=-1) nrows++;
  }
  if(nrows!=AliHLTTransform::GetNRows(fPatch))
    LOG(AliHLTLog::kError,"AliHLTTPCBeamTestMemHandler::RawData2Memory","nrows")
      <<AliHLTLog::kDec<<"Found Inconsistency "<<nrows<<" != "<<AliHLTTransform::GetNRows(fPatch)<<ENDLOG;

  //allocate memory
  Int_t size = sizeof(AliHLTDigitData)*ndigitcount
    + nrows*sizeof(AliHLTDigitRowData);
  LOG(AliHLTLog::kDebug,"AliHLTTPCBeamTestMemHandler::RawData2Memory","Digits")
    <<AliHLTLog::kDec<<"Found "<<ndigitcount<<" Digits on "<<nrows<<" rows"<<ENDLOG;

  data=(AliHLTDigitRowData*)Allocate(size);
  nrow = (UInt_t)nrows;
  //memset(data,1,size); //for debugging

  Int_t ndigitcounttest=0;
  AliHLTDigitRowData *tempPt = data;
  for(Int_t i=0;i<AliHLTTransform::GetNRows(fPatch);i++){
    Int_t slrow=i+fRowMin;
    
    if(slrow!=fRows[i].fRow){
      LOG(AliHLTLog::kFatal,"AliHLTTPCBeamTestMemHandler::RawData2Memory","Row Mismatch")
	<<AliHLTLog::kDec<<"Mismatch: slrow "<<slrow<<" row "
	<<fRows[i].fRow<<ENDLOG;
    }

    tempPt->fRow = slrow;
    tempPt->fNDigit = fRows[i].fNDigits;

    Int_t localcount=0;
    for(Int_t pad=0;pad<AliHLTTransform::GetNPads(slrow);pad++){
      Int_t pos=fRows[i].fPadPos[pad];
      if(pos==-1) continue; //no data on that pad;
      Int_t digmean=0;
#if 1
      for(Int_t timebin = fMinTimeBin ; timebin <= ntimebins ; timebin++){
        Int_t dig=fInputPtr[pos+timebin];
	digmean+=dig;
      }
      digmean/=(ntimebins-fMinTimeBin+1);
#else
    digmean = 40;
#endif
      for(Int_t timebin = fMinTimeBin ; timebin <= ntimebins ; timebin++){
        Int_t dig=fInputPtr[pos+timebin]-digmean;
    
	if(dig <= AliHLTTransform::GetZeroSup()) continue;
	if(dig >= AliHLTTransform::GetADCSat())
	    dig = AliHLTTransform::GetADCSat();

	//Exclude data outside cone:
	//AliHLTTransform::Raw2Local(xyz,sector,row,pad,time);
	//if(fParam->GetPadRowRadii(sector,row)<230./250.*fabs(xyz[2])) continue;

	tempPt->fDigitData[localcount].fCharge=(UShort_t)dig;
	tempPt->fDigitData[localcount].fPad=(UChar_t)pad;
	tempPt->fDigitData[localcount].fTime=(UShort_t)timebin-1;
	//cout << slrow << " " << pad << " " << timebin << " " << dig << endl;
#ifdef do_mc
	tempPt->fDigitData[localcount].fTrackID[0] = 0;
	tempPt->fDigitData[localcount].fTrackID[1] = 0;
	tempPt->fDigitData[localcount].fTrackID[2] = 0;
#endif
	localcount++;
	ndigitcounttest++;
      } //time
    } //pad 

    if(localcount != fRows[i].fNDigits)
      LOG(AliHLTLog::kFatal,"AliHLTTPCBeamTestMemHandler::RawData2Memory","Memory")
	<<AliHLTLog::kDec<<"Mismatch: localcount "<<localcount<<" ndigits "
	<<fRows[i].fNDigits<<ENDLOG;

    Byte_t *tmp = (Byte_t*)tempPt;
    Int_t size = sizeof(AliHLTDigitRowData)
      + localcount*sizeof(AliHLTDigitData);
    tmp += size;
    tempPt = (AliHLTDigitRowData*)tmp;
  }//row

  if(ndigitcount!=ndigitcounttest)
    LOG(AliHLTLog::kError,"AliHLTTPCBeamTestMemHandler::RawData2Memory","Digits")
      <<AliHLTLog::kDec<<"Found Inconsistency "<<ndigitcount<<" != "<<ndigitcounttest<<ENDLOG;

  for(Int_t i=0;i<nrowsdummy;i++){
    delete[] fRows[i].fPadPos;
  }
  delete[] fRows;
  return data;
}

Bool_t AliHLTTPCBeamTestMemHandler::RawData2CompBinary(Int_t event)
{ 
  //raw data to memory
  Bool_t out = kTRUE;
  UInt_t ndigits=0;
  AliHLTDigitRowData *digits=0;
  digits = RawData2Memory(ndigits,event);
  out = Memory2CompBinary(ndigits,digits);
  Free();
  return out;
}

