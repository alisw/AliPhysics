#include "AliMUONTrackerDDLDecoder.h"
#include "AliMUONTrackerDDLDecoderEventHandler.h"

#include "hlt_structs.C"


class DiMuonTrackerCustomHandler : public AliMUONTrackerDDLDecoderEventHandler
{
public:

//   void OnData(UInt_t data) // Old type
  void OnNewBusPatch(const AliMUONBusPatchHeaderStruct* header, const void* data)
  {
    fBusPatchId = header->fBusPatchId;
  }

  void OnData(const UInt_t data, bool isTrue)
  {

    fData[fDataCount].fDataId &= 0x0;
    fData[fDataCount].fDataId = (fData[fDataCount].fDataId|fBusPatchId)<<17;  
    fData[fDataCount].fDataId |= ((UInt_t)(data>>12) & 0x1FFFF) ;

    fData[fDataCount].fADC = (UShort_t)(data) & 0xFFF;
    
    fDataCount++;
    if(fDataCount>(27000-1))
      cerr<<"Running Out of memory "<<endl;

  };
  
  void OnError(ErrorCode code, const void* location)
  {
    
    printf("Errorcode : %d, %s, %s\n",
	   code,ErrorCodeToString(code),ErrorCodeToMessage(code));
  };
  
  void ResetDataCounter(){
    fDataCount = 0;
  }

  Int_t GetDataSize(){return fDataCount;};

  AliHLTMUONTrackerRawData GetData(Int_t iData){return fData[iData] ;};

private:
  Int_t fDataCount;
  Int_t fMaxDataCount;
  UInt_t fBusPatchId;
  AliHLTMUONTrackerRawData fData[27000];
};

class DiMuonTriggerDDLDecoder : public TObject
{
public : 

  DiMuonTriggerDDLDecoder(){;}
  virtual ~DiMuonTriggerDDLDecoder(void){;}
  void SetTriggerMappingData(AliHLTMUONTriggerMappingData* mappingData){ fMapData = mappingData;}
  bool Decode(int* rawData);

  void ResetDataCounter(){
    fDataCount = 0;
  }
  Int_t GetDataSize(){return fDataCount;}

  AliHLTMUONTriggerPointData  GetData(Int_t iData){return fData[iData] ;}
  
  
private:
  AliHLTMUONTriggerMappingData* fMapData;
  Int_t fDataCount;
  Int_t fMaxDataCount;
  AliHLTMUONTriggerPointData fData[27000];
};

bool DiMuonTriggerDDLDecoder::Decode(int* rawData)
{

  int nofTrigRec = 0;
  int index = 0;
  int reg_output, reg_phys_trig_occur;
  int iLocIndex,loc,locDec,triggY,sign,loDev,triggX;
  short pattern[2][4]; // 2 stands for two cathode planes and 4 stands for 4 chambers
	
  int phys_trig_occur = (rawData[index]>>30)&0x1; // 1 for physics trigger, 0 for software trigger
	
  if (not phys_trig_occur) // for software trigger
    index += 8 ;// corresponding to scalar words
	
  index += 1 ; // To skip the separator 0xDEADFACE
  
  index += 4 ; // corresponding to global input
  
  index += 1 ; // reaches to global output
  
  if (not phys_trig_occur) index += 10; // corresponds to scalar words
	
  index += 1; // separator 0xDEADBEEF 
	
  for (int iReg = 0; iReg < 8; iReg++){

    index += 1; // DARC Status Word
    index += 1; // Regeional Word
    reg_output = rawData[index] & 0xFF;
    reg_phys_trig_occur = ( rawData[index] >> 31) & 0x1;
//     cout<<" reg_phys_trig_occur : "<<reg_phys_trig_occur<<endl;
    index += 2; // 2 words for regional input
    
    index += 1; // L0 counter
    
    if (not reg_phys_trig_occur) index += 10;
    
    index += 1; // end of Regeonal header 0xBEEFFACE
    
    for(int iLoc = 0; iLoc < 16 ; iLoc++){
      
      iLocIndex = index;
      
      loc = (rawData[index+5] >> 19) &  0xF ;

      locDec = (rawData[index+5] >> 15) & 0xF;
      triggY = (rawData[index+5] >> 14) & 0x1;
      sign = (rawData[index+5] >> 9) & 0x1;
      loDev = (rawData[index+5] >> 5) & 0xF ;
      triggX = (loDev >> 4 & 0x1 ) && !(loDev & 0xF);
      
      if( locDec != 0x9 ){ // check for Dec
	  
	index += 1;
	pattern[0][0] = rawData[index] & 0xFFFF; // x-strip pattern for chamber 0 
	pattern[0][1] = (rawData[index] >> 16) & 0xFFFF; // x-strip pattern for chamber 1
	index += 1; 
	pattern[0][2] = rawData[index] & 0xFFFF; 
	pattern[0][3] = (rawData[index] >> 16) & 0xFFFF; 
	
	index += 1;
	pattern[1][0] = rawData[index] & 0xFFFF; // y-strip pattern for chamber 0
	pattern[1][1] = (rawData[index] >> 16) & 0xFFFF; // y-strip pattern for chamber 0 
	index += 1; 
	pattern[1][2] = rawData[index] & 0xFFFF; 
	pattern[1][3] = (rawData[index] >> 16) & 0xFFFF; 
	

	if (pattern[0][0] || pattern[0][1] || pattern[0][2] || pattern[0][3]
	    || pattern[1][0] || pattern[1][1] || pattern[1][2] || pattern[1][3]
	    ){
	      

	  bool Xset[4] = {false, false, false, false};
	  bool Yset[4] = {false, false, false, false};
	  AliHLTMUONTriggerPointData fHit[4];

	  for (int iChamber = 0; iChamber < 4 ; iChamber++){ //4 chambers per DDL 
	    for (int iPlane = 0; iPlane < 2 ; iPlane++){ // 2 cathode plane
	      for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy){
		
		if (((pattern[iPlane][iChamber] >> ibitxy) & 0x1) != 0x1)
		  continue;
			
		if (iPlane == 1){
		  if((fMapData->fLut[iReg][iLoc][iChamber][iPlane][ibitxy]).fX != 0){

		    fHit[iChamber].fX = (fMapData->fLut[iReg][iLoc][iChamber][iPlane][ibitxy]).fX ;
		    fHit[iChamber].fDetElemId = (fMapData->fLut[iReg][iLoc][iChamber][iPlane][ibitxy]).fDetElemId ;
		    
		    Xset[iChamber] = true ;
		  }
		}
		else{
		  if((fMapData->fLut[iReg][iLoc][iChamber][iPlane][ibitxy]).fY != 0){
		    
		    fHit[iChamber].fY = (fMapData->fLut[iReg][iLoc][iChamber][iPlane][ibitxy]).fY ;
		    fHit[iChamber].fZ = (fMapData->fLut[iReg][iLoc][iChamber][iPlane][ibitxy]).fZ ;
		    
// 		    cout<<"2: X : "<<fMapData[iReg][iLoc][iChamber][iPlane][ibitxy].fX<<endl;
// 		    cout<<"2: Y : "<<fMapData[iReg][iLoc][iChamber][iPlane][ibitxy].fY<<endl;
// 		    cout<<"2: Z : "<<fMapData[iReg][iLoc][iChamber][iPlane][ibitxy].fZ<<endl;
		    
		    Yset[iChamber] = true ;
		  }
		}

	      }// loop of ibitxy
	    }// iplane
	    
// 	    cout<<endl;
	      
	    }// ichamber
	  
	  for (int iChamber = 0; iChamber < 4 ; iChamber++){ //4 chambers per DDL 
	    if(Yset[iChamber] and Xset[iChamber]){
// 	      cout<<"chamber : "<<iChamber
// 		  <<", (X, Y, Z) : ("<<fHit[iChamber].fX
// 		  <<", "<<fHit[iChamber].fY
// 		  <<", "<<fHit[iChamber].fZ
// 		  <<") "<<endl;
	      fData[fDataCount].fDetElemId = fHit[iChamber].fDetElemId ;
	      fData[fDataCount].fX = fHit[iChamber].fX ;
	      fData[fDataCount].fY = fHit[iChamber].fY ;
	      fData[fDataCount].fZ = fHit[iChamber].fZ ;
	      fDataCount++;

	      if(fDataCount>(27000-1)){
		cerr<<"Running Out of memory "<<endl;
		return false;
	      }// if data overflow quit
	      
	    }// of both chamber
	  }// chamber loop
	      
	  nofTrigRec++;
	      
	      
	}// if any non zero pattern found
	
	index += 1 ; // the last word, important one
      }// Dec Condn
      
      if (not reg_phys_trig_occur)
	index += 45;
      
      index += 1; // end of local Data 0xCAFEFADE
      
      index = iLocIndex + 6; //important to reset the index counter for fake locids like 235 
    }// iLoc loop
  }// iReg Loop

  return true;
}
