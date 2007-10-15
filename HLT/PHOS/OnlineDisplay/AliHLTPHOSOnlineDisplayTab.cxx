#include "AliHLTPHOSOnlineDisplayTab.h"
#include "HOMERReader.h"
#include <iostream>
#include "AliHLTDataTypes.h"


using namespace std;


AliHLTPHOSOnlineDisplayTab::AliHLTPHOSOnlineDisplayTab():AliHLTPHOSBase(),
							 fgSyncronize(kFALSE)
{

}

AliHLTPHOSOnlineDisplayTab::~AliHLTPHOSOnlineDisplayTab()
{

}


void 
AliHLTPHOSOnlineDisplayTab::PrintBlockInfo(HOMERReader *homeReaderPtr, int i)
{
  char tmp1[9], tmp2[5];
  memset( tmp1, 0, 9 );
  memset( tmp2, 0, 5);
  void *tmp11 = tmp1;
  ULong64_t* tmp12 = (ULong64_t*)tmp11;
  *tmp12 =homeReaderPtr->GetBlockDataType( i );
  void *tmp21 = tmp2;
  ULong_t* tmp22 = (ULong_t*)tmp21;
  *tmp22 = homeReaderPtr->GetBlockDataOrigin( i );
  cout << "Dataype for block:  "<< i<<"  is:  "<< tmp1<<tmp2 <<endl;
}


int 
AliHLTPHOSOnlineDisplayTab::GetEventInfo(HOMERReader *homeReaderPtr, int i)
{
  int ret = 0;
  ret =homeReaderPtr->ReadNextEvent();  
  if( ret ) 
    {
      int ndx = homeReaderPtr->GetErrorConnectionNdx();
      printf( "------------ TRY AGAIN --------------->Error reading event from source %d: %s (%d)\n", ndx, strerror(ret), ret );
      cout << "HOMER getconncetioNdx  status = " << ndx << endl;
      return -1; 
    }
  else
    {
      unsigned long blockCnt = homeReaderPtr->GetBlockCnt();
      cout << "AliHLTPHOSOnlineDisplayEventTab::GetNextEvent:  blockCnt ="  << blockCnt <<endl;
      return blockCnt;
    }
}


int 
AliHLTPHOSOnlineDisplayTab::DoGetNextEvent()
{
  HOMERReader* CurrentReaderPtr;
  int ret = 0;
  unsigned long ndx;
  const AliHLTComponentBlockData* iter = NULL;   
  Bool_t nextSwitch=kTRUE; 
  int nLoops=0;

  if(fgSyncronize == kTRUE)
    {
      nLoops = 1;
    }
  else
    {
      nLoops =  fgNHosts; 
    }

  for(int reader = 0; reader <  nLoops; reader ++)
    {
      if(fgSyncronize == kTRUE)
	{
	  CurrentReaderPtr =fgHomerReaderPtr;
	}
      else
	{
	  CurrentReaderPtr =fgHomerReadersPtr[reader];
	}

      ret = GetEventInfo(CurrentReaderPtr, reader);

      if( ret < 0) 
      	{
	  return ret; 
      	}
      else
	{
	  unsigned long blockCnt = ret;
	  for ( unsigned long i = 0; i < blockCnt; i++ ) 
	    {
	      PrintBlockInfo(CurrentReaderPtr, i);
	    }
	  ReadBlockData(CurrentReaderPtr);
	}
    }
}
