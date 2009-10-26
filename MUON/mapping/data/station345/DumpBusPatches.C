#include "AliMpDataProcessor.h"
#include "Riostream.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpDDLStore.h"
#include "AliMpBusPatch.h"

//#define RESET		0
//#define BRIGHT 		1
//#define DIM		2
//#define UNDERLINE 	3
//#define BLINK		4
//#define REVERSE		7
//#define HIDDEN		8
//#define BLACK 		0
//#define RED		1
//#define GREEN		2
//#define YELLOW		3
//#define BLUE		4
//#define MAGENTA		5
//#define CYAN		6
//#define	WHITE		7
//void textcolor(int attr, int fg, int bg)
//{	
//  char command[13];
//  
//	/* Command is the control command to the terminal */
//	sprintf(command, "%c[%d;%d;%dm", 0x1B, attr, fg + 30, bg + 40);
//	printf("%s", command);
//}
//
void Green(ostream& out)
{
  out << "GREEN ";
//  out << Form("%c[0;32;47m",0x1B);
}

void Blue(ostream& out)
{
  out << "BLUE  ";
//  out << Form("%c[0;34;47m",0x1B);
}

void Black(ostream& /*out*/)
{
//  out << Form("%c[0;30;47m",0x1B);
}

void DumpBusPatches()
{
  /// Dump the list of manus for all the bus patches
  
  AliMpDataProcessor mp;
  AliMpDataMap* map = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(map);
  AliMpDDLStore::ReadData(dataStreams);

  TIter nextBP(AliMpDDLStore::Instance()->CreateBusPatchIterator());
  AliMpBusPatch* bp(0x0);

  while ( ( bp = static_cast<AliMpBusPatch*>(nextBP())) )
  {
    Int_t busPatchId = bp->GetId();

    Int_t firstManu = bp->GetManuId(0);
    
    if ( firstManu > 1024 ) Blue(cout);
    else Green(cout);
    
    cout << Form("BP %5d N PATCH MODULES %2d (",busPatchId,bp->GetNofPatchModules());
    
    for ( Int_t i = 0;  i < bp->GetNofPatchModules(); ++i ) 
    {
      cout << Form(" %2d",bp->GetNofManusPerModule(i));
    }
    
    for ( Int_t i = bp->GetNofPatchModules(); i < 4; ++i ) 
    {
      cout << "   ";
    }
    
    cout << " ) MANUS ";
    
    for ( Int_t imanu = 0; imanu < bp->GetNofManus(); ++imanu )
    {
      Int_t manuId = bp->GetManuId(imanu);
      cout << Form(" %5d",manuId);
    }
    cout << endl;
  }
  
  Black(cout);
}