

//Root includes

//AliRoot includes
#include "AliMUONLoader.h"
#include "AliMUONConstants.h"

ClassImp(AliMUONLoader)
//___________________________________________________________________
AliMUONLoader::AliMUONLoader():AliLoader()
{
//default constructor
}
//_______________________________________________________________________________
AliMUONLoader::AliMUONLoader(const Char_t* detname,const Char_t* eventfoldername):AliLoader(detname,eventfoldername)
{
}
//_______________________________________________________________________________
AliMUONLoader::AliMUONLoader(const Char_t * detname,TFolder* eventfolder):AliLoader(detname,eventfolder)
{
//constructor
}
//_______________________________________________________________________________
AliMUONLoader::~AliMUONLoader()
{
//detructor 
}
//_______________________________________________________________________________
