

//Root includes

//AliRoot includes
#include "AliMUONLoader.h"
#include "AliMUONConstants.h"

ClassImp(AliMUONLoader)
//___________________________________________________________________
AliMUONLoader::AliMUONLoader()
  : AliLoader(),
    fMUONData(0)
{
//default constructor
}
//_______________________________________________________________________________
AliMUONLoader::AliMUONLoader(const Char_t* detname,const Char_t* eventfoldername)
  : AliLoader(detname,eventfoldername),
    fMUONData(0)
{
}
//_______________________________________________________________________________
AliMUONLoader::AliMUONLoader(const Char_t * detname,TFolder* eventfolder)
  : AliLoader(detname,eventfolder),
  fMUONData(0)
{
//constructor
}
//_______________________________________________________________________________
AliMUONLoader::~AliMUONLoader()
{
//detructor 
}
//_______________________________________________________________________________


