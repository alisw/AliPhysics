

//Root includes

//AliRoot includes
#include "AliMUONLoader.h"

ClassImp(AliMUONLoader)
//___________________________________________________________________

 
/*****************************************************************************/ 
AliMUONLoader::AliMUONLoader():AliLoader()
{
//default constructor

}
/******************************************************************/

AliMUONLoader::AliMUONLoader(const Char_t* detname,const Char_t* eventfoldername):AliLoader(detname,eventfoldername)
{
    
}
/*****************************************************************************/ 

AliMUONLoader::AliMUONLoader(const Char_t * detname,TFolder* eventfolder):AliLoader(detname,eventfolder)
{
//constructor
  
}
/*****************************************************************************/ 

AliMUONLoader::~AliMUONLoader()
{
//detructor
 
}
/*****************************************************************************/ 

