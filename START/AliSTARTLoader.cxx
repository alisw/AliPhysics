//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// Loader for START digits and RecPoints inherit from AliBaseDataLoader     //
// STARTDigits is TObject  consists time flight signal from each PMTtube,   //
// mean time right and left array (start signal) and time differnce         //
// (vertex z position)                                                      // 
// START RecPoints is TObject with mean time (start signal)                 // 
// and evrtex position  (cm)                                                //
/////////////////////////////////////////////////////////////////////////////
#include "AliSTARTLoader.h"
#include <AliRunLoader.h>
#include "AliDataLoader.h"

#include <TTree.h>
#include <TFile.h>

#include <stdlib.h>
#include <Riostream.h>
#include <Riostream.h>


ClassImp(AliSTARTLoader)

/*****************************************************************************/ 
AliSTARTLoader::AliSTARTLoader(const Char_t *name,const Char_t *topfoldername):
  AliLoader(name,topfoldername)

{
//ctor   
  cout<<" AliSTARTLoader "<<endl;


   
}
/*****************************************************************************/ 

AliSTARTLoader::AliSTARTLoader(const Char_t *name,TFolder *topfolder):
 AliLoader(name,topfolder)
{
//ctor   

  cout<<"  My AliSTARTLoader!!!!! "<<endl;

   
}
/*****************************************************************************/ 
AliSTARTLoader::~AliSTARTLoader()
{
 //destructor
  UnloadDigits();
  fDataLoaders->Remove(&fDigitsDataLoader);

  UnloadRecPoints();
  fDataLoaders->Remove(&fVertexDataLoader);

 
}

