
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliFHistBrowser                                                      //
//                                                                      //
// helper class to browse AliFast Makers histograms.                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TBrowser.h>
#include "AliFast.h"
#include "AliFMaker.h"
#include "AliFHistBrowser.h"

ClassImp(AliFHistBrowser)



//_____________________________________________________________________________
AliFHistBrowser::AliFHistBrowser() 
                : TNamed("Histograms","ALIfast Histograms browser")
{

}

//_____________________________________________________________________________
void AliFHistBrowser::Browse(TBrowser *b)
{

  TIter next(gAliFast->Makers());
  AliFMaker *maker;
  while ((maker = (AliFMaker*)next())) {
     b->Add(maker->Histograms(),maker->GetName());
   }

}
