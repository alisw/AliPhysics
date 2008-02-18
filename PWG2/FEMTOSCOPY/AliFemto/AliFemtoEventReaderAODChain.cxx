////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderAODChain - the reader class for the Alice AOD from Chain//
// Reads in AOD information and converts it into internal AliFemtoEvent       //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoEventReaderAODChain.h"

#include "TFile.h"
#include "TTree.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"

#include "AliFmPhysicalHelixD.h"
#include "AliFmThreeVectorF.h"

#include "SystemOfUnits.h"

#include "AliFemtoEvent.h"
#include "AliFemtoModelHiddenInfo.h"

ClassImp(AliFemtoEventReaderAODChain)

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
//constructor with 0 parameters , look at default settings 
AliFemtoEventReaderAODChain::AliFemtoEventReaderAODChain():
  AliFemtoEventReaderAOD()
{
  // default constructor
}

AliFemtoEventReaderAODChain::AliFemtoEventReaderAODChain(const AliFemtoEventReaderAODChain &aReader) :
  AliFemtoEventReaderAOD(aReader)
{
  // copy constructor
}
//__________________
//Destructor
AliFemtoEventReaderAODChain::~AliFemtoEventReaderAODChain()
{
  // destructor
}

//__________________
AliFemtoEventReaderAODChain& AliFemtoEventReaderAODChain::operator=(const AliFemtoEventReaderAODChain& aReader)
{
  // assignment operator
  if (this == &aReader)
    return *this;

  *this = aReader;

  return *this;
}
//__________________
AliFemtoString AliFemtoEventReaderAODChain::Report()
{
  // create reader report
  AliFemtoString temp = "\n This is the AliFemtoEventReaderAODChain\n";
  return temp;
}

AliFemtoEvent* AliFemtoEventReaderAODChain::ReturnHbtEvent()
{
  // read in a next hbt event from the chain
  // convert it to AliFemtoEvent and return
  // for further analysis
  if (!fEvent) return 0;

  AliFemtoEvent *hbtEvent = 0;

  // Get the PWG2 specific information if it exists
  fPWG2AODTracks = (TClonesArray *) fEvent->GetList()->FindObject("pwg2aodtracks");
  
  if (fPWG2AODTracks) {
    cout << "Found additional PWG2 specific information in the AOD!" << endl;
    cout << "Reading only tracks with the additional information" << endl;
  }

  cout<<"starting to read event "<<fCurEvent<<endl;
	
  hbtEvent = new AliFemtoEvent;

  CopyAODtoFemtoEvent(hbtEvent);

  fCurEvent++;	
  return hbtEvent; 
}

//___________________
void AliFemtoEventReaderAODChain::SetAODSource(AliAODEvent *aAOD)
{
  // The chain loads the AOD for us
  // You must provide the address where it can be found
  fEvent = aAOD;
}




