////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoEventReaderNanoAODChain - the reader class for the Alice AOD from Chain//
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoEventReaderNanoAODChain.h"

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


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoEventReaderNanoAODChain);
  /// \endcond
#endif

#if !(ST_NO_NAMESPACES)
  using namespace units;
#endif

using namespace std;
//____________________________
//constructor with 0 parameters , look at default settings
AliFemtoEventReaderNanoAODChain::AliFemtoEventReaderNanoAODChain():
  AliFemtoEventReaderNanoAOD()
{
  // default constructor
}

AliFemtoEventReaderNanoAODChain::AliFemtoEventReaderNanoAODChain(const AliFemtoEventReaderNanoAODChain &aReader) :
  AliFemtoEventReaderNanoAOD(aReader)
{
  // copy constructor
}
//__________________
//Destructor
AliFemtoEventReaderNanoAODChain::~AliFemtoEventReaderNanoAODChain()
{
  // destructor
}

//__________________
AliFemtoEventReaderNanoAODChain& AliFemtoEventReaderNanoAODChain::operator=(const AliFemtoEventReaderNanoAODChain& aReader)
{
  // assignment operator
  if (this != &aReader) {
    AliFemtoEventReaderNanoAOD::operator=(aReader);
  }

  return *this;
}
//__________________
AliFemtoString AliFemtoEventReaderNanoAODChain::Report()
{
  // create reader report
  AliFemtoString temp = "\n This is the AliFemtoEventReaderNanoAODChain\n";
  return temp;
}

AliFemtoEvent* AliFemtoEventReaderNanoAODChain::ReturnHbtEvent()
{
  // read in a next hbt event from the chain
  // convert it to AliFemtoEvent and return
  // for further analysis
  if (!fEvent) {
    return nullptr;
  }

  AliFemtoEvent *hbtEvent = 0;

  // Get the PWG2 specific information if it exists
//   fPWG2AODTracks = (TClonesArray *) fEvent->GetList()->FindObject("pwg2aodtracks");

//   if (fPWG2AODTracks) {
//     cout << "Found additional PWG2 specific information in the AOD!" << endl;
//     cout << "Reading only tracks with the additional information" << endl;
//   }

  // cout<<"starting to read event "<<fCurEvent<<endl;

  //hbtEvent = new AliFemtoEvent;
  hbtEvent = CopyAODtoFemtoEvent();

  fCurEvent++;
  return hbtEvent;
}

//___________________
void AliFemtoEventReaderNanoAODChain::SetAODSource(AliAODEvent *aAOD)
{
  // The chain loads the AOD for us
  // You must provide the address where it can be found
  fEvent = aAOD;
}
