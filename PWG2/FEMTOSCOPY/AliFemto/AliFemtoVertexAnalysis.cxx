////////////////////////////////////////////////////////////////////////////
//                                                                        //
// AliFemtoVertexAnalysis - Femtoscopic analysis which mixes events       //
// with respect to the z position of the primary vertex                   //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id$
 *
 * Author: Frank Laue, Ohio State, laue@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *      This is the Class for Analysis objects.  Each of the simultaneous
 *      Analyses running should have one of these instantiated.  They link
 *      into the Manager in an Analysis Collection.
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.2.2.2  2007/10/12 14:28:37  akisiel
 * New wave of cleanup and rule conformance
 *
 * Revision 1.2.2.1  2007/10/05 09:38:17  akisiel
 * Fix stray colons
 *
 * Revision 1.2  2007/07/09 16:17:11  mlisa
 * several files changed to change name of AliFemtoAnalysis to AliFemtoSimpleAnalysis and AliFemtoBaseAnalysis to AliFemtoAnalysis.  Also removed some hard-coded cuts of Marek
 *
 * Revision 1.1  2007/05/16 10:22:12  akisiel
 * Making the directory structure of AliFemto flat. All files go into one common directory
 *
 * Revision 1.2  2007/05/03 09:39:37  akisiel
 * Fixing Effective C++ warnings
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.5  2001/05/25 23:24:00  lisa
 * Added in AliFemtoKink stuff
 *
 * Revision 1.4  2000/08/31 22:31:32  laue
 * AliFemtoSimpleAnalysis: output changed (a little bit less)
 * AliFemtoEvent: new version, members for reference mult added
 * AliFemtoIOBinary: new IO for new AliFemtoEvent version
 * AliFemtoTypes: TTree typedef to AliFemtoTTree added
 * AliFemtoVertexAnalysis: overflow and underflow added
 *
 * Revision 1.1  2000/07/16 21:44:11  laue
 * Collection and analysis for vertex dependent event mixing
 *
 *
 **************************************************************************/

#include "AliFemtoVertexAnalysis.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoKinkCut.h"
#include "AliFemtoPicoEventCollectionVector.h"
#include "AliFemtoPicoEventCollectionVectorHideAway.h"


#ifdef __ROOT__ 
ClassImp(AliFemtoVertexAnalysis)
#endif

extern void FillHbtParticleCollection(AliFemtoParticleCut*         partCut,
				     AliFemtoEvent*               hbtEvent,
				     AliFemtoParticleCollection*  partCollection);


//____________________________
AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(unsigned int bins, double min, double max):
  fVertexBins(0),
  fOverFlow(0),  
  fUnderFlow(0)
{
  //  mControlSwitch     = 0;
  fEventCut          = 0;
  fFirstParticleCut  = 0;
  fSecondParticleCut = 0;
  fPairCut           = 0;
  fCorrFctnCollection= 0;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fVertexBins = bins;
  fVertexZ[0] = min;
  fVertexZ[1] = max;
  fUnderFlow = 0; 
  fOverFlow = 0; 
  if (fMixingBuffer) delete fMixingBuffer;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexBins,fVertexZ[0],fVertexZ[1]);
}
//____________________________

AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a) : 
  AliFemtoSimpleAnalysis(),
  fVertexBins(0),
  fOverFlow(0),  
  fUnderFlow(0)
{
  //AliFemtoVertexAnalysis();
  fEventCut          = 0;
  fFirstParticleCut  = 0;
  fSecondParticleCut = 0;
  fPairCut           = 0;
  fCorrFctnCollection= 0;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fVertexBins = a.fVertexBins; 
  fVertexZ[0] = a.fVertexZ[0]; 
  fVertexZ[1] = a.fVertexZ[1];
  fUnderFlow = 0; 
  fOverFlow = 0; 
  if (fMixingBuffer) delete fMixingBuffer;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexBins,fVertexZ[0],fVertexZ[1]);

  // find the right event cut
  fEventCut = a.fEventCut->Clone();
  // find the right first particle cut
  fFirstParticleCut = a.fFirstParticleCut->Clone();
  // find the right second particle cut
  if (a.fFirstParticleCut==a.fSecondParticleCut) 
    SetSecondParticleCut(fFirstParticleCut); // identical particle hbt
  else
  fSecondParticleCut = a.fSecondParticleCut->Clone();

  fPairCut = a.fPairCut->Clone();
  
  if ( fEventCut ) {
      SetEventCut(fEventCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a) - event cut set " << endl;
  }
  if ( fFirstParticleCut ) {
      SetFirstParticleCut(fFirstParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a) - first particle cut set " << endl;
  }
  if ( fSecondParticleCut ) {
      SetSecondParticleCut(fSecondParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a) - second particle cut set " << endl;
  }  if ( fPairCut ) {
      SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a) - pair cut set " << endl;
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=a.fCorrFctnCollection->begin(); iter!=a.fCorrFctnCollection->end();iter++){
    cout << " AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a) - looking for correlation functions " << endl;
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
    else cout << " AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a) - correlation function not found " << endl;
  }

  fNumEventsToMix = a.fNumEventsToMix;

  cout << " AliFemtoVertexAnalysis::AliFemtoVertexAnalysis(const AliFemtoVertexAnalysis& a) - analysis copied " << endl;

}
//____________________________
AliFemtoVertexAnalysis::~AliFemtoVertexAnalysis(){
  // now delete every PicoEvent in the EventMixingBuffer and then the Buffer itself
  delete fPicoEventCollectionVectorHideAway;
}

//____________________________
AliFemtoString AliFemtoVertexAnalysis::Report()
{
  // prepare report fromt the execution
  cout << "AliFemtoVertexAnalysis - constructing Report..."<<endl;
  char ctemp[200];
  AliFemtoString temp = "-----------\nHbt AliFemtoVertexAnalysis Report:\n";
  sprintf(ctemp,"Events are mixed in %d bins in the range %E cm to %E cm.\n",fVertexBins,fVertexZ[0],fVertexZ[1]);
  temp += ctemp;
  sprintf(ctemp,"Events underflowing: %d\n",fUnderFlow);
  temp += ctemp;
  sprintf(ctemp,"Events overflowing: %d\n",fOverFlow);
  temp += ctemp;
  sprintf(ctemp,"Now adding AliFemtoSimpleAnalysis(base) Report\n");
  temp += ctemp;
  temp += AliFemtoSimpleAnalysis::Report();
  AliFemtoString returnThis=temp;
  return returnThis;
}
//_________________________
void AliFemtoVertexAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent) {
  // Put the event though the analysis
  cout << " AliFemtoVertexAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent) " << endl;
  // get right mixing buffer
  double vertexZ = hbtEvent->PrimVertPos().z();
  fMixingBuffer = fPicoEventCollectionVectorHideAway->PicoEventCollection(vertexZ); 
  if (!fMixingBuffer) {
    if ( vertexZ < fVertexZ[0] ) fUnderFlow++;
    if ( vertexZ > fVertexZ[1] ) fOverFlow++;
    return;
  }
  // call ProcessEvent() from AliFemtoSimpleAnalysis-base
  AliFemtoSimpleAnalysis::ProcessEvent(hbtEvent);
}
