#include "AliFemtoXiTrackCut.h"
#include "AliESDtrack.h"
#include <cstdio>

#ifdef __ROOT__ 
ClassImp(AliFemtoXiTrackCut)
#endif



AliFemtoXiTrackCut::AliFemtoXiTrackCut() : AliFemtoV0TrackCut()
{
  // Default constructor
 }
 //------------------------------
AliFemtoXiTrackCut::~AliFemtoXiTrackCut(){
  /* noop */
}
//------------------------------
bool AliFemtoXiTrackCut::Pass(const AliFemtoXi* aXi)
{
  // test the particle and return 
  // true if it meets all the criteria
  // false if it doesn't meet at least one of the criteria
  
  
  if(!AliFemtoV0TrackCut::Pass(aXi))
    return false;
  
  
  return true;
    
    
}
//------------------------------
AliFemtoString AliFemtoXiTrackCut::Report()
{
  // Prepare report from the execution
  string tStemp;
  char tCtemp[100];


  AliFemtoString returnThis = tStemp;
  return returnThis;
}
TList *AliFemtoXiTrackCut::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
 
  return tListSetttings;
}
