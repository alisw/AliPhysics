#include <iostream>

#include "AliPHOSGetter.h"

void Go(){
  cout << "AliPHOS:> Single File default reconstruction analysing" << endl ;
  AliPHOSGetter* gime=AliPHOSGetter::GetInstance("galice.root") ;  
  gime->Event(0,"SDR") ;
  if((gime->EmcRecPoints()==0)||(gime->EmcRecPoints()->At(0)==0)){
    cout << "        No EmcRecPoint  !!!!! " << endl ;
  }    
  if((gime->CpvRecPoints()==0)||(gime->CpvRecPoints()->At(0)==0)){
    cout << "        No CpvRecPoint  !!!!! " << endl ;
  }    
  if((gime->TrackSegments()==0)||(gime->TrackSegments()->At(0)==0)){
    cout << "        No TrackSegments !!!! " << endl ;
  }    
  if((gime->RecParticles()==0)||(gime->RecParticles()->At(0)==0)){
    cout << "        No RecParticles !!!!! " << endl ;
  }    
  cout << "AliPHOS:> Single File default analysing finished" << endl ;
  
  
  cout << "AliPHOS:> Single File branch TEST analyzing started" << endl ;
  gime=AliPHOSGetter::GetInstance("galice.root","test") ;  
  gime->Event(0,"SDR") ;
  if((gime->EmcRecPoints()==0)||(gime->EmcRecPoints()->At(0)==0)){
    cout << "No EmcRecPoint " << endl ;
  }    
  if((gime->CpvRecPoints()==0)||(gime->CpvRecPoints()->At(0)==0)){
    cout << "No CpvRecPoint " << endl ;
  }    
  if((gime->TrackSegments()==0)||(gime->TrackSegments()->At(0)==0)){
    cout << "No TrackSegments " << endl ;
  }    
  if((gime->RecParticles()==0)||(gime->RecParticles()->At(0)==0)){
    cout << "No RecParticles " << endl ;
  }    
  cout << "AliPHOS:> Single File branch TEST reconstruction ended" << endl ;
  
  
  cout << "AliPHOS:> Split File default reconstruction started" << endl ;
  gime=AliPHOSGetter::GetInstance("galice.root","Default",kTRUE) ;  
  gime->Event(0,"SDR") ;
  if((gime->EmcRecPoints()==0)||(gime->EmcRecPoints()->At(0)==0)){
    cout << "No EmcRecPoint " << endl ;
  }    
  if((gime->CpvRecPoints()==0)||(gime->CpvRecPoints()->At(0)==0)){
    cout << "No CpvRecPoint " << endl ;
  }    
  if((gime->TrackSegments()==0)||(gime->TrackSegments()->At(0)==0)){
    cout << "No TrackSegments " << endl ;
  }    
  if((gime->RecParticles()==0)||(gime->RecParticles()->At(0)==0)){
    cout << "No RecParticles " << endl ;
  }    
  cout << "AliPHOS:> Split File default reconstruction ended" << endl ;
  cout << "--------AliPHOS:> Reconstruction OK------------------"<< endl ;
  
}
