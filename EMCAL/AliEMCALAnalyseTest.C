#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>

#include "EMCAL/AliEMCALGetter.h"
#endif


void Go(){

  cout << "AliEMCAL:> Single File default reconstruction analysing" << endl ;
  
  AliEMCALGetter* gime=AliEMCALGetter::GetInstance("galice.root") ;  
  gime->Event(0,"SDR") ;
  
  if((gime->TowerRecPoints()==0)||(gime->TowerRecPoints()->At(0)==0)){
    cout << "        No TowerRecPoint  !!!!! " << endl ;
  }    

  if((gime->PreShowerRecPoints()==0)||(gime->PreShowerRecPoints()->At(0)==0)){
    cout << "        No PreShowerRecPoint  !!!!! " << endl ;
  } 
   
//   if((gime->TrackSegments()==0)||(gime->TrackSegments()->At(0)==0)){
//     cout << "        No TrackSegments !!!! " << endl ;
//   }    

//   if((gime->RecParticles()==0)||(gime->RecParticles()->At(0)==0)){
//     cout << "        No RecParticles !!!!! " << endl ;
//   }    

  cout << "AliEMCAL:> Single File default analysing finished" << endl ;


  cout << "AliEMCAL:> Single File branch TEST analyzing started" << endl ;

  gime=AliEMCALGetter::GetInstance("galice.root","test") ;  
  gime->Event(0,"SDR") ;

  if((gime->TowerRecPoints()==0)||(gime->TowerRecPoints()->At(0)==0)){
    cout << "No TowerRecPoint " << endl ;
  }    
  
  if((gime->PreShowerRecPoints()==0)||(gime->PreShowerRecPoints()->At(0)==0)){
    cout << "No PreShowerRecPoint " << endl ;
  }   
  
//   if((gime->TrackSegments()==0)||(gime->TrackSegments()->At(0)==0)){
//     cout << "No TrackSegments " << endl ;
//   }    

//   if((gime->RecParticles()==0)||(gime->RecParticles()->At(0)==0)){
//     cout << "No RecParticles " << endl ;
//   }   
 
  cout << "AliEMCAL:> Single File branch TEST reconstruction ended" << endl ;

  
  cout << "AliEMCAL:> Split File default reconstruction started" << endl ;
 
  gime=AliEMCALGetter::GetInstance("galice.root","Default",kTRUE) ;  
  gime->Event(0,"SDR") ;
 
  if((gime->TowerRecPoints()==0)||(gime->TowerRecPoints()->At(0)==0)){
    cout << "No TowerRecPoint " << endl ;
  }    
  if((gime->PreShowerRecPoints()==0)||(gime->PreShowerRecPoints()->At(0)==0)){
    cout << "No PreShowerRecPoint " << endl ;
  }    
  
  //  if((gime->TrackSegments()==0)||(gime->TrackSegments()->At(0)==0)){
//     cout << "No TrackSegments " << endl ;
//   }   
 
//   if((gime->RecParticles()==0)||(gime->RecParticles()->At(0)==0)){
//     cout << "No RecParticles " << endl ;
//  }   
 
  cout << "AliEMCAL:> Split File default reconstruction ended" << endl ;
  
  cout << "--------AliEMCAL:> Reconstruction OK------------------"<< endl ;
  
}
