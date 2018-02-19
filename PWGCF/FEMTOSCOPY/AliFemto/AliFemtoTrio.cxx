///
/// \file  AliFemtoTrio.h
/// \author Jeremi Niedziela

#include <TMath.h>
#include "AliFemtoTrio.h"
#include "AliFemtoLorentzVector.h"

AliFemtoTrio::AliFemtoTrio():
fTrack1(nullptr),
fTrack2(nullptr),
fTrack3(nullptr),
fTrack1type(kUnknown),
fTrack2type(kUnknown),
fTrack3type(kUnknown)
{
  
}

AliFemtoTrio::~AliFemtoTrio()
{

}

double AliFemtoTrio::MInv()
{
  if(!fTrack1 || !fTrack2 || !fTrack3){
    cout<<"W - AliFemtoTrio::MInv - track missing in a trio"<<endl;
    cout<<"track1:"<<fTrack1<<endl;
    cout<<"track2:"<<fTrack2<<endl;
    cout<<"track3:"<<fTrack3<<endl;
    return -1.0;
  }
  AliFemtoLorentzVector p1 = fTrack1->FourMomentum();
  AliFemtoLorentzVector p2 = fTrack2->FourMomentum();
  AliFemtoLorentzVector p3 = fTrack3->FourMomentum();
  
  double E  = p1.e()  + p2.e()  + p3.e();
  double px = p1.px() + p2.px() + p3.px();
  double py = p1.py() + p2.py() + p3.py();
  double pz = p1.pz() + p2.pz() + p3.pz();
  
//  double m_inv = abs(fTrack1->FourMomentum() + fTrack2->FourMomentum() + fTrack3->FourMomentum());
  return sqrt(E*E-(px*px+py*py+pz*pz));
}
