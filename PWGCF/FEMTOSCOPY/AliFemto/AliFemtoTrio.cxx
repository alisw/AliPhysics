///
/// \file  AliFemtoTrio.h
/// \author Jeremi Niedziela

#include <TMath.h>
#include "AliFemtoTrio.h"

AliFemtoTrio::AliFemtoTrio():
fTrack1(nullptr),
fTrack2(nullptr),
fTrack3(nullptr)
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
  double m_inv = abs(fTrack1->FourMomentum() + fTrack2->FourMomentum() + fTrack3->FourMomentum());
  return m_inv;
}
