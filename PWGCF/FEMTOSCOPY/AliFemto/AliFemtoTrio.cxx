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

AliFemtoTrio::AliFemtoTrio(const AliFemtoTrio& trio):
fTrack1(trio.fTrack1),
fTrack2(trio.fTrack2),
fTrack3(trio.fTrack3),
fTrack1type(trio.fTrack1type),
fTrack2type(trio.fTrack2type),
fTrack3type(trio.fTrack3type)
{
  
}


AliFemtoTrio::~AliFemtoTrio()
{

}

AliFemtoTrio& AliFemtoTrio::operator=(const AliFemtoTrio& trio)
{
  fTrack1 = trio.fTrack1;
  fTrack2 = trio.fTrack2;
  fTrack3 = trio.fTrack3;
  fTrack1type = trio.fTrack1type;
  fTrack2type = trio.fTrack2type;
  fTrack3type = trio.fTrack3type;

  return *this;
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

double AliFemtoTrio::MInv12()
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
  
  double E  = p1.e()  + p2.e();
  double px = p1.px() + p2.px();
  double py = p1.py() + p2.py();
  double pz = p1.pz() + p2.pz();
  return sqrt(E*E-(px*px+py*py+pz*pz));
}

double AliFemtoTrio::MInv23()
{
  if(!fTrack1 || !fTrack2 || !fTrack3){
    cout<<"W - AliFemtoTrio::MInv - track missing in a trio"<<endl;
    cout<<"track1:"<<fTrack1<<endl;
    cout<<"track2:"<<fTrack2<<endl;
    cout<<"track3:"<<fTrack3<<endl;
    return -1.0;
  }
  AliFemtoLorentzVector p2 = fTrack2->FourMomentum();
  AliFemtoLorentzVector p3 = fTrack3->FourMomentum();
  
  double E  = p2.e()  + p3.e();
  double px = p2.px() + p3.px();
  double py = p2.py() + p3.py();
  double pz = p2.pz() + p3.pz();
  return sqrt(E*E-(px*px+py*py+pz*pz));
}

double AliFemtoTrio::MInv31()
{
  if(!fTrack1 || !fTrack2 || !fTrack3){
    cout<<"W - AliFemtoTrio::MInv - track missing in a trio"<<endl;
    cout<<"track1:"<<fTrack1<<endl;
    cout<<"track2:"<<fTrack2<<endl;
    cout<<"track3:"<<fTrack3<<endl;
    return -1.0;
  }
  AliFemtoLorentzVector p3 = fTrack3->FourMomentum();
  AliFemtoLorentzVector p1 = fTrack1->FourMomentum();
  
  double E  = p3.e()  + p1.e();
  double px = p3.px() + p1.px();
  double py = p3.py() + p1.py();
  double pz = p3.pz() + p1.pz();
  return sqrt(E*E-(px*px+py*py+pz*pz));
}

double AliFemtoTrio::GetTheta12()
{
  AliFemtoLorentzVector a = fTrack1->FourMomentum();
  AliFemtoLorentzVector b = fTrack2->FourMomentum();
  AliFemtoLorentzVector c = fTrack3->FourMomentum();
  
  AliFemtoLorentzVector r = a+b;
  AliFemtoLorentzVector je = r+c;
  
  je.boost(r);
  return r.Theta()-je.Theta();
}

double AliFemtoTrio::GetTheta23()
{
  AliFemtoLorentzVector a = fTrack2->FourMomentum();
  AliFemtoLorentzVector b = fTrack3->FourMomentum();
  AliFemtoLorentzVector c = fTrack1->FourMomentum();
  
  AliFemtoLorentzVector r = a+b;
  AliFemtoLorentzVector je = r+c;
  
  je.boost(r);
  return r.Theta()-je.Theta();
}

double AliFemtoTrio::GetTheta31()
{
  AliFemtoLorentzVector a = fTrack3->FourMomentum();
  AliFemtoLorentzVector b = fTrack1->FourMomentum();
  AliFemtoLorentzVector c = fTrack2->FourMomentum();
  
  AliFemtoLorentzVector r = a+b;
  AliFemtoLorentzVector je = r+c;
  
  je.boost(r);
  return r.Theta()-je.Theta();
}

double AliFemtoTrio::GetTheta1()
{
  AliFemtoLorentzVector a = fTrack2->FourMomentum();
  AliFemtoLorentzVector b = fTrack3->FourMomentum();
  AliFemtoLorentzVector c = fTrack1->FourMomentum();
  
  AliFemtoLorentzVector r = a+b;
  AliFemtoLorentzVector je = r+c;
  
  je.boost(c);
  return c.Theta()-je.Theta();
}

double AliFemtoTrio::GetTheta2()
{
  AliFemtoLorentzVector a = fTrack3->FourMomentum();
  AliFemtoLorentzVector b = fTrack1->FourMomentum();
  AliFemtoLorentzVector c = fTrack2->FourMomentum();
  
  AliFemtoLorentzVector r = a+b;
  AliFemtoLorentzVector je = r+c;
  
  je.boost(c);
  return c.Theta()-je.Theta();
}

double AliFemtoTrio::GetTheta3()
{
  AliFemtoLorentzVector a = fTrack1->FourMomentum();
  AliFemtoLorentzVector b = fTrack2->FourMomentum();
  AliFemtoLorentzVector c = fTrack3->FourMomentum();
  
  AliFemtoLorentzVector r = a+b;
  AliFemtoLorentzVector je = r+c;
  
  je.boost(c);
  return c.Theta()-je.Theta();
}
