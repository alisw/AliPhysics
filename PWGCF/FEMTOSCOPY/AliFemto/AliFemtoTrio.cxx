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
  AliFemtoLorentzVector JeFrame = r+c;
  
  r.boost(JeFrame);
  c.boost(JeFrame);
  
//  cout<<"r\tpx:"<<r.px()<<"\t"<<r.py()<<"\t"<<r.pz()<<"\te:"<<r.e()<<endl;
//  cout<<"c\tpx:"<<c.px()<<"\t"<<c.py()<<"\t"<<c.pz()<<"\te:"<<c.e()<<endl;
  
  double modR = sqrt((r.px()*r.px()+r.py()*r.py()+r.pz()*r.pz()));
  double modC = sqrt((c.px()*c.px()+c.py()*c.py()+c.pz()*c.pz()));
  
  if(fabs(modR) < 0.0000001 || fabs(modC) < 0.0000001) return 0.0;
  
  double cosTheta = (r.px()*c.px()+r.py()*c.py()+r.pz()*c.pz());
//  cout<<"num:"<<cosTheta<<"\tden:"<<(modR*modC)<<endl;
  cosTheta /= (modR*modC);
  return acos(cosTheta);
}

double AliFemtoTrio::GetTheta23()
{
  AliFemtoLorentzVector a = fTrack2->FourMomentum();
  AliFemtoLorentzVector b = fTrack3->FourMomentum();
  AliFemtoLorentzVector c = fTrack1->FourMomentum();
  
  AliFemtoLorentzVector r = a+b;
  AliFemtoLorentzVector JeFrame = r+c;
  
  r.boost(JeFrame);
  c.boost(JeFrame);
  
  double modR = sqrt((r.px()*r.px()+r.py()*r.py()+r.pz()*r.pz()));
  double modC = sqrt((c.px()*c.px()+c.py()*c.py()+c.pz()*c.pz()));
  
  if(fabs(modR) < 0.0000001 || fabs(modC) < 0.0000001) return 0.0;
  
  double cosTheta = (r.px()*c.px()+r.py()*c.py()+r.pz()*c.pz());
  cosTheta /= (modR*modC);
  return acos(cosTheta);
}

double AliFemtoTrio::GetTheta31()
{
  AliFemtoLorentzVector a = fTrack3->FourMomentum();
  AliFemtoLorentzVector b = fTrack1->FourMomentum();
  AliFemtoLorentzVector c = fTrack2->FourMomentum();
  
  AliFemtoLorentzVector r = a+b;
  AliFemtoLorentzVector JeFrame = r+c;
  
  r.boost(JeFrame);
  c.boost(JeFrame);
  
  double modR = sqrt((r.px()*r.px()+r.py()*r.py()+r.pz()*r.pz()));
  double modC = sqrt((c.px()*c.px()+c.py()*c.py()+c.pz()*c.pz()));
  
  if(fabs(modR) < 0.0000001 || fabs(modC) < 0.0000001) return 0.0;
  
  double cosTheta = (r.px()*c.px()+r.py()*c.py()+r.pz()*c.pz());
  cosTheta /= (modR*modC);
  return acos(cosTheta);
}




