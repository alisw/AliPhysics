// $Id$

#include <Riostream.h>
#include <vector>

#include <TParticle.h>
#include <TMath.h>

#include "AliJFCluster.h"

ClassImp(AliJFCluster)

AliJFCluster::AliJFCluster(Int_t n) : fStatus(0),fNMerge(0),fPx(0),fPy(0),fPz(0),fE(0),
                                fY(0),fPhi(0),fPt2(0),fPt2dD(0),fList(n)
{
}

AliJFCluster::AliJFCluster(const AliJFCluster &copy)
{
  fStatus=copy.GetStatus();
  fNMerge=copy.GetNMerge();
  fPx=copy.GetPx();
  fPy=copy.GetPy();
  fPz=copy.GetPz();
  fE=copy.GetE();
  fY=copy.GetY();
  fPhi=copy.GetPhi();
  fPt2=copy.GetPt2();
  fPt2dD=copy.GetPt2D();
  fList=*copy.GetClusterList();
}

AliJFCluster::AliJFCluster(AliJFPreCluster &copy)
{
  fStatus=1; //valid
  fNMerge=0;
  fPx=copy.GetPx();
  fPy=copy.GetPy();
  fPz=copy.GetPz();
  fE=copy.GetE();

  SetValues();
  fList.push_back(&copy);
}

AliJFCluster::AliJFCluster(AliJFPreCluster *precluster)
{
  fStatus=1; //valid
  fNMerge=0;
  SetValues(precluster->GetPx(),fPy=precluster->GetPy(),fPz=precluster->GetPz(),fE=precluster->GetE());
  fList.push_back(precluster);
}

AliJFCluster::~AliJFCluster()
{
  fList.erase(fList.begin(),fList.end());
}

AliJFCluster& AliJFCluster::operator=(const AliJFCluster &rhs)
{
  fStatus=rhs.GetStatus();
  fNMerge=rhs.GetNMerge();
  fPx=rhs.GetPx();
  fPy=rhs.GetPy();
  fPz=rhs.GetPz();
  fE=rhs.GetE();
  fY=rhs.GetY();
  fPhi=rhs.GetPhi();
  fPt2=rhs.GetPt2();
  fPt2dD=rhs.GetPt2D();
  fList=*rhs.GetClusterList();
  return *this;
}

AliJFCluster& AliJFCluster::operator=(AliJFPreCluster &rhs)
{
  fStatus=1; //valid
  fNMerge=0;
  fPx=rhs.GetPx();
  fPy=rhs.GetPy();
  fPz=rhs.GetPz();
  fE=rhs.GetE();

  SetValues();
  fList.push_back(&rhs);
  return *this;
}

AliJFCluster& AliJFCluster::operator+=(AliJFCluster &rhs) //mark rhs as being merged!
{ 
  if(!rhs.IsValid()){
    cerr << "Cluster cannot be combined with invalid cluster:" << endl;
    cerr << *this << endl;
    cerr << rhs << endl;    
    return *this;
  }

  AddValues(rhs.GetPx(),fPy=rhs.GetPy(),fPz=rhs.GetPz(),fE=rhs.GetE());
  vector<AliJFPreCluster*> tmp=*rhs.GetClusterList();
  for(vector<AliJFPreCluster*>::iterator i=tmp.begin();i!=tmp.end();i++){
    fList.push_back(*i);
  }

  rhs.MarkIsMerged(); //change even rhs!

  return *this;
}

AliJFCluster& AliJFCluster::operator+=(AliJFPreCluster &rhs)
{
  AddValues(rhs.GetPx(),rhs.GetPy(),rhs.GetPz(),rhs.GetE());
  fList.push_back(&rhs);

  return *this;
}

ostream& operator<<(ostream& o, const AliJFCluster &c)
{
  o << c.GetStatus() << ": " << c.GetPx() << " " << c.GetPy() << " " << c.GetPz() << " " << c.GetE();

  return o;
}

void AliJFCluster::CombineCluster(AliJFCluster &rhs) //mark rhs as being merged!
{ 
  if(!rhs.IsValid()){
    cerr << "Error AliJFCluster: Cluster cannot be combined with invalid cluster:" << endl;
    cerr << *this << endl;
    cerr << rhs << endl;    
    return;
  }
  AddValues(rhs.GetPx(),fPy=rhs.GetPy(),rhs.GetPz(),rhs.GetE());
  vector<AliJFPreCluster*> tmp=*rhs.GetClusterList();
  for(vector<AliJFPreCluster*>::iterator i=tmp.begin();i!=tmp.end();i++){
    fList.push_back(*i);
  }

  rhs.MarkIsMerged(); //change even rhs!
}

void AliJFCluster::Print()
{
  cout << "Cluster " << " " << *this << endl;
  Int_t n=0;
  for(vector<AliJFPreCluster*>::iterator i=fList.begin();i!=fList.end();i++){
    n++;
    cout << "PreCluster " << n << ": " << *(*i) << endl;
  }
}

void AliJFCluster::SetValues(Float_t px, Float_t py, Float_t pz, Float_t E)
{
  fPx=px;
  fPy=py;
  fPz=pz;
  fE=E;

  SetValues();
}

void AliJFCluster::SetValues()
{
  if(IsValid()){
    fPt2=fPx*fPx+fPy*fPy;
    if(fE<0) fE=TMath::Sqrt(fPt2+fPz*fPz);
    fPt2dD=fPt2/D2;
    fPhi=TMath::Pi()+TMath::ATan2(-fPy,-fPx);
    fY=0.5*TMath::Log((fE+fPz)/(fE-fPz));
    if(fY>10) fY=10;
    else if(fY<-10) fY=-10;
  } else {
    fPt2=fPt2dD=fPhi=fY=0;
  }
}

void AliJFCluster::AddValues(Float_t px, Float_t py, Float_t pz, Float_t E)
{
  fNMerge++;
  fPx+=px;
  fPy+=py;
  fPz+=pz;
  fE+=E;

  SetValues();
}

Float_t AliJFCluster::D2=1.0;
