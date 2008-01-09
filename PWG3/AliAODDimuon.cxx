/* AliAODDimuon: a class for AODs for the MUON Arm of the ALICE Experiment
 * Author: P. Cortese, Universita' del Piemonte Orientale in Alessandria and
 * INFN of Torino - Italy
 */

#include "AliAODDimuon.h"
#define AliAODDimuon_CXX

ClassImp(AliAODDimuon)

AliAODDimuon::AliAODDimuon():ei(0),p(0),MProton(0.93827231)
{
  // default constructor
  //mu[0]=0;
  //mu[1]=0;
}

AliAODDimuon::AliAODDimuon(const AliAODDimuon& dimu):p(0),MProton(0.93827231)
{
  // default constructor
  mu[0]=dimu.mu[0];
  mu[1]=dimu.mu[1];
  ei=dimu.ei;
}

AliAODDimuon::AliAODDimuon(TObject *mu0, TObject *mu1, TObject *eipoint):p(0),MProton(0.93827231)
{
  ///printf("Creating dimuon from %p %p\n",mu0,mu1);
  mu[0]=mu0;
  mu[1]=mu1;
  ei=eipoint;
}

//______________________________________________________________________________
AliAODDimuon::~AliAODDimuon()
{
  // destructor
  if(p)delete p;
  p=0;
}

void AliAODDimuon::BookP(){
  static UInt_t UnID[2]={0,0};
  if(!p){
    p=new TLorentzVector(Px(),Py(),Pz(),E());
    UnID[0]=mu[0].GetUniqueID();
    UnID[1]=mu[1].GetUniqueID();
  }
  // For efficiency reasons
  if((UnID[0]!=mu[0].GetUniqueID())||(UnID[1]!=mu[1].GetUniqueID())){
    p->SetPxPyPzE(Px(),Py(),Pz(),E());
    UnID[0]=mu[0].GetUniqueID();
    UnID[1]=mu[1].GetUniqueID();
  }
}

//______________________________________________________________________________
Double_t AliAODDimuon::Px() const {
  if(this->CheckPointers())return -999999999;
  return ((AliAODTrack*)mu[0].GetObject())->Px()+((AliAODTrack*)mu[1].GetObject())->Px();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Py() const {
  if(this->CheckPointers())return -999999999;
  return ((AliAODTrack*)mu[0].GetObject())->Py()+((AliAODTrack*)mu[1].GetObject())->Py();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Pz() const {
  if(this->CheckPointers())return -999999999;
  return ((AliAODTrack*)mu[0].GetObject())->Pz()+((AliAODTrack*)mu[1].GetObject())->Pz();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Pt() const {
  if(this->CheckPointers())return -999999999;
  Double_t px=Px();
  Double_t py=Py();
  return TMath::Sqrt(px*px+py*py);
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::E() const {
  if(this->CheckPointers())return -999999999;
  return ((AliAODTrack*)mu[0].GetObject())->E()+((AliAODTrack*)mu[1].GetObject())->E();
}

//______________________________________________________________________________
Double_t AliAODDimuon::P() const {
  printf("You should never call: Double_t AliAODDimuon::P() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::P() {
  if(this->CheckPointers())return -999999999;
  BookP();
  return p->P();
}

//______________________________________________________________________________
Double_t AliAODDimuon::M() const {
  printf("You should never call: Double_t AliAODDimuon::M() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::M() {
  if(this->CheckPointers())return -999999999;
  BookP();
  return p->M();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Eta() const {
  printf("You should never call: Double_t AliAODDimuon::Eta() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::Eta() {
  if(this->CheckPointers())return -999999999;
  BookP();
  return p->Eta();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Phi() const {
  printf("You should never call: Double_t AliAODDimuon::Phi() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::Phi() {
  if(this->CheckPointers())return -999999999;
  BookP();
  return p->Phi();
}
//______________________________________________________________________________
Double_t AliAODDimuon::Theta() const {
  printf("You should never call: Double_t AliAODDimuon::Theta() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::Theta() {
  if(this->CheckPointers())return -999999999;
  BookP();
  return p->Theta();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Y() const {
  printf("You should never call: Double_t AliAODDimuon::Y() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::Y() {
  if(this->CheckPointers())return -999999999;
  BookP();
  return p->Rapidity();
}

//______________________________________________________________________________
Short_t AliAODDimuon::Charge() const
{
  if(this->CheckPointers())return -999;
  return ((AliAODTrack*)mu[0].GetObject())->Charge()+((AliAODTrack*)mu[1].GetObject())->Charge();
}

//______________________________________________________________________________
Int_t AliAODDimuon::CheckPointers() const{
  if(mu[0]==0||mu[1]==0){
    printf("Dimuon not initialized\n");
    return -999;
  }
  if((mu[0].GetObject())==0||(mu[1].GetObject())==0){
    printf("Can not get objects. Got: %p %p\n",mu[0].GetObject(),mu[1].GetObject());
    return -999;
  }
  return 0;
}

//______________________________________________________________________________
Double_t AliAODDimuon::xf() {
  Double_t EBeam=((AliAODEventInfo*)ei.GetObject())->EBeam();
  if(EBeam<=0){
    printf("AliAODDimuon::xf: can not compute xf with EBeam=%f\n",EBeam);
    return -999999999;
  }
  if(this->CheckPointers())return -999999999;
  BookP();
  Double_t MDimu=M();
  Double_t PMax=TMath::Sqrt(EBeam*EBeam-MDimu*MDimu);
  return Pz()/PMax;
}

//______________________________________________________________________________
// Calculation the Collins-Soper angle (adapted from code by R. Arnaldi)
Double_t AliAODDimuon::CostCS(){
  if(CheckPointers())return -999999999;
  if(ei==0){
    printf("Pointer to MuonHeader not initialized\n");
    return -999999999;
  }
  if(ei.GetObject()==0){
    printf("Can not get MuonHeader object\n");
    return -999999999;
  }
  Double_t EBeam=((AliAODEventInfo*)ei.GetObject())->EBeam();
  if(EBeam<=0){
    printf("Can not compute costCS with EBeam=%f\n",EBeam);
    return -999999999;
  }
  Double_t mp=MProton;
  Double_t PBeam=TMath::Sqrt(EBeam*EBeam-mp*mp);
  Double_t pla10=((AliAODTrack*)mu[0].GetObject())->Px();
  Double_t pla11=((AliAODTrack*)mu[0].GetObject())->Py();
  Double_t pla12=((AliAODTrack*)mu[0].GetObject())->Pz();
  Double_t e1=((AliAODTrack*)mu[0].GetObject())->E();
  Double_t Mu1Charge=((AliAODTrack*)mu[0].GetObject())->Charge();
  Double_t pla20=((AliAODTrack*)mu[1].GetObject())->Px();
  Double_t pla21=((AliAODTrack*)mu[1].GetObject())->Py();
  Double_t pla22=((AliAODTrack*)mu[1].GetObject())->Pz();
  Double_t e2=((AliAODTrack*)mu[1].GetObject())->E();
  Double_t Mu2Charge=((AliAODTrack*)mu[1].GetObject())->Charge();

  // Fill the Lorentz vector for projectile and target
  // For the moment we consider no crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector PProjLab(0.,0.,-PBeam,EBeam); // projectile
  TLorentzVector PTargLab(0.,0., PBeam,EBeam); // target
  //
  // --- Get the muons parameters in the LAB frame
  //
  TLorentzVector PMu1Lab(pla10,pla11,pla12,e1);
  TLorentzVector PMu2Lab(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the LAB frame
  //
  TLorentzVector PDimuLab=PMu1Lab+PMu2Lab;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta=(-1./PDimuLab.E())*PDimuLab.Vect();
  TLorentzVector PMu1Dimu=PMu1Lab;
  TLorentzVector PMu2Dimu=PMu2Lab;
  TLorentzVector PProjDimu=PProjLab;
  TLorentzVector PTargDimu=PTargLab;
  PMu1Dimu.Boost(beta);
  PMu2Dimu.Boost(beta);
  PProjDimu.Boost(beta);
  PTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle 
  //
  TVector3 zaxisCS=(((PProjDimu.Vect()).Unit())-((PTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  Double_t cost;
  if(Mu1Charge > 0) {
    cost = zaxisCS.Dot((PMu1Dimu.Vect()).Unit());
    // Theta CS is not properly defined for Like-Sign muons
    if(Mu2Charge > 0 && cost<0) cost=-cost;
  } else { 
    // Theta CS is not properly defined for Like-Sign muons
    cost = zaxisCS.Dot((PMu2Dimu.Vect()).Unit());
    if(Mu2Charge < 0 && cost<0) cost=-cost;
  }
  return cost;
}

//______________________________________________________________________________
// Calculation the Helicity polarization angle (adapted from code by R. Arnaldi)
Double_t AliAODDimuon::CostHe(){
  if(CheckPointers())return -999999999;
  if(ei==0){
    printf("Pointer to MuonHeader not initialized\n");
    return -999999999;
  }
  if(ei.GetObject()==0){
    printf("Can not get MuonHeader object\n");
    return -999999999;
  }
  Double_t EBeam=((AliAODEventInfo*)ei.GetObject())->EBeam();
  if(EBeam<=0){
    printf("Can not compute costCS with EBeam=%f\n",EBeam);
    return -999999999;
  }
  Double_t PBeam=TMath::Sqrt(EBeam*EBeam-MProton*MProton);
  Double_t pla10=((AliAODTrack*)mu[0].GetObject())->Px();
  Double_t pla11=((AliAODTrack*)mu[0].GetObject())->Py();
  Double_t pla12=((AliAODTrack*)mu[0].GetObject())->Pz();
  Double_t e1=((AliAODTrack*)mu[0].GetObject())->E();
  Double_t Mu1Charge=((AliAODTrack*)mu[0].GetObject())->Charge();
  Double_t pla20=((AliAODTrack*)mu[1].GetObject())->Px();
  Double_t pla21=((AliAODTrack*)mu[1].GetObject())->Py();
  Double_t pla22=((AliAODTrack*)mu[1].GetObject())->Pz();
  Double_t e2=((AliAODTrack*)mu[1].GetObject())->E();
  Double_t Mu2Charge=((AliAODTrack*)mu[1].GetObject())->Charge();

  // Fill the Lorentz vector for projectile and target
  // For the moment we consider no crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector PProjLab(0.,0.,-PBeam,EBeam); // projectile
  TLorentzVector PTargLab(0.,0., PBeam,EBeam); // target
  //
  // --- Get the muons parameters in the LAB frame
  //
  TLorentzVector PMu1Lab(pla10,pla11,pla12,e1);
  TLorentzVector PMu2Lab(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the LAB frame
  //
  TLorentzVector PDimuLab=PMu1Lab+PMu2Lab;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta=(-1./PDimuLab.E())*PDimuLab.Vect();
  TLorentzVector PMu1Dimu=PMu1Lab;
  TLorentzVector PMu2Dimu=PMu2Lab;
  PMu1Dimu.Boost(beta);
  PMu2Dimu.Boost(beta);
  //
  // --- Translate the dimuon parameters in the CM frame
  //
  TLorentzVector PDimuCM; //CM frame
  TVector3 beta2;
  beta2=(-1./(MProton+PProjLab.E()))*PProjLab.Vect();
  PDimuCM=PDimuLab;
  PDimuCM.Boost(beta2);
  //
  // --- Determine the z axis for the calculation of the polarization angle (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis;
  zaxis=(PDimuCM.Vect()).Unit();
  //
  // --- Calculation of the polarization angle (Kharzeev) (angle between mu+ and the z axis defined above)
  //
  Double_t cost;
  if(Mu1Charge > 0) {
    cost = zaxis.Dot((PMu1Dimu.Vect()).Unit());
    // Theta Kharzeev is not properly defined for Like-Sign muons
    if(Mu2Charge > 0 && cost<0) cost=-cost;
  } else { 
    cost = zaxis.Dot((PMu2Dimu.Vect()).Unit());
    // Theta Kharzeev is not properly defined for Like-Sign muons
    if(Mu2Charge < 0 && cost<0) cost=-cost;
  }  
  return cost;
}

Int_t AliAODDimuon::AnyPt(){
  if(this->CheckPointers())return 0;
  return (((AliAODTrack*)mu[0].GetObject())->MatchTriggerAnyPt())&&
         (((AliAODTrack*)mu[0].GetObject())->MatchTriggerAnyPt());
}

Int_t AliAODDimuon::LowPt(){
  if(this->CheckPointers())return 0;
  return (((AliAODTrack*)mu[0].GetObject())->MatchTriggerLowPt())&&
         (((AliAODTrack*)mu[0].GetObject())->MatchTriggerLowPt());
}

Int_t AliAODDimuon::HighPt(){
  if(this->CheckPointers())return 0;
  return (((AliAODTrack*)mu[0].GetObject())->MatchTriggerHighPt())&&
         (((AliAODTrack*)mu[0].GetObject())->MatchTriggerHighPt());
}

Double_t AliAODDimuon::MaxChi2Match(){
  if(this->CheckPointers())return -999999999;
  return TMath::Max((((AliAODTrack*)mu[0].GetObject())->GetChi2MatchTrigger()),
                    (((AliAODTrack*)mu[0].GetObject())->GetChi2MatchTrigger()));
}
