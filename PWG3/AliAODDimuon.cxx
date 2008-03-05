// AliAODDimuon: a class for AODs for the MUON Arm of the ALICE Experiment
// Author: P. Cortese, Universita' del Piemonte Orientale in Alessandria and
// INFN of Torino - Italy
//
// The class defines a dimuon pair object from two AliAODTrack objects.
// AliAODDimuon objects are supposed to be added to the AliAODEvent structure
// during analysis. They would then allow to calculate the dimuon-related
// kinematic variables with a minimal disk occupancy.
// The payload of the class has been reduced to two pointers to the two
// tracks with the addition of a pointer to the AliAODEventInfo. An instance of
// this class has also to be added to the AliAODEvent structure to provide
// additional information that is specific to MUON and therefore has not been
// included into the AOD header.
// Two transient data members are not stored on file as they can be recomputed
// at runtime.
//

#include "AliAODDimuon.h"
#include "TLorentzVector.h"
#define AliAODDimuon_CXX

ClassImp(AliAODDimuon)

//______________________________________________________________________________
AliAODDimuon::AliAODDimuon():fEi(0),fP(0),fMProton(0.93827231)
{
  // default constructor
  fMu[0]=0;
  fMu[1]=0;
}

//______________________________________________________________________________
AliAODDimuon::AliAODDimuon(const AliAODDimuon& dimu):fP(0),fMProton(0.93827231)
{
  // copy constructor
  fMu[0]=dimu.Mu(0);
  fMu[1]=dimu.Mu(1);
  fEi=dimu.Ei();
}

//______________________________________________________________________________
AliAODDimuon &AliAODDimuon::operator=(const AliAODDimuon& dimu)
{
  // assignment operator
  fP=0;
  fMProton=0.93827231;
  if(&dimu != this){
    fMu[0]=dimu.Mu(0);
    fMu[1]=dimu.Mu(1);
    fEi=dimu.Ei();
  }
  return *this;
}

//______________________________________________________________________________
AliAODDimuon::AliAODDimuon(TObject *mu0, TObject *mu1, TObject *ei):
  fP(0),fMProton(0.93827231)
{
  // Creates a dimuon pair from two tracks and the EventInfo
  
  //printf("Creating dimuon from %p %p\n",mu0,mu1);
  fMu[0]=mu0;
  fMu[1]=mu1;
  fEi=ei;
}

//______________________________________________________________________________
AliAODDimuon::~AliAODDimuon()
{
  // destructor
  if(fP)delete fP;
  fP=0;
}

//______________________________________________________________________________
void AliAODDimuon::BookP(){
  // Fills the dimuon momentum if not filled yet
  static UInt_t unID[2]={0,0};
  if(!fP){
    fP=new TLorentzVector(Px(),Py(),Pz(),E());
    unID[0]=fMu[0].GetUniqueID();
    unID[1]=fMu[1].GetUniqueID();
  }
  // For efficiency reasons
  if((unID[0]!=fMu[0].GetUniqueID())||(unID[1]!=fMu[1].GetUniqueID())){
    fP->SetPxPyPzE(Px(),Py(),Pz(),E());
    unID[0]=fMu[0].GetUniqueID();
    unID[1]=fMu[1].GetUniqueID();
  }
}

//______________________________________________________________________________
Double_t AliAODDimuon::Px() const {
  // Px of the dimuon
  if(this->CheckPointers())return -999999999;
  return ((AliAODTrack*)fMu[0].GetObject())->Px()+
         ((AliAODTrack*)fMu[1].GetObject())->Px();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Py() const {
  // Py of the dimuon
  if(this->CheckPointers())return -999999999;
  return ((AliAODTrack*)fMu[0].GetObject())->Py()+
         ((AliAODTrack*)fMu[1].GetObject())->Py();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Pz() const {
  // Pz of the dimuon
  if(this->CheckPointers())return -999999999;
  return ((AliAODTrack*)fMu[0].GetObject())->Pz()+
         ((AliAODTrack*)fMu[1].GetObject())->Pz();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Pt() const {
  // Pt of the dimuon
  if(this->CheckPointers())return -999999999;
  Double_t px=Px();
  Double_t py=Py();
  return TMath::Sqrt(px*px+py*py);
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::E() const {
  // Dimuon energy
  if(this->CheckPointers())return -999999999;
  return ((AliAODTrack*)fMu[0].GetObject())->E()+
         ((AliAODTrack*)fMu[1].GetObject())->E();
}

//______________________________________________________________________________
Double_t AliAODDimuon::P() const {
  // This is just to override the virtual function
  printf("You should never call: Double_t AliAODDimuon::P() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::P() {
  // Dimuon momentum
  if(this->CheckPointers())return -999999999;
  BookP();
  return fP->P();
}

//______________________________________________________________________________
Double_t AliAODDimuon::M() const {
  // This is just to override the virtual function
  printf("You should never call: Double_t AliAODDimuon::M() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::M() {
  // Dimuon invariant mass
  if(this->CheckPointers())return -999999999;
  BookP();
  return fP->M();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Mass() {
  // Dimuon invariant mass
  if(this->CheckPointers())return -999999999;
  BookP();
  return fP->M();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Eta() const {
  // This is just to override the virtual function
  printf("You should never call: Double_t AliAODDimuon::Eta() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::Eta() {
  // Dimuon pseudorapidity
  if(this->CheckPointers())return -999999999;
  BookP();
  return fP->Eta();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Phi() const {
  // This is just to override the virtual function
  printf("You should never call: Double_t AliAODDimuon::Phi() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::Phi() {
  // Dimuon asimuthal angle
  if(this->CheckPointers())return -999999999;
  BookP();
  return fP->Phi();
}
//______________________________________________________________________________
Double_t AliAODDimuon::Theta() const {
  // This is just to override the virtual function
  printf("You should never call: Double_t AliAODDimuon::Theta() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::Theta() {
  // Dimuon polar angle
  if(this->CheckPointers())return -999999999;
  BookP();
  return fP->Theta();
}

//______________________________________________________________________________
Double_t AliAODDimuon::Y() const {
  // This is just to override the virtual function
  printf("You should never call: Double_t AliAODDimuon::Y() const\n");
  return -999999999;
}

//______________________________________________________________________________
Double_t AliAODDimuon::Y() {
  // Dimuon rapidity
  if(this->CheckPointers())return -999999999;
  BookP();
  return fP->Rapidity();
}

//______________________________________________________________________________
Short_t AliAODDimuon::Charge() const {
  // Dimuon charge
  if(this->CheckPointers())return -999;
  return ((AliAODTrack*)fMu[0].GetObject())->Charge()+
         ((AliAODTrack*)fMu[1].GetObject())->Charge();
}

//______________________________________________________________________________
Int_t AliAODDimuon::CheckPointers() const{
  // Checks if the track pointers have been initialized
  if(fMu[0]==0||fMu[1]==0){
    printf("Dimuon not initialized\n");
    return -999;
  }
  if((fMu[0].GetObject())==0||(fMu[1].GetObject())==0){
    printf("Can not get objects. Got: %p %p\n",fMu[0].GetObject(),fMu[1].GetObject());
    return -999;
  }
  return 0;
}

//______________________________________________________________________________
void AliAODDimuon::SetMu(Int_t imu, AliAODTrack *mu){
  // Assign a track pointer
  if (imu==0||imu==1){
    fMu[imu]=mu;
  }
}

//______________________________________________________________________________
void AliAODDimuon::SetMuons(AliAODTrack *mu0, AliAODTrack *mu1){
  // Assign the track pointers
  fMu[0]=mu0;
  fMu[1]=mu1;
}

//______________________________________________________________________________
Double_t AliAODDimuon::XF() {
  // Dimuon Feynman x
  Double_t ebeam=((AliAODEventInfo*)fEi.GetObject())->EBeam();
  if(ebeam<=0){
    printf("AliAODDimuon::xf: can not compute xf with EBeam=%f\n",ebeam);
    return -999999999;
  }
  if(this->CheckPointers())return -999999999;
  BookP();
  Double_t mDimu=M();
  Double_t pMax=TMath::Sqrt(ebeam*ebeam-mDimu*mDimu);
  return Pz()/pMax;
}

//______________________________________________________________________________
// Calculation the Collins-Soper angle (adapted from code by R. Arnaldi)
Double_t AliAODDimuon::CostCS(){
  // Cosinus of the Collins-Soper polar decay angle
  if(CheckPointers())return -999999999;
  if(fEi==0){
    printf("Pointer to MuonHeader not initialized\n");
    return -999999999;
  }
  if(fEi.GetObject()==0){
    printf("Can not get MuonHeader object\n");
    return -999999999;
  }
  Double_t ebeam=((AliAODEventInfo*)fEi.GetObject())->EBeam();
  if(ebeam<=0){
    printf("Can not compute costCS with EBeam=%f\n",ebeam);
    return -999999999;
  }
  Double_t mp=fMProton;
  Double_t pbeam=TMath::Sqrt(ebeam*ebeam-mp*mp);
  Double_t pla10=((AliAODTrack*)fMu[0].GetObject())->Px();
  Double_t pla11=((AliAODTrack*)fMu[0].GetObject())->Py();
  Double_t pla12=((AliAODTrack*)fMu[0].GetObject())->Pz();
  Double_t e1=((AliAODTrack*)fMu[0].GetObject())->E();
  Double_t mu1Charge=((AliAODTrack*)fMu[0].GetObject())->Charge();
  Double_t pla20=((AliAODTrack*)fMu[1].GetObject())->Px();
  Double_t pla21=((AliAODTrack*)fMu[1].GetObject())->Py();
  Double_t pla22=((AliAODTrack*)fMu[1].GetObject())->Pz();
  Double_t e2=((AliAODTrack*)fMu[1].GetObject())->E();
  Double_t mu2Charge=((AliAODTrack*)fMu[1].GetObject())->Charge();

  // Fill the Lorentz vector for projectile and target
  // For the moment we do not consider the crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjLab(0.,0.,-pbeam,ebeam); // projectile
  TLorentzVector pTargLab(0.,0., pbeam,ebeam); // target
  //
  // --- Get the muons parameters in the LAB frame
  //
  TLorentzVector pMu1Lab(pla10,pla11,pla12,e1);
  TLorentzVector pMu2Lab(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the LAB frame
  //
  TLorentzVector pDimuLab=pMu1Lab+pMu2Lab;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta=(-1./pDimuLab.E())*pDimuLab.Vect();
  TLorentzVector pMu1Dimu=pMu1Lab;
  TLorentzVector pMu2Dimu=pMu2Lab;
  TLorentzVector pProjDimu=pProjLab;
  TLorentzVector pTargDimu=pTargLab;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle 
  //
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  Double_t cost;
  if(mu1Charge > 0) {
    cost = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
    // Theta CS is not properly defined for Like-Sign muons
    if(mu2Charge > 0 && cost<0) cost=-cost;
  } else { 
    // Theta CS is not properly defined for Like-Sign muons
    cost = zaxisCS.Dot((pMu2Dimu.Vect()).Unit());
    if(mu2Charge < 0 && cost<0) cost=-cost;
  }
  return cost;
}

//______________________________________________________________________________
// Calculation the Helicity polarization angle (adapted from code by R. Arnaldi)
Double_t AliAODDimuon::CostHe(){
  // Cosinus of the polar decay angle in the Helicity reference frame
  if(CheckPointers())return -999999999;
  if(fEi==0){
    printf("Pointer to MuonHeader not initialized\n");
    return -999999999;
  }
  if(fEi.GetObject()==0){
    printf("Can not get MuonHeader object\n");
    return -999999999;
  }
  Double_t ebeam=((AliAODEventInfo*)fEi.GetObject())->EBeam();
  if(ebeam<=0){
    printf("Can not compute costCS with EBeam=%f\n",ebeam);
    return -999999999;
  }
  Double_t pbeam=TMath::Sqrt(ebeam*ebeam-fMProton*fMProton);
  Double_t pla10=((AliAODTrack*)fMu[0].GetObject())->Px();
  Double_t pla11=((AliAODTrack*)fMu[0].GetObject())->Py();
  Double_t pla12=((AliAODTrack*)fMu[0].GetObject())->Pz();
  Double_t e1=((AliAODTrack*)fMu[0].GetObject())->E();
  Double_t mu1Charge=((AliAODTrack*)fMu[0].GetObject())->Charge();
  Double_t pla20=((AliAODTrack*)fMu[1].GetObject())->Px();
  Double_t pla21=((AliAODTrack*)fMu[1].GetObject())->Py();
  Double_t pla22=((AliAODTrack*)fMu[1].GetObject())->Pz();
  Double_t e2=((AliAODTrack*)fMu[1].GetObject())->E();
  Double_t mu2Charge=((AliAODTrack*)fMu[1].GetObject())->Charge();

  // Fill the Lorentz vector for projectile and target
  // For the moment we consider no crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjLab(0.,0.,-pbeam,ebeam); // projectile
  TLorentzVector pTargLab(0.,0., pbeam,ebeam); // target
  //
  // --- Get the muons parameters in the LAB frame
  //
  TLorentzVector pMu1Lab(pla10,pla11,pla12,e1);
  TLorentzVector pMu2Lab(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the LAB frame
  //
  TLorentzVector pDimuLab=pMu1Lab+pMu2Lab;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta=(-1./pDimuLab.E())*pDimuLab.Vect();
  TLorentzVector pMu1Dimu=pMu1Lab;
  TLorentzVector pMu2Dimu=pMu2Lab;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  //
  // --- Translate the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM; //CM frame
  TVector3 beta2;
  beta2=(-1./(fMProton+pProjLab.E()))*pProjLab.Vect();
  pDimuCM=pDimuLab;
  pDimuCM.Boost(beta2);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis;
  zaxis=(pDimuCM.Vect()).Unit();
  //
  // --- Calculation of the polarization angle (Helicity)
  // (angle between mu+ and the z axis defined above)
  //
  Double_t cost;
  if(mu1Charge > 0) {
    cost = zaxis.Dot((pMu1Dimu.Vect()).Unit());
    // Theta Helicity is not properly defined for Like-Sign muons
    if(mu2Charge > 0 && cost<0) cost=-cost;
  } else { 
    cost = zaxis.Dot((pMu2Dimu.Vect()).Unit());
    // Theta Helicity is not properly defined for Like-Sign muons
    if(mu2Charge < 0 && cost<0) cost=-cost;
  }  
  return cost;
}

//______________________________________________________________________________
Int_t AliAODDimuon::AnyPt(){
  // Test if the two muons match two trigger tracks
  if(this->CheckPointers())return 0;
  return (((AliAODTrack*)fMu[0].GetObject())->MatchTriggerAnyPt())&&
         (((AliAODTrack*)fMu[0].GetObject())->MatchTriggerAnyPt());
}

//______________________________________________________________________________
Int_t AliAODDimuon::LowPt(){
  // Test if the two muons match two trigger tracks with a "Low Pt" cut
  if(this->CheckPointers())return 0;
  return (((AliAODTrack*)fMu[0].GetObject())->MatchTriggerLowPt())&&
         (((AliAODTrack*)fMu[0].GetObject())->MatchTriggerLowPt());
}

//______________________________________________________________________________
Int_t AliAODDimuon::HighPt(){
  // Test if the two muons match two trigger tracks with a "High Pt" cut
  if(this->CheckPointers())return 0;
  return (((AliAODTrack*)fMu[0].GetObject())->MatchTriggerHighPt())&&
         (((AliAODTrack*)fMu[0].GetObject())->MatchTriggerHighPt());
}

//______________________________________________________________________________
Double_t AliAODDimuon::MaxChi2Match(){
  // Maximum matching Chi2 between track and trigger track
  if(this->CheckPointers())return -999999999;
  return TMath::Max((((AliAODTrack*)fMu[0].GetObject())->GetChi2MatchTrigger()),
                    (((AliAODTrack*)fMu[0].GetObject())->GetChi2MatchTrigger()));
}
