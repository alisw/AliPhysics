// $Id$

//__________________________________________________________
///////////////////////////////////////////////////////////////////
//
// class AliJetEventParticles
//
// loizides@ikf.uni-frankfurt.de
//
///////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include "AliJetParticle.h"
#include "AliJetEventParticles.h"

ClassImp(AliJetEventParticles)

AliJetEventParticles::AliJetEventParticles(Int_t size) :
  TObject(), 
  fHeader(),
  fNParticles(0),
  fParticles(new TClonesArray("AliJetParticle",size)),
  fVertexX(0.),
  fVertexY(0.),
  fVertexZ(0.),
  fTrials(0),
  fNJets(0),
  fNUQJets(0),
  fXJet(-1),
  fYJet(-1),
  fImpact(0.),
  fNHardScatters(0),
  fNwNwColl(0),
  fEventNr(0)
{
  // Default Constructor
  for (Int_t i = 0; i < 4; i++) fZquench[i] = 0.;
  for (Int_t i = 0; i < 10; i++) 
    for (Int_t j = 0; j < 4; j++) {
      fJets[j][i]=0;    // Trigger jets
      fUQJets[j][i]=0;  // Unquenched trigger jets
    }
  for (Int_t i = 0; i < 5; i++){
    fHard[i][0]=0;
    fHard[i][1]=0;
  }
}

AliJetEventParticles::AliJetEventParticles(const AliJetEventParticles& source) :
  TObject(source), 
  fHeader(),
  fNParticles(source.fNParticles),
  fParticles(new TClonesArray("AliJetParticle",source.fNParticles)),
  fVertexX(source.GetVertexX()),
  fVertexY(source.GetVertexY()),
  fVertexZ(source.GetVertexZ()),
  fTrials(source.Trials()),
  fNJets(source.NTriggerJets()),
  fNUQJets(source.NUQTriggerJets()),
  fXJet(source.GetXJet()),
  fYJet(source.GetXJet()),
  fImpact(source.GetImpact()),
  fNHardScatters(source.GetNhard()),
  fNwNwColl(source.GetNpart()),
  fEventNr(source.GetEventNr())
{
  //copy constructor
  for(Int_t i =0; i<fNParticles; i++)
    {
      const AliJetParticle *kjp=(const AliJetParticle *)source.fParticles->At(i);
      new((*fParticles)[i]) AliJetParticle(*(kjp));
    }
  for (Int_t i = 0; i < 4; i++) fZquench[i] = 0.;
  for (Int_t i = 0; i < 10; i++) 
    for (Int_t j = 0; j < 4; j++) {
      fJets[j][i]=0;    // Trigger jets
      fUQJets[j][i]=0;  // Unquenched trigger jets
    }
  source.GetZQuench(fZquench);
  for (Int_t i = 0; i < NTriggerJets(); i++){
    source.TriggerJet(i,fJets[0][i],fJets[1][i],fJets[2][i],fJets[3][i]);
  }
  for (Int_t i = 0; i < NUQTriggerJets(); i++){
    source.UQJet(i,fUQJets[0][i],fUQJets[1][i],fUQJets[2][i],fUQJets[3][i]);
  }

  source.Hard(0,fHard[0][0],fHard[1][0],fHard[2][0],fHard[3][0],fHard[4][0]);
  source.Hard(1,fHard[0][1],fHard[1][1],fHard[2][1],fHard[3][1],fHard[4][1]);
}

void AliJetEventParticles::Set(const AliJetEventParticles& source)
{
  source.Copy(*this); 
  fHeader=source.GetHeader();
  fNParticles=source.fNParticles;
  if(fParticles) delete fParticles;
  fParticles=new TClonesArray("AliJetParticle",fNParticles);
  fVertexX=source.GetVertexX();
  fVertexY=source.GetVertexY();
  fVertexZ=source.GetVertexZ();
  fTrials=source.Trials();
  fNJets=source.NTriggerJets();
  fNUQJets=source.NUQTriggerJets();
  fXJet=source.GetXJet();
  fYJet=source.GetXJet();
  fImpact=source.GetImpact();
  fNHardScatters=source.GetNhard();
  fNwNwColl=source.GetNpart();
  fEventNr=source.GetEventNr();

  for(Int_t i =0; i<fNParticles; i++)
    {
      const AliJetParticle *kjp=(const AliJetParticle *)source.fParticles->At(i);
      new((*fParticles)[i]) AliJetParticle(*(kjp));
    }
  for (Int_t i = 0; i < 4; i++) fZquench[i] = 0.;
  for (Int_t i = 0; i < 10; i++) 
    for (Int_t j = 0; j < 4; j++) {
      fJets[j][i]=0;    // Trigger jets
      fUQJets[j][i]=0;  // Unquenched trigger jets
    }
  source.GetZQuench(fZquench);
  for (Int_t i = 0; i < NTriggerJets(); i++){
    source.TriggerJet(i,fJets[0][i],fJets[1][i],fJets[2][i],fJets[3][i]);
  }
  for (Int_t i = 0; i < NUQTriggerJets(); i++){
    source.UQJet(i,fUQJets[0][i],fUQJets[1][i],fUQJets[2][i],fUQJets[3][i]);
  }

  source.Hard(0,fHard[0][0],fHard[1][0],fHard[2][0],fHard[3][0],fHard[4][0]);
  source.Hard(1,fHard[0][1],fHard[1][1],fHard[2][1],fHard[3][1],fHard[4][1]);
}

AliJetEventParticles::~AliJetEventParticles()
{
  //destructor   
  Reset();
  delete fParticles;
}

void  AliJetEventParticles::Reset(Int_t size)
{
  fHeader="";
  //deletes all particles from the event
  if(fParticles) fParticles->Clear();
  if(size>=0) fParticles->Expand(size);
  fNParticles = 0;

  fVertexX=0.;
  fVertexY=0.;
  fVertexZ=0.;
  fTrials=0;
  fNJets=0;
  fNUQJets=0;
  fXJet=-1;
  fYJet=-1;
  fImpact=0.;
  fNHardScatters=0;
  fNwNwColl=0;
  fEventNr=0;
  for (Int_t i = 0; i < 4; i++) fZquench[i] = 0.;
  for (Int_t i = 0; i < 10; i++) 
    for (Int_t j = 0; j < 4; j++) {
      fJets[j][i]=0;    // Trigger jets
      fUQJets[j][i]=0;  // Unquenched trigger jets
    }
  for (Int_t i = 0; i < 5; i++){
    fHard[i][0]=0;
    fHard[i][1]=0;
  }
} 

void AliJetEventParticles::AddSignal(const AliJetEventParticles& source)
{ //mark signal particles and add them to TClonesArray
  //note that fNParticles still keeps only background particles

  Int_t nSignalParts=source.GetNParticles();
  for(Int_t i=0; i<nSignalParts; i++)
    {
      const AliJetParticle *kjp=source.GetParticle(i);

      AliJetParticle *ap=new((*fParticles)[fNParticles+i]) AliJetParticle(*(kjp));
      ap->SetType(-123); //mark pythia particle
    }
  for(Int_t i=nSignalParts+fNParticles;i<fParticles->GetEntriesFast();i++)
    fParticles->RemoveAt(i);
  //cout << fParticles->GetEntries() << " " << fNParticles << " " << nSignalParts << endl;

  /* should we transform the vertex???
  fVertexX=source.GetVertexX();
  fVertexY=source.GetVertexY();
  fVertexZ=source.GetVertexZ();
  */
  fTrials=source.Trials();
  fNJets=source.NTriggerJets();
  fNUQJets=source.NUQTriggerJets();
  fXJet=source.GetXJet();
  fYJet=source.GetXJet();
  fEventNr=source.GetEventNr();

  for (Int_t i = 0; i < 4; i++) fZquench[i] = 0.;
  for (Int_t i = 0; i < 10; i++) 
    for (Int_t j = 0; j < 4; j++) {
      fJets[j][i]=0;    // Trigger jets
      fUQJets[j][i]=0;  // Unquenched trigger jets
    }
  source.GetZQuench(fZquench);
  for (Int_t i = 0; i < NTriggerJets(); i++){
    source.TriggerJet(i,fJets[0][i],fJets[1][i],fJets[2][i],fJets[3][i]);
  }
  for (Int_t i = 0; i < NUQTriggerJets(); i++){
    source.UQJet(i,fUQJets[0][i],fUQJets[1][i],fUQJets[2][i],fUQJets[3][i]);
  }
  source.Hard(0,fHard[0][0],fHard[1][0],fHard[2][0],fHard[3][0],fHard[4][0]);
  source.Hard(1,fHard[0][1],fHard[1][1],fHard[2][1],fHard[3][1],fHard[4][1]);
}

void AliJetEventParticles::AddParticle(AliJetParticle* part)
{
  //Adds new particle to the event
  fParticles->AddAt(part,fNParticles++);
}

void AliJetEventParticles::AddParticle(const AliJetParticle* part)
{
  //Adds new particle to the event
  new((*fParticles)[fNParticles++]) AliJetParticle(*part);
}

void AliJetEventParticles::AddParticle(const TParticle* part,Int_t idx, Int_t l, Int_t ncl)
{
  //Adds new particle to the event
  new((*fParticles)[fNParticles++]) AliJetParticle(part,idx,l,ncl);
}

void AliJetEventParticles::AddParticle(Float_t px, Float_t py, Float_t pz, 
                              Float_t etot, Int_t idx, Int_t l, Int_t ncl)
{
  //Adds new particle to the event
  new((*fParticles)[fNParticles++]) AliJetParticle(px,py,pz,etot,idx,l,ncl); 
}

void AliJetEventParticles::AddParticle(Float_t px, Float_t py, Float_t pz, Float_t etot, Int_t idx, Int_t l,
		              Int_t ncl, Float_t pt, Float_t phi, Float_t eta)
{
  //Adds new particle to the event
  new((*fParticles)[fNParticles++]) AliJetParticle(px,py,pz,etot,idx,l,ncl,pt,phi,eta); 
}

const AliJetParticle* AliJetEventParticles::GetParticleSafely(Int_t n)
{
  //returns nth particle with range check
  if( (n<0) || (fNParticles<=n) ) return 0;
  return (const AliJetParticle*)fParticles->At(n);
}

void AliJetEventParticles::AddJet(Float_t px, Float_t py, Float_t pz, Float_t e)
{
  //
  //  Add a jet 
  //
  if (fNJets < 10) {
    fJets[0][fNJets] = px;
    fJets[1][fNJets] = py;
    fJets[2][fNJets] = pz;
    fJets[3][fNJets] = e;
    fNJets++;
  } else {
    printf("\nWarning: More than 10 jets triggered !!\n");
  }
}

void AliJetEventParticles::AddJet(Float_t p[4])
{
  //
  //  Add a jet 
  //
  if (fNJets < 10) {
    fJets[0][fNJets] = p[0];
    fJets[1][fNJets] = p[1];
    fJets[2][fNJets] = p[2];
    fJets[3][fNJets] = p[3];
    fNJets++;
  } else {
    printf("\nWarning: More than 10 jets triggered !!\n");
  }
}

void AliJetEventParticles::AddUQJet(Float_t px, Float_t py, Float_t pz, Float_t e)
{
  //
  //  Add a jet 
  //
  if (fNUQJets < 10) {
    fUQJets[0][fNUQJets] = px;
    fUQJets[1][fNUQJets] = py;
    fUQJets[2][fNUQJets] = pz;
    fUQJets[3][fNUQJets] = e;
    fNUQJets++;
  } else {
    printf("\nWarning: More than 10 jets triggered !!\n");
  }
}

void AliJetEventParticles::AddUQJet(Float_t p[4])
{
  //
  //  Add a jet 
  //
  if (fNUQJets < 10) {
    fUQJets[0][fNUQJets] = p[0];
    fUQJets[1][fNUQJets] = p[1];
    fUQJets[2][fNUQJets] = p[2];
    fUQJets[3][fNUQJets] = p[3];
    fNUQJets++;
  } else {
    printf("\nWarning: More than 10 jets triggered !!\n");
  }
}

void AliJetEventParticles::AddHard(Int_t i,Float_t px, Float_t py, Float_t pz, Float_t e, Float_t type)
{
  //
  //  Add a had parton
  //
  if (i < 2) {
    fHard[0][i] = px;
    fHard[1][i] = py;
    fHard[2][i] = pz;
    fHard[3][i] = e;
    fHard[4][i] = type;
  } else {
    printf("\nWarning: More than 2 partons !!\n");
  }
}

void AliJetEventParticles::SetZQuench(Double_t z[4])
{
  //
  // Set quenching fraction
  //
  for (Int_t i = 0; i < 4; i++) fZquench[i] = z[i];
}

void AliJetEventParticles::GetZQuench(Double_t z[4]) const
{
  //
  // Get quenching fraction
  //
  for (Int_t i = 0; i < 4; i++) z[i] = fZquench[i];
}

void AliJetEventParticles::TriggerJet(Int_t i, Float_t p[4]) const
{
  //
  // Give back jet #i
  //
  if (i >= fNJets) {
    printf("\nWarning: TriggerJet, index out of Range!!\n");
  } else {
    p[0] = fJets[0][i];
    p[1] = fJets[1][i];
    p[2] = fJets[2][i];
    p[3] = fJets[3][i];
  }
}

void AliJetEventParticles::TriggerJet(Int_t i, Float_t &p1, Float_t &p2, Float_t &p3, Float_t &E) const
{
  //
  // Give back jet #i
  //
  if (i >= fNJets) {
    printf("\nWarning: TriggerJet, index out of Range!!\n");
  } else {
    p1   = fJets[0][i];
    p2   = fJets[1][i];
    p3   = fJets[2][i];
    E    = fJets[3][i];
  }
}

void AliJetEventParticles::UQJet(Int_t i, Float_t p[4]) const
{
  //
  // Give back jet #i
  //
  if (i >= fNUQJets) {
    printf("\nWarning: Unquenched Jets, index out of Range!!\n");
  } else {
    p[0] = fUQJets[0][i];
    p[1] = fUQJets[1][i];
    p[2] = fUQJets[2][i];
    p[3] = fUQJets[3][i];
  }
}

void AliJetEventParticles::UQJet(Int_t i, Float_t &p1, Float_t &p2, Float_t &p3, Float_t &E) const
{
  //
  // Give back jet #i
  //
  if (i >= fNUQJets) {
    printf("\nWarning: Unquenched Jets, index out of Range!!\n");
  } else {
    p1   = fUQJets[0][i];
    p2   = fUQJets[1][i];
    p3   = fUQJets[2][i];
    E    = fUQJets[3][i];
  }
}

void AliJetEventParticles::Hard(Int_t i, Float_t &p1, Float_t &p2, Float_t &p3, Float_t &E, Float_t &type) const
{
  //
  // Give back jet #i
  //
  if (i >= 2) {
    printf("\nWarning: Hard partons, index out of Range!!\n");
  } else {
    p1   = fHard[0][i];
    p2   = fHard[1][i];
    p3   = fHard[2][i];
    E    = fHard[3][i];
    type = fHard[4][i];
  }
}

void AliJetEventParticles::Hard(Int_t i, Float_t p[4], Float_t &type) const
{
  //
  // Give back jet #i
  //
  if (i >= 2) {
    printf("\nWarning: Hard partons, index out of Range!!\n");
  } else {
    p[0]   = fHard[0][i];
    p[1]   = fHard[1][i];
    p[2]   = fHard[2][i];
    p[3]   = fHard[3][i];
    type = fHard[4][i];
  }
}

void AliJetEventParticles::SetXYJet(Double_t x, Double_t y)
{
  //
  //  Add jet production point
  //
  fXJet = x; 
  fYJet = y; 
}

void AliJetEventParticles::Print(Option_t* /*t*/) const
{
  cout << "--- AliJetEventParticles ---" << endl;
  if(fHeader.Length()) cout << fHeader.Data() << endl;
  cout << "Event Number: " << fEventNr << endl;
  cout << "Particles in Event: " << fNParticles << endl;
  if(fNUQJets){
    cout << "Unquenched Jets: " << fNUQJets << endl;
    for(Int_t i = 0; i<fNUQJets; i++){
      Float_t x=fUQJets[0][i];
      Float_t y=fUQJets[1][i];
      Float_t z=fUQJets[2][i];
      Float_t e=fUQJets[3][i];
      Float_t ptj=TMath::Sqrt(x*x+y*y);
      Float_t thj=TMath::ATan2(ptj,z);
      Float_t etaj=-TMath::Log(TMath::Tan(thj/2));
      Float_t phj=TMath::Pi()+TMath::ATan2(-y,-x);
      Float_t et=e*TMath::Sin(thj);
      cout << i << " " << et << " " << etaj << " " << phj << endl;
    }
  }
  if(fNJets){
    cout << "Triggered Jets: " << fNJets << endl;
    for(Int_t i = 0; i<fNJets; i++){
      Float_t x=fJets[0][i];
      Float_t y=fJets[1][i];
      Float_t z=fJets[2][i];
      Float_t e=fJets[3][i];
      Float_t ptj=TMath::Sqrt(x*x+y*y);
      Float_t thj=TMath::ATan2(ptj,z);
      Float_t etaj=-TMath::Log(TMath::Tan(thj/2));
      Float_t phj=TMath::Pi()+TMath::ATan2(-y,-x);
      Float_t et=e*TMath::Sin(thj);
      cout << i << " " << et << " " << etaj << " " << phj << endl;
    }
  }
}
