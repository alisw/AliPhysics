/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// Class to generate particles from using paramtrized pT and y distributions.
// Distributions are obtained from pointer to object of type AliGenLib.
// (For example AliGenMUONlib)
// Decays are performed using Pythia.
// andreas.morsch@cern.ch

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TROOT.h>
#include <TVirtualMC.h>

#include "AliDecayer.h"
#include "AliGenMUONlib.h"
#include "AliGenParam.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliGenEventHeader.h"

ClassImp(AliGenParam)

//------------------------------------------------------------

//Begin_Html
/*
  <img src="picts/AliGenParam.gif">
*/
//End_Html

//____________________________________________________________
AliGenParam::AliGenParam()
: fPtParaFunc(0),
  fYParaFunc(0),
  fIpParaFunc(0),
  fV2ParaFunc(0),
  fPtPara(0),
  fYPara(0),
  fV2Para(0),
  fdNdPhi(0),
  fParam(0),
  fdNdy0(0.),
  fYWgt(0.),
  fPtWgt(0.),
  fBias(0.),
  fTrials(0),
  fDeltaPt(0.01),
  fSelectAll(kFALSE),
  fDecayer(0),
  fForceConv(kFALSE),
  fKeepParent(kFALSE),
  fKeepIfOneChildSelected(kFALSE),
  fPreserveFullDecayChain(kFALSE)
{
  // Default constructor
}
//____________________________________________________________
AliGenParam::AliGenParam(Int_t npart, const AliGenLib * Library,  Int_t param, const char* tname)
  :AliGenMC(npart),
   fPtParaFunc(Library->GetPt(param, tname)),
   fYParaFunc (Library->GetY (param, tname)),
   fIpParaFunc(Library->GetIp(param, tname)),
   fV2ParaFunc(Library->GetV2(param, tname)),
   fPtPara(0),
   fYPara(0),
   fV2Para(0),
   fdNdPhi(0),
   fParam(param),
   fdNdy0(0.),
   fYWgt(0.),
   fPtWgt(0.),
   fBias(0.),
   fTrials(0),
   fDeltaPt(0.01),
   fSelectAll(kFALSE),
   fDecayer(0),
   fForceConv(kFALSE),
   fKeepParent(kFALSE),
   fKeepIfOneChildSelected(kFALSE),
   fPreserveFullDecayChain(kFALSE)
{
  // Constructor using number of particles parameterisation id and library
  fName     = "Param";
  fTitle    = "Particle Generator using pT and y parameterisation";
  fAnalog   = kAnalog;
  SetForceDecay();
}
//____________________________________________________________
AliGenParam::AliGenParam(Int_t npart, Int_t param, const char* tname, const char* name):
  AliGenMC(npart),
  fPtParaFunc(0),
  fYParaFunc (0),
  fIpParaFunc(0),
  fV2ParaFunc(0),
  fPtPara(0),
  fYPara(0),
  fV2Para(0),
  fdNdPhi(0),
  fParam(param),
  fdNdy0(0.),
  fYWgt(0.),
  fPtWgt(0.),
  fBias(0.),
  fTrials(0),
  fDeltaPt(0.01),
  fSelectAll(kFALSE),
  fDecayer(0),
  fForceConv(kFALSE),
  fKeepParent(kFALSE),
  fKeepIfOneChildSelected(kFALSE),
  fPreserveFullDecayChain(kFALSE)
{
  // Constructor using parameterisation id and number of particles
  //
  fName       = name;
  fTitle      = "Particle Generator using pT and y parameterisation";

  AliGenLib* pLibrary = new AliGenMUONlib();
  fPtParaFunc = pLibrary->GetPt(param, tname);
  fYParaFunc  = pLibrary->GetY (param, tname);
  fIpParaFunc = pLibrary->GetIp(param, tname);
  fV2ParaFunc = pLibrary->GetV2(param, tname);

  fAnalog     = kAnalog;
  fChildSelect.Set(5);
  for (Int_t i=0; i<5; i++) fChildSelect[i]=0;
  SetForceDecay();
  SetCutOnChild();
  SetChildMomentumRange();
  SetChildPtRange();
  SetChildPhiRange();
  SetChildThetaRange();
}
//____________________________________________________________

AliGenParam::AliGenParam(Int_t npart, Int_t param,
                         Double_t (*PtPara) (const Double_t*, const Double_t*),
                         Double_t (*YPara ) (const Double_t* ,const Double_t*),
                         Double_t (*V2Para) (const Double_t* ,const Double_t*),
		         Int_t    (*IpPara) (TRandom *))
  :AliGenMC(npart),

   fPtParaFunc(PtPara),
   fYParaFunc(YPara),
   fIpParaFunc(IpPara),
   fV2ParaFunc(V2Para),
   fPtPara(0),
   fYPara(0),
   fV2Para(0),
   fdNdPhi(0),
   fParam(param),
   fdNdy0(0.),
   fYWgt(0.),
   fPtWgt(0.),
   fBias(0.),
   fTrials(0),
   fDeltaPt(0.01),
   fSelectAll(kFALSE),
   fDecayer(0),
   fForceConv(kFALSE),
   fKeepParent(kFALSE),
   fKeepIfOneChildSelected(kFALSE)
{
  // Constructor
  // Gines Martinez 1/10/99
  fName   = "Param";
  fTitle  = "Particle Generator using pT and y parameterisation";

  fAnalog = kAnalog;
  fChildSelect.Set(5);
  for (Int_t i=0; i<5; i++) fChildSelect[i]=0;
  SetForceDecay();
  SetCutOnChild();
  SetChildMomentumRange();
  SetChildPtRange();
  SetChildPhiRange();
  SetChildThetaRange();
}

//____________________________________________________________
AliGenParam::~AliGenParam()
{
  // Destructor
  delete  fPtPara;
  delete  fYPara;
  delete  fV2Para;
  delete  fdNdPhi;
}

//-------------------------------------------------------------------
TVector3 AliGenParam::OrthogonalVector(TVector3 &inVec){
  double abc[]={inVec.x(), inVec.y(), inVec.z()};
  double xyz[]={1,1,1};
  int solvDim=0;
  double tmp=abc[0];
  for(int i=0; i<3; i++)
    if(fabs(abc[i])>tmp){
      solvDim=i;
      tmp=fabs(abc[i]);
    }
  xyz[solvDim]=(-abc[(1+solvDim)%3]-abc[(2+solvDim)%3])/abc[(0+solvDim)%3];

  TVector3 res(xyz[0],xyz[1],xyz[2]);
  return res;
}

void AliGenParam::RotateVector(Double_t *pin, Double_t *pout, Double_t costheta, Double_t sintheta,
			       Double_t cosphi, Double_t sinphi)
{
  // Perform rotation
  pout[0] = pin[0]*costheta*cosphi-pin[1]*sinphi+pin[2]*sintheta*cosphi;
  pout[1] = pin[0]*costheta*sinphi+pin[1]*cosphi+pin[2]*sintheta*sinphi;
  pout[2] = -1.0  * pin[0] * sintheta + pin[2] * costheta;
  return;
}

double AliGenParam::ScreenFunction1(double screenVariable){
  if(screenVariable>1)
    return 42.24 - 8.368 * log(screenVariable + 0.952);
  else
    return 42.392 - screenVariable * (7.796 - 1.961 * screenVariable);
}

double AliGenParam::ScreenFunction2(double screenVariable){
  if(screenVariable>1)
    return 42.24 - 8.368 * log(screenVariable + 0.952);
  else
    return 41.405 - screenVariable * (5.828 - 0.8945 * screenVariable);
}

double AliGenParam::RandomEnergyFraction(double Z, double photonEnergy){
  double aZ=Z/137.036;
  double epsilon ;
  double epsilon0Local = 0.000511 / photonEnergy ;

  // Do it fast if photon energy < 2. MeV
  if (photonEnergy < 0.002 ){
    epsilon = epsilon0Local + (0.5 - epsilon0Local) * fRandom->Rndm();
  } else {
    double fZ = 8*log(Z)/3;
    double fcZ=(aZ*aZ)*(1/(1+aZ*aZ)+0.20206-0.0368*aZ*aZ+0.0083*aZ*aZ*aZ);
    if (photonEnergy > 0.050) fZ += 8*fcZ;

    // Limits of the screening variable
    double screenFactor = 136. * epsilon0Local / pow (Z,1/3);
    double screenMax = exp ((42.24 - fZ)/8.368) - 0.952 ;
    double screenMin = std::min(4.*screenFactor,screenMax) ;

    // Limits of the energy sampling
    double epsilon1 = 0.5 - 0.5 * sqrt(1. - screenMin / screenMax) ;
    double epsilonMin = std::max(epsilon0Local,epsilon1);
    double epsilonRange = 0.5 - epsilonMin ;

    // Sample the energy rate of the created electron (or positron)
    double screen;
    double gReject ;

    double f10 = ScreenFunction1(screenMin) - fZ;
    double f20 = ScreenFunction2(screenMin) - fZ;
    double normF1 = std::max(f10 * epsilonRange * epsilonRange,0.);
    double normF2 = std::max(1.5 * f20,0.);

    do {
      if (normF1 / (normF1 + normF2) > fRandom->Rndm() ){
          epsilon = 0.5 - epsilonRange * pow(fRandom->Rndm(), 0.333333) ;
          screen = screenFactor / (epsilon * (1. - epsilon));
          gReject = (ScreenFunction1(screen) - fZ) / f10 ;
      } else {
          epsilon = epsilonMin + epsilonRange * fRandom->Rndm();
          screen = screenFactor / (epsilon * (1 - epsilon));
          gReject = (ScreenFunction2(screen) - fZ) / f20 ;
      }
    } while ( gReject < fRandom->Rndm() );
  }   //  End of epsilon sampling
  return epsilon;
}

double AliGenParam::RandomPolarAngle(){
  double u;
  const double a1 = 0.625;
  double a2 = 3. * a1;
  //  double d = 27. ;

  //  if (9. / (9. + d) > fRandom->Rndm())
  if (0.25 > fRandom->Rndm()) {
    u = - log(fRandom->Rndm() * fRandom->Rndm()) / a1 ;
  } else {
    u = - log(fRandom->Rndm() * fRandom->Rndm()) / a2 ;
  }
  return u*0.000511;
}

Double_t AliGenParam::RandomMass(Double_t mh){
  while(true){
    double y=fRandom->Rndm();
    double mee=2*0.000511*TMath::Power(2*0.000511/mh,-y); //inverse of the enveloping cumulative distribution
    double apxkw=2.0/3.0/137.036/TMath::Pi()/mee; //enveloping probability density
    double val=fRandom->Uniform(0,apxkw);
    double kw=apxkw*sqrt(1-4*0.000511*0.000511/mee/mee)*(1+2*0.000511*0.000511/mee/mee)*1*1*TMath::Power(1-mee*mee/mh/mh,3);
    if(val<kw)
      return mee;
  }
}

Int_t AliGenParam::VirtualGammaPairProduction(TClonesArray *particles, Int_t nPart)
{
  Int_t nPartNew=nPart;
  for(int iPart=0; iPart<nPart; iPart++){
    TParticle *gamma = (TParticle *) particles->At(iPart);
    if(gamma->GetPdgCode()!=220001) continue;
    if(gamma->Pt()<0.002941) continue;  //approximation of kw in AliGenEMlib is 0 below 0.002941
    double mass=RandomMass(gamma->Pt());

    // lepton pair kinematics in virtual photon rest frame
    double Ee=mass/2;
    double Pe=TMath::Sqrt((Ee+0.000511)*(Ee-0.000511));

    double costheta = (2.0 * gRandom->Rndm()) - 1.;
    double sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
    double phi      = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
    double sinphi   = TMath::Sin(phi);
    double cosphi   = TMath::Cos(phi);

    // momentum vectors of leptons in virtual photon rest frame
    Double_t pProd1[3] = {Pe * sintheta * cosphi,
                          Pe * sintheta * sinphi,
                          Pe * costheta};

    Double_t pProd2[3] = {-1.0 * Pe * sintheta * cosphi,
                          -1.0 * Pe * sintheta * sinphi,
                          -1.0 * Pe * costheta};

    // lepton 4-vectors in properly rotated virtual photon rest frame
    Double_t pRot1[3] = {0.};
    RotateVector(pProd1, pRot1, costheta, -sintheta, -cosphi, -sinphi);
    Double_t pRot2[3] = {0.};
    RotateVector(pProd2, pRot2, costheta, -sintheta, -cosphi, -sinphi);

    TLorentzVector e1V4(pRot1[0],pRot1[1],pRot1[2],Ee);
    TLorentzVector e2V4(pRot2[0],pRot2[1],pRot2[2],Ee);

    TVector3 boost(gamma->Px(),gamma->Py(),gamma->Pz());
    boost*=1/sqrt(gamma->P()*gamma->P()+mass*mass);
    e1V4.Boost(boost);
    e2V4.Boost(boost);

    TLorentzVector vtx;
    gamma->ProductionVertex(vtx);
    new((*particles)[nPartNew]) TParticle(11, gamma->GetStatusCode(), iPart+1, -1, 0, 0, e1V4, vtx);
    nPartNew++;
    new((*particles)[nPartNew]) TParticle(-11, gamma->GetStatusCode(), iPart+1, -1, 0, 0, e2V4, vtx);
    nPartNew++;
  }
  return nPartNew;
}

Int_t AliGenParam::ForceGammaConversion(TClonesArray *particles, Int_t nPart)
{
  //based on: http://geant4.cern.ch/G4UsersDocuments/UsersGuides/PhysicsReferenceManual/html/node27.html
  //     and: http://geant4.cern.ch/G4UsersDocuments/UsersGuides/PhysicsReferenceManual/html/node58.html
  //     and: G4LivermoreGammaConversionModel.cc
  Int_t nPartNew=nPart;
  for(int iPart=0; iPart<nPart; iPart++){
    TParticle *gamma = (TParticle *) particles->At(iPart);
    if(gamma->GetPdgCode()!=22 & gamma->GetPdgCode()!=220000) continue;
    if(gamma->Energy()<=0.001022) continue;
    TVector3 gammaV3(gamma->Px(),gamma->Py(),gamma->Pz());
    double frac=RandomEnergyFraction(1,gamma->Energy());
    double Ee1=frac*gamma->Energy();
    double Ee2=(1-frac)*gamma->Energy();
    double Pe1=sqrt((Ee1+0.000511)*(Ee1-0.000511));
    double Pe2=sqrt((Ee2+0.000511)*(Ee2-0.000511));

    TVector3 rotAxis(OrthogonalVector(gammaV3));
    Float_t az=fRandom->Uniform(TMath::Pi()*2);
    rotAxis.Rotate(az,gammaV3);
    TVector3 e1V3(gammaV3);
    double u=RandomPolarAngle();
    e1V3.Rotate(u/Ee1,rotAxis);
    e1V3=e1V3.Unit();
    e1V3*=Pe1;
    TVector3 e2V3(gammaV3);
    e2V3.Rotate(-u/Ee2,rotAxis);
    e2V3=e2V3.Unit();
    e2V3*=Pe2;
    // gamma = new TParticle(*gamma);
    // particles->RemoveAt(iPart);
    gamma->SetFirstDaughter(nPartNew+1);
    gamma->SetLastDaughter(nPartNew+2);
    // new((*particles)[iPart]) TParticle(*gamma);
    // delete gamma;

    // conversion probability per atom
    // fitted G4EMLOW6.35/pair/pp-cs-8.dat, fit is great for E>20MeV
    double convProb = 1/(exp(-log(28.44*(gamma->Energy()-0.001022))*(0.775+0.0271*log(gamma->Energy()+1)))+1);

    // radiation length is not considered here, so you have to normalize yourself in after-production from infinite radiation length to whatever you want
    // double meanExessPathlength=0.5*(exp(0.9)-exp(-0.9))/0.9; double scale=(1-exp(-7.0/9.0*radLength*meanExessPathlength))/(1-0);

    TLorentzVector vtx;
    gamma->ProductionVertex(vtx);
    TParticle *currPart;
    Int_t sign = (gRandom->Rndm() < 0.5) ? 1 : -1;
    currPart = new((*particles)[nPartNew]) TParticle(sign * 220011, gamma->GetStatusCode(), iPart+1, -1, 0, 0, TLorentzVector(e1V3,Ee1), vtx);
    currPart->SetWeight(convProb);
    nPartNew++;
    currPart = new((*particles)[nPartNew]) TParticle(- sign * 220011, gamma->GetStatusCode(), iPart+1, -1, 0, 0, TLorentzVector(e2V3,Ee2), vtx);
    currPart->SetWeight(convProb);
    nPartNew++;
  }
  return nPartNew;
}

//____________________________________________________________
void AliGenParam::Init()
{
  // Initialisation

  if (TVirtualMC::GetMC()) fDecayer = TVirtualMC::GetMC()->GetDecayer();
  //Begin_Html
  /*
    <img src="picts/AliGenParam.gif">
  */
  //End_Html
  char name[256];
  snprintf(name, 256, "pt-parameterisation for %s", GetName());

  if (fPtPara) fPtPara->Delete();
  fPtPara = new TF1(name, fPtParaFunc, fPtMin, fPtMax,0);
  gROOT->GetListOfFunctions()->Remove(fPtPara);
  //  Set representation precision to 10 MeV
  Int_t npx= Int_t((fPtMax - fPtMin) / fDeltaPt);

  fPtPara->SetNpx(npx);

  snprintf(name, 256, "y-parameterisation  for %s", GetName());
  if (fYPara) fYPara->Delete();
  fYPara  = new TF1(name, fYParaFunc, fYMin, fYMax, 0);
  gROOT->GetListOfFunctions()->Remove(fYPara);

  snprintf(name, 256, "v2-parameterisation for %s", GetName());
  if (fV2Para) fV2Para->Delete();
  fV2Para  = new TF1(name, fV2ParaFunc, fPtMin, fPtMax, 0);
  // fV2Para  = new TF1(name, "2*[0]/(1+TMath::Exp([1]*([2]-x)))-[0]", fPtMin, fPtMax);
  // fV2Para->SetParameter(0, 0.236910);
  // fV2Para->SetParameter(1, 1.71122);
  // fV2Para->SetParameter(2, 0.0827617);
  //gROOT->GetListOfFunctions()->Remove(fV2Para);  //TR: necessary?

  snprintf(name, 256, "dNdPhi for %s", GetName());
  if (fdNdPhi) fdNdPhi->Delete();
  fdNdPhi  = new TF1(name, "1+2*[0]*TMath::Cos(2*(x-[1]))", fPhiMin, fPhiMax);
  //gROOT->GetListOfFunctions()->Remove(fdNdPhi);  //TR: necessary?

  snprintf(name, 256, "pt-for-%s", GetName());
  TF1 ptPara(name ,fPtParaFunc, 0, 15, 0);
  snprintf(name, 256, "y-for-%s", GetName());
  TF1 yPara(name, fYParaFunc, -6, 6, 0);

  //
  // dN/dy| y=0
  Double_t y1=0;
  Double_t y2=0;

  fdNdy0=fYParaFunc(&y1,&y2);
  //
  // Integral over generation region
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
  Float_t intYS  = yPara.Integral(fYMin, fYMax,(Double_t*) 0x0,1.e-6);
  Float_t intPt0 = ptPara.Integral(0,15,(Double_t *) 0x0,1.e-6);
  Float_t intPtS = ptPara.Integral(fPtMin,fPtMax,(Double_t*) 0x0,1.e-6);
#else
  Float_t intYS  = yPara.Integral(fYMin, fYMax,1.e-6);
  Float_t intPt0 = ptPara.Integral(0,15,1.e-6);
  Float_t intPtS = ptPara.Integral(fPtMin,fPtMax,1.e-6);
#endif
  Float_t phiWgt=(fPhiMax-fPhiMin)/2./TMath::Pi();    //TR: should probably be done differently in case of anisotropic phi...
  if (fAnalog == kAnalog) {
    fYWgt  = intYS/fdNdy0;
    fPtWgt = intPtS/intPt0;
    fParentWeight = fYWgt*fPtWgt*phiWgt/fNpart;
  } else {
    fYWgt = intYS/fdNdy0;
    fPtWgt = (fPtMax-fPtMin)/intPt0;
    fParentWeight = fYWgt*fPtWgt*phiWgt/fNpart;
  }
  //
  // particle decay related initialization
  fDecayer->SetForceDecay(fForceDecay);
  fDecayer->Init();

  //
  AliGenMC::Init();
}

//____________________________________________________________
void AliGenParam::Generate()
{
  //
  // Generate 1 event (see Generate(Int_t ntimes) for details
  //
  GenerateN(1);
}
//____________________________________________________________
void AliGenParam::GenerateN(Int_t ntimes)
{
  //
  // Generate ntimes*'npart' light and heavy mesons (J/Psi, upsilon or phi, Pion,
  // Kaons, Etas, Omegas) and Baryons (proton, antiprotons, neutrons and
  // antineutrons in the the desired theta, phi and momentum windows;
  // Gaussian smearing on the vertex is done if selected.
  // The decay of heavy mesons is done using lujet,
  //    and the childern particle are tracked by GEANT
  // However, light mesons are directly tracked by GEANT
  // setting fForceDecay = nodecay (SetForceDecay(nodecay))
  //
  //
  //  Reinitialize decayer
  fDecayer->SetForceDecay(fForceDecay);
  fDecayer->Init();

  //
  Float_t polar[3]= {0,0,0};  // Polarisation of the parent particle (for GEANT tracking)
  Double_t origin0[3];         // Origin of the generated parent particle (for GEANT tracking)
  Double_t time0;              // Time0 of the generated parent particle
  Double_t pt, pl, ptot;       // Transverse, logitudinal and total momenta of the parent particle
  Double_t phi, theta;         // Phi and theta spherical angles of the parent particle momentum
  Double_t p[3], pc[3], och[3];// Momentum, polarisation and origin of the children particles from lujet
  Double_t ty, xmt;
  Int_t nt, i, j;
  Double_t energy;
  Float_t  wgtp, wgtch;
  Double_t dummy;
  static TClonesArray *particles;
  //
  if(!particles) particles = new TClonesArray("TParticle",1000);

  TDatabasePDG *pDataBase = TDatabasePDG::Instance();
  //
  Float_t random[6];

  // Calculating vertex position per event
  for (j=0;j<3;j++) origin0[j]=fOrigin[j];
  time0 = fTimeOrigin;
  if(fVertexSmear==kPerEvent) {
    Vertex();
    for (j=0;j<3;j++) origin0[j]=fVertex[j];
    time0 = fTime;
  }

  Int_t ipa=0;

  // Generating fNpart particles
  fNprimaries = 0;

  Int_t nGen = fNpart*ntimes;
  while (ipa<nGen) {
    while(1) {
      //
      // particle type
      Int_t iPart = fIpParaFunc(fRandom);
      Int_t iTemp = iPart;

      // custom pdg codes to destinguish direct photons
      if(iPart>=220000 & iPart<=220001) iPart=22;

      fChildWeight=(fDecayer->GetPartialBranchingRatio(iPart))*fParentWeight;
      TParticlePDG *particle = pDataBase->GetParticle(iPart);
      Float_t am = particle->Mass();

      Rndm(random,2);

      // --- For Exodus -------------------------------
      Double_t awidth = particle->Width();
      if(awidth>0){
        TF1 rbw("rbw","pow([1],2)*pow([0],2)/(pow(x*x-[0]*[0],2)+pow(x*x*[1]/[0],2))",am-5*awidth,am+5*awidth);
        rbw.SetParameter(0,am);
        rbw.SetParameter(1,awidth);
        am = rbw.GetRandom();
      }
      // -----------------------------------------------//

      //
      // y
      ty = TMath::TanH(fYPara->GetRandom());

      //
      // pT
      if (fAnalog == kAnalog) {
        pt=fPtPara->GetRandom();
        wgtp=fParentWeight;
        wgtch=fChildWeight;
      } else {
        pt=fPtMin+random[1]*(fPtMax-fPtMin);
        Double_t ptd=pt;
        wgtp=fParentWeight*fPtParaFunc(& ptd, &dummy);
        wgtch=fChildWeight*fPtParaFunc(& ptd, &dummy);
      }
      xmt=sqrt(pt*pt+am*am);
      if (TMath::Abs(ty)==1.) {
        ty=0.;
        Fatal("AliGenParam",
              "Division by 0: Please check you rapidity range !");
      }
      //
      // phi
      //      if(!ipa)
      //phi=fEvPlane; //align first particle of each event with event plane
      //else{
      double v2 = fV2Para->Eval(pt);
      fdNdPhi->SetParameter(0,v2);
      fdNdPhi->SetParameter(1,fEvPlane);
      phi=fdNdPhi->GetRandom();
      //     }

      pl=xmt*ty/sqrt((1.-ty)*(1.+ty));
      theta=TMath::ATan2(pt,pl);
      // Cut on theta
      if(theta<fThetaMin || theta>fThetaMax) continue;
      ptot=TMath::Sqrt(pt*pt+pl*pl);
      // Cut on momentum
      if(ptot<fPMin || ptot>fPMax) continue;
      //
      p[0]=pt*TMath::Cos(phi);
      p[1]=pt*TMath::Sin(phi);
      p[2]=pl;
      energy=TMath::Sqrt(ptot*ptot+am*am);

      if(fVertexSmear==kPerTrack) {
        Rndm(random,6);
        for (j=0;j<3;j++) {
          origin0[j]= fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
                      TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
        }
        Rndm(random,2);
        time0 = fTimeOrigin + fOsigma[2]/TMath::Ccgs()*
                TMath::Cos(2*random[0]*TMath::Pi())*
                TMath::Sqrt(-2*TMath::Log(random[1]));
      }

      // Looking at fForceDecay :
      // if fForceDecay != none Primary particle decays using
      // AliPythia and children are tracked by GEANT
      //
      // if fForceDecay == none Primary particle is tracked by GEANT
      // (In the latest, make sure that GEANT actually does all the decays you want)
      //
      Bool_t decayed = kFALSE;


      if (fForceDecay != kNoDecay) {
        // Using lujet to decay particle
        TLorentzVector pmom(p[0], p[1], p[2], energy);
        fDecayer->Decay(iPart,&pmom);
        //
        // select decay particles
        Int_t np=fDecayer->ImportParticles(particles);

        iPart=iTemp;
        if(iPart>=220000 & iPart<=220001){
          TParticle *gamma = (TParticle *)particles->At(0);
          gamma->SetPdgCode(iPart);
          np=VirtualGammaPairProduction(particles,np);
        }
        if(fForceConv) np=ForceGammaConversion(particles,np);

        //  Selecting  GeometryAcceptance for particles fPdgCodeParticleforAcceptanceCut;
        if (fGeometryAcceptance)
          if (!CheckAcceptanceGeometry(np,particles)) continue;
        Int_t ncsel=0;
        Int_t pFlag    [np];
        Int_t pParent  [np];
        Int_t pSelected[np];
        Int_t trackIt  [np];

        for (i=0; i<np; i++) {
          pFlag[i]     =  0;
          pSelected[i] =  0;
          pParent[i]   = -1;
        }

        if (np >1) {
          decayed = kTRUE;
          TParticle* iparticle =  0;
          Int_t ipF, ipL;
          for (i = 1; i<np ; i++) {
            trackIt[i] = 1;
            iparticle = (TParticle *) particles->At(i);
            Int_t kf = iparticle->GetPdgCode();
            Int_t ks = iparticle->GetStatusCode();
            // flagged particle

            if(!fPreserveFullDecayChain){
              if (pFlag[i] == 1) {
                ipF = iparticle->GetFirstDaughter();
                ipL = iparticle->GetLastDaughter();
                if (ipF > 0) for (j=ipF-1; j<ipL; j++) pFlag[j]=1;
                continue;
              }
            }
            // flag decay products of particles with long life-time (c tau > .3 mum)

            if (ks != 1) {
              //  TParticlePDG *particle = pDataBase->GetParticle(kf);

              Double_t lifeTime = fDecayer->GetLifetime(kf);
              //  Double_t mass     = particle->Mass();
              //  Double_t width    = particle->Width();
              if (lifeTime > (Double_t) fMaxLifeTime) {
                ipF = iparticle->GetFirstDaughter();
                ipL = iparticle->GetLastDaughter();
                if (ipF > 0) for (j=ipF-1; j<ipL; j++) pFlag[j]=1;
                    } else{
                trackIt[i]     = 0;
                pSelected[i]   = 1;
              }
            } // ks==1 ?
            //
            // children

            if ((ChildSelected(TMath::Abs(kf)) || fForceDecay == kAll || fSelectAll) && trackIt[i]){
              if (fCutOnChild) {
                pc[0]=iparticle->Px();
                pc[1]=iparticle->Py();
                pc[2]=iparticle->Pz();
                Bool_t  childok = KinematicSelection(iparticle, 1);
                if(childok) {
                  pSelected[i]  = 1;
                  ncsel++;
                } else {
                  if(!fKeepIfOneChildSelected){
                    ncsel=-1;
                    break;
                  }
                } // child kine cuts
              } else {
                pSelected[i]  = 1;
                ncsel++;
              } // if child selection
            } // select muon
          } // decay particle loop
        } // if decay products

        Int_t iparent;

        if (fKeepParent || (fCutOnChild && ncsel >0) || !fCutOnChild){
          //
          // Parent

          // --- For Exodus --------------------------------//
              //PushTrack(0, -1, iPart, p, origin0, polar, time0, kPPrimary, nt, wgtp, ((decayed)? 11 : 1));
          PushTrack(0, -1, iPart, p[0],p[1],p[2],energy,origin0[0],origin0[1],origin0[2],time0,polar[0],polar[1],polar[2],kPPrimary, nt, wgtp, ((decayed)? 11 : 1));
          // -----------------------------------------------//

          pParent[0] = nt;
          KeepTrack(nt);
          fNprimaries++;

          //but count is as "generated" particle" only if it produced child(s) within cut
          if ((fCutOnChild && ncsel >0) || !fCutOnChild) {
            ipa++;
          }

          //
          // Decay Products
          //
          for (i = 1; i < np; i++) {
            if (pSelected[i]) {
              TParticle* iparticle = (TParticle *) particles->At(i);
              Int_t kf   = iparticle->GetPdgCode();
              Int_t ksc  = iparticle->GetStatusCode();
              Int_t jpa  = iparticle->GetFirstMother()-1;
              Double_t weight = iparticle->GetWeight();
              och[0] = origin0[0]+iparticle->Vx();
              och[1] = origin0[1]+iparticle->Vy();
              och[2] = origin0[2]+iparticle->Vz();
              pc[0]  = iparticle->Px();
              pc[1]  = iparticle->Py();
              pc[2]  = iparticle->Pz();
              Double_t ec   = iparticle->Energy();

              if (jpa > -1) {
                iparent = pParent[jpa];
              } else {
                iparent = -1;
              }

              PushTrack(fTrackIt * trackIt[i], iparent, kf, pc[0], pc[1], pc[2], ec,
                        och[0], och[1], och[2], time0 + iparticle->T(),
                        polar[0], polar[1], polar[2], kPDecay, nt, weight*wgtch, ksc);

              //	      PushTrack(fTrackIt * trackIt[i], iparent, kf,
              //		pc, och, polar,
              //	time0 + iparticle->T(), kPDecay, nt, weight*wgtch, ksc);

              pParent[i] = nt;
              KeepTrack(nt);
              fNprimaries++;
            } // Selected
          } // Particle loop
        }  // Decays by Lujet
        particles->Clear();
      // kinematic selection
      } else { // nodecay option, so parent will be tracked by GEANT (pions, kaons, eta, omegas, baryons)
        PushTrack(fTrackIt, -1, iPart, p[0], p[1], p[2], energy,
                  origin0[0], origin0[1], origin0[2], time0,
                  polar[0], polar[1], polar[2],
                  kPPrimary, nt, wgtp, 1);

          // gAlice->GetMCApp()->
          // PushTrack(fTrackIt,-1,iPart,p,origin0,polar,time0,kPPrimary,nt,wgtp, 1);
        ipa++;
        fNprimaries++;
      }
      break;
    } // while
  } // event loop

  SetHighWaterMark(nt);

  AliGenEventHeader* header = new AliGenEventHeader("PARAM");
  header->SetPrimaryVertex(fVertex);
  header->SetInteractionTime(fTime);
  header->SetNProduced(fNprimaries);
  AddHeader(header);
}
//____________________________________________________________________________________
Float_t AliGenParam::GetRelativeArea(Float_t ptMin, Float_t ptMax, Float_t yMin, Float_t yMax, Float_t phiMin, Float_t phiMax)
{
  //
  // Normalisation for selected kinematic region
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(5,99,0)
  Float_t ratio =
    fPtPara->Integral(ptMin,ptMax,(Double_t *)0,1.e-6) / fPtPara->Integral( fPtPara->GetXmin(), fPtPara->GetXmax(),(Double_t *)0,1.e-6) *
    fYPara->Integral(yMin,yMax,(Double_t *)0,1.e-6)/fYPara->Integral(fYPara->GetXmin(),fYPara->GetXmax(),(Double_t *)0,1.e-6)   *
    (phiMax-phiMin)/360.;
#else
  Float_t ratio =
    fPtPara->Integral(ptMin,ptMax,1.e-6) / fPtPara->Integral( fPtPara->GetXmin(), fPtPara->GetXmax(),1.e-6) *
    fYPara->Integral(yMin,yMax,1.e-6)/fYPara->Integral(fYPara->GetXmin(),fYPara->GetXmax(),1.e-6)   *
    (phiMax-phiMin)/360.;
#endif
  return TMath::Abs(ratio);
}

//____________________________________________________________________________________

void AliGenParam::Draw( const char * /*opt*/)
{
  //
  // Draw the pT and y Distributions
  //
  TCanvas *c0 = new TCanvas("c0","Canvas 0",400,10,600,700);
  c0->Divide(2,1);
  c0->cd(1);
  fPtPara->Draw();
  fPtPara->GetHistogram()->SetXTitle("p_{T} (GeV)");
  c0->cd(2);
  fYPara->Draw();
  fYPara->GetHistogram()->SetXTitle("y");
}
