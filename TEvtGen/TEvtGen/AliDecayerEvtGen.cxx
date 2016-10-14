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
///////////////////////////////////////////////////////////////////////////
//  Implementation of AliDecayer using EvtGen package.                   //
//                                                                       //
//  Giuseppe E. Bruno            &    Fiorella Fionda                    //
//  (Giuseppe.Bruno@ba.infn.it)       (Fiorella.Fionda@ba.infn.it)       //
///////////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "EvtGenBase/EvtStdHep.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtStdlibRandomEngine.hh" 
#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "AliDecayerEvtGen.h"
#include "EvtGenExternal/EvtExternalGenList.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "AliLog.h"

ClassImp(AliDecayerEvtGen)
//____________________________________________________________
AliDecayerEvtGen::AliDecayerEvtGen():
  fRandomEngine(0x0),
  fRadCorrEngine(0x0),
  fGenerator(0x0),
  fEvtstdhep(0x0),
  fDecayTablePath(0x0),
  fParticleTablePath(0x0),
  fDecay(kAll)
  {
  // Default constructor
  fEvtstdhep = new EvtStdHep();
  fDecayTablePath = gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DECAY.DEC"); //default decay table
  fParticleTablePath = gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/evt.pdl"); //particle table
  }
//_________________________________________________________________
AliDecayerEvtGen::AliDecayerEvtGen(const AliDecayerEvtGen &decayer):
  AliDecayer(decayer),
  fRandomEngine(decayer.fRandomEngine),
  fGenerator(decayer.fGenerator),
  fEvtstdhep(decayer.fEvtstdhep),
  fDecayTablePath(decayer.fDecayTablePath),
  fParticleTablePath(decayer.fParticleTablePath),
  fDecay(decayer.fDecay)
  {
  // Copy Constructor
  decayer.Copy(*this);
  }
//____________________________________________________________
AliDecayerEvtGen::~AliDecayerEvtGen()
  {
  // Destructor
  if(fRandomEngine) {delete fRandomEngine;}
  fRandomEngine = 0;
  if(fRadCorrEngine) {delete fRadCorrEngine;}
  fRadCorrEngine = 0;
  if(fGenerator) {delete fGenerator;}
  fGenerator = 0;
  if(fEvtstdhep) {delete fEvtstdhep;}
  fEvtstdhep = 0;
  if(fDecayTablePath) {delete fDecayTablePath;}
  fDecayTablePath = 0;
  if(fParticleTablePath) {delete fParticleTablePath;}
  fParticleTablePath = 0;
  }

//___________________________________________________________
void AliDecayerEvtGen::Init()
  {
  //Standard AliDecayerEvtGen initializer:
  //initialize EvtGen with default decay table (DECAY.DEC), particle table (evt.pdl)
  //and fRandomEngine for generation of random numbers
  //
  if(fGenerator){
  AliWarning(" AliDecayerEvtGen already initialized!!!!\n");
  return;
  }
  fRandomEngine = new EvtStdlibRandomEngine();
  std::list<EvtDecayBase*> extraModels;

  EvtExternalGenList genList;
  fRadCorrEngine = genList.getPhotosModel();
  extraModels = genList.getListOfModels();
  
  fGenerator=new EvtGen(fDecayTablePath,fParticleTablePath,fRandomEngine,fRadCorrEngine,&extraModels);
  }
//____________________________________________________________
void AliDecayerEvtGen::Decay(Int_t ipart, TLorentzVector *p)
  {
  //
  //Decay a particle
  //input: pdg code and momentum of the particle to be decayed  
  //all informations about decay products are stored in fEvtstdhep 
  //
  EvtId IPART=EvtPDL::evtIdFromStdHep(ipart);
  EvtVector4R p_init(p->E(),p->Px(),p->Py(),p->Pz());
  EvtParticle *froot_part=EvtParticleFactory::particleFactory(IPART,p_init);
  fGenerator->generateDecay(froot_part);
  fEvtstdhep->init();
  froot_part->makeStdHep(*fEvtstdhep);
  if(AliLog::GetDebugLevel("TEvtGen","AliDecayerEvtGen") > 0) 
    froot_part->printTree(); //to print the decay chain 
  froot_part->deleteTree();
  }

//____________________________________________________________
Int_t AliDecayerEvtGen::ImportParticles(TClonesArray *particles)
  {
  //
  //Input: pointer to a TClonesArray - Output(Int_t): number of decay products    
  //Put all the informations about the decay products in the 
  //TClonesArray particles
  //
  if (particles == 0) return 0;
  TClonesArray &clonesParticles = *particles;
  clonesParticles.Clear();

  int j;
  int istat;
  int partnum;
  double px,py,pz,e;
  double x,y,z,t;
  EvtVector4R p4,x4;

  Int_t npart=fEvtstdhep->getNPart();
  for(int i=0;i<fEvtstdhep->getNPart();i++){
  j=i+1;
  int jmotherfirst=fEvtstdhep->getFirstMother(i)+1;
  int jmotherlast=fEvtstdhep->getLastMother(i)+1;
  int jdaugfirst=fEvtstdhep->getFirstDaughter(i)+1;
  int jdauglast=fEvtstdhep->getLastDaughter(i)+1;

  partnum=fEvtstdhep->getStdHepID(i);
  
  //verify if all particles of decay chain are in the TDatabasePDG
  TParticlePDG *partPDG = TDatabasePDG::Instance()->GetParticle(partnum);
  if(!partPDG)
    {
    AliWarning("Particle code non known in TDatabasePDG - set pdg = 89");
    partnum=89; //internal use for unspecified resonance data
    }

  istat=fEvtstdhep->getIStat(i);

  if(istat!=1 && istat!=2) Info("ImportParticles","Attention: unknown status code!");
  if(istat == 2) istat = 11; //status decayed

  p4=fEvtstdhep->getP4(i);
  x4=fEvtstdhep->getX4(i);
  px=p4.get(1);
  py=p4.get(2);
  pz=p4.get(3);
  e=p4.get(0);
  const Float_t kconvT=0.001/2.999792458e8; // mm/c to seconds conversion
  const Float_t kconvL=1./10; // mm to cm conversion
  x=x4.get(1)*kconvL;//[cm]
  y=x4.get(2)*kconvL;//[cm]
  z=x4.get(3)*kconvL;//[cm]
  t=x4.get(0)*kconvT;//[s]

  AliDebug(1,Form("partnum = %d istat = %d primaMadre = %d ultimaMadre = %d primaF = %d ultimaF=%d x=%f y=%f z=%f t=%f e=%f px=%f \n",partnum,istat,jmotherfirst,jmotherlast,jdaugfirst,jdauglast,x,y,z,t,e,px));

  new(clonesParticles[i]) TParticle(partnum,istat,jmotherfirst,-1,jdaugfirst,jdauglast,px,py,pz,e,x,y,z,t);

  //set polarization!!!
  }

  return npart; 

  }

void  AliDecayerEvtGen::Copy(TObject &) const
  {
  //
  // Copy *this onto AliDecayerEvtGen -- not implemented
  //
   Fatal("Copy","Not implemented!\n");
  }

void AliDecayerEvtGen::ForceDecay()
  {
  //
  // Intupt: none - Output: none
  // Set the decay mode to decay particles: for each case is read a 
  // different decay table. case kAll read the default decay table only   
  //
  Decay_t decay = fDecay;
  switch(decay)
    {
     case kAll: // particles decayed "naturally" according to $ALICE_ROOT/TEvtGen/EvtGen/DECAY.DEC
      break;
     case kBJpsiDiElectron:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BTOJPSITOELE.DEC"));
      break;
     case kBJpsi:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BTOJPSI.DEC"));
      break;
     case kBJpsiDiMuon:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BTOJPSITOMU.DEC"));
      break;
     case kBSemiElectronic:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BTOELE.DEC"));
      break;
     case kHadronicD:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/HADRONICD.DEC"));
      break;
     case kHadronicDWithout4Bodies:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/HADRONICDWITHOUT4BODIES.DEC"));
      break;
     case kChiToJpsiGammaToElectronElectron:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/CHICTOJPSITOELE.DEC"));
      break;
     case kChiToJpsiGammaToMuonMuon:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/CHICTOJPSITOMUON.DEC"));
      break;
     case kSemiElectronic:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BANDCTOELE.DEC"));
      break;
     case kBSemiMuonic:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BTOMU.DEC"));
      break;
     case kSemiMuonic:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BANDCTOMU.DEC"));
      break;
     case kDiElectron:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/DIELECTRON.DEC"));
      break;
     case kDiMuon:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/DIMUON.DEC"));
      break;
     case kBPsiPrimeDiMuon:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BTOPSIPRIMETODIMUON.DEC"));
      break;
     case kBPsiPrimeDiElectron:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BTOPSIPRIMETODIELECTRON.DEC"));
      break;
     case kJpsiDiMuon:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/JPSIDIMUON.DEC"));
      break;
     case kPsiPrimeJpsiDiElectron:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/PSIPRIMETOJPSITOMU.DEC"));
      break;
     case kPhiKK:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/PHITOK.DEC"));
      break;
     case kOmega:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/OMEGATOLAMBDAK.DEC"));
      break;
     case kLambda:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/LAMBDATOPROTPI.DEC"));
      break;
     case kHardMuons:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/HARDMUONS.DEC"));
      break;
     case kElectronEM:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/ELECTRONEM.DEC"));
      break;
     case kDiElectronEM:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/DIELECTRONEM.DEC"));
      break;
     case kGammaEM:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/GAMMAEM.DEC"));
      break;
     case kBeautyUpgrade:
      SetDecayTablePath(gSystem->ExpandPathName("$ALICE_ROOT/TEvtGen/EvtGen/DecayTable/BEAUTYUPGRADE.DEC"));
      break;     
     case kPiToMu:
     case kKaToMu:
     case kAllMuonic:
     case kWToMuon:
     case kWToCharm:
     case kWToCharmToMuon:
     case kZDiMuon:
     case kZDiElectron:
     case kNoDecay:
     case kNoDecayHeavy:
     case kNeutralPion:
     case kBJpsiUndecayed:
     case kAllEM:
     case kNoDecayBeauty:
      AliWarning(Form("Warning: case %d not implemented for this class!",(int)decay));
     break;
     }
     ReadDecayTable();
  }

Float_t AliDecayerEvtGen::GetPartialBranchingRatio(Int_t)
  {
   // This method is dummy
   return  1.;
  }

Float_t AliDecayerEvtGen::GetLifetime(Int_t kf)
  {    
  //
  //Input: pdg code of a particle   
  //return lifetime in sec for a particle with particle code kf
  //
  EvtId IdPart=EvtPDL::evtIdFromStdHep(kf); 
  Double_t lifetime = EvtPDL::getctau(IdPart); //c*tau (mm)
  AliDebug(1,Form("lifetime is %f (mum) , particle id= %d",lifetime*1000,kf));
  return lifetime*kconv; //tau (sec)
  }

void AliDecayerEvtGen::ReadDecayTable()
    {
     //Input none - Output none 
     //Read the decay table that correspond to the path 
     //fDecayTablePath
     // 
     TString temp = fDecayTablePath;
     if(!temp.EndsWith("DECAY.DEC"))
     fGenerator->readUDecay(fDecayTablePath);
    }
////////////////////////////////////////////////////////////
Bool_t AliDecayerEvtGen::SetDecayTablePath(Char_t *path)
  {
  // 
  //Set the path of the decay table read to force particle decays 
  //
  if(gSystem->AccessPathName(path)) 
  {
   AliWarning("Attention: This path not exist!\n");
   return kFALSE;
  }
  fDecayTablePath = path;
  return kTRUE;
  } 



