/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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


#include "AliAODPWG4Particle.h"

/// \cond CLASSIMP
ClassImp(AliAODPWG4Particle)
/// \endcond


//______________________________________________________________________________
///
/// Constructor.
///
AliAODPWG4Particle::AliAODPWG4Particle() :
AliVParticle(),
fMomentum(0),fPdg(-1), fTag(0), fLabel(-1),
fCaloLabel(), fTrackLabel(), fDetectorTag(-1),
fBadDist(0), fNLM(0), fM02(0), fM20(0),
fTime(0),fNCells(0),fSuperModule(0),
fDecayTag(0),fIsolated(0), fLeadingParticle(0),
fIsoConePtLead(), fIsoConeSumPt(),
fDisp(0), fTof(0), fCharged(0),
fTagged(0), fFidArea(0), fInputFileIndex(0),fBtag(0)
{
  fCaloLabel [0] = -1;
  fCaloLabel [1] = -1;
  fTrackLabel[0] = -1;
  fTrackLabel[1] = -1;
  fTrackLabel[2] = -1;
  fTrackLabel[3] = -1;
  
  fIsoConePtLead[0] = 0.;
  fIsoConeSumPt [0] = 0.;
  fIsoConePtLead[1] = 0.;
  fIsoConeSumPt [1] = 0.;
}

//______________________________________________________________________________
///
/// Constructor.
///
/// \param px particle momentum in x
/// \param py particle momentum in y
/// \param pz particle momentum in z
/// \param e particle energy
///
/// particle: cluster or track
///
AliAODPWG4Particle::AliAODPWG4Particle(Double_t px, Double_t py, Double_t pz, Double_t e):
  AliVParticle(),
  fMomentum(0),fPdg(-1), fTag(0), fLabel(-1),
  fCaloLabel(), fTrackLabel(), fDetectorTag(-1),
  fBadDist(0), fNLM(0), fM02(0), fM20(0),
  fTime(0),fNCells(0),fSuperModule(0),
  fDecayTag(0),fIsolated(0), fLeadingParticle(0),
  fIsoConePtLead(), fIsoConeSumPt(),
  fDisp(0), fTof(0), fCharged(0),
  fTagged(0), fFidArea(0), fInputFileIndex(0),fBtag(0)
{
  fMomentum = new TLorentzVector(px, py, pz, e);
  
  fCaloLabel [0] = -1;
  fCaloLabel [1] = -1;
  fTrackLabel[0] = -1;
  fTrackLabel[1] = -1;	
  fTrackLabel[2] = -1;
  fTrackLabel[3] = -1;	
  
  fIsoConePtLead[0] = 0.;
  fIsoConeSumPt [0] = 0.;
  fIsoConePtLead[1] = 0.;
  fIsoConeSumPt [1] = 0.;
}

//______________________________________________________________________________
///
/// Constructor.
///
/// \param p: TLorentzVector of particle kinematics.
///
/// particle: cluster or track
///
AliAODPWG4Particle::AliAODPWG4Particle(TLorentzVector & p):
  AliVParticle(),
  fMomentum(0),fPdg(-1), fTag(0), fLabel(-1),
  fCaloLabel(), fTrackLabel(),fDetectorTag(-1),
  fBadDist(0), fNLM(0), fM02(0), fM20(0),
  fTime(0),fNCells(0),fSuperModule(0),
  fDecayTag(0),fIsolated(0), fLeadingParticle(0),
  fIsoConePtLead(), fIsoConeSumPt(),
  fDisp(0), fTof(0), fCharged(0),
  fTagged(0), fFidArea(0), fInputFileIndex(0),fBtag(0)
{
  fMomentum = new TLorentzVector(p);
  
  fCaloLabel [0] = -1;
  fCaloLabel [1] = -1;
  fTrackLabel[0] = -1;
  fTrackLabel[1] = -1;
  fTrackLabel[2] = -1;
  fTrackLabel[3] = -1;
  
  fIsoConePtLead[0] = 0.;
  fIsoConeSumPt [0] = 0.;
  fIsoConePtLead[1] = 0.;
  fIsoConeSumPt [1] = 0.;  
}

//______________________________________________________________________________
///
/// Destructor.
///
AliAODPWG4Particle::~AliAODPWG4Particle() 
{
    delete fMomentum;
}

//______________________________________________________________________________
///
/// Clear pointers.
///
void AliAODPWG4Particle::Clear(const Option_t* /*opt*/) 
{
  delete fMomentum;
}

//______________________________________________________________________________
///
///  Copy constructor.
///
AliAODPWG4Particle::AliAODPWG4Particle(const AliAODPWG4Particle& part) :
  AliVParticle(part),
  fMomentum(0), fPdg(part.fPdg), fTag(part.fTag), fLabel(part.fLabel),
  fCaloLabel(), fTrackLabel(), fDetectorTag(part.fDetectorTag),
  fBadDist(part.fBadDist),fNLM(part.fNLM), fM02(part.fM02), fM20(part.fM20),
  fTime(part.fTime),fNCells(part.fNCells),fSuperModule(part.fSuperModule),
  fDecayTag(part.fDecayTag),fIsolated(part.fIsolated), fLeadingParticle(part.fLeadingParticle),
  fDisp(part.fDisp), fTof(part.fTof), fCharged(part.fCharged),
  fTagged(part.fTagged), fFidArea(part.fFidArea), fInputFileIndex(part.fInputFileIndex),fBtag(part.fBtag)
{
  fMomentum = new TLorentzVector(*part.fMomentum);
  
  fCaloLabel [0] = part.fCaloLabel[0];
  fCaloLabel [1] = part.fCaloLabel[1];
  fTrackLabel[0] = part.fTrackLabel[0];
  fTrackLabel[1] = part.fTrackLabel[1];
  fTrackLabel[2] = part.fTrackLabel[2];
  fTrackLabel[3] = part.fTrackLabel[3];
  
  fIsoConePtLead[0] = part.fIsoConePtLead[0];
  fIsoConeSumPt [0] = part.fIsoConeSumPt [0];
  fIsoConePtLead[1] = part.fIsoConePtLead[1];
  fIsoConeSumPt [1] = part.fIsoConeSumPt [1];
}

//________________________________________________________________________________
///
/// Assignment operator.
///
AliAODPWG4Particle& AliAODPWG4Particle::operator=(const AliAODPWG4Particle & part)
{
  if(this!=&part)
  {
    fPdg   = part.fPdg;
    fTag   = part.fTag;
    fLabel = part.fLabel;
    
    fCaloLabel [0] = part.fCaloLabel[0];
    fCaloLabel [1] = part.fCaloLabel[1];
    fTrackLabel[0] = part.fTrackLabel[0];
    fTrackLabel[1] = part.fTrackLabel[1];
    
    fIsoConePtLead[0] = part.fIsoConePtLead[0];
    fIsoConeSumPt [0] = part.fIsoConeSumPt [0];
    fIsoConePtLead[1] = part.fIsoConePtLead[1];
    fIsoConeSumPt [1] = part.fIsoConeSumPt [1];
    
    fDetectorTag = part.fDetectorTag;
    fDisp     = part.fDisp;
    fTof      = part.fTof;
    fCharged  = part.fCharged;
    fBadDist  = part.fBadDist;
    fDecayTag = part.fDecayTag;
    
    fNLM      = part.fNLM;
    fM02      = part.fM02;
    fM20      = part.fM20;
    fIsolated = part.fIsolated;
    fLeadingParticle =part.fLeadingParticle;

    fBtag     = part.fBtag;
    fFidArea  = part.fFidArea;
    fTagged   = part.fTagged;
    fInputFileIndex =  part.fInputFileIndex;

    if (fMomentum ) delete fMomentum;
    fMomentum = new TLorentzVector(*part.fMomentum);
  }
  
  return *this;
}

//_______________________________________________________________
///
/// \return true if particle satisfies given PID criterium
///
Bool_t AliAODPWG4Particle::IsPIDOK(Int_t ipid, Int_t pdgwanted) const
{
	switch(ipid)
  {
    case 0: return kTRUE ; //No PID at all
    case 1:
	  {
	    if (fPdg == pdgwanted) return kTRUE;
	    else return kFALSE; //Overall PID calculated with bayesian methods.
	  }
    case 2: return fDisp ;   //only dispersion cut
    case 3: return fTof ;    //Only TOF cut
    case 4: return fCharged ;    //Only Charged cut
    case 5: return fDisp && fTof ;  //Dispersion and TOF
    case 6: return fDisp && fCharged ;  //Dispersion and Charged
    case 7: return fTof  && fCharged ;  //TOF and Charged
    case 8: return fDisp && fTof && fCharged ; // all 3 cuts
    default: return kFALSE ; //Not known combination
	}
}

//_________________________________________________________
///
/// Print information of all data members.
///
void AliAODPWG4Particle::Print(Option_t* /*option*/) const 
{  
  printf("Particle 4-vector:\n");
  printf("     E  = %13.3f", E() );
  printf("     Px = %13.3f", Px());
  printf("     Py = %13.3f", Py());
  printf("     Pz = %13.3f\n", Pz());
  printf("Id PDG     : %d\n",fPdg);
  printf("MC Tag     : %d\n",fTag);
  printf("Dist. to bad channel : %d\n",fBadDist);

  printf("Detector  : %d, Labels:\n",fDetectorTag);
  printf("      Calo: %d, %d \n",fCaloLabel[0],fCaloLabel[1]);
  printf("      Track: %d, %d, %d, %d \n",fTrackLabel[0],fTrackLabel[1],fTrackLabel[2],fTrackLabel[3]);
  
  if(fDetectorTag!=2) // Avoid tracks, AliFiducialCut::kCTS
  {
    printf("Calo param: \n");
    printf("      M02: %2.2f\n",fM02);
    printf("      M20: %2.2f\n",fM20);
    printf("      NCell: %d\n",fNCells);
    printf("      Time: %2.3f\n",fTime);
    printf("      SModule: %d\n",fSuperModule);
  }
  
  printf("Tags: \n");
//  printf("Btag      : %d\n",fBtag);
  printf("     Pi0 Tag   : %d\n",fDecayTag);
  if(fIsolated)        printf("      Isolated! \n");
  if(fLeadingParticle) printf("      Leading! \n");
  
  printf("Isolation cone: \n");
  printf("\t charged: pT Max %2.2f, Sum pT %2.2f",fIsoConePtLead[0],fIsoConeSumPt[0]);
  printf("\t neutral: pT Max %2.2f, Sum pT %2.2f",fIsoConePtLead[1],fIsoConeSumPt[1]);
  
  printf("PID bits :\n");
  printf("     TOF        : %d",fTof);
  printf("     Charged    : %d",fCharged);
  printf("     Dispersion : %d\n",fDisp);

  //  printf("Fid Area  : %d\n",fFidArea);
  //  printf("Input File Index : %d\n",fInputFileIndex);
}
