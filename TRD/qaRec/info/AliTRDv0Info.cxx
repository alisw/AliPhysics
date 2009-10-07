#include "TMath.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliTRDv0Info.h"
#include "AliTRDtrackInfo.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliTRDtrackInfo.h"
#include "AliLog.h"

//Gathers all information necessary for reference data selection about the
//track and (in case) its corresponding V0.
//Carries out the selection of electrons (from gamma conversions), pions
//(from K0s decays) and protons (from Lambda and Anti-Lambda decays) by
//cuts specific for the respective decay and particle species.
//(M.Heide, 2009/10/06)

// Authors:
//   Alex Bercuci <A.Bercuci@gsi.de>
//   Alex Wilk    <wilka@uni-muenster.de>
//   Markus Heide <mheide@uni-muenster.de>
// 


ClassImp(AliTRDv0Info)

//_________________________________________________
AliTRDv0Info::AliTRDv0Info()
  : TObject()
  ,fESD(0x0)
  ,fHasV0(0)      
  ,fQuality(0)
    
  ,fMomentum(0)
  ,fDCA(10)
  ,fPointingAngle(10)
  ,fOpenAngle(10)
  ,fPsiPair(99)
  ,fMagField(0)
  ,fRadius(0)
    
  ,fTrackID(0)
  ,fV0Momentum(0)


  ,fTrackP(0x0)
  ,fTrackN(0x0)
  ,fTrack(0x0)
  ,fNindex(0)
  ,fPindex(0)
  
{
  memset(fPplus, 0, 2*kNlayer*sizeof(Float_t));
  memset(fPminus, 0, 2*kNlayer*sizeof(Float_t));
  memset(fDetPID, 0, 2*kNDaughters*kNDetectors*AliPID::kSPECIES*sizeof(Float_t));
  memset(fInvMass, 0, kNMomBins*kNDecays*sizeof(Double_t));

  /////////////////////////////////////////////////////////////////////////////
  //Set Cut values: First specify decay in brackets, then the actual cut value!
  ///////////////////////////////////////////////////////////////////////////// 

  //Upper limit for distance of closest approach of two daughter tracks :
  fUpDCA[kGamma] = 0.25;
  fUpDCA[kK0s] = 0.25;
  fUpDCA[kLambda] = 0.25;
  fUpDCA[kAntiLambda] = 0.25;

  //Upper limit for pointing angle (= angle between between vector from primary to secondary vertex and reconstructed momentum of V0 mother particle) :
  fUpPointingAngle[kGamma] = 0.03;
  fUpPointingAngle[kK0s] = 0.03;
  fUpPointingAngle[kLambda] = 0.03;
  fUpPointingAngle[kAntiLambda] = 0.03;

  //Upper limit for invariant mass of V0 mother :
  fUpInvMass[kGamma][0] = 0.04;// second pair of brackets is for momentum bin: 0: below mother momentm of 2.5 GeV
  fUpInvMass[kGamma][1] = 0.07;//1: above 2.5 GeV
  fUpInvMass[kK0s][0] = fUpInvMass[kK0s][1] = 0.51;
  fUpInvMass[kLambda][0] = fUpInvMass[kLambda][1] = 1.22;
  fUpInvMass[kAntiLambda][0] = fUpInvMass[kAntiLambda][1] = 1.22;

  //Lower limit for invariant mass of V0 mother :
  fDownInvMass[kGamma] = -1.;
  fDownInvMass[kK0s] = 0.49;
  fDownInvMass[kLambda] = 1.;
  fDownInvMass[kAntiLambda] = 1.;

  //Lower limit for distance from secondary vertex to primary vertex in x-y plane :
  fDownRadius[kGamma] = 5.7;
  fDownRadius[kK0s] = 0.;
  fDownRadius[kLambda] = 10.;
  fDownRadius[kAntiLambda] = 10.;

  //Upper limit for distance from secondary vertex to primary vertex in x-y plane :
  fUpRadius[kGamma] = 1000.;
  fUpRadius[kK0s] = 1000.;
  fUpRadius[kLambda] = 1000.;
  fUpRadius[kAntiLambda] = 1000.;

  //Upper limit for opening angle between two daughter tracks (characteristically near zero for conversions) :
  fUpOpenAngle[kGamma] = 0.1;
  fUpOpenAngle[kK0s] = 3.15;
  fUpOpenAngle[kLambda] = 3.15;
  fUpOpenAngle[kAntiLambda] = 3.15;

  //Upper limit for angle between daughter momentum plane and plane perpendicular to magnetic field (characteristically around zero for conversions) :
  fUpPsiPair[kGamma] = 0.1;
  fUpPsiPair[kK0s] = 1.6;
  fUpPsiPair[kLambda] = 1.6;
  fUpPsiPair[kAntiLambda] = 1.6;

  //Lower limit for likelihood value of TPC PID :
  fDownTPCPIDneg[AliPID::kElectron] = 0.21;
  fDownTPCPIDpos[AliPID::kElectron] = 0.21;

  fDownTPCPIDneg[AliPID::kMuon] = 0.21;
  fDownTPCPIDpos[AliPID::kMuon] = 0.21;

  fDownTPCPIDneg[AliPID::kPion] = 0.21;
  fDownTPCPIDpos[AliPID::kPion] = 0.21;

  fDownTPCPIDneg[AliPID::kKaon] = 0.21;
  fDownTPCPIDpos[AliPID::kKaon] = 0.21;

  fDownTPCPIDneg[AliPID::kProton] = 0.21;
  fDownTPCPIDpos[AliPID::kProton] = 0.21;
  //////////////////////////////////////////////////////////////////////////////////

}

//_________________________________________________
void AliTRDv0Info::GetESDv0Info(AliESDv0 *esdv0)
{//Gets values of ESDv0 and daughter track properties
  //See header file for description of variables

  Int_t part1 = -1;
  Int_t part2 = -1;

  fQuality = Quality(esdv0);//Attributes an Int_t to the V0 due to quality cuts (= 1 if V0 is accepted, other integers depending on cut which excludes the vertex)    

  fRadius = Radius(esdv0);//distance from secondary vertex to primary vertex in x-y plane
      
  fDCA = esdv0->GetDcaV0Daughters();//distance of closest approach of two daughter tracks
      
  fPointingAngle = TMath::ACos(esdv0->GetV0CosineOfPointingAngle());// pointing angle (= angle between between vector from primary to secondary vertex and reconstructed momentum of V0 mother particle)
      
  fOpenAngle = OpenAngle(esdv0);//Opening angle between two daughter tracks
      
  fPsiPair = PsiPair(esdv0);//Angle between daughter momentum plane and plane perpendicular to magnetic field

  fV0Momentum = V0Momentum(esdv0);//Reconstructed momentum of the mother particle
      
  for(Int_t idecay = 0; idecay < kNDecays; idecay++)//4 decay types : conversions, K0s, Lambda, Anti-Lambda 
    //five particle types: electrons, muons, pions, kaons, protons (muons and kaons not involved)
    {
      if(idecay == kLambda)//protons and pions from Lambda
  {
    part1 = AliPID::kProton;
    part2 = AliPID::kPion;
  }
      else if(idecay == kAntiLambda)//antiprotons and pions from Anti-Lambda
  {
    part1 = AliPID::kPion;
    part2 = AliPID::kProton;
  }
      else if(idecay == kK0s)//pions from K0s
  part1 = part2 = AliPID::kPion;
      else if(idecay == kGamma)//electrons from conversions
  part1 = part2 = AliPID::kElectron;
    
      fInvMass[idecay] = InvMass(part1, part2, esdv0);//Calculate invariant mass for all of our four supposed decays
    }
  GetDetectorPID();//Gets all likelihood values from TPC, TOF and ITS PID for the fDetPID[kNDaughters][kNDetectors][AliPID::kSPECIES] array

    
}
//_________________________________________________
Float_t  AliTRDv0Info::V0Momentum(AliESDv0 *esdv0)
{//Reconstructed momentum of V0 mother particle
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};


  esdv0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter; 
  esdv0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;
  
  
  return TMath::Sqrt((mn[0]+mp[0])*(mn[0]+mp[0]) + (mn[1]+mp[1])*(mn[1]+mp[1])+(mn[2]+mp[2])*(mn[2]+mp[2]));
}

//_________________________________________________
Double_t AliTRDv0Info::InvMass(Int_t part1, Int_t part2, AliESDv0 *esdv0)
{//Invariant mass of reconstructed V0 mother

  const Double_t kpmass[5] = {AliPID::ParticleMass(AliPID::kElectron),AliPID::ParticleMass(AliPID::kMuon),AliPID::ParticleMass(AliPID::kPion),AliPID::ParticleMass(AliPID::kKaon),AliPID::ParticleMass(AliPID::kProton)};
  //Masses of electrons, muons, pions, kaons and protons, as implemented in ROOT


  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};  

  esdv0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  esdv0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;
  
  Double_t mass1 = kpmass[part1];//sets supposed rest masses for both daughters
  Double_t mass2 = kpmass[part2];   

  //Calculate daughters' energies :
  Double_t e1    = TMath::Sqrt(mass1*mass1+
            mp[0]*mp[0]+
            mp[1]*mp[1]+
            mp[2]*mp[2]);
  Double_t e2    = TMath::Sqrt(mass2*mass2+
            mn[0]*mn[0]+
            mn[1]*mn[1]+
            mn[2]*mn[2]);  

  //Sum of daughter momenta :   
  Double_t momsum =  
    (mn[0]+mp[0])*(mn[0]+mp[0])+
    (mn[1]+mp[1])*(mn[1]+mp[1])+
    (mn[2]+mp[2])*(mn[2]+mp[2]);

  //invariant mass :	  	     
  Double_t InvMass = TMath::Sqrt((e1+e2)*(e1+e2)-momsum);

  return InvMass;
  
}
//_________________________________________________
Float_t AliTRDv0Info::OpenAngle(AliESDv0 *esdv0)
{//Opening angle between two daughter tracks
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
    

  esdv0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  esdv0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter;

  
  fOpenAngle = TMath::ACos((mp[0]*mn[0] + mp[1]*mn[1] + mp[2]*mn[2])/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1] + mp[2]*mp[2])*TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1] + mn[2]*mn[2])));
  
  return fOpenAngle;
}

//_________________________________________________
Float_t AliTRDv0Info::PsiPair(AliESDv0 *esdv0)
{//Angle between daughter momentum plane and plane perpendicular to magnetic field
  Double_t x, y, z;
  esdv0->GetXYZ(x,y,z);//Reconstructed coordinates of V0; to be replaced by Markus Rammler's method in case of conversions!
  
  Double_t mn[3] = {0,0,0};
  Double_t mp[3] = {0,0,0};
  

  esdv0->GetNPxPyPz(mn[0],mn[1],mn[2]);//reconstructed cartesian momentum components of negative daughter;
  esdv0->GetPPxPyPz(mp[0],mp[1],mp[2]);//reconstructed cartesian momentum components of positive daughter; 


  Double_t deltat = 1.;
  deltat = TMath::ATan(mp[2]/(TMath::Sqrt(mp[0]*mp[0] + mp[1]*mp[1])+1.e-13)) -  TMath::ATan(mn[2]/(TMath::Sqrt(mn[0]*mn[0] + mn[1]*mn[1])+1.e-13));//difference of angles of the two daughter tracks with z-axis

  Double_t radiussum = TMath::Sqrt(x*x + y*y) + 50;//radius to which tracks shall be propagated

  Double_t MomPosProp[3];
  Double_t MomNegProp[3];
    
  AliExternalTrackParam nt(*fTrackN), pt(*fTrackP);
    
  fPsiPair = 4.;

  if(nt.PropagateTo(radiussum,fMagField) == 0)//propagate tracks to the outside
    fPsiPair =  -5.;
  if(pt.PropagateTo(radiussum,fMagField) == 0)
    fPsiPair = -5.;
  pt.GetPxPyPz(MomPosProp);//Get momentum vectors of tracks after propagation
  nt.GetPxPyPz(MomNegProp);
  
  Double_t p_ele =
    TMath::Sqrt(MomNegProp[0]*MomNegProp[0]+MomNegProp[1]*MomNegProp[1]+MomNegProp[2]*MomNegProp[2]);//absolute momentum value of negative daughter
  Double_t p_pos =
    TMath::Sqrt(MomPosProp[0]*MomPosProp[0]+MomPosProp[1]*MomPosProp[1]+MomPosProp[2]*MomPosProp[2]);//absolute momentum value of positive daughter
    
  Double_t scalarproduct =
    MomPosProp[0]*MomNegProp[0]+MomPosProp[1]*MomNegProp[1]+MomPosProp[2]*MomNegProp[2];//scalar product of propagated positive and negative daughters' momenta
    
  Double_t chipair = TMath::ACos(scalarproduct/(p_ele*p_pos));//Angle between propagated daughter tracks

  fPsiPair =  TMath::Abs(TMath::ASin(deltat/chipair));  

  return fPsiPair; 

}
//_________________________________________________
void AliTRDv0Info::V0fromTrack(AliTRDtrackInfo *track, Int_t ivertex)
{//Checks if track is a secondary vertex daughter (due to V0 finder)
  
  fMagField = fESD->GetMagneticField();

  fTrackID =  track->GetTrackId();//index of the track

  fTrack = fESD->GetTrack(fTrackID);//sets track information

  fHasV0 = 0;

  //comparing index of track with indices of pos./neg. V0 daughter :
  AliESDv0 * esdv0 = fESD->GetV0(ivertex);
  if((esdv0->GetIndex(0) == fTrackID)||(esdv0->GetIndex(1) == fTrackID))
    {
      fHasV0 = 1;//track belongs to vertex found by V0 finder!
      fNindex = esdv0->GetIndex(0);
      fPindex = esdv0->GetIndex(1);
      fTrackN = fESD->GetTrack(esdv0->GetIndex(0));//providing information about the other of the two daughter tracks 
      fTrackP = fESD->GetTrack(esdv0->GetIndex(1));
      GetESDv0Info(esdv0);//gets all the relevant information about our V0
    }
}
//_________________________________________________
void AliTRDv0Info::GetDetectorPID()
{//PID likelihoods from TPC, TOF, and ITS, for all particle species

  fTrackN->GetTPCpid(fDetPID[kNeg][kTPC]);
  fTrackP->GetTPCpid(fDetPID[kPos][kTPC]);
  fTrackN->GetTOFpid(fDetPID[kNeg][kTOF]);
  fTrackP->GetTOFpid(fDetPID[kPos][kTOF]);
  fTrackN->GetITSpid(fDetPID[kNeg][kITS]);
  fTrackP->GetITSpid(fDetPID[kPos][kITS]);

}

//_________________________________________________
Float_t AliTRDv0Info::Radius(AliESDv0 *esdv0)
{//distance from secondary vertex to primary vertex in x-y plane
  Double_t x, y, z;
  esdv0->GetXYZ(x,y,z); //Reconstructed coordinates of V0; to be replaced by Markus Rammler's method in case of conversions!
  fRadius = TMath::Sqrt(x*x + y*y);
  return fRadius;

}
//_________________________________________________
Int_t AliTRDv0Info::Quality(AliESDv0 *esdv0)
{ //Checking track and V0 quality status in order to exclude vertices based on poor information
  Float_t NclsN;
  NclsN = fTrackN->GetTPCNcls();//number of found clusters in TPC for negative track
  Float_t NclsFN;
  NclsFN = fTrackN->GetTPCNclsF();//number of findable clusters in TPC for negative track
  Float_t NclsP;
  NclsP = fTrackP->GetTPCNcls();//number of found clusters in TPC for positive track
  Float_t NclsFP;
  NclsFP = fTrackP->GetTPCNclsF();//number of findable clusters in TPC for positive track
  
  fQuality = 0;


  Float_t ClsRatioN; 
  Float_t ClsRatioP;

  if((NclsFN == 0) || (NclsFP == 0))
    return 2;
    
  ClsRatioN = NclsN/NclsFN; //ratios of found to findable clusters in TPC 
  ClsRatioP = NclsP/NclsFP;
    
  if (!(esdv0->GetOnFlyStatus()))//accept only vertices from online V0 finder
    return 3;
  if (!((fTrackP->GetStatus() &
  AliESDtrack::kTPCrefit)))//accept only vertices in which both tracks have TPC refit
    return 4;
  if (!((fTrackN->GetStatus() &
  AliESDtrack::kTPCrefit)))
    return 5;	
  if (fTrackP->GetKinkIndex(0)>0  ||
      fTrackN->GetKinkIndex(0)>0 )//exclude tracks with kinks
    return 6;
  if((ClsRatioN < 0.6)||(ClsRatioP < 0.6))//exclude tracks with low ratio of found to findable TPC clusters
    return 7;
  fQuality = 1;
  return fQuality;
}
//_________________________________________________
Bool_t AliTRDv0Info::GetV0PID(Int_t ipart, AliTRDtrackInfo *track)
{//decides if track is accepted for one of the reference data samples
  
  Int_t iDecay = -1;

  //decide which decay has to be considered for which particle sample (Anti-Lambda will be treated separately)
  if(ipart == AliPID::kElectron)
    iDecay = kGamma;
  else if(ipart == AliPID::kPion)
    iDecay = kK0s;
  else if(ipart == AliPID::kProton)
    iDecay = kLambda;

  Int_t iPSlot;//Mother momentum slots above/below 2.5 GeV

  
  Bool_t pid = 0;//return value for V0 PID decision

  if(!(track))
    {
      AliError("AliTRDv0Info::GetV0PID(Int_t ipart, AliTRDtrackInfo *track) : No track info found.\n");
      return 0;
    }

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!esdH)
    {
      AliError("AliTRDv0Info::GetV0PID(Int_t ipart, AliTRDtrackInfo *track) : ERROR - ESD input handler not found");
      return 0;
    } 
      
  
  fESD = esdH->GetEvent();
  
  for(Int_t ivertex=0; ivertex<fESD->GetNumberOfV0s(); ivertex++)
    {
    
      if(pid == 0)
  {     
    V0fromTrack(track, ivertex);//Get the V0 corresponding to the track (if there is a V0)
    
    if(fV0Momentum > 2.5)//divide into slots according to reconstructed momentum of the mother particle
      {iPSlot = 1;}
    else
      {iPSlot = 0;}
    //Accept track for a sample only if...

    if(!(fHasV0))//... there is a V0 found for it
      continue;
    if(!(fQuality == 1))//... it fulfills our quality criteria
      continue;
    if(!(fDCA < fUpDCA[iDecay]))//... distance of closest approach between daughters is reasonably small
      continue;
    if(!(fPointingAngle < fUpPointingAngle[iDecay]))//... pointing angle between momentum of mother particle and vector from prim. to sec. vertex is small
      continue;				  
    if(!(fRadius > fDownRadius[iDecay]))//... x-y plane distance of decay point to prim. vertex is bigger than a certain minimum value (for conversions)
      continue;
    if(!(fOpenAngle < fUpOpenAngle[iDecay]))//... opening angle is close enough to zero (for conversions)
      continue;
    if(!(TMath::Abs(fPsiPair) < fUpPsiPair[iDecay]))//... Psi-pair angle is close enough to zero(for conversions)
      continue;

    //specific cut criteria :
    if(ipart == AliPID::kProton)
      {//for proton sample: separate treatment of Lamba and Anti-Lambda decays:
        //for Anti-Lambda:
        //TPC PID likelihoods high enough for pi+ and anti-proton ; invariant mass calculated postulating these two particle species...
        if((fDetPID[kNeg][kTPC][AliPID::kProton] > fDownTPCPIDneg[AliPID::kProton]) && (fDetPID[kPos][kTPC][AliPID::kPion] > fDownTPCPIDpos[AliPID::kPion]))
    {
      if(fNindex == fTrackID)
        {
          if((fInvMass[kAntiLambda] < fUpInvMass[kAntiLambda][iPSlot]) && (fInvMass[kAntiLambda] > fDownInvMass[kAntiLambda]))
      {
        pid = 1;
      }
        }
    }
        //for Lambda:
        //TPC PID likelihoods high enough for pi- and proton ; invariant mass calculated accordingly
        if((fDetPID[kNeg][kTPC][AliPID::kPion] > fDownTPCPIDneg[AliPID::kPion]) && (fDetPID[kPos][kTPC][AliPID::kProton] > fDownTPCPIDpos[AliPID::kProton]))
    {
      if(fPindex == fTrackID)
        {
          if((fInvMass[kLambda] < fUpInvMass[kLambda][iPSlot]) && (fInvMass[kLambda] > fDownInvMass[kLambda]))
      {
        pid = 1;
      }
        }
    }
      }
    //for photon and K0s decays: equal TPC PID likelihood criteria for both daughters ; invariant mass calculated postulating two electrons/two pions
    else 		  
      if((fDetPID[kNeg][kTPC][ipart] > fDownTPCPIDneg[ipart]) && (fDetPID[kPos][kTPC][ipart] > fDownTPCPIDpos[ipart]))
        {
    if((fInvMass[iDecay] < fUpInvMass[iDecay][iPSlot]) && (fInvMass[iDecay] > fDownInvMass[iDecay]))
      {
        pid = 1;						  
      }
        }
  }
    }

  return pid;
  
}
//_________________________________________________
void AliTRDv0Info::Print(Option_t */*opt*/) const
{

}
