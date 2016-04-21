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


//-------------------------------------------------------------------------
//     AOD special track implementation of AliVTrack
//     Author: Michele Floris, CERN
//     michele.floris@cern.ch
//-------------------------------------------------------------------------

#include <TVector3.h>
#include "AliLog.h"
#include "AliExternalTrackParam.h"
#include "AliVVertex.h"
#include "AliDetectorPID.h"
#include "AliAODEvent.h"
#include "AliAODHMPIDrings.h"

#include "AliNanoAODTrack.h"
#include "AliNanoAODTrackMapping.h"

ClassImp(AliNanoAODTrack)


//______________________________________________________________________________
AliNanoAODTrack::AliNanoAODTrack() : 
  AliVTrack(),
  AliNanoAODStorage(),
  fLabel(0),
  fProdVertex(0),
  fCharge(0),
  fAODEvent(NULL)
{
  // default constructor
  // The default constructor should not allocate memory! You risk an infinite loop here.
  //  AllocateInternalStorage(AliNanoAODTrackMapping::GetInstance()->GetSize());

  //  AllocateInternalStorage("pt,theta,phi");
  // FIXME: TO BE REIMPLEMENTED
  // SetPosition((Float_t*)NULL);
  // SetXYAtDCA(-999., -999.);
  // SetPxPyPzAtDCA(-999., -999., -999.);
  // SetPID((Float_t*)NULL);
  // for (Int_t i = 0; i < 3; i++) {fTOFLabel[i] = -1;}
}

//______________________________________________________________________________
AliNanoAODTrack::AliNanoAODTrack(AliAODTrack * aodTrack, const char * vars) :
  AliVTrack(), 
  AliNanoAODStorage(),
  fLabel(0),
  fProdVertex(0),
  fCharge(0),
  fAODEvent(NULL)
{
  // constructor

  Double_t position[3];
  Bool_t isPosAvailable = aodTrack->GetPosition(position);
  AliNanoAODTrackMapping::GetInstance(vars);

  // Create internal structure
  AllocateInternalStorage(AliNanoAODTrackMapping::GetInstance()->GetSize());

  for (Int_t index = 0; index<AliNanoAODTrackMapping::GetInstance()->GetSize(); index++) {
    TString varString = AliNanoAODTrackMapping::GetInstance()->GetVarName(index);

    if     (varString == "pt"                     ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPt()               , aodTrack->Pt()                      );
    else if(varString == "phi"                    ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPhi()              , aodTrack->Phi()                     );
    else if(varString == "theta"                  ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTheta()            , aodTrack->Theta()                   );
    else if(varString == "chi2perNDF"             ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetChi2PerNDF()       , aodTrack->Chi2perNDF()              );  
    else if(varString == "posx" && isPosAvailable ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX()             , position[0]                         );
    else if(varString == "posy" && isPosAvailable ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY()             , position[1]                         );
    else if(varString == "posz" && isPosAvailable ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPosZ()             , position[2]                         );
    else if(varString == "posDCAx"                ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPosDCAx()          , aodTrack->XAtDCA()                  );
    else if(varString == "posDCAy"                ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPosDCAy()          , aodTrack->YAtDCA()                  );
    else if(varString == "pDCAx"                  ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPDCAX()            , aodTrack->PxAtDCA()                 );
    else if(varString == "pDCAy"                  ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPDCAY()            , aodTrack->PyAtDCA()                 );
    else if(varString == "pDCAz"                  ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetPDCAZ()            , aodTrack->PzAtDCA()                 );
    else if(varString == "RAtAbsorberEnd"         ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetRAtAbsorberEnd()   , aodTrack->GetRAtAbsorberEnd()       );
    else if(varString == "TPCncls"                ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCncls()          , aodTrack->GetTPCNcls()              );
    else if(varString == "id"                     ) SetVar(AliNanoAODTrackMapping::GetInstance()->Getid()               , aodTrack->GetID()                   );
    else if(varString == "TPCnclsF"               ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCnclsF()         , aodTrack->GetTPCNclsF()             );
    else if(varString == "TPCNCrossedRows"        ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCNCrossedRows()  , aodTrack->GetTPCNCrossedRows()      );
    else if(varString == "TrackPhiOnEMCal"        ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTrackPhiOnEMCal()  , aodTrack->GetTrackPhiOnEMCal()      );
    else if(varString == "TrackEtaOnEMCal"        ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTrackEtaOnEMCal()  , aodTrack->GetTrackEtaOnEMCal()      );
    else if(varString == "TrackPtOnEMCal"         ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTrackPtOnEMCal()   , aodTrack->GetTrackPtOnEMCal()       );
    else if(varString == "ITSsignal"              ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetITSsignal()        , aodTrack->GetITSsignal()            );
    else if(varString == "TPCsignal"              ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCsignal()        , aodTrack->GetTPCsignal()            );
    else if(varString == "TPCsignalTuned"         ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCsignalTuned()   , aodTrack->GetTPCsignalTunedOnData() );
    else if(varString == "TPCsignalN"             ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCsignalN()       , aodTrack->GetTPCsignalN()           );
    else if(varString == "TPCmomentum"            ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCmomentum()      , aodTrack->GetTPCmomentum()          );
    else if(varString == "TPCTgl"                 ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTPCTgl()           , aodTrack->GetTPCTgl()               );
    else if(varString == "TOFsignal"              ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTOFsignal()        , aodTrack->GetTOFsignal()            );
    else if(varString == "integratedLength"       ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetintegratedLenght() , aodTrack->GetIntegratedLength()     );
    else if(varString == "TOFsignalTuned"         ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTOFsignalTuned()   , aodTrack->GetTOFsignalTunedOnData() );
    else if(varString == "HMPIDsignal"            ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetHMPIDsignal()      , aodTrack->GetHMPIDsignal()          );
    else if(varString == "HMPIDoccupancy"         ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetHMPIDoccupancy()   , aodTrack->GetHMPIDoccupancy()       );
    else if(varString == "TRDsignal"              ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTRDsignal()        , aodTrack->GetTRDsignal()            );
    else if(varString == "TRDChi2"                ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTRDChi2()          , aodTrack->GetTRDchi2()              );
    else if(varString == "TRDnSlices"             ) SetVar(AliNanoAODTrackMapping::GetInstance()->GetTRDnSlices()       , aodTrack->GetNumberOfTRDslices()    );  
    else if(varString == "covmat"                 ) AliFatal("cov matrix To be implemented"                            );
  }


  fLabel = aodTrack->GetLabel();
  fCharge = aodTrack->Charge();
  fProdVertex = aodTrack->GetProdVertex();
  // SetUsedForVtxFit(usedForVtxFit);// FIXME: what is this
  // SetUsedForPrimVtxFit(usedForPrimVtxFit);// FIXME: what is this
  // //  if(covMatrix) SetCovMatrix(covMatrix);// FIXME: 
  // for (Int_t i=0;i<3;i++) {fTOFLabel[i]=-1;}

}

//______________________________________________________________________________
AliNanoAODTrack::AliNanoAODTrack(AliESDTrack * /*esdTrack*/, const char * /*vars*/) : 
  AliVTrack(), 
  AliNanoAODStorage(),
  fLabel(0),
  fProdVertex(0),
  fCharge(0),
  fAODEvent(NULL)
{
  // ctor: Creates a special track by copying the requested variables from an ESD track
  AliFatal("To be Implemented");
}


AliNanoAODTrack::AliNanoAODTrack(const char * vars) :
  AliVTrack(),
  AliNanoAODStorage(),
  fLabel(0),
  fProdVertex(0),
  fCharge(0),
  fAODEvent(NULL)
{
   // ctor: Creates a special track simply allocating the required variables
  AliNanoAODTrackMapping::GetInstance(vars);

  // Create internal structure
  AllocateInternalStorage(AliNanoAODTrackMapping::GetInstance()->GetSize());


}

//______________________________________________________________________________
AliNanoAODTrack::~AliNanoAODTrack() 
{
  //  std::cout << "1 " << this << " "  << fKinVars << " " << fNKinVars << std::endl;
  //  this->Print();
  // destructor
  // if(fKinVars)		 {
  //   std::cout << "2" << std::endl;
  //   delete [] fKinVars;
  //   std::cout << "3" << std::endl;
  //   fKinVars = 0;
  //   std::cout << "4" << std::endl;
  // }

  // if(fCovMatrix) {
  //   delete fCovMatrix;
  //   fCovMatrix = 0;
  // }
}


//______________________________________________________________________________
AliNanoAODTrack::AliNanoAODTrack(const AliNanoAODTrack& trk) :
  AliVTrack(),
  AliNanoAODStorage(),
  fLabel(trk.fLabel),
  fProdVertex(trk.fProdVertex),
  fCharge(trk.fCharge),
  fAODEvent(trk.fAODEvent)
{
  // Copy constructor
  // std::cout << "Copy Ctor" << std::endl;
  
  AllocateInternalStorage(AliNanoAODTrackMapping::GetInstance()->GetSize());
  for (Int_t isize = 0; isize<AliNanoAODTrackMapping::GetInstance()->GetSize(); isize++) {
    SetVar(isize, trk.GetVar(isize));    
  }


}

//______________________________________________________________________________
AliNanoAODTrack& AliNanoAODTrack::operator=(const AliNanoAODTrack& trk)
{
  // Assignment operator
  if(this!=&trk) {

    AliVTrack::operator=(trk); // FIXME: I think I should overload this...
    AliNanoAODStorage::operator=(trk);

    fLabel      = trk.fLabel;
    fProdVertex = trk.fProdVertex;
    fCharge     = trk.fCharge;
    fAODEvent   = trk.fAODEvent;
    
  }

  return *this;
}



//______________________________________________________________________________
Double_t AliNanoAODTrack::M(AliAODTrack::AODTrkPID_t pid) const
{
  // Returns the mass.
  // Masses for nuclei don't exist in the PDG tables, therefore they were put by hand.

  switch (pid) {

  case AliAODTrack::kElectron :
    return 0.000510999; //TDatabasePDG::Instance()->GetParticle(11/*::kElectron*/)->Mass();
    break;

  case AliAODTrack::kMuon :
    return 0.1056584; //TDatabasePDG::Instance()->GetParticle(13/*::kMuonMinus*/)->Mass();
    break;

  case AliAODTrack::kPion :
    return 0.13957; //TDatabasePDG::Instance()->GetParticle(211/*::kPiPlus*/)->Mass();
    break;

  case AliAODTrack::kKaon :
    return 0.4937; //TDatabasePDG::Instance()->GetParticle(321/*::kKPlus*/)->Mass();
    break;

  case AliAODTrack::kProton :
    return 0.9382720; //TDatabasePDG::Instance()->GetParticle(2212/*::kProton*/)->Mass();
    break;

  case AliAODTrack::kDeuteron :
    return 1.8756; //TDatabasePDG::Instance()->GetParticle(1000010020)->Mass();
    break;

  case AliAODTrack::kTriton :
    return 2.8089; //TDatabasePDG::Instance()->GetParticle(1000010030)->Mass();
    break;

  case AliAODTrack::kHelium3 :
    return 2.8084; //TDatabasePDG::Instance()->GetParticle(1000020030)->Mass();
    break;

  case AliAODTrack::kAlpha :
    return 3.7274; //TDatabasePDG::Instance()->GetParticle(1000020040)->Mass();
    break;

  case AliAODTrack::kUnknown :
    return -999.;
    break;

  default :
    return -999.;
  }
}

//______________________________________________________________________________
Double_t AliNanoAODTrack::E(AliAODTrack::AODTrkPID_t pid) const
{
  // Returns the energy of the particle of a given pid.
  
  if (pid != AliAODTrack::kUnknown) { // particle was identified
    Double_t m = M(pid);
    return TMath::Sqrt(P()*P() + m*m);
  } else { // pid unknown
    return -999.;
  }
}

//______________________________________________________________________________
Double_t AliNanoAODTrack::Y(AliAODTrack::AODTrkPID_t pid) const
{
  // Returns the rapidity of a particle of a given pid.
  
  if (pid != AliAODTrack::kUnknown) { // particle was identified
    Double_t e = E(pid);
    Double_t pz = Pz();
    if (e>=0 && e!=pz) { // energy was positive (e.g. not -999.) and not equal to pz
      return 0.5*TMath::Log((e+pz)/(e-pz));
    } else { // energy not known or equal to pz
      return -999.;
    }
  } else { // pid unknown
    return -999.;
  }
}

//______________________________________________________________________________
Double_t AliNanoAODTrack::Y(Double_t m) const
{
  // Returns the rapidity of a particle of a given mass.
  
  if (m >= 0.) { // mass makes sense
    Double_t e = E(m);
    Double_t pz = Pz();
    if (e>=0 && e!=pz) { // energy was positive (e.g. not -999.) and not equal to pz
      return 0.5*TMath::Log((e+pz)/(e-pz));
    } else { // energy not known or equal to pz
      return -999.;
    }
  } else { // pid unknown
    return -999.;
  }
}



//______________________________________________________________________________
template <typename T> void AliNanoAODTrack::SetP(const T *p, const Bool_t cartesian) 
{
  // Set the momentum

  if (p) {
    if (cartesian) {
      // This is inherited from AliAODtrack... I don't think we want/need this in the special track
      AliFatal("Not implemented");
    } else {
      SetVar(AliNanoAODTrackMapping::GetInstance()->GetPt()      , p[0]);  
      SetVar(AliNanoAODTrackMapping::GetInstance()->GetPhi()     , p[1]);  
      SetVar(AliNanoAODTrackMapping::GetInstance()->GetTheta()   , p[2]);  
    }
  } else {
      SetVar(AliNanoAODTrackMapping::GetInstance()->GetPt()      , p[0]);  
      SetVar(AliNanoAODTrackMapping::GetInstance()->GetPhi()     , p[1]);  
      SetVar(AliNanoAODTrackMapping::GetInstance()->GetTheta()   , p[2]);  
  }
}

/*
//______________________________________________________________________________
template <typename T> void AliNanoAODTrack::SetPosition(const T *x, const Bool_t dca) 
{
  // set the position

  if (x) {
    if (!dca) {
      ResetBit(kIsDCA);

      fPosition[0] = x[0];
      fPosition[1] = x[1];
      fPosition[2] = x[2];
    } else {
      SetBit(kIsDCA);
      // don't know any better yet
      fPosition[0] = -999.;
      fPosition[1] = -999.;
      fPosition[2] = -999.;
    }
  } else {
    ResetBit(kIsDCA);

    fPosition[0] = -999.;
    fPosition[1] = -999.;
    fPosition[2] = -999.;
  }
}
*/
//______________________________________________________________________________
void AliNanoAODTrack::SetDCA(Double_t d, Double_t z) 
{
  // set the dca 

  // FIXME: this is a hack which was taken over from the AliAODtrack,
  // where the same variable is used to store DCA or position,
  // according to the value of the bit kIsDCA. We can probably get rid
  // of this in the special track.
  SetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX(), d);
  SetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY(), z);
  SetVar(AliNanoAODTrackMapping::GetInstance()->GetPosZ(), 0);
  SetBit(AliAODTrack::kIsDCA);
}

//______________________________________________________________________________
void AliNanoAODTrack::Print(Option_t* /* option */) const
{
  // prints information about AliNanoAODTrack
  //  std::cout << "Size: " << AliNanoAODTrackMapping::GetInstance()->GetSize() << std::endl;
  AliNanoAODTrackMapping::GetInstance()->Print();

  for (Int_t index = 0; index<AliNanoAODTrackMapping::GetInstance()->GetSize(); index++) {
    printf(" - [%2.2d] %-10s : %f\n", index, AliNanoAODTrackMapping::GetInstance()->GetVarName(index), GetVar(index));    
  }
  std::cout << "" << std::endl;  

}



//______________________________________________________________________________
Bool_t AliNanoAODTrack::PropagateToDCA(const AliVVertex *vtx, 
    Double_t b, Double_t maxd, Double_t dz[2], Double_t covar[3])
{
  // compute impact parameters to the vertex vtx and their covariance matrix
  // b is the Bz, needed to propagate correctly the track to vertex 
  // only the track parameters are update after the propagation (pos and mom),
  // not the covariance matrix. This is OK for propagation over short distance
  // inside the beam pipe.
  // return kFALSE is something went wrong

  // allowed only for tracks inside the beam pipe
  Float_t xstart2 = GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX())*GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX())+GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY())*GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY());
  
  if(xstart2 > 3.*3.) { // outside beampipe radius
    AliError("This method can be used only for propagation inside the beam pipe");
    return kFALSE; 
  }

  // convert to AliExternalTrackParam
  AliExternalTrackParam etp; etp.CopyFromVTrack(this);  

  // propagate
  if(!etp.PropagateToDCA(vtx,b,maxd,dz,covar)) return kFALSE;

  // update track position and momentum
  Double_t mom[3];
  etp.GetPxPyPz(mom);
  SetP(mom,kTRUE);
  etp.GetXYZ(mom);
  SetPosition(mom,kFALSE);


  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliNanoAODTrack::GetPxPyPz(Double_t p[3]) const 
{
    //---------------------------------------------------------------------
    // This function returns the global track momentum components
    //---------------------------------------------------------------------
  p[0]=Px(); p[1]=Py(); p[2]=Pz();
  return kTRUE;
}



//_____________________________________________________________________________
//_____________________________________________________________________________

Bool_t AliNanoAODTrack::GetXYZAt(Double_t x, Double_t b, Double_t *r) const
{
  //---------------------------------------------------------------------
  // This function returns the global track position extrapolated to
  // the radial position "x" (cm) in the magnetic field "b" (kG)
  //---------------------------------------------------------------------

  //conversion of track parameter representation is
  //based on the implementation of AliExternalTrackParam::Set(...)
  //maybe some of this code can be moved to AliVTrack to avoid code duplication
  const double kSafe = 1e-5;
  Double_t alpha=0.0;
  Double_t radPos2 = GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX())*GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX())+GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY())*GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY());
  Double_t radMax  = 45.; // approximately ITS outer radius
  if (radPos2 < radMax*radMax) { // inside the ITS     
    alpha = TMath::ATan2(Py(),Px());
  } else { // outside the ITS
    Float_t phiPos = TMath::Pi()+TMath::ATan2(-GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY()), -GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX()));
     alpha = 
     TMath::DegToRad()*(20*((((Int_t)(phiPos*TMath::RadToDeg()))/20))+10);
  }
  //
  Double_t cs=TMath::Cos(alpha), sn=TMath::Sin(alpha);
  // protection:  avoid alpha being too close to 0 or +-pi/2
  if (TMath::Abs(sn)<kSafe) {
    alpha = kSafe;
    cs=TMath::Cos(alpha);
    sn=TMath::Sin(alpha);
  }
  else if (cs<kSafe) {
    alpha -= TMath::Sign(kSafe, alpha);
    cs=TMath::Cos(alpha);
    sn=TMath::Sin(alpha);    
  }
  
  // Get the vertex of origin and the momentum
  TVector3 ver(GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosX()), GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosY()), GetVar(AliNanoAODTrackMapping::GetInstance()->GetPosZ()));
  TVector3 mom(Px(),Py(),Pz());
  //
  // avoid momenta along axis
  if (TMath::Abs(mom[0])<kSafe) mom[0] = TMath::Sign(kSafe*TMath::Abs(mom[1]), mom[0]);
  if (TMath::Abs(mom[1])<kSafe) mom[1] = TMath::Sign(kSafe*TMath::Abs(mom[0]), mom[1]);

  // Rotate to the local coordinate system
  ver.RotateZ(-alpha);
  mom.RotateZ(-alpha);

  Double_t param0 = ver.Y();
  Double_t param1 = ver.Z();
  Double_t param2 = TMath::Sin(mom.Phi());
  Double_t param3 = mom.Pz()/mom.Pt();
  Double_t param4 = TMath::Sign(1/mom.Pt(),(Double_t)fCharge);

  //calculate the propagated coordinates
  //this is based on AliExternalTrackParam::GetXYZAt(Double_t x, Double_t b, Double_t *r)
  Double_t dx=x-ver.X();
  if(TMath::Abs(dx)<=kAlmost0) return GetXYZ(r);

  Double_t f1=param2;
  Double_t f2=f1 + dx*param4*b*kB2C;

  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;
  
  Double_t r1=TMath::Sqrt((1.-f1)*(1.+f1)), r2=TMath::Sqrt((1.-f2)*(1.+f2));
  r[0] = x;
  r[1] = param0 + dx*(f1+f2)/(r1+r2);
  r[2] = param1 + dx*(r2 + f2*(f1+f2)/(r1+r2))*param3;//Thanks to Andrea & Peter

  return Local2GlobalPosition(r,alpha);
}


//_______________________________________________________

void  AliNanoAODTrack::Clear(Option_t * /*opt*/) {
  // empty storage
  fVars.clear();
  fNVars = 0;
}
