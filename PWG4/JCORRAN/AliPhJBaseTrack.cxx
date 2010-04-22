#include "AliPhJBaseTrack.h"



// $Id: AliPhJBaseTrack.cxx,v 1.5 2008/05/08 15:19:52 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliPhJBaseTrack.cxx
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.5 $
  \date $Date: 2008/05/08 15:19:52 $
*/
////////////////////////////////////////////////////


//______________________________________________________________________________
AliPhJBaseTrack::AliPhJBaseTrack():
	fBasePt(-999), 
	fBaseTheta(-999), 
	fBaseEta(-999), 
	fBasePhi(-999), 
	fBaseE(-999),
        fBaseID(-999),
	fBaseCharge(-999),
        fBaseFlavor(kNone), 
        fBaseTof(-999),
        fBasePtot(-999)
{
  // constructor
}

//_____________________________________________________________
AliPhJBaseTrack::AliPhJBaseTrack(float pt,float theta, float eta, float phi, float e, short charge, float tof, float ptot):
	fBasePt(pt), 
	fBaseTheta(theta), 
	fBaseEta(eta), 
	fBasePhi(phi), 
	fBaseE(e),
        fBaseID(-999),
	fBaseCharge(charge),
        fBaseFlavor(kNone), 
        fBaseTof(tof),
        fBasePtot(ptot)
{
  // constructor
}

//_____________________________________________________________
AliPhJBaseTrack::AliPhJBaseTrack(const AliPhJBaseTrack& a):
  TObject(a),
  fBasePt(a.fBasePt),
  fBaseTheta(a.fBaseTheta),
  fBaseEta(a.fBaseEta),
  fBasePhi(a.fBasePhi),
  fBaseE(a.fBaseE),
  fBaseID(a.fBaseID),
  fBaseCharge(a.fBaseCharge),
  fBaseFlavor(a.fBaseFlavor),
  fBaseTof(a.fBaseTof),
  fBasePtot(a.fBasePtot)
{
  //copy constructor
}

//_____________________________________________________________

void AliPhJBaseTrack::PrintOut(const char *message=" " ) const{
    //object print out
    std::cout<<message
    <<" type="<<fBaseFlavor<<" "
    <<" pT="<<fBasePt
    <<" Theta="<<fBaseTheta
    <<" Eta="<<fBaseEta
    <<" Phi="<<fBasePhi
    <<" E="<<fBaseE
    <<" charge="<<fBaseCharge
    <<" tof="<<fBaseTof
    <<" ptot="<<fBasePtot<<std::endl;
}
//_____________________________________________________________
AliPhJBaseTrack& AliPhJBaseTrack::operator=(const AliPhJBaseTrack& trk){
  //operator =  
  if(this != &trk){
    TObject::operator=(trk);
     
    fBasePt     = trk.fBasePt;
    fBaseTheta  = trk.fBaseTheta;
    fBaseEta    = trk.fBaseEta;
    fBasePhi    = trk.fBasePhi;
    fBaseE      = trk.fBaseE;
    fBaseID     = trk.fBaseID;
    fBaseCharge = trk.fBaseCharge;
    fBaseFlavor = trk.fBaseFlavor;
    fBaseTof    = trk.fBaseTof;
    fBasePtot   = trk.fBasePtot;
  }
  return *this;
}

//_____________________________________________________________
void AliPhJBaseTrack::PrintJetInput(char *message=" " ) const{
    std::cout<<message
    <<"\t"<< GetPx() <<"\t"<< GetPy() <<"\t"<< GetPz()<<"\t"<<GetE() << std::endl;
}


ClassImp(AliPhJBaseTrack)

