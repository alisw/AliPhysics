//
// Author: Francesca Bellini (fbellini@cern.ch)
//
//

#include <Riostream.h>

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliRsnCutSetDaughterParticle.h"

class AliRsnCutPIDNSigma;
class AliRsnPIDRange;

ClassImp(AliRsnCutSetDaughterParticle)

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle() :
  AliRsnCutSet("AliRsnCutSetDaughterParticle", AliRsnTarget::kDaughter),
  fPID(AliPID::kPion),
  fAppliedCutSetID(AliRsnCutSetDaughterParticle::kNDaughterCuts),
  fNsigmaTPC(1E20),
  fNsigmaTOF(1E20)
{
   //
   // Default constructor
}

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const char *name, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID, AliPID::EParticleType pid, Float_t nSigmaFast = -1.0) :
  AliRsnCutSet(name, AliRsnTarget::kDaughter),
  fPID(pid),
  fAppliedCutSetID(cutSetID),
  fNsigmaTPC(1E20),
  fNsigmaTOF(1E20)
{
   //
   // Constructor
   //
  if ( (nSigmaFast<=0) && 
       ((cutSetID == AliRsnCutSetDaughterParticle::kFastTPCpidNsigma) || (cutSetID == AliRsnCutSetDaughterParticle::kFastTOFpidNsigma)) ) {
    AliError("Requested fast n-sigma PID with invalid value for n. Setting n = 1E20");
  } else {
    if (cutSetID == AliRsnCutSetDaughterParticle::kFastTPCpidNsigma){
      fNsigmaTPC = nSigmaFast;
    }
    if ( (cutSetID == AliRsnCutSetDaughterParticle::kFastTOFpidNsigma) ||
	(cutSetID == AliRsnCutSetDaughterParticle::kTOFpidKstarPbPb2010) ) { 
      fNsigmaTOF = nSigmaFast;
    }
  }
  
  Init();
}

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const AliRsnCutSetDaughterParticle &copy) :
   AliRsnCutSet(copy),
   fPID(copy.fPID),
   fAppliedCutSetID(copy.fAppliedCutSetID),
   fNsigmaTPC(copy.fNsigmaTPC),
   fNsigmaTOF(copy.fNsigmaTOF)
{
   //
   // copy constructor
}

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle &AliRsnCutSetDaughterParticle::operator=(const AliRsnCutSetDaughterParticle &copy)
 {
   //
   // operator =
  //
  AliRsnCutSet::operator=(copy);
  if (this == &copy)
    return *this;
  fPID=copy.fPID;
  fAppliedCutSetID=copy.fAppliedCutSetID;
  fNsigmaTPC=copy.fNsigmaTPC;
  fNsigmaTOF=copy.fNsigmaTOF;
  return (*this);
}

//----------------------------------------------------------------------------
void AliRsnCutSetDaughterParticle::Init() 
{ 
  //
  // init cut sets by setting variable params
  //  
      AliRsnCutTrackQuality * iCutQualityStd2010 = new AliRsnCutTrackQuality("CutQualityStd2010");
      iCutQualityStd2010->SetDefaults2010();
      /*
	requires ITS refit ,TPC refit, TPC in
	iCutQualityStd2010->SetPtRange(0.15, 1E+20);
	iCutQualityStd2010->SetEtaRange(-0.8, 0.8);
	iCutQualityStd2010->SetDCARPtFormula("0.0182+0.0350/pt^1.01");
	iCutQualityStd2010->SetDCAZmax(2.0);
	iCutQualityStd2010->SetSPDminNClusters(1);
	iCutQualityStd2010->SetITSminNClusters(0);
	iCutQualityStd2010->SetITSmaxChi2(36);
	iCutQualityStd2010->SetTPCminNClusters(70);
	iCutQualityStd2010->SetTPCmaxChi2(4.0);
	iCutQualityStd2010->SetRejectKinkDaughters();
	iCutQualityStd2010->SetAODTestFilterBit(5); 
      */	

      AliRsnCutTOFMatch  *iCutTOFMatch     = new AliRsnCutTOFMatch("CutTOFMatch");
      AliRsnCutPIDNSigma *iCutTPCNSigma    = new AliRsnCutPIDNSigma("CutTPCNSigma", fPID, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
      AliRsnCutPIDNSigma *iCutTPCTOFNSigma = new AliRsnCutPIDNSigma("CutTPCTOFNSigma", fPID, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
      AliRsnCutPIDNSigma *iCutTOFNSigma    = new AliRsnCutPIDNSigma("CutTOFNSigma", fPID, AliRsnCutPIDNSigma::kTOF);//, AliRsnCutPIDNSigma::kP );

      switch (fAppliedCutSetID) 
	{
	case AliRsnCutSetDaughterParticle::kNoCuts :
	  AliInfo("No cuts applied to daughter particle");
	  break;
	  
	case AliRsnCutSetDaughterParticle::kQualityStd2010 :	
	  AddCut(iCutQualityStd2010);
	  SetCutScheme(iCutQualityStd2010->GetName());
	  break;

	case AliRsnCutSetDaughterParticle::kQualityStd2011:	 
	  //these are also golden cuts for PbPb 2010 
	  //AliESDtrackCuts::GetStandardITSTPCTrackCuts2011()
	  // select golden cut for ESD -> to be implemented 
	  //select golden cut for AOD
	  iCutQualityStd2010->SetAODTestFilterBit(10); 	  //1024
	  AddCut(iCutQualityStd2010);
	  SetCutScheme(iCutQualityStd2010->GetName());
	  break;
	  
	case AliRsnCutSetDaughterParticle::kTOFMatch :
	  AddCut(iCutQualityStd2010);
	  AddCut(iCutTOFMatch);
	  SetCutScheme( Form("%s&(%s)",iCutQualityStd2010->GetName(), iCutTOFMatch->GetName()) );
	  break;
	  
	case    AliRsnCutSetDaughterParticle::kFastTPCpidNsigma :
	  if (fNsigmaTPC <= 0.0) {
	    AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTPCNSigma->GetName()));
	    SetNsigmaForFastTPCpid(10.0);
	  }
	  iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
	  AddCut(iCutQualityStd2010);
	  AddCut(iCutTPCNSigma);
	  SetCutScheme( Form("%s&%s",iCutQualityStd2010->GetName(), iCutTPCNSigma->GetName()) );
	  break;
      
	case    AliRsnCutSetDaughterParticle::kFastTOFpidNsigma :
	  if (fNsigmaTOF <= 0.0) {
	    AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
	    SetNsigmaForFastTOFpid(10.0);
	  }
	  iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
	  //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
	  AddCut(iCutQualityStd2010);
	  AddCut(iCutTOFNSigma);
	  SetCutScheme( Form("%s&%s",iCutQualityStd2010->GetName(), iCutTOFNSigma->GetName()) );
	  break;
	  
	case    AliRsnCutSetDaughterParticle::kTPCTOFpidKstarPP2010 :	  
	  /* Set TPC  PID (if no TOF)*/
	  // all   below   350         MeV: 5sigma
	  // all   between 350 and 500 MeV: 3sigma
	  // pios above   500         MeV: 2sigma
	  // kaons between 500 and 700 MeV: 2sigma
	  // kaons above   700         MeV: rejected
	  iCutTPCNSigma->AddPIDRange(5.0, 0.0, 0.35); 
	  iCutTPCNSigma->AddPIDRange(3.0, 0.35, 0.5); 
	  if (fPID==AliPID::kPion)
	    iCutTPCNSigma->AddPIDRange(2.0, 0.5, 1E20);
	  if (fPID==AliPID::kKaon)
	    iCutTPCNSigma->AddPIDRange(2.0, 0.5, 0.7);

	  AddCut(iCutQualityStd2010);
	  AddCut(iCutTOFMatch);
	  AddCut(iCutTPCNSigma);
	  
	  /* set TPC+TOF PID*/
	  iCutTPCTOFNSigma->SinglePIDRange(5.0);
	  iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1.5);
	  iCutTOFNSigma->AddPIDRange(2.0, 1.5, 1E20);

	  AddCut(iCutTPCTOFNSigma);
      	  AddCut(iCutTOFNSigma);
      
	  // scheme:
	  // quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]      
	  SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",iCutQualityStd2010->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
	  break;
	  
	case AliRsnCutSetDaughterParticle::kTOFpidKstarPbPb2010:
	  if (fNsigmaTOF <= 0.0) {
	    AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
	    SetNsigmaForFastTOFpid(10.0);
	  }
	  iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
	  //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
	  iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal
	  /* requires ITS refit ,TPC refit, TPC in*/
	  iCutQualityStd2010->SetPtRange(0.15, 1E+20);
	  iCutQualityStd2010->SetEtaRange(-0.8, 0.8);
	  iCutQualityStd2010->SetDCARPtFormula("0.0182+0.0350/pt^1.01");
	  iCutQualityStd2010->SetDCAZmax(2.0);
	  iCutQualityStd2010->SetSPDminNClusters(1);
	  iCutQualityStd2010->SetITSminNClusters(0);
	  iCutQualityStd2010->SetITSmaxChi2(36);
	  iCutQualityStd2010->SetTPCminNClusters(70);
	  iCutQualityStd2010->SetTPCmaxChi2(4.0);
	  iCutQualityStd2010->SetRejectKinkDaughters();
	  //iCutQualityStd2010->SetAODTestFilterBit(10); //AOD086
	  iCutQualityStd2010->SetAODTestFilterBit(5); //AOD049 no golden cuts  
	  AddCut(iCutQualityStd2010);
	  AddCut(iCutTOFNSigma);
	  AddCut(iCutTPCTOFNSigma);
	  SetCutScheme( Form("%s&%s&%s",iCutQualityStd2010->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName()) );
	  break;

	default :  
	  break;
	}
      
}



