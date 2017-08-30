/**************************************************************************
 * Copyright(c) 2007-2011, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////
//                                                               //
// Class to compute variables for correction framework           //  
// for 3-body decays of D mesons (D+, Ds, Lc)                    //
// in bins of cut variables                                      //
// Origin:       Francesco Prino (prino@to.infn.it)              //
//               Renu Bala       (bala@to.infn.it)               //
//               Davide Caffarri (cafarri@pd.infn.it)            //
///////////////////////////////////////////////////////////////////

#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "TClonesArray.h"
#include "AliCFVertexingHF.h"
#include "AliESDtrack.h"
#include "TDatabasePDG.h"

#include "AliCFVertexingHF3Prong.h"
#include "AliCFContainer.h"
#include "AliCFTaskVertexingHF.h"

/// \cond CLASSIMP
ClassImp(AliCFVertexingHF3Prong);
/// \endcond

//_________________________________________
AliCFVertexingHF3Prong::AliCFVertexingHF3Prong(Int_t decay, UInt_t resonantDecay):
AliCFVertexingHF(),
  fDecay(decay),
  fGenDsOption(kCountResonant),
  fResonantDecay(resonantDecay)
 {
  // 
  SetNProngs(3);

  fPtAccCut=new Float_t[fProngs];
  fEtaAccCut=new Float_t[fProngs];
  for(Int_t iP=0; iP<fProngs; iP++){
	  fPtAccCut[iP]=0.1;
	  fEtaAccCut[iP]=0.9;
  }

}
//_________________________________________
AliCFVertexingHF3Prong::AliCFVertexingHF3Prong(Int_t decay):
AliCFVertexingHF(),
  fDecay(decay),
  fGenDsOption(kCountResonant),
  fResonantDecay(0)
 {
  //
  SetNProngs(3);

  fPtAccCut=new Float_t[fProngs];
  fEtaAccCut=new Float_t[fProngs];
  for(Int_t iP=0; iP<fProngs; iP++){
          fPtAccCut[iP]=0.1;
          fEtaAccCut[iP]=0.9;
  }

}
//_________________________________________
AliCFVertexingHF3Prong::AliCFVertexingHF3Prong(TClonesArray *mcArray, UShort_t originDselection, Int_t decay, UInt_t resonantDecay):
  AliCFVertexingHF(mcArray, originDselection),
  fDecay(decay),
  fGenDsOption(kCountResonant),
  fResonantDecay(resonantDecay)
{
  //
  SetNProngs(3);
  fPtAccCut=new Float_t[fProngs];
  fEtaAccCut=new Float_t[fProngs];
  for(Int_t iP=0; iP<fProngs; iP++){
	  fPtAccCut[iP]=0.1;
	  fEtaAccCut[iP]=0.9;
  }
}

//_________________________________________
AliCFVertexingHF3Prong::AliCFVertexingHF3Prong(TClonesArray *mcArray, UShort_t originDselection, Int_t decay):
  AliCFVertexingHF(mcArray, originDselection),
  fDecay(decay),
  fGenDsOption(kCountResonant),
  fResonantDecay(0)
{
  //
  SetNProngs(3);
  fPtAccCut=new Float_t[fProngs];
  fEtaAccCut=new Float_t[fProngs];
  for(Int_t iP=0; iP<fProngs; iP++){
          fPtAccCut[iP]=0.1;
          fEtaAccCut[iP]=0.9;
  }
}
//_____________________________________
AliCFVertexingHF3Prong& AliCFVertexingHF3Prong::operator=(const AliCFVertexingHF3Prong& c){
  //
  if  (this != &c) {

    AliCFVertexingHF::operator=(c);
   
  }
    return *this;
}

//__________________________________________
Bool_t AliCFVertexingHF3Prong::SetRecoCandidateParam(AliAODRecoDecayHF *recoCand){
  // Checks if candidate is signal and D meson is present in MC array
  
  Bool_t bSignAssoc = kFALSE;
  fRecoCandidate = recoCand;

  if (!fRecoCandidate) {
    AliError("fRecoCandidate not found, problem in assignement\n");
    return bSignAssoc;
  }
  
  Int_t pdgCand = -1;
  Int_t pdgDaughter[3]={-1,-1,-1};
  if(fDecay==kDplustoKpipi){
    pdgCand=411;
    pdgDaughter[0]=321;
    pdgDaughter[1]=211;
    pdgDaughter[2]=211;
  }else if(fDecay==kDstoKKpi){
    pdgCand=431;
    pdgDaughter[0]=321;
    pdgDaughter[1]=321;
    pdgDaughter[2]=211;
  }else if(fDecay==kLctopKpi){
	  pdgCand=4122;
	  pdgDaughter[0]=2212;
	  pdgDaughter[1]=321;
	  pdgDaughter[2]=211;	  
  }else{
    AliError("WRONG DECAY SETTING");
    return bSignAssoc;    
  }

  Int_t mcLabel = fRecoCandidate->MatchToMC(pdgCand,fmcArray,3,pdgDaughter);  
  if (mcLabel == -1) return bSignAssoc;

  if (fRecoCandidate->NumberOfFakeDaughters()>0){
	  fFake = 0;    // fake candidate
	  if (fFakeSelection==1) return bSignAssoc;
  }
  if (fRecoCandidate->NumberOfFakeDaughters()==0){
	  fFake = 2;    // non-fake candidate
	  if (fFakeSelection==2) return bSignAssoc;
  }
  
  SetMCLabel(mcLabel);
  fmcPartCandidate = dynamic_cast<AliAODMCParticle*>(fmcArray->At(fmcLabel));
    
  if (!fmcPartCandidate){
    AliDebug(3,"No part candidate");
    return bSignAssoc;
  }

  if(fDecay==kDstoKKpi && fGenDsOption!=kCountAllDsKKpi){
    if(!CheckMCChannelDecay()){
      AliDebug(3,"Ds not from the selected resonant channel");
      return bSignAssoc;
    }
  }

  if (fDecay==kLctopKpi && fResonantDecay != AliCFTaskVertexingHF::kAll) {
    if (!CheckLc3Prong()) return bSignAssoc;
  }

  bSignAssoc = kTRUE;
  return bSignAssoc;
}

//______________________________________________
Bool_t AliCFVertexingHF3Prong::GetGeneratedValuesFromMCParticle(Double_t* vectorMC) {
	// 
	// collecting all the necessary info from MC particle and fill vectorMC: 12 variables
	// pt_D
	// y_D
	// phi_D
	// ctau
	// cos point
	// pt_1
	// pt_2
	// pt_3
	// d0_1
	// d0_2
	// d0_3
	// zPrimVert
	// centrality
	
	Bool_t bGenValues = kFALSE;
	
	Int_t pdgCand = -1;
	if(fDecay==kDplustoKpipi){
		pdgCand=411;
	}else if(fDecay==kDstoKKpi){
		pdgCand=431;
	}else if(fDecay==kLctopKpi){
		pdgCand=4122;
	}else{
		AliError("WRONG DECAY SETTING");
		return bGenValues;
	}
	
	Double_t vertD[3] = {0,0,0};   // D origin
	fmcPartCandidate->XvYvZv(vertD);  // cm
	
	Int_t nprongs = 3;
	Int_t daughter[3];
	Short_t charge = fmcPartCandidate->Charge();
	
	// order the daughters as LS,OS,LS, e.g. D+ -> pi+ K- pi+
	// the 2 LS are ordered so that in pos. 0 there is the one with lower label value
	Int_t index=0;
	Int_t nDauLS=0;
	Int_t nDauOS=0;
	
	
	Int_t nDau=fmcPartCandidate->GetNDaughters();
	Int_t labelFirstDau = fmcPartCandidate->GetDaughter(0); 
	if(nDau==3){
		for(Int_t iDau=0; iDau<3; iDau++){
			Int_t ind = labelFirstDau+iDau;
			AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(fmcArray->At(ind));
			if(!part){
				AliError("Daughter particle not found in MC array");
				return bGenValues;
			}
			Short_t signDau=part->Charge();
			if(signDau==charge){
				nDauLS++;
				daughter[index] = ind;
				index=2;
			}else{
				daughter[1] = ind;
				nDauOS++;
			}
		}
	}else if(nDau==2){
		for(Int_t iDau=0; iDau<2; iDau++){
			Int_t ind = labelFirstDau+iDau;
			AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(fmcArray->At(ind));
			if(!part){
				AliError("Daughter particle not found in MC array");
				return bGenValues;
			}
			Int_t pdgCode=TMath::Abs(part->GetPdgCode());
			if(pdgCode==211 || pdgCode==321 || pdgCode==2212){
				Short_t signDau=part->Charge();
				if(signDau==charge){
					nDauLS++;
					daughter[index] = ind;
					index=2;
				}else{
					daughter[1] = ind;
					nDauOS++;
				}
			}else{
				Int_t nDauRes=part->GetNDaughters();
				if(nDauRes!=2){
					AliError("Wrong resonant decay");
					return bGenValues;
				}
				Int_t labelFirstDauRes = part->GetDaughter(0); 	
				for(Int_t iDauRes=0; iDauRes<2; iDauRes++){
					Int_t indDR = labelFirstDauRes+iDauRes;
					AliAODMCParticle* partDR = dynamic_cast<AliAODMCParticle*>(fmcArray->At(indDR));
					if(!partDR){
						AliError("Daughter particle not found in MC array");
						return bGenValues;
					}
					Short_t signDau=partDR->Charge();
					if(signDau==charge){
						nDauLS++;
						daughter[index] = ind;
						index=2;
					}else{
						daughter[1] = ind;
						nDauOS++;
					}
				}
			}
		}
	}else{
		AliError(Form("Wrong number of daughters %d",nDau));
		return bGenValues;
	}
	
	if(nDauLS!=2 || nDauOS!=1){
		AliError(Form("Wrong decay channel: LS and OS daughters not OK: %d %d",nDauLS,nDauOS));
		return bGenValues;
	}
	if(daughter[0]>daughter[2]){
		Int_t tmp=daughter[0];
		daughter[0]=daughter[2];
		daughter[2]=tmp;
	}
	
	// getting the momentum from the daughters and decay vertex
	Double_t px[3],py[3],pz[3],pt[3];
	Double_t vertDec[3] = {0,0,0};   // decay vertex		
	for(Int_t iDau=0; iDau<3; iDau++){
		AliAODMCParticle* part=dynamic_cast<AliAODMCParticle*>(fmcArray->At(daughter[iDau]));
		if(!part){
			AliError("Daughter particle not found in MC array");
			return bGenValues;
		}
		px[iDau]=part->Px();
		py[iDau]=part->Py();
		pz[iDau]=part->Pz();
		pt[iDau]=part->Pt();
		if(iDau==0) part->XvYvZv(vertDec);
	}
	
	Double_t d0[3] = {0.,0.,0.}; // dummy values!!!!
	
	AliAODRecoDecayHF* decay = new AliAODRecoDecayHF(vertD,vertDec,nprongs,charge,px,py,pz,d0);
	Double_t cT = decay->Ct(pdgCand);
	
	switch (fConfiguration){
	case AliCFTaskVertexingHF::kSnail:
		vectorMC[0] = fmcPartCandidate->Pt();
		vectorMC[1] = fmcPartCandidate->Y() ;
		vectorMC[2] = fmcPartCandidate->Phi();
		vectorMC[3] = cT*1.E4 ;  // in micron
		vectorMC[4] = 1.01;    // cos pointing angle, dummy value, meaningless in MC
		vectorMC[5] = pt[0];
		vectorMC[6] = pt[1];
		vectorMC[7] = pt[2];
		vectorMC[8] = fzMCVertex;    // z of reconstructed of primary vertex
		vectorMC[9] = fCentValue; // reconstructed centrality value 
		vectorMC[10] = 1.;           // fake: always filling with 1 at MC level 
		vectorMC[11] = 1.01; // dummy value for cosPointingXY  multiplicity
		vectorMC[12] = 0.; // dummy value for NormalizedDecayLengthXY multiplicity
		vectorMC[13] = fMultiplicity; // reconstructed multiplicity
		
		if (fDecay==kLctopKpi){
			vectorMC[11] = 0.; //dist12
			vectorMC[12] = 0.; //dist23
			vectorMC[13] = 0.; //sigmaVtx
			vectorMC[14] = 0.; //sumd02
			vectorMC[15] = 1.01; // dummy value for cosPointingXY  multiplicity
			vectorMC[16] = 0.; // dummy value for NormalizedDecayLengthXY multiplicity
			vectorMC[17] = fMultiplicity; // reconstructed multiplicity
		}
		break;
		
	case AliCFTaskVertexingHF::kCheetah:
		vectorMC[0] = fmcPartCandidate->Pt();
		vectorMC[1] = fmcPartCandidate->Y() ;
		vectorMC[2] = cT*1.E4; // in micron
		vectorMC[3] = fmcPartCandidate->Phi();
		vectorMC[4] = fzMCVertex;
		vectorMC[5] = fCentValue;   // dummy value for dca, meaningless in MC
		vectorMC[6] = 1. ;  // fake: always filling with 1 at MC level 
		vectorMC[7] = fMultiplicity;   // dummy value for d0pi, meaningless in MC, in micron
		break;
	}
	delete decay;	
	bGenValues = kTRUE;
	return bGenValues;
}


//____________________________________________
Bool_t AliCFVertexingHF3Prong::GetRecoValuesFromCandidate(Double_t *vectorReco) const
{ 
	// Fill vector (see above) with reconstructed quantities
	Bool_t bFillRecoValues=kFALSE;
	
	Int_t pdgCand = -1;
	if(fDecay==kDplustoKpipi){
		pdgCand=411;
	}else if(fDecay==kDstoKKpi){
		pdgCand=431;
	}else if(fDecay==kLctopKpi){
		pdgCand=4122;
		// AliError("LambdaC not yet implemented");
		// return bFillRecoValues;
	}else{
		AliError("WRONG DECAY SETTING");
		return bFillRecoValues;
	}
	
	AliAODRecoDecayHF3Prong *decay3 = (AliAODRecoDecayHF3Prong*)fRecoCandidate;
	Short_t charge=decay3->Charge();
	Double_t rapidity=decay3->Y(pdgCand);
	Double_t cT=decay3->Ct(pdgCand); 
	Double_t pt = decay3->Pt();
	Double_t cosPointingAngle = decay3->CosPointingAngle();
	Double_t phi = decay3->Phi();
	Double_t dist12= decay3->GetDist12toPrim();
	Double_t dist23 = decay3->GetDist23toPrim();
	Double_t sigmVert = decay3->GetSigmaVert();
	Double_t cosPointingAngleXY = decay3->CosPointingAngleXY();
	Double_t normDecayLengthXY = decay3->NormalizedDecayLengthXY();
		
	Int_t daughtSorted[3];
	Int_t tmpIndex=0;
	Int_t nDauLS=0;
	Int_t nDauOS=0;
	for(Int_t iDau=0; iDau<3; iDau++){
		AliAODTrack *trk = (AliAODTrack*)decay3->GetDaughter(iDau);
		Short_t chargedau=trk->Charge();
		if(chargedau==charge){
		  daughtSorted[tmpIndex]=iDau;
			tmpIndex=2;
			nDauLS++;
		}else{
		  daughtSorted[1]=iDau;
			nDauOS++;
		}
	}
	
	if(nDauLS!=2 || nDauOS!=1){
		AliError("Wrong decay channel: number of OS and LS tracks not OK");
		return bFillRecoValues;
	}
	
	if(daughtSorted[0]>daughtSorted[2]){
		Int_t tmp=daughtSorted[0];
		daughtSorted[0]=daughtSorted[2];
		daughtSorted[2]=tmp;
	}
	
	
	switch (fConfiguration){
	case AliCFTaskVertexingHF::kSnail:
		vectorReco[0] = pt;
		vectorReco[1] = rapidity;
		vectorReco[2] = phi;
		vectorReco[3] = cT*1.E4;  // in micron
		vectorReco[4] = cosPointingAngle;  // in micron
		vectorReco[5] = decay3->PtProng(daughtSorted[0]);
		vectorReco[6] = decay3->PtProng(daughtSorted[1]);
		vectorReco[7] = decay3->PtProng(daughtSorted[2]);
		vectorReco[8] = fzPrimVertex;    // z of reconstructed of primary vertex
		vectorReco[9] = fCentValue; //reconstructed centrality value
		vectorReco[10] = fFake;      // whether the reconstructed candidate was a fake (fFake = 0) or not (fFake = 2) 
		vectorReco[11] = cosPointingAngleXY; 
		vectorReco[12] = normDecayLengthXY; // in cm
		vectorReco[13] = fMultiplicity; // reconstructed multiplicity
		
		if(fDecay==kLctopKpi){  
		        Double_t d0prong0 = decay3->Getd0Prong(daughtSorted[0]);
		        Double_t d0prong1 = decay3->Getd0Prong(daughtSorted[1]);
		        Double_t d0prong2 = decay3->Getd0Prong(daughtSorted[2]);
			Double_t sumd02 =(d0prong0*d0prong0 + d0prong1*d0prong1 + d0prong2*d0prong2); 
			vectorReco[11] = dist12*1.E4;
			vectorReco[12] = dist23*1.E4;
			vectorReco[13] = sigmVert*1.E4;
			vectorReco[14] = sumd02*1.E8;
			vectorReco[15] = cosPointingAngleXY; 
			vectorReco[16] = normDecayLengthXY; // in cm
			vectorReco[17] = fMultiplicity; // reconstructed multiplicity
		}
		break;
	case AliCFTaskVertexingHF::kCheetah:
		vectorReco[0] = pt;
		vectorReco[1] = rapidity ;
		vectorReco[2] = cT*1.E4; // in micron
		vectorReco[3] = phi; 
		vectorReco[4] = fzPrimVertex;
		vectorReco[5] = fCentValue;   
		vectorReco[6] = fFake ; 
		vectorReco[7] = fMultiplicity;  
		break;
	}

	bFillRecoValues = kTRUE;
	return bFillRecoValues;
}


//_____________________________________________________________
Bool_t AliCFVertexingHF3Prong::CheckMCChannelDecay() const
{ 
  // Check the pdg codes of the daughters
  Bool_t checkCD = kFALSE;

  Int_t pdgCand = -1;
  Int_t pdgDaughter[3]={-1,-1,-1};
  if(fDecay==kDplustoKpipi){
    pdgCand=411;
    pdgDaughter[0]=321;
    pdgDaughter[1]=211;
    pdgDaughter[2]=211;
  }else if(fDecay==kDstoKKpi){
    pdgCand=431;
    pdgDaughter[0]=321;
    pdgDaughter[1]=321;
    pdgDaughter[2]=211;
  }else if(fDecay==kLctopKpi){
	  pdgCand=4122;
	  pdgDaughter[0]=2212;
	  pdgDaughter[1]=321;
	  pdgDaughter[2]=211;
	  
  //  AliError("LambdaC not yet implemented");
  //  return checkCD;
  }else{
    AliError("WRONG DECAY SETTING");
    return checkCD;    
  }


  Int_t daughter[3];
  Double_t sumPxDau=0.;
  Double_t sumPyDau=0.;
  Double_t sumPzDau=0.;

  Int_t nDau=fmcPartCandidate->GetNDaughters();
  Int_t labelFirstDau = fmcPartCandidate->GetDaughter(0);

  if(fDecay==kLctopKpi && fResonantDecay!=AliCFTaskVertexingHF::kAll){if(!CheckLc3Prong()) return checkCD;}
 
  if(nDau==3){
    if(fDecay==kDstoKKpi && !(fGenDsOption==kCountAllDsKKpi || fGenDsOption==kCountNonResonant)){
      AliDebug(3,"Decay channel in direct KKpi, should be skipped");
      return checkCD;
    }
    for(Int_t iDau=0; iDau<3; iDau++){
      Int_t ind = labelFirstDau+iDau;
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(fmcArray->At(ind));
      if(!part){
	AliError("Daughter particle not found in MC array");
	return checkCD;
      }
      daughter[iDau]=TMath::Abs(part->GetPdgCode());
      sumPxDau+=part->Px();
      sumPyDau+=part->Py();
      sumPzDau+=part->Pz();
    }
  }else if(nDau==2){
    if(fDecay==kDstoKKpi && fGenDsOption==kCountNonResonant) return checkCD;
    Int_t nDauFound=0;
    for(Int_t iDau=0; iDau<2; iDau++){
      Int_t ind = labelFirstDau+iDau;
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(fmcArray->At(ind));
      if(!part){
	AliError("Daughter particle not found in MC array");
	return checkCD;
      }
      Int_t pdgCode=TMath::Abs(part->GetPdgCode());
      if(pdgCode==211 || pdgCode==321 || pdgCode==2212){
	if(nDauFound>=3) return checkCD;
	daughter[nDauFound]=pdgCode;
	sumPxDau+=part->Px();
	sumPyDau+=part->Py();
	sumPzDau+=part->Pz();
	nDauFound++;
      }else{
	if(fDecay==kDstoKKpi && fGenDsOption!=3){
	  Int_t pdgCodeRes=TMath::Abs(part->GetPdgCode());
	  if(fGenDsOption==kCountPhipi && pdgCodeRes!=333) return checkCD;
	  else if(fGenDsOption==kCountK0stK && pdgCodeRes!=313) return checkCD;
	}
	Int_t nDauRes=part->GetNDaughters();
	if(nDauRes!=2) return checkCD;
	Int_t labelFirstDauRes = part->GetDaughter(0); 	
	for(Int_t iDauRes=0; iDauRes<2; iDauRes++){
	  Int_t indDR = labelFirstDauRes+iDauRes;
	  AliAODMCParticle* partDR = dynamic_cast<AliAODMCParticle*>(fmcArray->At(indDR));
	  if(!partDR){
	    AliError("Daughter particle not found in MC array");
	    return checkCD;
	  }
	  Int_t pdgCodeDR=TMath::Abs(partDR->GetPdgCode());
	  if(nDauFound>=3) return checkCD;
	  daughter[nDauFound]=pdgCodeDR;
	  sumPxDau+=partDR->Px();
	  sumPyDau+=partDR->Py();
	  sumPzDau+=partDR->Pz();
	  nDauFound++;
	}
      }
    }
  }else{
    return checkCD;
  }
  for(Int_t iDau1=0; iDau1<3; iDau1++){
    for(Int_t iDau2=iDau1; iDau2<3; iDau2++){
      if(daughter[iDau1]<daughter[iDau2]){
	Int_t tmp=daughter[iDau1];
	daughter[iDau1]=daughter[iDau2];
	daughter[iDau2]=tmp;
      }
    }
  }
  for(Int_t iDau=0; iDau<3; iDau++){
    if(daughter[iDau]!=pdgDaughter[iDau]){
      AliDebug(2, "Wrong decay channel from MC, skipping!!");
      return checkCD;  
    }
  }
  
  Double_t pxMother=fmcPartCandidate->Px();
  Double_t pyMother=fmcPartCandidate->Py();
  Double_t pzMother=fmcPartCandidate->Pz();
  if(TMath::Abs(pxMother-sumPxDau)/(TMath::Abs(pxMother)+1.e-13)>0.00001 ||
     TMath::Abs(pyMother-sumPyDau)/(TMath::Abs(pyMother)+1.e-13)>0.00001 ||
     TMath::Abs(pzMother-sumPzDau)/(TMath::Abs(pzMother)+1.e-13)>0.00001){
    AliDebug(2, "Momentum conservation violated, skipping!!");
    return checkCD;  
  }

  checkCD = kTRUE;
  return checkCD;
  
}
//---------------------------------
Bool_t AliCFVertexingHF3Prong::CheckLc3Prong() const
{
  Int_t numberOfLambdac=0;
  if(TMath::Abs(fmcPartCandidate->GetPdgCode())!=4122) return kFALSE;
  Int_t nDaugh = (Int_t)fmcPartCandidate->GetNDaughters();
  if(nDaugh<2) return kFALSE;
  if(nDaugh>3) return kFALSE;
  AliAODMCParticle* pdaugh1 = (AliAODMCParticle*)fmcArray->At(fmcPartCandidate->GetDaughter(0));
  if(!pdaugh1) {return kFALSE;}
  Int_t number1 = TMath::Abs(pdaugh1->GetPdgCode());
  AliAODMCParticle* pdaugh2 = (AliAODMCParticle*)fmcArray->At(fmcPartCandidate->GetDaughter(1));
  if(!pdaugh2) {return kFALSE;}
  Int_t number2 = TMath::Abs(pdaugh2->GetPdgCode());
  if(nDaugh==3){
    if(fResonantDecay!=AliCFTaskVertexingHF::kNonResonant && fResonantDecay!=AliCFTaskVertexingHF::kAll)return kFALSE; 
    Int_t thirdDaugh=fmcPartCandidate->GetDaughter(1)-1;
    AliAODMCParticle* pdaugh3 = (AliAODMCParticle*)fmcArray->At(thirdDaugh);
    if(!pdaugh3) return kFALSE;
    Int_t number3 = TMath::Abs(pdaugh3->GetPdgCode());
    if((number1==321 && number2==211 && number3==2212) || (number1==211 && number2==321 && number3==2212) || (number1==211 && number2==2212 && number3==321) || (number1==321 && number2==2212 && number3==211) || (number1==2212 && number2==321 && number3==211) || (number1==2212 && number2==211 && number3==321)) numberOfLambdac++;
  }
  if(nDaugh==2){
    if(fResonantDecay==AliCFTaskVertexingHF::kNonResonant)return kFALSE; 
    Int_t nfiglieK=0;
    if((number1==2212 && number2==313)){
      if(fResonantDecay!=AliCFTaskVertexingHF::kKstar && fResonantDecay!=AliCFTaskVertexingHF::kAll)return kFALSE; 
      nfiglieK=pdaugh2->GetNDaughters();
      if(nfiglieK!=2) return kFALSE;
      AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)fmcArray->At(pdaugh2->GetDaughter(0));
      AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)fmcArray->At(pdaugh2->GetDaughter(1));
      if(!pdaughK1) return kFALSE;
      if(!pdaughK2) return kFALSE;
      if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac++;
    }
    if((number1==313 && number2==2212)){
      if(fResonantDecay!=AliCFTaskVertexingHF::kKstar && fResonantDecay!=AliCFTaskVertexingHF::kAll)return kFALSE; 
      nfiglieK=pdaugh1->GetNDaughters();
      if(nfiglieK!=2) return kFALSE;
      AliAODMCParticle* pdaughK1 = (AliAODMCParticle*)fmcArray->At(pdaugh1->GetDaughter(0));
      AliAODMCParticle* pdaughK2 = (AliAODMCParticle*)fmcArray->At(pdaugh1->GetDaughter(1));
      if(!pdaughK1) return kFALSE;
      if(!pdaughK2) return kFALSE;
      if((TMath::Abs(pdaughK1->GetPdgCode())==211 && TMath::Abs(pdaughK2->GetPdgCode())==321) || (TMath::Abs(pdaughK1->GetPdgCode())==321 && TMath::Abs(pdaughK2->GetPdgCode())==211)) numberOfLambdac++;
    }
    Int_t nfiglieDelta=0;
    if(number1==321 && number2==2224){
      if(fResonantDecay!=AliCFTaskVertexingHF::kDelta && fResonantDecay!=AliCFTaskVertexingHF::kAll)return kFALSE; 
      nfiglieDelta=pdaugh2->GetNDaughters();
      if(nfiglieDelta!=2) return kFALSE;
      AliAODMCParticle *pdaughD1=(AliAODMCParticle*)fmcArray->At(pdaugh2->GetDaughter(0));
      AliAODMCParticle *pdaughD2=(AliAODMCParticle*)fmcArray->At(pdaugh2->GetDaughter(1));
      if(!pdaughD1) return kFALSE;
      if(!pdaughD2) return kFALSE;
      if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac++;
    }
    if(number1==2224 && number2==321){
      if(fResonantDecay!=AliCFTaskVertexingHF::kDelta && fResonantDecay!=AliCFTaskVertexingHF::kAll)return kFALSE; 
      nfiglieDelta=pdaugh1->GetNDaughters();
      if(nfiglieDelta!=2) return kFALSE;
      AliAODMCParticle* pdaughD1 = (AliAODMCParticle*)fmcArray->At(pdaugh1->GetDaughter(0));
      AliAODMCParticle* pdaughD2 = (AliAODMCParticle*)fmcArray->At(pdaugh1->GetDaughter(1)); 
      if(!pdaughD1) return kFALSE;
      if(!pdaughD2) return kFALSE;
      if((TMath::Abs(pdaughD1->GetPdgCode())==211 && TMath::Abs(pdaughD2->GetPdgCode())==2212) || (TMath::Abs(pdaughD1->GetPdgCode())==2212 && TMath::Abs(pdaughD2->GetPdgCode())==211)) numberOfLambdac++;
    }
    
    Int_t nfiglieLa=0;
    if(number1==3124 && number2==211){
      if(fResonantDecay!=AliCFTaskVertexingHF::kL1520 && fResonantDecay!=AliCFTaskVertexingHF::kAll)return kFALSE; 
      nfiglieLa=pdaugh1->GetNDaughters();
      if(nfiglieLa!=2) return kFALSE;
      AliAODMCParticle *pdaughL1=(AliAODMCParticle*)fmcArray->At(pdaugh1->GetDaughter(0));
      AliAODMCParticle *pdaughL2=(AliAODMCParticle*)fmcArray->At(pdaugh1->GetDaughter(1));
      if(!pdaughL1) return kFALSE;
      if(!pdaughL2) return kFALSE;
      if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac++;
    }
    if(number1==211 && number2==3124){
      if(fResonantDecay!=AliCFTaskVertexingHF::kL1520 && fResonantDecay!=AliCFTaskVertexingHF::kAll)return kFALSE; 
      nfiglieLa=pdaugh2->GetNDaughters();
      if(nfiglieLa!=2) return kFALSE;
      AliAODMCParticle *pdaughL1=(AliAODMCParticle*)fmcArray->At(pdaugh2->GetDaughter(0));
      AliAODMCParticle *pdaughL2=(AliAODMCParticle*)fmcArray->At(pdaugh2->GetDaughter(1));
      if(!pdaughL1) return kFALSE;
      if(!pdaughL2) return kFALSE;
      if((TMath::Abs(pdaughL1->GetPdgCode())==321 && TMath::Abs(pdaughL2->GetPdgCode())==2212) || (TMath::Abs(pdaughL1->GetPdgCode())==2212 && TMath::Abs(pdaughL2->GetPdgCode())==321)) numberOfLambdac++;
      
    }
  }
  if(numberOfLambdac>0) return kTRUE;
  return kFALSE;
}
