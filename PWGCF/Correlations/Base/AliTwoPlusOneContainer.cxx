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

//
//
// encapsulates several AliUEHist objects for a full UE analysis plus additional control histograms
// the method FillCorrelations is the only method which fills the container
// associated particles are filled once for each near side and once for each away side
// the pT 2 for the near side particles is set to the highest pT of all accepted away side pT 2
// as consequence the pT 2 axis can only be cut for the minimum
// if the maximum of the axis should be changed the full analysis needs to rerun
//
//
// Author: Markus Zimmermann
// 


#include "AliTwoPlusOneContainer.h"

#include "TList.h"

#include "AliCFContainer.h"
#include "AliVParticle.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

ClassImp(AliTwoPlusOneContainer)

AliTwoPlusOneContainer::AliTwoPlusOneContainer(const char* name, const char* uEHist_name, const char* binning, Double_t alpha) : 
  TNamed(name, name),
  fTwoPlusOne(0),
  fAsymmetry(0),
  fAsymmetryMixed(0),
  fTriggerPt(0),
  fTriggerPt1Min(0),
  fTriggerPt1Max(0),
  fTriggerPt2Min(0),
  fTriggerPt2Max(0),
  fPtAssocMin(0),
  fPtAssocMax(0),
  fAlpha(alpha),
  fUseLeadingPt(1),
  fUseAllT1(1),
  fUseBackgroundSameOneSide(0),
  fUseSmallerPtAssoc(0),
  fEfficiencyCorrection(0),
  fMergeCount(0)
{
  // Constructor
  //
  TString defaultBinningStr;
  defaultBinningStr = "eta: -1.0, 1.0\n"//needs to be defined for AliUEHist
    "p_t_assoc: 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0\n"
    "p_t_leading: 6.0, 8.0, 10.0, 12.0, 14.0\n"
    "p_t_leading_course: 4.0, 6.0, 8.0, 10.0\n"
    "p_t_eff: 0.5, 8.0\n"//needs to be defined for AliUEHist
    "vertex_eff: -10, 10\n"//needs to be defined for AliUEHist
    "multiplicity: 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1\n"
    "vertex: -7, -5, -3, -1, 1, 3, 5, 7\n"
    "delta_phi: -1.570796, -1.396263, -1.221730, -1.047197, -0.872664, -0.698131, -0.523598, -0.349065, -0.174532, 0.0, 0.174532, 0.349065, 0.523598, 0.698131, 0.872664, 1.047197, 1.221730, 1.396263, 1.570796, 1.745329, 1.919862, 2.094395, 2.268928, 2.443460, 2.617993, 2.792526, 2.967059, 3.141592, 3.316125, 3.490658, 3.665191, 3.839724, 4.014257, 4.188790, 4.363323, 4.537856, 4.712388\n" // this binning starts at -pi/2 and it ends at +3pi/2, it has 36 bins
    "delta_eta: -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8\n";

  // combine customBinning with defaultBinningStr -> use customBinning where available and otherwise defaultBinningStr
  TString binningStr = AliUEHist::CombineBinning(defaultBinningStr, TString(binning));

  fTwoPlusOne = new AliUEHist(uEHist_name, binningStr);

  //set minimum and maximum trigger pt values
  fTriggerPt1Min = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(2)->GetXmin();
  fTriggerPt1Max = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(2)->GetXmax();
  fTriggerPt2Min = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(6)->GetXmin();
  fTriggerPt2Max = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(6)->GetXmax();
  fPtAssocMin = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(1)->GetXmin();
  fPtAssocMax = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(1)->GetXmax();
  
  fAsymmetry  = new TH1F("fAsymmetry", ";A;dN", 50, 0, 1);
  fAsymmetryMixed  = new TH1F("fAsymmetryMixed", ";A;dN", 50, 0, 1);
  Int_t pt1_bins = (fTriggerPt1Max - fTriggerPt1Min)/0.1;
  Int_t pt2_bins = (fTriggerPt2Max - fTriggerPt2Min)/0.1;
  fTriggerPt = new TH2F("fTriggerPt", ";p_{T,1};p_{T,2}", pt1_bins, fTriggerPt1Min, fTriggerPt1Max, pt2_bins, fTriggerPt2Min, fTriggerPt2Max);
  
  TH1::AddDirectory();
}

//_____________________________________________________________________________
AliTwoPlusOneContainer::AliTwoPlusOneContainer(const AliTwoPlusOneContainer &c) :
  TNamed(c),
  fTwoPlusOne(0),
  fAsymmetry(0),
  fAsymmetryMixed(0),
  fTriggerPt(0),
  fTriggerPt1Min(0),
  fTriggerPt1Max(0),
  fTriggerPt2Min(0),
  fTriggerPt2Max(0),
  fPtAssocMin(0),
  fPtAssocMax(0),
  fAlpha(0.2),
  fUseLeadingPt(1),
  fUseAllT1(1),
  fUseBackgroundSameOneSide(0),
  fUseSmallerPtAssoc(0),
  fEfficiencyCorrection(0),
  fMergeCount(0)
{
  //
  // AliTwoPlusOneContainer copy constructor
  //

  ((AliTwoPlusOneContainer &) c).Copy(*this);
  fTriggerPt1Min = ((AliTwoPlusOneContainer &) c).getTriggerPt1Min();
  fTriggerPt1Max = ((AliTwoPlusOneContainer &) c).getTriggerPt1Max();
  fTriggerPt2Min = ((AliTwoPlusOneContainer &) c).getTriggerPt2Min();
  fTriggerPt2Max = ((AliTwoPlusOneContainer &) c).getTriggerPt2Max();
  fPtAssocMin = ((AliTwoPlusOneContainer &) c).getPtAssocMin();
  fPtAssocMax = ((AliTwoPlusOneContainer &) c).getPtAssocMax();
  fAlpha = ((AliTwoPlusOneContainer &) c).fAlpha;
  fMergeCount = ((AliTwoPlusOneContainer &) c).fMergeCount;
}

//____________________________________________________________________
AliTwoPlusOneContainer::~AliTwoPlusOneContainer()
{
  // Destructor
  
  DeleteContainers();
}

//____________________________________________________________________
void AliTwoPlusOneContainer::DeleteContainers()
{
  if (fTwoPlusOne)
  {
    delete fTwoPlusOne;
    fTwoPlusOne = 0;
  }
  
  if(fAsymmetry)
  {
    delete fAsymmetry;
    fAsymmetry = 0;
  }

  if(fAsymmetryMixed)
  {
    delete fAsymmetryMixed;
    fAsymmetryMixed = 0;
  }

  if(fTriggerPt)
  {
    delete fTriggerPt;
    fTriggerPt = 0;
  }

  if (fEfficiencyCorrection)
  {
    delete fEfficiencyCorrection;
    fEfficiencyCorrection = 0;
  }
}


//____________________________________________________________________
void AliTwoPlusOneContainer::FillParticleDist(Double_t centrality, Float_t zVtx, TObjArray* particleDist, Double_t weight, Bool_t applyEfficiency)
{
  AliCFContainer* track_hist = fTwoPlusOne->GetTrackHist(AliUEHist::kToward);
  AliCFContainer* event_hist = fTwoPlusOne->GetEventHist();
  AliUEHist::CFStep stepUEHist = static_cast<AliUEHist::CFStep>(AliTwoPlusOneContainer::kParticleDist);


  for (Int_t i=0; i<particleDist->GetEntriesFast(); i++){

    AliVParticle* part = (AliVParticle*) particleDist->UncheckedAt(i);
    
    Double_t part_pt = part->Pt();

    if(part_pt<fPtAssocMin||part_pt>fPtAssocMax)
      continue;

    Double_t part_eta = part->Eta();
    Double_t part_phi = part->Phi();

    //have to bring the phi angle on the same range as the Delta phi, because it is saved in the delta phi axis
    if(part_phi>1.5*TMath::Pi())  part_phi-= TMath::TwoPi();
    else if(part_phi<-0.5*TMath::Pi())  part_phi+= TMath::TwoPi();

    Double_t vars[7];
    vars[0] = part_eta;
    vars[1] = part_pt;//the particle pT is saved in the associated pT in the histograms because the axis has the correct range
    vars[2] = (fTriggerPt1Max+fTriggerPt1Min)/2;;
    vars[3] = centrality;
    vars[4] = part_phi;
    vars[5] = zVtx;
    vars[6] = (fTriggerPt2Max+fTriggerPt2Min)/2;
    
    Double_t efficiency = 1.0;
    if(applyEfficiency)
      efficiency = getEfficiency(part_pt, part_eta, centrality, zVtx);

    track_hist->Fill(vars, stepUEHist, weight*efficiency);

  }

}


//____________________________________________________________________
Int_t AliTwoPlusOneContainer::FillCorrelations(Double_t centrality, Float_t zVtx, AliTwoPlusOneContainer::PlotKind step, TObjArray* triggerNear, TObjArray* triggerAway, TObjArray* assocNear, TObjArray* assocAway, Double_t weight, Bool_t is1plus1, Bool_t isBackgroundSame, Bool_t applyEfficiency)
{
  //Fill Correlations fills the UEHist fTwoPlusOne with the 2+1 correlation
  //the input variables centrality and zVtx are the centrality and the z vertex of the event
  //the 4 lists triggerNear, triggerAway, assocNear and assocAway are four lists of particles. For the same event analysis all for lists are identical. For the mixed event analysis assocNear and assocAway are from different events as the trigger lists and for the mixed combinatorics triggerNear and assocNear are from one event and triggerAway and assocAway are from another event.
  //several booleans in this container change the behaviour of the container:
  //fUseLeadingPt: decides if a particle is only accepted as trigger particle if it has the highest pT within an circle with the radius of alpha
  //fUseAllT1: in case multiple trigger 2 are accepted, the near side yield is filled for each found away side yield
  AliCFContainer* track_hist = fTwoPlusOne->GetTrackHist(AliUEHist::kToward);
  AliCFContainer* event_hist = fTwoPlusOne->GetEventHist();
  AliUEHist::CFStep stepUEHist = static_cast<AliUEHist::CFStep>(step);

  //in case of the computation of the background in the same event there are two possible positions: delta phi = +/- pi/2
  //both positions are used so the results could only be weighted with 0.5*weight
  if(isBackgroundSame && !fUseBackgroundSameOneSide)
     fAlpha*= 0.5;

  //return variable: found triggers
  Int_t found_triggers = 0;

  //efficiency value which is multiplied with the weight
  //this value is always adjusted to the correct value just before the usage
  Double_t efficiency = 1.;

  for (Int_t i=0; i<triggerNear->GetEntriesFast(); i++){
    AliVParticle* part = (AliVParticle*) triggerNear->UncheckedAt(i);
    
    Double_t part_pt = part->Pt();

    if(part_pt<fTriggerPt1Min || part_pt>fTriggerPt1Max)
      continue;

    Double_t part_eta = part->Eta();
    Double_t part_phi = part->Phi();

    //search for second trigger particle
    //entries is the number of particles in the event, so there can only be less particles within the trigger 2 conditions
    Int_t ind_found = 0;
    Int_t ind_max_found_pt = -1;
    Int_t triggerAway_entries = triggerAway->GetEntriesFast();
    AliVParticle* found_particle[triggerAway_entries];
    Double_t found_particle_efficiency[triggerAway_entries];

    Bool_t do_not_use_T1 = false;

    //in case only the leading pt of a jet should be used, check every particle on the trigger near side if it's closer than alpha and if it has a higher pt than trigger 1
    if(fUseLeadingPt){
      for (Int_t i2=0; i2<triggerNear->GetEntriesFast(); i2++){
	if(i==i2)
	  continue;

	AliVParticle* part_i2 = (AliVParticle*) triggerNear->UncheckedAt(i2);

	if(part_i2->Pt()<=part_pt)
	  continue;
    
	Double_t dphi_check = part_phi-part_i2->Phi(); 
	if(dphi_check>1.5*TMath::Pi()) dphi_check -= TMath::TwoPi();
	else if(dphi_check<-0.5*TMath::Pi()) dphi_check += TMath::TwoPi();

	if(TMath::Abs(dphi_check)<fAlpha){
	  do_not_use_T1 = true;
	  break;
	}
      }
    }

    //if there is a particle with higher energy than T1 closer than alpha to T1, do not use this T1
     if(do_not_use_T1)
      continue;

    //have to fake the away side triggers for the 1+1 analysis
    if(is1plus1){
      found_particle[ind_found] = part;//in 1plus1 use first trigger particle also as pseudo second trigger particle
      found_particle_efficiency[ind_found] = 1.0;
      ind_max_found_pt = ind_found;
      ind_found = 1;
    }else{
      //normal 2+1 analysis
      for (Int_t j=0; j<triggerAway_entries; j++){
	AliVParticle* part2 = (AliVParticle*) triggerAway->UncheckedAt(j);

	Double_t part2_pt = part2->Pt();
	//check if pT of trigger 2 has enough energy to be a trigger
	//maximum energy is checked later after checking this particle may have to much energy for trigger 1 or 2
	if(part2_pt<fTriggerPt2Min)
	  continue;

	// don't use the same particle (is in any case impossible because the Delta phi angle will be 0)
	if(part==part2){
	  continue;
	}
	
	Double_t dphi_triggers = part_phi-part2->Phi();
	if(dphi_triggers>1.5*TMath::Pi()) dphi_triggers -= TMath::TwoPi();
	else if(dphi_triggers<-0.5*TMath::Pi()) dphi_triggers += TMath::TwoPi();
      
	//if 2+1 analysis check if trigger particles have a delta phi = pi +/- alpha
	if(!isBackgroundSame)
	  dphi_triggers -= TMath::Pi();
	else if(isBackgroundSame){
	  //shift defined area of delta phi
	  if(dphi_triggers>TMath::Pi()) dphi_triggers -= TMath::TwoPi();
	
	  if(!fUseBackgroundSameOneSide){
	    //look at delta phi = +/- pi/2
	    if(dphi_triggers<0)
	      dphi_triggers += 0.5*TMath::Pi();
	    else if(dphi_triggers>0)
	      dphi_triggers -= 0.5*TMath::Pi();
	  }else if(fUseBackgroundSameOneSide){
	    dphi_triggers -= 0.5*TMath::Pi();
	  }
	}
	if(TMath::Abs(dphi_triggers)>fAlpha)
	  continue;

	//check if pT of trigger 2 is too high
	if(part2_pt>fTriggerPt2Max || part2_pt>=part_pt){
	  //pt of trigger 2 needs to be smaller than the pt of trigger 1 (to have an ordering if both pt are close to each other)
	  if(fUseLeadingPt){
	    do_not_use_T1 = true;
	    break;
	  }else
	    continue;
	}

	found_particle[ind_found] = part2;
	if(applyEfficiency)
	  found_particle_efficiency[ind_found] = getEfficiency(part2_pt, part2->Eta(), centrality, zVtx);
	if(ind_max_found_pt==-1 || part2_pt>found_particle[ind_max_found_pt]->Pt()) ind_max_found_pt = ind_found;
	ind_found++;
      }//end loop to search for the second trigger particle
    }//end if for 1+1 

    //if there is a particle with higher energy than T1 or max(T2) within Delta phi = pi +/- alpha to T1, do not use this T1
    if(do_not_use_T1)
      continue;

    //if no second trigger particle was found continue to search for the next first trigger particle
    if(ind_found==0)
      continue;
    else{
      found_triggers += ind_found;
    }

    //use only the highest energetic particle on the away side, if there is only 1 away side trigger this is already the case
    if(fUseLeadingPt && ind_found>1){
      found_particle[0] = found_particle[ind_max_found_pt];
      found_particle_efficiency[0] = found_particle_efficiency[ind_max_found_pt];
      ind_found=1;
      ind_max_found_pt = 0;
    }

    //define efficiency for the found trigger 1 particle
    Double_t part_efficiency = 1.0;
    if(applyEfficiency)
      part_efficiency = getEfficiency(part_pt, part_eta, centrality, zVtx);
  
    //the energy of the second trigger particle is set for the near side to the maximum energy of all trigger 2 particles on the away side
    // this leads to the fact that the number of accepted trigger combinations can be artificial smaller than the real number if there is a cut on the pT 2 energy from the top; cutting away the smallest energy of pT 2 is still save; this is the reason why it is not allowed to use a cut on the top pt of trigger particle 2
    //fill trigger particles
    if(ind_found>0){
      Double_t vars[4];
      vars[0] = part_pt;
      vars[1] = centrality;
      vars[2] = zVtx;
      vars[3] = found_particle[ind_max_found_pt]->Pt();
      if(is1plus1)
	vars[3] = (fTriggerPt2Max+fTriggerPt2Min)/2;

      if(is1plus1){
	if(applyEfficiency)
	  efficiency = part_efficiency;
	event_hist->Fill(vars, stepUEHist, weight*efficiency);//near side (one times)

      }else if(!is1plus1){
	if(!fUseAllT1){
	  if(applyEfficiency)
	    efficiency = part_efficiency*found_particle_efficiency[ind_max_found_pt];
	  event_hist->Fill(vars, stepUEHist, weight*efficiency);//near side (one times)
	}

	for(Int_t k=0; k< ind_found; k++){
	  vars[3] = found_particle[k]->Pt();
	  if(applyEfficiency)
	    efficiency = part_efficiency*found_particle_efficiency[k];

	  event_hist->Fill(vars, stepUEHist+1, weight*efficiency);//away side

	  if(fUseAllT1)
	    event_hist->Fill(vars, stepUEHist, weight*efficiency);//near side (for every away side)
	    
	}
      }

      //fill fTriggerPt only once, choosed kSameNS
      if(step==AliTwoPlusOneContainer::kSameNS)
	for(Int_t k=0; k< ind_found; k++)
	  fTriggerPt->Fill(part_pt, found_particle[k]->Pt());

      //fill asymmetry only for kSameNS and kMixedNS
      if(step==AliTwoPlusOneContainer::kSameNS||step==AliTwoPlusOneContainer::kMixedNS){
	for(Int_t k=0; k< ind_found; k++){
	  Float_t asymmetry = (part_pt-found_particle[k]->Pt())/(part_pt+found_particle[k]->Pt());
	  if(step==AliTwoPlusOneContainer::kSameNS){
	    fAsymmetry->Fill(asymmetry);
	  }else{
	    fAsymmetryMixed->Fill(asymmetry);
	  }
	}
      }
    }
    
    //add correlated particles on the near side
    for (Int_t k=0; k<assocNear->GetEntriesFast(); k++){
      AliVParticle* part3 = (AliVParticle*) assocNear->UncheckedAt(k);

      Double_t part3_pt = part3->Pt();
      if(part3_pt<fPtAssocMin || part3_pt>fPtAssocMax)
	continue;

      //do not add the trigger 1 particle
      if(part==part3)
	continue;

      //use only pT,assoc which is samller than the trigger pT
      if(fUseSmallerPtAssoc && part3_pt>=part_pt)
	continue;

      Double_t dphi_near = part_phi-part3->Phi(); 
      if(dphi_near>1.5*TMath::Pi()) dphi_near -= TMath::TwoPi();
      else if(dphi_near<-0.5*TMath::Pi()) dphi_near += TMath::TwoPi();

      Double_t deta_near = part_eta-part3->Eta();
      
      Double_t vars[7];
      vars[0] = deta_near;
      vars[1] = part3_pt;
      vars[2] = part_pt;
      vars[3] = centrality;
      vars[4] = dphi_near;
      vars[5] = zVtx;
      vars[6] = found_particle[ind_max_found_pt]->Pt();
      if(is1plus1)
	vars[6] = (fTriggerPt2Max+fTriggerPt2Min)/2;

      Double_t part3_efficiency = 1.;
      if(applyEfficiency)
	part3_efficiency = getEfficiency(part3_pt, part3->Eta(), centrality, zVtx);

      if(is1plus1){
	if(applyEfficiency)
	  efficiency = part_efficiency*part3_efficiency;

	track_hist->Fill(vars, stepUEHist, weight*efficiency);
      }else if(!is1plus1){
	if(!fUseAllT1){
	  //do not add the trigger 2 particle with the highest pT
	  if(found_particle[ind_max_found_pt]==part3)
	    continue;

	  if(applyEfficiency)
	    efficiency = part_efficiency*found_particle_efficiency[ind_max_found_pt]*part3_efficiency;
	  
	  track_hist->Fill(vars, stepUEHist, weight*efficiency);
	}else
	  for(int l=0; l<ind_found; l++){
	    //do not add the trigger 2 particle
	    if(found_particle[l]==part3)
	      continue;
	    
	    vars[6] = found_particle[l]->Pt();
	    if(applyEfficiency)
	      efficiency = part_efficiency*found_particle_efficiency[l]*part3_efficiency;

	    track_hist->Fill(vars, stepUEHist, weight*efficiency);//fill NS for all AS triggers
	  }
      }
    }

    //search only for the distribution of the 2nd trigger particle
    if(is1plus1)
      continue;

    //add correlated particles on the away side
    for (Int_t k=0; k<assocAway->GetEntriesFast(); k++){
      AliVParticle* part3 = (AliVParticle*) assocAway->UncheckedAt(k);
	
      Double_t part3_pt = part3->Pt();
      if(part3_pt<fPtAssocMin || part3_pt>fPtAssocMax)
	continue;

      //do not add the trigger 1 particle
      if(part==part3)
	continue;

      Double_t part3_efficiency = 1.;
      if(applyEfficiency)
	part3_efficiency = getEfficiency(part3_pt, part3->Eta(), centrality, zVtx);

      for(int l=0; l<ind_found; l++){
	//do not add the trigger 2 particle
	if(found_particle[l]==part3)
	  continue;

	//use only pT,assoc which is samller than the trigger pT
	if(fUseSmallerPtAssoc && part3_pt>=found_particle[l]->Pt())
	  continue;

	Double_t dphi_away = found_particle[l]->Phi()-part3->Phi();
	if(dphi_away>1.5*TMath::Pi()) dphi_away -= TMath::TwoPi();
	else if(dphi_away<-0.5*TMath::Pi()) dphi_away += TMath::TwoPi();
	
	Double_t deta_away = found_particle[l]->Eta()-part3->Eta();
      
	Double_t vars[7];
	vars[0] = deta_away;
	vars[1] = part3_pt;
	vars[2] = part_pt;
	vars[3] = centrality;
	vars[4] = dphi_away;
	vars[5] = zVtx;
	vars[6] = found_particle[l]->Pt();

	if(applyEfficiency)
	  efficiency = part_efficiency*found_particle_efficiency[l]*part3_efficiency;

	track_hist->Fill(vars, stepUEHist+1, weight*efficiency);//step +1 is the AS to the NS plot of step
      }
    }
  }//end loop to search for the first trigger particle

  //put fAlpha back on the old value in case this is for background same
  if(isBackgroundSame && !fUseBackgroundSameOneSide)
     fAlpha*= 2;

  return found_triggers;
}

//____________________________________________________________________
AliTwoPlusOneContainer &AliTwoPlusOneContainer::operator=(const AliTwoPlusOneContainer &c)
{
  // assigment operator

  DeleteContainers();

  if (this != &c)
    ((AliTwoPlusOneContainer &) c).Copy(*this);

  return *this;
}

//____________________________________________________________________
void AliTwoPlusOneContainer::Copy(TObject& c) const
{
  // copy function

  AliTwoPlusOneContainer& target = (AliTwoPlusOneContainer &) c;

  if (fTwoPlusOne)
    target.fTwoPlusOne = dynamic_cast<AliUEHist*> (fTwoPlusOne->Clone());
  
  if (fAsymmetry)
    target.fAsymmetry = dynamic_cast<TH1F*> (fAsymmetry->Clone());

  if (fAsymmetryMixed)
    target.fAsymmetryMixed = dynamic_cast<TH1F*> (fAsymmetryMixed->Clone());

  if (fTriggerPt)
    target.fTriggerPt = dynamic_cast<TH2F*> (fTriggerPt->Clone());
  
  if (fEfficiencyCorrection)
    target.fEfficiencyCorrection = dynamic_cast<THnF*> (fEfficiencyCorrection->Clone());
  
}

//____________________________________________________________________
Long64_t AliTwoPlusOneContainer::Merge(TCollection* list)
{
  // Merge a list of AliTwoPlusOneContainer objects with this. 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;
  
  // collections of objects
  const Int_t kMaxLists = 4;
  TList* lists[kMaxLists];

  for (Int_t i=0; i<kMaxLists; i++)
    lists[i] = new TList;

  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliTwoPlusOneContainer* entry = dynamic_cast<AliTwoPlusOneContainer*> (obj);
    if (entry == 0) 
      continue;

    lists[0]->Add(entry->fTwoPlusOne);
    lists[1]->Add(entry->fAsymmetry);
    lists[2]->Add(entry->fAsymmetryMixed);
    lists[3]->Add(entry->fTriggerPt);

    fMergeCount += entry->fMergeCount;
    count++;
  }
  
  fTwoPlusOne->Merge(lists[0]);
  fAsymmetry->Merge(lists[1]);
  fAsymmetryMixed->Merge(lists[2]);
  fTriggerPt->Merge(lists[3]);

  for (Int_t i=0; i<kMaxLists; i++)
  delete lists[i];

  //  delete lists;
  return count+1;

}

//____________________________________________________________________
Double_t AliTwoPlusOneContainer::getEfficiency(Double_t pt, Double_t eta, Double_t centrality, Double_t zVtx){
  //returns the efficency correction for this particle
  
    Int_t effVars[4];
    effVars[0] = fEfficiencyCorrection->GetAxis(0)->FindBin(eta);
    effVars[1] = fEfficiencyCorrection->GetAxis(1)->FindBin(pt); //pt
    effVars[2] = fEfficiencyCorrection->GetAxis(2)->FindBin(centrality); //centrality
    effVars[3] = fEfficiencyCorrection->GetAxis(3)->FindBin(zVtx); //zVtx

    return fEfficiencyCorrection->GetBinContent(effVars);
  
  return 1.;
}
