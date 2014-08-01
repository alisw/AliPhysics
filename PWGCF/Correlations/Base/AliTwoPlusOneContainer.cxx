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
#include "TMath.h"

ClassImp(AliTwoPlusOneContainer)

AliTwoPlusOneContainer::AliTwoPlusOneContainer(const char* name, const char* binning, Double_t alpha) : 
  TNamed(name, name),
  fTwoPlusOne(0),
  fTriggerPt1Min(0),
  fTriggerPt1Max(0),
  fTriggerPt2Min(0),
  fTriggerPt2Max(0),
  fPtAssocMin(0),
  fPtAssocMax(0),
  fAlpha(alpha),
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

  fTwoPlusOne = new AliUEHist("TwoPlusOne", binningStr);

  //set minimum and maximum trigger pt values
  fTriggerPt1Min = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(2)->GetXmin();
  fTriggerPt1Max = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(2)->GetXmax();
  fTriggerPt2Min = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(6)->GetXmin();
  fTriggerPt2Max = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(6)->GetXmax();
  fPtAssocMin = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(1)->GetXmin();
  fPtAssocMax = fTwoPlusOne->GetTrackHist(AliUEHist::kToward)->GetGrid((AliUEHist::CFStep) AliTwoPlusOneContainer::kSameNS)->GetAxis(1)->GetXmax();

  TH1::AddDirectory();
}

//_____________________________________________________________________________
AliTwoPlusOneContainer::AliTwoPlusOneContainer(const AliTwoPlusOneContainer &c) :
  TNamed(fName, fTitle),
  fTwoPlusOne(0),
  fTriggerPt1Min(0),
  fTriggerPt1Max(0),
  fTriggerPt2Min(0),
  fTriggerPt2Max(0),
  fPtAssocMin(0),
  fPtAssocMax(0),
  fAlpha(0.2),
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
}


//____________________________________________________________________
void AliTwoPlusOneContainer::FillCorrelations(Double_t centrality, Float_t zVtx, AliTwoPlusOneContainer::PlotKind step, TObjArray* triggerNear, TObjArray* triggerAway, TObjArray* assocNear, TObjArray* assocAway, Double_t weight, Bool_t is1plus1, Bool_t isBackgroundSame)
{
  //Fill Correlations fills the UEHist fTwoPlusOne with the 2+1 correlation
  //the input variables centrality and zVtx are the centrality and the z vertex of the event
  //the 4 lists triggerNear, triggerAway, assocNear and assocAway are four lists of particles. For the same event analysis all for lists are identical. For the mixed event analysis assocNear and assocAway are from different events as the trigger lists and for the mixed combinatorics triggerNear and assocNear are from one event and triggerAway and assocAway are from another event.

  AliCFContainer* track_hist = fTwoPlusOne->GetTrackHist(AliUEHist::kToward);
  AliCFContainer* event_hist = fTwoPlusOne->GetEventHist();
  AliUEHist::CFStep stepUEHist = static_cast<AliUEHist::CFStep>(step);

  //in case of the computation of the background in the same event there are two possible positions: delta phi = +/- pi/2
  //both positions are used so the results could only be weighted with 0.5*weight
  if(isBackgroundSame)
    weight *= 0.5;

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

    for (Int_t j=0; j<triggerAway_entries; j++){
      AliVParticle* part2 = (AliVParticle*) triggerAway->UncheckedAt(j);

      Double_t part2_pt = part2->Pt();
      //check if pT of trigger 2 is within the trigger range
      //also pt of trigger 2 needs to be smaller than the pt of trigger 1 (to have an ordering if both pt are close to each other)
      if(part2_pt<fTriggerPt2Min || part2_pt>fTriggerPt2Max || part2_pt>part_pt)
	continue;
      
       // don't use the same particle (is in any case impossible because the Delta phi angle will be 0)
      if(part==part2){
	continue;
      }
	
      Double_t dphi_triggers = part_phi-part2->Phi();

      if(dphi_triggers>1.5*TMath::Pi()) dphi_triggers -= TMath::TwoPi();
      else if(dphi_triggers<-0.5*TMath::Pi()) dphi_triggers += TMath::TwoPi();
      
      //if 2+1 analysis check if trigger particles have a delta phi = pi +/- alpha
      if(!is1plus1 && !isBackgroundSame)
	dphi_triggers -= TMath::Pi();
      else if(!is1plus1 && isBackgroundSame){
	//shift defined area of delta phi
	if(dphi_triggers>TMath::Pi()) dphi_triggers -= TMath::TwoPi();

	//look at delta phi = +/- pi/2
	if(TMath::Abs(dphi_triggers)<0)
	  dphi_triggers += 0.5*TMath::Pi();
	else if(TMath::Abs(dphi_triggers)>0)
	  dphi_triggers -= 0.5*TMath::Pi();
      }
      if(!is1plus1 && TMath::Abs(dphi_triggers)>fAlpha)
	continue;

      found_particle[ind_found] = part2;
      if(ind_max_found_pt==-1 || part2_pt>found_particle[ind_max_found_pt]->Pt()) ind_max_found_pt = ind_found;
      ind_found++;

      if(is1plus1){
	Double_t vars[7];
	vars[0] = part_eta-part2->Eta();
	vars[1] = part2_pt;
	vars[2] = part_pt;
	vars[3] = centrality;
	vars[4] = dphi_triggers;//in case of is1plus1 this did not change
	vars[5] = zVtx;
	vars[6] = part2_pt;
	
	track_hist->Fill(vars, stepUEHist, weight);
      }

    }//end loop to search for the second trigger particle

    //if no second trigger particle was found continue to search for the next first trigger particle
    if(ind_found==0)
      continue;

    //the energy of the second trigger particle is set for the near side to the maximum energy of all trigger 2 particles on the away side
    // this leads to the fact that the number of accepted trigger combinations can be artificial smaller than the real number if there is a cut on the pT 2 energy from the top; cutting away the smallest energy of pT 2 is still save; this is the reason why it is not allowed to use a cut on the top pt of trigger particle 2
    //fill trigger particles
    if(ind_found>0){
	Double_t vars[4];
	vars[0] = part_pt;
	vars[1] = centrality;
	vars[2] = zVtx;
	vars[3] = found_particle[ind_max_found_pt]->Pt();

	event_hist->Fill(vars, stepUEHist, weight);//near side

	if(!is1plus1)
	  for(Int_t k=0; k< ind_found; k++){
	    vars[3] = found_particle[k]->Pt();
	    event_hist->Fill(vars, stepUEHist+1, weight);//away side
	  }
    }

    //search only for the distribution of the 2nd trigger particle
    if(is1plus1)
      continue;
  
    //add correlated particles on the near side
    for (Int_t k=0; k<assocNear->GetEntriesFast(); k++){
      AliVParticle* part3 = (AliVParticle*) assocNear->UncheckedAt(k);

      Double_t part3_pt = part3->Pt();
      if(part3_pt<fPtAssocMin || part3_pt>fPtAssocMax)
	continue;

      //do not add the trigger 1 particle
      if(part==part3)
	continue;

      //do not add the trigger 2 particle
      if(found_particle[ind_max_found_pt]==part3)
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

      track_hist->Fill(vars, stepUEHist, weight);
    }

    //add correlated particles on the away side
    for (Int_t k=0; k<assocAway->GetEntriesFast(); k++){
      AliVParticle* part3 = (AliVParticle*) assocAway->UncheckedAt(k);
	
      Double_t part3_pt = part3->Pt();
      if(part3_pt<fPtAssocMin || part3_pt>fPtAssocMax)
	continue;

      //do not add the trigger 1 particle
      if(part==part3)
	continue;

      for(int l=0; l<ind_found; l++){
	//do not add the trigger 2 particle
	if(found_particle[l]==part3)
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
	
	track_hist->Fill(vars, stepUEHist+1, weight);//step +1 is the AS to the NS plot of step
      }
    }
  }//end loop to search for the first trigger particle
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
  TList* lists = new TList;

  Int_t count = 0;
  while ((obj = iter->Next())) {
    
    AliTwoPlusOneContainer* entry = dynamic_cast<AliTwoPlusOneContainer*> (obj);
    if (entry == 0) 
      continue;

    lists->Add(entry->fTwoPlusOne);

    fMergeCount += entry->fMergeCount;
    count++;
  }
  
  fTwoPlusOne->Merge(lists);

  delete lists;
  return count+1;

}
