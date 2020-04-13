#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODConversionPhoton.h"
#include "AliKFConversionPhoton.h"

using namespace std;

ClassImp(AliAODConversionPhoton)

AliAODConversionPhoton::AliAODConversionPhoton() :
AliAODConversionParticle(),
AliConversionPhotonBase(),
  fLeadingNeutralPionIndex(-1),
  fLeadingNeutralPionDaughterIndex(-1),
  fCaloClusterRef(-1),
  fDCArPrimVtx(0),
  fDCAzPrimVtx(0),
  fInvMassPair(0),
  fNCaloPhotonMCLabels(0),
  fNCaloPhotonMotherMCLabels(0),
  fNNeutralPionLabels(0),
  fCaloPhotonMCFlags(0),
  fPairedId(-1),
  fCaloPhoton(0),
  fUseForMesonPair(kTRUE)
{
  // initialize calo photon MC labels
  for (Int_t i =0; i<50; i++){
    fCaloPhotonMCLabels[i]=-1;
  }
  for (Int_t i =0; i<20; i++){
    fCaloPhotonMotherMCLabels[i]=-1;
    fNeutralPionLabels[i]=-1;
    fNeutralPionEnergyFraction[i]=0;
  }
  //Standard constructor
}

AliAODConversionPhoton::AliAODConversionPhoton(AliKFConversionPhoton *kfphoton) :
AliAODConversionParticle(kfphoton),
AliConversionPhotonBase(*((AliConversionPhotonBase*)kfphoton)),
  fLeadingNeutralPionIndex(-1),
  fLeadingNeutralPionDaughterIndex(-1),
  fCaloClusterRef(-1),
  fDCArPrimVtx(0),
  fDCAzPrimVtx(0),
  fInvMassPair(0),
  fNCaloPhotonMCLabels(0),
  fNCaloPhotonMotherMCLabels(0),
  fNNeutralPionLabels(0),
  fCaloPhotonMCFlags(0),
  fPairedId(-1),
  fCaloPhoton(0),
  fUseForMesonPair(kTRUE)
{
  //Constructor from kfphoton
  // puts the mass to zero and store dilepton mass
  SetMass(kfphoton->M());

  //SetE(P());
  // initialize calo photon MC labels
  for (Int_t i =0; i<50; i++){
    fCaloPhotonMCLabels[i]=-1;
  }
  for (Int_t i =0; i<20; i++){
    fCaloPhotonMotherMCLabels[i]=-1;
    fNeutralPionLabels[i]=-1;
    fNeutralPionEnergyFraction[i]=0;
  }

}

AliAODConversionPhoton::AliAODConversionPhoton(TLorentzVector *vec) :
AliAODConversionParticle(vec),
AliConversionPhotonBase(),
  fLeadingNeutralPionIndex(-1),
  fLeadingNeutralPionDaughterIndex(-1),
  fCaloClusterRef(-1),
  fDCArPrimVtx(0),
  fDCAzPrimVtx(0),
  fInvMassPair(0),
  fNCaloPhotonMCLabels(0),
  fNCaloPhotonMotherMCLabels(0),
  fNNeutralPionLabels(0),
  fCaloPhotonMCFlags(0),
  fPairedId(-1),
  fCaloPhoton(kFALSE),
  fUseForMesonPair(kTRUE)
{
  //Constructor from TLorentzVector

  // initialize calo photon MC labels
  for (Int_t i =0; i<50; i++){
    fCaloPhotonMCLabels[i]=-1;
  }
  for (Int_t i =0; i<20; i++){
    fCaloPhotonMotherMCLabels[i]=-1;
    fNeutralPionLabels[i]=-1;
    fNeutralPionEnergyFraction[i]=0;
  }
}



AliAODConversionPhoton::AliAODConversionPhoton(const AliAODConversionPhoton & original) :
AliAODConversionParticle(original),
AliConversionPhotonBase(original),
  fLeadingNeutralPionIndex(original.fLeadingNeutralPionIndex),
  fLeadingNeutralPionDaughterIndex(original.fLeadingNeutralPionDaughterIndex),
  fCaloClusterRef(original.fCaloClusterRef),
  fDCArPrimVtx(original.fDCArPrimVtx),
  fDCAzPrimVtx(original.fDCAzPrimVtx),
  fInvMassPair(original.fInvMassPair),
  fNCaloPhotonMCLabels(original.fNCaloPhotonMCLabels),
  fNCaloPhotonMotherMCLabels(original.fNCaloPhotonMotherMCLabels),
  fNNeutralPionLabels(original.fNNeutralPionLabels),
  fCaloPhotonMCFlags(original.fCaloPhotonMCFlags),
  fPairedId(original.fPairedId),
  fCaloPhoton(original.fCaloPhoton),
  fUseForMesonPair(original.fUseForMesonPair)
{
  //Copy constructor

  // initialize calo photon MC labels
  for (Int_t i =0; i<50; i++){
    fCaloPhotonMCLabels[i]=original.fCaloPhotonMCLabels[i];
  }
  for (Int_t i =0; i<20; i++){
      fCaloPhotonMotherMCLabels[i]=original.fCaloPhotonMotherMCLabels[i];
      fNeutralPionLabels[i]=original.fNeutralPionLabels[i];
      fNeutralPionEnergyFraction[i]=original.fNeutralPionEnergyFraction[i];
  }
}

AliAODConversionPhoton::~AliAODConversionPhoton()
{
  // empty standard destructor
}

AliAODConversionPhoton & AliAODConversionPhoton::operator = (const AliAODConversionPhoton & /*source*/)
{
  // assignment operator
  return *this;
}

///________________________________________________________________________
void AliAODConversionPhoton::CalculateDistanceOfClossetApproachToPrimVtx(const AliVVertex* primVertex ){

  Double_t primCo[3] = {primVertex->GetX(),primVertex->GetY(),primVertex->GetZ()};

  Double_t absoluteP = TMath::Sqrt(TMath::Power(GetPx(),2) + TMath::Power(GetPy(),2) + TMath::Power(GetPz(),2));
  Double_t p[3] = {GetPx()/absoluteP,GetPy()/absoluteP,GetPz()/absoluteP};
  Double_t CP[3];

  CP[0] =  fConversionPoint[0] - primCo[0];
  CP[1] =  fConversionPoint[1] - primCo[1];
  CP[2] =  fConversionPoint[2] - primCo[2];

  Double_t Lambda = - (CP[0]*p[0]+CP[1]*p[1]+CP[2]*p[2])/(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);

  Double_t S[3];
  S[0] = fConversionPoint[0] + p[0]*Lambda;
  S[1] = fConversionPoint[1] + p[1]*Lambda;
  S[2] = fConversionPoint[2] + p[2]*Lambda;

  fDCArPrimVtx = TMath::Sqrt( TMath::Power(primCo[0]-S[0],2) + TMath::Power(primCo[1]-S[1],2));
  fDCAzPrimVtx = primCo[2]-S[2];

  return;
}


void AliAODConversionPhoton::SetCaloPhotonMCFlags(AliMCEvent *mcEvent, Bool_t enableSort, Bool_t mergedAnalysis, AliVCluster* cluster){


  TParticle* PhotonDummyMerged;
  TParticle* PhotonDummyMergedMother;
  TParticle* PhotonDummyMergedGrandMother;
  TParticle* PhotonDummyMergedGGMother;
  Int_t photonDummyMergedPDG= -1;
  Int_t photonDummyMergedMotherPDG= -1;
  Int_t photonDummyMergedGrandMotherPDG= -1;
  Int_t photonDummyMergedGGMotherPDG= -1;
  Bool_t foundNeutralPion = kFALSE;
  Int_t neutralPionLabel = -1;
  if(mergedAnalysis){
    for (Int_t j = 0; j< fNCaloPhotonMCLabels; j++){
      neutralPionLabel = -1;
      foundNeutralPion = kFALSE;
      PhotonDummyMerged        = mcEvent->Particle(GetCaloPhotonMCLabel(j)); // main particle
      photonDummyMergedPDG = PhotonDummyMerged->GetPdgCode();
      if(TMath::Abs(photonDummyMergedPDG)==111){
        foundNeutralPion = kTRUE;
        neutralPionLabel = GetCaloPhotonMCLabel(j);
      }
      if(PhotonDummyMerged->GetMother(0)>-1 && !foundNeutralPion){
        PhotonDummyMergedMother  = mcEvent->Particle(PhotonDummyMerged->GetMother(0)); // mother
        photonDummyMergedMotherPDG = PhotonDummyMergedMother->GetPdgCode();
        if(TMath::Abs(photonDummyMergedMotherPDG)==111){
          foundNeutralPion = kTRUE;
          neutralPionLabel = PhotonDummyMerged->GetMother(0);
        }
        if(PhotonDummyMergedMother->GetMother(0)>-1 && !foundNeutralPion){
          PhotonDummyMergedGrandMother  = mcEvent->Particle(PhotonDummyMergedMother->GetMother(0)); // grandmother
          photonDummyMergedGrandMotherPDG = PhotonDummyMergedGrandMother->GetPdgCode();
          if(TMath::Abs(photonDummyMergedGrandMotherPDG)==111){
            foundNeutralPion = kTRUE;
            neutralPionLabel = PhotonDummyMergedMother->GetMother(0);
          }
          if(PhotonDummyMergedGrandMother->GetMother(0)>-1 && !foundNeutralPion){
            PhotonDummyMergedGGMother  = mcEvent->Particle(PhotonDummyMergedGrandMother->GetMother(0)); // grand-grandmother
            photonDummyMergedGGMotherPDG = PhotonDummyMergedGGMother->GetPdgCode();
            if(TMath::Abs(photonDummyMergedGGMotherPDG)==111){
              foundNeutralPion = kTRUE;
              neutralPionLabel = PhotonDummyMergedGrandMother->GetMother(0);
            }
          }
        }
      }
       // check if current found pi0 was already found before in the cluster
      // if it was found before, add the cluster energy fraction of the current label to this pi0
      if(foundNeutralPion){
        Bool_t newNeutralPion                    = true;
        for(Int_t k=0; k<fNNeutralPionLabels; k++){
          if (fNeutralPionLabels[k] == neutralPionLabel){ //check if mother is already contained in fNeutralPionLabels
            newNeutralPion                       = false;
            fNeutralPionEnergyFraction[k] += cluster->GetClusterMCEdepFraction(j);
          }
        }
        // if the pi0 is new, then increase pi0 count by one and set E-fraction and label to array
        if (newNeutralPion){
          fNeutralPionLabels[fNNeutralPionLabels]    = neutralPionLabel;
          fNeutralPionEnergyFraction[fNNeutralPionLabels] = cluster->GetClusterMCEdepFraction(j);
          fNNeutralPionLabels++; //only if particle label is not yet contained in array, count up fNeutralPionLabels
        }
      }
    }

    // search for the largest contributing pi0 in the cluster
    Double_t largestEfraction = 0;
    for(Int_t k=0; k<fNNeutralPionLabels; k++){
      if(fNeutralPionEnergyFraction[k]>largestEfraction){
        largestEfraction = fNeutralPionEnergyFraction[k];
        fLeadingNeutralPionIndex = k;
      }
    }
    // if the cluster contains a pi0 contribution, find the MC label/index of the largest contributing daughter
    // for this, go upwards three mothers from each MC label to find leading daughter index
    if(fLeadingNeutralPionIndex>-1){
      for(Int_t l=0; l<fNCaloPhotonMCLabels; l++){
        if(fLeadingNeutralPionDaughterIndex>-1)
          continue;
        PhotonDummyMerged        = mcEvent->Particle(GetCaloPhotonMCLabel(l)); // main particle
        if(PhotonDummyMerged->GetMother(0)>-1){
          if(PhotonDummyMerged->GetMother(0)==fNeutralPionLabels[fLeadingNeutralPionIndex]){
            fLeadingNeutralPionDaughterIndex = l;
          }else{
            PhotonDummyMergedMother        = mcEvent->Particle(PhotonDummyMerged->GetMother(0));
            if(PhotonDummyMergedMother->GetMother(0)>-1){
              if(PhotonDummyMergedMother->GetMother(0)==fNeutralPionLabels[fLeadingNeutralPionIndex]){
                fLeadingNeutralPionDaughterIndex = l;
              }else{
                PhotonDummyMergedGrandMother        = mcEvent->Particle(PhotonDummyMergedMother->GetMother(0));
                if(PhotonDummyMergedGrandMother->GetMother(0)>-1){
                  if(PhotonDummyMergedGrandMother->GetMother(0)==fNeutralPionLabels[fLeadingNeutralPionIndex]){
                    fLeadingNeutralPionDaughterIndex = l;
                  }else{
                    PhotonDummyMergedGGMother        = mcEvent->Particle(PhotonDummyMergedGrandMother->GetMother(0));
                    if(PhotonDummyMergedGGMother->GetMother(0)>-1){
                      if(PhotonDummyMergedGGMother->GetMother(0)==fNeutralPionLabels[fLeadingNeutralPionIndex]){
                        fLeadingNeutralPionDaughterIndex = l;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  Bool_t isPhoton                   = kFALSE; // largest contribution to cluster is photon
  Bool_t isElectron                 = kFALSE; // largest contribution to cluster is electron
  Bool_t isConversion               = kFALSE; // largest contribution to cluster is converted electron
  Bool_t isConversionFullyContained = kFALSE; // largest contribution to cluster is converted electron, second electron has been found in same cluster
  Bool_t isMerged                   = kFALSE; // cluster contains more than one particle from the same decay
  Bool_t isMergedPartConv           = kFALSE; // cluster contains more than one particle from the same decay and at least one of the particles came from a conversion
  Bool_t isDalitz                   = kFALSE; // this cluster was created by a particle stemming from a dalitz decay
  Bool_t isDalitzMerged             = kFALSE; // this cluster was created by a particle stemming from a dalitz decay and more than one particle of the dalitz decay is contained in the cluster
  Bool_t isPhotonWithElecMother     = kFALSE; // this cluster is from a photon with an electron as mother
  Bool_t isShower                   = kFALSE; // this cluster contains as a largest contribution a particle from a shower or radiative process
  Bool_t isSubLeadingEM             = kFALSE; // cluster contains at least one electron or photon from a pi0, eta or eta_prime in subleading contribution
  Bool_t isElectronFromFragPhoton   = kFALSE; // largest contribution to cluster is from converted electron, but photon stems from fragmentation photon ( q -> q gamma)

  TParticle* Photon = 0x0;
  TParticle* PhotonMother = 0x0;
  TParticle* PhotonGrandMother = 0x0;
  if (fNCaloPhotonMCLabels==0) return;

  if (enableSort){
    // sort the array according to the energy of contributing particles
    if (fNCaloPhotonMCLabels>1){
      // cout << "start sorting" << endl;
      Int_t* sortIdx            = new Int_t[fNCaloPhotonMCLabels];
      Double_t* energyPerPart   = new Double_t[fNCaloPhotonMCLabels];
      Long_t* orginalContrib    = new Long_t[fNCaloPhotonMCLabels];
      for (Int_t i = 0; i < fNCaloPhotonMCLabels; i++){
        orginalContrib[i]       = fCaloPhotonMCLabels[i];
        if (fCaloPhotonMCLabels[i]> -1){
          TParticle* dummy  = mcEvent->Particle(fCaloPhotonMCLabels[i]);
          if (dummy){
            energyPerPart[i]  = dummy->Energy();
            // suppress energy of hadrons !!! DIRTY hack !!!
            if (!(TMath::Abs(dummy->GetPdgCode())== 11 || TMath::Abs(dummy->GetPdgCode())== 22)){
              // cout << "suppressed hadron energy for:" << dummy->GetPdgCode() << endl;
              energyPerPart[i]= 0.25; // put energy to mip
            }
          }
        } else {
          energyPerPart[i]  = 0;
        }
      }

      TMath::Sort(fNCaloPhotonMCLabels,energyPerPart,sortIdx);
      for(Int_t index = 0; index < fNCaloPhotonMCLabels; index++) {
        fCaloPhotonMCLabels[index] = orginalContrib  [sortIdx[index]] ;

      }
      delete [] sortIdx;
      delete [] energyPerPart;
      delete [] orginalContrib;
    }
  }

  if(mergedAnalysis && fNNeutralPionLabels>0 && fLeadingNeutralPionDaughterIndex!=0 && fNeutralPionEnergyFraction[fLeadingNeutralPionIndex]>cluster->GetClusterMCEdepFraction(0)){
    // check if leading pi0 comes not from label 0 in cluster
    // for this do:
    // -> check if neutral pions were found in cluster
    // -> if the leading daughter index is not 0
    // -> the leading neutral pion has a larger cluster energy fraction than the cluster label 0
    // load particle corresponding to largest daughter of leading pi0
    Photon                            = mcEvent->Particle(GetCaloPhotonMCLabel(fLeadingNeutralPionDaughterIndex));
  } else {
    // load particle corresponding to MC label 0 in cluster
    if(GetCaloPhotonMCLabel(0)>-1)
      Photon                            = mcEvent->Particle(GetCaloPhotonMCLabel(0));
  }


  if(Photon == NULL){
    return;
  }

  Int_t particleMotherLabel           = Photon->GetMother(0);
  Int_t particleGrandMotherLabel      = -1;
  Int_t particleGrandMotherX2Label    = -1;
  Int_t particleGrandMotherX3Label    = -1;
  Int_t particleGrandMotherX4Label    = -1;
  Int_t particleGrandMotherX5Label    = -1;
  Int_t particleGrandMotherX6Label    = -1;
  Int_t particleGrandMotherX7Label    = -1;

  Int_t particleMotherPDG             = -1;
  Int_t particleGrandMotherPDG        = -1;
  Int_t particleGrandMotherX2PDG      = -1;
  Int_t particleGrandMotherX3PDG      = -1;
  Int_t particleGrandMotherX4PDG      = -1;
  Int_t particleGrandMotherX5PDG      = -1;
  Int_t particleGrandMotherX6PDG      = -1;

  Int_t particleMotherNDaugthers      = 0;
  Int_t particleGrandMotherNDaugthers = 0;

  if (particleMotherLabel > -1){
    PhotonMother                      = mcEvent->Particle(particleMotherLabel);
    particleMotherNDaugthers          = PhotonMother->GetNDaughters();
    particleGrandMotherLabel          = PhotonMother->GetMother(0);
    particleMotherPDG                 = PhotonMother->GetPdgCode();
    if (particleGrandMotherLabel > -1){
      PhotonGrandMother               = mcEvent->Particle(particleGrandMotherLabel);
      particleGrandMotherPDG          = PhotonGrandMother->GetPdgCode();
      particleGrandMotherNDaugthers   = PhotonGrandMother->GetNDaughters();
      particleGrandMotherX2Label      = PhotonGrandMother->GetMother(0);
      TParticle* dummyGMM             = 0x0;
      if (particleGrandMotherX2Label > -1){
        dummyGMM                        = mcEvent->Particle(particleGrandMotherX2Label);
        particleGrandMotherX2PDG        = dummyGMM->GetPdgCode();
        particleGrandMotherX3Label      = dummyGMM->GetMother(0);
        if (particleGrandMotherX3Label > -1){
          dummyGMM                        = mcEvent->Particle(particleGrandMotherX3Label);
          particleGrandMotherX3PDG        = dummyGMM->GetPdgCode();
          particleGrandMotherX4Label      = dummyGMM->GetMother(0);
          if (particleGrandMotherX4Label > -1){
            dummyGMM                        = mcEvent->Particle(particleGrandMotherX4Label);
            particleGrandMotherX4PDG        = dummyGMM->GetPdgCode();
            particleGrandMotherX5Label      = dummyGMM->GetMother(0);
            if (particleGrandMotherX5Label > -1){
              dummyGMM                        = mcEvent->Particle(particleGrandMotherX5Label);
              particleGrandMotherX5PDG        = dummyGMM->GetPdgCode();
              particleGrandMotherX6Label      = dummyGMM->GetMother(0);
              if (particleGrandMotherX6Label > -1){
                dummyGMM                        = mcEvent->Particle(particleGrandMotherX6Label);
                particleGrandMotherX6PDG        = dummyGMM->GetPdgCode();
                particleGrandMotherX7Label      = dummyGMM->GetMother(0);
              }
            }
          }
        }
      }
    }
  }

  //determine mother/grandmother of leading particle and if it is pion/eta/eta_prime: fill array fCaloPhotonMotherMCLabels at position 0
  if (particleMotherLabel > -1){
    if( TMath::Abs(particleMotherPDG) == 111 || TMath::Abs(particleMotherPDG) == 221 || TMath::Abs(particleMotherPDG) == 331 ){
      fCaloPhotonMotherMCLabels[0]    = particleMotherLabel;
      fNCaloPhotonMotherMCLabels++;
    } else if (particleGrandMotherLabel > -1 &&  TMath::Abs(particleMotherPDG) == 22 &&
                (TMath::Abs(particleGrandMotherPDG) == 111 || TMath::Abs(particleGrandMotherPDG) == 221 || TMath::Abs(particleGrandMotherPDG) == 331) ){
      fCaloPhotonMotherMCLabels[0]  = particleGrandMotherLabel;
      fNCaloPhotonMotherMCLabels++;
    } else {
      fCaloPhotonMotherMCLabels[0]    = particleMotherLabel;
      fNCaloPhotonMotherMCLabels++;
    }
  }

  // Check whether the first contribution was photon
  if(TMath::Abs(Photon->GetPdgCode()) == 22){
    isPhoton                          = kTRUE;
    // did it decay via the dalitz channel
    if (particleMotherLabel > -1 && particleMotherNDaugthers == 3)
      isDalitz                        = kTRUE;
    // Test whether particle stems from a shower or radiation
    if (TMath::Abs(particleMotherPDG) == 11){                    // check whether photon stems from electron
      isPhotonWithElecMother          = kTRUE;
      if (particleGrandMotherLabel > -1){                 // test whether first particle has a grandmother
        if (TMath::Abs(particleGrandMotherPDG) == 22 )
          isShower                    = kTRUE;            // check whether grandmother is a photon (meaning this is most likely a shower)
      }
    }
  }

  // Check whether the first contribution was electron
  if( TMath::Abs(Photon->GetPdgCode()) == 11 ){
    isElectron                        = kTRUE;
    if (particleMotherLabel > -1) {
      // was it a conversion
      if (TMath::Abs(particleMotherPDG) == 22)
        isConversion                  = kTRUE;
      // did it decay via the dalitz channel
      if (particleGrandMotherLabel > -1 && particleGrandMotherNDaugthers == 3 )
        isDalitz                      = kTRUE;
    }
    if (particleGrandMotherLabel > -1){                     // check whether electron has a grandmother
      if (TMath::Abs(particleGrandMotherPDG) == 11 ||  TMath::Abs(particleGrandMotherPDG) == 22){ // test whether electron has photon or electron as grandmother (meaning will most likely be a shower)
        isShower                      = kTRUE;
      }
    }
    // consider the rare case, where the conversion electron stems from photon, which stems from electron, which stems from photon, ...
    if (isConversion && TMath::Abs(particleGrandMotherPDG) == 11   && particleGrandMotherX2PDG == 22){
      SetCaloPhotonMCLabel(0,particleGrandMotherLabel);
      fCaloPhotonMotherMCLabels[0] = particleGrandMotherX3Label;
      if (TMath::Abs(particleGrandMotherX3PDG) == 11 && particleGrandMotherX4PDG == 22  ){
        SetCaloPhotonMCLabel(0,particleGrandMotherX3Label);
        fCaloPhotonMotherMCLabels[0] = particleGrandMotherX5Label;
        if (TMath::Abs(particleGrandMotherX5PDG) == 11 && particleGrandMotherX6PDG == 22  ){
          SetCaloPhotonMCLabel(0,particleGrandMotherX5Label);
          fCaloPhotonMotherMCLabels[0] = particleGrandMotherX7Label;
        }
      }
    }

    // consider the case, where a photon stems from a photon which stems from a photon etc...which is common for frag. photons
    if (particleMotherLabel > -1){
      TParticle *dummyMother = mcEvent->Particle(particleMotherLabel);
      Bool_t originReached   = kFALSE;
      if (dummyMother){
        while (dummyMother->GetPdgCode() == 22 && !originReached){ // follow conversion photon's history, as long as the mother is a photon
          if (dummyMother->GetMother(0) > -1){
            dummyMother = mcEvent->Particle(dummyMother->GetMother(0));
            if ((TMath::Abs(dummyMother->GetPdgCode()) == 11) || (TMath::Abs(dummyMother->GetPdgCode()) == 22)){ // in case of additional conversion skip to photon's grandma, which should be a photon
              if (dummyMother->GetMother(0) > -1){                                                               // also mother of photon could be a photon (can happen for fragmentation)
                dummyMother   = mcEvent->Particle(dummyMother->GetMother(0));
              } else {
                originReached = kTRUE;
              }
            } else {
              originReached = kTRUE;
            }
            isElectronFromFragPhoton = (TMath::Abs(dummyMother->GetPdgCode()) < 6);// photon stems from quark = fragmentation photon
          } else {
            originReached = kTRUE;
          }
        }
      }
    }
  }

  Bool_t enablePrintOuts            = kFALSE;
  // if (fNCaloPhotonMCLabels>2)
  //   enablePrintOuts                 = kTRUE;

  // check whether there were other contributions to the cluster
  if (fNCaloPhotonMCLabels>1){
    TParticle* dummyPart            =NULL;
    TParticle* dummyPartMother      =NULL;
    TParticle* dummyPartGrandMother =NULL;

    for (Int_t i = 1; i< fNCaloPhotonMCLabels; i++){
      if (i > 49) continue;													// abort if more than 50 entries to the cluster have been checked (more are not stored in these objects)
      if (enablePrintOuts) cout << "checking particle: " <<  i << endl;
      if (GetCaloPhotonMCLabel(i) < 0) continue;
      dummyPart = mcEvent->Particle(GetCaloPhotonMCLabel(i));
      Int_t dummyPartMotherLabel      = dummyPart->GetMother(0);
      Int_t dummyPartGrandMotherLabel = -1;
      Int_t dummyPartMotherPDG        = -1;
      Int_t dummyPartGrandMotherPDG   = -1;
      // check whether this particle has a mother & obtain the pdg code
      if (dummyPartMotherLabel > -1){
        dummyPartMother               = mcEvent->Particle(dummyPartMotherLabel);
        dummyPartGrandMotherLabel     = dummyPartMother->GetMother(0);
        dummyPartMotherPDG            = dummyPartMother->GetPdgCode();
        // check whether this particle has a grandmother & obtain its pdg code
        if (dummyPartGrandMotherLabel > -1){
          dummyPartGrandMother        = mcEvent->Particle(dummyPartGrandMotherLabel);
          dummyPartGrandMotherPDG     = dummyPartGrandMother->GetPdgCode();
        }
      }
      // largest contribution was from photon and is not from shower or electron mother
      if (isPhoton && (!isShower || !isPhotonWithElecMother )){
        if (enablePrintOuts) cout << "lead gamma" << endl;
        if (dummyPartMotherLabel > -1 && particleMotherLabel > -1){   // test whether first particle has a mother
          if (dummyPartMotherLabel == particleMotherLabel)
            isMerged                  = kTRUE;                        // test whether current and first particle have the same mother => i.e. other gamma from decay or dalitz electron
          if (dummyPartGrandMotherLabel > -1){                        // test whether first particle has a grandmother
            // check whether particle is an electron from a conversion of a photon from the original mother
            if (TMath::Abs(dummyPart->GetPdgCode()) == 11 && dummyPartGrandMotherLabel == particleMotherLabel )
              isMergedPartConv        = kTRUE;
            // check whether particle is an electron from a dalitz decay from the original mother
            if (TMath::Abs(dummyPart->GetPdgCode()) == 11 && dummyPartMotherLabel == particleMotherLabel )
              isDalitzMerged          = kTRUE;
          }
        }
      }

      // largest contribution was from electron & not a from a shower
      if (isElectron && !isShower){
        if (enablePrintOuts) cout << "lead electron" << endl;
        if (dummyPartMotherLabel > -1 && particleMotherLabel > -1){     // test whether first particle has a mother
          if ( TMath::Abs(dummyPart->GetPdgCode()) == 11 && isConversion && dummyPartMotherLabel == particleMotherLabel  ) {
            isConversionFullyContained  = kTRUE;                        // test whether conversion is fully contained in cluster
            if (enablePrintOuts) cout << "found full conversion" << endl;
          } else if (TMath::Abs(dummyPart->GetPdgCode()) == 11) { // current particle is an electron
            if (enablePrintOuts) cout << "current is electron" << endl;
            // check whether current particle is a shower
            Bool_t isShowerCurrent    = kFALSE;
            if (TMath::Abs(dummyPartGrandMotherPDG) == 11)
              isShowerCurrent         = kTRUE;
            if (TMath::Abs(dummyPartGrandMotherPDG) == 22)
              isShowerCurrent         = kTRUE;

            // check whether current particle is a conversion
            Bool_t isConvCurrent      = kFALSE;
            if (TMath::Abs(dummyPartMotherPDG) == 22 && !isShowerCurrent)
              isConvCurrent           = kTRUE;

            if (enablePrintOuts) cout << "conv current: " << isConvCurrent << "\t shower current: " << isShowerCurrent << endl;
            // check whether orginal electron and this electron stem from the same particle and electron originated in dalitz decay
            if( dummyPartMotherLabel == particleMotherLabel && TMath::Abs(particleMotherPDG) != 22 &&  !isShowerCurrent   ) {
              isDalitzMerged          = kTRUE;
              isMerged                = kTRUE;
            // both particles are from conversions
            } else if (isConversion && isConvCurrent){
              if (particleGrandMotherLabel > -1 && dummyPartGrandMotherLabel > -1){ // test whether first particle has a grandmother and current particle has one as well
                // check whether orginal electron and this electron stem from the same particle and electron stems from conversion
                if( dummyPartGrandMotherLabel == particleGrandMotherLabel &&  !isShowerCurrent   )
                  isMergedPartConv        = kTRUE;
              }
            // check whether current electron is from conversion & orginal e are from same mother
            } else if (dummyPartGrandMotherLabel > -1 && isConvCurrent){
              if (dummyPartGrandMotherLabel == particleMotherLabel && (TMath::Abs(particleMotherPDG) == 111 || TMath::Abs(particleMotherPDG) == 221 || TMath::Abs(particleMotherPDG) == 331) ){
                isDalitzMerged          = kTRUE;
                isMerged                = kTRUE;
              }
            // check whether current electron & orginal e is from conversion are from same mother
            } else if (isConversion && particleGrandMotherLabel > -1){
              if (dummyPartMotherLabel == particleGrandMotherLabel && (TMath::Abs(dummyPartMotherPDG) == 111 || TMath::Abs(dummyPartMotherPDG) == 221 || TMath::Abs(dummyPartMotherPDG) == 331) ){
                isDalitzMerged          = kTRUE;
                isMerged                = kTRUE;
              }
            }
          } else if (TMath::Abs(dummyPart->GetPdgCode()) == 22){                 // test whether this particle is a photon
            // check whether current particle is a shower
            Bool_t isShowerCurrent    = kFALSE;
            if (TMath::Abs(dummyPartMotherPDG) == 11)
              isShowerCurrent         = kTRUE;
            if (TMath::Abs(dummyPartMotherPDG) == 22)
              isShowerCurrent         = kTRUE;

            // check whether orginal electron and this photon stem from the same particle and electron originated in dalitz
            if( dummyPartMotherLabel == particleMotherLabel && !isShowerCurrent ) {
              isDalitzMerged          = kTRUE;
              isMerged                = kTRUE;
            } else if (particleGrandMotherLabel > -1){ // test whether first particle has a grandmother
              // check whether orginal electron and this photon stem from the same particle and electron stems from conversion
              if( dummyPartMotherLabel == particleGrandMotherLabel && !isShowerCurrent )
                isMergedPartConv        = kTRUE;
            }
          }
        }
      }


      if (dummyPartMotherLabel > -1){ // test whether particle has a mother
        if (TMath::Abs(dummyPart->GetPdgCode()) == 22){ // test whether particle is a photon
          //check if photon directly comes from a pion/eta/eta_prime decay
          Bool_t helpN                    = true;
          for(Int_t j=0; j<fNCaloPhotonMotherMCLabels; j++){
            if (fCaloPhotonMotherMCLabels[j] == dummyPartMotherLabel){ //check if mother is already contained in fCaloPhotonMotherMCLabels
              helpN                       = false;
            }
          }
          if (helpN){
            fCaloPhotonMotherMCLabels[fNCaloPhotonMotherMCLabels]    = dummyPartMotherLabel;
            fNCaloPhotonMotherMCLabels++; //only if particle label is not yet contained in array, count up fNCaloPhotonMotherMCLabels
          }
          if ((TMath::Abs(dummyPartMotherPDG) == 111 || TMath::Abs(dummyPartMotherPDG) == 221 || TMath::Abs(dummyPartMotherPDG) == 331) && !isPhoton && !isElectron)
            isSubLeadingEM                = kTRUE;
        } else if (TMath::Abs(dummyPart->GetPdgCode()) == 11){ //test whether particle is an electron
          //check if electron comes from a pion decay
          if ( TMath::Abs(dummyPartMotherPDG) != 22 ){
            Bool_t helpN                    = true;
            for(Int_t j=0; j<fNCaloPhotonMotherMCLabels; j++){
              if (fCaloPhotonMotherMCLabels[j] == dummyPartMotherLabel){ //check if mother is already contained in fCaloPhotonMotherMCLabels
                helpN                       = false;
              }
            }
            if (helpN){
              fCaloPhotonMotherMCLabels[fNCaloPhotonMotherMCLabels]    = dummyPartMotherLabel;
              fNCaloPhotonMotherMCLabels++; //only if particle label is not yet contained in array, count up fNCaloPhotonMotherMCLabels
            }
            if ((TMath::Abs(dummyPartMotherPDG) == 111 || TMath::Abs(dummyPartMotherPDG) == 221 || TMath::Abs(dummyPartMotherPDG) == 331) && !isPhoton && !isElectron)
              isSubLeadingEM                = kTRUE;
          } else if (dummyPartGrandMotherLabel > -1){ //if it is not a dalitz decay, test whether particle has a grandmother
            //check if it is a conversion electron that has pion/eta/eta_prime as grandmother
            if ( (TMath::Abs(dummyPartGrandMotherPDG) == 111 || TMath::Abs(dummyPartGrandMotherPDG) == 221 || TMath::Abs(dummyPartGrandMotherPDG) == 331)){

              Bool_t helpN                  = true;
              for(Int_t j=0; j<fNCaloPhotonMotherMCLabels; j++){
                if (fCaloPhotonMotherMCLabels[j] == dummyPartGrandMotherLabel){ //check if grandmother is already contained in fCaloPhotonMotherMCLabels
                  helpN                     = false;
                }
              }
              if (helpN){
                fCaloPhotonMotherMCLabels[fNCaloPhotonMotherMCLabels]  = dummyPartGrandMotherLabel;
                fNCaloPhotonMotherMCLabels++; //only if particle label is not yet contained in array, count up fNCaloPhotonMotherMCLabels
              }
              if (!isPhoton && !isElectron)
                isSubLeadingEM              = kTRUE;
            }
          }
        }
      }
    }
  }
  fCaloPhotonMCFlags = isPhoton *1 + isElectron *2 + isConversion*4+ isConversionFullyContained *8 + isMerged *16 + isMergedPartConv*32 + isDalitz *64 + isDalitzMerged *128 + isPhotonWithElecMother *256 + isShower * 512 + isSubLeadingEM * 1024 + isElectronFromFragPhoton * 2048;
}

void AliAODConversionPhoton::SetCaloPhotonMCFlagsAOD(AliVEvent* event, Bool_t enableSort, Bool_t mergedAnalysis, AliVCluster* cluster){

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!AODMCTrackArray) return;

  AliAODMCParticle* PhotonDummyMerged;
  AliAODMCParticle* PhotonDummyMergedMother;
  AliAODMCParticle* PhotonDummyMergedGrandMother;
  AliAODMCParticle* PhotonDummyMergedGGMother;
  Int_t photonDummyMergedPDG= -1;
  Int_t photonDummyMergedMotherPDG= -1;
  Int_t photonDummyMergedGrandMotherPDG= -1;
  Int_t photonDummyMergedGGMotherPDG= -1;
  Bool_t foundNeutralPion = kFALSE;
  Int_t neutralPionLabel = -1;
  if(mergedAnalysis){
    // this part searches the MC labels of the cluster/photon candidate for all contributions from pi0s and save all pi0s that contributed
    for (Int_t j = 0; j< fNCaloPhotonMCLabels; j++){
      neutralPionLabel = -1;
      foundNeutralPion = kFALSE;
      PhotonDummyMerged        = (AliAODMCParticle*) AODMCTrackArray->At(GetCaloPhotonMCLabel(j)); // main particle
      photonDummyMergedPDG = PhotonDummyMerged->GetPdgCode();
      if(TMath::Abs(photonDummyMergedPDG)==111){
        foundNeutralPion = kTRUE;
        neutralPionLabel = GetCaloPhotonMCLabel(j);
      }
      if(PhotonDummyMerged->GetMother()>-1 && !foundNeutralPion){
        PhotonDummyMergedMother  = (AliAODMCParticle*) AODMCTrackArray->At(PhotonDummyMerged->GetMother()); // mother
        photonDummyMergedMotherPDG = PhotonDummyMergedMother->GetPdgCode();
        if(TMath::Abs(photonDummyMergedMotherPDG)==111){
          foundNeutralPion = kTRUE;
          neutralPionLabel = PhotonDummyMerged->GetMother();
        }
        if(PhotonDummyMergedMother->GetMother()>-1 && !foundNeutralPion){
          PhotonDummyMergedGrandMother  = (AliAODMCParticle*) AODMCTrackArray->At(PhotonDummyMergedMother->GetMother()); // grandmother
          photonDummyMergedGrandMotherPDG = PhotonDummyMergedGrandMother->GetPdgCode();
          if(TMath::Abs(photonDummyMergedGrandMotherPDG)==111){
            foundNeutralPion = kTRUE;
            neutralPionLabel = PhotonDummyMergedMother->GetMother();
          }
         if(PhotonDummyMergedGrandMother->GetMother()>-1 && !foundNeutralPion){
            PhotonDummyMergedGGMother  = (AliAODMCParticle*) AODMCTrackArray->At(PhotonDummyMergedGrandMother->GetMother()); // grand-grandmother
            photonDummyMergedGGMotherPDG = PhotonDummyMergedGGMother->GetPdgCode();
            if(TMath::Abs(photonDummyMergedGGMotherPDG)==111){
              foundNeutralPion = kTRUE;
              neutralPionLabel = PhotonDummyMergedGrandMother->GetMother();
            }
          }
        }
      }
      // check if current found pi0 was already found before in the cluster
      // if it was found before, add the cluster energy fraction of the current label to this pi0
      if(foundNeutralPion){
        Bool_t newNeutralPion                    = true;
        for(Int_t k=0; k<fNNeutralPionLabels; k++){
          if (fNeutralPionLabels[k] == neutralPionLabel){ //check if mother is already contained in fNeutralPionLabels
            newNeutralPion                       = false;
            fNeutralPionEnergyFraction[k] += cluster->GetClusterMCEdepFraction(j);
          }
        }
        // if the pi0 is new, then increase pi0 count by one and set E-fraction and label to array
        if (newNeutralPion){
          fNeutralPionLabels[fNNeutralPionLabels]    = neutralPionLabel;
          fNeutralPionEnergyFraction[fNNeutralPionLabels] = cluster->GetClusterMCEdepFraction(j);
          fNNeutralPionLabels++; //only if particle label is not yet contained in array, count up fNeutralPionLabels
        }
      }
    }

    // search for the largest contributing pi0 in the cluster
    Double_t largestEfraction = 0;
    for(Int_t k=0; k<fNNeutralPionLabels; k++){
      if(fNeutralPionEnergyFraction[k]>largestEfraction){
        largestEfraction = fNeutralPionEnergyFraction[k];
        fLeadingNeutralPionIndex = k;
      }
    }
    // if the cluster contains a pi0 contribution, find the MC label/index of the largest contributing daughter
    // for this, go upwards three mothers from each MC label to find leading daughter index
    if(fLeadingNeutralPionIndex>-1){
      for(Int_t l=0; l<fNCaloPhotonMCLabels; l++){
        if(fLeadingNeutralPionDaughterIndex>-1)
          continue;
        PhotonDummyMerged        = (AliAODMCParticle*) AODMCTrackArray->At(GetCaloPhotonMCLabel(l)); // main particle
        if(PhotonDummyMerged->GetMother()>-1){
          if(PhotonDummyMerged->GetMother()==fNeutralPionLabels[fLeadingNeutralPionIndex]){
            fLeadingNeutralPionDaughterIndex = l;
          }else{
            PhotonDummyMergedMother        = (AliAODMCParticle*) AODMCTrackArray->At(PhotonDummyMerged->GetMother());
            if(PhotonDummyMergedMother->GetMother()>-1){
              if(PhotonDummyMergedMother->GetMother()==fNeutralPionLabels[fLeadingNeutralPionIndex]){
                fLeadingNeutralPionDaughterIndex = l;
              }else{
                PhotonDummyMergedGrandMother        = (AliAODMCParticle*) AODMCTrackArray->At(PhotonDummyMergedMother->GetMother());
                if(PhotonDummyMergedGrandMother->GetMother()>-1){
                  if(PhotonDummyMergedGrandMother->GetMother()==fNeutralPionLabels[fLeadingNeutralPionIndex]){
                    fLeadingNeutralPionDaughterIndex = l;
                  }else{
                    PhotonDummyMergedGGMother        = (AliAODMCParticle*) AODMCTrackArray->At(PhotonDummyMergedGrandMother->GetMother());
                    if(PhotonDummyMergedGGMother->GetMother()>-1){
                      if(PhotonDummyMergedGGMother->GetMother()==fNeutralPionLabels[fLeadingNeutralPionIndex]){
                        fLeadingNeutralPionDaughterIndex = l;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  Bool_t isPhoton                   = kFALSE; // largest contribution to cluster is photon
  Bool_t isElectron                 = kFALSE; // largest contribution to cluster is electron
  Bool_t isConversion               = kFALSE; // largest contribution to cluster is converted electron
  Bool_t isConversionFullyContained = kFALSE; // largest contribution to cluster is converted electron, second electron has been found in same cluster
  Bool_t isMerged                   = kFALSE; // largest contribution to cluster is photon, second photon or electron from dalitz has been found in same cluster
  Bool_t isMergedPartConv           = kFALSE; // cluster contains more than one particle from the same decay and at least one of the particles came from a conversion
  Bool_t isDalitz                   = kFALSE; // this cluster was created by a particle stemming from a dality decay
  Bool_t isDalitzMerged             = kFALSE; // this cluster was created by a particle stemming from a dality decay and more than one particle of the dalitz decay is contained in the cluster
  Bool_t isPhotonWithElecMother     = kFALSE; // this cluster is from a photon with an electron as mother
  Bool_t isShower                   = kFALSE; // this cluster contains as a largest contribution a particle from a shower or radiative process
  Bool_t isSubLeadingEM             = kFALSE; // cluster contains at least one electron or photon from a pi0 or eta in subleading contribution
  Bool_t isElectronFromFragPhoton   = kFALSE; // largest contribution to cluster is from converted electron, but photon stems from fragmentation photon ( q -> q gamma)

  AliAODMCParticle* Photon = NULL;
  AliAODMCParticle* PhotonMother = NULL;
  AliAODMCParticle* PhotonGrandMother = NULL;

  if (fNCaloPhotonMCLabels==0) return;

  if (enableSort){
    // sort the array according to the energy of contributing particles
    if (fNCaloPhotonMCLabels>1){
      // cout << "start sorting" << endl;
      Int_t* sortIdx            = new Int_t[fNCaloPhotonMCLabels];
      Double_t* energyPerPart   = new Double_t[fNCaloPhotonMCLabels];
      Long_t* orginalContrib    = new Long_t[fNCaloPhotonMCLabels];
      for (Int_t i = 0; i < fNCaloPhotonMCLabels; i++){
        orginalContrib[i]       = fCaloPhotonMCLabels[i];
        if (fCaloPhotonMCLabels[i]> -1){
          AliAODMCParticle* dummy  = (AliAODMCParticle*) AODMCTrackArray->At(GetCaloPhotonMCLabel(i));
          energyPerPart[i]  = dummy->E();
          // suppress energy of hadrons !!! DIRTY hack !!!
          if (!(TMath::Abs(dummy->GetPdgCode())== 11 || TMath::Abs(dummy->GetPdgCode())== 22)){
            // cout << "suppressed hadron energy" << endl;
            energyPerPart[i]= 0.25; // put energy to mip
          }
        } else {
          energyPerPart[i]  = 0;
        }
      }

      TMath::Sort(fNCaloPhotonMCLabels,energyPerPart,sortIdx);
      for(Int_t index = 0; index < fNCaloPhotonMCLabels; index++) {
        fCaloPhotonMCLabels[index] = orginalContrib  [sortIdx[index]] ;

      }
      delete [] sortIdx;
      delete [] energyPerPart;
      delete [] orginalContrib;
    }
  }

  if(mergedAnalysis && fNNeutralPionLabels>0 && fLeadingNeutralPionDaughterIndex!=0 && fNeutralPionEnergyFraction[fLeadingNeutralPionIndex]>cluster->GetClusterMCEdepFraction(0)){
    // check if leading pi0 comes not from label 0 in cluster
    // for this do:
    // -> check if neutral pions were found in cluster
    // -> if the leading daughter index is not 0
    // -> the leading neutral pion has a larger cluster energy fraction than the cluster label 0
    // load particle corresponding to largest daughter of leading pi0
    Photon                            = (AliAODMCParticle*) AODMCTrackArray->At(GetCaloPhotonMCLabel(fLeadingNeutralPionDaughterIndex));
  } else {
    // load particle corresponding to MC label 0 in cluster
    if(GetCaloPhotonMCLabel(0)>-1){
      Photon                            = (AliAODMCParticle*) AODMCTrackArray->At(GetCaloPhotonMCLabel(0));
    }
  }



  if(Photon == NULL){
    return;
  }

  Int_t particleMotherLabel           = Photon->GetMother();
  Int_t particleGrandMotherLabel      = -1;
  Int_t particleMotherPDG             = -1;
  Int_t particleGrandMotherPDG        = -1;
  Int_t particleMotherNDaugthers      = 0;
  Int_t particleGrandMotherNDaugthers = 0;
  Int_t particleGrandMotherX2Label    = -1;
  Int_t particleGrandMotherX3Label    = -1;
  Int_t particleGrandMotherX4Label    = -1;
  Int_t particleGrandMotherX5Label    = -1;
  Int_t particleGrandMotherX6Label    = -1;
  Int_t particleGrandMotherX7Label    = -1;
  Int_t particleGrandMotherX2PDG      = -1;
  Int_t particleGrandMotherX3PDG      = -1;
  Int_t particleGrandMotherX4PDG      = -1;
  Int_t particleGrandMotherX5PDG      = -1;
  Int_t particleGrandMotherX6PDG      = -1;

  if (particleMotherLabel > -1){
    PhotonMother                      = (AliAODMCParticle*) AODMCTrackArray->At(Photon->GetMother());
    particleMotherNDaugthers          = PhotonMother->GetNDaughters();
    particleGrandMotherLabel          = PhotonMother->GetMother();
    particleMotherPDG                 = PhotonMother->GetPdgCode();
    if (particleGrandMotherLabel > -1){
      PhotonGrandMother               = (AliAODMCParticle*) AODMCTrackArray->At(PhotonMother->GetMother());
      particleGrandMotherPDG          = PhotonGrandMother->GetPdgCode();
      particleGrandMotherNDaugthers   = PhotonGrandMother->GetNDaughters();
      AliAODMCParticle* dummyGMM      = NULL;

      particleGrandMotherX2Label      = PhotonGrandMother->GetMother();
      if (particleGrandMotherX2Label > -1){
        dummyGMM                        = (AliAODMCParticle*) AODMCTrackArray->At(particleGrandMotherX2Label);
        particleGrandMotherX2PDG        = dummyGMM->GetPdgCode();
        particleGrandMotherX3Label      = dummyGMM->GetMother();
        if (particleGrandMotherX3Label > -1){
          dummyGMM                        = (AliAODMCParticle*) AODMCTrackArray->At(particleGrandMotherX3Label);
          particleGrandMotherX3PDG        = dummyGMM->GetPdgCode();
          particleGrandMotherX4Label      = dummyGMM->GetMother();
          if (particleGrandMotherX4Label > -1){
            dummyGMM                        = (AliAODMCParticle*) AODMCTrackArray->At(particleGrandMotherX4Label);
            particleGrandMotherX4PDG        = dummyGMM->GetPdgCode();
            particleGrandMotherX5Label      = dummyGMM->GetMother();
            if (particleGrandMotherX5Label > -1){
              dummyGMM                        = (AliAODMCParticle*) AODMCTrackArray->At(particleGrandMotherX5Label);
              particleGrandMotherX5PDG        = dummyGMM->GetPdgCode();
              particleGrandMotherX6Label      = dummyGMM->GetMother();
              if (particleGrandMotherX6Label > -1){
                dummyGMM                        = (AliAODMCParticle*) AODMCTrackArray->At(particleGrandMotherX6Label);
                particleGrandMotherX6PDG        = dummyGMM->GetPdgCode();
                particleGrandMotherX7Label      = dummyGMM->GetMother();
              }
            }
          }
        }
      }
    }
  }

  //determine mother/grandmother of leading particle and if it is pion/eta/eta_prime: fill array fCaloPhotonMotherMCLabels at position 0
  if (particleMotherLabel > -1){
    if( TMath::Abs(particleMotherPDG) == 111 || TMath::Abs(particleMotherPDG) == 221 || TMath::Abs(particleMotherPDG) == 331 ){
      fCaloPhotonMotherMCLabels[0]    = particleMotherLabel;
      fNCaloPhotonMotherMCLabels++;
    } else if (TMath::Abs(particleMotherPDG) == 22 && particleGrandMotherLabel > -1){
      if ( TMath::Abs(particleMotherPDG) == 22 && (TMath::Abs(particleGrandMotherPDG) == 111 || TMath::Abs(particleGrandMotherPDG) == 221 || TMath::Abs(particleGrandMotherPDG) == 331) ){
        fCaloPhotonMotherMCLabels[0]  = particleGrandMotherLabel;
        fNCaloPhotonMotherMCLabels++;
      }
    } else {
      fCaloPhotonMotherMCLabels[0]    = particleMotherLabel;
      fNCaloPhotonMotherMCLabels++;
    }
  }

  // Check whether the first contribution was photon
  if(TMath::Abs(Photon->GetPdgCode()) == 22){
    isPhoton                          = kTRUE;
    // did it decay via the dalitz channel
    if (particleMotherLabel > -1 && particleMotherNDaugthers == 3)
      isDalitz                        = kTRUE;
    // Test whether particle stems from a shower or radiation
    if (TMath::Abs(particleMotherPDG) == 11){            // check whether photon stems from electron
      isPhotonWithElecMother          = kTRUE;
      if (particleGrandMotherLabel > -1){         // test whether first particle has a grandmother
        if (TMath::Abs(particleGrandMotherPDG) == 22 )
          isShower                    = kTRUE;    // check whether grandmother is a photon (meaning this is most likely a shower)
      }
    }
  }
  // Check whether the first contribution was electron
  if(TMath::Abs(Photon->GetPdgCode()) == 11 ){
    isElectron                        = kTRUE;
    if (particleMotherLabel > -1) {
      // was it a conversion
      if (TMath::Abs(particleMotherPDG) == 22)
        isConversion                  = kTRUE;
      // did it decay via the dalitz channel
      if (particleGrandMotherLabel > -1 && particleGrandMotherNDaugthers == 3 )
        isDalitz                      = kTRUE;
    }
    if (particleGrandMotherLabel > -1){           // check whether electron has a grandmother
      if (TMath::Abs(particleGrandMotherPDG) == 11 ||  TMath::Abs(particleGrandMotherPDG) == 22){	// test whether electron has photon or electron as grandmother (meaning will most likely be a shower)
        isShower                      = kTRUE;
      }
    }
    // consider the rare case, where the conversion electron stems from photon, which stems from electron, which stems from photon, ...
    if (isConversion && TMath::Abs(particleGrandMotherPDG) == 11   && particleGrandMotherX2PDG == 22){
      SetCaloPhotonMCLabel(0,particleGrandMotherLabel);
      fCaloPhotonMotherMCLabels[0] = particleGrandMotherX3Label;
      if (TMath::Abs(particleGrandMotherX3PDG) == 11 && particleGrandMotherX4PDG == 22  ){
        SetCaloPhotonMCLabel(0,particleGrandMotherX3Label);
        fCaloPhotonMotherMCLabels[0] = particleGrandMotherX5Label;
        if (TMath::Abs(particleGrandMotherX5PDG) == 11 && particleGrandMotherX6PDG == 22  ){
          SetCaloPhotonMCLabel(0,particleGrandMotherX5Label);
          fCaloPhotonMotherMCLabels[0] = particleGrandMotherX7Label;
        }
      }
    }

    // consider the case, where a photon stems from a photon which stems from a photon etc...which is common for frag. photons
    if (particleMotherLabel > -1){
      AliAODMCParticle* dummyMother = (AliAODMCParticle*) AODMCTrackArray->At(Photon->GetMother());
      if (dummyMother){
        Bool_t originReached   = kFALSE;
        while (dummyMother->GetPdgCode() == 22 && !originReached){ // follow conversion photon's history, as long as the mother is a photon
          if (dummyMother->GetMother() > -1){
            dummyMother = (AliAODMCParticle*) AODMCTrackArray->At(dummyMother->GetMother());
            if ((TMath::Abs(dummyMother->GetPdgCode()) == 11) || (TMath::Abs(dummyMother->GetPdgCode()) == 22)){  // in case of additional conversion skip to photon's grandma, which should be a photon
              if (dummyMother->GetMother() > -1){
                dummyMother = (AliAODMCParticle*) AODMCTrackArray->At(dummyMother->GetMother());
              } else {
                  originReached = kTRUE;
              }
            } else {
              originReached = kTRUE;
            }
            isElectronFromFragPhoton = (TMath::Abs(dummyMother->GetPdgCode()) < 6);// photon stems from quark = fragmentation photon
          } else {
            originReached = kTRUE;
          }
        }
      }
    }
  }



  // check whether there were other contributions to the cluster
  if (fNCaloPhotonMCLabels>1){
    AliAODMCParticle* dummyPart             = NULL;
    AliAODMCParticle* dummyPartMother       = NULL;
    AliAODMCParticle* dummyPartGrandMother  = NULL;
    for (Int_t i = 1; i< fNCaloPhotonMCLabels; i++){
      if (i > 49) continue;                       // abort if more than 50 entries to the cluster have been checked (more are not stored in these objects)
      dummyPart = (AliAODMCParticle*) AODMCTrackArray->At(GetCaloPhotonMCLabel(i));
      Int_t dummyPartMotherLabel            = dummyPart->GetMother();
      Int_t dummyPartGrandMotherLabel       = -1;
      Int_t dummyPartMotherPDG              = -1;
      Int_t dummyPartGrandMotherPDG         = -1;

      // check whether this particle has a mother & obtain the pdg code
      if (dummyPartMotherLabel > -1){
        dummyPartMother                     = (AliAODMCParticle*) AODMCTrackArray->At(dummyPart->GetMother());
        dummyPartGrandMotherLabel           = dummyPartMother->GetMother();
        dummyPartMotherPDG                  = dummyPartMother->GetPdgCode();
        // check whether this particle has a grandmother & obtain its pdg code
        if (dummyPartGrandMotherLabel > -1){
          dummyPartGrandMother              = (AliAODMCParticle*) AODMCTrackArray->At(dummyPartMother->GetMother());
          dummyPartGrandMotherPDG           = dummyPartGrandMother->GetPdgCode();
        }
      }


            // largest contribution was from photon and is not from shower or electron mother
      if (isPhoton && (!isShower || !isPhotonWithElecMother )){
        if (dummyPartMotherLabel > -1 && particleMotherLabel > -1){   // test whether first particle has a mother
          if (dummyPartMotherLabel == particleMotherLabel)
            isMerged                  = kTRUE;                        // test whether current and first particle have the same mother => i.e. other gamma from decay or dalitz electron
          if (dummyPartGrandMotherLabel > -1){                        // test whether first particle has a grandmother
            // check whether particle is an electron from a conversion of a photon from the original mother
            if (TMath::Abs(dummyPart->GetPdgCode()) == 11 && dummyPartGrandMotherLabel == particleMotherLabel )
              isMergedPartConv        = kTRUE;
            // check whether particle is an electron from a dalitz decay from the original mother
            if (TMath::Abs(dummyPart->GetPdgCode()) == 11 && dummyPartMotherLabel == particleMotherLabel )
              isDalitzMerged          = kTRUE;
          }
        }
      }

      // largest contribution was from electron & not a from a shower
      if (isElectron && !isShower){
        if (dummyPartMotherLabel > -1 && particleMotherLabel > -1){     // test whether first particle has a mother
          if ( TMath::Abs(dummyPart->GetPdgCode()) == 11 && isConversion && dummyPartMotherLabel == particleMotherLabel  ) {
            isConversionFullyContained  = kTRUE;                        // test whether conversion is fully contained in cluster
          } else if (TMath::Abs(dummyPart->GetPdgCode()) == 11) { // current particle is an electron

            // check whether current particle is a shower
            Bool_t isShowerCurrent    = kFALSE;
            if (TMath::Abs(dummyPartGrandMotherPDG) == 11)
              isShowerCurrent         = kTRUE;
            if (TMath::Abs(dummyPartGrandMotherPDG) == 22)
              isShowerCurrent         = kTRUE;

            // check whether current particle is a conversion
            Bool_t isConvCurrent      = kFALSE;
            if (TMath::Abs(dummyPartMotherPDG) == 22 && !isShowerCurrent)
              isConvCurrent           = kTRUE;

            // check whether orginal electron and this electron stem from the same particle and electron originated in dalitz decay
            if( dummyPartMotherLabel == particleMotherLabel && TMath::Abs(particleMotherPDG) != 22 &&  !isShowerCurrent   ) {
              isDalitzMerged          = kTRUE;
              isMerged                = kTRUE;
            // both particles are from conversions
            } else if (isConversion && isConvCurrent){
              if (particleGrandMotherLabel > -1 && dummyPartGrandMotherLabel > -1){ // test whether first particle has a grandmother and current particle has one as well
                // check whether orginal electron and this electron stem from the same particle and electron stems from conversion
                if( dummyPartGrandMotherLabel == particleGrandMotherLabel &&  !isShowerCurrent   )
                  isMergedPartConv        = kTRUE;
              }
            // check whether current electron is from conversion & orginal e are from same mother
            } else if (dummyPartGrandMotherLabel > -1 && isConvCurrent){
              if (dummyPartGrandMotherLabel == particleMotherLabel && (TMath::Abs(particleMotherPDG) == 111 || TMath::Abs(particleMotherPDG) == 221 || TMath::Abs(particleMotherPDG) == 331) ){
                isDalitzMerged          = kTRUE;
                isMerged                = kTRUE;
              }
            // check whether current electron & orginal e is from conversion are from same mother
            } else if (isConversion && particleGrandMotherLabel > -1){
              if (dummyPartMotherLabel == particleGrandMotherLabel && (TMath::Abs(dummyPartMotherPDG) == 111 || TMath::Abs(dummyPartMotherPDG) == 221 || TMath::Abs(dummyPartMotherPDG) == 331) ){
                isDalitzMerged          = kTRUE;
                isMerged                = kTRUE;
              }
            }
          } else if (TMath::Abs(dummyPart->GetPdgCode()) == 22){                 // test whether this particle is a photon
            // check whether current particle is a shower
            Bool_t isShowerCurrent    = kFALSE;
            if (TMath::Abs(dummyPartMotherPDG) == 11)
              isShowerCurrent         = kTRUE;
            if (TMath::Abs(dummyPartMotherPDG) == 22)
              isShowerCurrent         = kTRUE;

            // check whether orginal electron and this photon stem from the same particle and electron originated in dalitz
            if( dummyPartMotherLabel == particleMotherLabel && !isShowerCurrent ) {
              isDalitzMerged          = kTRUE;
              isMerged                = kTRUE;
            } else if (particleGrandMotherLabel > -1){ // test whether first particle has a grandmother
              // check whether orginal electron and this photon stem from the same particle and electron stems from conversion
              if( dummyPartMotherLabel == particleGrandMotherLabel && !isShowerCurrent )
                isMergedPartConv        = kTRUE;
            }
          }
        }
      }

      if (dummyPartMotherLabel > -1){ // test whether particle has a mother
        if (TMath::Abs(dummyPart->GetPdgCode()) == 22){ // test whether particle is a photon
          //check if photon directly comes from a pion/eta/eta_prime decay
          Bool_t helpN                    = true;
          for(Int_t j=0; j<fNCaloPhotonMotherMCLabels; j++){
            if (fCaloPhotonMotherMCLabels[j] == dummyPartMotherLabel){ //check if mother is already contained in fCaloPhotonMotherMCLabels
              helpN                       = false;
            }
          }
          if (helpN){
            fCaloPhotonMotherMCLabels[fNCaloPhotonMotherMCLabels]    = dummyPartMotherLabel;
            fNCaloPhotonMotherMCLabels++; //only if particle label is not yet contained in array, count up fNCaloPhotonMotherMCLabels
          }
          if ((TMath::Abs(dummyPartMotherPDG) == 111 || TMath::Abs(dummyPartMotherPDG) == 221 || TMath::Abs(dummyPartMotherPDG) == 331) && !isPhoton && !isElectron)
            isSubLeadingEM                = kTRUE;
        } else if (TMath::Abs(dummyPart->GetPdgCode()) == 11){ //test whether particle is an electron
          //check if electron comes from a pion decay
          if ( TMath::Abs(dummyPartMotherPDG) != 22 ){
            Bool_t helpN                    = true;
            for(Int_t j=0; j<fNCaloPhotonMotherMCLabels; j++){
              if (fCaloPhotonMotherMCLabels[j] == dummyPartMotherLabel){ //check if mother is already contained in fCaloPhotonMotherMCLabels
                helpN                       = false;
              }
            }
            if (helpN){
              fCaloPhotonMotherMCLabels[fNCaloPhotonMotherMCLabels]    = dummyPartMotherLabel;
              fNCaloPhotonMotherMCLabels++; //only if particle label is not yet contained in array, count up fNCaloPhotonMotherMCLabels
            }
            if ((TMath::Abs(dummyPartMotherPDG) == 111 || TMath::Abs(dummyPartMotherPDG) == 221 || TMath::Abs(dummyPartMotherPDG) == 331) && !isPhoton && !isElectron)
              isSubLeadingEM                = kTRUE;
          } else if (dummyPartGrandMotherLabel > -1){ //if it is not a dalitz decay, test whether particle has a grandmother
            //check if it is a conversion electron that has pion/eta/eta_prime as grandmother
            if ( (TMath::Abs(dummyPartGrandMotherPDG) == 111 || TMath::Abs(dummyPartGrandMotherPDG) == 221 || TMath::Abs(dummyPartGrandMotherPDG) == 331)){

              Bool_t helpN                  = true;
              for(Int_t j=0; j<fNCaloPhotonMotherMCLabels; j++){
                if (fCaloPhotonMotherMCLabels[j] == dummyPartGrandMotherLabel){ //check if grandmother is already contained in fCaloPhotonMotherMCLabels
                  helpN                     = false;
                }
              }
              if (helpN){
                fCaloPhotonMotherMCLabels[fNCaloPhotonMotherMCLabels]  = dummyPartGrandMotherLabel;
                fNCaloPhotonMotherMCLabels++; //only if particle label is not yet contained in array, count up fNCaloPhotonMotherMCLabels
              }
              if (!isPhoton && !isElectron)
                isSubLeadingEM              = kTRUE;
            }
          }
        }
      }
    }
  }
  fCaloPhotonMCFlags = isPhoton *1 + isElectron *2 + isConversion*4+ isConversionFullyContained *8 + isMerged *16 + isMergedPartConv*32 + isDalitz *64 + isDalitzMerged *128 + isPhotonWithElecMother *256 + isShower * 512 + isSubLeadingEM * 1024 + isElectronFromFragPhoton * 2048;;
}

//_____________________________________________________________________________
// prints information of MC particles contributing to cluster
//_____________________________________________________________________________
void AliAODConversionPhoton::PrintCaloMCLabelsAndInfo(AliMCEvent *mcEvent){
  cout << endl << endl << "particles contributing: " << endl;
  for (Int_t i =0 ; i < fNCaloPhotonMCLabels; i++ ){
    TParticle *dummy    = NULL;
    if (fCaloPhotonMCLabels[i]>0){
      dummy             = (TParticle*)mcEvent->Particle(fCaloPhotonMCLabels[i]);
      cout << i << "\t"<< fCaloPhotonMCLabels[i] << "\t pdg code: " <<dummy->GetPdgCode() << "\t prod radius: "<< dummy->R() << "\t energy: " << dummy->Energy() ;
      if (dummy->GetMother(0) > -1){
          TParticle* dummyMother   = (TParticle*)mcEvent->Particle(dummy->GetMother(0));
          cout << "\t mother part: " << dummy->GetMother(0) << "\t mother pdg code: " << dummyMother->GetPdgCode() << "\t energy: " << dummyMother->Energy() << endl;
          if (dummyMother->GetMother(0) > -1){
            TParticle* dummyGrandMother   = (TParticle*)mcEvent->Particle(dummyMother->GetMother(0));
            cout << "\t grandmother part: " << dummyMother->GetMother(0) << "\t grandmother pdg code: " << dummyGrandMother->GetPdgCode() << "\t energy: " << dummyGrandMother->Energy() << endl;
          } else {
              cout << endl;
          }
      } else {
        cout << endl;
      }
    }
//     if (dummy) delete dummy;
  }
  cout << "mothers contributing: " << endl;
  for (Int_t i =0 ; i < fNCaloPhotonMotherMCLabels; i++ ){
    TParticle *dummy    = NULL;
    if (fCaloPhotonMotherMCLabels[i]>0){
      dummy             = (TParticle*)mcEvent->Particle(fCaloPhotonMotherMCLabels[i]);
      cout << i << "\t"<< fCaloPhotonMotherMCLabels[i] << "\t pdg code: " <<dummy->GetPdgCode() << "\t prod radius: "<< dummy->R() << "\t energy: " << dummy->Energy() << endl;


    }
//     if (dummy) delete dummy;
  }
}

//_____________________________________________________________________________
// prints information of cluster as it is flagged
//_____________________________________________________________________________
void AliAODConversionPhoton::PrintCaloMCFlags(){
  cout << fCaloPhotonMCFlags << "\t photon: " << IsLargestComponentPhoton() << "\t electron: " << IsLargestComponentElectron() << "\t conversion: " << IsConversion() << "\t conversion full: "
       << IsConversionFullyContained() << "\t merged: " << IsMerged() << "\t merged p.conv: " << IsMergedPartConv() << "\t Dalitz: " << IsDalitz() << "\t Dalitz merged: " << IsDalitzMerged()
       << "\t rad: " << IsPhotonWithElecMother() << "\t shower: " << IsShower() << "\t no EM lead: " << IsEMNonLeading() << "\t EM sub lead: " << IsSubLeadingEM() << endl;
}
