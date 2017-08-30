#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODConversionPhoton.h"
#include "AliKFConversionPhoton.h"

using namespace std;

ClassImp(AliAODConversionPhoton)

AliAODConversionPhoton::AliAODConversionPhoton() :
AliAODConversionParticle(),
AliConversionPhotonBase(),
fDCArPrimVtx(0),
fDCAzPrimVtx(0),
fInvMassPair(0),
fCaloPhoton(0),
fCaloClusterRef(-1),
fNCaloPhotonMCLabels(0),
fNCaloPhotonMotherMCLabels(0),
fCaloPhotonMCFlags(0)
{
  // initialize calo photon MC labels
  for (Int_t i =0; i<50; i++){
    fCaloPhotonMCLabels[i]=-1;		
  }
  for (Int_t i =0; i<20; i++){
    fCaloPhotonMotherMCLabels[i]=-1;
  }
  //Standard constructor
}

AliAODConversionPhoton::AliAODConversionPhoton(AliKFConversionPhoton *kfphoton) :
AliAODConversionParticle(kfphoton),
AliConversionPhotonBase(*((AliConversionPhotonBase*)kfphoton)),
fDCArPrimVtx(0),
fDCAzPrimVtx(0),
fInvMassPair(0),
fCaloPhoton(0),
fCaloClusterRef(-1),
fNCaloPhotonMCLabels(0),
fNCaloPhotonMotherMCLabels(0),
fCaloPhotonMCFlags(0)
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
  }

}

AliAODConversionPhoton::AliAODConversionPhoton(TLorentzVector *vec) :
AliAODConversionParticle(vec),
AliConversionPhotonBase(),
fDCArPrimVtx(0),
fDCAzPrimVtx(0),
fInvMassPair(0),
fCaloPhoton(0),
fCaloClusterRef(-1),
fNCaloPhotonMCLabels(0),
fNCaloPhotonMotherMCLabels(0),
fCaloPhotonMCFlags(0)
{
  //Constructor from TLorentzVector

  // initialize calo photon MC labels
  for (Int_t i =0; i<50; i++){
    fCaloPhotonMCLabels[i]=-1;		
  }
  for (Int_t i =0; i<20; i++){
    fCaloPhotonMotherMCLabels[i]=-1;
  }
}



AliAODConversionPhoton::AliAODConversionPhoton(const AliAODConversionPhoton & original) :
AliAODConversionParticle(original),
AliConversionPhotonBase(original),
fDCArPrimVtx(original.fDCArPrimVtx),
fDCAzPrimVtx(original.fDCAzPrimVtx),
fInvMassPair(original.fInvMassPair),
fCaloPhoton(original.fCaloPhoton),
fCaloClusterRef(original.fCaloClusterRef),
fNCaloPhotonMCLabels(original.fNCaloPhotonMCLabels),
fNCaloPhotonMotherMCLabels(original.fNCaloPhotonMotherMCLabels),
fCaloPhotonMCFlags(original.fCaloPhotonMCFlags)
{
  //Copy constructor

  // initialize calo photon MC labels
  for (Int_t i =0; i<50; i++){
    fCaloPhotonMCLabels[i]=original.fCaloPhotonMCLabels[i];		
  }
  for (Int_t i =0; i<20; i++){
      fCaloPhotonMotherMCLabels[i]=original.fCaloPhotonMotherMCLabels[i];
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


void AliAODConversionPhoton::SetCaloPhotonMCFlags(AliMCEvent *mcEvent, Bool_t enableSort){
  
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
  
  
  TParticle* Photon = 0x0;
  if (fNCaloPhotonMCLabels==0) return;
  
  if (enableSort){
    // sort the array according to the energy of contributing particles
    if (fNCaloPhotonMCLabels>1){
  //     cout << "start sorting" << endl;
      Int_t* sortIdx            = new Int_t[fNCaloPhotonMCLabels];
      Double_t* energyPerPart   = new Double_t[fNCaloPhotonMCLabels];
      Long_t* orginalContrib    = new Long_t[fNCaloPhotonMCLabels];
      for (Int_t i = 0; i < fNCaloPhotonMCLabels; i++){
        orginalContrib[i]       = fCaloPhotonMCLabels[i];
        if (fCaloPhotonMCLabels[i]> -1){
          TParticle* dummy  = mcEvent->Particle(fCaloPhotonMCLabels[i]);
          energyPerPart[i]  = dummy->Energy();
          // suppress energy of hadrons !!! DIRTY hack !!!
          if (!(TMath::Abs(dummy->GetPdgCode())== 11 || TMath::Abs(dummy->GetPdgCode())== 22)){
  //           cout << "suppressed hadron energy for:" << dummy->GetPdgCode() << endl;
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
  
  if(GetCaloPhotonMCLabel(0)>-1) Photon = mcEvent->Particle(GetCaloPhotonMCLabel(0));
  
  if(Photon == NULL){
    return;
  }

  Int_t particleMotherLabel           = Photon->GetMother(0);
  Int_t particleGrandMotherLabel      = -1; 
  Int_t particleMotherPDG             = -1; 
  Int_t particleGrandMotherPDG        = -1; 
  Int_t particleMotherNDaugthers      = 0;
  Int_t particleGrandMotherNDaugthers = 0;
  if (particleMotherLabel > -1){
    particleMotherNDaugthers          = mcEvent->Particle(Photon->GetMother(0))->GetNDaughters();
    particleGrandMotherLabel          = mcEvent->Particle(Photon->GetMother(0))->GetMother(0);
    particleMotherPDG = mcEvent->Particle(Photon->GetMother(0))->GetPdgCode();
    if (particleGrandMotherLabel > -1){
      particleGrandMotherPDG          = mcEvent->Particle(mcEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode();
      particleGrandMotherNDaugthers   = mcEvent->Particle(mcEvent->Particle(Photon->GetMother(0))->GetMother(0))->GetNDaughters();
    }	
  }

  //determine mother/grandmother of leading particle and if it is pion/eta/eta_prime: fill array fCaloPhotonMotherMCLabels at position 0
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
  if(TMath::Abs(mcEvent->Particle(GetCaloPhotonMCLabel(0))->GetPdgCode()) == 22){
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
  if( TMath::Abs(mcEvent->Particle(GetCaloPhotonMCLabel(0))->GetPdgCode()) == 11 ){
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
  }

  Bool_t enablePrintOuts            = kFALSE;
//   if (fNCaloPhotonMCLabels>2)       
//     enablePrintOuts                 = kTRUE;

  // check whether there were other contributions to the cluster
  if (fNCaloPhotonMCLabels>1){
    TParticle* dummyPart =NULL;
    
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
        dummyPartGrandMotherLabel     = mcEvent->Particle(dummyPart->GetMother(0))->GetMother(0);
        dummyPartMotherPDG = mcEvent->Particle(dummyPart->GetMother(0))->GetPdgCode();
        // check whether this particle has a grandmother & obtain its pdg code
        if (dummyPartGrandMotherLabel > -1){
          dummyPartGrandMotherPDG     = mcEvent->Particle(mcEvent->Particle(dummyPart->GetMother(0))->GetMother(0))->GetPdgCode();
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
  fCaloPhotonMCFlags = isPhoton *1 + isElectron *2 + isConversion*4+ isConversionFullyContained *8 + isMerged *16 + isMergedPartConv*32 + isDalitz *64 + isDalitzMerged *128 + isPhotonWithElecMother *256 + isShower * 512 + isSubLeadingEM * 1024;
}

void AliAODConversionPhoton::SetCaloPhotonMCFlagsAOD(AliVEvent* event, Bool_t enableSort){
  
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!AODMCTrackArray) return;
  
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
  
  AliAODMCParticle* Photon;
  AliAODMCParticle* PhotonMother;
  AliAODMCParticle* PhotonGrandMother;

  if (fNCaloPhotonMCLabels==0) return;
  Photon                            = (AliAODMCParticle*) AODMCTrackArray->At(GetCaloPhotonMCLabel(0));

  if (enableSort){
    // sort the array according to the energy of contributing particles
    if (fNCaloPhotonMCLabels>1){
  //     cout << "start sorting" << endl;
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
  //           cout << "suppressed hadron energy" << endl;
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
  
  if(Photon == NULL){
    return;
  }

  Int_t particleMotherLabel           = Photon->GetMother();
  Int_t particleGrandMotherLabel      = -1; 
  Int_t particleMotherPDG             = -1; 
  Int_t particleGrandMotherPDG        = -1; 
  Int_t particleMotherNDaugthers      = 0;
  Int_t particleGrandMotherNDaugthers = 0;
  if (particleMotherLabel > -1){
    PhotonMother                      = (AliAODMCParticle*) AODMCTrackArray->At(Photon->GetMother());
    particleMotherNDaugthers          = PhotonMother->GetNDaughters();
    particleGrandMotherLabel          = PhotonMother->GetMother();
    particleMotherPDG                 = PhotonMother->GetPdgCode();
    if (particleGrandMotherLabel > -1){
      PhotonGrandMother               = (AliAODMCParticle*) AODMCTrackArray->At(PhotonMother->GetMother());
      particleGrandMotherPDG          = PhotonGrandMother->GetPdgCode();
      particleGrandMotherNDaugthers   = PhotonGrandMother->GetNDaughters();
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
  fCaloPhotonMCFlags = isPhoton *1 + isElectron *2 + isConversion*4+ isConversionFullyContained *8 + isMerged *16 + isMergedPartConv*32 + isDalitz *64 + isDalitzMerged *128 + isPhotonWithElecMother *256 + isShower * 512 + isSubLeadingEM * 1024;
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
