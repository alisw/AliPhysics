//////////////////////////////////////////////////////////////////////
//                                                                  //
// The macro allows to make all needed cuts on V0 particles and     //
// to remove from the V0Collection particles which share            //
// at least one of their daughters with accepted ones using lower   //
// DCA to primary vertex of V0 as the criterion.                    //
//                                                                  //
// Author: Dominik Arominski                                        //
// Email: dominik.arominski@cern.ch                                 //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include "AliFemtoV0SharedDaughterCut.h"

AliFemtoV0SharedDaughterCut::AliFemtoV0SharedDaughterCut() { }
AliFemtoV0SharedDaughterCut::~AliFemtoV0SharedDaughterCut() { }

AliFemtoV0Collection AliFemtoV0SharedDaughterCut::AliFemtoV0SharedDaughterCutCollection(AliFemtoV0Collection *V0Collection,
                                                                                        AliFemtoV0Cut *pCut) {

  Int_t collectionSize = V0Collection->size();         //determination of collection size
  Int_t positionInCollection=0;                        //used to iterate inside the collection
  Int_t count_pass=0;                                  //used to indicate how many V0s are accepted at the moment
  Int_t *IdPosArray = new Int_t[collectionSize]; //used to determine if there is more than one pretendent for one daughter ID
  Int_t *IdNegArray = new Int_t[collectionSize]; //same as above
  Double_t *dcaToPrimVertex = new Double_t[collectionSize]; //used to discriminate between pretendents using DCA to primary vertex
  bool acceptV0 = false;                                            //used to determine if particle should be pushed to the collection
  AliFemtoV0Collection V0CorrectedCollection;               //collection to be returned

  AliFemtoV0* pParticle;
  AliFemtoV0Iterator iterator;
  AliFemtoV0Iterator start = V0Collection->begin();
  AliFemtoV0Iterator end = V0Collection->end();

  for (iterator=start;iterator!=end;iterator++){
    pParticle = *iterator;
    bool tmpPassV0 = pCut->Pass(pParticle);
    // it is not certain if particle which passed all the cuts will stay in the collection but if it fails - it fails for sure
    pCut->FillCutMonitor(pParticle,tmpPassV0);
    if(tmpPassV0) {
      acceptV0=true;                                                //it stays true if there are no conflicts
      IdPosArray[count_pass] = pParticle->IdPos();                  //id array of positive V0's daughters
      IdNegArray[count_pass] = pParticle->IdNeg();                  //id array of negative V0's daughters
      dcaToPrimVertex[count_pass] = pParticle->DcaV0ToPrimVertex(); //array of corresponding values of V0's DCA to primary vertex

      for(Int_t ii=0; ii<count_pass; ii++) { //used to check if in collection are already V0s with these daughters' IDs
        if( (IdPosArray[count_pass]==IdPosArray[ii]) || (IdNegArray[count_pass]==IdNegArray[ii]) ) {
          //checking if any of the new V0 daughters doesn't have the same ID value as particles already in the set
          if(dcaToPrimVertex[count_pass] < dcaToPrimVertex[ii]) { //checking if new V0 is better than the old one
            //if true there is a need to remove worse particle from the collection
            for(AliFemtoV0Iterator iter=V0CorrectedCollection.begin(); iter!=V0CorrectedCollection.end(); iter++ ) {
              if(positionInCollection==ii) {
                V0CorrectedCollection.insert(iter, 1, pParticle);  //inserting new better V0
                V0CorrectedCollection.erase(iter);                 //removing old V0
                IdPosArray[ii] = IdPosArray[count_pass];           //update of positive daughter's ID
                IdNegArray[ii] = IdNegArray[count_pass];           //update of negative daughter's ID
                dcaToPrimVertex[ii] = dcaToPrimVertex[count_pass]; //update of V0's DCA to primary vertex
                break;
              }
              positionInCollection++;
            }
            positionInCollection=0;
          }
          acceptV0=false; //conflict found - we don't need to change anything(worse DCA) or change has been already done
          break;
        }
      }
      if(acceptV0) {                                //if there are no problems with the particle it is added to the collection
        V0CorrectedCollection.push_back(pParticle); //finally particle is added to collection
        count_pass++;                               //particle accepted so number of particles increases
      }
    }
  }


  delete [] IdPosArray;
  delete [] IdNegArray;
  delete [] dcaToPrimVertex;

  return V0CorrectedCollection;
}
