///
/// \file AliFemtoXiSharedDaughterCut.cxx
///



#include "AliFemtoXiSharedDaughterCut.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoXiSharedDaughterCut);
  /// \endcond
#endif

AliFemtoXiSharedDaughterCut::AliFemtoXiSharedDaughterCut() {}
AliFemtoXiSharedDaughterCut::~AliFemtoXiSharedDaughterCut() {}

AliFemtoXiCollection AliFemtoXiSharedDaughterCut::AliFemtoXiSharedDaughterCutCollection(AliFemtoXiCollection *XiCollection, AliFemtoXiTrackCut *pCut)
{
  Int_t collectionSize = XiCollection->size();	             //determination of collection size
  Int_t positionInCollection = 0;                            //used to iterate inside the collection
  Int_t count_pass=0;                                        //used to indicate how many V0s are accepted at the moment
  Int_t *IdPosArray = new Int_t[collectionSize];             //ID array of positive V0's daughters
  Int_t *IdNegArray = new Int_t[collectionSize];             //ID array of negative V0's daughters
  Int_t *IdBacArray = new Int_t[collectionSize];             //ID array of bacherlor pions
  Double_t *dcaXiToPrimVertex = new Double_t[collectionSize];  //used to discriminate between Xi's sharing daughters
  bool acceptXi = false;                                     //used to see if Xi should be pushed to the collection
  AliFemtoXiCollection XiCorrectedCollection;                //collection to be returned

  AliFemtoXi* tXi;
  AliFemtoXiIterator XiIterator;
  AliFemtoXiIterator XiStart = XiCollection->begin();
  AliFemtoXiIterator XiEnd = XiCollection->end();

  for(XiIterator=XiStart; XiIterator!=XiEnd; XiIterator++)
  {
    tXi = *XiIterator;
    bool tmpPassXi = pCut->Pass(tXi);
    // it is not certain if particle which passed all the cuts will stay in the collection but if it fails - it fails for sure
    pCut->FillCutMonitor(tXi,tmpPassXi);
    if(tmpPassXi)
    {
      acceptXi = true;                                              //acceptXi stays true if there are no conflicts
      IdPosArray[count_pass] = tXi->IdPos();
      IdNegArray[count_pass] = tXi->IdNeg();
      IdBacArray[count_pass] = tXi->IdBac();
      dcaXiToPrimVertex[count_pass] = tXi->DcaXiToPrimVertex();

      for(Int_t ii=0; ii<count_pass; ii++) //used to check if in collection are already Xis with these daughters' IDs
      {
        if(IdPosArray[count_pass]==IdPosArray[ii] || IdNegArray[count_pass]==IdNegArray[ii] || IdBacArray[count_pass]==IdBacArray[ii])
        {
          if(dcaXiToPrimVertex[count_pass] < dcaXiToPrimVertex[ii])  //is new Xi "better" than old one?
          {
            //If true, need to remove old one from the collection
            for(AliFemtoXiIterator iter=XiCorrectedCollection.begin(); iter!=XiCorrectedCollection.end(); iter++)
            {
              if(positionInCollection==ii)
              {
                XiCorrectedCollection.insert(iter, 1, tXi);             //inserting new better Xi
                XiCorrectedCollection.erase(iter);                      //removing old Xi
                IdPosArray[ii] = IdPosArray[count_pass];                //update of positive V0 daughter's ID
                IdNegArray[ii] = IdNegArray[count_pass];                //update of negative V0 daughter's ID
                IdBacArray[ii] = IdBacArray[count_pass];                //update of bachelor pion's ID
                dcaXiToPrimVertex[ii] = dcaXiToPrimVertex[count_pass];  //update of Xi's DCA to primary vertex
                break;
              }
              positionInCollection++;
            }
            positionInCollection=0;
          }
          acceptXi=false; //conflict found - we don't need to change anything(worse DCA) or change has been already done
          break;
        }
      }
      if(acceptXi)  //If there are no conflicts found, add the Xi to the collection and increase the count_pass
      {
        XiCorrectedCollection.push_back(tXi);
        count_pass++;
      }
    }
  }

  delete [] IdPosArray;
  delete [] IdNegArray;
  delete [] IdBacArray;
  delete [] dcaXiToPrimVertex;

  return XiCorrectedCollection;
}
