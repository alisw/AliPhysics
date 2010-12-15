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

/* $Id$ */

//.
// ITSupgrade base class to reconstruct an event
//.
//.
//.
#include "TObjArray.h"
#include "TTree.h"
#include "AliITSRecPoint.h"
#include "AliITSReconstructor.h"
#include "AliITSupgrade.h"
#include "AliITSUpgradeReconstructor.h" //class header
#include "AliITSDetTypeRec.h"
#include "AliITS.h"              //Reconstruct() 
#include "AliCDBEntry.h"           //ctor
#include "AliCDBManager.h"         //ctor
#include "AliESDEvent.h"           //FillEsd()
#include "AliRawReader.h"          //Reconstruct() for raw digits
#include "AliRun.h"
#include "AliLog.h"                //
#include "AliITSRawStream.h"     //ConvertDigits()
#include "AliRunLoader.h" 
#include "AliDataLoader.h"
#include "AliITSLoader.h"
#include "AliITSsegmentationUpgrade.h"
#include "AliITSUpgradeClusterFinder.h"
#include "AliITStrackerUpgrade.h"
#include "AliStack.h"
#include "TFile.h"
#include "TNtupleD.h"
ClassImp(AliITSUpgradeReconstructor)

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  AliITSUpgradeReconstructor::AliITSUpgradeReconstructor():
    AliITSReconstructor(), 
    fRecPoints(0),
    fNRecPoints(0),
    fDigits(0)

{
  //
  //ctor
  //
  fRecPoints = new TClonesArray("AliITSRecPoint",3000);
  fNRecPoints = 0;
  fDigits = new TObjArray(6);




}//AliITSReconstructor
//-----------------------------------------------------------------------
AliITSUpgradeReconstructor::~AliITSUpgradeReconstructor(){
  //Destructor
  if(fDigits){
    fDigits->Delete();
    delete fDigits;
    fDigits=0;
  }
  if(fRecPoints){
    fRecPoints->Delete();
    delete fRecPoints;
    fRecPoints=0;
  }
}

//_________________________________________________________________
void AliITSUpgradeReconstructor::Init() {
  // Initalize this constructor bet getting/creating the objects
  // nesseary for a proper ITS reconstruction.
  // Inputs:
  //   none.
  // Output:
  //   none.
  // Return:
  //   none.
  return;
}

//___________________________________________________________________________
void AliITSUpgradeReconstructor::SetTreeAddressD(TTree* const treeD){

  // Set branch address for the tree of digits.

  TBranch *branch;
  Int_t i;
  char branchname[30];
  if(!treeD) return;
  if (fDigits == 0x0) fDigits = new TObjArray(6);
  for (i=0; i<6; i++) {
    if(!(fDigits->At(i))) {
      fDigits->AddAt(new TClonesArray("AliITSDigitUpgrade",1000),i);
    }
    else{
      ResetDigits(i);
    }
    sprintf(branchname,"ITSDigits%d",i+1); 
    if (fDigits) {
      branch = treeD->GetBranch(branchname);
      if (branch) branch->SetAddress(&((*fDigits)[i]));
    }
  }
}
//__________________________________________________________________
void AliITSUpgradeReconstructor::SetTreeAddressR(TTree* const treeR){
  // Set branch address for the Reconstructed points Trees.
  // Inputs:
  //      TTree *treeR   Tree containing the RecPoints.
  // Outputs:
  //      none.
  // Return:
  char branchname[30];
  Char_t namedet[30]="ITS";
	
  if(!treeR) return;
  if(fRecPoints==0x0){
    fRecPoints = new TClonesArray("AliITSRecPoint",1000);
  }
	
  TBranch *branch;
  sprintf(branchname,"%sRecPoints",namedet);
  branch = treeR->GetBranch(branchname);
	
  if (branch) {
    branch->SetAddress(&fRecPoints);
  }
  else {
    sprintf(branchname,"%sRecPointsF",namedet); 
    branch = treeR->GetBranch(branchname);
    if (branch) {
      branch->SetAddress(&fRecPoints);
    }
  }

}
//_________________________________________________________________________
TBranch* AliITSUpgradeReconstructor::MakeBranchInTree(TTree* const tree,
						      const char* name, const char *classname,
						      void* address,Int_t size,Int_t splitlevel)
{
  //
  // Makes branch in given tree and diverts them to a separate file
  // 
  //
  //

  if (tree == 0x0) {
    Error("MakeBranchInTree","Making Branch %s Tree is NULL",name);
    return 0x0;
  }
  TBranch *branch = tree->GetBranch(name);
  if (branch) {
    return branch;
  }
  if (classname){
    branch = tree->Branch(name,classname,address,size,splitlevel);
  }
  else {
    branch = tree->Bronch(name, "TClonesArray", address, size, splitlevel);
  }

  return branch;
}
//________________________________________________________________________


void AliITSUpgradeReconstructor::MakeBranch(TTree* tree, Option_t* option){

  //Creates branches for clusters and recpoints
  Bool_t cR = (strstr(option,"R")!=0);
  Bool_t cRF = (strstr(option,"RF")!=0);

  if(cRF)cR = kFALSE;

  if(cR) MakeBranchR(tree);
  if(cRF) MakeBranchRF(tree);

}

//_________________________________________________________________
void AliITSUpgradeReconstructor::MakeBranchR(TTree *treeR, Option_t *opt){
  //Creates tree branches for recpoints
  // Inputs:
  //      cont char *file  File name where RecPoints branch is to be written
  //                       to. If blank it write the SDigits to the same
  //                       file in which the Hits were found.

  Int_t buffsz = 4000;
  char branchname[30];
  Bool_t oFast= (strstr(opt,"Fast")!=0);
  Char_t detname[30] = "ITS";

  if(oFast){
    sprintf(branchname,"%sRecPointsF",detname);
  } else {
    sprintf(branchname,"%sRecPoints",detname);
  }

  if(!fRecPoints)fRecPoints = new TClonesArray("AliITSRecPoint",1000);
  if (treeR)
    MakeBranchInTree(treeR,branchname,0,&fRecPoints,buffsz,99);
}

//_____________________________________________________________________
void AliITSUpgradeReconstructor::ResetDigits(){
  // Reset number of digits and the digits array for the ITS detector.

  if(!fDigits) return;
  for(Int_t i=0;i<6;i++){
    ResetDigits(i);
  }
}
//____________________________________________________________________
void AliITSUpgradeReconstructor::ResetDigits(Int_t branch){
  // Reset number of digits and the digits array for this branch.

  if(fDigits->At(branch)) ((TClonesArray*)fDigits->At(branch))->Clear();

}
//________________________________________________________________
void AliITSUpgradeReconstructor::AddRecPoint(const AliITSRecPoint &r){
  // Add a reconstructed space point to the list
  // Inputs:
  //      const AliITSRecPoint &r RecPoint class to be added to the tree
  //                              of reconstructed points TreeR.
  // Outputs:
  //      none.
  // Return:
  //      none.

  TClonesArray &lrecp = *fRecPoints;
  new(lrecp[fNRecPoints++]) AliITSRecPoint(r);
}

//___________________________________________________________________
void AliITSUpgradeReconstructor::DigitsToRecPoints(TTree *treeD,TTree *treeR){
  AliITSsegmentationUpgrade *segmentation2 = 0x0;
  AliITSRecPoint  recpnt; 
  if(!segmentation2)
  segmentation2 = new AliITSsegmentationUpgrade();
  AliITSUpgradeClusterFinder *clf = new AliITSUpgradeClusterFinder();
  TObjArray *digList= new TObjArray(6);
  for(Int_t il=0; il<6; il ++) digList->AddAt(new TClonesArray("AliITSDigitUpgrade"),il);
  for(Int_t ilayer=0; ilayer<6; ilayer ++) {
    treeD->SetBranchAddress(Form("Layer%d",ilayer),&(*digList)[ilayer]);
    //treeD->GetListOfLeaves()->Print();	
  }//loop layer per tbranch
  treeD->GetEntry(0);
		
  for(Int_t ilayer=0; ilayer < 6 ;ilayer ++){
    TClonesArray *pArrDig= (TClonesArray*)digList->At(ilayer);
    clf->StartEvent();
    for(Int_t ientr =0; ientr < pArrDig->GetEntries() ; ientr++){
      AliITSDigitUpgrade *dig = (AliITSDigitUpgrade*)pArrDig->At(ientr);
      Int_t colz=dig->GetzPixelNumber();
      Int_t rowx=dig->GetxPixelNumber();
      Double_t hitcharge= (dig->GetNelectrons());
      clf->ProcessHitOnline(ilayer,colz, rowx,(Short_t)hitcharge,dig->GetTracks());


    }//ientr
    clf->FinishEvent();
    for(UInt_t nClu = 0; nClu <  clf -> GetClusterCount(ilayer); nClu++){
      UShort_t charge = clf->GetCharge(ilayer, nClu);
      recpnt.SetQ(charge);
      recpnt.SetLayer(ilayer);
      Int_t *lab=clf->GetLabels(ilayer,nClu);
for(Int_t l=0; l<3; l++) {
                        recpnt.SetLabel(lab[l],l);
                        }
if(clf->GetClusterSize(ilayer,nClu)==2){
                        printf("\n\n CLUSTER LABELS (size %i) \n",clf->GetClusterSize(ilayer,nClu));
                        for(Int_t i=0; i<10; i++ ){
                         printf(" %i ",lab[i]);
                         }
                        printf("\n\n");
                       }


 			
      Bool_t check2;
      Double_t xcheck2=0.;
      Double_t ycheck2=0.;
      Double_t zcheck2=0.;
      Double_t xzl2[2]={0.,0.};
      Double_t XpixC2,ZpixC2=0.;

      XpixC2 = clf-> GetClusterMeanRow(ilayer, nClu);
      ZpixC2 = clf-> GetClusterMeanCol(ilayer, nClu);
      xzl2[0] = XpixC2*(segmentation2->GetCellSizeX(ilayer))+0.5*(segmentation2-> GetCellSizeX(ilayer));
      xzl2[1] = ZpixC2*(segmentation2->GetCellSizeZ(ilayer))+0.5*(segmentation2->GetCellSizeZ(ilayer))-(segmentation2->GetHalfLength(ilayer));
      check2 = segmentation2->DetToGlobal(ilayer,xzl2[0], xzl2[1],xcheck2,ycheck2,zcheck2);
      recpnt.SetType(clf->GetClusterType(ilayer,nClu ));
    // recpnt.SetLocalCoord(xzl2[0],xzl2[1]); //temporary solution (no LocalToTrack Matrix)
         //from global to tracking system coordinate
	// global coordinate -> local coordinate getting alpha angle of the recpoint
              Float_t xclg = xcheck2;//upgrade clusters global coordinate ( ITS official: GetX tracking coordinate)
              Float_t yclg = ycheck2;
              Float_t zclg = zcheck2;
              Double_t phiclu1rad, phiclu1deg;
              phiclu1rad=TMath::ATan2(yclg,xclg);//cluster phi angle (rad)
              if (phiclu1rad<0) phiclu1rad+=TMath::TwoPi();//from 0 to 360
              else if (phiclu1rad>=TMath::TwoPi()) phiclu1rad-=TMath::TwoPi();//

              phiclu1deg=180.*phiclu1rad/TMath::Pi();// in deg
              Int_t ladder;// virtual segmentation starting from the cluster phi

              ladder=(Int_t)phiclu1deg/18;// in which ladder the cluster is
              Double_t alpha= (ladder*18.+9.)*TMath::Pi()/180.;//angle at the center of the ladder (rad)

              //alpha rotation 
              Float_t xclu1_tr = xclg*TMath::Cos(alpha)-yclg*TMath::Sin(alpha);
              Float_t yclu1 = yclg*TMath::Cos(alpha)+xclg*TMath::Sin(alpha);
              Float_t xclu1 = TMath::Sqrt(xclu1_tr*xclu1_tr+yclu1*yclu1);
              Float_t zclu1 = zclg;
              Double_t phi_trk= (phiclu1rad-alpha);// cluster angle in the rotated system (rad)

              yclu1=xclu1*phi_trk; // tracking system coordinate: r*phi
      recpnt.SetX(0.);
      recpnt.SetY(yclu1);
      recpnt.SetZ(zclu1);
      //  
      Double_t xsize, zsize;
      segmentation2->GetSegmentation(ilayer,xsize, zsize);
      recpnt.SetSigmaY2(xsize/TMath::Sqrt(12)*xsize/TMath::Sqrt(12));
      recpnt.SetSigmaZ2(zsize/TMath::Sqrt(12)*zsize/TMath::Sqrt(12));
      TClonesArray &lrecp = *fRecPoints;
      new(lrecp[fNRecPoints++]) AliITSRecPoint(recpnt);
    }//cluster list entries
    treeR->Fill();
    ResetRecPoints();
  }//ilayer
}
//_______________________________________________________________________________________________________________
AliTracker* AliITSUpgradeReconstructor::CreateTracker() const
{
        // create a ITS tracker
   //     Int_t trackerOpt = GetRecoParam()->GetTracker();
    //    AliTracker* tracker;

      //          tracker = new AliITStrackerUpgrade();
                AliITStrackerUpgrade *trackUp = new AliITStrackerUpgrade();
                if(GetRecoParam()->GetTrackerSAOnly()) trackUp->SetSAFlag(kTRUE);
                if(trackUp->GetSAFlag())AliDebug(1,"Tracking Performed in ITS only\n");
                if(GetRecoParam()->GetInwardFindingSA()){
                        trackUp->SetInwardFinding();
                        trackUp->SetInnerStartLayer(GetRecoParam()->GetInnerStartLayerSA());
                }else{
                        trackUp->SetOutwardFinding();
                        trackUp->SetOuterStartLayer(GetRecoParam()->GetOuterStartLayerSA());
                }
                trackUp->SetMinNPoints(GetRecoParam()->GetMinNPointsSA());
                return trackUp;
}



