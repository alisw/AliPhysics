// This macro converts the networks that were saved by the AliTRDpidRefMaker to a root file (NN.root)
// which is readable by the AliTRDpidCDB.C macro to create the pid reference data bases.

#ifndef __CINT__
#include "TFile.h"
#include "TMultiLayerPerceptron.h"
#include "TEventList.h"
#include "TTree.h"
#include "TObject.h"

#include "AliTRDgeometry.h"
#include "AliTRDpidUtil.h"
#include "Cal/AliTRDCalPID.h"
#endif

Int_t ConvertMLPs(Int_t iDate = 20090129){

  const Int_t nCha = AliTRDgeometry::kNlayer;

  TFile *fOut = new TFile("NN.root","RECREATE");
    
  Float_t fv0pid[AliPID::kSPECIES], fdEdx[AliTRDpidUtil::kNNslices];

  for(Int_t iMomBin = 0; iMomBin <AliTRDCalPID::kNMom; iMomBin++){

// Load the TFile and the TTree to load the weights of the TMultiLayerPerceptron
    TFile *fIn = new TFile(Form("./%3.1fGeV/TRD.TaskPidRefMakerNN.root",AliTRDCalPID::GetMomentum(iMomBin)),"READ");
    TTree *tIn = (TTree *) fIn -> Get("NN");

// Create dummy TEventLists to load the weights of the TMultiLayerPerceptron
    TEventList* lTrainData = new TEventList();
    TEventList* lTestData = new TEventList();
    lTrainData -> Enter(1);
    lTestData -> Enter(1);
  
    tIn -> SetEventList(lTrainData);
    tIn -> SetEventList(lTestData);
    
    Char_t structure[1024];

    TMultiLayerPerceptron *mNet[nCha];
    
    for(Int_t iCha = 0; iCha < nCha; iCha++){
      sprintf(structure,"fdEdx[0],fdEdx[1],fdEdx[2],fdEdx[3],fdEdx[4],fdEdx[5],fdEdx[6],fdEdx[7]:15:7:fv0pid[0],fv0pid[1],fv0pid[2],fv0pid[3],fv0pid[4]!");      
      mNet[iCha] = new TMultiLayerPerceptron(structure,tIn,lTrainData,lTestData);

      printf("Network[%d] is ./%3.1fGeV/Network_%d/MomBin_%d/Net%d_49\n", iCha, AliTRDCalPID::GetMomentum(iMomBin), iDate, iMomBin, iCha);
      mNet[iCha] -> LoadWeights(Form("./%3.1fGeV/Network_%d/MomBin_%d/Net%d_49", AliTRDCalPID::GetMomentum(iMomBin), iDate, iMomBin, iCha));
    }
    
    fOut -> cd();
    for(Int_t iCha = 0; iCha < AliTRDgeometry::kNlayer; iCha++){
      mNet[iCha] -> Write(Form("NN_Mom%d_Plane%d",iMomBin,iCha));
    }
  }

  fOut -> Close();

  return 1;
    
}

