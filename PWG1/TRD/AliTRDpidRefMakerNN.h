#ifndef ALITRDPIDREFMAKERNN_H
#define ALITRDPIDREFMAKERNN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDpidRefMakerNN.h 27496 2008-07-22 08:35:45Z cblume $ */

//////////////////////////////////////////////////////
//
// Task to build PID reference tree for the training
// of neural networs for the TRD PID
//
// Authors: Alex Wilk    <wilka@uni-muenster.de>
//          Markus Heide <mheide@uni-muenster.de>
//
///////////////////////////////////////////////////////

#ifndef ALITRDPIDREFMAKER_H
#include "AliTRDpidRefMaker.h"
#endif

class TEventList;
class TMultiLayerPerceptron;
class AliTRDpidRefMakerNN : public AliTRDpidRefMaker
{

public:
  enum ETRDpidRefMakerNNgraph {
    kGraphTrain = 1
    ,kGraphTest = 2
  };

  enum ETRDpidRefMakerNNmoni {
    kMoniTrain = 50
  };

  enum ETRDpidRefMakerNNsteer{
    kMaxStat    = 20000 // maximum statistics/PID bin
   ,kMinStat    = 50     // minimum statistics/PID bin
  };
  AliTRDpidRefMakerNN();
  AliTRDpidRefMakerNN(const char *n);
  virtual ~AliTRDpidRefMakerNN();
  
  Int_t   GetEpochs() const {return fEpochs;};
  Int_t   GetMinTrain() const {return fMinTrain;};
  Int_t   GetTrainMomBin() const {return fTrainMomBin;};

  Bool_t  PostProcess();

  void    SetEpochs(Int_t epochs) {fEpochs = epochs;};
  void    SetMinTrain(Int_t mintrain) {fMinTrain = mintrain;};
  void    SetTrainMomBin(Int_t trainmombin) {fTrainMomBin = trainmombin;};
  void    SetDate(Int_t date) {fDate = date;};
  void    SetDoTraining(Bool_t train) {fDoTraining = train;};
  void    SetContinueTraining(Bool_t continTrain) {fContinueTraining = continTrain;};
  void    SetTrainPath(Int_t path) {fTrainPath = path;};
  Bool_t  LoadFile(const Char_t *InFileNN);
  void    SetScaledEdx(Float_t s) {fScale = s;};

  void    MakeTrainingLists(Int_t mombin = 0);                                 // build the training and the test list
  void    MonitorTraining(Int_t mombin);                       // monitor training process

protected:
  void MakeRefs(Int_t mombin);                         // train the neural networks for a given momentum bin

private:
  AliTRDpidRefMakerNN(const AliTRDpidRefMakerNN&);              // not implemented
  AliTRDpidRefMakerNN& operator=(const AliTRDpidRefMakerNN&);   // not implemented
  void    MakeTrainTestTrees();

  TEventList *fTrain[AliTRDCalPID::kNMom];          // Training list for each momentum
  TEventList *fTest[AliTRDCalPID::kNMom];           // Test list for each momentum
  TMultiLayerPerceptron *fNet; // artificial neural network

  Int_t         fTrainMomBin;              // momentum bin for the training
  Int_t         fEpochs;                   // Number of epochs for the training of the NNs
  Int_t         fMinTrain;                 // minimum of events needed for training
  Int_t         fDate;                     // date stamp for training of the NNs
  Bool_t        fDoTraining;               // checks if training will be done
  Bool_t        fContinueTraining;         // checks if training from an older run should be continued
  Int_t         fTrainPath;                // sets the path for continuing the training

  Float_t       fScale;
  Int_t         fLy;                       // TRD layer
  Int_t         fNtrkl;                    // No. tracklets
  TTree         *fTrainData[AliTRDCalPID::kNMom]; 
  TFile         *fRef;

  ClassDef(AliTRDpidRefMakerNN, 3); // TRD reference  maker for NN
};

#endif
