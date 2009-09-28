#ifndef ALITRDPIDREFMAKERNN_H
#define ALITRDPIDREFMAKERNN_H

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

class TMultiLayerPerceptron;
class AliTRDpidRefMakerNN : public AliTRDpidRefMaker
{
public:
  enum ETRDpidRefMakerNNgraph {
    kGraphTrain = 0
    ,kGraphTest = 1
  };

  enum ETRDpidRefMakerNNmoni {
    kMoniTrain = 50
  };

  AliTRDpidRefMakerNN();

  virtual ~AliTRDpidRefMakerNN();
  
  void    CreateOutputObjects();

  Int_t   GetEpochs() {return fEpochs;};
  Int_t   GetMinTrain() {return fMinTrain;};
  Int_t   GetTrainMomBin() {return fTrainMomBin;};

  Bool_t  PostProcess();

  void    SetEpochs(Int_t epochs) {fEpochs = epochs;};
  void    SetMinTrain(Int_t mintrain) {fMinTrain = mintrain;};
  void    SetTrainMomBin(Int_t trainmombin) {fTrainMomBin = trainmombin;};
  void    SetDate(Int_t date) {fDate = date;};
  void    SetDoTraining(Bool_t train) {fDoTraining = train;};
  void    SetContinueTraining(Bool_t continTrain) {fContinueTraining = continTrain;};
  void    SetTrainPath(Int_t path) {fTrainPath = path;};


  void    MonitorTraining(Int_t mombin);                       // monitor training process

protected:
  void MakeRefs(Int_t mombin);                         // train the neural networks for a given momentum bin

private:
  AliTRDpidRefMakerNN(const AliTRDpidRefMakerNN&);              // not implemented
  AliTRDpidRefMakerNN& operator=(const AliTRDpidRefMakerNN&);   // not implemented


  TMultiLayerPerceptron *fNet[AliTRDgeometry::kNlayer]; // artificial neural network

  Int_t         fTrainMomBin;              // momentum bin for the training
  Int_t         fEpochs;                   // Number of epochs for the training of the NNs
  Int_t         fMinTrain;                 // minimum of events needed for training
  Int_t         fDate;                     // date stamp for training of the NNs
  Bool_t        fDoTraining;               // checks if training will be done
  Bool_t        fContinueTraining;         // checks if training from an older run should be continued
  Int_t         fTrainPath;                // sets the path for continuing the training

  ClassDef(AliTRDpidRefMakerNN, 2); // TRD PID reference  maker for NN
};

#endif
