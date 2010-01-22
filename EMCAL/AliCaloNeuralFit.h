#ifndef ALICALONEURALFIT_H
#define ALICALONEURALFIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id: $ */

//_________________________________________________________________________
//  Utility Class for Neural Network fit
//
//  currently uses 5 input neurons
//  network configured via TMultiLayerPerceptron
//
//*-- Author: Paola La Rocca (Catania)
//

#include <TObject.h>

class AliCaloNeuralFit : public TObject
{ 
public:
   AliCaloNeuralFit(): fInput0(0.0), fInput1(0.0), fInput2(0.0), fInput3(0.0), fInput4(0.0) {}
   ~AliCaloNeuralFit() {}
   Double_t Value(int index, Double_t in0, Double_t in1, Double_t in2, Double_t in3, Double_t in4);
   Double_t Value(int index, Double_t* input) { return Value(index, input[0], input[1], input[2], input[3], input[4]); }
private:
   Double_t fInput0;   // neural network input neuron #1
   Double_t fInput1;   // neural network input neuron #2
   Double_t fInput2;   // neural network input neuron #3
   Double_t fInput3;   // neural network input neuron #4
   Double_t fInput4;   // neural network input neuron #5

   // private functions as exported by TMultiLayerPerceptron
   // containing the network data processing
   Double_t Neuron0x945cbe8() const;
   Double_t Neuron0x945cd78() const;
   Double_t Neuron0x945cf50() const;
   Double_t Neuron0x945d128() const;
   Double_t Neuron0x945d300() const;
   Double_t Input0x945d620() const;
   Double_t Neuron0x945d620() const;
   Double_t Input0x945d870() const;
   Double_t Neuron0x945d870() const;
   Double_t Input0x945db30() const;
   Double_t Neuron0x945db30() const;
   Double_t Input0x945ddf0() const;
   Double_t Neuron0x945ddf0() const;
   Double_t Input0x945e138() const;
   Double_t Neuron0x945e138() const;
   Double_t Input0x945e3b0() const;
   Double_t Neuron0x945e3b0() const;
   Double_t Input0x945e670() const;
   Double_t Neuron0x945e670() const;
   Double_t Input0x945e9e8() const;
   Double_t Neuron0x945e9e8() const;
   Double_t Input0x945eca8() const;
   Double_t Neuron0x945eca8() const;
   Double_t Input0x945ef68() const;
   Double_t Neuron0x945ef68() const;
   Double_t Input0x945d4f0() const;
   Double_t Neuron0x945d4f0() const;
   Double_t Input0x945f490() const;
   Double_t Neuron0x945f490() const;
   Double_t Synapse0x943edb8() const;
   Double_t Synapse0x945d7d0() const;
   Double_t Synapse0x945d7f8() const;
   Double_t Synapse0x945d820() const;
   Double_t Synapse0x945d848() const;
   Double_t Synapse0x945da68() const;
   Double_t Synapse0x945da90() const;
   Double_t Synapse0x945dab8() const;
   Double_t Synapse0x945dae0() const;
   Double_t Synapse0x945db08() const;
   Double_t Synapse0x945dd28() const;
   Double_t Synapse0x945dd50() const;
   Double_t Synapse0x945dd78() const;
   Double_t Synapse0x945dda0() const;
   Double_t Synapse0x945ddc8() const;
   Double_t Synapse0x945dfe8() const;
   Double_t Synapse0x945e010() const;
   Double_t Synapse0x945e0c0() const;
   Double_t Synapse0x945e0e8() const;
   Double_t Synapse0x945e110() const;
   Double_t Synapse0x945e2e8() const;
   Double_t Synapse0x945e310() const;
   Double_t Synapse0x945e338() const;
   Double_t Synapse0x945e360() const;
   Double_t Synapse0x945e388() const;
   Double_t Synapse0x945e5a8() const;
   Double_t Synapse0x945e5d0() const;
   Double_t Synapse0x945e5f8() const;
   Double_t Synapse0x945e620() const;
   Double_t Synapse0x945e648() const;
   Double_t Synapse0x945e868() const;
   Double_t Synapse0x945e890() const;
   Double_t Synapse0x945e8b8() const;
   Double_t Synapse0x945e038() const;
   Double_t Synapse0x945e060() const;
   Double_t Synapse0x945ebe0() const;
   Double_t Synapse0x945ec08() const;
   Double_t Synapse0x945ec30() const;
   Double_t Synapse0x945ec58() const;
   Double_t Synapse0x945ec80() const;
   Double_t Synapse0x945eea0() const;
   Double_t Synapse0x945eec8() const;
   Double_t Synapse0x945eef0() const;
   Double_t Synapse0x945ef18() const;
   Double_t Synapse0x945ef40() const;
   Double_t Synapse0x945f160() const;
   Double_t Synapse0x945f188() const;
   Double_t Synapse0x945f1b0() const;
   Double_t Synapse0x945f1d8() const;
   Double_t Synapse0x945f200() const;
   Double_t Synapse0x945f300() const;
   Double_t Synapse0x945f328() const;
   Double_t Synapse0x945f350() const;
   Double_t Synapse0x945f378() const;
   Double_t Synapse0x945f3a0() const;
   Double_t Synapse0x945f3c8() const;
   Double_t Synapse0x945f3f0() const;
   Double_t Synapse0x945f418() const;
   Double_t Synapse0x945f440() const;
   Double_t Synapse0x945f468() const;
   Double_t Synapse0x945f690() const;
   Double_t Synapse0x945f6b8() const;
   Double_t Synapse0x945f6e0() const;
   Double_t Synapse0x945f708() const;
   Double_t Synapse0x945f730() const;
   Double_t Synapse0x936a1f0() const;
   Double_t Synapse0x943ee18() const;
   Double_t Synapse0x945cb70() const;
   Double_t Synapse0x945cb98() const;
   Double_t Synapse0x945cbc0() const;

   ClassDef(AliCaloNeuralFit,1)
};

#endif // AliCaloNeuralFit_h

