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
   Double_t Value(int index, const Double_t* input) { return Value(index, input[0], input[1], input[2], input[3], input[4]); }
private:
   Double_t fInput0;   // neural network input neuron #1
   Double_t fInput1;   // neural network input neuron #2
   Double_t fInput2;   // neural network input neuron #3
   Double_t fInput3;   // neural network input neuron #4
   Double_t fInput4;   // neural network input neuron #5

   // private functions as exported by TMultiLayerPerceptron
   // containing the network data processing
   Double_t Neuron0x8ceb770() const;
   Double_t Neuron0x8ceb900() const;
   Double_t Neuron0x8cebad8() const;
   Double_t Neuron0x8cebcb0() const;
   Double_t Neuron0x93bc518() const;
   Double_t Input0x93bc838() const;
   Double_t Neuron0x93bc838() const;
   Double_t Input0x93bcab0() const;
   Double_t Neuron0x93bcab0() const;
   Double_t Input0x93bcd70() const;
   Double_t Neuron0x93bcd70() const;
   Double_t Input0x93bd030() const;
   Double_t Neuron0x93bd030() const;
   Double_t Input0x93bd378() const;
   Double_t Neuron0x93bd378() const;
   Double_t Input0x93bd5f0() const;
   Double_t Neuron0x93bd5f0() const;
   Double_t Input0x93bd8b0() const;
   Double_t Neuron0x93bd8b0() const;
   Double_t Input0x93bdc28() const;
   Double_t Neuron0x93bdc28() const;
   Double_t Input0x93bdee8() const;
   Double_t Neuron0x93bdee8() const;
   Double_t Input0x93be1a8() const;
   Double_t Neuron0x93be1a8() const;
   Double_t Input0x93bc708() const;
   Double_t Neuron0x93bc708() const;
   Double_t Input0x93be5d0() const;
   Double_t Neuron0x93be5d0() const;
   Double_t Synapse0x93bc9e8() const;
   Double_t Synapse0x93bca10() const;
   Double_t Synapse0x93bca38() const;
   Double_t Synapse0x93bca60() const;
   Double_t Synapse0x93bca88() const;
   Double_t Synapse0x93bcca8() const;
   Double_t Synapse0x93bccd0() const;
   Double_t Synapse0x93bccf8() const;
   Double_t Synapse0x93bcd20() const;
   Double_t Synapse0x93bcd48() const;
   Double_t Synapse0x93bcf68() const;
   Double_t Synapse0x93bcf90() const;
   Double_t Synapse0x93bcfb8() const;
   Double_t Synapse0x93bcfe0() const;
   Double_t Synapse0x93bd008() const;
   Double_t Synapse0x93bd228() const;
   Double_t Synapse0x93bd250() const;
   Double_t Synapse0x93bd300() const;
   Double_t Synapse0x93bd328() const;
   Double_t Synapse0x93bd350() const;
   Double_t Synapse0x93bd528() const;
   Double_t Synapse0x93bd550() const;
   Double_t Synapse0x93bd578() const;
   Double_t Synapse0x93bd5a0() const;
   Double_t Synapse0x93bd5c8() const;
   Double_t Synapse0x93bd7e8() const;
   Double_t Synapse0x93bd810() const;
   Double_t Synapse0x93bd838() const;
   Double_t Synapse0x93bd860() const;
   Double_t Synapse0x93bd888() const;
   Double_t Synapse0x93bdaa8() const;
   Double_t Synapse0x93bdad0() const;
   Double_t Synapse0x93bdaf8() const;
   Double_t Synapse0x8ce7098() const;
   Double_t Synapse0x8ce70c0() const;
   Double_t Synapse0x93bde20() const;
   Double_t Synapse0x93bde48() const;
   Double_t Synapse0x93bde70() const;
   Double_t Synapse0x93bde98() const;
   Double_t Synapse0x93bdec0() const;
   Double_t Synapse0x93be0e0() const;
   Double_t Synapse0x93be108() const;
   Double_t Synapse0x93be130() const;
   Double_t Synapse0x93be158() const;
   Double_t Synapse0x93be180() const;
   Double_t Synapse0x93be3a0() const;
   Double_t Synapse0x93be3c8() const;
   Double_t Synapse0x93be3f0() const;
   Double_t Synapse0x93be418() const;
   Double_t Synapse0x93be440() const;
   Double_t Synapse0x93bd2c0() const;
   Double_t Synapse0x93be468() const;
   Double_t Synapse0x93be490() const;
   Double_t Synapse0x93be4b8() const;
   Double_t Synapse0x93be4e0() const;
   Double_t Synapse0x93be508() const;
   Double_t Synapse0x93be530() const;
   Double_t Synapse0x93be558() const;
   Double_t Synapse0x93be580() const;
   Double_t Synapse0x93be5a8() const;
   Double_t Synapse0x93be7d0() const;
   Double_t Synapse0x93be7f8() const;
   Double_t Synapse0x93be820() const;
   Double_t Synapse0x93be848() const;
   Double_t Synapse0x93be870() const;
   Double_t Synapse0x934a7c8() const;
   Double_t Synapse0x93605e0() const;
   Double_t Synapse0x9360608() const;
   Double_t Synapse0x8ce6fe8() const;
   Double_t Synapse0x8ce7010() const;

	 ClassDef(AliCaloNeuralFit,1)
};

#endif // AliCaloNeuralFit_h

