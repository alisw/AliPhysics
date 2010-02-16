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
   Double_t Neuron0x9fd0208() const;
   Double_t Neuron0x9fd0398() const;
   Double_t Neuron0x9fd0570() const;
   Double_t Neuron0x9fd0748() const;
   Double_t Neuron0x9fdec20() const;
   Double_t Input0x9fdef28() const;
   Double_t Neuron0x9fdef28() const;
   Double_t Input0x9fdf1a0() const;
   Double_t Neuron0x9fdf1a0() const;
   Double_t Input0x9fdf460() const;
   Double_t Neuron0x9fdf460() const;
   Double_t Input0x9fdf720() const;
   Double_t Neuron0x9fdf720() const;
   Double_t Input0x9fdfa68() const;
   Double_t Neuron0x9fdfa68() const;
   Double_t Input0x9fdfce0() const;
   Double_t Neuron0x9fdfce0() const;
   Double_t Input0x9fdffa0() const;
   Double_t Neuron0x9fdffa0() const;
   Double_t Input0x9fe0318() const;
   Double_t Neuron0x9fe0318() const;
   Double_t Input0x9fe05b8() const;
   Double_t Neuron0x9fe05b8() const;
   Double_t Input0x9fe0878() const;
   Double_t Neuron0x9fe0878() const;
   Double_t Input0x9fdedf8() const;
   Double_t Neuron0x9fdedf8() const;
   Double_t Input0x9fe0da0() const;
   Double_t Neuron0x9fe0da0() const;
   Double_t Synapse0x9fdf0d8() const;
   Double_t Synapse0x9fdf100() const;
   Double_t Synapse0x9fdf128() const;
   Double_t Synapse0x9fdf150() const;
   Double_t Synapse0x9fdf178() const;
   Double_t Synapse0x9fdf398() const;
   Double_t Synapse0x9fdf3c0() const;
   Double_t Synapse0x9fdf3e8() const;
   Double_t Synapse0x9fdf410() const;
   Double_t Synapse0x9fdf438() const;
   Double_t Synapse0x9fdf658() const;
   Double_t Synapse0x9fdf680() const;
   Double_t Synapse0x9fdf6a8() const;
   Double_t Synapse0x9fdf6d0() const;
   Double_t Synapse0x9fdf6f8() const;
   Double_t Synapse0x9fdf918() const;
   Double_t Synapse0x9fdf940() const;
   Double_t Synapse0x9fdf9f0() const;
   Double_t Synapse0x9fdfa18() const;
   Double_t Synapse0x9fdfa40() const;
   Double_t Synapse0x9fdfc18() const;
   Double_t Synapse0x9fdfc40() const;
   Double_t Synapse0x9fdfc68() const;
   Double_t Synapse0x9fdfc90() const;
   Double_t Synapse0x9fdfcb8() const;
   Double_t Synapse0x9fdfed8() const;
   Double_t Synapse0x9fdff00() const;
   Double_t Synapse0x9fdff28() const;
   Double_t Synapse0x9fdff50() const;
   Double_t Synapse0x9fdff78() const;
   Double_t Synapse0x9fe0198() const;
   Double_t Synapse0x9fe01c0() const;
   Double_t Synapse0x9fe01e8() const;
   Double_t Synapse0x9fdf968() const;
   Double_t Synapse0x9fdf990() const;
   Double_t Synapse0x9fe04f0() const;
   Double_t Synapse0x9fe0518() const;
   Double_t Synapse0x9fe0540() const;
   Double_t Synapse0x9fe0568() const;
   Double_t Synapse0x9fe0590() const;
   Double_t Synapse0x9fe07b0() const;
   Double_t Synapse0x9fe07d8() const;
   Double_t Synapse0x9fe0800() const;
   Double_t Synapse0x9fe0828() const;
   Double_t Synapse0x9fe0850() const;
   Double_t Synapse0x9fe0a70() const;
   Double_t Synapse0x9fe0a98() const;
   Double_t Synapse0x9fe0ac0() const;
   Double_t Synapse0x9fe0ae8() const;
   Double_t Synapse0x9fe0b10() const;
   Double_t Synapse0x9fe0c10() const;
   Double_t Synapse0x9fe0c38() const;
   Double_t Synapse0x9fe0c60() const;
   Double_t Synapse0x9fe0c88() const;
   Double_t Synapse0x9fe0cb0() const;
   Double_t Synapse0x9fe0cd8() const;
   Double_t Synapse0x9fe0d00() const;
   Double_t Synapse0x9fe0d28() const;
   Double_t Synapse0x9fe0d50() const;
   Double_t Synapse0x9fe0d78() const;
   Double_t Synapse0x9fe0fa0() const;
   Double_t Synapse0x9fe0fc8() const;
   Double_t Synapse0x9fe0ff0() const;
   Double_t Synapse0x9fe1018() const;
   Double_t Synapse0x9fe1040() const;
   Double_t Synapse0x9882b78() const;
   Double_t Synapse0x9fd0158() const;
   Double_t Synapse0x9fd0180() const;
   Double_t Synapse0x9fd01a8() const;
   Double_t Synapse0x9fd01d0() const;

	 ClassDef(AliCaloNeuralFit,1)
};

#endif // AliCaloNeuralFit_h

