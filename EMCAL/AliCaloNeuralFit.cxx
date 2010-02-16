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

/* $Id: $ */

//_________________________________________________________________________
//  Utility Class for Neural Network fit
//
//  currently uses 5 input neurons
//  network configured via TMultiLayerPerceptron
//
//*-- Author: Paola La Rocca (Catania)
//

#include "AliCaloNeuralFit.h"
#include <cmath>


Double_t AliCaloNeuralFit::Value
(int index, Double_t in0, Double_t in1, Double_t in2, Double_t in3, Double_t in4)
{
//
// Compute the neural network answer,
// given the input values (taken from the signal TGraph)
//

   fInput0 = in0;
   fInput1 = in1;
   fInput2 = in2;
   fInput3 = in3;
   fInput4 = in4;
   switch(index) 
	 {
     case 0:
         return Neuron0x9fdedf8();
     case 1:
         return Neuron0x9fe0da0();
     default:
         return 0.;
   }
}

Double_t AliCaloNeuralFit::Neuron0x9fd0208() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput0;
}

Double_t AliCaloNeuralFit::Neuron0x9fd0398() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput1;
}

Double_t AliCaloNeuralFit::Neuron0x9fd0570() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput2;
}

Double_t AliCaloNeuralFit::Neuron0x9fd0748() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput3;
}

Double_t AliCaloNeuralFit::Neuron0x9fdec20() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput4;
}

Double_t AliCaloNeuralFit::Input0x9fdef28() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//

   Double_t input = 1.01367;
   input += Synapse0x9fdf0d8();
   input += Synapse0x9fdf100();
   input += Synapse0x9fdf128();
   input += Synapse0x9fdf150();
   input += Synapse0x9fdf178();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fdef28() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fdef28();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fdf1a0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.388335;
   input += Synapse0x9fdf398();
   input += Synapse0x9fdf3c0();
   input += Synapse0x9fdf3e8();
   input += Synapse0x9fdf410();
   input += Synapse0x9fdf438();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fdf1a0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fdf1a0();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fdf460() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.547781;
   input += Synapse0x9fdf658();
   input += Synapse0x9fdf680();
   input += Synapse0x9fdf6a8();
   input += Synapse0x9fdf6d0();
   input += Synapse0x9fdf6f8();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fdf460() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fdf460();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fdf720() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.525561;
   input += Synapse0x9fdf918();
   input += Synapse0x9fdf940();
   input += Synapse0x9fdf9f0();
   input += Synapse0x9fdfa18();
   input += Synapse0x9fdfa40();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fdf720() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fdf720();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fdfa68() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.360225;
   input += Synapse0x9fdfc18();
   input += Synapse0x9fdfc40();
   input += Synapse0x9fdfc68();
   input += Synapse0x9fdfc90();
   input += Synapse0x9fdfcb8();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fdfa68() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fdfa68();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fdfce0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.485774;
   input += Synapse0x9fdfed8();
   input += Synapse0x9fdff00();
   input += Synapse0x9fdff28();
   input += Synapse0x9fdff50();
   input += Synapse0x9fdff78();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fdfce0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fdfce0();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fdffa0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.467242;
   input += Synapse0x9fe0198();
   input += Synapse0x9fe01c0();
   input += Synapse0x9fe01e8();
   input += Synapse0x9fdf968();
   input += Synapse0x9fdf990();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fdffa0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fdffa0();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fe0318() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.566875;
   input += Synapse0x9fe04f0();
   input += Synapse0x9fe0518();
   input += Synapse0x9fe0540();
   input += Synapse0x9fe0568();
   input += Synapse0x9fe0590();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fe0318() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fe0318();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fe05b8() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.848776;
   input += Synapse0x9fe07b0();
   input += Synapse0x9fe07d8();
   input += Synapse0x9fe0800();
   input += Synapse0x9fe0828();
   input += Synapse0x9fe0850();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fe05b8() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fe05b8();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fe0878() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.408427;
   input += Synapse0x9fe0a70();
   input += Synapse0x9fe0a98();
   input += Synapse0x9fe0ac0();
   input += Synapse0x9fe0ae8();
   input += Synapse0x9fe0b10();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fe0878() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fe0878();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fdedf8() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.310322;
   input += Synapse0x9fe0c10();
   input += Synapse0x9fe0c38();
   input += Synapse0x9fe0c60();
   input += Synapse0x9fe0c88();
   input += Synapse0x9fe0cb0();
   input += Synapse0x9fe0cd8();
   input += Synapse0x9fe0d00();
   input += Synapse0x9fe0d28();
   input += Synapse0x9fe0d50();
   input += Synapse0x9fe0d78();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fdedf8() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fdedf8();
   return (input * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x9fe0da0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.0574773;
   input += Synapse0x9fe0fa0();
   input += Synapse0x9fe0fc8();
   input += Synapse0x9fe0ff0();
   input += Synapse0x9fe1018();
   input += Synapse0x9fe1040();
   input += Synapse0x9882b78();
   input += Synapse0x9fd0158();
   input += Synapse0x9fd0180();
   input += Synapse0x9fd01a8();
   input += Synapse0x9fd01d0();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x9fe0da0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x9fe0da0();
   return (input * 1)+0;
}

Double_t AliCaloNeuralFit::Synapse0x9fdf0d8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*1.53012);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf100()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*-0.316606);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf128()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*1.31047);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf150()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*0.31846);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf178()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*-1.43145);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf398()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*-0.199402);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf3c0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*0.0250046);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf3e8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*0.21622);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf410()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*0.0240984);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf438()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*0.492242);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf658()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*-2.9778);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf680()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*-0.598961);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf6a8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*-0.857305);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf6d0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*0.58472);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf6f8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*1.87975);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf918()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*0.334053);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf940()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*0.142395);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf9f0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*0.293019);
}

Double_t AliCaloNeuralFit::Synapse0x9fdfa18()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*-0.109163);
}

Double_t AliCaloNeuralFit::Synapse0x9fdfa40()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*0.482151);
}

Double_t AliCaloNeuralFit::Synapse0x9fdfc18()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*-0.270303);
}

Double_t AliCaloNeuralFit::Synapse0x9fdfc40()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*-0.212234);
}

Double_t AliCaloNeuralFit::Synapse0x9fdfc68()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*-0.136993);
}

Double_t AliCaloNeuralFit::Synapse0x9fdfc90()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*-0.322175);
}

Double_t AliCaloNeuralFit::Synapse0x9fdfcb8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*-0.137984);
}

Double_t AliCaloNeuralFit::Synapse0x9fdfed8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*0.00656827);
}

Double_t AliCaloNeuralFit::Synapse0x9fdff00()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*-0.0405237);
}

Double_t AliCaloNeuralFit::Synapse0x9fdff28()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*-0.437988);
}

Double_t AliCaloNeuralFit::Synapse0x9fdff50()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*-0.343747);
}

Double_t AliCaloNeuralFit::Synapse0x9fdff78()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*-0.168597);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0198()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*-0.47585);
}

Double_t AliCaloNeuralFit::Synapse0x9fe01c0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*1.7023);
}

Double_t AliCaloNeuralFit::Synapse0x9fe01e8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*0.193432);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf968()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*0.139333);
}

Double_t AliCaloNeuralFit::Synapse0x9fdf990()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*-0.400241);
}

Double_t AliCaloNeuralFit::Synapse0x9fe04f0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*0.757465);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0518()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*0.070475);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0540()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*0.412929);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0568()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*-0.161017);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0590()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*-0.168403);
}

Double_t AliCaloNeuralFit::Synapse0x9fe07b0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*0.490089);
}

Double_t AliCaloNeuralFit::Synapse0x9fe07d8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*-0.193417);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0800()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*-0.345063);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0828()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*-0.507424);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0850()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*-0.790068);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0a70()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0208()*1.717);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0a98()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0398()*3.29133);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0ac0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0570()*-0.60354);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0ae8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fd0748()*-0.553822);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0b10()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdec20()*-0.292983);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0c10()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdef28()*0.0177982);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0c38()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdf1a0()*0.088043);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0c60()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdf460()*0.0698223);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0c88()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdf720()*0.609632);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0cb0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdfa68()*-0.825672);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0cd8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdfce0()*-0.109339);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0d00()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdffa0()*-0.498954);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0d28()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fe0318()*0.346775);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0d50()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fe05b8()*-0.284703);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0d78()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fe0878()*0.0176923);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0fa0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdef28()*-1.65881);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0fc8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdf1a0()*0.0393971);
}

Double_t AliCaloNeuralFit::Synapse0x9fe0ff0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdf460()*1.71678);
}

Double_t AliCaloNeuralFit::Synapse0x9fe1018()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdf720()*0.408921);
}

Double_t AliCaloNeuralFit::Synapse0x9fe1040()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdfa68()*-0.508733);
}

Double_t AliCaloNeuralFit::Synapse0x9882b78()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdfce0()*-0.48213);
}

Double_t AliCaloNeuralFit::Synapse0x9fd0158()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fdffa0()*-1.36281);
}

Double_t AliCaloNeuralFit::Synapse0x9fd0180()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fe0318()*0.102217);
}

Double_t AliCaloNeuralFit::Synapse0x9fd01a8()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fe05b8()*-0.684593);
}

Double_t AliCaloNeuralFit::Synapse0x9fd01d0()  const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x9fe0878()*1.99419);
}
