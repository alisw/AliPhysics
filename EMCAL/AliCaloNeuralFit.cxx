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
         return Neuron0x93bc708();
     case 1:
         return Neuron0x93be5d0();
     default:
         return 0.;
   }
}

Double_t AliCaloNeuralFit::Neuron0x8ceb770() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput0;
}

Double_t AliCaloNeuralFit::Neuron0x8ceb900() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput1;
}

Double_t AliCaloNeuralFit::Neuron0x8cebad8() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput2;
}

Double_t AliCaloNeuralFit::Neuron0x8cebcb0() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput3;
}

Double_t AliCaloNeuralFit::Neuron0x93bc518() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput4;
}

Double_t AliCaloNeuralFit::Input0x93bc838() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.097207;
   input += Synapse0x93bc9e8();
   input += Synapse0x93bca10();
   input += Synapse0x93bca38();
   input += Synapse0x93bca60();
   input += Synapse0x93bca88();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bc838() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bc838();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93bcab0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.0466086;
   input += Synapse0x93bcca8();
   input += Synapse0x93bccd0();
   input += Synapse0x93bccf8();
   input += Synapse0x93bcd20();
   input += Synapse0x93bcd48();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bcab0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bcab0();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93bcd70() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.617042;
   input += Synapse0x93bcf68();
   input += Synapse0x93bcf90();
   input += Synapse0x93bcfb8();
   input += Synapse0x93bcfe0();
   input += Synapse0x93bd008();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bcd70() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bcd70();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93bd030() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.519529;
   input += Synapse0x93bd228();
   input += Synapse0x93bd250();
   input += Synapse0x93bd300();
   input += Synapse0x93bd328();
   input += Synapse0x93bd350();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bd030() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bd030();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93bd378() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.405186;
   input += Synapse0x93bd528();
   input += Synapse0x93bd550();
   input += Synapse0x93bd578();
   input += Synapse0x93bd5a0();
   input += Synapse0x93bd5c8();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bd378() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bd378();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93bd5f0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.253465;
   input += Synapse0x93bd7e8();
   input += Synapse0x93bd810();
   input += Synapse0x93bd838();
   input += Synapse0x93bd860();
   input += Synapse0x93bd888();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bd5f0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bd5f0();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93bd8b0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.0340672;
   input += Synapse0x93bdaa8();
   input += Synapse0x93bdad0();
   input += Synapse0x93bdaf8();
   input += Synapse0x8ce7098();
   input += Synapse0x8ce70c0();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bd8b0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bd8b0();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93bdc28() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.231599;
   input += Synapse0x93bde20();
   input += Synapse0x93bde48();
   input += Synapse0x93bde70();
   input += Synapse0x93bde98();
   input += Synapse0x93bdec0();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bdc28() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bdc28();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93bdee8() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.56831;
   input += Synapse0x93be0e0();
   input += Synapse0x93be108();
   input += Synapse0x93be130();
   input += Synapse0x93be158();
   input += Synapse0x93be180();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bdee8() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bdee8();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93be1a8() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.223889;
   input += Synapse0x93be3a0();
   input += Synapse0x93be3c8();
   input += Synapse0x93be3f0();
   input += Synapse0x93be418();
   input += Synapse0x93be440();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93be1a8() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93be1a8();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93bc708() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.565736;
   input += Synapse0x93bd2c0();
   input += Synapse0x93be468();
   input += Synapse0x93be490();
   input += Synapse0x93be4b8();
   input += Synapse0x93be4e0();
   input += Synapse0x93be508();
   input += Synapse0x93be530();
   input += Synapse0x93be558();
   input += Synapse0x93be580();
   input += Synapse0x93be5a8();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93bc708() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93bc708();
   return (input * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x93be5d0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.311772;
   input += Synapse0x93be7d0();
   input += Synapse0x93be7f8();
   input += Synapse0x93be820();
   input += Synapse0x93be848();
   input += Synapse0x93be870();
   input += Synapse0x934a7c8();
   input += Synapse0x93605e0();
   input += Synapse0x9360608();
   input += Synapse0x8ce6fe8();
   input += Synapse0x8ce7010();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x93be5d0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x93be5d0();
   return (input * 1)+0;
}

Double_t AliCaloNeuralFit::Synapse0x93bc9e8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*-0.173079);
}

Double_t AliCaloNeuralFit::Synapse0x93bca10() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*-0.356515);
}

Double_t AliCaloNeuralFit::Synapse0x93bca38() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*0.116333);
}

Double_t AliCaloNeuralFit::Synapse0x93bca60() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*-0.0647334);
}

Double_t AliCaloNeuralFit::Synapse0x93bca88() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*0.135181);
}

Double_t AliCaloNeuralFit::Synapse0x93bcca8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*-0.00703734);
}

Double_t AliCaloNeuralFit::Synapse0x93bccd0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*1.04245);
}

Double_t AliCaloNeuralFit::Synapse0x93bccf8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*0.70812);
}

Double_t AliCaloNeuralFit::Synapse0x93bcd20() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*-0.486738);
}

Double_t AliCaloNeuralFit::Synapse0x93bcd48() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*-0.623982);
}

Double_t AliCaloNeuralFit::Synapse0x93bcf68() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*0.0557636);
}

Double_t AliCaloNeuralFit::Synapse0x93bcf90() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*0.503374);
}

Double_t AliCaloNeuralFit::Synapse0x93bcfb8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*-0.0500216);
}

Double_t AliCaloNeuralFit::Synapse0x93bcfe0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*0.0561266);
}

Double_t AliCaloNeuralFit::Synapse0x93bd008() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*-0.403885);
}

Double_t AliCaloNeuralFit::Synapse0x93bd228() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*-0.0649038);
}

Double_t AliCaloNeuralFit::Synapse0x93bd250() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*-0.164437);
}

Double_t AliCaloNeuralFit::Synapse0x93bd300() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*-0.657412);
}

Double_t AliCaloNeuralFit::Synapse0x93bd328() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*0.175571);
}

Double_t AliCaloNeuralFit::Synapse0x93bd350() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*0.588896);
}

Double_t AliCaloNeuralFit::Synapse0x93bd528() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*-0.00361627);
}

Double_t AliCaloNeuralFit::Synapse0x93bd550() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*0.398821);
}

Double_t AliCaloNeuralFit::Synapse0x93bd578() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*0.38779);
}

Double_t AliCaloNeuralFit::Synapse0x93bd5a0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*0.341009);
}

Double_t AliCaloNeuralFit::Synapse0x93bd5c8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*0.290214);
}

Double_t AliCaloNeuralFit::Synapse0x93bd7e8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*-0.238609);
}

Double_t AliCaloNeuralFit::Synapse0x93bd810() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*-0.193207);
}

Double_t AliCaloNeuralFit::Synapse0x93bd838() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*0.304672);
}

Double_t AliCaloNeuralFit::Synapse0x93bd860() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*-0.327625);
}

Double_t AliCaloNeuralFit::Synapse0x93bd888() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*-0.268996);
}

Double_t AliCaloNeuralFit::Synapse0x93bdaa8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*0.725989);
}

Double_t AliCaloNeuralFit::Synapse0x93bdad0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*-0.217538);
}

Double_t AliCaloNeuralFit::Synapse0x93bdaf8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*-0.603975);
}

Double_t AliCaloNeuralFit::Synapse0x8ce7098() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*0.175658);
}

Double_t AliCaloNeuralFit::Synapse0x8ce70c0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*-0.262399);
}

Double_t AliCaloNeuralFit::Synapse0x93bde20() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*-1.08242);
}

Double_t AliCaloNeuralFit::Synapse0x93bde48() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*0.41456);
}

Double_t AliCaloNeuralFit::Synapse0x93bde70() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*0.117809);
}

Double_t AliCaloNeuralFit::Synapse0x93bde98() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*0.514938);
}

Double_t AliCaloNeuralFit::Synapse0x93bdec0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*-0.193678);
}

Double_t AliCaloNeuralFit::Synapse0x93be0e0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*0.580561);
}

Double_t AliCaloNeuralFit::Synapse0x93be108() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*0.610453);
}

Double_t AliCaloNeuralFit::Synapse0x93be130() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*0.4977);
}

Double_t AliCaloNeuralFit::Synapse0x93be158() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*-0.328315);
}

Double_t AliCaloNeuralFit::Synapse0x93be180() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*-0.729094);
}

Double_t AliCaloNeuralFit::Synapse0x93be3a0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb770()*0.172817);
}

Double_t AliCaloNeuralFit::Synapse0x93be3c8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8ceb900()*0.288833);
}

Double_t AliCaloNeuralFit::Synapse0x93be3f0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebad8()*0.0242409);
}

Double_t AliCaloNeuralFit::Synapse0x93be418() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x8cebcb0()*0.273568);
}

Double_t AliCaloNeuralFit::Synapse0x93be440() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc518()*0.261816);
}

Double_t AliCaloNeuralFit::Synapse0x93bd2c0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc838()*-0.332154);
}

Double_t AliCaloNeuralFit::Synapse0x93be468() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bcab0()*0.175458);
}

Double_t AliCaloNeuralFit::Synapse0x93be490() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bcd70()*0.211775);
}

Double_t AliCaloNeuralFit::Synapse0x93be4b8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bd030()*0.335443);
}

Double_t AliCaloNeuralFit::Synapse0x93be4e0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bd378()*0.341079);
}

Double_t AliCaloNeuralFit::Synapse0x93be508() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bd5f0()*-0.324381);
}

Double_t AliCaloNeuralFit::Synapse0x93be530() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bd8b0()*0.135666);
}

Double_t AliCaloNeuralFit::Synapse0x93be558() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bdc28()*-0.0213577);
}

Double_t AliCaloNeuralFit::Synapse0x93be580() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bdee8()*-0.598335);
}

Double_t AliCaloNeuralFit::Synapse0x93be5a8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93be1a8()*0.636813);
}

Double_t AliCaloNeuralFit::Synapse0x93be7d0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bc838()*0.173357);
}

Double_t AliCaloNeuralFit::Synapse0x93be7f8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bcab0()*-0.971616);
}

Double_t AliCaloNeuralFit::Synapse0x93be820() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bcd70()*-0.38099);
}

Double_t AliCaloNeuralFit::Synapse0x93be848() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bd030()*0.351755);
}

Double_t AliCaloNeuralFit::Synapse0x93be870() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bd378()*0.106307);
}

Double_t AliCaloNeuralFit::Synapse0x934a7c8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bd5f0()*0.118656);
}

Double_t AliCaloNeuralFit::Synapse0x93605e0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bd8b0()*-0.447119);
}

Double_t AliCaloNeuralFit::Synapse0x9360608() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bdc28()*0.677259);
}

Double_t AliCaloNeuralFit::Synapse0x8ce6fe8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93bdee8()*-0.246138);
}

Double_t AliCaloNeuralFit::Synapse0x8ce7010() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x93be1a8()*-0.117442);
}
