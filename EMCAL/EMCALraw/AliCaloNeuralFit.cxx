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

#include "AliCaloNeuralFit.h"
#include <cmath>

/// \cond CLASSIMP
ClassImp( AliCaloNeuralFit ) ;
/// \endcond

///
/// Compute the neural network answer,
/// given the input values (taken from the signal TGraph)
//______________________________________________________________________________
Double_t AliCaloNeuralFit::Value
(int index, Double_t in0, Double_t in1, Double_t in2, Double_t in3, Double_t in4)
{
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

///
/// Input neuron.
/// Just return activation value externally setted.
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x8ceb770() const
{
  return fInput0;
}

///
/// Input neuron.
/// Just return activation value externally setted.
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x8ceb900() const
{
  return fInput1;
}

///
/// Input neuron.
/// Just return activation value externally setted.
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x8cebad8() const
{
  return fInput2;
}

///
/// Input neuron.
/// Just return activation value externally setted.
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x8cebcb0() const
{
  return fInput3;
}

///
/// Input neuron.
/// Just return activation value externally setted.
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bc518() const
{
  return fInput4;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bc838() const
{
  Double_t input = -0.097207;
  input += Synapse0x93bc9e8();
  input += Synapse0x93bca10();
  input += Synapse0x93bca38();
  input += Synapse0x93bca60();
  input += Synapse0x93bca88();
  
  return input;
}

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bc838() const
{
  Double_t input = Input0x93bc838();
  return (tanh(input) * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bcab0() const
{
  Double_t input = 0.0466086;
  input += Synapse0x93bcca8();
  input += Synapse0x93bccd0();
  input += Synapse0x93bccf8();
  input += Synapse0x93bcd20();
  input += Synapse0x93bcd48();
  return input;
}

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bcab0() const
{
  Double_t input = Input0x93bcab0();
  return (tanh(input) * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bcd70() const
{
  Double_t input = 0.617042;
  input += Synapse0x93bcf68();
  input += Synapse0x93bcf90();
  input += Synapse0x93bcfb8();
  input += Synapse0x93bcfe0();
  input += Synapse0x93bd008();
  
  return input;
}

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bcd70() const
{
  Double_t input = Input0x93bcd70();
  return (tanh(input) * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bd030() const
{
  Double_t input = -0.519529;
  input += Synapse0x93bd228();
  input += Synapse0x93bd250();
  input += Synapse0x93bd300();
  input += Synapse0x93bd328();
  input += Synapse0x93bd350();
  return input;
}

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bd030() const
{
  Double_t input = Input0x93bd030();
  return (tanh(input) * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bd378() const
{
  Double_t input = -0.405186;
  input += Synapse0x93bd528();
  input += Synapse0x93bd550();
  input += Synapse0x93bd578();
  input += Synapse0x93bd5a0();
  input += Synapse0x93bd5c8();
  return input;
}

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bd378() const
{
  Double_t input = Input0x93bd378();
  return (tanh(input) * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bd5f0() const
{
  
  Double_t input = -0.253465;
  input += Synapse0x93bd7e8();
  input += Synapse0x93bd810();
  input += Synapse0x93bd838();
  input += Synapse0x93bd860();
  input += Synapse0x93bd888();
  return input;
}

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bd5f0() const
{
  Double_t input = Input0x93bd5f0();
  return (tanh(input) * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bd8b0() const
{
  Double_t input = 0.0340672;
  input += Synapse0x93bdaa8();
  input += Synapse0x93bdad0();
  input += Synapse0x93bdaf8();
  input += Synapse0x8ce7098();
  input += Synapse0x8ce70c0();
  return input;
}

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bd8b0() const
{
  Double_t input = Input0x93bd8b0();
  return (tanh(input) * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bdc28() const
{
  Double_t input = 0.231599;
  input += Synapse0x93bde20();
  input += Synapse0x93bde48();
  input += Synapse0x93bde70();
  input += Synapse0x93bde98();
  input += Synapse0x93bdec0();
  return input;
}

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bdc28() const
{
  Double_t input = Input0x93bdc28();
  return (tanh(input) * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bdee8() const
{
  Double_t input = 0.56831;
  input += Synapse0x93be0e0();
  input += Synapse0x93be108();
  input += Synapse0x93be130();
  input += Synapse0x93be158();
  input += Synapse0x93be180();
  return input;
}

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bdee8() const
{
  Double_t input = Input0x93bdee8();
  return (tanh(input) * 1)+0;
}

//________________________________________________
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

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93be1a8() const
{
  Double_t input = Input0x93be1a8();
  return (tanh(input) * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93bc708() const
{
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

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93bc708() const
{
  Double_t input = Input0x93bc708();
  return (input * 1)+0;
}

///
/// Hidden/Output neuron
/// Compute the activation from linear combination of
/// all neurons going into this, each one times its synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Input0x93be5d0() const
{
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

///
/// Hidden/Output neuron
/// Return computed activation
//________________________________________________
Double_t AliCaloNeuralFit::Neuron0x93be5d0() const
{
  Double_t input = Input0x93be5d0();
  return (input * 1)+0;
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bc9e8() const
{
  return (Neuron0x8ceb770()*-0.173079);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bca10() const
{
  return (Neuron0x8ceb900()*-0.356515);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bca38() const
{
  return (Neuron0x8cebad8()*0.116333);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bca60() const
{
  return (Neuron0x8cebcb0()*-0.0647334);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bca88() const
{
  return (Neuron0x93bc518()*0.135181);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bcca8() const
{
  return (Neuron0x8ceb770()*-0.00703734);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bccd0() const
{
  return (Neuron0x8ceb900()*1.04245);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bccf8() const
{
  return (Neuron0x8cebad8()*0.70812);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bcd20() const
{
  return (Neuron0x8cebcb0()*-0.486738);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bcd48() const
{
  return (Neuron0x93bc518()*-0.623982);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bcf68() const
{
  return (Neuron0x8ceb770()*0.0557636);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bcf90() const
{
  return (Neuron0x8ceb900()*0.503374);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bcfb8() const
{
  return (Neuron0x8cebad8()*-0.0500216);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bcfe0() const
{
  return (Neuron0x8cebcb0()*0.0561266);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd008() const
{
  return (Neuron0x93bc518()*-0.403885);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd228() const
{
  return (Neuron0x8ceb770()*-0.0649038);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd250() const
{
  return (Neuron0x8ceb900()*-0.164437);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd300() const
{
  return (Neuron0x8cebad8()*-0.657412);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd328() const
{
  return (Neuron0x8cebcb0()*0.175571);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd350() const
{
  return (Neuron0x93bc518()*0.588896);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd528() const
{
  return (Neuron0x8ceb770()*-0.00361627);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd550() const
{
  return (Neuron0x8ceb900()*0.398821);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd578() const
{
  return (Neuron0x8cebad8()*0.38779);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd5a0() const
{
  return (Neuron0x8cebcb0()*0.341009);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd5c8() const
{
  return (Neuron0x93bc518()*0.290214);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd7e8() const
{
  return (Neuron0x8ceb770()*-0.238609);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd810() const
{
  return (Neuron0x8ceb900()*-0.193207);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd838() const
{
  return (Neuron0x8cebad8()*0.304672);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd860() const
{
  return (Neuron0x8cebcb0()*-0.327625);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd888() const
{
  return (Neuron0x93bc518()*-0.268996);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bdaa8() const
{
  return (Neuron0x8ceb770()*0.725989);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bdad0() const
{
  return (Neuron0x8ceb900()*-0.217538);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bdaf8() const
{
  return (Neuron0x8cebad8()*-0.603975);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x8ce7098() const
{
  return (Neuron0x8cebcb0()*0.175658);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x8ce70c0() const
{
  return (Neuron0x93bc518()*-0.262399);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bde20() const
{
  return (Neuron0x8ceb770()*-1.08242);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bde48() const
{
  return (Neuron0x8ceb900()*0.41456);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bde70() const
{
  return (Neuron0x8cebad8()*0.117809);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bde98() const
{
  return (Neuron0x8cebcb0()*0.514938);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bdec0() const
{
  return (Neuron0x93bc518()*-0.193678);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be0e0() const
{
  return (Neuron0x8ceb770()*0.580561);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be108() const
{
  return (Neuron0x8ceb900()*0.610453);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be130() const
{
  return (Neuron0x8cebad8()*0.4977);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be158() const
{
  return (Neuron0x8cebcb0()*-0.328315);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be180() const
{
  return (Neuron0x93bc518()*-0.729094);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be3a0() const
{
  return (Neuron0x8ceb770()*0.172817);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be3c8() const
{
  return (Neuron0x8ceb900()*0.288833);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be3f0() const
{
  return (Neuron0x8cebad8()*0.0242409);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be418() const
{
  return (Neuron0x8cebcb0()*0.273568);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be440() const
{
  return (Neuron0x93bc518()*0.261816);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93bd2c0() const
{
  return (Neuron0x93bc838()*-0.332154);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be468() const
{
  return (Neuron0x93bcab0()*0.175458);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be490() const
{
  return (Neuron0x93bcd70()*0.211775);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be4b8() const
{
  return (Neuron0x93bd030()*0.335443);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be4e0() const
{
  return (Neuron0x93bd378()*0.341079);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be508() const
{
  return (Neuron0x93bd5f0()*-0.324381);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be530() const
{
  return (Neuron0x93bd8b0()*0.135666);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be558() const
{
  return (Neuron0x93bdc28()*-0.0213577);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be580() const
{
  return (Neuron0x93bdee8()*-0.598335);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be5a8() const
{
  return (Neuron0x93be1a8()*0.636813);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be7d0() const
{
  return (Neuron0x93bc838()*0.173357);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be7f8() const
{
  return (Neuron0x93bcab0()*-0.971616);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be820() const
{
  return (Neuron0x93bcd70()*-0.38099);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be848() const
{
  return (Neuron0x93bd030()*0.351755);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93be870() const
{
  return (Neuron0x93bd378()*0.106307);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x934a7c8() const
{
  return (Neuron0x93bd5f0()*0.118656);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x93605e0() const
{
  return (Neuron0x93bd8b0()*-0.447119);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x9360608() const
{
  return (Neuron0x93bdc28()*0.677259);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x8ce6fe8() const
{
  return (Neuron0x93bdee8()*-0.246138);
}

///
/// Synaptic connection
/// Multiplies input times synaptic weight
//________________________________________________
Double_t AliCaloNeuralFit::Synapse0x8ce7010() const
{
  return (Neuron0x93be1a8()*-0.117442);
}
