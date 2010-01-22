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
         return Neuron0x945d4f0();
     case 1:
         return Neuron0x945f490();
     default:
         return 0.;
   }
}

Double_t AliCaloNeuralFit::Neuron0x945cbe8() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput0;
}

Double_t AliCaloNeuralFit::Neuron0x945cd78() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput1;
}

Double_t AliCaloNeuralFit::Neuron0x945cf50() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput2;
}

Double_t AliCaloNeuralFit::Neuron0x945d128() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput3;
}

Double_t AliCaloNeuralFit::Neuron0x945d300() const
{
//
// Input neuron.
// Just return activation value externally setted.
//

   return fInput4;
}

Double_t AliCaloNeuralFit::Input0x945d620() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//

   Double_t input = -0.508174;
   input += Synapse0x943edb8();
   input += Synapse0x945d7d0();
   input += Synapse0x945d7f8();
   input += Synapse0x945d820();
   input += Synapse0x945d848();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945d620() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945d620();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945d870() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.29145;
   input += Synapse0x945da68();
   input += Synapse0x945da90();
   input += Synapse0x945dab8();
   input += Synapse0x945dae0();
   input += Synapse0x945db08();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945d870() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945d870();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945db30() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.132489;
   input += Synapse0x945dd28();
   input += Synapse0x945dd50();
   input += Synapse0x945dd78();
   input += Synapse0x945dda0();
   input += Synapse0x945ddc8();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945db30() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945db30();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945ddf0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -1.12891;
   input += Synapse0x945dfe8();
   input += Synapse0x945e010();
   input += Synapse0x945e0c0();
   input += Synapse0x945e0e8();
   input += Synapse0x945e110();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945ddf0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945ddf0();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945e138() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.576896;
   input += Synapse0x945e2e8();
   input += Synapse0x945e310();
   input += Synapse0x945e338();
   input += Synapse0x945e360();
   input += Synapse0x945e388();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945e138() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945e138();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945e3b0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.654194;
   input += Synapse0x945e5a8();
   input += Synapse0x945e5d0();
   input += Synapse0x945e5f8();
   input += Synapse0x945e620();
   input += Synapse0x945e648();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945e3b0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945e3b0();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945e670() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.356397;
   input += Synapse0x945e868();
   input += Synapse0x945e890();
   input += Synapse0x945e8b8();
   input += Synapse0x945e038();
   input += Synapse0x945e060();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945e670() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945e670();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945e9e8() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.798487;
   input += Synapse0x945ebe0();
   input += Synapse0x945ec08();
   input += Synapse0x945ec30();
   input += Synapse0x945ec58();
   input += Synapse0x945ec80();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945e9e8() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945e9e8();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945eca8() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.934985;
   input += Synapse0x945eea0();
   input += Synapse0x945eec8();
   input += Synapse0x945eef0();
   input += Synapse0x945ef18();
   input += Synapse0x945ef40();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945eca8() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945eca8();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945ef68() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.457775;
   input += Synapse0x945f160();
   input += Synapse0x945f188();
   input += Synapse0x945f1b0();
   input += Synapse0x945f1d8();
   input += Synapse0x945f200();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945ef68() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945ef68();
   return (tanh(input) * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945d4f0() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = 0.849942;
   input += Synapse0x945f300();
   input += Synapse0x945f328();
   input += Synapse0x945f350();
   input += Synapse0x945f378();
   input += Synapse0x945f3a0();
   input += Synapse0x945f3c8();
   input += Synapse0x945f3f0();
   input += Synapse0x945f418();
   input += Synapse0x945f440();
   input += Synapse0x945f468();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945d4f0() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945d4f0();
   return (input * 1)+0;
}

Double_t AliCaloNeuralFit::Input0x945f490() const
{
//
// Hidden/Output neuron
// Compute the activation from linear combination of
// all neurons going into this, each one times its synaptic weight
//
   Double_t input = -0.147416;
   input += Synapse0x945f690();
   input += Synapse0x945f6b8();
   input += Synapse0x945f6e0();
   input += Synapse0x945f708();
   input += Synapse0x945f730();
   input += Synapse0x936a1f0();
   input += Synapse0x943ee18();
   input += Synapse0x945cb70();
   input += Synapse0x945cb98();
   input += Synapse0x945cbc0();
   return input;
}

Double_t AliCaloNeuralFit::Neuron0x945f490() const
{
//
// Hidden/Output neuron
// Return computed activation
//
   Double_t input = Input0x945f490();
   return (input * 1)+0;
}

Double_t AliCaloNeuralFit::Synapse0x943edb8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*-0.104546);
}

Double_t AliCaloNeuralFit::Synapse0x945d7d0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*-0.0905177);
}

Double_t AliCaloNeuralFit::Synapse0x945d7f8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*-0.143637);
}

Double_t AliCaloNeuralFit::Synapse0x945d820() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*-0.413064);
}

Double_t AliCaloNeuralFit::Synapse0x945d848() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*0.883744);
}

Double_t AliCaloNeuralFit::Synapse0x945da68() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*-1.26724);
}

Double_t AliCaloNeuralFit::Synapse0x945da90() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*-0.14136);
}

Double_t AliCaloNeuralFit::Synapse0x945dab8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*0.27187);
}

Double_t AliCaloNeuralFit::Synapse0x945dae0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*0.563302);
}

Double_t AliCaloNeuralFit::Synapse0x945db08() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*1.38006);
}

Double_t AliCaloNeuralFit::Synapse0x945dd28() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*-0.235737);
}

Double_t AliCaloNeuralFit::Synapse0x945dd50() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*0.715314);
}

Double_t AliCaloNeuralFit::Synapse0x945dd78() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*-0.893506);
}

Double_t AliCaloNeuralFit::Synapse0x945dda0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*1.66689);
}

Double_t AliCaloNeuralFit::Synapse0x945ddc8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*0.433463);
}

Double_t AliCaloNeuralFit::Synapse0x945dfe8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*0.198835);
}

Double_t AliCaloNeuralFit::Synapse0x945e010() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*1.67429);
}

Double_t AliCaloNeuralFit::Synapse0x945e0c0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*-1.19328);
}

Double_t AliCaloNeuralFit::Synapse0x945e0e8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*2.5465);
}

Double_t AliCaloNeuralFit::Synapse0x945e110() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*0.153072);
}

Double_t AliCaloNeuralFit::Synapse0x945e2e8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*0.0815823);
}

Double_t AliCaloNeuralFit::Synapse0x945e310() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*0.0316826);
}

Double_t AliCaloNeuralFit::Synapse0x945e338() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*0.617448);
}

Double_t AliCaloNeuralFit::Synapse0x945e360() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*-0.749993);
}

Double_t AliCaloNeuralFit::Synapse0x945e388() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*-0.980764);
}

Double_t AliCaloNeuralFit::Synapse0x945e5a8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*-0.453657);
}

Double_t AliCaloNeuralFit::Synapse0x945e5d0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*0.146578);
}

Double_t AliCaloNeuralFit::Synapse0x945e5f8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*0.123041);
}

Double_t AliCaloNeuralFit::Synapse0x945e620() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*0.189871);
}

Double_t AliCaloNeuralFit::Synapse0x945e648() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*-0.153873);
}

Double_t AliCaloNeuralFit::Synapse0x945e868() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*-0.0569668);
}

Double_t AliCaloNeuralFit::Synapse0x945e890() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*-0.0208438);
}

Double_t AliCaloNeuralFit::Synapse0x945e8b8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*-0.662575);
}

Double_t AliCaloNeuralFit::Synapse0x945e038() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*-0.308952);
}

Double_t AliCaloNeuralFit::Synapse0x945e060() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*-0.0591419);
}

Double_t AliCaloNeuralFit::Synapse0x945ebe0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*0.203333);
}

Double_t AliCaloNeuralFit::Synapse0x945ec08() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*-0.210458);
}

Double_t AliCaloNeuralFit::Synapse0x945ec30() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*-0.46208);
}

Double_t AliCaloNeuralFit::Synapse0x945ec58() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*-0.213809);
}

Double_t AliCaloNeuralFit::Synapse0x945ec80() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*0.652572);
}

Double_t AliCaloNeuralFit::Synapse0x945eea0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*0.53005);
}

Double_t AliCaloNeuralFit::Synapse0x945eec8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*1.97055);
}

Double_t AliCaloNeuralFit::Synapse0x945eef0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*-0.934772);
}

Double_t AliCaloNeuralFit::Synapse0x945ef18() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*-0.253289);
}

Double_t AliCaloNeuralFit::Synapse0x945ef40() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*-0.190109);
}

Double_t AliCaloNeuralFit::Synapse0x945f160() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cbe8()*0.111492);
}

Double_t AliCaloNeuralFit::Synapse0x945f188() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cd78()*0.928076);
}

Double_t AliCaloNeuralFit::Synapse0x945f1b0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945cf50()*0.178153);
}

Double_t AliCaloNeuralFit::Synapse0x945f1d8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d128()*-0.750558);
}

Double_t AliCaloNeuralFit::Synapse0x945f200() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d300()*-1.40984);
}

Double_t AliCaloNeuralFit::Synapse0x945f300() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d620()*-0.838377);
}

Double_t AliCaloNeuralFit::Synapse0x945f328() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d870()*0.191143);
}

Double_t AliCaloNeuralFit::Synapse0x945f350() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945db30()*-0.453988);
}

Double_t AliCaloNeuralFit::Synapse0x945f378() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945ddf0()*-0.520562);
}

Double_t AliCaloNeuralFit::Synapse0x945f3a0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945e138()*-0.995398);
}

Double_t AliCaloNeuralFit::Synapse0x945f3c8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945e3b0()*-0.114216);
}

Double_t AliCaloNeuralFit::Synapse0x945f3f0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945e670()*-0.72899);
}

Double_t AliCaloNeuralFit::Synapse0x945f418() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945e9e8()*-0.453087);
}

Double_t AliCaloNeuralFit::Synapse0x945f440() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945eca8()*0.0891431);
}

Double_t AliCaloNeuralFit::Synapse0x945f468() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945ef68()*0.679937);
}

Double_t AliCaloNeuralFit::Synapse0x945f690() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d620()*0.806704);
}

Double_t AliCaloNeuralFit::Synapse0x945f6b8() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945d870()*-1.27447);
}

Double_t AliCaloNeuralFit::Synapse0x945f6e0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945db30()*1.0306);
}

Double_t AliCaloNeuralFit::Synapse0x945f708() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945ddf0()*2.09234);
}

Double_t AliCaloNeuralFit::Synapse0x945f730() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945e138()*0.0643316);
}

Double_t AliCaloNeuralFit::Synapse0x936a1f0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945e3b0()*-0.204933);
}

Double_t AliCaloNeuralFit::Synapse0x943ee18() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945e670()*0.423604);
}

Double_t AliCaloNeuralFit::Synapse0x945cb70() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945e9e8()*1.00527);
}

Double_t AliCaloNeuralFit::Synapse0x945cb98() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945eca8()*-1.54485);
}

Double_t AliCaloNeuralFit::Synapse0x945cbc0() const
{
//
// Synaptic connection
// Multiplies input times synaptic weight
//
   return (Neuron0x945ef68()*0.540381);
}

