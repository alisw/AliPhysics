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

/*
$Log$
Revision 1.1  2005/10/11 12:31:50  masera
Preprocessor classes for SPD (Paul Nilsson)

*/

///////////////////////////////////////////////////////////////////////////
// AliITSBadChannelsAuxSPD implementation by P. Nilsson 2005
// AUTHOR/CONTACT: Paul.Nilsson@cern.ch
//
// Auxiliary algorithms for the SPD
//
// This class contains converter methods and general algorithms for
// handling digit <-> channel convertions and methods for identifying
// changes (diff's) in arrays of bad channels, and for finding channels
// or digits in arrays of bad channels.
//
// The Diff algorithm can be used to check if there are any changes among
// the noisy channels. It returns two arrays, one with noisy channels
// are no longer visible (less likely to happen) and one with newly found
// noisy channels (more likely to happen).
//
// The Find algorithms looks for a given digit or channel in an array of
// known channels. It can be used by the clustering algorithms to check
// if a given digit that is about to be clustered, is in fact a known
// noisy channel. It should not be used in a normal cluster.
//
// Examples - Converters
//
// root [0] AliITSdigitSPD *d = new AliITSdigitSPD();
// root [1] d->SetCoord1(1);
// root [2] d->SetCoord2(2);
// root [3] d->SetSignal(1);
// root [4] AliITSBadChannelsAuxSPD *aux = new AliITSBadChannelsAuxSPD();
// root [5] AliITSChannelSPD *c = aux->CreateChannelFromDigit(d);
// root [6] cout << c->GetColumn() << endl;
// 1
// root [7] cout << c->GetRow() << endl;
// 2
// root [8] AliITSdigitSPD *d2 = aux->CreateDigitFromChannel(c);
// root [9] cout << d2->GetCoord1() << endl;
// 1
// root [10] cout << d2->GetCoord2() << endl;
// 2
// root [11] cout << d2->GetSignal() << endl;
// 1
// root [12] delete d2;
// root [13] delete d;
// root [14] delete c;
//
// The signal member of the digit is not a member of the channel class.
// It is artificially introduced by the CreateDigitFromChannel method and
// is per default set to 1.
//
// Modified by D. Elia, H. Tydesjo
// March 2006: Mixed up coordinates, bug fixed
//
///////////////////////////////////////////////////////////////////////////

#include "AliITSBadChannelsAuxSPD.h"

ClassImp(AliITSBadChannelsAuxSPD)

//__________________________________________________________________________
AliITSBadChannelsAuxSPD::AliITSBadChannelsAuxSPD(void)
{
  // Default constructor
}

//__________________________________________________________________________
Bool_t AliITSBadChannelsAuxSPD::Diff(TObjArray *&inputArray1, TObjArray *&inputArray2,
				     TObjArray *&outputArray1, TObjArray *&outputArray2) const
{
  // Make a diff between the input TObjArrays
  // 
  // Input: Two input TObjArrays of AliITSChannelSPD objects to be tested, two output TObjArrays corresponding
  //        to 
  // Output: Two output TObjArrays where outputArray1 contains the AliITSChannelSPD objects from inputArray1
  //         that are not in inputArray2, vice versa for outputArray2.
  // Return: kTRUE if the arrays differ

  Bool_t status = kFALSE;

  const Int_t kInputArray1Size = inputArray1->GetEntries();
  const Int_t kInputArray2Size = inputArray2->GetEntries();
  AliITSChannelSPD *ch1 = 0;
  AliITSChannelSPD *ch2 = 0;
  Bool_t found = kFALSE;

  Int_t i = 0, j = 0;

  // Pass 1
  Int_t lastFoundAtJ = -1;
  for (i = 0; i < kInputArray1Size; i++)
    {
      // Get the next channel from array 1
      ch1 = (AliITSChannelSPD *) inputArray1->At(i);

      // Is ch1 also in array 2?
      for (j = lastFoundAtJ + 1; j < kInputArray2Size; j++)
        {
          ch2 = (AliITSChannelSPD *) inputArray2->At(j);
          if (*ch1 == *ch2)
            {
              // Abort, go to next i
              found = kTRUE;
              lastFoundAtJ = j;
              break;
            }
        }

      // If ch1 was not found in array 2, store it
      if (!found)
        {
          outputArray1->Add(ch1);
        }
      else
        {
          found = kFALSE;
        }
    }

  // Pass 2
  lastFoundAtJ = -1;
  for (i = 0; i < kInputArray2Size; i++)
    {
      // Get the next channel from array 2
      ch2 = (AliITSChannelSPD *) inputArray2->At(i);

      // Is ch2 also in array 1?
      for (j = lastFoundAtJ + 1; j < kInputArray1Size; j++)
        {
          ch1 = (AliITSChannelSPD *) inputArray1->At(j);
          if (*ch1 == *ch2)
            {
              // Abort, go to next i
              found = kTRUE;
              lastFoundAtJ = j;
              break;
            }
        }

      // If ch1 was not found in array 1, store it
      if (!found)
        {
          outputArray2->Add(ch2);
        }
      else
        {
          found = kFALSE;
        }
    }

  if (outputArray1->GetEntries() > 0 || outputArray2->GetEntries() > 0) status = kTRUE;

  return status;
}

//__________________________________________________________________________
Bool_t AliITSBadChannelsAuxSPD::Find(AliITSChannelSPD *&channel, TObjArray *&array) const
{
  // Find the channel in the array
  //
  // Input: AliITSChannelSPD channel object, TObjArray of AliITSChannelSPD channel objects
  // Ouput: (none)
  // Return: kTRUE if channel is found in the array, kFALSE otherwise

  Bool_t status = kFALSE;

  // Loop over all channels in the array
  Int_t channelNr = 0;
  const Int_t kN = array->GetEntries();
  while (channelNr < kN)
    {
      if (*channel == *(AliITSChannelSPD *)array->At(channelNr))
	{
	  status = kTRUE;
	  break;
	}

      // Go to next channel
      channelNr++;
    }

  return status;
}

//__________________________________________________________________________
Bool_t AliITSBadChannelsAuxSPD::Find(AliITSdigitSPD *&digit, TObjArray *&array) const
{
  // Find the digit in the array
  //
  // WARNING: Using AliITSdigitSPD digits in this way is roughly 10% slower than to use AliITSChannelSPD channels
  //
  // Input: AliITSdigitSPD digit object, TObjArray of AliITSChannelSPD channel objects
  // Ouput: (none)
  // Return: kTRUE if digit is found in the array, kFALSE otherwise

  Bool_t status = kFALSE;

  AliITSChannelSPD *channel = 0;
  const Int_t kN = array->GetEntries();
  Int_t channelNr = 0;
  Int_t column = digit->GetCoord1();
  Int_t row = digit->GetCoord2();

  // Loop over all channels in the array
  while (channelNr < kN)
    {
      channel = (AliITSChannelSPD *)array->At(channelNr);
      if ( (channel->GetColumn() == column) && (channel->GetRow() == row) )
	{
	  status = kTRUE;
	  break;
	}

      // Go to next channel
      channelNr++;
    }

  return status;
}

//__________________________________________________________________________
AliITSdigitSPD* AliITSBadChannelsAuxSPD::CreateDigitFromChannel(const AliITSChannelSPD *&channel) const
{
  // Create a digit from a channel
  //
  // Input: AliITSChannelSPD object
  // Ouput: (none)
  // Return: AliITSdigitSPD object

  AliITSdigitSPD *digit = new AliITSdigitSPD();

  digit->SetCoord1(channel->GetColumn());
  digit->SetCoord2(channel->GetRow());
  digit->SetSignal(1);

  return digit;
}

//__________________________________________________________________________
AliITSChannelSPD* AliITSBadChannelsAuxSPD::CreateChannelFromDigit(const AliITSdigitSPD *&digit) const
{
  // Create a channel from a digit
  //
  // Input: AliITSdigitSPD object
  // Ouput: (none)
  // Return: AliITSChannelSPD object

  AliITSChannelSPD *channel = new AliITSChannelSPD();

  channel->SetColumn(digit->GetCoord1());
  channel->SetRow(digit->GetCoord2());

  return channel;
}

//__________________________________________________________________________
Int_t AliITSBadChannelsAuxSPD::GetNumberOfBadChannels(Int_t* &badChannelsArray, Int_t* &indexArray, Int_t size) const
{
  // Get the total number of bad channels

  Int_t n = 0;

  // Loop over all modules
  for (Int_t module = 0; module < size; module++)
    {
      // Get the module size (i.e. the number of bad channels)
      n += badChannelsArray[indexArray[module]];
    }

  return n;
}

//__________________________________________________________________________
Bool_t AliITSBadChannelsAuxSPD::CreateHTMLReport(char *name, Int_t* &badChannelsArray, Int_t* &indexArray,
						 Int_t indexArraySize, TString *buffer, Bool_t tags)
{
  // Create an HTML report from the bad channels array
  //
  // Input : file name, badChannelsArray, indexArray, size of indexArray (i.e. number of modules),
  //         tags boolean (if true, html tags will be added; if false, only formatted text will be created)
  // Output: TString (buffer) containing the html code/ASCII text
  // Return: kTRUE if a report has been created

  Bool_t status = kFALSE;

  Int_t totalNumberOfBadChannels = 0;
  Int_t numberOfModulesWithBadChannels = 0;

  if (tags)
    {
      buffer->Append("<html>");
      buffer->Append("<head><title>SPD bad channels</title></head>\n");
      buffer->Append("<body>\n");
    }

  buffer->Append("HTML report for file: ");
  buffer->Append(name);

  tags ? buffer->Append("<br>\n<br>\n") : buffer->Append("\n\n");

  char temp[10];

  // Loop over all modules
  for (Int_t module = 0; module < indexArraySize; module++)
    {
      // Get the start position of the data
      Int_t position = indexArray[module];

      // Get the module size (i.e. the number of bad channels)
      Int_t size = badChannelsArray[position++];

      // Only continue if there are bad channels in this module
      if (size > 0)
	{
	  // There are bad channels in this file
	  status = kTRUE;
	  numberOfModulesWithBadChannels++;
	  totalNumberOfBadChannels += size;

	  // Create report
	  buffer->Append("SPD module = ");
	  sprintf(temp,"%d",module);
	  buffer->Append(temp);
	  buffer->Append("<br>\n");
	  buffer->Append("Number of bad channels = ");
	  sprintf(temp,"%d",size);
	  buffer->Append(temp);

	  tags ? buffer->Append("<br>\n") : buffer->Append("\n");

	  buffer->Append("(column, row) = ");

	  // Get all bad channels
	  Int_t i = 0;
	  while (i < size)
	    {
	      // Create and add the current channel
	      buffer->Append("(");
	      sprintf(temp,"%d",badChannelsArray[position++]);
	      buffer->Append(temp);
	      buffer->Append(", ");
	      sprintf(temp,"%d",badChannelsArray[position++]);
	      buffer->Append(temp);
	      buffer->Append(")");

	      if (i < size - 1)
		{
		  buffer->Append(", ");
		}
	      else
		{
		  tags ? buffer->Append("<br>\n") : buffer->Append("\n");
		}

	      // Go to next bad channel
	      i++;
	    }

	  tags ? buffer->Append("<br>\n") : buffer->Append("\n");

	} // if size > 0
    } // end loop over modules

  if (!status)
    {
      buffer->Append("(Data does not contain any known bad channels)");
    }

  tags ? buffer->Append("<br>\n") : buffer->Append("\n");

  buffer->Append("Total number of bad channels = ");
  sprintf(temp,"%d",totalNumberOfBadChannels);
  buffer->Append(temp);

  tags ? buffer->Append("<br>\n") : buffer->Append("\n");

  buffer->Append("Number of modules with bad channels = ");
  sprintf(temp,"%d",numberOfModulesWithBadChannels);
  buffer->Append(temp);

  tags ? buffer->Append("<br>\n") : buffer->Append("\n");

  buffer->Append("Number of modules = ");
  sprintf(temp,"%d",indexArraySize);
  buffer->Append(temp);

  if (tags)
    {
      buffer->Append("<br>\n");
      buffer->Append("</body>\n");
      buffer->Append("</html>");
    }
  else
    {
      buffer->Append("\n");
    }

  return status;
}
