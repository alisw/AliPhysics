/* $Id$ */

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

/**
 * >> Flat structure representing a TPC seed <<
 *
 * To be used in the online and offline calibration schema.
 *
 *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli
 *
 **************************************************************************/

#include "AliFlatTPCseed.h"
#include "Riostream.h"

AliFlatTPCseed::AliFlatTPCseed()
  :
  fContentSize(0)
{
  // constructor
  fContent[0]=0;
}

void AliFlatTPCseed::Reset()
{
  // Reset
}

void AliFlatTPCseed::SetFromTPCseed( const AliTPCseed *p )
{
  // initialise from AliTPCseed

}
