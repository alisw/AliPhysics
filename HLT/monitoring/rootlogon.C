/**************************************************************************
 * Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
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
 * \file rootlogon.C
 * \date 9 Oct 2011
 * \author Artur Szostak <artursz@iafrica.com>
 * \brief Convenient rootlogon.C macro to setup AliRoot when the binary starts for using the MonitorSandbox.C macro.
 */
rootlogon()
{
	TString includePath = "-I${ALICE_ROOT}/include ";
	gSystem->SetIncludePath(includePath.Data());
	
	gSystem->Load("libCDB.so");
	gSystem->Load("libHLTbase.so");
	gSystem->Load("libHLTrec.so");
	gSystem->Load("libAliHLTHOMER.so");
	gSystem->Load("libAliHLTTrigger.so");
}
