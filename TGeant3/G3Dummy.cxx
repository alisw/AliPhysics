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
*/

//////////////////////////////////////////////////////
//  C++ dummy interface to Geant3 basic routines    //
//////////////////////////////////////////////////////

#include "AliCallf77.h"
 

extern "C" type_of_call {

  // This is only for ITS

  void type_of_call gsvolu_(){}
  void type_of_call gspos_(){}
  void type_of_call gsatt_(){}
  void type_of_call gsdvn_(){}
  void type_of_call gsposp_(){}
  void type_of_call atg_(){}
  void type_of_call sxsrot_(){}

  // All this for ZDC

  void type_of_call grndm_(){}
  void type_of_call gphits_(){}
  void type_of_call gfhits_() {}
  void type_of_call gsvert_(){}
  void type_of_call gskine_(){}
  void type_of_call gfpart_(){}
  void type_of_call lorenf_(){}
  void type_of_call gpvert_(){}
  void type_of_call gpkine_(){}
  void type_of_call rnorml_(){}
  void type_of_call gfpath_() {}
  void type_of_call uctoh_() {}
  void type_of_call glvolu_() {}
  void type_of_call gmtod_() {}
  void type_of_call gfkine_() {}
  void type_of_call vmod_() {}
  void type_of_call gsahit_() {}
  void type_of_call gschit_() {}
  void type_of_call gdtom_() {}
  void type_of_call rnpssn_() {}
  void type_of_call ucopy_() {}
}
