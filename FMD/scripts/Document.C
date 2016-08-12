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

/* $Id$ */

// Script to document the FMD code
/** @ingroup FMD_simple_script
 */
void
Document()
{
  gEnv->SetValue("Root.Html.SourceDir", "$(ALICE)/FMD");
  gEnv->SetValue("Root.Html.OutputDir", "$(ALICE)/FMD/html");

  gSystem->MakeDirectory("$(ALICE)/FMD/html");
  
  THtml* html = new THtml;
  html->MakeAll(kFALSE, "AliFMD*");
  html->Convert("$(ALICE)/FMD/Digitize.C", "Digitize", 
		"FMD/html/src");
  html->Convert("$(ALICE)/FMD/Reconstruct.C", "Reconstruct", 
		"FMD/html/src");
  html->Convert("$(ALICE)/FMD/Simulate.C", "Simulate", 
		"FMD/html/src");
  html->Convert("$(ALICE)/FMD/DrawFMD.C", "DrawFMD", 
		"FMD/html/src");
  html->Convert("$(ALICE)/FMD/ViewFMD.C", "ViewFMD", 
		"FMD/html/src");
  html->MakeIndex("AliFMD*");

  std::ofstream index("FMD/html/index.html");
  html->WriteHtmlHeader(index, "ALICE FMD Code - Index page");
  
  index << "<h1>ALICE FMD Code</h1>\n"
	<< "<ul>\n"
	<< "<li><a href=\"USER_Index.html\">Classes</a></li>\n"
        << "<li><a href=\"src/Digitize.C.html\">Digitize script</a></li>\n"
        << "<li><a href=\"src/Reconstruct.C.html\">Reconstruct script</a></li>\n"
        << "<li><a href=\"src/Simulate.C.html\">Simulate script</a></li>\n"
        << "<li><a href=\"src/DrawFMD.C.html\">DrawFMD script</a></li>\n"
        << "<li><a href=\"src/ViewFMD.C.html\">ViewFMD script</a></li>\n"
	<< "</ul>\n"
	<< std::endl;
  html->WriteHtmlFooter(index, "", "", "", "");
  index.close();
}

//
// EOF
//
