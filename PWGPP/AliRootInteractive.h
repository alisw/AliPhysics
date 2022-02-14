/*
 .L $AliPhysics_SRC/PWGPP/AliRootInteractive.h
*/
#ifndef AliRootInteractive_H
#define AliRootInteractive_H


namespace AliRootInteractive{
  void treeBokehDrawArray(const char * treeName, const char *query, const char *figureArray,  const char *widgets, const char *options, const char *htmlOut= nullptr);

  const char *importBokeh= "from Tools.aliTreePlayer import * \n"
                           "from anytree import * \n"
                           "from InteractiveDrawing.bokeh.bokehDrawSA import *\n"
                           "import numpy as np\n";
}



void AliRootInteractive::treeBokehDrawArray(const char * treeName, const char *query, const char *figureArray,  const char *widgets, const char *options, const char *htmlOut){
    //initPython();
    TString x=AliRootInteractive::importBokeh;
    x = x + "tree=ROOT.gROOT.GetGlobal(\""+treeName+"\")\n";
    if (htmlOut) x=x+"output_file('"+htmlOut +"')\n";
    x=x+"fig=bokehDrawSA.fromArray(tree ,'" + query + "'," + figureArray + ",'" + widgets + "'," + options+")";
    TPython::Exec(x);
}
#endif