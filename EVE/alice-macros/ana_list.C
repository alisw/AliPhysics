// How to use the list analyser (and especially this macro):

// -- Primary selection (to add e.g. tracklets, tracks, etc.):
// - Load the objects you want to analyse with a sufficient macro (e.g. for tracks and tracklets you can use ana_list_load_tracks.C)
// - Run this macro (ana_list.C)
// - In the tab "eve" in the browser select the list analyser (called "Analysis objects" in the standard case)
// - Select the "list" tab in the editor of this object.
// - Click the button "start"
// - Select the objects you want to add by left-clicking on them in the viewer (or as well in the browser (left panel))
// - When you have finished adding the desired objects, click the button "stop"
// Use the list analyser "as usual" (see class documentation)

// -- Secondary selection (to add e.g. single clusters (of a TEvePointSet) or digits (of a TEveQuadSet)):
// To add e.g. AliTrackPoints or tracks you can do the following:
// - Run e.g. the macro "esd_tracks.C" 
// - Select some track you want to analyse as follows: Hold "shift" und right-click on the track (sometimes you have to hold the right mouse button). The menu pops up
// -> Select "ImportClustersFromIndex"
// -> Do this for all tracks you want to analyse.
// - Run this macro (ana_list.C)
// - In the tab "eve" in the browser select the list analyser (called "Analysis objects" in the standard case)
// - Select the "list" tab in the editor of this object.
// - Click the button "start"
// - Select e.g. clusters by holding "ctrl"+"alt" (depending on your system, holding the "alt" only can also be fine) and left-clicking on the desired cluster
// - When you have finished adding the desired objects, click the button "stop"
// Use the list analyser "as usual" (see class documentation)

void ana_list(TEveElement *cont = 0)
{
  AliEveListAnalyser * objList = new AliEveListAnalyser("Analysis objects");
  
  objList->SetTitle("Analysis objects (0)");

  gEve->AddElement(objList, cont);

  gEve->Redraw3D();
}
