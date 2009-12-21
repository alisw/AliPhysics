// Start the list analyser. No objects will be loaded into the list! To add e.g. AliTrackPoints you can do the following:
// - Run the macro "esd_tracks.C" 
// - Select some track you want to analyse as follows: Hold "shift" und right-click on the track (sometimes you have to hold the right mouse button). The menu pops up
// -> Select "ImportClustersFromIndex"
// -> To this for all tracks you want to analyse.
// - Run this macro (ana_list.C)
// - In the tab "eve" in the browser select the list analyser (called "Analysis objects" in the standard case)
// - Select the "list" tab in the editor of this object.
// - Click the button "start"
// - Select e.g. clusters by holding "ctrl"+"alt" (depending on your system, holding the "alt" only can also be fine) and left-clicking on the desired cluster
// Use the list analyser "as usual" (see class documentation)

void ana_list(TEveElement *cont = 0)
{
  AliEveListAnalyser * objList = new AliEveListAnalyser("Analysis objects");
  
  objList->SetTitle("Analysis objects (0)");

  gEve->AddElement(objList, cont);

  gEve->Redraw3D();
}
