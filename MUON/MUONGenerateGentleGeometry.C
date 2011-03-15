/// \ingroup macros
/// \file MUONGenerateGentleGeometry.C
/// \brief Macro for generating a simplified geometry used by the
///  event display alieve. The generated file gentle_geo_muon.root
///  must be placed in EVE/alice-data/. The input is the sensitive
///  volumes list file svmap.dat and the geometry file geometry.root
///
/// To be run from aliroot:
///
/// .x MUONGenerateGentleGeometry.C
///
///
/// \author: M. Tadel, CERN and B. Vulpescu LPC, Clermont-Ferrand

void MUONGenerateGentleGeometry() {

  gSystem->Load("libGeom");

  TEveManager::Create();

  TGeoManager::Import("geometry.root");
  TGeoNode* tnode = gGeoManager->GetTopNode();
  TEveGeoTopNode* eve_tnode = new TEveGeoTopNode(gGeoManager, tnode);
  tnode->SetVisibility(kFALSE);
  eve_tnode->SetVisLevel(0);

  gEve->AddGlobalElement(eve_tnode);

  TString path;
  TObjArray *list;
  Int_t depth;

  Char_t line[256];
  ifstream in("data/svmap.dat", ios::in);

  while (!in.eof()) {
    
    in >> line;
    
    path = TString(line);
    if (!path.Contains("ALIC")) continue;

    list = path.Tokenize("/");
    depth = list->GetEntries();
    AddNodes(tnode,eve_tnode,depth,depth,list);
  }

  eve_tnode->SaveExtract("gentle_geo_muon.root", "Gentle MUON", kTRUE);

}

//_____________________________________________________________________________
void AddNodes(TGeoNode *node, TEveGeoNode *parent, Int_t depth, Int_t depthmax, TObjArray *list)
{  
  if (--depth <= 0)
    return;
  
  TObjString *nname = (TObjString*)list->At(depthmax-depth);
  TString sname = nname->GetString();

  TObjArray *nlist = node->GetVolume()->GetNodes();
  if (nlist == 0x0) return;
  Int_t     nNodes = nlist->GetEntries();
 
 for (Int_t in = 0; in < nNodes; in++)
 {
    TGeoNode *node = (TGeoNode*) nlist->At(in);
    TEveGeoNode *son;
    if (strcmp(node->GetName(),sname.Data()) == 0)
    {
      son = dynamic_cast<TEveGeoNode*>(parent->FindChild(nname->GetName()));
      if (!son)
      {
	son = new TEveGeoNode(node);
	parent->AddElement(son);
      }
    } else {
      continue;
    }
    AddNodes(node,son, depth, depthmax, list);
  }

}
