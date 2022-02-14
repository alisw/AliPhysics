#include "TXMLEngine.h"
#include <TString.h>
#include <TClonesArray.h>
#include <Riostream.h>
#include <vector>

#include "AliRDHFBDT.h"

//
// Macro to convert TMVA BDT model(.xml) to AliRDHFBDT object
// 
// Author: M. Cai cai.mengke@cern.ch
//

void ReadWeight(TXMLEngine &xml, XMLNodePointer_t WeightNode, AliRDHFBDT *objBDT);
void ReadDTNodes(TXMLEngine &xml, XMLNodePointer_t MNode_xml, TClonesArray &nodes);
void ReadAllAtt(TXMLEngine &xml, XMLNodePointer_t node);
void ReadInfo(TXMLEngine &xml, XMLNodePointer_t InfoNode, AliRDHFBDT *objBDT);
void ReadVarList(TXMLEngine &xml, XMLNodePointer_t VarListNode, AliRDHFBDT *objBDT);
void ReadDTNodes(TXMLEngine &xml, XMLNodePointer_t MNode_xml, TClonesArray &nodes);
void LinkDTNodes(TXMLEngine &xml, XMLNodePointer_t MNode_xml, TClonesArray &nodes);
AliRDHFDTNode *ConvertDTNode(TXMLEngine &xml, XMLNodePointer_t Node_xml);

TString ReadAttStr(TXMLEngine &xml, XMLNodePointer_t Node, const char* attName);
Int_t ReadAttInt(TXMLEngine &xml, XMLNodePointer_t Node, const char* attName);
Float_t ReadAttFlt(TXMLEngine &xml, XMLNodePointer_t Node, const char* attName);

// Main Function
void ReadTMVABDTxml(const char* infilename, const char* outfilename, const char* bdtname = "BDT", const char* bdttitle = "Example BDT")
{
   TXMLEngine XML;
   XMLDocPointer_t xmlfile = XML.ParseFile(infilename);
   if (!xmlfile) {std::cout << "ERROR: Failed to open input weightxml file!"<<endl; return;}

   XMLNodePointer_t BranchMain = XML.DocGetRootElement(xmlfile); // 'MethodSetup'
   // Find Node
   XMLNodePointer_t BranchInfo, BranchSetting, BranchVar, BranchBDT;
   XMLNodePointer_t branch = XML.GetChild(BranchMain);
   while (branch != 0) {
	  TString branchname = XML.GetNodeName(branch);
      if(branchname=="GeneralInfo") BranchInfo = branch;
      if(branchname=="Options") 	BranchSetting = branch;
      if(branchname=="Variables")	BranchVar = branch;
      if(branchname=="Weights") 	BranchBDT = branch;
      branch = XML.GetNext(branch);
   }
   if(!BranchInfo) 		std::cout<<"WARNING: General info cannot be read..."<<endl;
   if(!BranchSetting) 	std::cout<<"WARNING: Training setting cannot be read..."<<endl;
   if(!BranchVar) 		std::cout<<"ERROR: Variable list cannot be read..."<<endl;
   if(!BranchBDT) 		std::cout<<"ERROR: BDT weight cannot be read..."<<endl;
   // BranchInfo/BranchSetting/BranchVar/BranchBDT Ready
   
   TFile *out = new TFile(outfilename,"recreate");
   AliRDHFBDT *bdt = new AliRDHFBDT();
   bdt->SetName(bdtname);
   bdt->SetTitle(bdttitle);
   bdt->AddDescLine(bdttitle);
   
   //~ ReadAllAtt(XML,BranchInfo);
   //~ ReadAllAtt(XML,BranchSetting);
   //~ ReadAllAtt(XML,BranchVar);
   //~ ReadAllAtt(XML,BranchBDT);
   
   ReadInfo(XML,BranchInfo,bdt);
   ReadInfo(XML,BranchSetting,bdt);
   ReadVarList(XML,BranchVar,bdt);
   ReadWeight(XML,BranchBDT,bdt);
   
   //~ bdt->Description();
   //~ bdt->PrintAll();
   std::cout<<endl;
   bdt->Write();
   out->Close();
   //~ ReadNode(XML, mainnode);

   XML.FreeDoc(xmlfile);
}


void ReadWeight(TXMLEngine &xml, XMLNodePointer_t WeightNode, AliRDHFBDT *objBDT)
{
   TString infoName(xml.GetNodeName(WeightNode));
   if(!(infoName=="Weights")) return; // In case of wrong node input

   Int_t NTrees = ReadAttInt(xml,WeightNode,"NTrees");
   
   XMLNodePointer_t DT_xml = xml.GetChild(WeightNode);
   while(DT_xml){
	   TString DTName(xml.GetNodeName(DT_xml));
	   if(!(DTName=="BinaryTree")) continue; // In case of wrong node input
	   
	   Int_t 	iTree = ReadAttInt(xml,DT_xml,"itree");
	   Double_t boostw = (Double_t)ReadAttFlt(xml,DT_xml,"boostWeight");
	   AliRDHFDecisionTree *DT;
	   TClonesArray DTNodes("AliRDHFDTNode",20);
	   // Add root node
	   XMLNodePointer_t rootNode_xml = xml.GetChild(DT_xml);
	   xml.NewIntAttr(rootNode_xml,"NodeInd",0);
	   new(DTNodes[0]) AliRDHFDTNode(*ConvertDTNode(xml,rootNode_xml) );
	   // Add other nodes
	   ReadDTNodes(xml, rootNode_xml, DTNodes);
	   // Link all nodes
	   LinkDTNodes(xml, rootNode_xml, DTNodes);
	   // Create vHFBDT
	   //~ std::cout<< "INFO: Adding root node..."<<endl;
	   AliRDHFDTNode *rootNode = (AliRDHFDTNode*)DTNodes.UncheckedAt(0)->Clone();
	   DT = new AliRDHFDecisionTree(rootNode);
	   //~ rootNode->PrintAll();
	   for(Int_t i=1;i<DTNodes.GetEntriesFast();i++){
		   AliRDHFDTNode *tmpNode = (AliRDHFDTNode*)DTNodes.UncheckedAt(i)->Clone();
		   //~ std::cout<<"INFO: Adding "<<i<<"-th node..."<<endl;
		   //~ tmpNode->PrintAll();
		   DT->AddNode(tmpNode);
	   }
	   std::cout<<"INFO: Boostweight: "<<boostw<<" ..."<<endl;
	   //~ DT->PrintAll();
	   //~ if(DT->GetNNodes()<=2) {DT_xml = xml.GetNext(DT_xml);continue;}
	   objBDT->AddDecisionTree(DT,boostw);
	   DT_xml = xml.GetNext(DT_xml);
   }
   
   return;
}

void ReadDTNodes(TXMLEngine &xml, XMLNodePointer_t MNode_xml, TClonesArray &nodes)
{
   TString DTName(xml.GetNodeName(MNode_xml));
   if(!(DTName=="Node")) return 0x0; // In case of wrong node input
   
   XMLNodePointer_t dauNode_xml = xml.GetChild(MNode_xml);
   // First create new attribute 'index' for each node
   //~ xml.NewIntAttr(MNode_xml,"NodeInd",nodes.GetEntriesFast());
   while(dauNode_xml){	   
	   xml.NewIntAttr(dauNode_xml,"NodeInd",nodes.GetEntriesFast());
	   new(nodes[nodes.GetEntriesFast()]) AliRDHFDTNode(*ConvertDTNode(xml,dauNode_xml));
	   //~ std::cout<<"INFO: Add node at "<<nodes.GetEntriesFast()-1<<endl;
	   //~ ConvertDTNode(xml,dauNode_xml)->PrintAll();
	   // Recursively read all nodes(DFT)
	   ReadDTNodes(xml, dauNode_xml, nodes);	   
	   dauNode_xml = xml.GetNext(dauNode_xml);
   }
   return;
}

void LinkDTNodes(TXMLEngine &xml, XMLNodePointer_t MNode_xml, TClonesArray &nodes)
{
   TString DTName(xml.GetNodeName(MNode_xml));
   if(!(DTName=="Node")) return 0x0; // In case of wrong node input
   
   XMLNodePointer_t dauNode_xml = xml.GetChild(MNode_xml);
   while(dauNode_xml){
	   TString nodePos = ReadAttStr(xml, dauNode_xml, "pos");
	   Int_t nodeType = ReadAttInt(xml, dauNode_xml, "nType");
	   Int_t Ind_parent = ReadAttInt(xml, MNode_xml, "NodeInd");
	   Int_t Ind_this = ReadAttInt(xml, dauNode_xml, "NodeInd");
	   AliRDHFDTNode *Node_parent = (AliRDHFDTNode*)nodes.UncheckedAt(Ind_parent);
	   AliRDHFDTNode *Node_this = (AliRDHFDTNode*)nodes.UncheckedAt(Ind_this);
	   
	   Node_this->SetMNodeInd(Ind_parent);
	   if(nodePos=="r") Node_parent->SetRNodeInd(Ind_this);
	   else if (nodePos=="l") Node_parent->SetLNodeInd(Ind_this);
	   
	   if(nodeType) {Node_this->SetLNodeInd(-1); Node_this->SetRNodeInd(-1);}
	   else LinkDTNodes(xml, dauNode_xml, nodes);

	   dauNode_xml = xml.GetNext(dauNode_xml);
   }
   return;
}

void ReadAllAtt(TXMLEngine &xml, XMLNodePointer_t node)
{
   TString nodename = xml.GetNodeName(node);
   std::cout<<nodename<<":"<<endl;
   XMLAttrPointer_t attr = xml.GetFirstAttr(node);
   while(attr){
	  TString attName(xml.GetAttrName(attr));
	   printf(" attr.%s = %s\n", attName.Data(), xml.GetAttrValue(attr));
	   attr = xml.GetNextAttr(attr);
   }
   return;
}

void ReadInfo(TXMLEngine &xml, XMLNodePointer_t InfoNode, AliRDHFBDT *objBDT)
{
   TString infoName(xml.GetNodeName(InfoNode));
   if(!(infoName=="GeneralInfo"||infoName=="Options")) return; // In case of wrong node input
   
   XMLNodePointer_t dauNode = xml.GetChild(InfoNode);
   while(dauNode){
	  TString insert_this_line;
	  if(infoName=="GeneralInfo"){
		  TString nodeName(xml.GetAttr(dauNode,"name"));
		  if(nodeName=="Host"||nodeName=="Dir") {dauNode = xml.GetNext(dauNode);continue;}
		  insert_this_line = nodeName;
		  insert_this_line.Append(":  ");
		  insert_this_line.Append(xml.GetAttr(dauNode,"value"));
	  }
	  else if(infoName=="Options"){
		  TString nodeName(xml.GetAttr(dauNode,"name"));
		  //~ if(!(nodeName=="NTrees"||nodeName=="MaxDepth"||nodeName=="MinNodeSize"||nodeName=="nCuts"||nodeName=="BoostType"
		  //~ ||nodeName=="AdaBoostR2Loss"||nodeName=="UseBaggedBoost"||nodeName=="Shrinkage"||nodeName=="AdaBoostBeta"
		  //~ ||nodeName=="NodePurityLimit"||nodeName=="SeparationType"||nodeName=="RegressionLossFunctionBDTG"||nodeName=="HuberQuantile")) 
		  //~ {dauNode = xml.GetNext(dauNode);continue;}
		  insert_this_line = nodeName;
		  insert_this_line.Append(" = ");
		  insert_this_line.Append(xml.GetNodeContent(dauNode));
	  }
	  objBDT->AddDescLine(insert_this_line);
	  //~ printf(" attr.%s = %s\n", attName.Data(), xml.GetAttrValue(attr));
	  dauNode = xml.GetNext(dauNode);
   }
   return;   
}

void ReadVarList(TXMLEngine &xml, XMLNodePointer_t VarListNode, AliRDHFBDT *objBDT)
{
   TString varNodeName(xml.GetNodeName(VarListNode));
   if(!(varNodeName=="Variables")) return; // In case of wrong node input
   objBDT->SetNVars(ReadAttInt(xml,VarListNode,"NVar"));
   
   // Get variable names
   XMLNodePointer_t dauNode = xml.GetChild(VarListNode);
   while(dauNode){
	  Int_t ind(-1); TString varname;
	  ind = ReadAttInt(xml,dauNode,"VarIndex");
	  varname = ReadAttStr(xml,dauNode,"Expression");
	  if(ind>=0) objBDT->SetVarName(ind,varname);
	  dauNode = xml.GetNext(dauNode);
   }
   return; 
}

TString ReadAttStr(TXMLEngine &xml, XMLNodePointer_t Node, const char* attName)
{
	TString str_Att(xml.GetAttr(Node, attName));
	if(!str_Att) std::cout<<"ERROR: Cannot find the attribute named '"<<attName<<"'."<<endl;
	return str_Att;
}

Int_t ReadAttInt(TXMLEngine &xml, XMLNodePointer_t Node, const char* attName)
{
	TString str_Att(xml.GetAttr(Node, attName));
	if(!str_Att){
		std::cout<<"ERROR: Cannot find the attribute named '"<<attName<<"'."<<endl;
		return 0;
	}
	else return str_Att.Atoi();
}

Float_t ReadAttFlt(TXMLEngine &xml, XMLNodePointer_t Node, const char* attName)
{
	TString str_Att(xml.GetAttr(Node, attName));
	if(!str_Att){
		std::cout<<"ERROR: Cannot find the attribute named '"<<attName<<"'."<<endl;
		return 0;
	}
	else return str_Att.Atof();
}

AliRDHFDTNode *ConvertDTNode(TXMLEngine &xml, XMLNodePointer_t Node_xml)
{
   TString nodeName(xml.GetNodeName(Node_xml));
   if(!(nodeName=="Node")) return 0x0; // In case of wrong node input
   
   TString nPos;
   Int_t nType, nIvar, nCType;
   Double_t nCut, nRes, nPur;
   nPos		= ReadAttStr(xml, Node_xml, "pos");
   nType 	= ReadAttInt(xml, Node_xml, "nType");
   nIvar 	= ReadAttInt(xml, Node_xml, "IVar");
   nCut 	= (Double_t)ReadAttFlt(xml, Node_xml, "Cut");
   nCType 	= ReadAttInt(xml, Node_xml, "cType");
   nRes 	= (Double_t)ReadAttFlt(xml, Node_xml, "res");
   nPur 	= (Double_t)ReadAttFlt(xml, Node_xml, "purity");
   AliRDHFDTNode *node = new AliRDHFDTNode(-1,-1,-1,nIvar,nCut, nCType?AliRDHFDTNode::kGT:AliRDHFDTNode::kLT, AliRDHFDTNode::kNull, nPur, nRes);
   // For root node
   if(nPos=="s") node->SetNodeType(AliRDHFDTNode::kRoot);
   else if(!nType) node->SetNodeType(AliRDHFDTNode::kMed);
   else if(nPos=="l") node->SetNodeType(AliRDHFDTNode::kSignal);
   else if(nPos=="r") node->SetNodeType(AliRDHFDTNode::kBkg);
   
   return node;   
}
