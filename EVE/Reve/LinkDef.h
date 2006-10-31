#pragma link off all functions;
#pragma link off all globals;
#pragma link off all classes;

#pragma link C++ namespace Reve;
#pragma link C++ global   gReve; // In RGTopFrame ... should move.

//================================
// base/
//================================

// Reve
#pragma link C++ function Reve::SetupEnvironment;

#pragma link C++ function Reve::CheckMacro;
#pragma link C++ function Reve::AssertMacro;
#pragma link C++ function Reve::Macro;
#pragma link C++ function Reve::LoadMacro;

#pragma link C++ function Reve::PushPad;
#pragma link C++ function Reve::PopPad;
#pragma link C++ class Reve::Exc_t+;
#pragma link C++ class Reve::PadHolder+;
#pragma link C++ class Reve::GeoManagerHolder+;
#pragma link C++ class Reve::ReferenceCount+;

#pragma link C++ function Reve::ColorFromIdx;
#pragma link C++ class Reve::RGBAPalette;

// PODs
#pragma link C++ class Reve::Vector+;
#pragma link C++ class Reve::PathMark+;
#pragma link C++ class Reve::MCTrack+;
#pragma link C++ class Reve::MCTrackRef+;
#pragma link C++ class Reve::Hit+;
#pragma link C++ class Reve::Cluster+;
#pragma link C++ class Reve::RecTrack+;
#pragma link C++ class Reve::RecKink+;
#pragma link C++ class Reve::RecV0+;
#pragma link C++ class Reve::GenInfo+;

// ZTrans
#pragma link C++ class Reve::ZTrans-;

// EventBase, VSDEvent, VSD
#pragma link C++ class Reve::EventBase+;
#pragma link C++ class Reve::EvTree+;
#pragma link C++ class Reve::VSDEvent+;
#pragma link C++ class Reve::VSDTree+;
#pragma link C++ class Reve::VSD+;

// TTreeTools
#pragma link C++ class TSelectorToEventList+;
#pragma link C++ class TTreeQuery+;
#pragma link C++ class TPointSelectorConsumer+;
#pragma link C++ class TPointSelector+;

// RenderElement
#pragma link C++ class Reve::RenderElement+;
#pragma link C++ class Reve::RenderElement::ListTreeInfo+;
#pragma link C++ class Reve::RenderElementObjPtr+;
#pragma link C++ class Reve::RenderElementListBase+;
#pragma link C++ class Reve::RenderElementList+;

// Pad
#pragma link C++ class Reve::Pad+;

// VSDSelector
#pragma link C++ class Reve::VSDSelector+;

// RGBrowser
#pragma link C++ class Reve::RGBrowser+;

// RGEditor
#pragma link C++ class Reve::RGEditor+;

// RMacro
#pragma link C++ class Reve::RMacro+;

// RGTopFrame
#pragma link C++ class Reve::RGTopFrame+;

// RGValuators
#pragma link C++ class Reve::RGValuatorBase+;
#pragma link C++ class Reve::RGValuator+;
#pragma link C++ class Reve::RGDoubleValuator+;
#pragma link C++ class Reve::RGTriVecValuator+;

//================================
// g3d/
//================================

// Track
#pragma link C++ class Reve::Track+;
#pragma link C++ class Reve::TrackList+;
#pragma link C++ class Reve::TrackRnrStyle+;

// PointSet
#pragma link C++ class Reve::PointSet+;
#pragma link C++ class Reve::PointSetArray+;

// QuadSet
#pragma link C++ class Reve::Quad+;
#pragma link C++ class Reve::OldQuadSet+;
#pragma link C++ class Reve::QuadSet+;

// BoxSet
#pragma link C++ class Reve::Box+;
#pragma link C++ class Reve::BoxSet+;

// GeoNode
#pragma link C++ class Reve::GeoNodeRnrEl+;
#pragma link C++ class Reve::GeoTopNodeRnrEl+;

// TrianlgeSet
#pragma link C++ class Reve::TriangleSet;

//================================
// ged/
//================================

#pragma link C++ class Reve::ZTransSubEditor+;
#pragma link C++ class Reve::ZTransEditor+;

#pragma link C++ class Reve::RenderElementEditor+;
#pragma link C++ class Reve::TrackListEditor+;
#pragma link C++ class Reve::GeoNodeRnrElEditor+;
#pragma link C++ class Reve::GeoTopNodeRnrElEditor+;

#pragma link C++ class Reve::PointSetArrayEditor+;

#pragma link C++ class Reve::TriangleSetEditor;

//================================
// gl/
//================================

// ReveGLRenderers
#pragma link C++ class Reve::OldQuadSetGL+;
#pragma link C++ class Reve::BoxSetGL+;
#pragma link C++ class Reve::TriangleSetGL;
