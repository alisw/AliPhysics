#pragma link off all functions;
#pragma link off all globals;
#pragma link off all classes;

// Bertrand's new browser and helpers
#pragma link C++ class TGClassBrowser;
#pragma link C++ class TGCommandPlugin;
#pragma link C++ class TGFileBrowser;
#pragma link C++ class TGNewBrowser;


#pragma link C++ namespace Reve;
#pragma link C++ global   gReve; // In ReveManager ... should move.

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
#pragma link C++ class Reve::ReferenceBackPtr+;

#pragma link C++ function Reve::ColorFromIdx;
#pragma link C++ function Reve::FindColorVar;

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
#pragma link C++ class Reve::ZTransSubEditor+;
#pragma link C++ class Reve::ZTransEditor+;

// Stepper
#pragma link C++ class Reve::GridStepper+;
#pragma link C++ class Reve::GridStepperSubEditor+;
#pragma link C++ class Reve::GridStepperEditor+;

// RGBAPalette
#pragma link C++ class Reve::RGBAPalette+;
#pragma link C++ class Reve::RGBAPaletteEditor+;
#pragma link C++ class Reve::RGBAPaletteSubEditor+;

// Plexes
#pragma link C++ class Reve::VoidCPlex+;
#pragma link C++ class Reve::VoidCPlex::iterator-;

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
#pragma link C++ class Reve::RenderElementList+;
#pragma link C++ class Reve::RenderElementEditor+;

#pragma link C++ class std::list<Reve::RenderElement*>;
#pragma link C++ class std::list<Reve::RenderElement*>::iterator;
#pragma link C++ typedef Reve::RenderElement::List_t;
#pragma link C++ typedef Reve::RenderElement::List_i;

// GL-interface
#pragma link C++ class Reve::Scene+;
#pragma link C++ class Reve::SceneList+;
#pragma link C++ class Reve::SceneInfo+;
#pragma link C++ class Reve::Viewer+;
#pragma link C++ class Reve::ViewerList+;

// Pad
#pragma link C++ class Reve::Pad+;

// VSDSelector
#pragma link C++ class Reve::VSDSelector+;

// RGBrowser
#pragma link C++ class Reve::RGLTEFrame+;
#pragma link C++ class Reve::RGBrowser+;

// RGEditor
#pragma link C++ class Reve::RGEditor+;

// RMacro
#pragma link C++ class Reve::RMacro+;

// ReveManager
#pragma link C++ class Reve::ReveManager+;

// RGValuators
#pragma link C++ class Reve::RGValuatorBase+;
#pragma link C++ class Reve::RGValuator+;
#pragma link C++ class Reve::RGDoubleValuator+;
#pragma link C++ class Reve::RGTriVecValuator+;

// Non-linear transformations
#pragma link C++ class Reve::NLTProjectable+;
#pragma link C++ class Reve::NLTProjected+;
#pragma link C++ class Reve::NLTProjection+;
#pragma link C++ class Reve::RhoZ+;
#pragma link C++ class Reve::CircularFishEye+;

#pragma link C++ class Reve::NLTProjector+;
#pragma link C++ class Reve::NLTProjectorEditor+;
#pragma link C++ class Reve::NLTProjectorGL+;

#pragma link C++ class Reve::NLTTrack+;
#pragma link C++ class Reve::NLTTrackGL+;
#pragma link C++ class Reve::NLTTrackList+;

#pragma link C++ class Reve::NLTPolygon+;
#pragma link C++ class Reve::NLTPolygonSet+;
#pragma link C++ class Reve::NLTPolygonSetEditor+;
#pragma link C++ class Reve::NLTPolygonSetGL+;

#pragma link C++ class Reve::NLTPointSet+;
#pragma link C++ class Reve::NLTPointSetGL+;

//=====================================
// Graphical elements (with renderers)
//=====================================

// Track
#pragma link C++ class Reve::Track+;
#pragma link C++ class Reve::TrackGL+;
#pragma link C++ class Reve::TrackEditor+;
#pragma link C++ class Reve::TrackRnrStyle+;
#pragma link C++ class Reve::TrackList+;
#pragma link C++ class Reve::TrackListEditor+;
#pragma link C++ class Reve::TrackRnrStyleSubEditor+;
#pragma link C++ class Reve::TrackRnrStyleEditor+;
#pragma link C++ class Reve::TrackCounter+;
#pragma link C++ class Reve::TrackCounterEditor+;

// Cascade
#pragma link C++ class Reve::Cascade+;
#pragma link C++ class Reve::CascadeList+;
#pragma link C++ class Reve::CascadeListEditor+;

// V0
#pragma link C++ class Reve::V0+;
#pragma link C++ class Reve::V0List+;
#pragma link C++ class Reve::V0ListEditor+;

// PointSet
#pragma link C++ class Reve::PointSet+;
#pragma link C++ class Reve::PointSetArray+;
#pragma link C++ class Reve::PointSetArrayEditor+;

// Line
#pragma link C++ class Reve::Line+;
#pragma link C++ class Reve::LineEditor+;
#pragma link C++ class Reve::LineGL+;

// FrameBox
#pragma link C++ class Reve::FrameBox+;
#pragma link C++ class Reve::FrameBoxGL+;

// DigitSet
#pragma link C++ class Reve::DigitSet+;
#pragma link C++ class Reve::DigitSetEditor+;

// QuadSet
#pragma link C++ class Reve::Quad+;
#pragma link C++ class Reve::OldQuadSet+;
#pragma link C++ class Reve::OldQuadSetGL+;
#pragma link C++ class Reve::QuadSet+;
#pragma link C++ class Reve::QuadSetGL+;

// BoxSet
#pragma link C++ class Reve::BoxSet+;
#pragma link C++ class Reve::BoxSetGL+;

// GeoNode
#pragma link C++ class Reve::GeoRnrEl+;
#pragma link C++ class Reve::GeoNodeRnrEl+;
#pragma link C++ class Reve::GeoTopNodeRnrEl+;
#pragma link C++ class Reve::GeoNodeRnrElEditor+;
#pragma link C++ class Reve::GeoTopNodeRnrElEditor+;

#pragma link C++ class TGeoShapeExtract+;
#pragma link C++ class Reve::GeoShapeRnrEl+;

// TrianlgeSet
#pragma link C++ class Reve::TriangleSet+;
#pragma link C++ class Reve::TriangleSetEditor+;
#pragma link C++ class Reve::TriangleSetGL+;

// StraightLineSet
#pragma link C++ class Reve::StraightLineSet+;
#pragma link C++ class Reve::StraightLineSetGL+;
#pragma link C++ class Reve::StraightLineSetEditor+;
#pragma link C++ class Reve::NLTSLineSet+;
