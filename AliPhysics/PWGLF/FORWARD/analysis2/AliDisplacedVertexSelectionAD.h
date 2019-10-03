#ifndef ALIDISPLACEDVERTEXSELECTIONAD_H
#define ALIDISPLACEDVERTEXSELECTIONAD_H
#include <TObject.h>

class AliESDEvent;
class AliESDfriend;
class AliESDAD;
class AliESDADfriend;
class AliADCalibData;
class TList;
class TH2;
class TH1;

/** 
 * Check AD signal to see of this is a satelitte-main collision, and
 * give interaction point.
 */
class AliDisplacedVertexSelectionAD : public TObject
{
public:
  enum EEventType {
    kUnknown,
    kMain,
    kSatelliteA,
    kSatelliteC
  };
  /** 
   * Some default values 
   */
  enum {
    kInvalidTime = -9999,
    kInvalidVtxZ = 9999
  };
  /** 
   * Constructor 
   */
  AliDisplacedVertexSelectionAD()
    : fIPz(kInvalidVtxZ),
      fEventType(kUnknown),
      fAC(0),
      fSumDelta(0),
      fSumDeltaSatA(0),
      fSumDeltaSatC(0),
      fIPzAll(0),
      fIPzMain(0),
      fIPzSatA(0),
      fIPzSatC(0),
      fIPzDelta(0),
      fIPzSum(0),
      fIPzBunch(0),
      fSlewing(0),
      fCalib(0),
      fSpacing(2.5),
      fMaxBunch(15)
  {}
  /** 
   * Copy constructor 
   * 
   * @param o Object to copy from 
   */
  AliDisplacedVertexSelectionAD(const AliDisplacedVertexSelectionAD& o)
    : TObject(o),
      fIPz(kInvalidVtxZ),
      fEventType(kUnknown),
      fAC(0),
      fSumDelta(0),
      fSumDeltaSatA(0),
      fSumDeltaSatC(0),
      fIPzAll(0),
      fIPzMain(0),
      fIPzSatA(0),
      fIPzSatC(0),
      fIPzDelta(0),
      fIPzSum(0),
      fIPzBunch(0),
      fSlewing(0),
      fCalib(0),
      fSpacing(2.5),
      fMaxBunch(15)
  {}
  /** 
   * Destructor 
   */
  virtual ~AliDisplacedVertexSelectionAD() {}
  /** 
   * Assignement operator 
   * 
   * @param o Object to assign from 
   * 
   * @return Reference to this 
   */
  AliDisplacedVertexSelectionAD&
  operator=(const AliDisplacedVertexSelectionAD& o) { return *this; }
  /**
   * Check if this event is marked as a satellite-main interaction 
   *
   * @return true if the found vertex isn't invalid
   */
  Bool_t IsSatellite() const { return IsSatelliteA() || IsSatelliteC(); }
  /**
   * Check if this event is marked as a satellite-main interaction 
   *
   * @return true if the event is satellite-main from A-debunched beam
   */
  Bool_t IsSatelliteA() const { return fEventType == kSatelliteA; }
  /**
   * Check if this event is marked as a satellite-main interaction 
   *
   * @return true if the event is satellite-main from C-debunched beam
   */
  Bool_t IsSatelliteC() const { return fEventType == kSatelliteC; }
  /** 
   * Check if this is main-main event 
   * 
   * @return 
   */
  Bool_t IsMain() const { return fEventType == kMain; }
  /** 
   * Get the event type 
   * 
   * @return Event type 
   */
  EEventType GetEventType() const { return fEventType; }
  /** 
   * Get the interaction point Z-coordinate from ZDC timing. 
   * 
   * @return Interaction point Z-coordinate
   */
  Double_t GetIPz() const { return fIPz; }
  /** 
   * Define the output 
   * 
   * @param l     List to add output to
   * @param name  Name of the list 
   * @param mc    True if we're looking at MC data
   */
  void SetupForData(TList* l, const char* name=0, Bool_t mc=false);
  /** 
   * Print information 
   * 
   * @param option  Not used 
   */
  void Print(Option_t* option="") const;
  /** 
   * Process an ESD event to get the information 
   * 
   * @param esd ESD event 
   * 
   * @return true on success
   */
  Bool_t Process(const AliESDEvent* esd);
protected:
  /** 
   * Calculate the (weighted) mean time on either side 
   * 
   * @param adESD     AD object 
   * @param adFriend  AD friend object 
   * @param aNotC     A side if true, C otherwise 
   * 
   * @return The mean time 
   */
  Double_t MeanTime(AliESDAD*       adESD,
		    AliESDADfriend* adFriend,
		    Bool_t          aNotC) const;

  Double_t        fIPz;       //!
  EEventType      fEventType; // 
  TH2*            fAC;          //! 
  TH2*            fSumDelta;    //! 
  TH2*            fSumDeltaSatA;//! 
  TH2*            fSumDeltaSatC;//! 
  TH1*            fIPzAll;      //! 
  TH1*            fIPzMain;     //! 
  TH1*            fIPzSatA;     //! 
  TH1*            fIPzSatC;     //! 
  TH2*            fIPzDelta;    //! 
  TH2*            fIPzSum;      //!
  TH2*            fIPzBunch;    //!
  TList*          fSlewing;     //! 
  AliADCalibData* fCalib;       //!
  Double_t        fSpacing;
  Int_t           fMaxBunch;  // Maximum bunch number 
  ClassDef(AliDisplacedVertexSelectionAD,1);
};

#endif
// Local Variables:
//  mode: C++
// End:
