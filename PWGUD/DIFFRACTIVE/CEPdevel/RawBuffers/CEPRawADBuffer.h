#ifndef CEPRawADBuffer_H
#define CEPRawADBuffer_H

#include "TObject.h"
#include "AliESDAD.h"

class CEPRawADBuffer : public TObject 
{

  private:
    /// Number of AD cells
    UInt_t              fNCells;        // this is a fixed size (therefore hardcoded in constructor)
    /// Multiplicity in cell
    Float_t             fMult[16];
    /// ADC charge in cell 
    Float_t             fAdc[16];
    /// Leading time measured by TDC
    Float_t             fTime[16];
    /// Final decision for AD-A based on average time of channel 
    Int_t               fDecisionADA;
    /// Final decision for AD-C based on average time of channel
    Int_t               fDecisionADC;


  public:
    /// Constructor
                        CEPRawADBuffer();
    /// Destructor
    virtual             ~CEPRawADBuffer() {};

    // Modifiers
    /// Reset member variables
    void                Reset();

    // Setter functions
    /*! \brief Set the multiplicity in a certain AD cell. 
     *
     * The setters are currently straigt forward, no 
     * out_of_range checks are currently implemented.
     * Get the number cells by GetNCells().
     * @param mult the multiplicity in the cell i
     * @param i the cell index (0-15)
     */
    void                SetADMultiplicity  (Float_t mult,   UInt_t i)  { fMult[i] = mult;   }
    /*! \brief Set the charge in a certain AD cell. 
     *
     * @param charge the ADC charge in the cell i
     * @param i the cell index (0-15)
     */    
    void                SetADCharge        (Float_t adc,    UInt_t i)  { fAdc[i]  = adc;    }
    /*! \brief Set the the leading time measured by the TPC.
     * @param ADtime the time in the cell i
     * @param the cell index (0-15)
     */
    void                SetADTime          (Float_t ADtime, UInt_t i)  { fTime[i] = ADtime; }
    /*! \brief Set the final decision in AD-A.
     *
     * It is based on average time of channels.
     *
     * The key is set in AliVZERO:
     *      - kV0Invalid = -1
     *      - kV0EMpty = 0
     *      - kV0BB = 1
     *      - kV0BG = 2
     *      - kV0Fake = 3
     *
     * @param decisionA the decision for AD-A
     */
    void                SetADDecision_A    (Int_t decisionA)  { fDecisionADA = decisionA; }
    /*! \brief Set the final decision in AD-C.
     *
     * @see SetADDecision_A()
     * @param decisionC the decision for AD-C
     */
    void                SetADDecision_C    (Int_t decisionC)  { fDecisionADC = decisionC; }

    // Global Setter
    /*! \brief Global variable setter for AD.
     *
     * Calls all setter functions and fills the object.
     * @param ADobj the AliESDAD object which we get via GetADData()
     */
    void                SetADVariables(AliESDAD* ADobj);

    // Accessors
    /*! \brief Get the multiplicity in a single AD cell.
     *
     * This functions returns -999 if the integer 
     * provided goes above its boundary of 16 cells.
     * The number of cells is retrieved by GetNCells().
     * @see SetADMultiplicity()
     * @param i integer pointing to the cell.
     * @return The multiplicity in the cell i.
     */
    Float_t             GetADMultiplicity(Int_t i)  const;
     /*! \brief Get the ADC charge in a single AD cell.
      *
      * @see SetADCharge()
      * @param i integer pointing to the cell.
      * @return The ADC charge in the cell i.
      */
    Float_t             GetADCharge(Int_t i)        const;
     /*! \brief Get the time in a single AD cell.
      *
      * @see SetADTime()
      * @param i integer pointing to the cell.
      * @return The leading time measured by TDC in the cell i.
      */
    Float_t             GetADTime(Int_t i)          const;
     /*! \brief Get the total multiplicity of the AD (=sum of each cell)
      *
      * @return Sum of multiplicities in each cell
      */
    Float_t             GetADTotalMultiplicity()    const;
     /*! \brief Get the total time of the AD
      *
      * @return Sum of time in each cell
      */
    Float_t             GetADTotalTime()            const;
     /*! \brief Get the total charge of the AD
      *
      * @return Sum of charge in each cell
      */
    Float_t             GetADTotalCharge()          const;     
    /*! \brief Get the decision of AD-A.
      * 
      * @see SetADDecision_A()
      * @return The decision in AD-A
      */
    Int_t               GetADADecision_A()          const { return fDecisionADA; }
     /*! \brief Get the decision of AD-C.
      *
      * @see SetADDecision_C()
      * @return The decision in AD-C
      */
    Int_t               GetADADecision_C()          const { return fDecisionADC; }

    /// Return number of cells
    UInt_t              GetNCells()                 const { return fNCells; }


    ClassDef(CEPRawADBuffer, 1)     // CEP raw track buffer
};

#endif
