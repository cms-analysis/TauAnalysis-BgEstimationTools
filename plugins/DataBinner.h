#ifndef TauAnalysis_BgEstimationTools_DataBinner_h  
#define TauAnalysis_BgEstimationTools_DataBinner_h

/** \class DataBinner
 *
 * Extract observables to be binned from the event and fill them into DataBinning object,
 * using the bin-grid associated with the DataBinning object.
 *
 * Class can be used for Data and Monte Carlo, 
 * but is restricted to use reconstruction level information only
 * (cannot compute acceptance, purity, stability).
 * 
 * \author Christian Veelken, UC Davis
 *         (inspired by code written for H1 by Paul Laycock, University of Liverpool)
 *
 * \version $Revision: 1.1 $
 *
 * $Id: BinnerBase.h,v 1.1 2009/02/04 15:53:56 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/BinnerBase.h"

class DataBinner : public BinnerBase
{
 public: 
  explicit DataBinner(const edm::ParameterSet&);
  ~DataBinner();
};

#endif  


