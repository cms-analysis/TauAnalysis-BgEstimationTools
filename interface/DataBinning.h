#ifndef TauAnalysis_BgEstimationTools_DataBinning_h  
#define TauAnalysis_BgEstimationTools_DataBinning_h

/** \class DataBinning
 *
 * Store number of events passing selection in different bins of DataBinning object.
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
 * $Id: DataBinning.h,v 1.1 2009/02/04 15:53:56 veelken Exp $
 *
 */

#include "TauAnalysis/BgEstimationTools/interface/BinningBase.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

class DataBinning : public BinningBase
{
 public:
  explicit DataBinning(const edm::ParameterSet&);
  ~DataBinning();

  const std::string& name() const { return name_; }

  const BinGrid* binGrid() const { return binGrid_; }

  void bin(const std::vector<double>&, double = 1.);

  void print(std::ostream&) const;

 protected:
  typedef std::vector<double> vdouble;
  vdouble binContents_;
  vdouble binSumw2_;
  unsigned numBins_;
};

#endif  


