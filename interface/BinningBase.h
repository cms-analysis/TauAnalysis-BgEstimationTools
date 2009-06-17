#ifndef TauAnalysis_BgEstimationTools_BinningBase_h  
#define TauAnalysis_BgEstimationTools_BinningBase_h

/** \class BinningBase
 *
 * Pure virtual base-class for storing number of events passing selection in different bins
 * 
 * \author Christian Veelken, UC Davis
 *         (inspired by code written for H1 by Paul Laycock, University of Liverpool)
 *
 * \version $Revision: 1.1 $
 *
 * $Id: BinningBase.h,v 1.1 2009/06/11 07:23:28 veelken Exp $
 *
 */

#include "TauAnalysis/BgEstimationTools/interface/BinGrid.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>

class BinningBase
{
 public:
  BinningBase();
  BinningBase(const edm::ParameterSet&);
  virtual ~BinningBase();

  const std::string& name() const { return name_; }

  const BinGrid* binGrid() const { return binGrid_; }

  virtual void bin(const std::vector<double>&, double = 1.) = 0;

  virtual void print(std::ostream&) const = 0;

  friend std::vector<std::string>& operator<<(std::vector<std::string>&, const BinningBase&);
  friend std::vector<std::string>& operator>>(std::vector<std::string>&, BinningBase&);

 protected:
  virtual std::vector<std::string> encodeStringRep() const;
  virtual void decodeStringRep(std::vector<std::string>&);

  std::string name_;

  const BinGrid* binGrid_;
};

std::vector<std::string>& operator<<(std::vector<std::string>&, const BinningBase&);
std::vector<std::string>& operator>>(std::vector<std::string>&, BinningBase&);

#endif  


