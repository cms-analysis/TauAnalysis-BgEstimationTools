#ifndef TauAnalysis_BgEstimationTools_BgEstMedian_h  
#define TauAnalysis_BgEstimationTools_BgEstMedian_h

/** \class BgEstMedian
 *
 * Compute median of a list of values;
 * in case of values being a vector of dimension d,
 * the median of the vector is computed as vector of medians for each component
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: BgEstMedian.h,v 1.2 2009/08/07 11:49:11 veelken Exp $
 *
 */

#include <TVectorD.h>

#include <iostream>
#include <vector>
#include <string>

class BgEstMedian
{
 public:  
  explicit BgEstMedian(unsigned d);
  ~BgEstMedian() {}
  
  void update(const TVectorD& value);

  const TVectorD& operator()() const;
  double operator()(unsigned i) const { return (*this)()(i); }

  void print(std::ostream&, const std::vector<std::string>* = 0) const;

 private:
  typedef std::vector<double> vComponent;
  mutable std::vector<vComponent> values_;
  unsigned numVar_;
  mutable bool isSorted_;
  mutable TVectorD auxMedian_; 
};

#endif  


