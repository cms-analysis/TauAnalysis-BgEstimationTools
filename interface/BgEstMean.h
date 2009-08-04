#ifndef TauAnalysis_BgEstimationTools_BgEstMean_h  
#define TauAnalysis_BgEstimationTools_BgEstMean_h

/** \class BgEstMean
 *
 * Compute average of a list of values,
 * each value being a vector of dimension d
 *
 *    NOTE: algorithms for numerically stable computation of mean and covariance matrix 
 *          implemented in BgEstMean and BgEstCovMatrix classes 
 *          taken from:
 *            Ronald A. Thisted 
 *            "Elements of Statistical Computing: Numerical Computation", pp. 84-91
 *            Chapman & Hall, 1988 
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: BgEstMean.h,v 1.1 2009/06/11 07:23:29 veelken Exp $
 *
 */

#include <TVectorD.h>

class BgEstMean
{
 public:  
  explicit BgEstMean(unsigned d);
  ~BgEstMean() {}
  
  void update(const TVectorD& value);

  const TVectorD& operator()() const { return mean_; }
  double operator()(unsigned i) const { return mean_(i); }

 private:
  TVectorD mean_;
  unsigned numVar_;
  int iValue_;
  TVectorD auxDelta_;
};

#endif  


