#ifndef TauAnalysis_BgEstimationTools_BinGrid_h  
#define TauAnalysis_BgEstimationTools_BinGrid_h

/** \class BinGrid
 *
 * Store number of bins and bin-boundaries
 * used for binning event passing selection
 * (in any number of dimensions)
 * 
 * \author Christian Veelken, UC Davis
 *         (inspired by code written for H1 by Paul Laycock, University of Liverpool)
 *
 * \version $Revision: 1.1 $
 *
 * $Id: Binning.h,v 1.1 2009/02/04 15:53:56 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <string>

class BinGrid 
{
 public: 
  explicit BinGrid(const edm::ParameterSet&);
  virtual ~BinGrid();
  
  unsigned dimensions() const { return numDimensions_; }
  unsigned numBins() const { return numBinsTotal_; }
  
  const std::vector<std::string>& objVarNames() const { return objVarNames_; }

  virtual int binNumber(const std::vector<double>&) const;

  virtual std::vector<double> binCenter(unsigned) const;
  virtual double binVolume(unsigned) const;

 private:
//--- auxiliary functions
//    (for encoding/decoding of bin numbers)
  typedef std::vector<unsigned> vunsigned;
  unsigned encodeTotBin(const vunsigned&) const;
  vunsigned decodeTotBin(unsigned) const;
  unsigned getDimValue(unsigned) const;

//--- data-members
  unsigned numDimensions_;

  typedef std::vector<std::string> vstring;
  vstring objVarNames_;

  // list of "left" bin-edges (similar to ROOT's TH1 histogram class) for each dimension
  // (xMin_i = binEdges[i][0],..binBoundary_j = binEdges[i][j + 1].., xMax_i = binEdges[i][numBinsPerDimension])
  typedef std::vector<double> vdouble;  
  std::vector<vdouble> binEdges_;

  // number of bins per dimension
  vunsigned numBinsPerDimension_; 
  unsigned numBinsTotal_;

  // auxiliary data-structure for encoding/decoding
  // of bin numbers
  vunsigned dimValues_;
};

#endif  


