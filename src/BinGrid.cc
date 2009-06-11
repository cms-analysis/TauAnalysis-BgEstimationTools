#include "TauAnalysis/BgEstimationTools/interface/BinGrid.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

BinGrid::BinGrid(const edm::ParameterSet& cfg)
{
  numDimensions_ = 0;
  
  typedef std::vector<edm::ParameterSet> vParameterSet;
  vParameterSet cfgBinning = cfg.getParameter<vParameterSet>("config");
  for ( vParameterSet::const_iterator cfg1dObjVar = cfgBinning.begin(); 
	cfg1dObjVar != cfgBinning.end(); ++cfg1dObjVar ) {
    edm::ParameterSet cfg1dBinning = cfg1dObjVar->getParameter<edm::ParameterSet>("binning");
    
    double xMin = cfg1dBinning.getParameter<double>("min");
    double xMax = cfg1dBinning.getParameter<double>("max");

    typedef std::vector<double> vdouble;
    vdouble binBoundaries = cfg1dBinning.getParameter<vdouble>("boundaries");

    unsigned numBinBoundaries = binBoundaries.size();
    unsigned numBins_i = numBinBoundaries + 2;
    vdouble binEdges_i(numBins_i);
    binEdges_i[0] = xMin;
    for ( unsigned iBin = 0; iBin < numBinBoundaries; ++iBin ) {
      binEdges_i[iBin + 1] = binBoundaries[iBin];
    }
    binEdges_i[numBins_i - 1] = xMax;

    binEdges_.push_back(binEdges_i);
    numBinsPerDimension_.push_back(numBins_i);
    ++numDimensions_;
  }

//--- initialize "value" of each dimension i, defined by the product:
//
//      numBins(i + 1) * numBins(i + 2) * .. * numBins(numDimensions - 1)
//
  dimValues_.resize(numDimensions_);
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    dimValues_[iDimension] = getDimValue(iDimension);
  }

//--- compute total number of bins in multi-dimensional grid
  numBinsTotal_ = 1;
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    numBinsTotal_ *= numBinsPerDimension_[iDimension];
  }
}

BinGrid::~BinGrid()
{
//--- nothing to be done yet...
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

unsigned BinGrid::encodeTotBin(const vunsigned& binIndices) const
{
  assert(binIndices.size() == numDimensions_);

  unsigned totBin = 0;
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    totBin += binIndices[iDimension] * dimValues_[iDimension];
  }

  return totBin;
}

std::vector<unsigned> BinGrid::decodeTotBin(unsigned totBin) const
{
  std::cout << " totBin = " << totBin << std::endl;

  assert(totBin >= 0 && totBin < numBinsTotal_);

  vunsigned binIndices(numDimensions_);

  unsigned totBin_undecoded = totBin;
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    unsigned dimValue_i = dimValues_[iDimension];
    unsigned binIndex = totBin_undecoded / dimValue_i;
    binIndices[iDimension] = binIndex;
    totBin_undecoded -= (dimValue_i * binIndex);
  }

  std::cout << "totBin_undecoded = " << totBin_undecoded << std::endl;

  assert(totBin_undecoded == 0);

  return binIndices;
}

unsigned BinGrid::getDimValue(unsigned i) const
{
  assert(i >= 0 && i < numDimensions_);

  unsigned dimValue = 1;
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    dimValue *= numBinsPerDimension_[iDimension];
  }

  return dimValue;
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

int BinGrid::binNumber(const std::vector<double>& x) const
{
//--- check that vector x given as function argument
//    is of the same dimensionality as the bin-grid
  if ( x.size() != numDimensions_ ) {
    edm::LogError ("BinGrid::getBin") << " Invalid dimensionality = " << x.size() << " of vector x given as function argument," 
				      << " expected dimensionality = " << numDimensions_ << " !!";
    return -1;
  }

//--- compute bin value;
//    the bin value is computed as sum over values k_i calculated separately for each dimension i as:
//
//      k_i = bin_i * numBins(i + 1) * numBins(i + 2) * .. * numBins(numDimensions - 1)
//
//    where bin_i is taken as the index of the "left" bin-edge for which 
//
//      binEdges[i,k] < x[i] < binEdges[i,k + 1]
//
  vunsigned binIndices(numDimensions_);
  for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
    const vdouble& binEdges_i = binEdges_[iDimension];

    bool xInRange = false;

    unsigned numBins_i = numBinsPerDimension_[iDimension];

    for ( unsigned iBin = 0; iBin < numBins_i; ++iBin ) {
      if ( x[iDimension] >= binEdges_i[iBin] && x[iDimension] < binEdges_i[iBin + 1] ) {
	binIndices[iDimension] = iBin;
	xInRange = true;
	break;
      }
    }

    if ( !xInRange ) {
      edm::LogWarning ("BinGrid::getBin") << " Value x[" << iDimension << "] = " << x[iDimension] << " outside Range"
					  << " [xMin, xMax[ = [" << binEdges_i[0] << ", " << binEdges_i[numBins_i] << "[ !!";
      return -1;
    }
  }

  return encodeTotBin(binIndices);
}
 
//
//-----------------------------------------------------------------------------------------------------------------------
//

std::vector<double> BinGrid::binCenter(unsigned totBin) const
{
  if ( totBin >= 0 && totBin < numBinsTotal_ ) {
    std::vector<double> binCenter(numDimensions_);

    vunsigned binIndices = decodeTotBin(totBin);
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      const vdouble& binEdges_i = binEdges_[iDimension];
      unsigned binIndex_i = binIndices[iDimension];
      binCenter[iDimension] = 0.5*(binEdges_i[binIndex_i] + binEdges_i[binIndex_i + 1]);
    } 

    return binCenter; 
  } else {
    edm::LogError ("BinGrid::binCenter") << " Invalid bin = " << totBin << " !!";
    return std::vector<double>(numDimensions_);
  }
}

double BinGrid::binVolume(unsigned totBin) const
{
  if ( totBin >= 0 && totBin < numBinsTotal_ ) {
    double binVolume = 1.;

    vunsigned binIndices = decodeTotBin(totBin);
    for ( unsigned iDimension = 0; iDimension < numDimensions_; ++iDimension ) {
      const vdouble& binEdges_i = binEdges_[iDimension];
      unsigned binIndex_i = binIndices[iDimension];
      double binWidth = binEdges_i[binIndex_i + 1] - binEdges_i[binIndex_i];
      binVolume *= binWidth;
    }

    return binVolume; 
  } else {
    edm::LogError ("BinGrid::binVolume") << " Invalid bin = " << totBin << " !!";
    return 0.;
  }
}

