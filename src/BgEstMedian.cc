#include "TauAnalysis/BgEstimationTools/interface/BgEstMedian.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

BgEstMedian::BgEstMedian(unsigned d)
{
  numVar_ = d;
  isSorted_ = false;
}

void BgEstMedian::update(const TVectorD& value)
{
  if ( value.GetNoElements() != (int)numVar_ ) {
    edm::LogError("BgEstMedian::update") << "Given value has invalid dimension = " << value.GetNoElements() << "," 
					 << " expected = " << numVar_ << " --> mean value will NOT be updated !!";
    return;
  }

  for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
    values_[iVar].push_back(value(iVar));
  }

  isSorted_ = false;
}

const TVectorD& BgEstMedian::operator()() const
{
  if ( isSorted_ ) return auxMedian_;

  unsigned numValues = values_[0].size();
  if ( !numValues ) {
    edm::LogError("BgEstMedian::operator()") << " No values added yet; median cannot be computed !!";
    return auxMedian_;
  }

  for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
    vComponent& values_i = values_[iVar];
    assert(values_i.size() == numValues);

    std::sort(values_i.begin(), values_i.end());

    auxMedian_(iVar) = ( (numValues % 2) == 1 )  ?
      values_i[numValues / 2] : 0.5*(values_i[(numValues / 2) - 1] + values_i[numValues / 2]);
  }

  isSorted_ = true;

  return auxMedian_;
}

void BgEstMedian::print(std::ostream& outputStream, const std::vector<std::string>* varNames) const
{
  if ( varNames && varNames->size() != numVar_ ) {
    edm::LogError("BgEstMedian ::print") << "Given list of varable names of invalid size,"
					 << " expected lenght = " << numVar_ << " !!";
    return;
  }

  std::vector<std::string> varNames_temp;
  for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
    if ( varNames ) {
      varNames_temp.push_back(varNames->at(iVar));
    } else {
      std::ostringstream varName;
      varName << "var_" << iVar;
      varNames_temp.push_back(varName.str());
    }
  }

  for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
    outputStream << varNames_temp[iVar] << ": median = " << (*this)(iVar) << std::endl;
  }

  outputStream << std::endl;  
} 



