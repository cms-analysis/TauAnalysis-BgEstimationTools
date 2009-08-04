#include "TauAnalysis/BgEstimationTools/interface/BgEstCovMatrix.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TMath.h>

#include <iomanip>

BgEstCovMatrix::BgEstCovMatrix(unsigned d)
  : cov_(d, d), auxDelta_(d), auxMean_(d)
{
  numVar_ = d;

  iValue_ = 0;
}

void BgEstCovMatrix::update(const TVectorD& value)
{
  if ( value.GetNoElements() != (int)numVar_ ) {
    edm::LogError("BgEstCovMatrix::update") << "Given value has invalid dimension = " << value.GetNoElements() << "," 
					    << " expected = " << numVar_ << " --> mean value will NOT be updated !!";
    return;
  }

  if ( iValue_ == 0 ) {
    for ( unsigned iRow = 0; iRow < numVar_; ++iRow ) {
      for ( unsigned iColumn = 0; iColumn < numVar_; ++iColumn ) {
	cov_(iRow, iColumn) = 0.;
      }
    }

    for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
      auxMean_(iVar) = value(iVar);
    }
  } else {
    double weight = iValue_/(iValue_ + 1.);

    for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
      auxDelta_(iVar) = value(iVar) - auxMean_(iVar);
    }
    
    for ( unsigned iRow = 0; iRow < numVar_; ++iRow ) {
      for ( unsigned iColumn = 0; iColumn < numVar_; ++iColumn ) {
	cov_(iRow, iColumn) += weight*auxDelta_(iRow)*auxDelta_(iColumn);
      }
    }
    
    for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
      auxMean_(iVar) += auxDelta_(iVar)/(iValue_ + 1.);
    }
  }

  ++iValue_;
}

TMatrixD BgEstCovMatrix::operator()() const
{
//--- normalize covariance matrix
//    to number of values for which it has been computer

  TMatrixD cov_normalized(numVar_, numVar_);

  if ( iValue_ > 0 ) {
    for ( unsigned iRow = 0; iRow < numVar_; ++iRow ) {
      for ( unsigned iColumn = 0; iColumn < numVar_; ++iColumn ) {
	cov_normalized(iRow, iColumn) = cov_(iRow, iColumn)/iValue_;
      }
    }
  }
  
  return cov_normalized;
}

double BgEstCovMatrix::sigma(unsigned i) const
{
  if ( i >= 0 && i < numVar_ ) {
    return TMath::Sqrt(cov_(i, i));
  } else {
    edm::LogError("BgEstCovMatrix::sigma") << "Given index i = " << i << " out of bounds,"
					   << " expected range = 0.." << (numVar_ - 1) << " !!";
    return 0.;
  }
}

double BgEstCovMatrix::correlation(unsigned i, unsigned j) const
{
  if ( i >= 0 && i < numVar_ &&
       j >= 0 && j < numVar_ ) {
    double sigmaI = TMath::Sqrt(cov_(i, i));
    double sigmaJ = TMath::Sqrt(cov_(j, j));
    return cov_(i,j)/(sigmaI*sigmaJ);
  } else {
    edm::LogError("BgEstCovMatrix::correlation") << "Given indices i = " << i << ", j = " << j << " out of bounds,"
						 << " expected range = 0.." << (numVar_ - 1) << " !!";
    return 0.;
  }
}

void BgEstCovMatrix::print(std::ostream& outputStream, const std::vector<std::string>* varNames) const
{
  if ( varNames && varNames->size() != numVar_ ) {
    edm::LogError("BgEstCovMatrix::print") << "Given list of varable names of invalid size,"
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
    outputStream << varNames_temp[iVar] << ": mean = " << auxMean_(iVar) << ", sigma = " << sigma(iVar) << std::endl;
  }
  
  outputStream << "Covariance Matrix:" << std::endl;
  const unsigned widthColumn = 20;
  outputStream << std::setw(widthColumn) << std::left << " ";
  for ( unsigned iVar = 0; iVar < numVar_; ++iVar ) {
    outputStream << std::setw(widthColumn) << std::left << varNames_temp[iVar];
  }
  outputStream << std::endl;
  for ( unsigned iCharacter = 0; iCharacter < (widthColumn*(numVar_ + 1)); ++iCharacter ) {
    outputStream << "-";
  }
  outputStream << std::endl;
  TMatrixD cov_normalized = this->operator()();
  for ( unsigned iRow = 0; iRow < numVar_; ++iRow ) {
    outputStream << std::setw(widthColumn) << std::left << varNames_temp[iRow];
    for ( unsigned iColumn = 0; iColumn < numVar_; ++iColumn ) {
      outputStream << std::setw(4) << std::left << " ";
      outputStream << std::setw(widthColumn - 4) << std::setprecision(3) << std::fixed << std::left << cov_normalized(iRow, iColumn);
    }
    outputStream << std::endl;
  }
  for ( unsigned iCharacter = 0; iCharacter < (widthColumn*(numVar_ + 1)); ++iCharacter ) {
    outputStream << "-";
  }
  outputStream << std::endl;  
} 
